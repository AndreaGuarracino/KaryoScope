"""Lifecycle of the KMC query-names tee: nothing is left behind on any exit."""

from __future__ import annotations

import os
import subprocess
from pathlib import Path
from typing import ClassVar

import pytest

import karyoscope.core.io.kmc as kmc
from karyoscope.core.io.kmc import run_get_featureids


class _FakeProc:
    """Stands in for samtools / tee / the names consumer."""

    instances: ClassVar[list[_FakeProc]] = []

    def __init__(self, args, **kwargs):
        self.args = args
        self.killed = False
        self.waited = False
        self.returncode: int | None = None
        r, w = os.pipe()
        self.stdout = os.fdopen(r, "rb") if kwargs.get("stdout") is subprocess.PIPE else None
        self._w = os.fdopen(w, "wb")
        self.stdin = None
        self.stderr = None
        _FakeProc.instances.append(self)

    def poll(self):
        return self.returncode

    def kill(self):
        self.killed = True
        self.returncode = -9

    def wait(self, timeout=None):
        self.waited = True
        if self.returncode is None:
            self.returncode = 0
        self._w.close()
        return self.returncode

    def communicate(self):
        self.wait()
        return b"", b""


@pytest.fixture
def fake_children(monkeypatch: pytest.MonkeyPatch):
    _FakeProc.instances = []
    monkeypatch.setattr(kmc.subprocess, "Popen", _FakeProc)
    monkeypatch.setattr(kmc, "_kill_process_group", lambda proc: proc.kill())
    monkeypatch.setattr(kmc, "require_tool", lambda name, **kw: name)
    # The names sink resolves bgzip; the unit CI job has no htslib, and these
    # tests are about lifecycle, not compression.
    monkeypatch.setattr("karyoscope.core.io.bgzip.require_tool", lambda name, **kw: name)
    monkeypatch.setattr(kmc, "get_featureids_binary", lambda: "get_featureIDs")
    return _FakeProc.instances


def test_launch_failure_reaps_children_and_removes_fifo_and_sidecar(
    tmp_path: Path, fake_children, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Reproduces the review case: get_featureIDs cannot be executed. The
    consumer, samtools and tee were already started; all must be killed and
    reaped, the ks_names_*/decode.fifo directory removed, and no sidecar left."""

    def _permission_denied(*a, **kw):
        raise PermissionError(13, "Permission denied", "get_featureIDs")

    monkeypatch.setattr(kmc.subprocess, "run", _permission_denied)

    bam = tmp_path / "aln.bam"
    bam.write_bytes(b"")
    out = tmp_path / "out"
    out.mkdir()
    sidecar = out / "aln.query_names.txt.gz"
    sidecar.write_bytes(b"partial")

    with pytest.raises(PermissionError):
        run_get_featureids(
            db_path=tmp_path / "db",
            input_path=bam,
            output_dir=out,
            threads=1,
            prefix="aln.db",
            query_names_sidecar=sidecar,
        )

    assert len(fake_children) == 3, "consumer, samtools, tee were all launched"
    assert all(p.killed and p.waited for p in fake_children)
    assert not list(out.glob("ks_names_*")), "FIFO directory removed"
    assert not sidecar.exists(), "partial sidecar removed"


def test_interrupt_during_the_query_cleans_up_too(
    tmp_path: Path, fake_children, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(
        kmc.subprocess, "run", lambda *a, **kw: (_ for _ in ()).throw(KeyboardInterrupt())
    )
    bam = tmp_path / "aln.bam"
    bam.write_bytes(b"")
    out = tmp_path / "out"
    out.mkdir()
    with pytest.raises(KeyboardInterrupt):
        run_get_featureids(
            db_path=tmp_path / "db", input_path=bam, output_dir=out, threads=1, prefix="p",
            query_names_sidecar=out / "n.gz",
        )  # fmt: skip
    assert all(p.killed and p.waited for p in fake_children)
    assert not list(out.glob("ks_names_*"))


# --- real processes -------------------------------------------------------


def _group_is_empty(pgid: int, timeout: float = 3.0) -> bool:
    """True once no process in ``pgid`` remains (killed children are reaped by
    init within milliseconds; allow a little slack)."""
    import time

    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        try:
            os.killpg(pgid, 0)
        except ProcessLookupError:
            return True
        time.sleep(0.05)
    return False


def test_consumer_pipeline_children_die_with_it(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Real processes. The names consumer is `bash -c "a < fifo | b | c"`; if
    startup fails before tee ever opens the FIFO, its stages sit blocked on
    the FIFO. Killing only the bash would orphan them (three survivors were
    reproduced in review). They must all be gone when cleanup returns."""
    marker = f"KS_TEE_TEST_{os.getpid()}"

    def fake_sink(sidecar, threads=1, *, stdin_from=None):
        # Three real stages, like sed | sed | bgzip, none needing htslib. The
        # middle one carries a unique marker in its argv so a survivor can be
        # found by name whatever process group it ended up in.
        return f"cat < {stdin_from} | sed -n /{marker}/p | cat > {sidecar}"

    monkeypatch.setattr(kmc, "names_sink", fake_sink)
    fake_samtools = tmp_path / "samtools"
    fake_samtools.write_text("#!/bin/sh\nsleep 60\n")
    fake_samtools.chmod(0o755)

    def fake_require(name, **kw):
        if name == "samtools":
            return str(fake_samtools)
        if name == "tee":
            # Startup fails BEFORE tee opens the FIFO for writing -- but only
            # once the consumer has forked its stages, or there is nothing to
            # orphan and the test proves nothing (bash killed before its first
            # fork leaves no children). Wait for the marked stage to appear.
            assert _procs_matching(marker, settle=0.0, wait_for=True), "consumer stages up"
            raise PermissionError(13, "Permission denied", "tee")
        return name

    monkeypatch.setattr(kmc, "require_tool", fake_require)
    monkeypatch.setattr(kmc, "get_featureids_binary", lambda: "get_featureIDs")

    real_popen = subprocess.Popen
    leaders: list[int] = []

    def spy_popen(args, **kwargs):
        proc = real_popen(args, **kwargs)
        if args[0] == "bash":
            leaders.append(os.getpgid(proc.pid))
        return proc

    monkeypatch.setattr(kmc.subprocess, "Popen", spy_popen)

    bam = tmp_path / "aln.bam"
    bam.write_bytes(b"")
    out = tmp_path / "out"
    out.mkdir()
    try:
        with pytest.raises(PermissionError):
            run_get_featureids(
                db_path=tmp_path / "db", input_path=bam, output_dir=out, threads=1, prefix="p",
                query_names_sidecar=out / "n.gz",
            )  # fmt: skip

        survivors = _procs_matching(marker)
        assert survivors == [], f"pipeline stages orphaned by cleanup: {survivors}"
        assert len(leaders) == 1, "the consumer was started"
        assert leaders[0] != os.getpgid(0), "consumer ran in its own process group"
        assert _group_is_empty(leaders[0]), "consumer's whole process group is gone"
        assert not list(out.glob("ks_names_*"))
        assert not (out / "n.gz").exists()
    finally:
        # Never leave orphans behind on a failing run of this test.
        subprocess.run(["pkill", "-f", marker], check=False)
        subprocess.run(["pkill", "-f", str(fake_samtools)], check=False)


def _procs_matching(marker: str, settle: float = 0.5, *, wait_for: bool = False) -> list[str]:
    """Command lines of live processes whose argv contains ``marker``.

    ``wait_for`` polls up to a few seconds for at least one match (used to
    let the consumer fork its stages); otherwise ``settle`` gives killed
    children time to be reaped before the look.
    """
    import time

    def look() -> list[str]:
        r = subprocess.run(["pgrep", "-fl", marker], capture_output=True, text=True, check=False)
        return [line for line in r.stdout.splitlines() if line.strip()]

    if wait_for:
        deadline = time.monotonic() + 5.0
        while time.monotonic() < deadline:
            found = look()
            if found:
                return found
            time.sleep(0.02)
        return []
    time.sleep(settle)
    return look()
