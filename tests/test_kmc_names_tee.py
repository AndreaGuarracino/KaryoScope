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
    monkeypatch.setattr(kmc, "require_tool", lambda name, **kw: name)
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
