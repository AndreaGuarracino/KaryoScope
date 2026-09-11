"""Tests for the streaming ``bgzip`` writer in :mod:`karyoscope.core.io.bgzip`."""

from __future__ import annotations

import gzip
import shutil
import stat
from pathlib import Path

import pytest

import karyoscope.core.io.bgzip as bz
from karyoscope.core.external import ExternalToolError, ToolNotFoundError
from karyoscope.core.io.bgzip import bgzip_stage, open_bgzip_writer

requires_bgzip = pytest.mark.skipif(shutil.which("bgzip") is None, reason="bgzip not on PATH")


@requires_bgzip
def test_writer_round_trips_text_as_a_gzip_stream(tmp_path: Path) -> None:
    out = tmp_path / "x.bed.gz"
    with open_bgzip_writer(out, threads=2) as h:
        h.write("chr1\t0\t10\tA\n")
        h.write("chr1\t10\t20\tB\n")
    with gzip.open(out, "rt") as fh:
        assert fh.read() == "chr1\t0\t10\tA\nchr1\t10\t20\tB\n"
    # bgzip, not plain gzip: the BGZF extra subfield is in the first header
    assert out.read_bytes()[12:16] == b"BC\x02\x00"


@requires_bgzip
def test_writer_handles_an_empty_output(tmp_path: Path) -> None:
    out = tmp_path / "empty.gz"
    with open_bgzip_writer(out):
        pass
    with gzip.open(out, "rt") as fh:
        assert fh.read() == ""


def test_stage_adds_threads_only_above_one(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(bz, "require_tool", lambda name, **kw: name)
    assert bgzip_stage() == "bgzip"
    assert bgzip_stage(1) == "bgzip"
    assert bgzip_stage(6) == "bgzip -@ 6"


def test_missing_bgzip_is_a_tool_error(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    monkeypatch.setattr("karyoscope.core.external.shutil.which", lambda *_: None)
    with pytest.raises(ToolNotFoundError):
        open_bgzip_writer(tmp_path / "x.gz")


def test_failed_compressor_surfaces_at_close(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """A failing bgzip (full disk, killed) must not leave a silently bad file."""
    fake = tmp_path / "bgzip"
    fake.write_text("#!/bin/sh\ncat >/dev/null\necho 'disk full' >&2\nexit 3\n")
    fake.chmod(fake.stat().st_mode | stat.S_IEXEC)
    monkeypatch.setattr(bz, "require_tool", lambda name, **kw: str(fake))
    h = open_bgzip_writer(tmp_path / "out.gz")
    h.write("some text\n")
    with pytest.raises(ExternalToolError) as excinfo:
        h.close()
    assert "disk full" in str(excinfo.value)
