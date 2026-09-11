"""Tests for :mod:`karyoscope.core.io.query_names` -- the rank -> name sidecar."""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

import karyoscope.core.io.query_names as qn
from karyoscope.core.external import ExternalToolError
from karyoscope.core.io.query_names import (
    QUERY_NAMES_SUFFIX,
    _scan_pipeline,
    alignment_decode_cmd,
    names_sink,
    write_query_names_sidecar,
)
from karyoscope.exceptions import KaryoscopeError


def _names(sidecar: Path) -> list[str]:
    with gzip.open(sidecar, "rt") as fh:
        return fh.read().splitlines()


# --- real scan passes (awk + gzip) ---------------------------------------


def test_fastq_sidecar_lists_names_in_record_order(tmp_path: Path) -> None:
    """One name per record, header up to the first whitespace, rank order.

    The third record's quality line is all ``@`` (Phred 31), which is exactly
    why headers are picked by position and not by their leading character.
    """
    fq = tmp_path / "reads.fastq"
    fq.write_text(
        "@SRR123.1 1 length=8\nACGTACGT\n+\nIIIIIIII\n"
        "@SRR123.2 1 length=8\nTTTTAAAA\n+\nIIIIIIII\n"
        "@SRR123.3 1 length=8\nGGGGCCCC\n+\n@@@@@@@@\n"
    )
    sidecar = tmp_path / ("reads" + QUERY_NAMES_SUFFIX)
    assert write_query_names_sidecar(fq, sidecar) == sidecar
    assert _names(sidecar) == ["SRR123.1", "SRR123.2", "SRR123.3"]


def test_gzipped_fastq_and_multiline_fasta_are_handled(tmp_path: Path) -> None:
    fq_gz = tmp_path / "reads.fq.gz"
    with gzip.open(fq_gz, "wt") as fh:
        fh.write("@r1/1\nACGT\n+\nIIII\n@r1/2\nTTGG\n+\nIIII\n")
    fa = tmp_path / "asm.fa"
    fa.write_text(">contig_1 len=8\nACGT\nACGT\n>contig_2\nTT\n")

    s1 = tmp_path / "a.txt.gz"
    s2 = tmp_path / "b.txt.gz"
    write_query_names_sidecar(fq_gz, s1)
    write_query_names_sidecar(fa, s2)
    assert _names(s1) == ["r1/1", "r1/2"]
    assert _names(s2) == ["contig_1", "contig_2"]


def test_empty_input_yields_an_empty_sidecar_not_a_failure(tmp_path: Path) -> None:
    """awk rather than grep: zero records is empty output, not exit code 1."""
    fa = tmp_path / "empty.fa"
    fa.write_text("")
    sidecar = tmp_path / "empty.txt.gz"
    write_query_names_sidecar(fa, sidecar)
    assert _names(sidecar) == []


# --- pipeline shape ------------------------------------------------------


def test_scan_reads_a_plain_file_with_awk_alone() -> None:
    """No parser, no cat: awk opens the file itself."""
    assert _scan_pipeline(Path("asm.fa"), Path("n.gz")) == (
        "awk '/^>/ { print substr($1, 2) }' asm.fa | gzip > n.gz"
    )
    assert _scan_pipeline(Path("r.fastq"), Path("n.gz")) == (
        "awk 'NR % 4 == 1 { print substr($1, 2) }' r.fastq | gzip > n.gz"
    )


def test_scan_decompresses_gz_in_front_of_awk() -> None:
    assert _scan_pipeline(Path("r.fq.gz"), Path("n.gz")).startswith(
        "gzip -dc r.fq.gz | awk 'NR % 4 == 1"
    )
    assert _scan_pipeline(Path("asm.fasta.gz"), Path("n.gz")).startswith(
        "gzip -dc asm.fasta.gz | awk '/^>/"
    )


def test_scan_quotes_paths() -> None:
    pipeline = _scan_pipeline(Path("my reads.fq"), Path("out dir/n.gz"))
    assert "'my reads.fq'" in pipeline and "'out dir/n.gz'" in pipeline


def test_alignment_decode_matches_the_query_decode(monkeypatch: pytest.MonkeyPatch) -> None:
    """Names must line up with the ranks the query saw, so the samtools flags
    are the decode's flags: primary only, mate suffix on, the reference."""
    monkeypatch.setattr(qn, "require_tool", lambda name, **kw: f"/bin/{name}")
    cmd = alignment_decode_cmd(Path("tumor.cram"), reference=Path("ref.fa"), threads=4)
    assert cmd == [
        "/bin/samtools", "fasta", "-F", "0x900", "-N",
        "--reference", "ref.fa", "-@", "4", "tumor.cram",
    ]  # fmt: skip
    assert alignment_decode_cmd(Path("aln.bam")) == [
        "/bin/samtools", "fasta", "-F", "0x900", "-N", "aln.bam",
    ]  # fmt: skip


def test_cram_without_reference_is_refused(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(qn, "require_tool", lambda name, **kw: f"/bin/{name}")
    with pytest.raises(KaryoscopeError, match="--reference"):
        alignment_decode_cmd(Path("tumor.cram"))


def test_sink_is_awk_then_gzip_into_the_sidecar(tmp_path: Path) -> None:
    sidecar = tmp_path / "with space.txt.gz"
    sink = names_sink(sidecar)
    assert sink.startswith("awk '/^>/")
    assert sink.endswith(f"| gzip > '{sidecar}'")


def test_write_runs_one_linear_pipefail_pipeline(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    captured: dict[str, object] = {}

    class _Result:
        returncode = 0
        stderr = b""

    def _fake_run(cmd, **kwargs):
        captured["cmd"] = cmd
        return _Result()

    monkeypatch.setattr(qn.subprocess, "run", _fake_run)

    sidecar = tmp_path / "sub" / "names.txt.gz"
    write_query_names_sidecar(Path("reads.fq.gz"), sidecar, threads=8)
    cmd = captured["cmd"]
    assert cmd[:3] == ["bash", "-o", "pipefail"]
    pipeline = cmd[4]
    assert pipeline.startswith("gzip -dc reads.fq.gz | awk")
    assert pipeline.endswith(f"| gzip > {sidecar}")
    assert ">(" not in pipeline, "linear pipeline, no process substitution"
    assert sidecar.parent.is_dir(), "parent directory created up front"


def test_write_uses_samtools_for_an_alignment(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    captured: dict[str, object] = {}

    class _Result:
        returncode = 0
        stderr = b""

    monkeypatch.setattr(qn, "require_tool", lambda name, **kw: f"/bin/{name}")
    monkeypatch.setattr(
        qn.subprocess, "run", lambda cmd, **kw: captured.update(cmd=cmd) or _Result()
    )
    write_query_names_sidecar(Path("aln.bam"), tmp_path / "n.gz", threads=2)
    assert captured["cmd"][4].startswith("/bin/samtools fasta -F 0x900 -N -@ 2 aln.bam | awk")


def test_write_reports_a_failed_stage(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    class _Result:
        returncode = 1
        stderr = b"awk: boom"

    monkeypatch.setattr(qn.subprocess, "run", lambda *a, **kw: _Result())
    with pytest.raises(ExternalToolError) as excinfo:
        write_query_names_sidecar(Path("reads.fq"), tmp_path / "n.txt.gz", capture=True)
    assert "boom" in str(excinfo.value)
