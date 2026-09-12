"""Tests for :mod:`karyoscope.core.io.query_names` -- the rank -> name sidecar."""

from __future__ import annotations

import gzip
import shutil
from pathlib import Path

import pytest

import karyoscope.core.io.query_names as qn
from karyoscope.core.external import ExternalToolError
from karyoscope.core.io.query_names import (
    _NORMALISE,
    QUERY_NAMES_SUFFIX,
    _scan_pipeline,
    alignment_decode_cmd,
    names_sink,
    sniff_fastx,
    write_query_names_sidecar,
)
from karyoscope.exceptions import KaryoscopeError

requires_bgzip = pytest.mark.skipif(shutil.which("bgzip") is None, reason="bgzip not on PATH")


@pytest.fixture
def plain_bgzip(monkeypatch: pytest.MonkeyPatch) -> None:
    """Make the compressor stage read ``bgzip``, not an absolute path."""
    monkeypatch.setattr("karyoscope.core.io.bgzip.require_tool", lambda name, **kw: name)


def _names(sidecar: Path) -> list[str]:
    with gzip.open(sidecar, "rt") as fh:
        return fh.read().splitlines()


# --- real scan passes (sed/grep + cut + bgzip) ---------------------------


@requires_bgzip
def test_fastq_sidecar_lists_names_in_record_order(tmp_path: Path) -> None:
    """One name per record, header up to the first whitespace, rank order.

    The third record's quality line is all ``@`` (Phred 31), which is exactly
    why headers are picked by position and not by their leading character.
    """
    fq = tmp_path / "reads.fastq"
    fq.write_text(
        "@SRR123.1 1 length=8\nACGTACGT\n+\nIIIIIIII\n"
        "@SRR123.2\tBC:Z:ACGT 1\nTTTTAAAA\n+\nIIIIIIII\n"
        "@SRR123.3 1 length=8\nGGGGCCCC\n+\n@@@@@@@@\n"
    )
    sidecar = tmp_path / ("reads" + QUERY_NAMES_SUFFIX)
    assert write_query_names_sidecar(fq, sidecar) == sidecar
    assert _names(sidecar) == ["SRR123.1", "SRR123.2", "SRR123.3"]


@requires_bgzip
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


@requires_bgzip
def test_names_are_normalised_exactly_as_hks_normalises_them(tmp_path: Path) -> None:
    """hks: drop the first byte, then split_whitespace().next(). So leading
    whitespace is skipped, and CR, VT and FF end the name like space and tab."""
    fq = tmp_path / "odd.fastq"
    fq.write_bytes(
        b"@ read1 desc\nACGT\n+\nIIII\n"  # leading space: name is read1, not ""
        b"@read2\r\nACGT\n+\nIIII\n"  # CRLF: no trailing \r
        b"@read3\x0bdesc\nACGT\n+\nIIII\n"  # vertical tab separates
        b"@read4\tt1 t2\nACGT\n+\nIIII\n"  # tab separates
        b"@read5\x0cdesc\nACGT\n+\nIIII\n"  # form feed separates
    )
    sidecar = tmp_path / "odd.txt.gz"
    write_query_names_sidecar(fq, sidecar)
    assert _names(sidecar) == ["read1", "read2", "read3", "read4", "read5"]


@requires_bgzip
def test_format_and_compression_are_sniffed_not_inferred_from_the_name(tmp_path: Path) -> None:
    """hks reads the bytes, so an extensionless FASTQ or a gzip without .gz must
    get the same sidecar it would get under a conventional name."""
    plain_fastq_no_ext = tmp_path / "reads"
    plain_fastq_no_ext.write_text("@a 1\nACGT\n+\nIIII\n@b\nAC\n+\nII\n")
    gz_fasta_odd_ext = tmp_path / "asm.fq"  # says FASTQ, is gzipped FASTA
    with gzip.open(gz_fasta_odd_ext, "wt") as fh:
        fh.write(">c1 x\nACGT\n>c2\nAC\n")
    s1, s2 = tmp_path / "1.gz", tmp_path / "2.gz"
    write_query_names_sidecar(plain_fastq_no_ext, s1)
    write_query_names_sidecar(gz_fasta_odd_ext, s2)
    assert _names(s1) == ["a", "b"]
    assert _names(s2) == ["c1", "c2"]


def test_sniff_refuses_a_file_that_is_neither(tmp_path: Path) -> None:
    bad = tmp_path / "x.fa"
    bad.write_text("chr1\t0\t10\n")
    with pytest.raises(KaryoscopeError, match="not FASTA or FASTQ"):
        sniff_fastx(bad)
    empty = tmp_path / "e.fa"
    empty.write_bytes(b"")
    assert sniff_fastx(empty) == (False, False)


@requires_bgzip
def test_empty_input_yields_an_empty_sidecar_not_a_failure(tmp_path: Path) -> None:
    """Zero records is empty output under pipefail, not a failure."""
    fa = tmp_path / "empty.fa"
    fa.write_text("")
    sidecar = tmp_path / "empty.txt.gz"
    write_query_names_sidecar(fa, sidecar)
    assert _names(sidecar) == []


# --- pipeline shape ------------------------------------------------------

_FQ_SELECT = f"sed -n 'p;n;n;n' | sed '{_NORMALISE}'"
_FA_SELECT = f"sed -n '/^>/p' | sed '{_NORMALISE}'"


def test_scan_reads_a_plain_file_with_the_selector_alone(plain_bgzip, tmp_path: Path) -> None:
    """No parser, no cat: sed opens the file itself; the sink is one bgzip."""
    fa = tmp_path / "asm.fa"
    fa.write_text(">a\nACGT\n")
    fq = tmp_path / "r.fastq"
    fq.write_text("@a\nACGT\n+\nIIII\n")
    out = tmp_path / "n.gz"
    assert _scan_pipeline(fa, out) == f"sed -n '/^>/p' {fa} | sed '{_NORMALISE}' | bgzip > {out}"
    assert _scan_pipeline(fq, out) == f"sed -n 'p;n;n;n' {fq} | sed '{_NORMALISE}' | bgzip > {out}"


def test_scan_decompresses_gz_in_front_of_the_selector(plain_bgzip, tmp_path: Path) -> None:
    fq = tmp_path / "r.fq.gz"
    with gzip.open(fq, "wt") as fh:
        fh.write("@a\nACGT\n+\nIIII\n")
    out = tmp_path / "n.gz"
    assert _scan_pipeline(fq, out, threads=8) == (
        f"gzip -dc {fq} | {_FQ_SELECT} | bgzip -@ 8 > {out}"
    )


def test_scan_quotes_paths(plain_bgzip, tmp_path: Path) -> None:
    fq = tmp_path / "my reads.fq"
    fq.write_text("@a\nACGT\n+\nIIII\n")
    pipeline = _scan_pipeline(fq, tmp_path / "out dir" / "n.gz")
    assert f"'{fq}'" in pipeline and "out dir/n.gz'" in pipeline


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


def test_sink_selects_fasta_headers_then_bgzips(plain_bgzip, tmp_path: Path) -> None:
    sidecar = tmp_path / "with space.txt.gz"
    assert names_sink(sidecar, threads=3) == f"{_FA_SELECT} | bgzip -@ 3 > '{sidecar}'"
    fifo = tmp_path / "decode.fifo"
    assert names_sink(sidecar, stdin_from=fifo) == (
        f"sed -n '/^>/p' < {fifo} | sed '{_NORMALISE}' | bgzip > '{sidecar}'"
    )


def test_write_runs_one_linear_pipefail_pipeline(
    plain_bgzip, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    captured: dict[str, object] = {}

    class _Result:
        returncode = 0
        stderr = b""

    def _fake_run(cmd, **kwargs):
        captured["cmd"] = cmd
        return _Result()

    monkeypatch.setattr(qn.subprocess, "run", _fake_run)
    fq = tmp_path / "reads.fq.gz"
    with gzip.open(fq, "wt") as fh:
        fh.write("@a\nACGT\n+\nIIII\n")

    sidecar = tmp_path / "sub" / "names.txt.gz"
    write_query_names_sidecar(fq, sidecar, threads=8)
    cmd = captured["cmd"]
    assert cmd[:3] == ["bash", "-o", "pipefail"]
    pipeline = cmd[4]
    assert pipeline.startswith(f"gzip -dc {fq} | sed -n")
    assert pipeline.endswith(f"| bgzip -@ 8 > {sidecar}")
    assert ">(" not in pipeline, "linear pipeline, no process substitution"
    assert sidecar.parent.is_dir(), "parent directory created up front"


def test_write_uses_samtools_for_an_alignment(
    plain_bgzip, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
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
    assert captured["cmd"][4].startswith(
        f"/bin/samtools fasta -F 0x900 -N -@ 2 aln.bam | {_FA_SELECT}"
    )


def test_write_reports_a_failed_stage(
    plain_bgzip, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    class _Result:
        returncode = 1
        stderr = b"sed: boom"

    fq = tmp_path / "reads.fq"
    fq.write_text("@a\nACGT\n+\nIIII\n")
    monkeypatch.setattr(qn.subprocess, "run", lambda *a, **kw: _Result())
    with pytest.raises(ExternalToolError) as excinfo:
        write_query_names_sidecar(fq, tmp_path / "n.txt.gz", capture=True)
    assert "boom" in str(excinfo.value)


# --- against the real backend --------------------------------------------

requires_hks = pytest.mark.skipif(shutil.which("hks") is None, reason="hks binary not on PATH")


@requires_hks
@requires_bgzip
def test_sidecar_names_equal_the_names_hks_reports(tmp_path: Path) -> None:
    """The one test that pins the contract: for the same file, the sidecar's
    names are the names hks itself writes with --report-query-names, in order.
    Exercised on the awkward headers (leading space, CRLF, VT, FF, tab)."""
    import random

    from karyoscope.core.build import build_database
    from karyoscope.core.buildspec import BuildSpec, FeatureSetSpec
    from karyoscope.core.io.hks import run_hks_lookup

    rng = random.Random(3)
    genome_seq = "".join(rng.choice("ACGT") for _ in range(400))
    genome = tmp_path / "genome.fa"
    genome.write_text(f">chr1\n{genome_seq}\n")
    bed = tmp_path / "gene.bed"
    bed.write_text("chr1\t0\t200\texon\nchr1\t200\t400\tintron\n")
    spec = BuildSpec(
        id="HKS_names", version="1.0.0", sequence=genome, s=11, threads=1, mem_gigas=1,
        feature_sets=[FeatureSetSpec(name="gene", bed=bed)],
    )  # fmt: skip
    db_root = tmp_path / "db"
    build_database(spec, db_root, register=False, force=True)
    db_dir = db_root / "HKS_names"

    queries = tmp_path / "queries.fa"
    reads = [genome_seq[i : i + 40] for i in range(0, 200, 40)]
    queries.write_bytes(
        b"> lead desc\n" + reads[0].encode() + b"\n"
        b">crlf\r\n" + reads[1].encode() + b"\n"
        b">vt\x0bdesc\n" + reads[2].encode() + b"\n"
        b">tab\tt1 t2\n" + reads[3].encode() + b"\n"
        b">ff\x0cdesc\n" + reads[4].encode() + b"\n"
    )
    out = tmp_path / "lookup.bed"
    run_hks_lookup(
        base_path=db_dir / "index/features.hksb",
        feature_set_file=db_dir / "index/features.gene.hksf",
        k=11,
        input_path=queries,
        output_path=out,
        report_query_names=True,
        capture=True,
    )
    hks_names = list(dict.fromkeys(line.split("\t")[0] for line in out.read_text().splitlines()))
    sidecar = tmp_path / "queries.txt.gz"
    write_query_names_sidecar(queries, sidecar)
    assert _names(sidecar) == hks_names == ["lead", "crlf", "vt", "tab", "ff"]
