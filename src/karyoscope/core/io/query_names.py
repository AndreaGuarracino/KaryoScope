"""The rank -> record-name sidecar for rank-identified ``annotate`` output.

Read-level output is identified by ordinal rank, because ``hks`` buffers every
query name in memory before the first lookup and a human WGS sample's names
run to ~90 GB (see :func:`karyoscope.core.annotate._reject_query_names_for_reads`).
Rank N is the Nth record of the query file, so nothing is lost, but anything
that joins the output back to reads needs the mapping. ``--query-names-sidecar``
writes it next to the output.

This module owns the sidecar's content and how it is produced, so the file
looks the same whatever input format or backend made it: one name per line,
in the order the query tool assigns ranks (line N+1 is rank N); the header up
to its first whitespace, ``/1``/``/2`` mate suffix included for alignments;
bgzip-compressed like every other ``.gz`` KaryoScope writes.

There are two ways to produce it:

* **Teed off a decode that is happening anyway.** The HKS backend materialises
  a BAM/CRAM to a temp FASTA before querying; the caller builds that pipeline
  and appends :func:`names_sink` to it, so the names come off the same pass.
* **A scan pass of its own** (:func:`write_query_names_sidecar`). FASTA and
  FASTQ are read directly by both query tools, so there is no decode to tee
  off; and on the KMC backend an alignment is fed to ``get_featureIDs`` by a
  streaming pipe with no seekable copy left behind. Either way the names take
  one extra streaming read of the input -- cheap next to the lookup, which
  reads the input once per feature set, but for an alignment on KMC it is a
  full second decode, and the CLI help says so.

The FASTA/FASTQ scan is ``sed``/``grep`` over the (decompressed) file and
nothing else: no sequence parser, because ``hks`` itself imposes the shape
the scan relies on. Its reader (``jseqio``) takes a FASTQ record as exactly
four lines, so record N's header is line 4N+1 and a FASTQ the scan would
misread is one hks rejects; and it names a record by its header up to the
first whitespace, which is what the ``cut`` tail keeps.
"""

from __future__ import annotations

import logging
import shlex
import subprocess
from pathlib import Path

from karyoscope.core.external import ExternalToolError, require_tool
from karyoscope.core.io.bgzip import bgzip_stage
from karyoscope.exceptions import KaryoscopeError

logger = logging.getLogger(__name__)

#: Suffix appended to the input's basename to name its sidecar.
QUERY_NAMES_SUFFIX = ".query_names.txt.gz"

#: Alignment formats decoded with ``samtools fasta``; everything else is
#: FASTA or FASTQ and is scanned with ``awk`` directly.
_ALIGNMENT_EXTENSIONS: tuple[str, ...] = (".bam", ".cram")
_FASTQ_EXTENSIONS: tuple[str, ...] = (".fastq", ".fq", ".fastq.gz", ".fq.gz")

#: Header lines of a FASTA stream. ``{src}`` is `` <path>`` or empty (stdin).
#: ``grep`` exits 1 on no match, which under pipefail would fail an input
#: with zero records; the guard maps that one code to success and leaves a
#: real error (2) alone.
_FASTA_HEADERS = "(grep '^>'{src} || [ $? -eq 1 ])"

#: Header lines of a FASTQ stream, by POSITION: print a line, skip three.
#: Headers are not matched on their leading ``@`` because a quality line can
#: begin with one (Phred 31). Four lines per record is not an assumption of
#: ours -- it is how ``hks`` parses FASTQ, so this is the only shape whose
#: ranks exist to be mapped. sed rather than awk because it is the fast way
#: to say "every fourth line": 1.0 s against 10.6 s for 2 M reads with the
#: BSD awk on macOS, and no regex per line anywhere.
_FASTQ_HEADERS = "sed -n 'p;n;n;n'{src}"

#: Header line -> name: up to the first whitespace (tab or space, the same
#: split ``hks`` makes), minus the leading ``>``/``@``. This runs on the
#: header lines only, a small fraction of the input, so the three stages cost
#: nothing measurable.
_HEADER_TO_NAME = "cut -f1 | cut -d ' ' -f1 | cut -c2-"


def names_sink(sidecar: Path, threads: int = 1) -> str:
    """The pipeline tail that turns a FASTA stream into the sidecar at ``sidecar``.

    Returned as shell text (``grep ... | cut ... | bgzip > <sidecar>``) for the
    caller to append to a producer with ``|``. Keeping the tail here means the
    tee in :func:`karyoscope.core.io.hks.materialised_queries` and the scan
    pass in :func:`write_query_names_sidecar` cannot drift into two formats.
    """
    headers = _FASTA_HEADERS.format(src="")
    return f"{headers} | {_HEADER_TO_NAME} | {bgzip_stage(threads)} > {shlex.quote(str(sidecar))}"


def _scan_pipeline(input_path: Path, sidecar: Path, threads: int = 1) -> str:
    """Shell text producing ``sidecar`` from a FASTA/FASTQ ``input_path``.

    The header selector reads a plain file itself; a gzipped one is
    decompressed into it. Nothing is spent parsing sequence when only the
    header line of each record is wanted: measured on 2 M reads, the whole
    pipeline runs at ~650 MB/s of FASTQ on a laptop, under half the time of
    one ``hks lookup`` over the same file against a toy index, and a far
    smaller fraction against a human one.
    """
    name = input_path.name.lower()
    headers = _FASTQ_HEADERS if name.endswith(_FASTQ_EXTENSIONS) else _FASTA_HEADERS
    quoted_in = shlex.quote(str(input_path))
    tail = f"{_HEADER_TO_NAME} | {bgzip_stage(threads)} > {shlex.quote(str(sidecar))}"
    if name.endswith(".gz"):
        return f"gzip -dc {quoted_in} | {headers.format(src='')} | {tail}"
    return f"{headers.format(src=' ' + quoted_in)} | {tail}"


def alignment_decode_cmd(
    input_path: Path,
    *,
    reference: Path | None = None,
    threads: int = 0,
) -> list[str]:
    """The ``samtools fasta`` command that streams a BAM/CRAM's records as FASTA.

    The SAME invocation the query decode uses (``-F 0x900 -N``: primary
    records only, mate suffix forced on), so a sidecar produced from it lines
    up with the ranks ``hks`` or ``get_featureIDs`` assigned.
    """
    suffix = input_path.suffix.lower()
    if suffix in _ALIGNMENT_EXTENSIONS:
        if suffix == ".cram" and reference is None:
            raise KaryoscopeError(
                f"{input_path.name} is a CRAM, which stores bases as a diff "
                f"against the reference it was aligned to, so it cannot be "
                f"decoded without that reference. Pass --reference "
                f"<genome.fasta> (the same one used for alignment)."
            )
        samtools = require_tool(
            "samtools",
            install_hint=(
                "Install samtools to use BAM/CRAM inputs:\n  conda install -c bioconda samtools"
            ),
        )
        cmd = [samtools, "fasta", "-F", "0x900", "-N"]
        if reference is not None:
            # NOT -T: in `samtools fasta` that is the copy-tags-to-header
            # taglist, and a path there silently yields a tag-decorated FASTA
            # with no reference supplied at all.
            cmd += ["--reference", str(reference)]
        if threads > 0:
            cmd += ["-@", str(threads)]
        cmd.append(str(input_path))
        return cmd
    raise ValueError(f"not an alignment: {input_path.name}")


def write_query_names_sidecar(
    input_path: Path,
    sidecar: Path,
    *,
    reference: Path | None = None,
    threads: int = 0,
    capture: bool = False,
) -> Path:
    """Write ``input_path``'s rank -> name mapping to ``sidecar`` in one streaming pass.

    A BAM/CRAM is decoded with ``samtools fasta`` and the names taken off
    that stream; a FASTA/FASTQ is read by the header selector directly (see
    :func:`_scan_pipeline`). Either pipeline is LINEAR under ``pipefail``, so
    every stage's failure -- a full disk under the sidecar included -- is the
    run's failure, and every stage has finished when bash returns.
    ``capture`` collects stderr for the error message instead of letting it
    interleave with our own logging.
    """
    sidecar.parent.mkdir(parents=True, exist_ok=True)
    if input_path.suffix.lower() in _ALIGNMENT_EXTENSIONS:
        producer = alignment_decode_cmd(input_path, reference=reference, threads=threads)
        pipeline = f"{shlex.join(producer)} | {names_sink(sidecar, threads)}"
    else:
        producer = ["sed" if input_path.name.lower().endswith(_FASTQ_EXTENSIONS) else "grep"]
        pipeline = _scan_pipeline(input_path, sidecar, threads)
    logger.debug("scanning query names: %s", pipeline)
    result = subprocess.run(
        ["bash", "-o", "pipefail", "-c", pipeline],
        stderr=subprocess.PIPE if capture else None,
        check=False,
    )
    if result.returncode != 0:
        stderr = result.stderr.decode() if result.stderr else ""
        raise ExternalToolError(cmd=producer, returncode=result.returncode, stderr=stderr)
    logger.info("wrote query-name sidecar %s", sidecar)
    return sidecar
