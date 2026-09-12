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
* **Teed off a decode that is happening anyway.** The HKS backend materialises
  a BAM/CRAM to a temp FASTA before querying, and the KMC backend streams the
  decode into ``get_featureIDs``; both tee the names off that one decode
  (:func:`names_sink`), so an alignment is never decoded twice for this.
* **A scan pass of its own** (:func:`write_query_names_sidecar`). FASTA and
  FASTQ are read directly by both query tools, so there is no decode to tee
  off, and the names take one extra streaming read of the input. Against
  hks, which reads the input once per feature set, that is a small fraction;
  against ``get_featureIDs``, which reads it once for all feature sets, it is
  one more read of equal size.

The scan is ``sed`` over the (decompressed) file and nothing else: no
sequence parser, because ``hks`` itself imposes the shape the scan relies on.
Its reader (``jseqio``) takes a FASTQ record as exactly four lines, so record
N's header is line 4N+1 and a FASTQ the scan would misread is one hks
rejects; it tells FASTA from FASTQ, and gzip from plain, by content, which
the scan mirrors (:func:`sniff_fastx`); and it names a record by dropping the
first byte and taking the first whitespace-delimited token, which is what
``_NORMALISE`` does.
"""

from __future__ import annotations

import gzip
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
#: FASTA or FASTQ, told apart by CONTENT (see :func:`sniff_fastx`).
_ALIGNMENT_EXTENSIONS: tuple[str, ...] = (".bam", ".cram")

#: Header line -> name, exactly as ``hks`` does it (``load_seq_names``):
#: drop the record's first byte (``>``/``@``), then take the first
#: whitespace-delimited token -- ``split_whitespace().next()`` -- which skips
#: leading whitespace and stops at any whitespace, tab, CR, VT and FF
#: included, so a ``@ read1`` header is ``read1`` and a CRLF header loses its
#: ``\r``. ``[[:space:]]`` is the POSIX spelling of that set (ASCII; hks also
#: treats non-ASCII Unicode spaces as separators, which no header has).
_NORMALISE = "s/^.//;s/^[[:space:]]*//;s/[[:space:]].*//"

#: Header lines of a FASTQ stream, by POSITION: print a line, skip three.
#: Headers are not matched on their leading ``@`` because a quality line can
#: begin with one (Phred 31). Four lines per record is not an assumption of
#: ours -- it is how ``hks`` parses FASTQ, so this is the only shape whose
#: ranks exist to be mapped. sed rather than awk because it is the fast way
#: to say "every fourth line": 1.0 s against 10.6 s for 2 M reads with the
#: BSD awk on macOS. Normalisation is a SECOND sed, downstream, so it runs on
#: the header lines only and on another core: 1.4 s for the pair against
#: 2.1 s with the substitutions folded into the selector. ``{src}`` is
#: `` <path>`` or empty (stdin).
_FASTQ_SELECT = f"sed -n 'p;n;n;n'{{src}} | sed '{_NORMALISE}'"

#: Header lines of a FASTA stream: the ``>`` lines, then the same normaliser.
#: A sequence line fails the anchored match at its first byte, so a
#: multi-line assembly costs one comparison per line.
_FASTA_SELECT = f"sed -n '/^>/p'{{src}} | sed '{_NORMALISE}'"


def sniff_fastx(path: Path) -> tuple[bool, bool]:
    """``(gzipped, fastq)`` for ``path``, decided from its bytes like ``hks`` does.

    ``hks`` (jseqio) ignores the filename: gzip is the two magic bytes, and the
    format is the first byte of the decompressed stream, ``>`` or ``@``. A
    scan keyed on extensions would silently write an EMPTY sidecar for an
    extensionless FASTQ that hks annotated fine, so this looks at the same
    bytes. An empty file is FASTA with no records (an empty sidecar, which is
    right); any other first byte is refused.
    """
    with path.open("rb") as fh:
        magic = fh.read(2)
    gzipped = magic == b"\x1f\x8b"
    if gzipped:
        with gzip.open(path, "rb") as fh:
            first = fh.read(1)
    else:
        first = magic[:1]
    if first == b"@":
        return gzipped, True
    if first in (b">", b""):
        return gzipped, False
    raise KaryoscopeError(
        f"{path.name} is not FASTA or FASTQ (first byte {first!r}); cannot list its "
        f"record names for --query-names-sidecar."
    )


def names_sink(sidecar: Path, threads: int = 1, *, stdin_from: Path | None = None) -> str:
    """The pipeline tail that turns a FASTA stream into the sidecar at ``sidecar``.

    Returned as shell text (``sed ... | bgzip > <sidecar>``) for the caller to
    append to a producer with ``|``, or, with ``stdin_from``, to run on its own
    reading that path (a FIFO a ``tee`` writes into). Keeping the tail here
    means the tee in :func:`karyoscope.core.io.hks.materialised_queries`, the
    one in :mod:`karyoscope.core.io.kmc` and the scan pass in
    :func:`write_query_names_sidecar` cannot drift into two formats.
    """
    src = "" if stdin_from is None else f" < {shlex.quote(str(stdin_from))}"
    select = _FASTA_SELECT.format(src=src)
    return f"{select} | {bgzip_stage(threads)} > {shlex.quote(str(sidecar))}"


def _scan_pipeline(input_path: Path, sidecar: Path, threads: int = 1) -> str:
    """Shell text producing ``sidecar`` from a FASTA/FASTQ ``input_path``.

    The header selector reads a plain file itself; a gzipped one is
    decompressed into it. The whole input is read once -- there is no way
    to find the headers without reading the lines between them -- but
    nothing is spent parsing sequence when only the header line of each
    record is wanted: measured on 2 M reads, the pipeline runs at ~650 MB/s
    of FASTQ on a laptop.
    """
    gzipped, fastq = sniff_fastx(input_path)
    select = _FASTQ_SELECT if fastq else _FASTA_SELECT
    quoted_in = shlex.quote(str(input_path))
    tail = f"{bgzip_stage(threads)} > {shlex.quote(str(sidecar))}"
    if gzipped:
        return f"gzip -dc {quoted_in} | {select.format(src='')} | {tail}"
    return f"{select.format(src=' ' + quoted_in)} | {tail}"


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
        producer = ["sed"]
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
