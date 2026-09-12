"""``bgzip`` compression, the one compressor for every ``.gz`` KaryoScope writes.

Three shapes, one tool:

* :func:`bgzip_file` compresses a finished file in place (``annotate``,
  ``scaffold``, ``centromeres``, ``karyotype`` compress their outputs after
  writing them).
* :func:`open_bgzip_writer` streams text into ``bgzip`` as it is produced, for
  writers that build their output line by line (``bin``, ``remap-bed``, the
  FASTA/BED rewriters) -- no plain copy ever touches the disk.
* :func:`bgzip_stage` is the compressor as a shell pipeline stage, for the
  sidecar pipelines in :mod:`karyoscope.core.io.query_names`.

bgzip rather than gzip everywhere because its output *is* a gzip stream --
``zcat``, Python's ``gzip`` module and every other reader are unaffected --
while it compresses in parallel (``-@``) and is tabix-indexable. The price is
a slightly larger file from the 64 KB block boundaries.
"""

from __future__ import annotations

import io
import logging
import shlex
import subprocess
import time
from pathlib import Path

from karyoscope.core.external import ExternalToolError, require_tool, run_tool
from karyoscope.diskspace import format_bytes

logger = logging.getLogger(__name__)

_INSTALL_HINT = "Install htslib (`conda install -c bioconda htslib`) to write .gz output."


def _bgzip_cmd(threads: int) -> list[str]:
    """``bgzip`` reading stdin and writing stdout, with ``-@`` when it helps."""
    cmd = [require_tool("bgzip", install_hint=_INSTALL_HINT)]
    if threads > 1:
        cmd += ["-@", str(threads)]
    return cmd


def bgzip_stage(threads: int = 1) -> str:
    """Shell text for ``bgzip`` as a pipeline stage (stdin -> stdout)."""
    return shlex.join(_bgzip_cmd(threads))


class _BgzipWriter(io.TextIOWrapper):
    """A text handle whose bytes go through a ``bgzip`` child into ``path``.

    Closing it closes the child's stdin, waits for it to finish, and raises
    :class:`ExternalToolError` if it failed -- a full disk under the output
    surfaces as an exception at close, not as a silently truncated file.
    """

    def __init__(self, path: Path, threads: int):
        self._path = path
        self._sink = path.open("wb")
        try:
            self._proc = subprocess.Popen(
                _bgzip_cmd(threads),
                stdin=subprocess.PIPE,
                stdout=self._sink,
                stderr=subprocess.PIPE,
            )
        except BaseException:
            self._sink.close()
            raise
        assert self._proc.stdin is not None
        super().__init__(self._proc.stdin, encoding="utf-8", write_through=False)

    def close(self) -> None:
        if self.closed:
            return
        # A compressor that died early makes the final flush fail with
        # BrokenPipeError. That exception is the symptom; the child's exit
        # status and stderr are the cause, so reap the child first and let
        # its failure be the error the caller sees.
        flush_error: BaseException | None = None
        try:
            super().close()  # flushes and closes the child's stdin
        except (OSError, ValueError) as e:
            flush_error = e
        finally:
            # Not communicate(): it would try to flush the stdin just closed.
            assert self._proc.stderr is not None
            stderr = self._proc.stderr.read()
            self._proc.stderr.close()
            self._proc.wait()
            self._sink.close()
        if self._proc.returncode != 0:
            raise ExternalToolError(
                cmd=self._proc.args,
                returncode=self._proc.returncode,
                stderr=stderr.decode(errors="replace"),
            ) from flush_error
        if flush_error is not None:
            raise flush_error


def open_bgzip_writer(path: Path, threads: int = 1) -> io.TextIOWrapper:
    """Open ``path`` for text writing through a streaming ``bgzip``.

    A drop-in for ``gzip.open(path, "wt")``: use it as a context manager or
    close it explicitly. Compression runs in the child as the text is
    written, so the caller's peak disk usage is the compressed file alone.
    """
    return _BgzipWriter(path, threads)


def bgzip_file(path: Path, threads: int = 1) -> Path:
    """Compress ``path`` in-place with ``bgzip``, returning the new path.

    ``bgzip`` removes the source file by default (matches gzip's behaviour).
    Returns ``Path(str(path) + ".gz")``. Logs per-file start + completion
    at INFO so a long bgzip pass (12 files for a 6-feature-set human
    database) doesn't look like the pipeline has hung.

    ``threads`` is forwarded as ``bgzip -@``; the htslib bgzip compresses
    a single file in parallel when given more than one thread. We
    process files sequentially within the bgzip pass, so passing the
    user's full ``--threads`` here is the right call (no contention
    with concurrent file compressions). ``threads=1`` (the default)
    omits ``-@`` entirely for cleanest subprocess invocation.
    """
    bgzip = require_tool(
        "bgzip",
        install_hint="Install htslib (`conda install -c bioconda htslib`), "
        "or rerun with --no-bgzip to skip compression.",
    )
    orig_size = path.stat().st_size
    logger.info("bgzipping %s (%s, threads=%d)", path.name, format_bytes(orig_size), threads)
    t0 = time.perf_counter()
    cmd = [bgzip, "-f"]
    if threads > 1:
        cmd.extend(["-@", str(threads)])
    cmd.append(str(path))
    run_tool(cmd)
    out_path = Path(str(path) + ".gz")
    out_size = out_path.stat().st_size if out_path.is_file() else 0
    dt = time.perf_counter() - t0
    logger.info(
        "bgzipped %s (%s -> %s) in %.1fs",
        out_path.name,
        format_bytes(orig_size),
        format_bytes(out_size),
        dt,
    )
    return out_path
