"""Tests for the ``karyoscope build`` command wiring (no HKS binary needed).

The full construction path is covered by ``test_build_core.py`` (skipped when
``hks`` is absent); here we stub :func:`karyoscope.core.build.build_database` to
exercise spec assembly, flag/spec mutual exclusion, and error surfacing.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from click.testing import CliRunner

import karyoscope.commands.build as build_cmd
from karyoscope.cli import main
from karyoscope.core.build import BuildResult
from karyoscope.core.buildspec import BuildSpec


@pytest.fixture
def cli_runner() -> CliRunner:
    return CliRunner()


@pytest.fixture
def inputs(tmp_path: Path) -> tuple[Path, Path]:
    genome = tmp_path / "g.fa"
    genome.write_text(">chr1\nACGTACGTACGT\n")
    bed = tmp_path / "r.bed"
    bed.write_text("chr1\t0\t8\tLINE\n")
    return genome, bed


@pytest.fixture
def stub_build(monkeypatch: pytest.MonkeyPatch) -> dict:
    """Replace build_database with a recorder; return the captured spec/kwargs."""
    captured: dict = {}

    def _fake(spec: BuildSpec, db_root: Path, **kwargs: object) -> BuildResult:
        captured["spec"] = spec
        captured["db_root"] = db_root
        captured["kwargs"] = kwargs
        return BuildResult(
            db_id=spec.id, db_dir=db_root / spec.id, registered=kwargs.get("register") is not False
        )

    monkeypatch.setattr(build_cmd, "build_database", _fake)
    return captured


def test_simple_form_builds_spec(cli_runner: CliRunner, inputs, stub_build, tmp_path: Path) -> None:
    genome, bed = inputs
    result = cli_runner.invoke(
        main,
        [
            "build",
            "--id",
            "HKS_x",
            "--sequence",
            str(genome),
            "--feature-set",
            f"repeat={bed}",
            "--background",
            "repeat=nonrepeat",
            "-s",
            "15",
            "--db-root",
            str(tmp_path / "db"),
        ],
        catch_exceptions=False,
    )
    assert result.exit_code == 0, result.output
    spec = stub_build["spec"]
    assert spec.id == "HKS_x"
    assert spec.s == 15
    (fs,) = spec.feature_sets
    assert fs.name == "repeat" and fs.background == "nonrepeat"


def test_spec_and_flags_mutually_exclusive(
    cli_runner: CliRunner, inputs, stub_build, tmp_path: Path
) -> None:
    genome, bed = inputs
    spec_file = tmp_path / "spec.yaml"
    spec_file.write_text(
        f'id: X\nversion: "1"\nsequence: {genome}\nfeature_sets:\n  - name: r\n    bed: {bed}\n'
    )
    result = cli_runner.invoke(
        main,
        ["build", "--spec", str(spec_file), "--id", "HKS_x", "--feature-set", f"repeat={bed}"],
    )
    assert result.exit_code != 0
    assert "cannot be combined" in result.output


def test_missing_id_errors(cli_runner: CliRunner, inputs, stub_build) -> None:
    genome, bed = inputs
    result = cli_runner.invoke(
        main, ["build", "--sequence", str(genome), "--feature-set", f"r={bed}"]
    )
    assert result.exit_code != 0
    assert "--id is required" in result.output


def test_no_register_flag_passed_through(
    cli_runner: CliRunner, inputs, stub_build, tmp_path: Path
) -> None:
    genome, bed = inputs
    result = cli_runner.invoke(
        main,
        [
            "build",
            "--id",
            "HKS_x",
            "--sequence",
            str(genome),
            "--feature-set",
            f"repeat={bed}",
            "--no-register",
            "--db-root",
            str(tmp_path / "db"),
        ],
        catch_exceptions=False,
    )
    assert result.exit_code == 0, result.output
    assert stub_build["kwargs"]["register"] is False


def test_spec_form_reads_yaml(cli_runner: CliRunner, inputs, stub_build, tmp_path: Path) -> None:
    genome, bed = inputs
    spec_file = tmp_path / "spec.yaml"
    spec_file.write_text(
        f'id: HKS_spec\nversion: "3.1.0"\nsequence: {genome}\n'
        f"feature_sets:\n  - name: repeat\n    bed: {bed}\n"
    )
    result = cli_runner.invoke(
        main,
        ["build", "--spec", str(spec_file), "--db-root", str(tmp_path / "db")],
        catch_exceptions=False,
    )
    assert result.exit_code == 0, result.output
    assert stub_build["spec"].id == "HKS_spec"
    assert stub_build["spec"].version == "3.1.0"


def test_spec_form_keeps_spec_build_block_when_flags_absent(
    cli_runner: CliRunner, inputs, stub_build, tmp_path: Path
) -> None:
    """Click defaults (threads 4, mem_gigas 8) must not clobber the spec's build block."""
    genome, bed = inputs
    spec_file = tmp_path / "spec.yaml"
    spec_file.write_text(
        f'id: HKS_spec\nversion: "1.0.0"\nsequence: {genome}\n'
        f"build:\n  threads: 16\n  mem_gigas: 32\n  external_memory: {tmp_path / 'ext'}\n"
        "  forward_only: true\n"
        f"feature_sets:\n  - name: repeat\n    bed: {bed}\n"
    )
    result = cli_runner.invoke(
        main,
        ["build", "--spec", str(spec_file), "--db-root", str(tmp_path / "db")],
        catch_exceptions=False,
    )
    assert result.exit_code == 0, result.output
    spec = stub_build["spec"]
    assert spec.threads == 16
    assert spec.mem_gigas == 32
    assert spec.external_memory == tmp_path / "ext"
    assert spec.forward_only is True


def test_spec_form_tuning_flags_override_spec(
    cli_runner: CliRunner, inputs, stub_build, tmp_path: Path
) -> None:
    """--external-memory & co. given alongside --spec must reach the build, not be dropped.

    Regression: ``--spec build.yaml --external-memory DIR`` used to run the
    in-memory k-mer sort and get OOM-killed, because the flag never left the
    command layer.
    """
    genome, bed = inputs
    spec_file = tmp_path / "spec.yaml"
    spec_file.write_text(
        f'id: HKS_spec\nversion: "1.0.0"\nsequence: {genome}\n'
        f"build:\n  threads: 16\n  mem_gigas: 32\n"
        f"feature_sets:\n  - name: repeat\n    bed: {bed}\n"
    )
    result = cli_runner.invoke(
        main,
        [
            "build",
            "--spec",
            str(spec_file),
            "--db-root",
            str(tmp_path / "db"),
            "--threads",
            "3",
            "--external-memory",
            str(tmp_path / "scratch"),
            "--forward-only",
        ],
        catch_exceptions=False,
    )
    assert result.exit_code == 0, result.output
    spec = stub_build["spec"]
    assert spec.threads == 3
    assert spec.external_memory == tmp_path / "scratch"
    assert spec.forward_only is True
    # not passed on the command line -> the spec's value survives
    assert spec.mem_gigas == 32


def test_apply_tuning_overrides_without_click_context(tmp_path: Path) -> None:
    """Called from Python (no click context): the two nullable options override when set."""
    spec = BuildSpec(id="x", version="1.0.0", feature_sets=[], threads=16, mem_gigas=32)
    out = build_cmd._apply_tuning_overrides(
        spec, threads=4, mem_gigas=8, external_memory=tmp_path, forward_only=True
    )
    assert out.external_memory == tmp_path
    assert out.forward_only is True
    assert (out.threads, out.mem_gigas) == (16, 32)


@pytest.mark.parametrize(
    "flags,expected_threads,expected_memory",
    [
        (["--threads", "4"], 4, 32),
        (["--mem-gigas", "8"], 16, 8),
        (["--mem-gigas", "12"], 16, 12),
        (["-t", "4", "--mem-gigas", "8"], 4, 8),
    ],
)
def test_spec_explicit_tuning_defaults_override_yaml(
    cli_runner: CliRunner,
    inputs,
    stub_build,
    tmp_path: Path,
    flags: list[str],
    expected_threads: int,
    expected_memory: int,
) -> None:
    genome, bed = inputs
    spec_file = tmp_path / "spec.yaml"
    spec_file.write_text(
        f'id: HKS_spec\nversion: "1.0.0"\nsequence: {genome}\n'
        "build: {threads: 16, mem_gigas: 32, forward_only: true}\n"
        f"feature_sets:\n  - name: repeat\n    bed: {bed}\n"
    )
    result = cli_runner.invoke(
        main,
        ["build", "--spec", str(spec_file), "--db-root", str(tmp_path / "db"), *flags],
        catch_exceptions=False,
    )
    assert result.exit_code == 0, result.output
    spec = stub_build["spec"]
    assert (spec.threads, spec.mem_gigas) == (expected_threads, expected_memory)
    assert spec.forward_only is True


@pytest.mark.parametrize(
    "flag,value",
    [
        ("--id", "another"),
        ("--sequence", "{genome}"),
        ("--feature-set", "repeat={bed}"),
        ("--background", "repeat=other"),
        ("--exclude", "chr1"),
        ("--flatten-order", "repeat={bed}"),
        ("--hierarchy", "repeat={bed}"),
        ("--priority", "repeat={bed}"),
        ("--colors", "repeat={bed}"),
        ("--flatten", None),
        ("--variable-k", None),
        ("--s", "51"),
        ("--s", "31"),
        ("--db-version", "9.0.0"),
        ("--db-version", "1.0.0"),
    ],
)
def test_spec_rejects_explicit_database_definition_options(
    cli_runner: CliRunner,
    inputs,
    stub_build,
    tmp_path: Path,
    flag: str,
    value: str | None,
) -> None:
    genome, bed = inputs
    spec_file = tmp_path / "spec.yaml"
    spec_file.write_text(
        f'id: HKS_spec\nversion: "1.0.0"\nsequence: {genome}\n'
        f"feature_sets:\n  - name: repeat\n    bed: {bed}\n"
    )
    flags = [flag] if value is None else [flag, value.format(genome=genome, bed=bed)]
    result = cli_runner.invoke(
        main,
        ["build", "--spec", str(spec_file), "--db-root", str(tmp_path / "db"), *flags],
    )
    assert result.exit_code == 2, result.output
    assert f"--spec cannot be combined with {flag}" in result.output
    assert "spec" not in stub_build  # Reject before doing any expensive build work.
