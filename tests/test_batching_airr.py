"""Phase 4 AIRR streaming pipeline (#26): batch-size equivalence.

For the example AIRR input, running the CLI with ``--batch-size`` in
``{1, 2, 50, 10000}`` produces output that is structurally equal to
the existing golden the legacy path generated.
"""

from __future__ import annotations

import gzip
import json
import subprocess
from pathlib import Path

import pytest

from .test_cli_processing import (
    compare_consolidated_files,
    normalize_json,
    strip_volatile_fields,
)

EXAMPLE_DATA = Path(__file__).parent.parent / "example-data"
AIRR_INPUT = EXAMPLE_DATA / "airr" / "input-airr.json"
AIRR_GOLDEN = EXAMPLE_DATA / "airr" / "airr-olmsted-golden.json"
AIRR_GOLDEN_GZ = EXAMPLE_DATA / "airr" / "airr-olmsted-golden.json.gz"


def _run_airr(output, name, batch_size):
    cmd = [
        "olmsted",
        "process",
        "-f",
        "airr",
        "-i",
        str(AIRR_INPUT),
        "-o",
        str(output),
        "--seed",
        "42",
        "--name",
        name,
        "--batch-size",
        str(batch_size),
        "-q",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, (
        f"CLI failed (batch_size={batch_size}):\nSTDOUT: {result.stdout}\n"
        f"STDERR: {result.stderr}"
    )


@pytest.mark.parametrize("batch_size", [1, 2, 50, 10000])
def test_airr_batch_size_matches_golden(batch_size, tmp_path):
    out = tmp_path / f"out_bs{batch_size}.json"
    _run_airr(out, "airr-example", batch_size)
    ok, message = compare_consolidated_files(out, AIRR_GOLDEN)
    assert ok, f"batch_size={batch_size} diverges from golden:\n{message}"


def _run_airr_gzip(output, name, batch_size):
    cmd = [
        "olmsted", "process",
        "-f", "airr",
        "-i", str(AIRR_INPUT),
        "-o", str(output),
        "--seed", "42",
        "--name", name,
        "--batch-size", str(batch_size),
        "--json-format", "gzip",
        "-q",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, (
        f"CLI failed (batch_size={batch_size}):\nSTDOUT: {result.stdout}\n"
        f"STDERR: {result.stderr}"
    )


@pytest.mark.parametrize("batch_size", [1, 50])
def test_airr_streaming_gzip_matches_golden(batch_size, tmp_path):
    """Gzip output through ``write_olmsted_json_streaming`` must parse-equal the golden.

    Pre-fix, the only gzip integration test for AIRR ran at the default
    ``--batch-size 50``; this exercises multi-batch streaming (``=1``)
    plus the default to lock both code paths.
    """
    out_base = tmp_path / f"out_bs{batch_size}.json"
    _run_airr_gzip(out_base, "airr-example", batch_size)
    out_gz = Path(str(out_base) + ".gz")
    assert out_gz.exists()

    with gzip.open(out_gz, "rb") as fh:
        actual = json.loads(fh.read().decode("utf-8"))
    with gzip.open(AIRR_GOLDEN_GZ, "rb") as fh:
        expected = json.loads(fh.read().decode("utf-8"))
    assert normalize_json(strip_volatile_fields(actual)) == normalize_json(
        strip_volatile_fields(expected)
    )


def test_airr_streaming_gzip_is_deterministic(tmp_path):
    """Two streaming-gzip runs of the same input must produce byte-equal output.

    Locks the ``mtime=0`` + empty-filename header pinning the streaming
    writer inherits from ``data_io.write_olmsted_json``.
    """
    out_a = tmp_path / "a.json"
    out_b = tmp_path / "b.json"
    _run_airr_gzip(out_a, "airr-example", batch_size=1)
    _run_airr_gzip(out_b, "airr-example", batch_size=1)
    a_bytes = Path(str(out_a) + ".gz").read_bytes()
    b_bytes = Path(str(out_b) + ".gz").read_bytes()
    # The full payload also embeds non-deterministic metadata.created_at,
    # so check only the gzip header (first 10 bytes) for byte-equality —
    # same guarantee data_io.write_olmsted_json carries.
    assert a_bytes[:10] == b_bytes[:10], (
        "gzip header drift: streaming writer is not byte-deterministic"
    )


def test_streaming_default_path_uses_streaming_for_airr():
    """Catches regressions where _should_stream_airr silently flips back to
    legacy for the common case.
    """
    import argparse

    from olmsted_cli.process_data import _should_stream_airr

    args = argparse.Namespace(
        batch_size=50,
        split_files=None,
        mutations=None,
        validate=False,
    )
    assert _should_stream_airr(args) is True

    # Mutations: stays on streaming after phase 4 wired per-batch merge.
    args.mutations = "x.csv"
    assert _should_stream_airr(args) is True
    args.mutations = None

    # Validate: streaming bypassed until per-batch validation lands.
    args.validate = True
    assert _should_stream_airr(args) is False
    args.validate = False

    # Split-files: streaming bypassed.
    args.split_files = "/tmp/x"
    assert _should_stream_airr(args) is False
    args.split_files = None

    # batch_size 0: explicit opt-out.
    args.batch_size = 0
    assert _should_stream_airr(args) is False
