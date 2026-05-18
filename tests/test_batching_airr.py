"""Phase 4 AIRR streaming pipeline (#26): batch-size equivalence.

For the example AIRR input, running the CLI with ``--batch-size`` in
``{1, 2, 50, 10000}`` produces output that is structurally equal to
the existing golden the legacy path generated.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from .test_cli_processing import compare_consolidated_files

EXAMPLE_DATA = Path(__file__).parent.parent / "example-data"
AIRR_INPUT = EXAMPLE_DATA / "airr" / "input-airr.json"
AIRR_GOLDEN = EXAMPLE_DATA / "airr" / "airr-olmsted-golden.json"


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
