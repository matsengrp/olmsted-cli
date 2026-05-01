"""Tests for gzip JSON I/O across process / merge / tag.

Covers:
- writing ``.json.gz`` output (``--json-format gzip``) on all three commands
- reading a gzipped Olmsted JSON as input on ``tag`` and ``merge``

Comparisons of decompressed content go through
``compare_consolidated_files`` (volatile metadata stripped) — gzip output
is byte-deterministic at the compression layer (``mtime=0``, empty filename)
but the underlying JSON content has timestamp churn between runs.
"""

import gzip
import json
import shutil
import subprocess
from pathlib import Path

import pytest

from .test_cli_processing import compare_consolidated_files

REPO_ROOT = Path(__file__).parent.parent
EXAMPLE = REPO_ROOT / "example-data"


def _gz_output_path(out_path):
    """Path that --json-format gzip writes to (auto-appends .gz)."""
    return Path(str(out_path) + ".gz")


def _gzip_copy(src, dst):
    """Compress src (text) into dst (.gz). For building gzipped test inputs."""
    with open(src, "rb") as sfh, gzip.open(dst, "wb") as dfh:
        shutil.copyfileobj(sfh, dfh)


# --- write side ---


def test_process_pcp_gzip_output_matches_golden(tmp_path):
    """`process -f pcp --json-format gzip` decompresses to match the plain golden."""
    out = tmp_path / "out.json"
    subprocess.run(
        [
            "olmsted", "process", "-f", "pcp",
            "-i", str(EXAMPLE / "pcp" / "input-pcp.csv"),
            "-t", str(EXAMPLE / "pcp" / "input-trees.csv"),
            "-o", str(out),
            "--seed", "42", "--name", "pcp-example",
            "--json-format", "gzip", "-q",
        ],
        check=True, capture_output=True,
    )

    gz_out = _gz_output_path(out)
    assert gz_out.exists(), f"Expected gzipped output at {gz_out}"
    assert not out.exists(), "Plain .json should not exist when --json-format gzip"

    decompressed = tmp_path / "decompressed.json"
    with gzip.open(gz_out, "rt") as fh, open(decompressed, "w") as ofh:
        ofh.write(fh.read())

    match, message = compare_consolidated_files(
        str(EXAMPLE / "pcp" / "pcp-olmsted-golden.json"),
        str(decompressed),
    )
    assert match, f"Decompressed gzip output doesn't match golden:\n{message}"


def test_process_airr_gzip_output_matches_golden(tmp_path):
    """`process -f airr --json-format gzip` propagates json_format through airr_args.

    Regression: airr_args used to drop json_format, silently producing
    pretty JSON regardless of the user's --json-format flag.
    """
    out = tmp_path / "out.json"
    subprocess.run(
        [
            "olmsted", "process", "-f", "airr",
            "-i", str(EXAMPLE / "airr" / "input-airr.json"),
            "-o", str(out),
            "--seed", "42", "--name", "airr-example",
            "--json-format", "gzip", "-q",
        ],
        check=True, capture_output=True,
    )

    gz_out = _gz_output_path(out)
    assert gz_out.exists()
    assert not out.exists()

    decompressed = tmp_path / "decompressed.json"
    with gzip.open(gz_out, "rt") as fh, open(decompressed, "w") as ofh:
        ofh.write(fh.read())

    match, message = compare_consolidated_files(
        str(EXAMPLE / "airr" / "airr-olmsted-golden.json"),
        str(decompressed),
    )
    assert match, f"Decompressed gzip output doesn't match golden:\n{message}"


def test_tag_gzip_output(tmp_path):
    """`tag --json-format gzip` produces a valid decompressible output."""
    out = tmp_path / "tagged.json"
    subprocess.run(
        [
            "olmsted", "tag",
            "-i", str(EXAMPLE / "mutations" / "input-olmsted.json"),
            "-o", str(out),
            "--json-format", "gzip", "-q",
        ],
        check=True, capture_output=True,
    )

    gz_out = _gz_output_path(out)
    assert gz_out.exists()
    with gzip.open(gz_out, "rt") as fh:
        data = json.load(fh)
    assert "field_metadata" in data["datasets"][0]


def test_merge_gzip_output(tmp_path):
    """`merge --json-format gzip` produces a valid decompressible output."""
    out = tmp_path / "merged.json"
    subprocess.run(
        [
            "olmsted", "merge",
            "-i", str(EXAMPLE / "merge" / "input-olmsted.json"),
            "--mutations", str(EXAMPLE / "merge" / "input-mutations.csv"),
            "--mutations-use-depth",
            "-o", str(out),
            "--json-format", "gzip", "-q",
        ],
        check=True, capture_output=True,
    )

    gz_out = _gz_output_path(out)
    assert gz_out.exists()
    with gzip.open(gz_out, "rt") as fh:
        data = json.load(fh)
    assert len(data["trees"]) == 2


# --- read side ---


def test_tag_reads_gzip_input(tmp_path):
    """`tag` transparently reads a gzipped Olmsted JSON input."""
    src = EXAMPLE / "mutations" / "input-olmsted.json"
    gz_input = tmp_path / "input.json.gz"
    _gzip_copy(src, gz_input)

    out = tmp_path / "tagged.json"
    result = subprocess.run(
        ["olmsted", "tag", "-i", str(gz_input), "-o", str(out), "-q"],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, f"tag failed on .gz input: {result.stderr}"
    assert out.exists()
    with open(out) as fh:
        data = json.load(fh)
    assert "field_metadata" in data["datasets"][0]


def test_merge_reads_gzip_input(tmp_path):
    """`merge` transparently reads a gzipped Olmsted JSON input."""
    src = EXAMPLE / "merge" / "input-olmsted.json"
    gz_input = tmp_path / "input.json.gz"
    _gzip_copy(src, gz_input)

    out = tmp_path / "merged.json"
    result = subprocess.run(
        [
            "olmsted", "merge",
            "-i", str(gz_input),
            "--mutations", str(EXAMPLE / "merge" / "input-mutations.csv"),
            "--mutations-use-depth",
            "-o", str(out),
            "-q",
        ],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, f"merge failed on .gz input: {result.stderr}"
    assert out.exists()


# --- compressed inputs to `process` -----------------------------------------


def test_process_airr_reads_gz_input(tmp_path):
    """`process -f airr` accepts a gzipped AIRR JSON input (the gz gap that
    motivated this PR)."""
    src = EXAMPLE / "airr" / "input-airr.json"
    gz_input = tmp_path / "input.json.gz"
    _gzip_copy(src, gz_input)

    out = tmp_path / "airr_out.json"
    result = subprocess.run(
        [
            "olmsted", "process", "-f", "airr",
            "-i", str(gz_input),
            "-o", str(out),
            "--seed", "42", "--name", "airr-example",
            "--json-format", "pretty", "-q",
        ],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, f"process -f airr failed on .gz input: {result.stderr}"
    # Can't compare to golden directly — `metadata.source_files` records the
    # temp filename, not the golden's recorded one. Structural sanity instead.
    data = json.loads(out.read_text())
    assert "datasets" in data and "clones" in data and "trees" in data
    assert len(data["trees"]) > 0


def test_process_pcp_reads_gz_clones_input(tmp_path):
    """`process -f pcp -i pcp.csv.gz` accepts a gzipped PCP CSV."""
    src_clones = EXAMPLE / "pcp" / "input-pcp.csv"
    gz_clones = tmp_path / "input.csv.gz"
    _gzip_copy(src_clones, gz_clones)

    out = tmp_path / "pcp_out.json"
    result = subprocess.run(
        [
            "olmsted", "process", "-f", "pcp",
            "-i", str(gz_clones),
            "-t", str(EXAMPLE / "pcp" / "input-trees.csv"),
            "-o", str(out),
            "--seed", "42", "--name", "pcp-example",
            "--json-format", "pretty", "-q",
        ],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, f"process -f pcp failed on .gz clones: {result.stderr}"
    data = json.loads(out.read_text())
    assert "datasets" in data and "clones" in data and "trees" in data
    assert len(data["trees"]) > 0


def test_process_pcp_reads_gz_trees_input(tmp_path):
    """`process -f pcp -t trees.csv.gz` accepts a gzipped trees CSV."""
    src_trees = EXAMPLE / "pcp" / "input-trees.csv"
    gz_trees = tmp_path / "trees.csv.gz"
    _gzip_copy(src_trees, gz_trees)

    out = tmp_path / "pcp_out.json"
    result = subprocess.run(
        [
            "olmsted", "process", "-f", "pcp",
            "-i", str(EXAMPLE / "pcp" / "input-pcp.csv"),
            "-t", str(gz_trees),
            "-o", str(out),
            "--seed", "42", "--name", "pcp-example",
            "--json-format", "pretty", "-q",
        ],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, f"process -f pcp failed on .gz trees: {result.stderr}"
    data = json.loads(out.read_text())
    assert "datasets" in data and "clones" in data and "trees" in data
    assert len(data["trees"]) > 0


# --- determinism (gzip header layer only) ---


def test_gzip_header_is_deterministic(tmp_path):
    """Gzip header (mtime + filename) is pinned, so two writes of identical
    content produce byte-identical files. Useful when piping the same
    structure through gzip twice (the JSON content from process is not yet
    deterministic across runs because of `metadata.created_at`)."""
    from olmsted_cli.data_io import write_olmsted_json

    data = {"hello": "world", "n": [1, 2, 3]}
    a = tmp_path / "a.json"
    b = tmp_path / "b.json"
    write_olmsted_json(data, a, json_format="gzip")
    write_olmsted_json(data, b, json_format="gzip")
    assert (tmp_path / "a.json.gz").read_bytes() == (tmp_path / "b.json.gz").read_bytes()
