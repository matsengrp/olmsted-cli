"""Unit tests for the data_io module."""

import gzip
import json
from pathlib import Path

import pytest
import yaml

from olmsted_cli.constants import (
    FORMAT_AIRR,
    FORMAT_OLMSTED,
    FORMAT_PCP,
)
from olmsted_cli.data_io import (
    open_file,
    read_airr_json,
    read_csv_rows,
    read_olmsted_json,
    read_pcp_csv_rows,
    read_yaml_config,
    write_olmsted_json,
    write_file,
)

REPO_ROOT = Path(__file__).parent.parent
EXAMPLE = REPO_ROOT / "example-data"


# --- open_file -------------------------------------------------------------


def test_open_file_detects_olmsted():
    handle, fmt = open_file(EXAMPLE / "mutations" / "input-olmsted.json")
    handle.close()
    assert fmt == FORMAT_OLMSTED


def test_open_file_detects_airr():
    handle, fmt = open_file(EXAMPLE / "airr" / "input-airr.json")
    handle.close()
    assert fmt == FORMAT_AIRR


def test_open_file_detects_pcp():
    handle, fmt = open_file(EXAMPLE / "pcp" / "input-pcp.csv")
    handle.close()
    assert fmt == FORMAT_PCP


def test_open_file_handles_gz(tmp_path):
    """A gzipped Olmsted file detects as olmsted and opens for reading."""
    src = EXAMPLE / "mutations" / "input-olmsted.json"
    gz_path = tmp_path / "input.json.gz"
    with open(src, "rb") as src_fh, gzip.open(gz_path, "wb") as gz_fh:
        gz_fh.write(src_fh.read())

    handle, fmt = open_file(gz_path)
    try:
        # Confirm we can actually read decompressed content
        data = json.load(handle)
    finally:
        handle.close()
    assert fmt == FORMAT_OLMSTED
    assert "datasets" in data


def test_open_file_returns_unknown_when_no_expected(tmp_path):
    """Without expected_formats, an unrecognized file opens fine and returns
    'unknown' as the detected format. Callers that don't care about format
    (e.g., the validate command's per-record checks) rely on this."""
    bogus = tmp_path / "bogus.txt"
    bogus.write_text("not a recognized format")
    handle, fmt = open_file(bogus)
    handle.close()
    assert fmt == "unknown"


def test_open_file_rejects_unknown_when_expected_set(tmp_path):
    """With expected_formats given, unknown is treated as a mismatch."""
    bogus = tmp_path / "bogus.txt"
    bogus.write_text("not a recognized format")
    with pytest.raises(ValueError, match="detected 'unknown'"):
        open_file(bogus, expected_formats=("airr",))


def test_open_file_rejects_expected_mismatch():
    """Asking for olmsted on an airr file fails fast."""
    with pytest.raises(ValueError, match="Expected.*olmsted.*detected 'airr'"):
        open_file(EXAMPLE / "airr" / "input-airr.json", expected_formats=(FORMAT_OLMSTED,))


def test_open_file_accepts_when_in_expected_set():
    """Multi-format expected_formats works."""
    handle, fmt = open_file(
        EXAMPLE / "airr" / "input-airr.json",
        expected_formats=(FORMAT_AIRR, FORMAT_PCP),
    )
    handle.close()
    assert fmt == FORMAT_AIRR


# --- read_olmsted_json ------------------------------------------------------


def test_read_olmsted_json_happy():
    data = read_olmsted_json(EXAMPLE / "mutations" / "input-olmsted.json")
    assert "datasets" in data and "clones" in data and "trees" in data


def test_read_olmsted_json_rejects_airr_file():
    """Passing an AIRR file is a format mismatch — detected as airr, not olmsted."""
    with pytest.raises(ValueError, match="Expected.*olmsted"):
        read_olmsted_json(EXAMPLE / "airr" / "input-airr.json")


def test_read_olmsted_json_rejects_malformed_json(tmp_path):
    """A malformed JSON file fails at format detection (it can't be parsed
    well enough to find the format tag), so the user-facing message is
    'Could not infer format' rather than a JSON parse error. Either
    rejection is acceptable; we just want a fail-fast ValueError."""
    bad = tmp_path / "bad-olmsted.json"
    bad.write_text('{"metadata": {"format": "olmsted"}, "datasets": [BROKEN')
    with pytest.raises(ValueError):
        read_olmsted_json(bad)


def test_read_olmsted_json_rejects_missing_required_keys(tmp_path):
    """Detected as olmsted (format tag) but missing the required top-level keys."""
    skeletal = tmp_path / "skeletal.json"
    skeletal.write_text('{"metadata": {"format": "olmsted"}, "datasets": []}')
    with pytest.raises(ValueError, match="missing required Olmsted top-level keys"):
        read_olmsted_json(skeletal)


# --- read_airr_json ---------------------------------------------------------


def test_read_airr_json_happy():
    data = read_airr_json(EXAMPLE / "airr" / "input-airr.json")
    # AIRR file has clones; structural validation is the caller's job.
    assert "clones" in data


def test_read_airr_json_rejects_olmsted_file():
    with pytest.raises(ValueError, match="Expected.*airr"):
        read_airr_json(EXAMPLE / "mutations" / "input-olmsted.json")


# --- read_pcp_csv_rows ------------------------------------------------------


def test_read_pcp_csv_rows_yields_dicts():
    rows = list(read_pcp_csv_rows(EXAMPLE / "pcp" / "input-pcp.csv"))
    assert rows, "expected at least one row from PCP CSV"
    assert all(isinstance(r, dict) for r in rows)


def test_read_pcp_csv_rows_handles_gz(tmp_path):
    src = EXAMPLE / "pcp" / "input-pcp.csv"
    gz_path = tmp_path / "input.csv.gz"
    with open(src, "rb") as src_fh, gzip.open(gz_path, "wb") as gz_fh:
        gz_fh.write(src_fh.read())

    rows = list(read_pcp_csv_rows(gz_path))
    assert rows


def test_read_pcp_csv_rows_rejects_non_pcp_file():
    """A JSON file isn't PCP — should fail fast."""
    with pytest.raises(ValueError, match="Expected.*pcp"):
        list(read_pcp_csv_rows(EXAMPLE / "airr" / "input-airr.json"))


# --- read_csv_rows ----------------------------------------------------------


def test_read_csv_rows_iterates(tmp_path):
    """Generic CSV reader works on a hand-rolled CSV without format detection."""
    csv_path = tmp_path / "data.csv"
    csv_path.write_text("a,b,c\n1,2,3\n4,5,6\n")
    rows = list(read_csv_rows(csv_path))
    assert rows == [{"a": "1", "b": "2", "c": "3"}, {"a": "4", "b": "5", "c": "6"}]


def test_read_csv_rows_handles_gz(tmp_path):
    csv_path = tmp_path / "data.csv.gz"
    with gzip.open(csv_path, "wt") as fh:
        fh.write("a,b\n1,2\n")
    rows = list(read_csv_rows(csv_path))
    assert rows == [{"a": "1", "b": "2"}]


# --- read_yaml_config -------------------------------------------------------


def test_read_yaml_config_parses(tmp_path):
    cfg_path = tmp_path / "config.yaml"
    cfg_path.write_text("name: test\ncount: 3\n")
    cfg = read_yaml_config(cfg_path)
    assert cfg == {"name": "test", "count": 3}


def test_read_yaml_config_handles_gz(tmp_path):
    cfg_path = tmp_path / "config.yaml.gz"
    with gzip.open(cfg_path, "wt") as fh:
        yaml.safe_dump({"foo": "bar"}, fh)
    cfg = read_yaml_config(cfg_path)
    assert cfg == {"foo": "bar"}


# --- write_file dispatcher ------------------------------------------------


def test_write_file_routes_olmsted_json(tmp_path):
    out = tmp_path / "out.json"
    written = write_file({"datasets": [], "clones": {}, "trees": []}, out)
    assert written == str(out)
    assert json.loads(out.read_text()) == {"datasets": [], "clones": {}, "trees": []}


def test_write_file_passes_through_opts_to_write_olmsted_json(tmp_path):
    out = tmp_path / "out.json"
    written = write_file({"x": 1}, out, json_format="gzip")
    assert written.endswith(".gz")
    assert Path(written).exists()


def test_write_file_rejects_unknown_kind(tmp_path):
    with pytest.raises(ValueError, match="Unknown output_kind"):
        write_file({}, tmp_path / "x.json", output_kind="not_a_real_kind")


# --- write_olmsted_json (smoke; details covered by test_gzip_io) -----------


def test_write_olmsted_json_pretty_round_trip(tmp_path):
    out = tmp_path / "out.json"
    payload = {"hello": "world", "n": [1, 2, 3]}
    write_olmsted_json(payload, out, json_format="pretty")
    assert json.loads(out.read_text()) == payload
