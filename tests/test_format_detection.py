"""Tests for format detection including Olmsted JSON."""

import json
import os
import tempfile

import pytest

from olmsted_cli.format_detection import detect_file_format


class TestDetectFileFormat:
    def test_pcp_csv(self):
        assert detect_file_format("example-data/pcp/input-pcp.csv") == "pcp"

    def test_airr_json(self):
        assert detect_file_format("example-data/airr/input-airr.json") == "airr"

    def test_olmsted_json_with_format_tag(self):
        assert detect_file_format("example-data/mutations/input-olmsted.json") == "olmsted"

    def test_olmsted_json_consolidated(self):
        assert detect_file_format("example-data/pcp/pcp-olmsted-golden.json") == "olmsted"

    def test_olmsted_json_without_format_tag(self):
        """Heuristic detection: datasets + metadata keys → olmsted."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            json.dump({
                "metadata": {"schema_version": "2.0.0"},
                "datasets": [{"dataset_id": "test"}],
                "clones": {},
                "trees": [],
            }, f)
            path = f.name

        try:
            assert detect_file_format(path) == "olmsted"
        finally:
            os.unlink(path)

    def test_airr_json_not_misdetected_as_olmsted(self):
        """AIRR JSON (has clones but no datasets/metadata) → airr."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            json.dump({
                "dataset_id": "test",
                "ident": "abc",
                "clones": [],
                "subjects": [],
                "samples": [],
            }, f)
            path = f.name

        try:
            assert detect_file_format(path) == "airr"
        finally:
            os.unlink(path)

    def test_unknown_format(self):
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            f.write("this is not a data file")
            path = f.name

        try:
            assert detect_file_format(path) == "unknown"
        finally:
            os.unlink(path)
