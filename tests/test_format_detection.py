"""Tests for format detection including Olmsted JSON."""

import json
import os
import tempfile

from olmsted_cli.data_io import detect_file_format


class TestDetectFileFormat:
    def test_pcp_csv(self):
        assert detect_file_format("example-data/pcp/input-pcp.csv") == "pcp"

    def test_olmsted_json_with_format_tag(self):
        assert (
            detect_file_format("example-data/mutations/input-olmsted.json") == "olmsted"
        )

    def test_olmsted_json_consolidated(self):
        assert (
            detect_file_format("example-data/pcp/pcp-olmsted-golden.json") == "olmsted"
        )

    def test_olmsted_json_without_format_tag(self):
        """Heuristic detection: datasets + metadata keys → olmsted."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            json.dump(
                {
                    "metadata": {"schema_version": "2.0.0"},
                    "datasets": [{"dataset_id": "test"}],
                    "clones": {},
                    "trees": [],
                },
                f,
            )
            path = f.name

        try:
            assert detect_file_format(path) == "olmsted"
        finally:
            os.unlink(path)

    def test_airr_json_variants(self):
        """AIRR-C v2 Clone/Tree (top-level Clone + Rearrangement) → airr."""
        for variant in (
            "nocell-noinfo",
            "unpaired-noinfo",
            "paired-noinfo",
            "nocell-info",
            "unpaired-info",
            "paired-info",
        ):
            assert (
                detect_file_format(f"example-data/airr/input-{variant}.json") == "airr"
            )

    def test_minimal_clone_rearrangement_detects_as_airr(self):
        """A minimal {Clone, Rearrangement} object → airr."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            json.dump({"Clone": [], "Rearrangement": []}, f)
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
