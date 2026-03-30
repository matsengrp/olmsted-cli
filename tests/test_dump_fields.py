"""Tests for the dump-fields command."""

import json
import os
import subprocess
import tempfile

import pytest


class TestDumpFieldsOlmsted:
    def test_dump_olmsted_json(self):
        """dump-fields on Olmsted JSON produces valid YAML with field entries."""
        result = subprocess.run(
            ["olmsted", "dump-fields", "-i",
             "example_data/surprise/surprise_subset.json"],
            capture_output=True, text=True,
        )
        assert result.returncode == 0
        output = result.stdout
        assert "custom_fields:" in output
        assert "clone" in output
        assert "mutation" in output
        assert "surprise_mutsel" in output

    def test_dump_to_file(self):
        """dump-fields -o writes to file instead of stdout."""
        with tempfile.NamedTemporaryFile(suffix=".yaml", delete=False) as f:
            out_path = f.name

        try:
            result = subprocess.run(
                ["olmsted", "dump-fields", "-i",
                 "example_data/surprise/surprise_subset.json",
                 "-o", out_path],
                capture_output=True, text=True,
            )
            assert result.returncode == 0
            assert os.path.exists(out_path)
            with open(out_path) as f:
                content = f.read()
            assert "custom_fields:" in content
        finally:
            if os.path.exists(out_path):
                os.unlink(out_path)

    def test_dump_shows_ranges_for_continuous_mutation(self):
        """Continuous mutation fields show range comments."""
        result = subprocess.run(
            ["olmsted", "dump-fields", "-i",
             "example_data/surprise/surprise_subset.json"],
            capture_output=True, text=True,
        )
        assert "range in data:" in result.stdout


class TestDumpFieldsPcp:
    def test_dump_pcp_with_trees(self):
        """dump-fields on raw PCP CSV + tree CSV discovers fields."""
        result = subprocess.run(
            ["olmsted", "dump-fields", "-i",
             "example_data/pcp/pcp.csv",
             "-t", "example_data/pcp/trees.csv"],
            capture_output=True, text=True,
        )
        assert result.returncode == 0
        output = result.stdout
        assert "custom_fields:" in output
        assert "Clone level" in output
        assert "unique_seqs_count" in output

    def test_dump_pcp_extra_columns(self):
        """dump-fields on PCP with extra columns discovers them."""
        result = subprocess.run(
            ["olmsted", "dump-fields", "-i",
             "example_data/test-fields/pcp-test-fields.csv",
             "-t", "example_data/test-fields/trees-test-fields.csv"],
            capture_output=True, text=True,
        )
        assert result.returncode == 0
        output = result.stdout
        # Extra tree CSV columns → clone level
        assert "foobar_score" in output
        assert "foobar_category" in output
        # Extra PCP CSV columns → node level
        assert "foobar_weight" in output
        assert "foobar_class" in output


class TestDumpFieldsAirr:
    def test_dump_airr(self):
        """dump-fields on AIRR JSON discovers fields."""
        result = subprocess.run(
            ["olmsted", "dump-fields", "-i",
             "example_data/airr/airr.json"],
            capture_output=True, text=True,
        )
        assert result.returncode == 0
        output = result.stdout
        assert "custom_fields:" in output
        assert "v_call" in output or "V Gene" in output


class TestDumpFieldsFormatDetection:
    def test_detects_olmsted(self):
        """Correctly detects Olmsted JSON format."""
        result = subprocess.run(
            ["olmsted", "dump-fields", "-i",
             "example_data/surprise/surprise_subset.json"],
            capture_output=True, text=True,
        )
        assert "OLMSTED" in result.stderr

    def test_detects_pcp(self):
        """Correctly detects PCP format."""
        result = subprocess.run(
            ["olmsted", "dump-fields", "-i",
             "example_data/pcp/pcp.csv"],
            capture_output=True, text=True,
        )
        assert "PCP" in result.stderr

    def test_detects_airr(self):
        """Correctly detects AIRR format."""
        result = subprocess.run(
            ["olmsted", "dump-fields", "-i",
             "example_data/airr/airr.json"],
            capture_output=True, text=True,
        )
        assert "AIRR" in result.stderr
