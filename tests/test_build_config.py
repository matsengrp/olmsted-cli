"""Tests for the build-config command."""

import json
import os
import subprocess
import tempfile

import pytest


class TestBuildConfigOlmsted:
    def test_olmsted_json(self):
        """build-config on Olmsted JSON produces valid YAML with field entries."""
        result = subprocess.run(
            ["olmsted", "build-config", "-i",
             "example_data/surprise/surprise_subset.json"],
            capture_output=True, text=True,
        )
        assert result.returncode == 0
        output = result.stdout
        assert "custom_fields:" in output
        assert "clone" in output
        assert "mutation" in output
        assert "surprise_mutsel" in output

    def test_output_to_file(self):
        """build-config -o writes to file instead of stdout."""
        with tempfile.NamedTemporaryFile(suffix=".yaml", delete=False) as f:
            out_path = f.name

        try:
            result = subprocess.run(
                ["olmsted", "build-config", "-i",
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

    def test_shows_ranges_for_continuous_mutation(self):
        """Continuous mutation fields show range comments."""
        result = subprocess.run(
            ["olmsted", "build-config", "-i",
             "example_data/surprise/surprise_subset.json"],
            capture_output=True, text=True,
        )
        assert "range in data:" in result.stdout

    def test_includes_processing_options_template(self):
        """Output includes commented-out processing options."""
        result = subprocess.run(
            ["olmsted", "build-config", "-i",
             "example_data/surprise/surprise_subset.json"],
            capture_output=True, text=True,
        )
        assert "Processing Options" in result.stdout

    def test_includes_alias_reference(self):
        """Output includes cross-format alias reference."""
        result = subprocess.run(
            ["olmsted", "build-config", "-i",
             "example_data/surprise/surprise_subset.json"],
            capture_output=True, text=True,
        )
        assert "output_name" in result.stdout
        assert "v_call" in result.stdout

    def test_includes_output_name_example(self):
        """Output includes output_name usage example."""
        result = subprocess.run(
            ["olmsted", "build-config", "-i",
             "example_data/surprise/surprise_subset.json"],
            capture_output=True, text=True,
        )
        assert "rearrangement_count" in result.stdout
        assert "output_name: unique_seqs_count" in result.stdout


class TestBuildConfigPcp:
    def test_pcp_with_trees(self):
        """build-config on raw PCP CSV + tree CSV discovers fields."""
        result = subprocess.run(
            ["olmsted", "build-config", "-i",
             "example_data/pcp/pcp.csv",
             "-t", "example_data/pcp/trees.csv"],
            capture_output=True, text=True,
        )
        assert result.returncode == 0
        output = result.stdout
        assert "custom_fields:" in output
        assert "Clone level" in output
        assert "unique_seqs_count" in output

    def test_pcp_extra_columns(self):
        """build-config on PCP with extra columns discovers them."""
        result = subprocess.run(
            ["olmsted", "build-config", "-i",
             "example_data/test-fields/pcp-test-fields.csv",
             "-t", "example_data/test-fields/trees-test-fields.csv"],
            capture_output=True, text=True,
        )
        assert result.returncode == 0
        output = result.stdout
        assert "foobar_score" in output
        assert "foobar_category" in output
        assert "foobar_weight" in output
        assert "foobar_class" in output

    def test_pcp_shows_compute_metrics(self):
        """PCP config template includes compute_metrics option."""
        result = subprocess.run(
            ["olmsted", "build-config", "-i",
             "example_data/pcp/pcp.csv",
             "-t", "example_data/pcp/trees.csv"],
            capture_output=True, text=True,
        )
        assert "compute_metrics" in result.stdout


class TestBuildConfigAirr:
    def test_airr(self):
        """build-config on AIRR JSON discovers fields."""
        result = subprocess.run(
            ["olmsted", "build-config", "-i",
             "example_data/airr/airr.json"],
            capture_output=True, text=True,
        )
        assert result.returncode == 0
        output = result.stdout
        assert "custom_fields:" in output
        assert "v_call" in output or "V Gene" in output


class TestBuildConfigFormatDetection:
    def test_detects_olmsted(self):
        result = subprocess.run(
            ["olmsted", "build-config", "-i",
             "example_data/surprise/surprise_subset.json"],
            capture_output=True, text=True,
        )
        assert "OLMSTED" in result.stderr

    def test_detects_pcp(self):
        result = subprocess.run(
            ["olmsted", "build-config", "-i",
             "example_data/pcp/pcp.csv"],
            capture_output=True, text=True,
        )
        assert "PCP" in result.stderr

    def test_detects_airr(self):
        result = subprocess.run(
            ["olmsted", "build-config", "-i",
             "example_data/airr/airr.json"],
            capture_output=True, text=True,
        )
        assert "AIRR" in result.stderr
