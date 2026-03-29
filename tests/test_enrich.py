"""Tests for the enrich command."""

import json
import os
import tempfile

import pytest
import yaml


@pytest.fixture
def sample_olmsted_json():
    """Create a minimal Olmsted JSON file for testing."""
    data = {
        "metadata": {"format_version": "1.0", "schema_version": "2.0.0"},
        "datasets": [
            {"dataset_id": "test-ds", "name": "Test"}
        ],
        "clones": {
            "test-ds": [
                {
                    "clone_id": "c1",
                    "dataset_id": "test-ds",
                    "unique_seqs_count": 10,
                    "mean_mut_freq": 0.05,
                    "v_call": "IGHV3-48*01",
                    "j_call": "IGHJ4*02",
                    "sample_id": "s1",
                    "v_alignment_start": 0,
                    "v_alignment_end": 294,
                    "j_alignment_start": 300,
                    "j_alignment_end": 350,
                },
            ]
        },
        "trees": [
            {
                "ident": "tree-1",
                "clone_id": "c1",
                "newick": "(a:0.1)root;",
                "nodes": [
                    {
                        "sequence_id": "root",
                        "parent": None,
                        "type": "root",
                        "sequence_alignment": "ATCG",
                        "sequence_alignment_aa": "M",
                        "distance": 0.0,
                        "length": 0.0,
                        "multiplicity": 1,
                    },
                    {
                        "sequence_id": "a",
                        "parent": "root",
                        "type": "leaf",
                        "sequence_alignment": "ATCG",
                        "sequence_alignment_aa": "M",
                        "distance": 0.1,
                        "length": 0.1,
                        "multiplicity": 3,
                    },
                ],
            }
        ],
    }
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".json", delete=False
    ) as f:
        json.dump(data, f)
        return f.name


@pytest.fixture(autouse=True)
def cleanup_files(sample_olmsted_json):
    yield
    if os.path.exists(sample_olmsted_json):
        os.unlink(sample_olmsted_json)


class TestEnrichCommand:
    def test_enrich_adds_field_metadata(self, sample_olmsted_json):
        import subprocess

        with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as out:
            output_path = out.name

        try:
            result = subprocess.run(
                ["olmsted", "enrich", "-i", sample_olmsted_json, "-o", output_path],
                capture_output=True,
                text=True,
            )
            assert result.returncode == 0, f"stderr: {result.stderr}"

            with open(output_path) as f:
                data = json.load(f)

            ds = data["datasets"][0]
            assert "field_metadata" in ds
            fm = ds["field_metadata"]
            assert "clone" in fm
            assert "unique_seqs_count" in fm["clone"]
            assert fm["clone"]["unique_seqs_count"]["type"] == "continuous"
        finally:
            if os.path.exists(output_path):
                os.unlink(output_path)

    def test_enrich_preserves_data(self, sample_olmsted_json):
        import subprocess

        # Read original
        with open(sample_olmsted_json) as f:
            original = json.load(f)

        with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as out:
            output_path = out.name

        try:
            subprocess.run(
                ["olmsted", "enrich", "-i", sample_olmsted_json, "-o", output_path],
                capture_output=True,
                text=True,
            )

            with open(output_path) as f:
                enriched = json.load(f)

            # All original data should still be present
            # Metadata may have "format": "olmsted" added by enrich
            for key in original["metadata"]:
                assert enriched["metadata"][key] == original["metadata"][key]
            assert enriched["metadata"].get("format") == "olmsted"
            assert enriched["clones"] == original["clones"]
            assert enriched["trees"] == original["trees"]
            # Dataset should have field_metadata added but otherwise same
            assert enriched["datasets"][0]["dataset_id"] == original["datasets"][0]["dataset_id"]
        finally:
            if os.path.exists(output_path):
                os.unlink(output_path)

    def test_enrich_in_place(self, sample_olmsted_json):
        import subprocess

        result = subprocess.run(
            ["olmsted", "enrich", "-i", sample_olmsted_json, "--in-place"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0

        with open(sample_olmsted_json) as f:
            data = json.load(f)

        assert "field_metadata" in data["datasets"][0]

    def test_enrich_with_custom_fields(self, sample_olmsted_json):
        import subprocess

        config = {
            "custom_fields": [
                {
                    "name": "custom_score",
                    "level": "clone",
                    "type": "continuous",
                    "label": "Custom Score",
                },
            ]
        }

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".yaml", delete=False
        ) as cf:
            yaml.dump(config, cf)
            config_path = cf.name

        with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as out:
            output_path = out.name

        try:
            result = subprocess.run(
                [
                    "olmsted", "enrich",
                    "-i", sample_olmsted_json,
                    "-o", output_path,
                    "-c", config_path,
                ],
                capture_output=True,
                text=True,
            )
            assert result.returncode == 0, f"stderr: {result.stderr}"

            with open(output_path) as f:
                data = json.load(f)

            fm = data["datasets"][0]["field_metadata"]
            assert "custom_score" in fm["clone"]
            assert fm["clone"]["custom_score"]["label"] == "Custom Score"
        finally:
            for p in (output_path, config_path):
                if os.path.exists(p):
                    os.unlink(p)

    def test_enrich_node_level(self, sample_olmsted_json):
        import subprocess

        with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as out:
            output_path = out.name

        try:
            result = subprocess.run(
                ["olmsted", "enrich", "-i", sample_olmsted_json, "-o", output_path],
                capture_output=True,
                text=True,
            )
            assert result.returncode == 0

            with open(output_path) as f:
                data = json.load(f)

            fm = data["datasets"][0]["field_metadata"]
            # Nodes have multiplicity and distance
            assert "node" in fm
            assert "multiplicity" in fm["node"]
            assert "distance" in fm["node"]
        finally:
            if os.path.exists(output_path):
                os.unlink(output_path)
