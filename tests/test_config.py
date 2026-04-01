"""Tests for YAML config file support."""

import os
import tempfile

import pytest
import yaml

from olmsted_cli.process_data import load_config


# =============================================================================
# Fixtures
# =============================================================================


@pytest.fixture
def config_dir():
    """Create a temporary directory for config files."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield tmpdir


def write_config(config_dir, config_dict, filename="config.yaml"):
    """Write a config dict as YAML to a temp file."""
    path = os.path.join(config_dir, filename)
    with open(path, "w") as f:
        yaml.dump(config_dict, f)
    return path


# =============================================================================
# Tests: load_config
# =============================================================================


class TestLoadConfig:
    def test_load_valid_config(self, config_dir):
        path = write_config(config_dir, {
            "format": "pcp",
            "name": "Test Dataset",
            "seed": 42,
            "verbose": 2,
            "compute_metrics": True,
            "lbi_tau": 0.05,
        })
        config_dict, custom_fields = load_config(path)
        assert config_dict["format"] == "pcp"
        assert config_dict["name"] == "Test Dataset"
        assert config_dict["seed"] == 42
        assert config_dict["verbose"] == 2
        assert config_dict["compute_metrics"] is True
        assert config_dict["lbi_tau"] == 0.05
        assert custom_fields == []

    def test_missing_config_file(self):
        with pytest.raises(SystemExit):
            load_config("/nonexistent/config.yaml")

    def test_empty_config(self, config_dir):
        path = write_config(config_dir, None)
        config_dict, custom_fields = load_config(path)
        assert config_dict == {}
        assert custom_fields == []

    def test_unrecognized_keys_warn(self, config_dir, capsys):
        path = write_config(config_dir, {
            "format": "pcp",
            "typo_key": "value",
        })
        config_dict, _ = load_config(path)
        captured = capsys.readouterr()
        assert "Unrecognized config key 'typo_key'" in captured.err

    def test_relative_paths_resolved(self, config_dir):
        path = write_config(config_dir, {
            "inputs": ["data.csv"],
            "output": "output/result.json",
            "tree": "trees.csv",
        })
        config_dict, _ = load_config(path)
        assert config_dict["inputs"][0] == os.path.join(config_dir, "data.csv")
        assert config_dict["output"] == os.path.join(config_dir, "output/result.json")
        assert config_dict["tree"] == os.path.join(config_dir, "trees.csv")

    def test_absolute_paths_preserved(self, config_dir):
        path = write_config(config_dir, {
            "inputs": ["/absolute/path/data.csv"],
        })
        config_dict, _ = load_config(path)
        assert config_dict["inputs"][0] == "/absolute/path/data.csv"


# =============================================================================
# Tests: custom_fields parsing
# =============================================================================


class TestCustomFieldsParsing:
    def test_valid_custom_fields(self, config_dir):
        path = write_config(config_dir, {
            "custom_fields": [
                {
                    "name": "my_metric",
                    "level": "clone",
                    "type": "continuous",
                    "label": "My Metric",
                },
                {
                    "name": "my_category",
                    "level": "node",
                    "type": "categorical",
                    "label": "My Category",
                },
            ]
        })
        _, custom_fields = load_config(path)
        assert len(custom_fields) == 2
        assert custom_fields[0]["name"] == "my_metric"
        assert custom_fields[0]["level"] == "clone"
        assert custom_fields[1]["level"] == "node"

    def test_invalid_level_skipped(self, config_dir, capsys):
        path = write_config(config_dir, {
            "custom_fields": [
                {
                    "name": "bad",
                    "level": "invalid_level",
                    "type": "continuous",
                    "label": "Bad",
                },
            ]
        })
        _, custom_fields = load_config(path)
        assert len(custom_fields) == 0
        captured = capsys.readouterr()
        assert "invalid level" in captured.err

    def test_invalid_type_skipped(self, config_dir, capsys):
        path = write_config(config_dir, {
            "custom_fields": [
                {
                    "name": "bad",
                    "level": "clone",
                    "type": "unknown_type",
                    "label": "Bad",
                },
            ]
        })
        _, custom_fields = load_config(path)
        assert len(custom_fields) == 0
        captured = capsys.readouterr()
        assert "invalid type" in captured.err

    def test_missing_required_keys_skipped(self, config_dir, capsys):
        path = write_config(config_dir, {
            "custom_fields": [
                {"name": "incomplete"},  # missing level, type, label
            ]
        })
        _, custom_fields = load_config(path)
        assert len(custom_fields) == 0
        captured = capsys.readouterr()
        assert "missing required keys" in captured.err

    def test_custom_field_with_path(self, config_dir):
        path = write_config(config_dir, {
            "custom_fields": [
                {
                    "name": "surprise_mutsel",
                    "level": "mutation",
                    "type": "continuous",
                    "label": "Surprise Score",
                    "path": "nodes[].mutations[].surprise_mutsel",
                },
            ]
        })
        _, custom_fields = load_config(path)
        assert custom_fields[0]["path"] == "nodes[].mutations[].surprise_mutsel"
