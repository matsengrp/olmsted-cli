"""
Test that committed JSON schema files match their YAML sources.

Run `make schemas` to regenerate JSON files after editing YAML.
"""

import json
from pathlib import Path

import pytest
import yaml

SCHEMA_DIR = Path(__file__).parent.parent / "olmsted_cli" / "schemas"
YAML_FILES = sorted(SCHEMA_DIR.glob("*.schema.yaml"))
assert YAML_FILES, f"No *.schema.yaml files found in {SCHEMA_DIR} — check path"


@pytest.mark.parametrize("yaml_path", YAML_FILES, ids=lambda p: p.name)
def test_json_matches_yaml(yaml_path):
    """Committed JSON must match what would be generated from the YAML source."""
    json_path = yaml_path.with_suffix(".json")
    assert json_path.exists(), (
        f"Missing JSON schema: {json_path.name}. Run `make schemas` to generate it."
    )

    yaml_data = yaml.safe_load(yaml_path.read_text())
    json_data = json.loads(json_path.read_text())

    assert yaml_data == json_data, (
        f"{json_path.name} is stale. Run `make schemas` to regenerate it from "
        f"{yaml_path.name}."
    )
