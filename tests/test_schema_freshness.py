"""
Test that the committed JSON schema file matches its YAML source.

Run `make schemas` to regenerate the JSON file after editing the YAML.
"""

import json
from pathlib import Path

import yaml

SCHEMA_DIR = Path(__file__).parent.parent / "olmsted_cli" / "schemas"
YAML_PATH = SCHEMA_DIR / "olmsted-schema.yaml"
JSON_PATH = SCHEMA_DIR / "olmsted-schema.json"


def test_json_matches_yaml():
    """Committed olmsted-schema.json must match what would be generated from the YAML source."""
    assert YAML_PATH.exists(), f"Missing YAML schema: {YAML_PATH}"
    assert JSON_PATH.exists(), (
        f"Missing JSON schema: {JSON_PATH.name}. Run `make schemas` to generate it."
    )

    yaml_data = yaml.safe_load(YAML_PATH.read_text())
    json_data = json.loads(JSON_PATH.read_text())

    assert yaml_data == json_data, (
        f"{JSON_PATH.name} is stale. Run `make schemas` to regenerate it from "
        f"{YAML_PATH.name}."
    )
