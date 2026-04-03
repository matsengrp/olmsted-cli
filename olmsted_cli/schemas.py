#!/usr/bin/env python
"""
Unified schema definitions for Olmsted data structures.

Schemas are authored in YAML files under olmsted_cli/schemas/ and loaded at
import time. This makes the YAML the single authoritative source: the Python
dicts used for validation are identical to the published JSON artifacts.

Dynamic fragments (field_metadata.properties, field type/display enums) are
patched in after YAML loading because they are derived from constants.py.

NOTE: The AIRR schema components reference the official AIRR schema from
airr-standards/specs/airr-schema.yaml. The SCHEMA_VERSION constant corresponds
to the 'version' field in the Info section of that schema.
"""

from pathlib import Path

import yaml

from .constants import DISPLAY_MODES, FIELD_LEVELS, FIELD_TYPES

# Version Constants
# SCHEMA_VERSION corresponds to Info.version in airr-standards/specs/airr-schema.yaml
SCHEMA_VERSION = "2.0.0"

# Output display modes (skip means "not in output", so exclude it from schema)
_OUTPUT_DISPLAY_MODES = sorted(DISPLAY_MODES - {"skip"})

# Schema fragment for a single field_metadata entry.
# Built from constants so that adding a new FIELD_TYPE automatically updates validation.
_FIELD_ENTRY_SCHEMA = {
    "type": "object",
    "properties": {
        "type": {
            "type": "string",
            "enum": sorted(FIELD_TYPES),
        },
        "display": {
            "type": "string",
            "enum": _OUTPUT_DISPLAY_MODES,
        },
        "label": {"type": "string"},
        "range": {
            "type": "array",
            "items": {"type": "number"},
            "minItems": 2,
            "maxItems": 2,
        },
    },
    "required": ["type", "label"],
}

_SCHEMA_DIR = Path(__file__).parent / "schemas"


def _load(name):
    """Load a YAML schema file from the schemas directory."""
    return yaml.safe_load((_SCHEMA_DIR / name).read_text())


node_spec = _load("node.schema.yaml")
tree_spec = _load("tree.schema.yaml")
clone_spec = _load("clone.schema.yaml")
dataset_spec = _load("dataset.schema.yaml")

# Patch field_metadata.properties dynamically from FIELD_LEVELS so that adding
# a new level in constants.py automatically extends the schema.
dataset_spec["properties"]["field_metadata"]["properties"] = {
    level: {
        "description": f"{level.title()}-level field metadata",
        "type": "object",
        "additionalProperties": _FIELD_ENTRY_SCHEMA,
    }
    for level in sorted(FIELD_LEVELS)
}


# Helper functions for AIRR-specific schema generation
def id_spec(description="Identifier"):
    """Create an ID specification with custom description."""
    return {
        "description": description,
        "type": "string",
    }


def sequence_spec(description):
    """Create a sequence specification with custom description."""
    return {
        "description": description,
        "type": "string",
    }


def multiplicity_spec(description=None):
    """Create a multiplicity specification for AIRR fields."""
    return {
        "description": description
        or "Number of times sequence was observed in the sample. The presence of a given sequence in a clonal family may represent many identical such sequences in the original sample.",
        "type": ["integer", "null"],
        "minimum": 0,
    }


# AIRR-specific schemas (used by legacy code paths)
ident_spec = {
    "description": "UUID specific to the given object",
    "type": "string",
}

build_spec = {
    "description": "Information about how a dataset was built",
    "type": "object",
    "required": ["commit"],
    "title": "Build info",
    "properties": {
        "commit": {
            "description": "Commit sha of whatever build system you used to process the data",
            "type": "string",
        },
        "time": {
            "description": "Time at which build was initiated",
            "type": "string",
        },
    },
}

timepoint_multiplicity_spec = {
    "title": "Timepoint multiplicity",
    "description": "Multiplicity at a specific time",
    "type": "object",
    "required": ["multiplicity", "timepoint_id"],
    "properties": {
        "timepoint_id": {
            "description": "Id associated with the timepoint in question",
            "type": ["string", "null"],
        },
        "multiplicity": {
            "description": "Number of times sequence was observed at the given timepoint",
            "type": "integer",
            "minimum": 0,
        },
    },
    "additionalProperties": False,
}

sample_spec = {
    "title": "Sample",
    "description": "A sample is a collection of sequences",
    "type": "object",
    "required": ["locus"],
    "properties": {
        "ident": ident_spec,
        "sample_id": {
            "description": "Sample id",
            "type": "string",
        },
        "timepoint_id": {
            "description": 'Timepoint associated with this sample (may choose "merged" if data has been combined from multiple timepoints)',
            "type": "string",
        },
        "locus": {
            "description": "B-cell Locus",
            "type": "string",
        },
    },
}

subject_spec = {
    "title": "Subject",
    "description": "Subject from which the clonal family was sampled",
    "type": "object",
    "required": ["subject_id"],
    "properties": {
        "ident": ident_spec,
        "subject_id": {
            "description": "Subject id",
            "type": "string",
        },
    },
}

seed_spec = {
    "title": "Seed",
    "description": "A sequence of interest among other clonal family members",
    "type": ["object", "null"],
    "required": ["seed_id"],
    "properties": {
        "ident": ident_spec,
        "seed_id": {
            "description": "Seed id",
            "type": "string",
        },
    },
}

# Export all schemas, constants, and helper functions
__all__ = [
    # Version constants
    "SCHEMA_VERSION",
    # Main unified schemas
    "node_spec",
    "tree_spec",
    "clone_spec",
    "dataset_spec",
    # AIRR-specific schemas
    "ident_spec",
    "build_spec",
    "timepoint_multiplicity_spec",
    "sample_spec",
    "subject_spec",
    "seed_spec",
    # Helper functions
    "id_spec",
    "sequence_spec",
    "multiplicity_spec",
]
