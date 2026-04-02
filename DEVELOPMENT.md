# Development Guide

This guide covers development setup, testing, and contributing to olmsted-cli.

See also:
- **[README.md](./README.md)**: User documentation and examples
- **[ARCHITECTURE.md](./ARCHITECTURE.md)**: Data flow, processing pipelines, field_metadata system
- **[FORMATS.md](./FORMATS.md)**: Input/output format specs, field mapping, validation rules
- **[CLAUDE.md](./CLAUDE.md)**: Code quality rules, terminology, quick reference

## Table of Contents

- [Quick Start](#quick-start)
- [Prerequisites](#prerequisites)
- [Running Tests](#running-tests)
- [Project Structure](#project-structure)
- [Common Tasks](#common-tasks)
- [Linting](#linting)
- [Release](#release)

---

## Quick Start

```bash
git clone https://github.com/matsengrp/olmsted-cli.git
cd olmsted-cli
pip install -e ".[dev]"
pytest
```

## Prerequisites

- **Python**: 3.8+
- **pip**: For package installation
- **Dependencies**: Automatically installed (ete3, jsonschema, numpy, pyyaml, scipy, etc.)

## Running Tests

```bash
# Run all tests
pytest

# Verbose output
pytest -v

# Specific test file
pytest tests/test_field_metadata.py

# Specific test
pytest tests/test_field_metadata.py::TestInferFieldType::test_aa_values

# Pattern matching
pytest -k "pcp"

# Coverage
pytest --cov=olmsted_cli
```

### Test Organization

| File | What it tests |
|------|---------------|
| `test_cli_processing.py` | Full CLI integration (subprocess calls, golden data comparison) |
| `test_field_metadata.py` | Field registries, type inference, AA/DNA detection, ranges, custom fields |
| `test_build_config.py` | Config generation from all formats, alias suggestions, skip section |
| `test_config.py` | YAML config loading, custom_fields validation, path resolution |
| `test_tag.py` | Tag command (add metadata, preserve data, in-place, custom fields) |
| `test_pcp_extras.py` | Extra CSV columns, chain partitioning, column aliases, coercion |
| `test_format_detection.py` | Olmsted/AIRR/PCP format detection |
| `test_validation.py` | Schema validation for datasets, clones, trees |

### Golden Data

Integration tests compare CLI output against golden data in `example_data/`. Regenerate after output format changes:

```bash
olmsted process -f airr -i example_data/airr/airr.json \
  --split-files example_data/airr/split_golden_data --seed 42 --name airr-example -q
olmsted process -f airr -i example_data/airr/airr.json \
  -o example_data/airr/consolidated_golden_data.json --seed 42 --name airr-example -q
olmsted process -f pcp -i example_data/pcp/pcp.csv -t example_data/pcp/trees.csv \
  --split-files example_data/pcp/split_golden_data --seed 42 --name pcp-example -q
olmsted process -f pcp -i example_data/pcp/pcp.csv -t example_data/pcp/trees.csv \
  -o example_data/pcp/consolidated_golden_data.json --seed 42 --name pcp-example -q
```

## Project Structure

```
olmsted-cli/
├── olmsted_cli/              # Main package
│   ├── cli.py                # CLI entry point
│   ├── constants.py          # All constants, enums, registries
│   ├── types.py              # TypedDict definitions
│   ├── schemas.py            # JSON Schema definitions
│   ├── metrics.py            # LBI, LBR, scaled_affinity
│   ├── field_metadata.py     # Field metadata generation
│   ├── process_data.py       # Unified processor + YAML config
│   ├── process_pcp_data.py   # PCP CSV parsing and conversion
│   ├── process_airr_data.py  # AIRR JSON processing
│   ├── process_utils.py      # Shared utilities, output writing
│   ├── build_config.py       # build-config command
│   ├── tag.py                # tag command
│   ├── api.py                # Programmatic API
│   └── configs/              # Default YAML configs
├── tests/                    # Test suite
├── example_data/             # Example datasets and golden data
│   ├── airr/                 # AIRR format + golden data
│   ├── pcp/                  # PCP format + golden data
│   ├── pcp-paired/           # Paired heavy/light PCP
│   ├── surprise/             # Pre-built Olmsted JSON with mutation-level data
│   └── test-fields/          # Foobar test data (all types × levels)
├── ARCHITECTURE.md           # System architecture
├── DEVELOPMENT.md            # This file
├── CLAUDE.md                 # AI assistant guidance
└── README.md                 # User documentation
```

## Common Tasks

### Adding a New CLI Command

1. Create `olmsted_cli/mycommand.py` with `get_args()` and `main()`
2. Register in `cli.py`: import, add subparser, add routing, add wrapper function

### Adding a New Known Field

```python
# In constants.py:
KNOWN_CLONE_FIELDS["my_field"] = {"type": "continuous", "label": "My Field"}
```

### Adding a New Field Type

1. Add to `FIELD_TYPES` in `constants.py`
2. Schema auto-updates (generated from `FIELD_TYPES`)
3. Update `infer_field_type()` in `field_metadata.py` if auto-detection needed
4. Update type docs in `build_config.py`

### Adding a New Level

1. Add to `FIELD_LEVELS` in `constants.py`
2. If alias, add to `LEVEL_ALIASES`
3. Add generator in `field_metadata.py`
4. Add section in `build_config.py`'s `_build_yaml()`

### Modifying PCP Column Handling

- **Known columns**: `KNOWN_PCP_COLUMNS` in `constants.py`
- **Column aliases**: `CHAIN_COLUMN_ALIASES` in `constants.py`
- **Chain partitioning**: `_partition_chain_fields()` in `process_pcp_data.py`

## Linting

```bash
ruff check .
ruff format .
```

## Release

Version is in `pyproject.toml`.

```bash
pip install build twine
python -m build
twine upload dist/*
```

---

_Last updated: 2026-03-29_
