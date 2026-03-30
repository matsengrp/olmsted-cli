# Development Guide

This guide covers development setup, common tasks, and contributing to olmsted-cli.

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

# Run with verbose output
pytest -v

# Run a specific test file
pytest tests/test_field_metadata.py

# Run a specific test
pytest tests/test_field_metadata.py::TestInferFieldType::test_aa_values

# Run tests matching a pattern
pytest -k "pcp"

# Run with coverage
pytest --cov=olmsted_cli
```

### Test Organization

| File | What it tests |
|------|---------------|
| `test_cli_processing.py` | Full CLI integration (subprocess calls, golden data comparison) |
| `test_field_metadata.py` | Field registries, type inference, AA/DNA detection, ranges, custom fields |
| `test_build_config.py` | Config generation from all formats, alias suggestions, skip section |
| `test_config.py` | YAML config loading, custom_fields validation, path resolution |
| `test_enrich.py` | Enrich command (add metadata, preserve data, in-place, custom fields) |
| `test_pcp_extras.py` | Extra CSV columns, chain partitioning, column aliases, coercion |
| `test_format_detection.py` | Olmsted/AIRR/PCP format detection |
| `test_validation.py` | Schema validation for datasets, clones, trees |

### Golden Data

Integration tests compare CLI output against golden data in `example_data/`. When the output format changes (e.g., adding `field_metadata`), regenerate:

```bash
# AIRR golden data
olmsted process -f airr -i example_data/airr/airr.json \
  --split-files example_data/airr/split_golden_data --seed 42 --name airr-example -q
olmsted process -f airr -i example_data/airr/airr.json \
  -o example_data/airr/consolidated_golden_data.json --seed 42 --name airr-example -q

# PCP golden data
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
│   ├── metrics.py            # LBI, LBR, scaled_affinity computation
│   ├── field_metadata.py     # Field metadata generation
│   ├── process_data.py       # Unified processor + YAML config support
│   ├── process_pcp_data.py   # PCP CSV parsing and conversion
│   ├── process_airr_data.py  # AIRR JSON processing
│   ├── process_utils.py      # Shared utilities, output writing
│   ├── build_config.py       # build-config command
│   ├── enrich.py             # enrich command
│   ├── api.py                # Programmatic API (OlmstedData class)
│   ├── validate.py           # validate command
│   ├── summary.py            # summary command
│   ├── split.py              # split command
│   └── configs/              # Default YAML configs
├── tests/                    # Test suite
├── example_data/             # Example datasets and golden data
│   ├── airr/                 # AIRR format examples + golden data
│   ├── pcp/                  # PCP format examples + golden data
│   ├── pcp-paired/           # Paired heavy/light chain PCP
│   ├── pcp-light/            # Light chain only PCP
│   ├── pcp-byhand/           # Artificial test data
│   ├── surprise/             # DASM2 surprise analysis subset
│   └── test-fields/          # Foobar test data for all field types
├── pyproject.toml            # Package configuration
├── CLAUDE.md                 # Claude Code guidance
├── DEVELOPMENT.md            # This file
└── README.md                 # User documentation
```

## Common Development Tasks

### Adding a New CLI Command

1. Create `olmsted_cli/mycommand.py` with `get_args()` and `main()` functions
2. Register in `olmsted_cli/cli.py`:
   - Import the module
   - Add subparser entry
   - Add command routing in the `if/elif` chain
   - Add a `mycommand_command()` wrapper function

### Adding a New Known Field

Add to the appropriate registry in `constants.py`:

```python
KNOWN_CLONE_FIELDS = {
    ...
    "my_new_field": {"type": "continuous", "label": "My New Field"},
}
```

### Adding a New Field Type

1. Add to `FIELD_TYPES` in `constants.py`
2. Schema enum auto-updates (generated from `FIELD_TYPES` in `schemas.py`)
3. Update `infer_field_type()` in `field_metadata.py` if auto-detection needed
4. Update build-config type docs in `build_config.py`

### Adding a New Level

1. Add to `FIELD_LEVELS` in `constants.py`
2. If alias needed, add to `LEVEL_ALIASES`
3. Add generator function in `field_metadata.py`
4. Add section in `build_config.py`'s `_build_yaml()`
5. Schema auto-updates (generated from `FIELD_LEVELS`)

### Adding a Suggested Skip or Type Override

```python
# In constants.py:

# Field should be skipped by default in build-config:
SUGGESTED_SKIP_FIELDS.add("my_internal_field")

# Field gets wrong type from auto-inference:
SUGGESTED_FIELD_TYPES["my_field"] = "tooltip"
```

### Modifying PCP Column Handling

- **Known columns**: `KNOWN_PCP_COLUMNS` in `constants.py`
- **Column aliases**: `CHAIN_COLUMN_ALIASES` in `constants.py`
- **Chain partitioning**: `_partition_chain_fields()` in `process_pcp_data.py`
- Extra columns not in `KNOWN_PCP_COLUMNS` are automatically captured on nodes

## Linting

```bash
ruff check .
ruff format .
```

## Dependencies

Runtime dependencies are in `pyproject.toml` under `[project.dependencies]`. Dev dependencies under `[project.optional-dependencies.dev]`.

Key dependencies:
- **ete3**: Newick tree parsing
- **pyyaml**: YAML config support
- **jsonschema**: Output validation
- **numpy/scipy**: Metric computations
- **tqdm**: Progress bars

## Release Process

The package is published to PyPI as `olmsted-cli`. Version is in `pyproject.toml`.

```bash
pip install build twine
python -m build
twine upload dist/*
```
