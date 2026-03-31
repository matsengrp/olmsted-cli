# CLAUDE.md

This file provides guidance to Claude Code when working with the olmsted-cli repository.

## Project Overview

olmsted-cli is a Python CLI tool that processes immunological data (AIRR JSON and PCP CSV formats) into Olmsted JSON format for visualization in the [Olmsted web application](https://github.com/matsengrp/olmsted). It also generates field metadata describing available data fields for dynamic visualization controls.

### Related Repository

- **Olmsted web app**: https://github.com/matsengrp/olmsted (React/Redux/Vega)
- olmsted-cli produces JSON that the web app consumes
- `field_metadata` on datasets drives the web app's dropdown controls

## Design Documents

- **[README.md](./README.md)**: User documentation and examples
- **[ARCHITECTURE.md](./ARCHITECTURE.md)**: Data flow, processing pipelines, field_metadata system, config flow
- **[FORMATS.md](./FORMATS.md)**: Input/output format specs, field mapping, validation rules
- **[DEVELOPMENT.md](./DEVELOPMENT.md)**: Setup guide, testing, project structure, how-to guides

## Common Commands

| Command | Purpose |
|---------|---------|
| `pip install -e ".[dev]"` | Install with dev dependencies |
| `pytest` | Run all tests |
| `pytest tests/test_field_metadata.py -v` | Run specific test file |
| `ruff check .` | Lint |
| `ruff format .` | Format |

## Code Quality Rules

- **No inline TODOs** — file GitHub issues instead
- **Constants go in `constants.py`** — never hardcode enums, registries, or reference tables inline
- **Format processors don't cross-import** — shared logic goes in `metrics.py` or `process_utils.py`
- **Tests for every new feature** — unit tests for logic, integration tests for CLI commands
- `constants.py` and `types.py` have **no dependencies** on other project modules

## Terminology

| User-facing | Internal | Notes |
|-------------|----------|-------|
| `family` | `clone` | Clonal family level. Configs use "family", output JSON uses "clone" |
| `skip: true` | (field excluded) | Custom field keyword to exclude from metadata |
| `field_metadata` | `field_metadata` | Dict on each dataset describing available fields |

### Field Levels

- **family/clone**: Clonal family level (scatterplot axes, color, facet)
- **node**: Tree node level (node properties, tooltips)
- **branch**: Tree branch level (branch coloring, width)
- **mutation**: Per-mutation level (alignment coloring)

### Field Types

- **continuous**: Numeric values (axes, size, color scales)
- **categorical**: String/enum values (color, shape, facet)
- **tooltip**: Display-only (shown in tooltips)
- **aa**: Amino acid identity (uses full genetic alphabet)
- **dna**: Nucleotide identity (uses full genetic alphabet)

## Quick Reference: Modifying Constants

### Adding a Known Field

Add to the appropriate registry in `constants.py`:
```python
KNOWN_CLONE_FIELDS["my_field"] = {"type": "continuous", "label": "My Field"}
```

### Adding a Field Type or Level

1. Add to `FIELD_TYPES` or `FIELD_LEVELS` in `constants.py`
2. If level alias, add to `LEVEL_ALIASES`
3. Schema auto-updates from constants
4. Add generator in `field_metadata.py` and section in `build_config.py`
5. Regenerate golden data

### Adding Suggested Skip or Type Override

```python
# In constants.py:
SUGGESTED_SKIP_FIELDS.add("my_internal_field")
SUGGESTED_FIELD_TYPES["my_field"] = "tooltip"
```

### Regenerating Golden Data

```bash
olmsted process -f airr -i example_data/airr/airr.json --split-files example_data/airr/split_golden_data --seed 42 --name airr-example -q
olmsted process -f airr -i example_data/airr/airr.json -o example_data/airr/consolidated_golden_data.json --seed 42 --name airr-example -q
olmsted process -f pcp -i example_data/pcp/pcp.csv -t example_data/pcp/trees.csv --split-files example_data/pcp/split_golden_data --seed 42 --name pcp-example -q
olmsted process -f pcp -i example_data/pcp/pcp.csv -t example_data/pcp/trees.csv -o example_data/pcp/consolidated_golden_data.json --seed 42 --name pcp-example -q
```

---

_Last updated: 2026-03-31_
