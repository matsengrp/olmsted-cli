# CLAUDE.md

This file provides guidance to Claude Code when working with the olmsted-cli repository.

## Project Overview

olmsted-cli is a Python CLI tool that processes immunological data (AIRR JSON and PCP CSV formats) into Olmsted JSON format for visualization in the [Olmsted web application](https://github.com/matsengrp/olmsted). It also generates field metadata describing available data fields for dynamic visualization controls.

### Related Repository

- **Olmsted web app**: https://github.com/matsengrp/olmsted (React/Redux/Vega)
- olmsted-cli produces JSON that the web app consumes
- `field_metadata` on datasets drives the web app's dropdown controls

## Architecture

```
olmsted_cli/
├── cli.py                 # CLI entry point, subcommand routing
├── constants.py           # All shared constants, enums, registries, ref tables
├── types.py               # TypedDict definitions for all data structures
├── schemas.py             # JSON Schema definitions (uses constants.py)
├── metrics.py             # Shared phylogenetic metrics (LBI, LBR, etc.)
├── field_metadata.py      # Field metadata generation for all levels
├── process_data.py        # Unified processor: format detection, arg parsing, YAML config
├── process_pcp_data.py    # PCP CSV parsing and conversion
├── process_airr_data.py   # AIRR JSON processing
├── process_utils.py       # Shared utilities, output writing, validation
├── build_config.py        # build-config command: generate YAML from data
├── enrich.py              # enrich command: add field_metadata to existing files
├── api.py                 # High-level OlmstedData API
├── validate.py            # validate command
├── summary.py             # summary command
├── split.py               # split command
└── configs/               # Default YAML configs shipped with package
    ├── pcp.yaml
    ├── airr.yaml
    └── surprise.yaml
```

### Key Design Patterns

- **constants.py has no dependencies** on other modules — it's the bottom of the import hierarchy
- **types.py has no dependencies** on other modules
- **metrics.py** is format-agnostic — both PCP and AIRR import from it (never cross-import between format processors)
- **field_metadata.py** imports from constants.py only
- **build_config.py** is for config generation output only — its suggestion tables (FIELD_ALIASES, SUGGESTED_SKIP_FIELDS) are not used at runtime by the processing pipeline

### Data Flow

```
Input (PCP CSV / AIRR JSON / Olmsted JSON)
    │
    ├── process command ──→ parse → transform → compute metrics → generate field_metadata → write JSON
    ├── enrich command ───→ load Olmsted JSON → generate field_metadata → merge → write JSON
    └── build-config ─────→ parse → introspect fields → generate YAML config
```

## Terminology

| User-facing | Internal | Notes |
|-------------|----------|-------|
| `family` | `clone` | Clonal family level. Config uses "family", output JSON uses "clone" |
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

## Development

### Setup

```bash
git clone https://github.com/matsengrp/olmsted-cli.git
cd olmsted-cli
pip install -e ".[dev]"
```

### Common Commands

| Command | Purpose |
|---------|---------|
| `pytest` | Run all tests |
| `pytest tests/test_field_metadata.py -v` | Run specific test file |
| `ruff check .` | Lint |

### Testing

Tests are in `tests/`. Key test files:

| File | Tests |
|------|-------|
| `test_field_metadata.py` | Field registries, type inference, AA/DNA detection, ranges |
| `test_build_config.py` | Config generation, format detection, alias suggestions |
| `test_config.py` | YAML config loading, custom_fields parsing |
| `test_enrich.py` | Enrich command end-to-end |
| `test_pcp_extras.py` | Extra columns, chain partitioning, column aliases |
| `test_format_detection.py` | Olmsted/AIRR/PCP format detection |
| `test_cli_processing.py` | Full CLI integration against golden data |
| `test_validation.py` | Schema validation |

### Golden Data

`example_data/airr/` and `example_data/pcp/` contain golden data files used by integration tests. Regenerate after changes to output format:

```bash
olmsted process -f airr -i example_data/airr/airr.json --split-files example_data/airr/split_golden_data --seed 42 --name airr-example -q
olmsted process -f airr -i example_data/airr/airr.json -o example_data/airr/consolidated_golden_data.json --seed 42 --name airr-example -q
olmsted process -f pcp -i example_data/pcp/pcp.csv -t example_data/pcp/trees.csv --split-files example_data/pcp/split_golden_data --seed 42 --name pcp-example -q
olmsted process -f pcp -i example_data/pcp/pcp.csv -t example_data/pcp/trees.csv -o example_data/pcp/consolidated_golden_data.json --seed 42 --name pcp-example -q
```

### Adding New Field Types or Levels

1. Add to `FIELD_TYPES` or `FIELD_LEVELS` in `constants.py`
2. If level alias needed, add to `LEVEL_ALIASES`
3. Update `_FIELD_ENTRY_SCHEMA` enum in `schemas.py` (auto-generated from constants)
4. Update build-config output sections in `build_config.py`
5. Regenerate golden data

### Adding New Known Fields

Add to the appropriate registry in `constants.py`: `KNOWN_CLONE_FIELDS`, `KNOWN_NODE_FIELDS`, `KNOWN_BRANCH_FIELDS`, or `KNOWN_MUTATION_FIELDS`.

### Adding Suggested Skip or Type Overrides

- `SUGGESTED_SKIP_FIELDS` in `constants.py`: fields shown with `skip: true` in build-config
- `SUGGESTED_FIELD_TYPES` in `constants.py`: type overrides where auto-inference is wrong

## Code Quality

- **No inline TODOs** — file GitHub issues instead
- **Constants go in constants.py** — never hardcode enums, registries, or reference tables inline
- **Format processors don't cross-import** — shared logic goes in metrics.py or process_utils.py
- **Tests for every new feature** — unit tests for logic, integration tests for CLI commands
