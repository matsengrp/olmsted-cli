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
| `tree_name` | `tree_name` | Per-tree-within-clone identifier; drives the webapp's tree dropdown when one clone has multiple alternate reconstructions |
| `skip: true` | (field excluded) | Custom field keyword to exclude from metadata |
| `field_metadata` | `field_metadata` | Dict on each dataset describing available fields |

### Identifier Fields

| Field | Role |
|-------|------|
| `*_id` (e.g. `clone_id`, `tree_id`, `dataset_id`) | **Input-derived** identifier. When synthesis is unavoidable, use `{datatype}-{uuid}` — never format-origin prefixes like `pcp-`. |
| `ident` | **CLI-minted** primary key, always via `IdentMinter.mint(datatype)` producing `{datatype}-{uuid}`. Only on objects the webapp keys on it (today: `clone`, `tree`). |

See `ARCHITECTURE.md#identifier-conventions` for the full rules and per-field uniqueness scopes.

### Field Levels

- **family/clone**: Clonal family level (scatterplot axes, color, facet)
- **tree**: Per-tree-within-clone level (drives the tree dropdown's color/filter/sort controls; only populated when at least one clone has trees that disagree on the field's value)
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

All regen commands pass `--json-format pretty` explicitly. Pretty is the
current default, but pinning it locks in test-friendly formatting if the
default ever changes. The `gzip` variants below pin the gzip header
(`mtime=0`, empty filename) so the *compression layer* is deterministic;
note that the JSON content itself still varies slightly between runs
(`metadata.created_at`, some field-iteration ordering), so a tracked
`.json.gz` will show a small diff on every regeneration. Tests compare
decompressed content, not gzipped bytes.

Consolidated goldens (the canonical single-file output for each dataset):

```bash
olmsted process -f airr -i example-data/airr/input-airr.json -o example-data/airr/airr-olmsted-golden.json --seed 42 --name airr-example --json-format pretty -q
olmsted process -f pcp -i example-data/pcp/input-pcp.csv -t example-data/pcp/input-trees.csv -o example-data/pcp/pcp-olmsted-golden.json --seed 42 --name pcp-example --json-format pretty -q
olmsted process -f pcp -i example-data/pcp-byhand/input-pcp.csv -t example-data/pcp-byhand/input-trees.csv -o example-data/pcp-byhand/pcp-byhand-olmsted-golden.json --seed 42 --name pcp-byhand-example --json-format pretty -q
olmsted process -f pcp -i example-data/pcp-light/input-pcp.csv -t example-data/pcp-light/input-trees.csv -o example-data/pcp-light/pcp-light-olmsted-golden.json --seed 42 --name pcp-light-example --json-format pretty -q
olmsted process -f pcp -i example-data/pcp-paired/input-pcp.csv -t example-data/pcp-paired/input-trees.csv -o example-data/pcp-paired/pcp-paired-olmsted-golden.json --seed 42 --name pcp-paired-example --json-format pretty -q
```

Gzipped consolidated goldens (tracked alongside the plain JSON for `.json.gz` upload coverage):

```bash
olmsted process -f airr -i example-data/airr/input-airr.json -o example-data/airr/airr-olmsted-golden.json --seed 42 --name airr-example --json-format gzip -q
olmsted process -f pcp -i example-data/pcp/input-pcp.csv -t example-data/pcp/input-trees.csv -o example-data/pcp/pcp-olmsted-golden.json --seed 42 --name pcp-example --json-format gzip -q
```

Split-format goldens (legacy, pinned for integrity testing as long as `--split-files` is supported):

```bash
olmsted process -f airr -i example-data/airr/input-airr.json --split-files example-data/airr/split-golden-data --seed 42 --name airr-example --json-format pretty -q
olmsted process -f pcp -i example-data/pcp/input-pcp.csv -t example-data/pcp/input-trees.csv --split-files example-data/pcp/split-golden-data --seed 42 --name pcp-example --json-format pretty -q
```

Merge golden (post-merge source-of-truth):

```bash
olmsted merge -i example-data/merge/input-olmsted.json --mutations example-data/merge/input-mutations.csv --mutations-use-depth -o example-data/merge/merge-olmsted-golden.json --json-format pretty -q
```

---

_Last updated: 2026-03-31_
