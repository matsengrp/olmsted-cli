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
| `test_merge.py` | Merge command + mutations CSV utility (sequence diff derivation, CSV loading, end-to-end merge, unmatched reporting) |
| `test_pcp_extras.py` | Extra CSV columns, chain partitioning, column aliases, coercion |
| `test_format_detection.py` | Olmsted/AIRR/PCP format detection |
| `test_validation.py` | Schema validation for datasets, clones, trees |
| `test_streaming.py` | Streaming primitives in isolation (`FieldTypeEvidence`, `RangeEvidence`, `BatchAccumulator`, `BatchSpooler`, `write_olmsted_json_streaming`) |
| `test_clone_group_iterators.py` | Per-clone-group iterators (`iter_pcp_clone_groups`, `iter_airr_clones`) — batch-size invariance, alt-reconstruction co-emission |
| `test_batching_pcp.py` | End-to-end PCP streaming: every example dataset matches its golden across `--batch-size in {1, 2, 50, 10000}`; synthetic hoist-correctness test |
| `test_batching_airr.py` | End-to-end AIRR streaming: example matches golden across `--batch-size in {1, 2, 50, 10000}` |
| `test_batching_mutations.py` | `process --mutations` streaming output matches the legacy `--batch-size 0` run across batch sizes |

### Golden Data

Integration tests compare CLI output against golden data in `example-data/`. Regenerate after output format changes:

All regen commands pass `--json-format pretty` explicitly. Pretty is the
current default but pinning it locks in test-friendly formatting if the
default ever changes. The `gzip` variants below pin the gzip header
(`mtime=0`, empty filename) so the compression layer is deterministic;
the JSON content still varies between runs (`metadata.created_at`, some
field-iteration ordering), so a tracked `.json.gz` will show a small diff
on every regeneration. Tests compare decompressed content, not bytes.

PCP regen-output expectations (post-issue #23):

- `clones[]` is keyed on `(sample, family)` — so families that repeat
  across samples (e.g. the pilot PCP data) now emit a clone per
  `(sample, family)` pair rather than silently merging.
- `clone_id` is synthesized as `{sample}_{family}` for the same reason.
- `field_metadata.tree` appears only for datasets where at least one
  clone has multiple alternate reconstructions (AIRR's
  `downsampling_strategy`, or PCP inputs with a `tree_name` column).
- Each tree record carries a `tree_name` field — either the input
  tree-role value or the synthesized `tree-{clone_id}` fallback.

Consolidated goldens (one per dataset folder):

```bash
olmsted process -f airr -i example-data/airr/input-airr.json \
  -o example-data/airr/airr-olmsted-golden.json --seed 42 --name airr-example --json-format pretty -q
olmsted process -f pcp -i example-data/pcp/input-pcp.csv -t example-data/pcp/input-trees.csv \
  -o example-data/pcp/pcp-olmsted-golden.json --seed 42 --name pcp-example --json-format pretty -q
olmsted process -f pcp -i example-data/pcp-byhand/input-pcp.csv -t example-data/pcp-byhand/input-trees.csv \
  -o example-data/pcp-byhand/pcp-byhand-olmsted-golden.json --seed 42 --name pcp-byhand-example --json-format pretty -q
olmsted process -f pcp -i example-data/pcp-light/input-pcp.csv -t example-data/pcp-light/input-trees.csv \
  -o example-data/pcp-light/pcp-light-olmsted-golden.json --seed 42 --name pcp-light-example --json-format pretty -q
olmsted process -f pcp -i example-data/pcp-paired/input-pcp.csv -t example-data/pcp-paired/input-trees.csv \
  -o example-data/pcp-paired/pcp-paired-olmsted-golden.json --seed 42 --name pcp-paired-example --json-format pretty -q
```

Gzipped consolidated goldens (tracked alongside plain JSON for `.json.gz` upload coverage):

```bash
olmsted process -f airr -i example-data/airr/input-airr.json \
  -o example-data/airr/airr-olmsted-golden.json --seed 42 --name airr-example --json-format gzip -q
olmsted process -f pcp -i example-data/pcp/input-pcp.csv -t example-data/pcp/input-trees.csv \
  -o example-data/pcp/pcp-olmsted-golden.json --seed 42 --name pcp-example --json-format gzip -q
```

Split-format goldens (legacy, kept for integrity testing while `--split-files` is supported):

```bash
olmsted process -f airr -i example-data/airr/input-airr.json \
  --split-files example-data/airr/split-golden-data --seed 42 --name airr-example --json-format pretty -q
olmsted process -f pcp -i example-data/pcp/input-pcp.csv -t example-data/pcp/input-trees.csv \
  --split-files example-data/pcp/split-golden-data --seed 42 --name pcp-example --json-format pretty -q
```

Merge golden (post-merge drift coverage for `olmsted merge`):

```bash
olmsted merge -i example-data/merge/input-olmsted.json \
  --mutations example-data/merge/input-mutations.csv --mutations-use-depth \
  -o example-data/merge/merge-olmsted-golden.json --json-format pretty -q
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
│   ├── column_resolution.py  # Auto-detect sample/family/tree role columns in CSV inputs
│   ├── utils.py              # General-purpose utilities (no project deps)
│   ├── format_detection.py   # File format detection
│   ├── process_utils.py      # Processing utilities, output writing, validation
│   ├── merge_mutations.py    # Mutations CSV merge logic (shared by merge and process --mutations)
│   ├── streaming.py          # Streaming primitives: BatchAccumulator, BatchSpooler, evidence accumulators, write_olmsted_json_streaming
│   ├── build_config.py       # build-config command, generate_default_config()
│   ├── tag.py                # tag command
│   ├── merge.py              # merge command
│   ├── api.py                # Programmatic API
│   └── configs/              # Default YAML configs
├── tests/                    # Test suite
├── example-data/             # Example datasets and golden data
│   ├── airr/                 # AIRR format + golden data
│   ├── pcp/                  # PCP format + golden data
│   ├── pcp-paired/           # Paired heavy/light PCP
│   ├── mutations/            # Pre-built Olmsted JSON with mutation-level data
│   └── fields-config/        # Foobar test data (all types × levels) + config
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
- **Role-column auto-detection** (sample/family/tree): `column_resolution.py`
- **Chain partitioning**: `_partition_chain_fields()` in `process_pcp_data.py`

## Linting

```bash
ruff check .
ruff format .
```

## Release

Releases are **tag-driven**: the version comes from the git tag (via
`setuptools-scm`), and pushing a `vX.Y.Z` tag triggers a GitHub Actions
workflow that builds and publishes to PyPI via Trusted Publishing — no manual
`build`/`twine` and no stored tokens.

```bash
git tag v0.4.0 && git push origin v0.4.0   # publishes to PyPI
```

See **[RELEASING.md](RELEASING.md)** for the full flow, the TestPyPI dry-run
path, version-number guidance, and one-time trusted-publishing setup.

---

_Last updated: 2026-06-11_
