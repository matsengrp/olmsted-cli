# Olmsted CLI Architecture

See also:
- **[FORMATS.md](./FORMATS.md)**: Input/output format specs, field mapping, validation rules
- **[DEVELOPMENT.md](./DEVELOPMENT.md)**: Setup guide, testing, how-to guides
- **[CLAUDE.md](./CLAUDE.md)**: Code quality rules, terminology, quick reference
- **[README.md](./README.md)**: User documentation and examples

## Table of Contents

- [Overview](#overview)
- [Module Dependency Hierarchy](#module-dependency-hierarchy)
- [Processing Pipelines](#processing-pipelines)
- [Field Metadata System](#field-metadata-system)
- [YAML Config System](#yaml-config-system)
- [PCP Column Handling](#pcp-column-handling)
- [Phylogenetic Metrics](#phylogenetic-metrics)
- [Output Format](#output-format)
- [Key Files Reference](#key-files-reference)

---

## Overview

olmsted-cli converts immunological sequencing data (AIRR JSON, PCP CSV) into Olmsted JSON format for the [Olmsted web application](https://github.com/matsengrp/olmsted). The CLI generates `field_metadata` on each dataset describing available fields for dynamic visualization controls.

### Technology Stack

- **Python 3.8+**
- **ete3**: Newick tree parsing
- **pyyaml**: YAML config support
- **jsonschema**: Output validation
- **numpy/scipy**: Metric computations

---

## Module Dependency Hierarchy

```
constants.py            (no imports — bottom of hierarchy)
types.py                (no imports)
utils.py                (no project imports — general-purpose utilities)
    │
    ├── schemas.py          (imports constants)
    ├── metrics.py          (standalone)
    ├── field_metadata.py   (imports constants)
    ├── format_detection.py (imports constants)
    │
    ├── build_config.py     (imports constants, field_metadata, format_detection)
    │
    ├── process_utils.py    (imports build_config, field_metadata, schemas, utils)
    ├── process_pcp_data.py (imports constants, metrics, process_utils)
    ├── process_airr_data.py(imports constants, metrics, process_utils)
    │
    ├── process_data.py     (imports constants, format_detection, process_utils,
    │                         process_pcp_data, process_airr_data)
    ├── tag.py              (imports process_data, process_utils)
    │
    └── cli.py              (imports all command modules)
```

**Key rules**:
- `process_pcp_data.py` and `process_airr_data.py` never import from each other. Shared logic lives in `metrics.py` or `process_utils.py`.
- `utils.py`, `constants.py`, and `types.py` have **no project dependencies** — any module can import from them without creating cycles.
- `format_detection.py` is a leaf module (depends only on constants) so it can be imported by both `build_config.py` and `process_data.py` without cycles.

---

## Processing Pipelines

### `process` Command

```
CLI args + YAML config
    │
    ▼
get_args() ──→ build_parser() + load_config() + merge
    │
    ▼
detect_file_format()
    │
    ├── PCP: process_pcp_format(args)
    │       │
    │       ├── parse_pcp_csv()
    │       │   ├── _normalize_column_names()    (alias mapping)
    │       │   ├── Extra columns captured on nodes
    │       │   └── CSV values coerced (int/float/JSON/bool/string)
    │       │
    │       ├── parse_newick_csv()
    │       │   └── Extra columns captured as family-level fields
    │       │
    │       ├── process_pcp_to_olmsted()
    │       │   ├── merge_tree_topology_with_pcp()
    │       │   ├── Node processing (heavy + light chain partitioning)
    │       │   ├── _partition_chain_fields() for extras
    │       │   ├── compute_tree_metrics() if --compute-metrics
    │       │   ├── mean_mut_freq calculation
    │       │   └── tag_field_metadata()
    │       │
    │       └── create_consolidated_data() → write_out()
    │
    └── AIRR: process_airr_format(args)
            │
            ├── Read JSON, validate
            ├── process_dataset()
            │   ├── process_clone()   (position adjustment, sample lookup)
            │   ├── process_tree()    (tree parsing, node processing)
            │   ├── compute_tree_metrics() if --compute-metrics
            │   └── tag_field_metadata()
            │
            └── create_consolidated_data() → write_out()
```

### `tag` Command

```
Input Olmsted JSON
    │
    ▼
Load JSON → validate structure
    │
    ▼
Ensure metadata.format = "olmsted"
    │
    ▼
For each dataset:
    ├── Collect clones and matching trees
    ├── tag_field_metadata(clones, trees, custom_fields)
    │       ├── generate_default_config() if no custom_fields
    │       ├── unpack_encoded_mutations()
    │       └── generate_field_metadata()
    └── Merge with existing field_metadata (add mode) or replace (overwrite mode)
    │
    ▼
Write output
```

### `build-config` Command

```
Input (any format)
    │
    ▼
detect_file_format()
    │
    ├── Olmsted: load directly
    ├── PCP: process through pipeline in memory
    └── AIRR: process through pipeline in memory
    │
    ▼
_build_yaml()
    ├── Processing options template (inputs, output, format, etc.)
    ├── Field declarations header (types, aliases, syntax)
    │
    ├── Active fields by level:
    │   ├── Family level (clone) — from clone data
    │   ├── Node level — from tree node data
    │   ├── Branch level — from known branch fields
    │   └── Mutation level — from mutations arrays + derived AA
    │
    └── Skipped fields section (from SUGGESTED_SKIP_FIELDS)
```

---

## Field Metadata System

### How Fields Are Classified

```
Field on a clone/node/mutation dict
    │
    ▼
In EXCLUDED_*_FIELDS? ──yes──→ Hidden (never appears anywhere)
    │ no
    ▼
In KNOWN_*_FIELDS registry? ──yes──→ Use registry type + label
    │ no
    ▼
Infer type from values:
    ├── All numeric → "continuous"
    ├── Single-char AA-only → "aa"
    ├── Single-char DNA-only → "dna"
    ├── All strings → "categorical"
    ├── Booleans → "categorical"
    └── Mixed/complex → "tooltip"
    │
    ▼
In SUGGESTED_FIELD_TYPES? ──yes──→ Override inferred type (build-config only)
    │
    ▼
In SUGGESTED_SKIP_FIELDS? ──yes──→ Render with skip: true (build-config only)
    │
    ▼
Custom field declaration? ──yes──→ Override type/label/output_name
```

### Derived Mutation Fields

When nodes have `sequence_alignment_aa` but no `mutations` arrays, the field_metadata still declares:
- `child_aa` (type: `aa`) — derived by the web app during alignment rendering
- `parent_aa` (type: `tooltip`) — derived alongside child_aa

This is a "promise" that the web app will create these fields at render time.

### The `skip` Keyword

`skip: true` on a custom field entry tells `_apply_custom_fields()` to remove that field from the output metadata. It's separate from `type` — the field retains its type and label in the config for documentation, but is excluded from the output JSON.

### Level Aliasing

`family` is the user-facing name for the clonal family level. Internally, the output JSON uses `clone` for backward compatibility with the web app:

```
Config YAML: level: family
    │
    ▼ normalize_level()
Internal: level: clone
    │
    ▼
Output JSON: field_metadata.clone.{...}
```

---

## YAML Config System

### Precedence

```
argparse defaults  <  YAML config values  <  explicit CLI arguments
```

### Two-Pass Argument Parsing

1. **First pass**: Parse with all defaults set to `None` to detect which args the user explicitly provided on the command line
2. **Load config**: `load_config()` reads YAML, validates keys, resolves paths, parses `custom_fields`
3. **Merge**: For each arg still at `None` (not explicitly provided), substitute the config value
4. **Defaults**: Anything still unset gets the argparse default

### Config ↔ Command Compatibility

The same config file can be used with both `process` and `tag`. `process`-specific args (inputs, format, compute_metrics, etc.) are recognized by `load_config()` but silently ignored by `tag`, which only reads `custom_fields`.

---

## PCP Column Handling

### Standard vs Extra Columns

```
CSV columns
    │
    ├── In KNOWN_PCP_COLUMNS? → Parsed by standard PCP logic
    ├── Empty/None name? → Filtered out (unnamed index columns)
    └── Everything else → Captured as extra node-level fields
```

### Column Alias Mapping

`_normalize_column_names()` maps common alternative names to canonical PCP column names before parsing:

```
v_gene    → v_gene_heavy    (when v_gene_heavy not already present)
v_call    → v_gene_heavy
parent_seq → parent_heavy
```

### Chain Partitioning (Paired Data)

`_partition_chain_fields()` splits extra fields by suffix:

```
foobar_score       → shared (both heavy and light clones)
foobar_score_heavy → heavy clone only (suffix stripped → foobar_score)
foobar_score_light → light clone only (suffix stripped → foobar_score)
```

Chain-specific values override shared defaults when both are present.

### Tree CSV Extra Columns

Extra columns in the tree CSV become family-level (clone-level) fields on the output clone objects. Same chain partitioning applies for paired data.

---

## Phylogenetic Metrics

`metrics.py` provides format-agnostic tree metric computation:

| Metric | Function | What it needs |
|--------|----------|---------------|
| **LBI** (Local Branching Index) | `compute_lbi_for_tree()` | nodes_dict, edges, root_id, tau |
| **LBR** (Local Branching Ratio) | `compute_lbr_for_tree()` | nodes_dict, edges, root_id |
| **Scaled Affinity** | `compute_scaled_affinity()` | affinity_values dict |
| **All three** | `compute_tree_metrics()` | nodes_dict, edges, root_id, tau |

These work on any tree with branch lengths — both PCP and AIRR data.

---

## Output Format

### Consolidated Olmsted JSON

```json
{
  "metadata": {
    "format": "olmsted",
    "format_version": "1.0",
    "schema_version": "2.0.0",
    "source_format": "pcp",
    "created_at": "...",
    "generated_by": {"tool": "olmsted-cli", "version": "2.0.0"},
    "name": "My Dataset",
    "processing_info": {
      "datasets_count": 1,
      "total_clones_count": 8,
      "total_trees_count": 8,
      "total_leaf_nodes_count": 186
    }
  },
  "datasets": [
    {
      "dataset_id": "...",
      "field_metadata": {
        "clone": {
          "unique_seqs_count": {"type": "continuous", "label": "Unique Sequences Count"},
          "v_call": {"type": "categorical", "label": "V Gene"}
        },
        "node": {
          "lbi": {"type": "continuous", "label": "LBI"}
        },
        "mutation": {
          "child_aa": {"type": "aa", "label": "Child Amino Acid"},
          "selection_contribution": {"type": "continuous", "label": "Selection Contribution", "range": [-2.5, 5.1]}
        }
      }
    }
  ],
  "clones": {
    "dataset-id": [...]
  },
  "trees": [...]
}
```

### Format Detection

`detect_file_format()` identifies input format:

1. CSV extension → `pcp`
2. JSON with `metadata.format == "olmsted"` → `olmsted` (explicit tag)
3. JSON with `datasets` + `metadata` keys → `olmsted` (heuristic)
4. JSON with `dataset_id` or `clones` → `airr`
5. Otherwise → `unknown`

---

## Key Files Reference

| File | Purpose |
|------|---------|
| `constants.py` | All enums, registries, exclusion sets, alias tables |
| `types.py` | TypedDict definitions for OlmstedNode, OlmstedClone, etc. |
| `schemas.py` | JSON Schema for validation (auto-generated from constants) |
| `metrics.py` | LBI, LBR, scaled_affinity computation |
| `field_metadata.py` | generate_field_metadata() and helpers |
| `utils.py` | General-purpose utilities (VerbosePrinter, dict helpers, translate_dna_to_aa) |
| `format_detection.py` | File format detection (AIRR, PCP, Olmsted) |
| `process_data.py` | CLI entry for process, YAML config loading |
| `process_pcp_data.py` | PCP CSV parsing, column handling, clone assembly |
| `process_airr_data.py` | AIRR JSON processing |
| `process_utils.py` | tag_field_metadata(), create_consolidated_data(), write_out(), validation |
| `build_config.py` | build-config command, generate_default_config() |
| `tag.py` | tag command |
| `cli.py` | Subcommand routing |

---

_Last updated: 2026-04-02_
