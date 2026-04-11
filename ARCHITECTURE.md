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
- [Mutations CSV Merge](#mutations-csv-merge)
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
    ├── merge_mutations.py  (imports constants, process_utils, utils)
    │
    ├── process_data.py     (imports constants, format_detection, merge_mutations,
    │                         process_utils, process_pcp_data, process_airr_data)
    ├── tag.py              (imports process_data, process_utils)
    ├── merge.py            (imports merge_mutations, process_data, process_utils)
    │
    └── cli.py              (imports all command modules)
```

**Key rules**:
- `process_pcp_data.py` and `process_airr_data.py` never import from each other. Shared logic lives in `metrics.py` or `process_utils.py`.
- `utils.py`, `constants.py`, and `types.py` have **no project dependencies** — any module can import from them without creating cycles.
- `format_detection.py` is a leaf module (depends only on constants) so it can be imported by both `build_config.py` and `process_data.py` without cycles.
- `merge_mutations.py` is the single source of truth for the mutations-CSV merge logic. Both the `merge` command and `process --mutations` flag call into `apply_mutations_csv()` so behavior is identical regardless of entry point.

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
retag_datasets_field_metadata(datasets, clones_dict, trees, custom_fields, mode)
    ├── For each dataset:
    │   ├── Collect clones and matching trees
    │   ├── tag_field_metadata(clones, trees, custom_fields)
    │   │       ├── generate_default_config() if no custom_fields
    │   │       ├── unpack_encoded_mutations()
    │   │       └── generate_field_metadata()
    │   └── Merge with existing field_metadata (add) or replace (overwrite)
    │
    ▼
Write output
```

### `merge` Command

```
Input Olmsted JSON + mutations CSV
    │
    ▼
Load JSON → validate structure
    │
    ▼
apply_mutations_csv(path, datasets, clones_dict, trees, custom_fields)
    ├── load_mutations_csv()             (parse CSV, group by family, force site→int)
    ├── merge_mutations_into_trees()     (derive mutations from AA diffs, match, merge)
    ├── Print warnings for unmatched families / unmatched mutations
    └── retag_datasets_field_metadata()  (regenerate field_metadata with merge=add)
    │
    ▼
Refuse to overwrite if --in-place and zero trees matched
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

## Mutations CSV Merge

The `merge` command and `process --mutations` flag share a single implementation in `merge_mutations.py`. By the time the merge runs, the input data — whether it started as PCP CSV, AIRR JSON, or pre-built Olmsted JSON — has already been converted to the internal `(datasets, clones_dict, trees)` representation. The merge operates on that representation, so behavior is identical across all three entry points.

### Matching Strategy

For each tree node:

1. If the node already has a `mutations` array, use it as-is.
2. Otherwise, derive mutations by diffing `node.sequence_alignment_aa` against its parent's. Each differing position becomes `{site, parent_aa, child_aa}`. Gap characters (`-`, `.`, `X`, `*`, `?`) are skipped.
3. For each derived (or pre-existing) mutation, look up the join key in the CSV index for that tree's `clone_id`. On match, merge the CSV's score columns onto the mutation dict.

The CSV's `family` column is the join key against `tree.clone_id`. The base join key is `(site, parent_aa, child_aa)`, optionally extended by disambiguation columns (see below).

### Disambiguation Columns

Naive matching by `(site, parent_aa, child_aa)` alone causes **fan-out**: a single CSV row will match every node in the tree whose parent→child diff produces that exact substitution at that exact site. Convergent mutations on independent lineages, or repeated substitutions at the same site after a back-mutation, all collapse onto the same CSV row and receive identical enrichment data — which is incorrect if the upstream pipeline computed per-event scores.

To narrow matches, optional disambiguation columns extend the join key when present in the CSV:

| Column | What it contributes | How it's matched on the tree |
|--------|---------------------|------------------------------|
| `depth` | Edges from the nearest root to the child node | Computed at merge time via BFS in `_compute_node_depths()` |

If `depth` is present in the CSV, the join key becomes `(site, parent_aa, child_aa, depth)` and node-mutations are matched only when their computed depth matches the CSV row's depth. Depth is parsed as an integer at load time and excluded from the *enriched output*.

Disambiguation is auto-detected: if any loaded row carries a disambiguation column, it's used for the entire run. The active columns are reported at status verbosity and recorded in `MergeStats.disambiguation_columns_used`.

Even with depth, ambiguity can remain — two convergent substitutions at the same site at the same depth on different lineages will still broadcast. Broadcasts are tracked and warned about (see below).

### Excluded CSV Columns

These columns are recognized as structural/join keys and are **not** included in the merged output (see `MUTATIONS_CSV_KEY_COLUMNS` in `constants.py`):

```
family, sample_id, site, parent_aa, child_aa, pcp_index, depth
```

`site`, `parent_aa`, and `child_aa` are excluded from the *merged extras dict* but are still kept on the mutation record (they identify the substitution). `depth` is excluded entirely from the merged record but is retained on the loaded row dict so it can serve as a join key.

### Stats and Reporting

`merge_mutations_into_trees()` returns a `MergeStats` dataclass with:

| Field | Meaning |
|-------|---------|
| `trees_matched` | Number of trees whose `clone_id` appeared in the CSV |
| `nodes_enriched` | Nodes that received at least one CSV-sourced field on this run |
| `mutations_enriched` | Individual `(node, mutation)` pairs that received CSV data |
| `unmatched_families` | Sorted list of CSV families that had no matching tree |
| `unmatched_family_rows` | Total CSV rows belonging to those unmatched families |
| `unmatched_mutations` | CSV rows in matched families whose join key didn't match any derived mutation |
| `broadcast_csv_rows` | CSV rows that matched **more than one** node-mutation pair (ambiguous join) |
| `disambiguation_columns_used` | Optional disambiguation columns active on this run (e.g. `["depth"]`) |

Counts are scoped to the current run: `nodes_enriched` and `mutations_enriched` exclude pre-existing mutation arrays from upstream pipelines.

`apply_mutations_csv()` reports at status verbosity:

```
Loading mutations CSV: {path}
Loaded {N} CSV rows across {M} families
Disambiguation columns in CSV: {cols}             # only if disambiguation active
Enriched {X} mutations across {Y} nodes in {Z} trees
Unmatched: {U}/{N} CSV rows ({A} in {B} unmatched families, {C} with no node match)
Broadcast: {K} CSV rows matched multiple nodes (ambiguous join — same data applied to every match)
```

Warnings (at error level, exit 0):
- **Unmatched mutations** in matched families — the CSV references substitutions that don't appear in any derived diff
- **Broadcast rows** — the CSV row's enrichment was applied to multiple node-mutations; may not be correct for per-event scores

Per-family detail at `-v 2`: which keys went unmatched and which broadcast to multiple nodes.

### Safety: `--in-place` Guard

If `--in-place` is requested and `stats.trees_matched == 0`, the `merge` command exits with an error rather than overwriting the input file. This prevents a typo in the CSV's `family` column from silently destroying the input.

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
| `constants.py` | All enums, registries, exclusion sets, alias tables, `MUTATIONS_CSV_KEY_COLUMNS` |
| `types.py` | TypedDict definitions for OlmstedNode, OlmstedClone, etc. |
| `schemas.py` | JSON Schema for validation (auto-generated from constants) |
| `metrics.py` | LBI, LBR, scaled_affinity computation |
| `field_metadata.py` | generate_field_metadata() and helpers |
| `utils.py` | General-purpose utilities (VerbosePrinter, dict helpers, translate_dna_to_aa) |
| `format_detection.py` | File format detection (AIRR, PCP, Olmsted) |
| `process_data.py` | CLI entry for process, YAML config loading |
| `process_pcp_data.py` | PCP CSV parsing, column handling, clone assembly |
| `process_airr_data.py` | AIRR JSON processing |
| `process_utils.py` | tag_field_metadata(), retag_datasets_field_metadata(), coerce_csv_value(), create_consolidated_data(), write_out(), validation |
| `merge_mutations.py` | load_mutations_csv(), derive_node_mutations(), merge_mutations_into_trees(), apply_mutations_csv(), MergeStats |
| `build_config.py` | build-config command, generate_default_config() |
| `tag.py` | tag command |
| `merge.py` | merge command |
| `cli.py` | Subcommand routing |

---

_Last updated: 2026-04-10_
