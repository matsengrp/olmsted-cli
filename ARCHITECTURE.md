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
- [Streaming Pipeline](#streaming-pipeline)
- [Mutations CSV Merge](#mutations-csv-merge)
- [Field Metadata System](#field-metadata-system)
- [YAML Config System](#yaml-config-system)
- [PCP Column Handling](#pcp-column-handling)
- [Phylogenetic Metrics](#phylogenetic-metrics)
- [Output Format](#output-format)
- [Identifier Conventions](#identifier-conventions)
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
    │       ├── _should_stream_pcp(args)?
    │       │   ├── yes + n_families > batch_size →
    │       │   │     _process_pcp_streaming() — see Streaming Pipeline
    │       │   └── no (or single-batch fast path) →
    │       │         process_pcp_to_olmsted()
    │       │           ├── iter_pcp_clone_groups()
    │       │           ├── _process_clone_group() (heavy + light partitioning,
    │       │           │     compute_tree_metrics(), mean_mut_freq)
    │       │           ├── _hoist_clone_invariant_extras()
    │       │           └── tag_field_metadata()
    │       │         create_consolidated_data() → write_out()
    │
    └── AIRR: process_airr_format(args)
            │
            ├── _should_stream_airr(args)?
            │   ├── yes → _process_airr_streaming() — see Streaming Pipeline
            │   └── no → process_dataset() for each input file
            │             ├── iter_airr_clones()
            │             ├── _process_airr_clone() (positions, sample lookup,
            │             │     process_tree(), compute_tree_metrics())
            │             └── tag_field_metadata()
            │           create_consolidated_data() → write_out()
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
populate_branch_lengths_from_newick(tree) for each tree
    └── backfill per-node length/distance from the newick when the string
        carries branch lengths but the nodes don't (no-clobber). Lets a
        hand-built base JSON drive the webapp's "distance from naive" mode.
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

## Streaming Pipeline

`olmsted process` ships with always-on family-batched streaming so peak
memory tracks `--batch-size` (default 50 families) rather than total
dataset size. Each batch's clones and trees are spooled to JSONL temp
files; the final consolidated JSON is written by stream-stitching
`metadata` → `datasets` → `clones` → `trees` directly from disk.

```
parse_pcp_csv() / read_airr_json()
    │
    ▼
BatchAccumulator.register_dataset(dataset_id, hoist_tree_extras_to_clone=?)
    └── PCP: True   (tree-csv extras hoist to clone-level)
        AIRR: False (tree-level fields live on trees natively)
    │
    ▼
begin_merge(mutations_csv)            (optional, --mutations)
    │
    ▼
for batch in iter_pcp_clone_groups | iter_airr_clones(batch_size=N):
    ├── unpack_encoded_mutations()    (if user passed custom_fields with encoding)
    ├── apply_mutations_to_trees()    (if --mutations)
    ├── BatchAccumulator.observe_batch()
    │     └── folds in field-type/range evidence, ID-uniqueness sets,
    │         running totals, tree-level-variance keys
    └── BatchSpooler.write_batch()    (JSONL, one record per line)
    │
    ▼
finalize_merge(merge_ctx) + report_merge_stats()    (if --mutations)
    │
    ▼
apply_dataset_hoist(spooler, dataset_id, tree_level_keys)    (PCP only)
    └── rewrites spool JSONL: keys constant across a clone's trees
        AND not in tree_level_keys move from trees to the clone
    │
    ▼
BatchAccumulator.finalize_field_metadata(dataset_id, custom_fields)
    └── mirrors generate_field_metadata's shape; applies FIELD_ALIASES
        when custom_fields is None (rearrangement_count → unique_seqs_count, …)
    │
    ▼
create_consolidated_data() builds the metadata wrapper;
    accumulator.finalize_totals() overrides processing_info
    │
    ▼
write_olmsted_json_streaming()
    └── emits {metadata, datasets, clones, trees} in canonical key order,
        stream-stitching clones/trees from spooler. gzip output pins
        mtime=0 + empty filename (same determinism guarantee as the
        legacy data_io.write_olmsted_json).
```

### Fallback to the legacy in-memory path

`_should_stream_pcp` / `_should_stream_airr` route through
`process_pcp_to_olmsted` / `process_dataset` + `write_out` when:

| Condition | Why |
|---|---|
| `--batch-size 0` | Explicit opt-out. |
| `--split-files DIR` | Multi-file output has a different write shape. |
| `--validate` | Per-batch validation isn't wired yet; `validate_output_data` consumes the whole assembled output. |
| Single-batch fast path (PCP only): `n_families ≤ batch_size` | Spool round-trip would cost more than the in-memory pipeline; skipped automatically. |

### Per-batch correctness notes

- **Type inference** uses `FieldTypeEvidence` counters across all batches,
  so a value contradicting the inferred type later (e.g. a string
  arriving after 50 ints) still flips the result. The legacy
  sample-capped path can miss this.
- **Ranges** track running `(min, max, count)` via `RangeEvidence`; merged
  across batches.
- **Tree-level classification** is per-clone (variance within one
  clone's trees) — the iterator guarantees a clone's alt
  reconstructions co-emit in one batch, so per-batch classification is
  correct. The union across batches is the dataset-scope answer.
- **`--mutations`** loads the CSV once via `load_mutations_csv` and
  threads one `MergeContext` through every batch; `MergeStats` and the
  `unmatched_family_set` aggregate exactly like a one-shot run would
  aggregate across trees.

### Key files

| File | Role |
|---|---|
| `olmsted_cli/streaming.py` | `FieldTypeEvidence`, `RangeEvidence`, `BatchAccumulator`, `BatchSpooler`, `write_olmsted_json_streaming`, `apply_dataset_hoist` |
| `olmsted_cli/process_pcp_data.py` | `iter_pcp_clone_groups`, `_process_clone_group`, `_hoist_clone_invariant_extras` |
| `olmsted_cli/process_airr_data.py` | `iter_airr_clones`, `_process_airr_clone` |
| `olmsted_cli/process_data.py` | `_should_stream_pcp` / `_should_stream_airr`, `_process_pcp_streaming` / `_process_airr_streaming`, `_begin_mutations_merge` / `_finalize_mutations_merge` |
| `olmsted_cli/merge_mutations.py` | `MergeContext`, `begin_merge`, `apply_mutations_to_trees`, `finalize_merge`, `report_merge_stats` |

---

## Mutations CSV Merge

The `merge` command and `process --mutations` flag share a single implementation in `merge_mutations.py`. By the time the merge runs, the input data — whether it started as PCP CSV, AIRR JSON, or pre-built Olmsted JSON — has already been converted to the internal `(datasets, clones_dict, trees)` representation. The merge operates on that representation, so behavior is identical across all three entry points.

### Matching Strategy

For each tree node:

1. If the node already has a `mutations` array, use it as-is.
2. Otherwise, derive mutations by diffing `node.sequence_alignment_aa` against its parent's. Each differing position becomes `{site, parent_aa, child_aa}`. Gap characters (`-`, `.`, `X`, `*`, `?`) are skipped.
3. For each derived (or pre-existing) mutation, look up the join key in the CSV index for that tree's `clone_id`. On match, merge the CSV's score columns onto the mutation dict.

The CSV's `family` column is the join key against `tree.clone_id`. What happens next depends on which optional columns are present:

### Match Mode Selection

The merge picks a match mode based on what the CSV carries. The chosen mode is reported at status verbosity and recorded in `MergeStats.match_mode`.

| CSV has… | Mode | Join key | Notes |
|---|---|---|---|
| `node_name` or `child_name` | `name_site` | `(node_name, site)` | Fully disambiguating; no broadcast possible. `parent_aa`/`child_aa` are integrity checks; `depth` is too when `--mutations-use-depth` is set, otherwise ignored. |
| `depth` + `--mutations-use-depth` | `site_paa_caa_depth` | `(site, parent_aa, child_aa, depth)` | Narrows fan-out; still allows broadcast on convergent same-depth substitutions. |
| neither | `site_paa_caa` | `(site, parent_aa, child_aa)` | May broadcast (tracked). |

`node_name`/`child_name`: When the CSV has a node-name column (either alias — `node_name` wins if both are present), values are the `sequence_id` of the target node. The loader normalizes both aliases onto the canonical `node_name` key.

`depth`: Edges from the nearest root to the child node, computed at merge time via BFS in `_compute_node_depths()` (prefers naive as origin when connected; falls back to directed root BFS). Depth is **opt-in via `--mutations-use-depth`** because depth arithmetic depends on the upstream pipeline's rooting convention, which the CLI can't infer. Without the flag, a `depth` column in the CSV is ignored *entirely* — neither as a match-key participant nor as an integrity check. When the column is seen but the flag is absent, `apply_mutations_csv` logs a verbose note. Conversely, passing `--mutations-use-depth` when the CSV has **no** `depth` column raises `ValueError` — opting in with no data to opt into is a misuse signal.

### Integrity Checks and `--mutations-allow-mismatch`

In `name_site` mode, `parent_aa`/`child_aa` (and `depth`, only when `--mutations-use-depth` is set) aren't part of the join key — they're cross-checked against the tree's derived mutation at the identified `(node, site)`. Mismatched rows are **always skipped** (never attached, regardless of flags). The flag controls whether the command exits:
- **Default:** any integrity mismatch raises `ValueError` → the command exits non-zero. Callers can't accidentally ship a partially-wrong merge.
- **`--mutations-allow-mismatch`:** downgrade to a warning; skipped rows are still reported via `MergeStats.integrity_mismatches` but the command exits 0.

Rationale: attaching upstream scores to a mutation whose parent/child residues don't match what the CSV claimed would attach data to the wrong biological event. Skipping is always safer than attributing; failing loud by default surfaces CSV/tree drift rather than letting it pass silently. The flag exists as an explicit "I know, proceed anyway" escape hatch.

### Authoritative CSV: `--mutations-listed-only`

By default the merge enriches mutations matched by the CSV but still emits sequence-diff-derived mutations that have no CSV row (just without the extra columns). Pipelines that use the CSV as an *intentional filter* — e.g. dropping multi-nt-per-codon mutations because they're noisy under a downstream model — get those filtered events back as bare mutations on the tree.

Passing `--mutations-listed-only` makes the CSV authoritative: on every tree whose `clone_id` matches a family in the CSV, derived mutations that don't have a corresponding CSV row are removed. The number dropped is reported via `MergeStats.mutations_dropped` and at status verbosity. Trees whose family is absent from the CSV are not filtered (the user has no opinion on them), so the flag is safe to combine with a CSV that only covers a subset of families.

**Caveats**

- *The `mutations` array becomes the authoritative event log; `sequence_alignment_aa` is not rewritten.* The on-disk JSON is deliberately inconsistent after the flag runs: a fresh AA diff between a node and its parent may imply events that no longer appear in the node's `mutations` array. The webapp consumes the array (not the sequences) for mutation events, so the dropped events disappear from mutation views while the sequence display still carries the residue change. Pre-existing upstream `mutations` arrays on input nodes are filtered the same way as freshly-derived ones.
- *Re-derivation resurrection.* `_get_or_derive_mutations` re-derives from the AA diff whenever a node's `mutations` array is missing or empty. Chaining a second merge over the only-listed output (especially with a different CSV) would resurrect the dropped events as bare entries again. The flag is therefore intended for terminal/output-stage use; running it earlier in a multi-merge pipeline is unlikely to do what you want.
- *Interaction with integrity mismatches.* In `name_site` mode, a CSV row that resolves to a real `(node, site)` but fails the `parent_aa`/`child_aa` integrity check is skipped (its site never enters the listed set). Paired with `--mutations-allow-mismatch` the run continues — and `--mutations-listed-only` then drops the derived mutation at that site as well, since no row authoritatively claimed it. This cascade is intentional: a rejected CSV claim is treated as "no claim," not as evidence that the bare event should survive.

### Excluded CSV Columns

These columns are recognized as structural/join keys and are **not** included in the merged output (see `MUTATIONS_CSV_KEY_COLUMNS` in `constants.py`):

```
family, sample_id, site, parent_aa, child_aa, pcp_index, depth, node_name, child_name
```

`site`, `parent_aa`, and `child_aa` are excluded from the *merged extras dict* but are still kept on the mutation record (they identify the substitution). `depth` and `node_name`/`child_name` are excluded entirely from the merged record — they're retained on the loaded row dict only so they can serve as the join key / integrity check.

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
| `broadcast_csv_rows` | CSV rows that matched **more than one** node-mutation pair (ambiguous join). Only non-zero in `site_paa_caa[_depth]` modes. |
| `integrity_mismatches` | CSV rows that resolved to a real `(node, site)` in `name_site` mode but whose `parent_aa`/`child_aa`/`depth` disagreed with the tree's derived mutation. Not attached. |
| `disambiguation_columns_used` | Optional disambiguation / integrity-check columns active on this run (e.g. `["node_name"]` or `["depth"]`) |
| `match_mode` | Chosen match mode: `name_site`, `site_paa_caa_depth`, or `site_paa_caa` |
| `mutations_dropped` | Derived mutations removed under `--mutations-listed-only` because they had no matching CSV row. Always 0 without the flag. |

Counts are scoped to the current run: `nodes_enriched` and `mutations_enriched` exclude pre-existing mutation arrays from upstream pipelines.

`apply_mutations_csv()` reports at status verbosity:

```
Loading mutations CSV: {path}
Loaded {N} CSV rows across {M} families
Match mode: {mode}
Disambiguation columns in CSV: {cols}             # only if any disambig/integrity column is present
Enriched {X} mutations across {Y} nodes in {Z} trees
Unmatched: {U}/{N} CSV rows ({A} in {B} unmatched families, {C} with no node match)
Broadcast: {K} CSV rows matched multiple nodes (ambiguous join — same data applied to every match)
Integrity mismatches: {I} CSV rows matched a (node, site) but disagreed with the tree's derived mutation
```

Warnings (at error level, exit 0):
- **Unmatched mutations** in matched families — the CSV references substitutions / nodes that don't appear in any derived diff
- **Broadcast rows** — the CSV row's enrichment was applied to multiple node-mutations; may not be correct for per-event scores

Hard failures (exit non-zero):
- **Integrity mismatches** in `name_site` mode — the CSV's `parent_aa`/`child_aa`/`depth` disagreed with what the tree derived at that position. Mismatched rows are always skipped (never attached); the default is to also exit non-zero. `--mutations-allow-mismatch` downgrades this to a warning and keeps exit 0.

Per-family detail at `-v 2`: which keys went unmatched, which broadcast, and which failed integrity.

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

## Identifier Conventions

Objects in the output JSON carry up to two identifier fields with
distinct roles. The separation is load-bearing — mixing them produces
format-origin leaks (`"pcp-<uuid>"`) or silently-colliding primary
keys.

| Field | Role | Source | Format |
|---|---|---|---|
| `*_id` — `dataset_id`, `clone_id`, `tree_id`, `sample_id`, `subject_id`, `sequence_id`, `timepoint_id` | Semantic identifier carrying meaning from source data | **Input-derived.** When the input format lacks the field (e.g. PCP has no dataset concept), the CLI synthesizes a typed `{datatype}-{uuid}` value using the same minter as `ident`. | Input-supplied string, or synthesized `{datatype}-{uuid}` |
| `ident` | Primary key the webapp uses for cross-referencing (Redux state, Dexie DB) | **Always CLI-minted.** Never derived from input. | `{datatype}-{uuid}` |

### Rules

- **`*_id` is reserved for input-derived values.** When synthesis is
  unavoidable, use the same `{datatype}-{uuid}` shape as `ident`.
  Never use format-origin prefixes (`pcp-`, `cft-`) — those belong in
  `metadata.source_format` at the output root only.
- **`ident` is reserved for CLI-minted identifiers.** All minting goes
  through `IdentMinter.mint(datatype)` in `olmsted_cli/identifier.py`,
  which enforces the `{datatype}-{uuid}` shape at the signature level.
  Deterministic under `--seed`, random otherwise.
- **`ident` is minted on objects that are — or will become — webapp
  primary keys.** See the table below for per-object status.

| Object | `ident` minted? | Webapp status | Notes |
|---|---|---|---|
| `tree` | yes | **Dexie PK** (`trees.where("ident")`) and used in every lookup path | Load-bearing today |
| `clone` | yes | Redux PK (`clonalFamilies` state keyed on `ident`; `selectedFamily`, starred families flow through it); Dexie still uses compound `[dataset_id+clone_id]` | Half-migrated; webapp side still needs the Dexie PK swap |
| `dataset` | yes | Dexie PK is `dataset_id`; `ident` is written but not yet read | Minted in anticipation of the Dexie PK migration |
| `sample` (PCP only) | yes | Not its own Dexie store today; `sample_id` is indexed on the `clones` store | Minted in anticipation of a `samples` store |
| `subject` | no | No webapp presence | `IdentDatatype` registers the slot for future use |

### Uniqueness guarantees

Validated during processing — duplicates fail fast rather than
silently overwriting downstream.

| Field | Scope |
|---|---|
| `dataset_id` | unique across `datasets[]` |
| `clone_id` | unique within a dataset |
| `tree_id` | unique within a clone |
| `sample_id` | unique within `dataset.samples[]` |
| `subject_id` | unique within `dataset.subjects[]` |
| `sequence_id` | unique within a tree (Newick parser suffix-disambiguates) |

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

_Last updated: 2026-06-09_
