# File Format Reference

This document describes the input and output formats for olmsted-cli: what data is expected, where fields are found, how they map to the output, and what is enforced by validation.

See also:
- **[README.md](./README.md)**: Command usage and examples
- **[ARCHITECTURE.md](./ARCHITECTURE.md)**: Processing pipeline and data flow
- **[DEVELOPMENT.md](./DEVELOPMENT.md)**: Setup guide, testing, how-to guides
- **[CLAUDE.md](./CLAUDE.md)**: Code quality rules, terminology

---

## Table of Contents

- [PCP Input Format](#pcp-input-format)
- [AIRR Input Format](#airr-input-format)
- [Olmsted JSON Output Format](#olmsted-json-output-format)
- [Field Metadata](#field-metadata)
- [Validation](#validation)
- [Field Mapping: Input → Output](#field-mapping-input--output)

---

## PCP Input Format

PCP (Parent-Child Pair) format uses one or two CSV files.

### Main CSV (required)

Each row represents one parent-child edge in a phylogenetic tree.

**Required columns** (processing will fail without these):

| Column | Description |
|--------|-------------|
| `sample_id` | Sample/dataset identifier |
| `parent_name` | Parent node name |
| `child_name` | Child node name |

**Standard optional columns** (parsed by the standard pipeline):

| Column | Maps to | Level | Notes |
|--------|---------|-------|-------|
| `family` | `clone_id` | family | Defaults to `sample_id` if absent |
| `parent_heavy` / `child_heavy` | `sequence_alignment` | node | Heavy chain DNA sequences |
| `parent_light` / `child_light` | `sequence_alignment` (light clone) | node | Light chain DNA sequences (paired data) |
| `branch_length` or `edge_length` | `length` | branch | Branch length |
| `distance` | `distance` | node | Cumulative distance from root |
| `depth` | — | — | Depth in tree (not carried to output) |
| `sample_count` | `multiplicity` | node | Sequence abundance |
| `v_gene_heavy` | `v_call` | family | V gene assignment |
| `d_gene_heavy` | `d_call` | family | D gene assignment |
| `j_gene_heavy` | `j_call` | family | J gene assignment |
| `v_gene_light` / `j_gene_light` | `v_call_light` / `j_call_light` | family | Light chain gene calls |
| `cdr1_codon_start_heavy` / `_end` | `cdr1_alignment_start` / `_end` | family | CDR1 positions |
| `cdr2_codon_start_heavy` / `_end` | `cdr2_alignment_start` / `_end` | family | CDR2 positions |
| `cdr3_codon_start_heavy` / `_end` | `junction_start` / `junction_length` | family | CDR3/junction positions |
| `parent_is_naive` | Node `type: "root"` | node | Boolean |
| `child_is_leaf` | Node `type: "leaf"` | node | Boolean |
| `light_chain_type` | `light_chain_type` | family | `"kappa"` or `"lambda"` |

**Column aliases** (automatically mapped to canonical names):

| Alias | Maps to |
|-------|---------|
| `v_gene`, `v_call` | `v_gene_heavy` |
| `d_gene`, `d_call` | `d_gene_heavy` |
| `j_gene`, `j_call` | `j_gene_heavy` |
| `parent_seq`, `parent_sequence` | `parent_heavy` |
| `child_seq`, `child_sequence` | `child_heavy` |

**Extra columns**: Any column not listed above is captured as a custom node-level field on the child node. Values are auto-coerced: int → float → JSON → bool → string.

**Chain suffix convention** (paired data): Columns ending with `_heavy` are applied only to the heavy chain clone/nodes. Columns ending with `_light` are applied only to the light chain clone/nodes. Columns without a suffix are shared between both chains. Suffixes are stripped in the output (e.g., `score_heavy` → `score`).

### Trees CSV (optional)

Each row provides a Newick tree for one clonal family.

**Required columns** (one of each pair):

| Column | Description |
|--------|-------------|
| `family_name` or `family` | Clonal family ID (matches `family` in main CSV) |
| `newick_tree` or `newick` | Newick format tree string |

**Optional columns**:

| Column | Maps to | Level |
|--------|---------|-------|
| `sample_id` | Used as composite key with `family_name` | — |
| `rate_scale_heavy` | `rate_scale_heavy` | family |
| `rate_scale_light` | `rate_scale_light` | family |

**Extra columns**: Any column not listed above is captured as a family-level (clone-level) field. Same chain suffix convention applies.

### Minimum Viable PCP Data

To produce a valid tree, you need:
1. `sample_id`, `parent_name`, `child_name` columns
2. At least one parent-child edge with a root node (detected automatically as the node that never appears as a child)
3. Sequence data (`parent_heavy`/`child_heavy`) on at least the root node (required for `mean_mut_freq` calculation)

Missing gene calls, CDR positions, and alignment positions are handled gracefully (empty strings/zeros).

---

## AIRR Input Format

AIRR (Adaptive Immune Receptor Repertoire) format is a single JSON file following the [AIRR Community standards](https://docs.airr-community.org/).

### Top-level structure

```json
{
  "dataset_id": "...",           // Required
  "ident": "...",
  "subjects": [...],
  "samples": [...],
  "seeds": [...],
  "clones": [...]                // Required: array of clone objects
}
```

### Clone object

| Field | Required | Maps to | Notes |
|-------|----------|---------|-------|
| `clone_id` | Yes | `clone_id` | Unique family identifier |
| `sample_id` | Yes | `sample_id` | Links to samples array |
| `subject_id` | No | `subject_id` | Defaults to `"unknown"` |
| `v_call` | No | `v_call` | V gene assignment |
| `d_call` | No | `d_call` | D gene assignment |
| `j_call` | No | `j_call` | J gene assignment |
| `germline_alignment` | No | `germline_alignment` | Germline sequence |
| `unique_seqs_count` | Yes* | `unique_seqs_count` | Schema-required |
| `mean_mut_freq` | Yes* | `mean_mut_freq` | Schema-required |
| `v_alignment_start` | No | `v_alignment_start` | 1-based → 0-based conversion |
| `d_alignment_start` | No | `d_alignment_start` | 1-based → 0-based conversion |
| `j_alignment_start` | No | `j_alignment_start` | 1-based → 0-based conversion |
| `junction_start` | No | `junction_start` | 1-based → 0-based conversion |
| `junction_length` | No | `junction_length` | |
| `trees` | Yes | `trees` | Array of tree objects |

*Required by schema validation, but processing won't crash without them.

**Extra fields**: Any additional fields on clone objects are preserved in the output and auto-detected by field_metadata generation.

### Tree object (within clone)

| Field | Required | Notes |
|-------|----------|-------|
| `newick` | Yes | Newick tree string (schema-required) |
| `tree_id` | No | Tree identifier |
| `nodes` | Yes | Dict or array of node objects |

### Node object (within tree)

| Field | Required | Notes |
|-------|----------|-------|
| `sequence_id` | Yes | Node identifier (schema-required) |
| `sequence_alignment` | Yes | DNA sequence (schema-required) |
| `sequence_alignment_aa` | Yes | Amino acid sequence (schema-required) |
| `parent` | No | Parent node ID (null for root) |
| `type` | No | `"root"`, `"internal"`, or `"leaf"` |
| `length` | No | Branch length to parent |
| `distance` | No | Cumulative distance from root |
| `multiplicity` | No | Sequence abundance |
| `lbi`, `lbr` | No | Computed with `--compute-metrics` |
| `affinity` | No | |
| `timepoint_id` | No | Sampling timepoint |
| `mutations` | No | Array of per-mutation records |

**Extra fields**: Preserved in output and auto-detected.

### AIRR Position Convention

AIRR uses 1-based closed intervals. olmsted-cli converts `*_start` positions to 0-based (subtracting 1) during processing. Missing positions are skipped gracefully.

---

## Olmsted JSON Output Format

The consolidated output format produced by `olmsted process`.

### Top-level structure

```json
{
  "metadata": {
    "format": "olmsted",
    "format_version": "1.0",
    "schema_version": "2.0.0",
    "created_at": "ISO 8601 timestamp",
    "source_format": "pcp" | "airr",
    "source_files": ["filename.csv"],
    "processing_info": {
      "datasets_count": 1,
      "total_clones_count": 8,
      "total_trees_count": 8,
      "total_leaf_nodes_count": 186
    },
    "generated_by": {
      "tool": "olmsted-cli",
      "version": "0.2.0",
      "git_hash": "abc1234"
    },
    "name": "My Dataset",
    "processing_options": { ... }
  },
  "datasets": [ ... ],
  "clones": { "dataset_id": [ ... ] },
  "trees": [ ... ]
}
```

### Dataset object

| Field | Description |
|-------|-------------|
| `dataset_id` | Unique identifier (required) |
| `name` | User-provided name |
| `clone_count` | Number of clonal families |
| `field_metadata` | Describes available fields (see below) |
| `subjects`, `samples`, `timepoints` | Metadata arrays |

### Clone object (in `clones[dataset_id]`)

Contains all family-level data: gene calls, alignment positions, CDR positions, gene support probabilities, and any extra fields from input data.

### Tree object (in `trees[]`)

| Field | Description |
|-------|-------------|
| `ident` | Unique tree identifier |
| `clone_id` | Links to parent clone |
| `newick` | Newick tree string |
| `nodes` | Array of node objects |

### Node object (in tree `nodes[]`)

Contains sequence data, tree topology (`parent`, `type`), metrics (`lbi`, `lbr`, `affinity`), and any extra fields from input data.

---

## Field Metadata

The `field_metadata` object on each dataset describes available data fields for the web app's visualization controls.

### Structure

```json
{
  "field_metadata": {
    "clone": {
      "unique_seqs_count": {
        "type": "continuous",
        "display": "dropdown",
        "label": "Unique Sequences Count"
      },
      "v_call": {
        "type": "categorical",
        "display": "dropdown",
        "label": "V Gene"
      }
    },
    "node": { ... },
    "branch": { ... },
    "mutation": {
      "child_aa": {
        "type": "aa",
        "display": "dropdown",
        "label": "Child Amino Acid"
      },
      "parent_aa": {
        "type": "aa",
        "display": "tooltip",
        "label": "Parent Amino Acid"
      },
      "surprise_mutsel": {
        "type": "continuous",
        "display": "dropdown",
        "label": "Surprise (MutSel)",
        "range": [0.68, 13.03]
      }
    }
  }
}
```

### Entry fields

| Key | Values | Description |
|-----|--------|-------------|
| `type` | `continuous`, `categorical`, `aa`, `dna` | What the data is |
| `display` | `dropdown`, `tooltip` | How the web app uses it |
| `label` | String | Human-readable display name |
| `range` | `[min, max]` | Value bounds (continuous mutation fields) |

### How fields are classified

1. **Known field registries** (`KNOWN_CLONE_FIELDS`, etc. in `constants.py`) — matched by name, provides type/display/label
2. **Auto-inference** — values sampled from data, type inferred (numeric → continuous, string → categorical, single-char AA → aa, single-char DNA → dna)
3. **Suggestions** (`SUGGESTED_SKIP_FIELDS`, `SUGGESTED_DISPLAY_MODES`) — removes non-visualization fields, overrides display for context fields
4. **Custom fields** (from YAML config) — overrides everything above

### Levels

| Level | Internal key | Source data | Used for |
|-------|-------------|-------------|----------|
| **family** | `clone` | Clone/family objects | Scatterplot axes, color, shape, facet |
| **node** | `node` | Tree node objects | Tree node tooltips, properties |
| **branch** | `branch` | Branch length on nodes | Tree branch coloring, width |
| **mutation** | `mutation` | `mutations[]` on nodes, or derived by web app | Alignment mutation coloring |

### Derived fields

When nodes have `sequence_alignment_aa`, the mutation level includes `child_aa` and `parent_aa` even though they aren't in the data — the web app derives them at render time by diffing parent/child sequences.

---

## Validation

### Schema-enforced (required fields)

| Object | Required fields |
|--------|----------------|
| **Dataset** | `dataset_id` |
| **Clone** | `unique_seqs_count`, `mean_mut_freq` |
| **Tree** | `newick` |
| **Node** | `sequence_id`, `sequence_alignment`, `sequence_alignment_aa` |

All schemas allow `additionalProperties: true` — extra fields are preserved.

### Gracefully handled when missing

| Missing field | Behavior |
|---------------|----------|
| `v/d/j_alignment_start` | Skipped (not adjusted), notification at verbose ≥ 2 |
| `subject_id` | Defaults to `"unknown"` |
| `timepoint_id` (on samples) | Defaults to `"unknown"` |
| Sample not found for `sample_id` | Default sample created, notification at verbose ≥ 1 |
| Gene calls (`v_call`, `d_call`, `j_call`) | Empty string, locus defaults to `"igh"` |
| CDR/alignment positions | Zero values |
| Tree file (PCP) | Trees built from parent-child edges |

### Format detection

`detect_file_format()` identifies input format:

| Check | Result |
|-------|--------|
| `.csv` extension | `pcp` |
| JSON with `metadata.format == "olmsted"` | `olmsted` (explicit tag) |
| JSON with `datasets` + `metadata` keys | `olmsted` (heuristic) |
| JSON with `dataset_id` or `clones` key | `airr` |
| Otherwise | `unknown` |

---

## Field Mapping: Input → Output

### PCP → Olmsted JSON

| PCP Column | Output Location | Output Field |
|------------|-----------------|-------------|
| `sample_id` | clone | `sample_id` |
| `family` | clone | `clone_id` |
| `parent_heavy`/`child_heavy` | node | `sequence_alignment` |
| (translated) | node | `sequence_alignment_aa` |
| `branch_length` | node | `length` |
| `distance` | node | `distance` |
| `v_gene_heavy` | clone | `v_call` |
| `d_gene_heavy` | clone | `d_call` |
| `j_gene_heavy` | clone | `j_call` |
| `cdr1_codon_start_heavy` | clone | `cdr1_alignment_start` |
| `parent_is_naive` | node | `type: "root"` |
| `child_is_leaf` | node | `type: "leaf"` |
| (computed) | clone | `mean_mut_freq` |
| (computed) | clone | `unique_seqs_count` |
| (tree CSV extra cols) | clone | (field name preserved) |
| (PCP CSV extra cols) | node | (field name preserved) |

### AIRR → Olmsted JSON

AIRR fields are mostly passed through directly. Key transformations:

| Transformation | Description |
|---------------|-------------|
| `*_start` positions | 1-based → 0-based (subtract 1) |
| `dataset.samples` | Nested as `clone.sample` |
| `dataset` (minus clones) | Nested as `clone.dataset` |
| Tree nodes | Extracted from clones, stored in top-level `trees[]` |
| `clone.trees` | Reduced to metadata references (nodes removed) |

---

_Last updated: 2026-03-31_
