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
- [Mutations CSV Format](#mutations-csv-format)
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
| `tree_id` | `tree.tree_id` | tree |
| `reconstruction_method` | `tree.reconstruction_method` | tree |
| `rate_scale_heavy` | `rate_scale_heavy` | family |
| `rate_scale_light` | `rate_scale_light` | family |

**Extra columns**: Any column not listed above is captured as a family-level (clone-level) field. Same chain suffix convention applies.

**Multiple trees per family**: The same `(family_name, sample_id)` may appear on multiple rows, one per alternate phylogenetic reconstruction of the family. Each row becomes a separate entry in `clone.trees[]` and in the top-level `trees[]`. Supply a distinct `tree_id` on each row to label the alternate — duplicates within a `(family, sample_id)` pair cause the output uniqueness check to fail (see `--allow-duplicate-ids` to opt out).

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

## Mutations CSV Format

External mutation-level annotations consumed by the `merge` command and the `process --mutations` flag. Each row describes one substitution.

### Required columns

| Column | Description |
|--------|-------------|
| `family` | Clonal family identifier — joined against `clone_id` in the Olmsted JSON |
| `site` | Integer amino acid position (0-based, matching the `sequence_alignment_aa` index) |
| `parent_aa` | Single-character parent amino acid |
| `child_aa` | Single-character child amino acid |

If `site` is non-numeric, parsing fails with a clear `ValueError` pointing at the offending row.

### Recognized structural columns

These are read for context but **not** added as mutation-level fields on the output:

| Column | Purpose |
|--------|---------|
| `sample_id` | Optional sample identifier for cross-checking (not enforced) |
| `pcp_index` | Optional integer index into the source PCP CSV |
| `depth` | Optional tree depth where the mutation occurs |

### Score columns

Any column not listed above becomes a mutation-level field on matching nodes. Common examples produced by upstream pipelines:

| Column | Output type | Auto-detected label |
|--------|-------------|---------------------|
| `surprise_mutsel` | continuous | Surprise (MutSel) |
| `surprise_neutral` | continuous | Surprise (Neutral) |
| `surprise_mutsel_theoretical` | continuous | Surprise (MutSel, Theoretical) |
| `selection_contribution` | continuous | Selection Contribution |
| `log_selection_factor` | continuous | Log Selection Factor |
| `num_codon_changes` | continuous | Number of Codon Changes |

These are pre-registered in `KNOWN_MUTATION_FIELDS`. Any other score column will be auto-detected (continuous if numeric, categorical if string) and appear in `field_metadata.mutation` with a generated label.

### Matching semantics

For each tree whose `clone_id` matches a CSV `family`:

1. The CSV rows for that family are indexed by `(site, parent_aa, child_aa)`.
2. For each tree node, mutations are derived by diffing `node.sequence_alignment_aa` against its parent's (or read directly if a `mutations` array already exists). Gap characters (`-`, `.`, `X`, `*`, `?`) are skipped.
3. Each derived mutation is looked up in the CSV index. On match, the score columns are merged onto the mutation dict.

### Unmatched rows

Rows whose `family` doesn't appear in the JSON, or whose `(site, parent_aa, child_aa)` doesn't appear on any node in the matched tree, are reported as warnings. The merge still completes; warnings include counts at normal verbosity and per-family detail at `-v 2`.

### Example

```csv
family,site,parent_aa,child_aa,surprise_mutsel,selection_contribution,sample_id,depth
clone-abc,12,K,R,4.21,0.77,s1,3
clone-abc,57,A,T,3.06,1.31,s1,3
clone-xyz,9,G,D,5.21,0.94,s1,4
```

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
| `ident` | CLI-minted primary key (`tree-{uuid}`) |
| `tree_id` | Semantic identifier: from PCP trees.csv `tree_id` column if present, otherwise synthesized as `tree-{family_id}` (paired: `-heavy` / `-light` suffix). AIRR: passed through from input, or falls back to `ident`. |
| `clone_id` | Links to parent clone |
| `reconstruction_method` | *(optional)* Method used to build the tree (e.g. `"dnapars"`, `"raxml_ng"`). Only present when the input provided one; absent means unknown. |
| `newick` | Newick tree string |
| `nodes` | Array of node objects |

A clone can carry multiple alternate-reconstruction trees in its `clone.trees[]` list; each gets its own entry in the top-level `trees[]` with a full `nodes` array.

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
      "selection_contribution": {
        "type": "continuous",
        "display": "dropdown",
        "label": "Selection Contribution",
        "range": [-2.5, 5.1]
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
| `subject_id` | Left unset; webapp renders its own unknown marker |
| `timepoint_id` | Left unset; webapp renders its own unknown marker |
| Sample not found for `sample_id` (AIRR) | `clone.sample` left unset, notification at verbose ≥ 1 |
| Gene calls (`v_call`, `d_call`, `j_call`) | Empty string; locus is inferred from V-gene prefix when possible, else left unset |
| CDR/alignment positions | Zero values |
| Tree file (PCP) | Trees built from parent-child edges |
| `tree.reconstruction_method` | Left unset when not supplied by input |
| `tree.type` / `dataset.type` | Not synthesized — passed through from input only |

### Uniqueness enforcement

`*_id` fields that the webapp uses to cross-reference objects must be unique within their natural scope. `olmsted process`, `olmsted tag`, and `olmsted merge` all check this before writing output and fail fast on collisions:

| Scope | Field |
|---|---|
| Within output | `dataset.dataset_id` |
| Within a dataset | `clone.clone_id` |
| Within a clone | `tree.tree_id` |
| Within `dataset.samples[]` | `sample.sample_id` |
| Within `dataset.subjects[]` | `subject.subject_id` |

Pass `--allow-duplicate-ids` to downgrade these to warnings and let the data pass through unchanged. `sequence_id` uniqueness within a tree is always enforced upstream by the Newick parser (duplicate leaf names get `_1`, `_2` suffixes).

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
| Matching `dataset.samples[]` entry | Denormalized as `clone.sample` (webapp reads `clone.sample.locus` etc.) |
| Tree nodes | Extracted from clones, stored in top-level `trees[]` |
| `clone.trees` | Reduced to metadata references (nodes removed) |
| `tree.tree_id` | Passed through from input; when absent, filled with the CLI-minted `tree.ident` to satisfy the AIRR Community schema's required-field contract. |

---

_Last updated: 2026-04-22_
