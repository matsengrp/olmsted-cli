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

### Role columns (sample / family / tree)

Three columns identify the row's clonal-family membership. Each role
auto-detects across a small family of names; pass `--<role>-col` to
override.

| Role | Required? | Accepted column names | CLI override |
|------|-----------|----------------------|--------------|
| Sample | required (PCP CSV); optional (trees CSV) | `sample`, `sample_id`, `sample_name` | `--sample-col` |
| Family | required | `family`, `family_id`, `family_name` | `--family-col` |
| Tree | optional | `tree`, `tree_id`, `tree_name` | `--tree-col` |

Preference order when multiple variants are present: `_id` > bare >
`_name`. If two variants appear with conflicting values on the same
row, processing fails with a `Column conflict` error — pass the
override to disambiguate.

The output clone_id is synthesized as `{sample}_{family}` so families
that repeat across samples no longer silently merge.

The tree role enables multiple per-row alternate reconstructions of the
same clonal family. When absent, every `(sample, family)` collapses to a
single tree (today's common case). When present, rows are grouped by
`(sample, family, tree)` and each tree value becomes a separate entry
in `clone.trees[]`.

### Main CSV (required)

Each row represents one parent-child edge in a phylogenetic tree.

**Required columns** (processing will fail without these):

| Column | Description |
|--------|-------------|
| sample role (one of `sample` / `sample_id` / `sample_name`) | Sample identifier |
| family role (one of `family` / `family_id` / `family_name`) | Clonal family identifier |
| `parent_name` | Parent node name |
| `child_name` | Child node name |

**Standard optional columns** (parsed by the standard pipeline):

| Column | Maps to | Level | Notes |
|--------|---------|-------|-------|
| tree role (`tree` / `tree_id` / `tree_name`) | `tree.tree_name` | tree | Disambiguates rows that belong to alternate reconstructions of the same `(sample, family)`. Omit for single-tree-per-family data. |
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
| `cdr1_codon_start_heavy` / `_end` | `cdr1_alignment_start` / `_end`, `cdr1_length` | family | CDR1 positions + length; `cdr1_length` = end − start (nucleotides) |
| `cdr2_codon_start_heavy` / `_end` | `cdr2_alignment_start` / `_end`, `cdr2_length` | family | CDR2 positions + length; `cdr2_length` = end − start (nucleotides) |
| `cdr3_codon_start_heavy` / `_end` | `cdr3_alignment_start` / `_end`, `cdr3_length` | family | CDR3 (junction) positions + length (nucleotides). Despite the `_codon_` column name, the input values are nucleotide coordinates. |
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

Each row provides a Newick tree for one clonal family. Role columns
(sample / family / tree) follow the same auto-detection rules as the
main PCP CSV (see [Role columns](#role-columns-sample--family--tree)
above).

**Required columns**:

| Column | Description |
|--------|-------------|
| family role (`family` / `family_id` / `family_name`) | Clonal family ID (matches the main CSV's family role) |
| `newick_tree` or `newick` | Newick format tree string |

**Optional columns**:

| Column | Maps to | Level |
|--------|---------|-------|
| sample role (`sample` / `sample_id` / `sample_name`) | Composite-key component with family role | — |
| tree role (`tree` / `tree_id` / `tree_name`) | `tree.tree_name`, drives multi-tree-per-clone grouping | tree |
| `reconstruction_method` | `tree.reconstruction_method` | tree |
| `rate_scale_heavy` | `rate_scale_heavy` | family |
| `rate_scale_light` | `rate_scale_light` | family |

**Extra columns**: Any column not listed above is auto-classified by
intra-clone variance. A column whose values vary across the trees of at
least one clone is classified at the **tree** level (drives the tree
dropdown's color/filter/sort controls); a column that is constant
within every clone is classified at the **clone** level. Same chain
suffix convention applies for clone-level extras.

**Multiple trees per family**: Rows are grouped by composite key
`(sample, family, tree)`. Supply a distinct tree-role value
(`tree_name` recommended) on each row that belongs to an alternate
reconstruction of the same `(sample, family)`. Each composite key
becomes a separate entry in `clone.trees[]` and in the top-level
`trees[]`. Without a tree role column, every row attaching to the same
`(sample, family)` is still treated as an alternate reconstruction
(back-compatible with the older "list-per-key" trees CSV convention).

### Minimum Viable PCP Data

To produce a valid tree, you need:
1. `sample_id`, `parent_name`, `child_name` columns
2. At least one parent-child edge with a root node (detected automatically as the node that never appears as a child)
3. Sequence data (`parent_heavy`/`child_heavy`) on at least the root node (required for `mean_mut_freq` calculation)

Missing gene calls, CDR positions, and alignment positions are handled gracefully (empty strings/zeros).

---

## AIRR Input Format

`-f airr` (auto-detected) reads the **AIRR-C v2 Clone/Tree/Node/Cell schema**
— **AIRR Schema v2.0.0** (released 2026-06-05/08), documented at
[docs.airr-community.org/en/latest/datarep/clone.html](https://docs.airr-community.org/en/latest/datarep/clone.html)
— the format Dowser's `writeTreesJSON` emits. It is handled by
`process_airr_data.py` (entry point `process_airr_to_olmsted`).

This is a pragmatic mapping of the schema's concepts onto Olmsted's own JSON
shape, not a byte-for-byte implementation of it — most notably, `Clone.nodes`
here is a **list**, where the official schema's `Tree.nodes` is a **dict
keyed by `sequence_id`**. Don't assume this input round-trips through a
generic AIRR-schema validator unchanged.

(An earlier, Olmsted-flavored `-f airr` container predated this and never
corresponded to an official AIRR release at any version — see
[issue #47](https://github.com/matsengrp/olmsted-cli/issues/47) for the
research and the removal.)

**`mean_mut_freq` computation**: olmsted-cli always computes `mean_mut_freq`
for AIRR input — the Clone schema carries no such field to read or override,
unlike PCP where an upstream-supplied value could in principle disagree with
the recomputed one (see [issue 24](https://github.com/matsengrp/olmsted-cli/issues/24)).
Unlike PCP's multiplicity-weighted convention, AIRR's is necessarily
*unweighted* when the input has no real per-node multiplicity (the `noinfo`
shape, see below) — see [Mapping to Olmsted](#mapping-to-olmsted).

### Top-level structure

A single JSON object with two tables:

```json
{
  "Clone":        [ { "clone_id": "...", "clone_class": "...", "tree": "<newick>",
                      "inferred_ancestor": "...", "nodes": [ ... ] }, ... ],
  "Rearrangement":[ { "sequence_id": "...", "cell_id": "...",
                      "sequence_alignment": "...", "locus": "IGH", ... }, ... ]
}
```

Nodes do **not** carry sequences. Each node points into the `Rearrangement`
table by `sequence_id` (`clone_class: Rearrangement`) or `cell_id`
(`clone_class: Cell`). Inferred internal/germline nodes have their own
synthesized `Rearrangement` records (e.g. `sequence_id: Germline-10004`), so
every node — observed or ASR-inferred — resolves to a sequence.

### Node object (within `Clone.nodes`)

| Field | Role |
|-------|------|
| `node_id` | Node identifier; matches the label in `Clone.tree` (the Newick). |
| `sequence_id` / `cell_id` | Pointer into `Rearrangement` (by class). |
| `node_type` | `observed` or `inferred` (ASR ancestor) — carried onto the output node. |
| `node_class` | `Rearrangement` or `Cell`. |

### Mapping to Olmsted

- **Topology** comes from `Clone.tree` (Newick; its labels are the node ids),
  rerooted so `inferred_ancestor` (the germline) is the root — matching PCP's
  naive-at-root convention. Branch lengths are read from the Newick.
- **Sequences** are joined from `Rearrangement`; `sequence_alignment_aa` is
  translated. `mean_mut_freq` is computed from observed leaves vs the germline.
- **Chains**: `Rearrangement`-class and single-locus `Cell`-class clones → one
  clone + one tree. Paired `Cell` clones (IGH + IGK/IGL) → **two** clones + two
  trees sharing one topology, suffixed `-heavy` / `-light`, each node's sequence
  taken from the same-locus `Rearrangement`.
- **Dataset** is synthesized (like PCP): no subjects/seeds; one sample per
  `repertoire_id`; `--name` sets the dataset name.

Clone-level immunological fields (`v_call`, `j_call`, `cdr3_length`) are passed
through only when the input `Clone` supplies them (the Dowser `info` variant);
the clean v2 (`noinfo`) variant omits them and they are left unset rather than
fabricated.

When the input carries Dowser's `info` catchall (`dowser_fields=TRUE`, Dowser's
default — the `noinfo` variant is written with `dowser_fields=FALSE`), these
are additionally read (see `process_airr_data.py`, issue #45):

- `nodes[].info.tipdata.collapse_count` → node `multiplicity` (unset, not
  fabricated `1`, for nodes/inputs without it — e.g. every node in `noinfo`
  input, and inferred/ASR nodes even in `info` input, which carry `info: []`).
- `Clone.info.region` (a per-position IMGT region label array, indexed
  against the *ungapped* `Rearrangement.sequence`) → `cdr1_alignment_start`/
  `_end`, `cdr2_alignment_start`/`_end`, `cdr3_alignment_start`/`_end`, and
  the matching `cdr{1,2,3}_length` (0-based, half-open, nucleotide positions
  in the *gapped* `germline_alignment` — the same convention as PCP's
  `cdr*_alignment_start`/`_end`). `region`'s own `cdr3` span is
  **strict IMGT CDR3**, which excludes the 2 conserved anchor residues
  (V-gene 2nd-CYS, J-gene TRP/PHE) that "junction" includes — per the
  IMGT/AIRR Community convention, junction is exactly those 2 residues (1
  codon = 3 nucleotides each) longer on both ends. Since `cdr3_length` is a
  synonym for junction length everywhere else in this project (PCP,
  `schemas.py`), the derived `cdr3` span is normalized to the junction
  convention (extended by 1 anchor codon each side, in ungapped coordinates,
  before remapping to gapped ones) — resolving the inconsistency issue #46
  raised. When `junction_length` is also available, it's now expected to
  agree with the normalized `cdr3_length` exactly; a real disagreement (not
  the already-accounted-for CDR3/junction anchor difference) logs a warning
  and keeps the `region`-derived value. `cdr1`/`cdr2` have no such
  distinction and are used as `region` gives them.
  **Not yet handled for paired (heavy+light) clones**: `region` covers both
  chains concatenated, which doesn't match either chain's own
  `germline_alignment` length, so it's safely skipped rather than
  misattributed.

Clone-level `v_call`/`j_call` fall back to the germline/root node's own
`Rearrangement` record when the `Clone` doesn't supply them (cf. #24);
`d_call` — which the `Clone`-level `info` catchall never carries at all —
comes only from that fallback. This works for both the `info` and `noinfo`
variants, since it reads the `Rearrangement` table directly rather than the
`info` catchall.

Still deferred (#45): `program_origin`, arbitrary Dowser `trait=` columns in
per-node tipdata, and **arbitrary custom fields generally** — unlike PCP and
Olmsted JSON, this format's clone/node dicts are built from a fixed field
set, so there's no equivalent to those formats' `--capture-all`/custom-field
passthrough (see `example-data/fields-config/README.md`).

**Streaming** (`--batch-size`) is not supported for this format — it always
runs in-memory regardless of `--batch-size`, unlike PCP (see issue #36).

See `example-data/airr/` for `nocell`/`unpaired`/`paired` inputs + goldens,
in both the `noinfo` and `info` flavors.

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
| **tree** | `tree` | Per-tree refs in `clone.trees[]` | Tree-dropdown color/filter/sort (multi-tree clones) |
| **node** | `node` | Tree node objects | Tree node tooltips, properties |
| **branch** | `branch` | Branch length on nodes | Tree branch coloring, width |
| **mutation** | `mutation` | `mutations[]` on nodes, or derived by web app | Alignment mutation coloring |

A clone-level extra is auto-promoted to **tree** level when its value
varies across the trees of at least one clone — single-tree-per-clone
datasets classify everything as clone-level, preserving today's output
shape.

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
| Gene calls (`v_call`, `d_call`, `j_call`) | Empty string; locus is inferred from V-gene prefix when possible, else left unset |
| CDR/alignment positions | Zero values |
| Tree file (PCP) | Trees built from parent-child edges |
| `tree.reconstruction_method` | Left unset when not supplied by input |
| `tree.type` / `dataset.type` | Not synthesized — passed through from input only |
| Node `length` / `distance` | `olmsted merge` backfills these from the tree's `newick` branch lengths when absent (no-clobber); `olmsted validate` warns when a branch-length newick has nodes that lack them |

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
| JSON with top-level `Clone` + `Rearrangement` keys | `airr` |
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

| Source | Output Location | Output Field |
|--------|-----------------|-------------|
| `Clone.tree` (Newick, rerooted on `inferred_ancestor`) | tree | `newick`, node `parent`/`length`/`distance` |
| `Rearrangement.sequence_alignment` (joined via node) | node | `sequence_alignment` |
| (translated) | node | `sequence_alignment_aa` |
| `Rearrangement.v_call`/`d_call`/`j_call` (germline fallback) | clone | `v_call`/`d_call`/`j_call` |
| `nodes[].node_type` | node | `node_type` (`observed`/`inferred`) |
| (computed) | clone | `mean_mut_freq` |
| `Clone.info.tipdata.collapse_count` (when present) | node | `multiplicity` |
| `Clone.info.region` (when present, junction-normalized) | clone | `cdr{1,2,3}_alignment_start`/`_end`/`_length` |
| Paired `Cell` clones (IGH + IGK/IGL) | 2 clones + 2 trees | `-heavy`/`-light` suffix, per-locus sequences |
| (synthesized) | dataset | one sample per `repertoire_id` |

See [Mapping to Olmsted](#mapping-to-olmsted) in the AIRR section above for the full detail.

---

_Last updated: 2026-08-23_
