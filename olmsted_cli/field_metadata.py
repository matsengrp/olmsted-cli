"""
Field metadata generation for Olmsted datasets.

Generates field_metadata describing available data fields at each level
(clone, node, branch, mutation) with their types and human-readable labels.
This metadata drives dynamic dropdown construction in the Olmsted web app.

Field types:
    - "continuous": Numeric values suitable for axes, size, color scales
    - "categorical": String/enum values suitable for color, shape, facet
    - "tooltip": Display-only values shown in tooltips, not for encoding

Levels:
    - "clone": Clone/clonal family level (scatterplot axes, color, facet)
    - "node": Tree node level (node properties, tooltips)
    - "branch": Tree branch level (branch coloring, width)
    - "mutation": Per-mutation level (alignment coloring)
"""

from __future__ import annotations

import re
from typing import Any, Dict, List, Optional


# =============================================================================
# Static registries of known fields at each level
# =============================================================================

KNOWN_CLONE_FIELDS = {
    "unique_seqs_count": {"type": "continuous", "label": "Unique Sequences Count"},
    "total_read_count": {"type": "continuous", "label": "Total Read Count"},
    "mean_mut_freq": {"type": "continuous", "label": "Mean Mutation Frequency"},
    "junction_length": {"type": "continuous", "label": "Junction Length"},
    "clone_count": {"type": "continuous", "label": "Clone Count"},
    "v_call": {"type": "categorical", "label": "V Gene"},
    "d_call": {"type": "categorical", "label": "D Gene"},
    "j_call": {"type": "categorical", "label": "J Gene"},
    "locus": {"type": "categorical", "label": "Locus"},
    "subject_id": {"type": "categorical", "label": "Subject"},
    "sample_id": {"type": "categorical", "label": "Sample"},
    "has_seed": {"type": "categorical", "label": "Has Seed"},
    "is_paired": {"type": "categorical", "label": "Is Paired"},
    "light_chain_type": {"type": "categorical", "label": "Light Chain Type"},
    "v_call_light": {"type": "categorical", "label": "V Gene (Light)"},
    "j_call_light": {"type": "categorical", "label": "J Gene (Light)"},
    "rate_scale_heavy": {"type": "continuous", "label": "Rate Scale (Heavy)"},
    "rate_scale_light": {"type": "continuous", "label": "Rate Scale (Light)"},
}

KNOWN_NODE_FIELDS = {
    "lbi": {"type": "continuous", "label": "LBI"},
    "lbr": {"type": "continuous", "label": "LBR"},
    "multiplicity": {"type": "continuous", "label": "Multiplicity"},
    "cluster_multiplicity": {"type": "continuous", "label": "Cluster Multiplicity"},
    "affinity": {"type": "continuous", "label": "Affinity"},
    "scaled_affinity": {"type": "continuous", "label": "Scaled Affinity"},
    "relative_affinity": {"type": "continuous", "label": "Relative Affinity"},
    "distance": {"type": "continuous", "label": "Distance from Root"},
    "subtree_size": {"type": "continuous", "label": "Subtree Size"},
    "count": {"type": "continuous", "label": "Count"},
    "confidence": {"type": "continuous", "label": "Confidence"},
    "timepoint_id": {"type": "categorical", "label": "Timepoint"},
    "affinity_class": {"type": "categorical", "label": "Affinity Class"},
}

KNOWN_BRANCH_FIELDS = {
    "length": {"type": "continuous", "label": "Branch Length"},
    "branch_length": {"type": "continuous", "label": "Branch Length"},
}

KNOWN_MUTATION_FIELDS = {
    "surprise_mutsel": {"type": "continuous", "label": "Surprise (MutSel)"},
    "surprise_neutral": {"type": "continuous", "label": "Surprise (Neutral)"},
    "selection_contribution": {
        "type": "continuous",
        "label": "Selection Contribution",
    },
    "region": {"type": "categorical", "label": "Region"},
    "parent_aa": {"type": "tooltip", "label": "Parent Amino Acid"},
    "child_aa": {"type": "aa", "label": "Child Amino Acid"},
    "parent_nt": {"type": "tooltip", "label": "Parent Nucleotide"},
    "child_nt": {"type": "dna", "label": "Child Nucleotide"},
}

# Mapping from level name to its known fields registry
KNOWN_FIELDS_BY_LEVEL = {
    "clone": KNOWN_CLONE_FIELDS,
    "node": KNOWN_NODE_FIELDS,
    "branch": KNOWN_BRANCH_FIELDS,
    "mutation": KNOWN_MUTATION_FIELDS,
}


# =============================================================================
# Cross-format field aliases
# =============================================================================
#
# Maps alternative field names (from different formats/conventions) to their
# canonical Olmsted output name. Used during field_metadata generation and
# data renaming when output_name is not explicitly specified.
#
# Format: {input_name: canonical_output_name}

FIELD_ALIASES = {
    # Gene call naming: AIRR uses v_call, PCP uses v_gene_heavy
    "v_gene": "v_call",
    "v_gene_heavy": "v_call",
    "d_gene": "d_call",
    "d_gene_heavy": "d_call",
    "j_gene": "j_call",
    "j_gene_heavy": "j_call",
    "v_gene_light": "v_call_light",
    "j_gene_light": "j_call_light",
    # Sequence count naming
    "rearrangement_count": "unique_seqs_count",
    "sampled_seqs_count": "unique_seqs_count",
    # Size/count
    "size": "total_read_count",
    # Tree metrics
    "branch_length": "length",
    # Mutation fields
    "mut_to": "child_aa",
    "mut_from": "parent_aa",
}


# =============================================================================
# Exclusion lists — structural/internal fields not useful for visualization
# =============================================================================

EXCLUDED_CLONE_FIELDS = {
    # Identifiers and structural refs
    "ident",
    "clone_id",
    "dataset_id",
    "dataset",
    "sample",
    "trees",
    "repertoire_id",
    "pair_id",
    "seed_id",
    "unique_ids",
    "timepoint_ids",
    "type",
    "schema_version",
    # Sequences (long strings, not for dropdowns)
    "germline_alignment",
    "germline_sequence",
    "germline_alignment_light",
    "naive_sequence",
    "cdr3_sequence",
    # Alignment position fields
    "v_alignment_start",
    "v_alignment_end",
    "v_sequence_start",
    "v_sequence_end",
    "v_germline_start",
    "v_germline_end",
    "d_alignment_start",
    "d_alignment_end",
    "d_sequence_start",
    "d_sequence_end",
    "d_germline_start",
    "d_germline_end",
    "j_alignment_start",
    "j_alignment_end",
    "j_sequence_start",
    "j_sequence_end",
    "j_germline_start",
    "j_germline_end",
    "cdr1_alignment_start",
    "cdr1_alignment_end",
    "cdr2_alignment_start",
    "cdr2_alignment_end",
    "cdr1_alignment_start_light",
    "cdr1_alignment_end_light",
    "cdr2_alignment_start_light",
    "cdr2_alignment_end_light",
    "junction_start",
    "junction_end",
    "junction_start_light",
    "junction_length_light",
    # Gene support probability arrays
    "v_per_gene_support",
    "d_per_gene_support",
    "j_per_gene_support",
    # Build/metadata
    "build",
    # Tree reference metadata
    "trees_meta",
}

EXCLUDED_NODE_FIELDS = {
    # Identifiers and structural
    "sequence_id",
    "node_id",
    "parent",
    "is_root",
    "type",
    # Sequences (long strings)
    "sequence_alignment",
    "sequence_alignment_aa",
    "sequence_alignment_light",
    "sequence_alignment_light_aa",
    "aa_sequence",
    "junction",
    "junction_aa",
    # The surprise_mutations array itself (sub-fields are mutation-level)
    "surprise_mutations",
    # Multiplicity arrays (complex objects, not scalar)
    "timepoint_multiplicities",
    "cluster_timepoint_multiplicities",
    # Gene calls on nodes (usually same as clone-level)
    "v_call",
    "d_call",
    "j_call",
    # Sample/dataset refs
    "sample_id",
    "timepoint",
}

EXCLUDED_BRANCH_FIELDS = {
    "sequence_id",
    "node_id",
    "parent",
}

EXCLUDED_MUTATION_FIELDS = {
    "site",
}

EXCLUDED_FIELDS_BY_LEVEL = {
    "clone": EXCLUDED_CLONE_FIELDS,
    "node": EXCLUDED_NODE_FIELDS,
    "branch": EXCLUDED_BRANCH_FIELDS,
    "mutation": EXCLUDED_MUTATION_FIELDS,
}


# =============================================================================
# Abbreviation map for label generation
# =============================================================================

ABBREVIATION_MAP = {
    "lbi": "LBI",
    "lbr": "LBR",
    "cdr": "CDR",
    "cdr1": "CDR1",
    "cdr2": "CDR2",
    "cdr3": "CDR3",
    "shm": "SHM",
    "aa": "AA",
    "dna": "DNA",
    "id": "ID",
    "v": "V",
    "d": "D",
    "j": "J",
    "mut": "Mutation",
    "freq": "Frequency",
    "seq": "Sequence",
    "seqs": "Sequences",
    "mutsel": "MutSel",
}


# =============================================================================
# Utility functions
# =============================================================================


_AA_CHARS = set("ACDEFGHIKLMNPQRSTVWY*-X")
_DNA_CHARS = set("ACGTURYSWKMBDHVN-.")


def infer_field_type(values: List[Any]) -> str:
    """
    Infer field type from sample values.

    Args:
        values: Non-null sample values from the field.

    Returns:
        "continuous" if all values are numeric,
        "aa" if all values are single amino acid characters,
        "dna" if all values are single nucleotide characters,
        "categorical" if all values are strings,
        "tooltip" if mixed types or unclassifiable.
    """
    if not values:
        return "tooltip"

    numeric_count = 0
    string_count = 0
    string_values = []
    for v in values:
        if isinstance(v, bool):
            string_count += 1
        elif isinstance(v, (int, float)):
            numeric_count += 1
        elif isinstance(v, str):
            string_count += 1
            string_values.append(v)
        else:
            # Complex types (lists, dicts) are tooltip-only
            return "tooltip"

    if numeric_count > 0 and string_count == 0:
        return "continuous"
    if string_count > 0 and numeric_count == 0:
        # Check for single-character DNA or AA
        # DNA checked first: if values contain AA-only chars (e.g., D, E, F),
        # they can't be DNA. Pure ACGTU ambiguity is resolved by the known
        # fields registry (parent_aa/child_aa vs parent_nt/child_nt).
        if string_values and all(len(s) == 1 for s in string_values):
            upper_vals = {s.upper() for s in string_values}
            # If any char is AA-only (not in DNA alphabet), it's AA
            if upper_vals <= _AA_CHARS and not upper_vals <= _DNA_CHARS:
                return "aa"
            if upper_vals <= _DNA_CHARS:
                return "dna"
            if upper_vals <= _AA_CHARS:
                return "aa"
        return "categorical"
    return "tooltip"


def compute_range(dicts: List[Dict], field: str) -> Optional[List[float]]:
    """
    Compute [min, max] range for a numeric field across a list of dicts.

    Returns None if no numeric values are found.
    """
    values = []
    for d in dicts:
        v = d.get(field)
        if isinstance(v, (int, float)) and not isinstance(v, bool):
            values.append(v)
    if not values:
        return None
    return [min(values), max(values)]


def humanize_label(field_name: str) -> str:
    """
    Convert a snake_case field name to a human-readable label.

    Uses ABBREVIATION_MAP for known terms, title-cases the rest.

    Examples:
        "unique_seqs_count" -> "Unique Sequences Count"
        "lbi" -> "LBI"
        "mean_mut_freq" -> "Mean Mutation Frequency"
        "mean_surprise_mutsel" -> "Mean Surprise MutSel"
    """
    parts = field_name.split("_")
    result = []
    for part in parts:
        lower = part.lower()
        if lower in ABBREVIATION_MAP:
            result.append(ABBREVIATION_MAP[lower])
        else:
            result.append(part.capitalize())
    return " ".join(result)


def _sample_values(dicts: List[Dict], field: str, max_samples: int = 50) -> List[Any]:
    """Sample non-null values for a field across a list of dicts."""
    values = []
    for d in dicts:
        if field in d and d[field] is not None:
            values.append(d[field])
            if len(values) >= max_samples:
                break
    return values


def _collect_keys(dicts: List[Dict]) -> set:
    """Collect the union of all keys across a list of dicts."""
    keys = set()
    for d in dicts:
        keys.update(d.keys())
    return keys


def _apply_custom_fields(metadata, custom_fields, level, existing_metadata=None):
    """
    Apply custom field declarations to a metadata dict for a given level.

    Handles output_name renaming: if a custom field specifies output_name,
    the field is registered under that name in field_metadata (and the
    input name is noted for data renaming during processing).

    Args:
        metadata: The field metadata dict to update (modified in place).
        custom_fields: List of custom field declarations.
        level: The level to filter on ("clone", "node", etc.).
        existing_metadata: Optional existing metadata for range preservation.
    """
    if not custom_fields:
        return
    for cf in custom_fields:
        if cf.get("level") != level:
            continue
        output_key = cf.get("output_name", cf["name"])
        entry = {"type": cf["type"], "label": cf["label"]}
        # Preserve range from auto-detection or existing
        if "range" in cf:
            entry["range"] = cf["range"]
        elif output_key in metadata and "range" in metadata[output_key]:
            entry["range"] = metadata[output_key]["range"]
        elif cf["name"] in metadata and "range" in metadata[cf["name"]]:
            entry["range"] = metadata[cf["name"]]["range"]
        # If renaming, remove the original name from metadata
        if output_key != cf["name"] and cf["name"] in metadata:
            del metadata[cf["name"]]
        metadata[output_key] = entry


def _get_nested_value(d: Dict, path: str) -> Any:
    """
    Resolve a simple dot-path against a dict.

    Supports paths like "sample.locus" for nested field access.
    Does NOT support array indexing (e.g., "nodes[].field").
    """
    parts = path.split(".")
    current = d
    for part in parts:
        if isinstance(current, dict) and part in current:
            current = current[part]
        else:
            return None
    return current


# =============================================================================
# Per-level metadata generators
# =============================================================================


def generate_clone_metadata(
    clones: List[Dict],
    custom_fields: Optional[List[Dict]] = None,
) -> Dict[str, Dict[str, str]]:
    """
    Generate field metadata for clone-level fields.

    Introspects actual clone dicts, matches against known registry,
    infers type for unknown fields, and merges custom declarations.

    Args:
        clones: List of clone dictionaries.
        custom_fields: Optional custom field declarations (level=="clone" only).

    Returns:
        Dict mapping field_name -> {"type": ..., "label": ...}
    """
    if not clones:
        return {}

    metadata = {}
    all_keys = _collect_keys(clones)

    # Check for locus in nested sample object
    has_locus = any(
        _get_nested_value(c, "sample.locus") is not None for c in clones[:10]
    )
    if has_locus:
        all_keys.add("locus")

    # Filter out excluded fields
    candidate_keys = all_keys - EXCLUDED_CLONE_FIELDS

    for key in sorted(candidate_keys):
        if key in KNOWN_CLONE_FIELDS:
            # Verify the field actually has data
            if key == "locus":
                values = [
                    _get_nested_value(c, "sample.locus")
                    for c in clones[:50]
                    if _get_nested_value(c, "sample.locus") is not None
                ]
            else:
                values = _sample_values(clones, key)
            if values:
                metadata[key] = dict(KNOWN_CLONE_FIELDS[key])
        else:
            values = _sample_values(clones, key)
            if values:
                field_type = infer_field_type(values)
                metadata[key] = {
                    "type": field_type,
                    "label": humanize_label(key),
                }

    _apply_custom_fields(metadata, custom_fields, "clone")

    return metadata


def generate_node_metadata(
    trees: List[Dict],
    custom_fields: Optional[List[Dict]] = None,
) -> Dict[str, Dict[str, str]]:
    """
    Generate field metadata for node-level fields.

    Introspects nodes across all trees.

    Args:
        trees: List of tree dictionaries with "nodes" key.
        custom_fields: Optional custom field declarations (level=="node" only).

    Returns:
        Dict mapping field_name -> {"type": ..., "label": ...}
    """
    all_nodes = _collect_nodes(trees)
    if not all_nodes:
        return {}

    metadata = {}
    all_keys = _collect_keys(all_nodes) - EXCLUDED_NODE_FIELDS

    for key in sorted(all_keys):
        if key in KNOWN_NODE_FIELDS:
            values = _sample_values(all_nodes, key)
            if values:
                metadata[key] = dict(KNOWN_NODE_FIELDS[key])
        elif key in KNOWN_BRANCH_FIELDS:
            # Branch fields on nodes are handled at branch level
            continue
        else:
            values = _sample_values(all_nodes, key)
            if values:
                field_type = infer_field_type(values)
                metadata[key] = {
                    "type": field_type,
                    "label": humanize_label(key),
                }

    _apply_custom_fields(metadata, custom_fields, "node")

    return metadata


def generate_branch_metadata(
    trees: List[Dict],
    custom_fields: Optional[List[Dict]] = None,
) -> Dict[str, Dict[str, str]]:
    """
    Generate field metadata for branch-level fields.

    Branch fields are stored on node dicts (length, branch_length).

    Args:
        trees: List of tree dictionaries with "nodes" key.
        custom_fields: Optional custom field declarations (level=="branch" only).

    Returns:
        Dict mapping field_name -> {"type": ..., "label": ...}
    """
    all_nodes = _collect_nodes(trees)
    if not all_nodes:
        return {}

    metadata = {}
    all_keys = _collect_keys(all_nodes) - EXCLUDED_BRANCH_FIELDS

    for key in sorted(all_keys):
        if key in KNOWN_BRANCH_FIELDS:
            values = _sample_values(all_nodes, key)
            if values:
                metadata[key] = dict(KNOWN_BRANCH_FIELDS[key])

    _apply_custom_fields(metadata, custom_fields, "branch")

    return metadata


def generate_mutation_metadata(
    trees: List[Dict],
    custom_fields: Optional[List[Dict]] = None,
) -> Dict[str, Dict[str, str]]:
    """
    Generate field metadata for mutation-level fields.

    Mutation fields come from surprise_mutations arrays on nodes.
    Continuous fields include a "range" key with [min, max] for color scale domains.

    Args:
        trees: List of tree dictionaries with "nodes" key.
        custom_fields: Optional custom field declarations (level=="mutation" only).

    Returns:
        Dict mapping field_name -> {"type": ..., "label": ..., "range"?: [...]}
    """
    all_mutations = _collect_mutations(trees)

    # Check if nodes have AA sequence data — if so, the web app will derive
    # per-mutation child_aa/parent_aa fields during alignment rendering,
    # even if no surprise_mutations arrays exist in the data.
    all_nodes = _collect_nodes(trees, max_nodes=20)
    has_aa_sequences = any(
        n.get("sequence_alignment_aa") for n in all_nodes if isinstance(n, dict)
    )

    if not all_mutations:
        # No pre-computed mutation data; declare derived fields if sequences exist
        metadata = {}
        if has_aa_sequences:
            metadata["child_aa"] = {"type": "aa", "label": "Child Amino Acid"}
            metadata["parent_aa"] = {"type": "tooltip", "label": "Parent Amino Acid"}
        _apply_custom_fields(metadata, custom_fields, "mutation")
        return metadata

    # Collect ALL mutations (not sampled) for accurate range computation
    all_mutations_full = _collect_mutations(trees, max_mutations=100000)

    metadata = {}
    all_keys = _collect_keys(all_mutations) - EXCLUDED_MUTATION_FIELDS

    for key in sorted(all_keys):
        if key in KNOWN_MUTATION_FIELDS:
            values = _sample_values(all_mutations, key)
            if values:
                entry = dict(KNOWN_MUTATION_FIELDS[key])
                if entry["type"] == "continuous":
                    field_range = compute_range(all_mutations_full, key)
                    if field_range:
                        entry["range"] = field_range
                metadata[key] = entry
        else:
            values = _sample_values(all_mutations, key)
            if values:
                field_type = infer_field_type(values)
                entry = {
                    "type": field_type,
                    "label": humanize_label(key),
                }
                if field_type == "continuous":
                    field_range = compute_range(all_mutations_full, key)
                    if field_range:
                        entry["range"] = field_range
                metadata[key] = entry

    _apply_custom_fields(metadata, custom_fields, "mutation")

    return metadata


# =============================================================================
# Helper functions for collecting data across trees
# =============================================================================


def _collect_nodes(trees: List[Dict], max_nodes: int = 200) -> List[Dict]:
    """Collect node dicts from across trees (up to max_nodes for sampling)."""
    nodes = []
    for tree in trees:
        tree_nodes = tree.get("nodes", [])
        if isinstance(tree_nodes, dict):
            tree_nodes = list(tree_nodes.values())
        for node in tree_nodes:
            if isinstance(node, dict):
                nodes.append(node)
                if len(nodes) >= max_nodes:
                    return nodes
    return nodes


def _collect_mutations(trees: List[Dict], max_mutations: int = 200) -> List[Dict]:
    """Collect mutation dicts from surprise_mutations arrays across tree nodes."""
    mutations = []
    for tree in trees:
        tree_nodes = tree.get("nodes", [])
        if isinstance(tree_nodes, dict):
            tree_nodes = list(tree_nodes.values())
        for node in tree_nodes:
            if isinstance(node, dict):
                surprise = node.get("surprise_mutations")
                if isinstance(surprise, list):
                    for mut in surprise:
                        if isinstance(mut, dict):
                            mutations.append(mut)
                            if len(mutations) >= max_mutations:
                                return mutations
    return mutations


# =============================================================================
# Top-level generation function
# =============================================================================


def generate_field_metadata(
    clones: List[Dict],
    trees: Optional[List[Dict]] = None,
    custom_fields: Optional[List[Dict]] = None,
) -> Dict[str, Dict[str, Dict[str, str]]]:
    """
    Generate complete field_metadata for a dataset.

    Introspects clones and trees to discover available fields at each level,
    classifies them by type, and generates human-readable labels. Custom field
    declarations from YAML config override auto-introspection.

    Args:
        clones: List of clone dictionaries for the dataset.
        trees: List of tree dictionaries for the dataset (optional).
        custom_fields: List of custom field declaration dicts, each with
            keys: name, level, type, label, and optionally path.

    Returns:
        Dict with level keys ("clone", "node", "branch", "mutation"),
        each containing a dict of field_name -> {"type": ..., "label": ...}.
        Levels with no fields are omitted.
    """
    if trees is None:
        trees = []

    result = {}

    clone_meta = generate_clone_metadata(clones, custom_fields)
    if clone_meta:
        result["clone"] = clone_meta

    node_meta = generate_node_metadata(trees, custom_fields)
    if node_meta:
        result["node"] = node_meta

    branch_meta = generate_branch_metadata(trees, custom_fields)
    if branch_meta:
        result["branch"] = branch_meta

    mutation_meta = generate_mutation_metadata(trees, custom_fields)
    if mutation_meta:
        result["mutation"] = mutation_meta

    return result
