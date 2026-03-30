"""
Shared constants, enumerated types, field registries, and reference tables.

This module centralizes all configuration constants used across olmsted-cli.
It has NO dependencies on other project modules.
"""

# =============================================================================
# Enumerated Types
# =============================================================================

#: Valid field types for field_metadata entries.
#: - continuous: numeric values (axes, size, color scales)
#: - categorical: string/enum values (color, shape, facet)
#: - tooltip: display-only (shown in tooltips, not for encoding)
#: - aa: amino acid identity (uses full genetic alphabet)
#: - dna: nucleotide identity (uses full genetic alphabet)
#: - skip: exclude from output metadata (keeps entry in config for docs)
FIELD_TYPES = {"continuous", "categorical", "tooltip", "aa", "dna", "skip"}

#: Valid data levels for field_metadata.
#: - clone: clonal family level (scatterplot axes, color, facet)
#: - node: tree node level (node properties, tooltips)
#: - branch: tree branch level (branch coloring, width)
#: - mutation: per-mutation level (alignment coloring)
FIELD_LEVELS = {"clone", "node", "branch", "mutation"}


# =============================================================================
# Known Field Registries
# =============================================================================
#
# Maps field names to their default type and label. Used during field_metadata
# generation: if a field name matches a registry entry, its type/label are
# used instead of being inferred from values.

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

#: Mapping from level name to its known fields registry.
KNOWN_FIELDS_BY_LEVEL = {
    "clone": KNOWN_CLONE_FIELDS,
    "node": KNOWN_NODE_FIELDS,
    "branch": KNOWN_BRANCH_FIELDS,
    "mutation": KNOWN_MUTATION_FIELDS,
}


# =============================================================================
# Exclusion Sets
# =============================================================================
#
# Fields to never include in field_metadata — structural, positional, or
# internal fields not useful for visualization controls.

# Excluded fields are never shown anywhere — not in build-config output,
# not in field_metadata. These are fields whose values are structurally
# unpresentable (nested objects, long sequences, complex arrays).

EXCLUDED_CLONE_FIELDS = {
    # Nested objects/arrays (not scalar values)
    "dataset", "sample", "trees",
    # Sequences (long strings, never useful in dropdowns)
    "germline_alignment", "germline_sequence", "germline_alignment_light",
    "naive_sequence", "cdr3_sequence",
    # Gene support probability arrays (complex nested objects)
    "v_per_gene_support", "d_per_gene_support", "j_per_gene_support",
    # ID arrays
    "unique_ids", "timepoint_ids",
}

EXCLUDED_NODE_FIELDS = {
    # Sequences (long strings)
    "sequence_alignment", "sequence_alignment_aa",
    "sequence_alignment_light", "sequence_alignment_light_aa",
    "aa_sequence", "junction", "junction_aa",
    # The surprise_mutations array itself (sub-fields are mutation-level)
    "surprise_mutations",
    # Multiplicity arrays (complex objects, not scalar)
    "timepoint_multiplicities", "cluster_timepoint_multiplicities",
    # Structural (required for tree topology, not for viz encoding)
    "sequence_id", "node_id", "parent", "is_root",
    # Gene calls on nodes (redundant with clone-level)
    "v_call", "d_call", "j_call",
    # Refs (redundant with clone-level)
    "sample_id", "timepoint",
}

EXCLUDED_BRANCH_FIELDS = {
    "sequence_id", "node_id", "parent",
}

EXCLUDED_MUTATION_FIELDS = {
    "site",
}

#: Mapping from level name to its exclusion set.
EXCLUDED_FIELDS_BY_LEVEL = {
    "clone": EXCLUDED_CLONE_FIELDS,
    "node": EXCLUDED_NODE_FIELDS,
    "branch": EXCLUDED_BRANCH_FIELDS,
    "mutation": EXCLUDED_MUTATION_FIELDS,
}


# =============================================================================
# PCP CSV Column Sets
# =============================================================================

#: Known PCP CSV columns handled by the parser. Extra columns are captured
#: as custom node-level fields.
KNOWN_PCP_COLUMNS = {
    "sample_id", "family", "parent_name", "child_name",
    "parent_heavy", "child_heavy", "parent_light", "child_light",
    "branch_length", "edge_length", "depth", "distance", "sample_count",
    "v_gene_heavy", "d_gene_heavy", "j_gene_heavy",
    "v_gene_light", "d_gene_light", "j_gene_light",
    "v_gene_start_heavy", "v_gene_end_heavy",
    "d_gene_start_heavy", "d_gene_end_heavy",
    "j_gene_start_heavy", "j_gene_end_heavy",
    "v_gene_start_light", "v_gene_end_light",
    "d_gene_start_light", "d_gene_end_light",
    "j_gene_start_light", "j_gene_end_light",
    "cdr1_codon_start_heavy", "cdr1_codon_end_heavy",
    "cdr2_codon_start_heavy", "cdr2_codon_end_heavy",
    "cdr3_codon_start_heavy", "cdr3_codon_end_heavy",
    "cdr1_codon_start_light", "cdr1_codon_end_light",
    "cdr2_codon_start_light", "cdr2_codon_end_light",
    "cdr3_codon_start_light", "cdr3_codon_end_light",
    "parent_is_naive", "child_is_leaf",
    "light_chain_type",
}

#: Known tree CSV columns handled by the parser. Extra columns are captured
#: as clone-level fields.
KNOWN_TREE_COLUMNS = {
    "family_name", "family", "sample_id",
    "newick_tree", "newick",
    "rate_scale_heavy", "rate_scale_light",
}


# =============================================================================
# Reference / Alias Tables
# =============================================================================

#: Chain-specific column name aliases for PCP CSV parsing.
#: Maps common alternative names to canonical column names.
CHAIN_COLUMN_ALIASES = {
    "parent_seq": "parent_heavy",
    "child_seq": "child_heavy",
    "parent_sequence": "parent_heavy",
    "child_sequence": "child_heavy",
    "v_gene": "v_gene_heavy",
    "d_gene": "d_gene_heavy",
    "j_gene": "j_gene_heavy",
    "v_call": "v_gene_heavy",
    "d_call": "d_gene_heavy",
    "j_call": "j_gene_heavy",
    "cdr1_start": "cdr1_codon_start_heavy",
    "cdr1_end": "cdr1_codon_end_heavy",
    "cdr2_start": "cdr2_codon_start_heavy",
    "cdr2_end": "cdr2_codon_end_heavy",
    "cdr3_start": "cdr3_codon_start_heavy",
    "cdr3_end": "cdr3_codon_end_heavy",
}

#: Cross-format field aliases for output_name suggestions in build-config.
#: Maps input field names to suggested canonical output names.
#: NOT auto-applied during processing — user-facing suggestions only.
FIELD_ALIASES = {
    "v_gene": "v_call",
    "v_gene_heavy": "v_call",
    "d_gene": "d_call",
    "d_gene_heavy": "d_call",
    "j_gene": "j_call",
    "j_gene_heavy": "j_call",
    "v_gene_light": "v_call_light",
    "j_gene_light": "j_call_light",
    "rearrangement_count": "unique_seqs_count",
    "sampled_seqs_count": "unique_seqs_count",
    "size": "total_read_count",
    "branch_length": "length",
    "mut_to": "child_aa",
    "mut_from": "parent_aa",
}

#: Suggested type overrides for fields that auto-detection may not classify
#: correctly. Used by build-config to pre-fill the type field. These are
#: suggestions only — the user can change them in their config.
#:
#: Fields mapped to "skip" are non-visualization metadata that would
#: pollute the web app dropdowns if included.
SUGGESTED_FIELD_TYPES = {
    # --- Identifiers (useful for debugging, not for viz encoding) ---
    "ident": "skip",
    "clone_id": "skip",
    "dataset_id": "skip",
    "repertoire_id": "skip",
    "pair_id": "skip",
    "seed_id": "skip",
    "schema_version": "skip",
    "type": "skip",
    "build": "skip",
    "trees_meta": "skip",
    # --- Non-visualization metadata ---
    "partition": "skip",
    "path": "skip",
    "sorted_index": "skip",
    # --- Alignment positions (may be useful as tooltips) ---
    "v_alignment_start": "skip",
    "v_alignment_end": "skip",
    "v_sequence_start": "skip",
    "v_sequence_end": "skip",
    "v_germline_start": "skip",
    "v_germline_end": "skip",
    "d_alignment_start": "skip",
    "d_alignment_end": "skip",
    "d_sequence_start": "skip",
    "d_sequence_end": "skip",
    "d_germline_start": "skip",
    "d_germline_end": "skip",
    "j_alignment_start": "skip",
    "j_alignment_end": "skip",
    "j_sequence_start": "skip",
    "j_sequence_end": "skip",
    "j_germline_start": "skip",
    "j_germline_end": "skip",
    "cdr1_alignment_start": "skip",
    "cdr1_alignment_end": "skip",
    "cdr2_alignment_start": "skip",
    "cdr2_alignment_end": "skip",
    "cdr1_alignment_start_light": "skip",
    "cdr1_alignment_end_light": "skip",
    "cdr2_alignment_start_light": "skip",
    "cdr2_alignment_end_light": "skip",
    "junction_start": "skip",
    "junction_end": "skip",
    "junction_start_light": "skip",
    "junction_length_light": "skip",
}


#: Abbreviation map for humanize_label(): maps lowercase tokens to
#: their preferred display form.
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
    "nt": "Nucleotide",
}

#: Amino acid single-character alphabet (for type inference).
AA_CHARS = set("ACDEFGHIKLMNPQRSTVWY*-X")

#: DNA/RNA single-character alphabet (for type inference).
DNA_CHARS = set("ACGTURYSWKMBDHVN-.")
