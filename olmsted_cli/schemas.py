#!/usr/bin/env python
"""
Unified schema definitions for Olmsted data structures.

This module contains JSON Schema definitions for datasets, clones, trees, and nodes
used throughout the Olmsted CLI for both AIRR and PCP data processing.

The schemas are designed to be flexible enough to accommodate both formats while
maintaining consistency in the output structure.

NOTE: The AIRR schema components reference the official AIRR schema from
airr-standards/specs/airr-schema.yaml. The SCHEMA_VERSION constant corresponds
to the 'version' field in the Info section of that schema.
"""

from .constants import DISPLAY_MODES, FIELD_LEVELS, FIELD_TYPES

# Version Constants
# SCHEMA_VERSION corresponds to Info.version in airr-standards/specs/airr-schema.yaml
SCHEMA_VERSION = "2.0.0"

# Output display modes (skip means "not in output", so exclude it from schema)
_OUTPUT_DISPLAY_MODES = sorted(DISPLAY_MODES - {"skip"})

# Schema fragment for a single field_metadata entry
_FIELD_ENTRY_SCHEMA = {
    "type": "object",
    "properties": {
        "type": {
            "type": "string",
            "enum": sorted(FIELD_TYPES),
        },
        "display": {
            "type": "string",
            "enum": _OUTPUT_DISPLAY_MODES,
        },
        "label": {"type": "string"},
        "range": {
            "type": "array",
            "items": {"type": "number"},
            "minItems": 2,
            "maxItems": 2,
        },
    },
    "required": ["type", "label"],
}

# Timepoint multiplicity schema - for individual timepoint/multiplicity pairs
timepoint_multiplicity_spec = {
    "type": "object",
    "properties": {
        "multiplicity": {
            "description": "Number of sequences at this timepoint",
            "type": "integer",
        },
        "timepoint_id": {
            "description": "Timepoint identifier",
            "type": ["string", "null"],
        },
    },
    "required": ["multiplicity", "timepoint_id"],
    "additionalProperties": False,
}

# Node schema - unified for both AIRR and PCP
node_spec = {
    "title": "Node",
    "description": "Node in a phylogenetic tree",
    "type": "object",
    "required": ["sequence_id", "sequence_alignment", "sequence_alignment_aa"],
    "properties": {
        "sequence_id": {
            "description": "Identifier for this node that matches the id in the newick string",
            "type": "string",
        },
        "node_id": {
            "description": "Alternative identifier for this node",
            "type": "string",
        },
        "parent": {
            "description": "Parent node ID",
            "type": ["string", "null"],
        },
        "branch_length": {
            "description": "Branch length to parent",
            "type": ["number", "null"],
        },
        "length": {
            "description": "Branch length",
            "type": ["number", "null"],
        },
        "distance": {
            "description": "Distance from root",
            "type": ["number", "null"],
        },
        "sequence_alignment": {
            "description": "Nucleotide sequence alignment, including any indel corrections or spacers",
            "type": "string",
        },
        "sequence_alignment_aa": {
            "description": "Amino acid sequence alignment, including any indel corrections or spacers",
            "type": ["string", "null"],
        },
        "sample_id": {
            "description": "Sample identifier",
            "type": ["string", "null"],
        },
        "type": {
            "description": "Type of node (leaf, internal, node, or root)",
            "enum": ["leaf", "node", "internal", "root", None],
            "type": ["string", "null"],
        },
        "is_root": {
            "description": "Whether this node is the root",
            "type": ["boolean", "null"],
        },
        "subtree_size": {
            "description": "Number of descendants",
            "type": ["integer", "null"],
        },
        "timepoint": {
            "description": "Timepoint identifier",
            "type": ["string", "null"],
        },
        "timepoint_id": {
            "description": "Timepoint identifier (alternative field name)",
            "type": ["string", "null"],
        },
        "lbi": {
            "description": "Local branching index",
            "type": ["number", "null"],
        },
        "lbr": {
            "description": "Local branching ratio",
            "type": ["number", "null"],
        },
        "affinity": {
            "description": "Binding affinity",
            "type": ["number", "null"],
        },
        "relative_affinity": {
            "description": "Relative binding affinity",
            "type": ["number", "null"],
        },
        "affinity_class": {
            "description": "Affinity classification",
            "type": ["string", "null"],
        },
        "aa_sequence": {
            "description": "Full amino acid sequence",
            "type": ["string", "null"],
        },
        "junction": {
            "description": "CDR3 junction nucleotide sequence",
            "type": ["string", "null"],
        },
        "junction_aa": {
            "description": "CDR3 junction amino acid sequence",
            "type": ["string", "null"],
        },
        "v_call": {
            "description": "V gene assignment",
            "type": ["string", "null"],
        },
        "d_call": {
            "description": "D gene assignment",
            "type": ["string", "null"],
        },
        "j_call": {
            "description": "J gene assignment",
            "type": ["string", "null"],
        },
        "count": {
            "description": "Sequence count",
            "type": ["integer", "null"],
        },
        "confidence": {
            "description": "Confidence score or bootstrap value",
            "type": ["number", "null"],
        },
        "multiplicity": {
            "description": "Number of times sequence was observed",
            "type": ["integer", "null"],
        },
        "cluster_multiplicity": {
            "description": "Cumulative count if sequences were clustered",
            "type": ["integer", "null"],
        },
        "timepoint_multiplicities": {
            "description": "Multiplicities per timepoint",
            "type": ["array", "null"],
            "items": timepoint_multiplicity_spec,
        },
        "cluster_timepoint_multiplicities": {
            "description": "Cluster multiplicities per timepoint",
            "type": ["array", "null"],
            "items": timepoint_multiplicity_spec,
        },
        "scaled_affinity": {
            "description": "Scaled binding affinity (min-max normalized)",
            "type": ["number", "null"],
        },
    },
    "additionalProperties": True,
}

# Tree schema - unified for both AIRR and PCP
tree_spec = {
    "title": "Tree",
    "description": "Phylogenetic tree and possibly ancestral state reconstruction of sequences in a clonal family",
    "type": "object",
    "required": ["newick"],
    "properties": {
        "ident": {
            "description": "Tree identifier",
            "type": ["string", "null"],
        },
        "tree_id": {
            "description": "Unique identifier for the tree",
            "type": ["string", "null"],
        },
        "clone_id": {
            "description": "Identifier for the associated clone",
            "type": ["string", "null"],
        },
        "timepoint_ids": {
            "description": "Time points included in this tree",
            "type": ["array", "null"],
            "items": {"type": "string"},
        },
        "downsampling_strategy": {
            "description": "Method used to downsample sequences before tree inference",
            "type": ["string", "null"],
        },
        "diversity": {
            "description": "Mean distance from all tree nodes to their most recent common ancestor",
            "type": ["number", "null"],
        },
        "min_junction_length": {
            "description": "Minimum CDR3 junction nucleotide length",
            "type": ["number", "null"],
        },
        "max_junction_length": {
            "description": "Maximum CDR3 junction nucleotide length",
            "type": ["number", "null"],
        },
        "sample_id": {
            "description": "Sample identifier for associated sequences",
            "type": ["string", "null"],
        },
        "reconstruction_method": {
            "description": (
                "Method used to build this tree (e.g. 'dnapars', 'raxml_ng'). "
                "Only populated when supplied by the input; absent means the "
                "reconstruction method is unknown."
            ),
            "type": ["string", "null"],
        },
        "newick": {
            "description": "Tree in Newick format",
            "type": "string",
        },
        "nodes": {
            "description": "Tree nodes as array of node objects or dict keyed by node ID",
            "oneOf": [
                {
                    "type": "array",
                    "items": node_spec,
                },
                {
                    "type": "object",
                    "additionalProperties": node_spec,
                },
            ],
        },
    },
    "additionalProperties": True,
}

# Clone schema - unified for both AIRR and PCP
clone_spec = {
    "title": "Clone",
    "description": "Clonal family of sequences deriving from a particular rearrangement event",
    "type": "object",
    "required": [
        "unique_seqs_count",
        "mean_mut_freq",
    ],
    "properties": {
        "ident": {
            "description": "Clone identifier",
            "type": ["string", "null"],
        },
        "clone_id": {
            "description": "Unique identifier of the clone",
            "type": ["string", "null"],
        },
        "dataset_id": {
            "description": "Dataset identifier",
            "type": ["string", "null"],
        },
        "sample_id": {
            "description": "Sample identifier",
            "type": ["string", "null"],
        },
        "subject_id": {
            "description": "Subject/participant identifier",
            "type": ["string", "null"],
        },
        "unique_seqs_count": {
            "description": "Number of unique sequences in the clone",
            "type": "integer",
        },
        "total_read_count": {
            "description": "Total number of reads in the clone",
            "type": ["integer", "null"],
        },
        "clone_count": {
            "description": "Total number of sequences within the clone",
            "type": ["integer", "null"],
        },
        "mean_mut_freq": {
            "description": "Mean mutation frequency",
            "type": "number",
        },
        "naive_sequence": {
            "description": "Unmutated common ancestor sequence",
            "type": ["string", "null"],
        },
        "cdr3_sequence": {
            "description": "CDR3 nucleotide sequence",
            "type": ["string", "null"],
        },
        "v_call": {
            "description": "V gene assignment with allele",
            "type": ["string", "null"],
        },
        "d_call": {
            "description": "D gene assignment with allele",
            "type": ["string", "null"],
        },
        "j_call": {
            "description": "J gene assignment with allele",
            "type": ["string", "null"],
        },
        "v_sequence_start": {
            "description": "Start position of V gene alignment in query sequence",
            "type": ["integer", "null"],
        },
        "v_sequence_end": {
            "description": "End position of V gene alignment in query sequence",
            "type": ["integer", "null"],
        },
        "v_germline_start": {
            "description": "Start position of V gene alignment in germline sequence",
            "type": ["integer", "null"],
        },
        "v_germline_end": {
            "description": "End position of V gene alignment in germline sequence",
            "type": ["integer", "null"],
        },
        "v_alignment_start": {
            "description": "Start position in the V gene alignment",
            "type": "integer",
        },
        "v_alignment_end": {
            "description": "End position in the V gene alignment",
            "type": "integer",
        },
        "d_sequence_start": {
            "description": "Start position of D gene alignment in query sequence",
            "type": ["integer", "null"],
        },
        "d_sequence_end": {
            "description": "End position of D gene alignment in query sequence",
            "type": ["integer", "null"],
        },
        "d_germline_start": {
            "description": "Start position of D gene alignment in germline sequence",
            "type": ["integer", "null"],
        },
        "d_germline_end": {
            "description": "End position of D gene alignment in germline sequence",
            "type": ["integer", "null"],
        },
        "d_alignment_start": {
            "description": "Start position in the D gene alignment",
            "type": ["integer", "null"],
        },
        "d_alignment_end": {
            "description": "End position in the D gene alignment",
            "type": ["integer", "null"],
        },
        "j_sequence_start": {
            "description": "Start position of J gene alignment in query sequence",
            "type": ["integer", "null"],
        },
        "j_sequence_end": {
            "description": "End position of J gene alignment in query sequence",
            "type": ["integer", "null"],
        },
        "j_germline_start": {
            "description": "Start position of J gene alignment in germline sequence",
            "type": ["integer", "null"],
        },
        "j_germline_end": {
            "description": "End position of J gene alignment in germline sequence",
            "type": ["integer", "null"],
        },
        "j_alignment_start": {
            "description": "Start position in the J gene alignment",
            "type": "integer",
        },
        "j_alignment_end": {
            "description": "End position in the J gene alignment",
            "type": "integer",
        },
        "junction_start": {
            "description": "Start position of CDR3 junction",
            "type": ["integer", "null"],
        },
        "junction_end": {
            "description": "End position of CDR3 junction",
            "type": ["integer", "null"],
        },
        "junction_length": {
            "description": "Length of CDR3 junction",
            "type": ["integer", "null"],
        },
        "germline_alignment": {
            "description": "Assembled germline sequence aligned to query",
            "type": ["string", "null"],
        },
        "germline_sequence": {
            "description": "Full germline sequence",
            "type": ["string", "null"],
        },
        "has_seed": {
            "description": "Whether clone has seed sequence",
            "type": ["boolean", "null"],
        },
        "seed_id": {
            "description": "Seed sequence identifier",
            "type": ["string", "null"],
        },
        "trees": {
            "description": "Associated phylogenetic trees",
            "type": ["array", "null"],
            "items": tree_spec,
        },
        "unique_ids": {
            "description": "List of unique sequence identifiers",
            "type": ["array", "null"],
            "items": {"type": "string"},
        },
        "timepoint_ids": {
            "description": "Time points included",
            "type": ["array", "null"],
            "items": {"type": "string"},
        },
        "v_per_gene_support": {
            "description": "V gene assignment probabilities",
            "type": ["array", "null"],
            "items": {
                "type": "object",
                "properties": {
                    "gene": {"type": "string"},
                    "ident": {"type": "string"},
                    "prob": {"type": "number"},
                },
                "required": ["gene", "ident", "prob"],
                "additionalProperties": False,
            },
        },
        "j_per_gene_support": {
            "description": "J gene assignment probabilities",
            "type": ["array", "null"],
            "items": {
                "type": "object",
                "properties": {
                    "gene": {"type": "string"},
                    "ident": {"type": "string"},
                    "prob": {"type": "number"},
                },
                "required": ["gene", "ident", "prob"],
                "additionalProperties": False,
            },
        },
        "d_per_gene_support": {
            "description": "D gene assignment probabilities",
            "type": ["array", "null"],
            "items": {
                "type": "object",
                "properties": {
                    "gene": {"type": "string"},
                    "ident": {"type": "string"},
                    "prob": {"type": "number"},
                },
                "required": ["gene", "ident", "prob"],
                "additionalProperties": False,
            },
        },
        # Paired data fields
        "is_paired": {
            "description": "Whether this clone is part of a paired heavy/light chain dataset",
            "type": ["boolean", "null"],
        },
        "pair_id": {
            "description": "Identifier linking paired heavy and light chain clones",
            "type": ["string", "null"],
        },
    },
    "additionalProperties": True,
}

# Dataset schema - unified for both AIRR and PCP
dataset_spec = {
    "$schema": "https://json-schema.org/draft-07/schema#",
    "$id": "https://olmstedviz.org/input.schema.json",
    "title": "Olmsted Dataset",
    "description": "Olmsted dataset input file",
    "type": "object",
    "required": ["dataset_id"],
    "properties": {
        "schema_version": {
            "description": "Schema version",
            "type": ["string", "null"],
        },
        "ident": {
            "description": "Dataset identifier",
            "type": ["string", "null"],
        },
        "dataset_id": {
            "description": "Unique identifier for the dataset",
            "type": "string",
        },
        "type": {
            "description": (
                "Free-form dataset-type label. Passed through from input when "
                "present; never synthesized by olmsted-cli. Source format is "
                "recorded separately in metadata.source_format."
            ),
            "type": ["string", "null"],
        },
        "build": {
            "description": "Build information",
            "type": ["object", "null"],
            "properties": {
                "commit": {"type": "string"},
                "time": {"type": "string"},
            },
        },
        "study_id": {
            "description": "Study identifier",
            "type": ["string", "null"],
        },
        "study_title": {
            "description": "Study title",
            "type": ["string", "null"],
        },
        "study_description": {
            "description": "Study description",
            "type": ["string", "null"],
        },
        "study_contact": {
            "description": "Study contact information",
            "type": ["string", "null"],
        },
        "inclusion_exclusion_criteria": {
            "description": "Inclusion/exclusion criteria",
            "type": ["string", "null"],
        },
        "contributors": {
            "description": "List of contributors",
            "type": ["array", "null"],
            "items": {
                "type": "object",
                "properties": {
                    "contributor_id": {"type": "string"},
                    "name": {"type": "string"},
                    "orcid_id": {"type": "string"},
                    "affiliation": {"type": "string"},
                    "affiliation_ror_id": {"type": "string"},
                },
            },
        },
        "pub_ids": {
            "description": "Publication identifiers",
            "type": ["array", "null"],
            "items": {"type": "string"},
        },
        "grants": {
            "description": "Grant information",
            "type": ["array", "null"],
            "items": {"type": "string"},
        },
        "keywords_study": {
            "description": "Study keywords",
            "type": ["array", "null"],
            "items": {"type": "string"},
        },
        "subjects": {
            "description": "Subject information",
            "type": ["array", "null"],
            "items": {
                "type": "object",
                "properties": {
                    "ident": {"type": "string"},
                    "subject_id": {"type": "string"},
                    "subject_species": {"type": "string"},
                    "sex": {"type": "string"},
                    "age_min": {"type": "number"},
                    "age_max": {"type": "number"},
                    "age_unit": {"type": "string"},
                    "age_event": {"type": "string"},
                    "ancestry_population": {"type": "string"},
                    "ethnicity": {"type": "string"},
                    "race": {"type": "string"},
                    "strain_name": {"type": "string"},
                },
            },
        },
        "subject_id": {
            "description": "Subject/participant identifier (single subject datasets)",
            "type": ["string", "null"],
        },
        "subject_species": {
            "description": "Subject species",
            "type": ["string", "null"],
        },
        "sex": {
            "description": "Subject sex",
            "type": ["string", "null"],
        },
        "age_min": {
            "description": "Minimum age",
            "type": ["number", "null"],
        },
        "age_max": {
            "description": "Maximum age",
            "type": ["number", "null"],
        },
        "age_unit": {
            "description": "Age unit",
            "type": ["string", "null"],
        },
        "age_event": {
            "description": "Age event",
            "type": ["string", "null"],
        },
        "ancestry_population": {
            "description": "Ancestry population",
            "type": ["string", "null"],
        },
        "ethnicity": {
            "description": "Ethnicity",
            "type": ["string", "null"],
        },
        "race": {
            "description": "Race",
            "type": ["string", "null"],
        },
        "strain_name": {
            "description": "Strain name",
            "type": ["string", "null"],
        },
        "samples": {
            "description": "Sample information",
            "type": ["array", "null"],
            "items": {
                "type": "object",
                "properties": {
                    "ident": {"type": "string"},
                    "sample_id": {"type": "string"},
                    "subject_id": {"type": "string"},
                },
            },
        },
        "timepoints": {
            "description": "Time point information",
            "type": ["array", "null"],
            "items": {
                "type": "object",
                "properties": {
                    "timepoint_id": {"type": "string"},
                    "label": {"type": "string"},
                    "sample_id": {"type": "string"},
                },
            },
        },
        "field_metadata": {
            "description": "Metadata describing available data fields at each level",
            "type": ["object", "null"],
            "properties": {
                level: {
                    "description": f"{level.title()}-level field metadata",
                    "type": "object",
                    "additionalProperties": _FIELD_ENTRY_SCHEMA,
                }
                for level in sorted(FIELD_LEVELS)
            },
            "additionalProperties": False,
        },
        "clones": {
            "description": "Clonal families in the dataset",
            "type": ["array", "null"],
            "items": clone_spec,
        },
        "trees": {
            "description": "Phylogenetic trees in the dataset",
            "type": ["array", "null"],
            "items": tree_spec,
        },
    },
    "additionalProperties": True,
}


# Helper functions for AIRR-specific schema generation
def id_spec(description="Identifier"):
    """Create an ID specification with custom description."""
    return {
        "description": description,
        "type": "string",
    }


def sequence_spec(description):
    """Create a sequence specification with custom description."""
    return {
        "description": description,
        "type": "string",
    }


def multiplicity_spec(description=None):
    """Create a multiplicity specification for AIRR fields."""
    return {
        "description": description
        or "Number of times sequence was observed in the sample. The presence of a given sequence in a clonal family may represent many identical such sequences in the original sample.",
        "type": ["integer", "null"],
        "minimum": 0,
    }


# AIRR-specific schemas
ident_spec = {
    "description": "UUID specific to the given object",
    "type": "string",
}

build_spec = {
    "description": "Information about how a dataset was built",
    "type": "object",
    "required": ["commit"],
    "title": "Build info",
    "properties": {
        "commit": {
            "description": "Commit sha of whatever build system you used to process the data",
            "type": "string",
        },
        "time": {
            "description": "Time at which build was initiated",
            "type": "string",
        },
    },
}

timepoint_multiplicity_spec = {
    "title": "Timepoint multiplicity",
    "description": "Multiplicity at a specific time",
    "type": "object",
    "properties": {
        "timepoint_id": {
            "description": "Id associated with the timepoint in question",
            "type": "string",
        },
        "multiplicity": {
            "description": "Number of times sequence was observed at the given timepoint",
            "type": ["integer", "null"],
            "minimum": 0,
        },
    },
}

sample_spec = {
    "title": "Sample",
    "description": "A sample is a collection of sequences",
    "type": "object",
    "required": ["locus"],
    "properties": {
        "ident": ident_spec,
        "sample_id": {
            "description": "Sample id",
            "type": "string",
        },
        "timepoint_id": {
            "description": 'Timepoint associated with this sample (may choose "merged" if data has been combined from multiple timepoints)',
            "type": "string",
        },
        "locus": {
            "description": "B-cell Locus",
            "type": "string",
        },
    },
}

subject_spec = {
    "title": "Subject",
    "description": "Subject from which the clonal family was sampled",
    "type": "object",
    "required": ["subject_id"],
    "properties": {
        "ident": ident_spec,
        "subject_id": {
            "description": "Subject id",
            "type": "string",
        },
    },
}

seed_spec = {
    "title": "Seed",
    "description": "A sequence of interest among other clonal family members",
    "type": ["object", "null"],
    "required": ["seed_id"],
    "properties": {
        "ident": ident_spec,
        "seed_id": {
            "description": "Seed id",
            "type": "string",
        },
    },
}

# Export all schemas, constants, and helper functions
__all__ = [
    # Version constants
    "SCHEMA_VERSION",
    # Main unified schemas
    "node_spec",
    "tree_spec",
    "clone_spec",
    "dataset_spec",
    # AIRR-specific schemas
    "ident_spec",
    "build_spec",
    "timepoint_multiplicity_spec",
    "sample_spec",
    "subject_spec",
    "seed_spec",
    # Helper functions
    "id_spec",
    "sequence_spec",
    "multiplicity_spec",
]
