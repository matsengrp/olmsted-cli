"""
Typed data structures for Olmsted CLI.

This module provides TypedDict definitions for all Olmsted data structures,
enabling IDE autocompletion, static type checking, and clear documentation.

These types correspond to the JSON Schema definitions in schemas.py and
represent the structure of Olmsted output files.

This module has NO dependencies on other project modules, ensuring clean
architecture with types at the bottom of the dependency hierarchy.

Usage:
    from olmsted_cli.types import OlmstedNode, OlmstedClone, OlmstedTree

    # Type hints for function parameters
    def process_node(node: OlmstedNode) -> None:
        print(node["sequence_id"])  # IDE autocomplete works here

    # For the high-level API, use:
    from olmsted_cli.api import OlmstedData
"""

from __future__ import annotations

from typing import Dict, List, Literal, Optional, TypedDict


# =============================================================================
# Node Types
# =============================================================================


class TimepointMultiplicity(TypedDict, total=False):
    """Multiplicity at a specific timepoint."""

    timepoint_id: Optional[str]
    multiplicity: int


class OlmstedNode(TypedDict, total=False):
    """
    Node in a phylogenetic tree.

    Represents either an observed sequence (leaf), inferred ancestor (internal),
    or unmutated common ancestor (root/naive).
    """

    # Required fields
    sequence_id: str
    sequence_alignment: str
    sequence_alignment_aa: Optional[str]

    # Tree structure
    parent: Optional[str]
    type: Optional[Literal["root", "internal", "leaf", "node"]]
    is_root: Optional[bool]

    # Branch metrics
    length: Optional[float]  # Branch length to parent
    distance: Optional[float]  # Cumulative distance from root
    branch_length: Optional[float]  # Alternative field name

    # Multiplicity (observation counts)
    multiplicity: Optional[int]
    cluster_multiplicity: Optional[int]
    timepoint_multiplicities: Optional[List[TimepointMultiplicity]]
    cluster_timepoint_multiplicities: Optional[List[TimepointMultiplicity]]

    # Phylogenetic metrics
    lbi: Optional[float]  # Local Branching Index
    lbr: Optional[float]  # Local Branching Ratio
    affinity: Optional[float]  # Binding affinity (often = LBI)
    scaled_affinity: Optional[float]  # Min-max normalized affinity
    relative_affinity: Optional[float]
    affinity_class: Optional[str]
    subtree_size: Optional[int]

    # Gene annotations (node-level, if available)
    v_call: Optional[str]
    d_call: Optional[str]
    j_call: Optional[str]

    # Sequence variants
    junction: Optional[str]  # CDR3 junction nucleotide
    junction_aa: Optional[str]  # CDR3 junction amino acid
    aa_sequence: Optional[str]  # Full AA sequence

    # Sample/timepoint info
    sample_id: Optional[str]
    timepoint: Optional[str]
    timepoint_id: Optional[str]

    # Confidence
    confidence: Optional[float]
    count: Optional[int]

    # Light chain fields (paired heavy/light data)
    sequence_alignment_light: Optional[str]
    sequence_alignment_light_aa: Optional[str]


# =============================================================================
# Tree Types
# =============================================================================


class OlmstedTree(TypedDict, total=False):
    """
    Phylogenetic tree with ancestral state reconstruction.

    Contains a Newick string representation and detailed node information.
    """

    # Identifiers
    ident: Optional[str]
    tree_id: Optional[str]
    clone_id: Optional[str]

    # Tree structure
    newick: str  # Required: Newick format string
    nodes: List[OlmstedNode] | Dict[str, OlmstedNode]

    # Tree metadata
    type: Optional[str]  # e.g., "pcp.reconstruction", "dnapars"
    sample_id: Optional[str]
    timepoint_ids: Optional[List[str]]

    # Tree statistics
    diversity: Optional[float]  # Mean distance to MRCA
    min_junction_length: Optional[float]
    max_junction_length: Optional[float]

    # Processing info
    downsampling_strategy: Optional[str]


# =============================================================================
# Clone Types
# =============================================================================


class GeneSupport(TypedDict):
    """Per-gene probability support for gene assignment."""

    gene: str
    ident: str
    prob: float


class SampleInfo(TypedDict, total=False):
    """Sample information nested in clone."""

    ident: str
    sample_id: str
    locus: str
    timepoint_id: str


class DatasetInfo(TypedDict, total=False):
    """Dataset information nested in clone."""

    ident: str
    dataset_id: str


class OlmstedClone(TypedDict, total=False):
    """
    Clonal family of sequences from a single rearrangement event.

    Contains gene assignments, CDR positions, and associated trees.
    """

    # Identifiers
    ident: Optional[str]
    clone_id: Optional[str]
    dataset_id: Optional[str]
    sample_id: Optional[str]
    subject_id: Optional[str]

    # Sequence counts
    unique_seqs_count: int  # Required
    total_read_count: Optional[int]
    clone_count: Optional[int]

    # Mutation frequency
    mean_mut_freq: float  # Required

    # Heavy chain gene calls
    v_call: Optional[str]
    d_call: Optional[str]
    j_call: Optional[str]

    # Heavy chain V gene alignment positions
    v_alignment_start: int  # Required
    v_alignment_end: int  # Required
    v_sequence_start: Optional[int]
    v_sequence_end: Optional[int]
    v_germline_start: Optional[int]
    v_germline_end: Optional[int]

    # Heavy chain D gene alignment positions
    d_alignment_start: Optional[int]
    d_alignment_end: Optional[int]
    d_sequence_start: Optional[int]
    d_sequence_end: Optional[int]
    d_germline_start: Optional[int]
    d_germline_end: Optional[int]

    # Heavy chain J gene alignment positions
    j_alignment_start: int  # Required
    j_alignment_end: int  # Required
    j_sequence_start: Optional[int]
    j_sequence_end: Optional[int]
    j_germline_start: Optional[int]
    j_germline_end: Optional[int]

    # Heavy chain CDR positions
    cdr1_alignment_start: Optional[int]
    cdr1_alignment_end: Optional[int]
    cdr2_alignment_start: Optional[int]
    cdr2_alignment_end: Optional[int]

    # Heavy chain junction (CDR3)
    junction_start: Optional[int]
    junction_end: Optional[int]
    junction_length: Optional[int]

    # Germline/naive sequences
    germline_alignment: Optional[str]
    germline_sequence: Optional[str]
    naive_sequence: Optional[str]
    cdr3_sequence: Optional[str]

    # Seed sequence
    has_seed: Optional[bool]
    seed_id: Optional[str]

    # Associated trees
    trees: Optional[List[OlmstedTree]]

    # Unique sequence IDs
    unique_ids: Optional[List[str]]
    timepoint_ids: Optional[List[str]]

    # Gene assignment probabilities
    v_per_gene_support: Optional[List[GeneSupport]]
    d_per_gene_support: Optional[List[GeneSupport]]
    j_per_gene_support: Optional[List[GeneSupport]]

    # Nested references (for webapp compatibility)
    sample: Optional[SampleInfo]
    dataset: Optional[DatasetInfo]

    # === Light chain fields (paired heavy/light data) ===
    is_paired: Optional[bool]
    light_chain_type: Optional[Literal["kappa", "lambda"]]

    # Light chain gene calls
    v_call_light: Optional[str]
    j_call_light: Optional[str]

    # Light chain CDR positions
    cdr1_alignment_start_light: Optional[int]
    cdr1_alignment_end_light: Optional[int]
    cdr2_alignment_start_light: Optional[int]
    cdr2_alignment_end_light: Optional[int]

    # Light chain junction
    junction_start_light: Optional[int]
    junction_length_light: Optional[int]

    # Light chain germline
    germline_alignment_light: Optional[str]

    # Rate scaling (from paired trees)
    rate_scale_heavy: Optional[float]
    rate_scale_light: Optional[float]


# =============================================================================
# Dataset Types
# =============================================================================


class SubjectInfo(TypedDict, total=False):
    """Subject/participant information."""

    ident: str
    subject_id: str
    subject_species: str
    sex: str
    age_min: float
    age_max: float
    age_unit: str
    age_event: str
    ancestry_population: str
    ethnicity: str
    race: str
    strain_name: str


class SampleMetadata(TypedDict, total=False):
    """Sample metadata in dataset."""

    ident: str
    sample_id: str
    subject_id: str
    locus: str
    timepoint_id: str


class TimepointInfo(TypedDict, total=False):
    """Timepoint information."""

    timepoint_id: str
    label: str
    sample_id: str


class ContributorInfo(TypedDict, total=False):
    """Contributor information."""

    contributor_id: str
    name: str
    orcid_id: str
    affiliation: str
    affiliation_ror_id: str


class BuildInfo(TypedDict, total=False):
    """Build/processing information."""

    commit: str
    time: str


class OlmstedDataset(TypedDict, total=False):
    """
    Top-level dataset containing study metadata and sample information.

    This is the root container for an Olmsted data file.
    """

    # Identifiers
    ident: Optional[str]
    dataset_id: str  # Required

    # Schema info
    schema_version: Optional[str]
    type: Optional[str]  # e.g., "pcp.dataset", "airr.dataset"

    # Build info
    build: Optional[BuildInfo]

    # Study metadata
    study_id: Optional[str]
    study_title: Optional[str]
    study_description: Optional[str]
    study_contact: Optional[str]
    inclusion_exclusion_criteria: Optional[str]

    # Contributors and publications
    contributors: Optional[List[ContributorInfo]]
    pub_ids: Optional[List[str]]
    grants: Optional[List[str]]
    keywords_study: Optional[List[str]]

    # Subject info (single subject datasets)
    subject_id: Optional[str]
    subject_species: Optional[str]
    sex: Optional[str]
    age_min: Optional[float]
    age_max: Optional[float]
    age_unit: Optional[str]
    age_event: Optional[str]
    ancestry_population: Optional[str]
    ethnicity: Optional[str]
    race: Optional[str]
    strain_name: Optional[str]

    # Collections
    subjects: Optional[List[SubjectInfo]]
    samples: Optional[List[SampleMetadata]]
    timepoints: Optional[List[TimepointInfo]]

    # Embedded data (for single-file format)
    clones: Optional[List[OlmstedClone]]
    trees: Optional[List[OlmstedTree]]

    # Dataset name (user-provided)
    name: Optional[str]


# =============================================================================
# Consolidated Output Types
# =============================================================================


class OutputMetadata(TypedDict, total=False):
    """Metadata for consolidated output files."""

    format_version: str
    schema_version: str
    source_format: str  # "pcp" or "airr"
    generated_at: str
    generator: str
    generator_version: str


class OlmstedOutput(TypedDict, total=False):
    """
    Complete Olmsted output file structure.

    This represents the consolidated JSON output from olmsted-cli,
    containing all datasets, clones, and trees.

    Example:
        output: OlmstedOutput = {
            "metadata": {...},
            "datasets": [...],
            "clones": {"dataset-1": [...]},
            "trees": [...]
        }
    """

    metadata: Optional[OutputMetadata]
    datasets: List[OlmstedDataset]
    clones: Dict[str, List[OlmstedClone]]  # Keyed by dataset_id
    trees: List[OlmstedTree]


# =============================================================================
# Type aliases for convenience
# =============================================================================

# Node type literal
NodeType = Literal["root", "internal", "leaf", "node"]

# Light chain type literal
LightChainType = Literal["kappa", "lambda"]

# Clones dictionary type
ClonesDict = Dict[str, List[OlmstedClone]]


# =============================================================================
# Export all types
# =============================================================================

__all__ = [
    # Core types
    "OlmstedNode",
    "OlmstedTree",
    "OlmstedClone",
    "OlmstedDataset",
    "OlmstedOutput",
    # Supporting types
    "TimepointMultiplicity",
    "GeneSupport",
    "SampleInfo",
    "DatasetInfo",
    "SubjectInfo",
    "SampleMetadata",
    "TimepointInfo",
    "ContributorInfo",
    "BuildInfo",
    "OutputMetadata",
    # Type aliases
    "NodeType",
    "LightChainType",
    "ClonesDict",
]
