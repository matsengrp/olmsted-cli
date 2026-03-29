"""Tests for field_metadata generation."""

import pytest

from olmsted_cli.field_metadata import (
    generate_clone_metadata,
    generate_field_metadata,
    generate_mutation_metadata,
    generate_node_metadata,
    humanize_label,
    infer_field_type,
)


# =============================================================================
# Fixtures
# =============================================================================


@pytest.fixture
def pcp_clones():
    """Mock PCP-style clones with standard fields."""
    return [
        {
            "clone_id": "family-1",
            "ident": "clone-abc",
            "dataset_id": "ds-1",
            "unique_seqs_count": 10,
            "total_read_count": 25,
            "mean_mut_freq": 0.05,
            "v_call": "IGHV3-48*01",
            "d_call": "IGHD3-10*01",
            "j_call": "IGHJ4*02",
            "junction_length": 51,
            "v_alignment_start": 0,
            "v_alignment_end": 294,
            "j_alignment_start": 300,
            "j_alignment_end": 350,
            "d_alignment_start": 294,
            "d_alignment_end": 300,
            "germline_alignment": "ATCG" * 50,
            "has_seed": False,
            "trees": [{"ident": "tree-1"}],
            "sample": {"locus": "igh", "sample_id": "s1"},
            "dataset": {"dataset_id": "ds-1"},
        },
        {
            "clone_id": "family-2",
            "ident": "clone-def",
            "dataset_id": "ds-1",
            "unique_seqs_count": 5,
            "total_read_count": 12,
            "mean_mut_freq": 0.08,
            "v_call": "IGHV1-18*01",
            "d_call": "IGHD2-2*01",
            "j_call": "IGHJ6*02",
            "junction_length": 45,
            "v_alignment_start": 0,
            "v_alignment_end": 290,
            "j_alignment_start": 310,
            "j_alignment_end": 360,
            "d_alignment_start": 290,
            "d_alignment_end": 310,
            "germline_alignment": "GCTA" * 50,
            "has_seed": False,
            "trees": [{"ident": "tree-2"}],
            "sample": {"locus": "igh", "sample_id": "s1"},
            "dataset": {"dataset_id": "ds-1"},
        },
    ]


@pytest.fixture
def trees_with_nodes():
    """Mock trees with node-level data."""
    return [
        {
            "ident": "tree-1",
            "clone_id": "family-1",
            "newick": "(a:0.1,b:0.2)root;",
            "nodes": [
                {
                    "sequence_id": "root",
                    "parent": None,
                    "type": "root",
                    "sequence_alignment": "ATCG",
                    "sequence_alignment_aa": "M",
                    "distance": 0.0,
                    "length": 0.0,
                    "multiplicity": 1,
                    "lbi": 0.5,
                    "lbr": 0.3,
                },
                {
                    "sequence_id": "a",
                    "parent": "root",
                    "type": "leaf",
                    "sequence_alignment": "ATCG",
                    "sequence_alignment_aa": "M",
                    "distance": 0.1,
                    "length": 0.1,
                    "multiplicity": 3,
                    "lbi": 0.8,
                    "lbr": 1.6,
                    "timepoint_id": "day30",
                },
            ],
        }
    ]


@pytest.fixture
def trees_with_surprise():
    """Mock trees with surprise_mutations on nodes."""
    return [
        {
            "ident": "tree-1",
            "clone_id": "family-1",
            "newick": "(a:0.1)root;",
            "nodes": [
                {
                    "sequence_id": "a",
                    "parent": "root",
                    "type": "leaf",
                    "sequence_alignment": "ATCG",
                    "sequence_alignment_aa": "M",
                    "surprise_mutations": [
                        {
                            "site": 10,
                            "parent_aa": "A",
                            "child_aa": "V",
                            "surprise_mutsel": 4.5,
                            "surprise_neutral": 2.1,
                            "selection_contribution": 2.4,
                            "region": "CDR2",
                        },
                        {
                            "site": 25,
                            "parent_aa": "S",
                            "child_aa": "T",
                            "surprise_mutsel": 1.2,
                            "surprise_neutral": 0.8,
                            "selection_contribution": 0.4,
                            "region": "FWR3",
                        },
                    ],
                },
            ],
        }
    ]


# =============================================================================
# Tests: infer_field_type
# =============================================================================


class TestInferFieldType:
    def test_numeric_values(self):
        assert infer_field_type([1, 2, 3]) == "continuous"
        assert infer_field_type([1.0, 2.5, 3.7]) == "continuous"
        assert infer_field_type([1, 2.5, 3]) == "continuous"

    def test_string_values(self):
        assert infer_field_type(["hello", "world"]) == "categorical"
        assert infer_field_type(["CDR1", "FWR2", "CDR3"]) == "categorical"

    def test_aa_values(self):
        # Values with AA-only chars (not in DNA alphabet) → aa
        assert infer_field_type(["A", "V", "L", "M"]) == "aa"
        assert infer_field_type(["S", "R", "K", "D", "E"]) == "aa"
        assert infer_field_type(["*", "-"]) == "aa"

    def test_dna_values(self):
        # Pure ACGT/U → dna (ambiguous with AA, but DNA is preferred)
        assert infer_field_type(["A", "C", "G", "T"]) == "dna"
        assert infer_field_type(["A", "T", "G", "N"]) == "dna"

    def test_boolean_values(self):
        assert infer_field_type([True, False, True]) == "categorical"

    def test_mixed_values(self):
        assert infer_field_type([1, "a", 2]) == "tooltip"

    def test_empty_values(self):
        assert infer_field_type([]) == "tooltip"

    def test_complex_values(self):
        assert infer_field_type([{"a": 1}, {"b": 2}]) == "tooltip"
        assert infer_field_type([[1, 2], [3, 4]]) == "tooltip"


# =============================================================================
# Tests: humanize_label
# =============================================================================


class TestHumanizeLabel:
    def test_simple_field(self):
        assert humanize_label("distance") == "Distance"

    def test_multi_word(self):
        assert humanize_label("unique_seqs_count") == "Unique Sequences Count"

    def test_abbreviations(self):
        assert humanize_label("lbi") == "LBI"
        assert humanize_label("lbr") == "LBR"
        assert humanize_label("mean_mut_freq") == "Mean Mutation Frequency"

    def test_compound_abbreviation(self):
        assert humanize_label("mean_surprise_mutsel") == "Mean Surprise MutSel"

    def test_single_letter_abbreviations(self):
        assert humanize_label("v_call") == "V Call"
        assert humanize_label("d_call") == "D Call"


# =============================================================================
# Tests: generate_clone_metadata
# =============================================================================


class TestGenerateCloneMetadata:
    def test_standard_pcp_fields(self, pcp_clones):
        meta = generate_clone_metadata(pcp_clones)
        assert "unique_seqs_count" in meta
        assert meta["unique_seqs_count"]["type"] == "continuous"
        assert meta["v_call"]["type"] == "categorical"
        assert meta["mean_mut_freq"]["type"] == "continuous"
        assert meta["junction_length"]["type"] == "continuous"

    def test_locus_extracted_from_sample(self, pcp_clones):
        meta = generate_clone_metadata(pcp_clones)
        assert "locus" in meta
        assert meta["locus"]["type"] == "categorical"
        assert meta["locus"]["label"] == "Locus"

    def test_excluded_fields_not_present(self, pcp_clones):
        meta = generate_clone_metadata(pcp_clones)
        excluded = [
            "clone_id", "ident", "dataset_id", "dataset", "sample", "trees",
            "germline_alignment", "v_alignment_start", "v_alignment_end",
            "j_alignment_start", "j_alignment_end",
        ]
        for field in excluded:
            assert field not in meta, f"Excluded field '{field}' should not be in metadata"

    def test_custom_fields_override(self, pcp_clones):
        custom = [
            {"name": "unique_seqs_count", "level": "clone", "type": "tooltip", "label": "Custom Label"},
        ]
        meta = generate_clone_metadata(pcp_clones, custom_fields=custom)
        assert meta["unique_seqs_count"]["type"] == "tooltip"
        assert meta["unique_seqs_count"]["label"] == "Custom Label"

    def test_custom_fields_add_new(self, pcp_clones):
        custom = [
            {"name": "my_metric", "level": "clone", "type": "continuous", "label": "My Metric"},
        ]
        meta = generate_clone_metadata(pcp_clones, custom_fields=custom)
        assert "my_metric" in meta
        assert meta["my_metric"]["type"] == "continuous"

    def test_custom_fields_wrong_level_ignored(self, pcp_clones):
        custom = [
            {"name": "node_metric", "level": "node", "type": "continuous", "label": "Node Metric"},
        ]
        meta = generate_clone_metadata(pcp_clones, custom_fields=custom)
        assert "node_metric" not in meta

    def test_empty_clones(self):
        meta = generate_clone_metadata([])
        assert meta == {}

    def test_unknown_fields_inferred(self):
        clones = [
            {"clone_id": "c1", "my_number": 42, "my_string": "hello"},
        ]
        meta = generate_clone_metadata(clones)
        assert "my_number" in meta
        assert meta["my_number"]["type"] == "continuous"
        assert "my_string" in meta
        assert meta["my_string"]["type"] == "categorical"


# =============================================================================
# Tests: generate_node_metadata
# =============================================================================


class TestGenerateNodeMetadata:
    def test_standard_node_fields(self, trees_with_nodes):
        meta = generate_node_metadata(trees_with_nodes)
        assert "lbi" in meta
        assert meta["lbi"]["type"] == "continuous"
        assert meta["lbi"]["label"] == "LBI"
        assert "lbr" in meta
        assert "multiplicity" in meta
        assert "distance" in meta

    def test_excluded_node_fields(self, trees_with_nodes):
        meta = generate_node_metadata(trees_with_nodes)
        excluded = ["sequence_id", "parent", "type", "sequence_alignment", "sequence_alignment_aa"]
        for field in excluded:
            assert field not in meta

    def test_metrics_conditional(self):
        """LBI/LBR only appear when present on nodes."""
        trees_no_metrics = [
            {
                "nodes": [
                    {"sequence_id": "a", "multiplicity": 3, "distance": 0.1},
                ]
            }
        ]
        meta = generate_node_metadata(trees_no_metrics)
        assert "lbi" not in meta
        assert "lbr" not in meta
        assert "multiplicity" in meta

    def test_empty_trees(self):
        meta = generate_node_metadata([])
        assert meta == {}


# =============================================================================
# Tests: generate_mutation_metadata
# =============================================================================


class TestGenerateMutationMetadata:
    def test_surprise_fields(self, trees_with_surprise):
        meta = generate_mutation_metadata(trees_with_surprise)
        assert "surprise_mutsel" in meta
        assert meta["surprise_mutsel"]["type"] == "continuous"
        assert meta["surprise_mutsel"]["label"] == "Surprise (MutSel)"
        assert "surprise_neutral" in meta
        assert "selection_contribution" in meta
        assert "region" in meta
        assert meta["region"]["type"] == "categorical"

    def test_aa_fields(self, trees_with_surprise):
        meta = generate_mutation_metadata(trees_with_surprise)
        assert "parent_aa" in meta
        assert meta["parent_aa"]["type"] == "aa"
        assert meta["child_aa"]["type"] == "aa"

    def test_ranges_on_continuous_fields(self, trees_with_surprise):
        meta = generate_mutation_metadata(trees_with_surprise)
        assert "range" in meta["surprise_mutsel"]
        r = meta["surprise_mutsel"]["range"]
        assert len(r) == 2
        assert r[0] <= r[1]
        # region is categorical, no range
        assert "range" not in meta["region"]

    def test_excluded_mutation_fields(self, trees_with_surprise):
        meta = generate_mutation_metadata(trees_with_surprise)
        assert "site" not in meta

    def test_no_surprise_data(self, trees_with_nodes):
        meta = generate_mutation_metadata(trees_with_nodes)
        assert meta == {}

    def test_custom_mutation_fields_without_data(self):
        custom = [
            {"name": "custom_score", "level": "mutation", "type": "continuous", "label": "Custom Score"},
        ]
        meta = generate_mutation_metadata([], custom_fields=custom)
        assert "custom_score" in meta


# =============================================================================
# Tests: generate_field_metadata (top-level)
# =============================================================================


class TestGenerateFieldMetadata:
    def test_all_levels_present(self, pcp_clones, trees_with_surprise):
        meta = generate_field_metadata(pcp_clones, trees_with_surprise)
        assert "clone" in meta
        assert "mutation" in meta

    def test_empty_levels_omitted(self, pcp_clones):
        meta = generate_field_metadata(pcp_clones, trees=[])
        assert "clone" in meta
        assert "node" not in meta
        assert "branch" not in meta
        assert "mutation" not in meta

    def test_none_trees_handled(self, pcp_clones):
        meta = generate_field_metadata(pcp_clones, trees=None)
        assert "clone" in meta

    def test_custom_fields_distributed(self, pcp_clones, trees_with_nodes):
        custom = [
            {"name": "my_clone_field", "level": "clone", "type": "continuous", "label": "My Clone"},
            {"name": "my_node_field", "level": "node", "type": "categorical", "label": "My Node"},
        ]
        meta = generate_field_metadata(pcp_clones, trees_with_nodes, custom_fields=custom)
        assert "my_clone_field" in meta["clone"]
        assert "my_node_field" in meta["node"]
