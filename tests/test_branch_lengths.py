"""Tests for newick branch-length backfill and the validate advisory.

Covers ``populate_branch_lengths_from_newick`` (the merge-side backfill), the
shared ``_assign_branch_lengths`` core, and the always-on ``validate_tree``
warning that flags a newick carrying branch lengths whose nodes don't.

See issue #29.
"""

import pytest

from olmsted_cli.process_utils import (
    populate_branch_lengths_from_newick,
    validate_tree,
)


def _index(tree):
    return {n["sequence_id"]: n for n in tree["nodes"]}


def _named_internal_tree():
    """VRC26-shape tree: named internal node, branch lengths in the newick,
    bare nodes carrying only sequence_id/type/parent (no length/distance)."""
    return {
        "newick": "((A:0.1,B:0.2)anc1:0.3,naive:0.0)root;",
        "nodes": [
            {"sequence_id": "root", "type": "root", "parent": None},
            {"sequence_id": "anc1", "type": "node", "parent": "root"},
            {"sequence_id": "A", "type": "leaf", "parent": "anc1"},
            {"sequence_id": "B", "type": "leaf", "parent": "anc1"},
            {"sequence_id": "naive", "type": "leaf", "parent": "root"},
        ],
    }


def test_named_internal_round_trip():
    """Branch lengths and cumulative distances are populated from the newick."""
    tree = _named_internal_tree()
    warnings = populate_branch_lengths_from_newick(tree)
    assert warnings == []

    nodes = _index(tree)
    # Branch length to parent.
    assert nodes["A"]["length"] == pytest.approx(0.1, abs=1e-9)
    assert nodes["B"]["length"] == pytest.approx(0.2, abs=1e-9)
    assert nodes["anc1"]["length"] == pytest.approx(0.3, abs=1e-9)
    # Cumulative distance from root: A is anc1(0.3) + 0.1.
    assert nodes["A"]["distance"] == pytest.approx(0.4, abs=1e-9)
    assert nodes["B"]["distance"] == pytest.approx(0.5, abs=1e-9)
    assert nodes["anc1"]["distance"] == pytest.approx(0.3, abs=1e-9)
    # Root is anchored at zero.
    assert nodes["root"]["length"] == 0.0
    assert nodes["root"]["distance"] == 0.0


def test_unnamed_internal_graceful_skip():
    """An unnamed internal node can't be matched: warn and leave it unset."""
    tree = {
        "newick": "((A:0.1,B:0.2):0.3,naive:0.0)root;",
        "nodes": [
            {"sequence_id": "root", "type": "root", "parent": None},
            {"sequence_id": "", "type": "node", "parent": "root"},
            {"sequence_id": "A", "type": "leaf", "parent": ""},
        ],
    }
    warnings = populate_branch_lengths_from_newick(tree)
    assert len(warnings) == 1
    assert "no matching node" in warnings[0]

    nodes = _index(tree)
    # The unmatched node is left untouched rather than raising.
    assert "length" not in nodes[""]
    assert "distance" not in nodes[""]
    # Named nodes are still populated.
    assert nodes["A"]["length"] == pytest.approx(0.1, abs=1e-9)


def test_no_clobber():
    """Values already present on a node are not overwritten."""
    tree = {
        "newick": "(A:0.1,B:0.2)root;",
        "nodes": [
            {"sequence_id": "root", "type": "root", "parent": None},
            {"sequence_id": "A", "type": "leaf", "parent": "root",
             "length": 9.9, "distance": 9.9},
            {"sequence_id": "B", "type": "leaf", "parent": "root"},
        ],
    }
    populate_branch_lengths_from_newick(tree)
    nodes = _index(tree)
    # Pre-existing values survive.
    assert nodes["A"]["length"] == 9.9
    assert nodes["A"]["distance"] == 9.9
    # Missing values are filled.
    assert nodes["B"]["length"] == pytest.approx(0.2, abs=1e-9)


def test_no_branch_lengths_noop():
    """A newick without branch lengths produces no length/distance, no error."""
    tree = {
        "newick": "(A,B)root;",
        "nodes": [
            {"sequence_id": "root", "type": "root", "parent": None},
            {"sequence_id": "A", "type": "leaf", "parent": "root"},
            {"sequence_id": "B", "type": "leaf", "parent": "root"},
        ],
    }
    warnings = populate_branch_lengths_from_newick(tree)
    assert warnings == []
    for node in tree["nodes"]:
        assert "length" not in node
        assert "distance" not in node


def test_distance_monotonic_along_edges():
    """After backfill, every child's distance >= its parent's (additive)."""
    tree = _named_internal_tree()
    populate_branch_lengths_from_newick(tree)
    nodes = _index(tree)
    for node in tree["nodes"]:
        parent_id = node.get("parent")
        if parent_id is None:
            continue
        assert node["distance"] >= nodes[parent_id]["distance"] - 1e-12


def _valid_tree_missing_distances():
    """Schema-valid tree (nodes carry the required sequence fields) whose
    newick has branch lengths but whose nodes lack length/distance."""

    def node(seq_id, ntype, parent):
        return {
            "sequence_id": seq_id,
            "sequence_alignment": "ACGTACGT",
            "sequence_alignment_aa": "MR",
            "type": ntype,
            "parent": parent,
        }

    return {
        "newick": "((A:0.1,B:0.2)anc1:0.3,naive:0.0)root;",
        "nodes": [
            node("root", "root", None),
            node("anc1", "node", "root"),
            node("A", "leaf", "anc1"),
            node("B", "leaf", "anc1"),
            node("naive", "leaf", "root"),
        ],
    }


def test_validate_warns_then_clears():
    """validate_tree warns when the newick has branch lengths the nodes lack;
    the warning clears once they're populated."""
    tree = _valid_tree_missing_distances()

    result = validate_tree(tree)
    assert any("fall back to topological depth" in w for w in result.warnings)
    # Advisory only — it must not make an otherwise-valid tree invalid.
    assert result.ok, result.errors

    populate_branch_lengths_from_newick(tree)
    assert validate_tree(tree).warnings == []
