"""Tests for newick branch-length backfill and the validate advisory.

Covers ``populate_branch_lengths_from_newick`` (the merge-side backfill), the
shared ``_assign_branch_lengths`` core, and the always-on ``validate_tree``
warning that flags a newick carrying branch lengths whose nodes don't.

See issue #29.
"""

import pytest

from olmsted_cli.process_utils import (
    assign_branch_lengths,
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


def test_empty_sequence_id_graceful_skip():
    """A node whose sequence_id is "" can't be matched (ete skips empty names):
    warn and leave it unset rather than colliding on the empty string."""
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


def test_none_sequence_id_excluded():
    """A node with sequence_id=None is filtered out of the index entirely — no
    match attempt, so no warning and the node is left untouched."""
    tree = {
        "newick": "(A:0.1,B:0.2)root;",
        "nodes": [
            {"sequence_id": "root", "type": "root", "parent": None},
            {"sequence_id": None, "type": "node", "parent": "root"},
            {"sequence_id": "A", "type": "leaf", "parent": "root"},
        ],
    }
    warnings = populate_branch_lengths_from_newick(tree)
    assert warnings == []  # None-id node never enters node_index, so no warning

    none_node = next(n for n in tree["nodes"] if n["sequence_id"] is None)
    assert "length" not in none_node
    assert "distance" not in none_node


def test_reroot_measures_distance_from_naive():
    """assign_branch_lengths measures distance from an already-rerooted
    tree's root — the reroot-on-naive/germline path AIRR ingest uses (see
    process_airr_data._reroot_on_ancestor) before calling this shared core."""
    import ete3

    newick = "((A:0.1,B:0.2)inner:0.3,naive:0.05)root;"
    tree = ete3.PhyloTree(newick, format=1)

    # Reroot on "naive", mirroring process_airr_data._reroot_on_ancestor.
    naive = tree.search_nodes(name="naive")[0]
    tree.set_outgroup(naive)
    tree.remove_child(naive)
    naive.add_child(tree)
    naive.dist = 0
    tree = naive

    nodes = {n.name: {"sequence_id": n.name} for n in tree.traverse() if n.name}
    assign_branch_lengths(tree, nodes, overwrite=True)

    # Naive is the rerooted origin: anchored at zero.
    assert nodes["naive"]["length"] == 0.0
    assert nodes["naive"]["distance"] == 0.0

    # Every other node's distance is measured from naive and is monotonic down
    # the (rerooted) tree — i.e. distances are from-naive, not from-old-root.
    for ete_node in tree.traverse():
        if ete_node.up is None or not ete_node.name or not ete_node.up.name:
            continue
        assert (
            nodes[ete_node.name]["distance"]
            >= nodes[ete_node.up.name]["distance"] - 1e-9
        )


def test_no_clobber():
    """Values already present on a node are not overwritten."""
    tree = {
        "newick": "(A:0.1,B:0.2)root;",
        "nodes": [
            {"sequence_id": "root", "type": "root", "parent": None},
            {
                "sequence_id": "A",
                "type": "leaf",
                "parent": "root",
                "length": 9.9,
                "distance": 9.9,
            },
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
