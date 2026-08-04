"""Tests for PCP forest-input handling (olmsted-cli issue 35).

A PCP family is a "forest" when its edges have more than one node that
is a parent but never a child (more than one connected component). This
happens when the input dropped a parent edge for some internal node.
``--on-forest`` controls how olmsted-cli responds: reconcile / drop /
skip (default) / fail.
"""

import pytest

from olmsted_cli.process_pcp_data import (
    ForestTopologyError,
    parse_pcp_csv,
    process_pcp_to_olmsted,
)

# naive -> Node1 -> Leaf1 is the well-formed subtree.
# Node2 -> Leaf2 is an orphan subtree: Node2 is a parent but never a child.
FOREST_PCP_CSV = """sample_id,family,parent_name,child_name,parent_heavy,child_heavy,parent_is_naive,child_is_leaf,v_gene_heavy,j_gene_heavy
S1,F1,naive,Node1,GAATTCAAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGG,GAATTCAAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGG,true,false,IGHV1*01,IGHJ1*01
S1,F1,Node1,Leaf1,GAATTCAAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGG,GAATTCCAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGG,false,true,IGHV1*01,IGHJ1*01
S1,F1,Node2,Leaf2,GAATTCAAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGG,GAATTCAAATAGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGG,false,true,IGHV1*01,IGHJ1*01
"""


def _write_csv(tmp_path, text):
    path = tmp_path / "input-pcp.csv"
    path.write_text(text)
    return path


def _load_forest_families(tmp_path):
    pcp_path = _write_csv(tmp_path, FOREST_PCP_CSV)
    return parse_pcp_csv(str(pcp_path))


def test_reconcile_reattaches_orphan_under_naive(tmp_path):
    families = _load_forest_families(tmp_path)
    _, _, trees = process_pcp_to_olmsted(families, on_forest="reconcile")

    assert len(trees) == 1
    nodes = {n["sequence_id"]: n for n in trees[0]["nodes"]}
    assert set(nodes) == {"naive", "Node1", "Leaf1", "Node2", "Leaf2"}

    roots = [n for n in nodes.values() if n["parent"] is None]
    assert len(roots) == 1
    assert roots[0]["sequence_id"] == "naive"
    assert nodes["Node2"]["parent"] == "naive"


def test_drop_discards_orphan_subtree(tmp_path):
    families = _load_forest_families(tmp_path)
    _, _, trees = process_pcp_to_olmsted(families, on_forest="drop")

    assert len(trees) == 1
    nodes = {n["sequence_id"]: n for n in trees[0]["nodes"]}
    assert set(nodes) == {"naive", "Node1", "Leaf1"}

    roots = [n for n in nodes.values() if n["parent"] is None]
    assert len(roots) == 1
    assert roots[0]["sequence_id"] == "naive"


def test_skip_drops_the_whole_family(tmp_path):
    families = _load_forest_families(tmp_path)
    _, clones_dict, trees = process_pcp_to_olmsted(families, on_forest="skip")

    assert trees == []
    assert all(len(clones) == 0 for clones in clones_dict.values())


def test_skip_is_the_default(tmp_path):
    families = _load_forest_families(tmp_path)
    _, clones_dict, trees = process_pcp_to_olmsted(families)  # no on_forest kwarg

    assert trees == []
    assert all(len(clones) == 0 for clones in clones_dict.values())


def test_fail_raises(tmp_path):
    families = _load_forest_families(tmp_path)
    with pytest.raises(ForestTopologyError, match="Node2"):
        process_pcp_to_olmsted(families, on_forest="fail")


def test_well_formed_family_unaffected_by_on_forest(tmp_path):
    """A single-root family produces identical output under every mode."""
    well_formed_csv = FOREST_PCP_CSV.replace("S1,F1,Node2,Leaf2,", "S1,F1,Node1,Leaf2,")
    pcp_path = _write_csv(tmp_path, well_formed_csv)

    node_sets = {}
    for mode in ("reconcile", "drop", "skip", "fail"):
        families = parse_pcp_csv(str(pcp_path))
        _, _, trees = process_pcp_to_olmsted(families, on_forest=mode)
        node_sets[mode] = frozenset(n["sequence_id"] for n in trees[0]["nodes"])

    assert len(set(node_sets.values())) == 1
    assert node_sets["drop"] == {"naive", "Node1", "Leaf1", "Leaf2"}
