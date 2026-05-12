"""Multi-tree-per-family + flexible-column-naming integration tests.

Covers issue #23:
- ``parse_pcp_csv`` keys families on ``(sample, family, tree)``.
- ``process_pcp_to_olmsted`` groups composite keys by ``(sample, family)``
  to emit one clone with multiple trees.
- Column auto-detection across the ``sample/sample_id/sample_name``,
  ``family/family_id/family_name``, ``tree/tree_id/tree_name`` families.
- ``--sample-col`` / ``--family-col`` / ``--tree-col`` overrides.
- ``field_metadata.tree`` populated for fields that vary across trees
  within at least one clone.
"""

import json
import subprocess
from pathlib import Path

import pytest

from olmsted_cli.column_resolution import RoleColumnConflict
from olmsted_cli.identifier import IdentMinter
from olmsted_cli.process_pcp_data import (
    PCP_NO_TREE_SENTINEL,
    parse_newick_csv,
    parse_pcp_csv,
    process_pcp_to_olmsted,
)


# Minimal multi-tree PCP CSV: one sample, one family, two reconstruction
# methods. Topology overlap on tip names L1, L2; internal node names
# disjoint (mrca-* vs N#).
PCP_CSV_MULTI = """sample_id,family,tree_name,parent_name,child_name,parent_heavy,child_heavy,parent_is_naive,child_is_leaf,v_gene_heavy,j_gene_heavy
S1,F1,T_a,naive,mrca-1,GAATTCAAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGG,GAATTCAAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGG,true,false,IGHV1*01,IGHJ1*01
S1,F1,T_a,mrca-1,L1,GAATTCAAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGG,GAATTCCAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGG,false,true,IGHV1*01,IGHJ1*01
S1,F1,T_a,mrca-1,L2,GAATTCAAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGG,GAATTCAAATAGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGG,false,true,IGHV1*01,IGHJ1*01
S1,F1,T_b,naive,N1,GAATTCAAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGG,GAATTCAAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGG,true,false,IGHV1*01,IGHJ1*01
S1,F1,T_b,N1,L1,GAATTCAAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGG,GAATTCCAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGG,false,true,IGHV1*01,IGHJ1*01
S1,F1,T_b,N1,L2,GAATTCAAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGG,GAATTCAAATAGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGG,false,true,IGHV1*01,IGHJ1*01
"""

TREES_CSV_MULTI = """family_name,sample_id,tree_name,newick_tree,reconstruction_method
F1,S1,T_a,"((L1:0.2,L2:0.2)mrca-1:0.1)naive;",ground-truth
F1,S1,T_b,"((L1:0.2,L2:0.2)N1:0.1)naive;",bcrlarch-parsimony
"""


def _write(tmp_path: Path, name: str, content: str) -> Path:
    path = tmp_path / name
    path.write_text(content)
    return path


def _process_inline(tmp_path, pcp_text, trees_text, *, name="inline-test", **kwargs):
    """Run the full PCP → Olmsted pipeline in-process on inline CSV strings.

    Returns ``(datasets, clones_dict, trees)`` exactly as
    :func:`process_pcp_to_olmsted` returns them. Use when the test
    checks pipeline output shape, not CLI/argparse behavior — keeps
    runtime in-process rather than shelling out.
    """
    pcp = _write(tmp_path, "input-pcp.csv", pcp_text)
    trees = _write(tmp_path, "input-trees.csv", trees_text)
    pcp_families = parse_pcp_csv(str(pcp), **kwargs)
    newick_trees = parse_newick_csv(str(trees), **kwargs)
    return process_pcp_to_olmsted(
        pcp_families,
        newick_trees,
        minter=IdentMinter(seed=42),
        name=name,
        verbosity=0,
    )


# ---------------------------------------------------------------------------
# parse_pcp_csv composite-key behavior
# ---------------------------------------------------------------------------


def test_parse_pcp_csv_emits_composite_keys(tmp_path):
    csv_path = _write(tmp_path, "input-pcp.csv", PCP_CSV_MULTI)
    families = parse_pcp_csv(csv_path)
    keys = sorted(families.keys())
    assert keys == [("S1", "F1", "T_a"), ("S1", "F1", "T_b")]


def test_parse_pcp_csv_no_tree_column_uses_sentinel(tmp_path):
    csv = "\n".join(
        line for line in PCP_CSV_MULTI.splitlines()
        if line.strip()
    )
    # Strip the tree_name column out
    rows = [r.split(",") for r in csv.splitlines()]
    header_idx = rows[0].index("tree_name")
    stripped = "\n".join(",".join(c for i, c in enumerate(r) if i != header_idx) for r in rows) + "\n"
    csv_path = _write(tmp_path, "input-pcp.csv", stripped)
    families = parse_pcp_csv(csv_path)
    # All rows collapse to the same composite key with the sentinel.
    assert list(families.keys()) == [("S1", "F1", PCP_NO_TREE_SENTINEL)]


def test_parse_newick_csv_emits_composite_keys(tmp_path):
    csv_path = _write(tmp_path, "input-trees.csv", TREES_CSV_MULTI)
    trees = parse_newick_csv(csv_path)
    assert sorted(trees.keys()) == [("S1", "F1", "T_a"), ("S1", "F1", "T_b")]
    # Each composite key carries exactly one tree row.
    for v in trees.values():
        assert len(v) == 1


# ---------------------------------------------------------------------------
# end-to-end: multi-tree-per-family processing
# ---------------------------------------------------------------------------


def test_process_emits_one_clone_with_two_trees(tmp_path):
    datasets, clones_dict, trees = _process_inline(
        tmp_path, PCP_CSV_MULTI, TREES_CSV_MULTI, name="multi-tree-test"
    )
    all_clones = [c for clist in clones_dict.values() for c in clist]
    assert len(all_clones) == 1, f"expected 1 clone, got {len(all_clones)}"
    clone = all_clones[0]
    assert clone["sample_id"] == "S1"
    assert clone["clone_id"] == "S1_F1"
    assert len(clone["trees"]) == 2
    tree_names = sorted(t.get("tree_name") for t in clone["trees"])
    assert tree_names == ["T_a", "T_b"]

    # Top-level trees array carries both reconstructions, with full nodes.
    assert len(trees) == 2
    leaf_names_per_tree = {
        t["tree_name"]: sorted(
            n["sequence_id"]
            for n in t.get("nodes", [])
            if n.get("type") == "leaf"
        )
        for t in trees
    }
    # Tip names (L1, L2) match across both reconstructions.
    assert leaf_names_per_tree["T_a"] == ["L1", "L2"]
    assert leaf_names_per_tree["T_b"] == ["L1", "L2"]


def test_field_metadata_tree_populated_for_multi_tree(tmp_path):
    datasets, _, _ = _process_inline(
        tmp_path, PCP_CSV_MULTI, TREES_CSV_MULTI, name="multi-tree-test"
    )
    fm = datasets[0]["field_metadata"]
    assert "tree" in fm, f"expected field_metadata.tree, got levels {list(fm)}"
    # tree_name and reconstruction_method both differ between T_a/T_b.
    assert "tree_name" in fm["tree"]
    assert "reconstruction_method" in fm["tree"]


# ---------------------------------------------------------------------------
# column-name variants and override flags
# ---------------------------------------------------------------------------


def _swap_column_header(csv_text: str, old: str, new: str) -> str:
    """Replace `old` exactly in the header line."""
    lines = csv_text.splitlines()
    cols = lines[0].split(",")
    cols = [new if c == old else c for c in cols]
    lines[0] = ",".join(cols)
    return "\n".join(lines) + "\n"


def test_family_id_column_variant_auto_detected(tmp_path):
    """``family_id`` (instead of bare ``family``) is auto-detected."""
    _, clones_dict, _ = _process_inline(
        tmp_path,
        _swap_column_header(PCP_CSV_MULTI, "family", "family_id"),
        _swap_column_header(TREES_CSV_MULTI, "family_name", "family_id"),
        name="variant-test",
    )
    all_clones = [c for clist in clones_dict.values() for c in clist]
    assert len(all_clones) == 1
    assert all_clones[0]["clone_id"] == "S1_F1"


def test_conflicting_columns_fail_fast(tmp_path):
    """Two role variants present with disagreeing values raises."""
    # Inject a `family_id` column with a different value than `family`.
    rows = [r.split(",") for r in PCP_CSV_MULTI.strip().splitlines()]
    header = rows[0]
    family_idx = header.index("family")
    header.append("family_id")
    for r in rows[1:]:
        # different value for family_id vs family
        r.append("DIFFERENT")
    csv_text = "\n".join(",".join(r) for r in rows) + "\n"
    pcp = _write(tmp_path, "input-pcp.csv", csv_text)

    # parse_pcp_csv directly rather than the CLI, so the exception type
    # is exposed (CLI prints a friendly error and exits non-zero).
    with pytest.raises(RoleColumnConflict, match="family"):
        parse_pcp_csv(pcp)


def test_override_flag_supersedes_auto_detection(tmp_path):
    """``--family-col`` picks the explicitly-named column on both CSVs.

    Overrides apply uniformly to both the PCP CSV and the trees CSV — the
    user is expected to keep the role-column name consistent across the
    two files.
    """
    # Add `family_id` (with the same values as `family`) to BOTH the PCP
    # CSV and the trees CSV. Then force the override to point at the new
    # column on both.
    pcp_rows = [r.split(",") for r in PCP_CSV_MULTI.strip().splitlines()]
    pcp_header = pcp_rows[0]
    family_idx = pcp_header.index("family")
    pcp_header.append("family_id")
    for r in pcp_rows[1:]:
        r.append(r[family_idx])
    pcp_csv_text = "\n".join(",".join(r) for r in pcp_rows) + "\n"

    trees_rows = [r.split(",") for r in TREES_CSV_MULTI.strip().splitlines()]
    trees_header = trees_rows[0]
    fname_idx = trees_header.index("family_name")
    trees_header.append("family_id")
    for r in trees_rows[1:]:
        r.append(r[fname_idx])
    trees_csv_text = "\n".join(",".join(r) for r in trees_rows) + "\n"

    pcp = _write(tmp_path, "input-pcp.csv", pcp_csv_text)
    trees = _write(tmp_path, "input-trees.csv", trees_csv_text)
    out = tmp_path / "out.json"

    subprocess.run(
        [
            "olmsted", "process", "-f", "pcp",
            "-i", str(pcp), "-t", str(trees),
            "-o", str(out),
            "--family-col", "family_id",
            "--seed", "42", "--name", "override-test", "-q",
        ],
        check=True, capture_output=True,
    )

    data = json.loads(out.read_text())
    all_clones = [c for clist in data["clones"].values() for c in clist]
    # Successful processing under the override → clone_id reflects family_id value.
    assert len(all_clones) == 1
    assert all_clones[0]["clone_id"] == "S1_F1"
