"""Tests for multi-tree-per-family support in PCP processing.

The tree-CSV parser returns a list per (family, sample_id); multiple
entries become alternate reconstructions in clone.trees[]. Verifies:

- Multiple CSV rows per (family, sample_id) → N tree records.
- tree_id and reconstruction_method columns flow through to tree output.
- Distinct ident values per tree within the same clone.
- Duplicate tree_id within a clone caught by uniqueness enforcement.
- --allow-duplicate-ids downgrades tree_id collision to a warning.
"""

import json
import subprocess
from pathlib import Path

import pytest


# Minimal PCP inputs for these tests — small enough to stay readable.
_PCP_CSV = (
    "sample_id,family,parent_name,child_name,parent_heavy,child_heavy,"
    "branch_length,parent_is_naive,child_is_leaf\n"
    "s1,f1,naive,leaf1,ATCG,ATGG,0.01,True,True\n"
    "s1,f1,naive,leaf2,ATCG,ATGT,0.02,True,True\n"
)


def _write_trees_csv(path, rows):
    """Write a trees CSV with the fields we care about (newick_tree
    quoted to protect commas in the Newick)."""
    with open(path, "w") as f:
        f.write("family_name,sample_id,newick_tree,tree_id,reconstruction_method\n")
        for row in rows:
            f.write(
                f'{row["family_name"]},{row["sample_id"]},'
                f'"{row["newick_tree"]}",{row.get("tree_id", "")},'
                f'{row.get("reconstruction_method", "")}\n'
            )


def _run_process(tmp_path, pcp_path, trees_path, *extra_args, expect_success=True):
    out = tmp_path / "out.json"
    result = subprocess.run(
        [
            "olmsted", "process",
            "-f", "pcp",
            "-i", str(pcp_path),
            "-t", str(trees_path),
            "-o", str(out),
            "--seed", "42",
            "-q",
            *extra_args,
        ],
        capture_output=True, text=True,
    )
    if expect_success:
        assert result.returncode == 0, f"process failed: {result.stderr}"
        return json.loads(out.read_text())
    return result


@pytest.fixture
def pcp_csv(tmp_path):
    path = tmp_path / "pcp.csv"
    path.write_text(_PCP_CSV)
    return path


class TestMultiTreePipeline:
    def test_two_trees_per_family_produce_two_tree_records(self, tmp_path, pcp_csv):
        trees_csv = tmp_path / "trees.csv"
        _write_trees_csv(trees_csv, [
            {
                "family_name": "f1", "sample_id": "s1",
                "newick_tree": "(leaf1:0.01,leaf2:0.02)naive:0.0;",
                "tree_id": "tree-a", "reconstruction_method": "dnapars",
            },
            {
                "family_name": "f1", "sample_id": "s1",
                "newick_tree": "((leaf1:0.005,leaf2:0.015)n1:0.01)naive:0.0;",
                "tree_id": "tree-b", "reconstruction_method": "raxml_ng",
            },
        ])
        data = _run_process(tmp_path, pcp_csv, trees_csv)

        ds_id = data["datasets"][0]["dataset_id"]
        clones = data["clones"][ds_id]
        assert len(clones) == 1, "alternate reconstructions share a clone"

        tree_refs = clones[0]["trees"]
        assert len(tree_refs) == 2
        tree_ids = sorted(t["tree_id"] for t in tree_refs)
        assert tree_ids == ["tree-a", "tree-b"]

        methods = {t.get("reconstruction_method") for t in tree_refs}
        assert methods == {"dnapars", "raxml_ng"}

        # Each tree gets its own ident
        idents = {t["ident"] for t in tree_refs}
        assert len(idents) == 2

        # Top-level trees[] has one entry per tree, each with its own
        # nodes array (the topologies differ so node counts may differ)
        top_trees = data["trees"]
        assert len(top_trees) == 2
        for t in top_trees:
            assert isinstance(t.get("nodes"), list) and len(t["nodes"]) > 0

    def test_csv_tree_id_takes_precedence(self, tmp_path, pcp_csv):
        """When the CSV supplies tree_id, synthesized tree-{family_id}
        fallback is not used."""
        trees_csv = tmp_path / "trees.csv"
        _write_trees_csv(trees_csv, [
            {
                "family_name": "f1", "sample_id": "s1",
                "newick_tree": "(leaf1:0.01,leaf2:0.02)naive:0.0;",
                "tree_id": "my-explicit-tree",
            },
        ])
        data = _run_process(tmp_path, pcp_csv, trees_csv)
        tree_ids = [t["tree_id"] for t in data["trees"]]
        assert tree_ids == ["my-explicit-tree"]

    def test_reconstruction_method_optional(self, tmp_path, pcp_csv):
        """When reconstruction_method column is absent or empty, the
        output field is not present."""
        trees_csv = tmp_path / "trees.csv"
        _write_trees_csv(trees_csv, [
            {
                "family_name": "f1", "sample_id": "s1",
                "newick_tree": "(leaf1:0.01,leaf2:0.02)naive:0.0;",
                "tree_id": "t1",
                # reconstruction_method omitted → empty cell
            },
        ])
        data = _run_process(tmp_path, pcp_csv, trees_csv)
        for t in data["trees"]:
            assert "reconstruction_method" not in t or t["reconstruction_method"] is None

    def test_duplicate_tree_id_within_clone_fails(self, tmp_path, pcp_csv):
        """Two rows with the same tree_id in the same (family, sample_id)
        fail the uniqueness check."""
        trees_csv = tmp_path / "trees.csv"
        _write_trees_csv(trees_csv, [
            {
                "family_name": "f1", "sample_id": "s1",
                "newick_tree": "(leaf1:0.01,leaf2:0.02)naive:0.0;",
                "tree_id": "same-id", "reconstruction_method": "dnapars",
            },
            {
                "family_name": "f1", "sample_id": "s1",
                "newick_tree": "((leaf1:0.005,leaf2:0.015)n1:0.01)naive:0.0;",
                "tree_id": "same-id", "reconstruction_method": "raxml_ng",
            },
        ])
        result = _run_process(tmp_path, pcp_csv, trees_csv, expect_success=False)
        assert result.returncode != 0
        combined = result.stdout + result.stderr
        assert "tree_id" in combined
        assert "same-id" in combined

    def test_duplicate_tree_id_allowed_with_flag(self, tmp_path, pcp_csv):
        """--allow-duplicate-ids downgrades the collision to a warning."""
        trees_csv = tmp_path / "trees.csv"
        _write_trees_csv(trees_csv, [
            {
                "family_name": "f1", "sample_id": "s1",
                "newick_tree": "(leaf1:0.01,leaf2:0.02)naive:0.0;",
                "tree_id": "same-id",
            },
            {
                "family_name": "f1", "sample_id": "s1",
                "newick_tree": "((leaf1:0.005,leaf2:0.015)n1:0.01)naive:0.0;",
                "tree_id": "same-id",
            },
        ])
        data = _run_process(tmp_path, pcp_csv, trees_csv, "--allow-duplicate-ids")
        ds_id = data["datasets"][0]["dataset_id"]
        # Both trees still emitted, both still carry "same-id"
        tree_refs = data["clones"][ds_id][0]["trees"]
        assert len(tree_refs) == 2
        assert all(t["tree_id"] == "same-id" for t in tree_refs)
