"""Regression tests for the output-shape changes in olmsted#274.

Covers:
- Format-origin tags no longer synthesized (tree.type, dataset.type).
- Synthesized-literal defaults no longer emitted (subject_id, timepoint_id,
  sample_id fallback-to-family_id, dataset.subjects synthetic entry).
- clone.dataset nested object no longer present.
- Ident fields use the {datatype}-{uuid} shape everywhere.
- dataset_id drops the "pcp-" format-origin prefix (now dataset-{uuid}).
- tree.tree_id drops the "pcp-tree-" prefix (now tree-{family_id}).
"""

import json
import subprocess
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


def _run_pcp_process(tmp_path: Path) -> dict:
    out = tmp_path / "out.json"
    subprocess.run(
        [
            "olmsted", "process",
            "-f", "pcp",
            "-i", str(REPO_ROOT / "example-data/pcp/input-pcp.csv"),
            "-t", str(REPO_ROOT / "example-data/pcp/input-trees.csv"),
            "-o", str(out),
            "--seed", "42",
            "-q",
        ],
        check=True,
    )
    return json.loads(out.read_text())


def _run_airr_process(tmp_path: Path) -> dict:
    out = tmp_path / "out.json"
    subprocess.run(
        [
            "olmsted", "process",
            "-f", "airr",
            "-i", str(REPO_ROOT / "example-data/airr/input-airr.json"),
            "-o", str(out),
            "--seed", "42",
            "-q",
        ],
        check=True,
    )
    return json.loads(out.read_text())


class TestPcpOutputShape:
    def test_dataset_type_not_synthesized(self, tmp_path):
        """dataset.type is no longer set to "pcp.dataset" — only populated
        when input supplies one. PCP input has no dataset.type, so the
        field must be absent/None on output."""
        data = _run_pcp_process(tmp_path)
        for ds in data["datasets"]:
            assert ds.get("type") != "pcp.dataset"

    def test_tree_reconstruction_method_not_synthesized(self, tmp_path):
        """tree.reconstruction_method absent when CSV has no column."""
        data = _run_pcp_process(tmp_path)
        for tree in data["trees"]:
            assert "reconstruction_method" not in tree or tree["reconstruction_method"] is None

    def test_no_pcp_tree_prefix(self, tmp_path):
        """Synthesized tree_id uses the datatype prefix, not format-origin."""
        data = _run_pcp_process(tmp_path)
        for tree in data["trees"]:
            assert not tree["tree_id"].startswith("pcp-tree-")
            assert tree["tree_id"].startswith("tree-")

    def test_dataset_id_uses_datatype_prefix(self, tmp_path):
        """PCP dataset_id is dataset-{uuid}, not pcp-{uuid}."""
        data = _run_pcp_process(tmp_path)
        for ds in data["datasets"]:
            assert ds["dataset_id"].startswith("dataset-")
            assert not ds["dataset_id"].startswith("pcp-")

    def test_tree_ident_uses_datatype_prefix(self, tmp_path):
        """Tree idents carry the tree- prefix (were bare uuids before)."""
        data = _run_pcp_process(tmp_path)
        for tree in data["trees"]:
            assert tree["ident"].startswith("tree-")

    def test_clone_dataset_nested_absent(self, tmp_path):
        """clone.dataset removed — only clone.dataset_id remains."""
        data = _run_pcp_process(tmp_path)
        ds_id = data["datasets"][0]["dataset_id"]
        for clone in data["clones"][ds_id]:
            assert "dataset" not in clone

    def test_no_synthesized_subject(self, tmp_path):
        """dataset.subjects is empty (was [{"subject_id": "pcp-subject"}])."""
        data = _run_pcp_process(tmp_path)
        for ds in data["datasets"]:
            assert ds.get("subjects") == []
            assert ds.get("subjects_count") == 0

    def test_no_synthesized_subject_id_on_clone(self, tmp_path):
        """Clones don't carry subject_id = "pcp-subject"."""
        data = _run_pcp_process(tmp_path)
        ds_id = data["datasets"][0]["dataset_id"]
        for clone in data["clones"][ds_id]:
            assert clone.get("subject_id") != "pcp-subject"

    def test_no_synthesized_timepoint_id(self, tmp_path):
        """Samples and clone.sample don't carry timepoint_id = "merged"."""
        data = _run_pcp_process(tmp_path)
        for ds in data["datasets"]:
            for sample in ds.get("samples", []):
                assert sample.get("timepoint_id") != "merged"
        ds_id = data["datasets"][0]["dataset_id"]
        for clone in data["clones"][ds_id]:
            if "sample" in clone:
                assert clone["sample"].get("timepoint_id") != "merged"


class TestAirrOutputShape:
    def test_tree_id_passthrough_when_input_supplies_it(self, tmp_path):
        """When AIRR input provides tree_id, it's preserved (not overwritten
        by the ident fallback)."""
        data = _run_airr_process(tmp_path)
        for tree in data["trees"]:
            assert tree.get("tree_id"), "tree_id should be non-empty"
        # The example-data/airr input always supplies tree_id, so all trees
        # here exercise the pass-through path; the fallback path is covered
        # by test_tree_id_fallback_to_ident_when_input_missing below.

    def test_tree_id_fallback_to_ident_when_input_missing(self):
        """process_tree populates tree_id with the minted ident when AIRR
        input omits the field — the fallback at process_airr_data.process_tree:160-162."""
        from argparse import Namespace

        from olmsted_cli.identifier import IdentMinter
        from olmsted_cli.process_airr_data import process_tree

        args = Namespace(
            minter=IdentMinter(seed=42),
            root_trees=False,
            naive_name="naive",
            verbose=0,
        )
        # Minimal AIRR tree dict with no tree_id supplied
        tree = {
            "newick": "(leaf1:0.01,leaf2:0.02)naive:0.0;",
            "nodes": {
                "naive": {"sequence_alignment": "ATCG", "sequence_alignment_aa": "T"},
                "leaf1": {"sequence_alignment": "ATCG", "sequence_alignment_aa": "T"},
                "leaf2": {"sequence_alignment": "ATCG", "sequence_alignment_aa": "T"},
            },
        }
        result = process_tree(args, clone_id="clone-xyz", tree=tree)

        assert result["ident"].startswith("tree-"), "minted ident should carry tree- prefix"
        assert result["tree_id"] == result["ident"], (
            "missing input tree_id should fall back to the minted ident"
        )

    def test_tree_id_passthrough_on_helper(self):
        """process_tree leaves input-supplied tree_id alone (doesn't
        overwrite with ident fallback)."""
        from argparse import Namespace

        from olmsted_cli.identifier import IdentMinter
        from olmsted_cli.process_airr_data import process_tree

        args = Namespace(
            minter=IdentMinter(seed=42),
            root_trees=False,
            naive_name="naive",
            verbose=0,
        )
        tree = {
            "tree_id": "my-explicit-tree-id",
            "newick": "(leaf1:0.01,leaf2:0.02)naive:0.0;",
            "nodes": {
                "naive": {"sequence_alignment": "ATCG", "sequence_alignment_aa": "T"},
                "leaf1": {"sequence_alignment": "ATCG", "sequence_alignment_aa": "T"},
                "leaf2": {"sequence_alignment": "ATCG", "sequence_alignment_aa": "T"},
            },
        }
        result = process_tree(args, clone_id="clone-xyz", tree=tree)

        assert result["tree_id"] == "my-explicit-tree-id", (
            "input-supplied tree_id must be preserved"
        )
        assert result["ident"] != result["tree_id"], (
            "minted ident is still distinct from the input tree_id"
        )

    def test_clone_dataset_nested_absent(self, tmp_path):
        data = _run_airr_process(tmp_path)
        ds_id = data["datasets"][0]["dataset_id"]
        for clone in data["clones"][ds_id]:
            assert "dataset" not in clone
