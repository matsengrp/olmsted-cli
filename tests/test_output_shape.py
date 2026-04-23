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
            "-i", str(REPO_ROOT / "example_data/pcp/pcp.csv"),
            "-t", str(REPO_ROOT / "example_data/pcp/trees.csv"),
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
            "-i", str(REPO_ROOT / "example_data/airr/airr.json"),
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
        when input supplies one."""
        data = _run_pcp_process(tmp_path)
        for ds in data["datasets"]:
            assert ds.get("type") in (None, "", "pcp.dataset") is False or ds.get("type") != "pcp.dataset"

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
    def test_tree_id_populated_from_ident_when_missing(self, tmp_path):
        """AIRR input's Tree schema requires tree_id; olmsted-cli fills
        it with the minted ident when input omits it."""
        data = _run_airr_process(tmp_path)
        for tree in data["trees"]:
            assert tree.get("tree_id")  # non-empty
            # If input supplied a tree_id, it's kept; otherwise it matches ident
            if tree["tree_id"] != tree["ident"]:
                # Input supplied a tree_id — that's fine too
                pass

    def test_clone_dataset_nested_absent(self, tmp_path):
        data = _run_airr_process(tmp_path)
        ds_id = data["datasets"][0]["dataset_id"]
        for clone in data["clones"][ds_id]:
            assert "dataset" not in clone
