"""Tests for the merge command and merge_mutations utility."""

import json
import subprocess

import pytest

from olmsted_cli.merge_mutations import (
    derive_node_mutations,
    load_mutations_csv,
    merge_mutations_into_trees,
)


@pytest.fixture
def sample_olmsted_json():
    """Olmsted JSON with two clones, AA sequences on nodes, no mutations arrays."""
    return {
        "metadata": {
            "format": "olmsted",
            "format_version": "1.0",
            "schema_version": "2.0.0",
        },
        "datasets": [{"dataset_id": "test-ds", "name": "Test"}],
        "clones": {
            "test-ds": [
                {
                    "clone_id": "fam1",
                    "dataset_id": "test-ds",
                    "unique_seqs_count": 2,
                    "mean_mut_freq": 0.1,
                    "sample_id": "s1",
                },
                {
                    "clone_id": "fam2",
                    "dataset_id": "test-ds",
                    "unique_seqs_count": 2,
                    "mean_mut_freq": 0.0,
                    "sample_id": "s1",
                },
            ]
        },
        "trees": [
            {
                "ident": "tree-1",
                "clone_id": "fam1",
                "newick": "(child:0.1)root;",
                "nodes": [
                    {
                        "sequence_id": "root",
                        "parent": None,
                        "type": "root",
                        "sequence_alignment_aa": "MKTV",
                    },
                    {
                        "sequence_id": "child",
                        "parent": "root",
                        "type": "leaf",
                        "sequence_alignment_aa": "MRTV",
                    },
                ],
            },
            {
                "ident": "tree-2",
                "clone_id": "fam2",
                "newick": "(child:0.1)root;",
                "nodes": [
                    {
                        "sequence_id": "root",
                        "parent": None,
                        "type": "root",
                        "sequence_alignment_aa": "MKTV",
                    },
                ],
            },
        ],
    }


@pytest.fixture
def sample_csv():
    return (
        "family,site,parent_aa,child_aa,surprise_mutsel,log_selection_factor,sample_id,depth\n"
        "fam1,1,K,R,4.2,-0.5,s1,1\n"
        "fam1,99,X,Y,1.0,0.0,s1,1\n"  # site that doesn't exist on the node
        "fam99,1,A,B,9.9,0.0,s1,1\n"  # family that doesn't exist
    )


def test_derive_node_mutations_basic():
    parent = {"sequence_alignment_aa": "MKTV"}
    child = {"sequence_alignment_aa": "MRTV"}
    muts = derive_node_mutations(child, parent)
    assert muts == [{"site": 1, "parent_aa": "K", "child_aa": "R"}]


def test_derive_node_mutations_skips_gaps():
    parent = {"sequence_alignment_aa": "MK-V"}
    child = {"sequence_alignment_aa": "MRTV"}
    muts = derive_node_mutations(child, parent)
    # Position 1 (K→R) is a real mutation; position 2 (-→T) is a gap, skipped
    assert muts == [{"site": 1, "parent_aa": "K", "child_aa": "R"}]


def test_derive_node_mutations_no_parent():
    child = {"sequence_alignment_aa": "MRTV"}
    assert derive_node_mutations(child, None) == []


def test_load_mutations_csv(tmp_path, sample_csv):
    csv_path = tmp_path / "muts.csv"
    csv_path.write_text(sample_csv)
    by_family = load_mutations_csv(str(csv_path))

    assert set(by_family.keys()) == {"fam1", "fam99"}
    assert len(by_family["fam1"]) == 2

    # Key columns (sample_id, depth, family) should be stripped
    row = by_family["fam1"][0]
    assert row["site"] == 1
    assert row["parent_aa"] == "K"
    assert row["child_aa"] == "R"
    assert row["surprise_mutsel"] == 4.2
    assert row["log_selection_factor"] == -0.5
    assert "sample_id" not in row
    assert "depth" not in row
    assert "family" not in row


def test_load_mutations_csv_missing_required(tmp_path):
    csv_path = tmp_path / "bad.csv"
    csv_path.write_text("family,site\nfam1,1\n")
    with pytest.raises(ValueError, match="missing required columns"):
        load_mutations_csv(str(csv_path))


def test_merge_mutations_into_trees(sample_olmsted_json, sample_csv, tmp_path):
    csv_path = tmp_path / "muts.csv"
    csv_path.write_text(sample_csv)
    by_family = load_mutations_csv(str(csv_path))

    trees = sample_olmsted_json["trees"]
    matched, nodes_with, merged = merge_mutations_into_trees(trees, by_family)

    assert matched == 1  # Only fam1 has nodes that produce mutations
    assert merged == 1  # Only one (site=1, K, R) actually matched a derived mutation

    # Verify the child node now has the enriched mutation
    child_node = next(n for n in trees[0]["nodes"] if n["sequence_id"] == "child")
    assert "mutations" in child_node
    assert len(child_node["mutations"]) == 1
    mut = child_node["mutations"][0]
    assert mut["site"] == 1
    assert mut["parent_aa"] == "K"
    assert mut["child_aa"] == "R"
    assert mut["surprise_mutsel"] == 4.2
    assert mut["log_selection_factor"] == -0.5


def test_merge_command_end_to_end(sample_olmsted_json, sample_csv, tmp_path):
    """Run `olmsted merge` as a subprocess and verify the output JSON."""
    json_path = tmp_path / "input.json"
    csv_path = tmp_path / "muts.csv"
    out_path = tmp_path / "output.json"

    json_path.write_text(json.dumps(sample_olmsted_json))
    csv_path.write_text(sample_csv)

    result = subprocess.run(
        [
            "olmsted",
            "merge",
            "-i",
            str(json_path),
            "--mutations",
            str(csv_path),
            "-o",
            str(out_path),
            "-q",
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, f"merge failed: {result.stderr}"
    assert out_path.exists()

    out = json.loads(out_path.read_text())

    # field_metadata should now include the new mutation field
    fm = out["datasets"][0]["field_metadata"]["mutation"]
    assert "surprise_mutsel" in fm
    assert "log_selection_factor" in fm

    # The child node in fam1's tree should have the merged mutation
    fam1_tree = next(t for t in out["trees"] if t["clone_id"] == "fam1")
    child = next(n for n in fam1_tree["nodes"] if n["sequence_id"] == "child")
    assert child["mutations"][0]["surprise_mutsel"] == 4.2
    assert child["mutations"][0]["log_selection_factor"] == -0.5


def test_merge_command_in_place(sample_olmsted_json, sample_csv, tmp_path):
    json_path = tmp_path / "input.json"
    csv_path = tmp_path / "muts.csv"
    json_path.write_text(json.dumps(sample_olmsted_json))
    csv_path.write_text(sample_csv)

    result = subprocess.run(
        [
            "olmsted",
            "merge",
            "-i",
            str(json_path),
            "--mutations",
            str(csv_path),
            "--in-place",
            "-q",
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, f"merge failed: {result.stderr}"

    out = json.loads(json_path.read_text())
    fam1_tree = next(t for t in out["trees"] if t["clone_id"] == "fam1")
    child = next(n for n in fam1_tree["nodes"] if n["sequence_id"] == "child")
    assert child["mutations"][0]["surprise_mutsel"] == 4.2
