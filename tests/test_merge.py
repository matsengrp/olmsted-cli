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

    # Structural columns (sample_id, family) should be stripped from loaded rows.
    # depth is now retained as a join key (parsed to int) — it's excluded from
    # the *enriched output* in merge_mutations_into_trees, not at load time.
    row = by_family["fam1"][0]
    assert row["site"] == 1
    assert row["parent_aa"] == "K"
    assert row["child_aa"] == "R"
    assert row["surprise_mutsel"] == 4.2
    assert row["log_selection_factor"] == -0.5
    assert row["depth"] == 1  # parsed as int, used as disambiguation key
    assert "sample_id" not in row
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
    stats = merge_mutations_into_trees(trees, by_family)

    assert stats.trees_matched == 1  # Only fam1 has nodes that produce mutations
    assert stats.mutations_enriched == 1  # Only (site=1, K, R) matched a derived mutation
    assert stats.nodes_enriched == 1  # Only the child node received the merge
    # fam99 is in the CSV but not in the JSON → unmatched family
    assert stats.unmatched_families == ["fam99"]
    assert stats.unmatched_family_rows == 1  # one row in fam99
    # (site=99, X, Y) is in the CSV for fam1 but no matching derived mutation
    assert stats.unmatched_mutations == 1

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


def test_merge_command_reports_unmatched(sample_olmsted_json, sample_csv, tmp_path):
    """The merge command should surface unmatched families/mutations as errors."""
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
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0
    # Combined output (status goes to stderr in VerbosePrinter)
    combined = result.stdout + result.stderr
    assert "1 families in the mutations CSV had no matching clone" in combined
    assert "fam99" in combined
    assert "1 CSV mutation records" in combined


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


def test_merge_depth_disambiguation(tmp_path):
    """When the CSV has a `depth` column, the join key includes node depth.

    Build a tree with two independent K→R substitutions at site 1, one at
    depth 1 and one at depth 2. Without depth, both would receive the same
    CSV row's data (broadcast). With depth, each receives the matching row.
    """
    olmsted = {
        "metadata": {"format": "olmsted", "format_version": "1.0"},
        "datasets": [{"dataset_id": "ds", "name": "Depth Test"}],
        "clones": {
            "ds": [
                {
                    "clone_id": "fam1",
                    "dataset_id": "ds",
                    "unique_seqs_count": 3,
                    "mean_mut_freq": 0.0,
                    "sample_id": "s1",
                }
            ]
        },
        "trees": [
            {
                "ident": "tree-1",
                "clone_id": "fam1",
                "newick": "((leaf:0.1)inner:0.1)root;",
                "nodes": [
                    {
                        "sequence_id": "root",
                        "parent": None,
                        "type": "root",
                        "sequence_alignment_aa": "MKT",
                    },
                    {
                        # depth 1: K→R at site 1
                        "sequence_id": "inner",
                        "parent": "root",
                        "type": "internal",
                        "sequence_alignment_aa": "MRT",
                    },
                    {
                        # depth 2: same K→R at site 1 (back-mutation R→K then K→R)
                        # We construct this by setting parent residue back to K so
                        # the diff produces another K→R event.
                        "sequence_id": "leaf",
                        "parent": "back",
                        "type": "leaf",
                        "sequence_alignment_aa": "MRT",
                    },
                    {
                        # depth 2 intermediate that puts K back so the leaf shows K→R again
                        "sequence_id": "back",
                        "parent": "inner",
                        "type": "internal",
                        "sequence_alignment_aa": "MKT",
                    },
                ],
            }
        ],
    }
    # back is at depth 2, leaf is at depth 3 — derived diff: leaf's parent (back, MKT) → leaf (MRT) = K→R at depth 3
    # inner is at depth 1 — derived diff: root (MKT) → inner (MRT) = K→R at depth 1
    csv_text = (
        "family,site,parent_aa,child_aa,depth,score\n"
        "fam1,1,K,R,1,11.1\n"  # matches inner only
        "fam1,1,K,R,3,33.3\n"  # matches leaf only
    )
    json_path = tmp_path / "input.json"
    csv_path = tmp_path / "muts.csv"
    out_path = tmp_path / "out.json"
    json_path.write_text(json.dumps(olmsted))
    csv_path.write_text(csv_text)

    result = subprocess.run(
        [
            "olmsted",
            "merge",
            "-i",
            str(json_path),
            "--mutations",
            str(csv_path),
            "--mutations-use-depth",
            "-o",
            str(out_path),
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, f"merge failed: {result.stderr}"
    combined = result.stdout + result.stderr
    assert "Match mode: site_paa_caa_depth" in combined
    assert "Disambiguation columns in CSV: depth" in combined
    # Each CSV row matches exactly one node-mutation → no broadcasts
    assert "Enriched 2 mutations across 2 nodes" in combined
    assert "Broadcast" not in combined

    out = json.loads(out_path.read_text())
    tree = out["trees"][0]
    by_id = {n["sequence_id"]: n for n in tree["nodes"]}
    inner_mut = next(m for m in by_id["inner"]["mutations"] if m["site"] == 1)
    leaf_mut = next(m for m in by_id["leaf"]["mutations"] if m["site"] == 1)
    assert inner_mut["score"] == 11.1
    assert leaf_mut["score"] == 33.3
    # depth itself should NOT be in the enriched mutation record
    assert "depth" not in inner_mut
    assert "depth" not in leaf_mut


def test_merge_broadcast_detection(tmp_path):
    """Without depth disambiguation, identical substitutions broadcast and are flagged."""
    olmsted = {
        "metadata": {"format": "olmsted", "format_version": "1.0"},
        "datasets": [{"dataset_id": "ds", "name": "Broadcast Test"}],
        "clones": {
            "ds": [
                {
                    "clone_id": "fam1",
                    "dataset_id": "ds",
                    "unique_seqs_count": 2,
                    "mean_mut_freq": 0.0,
                    "sample_id": "s1",
                }
            ]
        },
        "trees": [
            {
                "ident": "tree-1",
                "clone_id": "fam1",
                "newick": "(a,b)root;",
                "nodes": [
                    {
                        "sequence_id": "root",
                        "parent": None,
                        "type": "root",
                        "sequence_alignment_aa": "MKT",
                    },
                    {
                        "sequence_id": "a",
                        "parent": "root",
                        "type": "leaf",
                        "sequence_alignment_aa": "MRT",
                    },
                    {
                        "sequence_id": "b",
                        "parent": "root",
                        "type": "leaf",
                        "sequence_alignment_aa": "MRT",
                    },
                ],
            }
        ],
    }
    # No `depth` column → fall back to (site, parent_aa, child_aa) only
    csv_text = "family,site,parent_aa,child_aa,score\nfam1,1,K,R,5.5\n"
    json_path = tmp_path / "input.json"
    csv_path = tmp_path / "muts.csv"
    out_path = tmp_path / "out.json"
    json_path.write_text(json.dumps(olmsted))
    csv_path.write_text(csv_text)

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
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0
    combined = result.stdout + result.stderr
    # No disambiguation column → not active
    assert "Disambiguation columns in CSV" not in combined
    # The single CSV row matches both 'a' and 'b' → 1 broadcast row
    assert "Enriched 2 mutations across 2 nodes" in combined
    assert "Broadcast: 1 CSV rows matched multiple nodes" in combined
    assert "broadcast to multiple node-mutations" in combined


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


def _name_keyed_fixture(tmp_path):
    """Minimal Olmsted JSON with two nodes that both have a K→R mutation at site 1.

    Without a name column, a CSV keyed on (site, K, R) would broadcast to both.
    With a name column, each CSV row can target exactly one node.
    """
    olmsted = {
        "metadata": {"format": "olmsted", "format_version": "1.0"},
        "datasets": [{"dataset_id": "ds", "name": "Name Test"}],
        "clones": {
            "ds": [
                {
                    "clone_id": "fam1",
                    "dataset_id": "ds",
                    "unique_seqs_count": 3,
                    "mean_mut_freq": 0.0,
                    "sample_id": "s1",
                }
            ]
        },
        "trees": [
            {
                "ident": "tree-1",
                "clone_id": "fam1",
                "newick": "((leaf_a:0.1)inner:0.1)root;",
                "nodes": [
                    {
                        "sequence_id": "root",
                        "parent": None,
                        "type": "root",
                        "sequence_alignment_aa": "MKT",
                    },
                    {
                        "sequence_id": "inner",
                        "parent": "root",
                        "type": "internal",
                        "sequence_alignment_aa": "MRT",  # K→R at site 1
                    },
                    {
                        "sequence_id": "leaf_a",
                        "parent": "root",
                        "type": "leaf",
                        "sequence_alignment_aa": "MRT",  # independent K→R at site 1
                    },
                ],
            }
        ],
    }
    json_path = tmp_path / "input.json"
    json_path.write_text(json.dumps(olmsted))
    return json_path


def test_merge_name_keyed_disambiguation(tmp_path):
    """With a `node_name` column, each CSV row targets one specific node."""
    json_path = _name_keyed_fixture(tmp_path)
    csv_path = tmp_path / "muts.csv"
    out_path = tmp_path / "out.json"
    csv_path.write_text(
        "family,node_name,site,parent_aa,child_aa,score\n"
        "fam1,inner,1,K,R,111\n"
        "fam1,leaf_a,1,K,R,222\n"
    )

    result = subprocess.run(
        ["olmsted", "merge", "-i", str(json_path), "--mutations", str(csv_path),
         "-o", str(out_path)],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, f"merge failed: {result.stderr}"
    combined = result.stdout + result.stderr
    assert "Match mode: name_site" in combined
    assert "Enriched 2 mutations across 2 nodes" in combined
    assert "Broadcast" not in combined
    assert "Integrity mismatches" not in combined

    out = json.loads(out_path.read_text())
    by_id = {n["sequence_id"]: n for n in out["trees"][0]["nodes"]}
    assert next(m for m in by_id["inner"]["mutations"] if m["site"] == 1)["score"] == 111
    assert next(m for m in by_id["leaf_a"]["mutations"] if m["site"] == 1)["score"] == 222
    # node_name is structural — must not leak into the enriched record
    for name in ("inner", "leaf_a"):
        for mut in by_id[name]["mutations"]:
            assert "node_name" not in mut


def test_merge_child_name_alias(tmp_path):
    """`child_name` is accepted as an alias for `node_name`."""
    json_path = _name_keyed_fixture(tmp_path)
    csv_path = tmp_path / "muts.csv"
    out_path = tmp_path / "out.json"
    csv_path.write_text(
        "family,child_name,site,parent_aa,child_aa,score\n"
        "fam1,inner,1,K,R,111\n"
    )

    result = subprocess.run(
        ["olmsted", "merge", "-i", str(json_path), "--mutations", str(csv_path),
         "-o", str(out_path)],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, f"merge failed: {result.stderr}"
    assert "Match mode: name_site" in result.stdout + result.stderr


def test_merge_integrity_mismatch_warns_and_skips(tmp_path):
    """Name+site match but parent_aa/child_aa disagreement is warned and skipped."""
    json_path = _name_keyed_fixture(tmp_path)
    csv_path = tmp_path / "muts.csv"
    out_path = tmp_path / "out.json"
    # Tree has K→R at (inner, site 1); CSV claims K→Q — a mismatch.
    csv_path.write_text(
        "family,node_name,site,parent_aa,child_aa,score\n"
        "fam1,inner,1,K,Q,111\n"
    )

    result = subprocess.run(
        ["olmsted", "merge", "-i", str(json_path), "--mutations", str(csv_path),
         "-o", str(out_path)],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, f"merge failed: {result.stderr}"
    combined = result.stdout + result.stderr
    assert "Integrity mismatches: 1" in combined
    assert "Enriched 0 mutations across 0 nodes" in combined

    out = json.loads(out_path.read_text())
    by_id = {n["sequence_id"]: n for n in out["trees"][0]["nodes"]}
    inner_mut = next(m for m in by_id["inner"]["mutations"] if m["site"] == 1)
    assert "score" not in inner_mut, "Enrichment must not attach on integrity mismatch"


def test_merge_strict_check_fails_on_mismatch(tmp_path):
    """--mutations-strict-check turns integrity mismatches into a hard failure."""
    json_path = _name_keyed_fixture(tmp_path)
    csv_path = tmp_path / "muts.csv"
    out_path = tmp_path / "out.json"
    csv_path.write_text(
        "family,node_name,site,parent_aa,child_aa,score\n"
        "fam1,inner,1,K,Q,111\n"
    )

    result = subprocess.run(
        ["olmsted", "merge", "-i", str(json_path), "--mutations", str(csv_path),
         "--mutations-strict-check", "-o", str(out_path)],
        capture_output=True, text=True,
    )
    assert result.returncode != 0
    combined = result.stdout + result.stderr
    assert "--mutations-strict-check" in combined
    assert "integrity mismatches" in combined.lower()


def test_merge_depth_ignored_without_flag(tmp_path):
    """A `depth` column in the CSV is silently ignored unless --mutations-use-depth is passed."""
    olmsted = {
        "metadata": {"format": "olmsted", "format_version": "1.0"},
        "datasets": [{"dataset_id": "ds", "name": "Depth Opt-in Test"}],
        "clones": {"ds": [{"clone_id": "fam1", "dataset_id": "ds", "unique_seqs_count": 2,
                           "mean_mut_freq": 0.0, "sample_id": "s1"}]},
        "trees": [
            {
                "ident": "tree-1",
                "clone_id": "fam1",
                "newick": "(leaf:0.1)root;",
                "nodes": [
                    {"sequence_id": "root", "parent": None, "type": "root",
                     "sequence_alignment_aa": "MKT"},
                    {"sequence_id": "leaf", "parent": "root", "type": "leaf",
                     "sequence_alignment_aa": "MRT"},  # K→R at site 1, depth 1
                ],
            }
        ],
    }
    json_path = tmp_path / "input.json"
    csv_path = tmp_path / "muts.csv"
    out_path = tmp_path / "out.json"
    json_path.write_text(json.dumps(olmsted))
    # CSV row has depth=99 which would not match any node if depth were used
    csv_path.write_text("family,site,parent_aa,child_aa,depth,score\nfam1,1,K,R,99,111\n")

    result = subprocess.run(
        ["olmsted", "merge", "-i", str(json_path), "--mutations", str(csv_path),
         "-o", str(out_path)],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, f"merge failed: {result.stderr}"
    combined = result.stdout + result.stderr
    # Depth ignored → (site, K, R) matches → enrichment happens
    assert "Match mode: site_paa_caa" in combined
    assert "Match mode: site_paa_caa_depth" not in combined
    assert "Enriched 1 mutations across 1 nodes" in combined

    out = json.loads(out_path.read_text())
    by_id = {n["sequence_id"]: n for n in out["trees"][0]["nodes"]}
    leaf_mut = next(m for m in by_id["leaf"]["mutations"] if m["site"] == 1)
    assert leaf_mut["score"] == 111


@pytest.mark.parametrize("flag", ["--mutations-use-depth", "--mutations-strict-check"])
def test_process_rejects_mutation_flags_without_mutations(tmp_path, flag):
    """`process` argparse rejects mutation-related flags when --mutations is absent.

    Catches a regression where the flags silently no-op. Uses a dummy
    input path — argparse validation should fire before file I/O.
    """
    result = subprocess.run(
        ["olmsted", "process", "-f", "pcp",
         "-i", str(tmp_path / "dummy.csv"),
         "-o", str(tmp_path / "out.json"),
         flag],
        capture_output=True, text=True,
    )
    assert result.returncode != 0
    # argparse errors go to stderr; message mentions the flag requirement
    combined = result.stdout + result.stderr
    assert "--mutations" in combined
    assert "require" in combined.lower()
