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
    # only_listed defaults to False → no derived mutations are dropped
    assert stats.mutations_dropped == 0

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


def test_merge_mutations_into_trees_only_listed_stat(tmp_path):
    """`mutations_dropped` counts derived mutations removed under only_listed.

    Unit-level coverage of the new stat: a leaf has two derived mutations
    (sites 1 and 2) but the CSV lists only site 1. With only_listed=True
    the second mutation is dropped and the count surfaces in stats.
    """
    trees = [
        {
            "ident": "tree-1",
            "clone_id": "fam1",
            "nodes": [
                {"sequence_id": "root", "parent": None,
                 "sequence_alignment_aa": "MQQ"},
                {"sequence_id": "leaf", "parent": "root",
                 "sequence_alignment_aa": "MKR"},  # Q→K at site 1, Q→R at site 2
            ],
        }
    ]
    csv_path = tmp_path / "muts.csv"
    csv_path.write_text("family,site,parent_aa,child_aa,score\nfam1,1,Q,K,9.9\n")
    by_family = load_mutations_csv(str(csv_path))

    stats = merge_mutations_into_trees(trees, by_family, only_listed=True)

    assert stats.mutations_enriched == 1
    assert stats.mutations_dropped == 1
    leaf = next(n for n in trees[0]["nodes"] if n["sequence_id"] == "leaf")
    assert [m["site"] for m in leaf["mutations"]] == [1]
    assert leaf["mutations"][0]["score"] == 9.9


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


def test_merge_integrity_mismatch_fails_by_default(tmp_path):
    """By default, any parent_aa/child_aa disagreement is a hard failure.

    The merge skips the row (never attaches wrong-looking data) AND exits
    non-zero so callers can't accidentally ship a partially-wrong merge.
    """
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
    assert result.returncode != 0
    combined = result.stdout + result.stderr
    assert "integrity mismatches" in combined.lower()
    assert "--mutations-allow-mismatch" in combined


def test_merge_allow_mismatch_downgrades_to_warning(tmp_path):
    """--mutations-allow-mismatch downgrades integrity mismatches to a warning.

    The row is still skipped — the flag never attaches wrong-looking data,
    it only stops the non-zero exit.
    """
    json_path = _name_keyed_fixture(tmp_path)
    csv_path = tmp_path / "muts.csv"
    out_path = tmp_path / "out.json"
    csv_path.write_text(
        "family,node_name,site,parent_aa,child_aa,score\n"
        "fam1,inner,1,K,Q,111\n"
    )

    result = subprocess.run(
        ["olmsted", "merge", "-i", str(json_path), "--mutations", str(csv_path),
         "--mutations-allow-mismatch", "-o", str(out_path)],
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


def test_merge_use_depth_flag_without_depth_column_fails(tmp_path):
    """--mutations-use-depth is a misuse signal when the CSV has no depth column."""
    json_path = _name_keyed_fixture(tmp_path)
    csv_path = tmp_path / "muts.csv"
    out_path = tmp_path / "out.json"
    # No 'depth' column at all
    csv_path.write_text(
        "family,node_name,site,parent_aa,child_aa,score\n"
        "fam1,inner,1,K,R,111\n"
    )

    result = subprocess.run(
        ["olmsted", "merge", "-i", str(json_path), "--mutations", str(csv_path),
         "--mutations-use-depth", "-o", str(out_path)],
        capture_output=True, text=True,
    )
    assert result.returncode != 0
    combined = result.stdout + result.stderr
    assert "--mutations-use-depth" in combined
    assert "no 'depth' column" in combined.lower()


def test_merge_name_keyed_depth_ignored_without_flag(tmp_path):
    """In name-keyed mode, a wrong `depth` value is ignored unless --mutations-use-depth is passed.

    Depth arithmetic depends on upstream rooting conventions the CLI can't
    infer with certainty, so depth is opt-in for BOTH match-key and
    integrity-check use. Without the flag, a CSV depth disagreement should
    not trigger an integrity mismatch.
    """
    json_path = _name_keyed_fixture(tmp_path)
    csv_path = tmp_path / "muts.csv"
    out_path = tmp_path / "out.json"
    # Tree has K→R at (inner, site 1) at depth 1. CSV claims depth=99 — would
    # be an integrity mismatch if depth were being checked.
    csv_path.write_text(
        "family,node_name,site,parent_aa,child_aa,depth,score\n"
        "fam1,inner,1,K,R,99,111\n"
    )

    result = subprocess.run(
        ["olmsted", "merge", "-i", str(json_path), "--mutations", str(csv_path),
         "-o", str(out_path)],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, f"merge failed: {result.stderr}"
    combined = result.stdout + result.stderr
    # Without the flag: depth ignored, enrichment proceeds, no integrity mismatch
    assert "Integrity mismatches" not in combined
    assert "Enriched 1 mutations across 1 nodes" in combined

    out = json.loads(out_path.read_text())
    by_id = {n["sequence_id"]: n for n in out["trees"][0]["nodes"]}
    inner_mut = next(m for m in by_id["inner"]["mutations"] if m["site"] == 1)
    assert inner_mut["score"] == 111


def test_merge_name_keyed_depth_check_with_flag(tmp_path):
    """With --mutations-use-depth, name-keyed mode checks depth as an integrity field.

    A depth disagreement triggers the same default-fail behavior as a
    parent_aa/child_aa disagreement.
    """
    json_path = _name_keyed_fixture(tmp_path)
    csv_path = tmp_path / "muts.csv"
    out_path = tmp_path / "out.json"
    csv_path.write_text(
        "family,node_name,site,parent_aa,child_aa,depth,score\n"
        "fam1,inner,1,K,R,99,111\n"  # depth=99 disagrees with tree depth 1
    )

    result = subprocess.run(
        ["olmsted", "merge", "-i", str(json_path), "--mutations", str(csv_path),
         "--mutations-use-depth", "-o", str(out_path)],
        capture_output=True, text=True,
    )
    # With the flag: depth checked, mismatch detected → default-fail
    assert result.returncode != 0
    combined = result.stdout + result.stderr
    assert "Integrity mismatches: 1" in combined
    assert "--mutations-allow-mismatch" in combined


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


def test_only_listed_drops_unlisted_derived_mutations(tmp_path):
    """--only-listed-mutations: derived mutations not in the CSV are dropped.

    Reproduces the scenario from issue #18: a leaf has two derived
    mutations (sites 1 and 2), but the CSV only lists site 1. Without
    the flag, both appear in the output (site 2 unannotated). With the
    flag, only site 1 survives.
    """
    olmsted = {
        "metadata": {"format": "olmsted", "format_version": "1.0"},
        "datasets": [{"dataset_id": "ds", "name": "Only-Listed Test"}],
        "clones": {"ds": [{"clone_id": "fam1", "dataset_id": "ds",
                           "unique_seqs_count": 2, "mean_mut_freq": 0.0,
                           "sample_id": "s1"}]},
        "trees": [
            {
                "ident": "tree-1",
                "clone_id": "fam1",
                "newick": "(leaf:0.1)root;",
                "nodes": [
                    {"sequence_id": "root", "parent": None, "type": "root",
                     "sequence_alignment_aa": "MQQ"},
                    # K at site 1 (Q→K), R at site 2 (Q→R) — two derived mutations
                    {"sequence_id": "leaf", "parent": "root", "type": "leaf",
                     "sequence_alignment_aa": "MKR"},
                ],
            }
        ],
    }
    json_path = tmp_path / "input.json"
    csv_path = tmp_path / "muts.csv"
    out_path_default = tmp_path / "out_default.json"
    out_path_filtered = tmp_path / "out_filtered.json"
    json_path.write_text(json.dumps(olmsted))
    # CSV lists only the site-1 mutation
    csv_path.write_text(
        "family,site,parent_aa,child_aa,score\nfam1,1,Q,K,9.9\n"
    )

    # Default behavior: site-2 derived mutation passes through unannotated
    result = subprocess.run(
        ["olmsted", "merge", "-i", str(json_path), "--mutations", str(csv_path),
         "-o", str(out_path_default)],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    out = json.loads(out_path_default.read_text())
    leaf = next(n for n in out["trees"][0]["nodes"] if n["sequence_id"] == "leaf")
    sites = sorted(m["site"] for m in leaf["mutations"])
    assert sites == [1, 2], "Default keeps both derived mutations"
    site2 = next(m for m in leaf["mutations"] if m["site"] == 2)
    assert "score" not in site2, "Unlisted mutation comes through unannotated"

    # With --only-listed-mutations: site-2 derived mutation is dropped
    result = subprocess.run(
        ["olmsted", "merge", "-i", str(json_path), "--mutations", str(csv_path),
         "--only-listed-mutations", "-o", str(out_path_filtered)],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    combined = result.stdout + result.stderr
    assert "Dropped 1 derived mutations" in combined
    out = json.loads(out_path_filtered.read_text())
    leaf = next(n for n in out["trees"][0]["nodes"] if n["sequence_id"] == "leaf")
    assert len(leaf["mutations"]) == 1
    assert leaf["mutations"][0]["site"] == 1
    assert leaf["mutations"][0]["score"] == 9.9


def test_only_listed_name_keyed(tmp_path):
    """--only-listed-mutations works in name-keyed mode too.

    The leaf node has two derived mutations (K→R at site 1, T→R at site 2)
    but the CSV only lists site 1. Under --only-listed-mutations, site 2
    is dropped even though it would have been derived and emitted bare
    by the default merge.
    """
    olmsted = {
        "metadata": {"format": "olmsted", "format_version": "1.0"},
        "datasets": [{"dataset_id": "ds", "name": "Name Only-Listed Test"}],
        "clones": {"ds": [{"clone_id": "fam1", "dataset_id": "ds",
                           "unique_seqs_count": 2, "mean_mut_freq": 0.0,
                           "sample_id": "s1"}]},
        "trees": [
            {
                "ident": "tree-1",
                "clone_id": "fam1",
                "newick": "(leaf:0.1)root;",
                "nodes": [
                    {"sequence_id": "root", "parent": None, "type": "root",
                     "sequence_alignment_aa": "MKT"},
                    # K→R at site 1, T→R at site 2 — two derived mutations
                    {"sequence_id": "leaf", "parent": "root", "type": "leaf",
                     "sequence_alignment_aa": "MRR"},
                ],
            }
        ],
    }
    json_path = tmp_path / "input.json"
    csv_path = tmp_path / "muts.csv"
    out_path = tmp_path / "out.json"
    json_path.write_text(json.dumps(olmsted))
    csv_path.write_text(
        "family,node_name,site,parent_aa,child_aa,score\n"
        "fam1,leaf,1,K,R,111\n"
    )

    result = subprocess.run(
        ["olmsted", "merge", "-i", str(json_path), "--mutations", str(csv_path),
         "--only-listed-mutations", "-o", str(out_path)],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    combined = result.stdout + result.stderr
    assert "Match mode: name_site" in combined
    assert "Dropped 1 derived mutations" in combined

    out = json.loads(out_path.read_text())
    leaf = next(n for n in out["trees"][0]["nodes"] if n["sequence_id"] == "leaf")
    assert len(leaf["mutations"]) == 1
    assert leaf["mutations"][0]["site"] == 1
    assert leaf["mutations"][0]["score"] == 111


def test_only_listed_leaves_unmatched_families_alone(tmp_path, sample_olmsted_json,
                                                     sample_csv):
    """Trees whose family is absent from the CSV pass through untouched.

    The CSV only mentions fam1 (and fam99, which has no tree). Pre-existing
    mutations on fam2's tree must survive --only-listed-mutations untouched
    — the filter is scoped to CSV-matched trees only.
    """
    # Pre-populate fam2's child node with mutations to confirm they survive.
    fam2_tree = next(t for t in sample_olmsted_json["trees"] if t["clone_id"] == "fam2")
    fam2_tree["nodes"].append({
        "sequence_id": "child2",
        "parent": "root",
        "type": "leaf",
        "sequence_alignment_aa": "MRTV",
        "mutations": [{"site": 1, "parent_aa": "K", "child_aa": "R"}],
    })

    json_path = tmp_path / "input.json"
    csv_path = tmp_path / "muts.csv"
    out_path = tmp_path / "out.json"
    json_path.write_text(json.dumps(sample_olmsted_json))
    csv_path.write_text(sample_csv)

    result = subprocess.run(
        ["olmsted", "merge", "-i", str(json_path), "--mutations", str(csv_path),
         "--only-listed-mutations", "-o", str(out_path)],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    out = json.loads(out_path.read_text())
    # fam1 IS in the CSV → its child mutation receives the CSV enrichment.
    # Verifying both halves: filtering ran on fam1, was a no-op on fam2.
    fam1_tree_out = next(t for t in out["trees"] if t["clone_id"] == "fam1")
    fam1_child = next(n for n in fam1_tree_out["nodes"] if n["sequence_id"] == "child")
    assert fam1_child["mutations"] == [
        {"site": 1, "parent_aa": "K", "child_aa": "R",
         "surprise_mutsel": 4.2, "log_selection_factor": -0.5}
    ]
    fam2_tree_out = next(t for t in out["trees"] if t["clone_id"] == "fam2")
    child2 = next(n for n in fam2_tree_out["nodes"] if n["sequence_id"] == "child2")
    # fam2 isn't in the CSV → its pre-existing mutations are not filtered
    assert child2["mutations"] == [{"site": 1, "parent_aa": "K", "child_aa": "R"}]


def test_only_listed_filters_preexisting_upstream_mutations(tmp_path):
    """Pre-existing upstream mutations on a CSV-matched family are filtered too.

    The flag's contract is "the CSV is authoritative for which mutations
    appear" — applies regardless of whether the array on the node came
    from sequence-diff derivation or from an upstream pipeline that
    pre-populated `mutations`. Counterpart to
    `test_only_listed_leaves_unmatched_families_alone`, which covers the
    *unmatched*-family case.
    """
    olmsted = {
        "metadata": {"format": "olmsted", "format_version": "1.0"},
        "datasets": [{"dataset_id": "ds", "name": "Pre-existing Filter Test"}],
        "clones": {"ds": [{"clone_id": "fam1", "dataset_id": "ds",
                           "unique_seqs_count": 2, "mean_mut_freq": 0.0,
                           "sample_id": "s1"}]},
        "trees": [
            {
                "ident": "tree-1",
                "clone_id": "fam1",
                "newick": "(leaf:0.1)root;",
                "nodes": [
                    {"sequence_id": "root", "parent": None, "type": "root",
                     "sequence_alignment_aa": "MQQQ"},
                    # AA sequence still shows three changes vs. parent, but
                    # the upstream pipeline pre-populated only two of them
                    # — exercising the "existing array, no derive" path.
                    {"sequence_id": "leaf", "parent": "root", "type": "leaf",
                     "sequence_alignment_aa": "MKRS",
                     "mutations": [
                         {"site": 1, "parent_aa": "Q", "child_aa": "K"},
                         {"site": 2, "parent_aa": "Q", "child_aa": "R"},
                     ]},
                ],
            }
        ],
    }
    json_path = tmp_path / "input.json"
    csv_path = tmp_path / "muts.csv"
    out_path = tmp_path / "out.json"
    json_path.write_text(json.dumps(olmsted))
    # CSV lists only site 1; site 2 is a pre-existing entry the user has
    # no opinion about and wants filtered out.
    csv_path.write_text(
        "family,site,parent_aa,child_aa,score\nfam1,1,Q,K,9.9\n"
    )

    result = subprocess.run(
        ["olmsted", "merge", "-i", str(json_path), "--mutations", str(csv_path),
         "--only-listed-mutations", "-o", str(out_path)],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    out = json.loads(out_path.read_text())
    leaf = next(n for n in out["trees"][0]["nodes"] if n["sequence_id"] == "leaf")
    # Site-1 enriched, site-2 pre-existing entry filtered out.
    assert len(leaf["mutations"]) == 1
    assert leaf["mutations"][0]["site"] == 1
    assert leaf["mutations"][0]["score"] == 9.9


@pytest.mark.parametrize("mode", ["site_keyed", "name_keyed"])
def test_only_listed_deletes_empty_mutations_array(tmp_path, mode):
    """When every derived mutation on a node is unlisted, drop the key entirely.

    Both merge modes should `del node["mutations"]` rather than leave an
    empty list behind, so a node that loses all its events looks the
    same as one that never had any.
    """
    olmsted = {
        "metadata": {"format": "olmsted", "format_version": "1.0"},
        "datasets": [{"dataset_id": "ds", "name": "Empty-Result Test"}],
        "clones": {"ds": [{"clone_id": "fam1", "dataset_id": "ds",
                           "unique_seqs_count": 2, "mean_mut_freq": 0.0,
                           "sample_id": "s1"}]},
        "trees": [
            {
                "ident": "tree-1",
                "clone_id": "fam1",
                "newick": "(leaf:0.1)root;",
                "nodes": [
                    {"sequence_id": "root", "parent": None, "type": "root",
                     "sequence_alignment_aa": "MQQ"},
                    {"sequence_id": "leaf", "parent": "root", "type": "leaf",
                     "sequence_alignment_aa": "MKR"},
                ],
            }
        ],
    }
    json_path = tmp_path / "input.json"
    csv_path = tmp_path / "muts.csv"
    out_path = tmp_path / "out.json"
    json_path.write_text(json.dumps(olmsted))
    # CSV is for fam1 (so the tree is matched and filtering runs) but
    # lists a site that isn't in the tree — every derived mutation is
    # unlisted. Mode is selected by whether `node_name` is in the header.
    if mode == "name_keyed":
        csv_path.write_text(
            "family,node_name,site,parent_aa,child_aa,score\n"
            "fam1,leaf,99,X,Y,1.0\n"
        )
    else:
        csv_path.write_text(
            "family,site,parent_aa,child_aa,score\nfam1,99,X,Y,1.0\n"
        )

    result = subprocess.run(
        ["olmsted", "merge", "-i", str(json_path), "--mutations", str(csv_path),
         "--only-listed-mutations", "-o", str(out_path)],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    out = json.loads(out_path.read_text())
    leaf = next(n for n in out["trees"][0]["nodes"] if n["sequence_id"] == "leaf")
    assert "mutations" not in leaf, (
        f"Expected `mutations` key to be deleted when all entries are dropped; "
        f"got: {leaf.get('mutations')!r}"
    )


def test_only_listed_drops_integrity_mismatched_sites(tmp_path):
    """Integrity mismatch + allow-mismatch + only-listed cascades into a drop.

    A name-keyed CSV row that resolves to a real (node, site) but
    disagrees with the tree's parent_aa/child_aa is skipped (its site
    never enters the listed set). Combined with --mutations-allow-mismatch
    the run continues, and --only-listed-mutations then drops the bare
    derived mutation at that site as well — the rejected CSV claim is
    treated as "no claim," not as evidence the bare event should survive.
    """
    olmsted = {
        "metadata": {"format": "olmsted", "format_version": "1.0"},
        "datasets": [{"dataset_id": "ds", "name": "Integrity Cascade Test"}],
        "clones": {"ds": [{"clone_id": "fam1", "dataset_id": "ds",
                           "unique_seqs_count": 2, "mean_mut_freq": 0.0,
                           "sample_id": "s1"}]},
        "trees": [
            {
                "ident": "tree-1",
                "clone_id": "fam1",
                "newick": "(leaf:0.1)root;",
                "nodes": [
                    {"sequence_id": "root", "parent": None, "type": "root",
                     "sequence_alignment_aa": "MQQ"},
                    # Tree-derived mutations: site 1 Q→K, site 2 Q→R.
                    {"sequence_id": "leaf", "parent": "root", "type": "leaf",
                     "sequence_alignment_aa": "MKR"},
                ],
            }
        ],
    }
    json_path = tmp_path / "input.json"
    csv_path = tmp_path / "muts.csv"
    out_path = tmp_path / "out.json"
    json_path.write_text(json.dumps(olmsted))
    # Site 1 matches; site 2 lies about child_aa (claims S, tree has R)
    # → integrity mismatch → row skipped → site 2 not in listed set.
    csv_path.write_text(
        "family,node_name,site,parent_aa,child_aa,score\n"
        "fam1,leaf,1,Q,K,9.9\n"
        "fam1,leaf,2,Q,S,7.7\n"
    )

    result = subprocess.run(
        ["olmsted", "merge", "-i", str(json_path), "--mutations", str(csv_path),
         "--mutations-allow-mismatch", "--only-listed-mutations",
         "-o", str(out_path)],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    out = json.loads(out_path.read_text())
    leaf = next(n for n in out["trees"][0]["nodes"] if n["sequence_id"] == "leaf")
    # Only the integrity-clean site survives. The mismatched site's
    # bare derived mutation is dropped under --only-listed-mutations
    # even though --mutations-allow-mismatch keeps the run alive.
    assert len(leaf["mutations"]) == 1
    assert leaf["mutations"][0]["site"] == 1
    assert leaf["mutations"][0]["score"] == 9.9


@pytest.mark.parametrize("flag", ["--mutations-use-depth", "--mutations-allow-mismatch",
                                  "--only-listed-mutations"])
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
