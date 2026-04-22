"""Integration test for `olmsted merge` against the real-data fixture in example_data/merge/.

The fixture is a 2-clone subset of a DASM2 surprise analysis (top20_olmsted.json) plus
the matching subset of its mutation scores CSV. It exercises both the happy path
(mutations successfully merged into derived parent/child diffs) and the warning path
(CSV rows that reference sites with no corresponding derived mutation).
"""

import json
import subprocess
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).parent.parent
FIXTURE_DIR = REPO_ROOT / "example_data" / "merge"
FIXTURE_JSON = FIXTURE_DIR / "olmsted_input.json"
FIXTURE_CSV = FIXTURE_DIR / "mutations.csv"


@pytest.fixture(scope="module")
def fixture_files_exist():
    assert FIXTURE_JSON.exists(), f"Missing fixture: {FIXTURE_JSON}"
    assert FIXTURE_CSV.exists(), f"Missing fixture: {FIXTURE_CSV}"


def test_merge_fixture_end_to_end(fixture_files_exist, tmp_path):
    """Run `olmsted merge` against the real-data fixture and verify expected counts."""
    out_path = tmp_path / "merged.json"

    result = subprocess.run(
        [
            "olmsted",
            "merge",
            "-i",
            str(FIXTURE_JSON),
            "--mutations",
            str(FIXTURE_CSV),
            "--mutations-use-depth",
            "-o",
            str(out_path),
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, f"merge failed: {result.stderr}"
    assert out_path.exists()

    combined = result.stdout + result.stderr

    # CSV load: 37 rows across 2 families
    assert "Loaded 37 CSV rows across 2 families" in combined
    # The fixture's CSV has a `depth` column → disambiguation is active
    assert "Disambiguation columns in CSV: depth" in combined
    # With depth disambiguation (naive-rooted): 33 enrichments across 18 nodes
    assert "Enriched 33 mutations across 18 nodes in 2 trees" in combined
    # 19 of 37 rows have no node-mutation match at the matching depth
    assert "Unmatched: 19/37 CSV rows" in combined
    # 9 CSV rows broadcast to multiple nodes (convergent mutations at same depth)
    assert "Broadcast: 9 CSV rows matched multiple nodes" in combined
    # And the corresponding warnings
    assert "had no corresponding derived mutation" in combined
    assert "broadcast to multiple node-mutations" in combined


def test_merge_fixture_output_structure(fixture_files_exist, tmp_path):
    """Verify the merged output has the expected enriched fields and metadata."""
    out_path = tmp_path / "merged.json"

    subprocess.run(
        [
            "olmsted",
            "merge",
            "-i",
            str(FIXTURE_JSON),
            "--mutations",
            str(FIXTURE_CSV),
            "--mutations-use-depth",
            "-o",
            str(out_path),
            "-q",
        ],
        check=True,
        capture_output=True,
    )

    out = json.loads(out_path.read_text())

    # Both clones preserved
    assert len(out["trees"]) == 2
    clone_ids = {t["clone_id"] for t in out["trees"]}
    assert clone_ids == {"192491-igk-192491", "203694-igk-203694"}

    # field_metadata should declare the merged mutation fields
    fm_mut = out["datasets"][0]["field_metadata"]["mutation"]
    for field in (
        "surprise_mutsel",
        "surprise_neutral",
        "selection_contribution",
        "log_selection_factor",
        "num_codon_changes",
        "surprise_mutsel_theoretical",
    ):
        assert field in fm_mut, f"Expected merged mutation field {field!r} in field_metadata"

    # Spot-check that at least one node has a mutation entry with the merged fields
    enriched_count = 0
    for tree in out["trees"]:
        for node in tree["nodes"]:
            for mut in node.get("mutations", []) or []:
                if "surprise_mutsel" in mut:
                    enriched_count += 1
                    # Sanity-check shape of enriched record
                    assert "site" in mut
                    assert "parent_aa" in mut
                    assert "child_aa" in mut
                    assert isinstance(mut["surprise_mutsel"], (int, float))
    assert enriched_count == 33, (
        f"Expected 33 enriched mutation records, found {enriched_count}"
    )


def test_merge_fixture_key_columns_excluded(fixture_files_exist, tmp_path):
    """Key/structural CSV columns should not leak into the merged mutation records."""
    out_path = tmp_path / "merged.json"

    subprocess.run(
        [
            "olmsted",
            "merge",
            "-i",
            str(FIXTURE_JSON),
            "--mutations",
            str(FIXTURE_CSV),
            "--mutations-use-depth",
            "-o",
            str(out_path),
            "-q",
        ],
        check=True,
        capture_output=True,
    )

    out = json.loads(out_path.read_text())
    excluded = {"family", "sample_id", "pcp_index", "depth"}
    for tree in out["trees"]:
        for node in tree["nodes"]:
            for mut in node.get("mutations", []) or []:
                if "surprise_mutsel" not in mut:
                    continue  # not an enriched mutation
                leaked = excluded & set(mut.keys())
                assert not leaked, f"Key columns leaked into merged mutation: {leaked}"
