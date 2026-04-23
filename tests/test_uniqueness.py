"""Tests for check_output_id_uniqueness and the --allow-duplicate-ids flag.

Exercises every scope the checker enforces:
- dataset_id across datasets[]
- clone_id within a dataset
- tree_id within a clone
- sample_id within dataset.samples[]
- subject_id within dataset.subjects[]
"""

import json
import subprocess
from pathlib import Path

import pytest

from olmsted_cli.process_utils import check_output_id_uniqueness


REPO_ROOT = Path(__file__).resolve().parents[1]


def _minimal_datasets(dataset_ids, *, samples=None, subjects=None):
    return [
        {"dataset_id": d, "samples": list(samples or []), "subjects": list(subjects or [])}
        for d in dataset_ids
    ]


class TestChecker:
    def test_passes_when_all_unique(self):
        datasets = _minimal_datasets(["ds-1"])
        clones = {"ds-1": [{"clone_id": "c1", "trees": []}]}
        # Should not raise
        check_output_id_uniqueness(datasets, clones, [])

    def test_detects_duplicate_dataset_id(self):
        datasets = _minimal_datasets(["ds-1", "ds-1"])
        with pytest.raises(ValueError, match="dataset_id"):
            check_output_id_uniqueness(datasets, {"ds-1": []}, [])

    def test_detects_duplicate_clone_id(self):
        datasets = _minimal_datasets(["ds-1"])
        clones = {
            "ds-1": [
                {"clone_id": "c1", "trees": []},
                {"clone_id": "c1", "trees": []},
            ]
        }
        with pytest.raises(ValueError, match="clone_id.*c1"):
            check_output_id_uniqueness(datasets, clones, [])

    def test_detects_duplicate_tree_id_within_clone(self):
        datasets = _minimal_datasets(["ds-1"])
        clones = {
            "ds-1": [
                {
                    "clone_id": "c1",
                    "trees": [{"tree_id": "t-a"}, {"tree_id": "t-a"}],
                }
            ]
        }
        with pytest.raises(ValueError, match="tree_id.*t-a"):
            check_output_id_uniqueness(datasets, clones, [])

    def test_tree_ids_may_repeat_across_clones(self):
        """The scope is within-a-clone, not across clones."""
        datasets = _minimal_datasets(["ds-1"])
        clones = {
            "ds-1": [
                {"clone_id": "c1", "trees": [{"tree_id": "t-a"}]},
                {"clone_id": "c2", "trees": [{"tree_id": "t-a"}]},
            ]
        }
        check_output_id_uniqueness(datasets, clones, [])  # no raise

    def test_detects_duplicate_sample_id(self):
        samples = [{"sample_id": "s1"}, {"sample_id": "s1"}]
        datasets = _minimal_datasets(["ds-1"], samples=samples)
        with pytest.raises(ValueError, match="sample_id"):
            check_output_id_uniqueness(datasets, {"ds-1": []}, [])

    def test_detects_duplicate_subject_id(self):
        subjects = [{"subject_id": "subj1"}, {"subject_id": "subj1"}]
        datasets = _minimal_datasets(["ds-1"], subjects=subjects)
        with pytest.raises(ValueError, match="subject_id"):
            check_output_id_uniqueness(datasets, {"ds-1": []}, [])

    def test_reports_all_violations_in_one_error(self):
        datasets = _minimal_datasets(["ds-1", "ds-1"])  # dataset_id dup
        clones = {
            "ds-1": [
                {"clone_id": "c1", "trees": []},
                {"clone_id": "c1", "trees": []},  # clone_id dup
            ]
        }
        with pytest.raises(ValueError) as excinfo:
            check_output_id_uniqueness(datasets, clones, [])
        msg = str(excinfo.value)
        # Both violations show up in the same error
        assert "dataset_id" in msg
        assert "clone_id" in msg

    def test_allow_duplicates_does_not_raise(self):
        datasets = _minimal_datasets(["ds-1", "ds-1"])
        # Should not raise when allow_duplicates=True
        check_output_id_uniqueness(datasets, {"ds-1": []}, [], allow_duplicates=True)

    def test_empty_or_absent_ids_not_counted(self):
        """Missing/empty *_id values don't collide with each other."""
        datasets = _minimal_datasets(["ds-1"])
        clones = {
            "ds-1": [
                {"clone_id": "", "trees": []},
                {"clone_id": None, "trees": []},
                {"trees": []},  # no clone_id key at all
            ]
        }
        # None of these should register as duplicates
        check_output_id_uniqueness(datasets, clones, [])


class TestCLIIntegration:
    """Smoke-test the --allow-duplicate-ids CLI flag via the `tag` command.

    `tag` is the simplest vehicle: we can inject duplicates into an
    Olmsted JSON file by editing a pre-processed golden input.
    """

    def _inject_duplicate_clone_id(self, path: Path) -> str:
        """Corrupt the golden by making two clones share the same clone_id.
        Returns the duplicated clone_id value."""
        data = json.loads(path.read_text())
        ds_id = data["datasets"][0]["dataset_id"]
        clones = data["clones"][ds_id]
        assert len(clones) >= 2, "need at least 2 clones to create a duplicate"
        dup_id = clones[0]["clone_id"]
        clones[1]["clone_id"] = dup_id
        path.write_text(json.dumps(data))
        return dup_id

    def test_tag_fails_on_duplicate_ids(self, tmp_path):
        """Default behavior: tag exits non-zero with a clear error."""
        # Produce a fresh Olmsted JSON to mutate
        source = tmp_path / "source.json"
        subprocess.run(
            [
                "olmsted", "process",
                "-f", "pcp",
                "-i", str(REPO_ROOT / "example_data/pcp/pcp.csv"),
                "-t", str(REPO_ROOT / "example_data/pcp/trees.csv"),
                "-o", str(source),
                "--seed", "42",
                "-q",
            ],
            check=True,
        )
        dup_id = self._inject_duplicate_clone_id(source)

        result = subprocess.run(
            ["olmsted", "tag", "-i", str(source), "-o", str(tmp_path / "out.json")],
            capture_output=True, text=True,
        )
        combined = result.stdout + result.stderr
        assert result.returncode != 0, f"expected failure, got: {combined}"
        assert "clone_id" in combined
        assert dup_id in combined

    def test_tag_passes_with_allow_duplicate_ids(self, tmp_path):
        """With the flag, duplicates warn but don't fail; output is
        written unchanged (duplicates preserved)."""
        source = tmp_path / "source.json"
        subprocess.run(
            [
                "olmsted", "process",
                "-f", "pcp",
                "-i", str(REPO_ROOT / "example_data/pcp/pcp.csv"),
                "-t", str(REPO_ROOT / "example_data/pcp/trees.csv"),
                "-o", str(source),
                "--seed", "42",
                "-q",
            ],
            check=True,
        )
        dup_id = self._inject_duplicate_clone_id(source)

        out = tmp_path / "out.json"
        result = subprocess.run(
            [
                "olmsted", "tag",
                "-i", str(source),
                "-o", str(out),
                "--allow-duplicate-ids",
            ],
            capture_output=True, text=True,
        )
        combined = result.stdout + result.stderr
        assert result.returncode == 0, f"unexpected failure: {combined}"
        # Warning reached one of the streams
        assert "clone_id" in combined
        # Output exists and still carries the duplicate — the flag
        # doesn't mutate data
        data = json.loads(out.read_text())
        ds_id = data["datasets"][0]["dataset_id"]
        clone_ids = [c["clone_id"] for c in data["clones"][ds_id]]
        assert clone_ids.count(dup_id) == 2
