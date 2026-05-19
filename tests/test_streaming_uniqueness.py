"""Streaming-path id-uniqueness + empty-dataset output shape (#26 review).

Regression coverage for two issues clean-code review flagged:

- C1: the streaming path skipped ``check_output_id_uniqueness``, letting
  duplicate ``sample_id`` / ``subject_id`` slip through.
- C2: a dataset with zero clones registered but no batches spooled
  produced ``"clones": {}`` instead of ``"clones": {"<id>": []}`` —
  diverging from the legacy output shape.
"""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pytest


def _run(cmd):
    return subprocess.run(cmd, capture_output=True, text=True)


def _write_airr(tmp_path: Path, dataset: dict) -> Path:
    path = tmp_path / "input.json"
    path.write_text(json.dumps(dataset))
    return path


def _airr_dataset_with_duplicate_samples():
    return {
        "dataset_id": "ds-dup-sample",
        "subjects": [{"subject_id": "subj1"}],
        "samples": [
            {"sample_id": "s1", "locus": "IGH"},
            {"sample_id": "s1", "locus": "IGK"},
        ],
        "seeds": [],
        "clones": [
            {
                "clone_id": "c1",
                "sample_id": "s1",
                "subject_id": "subj1",
                "unique_seqs_count": 1,
                "mean_mut_freq": 0.0,
                "trees": [
                    {
                        "newick": "(a:1);",
                        "nodes": {
                            "a": {
                                "sequence_id": "a",
                                "sequence_alignment": "ACGT",
                                "sequence_alignment_aa": "T",
                                "type": "leaf",
                            }
                        },
                    }
                ],
            }
        ],
    }


def _airr_dataset_with_duplicate_subjects():
    d = _airr_dataset_with_duplicate_samples()
    d["dataset_id"] = "ds-dup-subj"
    d["samples"] = [{"sample_id": "s1", "locus": "IGH"}]
    d["subjects"] = [
        {"subject_id": "subj1"},
        {"subject_id": "subj1"},
    ]
    return d


def _airr_dataset_empty_clones():
    return {
        "dataset_id": "ds-empty",
        "subjects": [{"subject_id": "subj1"}],
        "samples": [{"sample_id": "s1", "locus": "IGH"}],
        "seeds": [],
        "clones": [],
    }


@pytest.mark.parametrize(
    "dataset_factory,expected_scope",
    [
        (_airr_dataset_with_duplicate_samples, "sample_id"),
        (_airr_dataset_with_duplicate_subjects, "subject_id"),
    ],
)
def test_airr_streaming_rejects_duplicate_ids(
    dataset_factory, expected_scope, tmp_path
):
    """The streaming path must raise on ``sample_id`` / ``subject_id`` dups.

    Legacy ``process_airr_format`` calls ``check_output_id_uniqueness``
    just before write; the streaming refactor needs the same guarantee
    so the webapp's Redux-keyed-on-id store doesn't silently overwrite.
    """
    inp = _write_airr(tmp_path, dataset_factory())
    out = tmp_path / "out.json"
    result = _run([
        "olmsted", "process",
        "-f", "airr",
        "-i", str(inp),
        "-o", str(out),
        "--seed", "42",
        "--batch-size", "1",
        "-q",
    ])
    assert result.returncode != 0, (
        f"streaming path silently accepted duplicate {expected_scope}: "
        f"{result.stdout}\n{result.stderr}"
    )
    assert expected_scope in (result.stdout + result.stderr)


@pytest.mark.parametrize(
    "dataset_factory",
    [
        _airr_dataset_with_duplicate_samples,
        _airr_dataset_with_duplicate_subjects,
    ],
)
def test_airr_streaming_allow_duplicate_ids_downgrades(dataset_factory, tmp_path):
    """``--allow-duplicate-ids`` downgrades the failure to a warning."""
    inp = _write_airr(tmp_path, dataset_factory())
    out = tmp_path / "out.json"
    result = _run([
        "olmsted", "process",
        "-f", "airr",
        "-i", str(inp),
        "-o", str(out),
        "--seed", "42",
        "--batch-size", "1",
        "--allow-duplicate-ids",
    ])
    assert result.returncode == 0, (
        f"--allow-duplicate-ids should not exit non-zero: "
        f"{result.stdout}\n{result.stderr}"
    )
    assert out.exists()


def test_airr_streaming_empty_clones_emits_dataset_key(tmp_path):
    """Empty-clones dataset must still appear in the clones map.

    The legacy path always assigns ``clones_dict[dataset_id] = []`` so
    every dataset has a corresponding key; the webapp iterates
    ``datasets[]`` and looks up ``clones[dataset_id]`` — a missing key
    would crash.
    """
    inp = _write_airr(tmp_path, _airr_dataset_empty_clones())
    out = tmp_path / "out.json"
    result = _run([
        "olmsted", "process",
        "-f", "airr",
        "-i", str(inp),
        "-o", str(out),
        "--seed", "42",
        "--batch-size", "1",
        "-q",
    ])
    assert result.returncode == 0, (
        f"empty-clones AIRR input failed: {result.stdout}\n{result.stderr}"
    )

    data = json.loads(out.read_text())
    assert "clones" in data
    ds_ids = [d["dataset_id"] for d in data["datasets"]]
    for ds_id in ds_ids:
        assert ds_id in data["clones"], (
            f"streaming output missing clones[{ds_id!r}] for empty-clones dataset"
        )
        assert data["clones"][ds_id] == []
