"""Phase 2 streaming primitives (#26): evidence accumulators, spooler, writer.

These tests pin the behavior phase 3 will rely on:

- ``FieldTypeEvidence`` counters infer the same field types as the legacy
  ``field_metadata.infer_field_type`` for in-distribution data, and remain
  correct after merging contradicting batches (the failure mode the
  sample-capped path misses).
- ``RangeEvidence.merge`` preserves running min/max across batches.
- ``BatchAccumulator.finalize_field_metadata`` produces the same shape
  ``field_metadata.generate_field_metadata`` does on a representative
  dataset.
- ``BatchSpooler`` round-trips clones and trees written across multiple
  ``write_batch`` calls.
- ``write_olmsted_json_streaming`` emits canonical key order (``metadata``
  first) and produces a parsed structure equal to the all-in-one
  ``write_olmsted_json``.
"""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import pytest

from olmsted_cli.field_metadata import infer_field_type
from olmsted_cli.streaming import (
    BatchAccumulator,
    BatchSpooler,
    DuplicateIdError,
    FieldTypeEvidence,
    RangeEvidence,
    write_olmsted_json_streaming,
)

# =============================================================================
# FieldTypeEvidence
# =============================================================================


def _evidence_for(values):
    ev = FieldTypeEvidence()
    for v in values:
        ev.record(v)
    return ev


@pytest.mark.parametrize(
    "values,expected",
    [
        ([1, 2, 3], "continuous"),
        ([1.5, 2.5], "continuous"),
        (["J", "O", "Z"], "categorical"),  # letters outside both DNA and AA
        (["A", "C", "G", "T"], "dna"),
        (["E", "F", "I"], "aa"),  # AA-only (not in DNA's IUPAC alphabet)
        (["A", "C", "G", "T", "B"], "dna"),  # B is DNA IUPAC ambiguity
        ([True, False], "categorical"),
        ([[1, 2], [3]], "list"),
        ([{"a": 1}], "json"),
        ([1, "a"], "categorical"),
        ([[1], {"a": 1}], "categorical"),
        ([], "categorical"),
    ],
)
def test_field_type_evidence_matches_legacy_inference(values, expected):
    """Single-batch behavior must match ``infer_field_type`` for parity."""
    legacy = infer_field_type(list(values))
    assert legacy == expected, f"sanity: legacy gave {legacy!r}"
    assert _evidence_for(values).infer() == expected


def test_field_type_evidence_merge_int_then_string_flips_to_categorical():
    """Counter-based inference catches type drift that capped sampling misses."""
    a = _evidence_for([1, 2, 3, 4, 5] * 20)  # 100 ints
    b = _evidence_for(["surprise"])
    a.merge(b)
    assert a.infer() == "categorical"


def test_field_type_evidence_merge_dna_then_multi_char_flips_to_categorical():
    a = _evidence_for(["A", "C", "G", "T"])
    assert a.infer() == "dna"
    a.merge(_evidence_for(["AC"]))
    assert a.infer() == "categorical"


def test_field_type_evidence_merge_preserves_dna_when_all_single_char():
    a = _evidence_for(["A", "C"])
    b = _evidence_for(["G", "T", "U"])
    a.merge(b)
    assert a.infer() == "dna"


# =============================================================================
# RangeEvidence
# =============================================================================


def test_range_evidence_records_min_max():
    rng = RangeEvidence()
    for v in [3, 1, 2, 5, 4]:
        rng.record(v)
    assert rng.as_list() == [1, 5]


def test_range_evidence_ignores_booleans():
    rng = RangeEvidence()
    rng.record(True)
    rng.record(False)
    assert rng.as_list() is None


def test_range_evidence_merge():
    a = RangeEvidence()
    for v in [3.0, 5.0]:
        a.record(v)
    b = RangeEvidence()
    for v in [1.0, 10.0]:
        b.record(v)
    a.merge(b)
    assert a.as_list() == [1.0, 10.0]


def test_range_evidence_merge_into_empty():
    a = RangeEvidence()
    b = RangeEvidence()
    b.record(7)
    a.merge(b)
    assert a.as_list() == [7, 7]


# =============================================================================
# BatchAccumulator
# =============================================================================


def _toy_clone(clone_id, sample_id, extras=None, trees=None):
    clone = {
        "clone_id": clone_id,
        "sample_id": sample_id,
        "unique_seqs_count": 5,
        "mean_mut_freq": 0.1,
        "trees": trees or [],
    }
    if extras:
        clone.update(extras)
    return clone


def _toy_tree(tree_id, clone_id, nodes):
    return {
        "ident": f"tree-{tree_id}",
        "tree_id": tree_id,
        "clone_id": clone_id,
        "newick": "(A:1,B:1)naive;",
        "nodes": nodes,
    }


def _toy_node(seq_id, *, type_="leaf", parent="naive", extras=None, mutations=None):
    node = {
        "sequence_id": seq_id,
        "sequence_alignment": "ACGT",
        "sequence_alignment_aa": "MD",
        "type": type_,
        "parent": parent,
    }
    if extras:
        node.update(extras)
    if mutations:
        node["mutations"] = mutations
    return node


def test_accumulator_totals_and_uniqueness():
    acc = BatchAccumulator()
    acc.register_dataset("ds-1")

    clone_a = _toy_clone("c-a", "s1", trees=[{"tree_id": "t-a-1", "clone_id": "c-a"}])
    clone_b = _toy_clone("c-b", "s1", trees=[{"tree_id": "t-b-1", "clone_id": "c-b"}])
    tree_a = _toy_tree("t-a-1", "c-a", [_toy_node("L1"), _toy_node("L2")])
    tree_b = _toy_tree("t-b-1", "c-b", [_toy_node("L3")])

    acc.observe_batch("ds-1", [clone_a], [tree_a])
    acc.observe_batch("ds-1", [clone_b], [tree_b])

    totals = acc.finalize_totals()
    assert totals == {
        "datasets_count": 1,
        "total_clones_count": 2,
        "total_trees_count": 2,
        "total_leaf_nodes_count": 3,
    }


def test_accumulator_raises_on_duplicate_clone_id():
    acc = BatchAccumulator()
    acc.register_dataset("ds-1")
    acc.observe_batch("ds-1", [_toy_clone("c-a", "s1")], [])
    with pytest.raises(DuplicateIdError, match="duplicate clone_id"):
        acc.observe_batch("ds-1", [_toy_clone("c-a", "s1")], [])


def test_accumulator_allows_duplicate_ids_when_flagged():
    acc = BatchAccumulator(allow_duplicate_ids=True)
    acc.register_dataset("ds-1")
    acc.observe_batch("ds-1", [_toy_clone("c-a", "s1")], [])
    acc.observe_batch("ds-1", [_toy_clone("c-a", "s1")], [])
    assert any("clone_id" in msg for msg in acc.duplicate_warnings)


def test_accumulator_field_metadata_levels_and_ranges():
    acc = BatchAccumulator()
    acc.register_dataset("ds-1")

    nodes_a = [
        _toy_node(
            "L1",
            extras={"score": 1.0},
            mutations=[
                {"site": 1, "parent_aa": "A", "child_aa": "C", "surprise_mutsel": 2.0},
            ],
        ),
        _toy_node("L2", extras={"score": 3.0}),
    ]
    nodes_b = [
        _toy_node(
            "L3",
            extras={"score": 5.0},
            mutations=[
                {"site": 2, "parent_aa": "G", "child_aa": "D", "surprise_mutsel": -1.0},
            ],
        ),
    ]
    clone_a = _toy_clone(
        "c-a",
        "s1",
        trees=[{"tree_id": "t-a-1", "clone_id": "c-a"}],
    )
    clone_b = _toy_clone(
        "c-b",
        "s1",
        trees=[{"tree_id": "t-b-1", "clone_id": "c-b"}],
    )
    acc.observe_batch(
        "ds-1",
        [clone_a],
        [_toy_tree("t-a-1", "c-a", nodes_a)],
    )
    acc.observe_batch(
        "ds-1",
        [clone_b],
        [_toy_tree("t-b-1", "c-b", nodes_b)],
    )

    fm = acc.finalize_field_metadata("ds-1")

    # Clone-level numeric field has a range spanning both batches
    assert "clone" in fm
    assert fm["clone"]["unique_seqs_count"]["type"] == "continuous"

    # Node-level "score" is continuous but doesn't get a range key
    # (matches the legacy generate_node_metadata, which only emits range
    # at tree and mutation levels).
    assert fm["node"]["score"]["type"] == "continuous"
    assert "range" not in fm["node"]["score"]

    # Mutation-level surprise_mutsel range spans the two batches
    assert fm["mutation"]["surprise_mutsel"]["type"] == "continuous"
    assert fm["mutation"]["surprise_mutsel"]["range"] == [-1.0, 2.0]


def test_accumulator_derives_aa_fields_without_mutations():
    """No mutation arrays + AA sequences on nodes → derived child_aa/parent_aa."""
    acc = BatchAccumulator()
    acc.register_dataset("ds-1")
    clone = _toy_clone("c-a", "s1", trees=[{"tree_id": "t-a-1", "clone_id": "c-a"}])
    acc.observe_batch(
        "ds-1",
        [clone],
        [_toy_tree("t-a-1", "c-a", [_toy_node("L1"), _toy_node("L2")])],
    )

    fm = acc.finalize_field_metadata("ds-1")
    assert fm["mutation"]["child_aa"]["type"] == "aa"
    assert fm["mutation"]["parent_aa"]["display"] == "tooltip"


def test_accumulator_promotes_tree_level_field_when_varies_in_clone():
    acc = BatchAccumulator()
    acc.register_dataset("ds-1")

    clone = _toy_clone(
        "c-a",
        "s1",
        trees=[
            {"tree_id": "t-a-1", "clone_id": "c-a", "method_score": 0.1},
            {"tree_id": "t-a-2", "clone_id": "c-a", "method_score": 0.9},
        ],
    )
    acc.observe_batch(
        "ds-1",
        [clone],
        [
            _toy_tree("t-a-1", "c-a", [_toy_node("L1")]),
            _toy_tree("t-a-2", "c-a", [_toy_node("L1")]),
        ],
    )

    fm = acc.finalize_field_metadata("ds-1")
    assert "tree" in fm
    assert fm["tree"]["method_score"]["type"] == "continuous"
    assert "method_score" not in fm.get("clone", {})


def test_accumulator_samples_dedup_by_sample_id():
    acc = BatchAccumulator()
    acc.register_dataset("ds-1")
    acc.add_samples("ds-1", [{"sample_id": "s1", "locus": "IGH"}])
    acc.add_samples("ds-1", [{"sample_id": "s1", "locus": "IGH"}, {"sample_id": "s2"}])
    samples = acc.samples_for("ds-1")
    assert [s["sample_id"] for s in samples] == ["s1", "s2"]


# =============================================================================
# BatchSpooler
# =============================================================================


def test_batch_spooler_round_trip():
    with BatchSpooler() as spooler:
        spooler.write_batch(
            "ds-1",
            [{"clone_id": "c-a"}, {"clone_id": "c-b"}],
            [{"tree_id": "t-1"}],
        )
        spooler.write_batch(
            "ds-1",
            [{"clone_id": "c-c"}],
            [{"tree_id": "t-2"}, {"tree_id": "t-3"}],
        )
        clones = list(spooler.iter_clones("ds-1"))
        trees = list(spooler.iter_trees("ds-1"))

    assert [c["clone_id"] for c in clones] == ["c-a", "c-b", "c-c"]
    assert [t["tree_id"] for t in trees] == ["t-1", "t-2", "t-3"]


def test_batch_spooler_keeps_datasets_separate():
    with BatchSpooler() as spooler:
        spooler.write_batch("ds-1", [{"clone_id": "x"}], [])
        spooler.write_batch("ds-2", [{"clone_id": "y"}], [])
        assert {ds for ds in spooler.dataset_ids()} == {"ds-1", "ds-2"}
        assert list(spooler.iter_clones("ds-1")) == [{"clone_id": "x"}]
        assert list(spooler.iter_clones("ds-2")) == [{"clone_id": "y"}]


def test_batch_spooler_requires_context_manager():
    with pytest.raises(RuntimeError, match="context manager"):
        BatchSpooler().write_batch("ds-1", [], [])


# =============================================================================
# write_olmsted_json_streaming
# =============================================================================


def _build_streamed(tmp_path: Path, json_format: str = "pretty") -> Path:
    metadata = {
        "format": "olmsted",
        "format_version": "1.0",
        "processing_info": {
            "datasets_count": 1,
            "total_clones_count": 2,
            "total_trees_count": 1,
            "total_leaf_nodes_count": 2,
        },
    }
    datasets = [{"dataset_id": "ds-1", "field_metadata": {}}]
    out = tmp_path / "stream.json"
    with BatchSpooler() as spooler:
        spooler.write_batch(
            "ds-1",
            [{"clone_id": "c-a"}, {"clone_id": "c-b"}],
            [{"tree_id": "t-1", "newick": "(A:1);"}],
        )
        written = write_olmsted_json_streaming(
            metadata, datasets, spooler, str(out), json_format=json_format
        )
    return Path(written)


def test_streaming_writer_emits_metadata_first_pretty(tmp_path):
    path = _build_streamed(tmp_path, "pretty")
    text = path.read_text()

    # object_pairs_hook surfaces the actual on-disk key order
    parsed_pairs = json.loads(text, object_pairs_hook=list)
    top_keys = [k for k, _ in parsed_pairs]
    assert top_keys == ["metadata", "datasets", "clones", "trees"]


def test_streaming_writer_compact_format_parses_equivalently(tmp_path):
    path = _build_streamed(tmp_path, "compact")
    parsed = json.loads(path.read_text())

    assert parsed["metadata"]["format"] == "olmsted"
    assert [c["clone_id"] for c in parsed["clones"]["ds-1"]] == ["c-a", "c-b"]
    assert [t["tree_id"] for t in parsed["trees"]] == ["t-1"]


def test_streaming_writer_gzip_is_deterministic(tmp_path):
    """Gzip header pinning matches data_io.write_olmsted_json's guarantee."""
    metadata = {"format": "olmsted"}
    datasets = [{"dataset_id": "ds-1", "field_metadata": {}}]

    def _run(out_name: str) -> bytes:
        out = tmp_path / out_name
        with BatchSpooler() as spooler:
            spooler.write_batch("ds-1", [{"clone_id": "c-a"}], [])
            write_olmsted_json_streaming(
                metadata, datasets, spooler, str(out), json_format="gzip"
            )
        return (tmp_path / (out_name + ".gz")).read_bytes()

    first = _run("a.json")
    second = _run("b.json")
    # Header bytes (first 10) include the timestamp / filename fields we pin.
    assert first[:10] == second[:10]

    payload_a = gzip.decompress(first).decode("utf-8")
    parsed = json.loads(payload_a)
    assert parsed["metadata"]["format"] == "olmsted"
    assert parsed["clones"]["ds-1"][0]["clone_id"] == "c-a"


def test_streaming_writer_handles_empty_trees(tmp_path):
    metadata = {"format": "olmsted"}
    datasets = [{"dataset_id": "ds-1", "field_metadata": {}}]
    out = tmp_path / "empty.json"
    with BatchSpooler() as spooler:
        spooler.write_batch("ds-1", [{"clone_id": "only"}], [])
        write_olmsted_json_streaming(
            metadata, datasets, spooler, str(out), json_format="compact"
        )
    parsed = json.loads(out.read_text())
    assert parsed["trees"] == []
    assert parsed["clones"]["ds-1"] == [{"clone_id": "only"}]
