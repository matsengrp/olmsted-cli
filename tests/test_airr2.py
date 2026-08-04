#!/usr/bin/env python3
"""Unit tests for the AIRR-C v2 Clone/Tree ("airr2") processor.

Covers the ``process_airr2_data`` module: Rearrangement indexing, Newick-driven
topology with reroot-on-germline, node → sequence joins (incl. inferred
germline), ``node_type`` passthrough, the paired heavy/light split, dataset
synthesis, and graceful handling of a missing sequence.
"""

import json
from pathlib import Path

import pytest

from olmsted_cli.identifier import IdentMinter
from olmsted_cli.process_airr2_data import (
    _locus_chain,
    _mean_mutation_frequency,
    index_rearrangements,
    process_airr2_to_olmsted,
)

EXAMPLE_DIR = Path(__file__).parent.parent / "example-data" / "airr2"


def _load(variant):
    with open(EXAMPLE_DIR / f"input-{variant}.json") as f:
        return json.load(f)


def _run(variant, **kwargs):
    data = _load(variant)
    return process_airr2_to_olmsted(
        data["Clone"],
        data["Rearrangement"],
        minter=IdentMinter(seed=42),
        verbosity=0,
        **kwargs,
    )


def _clones(clones_dict):
    return [c for cl in clones_dict.values() for c in cl]


@pytest.mark.airr2
class TestLocusChain:
    def test_heavy(self):
        assert _locus_chain("IGH") == "heavy"

    def test_light_kappa_lambda(self):
        assert _locus_chain("IGK") == "light"
        assert _locus_chain("IGL") == "light"

    def test_unknown_and_none(self):
        assert _locus_chain("TRB") is None
        assert _locus_chain(None) is None


@pytest.mark.airr2
class TestIndexRearrangements:
    def test_by_sequence_and_cell(self):
        data = _load("paired")
        by_seq, by_cell = index_rearrangements(data["Rearrangement"])
        # Rearrangement class join key
        assert all(isinstance(v, dict) for v in by_seq.values())
        # Paired cells carry one record per locus.
        counts = {len(v) for v in by_cell.values()}
        assert counts == {2}


@pytest.mark.airr2
class TestNocell:
    def test_clone_and_tree_counts(self):
        _datasets, clones_dict, trees = _run("nocell")
        clones = _clones(clones_dict)
        assert len(clones) == 2
        assert len(trees) == 2
        assert {c["clone_class"] for c in clones} == {"Rearrangement"}

    def test_topology_root_is_germline_inferred(self):
        _datasets, _clones_dict, trees = _run("nocell")
        for tree in trees:
            roots = [n for n in tree["nodes"] if n["type"] == "root"]
            assert len(roots) == 1
            root = roots[0]
            # Germline is the rerooted root and is inferred (ASR), with sequence.
            assert root["sequence_id"].startswith("Germline-")
            assert root["node_type"] == "inferred"
            assert root["sequence_alignment"]
            assert root["parent"] is None

    def test_every_node_has_sequence_and_aa(self):
        _datasets, _clones_dict, trees = _run("nocell")
        for tree in trees:
            for node in tree["nodes"]:
                assert node["sequence_alignment"], node["sequence_id"]
                assert node["sequence_alignment_aa"]
                assert node["node_type"] in ("observed", "inferred")

    def test_single_root_and_parents_resolve(self):
        _datasets, _clones_dict, trees = _run("nocell")
        for tree in trees:
            ids = {n["sequence_id"] for n in tree["nodes"]}
            for node in tree["nodes"]:
                if node["type"] != "root":
                    assert node["parent"] in ids


@pytest.mark.airr2
class TestUnpairedCell:
    def test_cell_class_single_chain(self):
        _datasets, clones_dict, trees = _run("unpaired")
        clones = _clones(clones_dict)
        assert len(clones) == 2
        assert {c["clone_class"] for c in clones} == {"Cell"}
        # Single-chain cell → no -heavy/-light suffix.
        assert all("-heavy" not in c["clone_id"] for c in clones)


@pytest.mark.airr2
@pytest.mark.paired
class TestPairedCell:
    def test_splits_into_heavy_and_light(self):
        _datasets, clones_dict, trees = _run("paired")
        clones = _clones(clones_dict)
        # Two input clones → four output clones (heavy + light each).
        assert len(clones) == 4
        assert len(trees) == 4
        ids = {c["clone_id"] for c in clones}
        assert ids == {"10004-heavy", "10004-light", "8232-heavy", "8232-light"}

    def test_pair_id_links_chains(self):
        _datasets, clones_dict, _trees = _run("paired")
        by_id = {c["clone_id"]: c for c in _clones(clones_dict)}
        assert by_id["10004-heavy"]["pair_id"] == by_id["10004-light"]["pair_id"]
        assert by_id["10004-heavy"]["is_paired"] is True

    def test_per_locus_sequences_differ(self):
        """Heavy and light trees share topology but carry different sequences."""
        _datasets, _clones_dict, trees = _run("paired")
        by_ident = {t["ident"]: t for t in trees}
        heavy = next(t for t in trees if t["ident"].endswith("-heavy"))
        light = next(
            t
            for t in trees
            if t["ident"].endswith("-light")
            and t["clone_id"].split("-")[0] == heavy["clone_id"].split("-")[0]
        )
        assert by_ident  # sanity
        # Same node ids (shared topology)...
        assert {n["sequence_id"] for n in heavy["nodes"]} == {
            n["sequence_id"] for n in light["nodes"]
        }
        # ...but the germline sequences differ (different loci).
        h_root = next(n for n in heavy["nodes"] if n["type"] == "root")
        l_root = next(n for n in light["nodes"] if n["type"] == "root")
        assert h_root["sequence_alignment"] != l_root["sequence_alignment"]
        assert h_root["locus"] == "IGH"
        assert l_root["locus"] in ("IGK", "IGL")

    def test_sample_loci(self):
        datasets, clones_dict, _trees = _run("paired")
        loci = sorted(c["sample"]["locus"] for c in _clones(clones_dict))
        assert loci == ["igh", "igh", "igk", "igl"]
        # Dataset samples dedup by sample_id (all share repertoire_id "sample").
        assert len(datasets[0]["samples"]) == 1


@pytest.mark.airr2
class TestDatasetSynthesis:
    def test_dataset_shape(self):
        datasets, clones_dict, _trees = _run("nocell", name="my-airr2")
        assert len(datasets) == 1
        ds = datasets[0]
        assert ds["name"] == "my-airr2"
        assert ds["dataset_id"] in clones_dict
        assert ds["clone_count"] == 2
        assert ds["subjects"] == []
        assert "field_metadata" in ds

    def test_required_clone_fields_present(self):
        _datasets, clones_dict, _trees = _run("nocell")
        for clone in _clones(clones_dict):
            assert clone["unique_seqs_count"] > 0
            assert "mean_mut_freq" in clone

    def test_node_type_label_does_not_collide_with_topology_type(self):
        """node_type must not share the "Node Type" label the webapp uses for
        the built-in topological ``type`` field.

        Two node fields with the same tooltip label produce a duplicate object
        key -> invalid Vega expression -> the webapp fails to render the tree.
        """
        datasets, _clones_dict, _trees = _run("paired")
        node_meta = datasets[0]["field_metadata"]["node"]
        assert node_meta["node_type"]["label"] != "Node Type"
        # No two node fields share a label, and none collide with the webapp's
        # built-in node-field labels.
        builtin_labels = {
            "Sequence ID",
            "Parent ID",
            "Node Type",
            "Distance",
            "Depth",
        }
        labels = [m["label"] for m in node_meta.values()]
        assert len(labels) == len(set(labels)), "duplicate field_metadata labels"
        assert not (set(labels) & builtin_labels), "collides with a built-in label"


@pytest.mark.airr2
class TestMetrics:
    def test_compute_metrics_populates_lbi(self):
        _datasets, _clones_dict, trees = _run("nocell", compute_metrics=True)
        # At least some nodes get a numeric LBI when metrics are requested.
        lbis = [n["lbi"] for t in trees for n in t["nodes"]]
        assert any(v is not None for v in lbis)


@pytest.mark.airr2
class TestMeanMutationFrequency:
    def test_zero_without_germline(self):
        assert _mean_mutation_frequency([], None) == 0.0

    def test_counts_only_leaves(self):
        germline = "AAAA"
        nodes = [
            {"type": "root", "sequence_alignment": "AAAA"},
            {"type": "leaf", "sequence_alignment": "AAAT"},  # 1/4
            {"type": "leaf", "sequence_alignment": "AATT"},  # 2/4
        ]
        # (0.25 + 0.5) / 2 leaves
        assert _mean_mutation_frequency(nodes, germline) == pytest.approx(0.375)

    def test_skips_gaps_and_n(self):
        germline = "AAAA"
        nodes = [{"type": "leaf", "sequence_alignment": "A.NT"}]
        # Only positions 0 and 3 comparable; 1 mismatch of 2 → 0.5
        assert _mean_mutation_frequency(nodes, germline) == pytest.approx(0.5)


@pytest.mark.airr2
class TestMissingSequenceGraceful:
    def test_missing_rearrangement_yields_empty_sequence(self):
        """A node with no matching Rearrangement gets an empty sequence, not a crash."""
        data = _load("nocell")
        # Drop one observed leaf's Rearrangement record.
        clone = data["Clone"][0]
        leaf = next(n for n in clone["nodes"] if n["node_type"] == "observed")
        dropped_id = leaf["sequence_id"]
        rearrangements = [
            r for r in data["Rearrangement"] if r.get("sequence_id") != dropped_id
        ]
        _datasets, _clones_dict, trees = process_airr2_to_olmsted(
            [clone], rearrangements, minter=IdentMinter(seed=42), verbosity=0
        )
        node = next(
            n for t in trees for n in t["nodes"] if n["sequence_id"] == dropped_id
        )
        assert node["sequence_alignment"] == ""
        assert node["sequence_alignment_aa"] == ""
