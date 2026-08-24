#!/usr/bin/env python3
"""Unit tests for the AIRR-C v2 Clone/Tree processor.

Covers the ``process_airr_data`` module: Rearrangement indexing, Newick-driven
topology with reroot-on-germline, node → sequence joins (incl. inferred
germline), ``node_type`` passthrough, the paired heavy/light split, dataset
synthesis, and graceful handling of a missing sequence.
"""

import json
from pathlib import Path

import pytest

from olmsted_cli.identifier import IdentMinter
from olmsted_cli.process_airr_data import (
    _cdr_boundaries_from_region,
    _locus_chain,
    _mean_mutation_frequency,
    _node_multiplicity,
    index_rearrangements,
    process_airr_to_olmsted,
)

EXAMPLE_DIR = Path(__file__).parent.parent / "example-data" / "airr"


def _load(variant):
    with open(EXAMPLE_DIR / f"input-{variant}.json") as f:
        return json.load(f)


def _run(variant, **kwargs):
    data = _load(variant)
    return process_airr_to_olmsted(
        data["Clone"],
        data["Rearrangement"],
        minter=IdentMinter(seed=42),
        verbosity=0,
        **kwargs,
    )


def _clones(clones_dict):
    return [c for cl in clones_dict.values() for c in cl]


@pytest.mark.airr
class TestLocusChain:
    def test_heavy(self):
        assert _locus_chain("IGH") == "heavy"

    def test_light_kappa_lambda(self):
        assert _locus_chain("IGK") == "light"
        assert _locus_chain("IGL") == "light"

    def test_unknown_and_none(self):
        assert _locus_chain("TRB") is None
        assert _locus_chain(None) is None


@pytest.mark.airr
class TestIndexRearrangements:
    def test_by_sequence_and_cell(self):
        data = _load("paired-noinfo")
        by_seq, by_cell = index_rearrangements(data["Rearrangement"])
        # Rearrangement class join key
        assert all(isinstance(v, dict) for v in by_seq.values())
        # Paired cells carry one record per locus.
        counts = {len(v) for v in by_cell.values()}
        assert counts == {2}


@pytest.mark.airr
class TestNocell:
    def test_clone_and_tree_counts(self):
        _datasets, clones_dict, trees = _run("nocell-noinfo")
        clones = _clones(clones_dict)
        assert len(clones) == 2
        assert len(trees) == 2
        assert {c["clone_class"] for c in clones} == {"Rearrangement"}

    def test_topology_root_is_germline_inferred(self):
        _datasets, _clones_dict, trees = _run("nocell-noinfo")
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
        _datasets, _clones_dict, trees = _run("nocell-noinfo")
        for tree in trees:
            for node in tree["nodes"]:
                assert node["sequence_alignment"], node["sequence_id"]
                assert node["sequence_alignment_aa"]
                assert node["node_type"] in ("observed", "inferred")

    def test_single_root_and_parents_resolve(self):
        _datasets, _clones_dict, trees = _run("nocell-noinfo")
        for tree in trees:
            ids = {n["sequence_id"] for n in tree["nodes"]}
            for node in tree["nodes"]:
                if node["type"] != "root":
                    assert node["parent"] in ids


@pytest.mark.airr
class TestUnpairedCell:
    def test_cell_class_single_chain(self):
        _datasets, clones_dict, trees = _run("unpaired-noinfo")
        clones = _clones(clones_dict)
        assert len(clones) == 2
        assert {c["clone_class"] for c in clones} == {"Cell"}
        # Single-chain cell → no -heavy/-light suffix.
        assert all("-heavy" not in c["clone_id"] for c in clones)


@pytest.mark.airr
@pytest.mark.paired
class TestPairedCell:
    def test_splits_into_heavy_and_light(self):
        _datasets, clones_dict, trees = _run("paired-noinfo")
        clones = _clones(clones_dict)
        # Two input clones → four output clones (heavy + light each).
        assert len(clones) == 4
        assert len(trees) == 4
        ids = {c["clone_id"] for c in clones}
        assert ids == {"10004-heavy", "10004-light", "8232-heavy", "8232-light"}

    def test_pair_id_links_chains(self):
        _datasets, clones_dict, _trees = _run("paired-noinfo")
        by_id = {c["clone_id"]: c for c in _clones(clones_dict)}
        assert by_id["10004-heavy"]["pair_id"] == by_id["10004-light"]["pair_id"]
        assert by_id["10004-heavy"]["is_paired"] is True

    def test_per_locus_sequences_differ(self):
        """Heavy and light trees share topology but carry different sequences."""
        _datasets, _clones_dict, trees = _run("paired-noinfo")
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
        datasets, clones_dict, _trees = _run("paired-noinfo")
        loci = sorted(c["sample"]["locus"] for c in _clones(clones_dict))
        assert loci == ["igh", "igh", "igk", "igl"]
        # Dataset samples dedup by sample_id (all share repertoire_id "sample").
        assert len(datasets[0]["samples"]) == 1


@pytest.mark.airr
class TestDatasetSynthesis:
    def test_dataset_shape(self):
        datasets, clones_dict, _trees = _run("nocell-noinfo", name="my-airr")
        assert len(datasets) == 1
        ds = datasets[0]
        assert ds["name"] == "my-airr"
        assert ds["dataset_id"] in clones_dict
        assert ds["clone_count"] == 2
        assert ds["subjects"] == []
        assert "field_metadata" in ds

    def test_required_clone_fields_present(self):
        _datasets, clones_dict, _trees = _run("nocell-noinfo")
        for clone in _clones(clones_dict):
            assert clone["unique_seqs_count"] > 0
            assert "mean_mut_freq" in clone

    def test_node_type_label_does_not_collide_with_topology_type(self):
        """node_type must not share the "Node Type" label the webapp uses for
        the built-in topological ``type`` field.

        Two node fields with the same tooltip label produce a duplicate object
        key -> invalid Vega expression -> the webapp fails to render the tree.
        """
        datasets, _clones_dict, _trees = _run("paired-noinfo")
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


@pytest.mark.airr
class TestMetrics:
    def test_compute_metrics_populates_lbi(self):
        _datasets, _clones_dict, trees = _run("nocell-noinfo", compute_metrics=True)
        # At least some nodes get a numeric LBI when metrics are requested.
        lbis = [n["lbi"] for t in trees for n in t["nodes"]]
        assert any(v is not None for v in lbis)


@pytest.mark.airr
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

    def test_weights_by_real_multiplicity_when_present(self):
        """A node carrying Dowser's collapse_count (#45) is weighted by it,
        not treated as count-1, unlike the noinfo (unset multiplicity) case."""
        germline = "AAAA"
        nodes = [
            {"type": "leaf", "sequence_alignment": "AAAT", "multiplicity": 3},  # 1/4
            {"type": "leaf", "sequence_alignment": "AATT", "multiplicity": 1},  # 2/4
        ]
        # (0.25*3 + 0.5*1) / (3+1) = 1.25/4
        assert _mean_mutation_frequency(nodes, germline) == pytest.approx(0.3125)

    def test_unset_multiplicity_falls_back_to_one(self):
        """The noinfo schema (no per-node multiplicity) still gets an
        unweighted mean, matching pre-#45 behavior exactly."""
        germline = "AAAA"
        nodes = [
            {"type": "leaf", "sequence_alignment": "AAAT", "multiplicity": None},
            {"type": "leaf", "sequence_alignment": "AATT"},  # key absent entirely
        ]
        assert _mean_mutation_frequency(nodes, germline) == pytest.approx(0.375)


@pytest.mark.airr
class TestNodeMultiplicity:
    """Reading Dowser's collapse_count from the info catchall (#45)."""

    def test_reads_collapse_count_from_tipdata(self):
        node_record = {"info": {"tipdata": {"collapse_count": 5, "tip_order": 1}}}
        assert _node_multiplicity(node_record) == 5

    def test_none_when_info_absent(self):
        """The noinfo schema never has an `info` key at all."""
        assert _node_multiplicity({}) is None

    def test_none_when_info_is_empty_list(self):
        """Inferred/ASR nodes carry `info: []`, not a dict — a real shape
        gotcha in Dowser's output that must not raise or misread."""
        assert _node_multiplicity({"info": []}) is None

    def test_none_when_tipdata_missing(self):
        assert _node_multiplicity({"info": {}}) is None

    def test_none_when_collapse_count_missing(self):
        assert _node_multiplicity({"info": {"tipdata": {"tip_order": 1}}}) is None


@pytest.mark.airr
class TestCdrBoundariesFromRegion:
    """Deriving cdr1/cdr2/cdr3 alignment boundaries from Clone.info.region (#45)."""

    def test_no_gaps_maps_directly(self):
        region = ["fwr1"] * 3 + ["cdr1"] * 2 + ["fwr2"] * 3
        germline_alignment = "AAAAAAAA"  # 8 chars, no gaps, matches len(region)
        boundaries = _cdr_boundaries_from_region(region, germline_alignment)
        assert boundaries == {"cdr1": (3, 5)}

    def test_remaps_through_gap_padding(self):
        """Interior '.' padding in germline_alignment shifts boundaries, and
        a gap run is attributed to whichever region it falls within."""
        # Ungapped (real chars only): fwr1(3) cdr1(2) fwr2(3) -- matches
        # region below, one entry per real nucleotide.
        region = ["fwr1"] * 3 + ["cdr1"] * 2 + ["fwr2"] * 3
        # Gapped: "AAAA..AAAA" -- 2 gap chars inserted between the two real
        # cdr1 characters (positions 3 and 6).
        germline_alignment = "AAA" + "A..A" + "AAA"  # len 10
        boundaries = _cdr_boundaries_from_region(region, germline_alignment)
        # cdr1's real characters land at gapped positions 3 and 6; the next
        # region (fwr2) starts at gapped position 7 -- so the interior gap
        # pair (positions 4, 5) is attributed to cdr1, giving span [3, 7).
        assert boundaries["cdr1"] == (3, 7)

    def test_cdr3_extended_by_one_anchor_codon_each_side(self):
        """cdr3 is normalized to the junction convention (#46): region's
        strict IMGT cdr3 span gets +1 codon (3nt) on each side for the
        conserved anchor residues (V-gene Cys, J-gene Phe/Trp) junction
        includes and CDR3 excludes. cdr1/cdr2 have no such distinction."""
        region = ["fwr3"] * 6 + ["cdr3"] * 4 + ["fwr4"] * 6
        germline_alignment = "A" * 16
        boundaries = _cdr_boundaries_from_region(region, germline_alignment)
        # strict cdr3 = [6, 10); junction = [6-3, 10+3) = [3, 13)
        assert boundaries["cdr3"] == (3, 13)

    def test_cdr3_extension_clamps_at_sequence_boundaries(self):
        """A pathological cdr3 abutting the start/end of the sequence
        clamps rather than underflows/overflows."""
        region = ["cdr3"] * 2 + ["fwr4"] * 6
        germline_alignment = "A" * 8
        boundaries = _cdr_boundaries_from_region(region, germline_alignment)
        assert boundaries["cdr3"][0] == 0

    def test_matches_real_dowser_fixture(self):
        """Cross-check against the real nocell-info fixture's clone 10004,
        whose boundaries were independently hand-computed. cdr3 is the
        junction-normalized span (#46) -- strict IMGT cdr3 is (312, 363);
        junction extends it by 1 anchor codon (3nt) each side, matching this
        clone's junction_length of 57 (= 363-312 + 2*3) exactly."""
        data = _load("nocell-info")
        clone = data["Clone"][0]
        assert clone["clone_id"] == "10004"
        assert clone["junction_length"] == 57
        region = clone["info"]["region"]
        rear_by_id = {r.get("sequence_id"): r for r in data["Rearrangement"]}
        germline = rear_by_id[clone["inferred_ancestor"]]
        boundaries = _cdr_boundaries_from_region(region, germline["sequence_alignment"])
        assert boundaries == {
            "cdr1": (78, 114),
            "cdr2": (165, 195),
            "cdr3": (309, 366),
        }

    def test_only_cdr_regions_returned(self):
        region = ["fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4"]
        boundaries = _cdr_boundaries_from_region(region, "A" * 7)
        assert set(boundaries.keys()) == {"cdr1", "cdr2", "cdr3"}

    def test_empty_when_region_missing(self):
        assert _cdr_boundaries_from_region(None, "AAAA") == {}
        assert _cdr_boundaries_from_region([], "AAAA") == {}

    def test_empty_when_germline_alignment_missing(self):
        assert _cdr_boundaries_from_region(["fwr1"], None) == {}
        assert _cdr_boundaries_from_region(["fwr1"], "") == {}

    def test_empty_on_length_mismatch(self):
        """The paired heavy+light case (#45 known gap): region covers both
        chains concatenated, so it won't match either chain's alignment
        alone -- must skip safely, not misattribute boundaries."""
        region = ["fwr1"] * 10
        germline_alignment = "A" * 5  # deliberately wrong length
        assert _cdr_boundaries_from_region(region, germline_alignment) == {}


@pytest.mark.airr
class TestInfoCatchall:
    """End-to-end: Dowser's info catchall flows through to real output (#45)."""

    def test_collapse_count_becomes_node_multiplicity(self):
        data = _load("nocell-info")
        clone = data["Clone"][0]
        observed = next(n for n in clone["nodes"] if n["node_type"] == "observed")
        observed["info"]["tipdata"]["collapse_count"] = 7
        target_id = observed["sequence_id"]

        _datasets, _clones_dict, trees = process_airr_to_olmsted(
            [clone], data["Rearrangement"], minter=IdentMinter(seed=42), verbosity=0
        )
        node = next(
            n for t in trees for n in t["nodes"] if n["sequence_id"] == target_id
        )
        assert node["multiplicity"] == 7

    def test_inferred_node_multiplicity_stays_unset(self):
        """Inferred/ASR nodes have info: [] even in the info variant — no
        collapse_count to read, so multiplicity stays None, not fabricated."""
        _datasets, _clones_dict, trees = _run("nocell-info")
        inferred_nodes = [
            n for t in trees for n in t["nodes"] if n["type"] in ("root", "internal")
        ]
        assert inferred_nodes  # sanity: the tree actually has some
        assert all(n["multiplicity"] is None for n in inferred_nodes)

    def test_noinfo_variant_multiplicity_still_unset(self):
        """Unchanged behavior for the noinfo schema (no info key at all)."""
        _datasets, _clones_dict, trees = _run("nocell-noinfo")
        assert all(n["multiplicity"] is None for t in trees for n in t["nodes"])

    def test_cdr_boundaries_appear_on_output_clone(self):
        """Clone.info.region flows through to cdr1/cdr2/cdr3 alignment
        fields on unpaired (single-chain) clones."""
        _datasets, clones_dict, _trees = _run("nocell-info")
        clone = next(c for c in _clones(clones_dict) if c["clone_id"] == "10004")
        assert clone["cdr1_alignment_start"] == 78
        assert clone["cdr1_alignment_end"] == 114
        assert clone["cdr1_length"] == 36
        assert clone["cdr2_alignment_start"] == 165
        assert clone["cdr2_alignment_end"] == 195
        # cdr3 is junction-normalized (#46): region's strict IMGT cdr3 is
        # (312, 363)/51nt; +1 anchor codon (3nt) each side matches this
        # clone's junction_length of 57 exactly.
        assert clone["cdr3_alignment_start"] == 309
        assert clone["cdr3_alignment_end"] == 366
        assert clone["cdr3_length"] == 57
        # matches this clone's input junction_length exactly (57), by design
        input_clone = _load("nocell-info")["Clone"][0]
        assert input_clone["junction_length"] == 57
        assert clone["cdr3_length"] == input_clone["junction_length"]

    def test_warns_when_junction_normalized_cdr3_disagrees_with_junction_length(
        self, capsys
    ):
        """A genuine mismatch (not the already-accounted-for CDR3/junction
        anchor-codon difference) prints a warning naming the clone and both
        values, and the region-derived value still wins."""
        data = _load("nocell-info")
        clone = data["Clone"][0]
        clone["junction_length"] = 999  # deliberately wrong

        _datasets, clones_dict, _trees = process_airr_to_olmsted(
            [clone], data["Rearrangement"], minter=IdentMinter(seed=42), verbosity=1
        )
        out_clone = _clones(clones_dict)[0]
        assert out_clone["cdr3_length"] == 57  # region-derived value kept

        captured = capsys.readouterr()
        assert "10004" in captured.out
        assert "57" in captured.out
        assert "999" in captured.out

    def test_no_warning_when_junction_normalized_cdr3_matches(self, capsys):
        data = _load("nocell-info")
        process_airr_to_olmsted(
            data["Clone"],
            data["Rearrangement"],
            minter=IdentMinter(seed=42),
            verbosity=1,
        )
        captured = capsys.readouterr()
        assert "disagrees with junction_length" not in captured.out

    def test_noinfo_variant_has_no_cdr_alignment_fields(self):
        """No Clone.info.region in the noinfo schema -> no alignment fields,
        only the pre-existing junction_length-derived cdr3_length."""
        _datasets, clones_dict, _trees = _run("nocell-noinfo")
        for clone in _clones(clones_dict):
            assert "cdr1_alignment_start" not in clone
            assert "cdr2_alignment_start" not in clone
            assert "cdr3_alignment_start" not in clone

    def test_paired_clones_get_no_cdr_boundaries_yet(self):
        """Known gap (#45): a paired clone's region spans both chains
        concatenated, which doesn't match either chain's own
        germline_alignment length, so the safe (skip, not guess) behavior
        kicks in -- no cdr*_alignment fields on paired output, but the
        existing junction_length-derived cdr3_length is unaffected."""
        _datasets, clones_dict, _trees = _run("paired-info")
        for clone in _clones(clones_dict):
            assert "cdr1_alignment_start" not in clone
            assert "cdr2_alignment_start" not in clone
            assert "cdr3_alignment_start" not in clone
            assert "cdr3_length" in clone


@pytest.mark.airr
class TestGermlineGeneCallFallback:
    """v_call/d_call/j_call derived from the germline Rearrangement record
    (cf. #24/#45) — d_call's only source, and a fallback for v_call/j_call
    when the Clone-level info catchall doesn't supply them."""

    def _patch_germline_record(self, variant, **calls):
        """Load `variant`, inject `calls` onto the germline's own
        Rearrangement record, and return (clone, rearrangements)."""
        data = _load(variant)
        clone = data["Clone"][0]
        germline_id = clone["inferred_ancestor"]
        record = next(
            r for r in data["Rearrangement"] if r.get("sequence_id") == germline_id
        )
        record.update(calls)
        return clone, data["Rearrangement"]

    def test_d_call_has_no_clone_level_source_and_comes_from_germline_record(self):
        """d_call is never present on Clone (even in the info variant) — the
        germline Rearrangement record is the only possible source."""
        clone, rearrangements = self._patch_germline_record(
            "nocell-info", d_call="IGHD3-10*01"
        )
        assert clone.get("d_call") is None  # sanity: truly absent from Clone
        _datasets, clones_dict, _trees = process_airr_to_olmsted(
            [clone], rearrangements, minter=IdentMinter(seed=42), verbosity=0
        )
        assert _clones(clones_dict)[0]["d_call"] == "IGHD3-10*01"

    def test_v_call_j_call_fall_back_to_germline_record_when_clone_omits_them(self):
        """The noinfo variant's Clone has no v_call/j_call at all -- the
        germline record fills the gap when it has them."""
        clone, rearrangements = self._patch_germline_record(
            "nocell-noinfo", v_call="IGHV1-2*02", j_call="IGHJ4*02"
        )
        _datasets, clones_dict, _trees = process_airr_to_olmsted(
            [clone], rearrangements, minter=IdentMinter(seed=42), verbosity=0
        )
        out_clone = _clones(clones_dict)[0]
        assert out_clone["v_call"] == "IGHV1-2*02"
        assert out_clone["j_call"] == "IGHJ4*02"

    def test_clone_level_v_call_j_call_take_priority_over_germline_record(self):
        """When the info variant's Clone already supplies v_call/j_call,
        that value wins even if the germline record disagrees."""
        clone, rearrangements = self._patch_germline_record(
            "nocell-info", v_call="IGHV1-2*02", j_call="IGHJ4*02"
        )
        assert clone["v_call"] == "IGHV3-23"  # sanity: Clone-level value present
        _datasets, clones_dict, _trees = process_airr_to_olmsted(
            [clone], rearrangements, minter=IdentMinter(seed=42), verbosity=0
        )
        out_clone = _clones(clones_dict)[0]
        assert out_clone["v_call"] == "IGHV3-23"
        assert out_clone["j_call"] == "IGHJ3"

    def test_no_gene_calls_when_germline_record_has_none(self):
        """Unchanged behavior for the tracked fixtures as-is (germline
        records carry null v_call/d_call/j_call in the sample data)."""
        _datasets, clones_dict, _trees = _run("nocell-noinfo")
        out_clone = _clones(clones_dict)[0]
        assert out_clone.get("v_call") is None
        assert out_clone.get("d_call") is None
        assert out_clone.get("j_call") is None


@pytest.mark.airr
class TestMissingSequenceGraceful:
    def test_missing_rearrangement_yields_empty_sequence(self):
        """A node with no matching Rearrangement gets an empty sequence, not a crash."""
        data = _load("nocell-noinfo")
        # Drop one observed leaf's Rearrangement record.
        clone = data["Clone"][0]
        leaf = next(n for n in clone["nodes"] if n["node_type"] == "observed")
        dropped_id = leaf["sequence_id"]
        rearrangements = [
            r for r in data["Rearrangement"] if r.get("sequence_id") != dropped_id
        ]
        _datasets, _clones_dict, trees = process_airr_to_olmsted(
            [clone], rearrangements, minter=IdentMinter(seed=42), verbosity=0
        )
        node = next(
            n for t in trees for n in t["nodes"] if n["sequence_id"] == dropped_id
        )
        assert node["sequence_alignment"] == ""
        assert node["sequence_alignment_aa"] == ""
