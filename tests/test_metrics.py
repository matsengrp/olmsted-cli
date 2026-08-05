"""Tests for olmsted_cli.metrics.compute_mean_mut_freq (issue 24).

``compute_mean_mut_freq`` is the single source of truth for ``mean_mut_freq``
across both the PCP and AIRR ingest paths — these tests cover the pure
function directly (known-answer, gap handling, boundary cases) plus one
end-to-end test that both formats agree when fed equivalent biological data.
"""

from __future__ import annotations

import argparse

import pytest

from olmsted_cli.identifier import IdentMinter
from olmsted_cli.metrics import compute_mean_mut_freq
from olmsted_cli.process_airr_data import _process_airr_clone
from olmsted_cli.process_pcp_data import parse_pcp_csv, process_pcp_to_olmsted
from olmsted_cli.utils import set_verbosity

TOLERANCE = 1e-12


def _leaf(sequence_id, sequence, multiplicity, node_type="leaf"):
    return {
        "sequence_id": sequence_id,
        "type": node_type,
        "multiplicity": multiplicity,
        "sequence_alignment": sequence,
    }


def test_known_answer():
    """germline=AAAA, leaves [AATA mult=2, ATAA mult=3, AAAA mult=5]."""
    nodes = [
        _leaf("L1", "AATA", 2),
        _leaf("L2", "ATAA", 3),
        _leaf("L3", "AAAA", 5),
    ]
    mean_mut_freq, debug_info, skipped_nodes = compute_mean_mut_freq("AAAA", nodes)

    expected = (0.25 * 2 + 0.25 * 3 + 0.0 * 5) / 10
    assert mean_mut_freq == pytest.approx(expected, abs=TOLERANCE)
    assert len(debug_info) == 3
    assert skipped_nodes == []


def test_gap_position_symmetry():
    """A gap-only difference (no other mutations) gives 0.0 regardless of
    which side the gap is on, or whether the position is removed entirely.

    "." is the gap sentinel this codebase's alignment convention uses (see
    align_and_calculate_mutations's "pad" mode); "-" is not treated as a
    gap and would count as a real mismatch, so it's deliberately not used
    here — that's existing, preserved behavior, not something this issue
    changes.
    """
    gap_on_leaf, _, _ = compute_mean_mut_freq("AAAA", [_leaf("L", "A.AA", 1)])
    gap_on_germline, _, _ = compute_mean_mut_freq("A.AA", [_leaf("L", "AAAA", 1)])
    position_removed, _, _ = compute_mean_mut_freq("AAA", [_leaf("L", "AAA", 1)])

    assert gap_on_leaf == pytest.approx(0.0, abs=TOLERANCE)
    assert gap_on_germline == pytest.approx(0.0, abs=TOLERANCE)
    assert position_removed == pytest.approx(0.0, abs=TOLERANCE)


def test_identical_leaves_give_zero():
    nodes = [_leaf("L1", "AAAA", 1), _leaf("L2", "AAAA", 3)]
    mean_mut_freq, _, _ = compute_mean_mut_freq("AAAA", nodes)
    assert mean_mut_freq == pytest.approx(0.0, abs=TOLERANCE)


def test_fully_diverged_leaves_give_one():
    nodes = [_leaf("L1", "TTTT", 1), _leaf("L2", "TTTT", 4)]
    mean_mut_freq, _, _ = compute_mean_mut_freq("AAAA", nodes)
    assert mean_mut_freq == pytest.approx(1.0, abs=TOLERANCE)


def test_no_leaves_gives_zero_and_does_not_raise():
    mean_mut_freq, debug_info, skipped_nodes = compute_mean_mut_freq("AAAA", [])
    assert mean_mut_freq == 0.0
    assert debug_info == []
    assert skipped_nodes == []


def test_empty_germline_gives_zero_and_does_not_raise():
    nodes = [_leaf("L1", "AATA", 2)]
    mean_mut_freq, debug_info, skipped_nodes = compute_mean_mut_freq("", nodes)
    assert mean_mut_freq == 0.0
    assert debug_info == []
    assert len(skipped_nodes) == 1


def test_all_zero_multiplicity_gives_zero_and_does_not_raise():
    nodes = [_leaf("L1", "AATA", 0), _leaf("L2", "ATAA", 0)]
    mean_mut_freq, debug_info, skipped_nodes = compute_mean_mut_freq("AAAA", nodes)
    assert mean_mut_freq == 0.0
    assert debug_info == []
    assert len(skipped_nodes) == 2


def test_non_leaf_nodes_are_ignored():
    nodes = [
        _leaf("naive", "AAAA", 0, node_type="root"),
        _leaf("N1", "AAAA", 0, node_type="internal"),
        _leaf("L1", "AATA", 2),
    ]
    mean_mut_freq, debug_info, skipped_nodes = compute_mean_mut_freq("AAAA", nodes)
    assert mean_mut_freq == pytest.approx(0.25, abs=TOLERANCE)
    assert len(debug_info) == 1
    assert len(skipped_nodes) == 2


# ---------------------------------------------------------------------------
# Cross-format equivalence: same biological data, through the full ingest
# pipeline for both PCP and AIRR, must agree to within 1e-9.
# ---------------------------------------------------------------------------

# naive=AAAA; three leaves with the same sequences/multiplicities as the
# known-answer test above (expected mean_mut_freq = 0.125).
PCP_CSV = """sample_id,family,parent_name,child_name,parent_heavy,child_heavy,parent_is_naive,child_is_leaf,sample_count,v_gene_heavy,j_gene_heavy
S1,F1,naive,L1,AAAA,AATA,true,true,2,IGHV1*01,IGHJ1*01
S1,F1,naive,L2,AAAA,ATAA,true,true,3,IGHV1*01,IGHJ1*01
S1,F1,naive,L3,AAAA,AAAA,true,true,5,IGHV1*01,IGHJ1*01
"""


def _airr_args(seed=42):
    args = argparse.Namespace()
    args.minter = IdentMinter(seed=seed)
    args.verbose = 0
    args.compute_metrics = False
    args.lbi_tau = 0.0125
    args.naive_name = "naive"
    args.root_trees = False
    args.custom_fields = None
    return args


def _airr_dataset():
    return {
        "dataset_id": "test-ds",
        "samples": [{"sample_id": "S1"}],
        "clones": [
            {
                "clone_id": "clust-test",
                "sample_id": "S1",
                "germline_alignment": "AAAA",
                "trees": [
                    {
                        "newick": "(L1:0.1,L2:0.1,L3:0.1)naive;",
                        "nodes": {
                            "naive": {"sequence_alignment": "AAAA", "multiplicity": 0},
                            "L1": {"sequence_alignment": "AATA", "multiplicity": 2},
                            "L2": {"sequence_alignment": "ATAA", "multiplicity": 3},
                            "L3": {"sequence_alignment": "AAAA", "multiplicity": 5},
                        },
                    }
                ],
            }
        ],
    }


def test_cross_format_equivalence(tmp_path):
    expected = 0.125

    pcp_path = tmp_path / "input-pcp.csv"
    pcp_path.write_text(PCP_CSV)
    families = parse_pcp_csv(str(pcp_path))
    _, pcp_clones_dict, _ = process_pcp_to_olmsted(families, verbosity=0)
    pcp_clone = next(c for clist in pcp_clones_dict.values() for c in clist)

    dataset = _airr_dataset()
    airr_clone, _ = _process_airr_clone(_airr_args(), dataset, dataset["clones"][0])

    assert pcp_clone["mean_mut_freq"] == pytest.approx(expected, abs=1e-9)
    assert airr_clone["mean_mut_freq"] == pytest.approx(expected, abs=1e-9)
    assert pcp_clone["mean_mut_freq"] == pytest.approx(
        airr_clone["mean_mut_freq"], abs=1e-9
    )


def test_airr_overwrites_input_mean_mut_freq():
    """Any mean_mut_freq present on the input clone is ignored/overwritten."""
    dataset = _airr_dataset()
    dataset["clones"][0]["mean_mut_freq"] = 0.999  # bogus upstream value
    airr_clone, _ = _process_airr_clone(_airr_args(), dataset, dataset["clones"][0])
    assert airr_clone["mean_mut_freq"] == pytest.approx(0.125, abs=1e-9)


def test_airr_warns_on_input_mean_mut_freq_divergence(capsys):
    """A real convention mismatch (weighted vs unweighted) prints a warning
    at default verbosity, naming the clone and both values."""
    set_verbosity(1)
    dataset = _airr_dataset()
    dataset["clones"][0]["mean_mut_freq"] = 0.5  # far from the computed 0.125
    _process_airr_clone(_airr_args(), dataset, dataset["clones"][0])

    captured = capsys.readouterr()
    assert "clust-test" in captured.out
    assert "0.500000" in captured.out
    assert "0.125000" in captured.out


def test_airr_no_warning_when_input_matches_computed(capsys):
    set_verbosity(1)
    dataset = _airr_dataset()
    dataset["clones"][0]["mean_mut_freq"] = 0.125  # already the computed value
    _process_airr_clone(_airr_args(), dataset, dataset["clones"][0])

    captured = capsys.readouterr()
    assert "differs from the recomputed value" not in captured.out


def test_airr_no_warning_when_input_has_no_mean_mut_freq(capsys):
    set_verbosity(1)
    dataset = _airr_dataset()
    assert "mean_mut_freq" not in dataset["clones"][0]
    _process_airr_clone(_airr_args(), dataset, dataset["clones"][0])

    captured = capsys.readouterr()
    assert "differs from the recomputed value" not in captured.out


def test_airr_no_usable_leaves_gives_zero():
    dataset = {
        "dataset_id": "test-ds",
        "samples": [{"sample_id": "S1"}],
        "clones": [
            {
                "clone_id": "clust-empty",
                "sample_id": "S1",
                "germline_alignment": "AAAA",
                "trees": [
                    {
                        "newick": "naive;",
                        "nodes": {
                            "naive": {"sequence_alignment": "AAAA", "multiplicity": 0}
                        },
                    }
                ],
            }
        ],
    }
    airr_clone, _ = _process_airr_clone(_airr_args(), dataset, dataset["clones"][0])
    assert airr_clone["mean_mut_freq"] == 0.0
