"""Tests for olmsted_cli.metrics.compute_mean_mut_freq (issue 24).

``compute_mean_mut_freq`` is the single source of truth for ``mean_mut_freq``
across both the PCP and AIRR ingest paths — these tests cover the pure
function directly (known-answer, gap handling, boundary cases) plus one
end-to-end test that both formats agree when fed equivalent biological data.
"""

from __future__ import annotations

import pytest

from olmsted_cli.metrics import compute_mean_mut_freq
from olmsted_cli.process_airr_data import _mean_mutation_frequency
from olmsted_cli.process_pcp_data import parse_pcp_csv, process_pcp_to_olmsted

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


def test_cross_format_equivalence(tmp_path):
    """PCP and AIRR ingest paths agree on mean_mut_freq for equivalent data.

    AIRR's ``_mean_mutation_frequency`` has no input ``mean_mut_freq`` to
    override or warn about — the Clone/Rearrangement schema doesn't carry
    that field at all, unlike the legacy Olmsted-flavored AIRR container
    this format replaced (see issue #47) — so this only needs to check the
    computed value agrees with PCP's, not any override/warning behavior.
    """
    expected = 0.125

    pcp_path = tmp_path / "input-pcp.csv"
    pcp_path.write_text(PCP_CSV)
    families = parse_pcp_csv(str(pcp_path))
    _, pcp_clones_dict, _ = process_pcp_to_olmsted(families, verbosity=0)
    pcp_clone = next(c for clist in pcp_clones_dict.values() for c in clist)

    nodes = [
        _leaf("L1", "AATA", 2),
        _leaf("L2", "ATAA", 3),
        _leaf("L3", "AAAA", 5),
    ]
    airr_mean_mut_freq = _mean_mutation_frequency(nodes, "AAAA")

    assert pcp_clone["mean_mut_freq"] == pytest.approx(expected, abs=1e-9)
    assert airr_mean_mut_freq == pytest.approx(expected, abs=1e-9)
    assert pcp_clone["mean_mut_freq"] == pytest.approx(airr_mean_mut_freq, abs=1e-9)


def test_airr_no_usable_leaves_gives_zero():
    assert _mean_mutation_frequency([], "AAAA") == 0.0
