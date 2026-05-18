"""Phase 4 per-batch ``--mutations`` (#26): equivalence across batch sizes.

The streaming pipeline loads the mutations CSV once and threads a single
:class:`MergeContext` through every batch.  Stats and the unmatched-family
set aggregate across calls so the final summary matches the one-shot
legacy path.  These tests verify the end-to-end output is byte-equivalent
regardless of ``--batch-size``.
"""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pytest

from .test_cli_processing import (
    normalize_json,
    strip_volatile_fields,
)

# Two PCP families, each with one tree that produces a K→R substitution
# at AA site 0 (M→V in the actual translation, but we just need a
# mutation to attach scores to).  The CSV provides a surprise_mutsel
# score keyed by ``(family, site, parent_aa, child_aa)`` so we exercise
# the site-keyed match mode.

PCP_CSV = """sample_id,family,parent_name,child_name,parent_heavy,child_heavy,parent_is_naive,child_is_leaf,v_gene_heavy,j_gene_heavy
S1,F1,naive,L1,ATGAAA,GTGAAA,true,true,IGHV1*01,IGHJ1*01
S1,F2,naive,L2,ATGCCC,GTGCCC,true,true,IGHV2*01,IGHJ2*01
"""

# Diff of ATGAAA → GTGAAA in AA: site 0 M (ATG) → V (GTG); rest matches.
MUTATIONS_CSV = """family,site,parent_aa,child_aa,surprise_mutsel
S1_F1,0,M,V,2.5
S1_F2,0,M,V,4.2
"""


def _write_inputs(tmp_path: Path) -> tuple[Path, Path]:
    pcp_path = tmp_path / "input-pcp.csv"
    pcp_path.write_text(PCP_CSV)
    mutations_path = tmp_path / "mutations.csv"
    mutations_path.write_text(MUTATIONS_CSV)
    return pcp_path, mutations_path


def _run(pcp_path, mutations_path, output, batch_size):
    cmd = [
        "olmsted",
        "process",
        "-f",
        "pcp",
        "-i",
        str(pcp_path),
        "--mutations",
        str(mutations_path),
        "-o",
        str(output),
        "--seed",
        "42",
        "--name",
        "mut-batching",
        "--batch-size",
        str(batch_size),
        "-q",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, (
        f"CLI failed (batch_size={batch_size}):\n"
        f"STDOUT: {result.stdout}\nSTDERR: {result.stderr}"
    )


def _load_stripped(path: Path):
    with path.open() as fh:
        return normalize_json(strip_volatile_fields(json.load(fh)))


@pytest.mark.parametrize("batch_size", [1, 2, 10000])
def test_process_mutations_streaming_matches_legacy(batch_size, tmp_path):
    """``process --mutations`` at ``--batch-size N`` matches the legacy run.

    ``--batch-size 0`` routes to the legacy in-memory pipeline; every
    other size goes through the streaming path with per-batch
    :func:`apply_mutations_to_trees`.  Outputs must be structurally equal.
    """
    pcp_path, mutations_path = _write_inputs(tmp_path)

    legacy_out = tmp_path / "legacy.json"
    streaming_out = tmp_path / f"streaming_bs{batch_size}.json"
    _run(pcp_path, mutations_path, legacy_out, 0)
    _run(pcp_path, mutations_path, streaming_out, batch_size)

    legacy = _load_stripped(legacy_out)
    streaming = _load_stripped(streaming_out)
    assert streaming == legacy, (
        f"batch_size={batch_size} diverges from legacy --batch-size 0"
    )


def test_process_mutations_streaming_attaches_score(tmp_path):
    """The streaming path actually enriches mutations — not just byte-equivalent."""
    pcp_path, mutations_path = _write_inputs(tmp_path)
    out = tmp_path / "out.json"
    _run(pcp_path, mutations_path, out, batch_size=1)

    data = json.loads(out.read_text())

    enriched_scores = []
    for tree in data["trees"]:
        for node in tree.get("nodes", []) or []:
            for mut in node.get("mutations", []) or []:
                if "surprise_mutsel" in mut:
                    enriched_scores.append((tree["clone_id"], mut["surprise_mutsel"]))

    enriched_scores.sort()
    assert enriched_scores == [("S1_F1", 2.5), ("S1_F2", 4.2)]

    # field_metadata reflects the new mutation field
    fm = data["datasets"][0]["field_metadata"]["mutation"]
    assert "surprise_mutsel" in fm


def test_process_mutations_streaming_reports_unmatched_family(tmp_path):
    """Unmatched-family warnings accumulate across batches and print at finalize."""
    pcp_path = tmp_path / "input-pcp.csv"
    pcp_path.write_text(PCP_CSV)

    mutations_path = tmp_path / "mutations.csv"
    mutations_path.write_text(
        "family,site,parent_aa,child_aa,surprise_mutsel\n"
        "S1_F1,0,M,V,2.5\n"
        "S1_GHOST,0,K,R,9.9\n"  # family with no matching clone
    )

    out = tmp_path / "out.json"
    cmd = [
        "olmsted",
        "process",
        "-f",
        "pcp",
        "-i",
        str(pcp_path),
        "--mutations",
        str(mutations_path),
        "-o",
        str(out),
        "--seed",
        "42",
        "--name",
        "mut-unmatched",
        "--batch-size",
        "1",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr

    combined = result.stdout + result.stderr
    assert "S1_GHOST" in combined or "had no matching clone" in combined
