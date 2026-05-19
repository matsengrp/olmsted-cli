"""Phase 3 PCP streaming pipeline (#26): batch-size equivalence and hoist correctness.

These tests exercise ``process_pcp_format`` through the streaming path
end-to-end:

- For every PCP example dataset, running the CLI with ``--batch-size`` in
  ``{1, 2, 50, 10000}`` produces output that is structurally equal to
  the existing golden (which the legacy path generated).  Batching
  shouldn't perturb data shape, only how it's assembled.
- A synthetic mixed-variance input verifies that a tree-csv extra
  constant across most clones' trees but varying within one multi-tree
  clone correctly classifies as tree-level dataset-wide, even when the
  varying clone lands in a separate batch from the constant ones.
"""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pytest

from .test_cli_processing import (
    compare_consolidated_files,
)

EXAMPLE_DATA = Path(__file__).parent.parent / "example-data"


PCP_DATASETS = [
    pytest.param(
        EXAMPLE_DATA / "pcp" / "input-pcp.csv",
        EXAMPLE_DATA / "pcp" / "input-trees.csv",
        EXAMPLE_DATA / "pcp" / "pcp-olmsted-golden.json",
        "pcp-example",
        id="pcp",
    ),
    pytest.param(
        EXAMPLE_DATA / "pcp-byhand" / "input-pcp.csv",
        EXAMPLE_DATA / "pcp-byhand" / "input-trees.csv",
        EXAMPLE_DATA / "pcp-byhand" / "pcp-byhand-olmsted-golden.json",
        "pcp-byhand-example",
        id="pcp-byhand",
    ),
    pytest.param(
        EXAMPLE_DATA / "pcp-light" / "input-pcp.csv",
        EXAMPLE_DATA / "pcp-light" / "input-trees.csv",
        EXAMPLE_DATA / "pcp-light" / "pcp-light-olmsted-golden.json",
        "pcp-light-example",
        id="pcp-light",
    ),
    pytest.param(
        EXAMPLE_DATA / "pcp-paired" / "input-pcp.csv",
        EXAMPLE_DATA / "pcp-paired" / "input-trees.csv",
        EXAMPLE_DATA / "pcp-paired" / "pcp-paired-olmsted-golden.json",
        "pcp-paired-example",
        id="pcp-paired",
    ),
]


def _run_cli(inputs, trees, output, name, batch_size):
    cmd = [
        "olmsted",
        "process",
        "-f",
        "pcp",
        "-i",
        str(inputs),
        "-t",
        str(trees),
        "-o",
        str(output),
        "--seed",
        "42",
        "--name",
        name,
        "--batch-size",
        str(batch_size),
        "-q",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, (
        f"CLI failed (batch_size={batch_size}):\nSTDOUT: {result.stdout}\n"
        f"STDERR: {result.stderr}"
    )


@pytest.mark.parametrize("pcp_input,trees_input,golden,name", PCP_DATASETS)
@pytest.mark.parametrize("batch_size", [1, 2, 50, 10000])
def test_pcp_batch_size_matches_golden(
    pcp_input, trees_input, golden, name, batch_size, tmp_path
):
    out = tmp_path / f"out_bs{batch_size}.json"
    _run_cli(pcp_input, trees_input, out, name, batch_size)
    ok, message = compare_consolidated_files(out, golden)
    assert ok, f"batch_size={batch_size} diverges from golden:\n{message}"


# Synthetic PCP input where method_score is constant across alt
# reconstructions on clones 1..9 but varies within clone 10's two trees.
# The dataset-scope rule says method_score is tree-level for the whole
# dataset; verify that classification holds even when batches split
# the varying and non-varying clones apart.


def _pcp_block(family: str, tree_name: str, method_score: float) -> str:
    seq = "GAATTCAAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGG"
    return (
        f"S1,{family},{tree_name},naive,mrca,{seq},{seq},true,false,IGHV1*01,IGHJ1*01,{method_score}\n"
        f"S1,{family},{tree_name},mrca,L1,{seq},{seq[:5]}CCAAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGG,false,true,IGHV1*01,IGHJ1*01,{method_score}\n"
    )


def _build_synthetic_pcp(tmp_path: Path) -> tuple[Path, Path]:
    header = (
        "sample_id,family,tree_name,parent_name,child_name,parent_heavy,"
        "child_heavy,parent_is_naive,child_is_leaf,v_gene_heavy,j_gene_heavy,"
        "method_score\n"
    )
    pcp_lines = [header]
    trees_header = "family_name,sample_id,tree_name,newick_tree,method_score\n"
    trees_lines = [trees_header]

    # Clones F1..F9: two alt reconstructions per clone with the same method_score
    for i in range(1, 10):
        family = f"F{i}"
        for tree_name, method_score in [("T_a", 0.5), ("T_b", 0.5)]:
            pcp_lines.append(_pcp_block(family, tree_name, method_score))
            trees_lines.append(
                f'{family},S1,{tree_name},"((L1:0.2)mrca:0.1)naive;",{method_score}\n'
            )

    # Clone F10: two alt reconstructions with DIFFERENT method_score values
    for tree_name, method_score in [("T_a", 0.1), ("T_b", 0.9)]:
        pcp_lines.append(_pcp_block("F10", tree_name, method_score))
        trees_lines.append(
            f'F10,S1,{tree_name},"((L1:0.2)mrca:0.1)naive;",{method_score}\n'
        )

    pcp_path = tmp_path / "input-pcp.csv"
    pcp_path.write_text("".join(pcp_lines))
    trees_path = tmp_path / "input-trees.csv"
    trees_path.write_text("".join(trees_lines))
    return pcp_path, trees_path


@pytest.mark.parametrize("batch_size", [1, 2, 3, 5, 10000])
def test_streaming_hoist_respects_dataset_scope_variance(batch_size, tmp_path):
    """``method_score`` varies in one clone — it must classify as tree-level
    dataset-wide, regardless of which batch boundary the varying clone falls on.
    """
    pcp_path, trees_path = _build_synthetic_pcp(tmp_path)
    out = tmp_path / f"out_bs{batch_size}.json"
    _run_cli(pcp_path, trees_path, out, "hoist-test", batch_size)

    with out.open() as fh:
        data = json.load(fh)

    field_metadata = data["datasets"][0]["field_metadata"]

    assert "tree" in field_metadata, "expected tree-level field_metadata"
    assert "method_score" in field_metadata["tree"], (
        f"method_score must be tree-level (batch_size={batch_size}); got "
        f"tree fields {list(field_metadata['tree'].keys())}"
    )
    assert "method_score" not in field_metadata.get("clone", {}), (
        f"method_score must not be hoisted to clone level (batch_size={batch_size})"
    )

    # No clone should carry method_score (it stays on tree records).
    for clones in data["clones"].values():
        for clone in clones:
            assert "method_score" not in clone, (
                f"clone {clone.get('clone_id')} unexpectedly has method_score "
                f"hoisted (batch_size={batch_size})"
            )

    # Every tree should still carry method_score.
    for tree in data["trees"]:
        assert "method_score" in tree, (
            f"tree {tree.get('tree_id')} missing method_score (batch_size={batch_size})"
        )


def test_streaming_default_path_uses_streaming_for_pcp(tmp_path):
    """At the default batch_size the streaming code path is what runs.

    Catches regressions where _should_stream_pcp silently flips back to
    legacy for the common case.
    """
    import argparse

    from olmsted_cli.process_data import _should_stream_pcp

    args = argparse.Namespace(
        batch_size=50,
        split_files=None,
        mutations=None,
        validate=False,
    )
    assert _should_stream_pcp(args) is True

    # Mutations: stays on streaming after phase 4 wired per-batch merge.
    args.mutations = "x.csv"
    assert _should_stream_pcp(args) is True
    args.mutations = None

    # Validate: streaming bypassed until per-batch validation lands.
    args.validate = True
    assert _should_stream_pcp(args) is False
    args.validate = False

    # Split-files: streaming bypassed.
    args.split_files = "/tmp/x"
    assert _should_stream_pcp(args) is False
    args.split_files = None

    # batch_size 0: explicit opt-out.
    args.batch_size = 0
    assert _should_stream_pcp(args) is False
