"""Round-trip tests for the per-clone-group iterators.

Phase 1 of the streaming-batch refactor introduces ``iter_pcp_clone_groups``
(PCP). These tests assert the property the streaming pipeline depends on:
iterating with any ``batch_size`` produces the same concatenated clones and
trees as iterating with ``batch_size=None``. (AIRR has no equivalent
iterator — it has no streaming/batching support at all, see issue #36.)
"""

from __future__ import annotations

import copy

from olmsted_cli.identifier import IdentMinter
from olmsted_cli.process_pcp_data import (
    TreeProcessingConfig,
    iter_pcp_clone_groups,
    parse_newick_csv,
    parse_pcp_csv,
)

# Two-tree-per-family, two-family PCP input gives the iterator enough material
# to slice at every meaningful boundary (1, 2, all-at-once).
PCP_CSV = """sample_id,family,tree_name,parent_name,child_name,parent_heavy,child_heavy,parent_is_naive,child_is_leaf,v_gene_heavy,j_gene_heavy
S1,F1,T_a,naive,mrca-1,GAATTCAAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGG,GAATTCAAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGG,true,false,IGHV1*01,IGHJ1*01
S1,F1,T_a,mrca-1,L1,GAATTCAAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGG,GAATTCCAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGG,false,true,IGHV1*01,IGHJ1*01
S1,F1,T_a,mrca-1,L2,GAATTCAAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGG,GAATTCAAATAGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGG,false,true,IGHV1*01,IGHJ1*01
S1,F1,T_b,naive,N1,GAATTCAAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGG,GAATTCAAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGG,true,false,IGHV1*01,IGHJ1*01
S1,F1,T_b,N1,L1,GAATTCAAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGG,GAATTCCAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGG,false,true,IGHV1*01,IGHJ1*01
S1,F1,T_b,N1,L2,GAATTCAAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGG,GAATTCAAATAGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGG,false,true,IGHV1*01,IGHJ1*01
S1,F2,T_a,naive,mrca-2,GAATTCAAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGG,GAATTCAAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGG,true,false,IGHV2*01,IGHJ2*01
S1,F2,T_a,mrca-2,L3,GAATTCAAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGG,GAATTCCAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGG,false,true,IGHV2*01,IGHJ2*01
S1,F2,T_a,mrca-2,L4,GAATTCAAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGG,GAATTCAAATAGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGG,false,true,IGHV2*01,IGHJ2*01
"""

TREES_CSV = """family_name,sample_id,tree_name,newick_tree,reconstruction_method
F1,S1,T_a,"((L1:0.2,L2:0.2)mrca-1:0.1)naive;",ground-truth
F1,S1,T_b,"((L1:0.2,L2:0.2)N1:0.1)naive;",bcrlarch-parsimony
F2,S1,T_a,"((L3:0.2,L4:0.2)mrca-2:0.1)naive;",ground-truth
"""


def _drain(iterator):
    """Consume an iterator, returning concatenated clones, trees, and yield count."""
    clones, trees = [], []
    yields = 0
    for batch_clones, batch_trees in iterator:
        clones.extend(batch_clones)
        trees.extend(batch_trees)
        yields += 1
    return clones, trees, yields


def _run_pcp(pcp_families, newick_trees, batch_size):
    """Drive iter_pcp_clone_groups with a fresh deterministic minter.

    Mints the two dataset idents that the legacy ``process_pcp_to_olmsted``
    consumes before entering its clone loop, so subsequent clone/sample/tree
    idents land at the same minter positions as the real pipeline.
    """
    minter = IdentMinter(seed=42)
    dataset_id = minter.mint("dataset")
    _ = minter.mint("dataset")
    config = TreeProcessingConfig()
    samples: list = []

    clones, trees, yields = _drain(
        iter_pcp_clone_groups(
            pcp_families,
            newick_trees,
            minter,
            dataset_id,
            samples,
            config,
            batch_size=batch_size,
            progress=False,
        )
    )
    return clones, trees, samples, yields


def test_iter_pcp_clone_groups_concatenation_invariant(tmp_path):
    pcp = tmp_path / "input-pcp.csv"
    pcp.write_text(PCP_CSV)
    trees = tmp_path / "input-trees.csv"
    trees.write_text(TREES_CSV)

    pcp_families = parse_pcp_csv(str(pcp))
    newick_trees = parse_newick_csv(str(trees))

    baseline = _run_pcp(copy.deepcopy(pcp_families), copy.deepcopy(newick_trees), None)
    base_clones, base_trees, base_samples, base_yields = baseline
    assert base_yields == 1

    for size in (1, 2, 10):
        run = _run_pcp(copy.deepcopy(pcp_families), copy.deepcopy(newick_trees), size)
        clones, trees_out, samples, yields = run

        assert clones == base_clones, f"batch_size={size}: clones differ"
        assert trees_out == base_trees, f"batch_size={size}: trees differ"
        assert samples == base_samples, f"batch_size={size}: samples differ"
        assert yields >= 1


def test_iter_pcp_clone_groups_groups_alt_reconstructions_in_same_yield(tmp_path):
    """A clone with multiple tree reconstructions must co-emit them in one yield."""
    pcp = tmp_path / "input-pcp.csv"
    pcp.write_text(PCP_CSV)
    trees = tmp_path / "input-trees.csv"
    trees.write_text(TREES_CSV)

    pcp_families = parse_pcp_csv(str(pcp))
    newick_trees = parse_newick_csv(str(trees))

    minter = IdentMinter(seed=42)
    dataset_id = minter.mint("dataset")
    _ = minter.mint("dataset")
    config = TreeProcessingConfig()
    samples: list = []

    # batch_size=1 → one clone per yield. F1 has two trees; both must
    # appear in the same yield as the F1 clone.
    yielded_batches = list(
        iter_pcp_clone_groups(
            pcp_families,
            newick_trees,
            minter,
            dataset_id,
            samples,
            config,
            batch_size=1,
            progress=False,
        )
    )
    assert len(yielded_batches) == 2  # F1 and F2

    f1_clones, f1_trees = yielded_batches[0]
    assert len(f1_clones) == 1
    assert f1_clones[0]["clone_id"] == "S1_F1"
    assert len(f1_trees) == 2  # T_a and T_b

    f2_clones, f2_trees = yielded_batches[1]
    assert len(f2_clones) == 1
    assert f2_clones[0]["clone_id"] == "S1_F2"
    assert len(f2_trees) == 1
