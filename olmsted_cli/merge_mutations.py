"""
Shared utilities for merging external mutation-level CSV data into Olmsted JSON.

Used by:
- The `merge` command (post-hoc merging into existing Olmsted JSON)
- The `process` command's `--mutations` flag (during initial processing)

The CSV is expected to have at minimum these columns:
    family, site, parent_aa, child_aa
plus any number of score/annotation columns.

Mutations are matched to specific tree nodes by deriving each node's mutations
from its AA sequence diff against its parent, then matching CSV rows by
(clone_id, site, parent_aa, child_aa).
"""

import csv
import gzip
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from .constants import MUTATIONS_CSV_KEY_COLUMNS
from .process_pcp_data import _coerce_csv_value
from .utils import vprint

# Characters in AA sequences that should not be treated as a mutation event.
_GAP_AA_CHARS = {"-", ".", "X", "*", "?"}


def load_mutations_csv(csv_path: str) -> Dict[str, List[Dict[str, Any]]]:
    """Load a mutations CSV and group rows by family.

    Args:
        csv_path: Path to the mutations CSV (may be gzipped).

    Returns:
        Dict mapping family/clone_id -> list of row dicts. Each row dict
        contains the mutation site keys (site, parent_aa, child_aa) plus all
        score/annotation columns. Structural join columns (family, sample_id,
        pcp_index, depth) are stripped.
    """
    path = Path(csv_path)
    if not path.exists():
        raise FileNotFoundError(f"Mutations CSV not found: {csv_path}")

    if str(path).endswith(".gz"):
        file_handle = gzip.open(path, "rt")
    else:
        file_handle = open(path, "r")

    mutations_by_family: Dict[str, List[Dict[str, Any]]] = {}

    with file_handle:
        reader = csv.DictReader(file_handle)
        if not reader.fieldnames:
            raise ValueError(f"Mutations CSV has no header: {csv_path}")

        required = {"family", "site", "parent_aa", "child_aa"}
        missing = required - set(reader.fieldnames)
        if missing:
            raise ValueError(
                f"Mutations CSV missing required columns: {sorted(missing)}"
            )

        for row in reader:
            family = row.get("family")
            if not family:
                continue

            mutation: Dict[str, Any] = {}
            for key, val in row.items():
                if key is None or val == "" or val is None:
                    continue
                # Always coerce site/parent_aa/child_aa for matching
                if key == "site":
                    mutation["site"] = _coerce_csv_value(val)
                elif key in ("parent_aa", "child_aa"):
                    mutation[key] = val
                elif key in MUTATIONS_CSV_KEY_COLUMNS:
                    # Skip other key/structural columns
                    continue
                else:
                    mutation[key] = _coerce_csv_value(val)

            mutations_by_family.setdefault(family, []).append(mutation)

    return mutations_by_family


def derive_node_mutations(
    node: Dict[str, Any], parent_node: Optional[Dict[str, Any]]
) -> List[Dict[str, Any]]:
    """Derive mutations for a node by diffing its AA sequence against its parent's.

    Args:
        node: Node dict (must have ``sequence_alignment_aa``).
        parent_node: Parent node dict (must have ``sequence_alignment_aa``).
            If None, returns an empty list.

    Returns:
        List of mutation dicts, each with ``site`` (0-based position),
        ``parent_aa``, and ``child_aa``. Skips positions where either
        residue is a gap character.
    """
    if parent_node is None:
        return []

    child_seq = node.get("sequence_alignment_aa") or ""
    parent_seq = parent_node.get("sequence_alignment_aa") or ""
    if not child_seq or not parent_seq:
        return []

    mutations = []
    for i, (p, c) in enumerate(zip(parent_seq, child_seq)):
        if p == c:
            continue
        if p in _GAP_AA_CHARS or c in _GAP_AA_CHARS:
            continue
        mutations.append({"site": i, "parent_aa": p, "child_aa": c})
    return mutations


def _build_parent_lookup(nodes: List[Dict[str, Any]]) -> Dict[str, Dict[str, Any]]:
    """Build a sequence_id -> node lookup for parent resolution."""
    return {
        n["sequence_id"]: n for n in nodes if isinstance(n, dict) and "sequence_id" in n
    }


def _normalize_nodes(tree: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Return tree nodes as a list (handling both dict and list representations)."""
    nodes = tree.get("nodes", [])
    if isinstance(nodes, dict):
        return list(nodes.values())
    return nodes


def merge_mutations_into_trees(
    trees: List[Dict[str, Any]],
    mutations_by_family: Dict[str, List[Dict[str, Any]]],
) -> Tuple[int, int, int]:
    """Merge mutation-level CSV data into tree nodes (in place).

    For each tree whose ``clone_id`` matches a family in the CSV:
      1. Index the family's CSV rows by (site, parent_aa, child_aa)
      2. For each node, derive (or read existing) mutation records
      3. Match each mutation against the CSV index and merge any extra fields

    Args:
        trees: List of tree dicts (modified in place).
        mutations_by_family: Output of ``load_mutations_csv``.

    Returns:
        Tuple of (trees_matched, nodes_with_mutations, mutations_merged).
    """
    trees_matched = 0
    nodes_with_mutations = 0
    mutations_merged = 0

    unmatched_families = set(mutations_by_family.keys())

    for tree in trees:
        clone_id = tree.get("clone_id")
        if not clone_id or clone_id not in mutations_by_family:
            continue

        unmatched_families.discard(clone_id)
        trees_matched += 1

        # Build CSV index keyed by (site, parent_aa, child_aa)
        csv_index: Dict[Tuple[Any, str, str], Dict[str, Any]] = {}
        for row in mutations_by_family[clone_id]:
            site = row.get("site")
            paa = row.get("parent_aa")
            caa = row.get("child_aa")
            if site is None or paa is None or caa is None:
                continue
            extras = {
                k: v
                for k, v in row.items()
                if k not in ("site", "parent_aa", "child_aa")
            }
            csv_index[(site, paa, caa)] = extras

        nodes = _normalize_nodes(tree)
        parent_lookup = _build_parent_lookup(nodes)

        for node in nodes:
            if not isinstance(node, dict):
                continue

            existing = node.get("mutations")
            if isinstance(existing, list) and existing:
                node_mutations = existing
            else:
                parent_id = node.get("parent")
                parent_node = parent_lookup.get(parent_id) if parent_id else None
                node_mutations = derive_node_mutations(node, parent_node)
                if not node_mutations:
                    continue

            any_merged = False
            for mut in node_mutations:
                key = (mut.get("site"), mut.get("parent_aa"), mut.get("child_aa"))
                extras = csv_index.get(key)
                if extras:
                    mut.update(extras)
                    mutations_merged += 1
                    any_merged = True

            # Write back if we derived (or already had) mutations
            node["mutations"] = sorted(node_mutations, key=lambda m: m.get("site", 0))
            if any_merged or existing:
                nodes_with_mutations += 1

    if unmatched_families:
        vprint.verbose(
            f"  Mutations CSV had {len(unmatched_families)} families with no matching clone "
            f"(e.g., {sorted(unmatched_families)[:3]})"
        )

    return trees_matched, nodes_with_mutations, mutations_merged
