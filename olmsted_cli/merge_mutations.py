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
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from .constants import MUTATIONS_CSV_KEY_COLUMNS
from .process_utils import coerce_csv_value
from .utils import vprint

# Characters in AA sequences that should not be treated as a mutation event.
_GAP_AA_CHARS = {"-", ".", "X", "*", "?"}


@dataclass
class MergeStats:
    """Stats from merging a mutations CSV into a set of trees."""

    trees_matched: int = 0
    nodes_with_mutations: int = 0
    mutations_merged: int = 0
    # Families present in the CSV that have no matching tree clone_id.
    unmatched_families: List[str] = field(default_factory=list)
    # Count of CSV rows whose (family, site, parent_aa, child_aa) key had no
    # matching derived mutation in the corresponding tree.
    unmatched_mutations: int = 0


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
                if key == "site":
                    # site must be an integer index to match derived mutations
                    try:
                        mutation["site"] = int(val)
                    except (ValueError, TypeError) as e:
                        raise ValueError(
                            f"Mutations CSV row {reader.line_num}: 'site' must be "
                            f"an integer, got {val!r}"
                        ) from e
                elif key in ("parent_aa", "child_aa"):
                    mutation[key] = val
                elif key in MUTATIONS_CSV_KEY_COLUMNS:
                    # Skip other key/structural columns (family, sample_id, pcp_index, depth)
                    continue
                else:
                    mutation[key] = coerce_csv_value(val)

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
) -> MergeStats:
    """Merge mutation-level CSV data into tree nodes (in place).

    For each tree whose ``clone_id`` matches a family in the CSV:
      1. Index the family's CSV rows by (site, parent_aa, child_aa)
      2. For each node, derive (or read existing) mutation records
      3. Match each mutation against the CSV index and merge any extra fields

    Args:
        trees: List of tree dicts (modified in place).
        mutations_by_family: Output of ``load_mutations_csv``.

    Returns:
        A ``MergeStats`` describing what was merged and what went unmatched.
    """
    stats = MergeStats()
    unmatched_family_set = set(mutations_by_family.keys())

    for tree in trees:
        clone_id = tree.get("clone_id")
        if not clone_id or clone_id not in mutations_by_family:
            continue

        unmatched_family_set.discard(clone_id)
        stats.trees_matched += 1

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

        # Track which CSV keys were actually used — remaining keys are
        # unmatched mutations for this family.
        matched_keys: set = set()

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
                    stats.mutations_merged += 1
                    matched_keys.add(key)
                    any_merged = True

            # Write back if we derived (or already had) mutations
            node["mutations"] = sorted(node_mutations, key=lambda m: m.get("site", 0))
            if any_merged or existing:
                stats.nodes_with_mutations += 1

        # Account for CSV rows that never matched a derived mutation in this tree
        unmatched_in_family = [k for k in csv_index if k not in matched_keys]
        if unmatched_in_family:
            stats.unmatched_mutations += len(unmatched_in_family)
            sample = unmatched_in_family[:3]
            vprint.verbose(
                f"  {clone_id}: {len(unmatched_in_family)} CSV mutations had no "
                f"matching node (e.g., {sample})"
            )

    stats.unmatched_families = sorted(unmatched_family_set)
    return stats


def apply_mutations_csv(
    mutations_path: Optional[str],
    datasets: List[Dict[str, Any]],
    clones_dict: Dict[str, List[Dict[str, Any]]],
    trees: List[Dict[str, Any]],
    custom_fields: Optional[List[Dict[str, Any]]] = None,
) -> Optional[MergeStats]:
    """High-level entry point: load a mutations CSV, merge, warn, retag.

    Used by both the ``merge`` command and ``process --mutations``. Modifies
    ``datasets`` and ``trees`` in place. Returns ``None`` if ``mutations_path``
    is falsy (no-op), otherwise the ``MergeStats``.

    Warnings for unmatched families and unmatched mutations are emitted via
    ``vprint.error`` at normal verbosity. Per-family detail is at -v 2.
    """
    if not mutations_path:
        return None

    # Import here to avoid a circular import: process_utils imports from
    # field_metadata which is pulled in via tag_field_metadata.
    from .process_utils import retag_datasets_field_metadata

    vprint.status(f"Loading mutations CSV: {mutations_path}")
    mutations_by_family = load_mutations_csv(mutations_path)
    total = sum(len(rows) for rows in mutations_by_family.values())
    vprint.status(
        f"Loaded {total} mutation records across {len(mutations_by_family)} families"
    )

    stats = merge_mutations_into_trees(trees, mutations_by_family)
    vprint.status(
        f"Merged {stats.mutations_merged} mutation records into "
        f"{stats.nodes_with_mutations} nodes across {stats.trees_matched} trees"
    )

    if stats.trees_matched == 0:
        vprint.error(
            "Warning: No trees matched the families in the mutations CSV. "
            "Check that the CSV 'family' column matches clone_id values."
        )
    if stats.unmatched_families:
        sample = stats.unmatched_families[:5]
        vprint.error(
            f"Warning: {len(stats.unmatched_families)} families in the mutations CSV "
            f"had no matching clone (e.g., {sample})"
        )
    if stats.unmatched_mutations:
        vprint.error(
            f"Warning: {stats.unmatched_mutations} CSV mutation records in matched "
            f"families had no corresponding derived mutation in any node. "
            f"Run with -v 2 to see per-family details."
        )

    retag_datasets_field_metadata(
        datasets, clones_dict, trees, custom_fields=custom_fields
    )
    return stats
