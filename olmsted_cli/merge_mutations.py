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
from collections import deque
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from .constants import MUTATIONS_CSV_KEY_COLUMNS
from .process_utils import coerce_csv_value
from .utils import vprint

# Characters in AA sequences that should not be treated as a mutation event.
_GAP_AA_CHARS = {"-", ".", "X", "*", "?"}

# Optional CSV columns that, when present, narrow node-mutation matching beyond
# the (site, parent_aa, child_aa) base key. Each column also has to be derivable
# (or already present) on tree nodes for the disambiguation to work.
_DISAMBIGUATION_COLUMNS = ("depth",)


@dataclass
class MergeStats:
    """Stats from merging a mutations CSV into a set of trees.

    Counts are per-merge-run: ``nodes_enriched`` and ``mutations_enriched``
    only count nodes/mutations that received CSV-sourced data on this run.
    Pre-existing mutation arrays from upstream pipelines are *not* counted.
    """

    # Tree-side counts: what was actually changed in the JSON this run.
    trees_matched: int = 0
    nodes_enriched: int = 0
    mutations_enriched: int = 0

    # CSV-side counts: what didn't make it through.
    # Families present in the CSV that have no matching tree clone_id.
    unmatched_families: List[str] = field(default_factory=list)
    # Total CSV rows belonging to unmatched families.
    unmatched_family_rows: int = 0
    # Count of CSV rows whose match key had no corresponding derived mutation
    # in the matched tree (within families that did match a tree).
    unmatched_mutations: int = 0
    # Count of CSV rows that matched more than one node-mutation pair. This
    # indicates a potentially-ambiguous join: the same residue substitution
    # at the same site (and depth, if present) occurred on multiple lineages.
    # The same enrichment data is broadcast to every matching node.
    broadcast_csv_rows: int = 0
    # Names of optional disambiguation columns that were present in the CSV
    # and used as part of the match key.
    disambiguation_columns_used: List[str] = field(default_factory=list)


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
                elif key == "depth":
                    # depth, when present, is used as a disambiguation key —
                    # must be an integer (edges from root)
                    try:
                        mutation["depth"] = int(val)
                    except (ValueError, TypeError) as e:
                        raise ValueError(
                            f"Mutations CSV row {reader.line_num}: 'depth' must "
                            f"be an integer, got {val!r}"
                        ) from e
                elif key in ("parent_aa", "child_aa"):
                    mutation[key] = val
                elif key in MUTATIONS_CSV_KEY_COLUMNS:
                    # Skip other key/structural columns (family, sample_id, pcp_index)
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


def _compute_node_depths(nodes: List[Dict[str, Any]]) -> Dict[str, int]:
    """Compute each node's depth (edges from the nearest root) via BFS.

    A node is treated as a root if its ``parent`` field is missing/None or
    refers to a sequence_id not in the tree.
    """
    sid_set = {n["sequence_id"] for n in nodes if isinstance(n, dict) and "sequence_id" in n}
    children_by_parent: Dict[str, List[str]] = {}
    for n in nodes:
        if not isinstance(n, dict):
            continue
        sid = n.get("sequence_id")
        if not sid:
            continue
        parent = n.get("parent")
        if parent and parent in sid_set:
            children_by_parent.setdefault(parent, []).append(sid)

    depths: Dict[str, int] = {}
    queue: deque = deque()
    for n in nodes:
        if not isinstance(n, dict):
            continue
        sid = n.get("sequence_id")
        if not sid:
            continue
        parent = n.get("parent")
        if not parent or parent not in sid_set:
            depths[sid] = 0
            queue.append(sid)

    while queue:
        sid = queue.popleft()
        for child_sid in children_by_parent.get(sid, []):
            if child_sid not in depths:
                depths[child_sid] = depths[sid] + 1
                queue.append(child_sid)
    return depths


def _detect_disambiguation_columns(
    mutations_by_family: Dict[str, List[Dict[str, Any]]],
) -> List[str]:
    """Return the list of optional disambiguation columns present in the CSV.

    A column counts as "present" if at least one loaded row carries it. We
    take the union across families because the loader populates fields
    sparsely (empty cells are skipped).
    """
    present = set()
    for rows in mutations_by_family.values():
        for row in rows:
            for col in _DISAMBIGUATION_COLUMNS:
                if col in row:
                    present.add(col)
        if len(present) == len(_DISAMBIGUATION_COLUMNS):
            break
    return [c for c in _DISAMBIGUATION_COLUMNS if c in present]


def merge_mutations_into_trees(
    trees: List[Dict[str, Any]],
    mutations_by_family: Dict[str, List[Dict[str, Any]]],
) -> MergeStats:
    """Merge mutation-level CSV data into tree nodes (in place).

    For each tree whose ``clone_id`` matches a family in the CSV:
      1. Detect optional disambiguation columns (e.g. ``depth``) in the CSV
      2. Index the family's CSV rows by (site, parent_aa, child_aa[, depth])
      3. For each node, derive (or read existing) mutation records
      4. Match each mutation against the CSV index and merge any extra fields

    When optional disambiguation columns are present, they tighten the join
    key and reduce false fan-out. The number of CSV rows that still match
    multiple nodes is tracked in ``stats.broadcast_csv_rows``.

    Args:
        trees: List of tree dicts (modified in place).
        mutations_by_family: Output of ``load_mutations_csv``.

    Returns:
        A ``MergeStats`` describing what was merged and what went unmatched.
    """
    stats = MergeStats()
    stats.disambiguation_columns_used = _detect_disambiguation_columns(mutations_by_family)
    use_depth = "depth" in stats.disambiguation_columns_used
    unmatched_family_set = set(mutations_by_family.keys())

    # Columns that are part of the join key and should NOT be enriched onto
    # the output mutation record.
    extras_excluded = {"site", "parent_aa", "child_aa", *stats.disambiguation_columns_used}

    for tree in trees:
        clone_id = tree.get("clone_id")
        if not clone_id or clone_id not in mutations_by_family:
            continue

        unmatched_family_set.discard(clone_id)
        stats.trees_matched += 1

        nodes = _normalize_nodes(tree)
        parent_lookup = _build_parent_lookup(nodes)
        node_depths = _compute_node_depths(nodes) if use_depth else {}

        # Build CSV index keyed by (site, parent_aa, child_aa[, depth])
        csv_index: Dict[Tuple, Dict[str, Any]] = {}
        csv_match_counts: Dict[Tuple, int] = {}
        for row in mutations_by_family[clone_id]:
            site = row.get("site")
            paa = row.get("parent_aa")
            caa = row.get("child_aa")
            if site is None or paa is None or caa is None:
                continue
            if use_depth:
                depth = row.get("depth")
                if depth is None:
                    # CSV declares depth but this row is missing it; skip.
                    continue
                key: Tuple = (site, paa, caa, depth)
            else:
                key = (site, paa, caa)
            extras = {k: v for k, v in row.items() if k not in extras_excluded}
            csv_index[key] = extras
            csv_match_counts[key] = 0

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

            node_depth = node_depths.get(node.get("sequence_id")) if use_depth else None

            any_merged = False
            for mut in node_mutations:
                if use_depth:
                    key = (
                        mut.get("site"),
                        mut.get("parent_aa"),
                        mut.get("child_aa"),
                        node_depth,
                    )
                else:
                    key = (mut.get("site"), mut.get("parent_aa"), mut.get("child_aa"))
                extras = csv_index.get(key)
                if extras:
                    mut.update(extras)
                    stats.mutations_enriched += 1
                    csv_match_counts[key] = csv_match_counts.get(key, 0) + 1
                    any_merged = True

            # Write back if we derived (or already had) mutations
            node["mutations"] = sorted(node_mutations, key=lambda m: m.get("site", 0))
            if any_merged:
                stats.nodes_enriched += 1

        # Account for CSV rows that never matched a derived mutation in this tree
        unmatched_in_family = [k for k, c in csv_match_counts.items() if c == 0]
        if unmatched_in_family:
            stats.unmatched_mutations += len(unmatched_in_family)
            sample = unmatched_in_family[:3]
            vprint.verbose(
                f"  {clone_id}: {len(unmatched_in_family)} CSV mutations had no "
                f"matching node (e.g., {sample})"
            )

        # Account for CSV rows that broadcast to multiple nodes (ambiguous join)
        broadcast_in_family = [(k, c) for k, c in csv_match_counts.items() if c > 1]
        if broadcast_in_family:
            stats.broadcast_csv_rows += len(broadcast_in_family)
            sample = broadcast_in_family[:3]
            vprint.verbose(
                f"  {clone_id}: {len(broadcast_in_family)} CSV rows broadcast to "
                f"multiple nodes (e.g., {sample})"
            )

    stats.unmatched_families = sorted(unmatched_family_set)
    stats.unmatched_family_rows = sum(
        len(mutations_by_family[fam]) for fam in unmatched_family_set
    )
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
    total_csv_rows = sum(len(rows) for rows in mutations_by_family.values())
    vprint.status(
        f"Loaded {total_csv_rows} CSV rows across {len(mutations_by_family)} families"
    )

    stats = merge_mutations_into_trees(trees, mutations_by_family)

    if stats.disambiguation_columns_used:
        vprint.status(
            f"Disambiguation columns in CSV: "
            f"{', '.join(stats.disambiguation_columns_used)}"
        )

    vprint.status(
        f"Enriched {stats.mutations_enriched} mutations across "
        f"{stats.nodes_enriched} nodes in {stats.trees_matched} trees"
    )

    total_unmatched_rows = stats.unmatched_family_rows + stats.unmatched_mutations
    if total_unmatched_rows or stats.unmatched_families:
        vprint.status(
            f"Unmatched: {total_unmatched_rows}/{total_csv_rows} CSV rows "
            f"({stats.unmatched_family_rows} in {len(stats.unmatched_families)} "
            f"unmatched families, {stats.unmatched_mutations} with no node match)"
        )
    if stats.broadcast_csv_rows:
        vprint.status(
            f"Broadcast: {stats.broadcast_csv_rows} CSV rows matched multiple nodes "
            f"(ambiguous join — same data applied to every match)"
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
    if stats.broadcast_csv_rows:
        vprint.error(
            f"Warning: {stats.broadcast_csv_rows} CSV rows broadcast to multiple "
            f"node-mutations. The same enrichment data was applied to every "
            f"matching node, which may not be correct if the upstream pipeline "
            f"computed per-event scores. Run with -v 2 to see per-family details."
        )

    retag_datasets_field_metadata(
        datasets, clones_dict, trees, custom_fields=custom_fields
    )
    return stats
