"""
Shared utilities for merging external mutation-level CSV data into Olmsted JSON.

Used by:
- The `merge` command (post-hoc merging into existing Olmsted JSON)
- The `process` command's `--mutations` flag (during initial processing)

The CSV is expected to have at minimum these columns:
    family, site, parent_aa, child_aa
plus any number of score/annotation columns.

Matching strategy (in order of preference):
  1. If a node-name column is present (``node_name`` or ``child_name``), the
     join key is ``(node_name, site)`` — fully disambiguating, no fan-out
     possible. ``parent_aa``/``child_aa`` become integrity checks against the
     tree's derived mutation at that (node, site).
  2. Otherwise the join key is ``(site, parent_aa, child_aa)``, which may
     fan out to multiple nodes on convergent lineages (counted as broadcast).
  3. Depth narrowing (``--mutations-use-depth``): extends the fallback key
     to ``(site, parent_aa, child_aa, depth)``. Opt-in because depth
     arithmetic depends on upstream rooting conventions the CLI can't infer.
"""

import csv
import gzip
from collections import deque
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from .constants import MUTATIONS_CSV_KEY_COLUMNS, MUTATIONS_CSV_NAME_ALIASES
from .process_utils import coerce_csv_value
from .utils import vprint

# Characters in AA sequences that should not be treated as a mutation event.
_GAP_AA_CHARS = {"-", ".", "X", "*", "?"}


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
    # Only non-zero when name-based keying is NOT in use.
    broadcast_csv_rows: int = 0
    # Count of CSV rows whose (node_name, site) resolved to a real node + site
    # but where parent_aa/child_aa or depth disagreed with the tree's derived
    # mutation at that position. The CSV data is NOT attached on mismatch,
    # and (under --mutations-strict-check) the command exits non-zero.
    integrity_mismatches: int = 0
    # Names of optional disambiguation / structural columns that were present
    # in the CSV and used as part of the match key or as integrity checks.
    disambiguation_columns_used: List[str] = field(default_factory=list)
    # The chosen matching mode: "name_site" (deterministic) or
    # "site_paa_caa[_depth]" (may broadcast).
    match_mode: str = ""


def _detect_name_column(fieldnames: List[str]) -> Optional[str]:
    """Return the first node-name alias present in the CSV header, or None.

    Searches ``MUTATIONS_CSV_NAME_ALIASES`` in order — ``node_name`` wins
    over ``child_name`` if both are present.
    """
    for alias in MUTATIONS_CSV_NAME_ALIASES:
        if alias in fieldnames:
            return alias
    return None


def load_mutations_csv(csv_path: str) -> Dict[str, List[Dict[str, Any]]]:
    """Load a mutations CSV and group rows by family.

    Args:
        csv_path: Path to the mutations CSV (may be gzipped).

    Returns:
        Dict mapping family/clone_id -> list of row dicts. Each row dict
        contains the mutation keys (site, parent_aa, child_aa; optionally
        node_name and/or depth) plus all score/annotation columns.
        Structural join columns not needed downstream (family, sample_id,
        pcp_index) are stripped. If the CSV has a ``child_name`` column,
        values are normalized onto the ``node_name`` key for downstream use.
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

        name_col = _detect_name_column(list(reader.fieldnames))

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
                    # depth, when present, may be used as a disambiguation key
                    # or integrity check — must be an integer.
                    try:
                        mutation["depth"] = int(val)
                    except (ValueError, TypeError) as e:
                        raise ValueError(
                            f"Mutations CSV row {reader.line_num}: 'depth' must "
                            f"be an integer, got {val!r}"
                        ) from e
                elif key in ("parent_aa", "child_aa"):
                    mutation[key] = val
                elif key == name_col:
                    # Normalize onto the canonical "node_name" key so downstream
                    # code doesn't need to care which alias the CSV used.
                    mutation["node_name"] = val
                elif key in MUTATIONS_CSV_KEY_COLUMNS:
                    # Skip other key/structural columns (family, sample_id,
                    # pcp_index, and the losing name alias if both are present).
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


def _compute_node_depths(
    nodes: List[Dict[str, Any]], naive_name: str = "naive"
) -> Dict[str, int]:
    """Compute each node's depth via BFS, measuring from ``naive`` when present.

    If a node with ``sequence_id == naive_name`` exists, depths are measured
    as the **undirected** graph distance from that node (treating parent→child
    edges as bidirectional). This matches the depth convention used by
    upstream mutation/surprise pipelines that treat the naive/germline
    sequence as the root reference point — even though olmsted trees
    typically place naive as a sibling leaf next to a synthetic root.

    Falls back to directed BFS from the tree root (nodes with no parent) when
    no naive node is found.
    """
    sid_set = {n["sequence_id"] for n in nodes if isinstance(n, dict) and "sequence_id" in n}

    # Build undirected adjacency for naive-rooted BFS
    adj: Dict[str, List[str]] = {}
    for n in nodes:
        if not isinstance(n, dict):
            continue
        sid = n.get("sequence_id")
        parent = n.get("parent")
        if sid and parent and parent in sid_set:
            adj.setdefault(sid, []).append(parent)
            adj.setdefault(parent, []).append(sid)

    # Prefer naive as the reference origin (depth 0), but only if it's
    # connected to the rest of the tree. In some tree representations,
    # naive is an isolated root in a forest of subtrees — in that case
    # BFS from naive can't reach anything useful, so we fall back.
    if naive_name in sid_set and adj.get(naive_name):
        depths: Dict[str, int] = {naive_name: 0}
        queue: deque = deque([naive_name])
        while queue:
            sid = queue.popleft()
            for nb in adj.get(sid, []):
                if nb not in depths:
                    depths[nb] = depths[sid] + 1
                    queue.append(nb)
        if len(depths) > 1:
            return depths
        # naive is present but disconnected — fall through to root-based BFS

    # Fallback: directed BFS from tree root(s)
    depths = {}
    queue = deque()
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

    children_by_parent: Dict[str, List[str]] = {}
    for n in nodes:
        if not isinstance(n, dict):
            continue
        sid = n.get("sequence_id")
        parent = n.get("parent")
        if sid and parent and parent in sid_set:
            children_by_parent.setdefault(parent, []).append(sid)

    while queue:
        sid = queue.popleft()
        for child_sid in children_by_parent.get(sid, []):
            if child_sid not in depths:
                depths[child_sid] = depths[sid] + 1
                queue.append(child_sid)
    return depths


def _has_name_column(mutations_by_family: Dict[str, List[Dict[str, Any]]]) -> bool:
    """Return True if any loaded row carries a canonical node_name value."""
    for rows in mutations_by_family.values():
        for row in rows:
            if "node_name" in row:
                return True
    return False


def _has_depth_column(mutations_by_family: Dict[str, List[Dict[str, Any]]]) -> bool:
    """Return True if any loaded row carries a depth value."""
    for rows in mutations_by_family.values():
        for row in rows:
            if "depth" in row:
                return True
    return False


def merge_mutations_into_trees(
    trees: List[Dict[str, Any]],
    mutations_by_family: Dict[str, List[Dict[str, Any]]],
    *,
    use_depth: bool = False,
) -> MergeStats:
    """Merge mutation-level CSV data into tree nodes (in place).

    For each tree whose ``clone_id`` matches a family in the CSV:
      1. Decide match mode:
         - ``name_site``: CSV has a node-name column → deterministic
           ``(node_name, site)`` join with ``parent_aa``/``child_aa``/``depth``
           (if present) as integrity checks.
         - ``site_paa_caa``: fallback; ``(site, parent_aa, child_aa)`` +
           optionally ``depth`` if ``use_depth=True``.
      2. For each node, derive (or read existing) mutation records and
         match against the family's CSV rows.
      3. Attach extra fields to matching mutations. Integrity mismatches
         (in name_site mode) warn and skip the attachment.

    Args:
        trees: List of tree dicts (modified in place).
        mutations_by_family: Output of ``load_mutations_csv``.
        use_depth: When True, and the CSV has a ``depth`` column, use it to
            extend the fallback join key. Ignored when a name column is
            present (depth is an integrity check in that mode).

    Returns:
        A ``MergeStats`` describing what was merged, what was skipped, and
        the chosen match mode.
    """
    stats = MergeStats()
    name_keyed = _has_name_column(mutations_by_family)
    depth_present = _has_depth_column(mutations_by_family)

    if name_keyed:
        stats.match_mode = "name_site"
        extend_with_depth = False  # depth is integrity-only in this mode
    elif use_depth and depth_present:
        stats.match_mode = "site_paa_caa_depth"
        extend_with_depth = True
    else:
        stats.match_mode = "site_paa_caa"
        extend_with_depth = False

    # Columns used for matching or integrity; surface them in stats/logs.
    disambig = []
    if name_keyed:
        disambig.append("node_name")
    if depth_present and (extend_with_depth or name_keyed):
        disambig.append("depth")
    stats.disambiguation_columns_used = disambig

    unmatched_family_set = set(mutations_by_family.keys())

    # Columns that are part of the join key / integrity, not payload.
    extras_excluded = {"site", "parent_aa", "child_aa", "node_name", "depth"}

    if name_keyed:
        _merge_name_keyed(
            trees, mutations_by_family, stats, unmatched_family_set, extras_excluded,
            check_depth=depth_present,
        )
    else:
        _merge_site_keyed(
            trees, mutations_by_family, stats, unmatched_family_set, extras_excluded,
            use_depth=extend_with_depth,
        )

    stats.unmatched_families = sorted(unmatched_family_set)
    stats.unmatched_family_rows = sum(
        len(mutations_by_family[fam]) for fam in unmatched_family_set
    )
    return stats


def _merge_name_keyed(
    trees: List[Dict[str, Any]],
    mutations_by_family: Dict[str, List[Dict[str, Any]]],
    stats: MergeStats,
    unmatched_family_set: set,
    extras_excluded: set,
    check_depth: bool,
) -> None:
    """Deterministic ``(node_name, site)`` merge with integrity checks."""
    for tree in trees:
        clone_id = tree.get("clone_id")
        if not clone_id or clone_id not in mutations_by_family:
            continue

        unmatched_family_set.discard(clone_id)
        stats.trees_matched += 1

        nodes = _normalize_nodes(tree)
        nodes_by_id = _build_parent_lookup(nodes)
        parent_lookup = nodes_by_id  # same dict; different semantic name
        node_depths = _compute_node_depths(nodes) if check_depth else {}

        unmatched = 0
        mismatched = 0
        enriched_nodes: set = set()

        for row in mutations_by_family[clone_id]:
            site = row.get("site")
            paa = row.get("parent_aa")
            caa = row.get("child_aa")
            name = row.get("node_name")
            if site is None or paa is None or caa is None or name is None:
                unmatched += 1
                continue

            node = nodes_by_id.get(name)
            if node is None:
                unmatched += 1
                continue

            # Derive or read the node's mutations
            existing = node.get("mutations")
            if isinstance(existing, list) and existing:
                node_mutations = existing
            else:
                parent_id = node.get("parent")
                parent_node = parent_lookup.get(parent_id) if parent_id else None
                node_mutations = derive_node_mutations(node, parent_node)
                node["mutations"] = node_mutations

            # Locate the mutation at this site
            target = next((m for m in node_mutations if m.get("site") == site), None)
            if target is None:
                unmatched += 1
                continue

            # Integrity: parent_aa / child_aa must agree with the derived mutation
            if target.get("parent_aa") != paa or target.get("child_aa") != caa:
                mismatched += 1
                vprint.verbose(
                    f"  {clone_id}: integrity mismatch on ({name}, site {site}): "
                    f"CSV says {paa}->{caa}, tree has "
                    f"{target.get('parent_aa')}->{target.get('child_aa')}; skipping"
                )
                continue

            # Integrity: depth (when the CSV has it) must agree with the tree
            if check_depth and "depth" in row:
                tree_depth = node_depths.get(name)
                if tree_depth is not None and row["depth"] != tree_depth:
                    mismatched += 1
                    vprint.verbose(
                        f"  {clone_id}: depth mismatch on ({name}, site {site}): "
                        f"CSV depth {row['depth']}, tree depth {tree_depth}; skipping"
                    )
                    continue

            extras = {k: v for k, v in row.items() if k not in extras_excluded}
            target.update(extras)
            stats.mutations_enriched += 1
            enriched_nodes.add(name)

        # Keep mutations arrays sorted where we touched them
        for name in enriched_nodes:
            node = nodes_by_id[name]
            node["mutations"] = sorted(
                node.get("mutations", []), key=lambda m: m.get("site", 0)
            )

        stats.nodes_enriched += len(enriched_nodes)
        stats.unmatched_mutations += unmatched
        stats.integrity_mismatches += mismatched

        if unmatched:
            vprint.verbose(
                f"  {clone_id}: {unmatched} CSV rows had no matching (node, site)"
            )
        if mismatched:
            vprint.verbose(
                f"  {clone_id}: {mismatched} CSV rows matched (node, site) but "
                f"failed integrity check"
            )


def _merge_site_keyed(
    trees: List[Dict[str, Any]],
    mutations_by_family: Dict[str, List[Dict[str, Any]]],
    stats: MergeStats,
    unmatched_family_set: set,
    extras_excluded: set,
    use_depth: bool,
) -> None:
    """Fallback ``(site, parent_aa, child_aa[, depth])`` merge. May broadcast."""
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
                    continue  # CSV declares depth but this row is missing it
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

            node["mutations"] = sorted(node_mutations, key=lambda m: m.get("site", 0))
            if any_merged:
                stats.nodes_enriched += 1

        unmatched_in_family = [k for k, c in csv_match_counts.items() if c == 0]
        if unmatched_in_family:
            stats.unmatched_mutations += len(unmatched_in_family)
            sample = unmatched_in_family[:3]
            vprint.verbose(
                f"  {clone_id}: {len(unmatched_in_family)} CSV mutations had no "
                f"matching node (e.g., {sample})"
            )

        broadcast_in_family = [(k, c) for k, c in csv_match_counts.items() if c > 1]
        if broadcast_in_family:
            stats.broadcast_csv_rows += len(broadcast_in_family)
            sample = broadcast_in_family[:3]
            vprint.verbose(
                f"  {clone_id}: {len(broadcast_in_family)} CSV rows broadcast to "
                f"multiple nodes (e.g., {sample})"
            )


def apply_mutations_csv(
    mutations_path: Optional[str],
    datasets: List[Dict[str, Any]],
    clones_dict: Dict[str, List[Dict[str, Any]]],
    trees: List[Dict[str, Any]],
    custom_fields: Optional[List[Dict[str, Any]]] = None,
    *,
    use_depth: bool = False,
    strict_check: bool = False,
) -> Optional[MergeStats]:
    """High-level entry point: load a mutations CSV, merge, warn, retag.

    Used by both the ``merge`` command and ``process --mutations``. Modifies
    ``datasets`` and ``trees`` in place. Returns ``None`` if ``mutations_path``
    is falsy (no-op), otherwise the ``MergeStats``.

    Args:
        use_depth: Enable the optional ``depth`` column as part of the
            fallback match key. No effect when the CSV has a node-name column.
        strict_check: Treat integrity mismatches (``parent_aa``/``child_aa``
            or ``depth`` disagreement when a row matches by ``(node_name,
            site)``) as a hard error; the caller is responsible for exiting.

    Warnings for unmatched families, unmatched mutations, and integrity
    mismatches are emitted via ``vprint.error`` at normal verbosity.
    Per-family detail is at -v 2.

    Raises:
        ValueError: When ``strict_check=True`` and at least one integrity
            mismatch was recorded. The merge has still been applied to
            the rows that matched cleanly.
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

    stats = merge_mutations_into_trees(
        trees, mutations_by_family, use_depth=use_depth
    )

    vprint.status(f"Match mode: {stats.match_mode}")
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
    if stats.integrity_mismatches:
        vprint.status(
            f"Integrity mismatches: {stats.integrity_mismatches} CSV rows matched "
            f"a (node, site) but disagreed with the tree's derived mutation"
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
    if stats.integrity_mismatches:
        vprint.error(
            f"Warning: {stats.integrity_mismatches} CSV rows resolved to a real "
            f"(node, site) but had parent_aa/child_aa/depth that disagreed with "
            f"the tree's derived mutation. The enrichment data was NOT attached "
            f"for those rows. Run with -v 2 to see per-family details."
        )

    retag_datasets_field_metadata(
        datasets, clones_dict, trees, custom_fields=custom_fields
    )

    if strict_check and stats.integrity_mismatches:
        raise ValueError(
            f"--mutations-strict-check: {stats.integrity_mismatches} integrity "
            f"mismatches between CSV rows and tree mutations. Re-run without "
            f"--mutations-strict-check to continue despite the mismatches."
        )

    return stats
