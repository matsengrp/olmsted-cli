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
from collections import deque
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Literal, Optional, Set, Tuple

from .constants import MUTATIONS_CSV_KEY_COLUMNS, MUTATIONS_CSV_NAME_ALIASES
from .data_io import open_file
from .process_utils import coerce_csv_value
from .utils import vprint

# Characters in AA sequences that should not be treated as a mutation event.
_GAP_AA_CHARS = {"-", ".", "X", "*", "?"}

# Match-mode identifiers for MergeStats.match_mode. The empty string is the
# pre-merge sentinel (before merge_mutations_into_trees has run).
MatchMode = Literal["", "name_site", "site_paa_caa", "site_paa_caa_depth"]
MATCH_MODE_NAME_SITE: MatchMode = "name_site"
MATCH_MODE_SITE_PAA_CAA: MatchMode = "site_paa_caa"
MATCH_MODE_SITE_PAA_CAA_DEPTH: MatchMode = "site_paa_caa_depth"


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
    # mutation at that position. The CSV data is NOT attached on mismatch.
    # By default any non-zero value causes the command to exit non-zero;
    # --mutations-allow-mismatch downgrades this to a warning.
    integrity_mismatches: int = 0
    # Names of optional disambiguation / structural columns that were present
    # in the CSV and used as part of the match key or as integrity checks.
    disambiguation_columns_used: List[str] = field(default_factory=list)
    # The chosen matching mode: name_site (deterministic) or
    # site_paa_caa[_depth] (may broadcast).
    match_mode: MatchMode = ""
    # Mutations dropped under --mutations-listed-only: derived mutations on
    # nodes in CSV-matched trees that didn't appear in the CSV. Always 0
    # when only_listed=False.
    mutations_dropped: int = 0


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

    mutations_by_family: Dict[str, List[Dict[str, Any]]] = {}

    # Mutations CSVs aren't auto-detectable as a known Olmsted format; format
    # detection sees the .csv extension and labels it 'pcp', which is a fine
    # passthrough — caller doesn't act on the detected format here.
    handle, _ = open_file(path)
    with handle as file_handle:
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
                raise ValueError(
                    f"Mutations CSV row {reader.line_num}: 'family' is empty. "
                    f"Every row must have a family/clone_id for routing."
                )

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


def _build_node_id_lookup(nodes: List[Dict[str, Any]]) -> Dict[str, Dict[str, Any]]:
    """Build a sequence_id -> node lookup used for both identity and parent resolution."""
    return {
        n["sequence_id"]: n for n in nodes if isinstance(n, dict) and "sequence_id" in n
    }


def _get_or_derive_mutations(
    node: Dict[str, Any],
    nodes_by_id: Dict[str, Dict[str, Any]],
) -> List[Dict[str, Any]]:
    """Return the node's mutations, deriving from parent-AA diff if needed.

    Prefers a pre-existing non-empty ``mutations`` list (upstream pipeline
    data). Otherwise derives via ``derive_node_mutations`` and writes the
    result back onto the node if non-empty. Empty results are NOT persisted
    onto the node to avoid polluting it with empty arrays.
    """
    existing = node.get("mutations")
    if isinstance(existing, list) and existing:
        return existing
    parent_id = node.get("parent")
    parent_node = nodes_by_id.get(parent_id) if parent_id else None
    derived = derive_node_mutations(node, parent_node)
    if derived:
        node["mutations"] = derived
    return derived


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
    sid_set = {
        n["sequence_id"] for n in nodes if isinstance(n, dict) and "sequence_id" in n
    }

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


@dataclass
class MergeContext:
    """Per-CSV state threaded across one or many ``apply_mutations_to_trees`` calls.

    The merge command and one-shot PCP/AIRR runs use the legacy
    ``merge_mutations_into_trees`` wrapper.  The streaming PCP/AIRR
    pipelines invoke the merge once per batch and share a single context
    so per-tree counts, integrity-mismatch tallies, and the running
    ``unmatched_family_set`` aggregate across batches the way a one-shot
    run would aggregate across trees.
    """

    mutations_by_family: Dict[str, List[Dict[str, Any]]]
    stats: MergeStats
    unmatched_family_set: Set[str]
    name_keyed: bool
    honor_depth: bool
    extend_with_depth: bool
    extras_excluded: Set[str]


def begin_merge(
    mutations_by_family: Dict[str, List[Dict[str, Any]]],
    *,
    use_depth: bool = False,
) -> MergeContext:
    """Initialize a ``MergeContext`` from a loaded mutations index.

    Performs the up-front mode-decision and flag validation that
    ``merge_mutations_into_trees`` used to do inline: chooses
    name-keyed vs site-keyed mode, validates ``--mutations-use-depth``
    against the CSV's columns, populates ``stats.match_mode`` and
    ``stats.disambiguation_columns_used``.
    """
    stats = MergeStats()
    name_keyed = _has_name_column(mutations_by_family)
    depth_present = _has_depth_column(mutations_by_family)

    if use_depth and not depth_present:
        raise ValueError(
            "--mutations-use-depth was passed but the mutations CSV has no "
            "'depth' column. Either add depth values to the CSV or drop the flag."
        )

    honor_depth = use_depth and depth_present

    if name_keyed:
        stats.match_mode = MATCH_MODE_NAME_SITE
        extend_with_depth = False
    elif honor_depth:
        stats.match_mode = MATCH_MODE_SITE_PAA_CAA_DEPTH
        extend_with_depth = True
    else:
        stats.match_mode = MATCH_MODE_SITE_PAA_CAA
        extend_with_depth = False

    disambig = []
    if name_keyed:
        disambig.append("node_name")
    if honor_depth:
        disambig.append("depth")
    stats.disambiguation_columns_used = disambig

    if depth_present and not use_depth:
        vprint.verbose(
            "Mutations CSV has a 'depth' column but --mutations-use-depth was "
            "not set; ignoring depth for both match-key and integrity checks."
        )

    return MergeContext(
        mutations_by_family=mutations_by_family,
        stats=stats,
        unmatched_family_set=set(mutations_by_family.keys()),
        name_keyed=name_keyed,
        honor_depth=honor_depth,
        extend_with_depth=extend_with_depth,
        extras_excluded={"site", "parent_aa", "child_aa", "node_name", "depth"},
    )


def apply_mutations_to_trees(
    ctx: MergeContext,
    trees: List[Dict[str, Any]],
    *,
    only_listed: bool = False,
) -> None:
    """Apply the merge to a slice of trees, folding counts into ``ctx.stats``.

    Safe to call repeatedly: the context's stats / unmatched-family set
    aggregate across calls.  Used by the streaming pipeline to merge one
    batch at a time without re-loading the CSV or re-deciding the
    match mode.
    """
    if ctx.name_keyed:
        _merge_name_keyed(
            trees,
            ctx.mutations_by_family,
            ctx.stats,
            ctx.unmatched_family_set,
            ctx.extras_excluded,
            check_depth=ctx.honor_depth,
            only_listed=only_listed,
        )
    else:
        _merge_site_keyed(
            trees,
            ctx.mutations_by_family,
            ctx.stats,
            ctx.unmatched_family_set,
            ctx.extras_excluded,
            use_depth=ctx.extend_with_depth,
            only_listed=only_listed,
        )


def finalize_merge(ctx: MergeContext) -> MergeStats:
    """Seal the running stats — populates ``unmatched_families`` / rows."""
    ctx.stats.unmatched_families = sorted(ctx.unmatched_family_set)
    ctx.stats.unmatched_family_rows = sum(
        len(ctx.mutations_by_family[fam]) for fam in ctx.unmatched_family_set
    )
    return ctx.stats


def merge_mutations_into_trees(
    trees: List[Dict[str, Any]],
    mutations_by_family: Dict[str, List[Dict[str, Any]]],
    *,
    use_depth: bool = False,
    only_listed: bool = False,
) -> MergeStats:
    """Merge mutation-level CSV data into tree nodes (in place).

    For each tree whose ``clone_id`` matches a family in the CSV:
      1. Decide match mode:
         - ``name_site``: CSV has a node-name column → deterministic
           ``(node_name, site)`` join with ``parent_aa``/``child_aa`` as
           integrity checks (``depth`` too when ``use_depth=True``).
         - ``site_paa_caa``: fallback; ``(site, parent_aa, child_aa)`` +
           optionally ``depth`` if ``use_depth=True``.
      2. For each node, derive (or read existing) mutation records and
         match against the family's CSV rows.
      3. Attach extra fields to matching mutations. Integrity mismatches
         (in name_site mode) warn and skip the attachment.
      4. If ``only_listed`` is True, drop any mutations on those nodes
         that did not match a CSV row. Trees whose family is absent
         from the CSV are passed through untouched.

    Args:
        trees: List of tree dicts (modified in place).
        mutations_by_family: Output of ``load_mutations_csv``.
        use_depth: When True, the ``depth`` column (if present in the CSV)
            is honored — as a match-key participant in site-keyed mode,
            or as an integrity check in name-keyed mode. When False, a
            ``depth`` column in the CSV is ignored entirely for both
            purposes. Opt-in because depth arithmetic depends on the
            upstream pipeline's rooting convention, which the CLI cannot
            infer with certainty.
        only_listed: When True, the CSV is treated as authoritative for
            which mutations should appear on each node of a CSV-matched
            tree. Sequence-diff-derived (or pre-existing upstream)
            mutations that don't have a corresponding CSV row are
            removed. Trees whose ``clone_id`` does not appear in the
            CSV are not filtered.

    Returns:
        A ``MergeStats`` describing what was merged, what was skipped, and
        the chosen match mode.
    """
    ctx = begin_merge(mutations_by_family, use_depth=use_depth)
    apply_mutations_to_trees(ctx, trees, only_listed=only_listed)
    return finalize_merge(ctx)


def _merge_name_keyed(
    trees: List[Dict[str, Any]],
    mutations_by_family: Dict[str, List[Dict[str, Any]]],
    stats: MergeStats,
    unmatched_family_set: set,
    extras_excluded: set,
    check_depth: bool,
    only_listed: bool = False,
) -> None:
    """Deterministic ``(node_name, site)`` merge with integrity checks."""
    for tree in trees:
        clone_id = tree.get("clone_id")
        if not clone_id or clone_id not in mutations_by_family:
            continue

        unmatched_family_set.discard(clone_id)
        stats.trees_matched += 1

        nodes = _normalize_nodes(tree)
        nodes_by_id = _build_node_id_lookup(nodes)
        node_depths = _compute_node_depths(nodes) if check_depth else {}

        unmatched = 0
        mismatched = 0
        # node_name -> set of sites enriched by the CSV. Drives both the
        # only-listed sweep and the per-tree node-count stat.
        enriched_sites_by_node: Dict[str, set] = {}

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

            node_mutations = _get_or_derive_mutations(node, nodes_by_id)
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
            enriched_sites_by_node.setdefault(name, set()).add(site)

        # Under --only-listed, drop mutations on every node of this matched
        # tree that did not receive a CSV row. Includes nodes never enriched.
        if only_listed:
            for node in nodes:
                if not isinstance(node, dict):
                    continue
                muts = node.get("mutations")
                if not isinstance(muts, list) or not muts:
                    continue
                node_name = node.get("sequence_id")
                kept_sites = enriched_sites_by_node.get(node_name, set())
                kept = [m for m in muts if m.get("site") in kept_sites]
                stats.mutations_dropped += len(muts) - len(kept)
                if kept:
                    node["mutations"] = kept
                else:
                    # Drop the empty list rather than leaving [] behind, so
                    # nodes look the same as ones that never had mutations.
                    del node["mutations"]

        # Keep mutations arrays sorted where we touched them. Nodes that
        # received at least one CSV match always retain their `mutations`
        # array even after the only-listed sweep (the sweep only deletes
        # the key when *no* sites matched — i.e., nodes never recorded
        # in enriched_sites_by_node), so no guard is needed here.
        for name in enriched_sites_by_node:
            node = nodes_by_id[name]
            node["mutations"] = sorted(
                node["mutations"], key=lambda m: m.get("site", 0)
            )

        stats.nodes_enriched += len(enriched_sites_by_node)
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
    only_listed: bool = False,
) -> None:
    """Fallback ``(site, parent_aa, child_aa[, depth])`` merge. May broadcast."""
    for tree in trees:
        clone_id = tree.get("clone_id")
        if not clone_id or clone_id not in mutations_by_family:
            continue

        unmatched_family_set.discard(clone_id)
        stats.trees_matched += 1

        nodes = _normalize_nodes(tree)
        nodes_by_id = _build_node_id_lookup(nodes)
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

            node_mutations = _get_or_derive_mutations(node, nodes_by_id)
            if not node_mutations:
                continue

            node_depth = node_depths.get(node.get("sequence_id")) if use_depth else None

            any_merged = False
            kept_mutations: List[Dict[str, Any]] = []
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
                    kept_mutations.append(mut)
                elif only_listed:
                    stats.mutations_dropped += 1
                else:
                    kept_mutations.append(mut)

            if kept_mutations:
                node["mutations"] = sorted(
                    kept_mutations, key=lambda m: m.get("site", 0)
                )
            elif "mutations" in node:
                del node["mutations"]
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


def report_merge_stats(stats: MergeStats, total_csv_rows: int) -> None:
    """Emit the user-facing summary lines + warnings for a merge run.

    Shared between the one-shot ``apply_mutations_csv`` and the streaming
    pipelines so the on-screen output is identical regardless of whether
    the merge ran in one pass or across many batches.  The caller is
    responsible for raising on integrity mismatches when
    ``--mutations-allow-mismatch`` isn't set.
    """
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
    if stats.mutations_dropped:
        vprint.status(
            f"Dropped {stats.mutations_dropped} derived mutations not listed in "
            f"the CSV (--mutations-listed-only)"
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


def apply_mutations_csv(
    mutations_path: Optional[str],
    trees: List[Dict[str, Any]],
    *,
    use_depth: bool = False,
    allow_mismatch: bool = False,
    only_listed: bool = False,
) -> Optional[MergeStats]:
    """Load a mutations CSV and merge it into ``trees`` (in place).

    Used by both the ``merge`` command and ``process --mutations``. Returns
    ``None`` when ``mutations_path`` is falsy (no-op), otherwise the
    ``MergeStats`` from the merge.

    Callers are responsible for regenerating ``field_metadata`` (via
    ``retag_datasets_field_metadata``) afterwards. This is deliberate: the
    merge doesn't know about dataset-level metadata, and keeping the two
    steps separate keeps the function testable with just a list of trees.

    Args:
        use_depth: Enable the optional ``depth`` column as part of the
            fallback match key. No effect when the CSV has a node-name column.
        allow_mismatch: Downgrade integrity mismatches (``parent_aa``/
            ``child_aa`` or ``depth`` disagreement when a row matches by
            ``(node_name, site)``) from a hard failure to a warning. Default
            behavior is to raise ``ValueError`` if any mismatch occurs, so
            that callers can't accidentally ship a partially-wrong merge.
            Mismatched rows are never attached regardless of this flag.
        only_listed: Treat the CSV as authoritative — drop derived
            mutations on CSV-matched trees that don't appear in the CSV.
            See ``merge_mutations_into_trees`` for details.

    Warnings for unmatched families, unmatched mutations, and broadcast
    rows are emitted via ``vprint.error`` at normal verbosity.
    Per-family detail is at -v 2.

    Raises:
        ValueError: When at least one integrity mismatch was recorded and
            ``allow_mismatch`` is False (the default). The merge has still
            been applied to the rows that matched cleanly; mismatched rows
            are skipped.
    """
    if not mutations_path:
        return None

    vprint.status(f"Loading mutations CSV: {mutations_path}")
    mutations_by_family = load_mutations_csv(mutations_path)
    total_csv_rows = sum(len(rows) for rows in mutations_by_family.values())
    vprint.status(
        f"Loaded {total_csv_rows} CSV rows across {len(mutations_by_family)} families"
    )

    stats = merge_mutations_into_trees(
        trees, mutations_by_family, use_depth=use_depth, only_listed=only_listed
    )
    report_merge_stats(stats, total_csv_rows)

    if stats.integrity_mismatches and not allow_mismatch:
        raise ValueError(
            f"{stats.integrity_mismatches} integrity mismatches between CSV "
            f"rows and tree mutations. Mismatched rows were skipped "
            f"(never attached). Re-run with --mutations-allow-mismatch to "
            f"proceed anyway — but investigate the CSV/tree disagreement first."
        )

    return stats
