"""
Shared phylogenetic metric computations.

These functions compute tree-based metrics (LBI, LBR, scaled affinity)
from tree topology and branch lengths. They are format-agnostic and can
be used with data from any input format (PCP, AIRR, Olmsted JSON).
"""

from __future__ import annotations

import math
from collections import defaultdict
from typing import Any, Dict, Iterable, List, Optional, Tuple


def compute_lbi_for_tree(
    nodes_dict: Dict[str, Any],
    edges: List[Tuple[str, str, float]],
    root_id: str,
    tau: float = 0.0125,
) -> Dict[str, float]:
    """
    Compute Local Branching Index (LBI) for all nodes in a tree.

    LBI measures the local branching structure around each node, capturing
    the rate of diversification in the recent evolutionary history.

    Reference: Neher & Bedford (2015) "nextflu: real-time tracking of seasonal
    influenza virus evolution in humans" Bioinformatics 31(21):3546-3548

    Args:
        nodes_dict: Dictionary of {node_id: node_data} with branch length data.
        edges: List of (parent, child, length) tuples.
        root_id: ID of the root node.
        tau: Time scale parameter (default: 0.0125).

    Returns:
        dict: {node_id: lbi_value}
    """
    children_map = defaultdict(list)
    parent_map = {}
    edge_length_map = {}

    for parent, child, length in edges:
        children_map[parent].append(child)
        parent_map[child] = parent
        edge_length_map[(parent, child)] = length

    up_polarizer = {node_id: 0.0 for node_id in nodes_dict.keys()}
    down_polarizer = {node_id: 0.0 for node_id in nodes_dict.keys()}

    def postorder(node):
        if node not in children_map or len(children_map[node]) == 0:
            up_polarizer[node] = 0.0
            return

        for child in children_map[node]:
            postorder(child)

        total = 0.0
        for child in children_map[node]:
            branch_length = edge_length_map.get((node, child), 0.0)
            weight = math.exp(-branch_length / tau)
            total += (branch_length + up_polarizer[child]) * weight

        up_polarizer[node] = total

    def preorder(node):
        if node not in children_map:
            return

        for child in children_map[node]:
            branch_length = edge_length_map.get((node, child), 0.0)
            weight = math.exp(-branch_length / tau)

            parent_contribution = down_polarizer[node]

            sibling_contribution = 0.0
            for sibling in children_map[node]:
                if sibling != child:
                    sibling_branch_length = edge_length_map.get((node, sibling), 0.0)
                    sibling_contribution += (
                        sibling_branch_length + up_polarizer[sibling]
                    )

            down_polarizer[child] = (
                branch_length + parent_contribution + sibling_contribution
            ) * weight

            preorder(child)

    postorder(root_id)
    down_polarizer[root_id] = 0.0
    preorder(root_id)

    lbi = {}
    for node_id in nodes_dict.keys():
        lbi[node_id] = up_polarizer[node_id] + down_polarizer[node_id]

    return lbi


def compute_lbr_for_tree(
    nodes_dict: Dict[str, Any],
    edges: List[Tuple[str, str, float]],
    root_id: str,
) -> Dict[str, float]:
    """
    Compute Local Branching Ratio (LBR) for all nodes in a tree.

    LBR = log(downstream_branches / upstream_branches).

    Args:
        nodes_dict: Dictionary of {node_id: node_data}.
        edges: List of (parent, child, length) tuples.
        root_id: ID of the root node.

    Returns:
        dict: {node_id: lbr_value}
    """
    children_map = defaultdict(list)
    parent_map = {}

    for parent, child, length in edges:
        children_map[parent].append(child)
        parent_map[child] = parent

    downstream_count = {}

    def count_downstream(node):
        if node not in children_map or len(children_map[node]) == 0:
            downstream_count[node] = 0
            return 0

        total = 0
        for child in children_map[node]:
            total += 1
            total += count_downstream(child)

        downstream_count[node] = total
        return total

    upstream_count = {}

    def count_upstream(node):
        if node == root_id:
            upstream_count[node] = 0
            return 0

        count = 0
        current = node
        while current != root_id and current in parent_map:
            count += 1
            current = parent_map[current]

        upstream_count[node] = count
        return count

    count_downstream(root_id)
    for node_id in nodes_dict.keys():
        count_upstream(node_id)

    lbr = {}
    for node_id in nodes_dict.keys():
        down = downstream_count.get(node_id, 0)
        up = upstream_count.get(node_id, 0)

        if up == 0 or down == 0:
            lbr[node_id] = 0.0
        else:
            lbr[node_id] = math.log(down / up)

    return lbr


def compute_scaled_affinity(
    affinity_values: Dict[str, Optional[float]],
) -> Dict[str, Optional[float]]:
    """
    Compute scaled affinity using min-max normalization.

    scaled_affinity = (affinity - min) / (max - min)

    Args:
        affinity_values: dict of {node_id: affinity} (None values allowed).

    Returns:
        dict: {node_id: scaled_affinity}
    """
    valid_affinities = {k: v for k, v in affinity_values.items() if v is not None}

    if not valid_affinities:
        return {k: None for k in affinity_values.keys()}

    min_affinity = min(valid_affinities.values())
    max_affinity = max(valid_affinities.values())

    if max_affinity == min_affinity:
        return {k: (0.5 if v is not None else None) for k, v in affinity_values.items()}

    scaled = {}
    for node_id, affinity in affinity_values.items():
        if affinity is None:
            scaled[node_id] = None
        else:
            scaled[node_id] = (affinity - min_affinity) / (max_affinity - min_affinity)

    return scaled


def compute_tree_metrics(
    nodes_dict: Dict[str, Any],
    edges: List[Tuple[str, str, float]],
    root_id: str,
    tau: float = 0.0125,
) -> None:
    """
    Compute all standard phylogenetic metrics on a tree in place.

    Modifies nodes_dict to add lbi, lbr, affinity, and scaled_affinity.

    Args:
        nodes_dict: Dictionary of {node_id: node_data}. Modified in place.
        edges: List of (parent, child, length) tuples.
        root_id: ID of the root node.
        tau: Time scale parameter for LBI (default: 0.0125).
    """
    lbi_values = compute_lbi_for_tree(nodes_dict, edges, root_id, tau=tau)
    for nid in nodes_dict:
        lbi = lbi_values.get(nid)
        nodes_dict[nid]["lbi"] = lbi
        nodes_dict[nid]["affinity"] = lbi

    affinity_values = {nid: nodes_dict[nid].get("affinity") for nid in nodes_dict}
    scaled = compute_scaled_affinity(affinity_values)
    for nid in nodes_dict:
        nodes_dict[nid]["scaled_affinity"] = scaled.get(nid)

    lbr_values = compute_lbr_for_tree(nodes_dict, edges, root_id)
    for nid in nodes_dict:
        nodes_dict[nid]["lbr"] = lbr_values.get(nid)


# Positions where either sequence has one of these characters don't count as
# a comparable (mutation-eligible) position: "." is the pad/gap sentinel this
# module's alignment convention uses, "-" is the standard alignment-gap
# character, and "N" is an ambiguous/undetermined base. Shared by
# align_and_calculate_mutations and compute_mean_mut_freq so every ingest
# path treats "no information at this position" the same way.
NON_COMPARABLE = frozenset({"", ".", "-", "N"})


def align_and_calculate_mutations(
    germline: str, leaf: str, alignment_method: str = "truncate"
) -> Tuple[int, int, str, str]:
    """
    Align sequences and count mutations between germline and leaf.

    This function consolidates the shared logic between truncate and pad alignment methods,
    eliminating ~90% code duplication.

    Args:
        germline: Germline sequence string
        leaf: Leaf sequence string
        alignment_method: Either "truncate" or "pad"
            - "truncate": Use min length, truncate longer sequence
            - "pad": Use max length, pad shorter sequence with "."

    Returns:
        tuple: (mutation_count, comparable_length, germline_aligned, leaf_aligned)
            - mutation_count: Number of mismatches, counted only at comparable
              positions (see NON_COMPARABLE)
            - comparable_length: Number of positions where neither sequence
              has a non-comparable character. This is the intended
              denominator for a mutation frequency, not the raw aligned
              length — a position with an unknown base carries no
              information and shouldn't dilute the frequency.
            - germline_aligned: Aligned germline sequence
            - leaf_aligned: Aligned leaf sequence

    Examples:
        >>> align_and_calculate_mutations("ATGC", "ATGT", "truncate")
        (1, 4, 'ATGC', 'ATGT')
        >>> align_and_calculate_mutations("ATG", "ATGCC", "truncate")
        (0, 3, 'ATG', 'ATG')
        >>> align_and_calculate_mutations("ATG", "ATGCC", "pad")
        (0, 3, 'ATG..', 'ATGCC')
        >>> align_and_calculate_mutations("ATGN", "ATGC", "truncate")
        (0, 3, 'ATGN', 'ATGC')
    """
    if alignment_method == "truncate":
        # Truncate to shorter sequence length
        length = min(len(germline), len(leaf))
        g_aligned = germline[:length]
        l_aligned = leaf[:length]
    else:  # pad
        # Pad to longer sequence length
        length = max(len(germline), len(leaf))
        g_aligned = germline.ljust(length, ".")
        l_aligned = leaf.ljust(length, ".")

    # Count mutations and comparable positions together: a non-comparable
    # character (gap or ambiguous base) on either side excludes the position
    # from both the mismatch count and the denominator.
    mutations = 0
    comparable = 0
    for g, leaf_base in zip(g_aligned, l_aligned):
        if g in NON_COMPARABLE or leaf_base in NON_COMPARABLE:
            continue
        comparable += 1
        if g != leaf_base:
            mutations += 1

    return mutations, comparable, g_aligned, l_aligned


def compute_mean_mut_freq(
    germline_alignment: str,
    nodes: Iterable[Dict[str, Any]],
    alignment_method: str = "truncate",
) -> Tuple[float, List[Dict[str, Any]], List[Dict[str, Any]]]:
    """Multiplicity-weighted mean mutation frequency across leaf nodes.

    ``mean_mut_freq = sum_over_leaves(num_mutations / comparable_length * multiplicity)
    / sum_over_leaves(multiplicity)``

    Leaves are nodes with ``type == "leaf"`` and ``multiplicity > 0``;
    non-comparable positions (``.``, ``-``, ``N`` — see :data:`NON_COMPARABLE`)
    on either side of a mismatch are skipped, matching
    :func:`align_and_calculate_mutations` semantics. This is the single
    source of truth for ``mean_mut_freq`` across input formats — PCP and
    AIRR both route through this helper so the same biological data
    produces the same value regardless of which pipeline ingested it.

    Format-agnostic: ``nodes`` can be any iterable of node dicts (a PCP
    processed-nodes dict's ``.values()``, an AIRR tree's ``nodes`` list,
    etc.) that carry ``sequence_alignment``, ``multiplicity``, and
    ``type``. Non-leaf nodes and leaves with ``multiplicity == 0`` are
    recorded in ``skipped_nodes`` and otherwise ignored — the full node
    collection is passed in, not a pre-filtered leaf list.

    Args:
        germline_alignment: The clone's germline (naive) DNA sequence.
        nodes: Iterable of node dicts (not pre-filtered to leaves).
        alignment_method: "truncate" or "pad" (see
            :func:`align_and_calculate_mutations`).

    Returns:
        ``(mean_mut_freq, debug_info, skipped_nodes)``. ``mean_mut_freq``
        is ``0.0`` when there are no counted leaves (empty germline, no
        leaves, or all leaves have ``multiplicity == 0``) — never raises.
        ``debug_info`` has one record per counted leaf (node id,
        num_mutations, mut_freq, multiplicity, weighted_contribution,
        per-position mutations, ...) for callers that want to log it.
        ``skipped_nodes`` is the matching diagnostic list for nodes that
        didn't contribute, with a human-readable ``reason``.
    """
    total_mut_freq = 0.0
    total_multiplicity = 0
    debug_info: List[Dict[str, Any]] = []
    skipped_nodes: List[Dict[str, Any]] = []

    for node_data in nodes:
        node_id = node_data.get("sequence_id")
        node_type = node_data.get("type")
        multiplicity = node_data.get("multiplicity", 0)

        if node_type == "leaf" and multiplicity > 0:
            leaf_sequence = node_data.get("sequence_alignment", "")

            if germline_alignment and leaf_sequence:
                num_mutations, seq_length, germline_aligned, leaf_aligned = (
                    align_and_calculate_mutations(
                        germline_alignment, leaf_sequence, alignment_method
                    )
                )
                mut_freq = num_mutations / seq_length if seq_length > 0 else 0.0
                total_mut_freq += mut_freq * multiplicity
                total_multiplicity += multiplicity

                mutation_positions = []
                for pos, (g, leaf_base) in enumerate(
                    zip(germline_aligned, leaf_aligned)
                ):
                    if (
                        g != leaf_base
                        and g not in NON_COMPARABLE
                        and leaf_base not in NON_COMPARABLE
                    ):
                        mutation_positions.append(
                            {"pos": pos, "germline": g, "leaf": leaf_base}
                        )

                debug_info.append(
                    {
                        "node": node_id,
                        "type": node_type,
                        "distance": node_data.get("distance", 0.0),
                        "num_mutations": num_mutations,
                        "seq_length": seq_length,
                        "original_leaf_len": len(leaf_sequence),
                        "original_germline_len": len(germline_alignment),
                        "mut_freq": mut_freq,
                        "multiplicity": multiplicity,
                        "weighted_contribution": mut_freq * multiplicity,
                        "germline_seq": germline_aligned,
                        "leaf_seq": leaf_aligned,
                        "mutations": mutation_positions,
                        "was_aligned": len(germline_alignment) != len(leaf_sequence),
                        "alignment_method": alignment_method,
                    }
                )
            else:
                # Track why we skipped this sequence - be specific
                if not leaf_sequence:
                    reason = "missing sequence (empty)"
                elif not germline_alignment:
                    reason = "missing germline sequence"
                elif len(leaf_sequence) != len(germline_alignment):
                    reason = f"sequence length mismatch (leaf={len(leaf_sequence)}, germline={len(germline_alignment)})"
                else:
                    reason = "unknown"

                skipped_nodes.append(
                    {
                        "node": node_id,
                        "type": node_type,
                        "multiplicity": multiplicity,
                        "reason": reason,
                        "has_sequence": bool(leaf_sequence),
                        "seq_len": len(leaf_sequence) if leaf_sequence else 0,
                        "germline_len": len(germline_alignment)
                        if germline_alignment
                        else 0,
                    }
                )
        else:
            # Track all non-leaf nodes or leaves with multiplicity 0
            reason = ""
            if node_type != "leaf":
                reason = f"not a leaf (type={node_type})"
            elif multiplicity == 0:
                reason = "leaf with multiplicity=0 (no observed sequences)"
            else:
                reason = "unknown"

            skipped_nodes.append(
                {
                    "node": node_id,
                    "type": node_type,
                    "multiplicity": multiplicity,
                    "reason": reason,
                }
            )

    mean_mut_freq = (
        total_mut_freq / total_multiplicity if total_multiplicity > 0 else 0.0
    )
    return mean_mut_freq, debug_info, skipped_nodes
