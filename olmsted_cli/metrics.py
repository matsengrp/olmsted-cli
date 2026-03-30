"""
Shared phylogenetic metric computations.

These functions compute tree-based metrics (LBI, LBR, scaled affinity)
from tree topology and branch lengths. They are format-agnostic and can
be used with data from any input format (PCP, AIRR, Olmsted JSON).
"""

from __future__ import annotations

import math
from collections import defaultdict
from typing import Any, Dict, List, Optional, Tuple


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
                    sibling_branch_length = edge_length_map.get(
                        (node, sibling), 0.0
                    )
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
        return {
            k: (0.5 if v is not None else None)
            for k, v in affinity_values.items()
        }

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
