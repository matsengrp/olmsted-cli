#!/usr/bin/env python

"""
Process Parent-Child Pair (PCP) CSV files and Newick trees for Olmsted visualization.

This script handles PCP CSV format with columns:
- sample_id: Sample identifier
- parent_name: Parent node name
- child_name: Child node name
- edge_length: Branch length
- sample_count: Number of sequences

And a CSV file containing Newick trees:
- family_name: Family identifier
- newick_tree: Newick format tree string
"""

import argparse
import csv
import gzip
import hashlib
import html
import os
import sys
import traceback
import uuid
from collections import defaultdict
from urllib.parse import parse_qs, parse_qsl

# Python 3.13+ compatibility: make cgi module available before ete3 import
try:
    import cgi  # noqa: F401
except ImportError:
    # Create a mock cgi module using our compatibility layer

    class CGIModule:
        """Mock cgi module for Python 3.13+ compatibility."""

        escape = html.escape

        # Add other cgi functions that might be needed by ete3
        def parse_qs(self, *args, **kwargs):
            return parse_qs(*args, **kwargs)

        def parse_qsl(self, *args, **kwargs):
            return parse_qsl(*args, **kwargs)

    # Make cgi available as a module
    sys.modules["cgi"] = CGIModule()

import ete3
from tqdm import tqdm

# Import shared utilities from process_data_utils
from .process_utils import (
    SCHEMA_VERSION,
    VerbosePrinter,
    create_consolidated_data,
    translate_dna_to_aa,
    validate_output_data,
    write_out,
)


def parse_pcp_csv(csv_path):
    """
    Parse PCP CSV file and return a dict of families with rich immunological data.

    Expected CSV format (required columns):
    sample_id,parent_name,child_name

    Optional columns for rich immunological annotations:
    - family: Family identifier (defaults to sample_id if not present)
    - parent_heavy, child_heavy: DNA sequences
    - branch_length or edge_length: Branch length
    - sample_count: Number of sequences
    - v_gene_heavy, d_gene_heavy, j_gene_heavy: Gene calls
    - v_gene_start_heavy, v_gene_end_heavy: V gene alignment positions
    - d_gene_start_heavy, d_gene_end_heavy: D gene alignment positions
    - j_gene_start_heavy, j_gene_end_heavy: J gene alignment positions
    - cdr1_codon_start_heavy, cdr1_codon_end_heavy: CDR1 positions
    - cdr2_codon_start_heavy, cdr2_codon_end_heavy: CDR2 positions
    - cdr3_codon_start_heavy, cdr3_codon_end_heavy: CDR3 positions
    - parent_is_naive, child_is_leaf: Boolean flags
    - distance: Distance from root

    Returns:
        dict: {family_id: {
            nodes: {node_id: node_data},
            edges: [(parent, child, length)],
            family_data: {v_gene, d_gene, j_gene, gene_positions, cdr_positions, etc.}
        }}
    """
    families = defaultdict(lambda: {"nodes": {}, "edges": [], "family_data": {}})

    # Determine if file is gzipped
    if csv_path.endswith(".gz"):
        file_handle = gzip.open(csv_path, "rt")
    else:
        file_handle = open(csv_path, "r")

    with file_handle:
        reader = csv.DictReader(file_handle)

        # Validate required columns (flexible format support)
        required_cols = {"sample_id", "parent_name", "child_name"}
        if not required_cols.issubset(reader.fieldnames):
            missing = required_cols - set(reader.fieldnames)
            raise ValueError(f"Missing required columns: {missing}")

        for row in reader:
            sample_id = row["sample_id"]
            family_id = row.get(
                "family", sample_id
            )  # Use family if available, fallback to sample_id
            parent = row["parent_name"]
            child = row["child_name"]

            # Handle different edge length column names
            edge_length = 0.0
            if "branch_length" in row:
                edge_length = float(row["branch_length"])
            elif "edge_length" in row:
                edge_length = float(row["edge_length"])

            # Handle sample count (default to 1 if not present)
            sample_count = 1
            if "sample_count" in row:
                sample_count = int(row["sample_count"])

            # Extract rich immunological fields
            parent_sequence = row.get("parent_heavy", "")
            child_sequence = row.get("child_heavy", "")
            v_gene = row.get("v_gene_heavy", "")
            d_gene = row.get("d_gene_heavy", "")
            j_gene = row.get("j_gene_heavy", "")
            parent_is_naive = row.get("parent_is_naive", "").lower() == "true"
            child_is_leaf = row.get("child_is_leaf", "").lower() == "true"

            # Extract distance/mutation data
            distance = float(row.get("distance", 0)) if row.get("distance") else 0.0
            branch_length = (
                float(row.get("branch_length", 0)) if row.get("branch_length") else 0.0
            )

            # Extract CDR position data
            cdr1_start = (
                int(row.get("cdr1_codon_start_heavy", 0))
                if row.get("cdr1_codon_start_heavy")
                else 0
            )
            cdr1_end = (
                int(row.get("cdr1_codon_end_heavy", 0))
                if row.get("cdr1_codon_end_heavy")
                else 0
            )
            cdr2_start = (
                int(row.get("cdr2_codon_start_heavy", 0))
                if row.get("cdr2_codon_start_heavy")
                else 0
            )
            cdr2_end = (
                int(row.get("cdr2_codon_end_heavy", 0))
                if row.get("cdr2_codon_end_heavy")
                else 0
            )
            cdr3_start = (
                int(row.get("cdr3_codon_start_heavy", 0))
                if row.get("cdr3_codon_start_heavy")
                else 0
            )
            cdr3_end = (
                int(row.get("cdr3_codon_end_heavy", 0))
                if row.get("cdr3_codon_end_heavy")
                else 0
            )

            # Extract gene position data (V, D, J gene start/end positions)
            v_gene_start = (
                int(row.get("v_gene_start_heavy", 0))
                if row.get("v_gene_start_heavy")
                else None
            )
            v_gene_end = (
                int(row.get("v_gene_end_heavy", 0))
                if row.get("v_gene_end_heavy")
                else None
            )
            d_gene_start = (
                int(row.get("d_gene_start_heavy", 0))
                if row.get("d_gene_start_heavy")
                else None
            )
            d_gene_end = (
                int(row.get("d_gene_end_heavy", 0))
                if row.get("d_gene_end_heavy")
                else None
            )
            j_gene_start = (
                int(row.get("j_gene_start_heavy", 0))
                if row.get("j_gene_start_heavy")
                else None
            )
            j_gene_end = (
                int(row.get("j_gene_end_heavy", 0))
                if row.get("j_gene_end_heavy")
                else None
            )

            # Store family-level data and sample_id for each family (will be same for all rows of same family)
            families[family_id]["family_data"] = {
                "sample_id": sample_id,  # Store original sample_id for reference
                "v_gene": v_gene,
                "d_gene": d_gene,
                "j_gene": j_gene,
                "v_gene_start": v_gene_start,
                "v_gene_end": v_gene_end,
                "d_gene_start": d_gene_start,
                "d_gene_end": d_gene_end,
                "j_gene_start": j_gene_start,
                "j_gene_end": j_gene_end,
                "cdr1_start": cdr1_start,
                "cdr1_end": cdr1_end,
                "cdr2_start": cdr2_start,
                "cdr2_end": cdr2_end,
                "cdr3_start": cdr3_start,
                "cdr3_end": cdr3_end,
            }

            # Add parent node if not already present
            if parent not in families[family_id]["nodes"]:
                families[family_id]["nodes"][parent] = {
                    "sequence_id": parent,
                    "multiplicity": 0,
                    "timepoint_multiplicities": [],
                    "sequence_alignment": parent_sequence,
                    "is_naive": parent_is_naive,
                    "is_leaf": False,
                    "distances": [],  # Will be updated if this node appears as child
                    "distance": 0.0
                    if parent_is_naive
                    else None,  # Will be set when node appears as child
                    "length": 0.0
                    if parent_is_naive
                    else None,  # Will be set when node appears as child
                }

            # Add child node if not already present
            if child not in families[family_id]["nodes"]:
                families[family_id]["nodes"][child] = {
                    "sequence_id": child,
                    "multiplicity": sample_count,
                    "timepoint_multiplicities": [],
                    "sequence_alignment": child_sequence,
                    "is_naive": False,
                    "is_leaf": child_is_leaf,
                    "distances": [distance] if distance > 0 else [],
                    "distance": distance,  # Distance from root
                    "length": branch_length,  # Branch length to this node
                }
            else:
                # Update multiplicity if node appears multiple times
                families[family_id]["nodes"][child]["multiplicity"] += sample_count

            # Update parent node distance/length if this parent appears as a child in another row
            if parent in families[family_id]["nodes"]:
                parent_node = families[family_id]["nodes"][parent]
                # Only update if not already set (and not naive)
                if parent_node["distance"] is None and not parent_node["is_naive"]:
                    # We need to find the row where this parent is a child to get its distance
                    # For now, we'll handle this in a post-processing step
                    pass
                # Add distance data
                if distance > 0:
                    families[family_id]["nodes"][child]["distances"].append(distance)
                # Update distance and length if this edge provides better data
                if distance > 0:
                    families[family_id]["nodes"][child]["distance"] = distance
                if branch_length > 0:
                    families[family_id]["nodes"][child]["length"] = branch_length

            # Add edge
            families[family_id]["edges"].append((parent, child, edge_length))

    # Post-process to ensure all nodes have correct distance/length values
    for family_id, family_data in families.items():
        _fix_node_distances_and_lengths(family_data)

    return dict(families)


def _fix_node_distances_and_lengths(family_data):
    """
    Post-process family data to ensure all nodes have correct distance and length values.

    Args:
        family_data: Dictionary containing nodes and edges for one family
    """
    nodes = family_data["nodes"]
    edges = family_data["edges"]

    # Create lookup for child -> (parent, edge_length, child_distance)
    child_info = {}
    for parent, child, edge_length in edges:
        if child in nodes:
            child_distance = nodes[child]["distance"]
            child_info[child] = (parent, edge_length, child_distance)

    # Update parent nodes that don't have distance/length set
    for node_id, node_data in nodes.items():
        if node_data["distance"] is None or node_data["length"] is None:
            # This node is a parent but we haven't seen it as a child yet
            # Check if it appears as a child in child_info
            if node_id in child_info:
                parent, edge_length, child_distance = child_info[node_id]
                if node_data["distance"] is None:
                    node_data["distance"] = child_distance
                if node_data["length"] is None:
                    node_data["length"] = edge_length
                # Add to distances list for mutation frequency calculation
                if child_distance > 0 and child_distance not in node_data["distances"]:
                    node_data["distances"].append(child_distance)


def parse_newick_tree(newick_string):
    """
    Parse a Newick string to extract complete tree topology using ETE3.

    This parser uses ETE3 for robust tree handling and ensures that if a
    'naive' node exists, it becomes the root (as it represents the unmutated
    ancestor in B cell lineage trees).

    Args:
        newick_string: Newick format tree string

    Returns:
        tuple: (nodes_dict, edges_list, root_node_id)
            - nodes_dict: {node_id: {"label": node_id, "branch_length": float}}
            - edges_list: [(parent, child, branch_length), ...]
            - root_node_id: ID of the root node
    """
    nodes = {}
    edges = []

    # Parse tree with ETE3
    try:
        tree = ete3.Tree(newick_string, format=1)  # format=1 means flexible with internal node names
    except Exception as e:
        print(f"Warning: Failed to parse Newick tree with ETE3: {e}")
        # Return empty structure on parse failure
        return {}, [], None

    # Don't try to reroot with ETE3 - we'll handle it after extracting edges

    # Build a mapping of ETE3 nodes to names for consistent naming
    # Handle duplicate names by adding suffixes
    node_name_map = {}
    name_counts = {}
    internal_node_counter = 1

    # First pass - count names to detect duplicates
    for node in tree.traverse():
        if node.name:
            name = node.name
            name_counts[name] = name_counts.get(name, 0) + 1

    # Second pass - assign unique names
    name_usage = {}
    for node in tree.traverse():
        if node.name:
            base_name = node.name
            # If name appears multiple times, add suffix
            if name_counts[base_name] > 1:
                usage_count = name_usage.get(base_name, 0)
                if usage_count > 0:
                    node_name_map[node] = f"{base_name}_dup{usage_count}"
                else:
                    node_name_map[node] = base_name
                name_usage[base_name] = usage_count + 1
            else:
                node_name_map[node] = base_name
        else:
            # Generate name for unnamed internal nodes
            while f"Node{internal_node_counter}" in name_counts:
                internal_node_counter += 1
            node_name_map[node] = f"Node{internal_node_counter}"
            internal_node_counter += 1

    # Extract nodes and edges from ETE3 tree
    for node in tree.traverse():
        node_name = node_name_map[node]

        # Add node to nodes dict
        nodes[node_name] = {
            "label": node_name,
            "branch_length": node.dist
        }

        # Add edge from parent to this node (if not root)
        if not node.is_root():
            parent = node.up
            parent_name = node_name_map[parent]
            edges.append((parent_name, node_name, node.dist))

    # Get root node
    root = tree.get_tree_root()
    root_node_id = node_name_map[root]

    # Special handling: if "naive" exists but is not the root, fix the edges
    if "naive" in nodes and root_node_id != "naive":
        # Find path from naive to current root and reverse those edges
        fixed_edges = []
        parent_map = {child: (parent, length) for parent, child, length in edges}

        # Find path from naive to root
        path_to_reverse = []
        current = "naive"
        while current in parent_map:
            parent, length = parent_map[current]
            path_to_reverse.append((parent, current, length))
            current = parent
            if current == root_node_id:
                break

        # Create set of edges to reverse
        edges_to_reverse = set((p, c) for p, c, _ in path_to_reverse)

        # Rebuild edge list with reversed edges
        for parent, child, length in edges:
            if (parent, child) in edges_to_reverse:
                # Reverse this edge
                fixed_edges.append((child, parent, length))
            else:
                fixed_edges.append((parent, child, length))

        edges = fixed_edges
        root_node_id = "naive"

    return nodes, edges, root_node_id


def merge_tree_topology_with_pcp(pcp_family_data, newick_string, warn_disagreements=False, family_id=None):
    """
    Merge complete tree topology from Newick with PCP data.

    This function ensures that the family topology includes ALL nodes from the
    Newick tree, even if some nodes are missing from the PCP data. This solves
    the issue where PCP files have incomplete topology.

    Args:
        pcp_family_data: Family data from parse_pcp_csv
        newick_string: Newick tree string for this family
        warn_disagreements: If True, print warnings when tree and PCP data disagree
        family_id: Family identifier for warning messages

    Returns:
        dict: Updated family data with complete topology
    """
    # Parse the Newick tree to get complete topology
    tree_nodes, tree_edges, tree_root = parse_newick_tree(newick_string)

    # Start with existing PCP data
    merged_nodes = pcp_family_data["nodes"].copy()
    merged_edges = pcp_family_data["edges"].copy()

    # Add any missing nodes from the tree
    for node_id, tree_node in tree_nodes.items():
        if node_id not in merged_nodes:
            # Create minimal node data for nodes not in PCP
            merged_nodes[node_id] = {
                "multiplicity": 0,  # No observed sequences
                "timepoint_multiplicities": [],
                "sequence_alignment": "",  # No sequence data
                "is_naive": node_id == "naive",
                "is_leaf": node_id not in {parent for parent, _, _ in tree_edges},
                "distances": [],
                "distance": None,  # Will be calculated from tree
                "length": tree_node["branch_length"],
            }

    # Build a set of nodes that are in the tree
    tree_node_ids = set(tree_nodes.keys())

    # Build dictionaries for easy edge lookup
    tree_edge_dict = {(parent, child): branch_length for parent, child, branch_length in tree_edges}
    pcp_edge_dict = {(parent, child): edge_length for parent, child, edge_length in pcp_family_data["edges"]}

    # Track disagreements if requested
    disagreements = []

    # Check for disagreements between tree and PCP edges
    if warn_disagreements:
        for (parent, child), tree_length in tree_edge_dict.items():
            if (parent, child) in pcp_edge_dict:
                pcp_length = pcp_edge_dict[(parent, child)]
                # Check if branch lengths disagree (with tolerance for floating point)
                if abs(tree_length - pcp_length) > 1e-10:
                    disagreements.append({
                        "type": "branch_length",
                        "edge": (parent, child),
                        "tree_value": tree_length,
                        "pcp_value": pcp_length
                    })

        # Check for edges that exist in PCP but not in tree (for shared nodes)
        for (parent, child), pcp_length in pcp_edge_dict.items():
            if parent in tree_node_ids and child in tree_node_ids:
                # Both nodes are in tree, but edge might not be
                if (parent, child) not in tree_edge_dict:
                    disagreements.append({
                        "type": "missing_edge",
                        "edge": (parent, child),
                        "source": "PCP has edge not in tree"
                    })

        # Check for edges that exist in tree but not in PCP (for shared nodes)
        for (parent, child), tree_length in tree_edge_dict.items():
            if parent in pcp_family_data["nodes"] and child in pcp_family_data["nodes"]:
                # Both nodes are in PCP, but edge might not be
                if (parent, child) not in pcp_edge_dict:
                    disagreements.append({
                        "type": "missing_edge",
                        "edge": (parent, child),
                        "source": "Tree has edge not in PCP"
                    })

        # Check for nodes that exist in tree but not in PCP
        pcp_node_ids = set(pcp_family_data["nodes"].keys())
        for node_id in tree_node_ids:
            if node_id not in pcp_node_ids:
                disagreements.append({
                    "type": "missing_node",
                    "node": node_id,
                    "source": "Tree has node not in PCP"
                })

        # Check for nodes that exist in PCP but not in tree
        for node_id in pcp_node_ids:
            if node_id not in tree_node_ids:
                disagreements.append({
                    "type": "missing_node",
                    "node": node_id,
                    "source": "PCP has node not in tree"
                })

    # Rebuild edges: use tree edges for tree nodes, keep PCP edges for others
    merged_edges = []

    # First, add all tree edges
    for parent, child, branch_length in tree_edges:
        merged_edges.append((parent, child, branch_length))

    # Then, add PCP edges for nodes not in the tree
    # (these are additional sequences that aren't in the phylogenetic reconstruction)
    pcp_edges_to_keep = []
    for parent, child, edge_length in pcp_family_data["edges"]:
        # Keep edge if either parent or child is not in the tree
        # (these represent additional observed sequences)
        if child not in tree_node_ids:
            # This is a PCP-only node, keep its edge
            pcp_edges_to_keep.append((parent, child, edge_length))

    # Add the PCP-only edges
    merged_edges.extend(pcp_edges_to_keep)

    # Print warnings if requested and disagreements were found
    if warn_disagreements and disagreements:
        display_id = family_id if family_id else pcp_family_data.get("family_data", {}).get("sample_id", "unknown")
        print(f"\nWarning: Found {len(disagreements)} disagreement(s) in family {display_id}:")
        for i, disagree in enumerate(disagreements, 1):
            if disagree["type"] == "branch_length":
                print(f"  {i}. Branch length mismatch for edge {disagree['edge'][0]} -> {disagree['edge'][1]}:")
                print(f"     Tree: {disagree['tree_value']:.10f}, PCP: {disagree['pcp_value']:.10f}")
            elif disagree["type"] == "missing_edge":
                if disagree["source"] == "PCP has edge not in tree":
                    print(f"  {i}. Edge {disagree['edge'][0]} -> {disagree['edge'][1]} exists in PCP but not in tree")
                elif disagree["source"] == "Tree has edge not in PCP":
                    print(f"  {i}. Edge {disagree['edge'][0]} -> {disagree['edge'][1]} exists in tree but not in PCP")
            elif disagree["type"] == "missing_node":
                if disagree["source"] == "PCP has node not in tree":
                    print(f"  {i}. Node '{disagree['node']}' exists in PCP but not in tree")
                elif disagree["source"] == "Tree has node not in PCP":
                    print(f"  {i}. Node '{disagree['node']}' exists in tree but not in PCP")

    # Update distances based on tree topology
    # Calculate distance from root for each node
    def calculate_distances_from_root(edges, root):
        """Calculate cumulative distance from root for each node."""
        distances = {root: 0.0}

        # Build adjacency list
        children_dict = {}
        edge_lengths = {}
        for parent, child, length in edges:
            if parent not in children_dict:
                children_dict[parent] = []
            children_dict[parent].append(child)
            edge_lengths[(parent, child)] = length

        # DFS to calculate distances
        def dfs(node, current_distance):
            distances[node] = current_distance
            if node in children_dict:
                for child in children_dict[node]:
                    edge_length = edge_lengths.get((node, child), 0.0)
                    dfs(child, current_distance + edge_length)

        dfs(root, 0.0)
        return distances

    # Calculate distances from root using all merged edges
    node_distances = calculate_distances_from_root(merged_edges, tree_root)

    # Update node distances and preserve PCP metadata
    for node_id in merged_nodes:
        merged_nodes[node_id]["distance"] = node_distances.get(node_id, 0.0)
        # Keep existing distance data if available
        if merged_nodes[node_id]["distances"]:
            # Use existing distances from PCP data
            pass
        else:
            # Use calculated distance
            if node_distances.get(node_id, 0.0) > 0:
                merged_nodes[node_id]["distances"] = [node_distances[node_id]]

    # Return updated family data
    return {
        "family_data": pcp_family_data["family_data"],
        "nodes": merged_nodes,
        "edges": merged_edges,
    }


def parse_newick_csv(csv_path):
    """
    Parse CSV file containing Newick trees.

    Expected CSV format:
    family_name,sample_id,newick_tree (or family_name,newick_tree for backwards compatibility)

    Returns:
        dict: {(family_name, sample_id): newick_string} if sample_id present,
              {family_name: newick_string} otherwise
    """
    newick_trees = {}

    # Determine if file is gzipped
    if csv_path.endswith(".gz"):
        file_handle = gzip.open(csv_path, "rt")
    else:
        file_handle = open(csv_path, "r")

    with file_handle:
        reader = csv.DictReader(file_handle)

        # Validate required columns
        required_cols = {"family_name", "newick_tree"}
        if not required_cols.issubset(reader.fieldnames):
            missing = required_cols - set(reader.fieldnames)
            raise ValueError(f"Missing required columns: {missing}")

        # Check if sample_id column exists
        has_sample_id = "sample_id" in reader.fieldnames

        for row in reader:
            family_name = row["family_name"]
            newick_tree = row["newick_tree"]

            if has_sample_id:
                sample_id = row["sample_id"]
                # Use composite key (family_name, sample_id) to handle multiple samples with same family ID
                newick_trees[(family_name, sample_id)] = newick_tree
            else:
                # Backwards compatibility: use just family_name
                newick_trees[family_name] = newick_tree

    return newick_trees


def compute_lbi_for_tree(nodes_dict, edges, root_id, tau=0.0125):
    """
    Compute Local Branching Index (LBI) for all nodes in a tree.

    LBI measures the local branching structure around each node, capturing
    the rate of diversification in the recent evolutionary history.

    Reference: Neher & Bedford (2015) "nextflu: real-time tracking of seasonal influenza
    virus evolution in humans" Bioinformatics 31(21):3546-3548

    Algorithm:
    1. For each node, calculate exponentially weighted sum of branches below it
    2. Weight each branch by exp(-distance/tau) where tau is time scale parameter
    3. LBI = sum of exponentially weighted branches in the subtree below

    Args:
        nodes_dict: Dictionary of {node_id: node_data} with 'distance' and 'length' fields
        edges: List of (parent, child, length) tuples
        root_id: ID of the root node
        tau: Time scale parameter (default: 0.0125)

    Returns:
        dict: {node_id: lbi_value}
    """
    import math
    from collections import defaultdict

    # Build data structures for tree traversal
    children_map = defaultdict(list)  # {parent: [children]}
    parent_map = {}  # {child: parent}
    edge_length_map = {}  # {(parent, child): length}

    for parent, child, length in edges:
        children_map[parent].append(child)
        parent_map[child] = parent
        edge_length_map[(parent, child)] = length

    # Initialize up and down polarizers
    up_polarizer = {node_id: 0.0 for node_id in nodes_dict.keys()}
    down_polarizer = {node_id: 0.0 for node_id in nodes_dict.keys()}

    # Postorder traversal: Calculate up_polarizer (tree length below each node)
    def postorder(node):
        """
        Calculate up_polarizer for a node and all its descendants.
        up_polarizer[node] = sum over children of exponentially weighted (branch_length + subtree_length)
        """
        if node not in children_map or len(children_map[node]) == 0:
            # Leaf node
            up_polarizer[node] = 0.0
            return

        # Process children first (postorder)
        for child in children_map[node]:
            postorder(child)

        # Calculate up_polarizer for this node
        total = 0.0
        for child in children_map[node]:
            branch_length = edge_length_map.get((node, child), 0.0)
            # Add exponentially weighted contribution from child subtree
            # Weight decreases with distance: exp(-branch_length/tau)
            weight = math.exp(-branch_length / tau)
            total += (branch_length + up_polarizer[child]) * weight

        up_polarizer[node] = total

    # Preorder traversal: Calculate down_polarizer (tree length above each node)
    def preorder(node):
        """
        Calculate down_polarizer for a node and all its descendants.
        down_polarizer[child] = exponentially weighted (down_polarizer[parent] + sibling contributions)
        """
        if node not in children_map:
            return

        # Process each child
        for child in children_map[node]:
            branch_length = edge_length_map.get((node, child), 0.0)
            weight = math.exp(-branch_length / tau)

            # Calculate contribution from parent's down_polarizer and siblings
            parent_contribution = down_polarizer[node]

            # Add contributions from sibling subtrees
            sibling_contribution = 0.0
            for sibling in children_map[node]:
                if sibling != child:
                    sibling_branch_length = edge_length_map.get((node, sibling), 0.0)
                    sibling_contribution += sibling_branch_length + up_polarizer[sibling]

            # Set down_polarizer for child
            down_polarizer[child] = (branch_length + parent_contribution + sibling_contribution) * weight

            # Recursively process child's descendants
            preorder(child)

    # Execute traversals
    postorder(root_id)
    down_polarizer[root_id] = 0.0  # Root has no ancestors
    preorder(root_id)

    # Calculate LBI = up_polarizer + down_polarizer
    lbi = {}
    for node_id in nodes_dict.keys():
        lbi[node_id] = up_polarizer[node_id] + down_polarizer[node_id]

    return lbi


def compute_lbr_for_tree(nodes_dict, edges, root_id):
    """
    Compute Local Branching Ratio (LBR) for all nodes in a tree.

    LBR measures the relative branching rate comparing upstream vs downstream lineages.

    Algorithm:
    1. For each node, count number of descendant branches (downstream)
    2. Count number of ancestral branches (upstream)
    3. LBR = log(downstream_branches / upstream_branches)

    Args:
        nodes_dict: Dictionary of {node_id: node_data}
        edges: List of (parent, child, length) tuples
        root_id: ID of the root node

    Returns:
        dict: {node_id: lbr_value}
    """
    import math
    from collections import defaultdict

    # Build data structures for tree traversal
    children_map = defaultdict(list)  # {parent: [children]}
    parent_map = {}  # {child: parent}

    for parent, child, length in edges:
        children_map[parent].append(child)
        parent_map[child] = parent

    # Count downstream branches for each node (postorder traversal)
    downstream_count = {}

    def count_downstream(node):
        """
        Count the number of branches (edges) in the subtree rooted at node.
        """
        if node not in children_map or len(children_map[node]) == 0:
            # Leaf node has 0 downstream branches
            downstream_count[node] = 0
            return 0

        # Count branches to children plus all branches in their subtrees
        total = 0
        for child in children_map[node]:
            total += 1  # Edge to this child
            total += count_downstream(child)  # All branches in child's subtree

        downstream_count[node] = total
        return total

    # Count upstream branches for each node (distance to root in edges)
    upstream_count = {}

    def count_upstream(node):
        """
        Count the number of branches (edges) from node to the root.
        """
        if node == root_id:
            upstream_count[node] = 0
            return 0

        # Count edges along path to root
        count = 0
        current = node
        while current != root_id and current in parent_map:
            count += 1
            current = parent_map[current]

        upstream_count[node] = count
        return count

    # Calculate counts for all nodes
    count_downstream(root_id)
    for node_id in nodes_dict.keys():
        count_upstream(node_id)

    # Calculate LBR = log(downstream / upstream)
    lbr = {}
    for node_id in nodes_dict.keys():
        down = downstream_count.get(node_id, 0)
        up = upstream_count.get(node_id, 0)

        # Handle edge cases - default to 0.0 for undefined values
        if up == 0 and down == 0:
            # Root node with no children (degenerate tree)
            lbr[node_id] = 0.0
        elif up == 0:
            # Root node or nodes with no ancestors
            # LBR is undefined (log(down/0))
            # Default to 0.0 for root nodes
            lbr[node_id] = 0.0
        elif down == 0:
            # Leaf node with no descendants
            # LBR = log(0/up) is undefined
            # Default to 0.0 for leaf nodes
            lbr[node_id] = 0.0
        else:
            # Standard case: LBR = log(downstream / upstream)
            lbr[node_id] = math.log(down / up)

    return lbr


def compute_cluster_multiplicity_for_tree(nodes_dict, edges, root_id):
    """
    Compute cluster multiplicity for all nodes in a tree.

    Cluster multiplicity is the sum of multiplicities of all descendant leaves.
    For leaf nodes, cluster_multiplicity equals their own multiplicity.
    For internal nodes, it's the sum of all leaf multiplicities in the subtree.

    Args:
        nodes_dict: Dictionary of {node_id: node_data} with 'multiplicity' field
        edges: List of (parent, child, length) tuples
        root_id: ID of the root node

    Returns:
        dict: {node_id: cluster_multiplicity}
    """
    from collections import defaultdict

    # Build children map
    children_map = defaultdict(list)
    for parent, child, length in edges:
        children_map[parent].append(child)

    # Postorder traversal to compute cluster multiplicities
    cluster_mult = {}

    def postorder(node):
        """
        Compute cluster_multiplicity for node and all descendants.
        """
        if node not in children_map or len(children_map[node]) == 0:
            # Leaf node: cluster_multiplicity = own multiplicity
            cluster_mult[node] = nodes_dict[node].get("multiplicity", 0)
            return cluster_mult[node]

        # Internal node: sum of all descendant leaf multiplicities
        total = 0
        for child in children_map[node]:
            total += postorder(child)

        cluster_mult[node] = total
        return total

    # Compute for entire tree
    postorder(root_id)

    return cluster_mult


def compute_scaled_affinity(affinity_values):
    """
    Compute scaled affinity using min-max normalization.

    scaled_affinity = (affinity - min) / (max - min)

    Nodes with None affinity will have None scaled_affinity.

    Args:
        affinity_values: dict of {node_id: affinity}

    Returns:
        dict: {node_id: scaled_affinity}
    """
    # Get non-None affinity values
    valid_affinities = {k: v for k, v in affinity_values.items() if v is not None}

    if not valid_affinities:
        # No valid affinity values, return all None
        return {k: None for k in affinity_values.keys()}

    # Find min and max
    min_affinity = min(valid_affinities.values())
    max_affinity = max(valid_affinities.values())

    # Handle case where all affinities are the same
    if max_affinity == min_affinity:
        # All nodes have same affinity, scale to 0.5 or 1.0
        return {k: (0.5 if v is not None else None) for k, v in affinity_values.items()}

    # Min-max normalization
    scaled = {}
    for node_id, affinity in affinity_values.items():
        if affinity is None:
            scaled[node_id] = None
        else:
            scaled[node_id] = (affinity - min_affinity) / (max_affinity - min_affinity)

    return scaled


def build_newick_from_edges(nodes, edges):
    """
    Build a Newick string from parent-child edges.

    Args:
        nodes: dict of {node_id: node_data}
        edges: list of (parent, child, edge_length) tuples

    Returns:
        str: Newick format tree string
    """
    # Build adjacency list
    children = defaultdict(list)
    edge_lengths = {}

    for parent, child, length in edges:
        children[parent].append(child)
        edge_lengths[(parent, child)] = length

    # Find root (node with no parent)
    all_children = {child for _, child, _ in edges}
    all_parents = {parent for parent, _, _ in edges}
    roots = all_parents - all_children

    if len(roots) != 1:
        raise ValueError(f"Expected exactly one root, found {len(roots)}: {roots}")

    root = roots.pop()

    def build_subtree(node, parent_node=None):
        """Recursively build Newick subtree."""
        if node not in children:
            # Leaf node
            edge_key = (parent_node, node) if parent_node else (node, node)
            edge_len = edge_lengths.get(edge_key, 0.0)
            return f"{node}:{edge_len}"

        # Internal node
        subtrees = []
        for child in children[node]:
            subtrees.append(build_subtree(child, node))

        edge_key = (parent_node, node) if parent_node else (node, node)
        edge_len = edge_lengths.get(edge_key, 0.0)
        return f"({','.join(subtrees)}){node}:{edge_len}"

    # Build the tree starting from root
    # Root doesn't have a parent edge, so handle specially
    if root not in children:
        return f"{root}:0.0;"

    subtrees = []
    for child in children[root]:
        subtrees.append(build_subtree(child, root))

    return f"({','.join(subtrees)}){root}:0.0;"


def process_pcp_to_olmsted(pcp_families, newick_trees=None, uuid_generator=None, warn_disagreements=False, compute_metrics=False, lbi_tau=0.0125, standardize_names=False, alignment_method="truncate", name=None, verbosity=1):
    """
    Convert PCP format data to Olmsted format.

    Args:
        pcp_families: dict from parse_pcp_csv
        newick_trees: dict from parse_newick_csv (optional)
        uuid_generator: Function to generate UUIDs (defaults to random)
        warn_disagreements: If True, print warnings when tree and PCP data disagree
        compute_metrics: If True, compute all phylogenetic metrics (LBI, LBR, affinity, scaled_affinity, mean_mut_freq)
        lbi_tau: Time scale parameter for LBI calculation (default: 0.0125)
        standardize_names: If True, rename nodes to naive/Node1/Node2.../Leaf1/Leaf2...
        alignment_method: Method for sequence alignment ("truncate" or "pad", default: "truncate")
        name: Optional name for the dataset (default: None)
        verbosity: Verbosity level (0=quiet, 1=normal, 2=verbose, 3=debug)

    Returns:
        tuple: (datasets, clones_dict, trees)
    """
    # Create verbosity printer
    vprint = VerbosePrinter(verbosity)

    if uuid_generator is None:
        uuid_generator = lambda prefix="": f"{prefix}{str(uuid.uuid4())}" if prefix else str(uuid.uuid4())

    dataset_id = f"pcp-{uuid_generator()}"
    dataset_ident = uuid_generator("dataset-")

    datasets = []
    clones_dict = {dataset_id: []}
    trees = []

    # Create dataset
    dataset = {
        "ident": dataset_ident,
        "dataset_id": dataset_id,
        "schema_version": SCHEMA_VERSION,
        "type": "pcp.dataset",
        "build": {"commit": "pcp-import", "time": ""},
        "subjects": [{"ident": uuid_generator("subject-"), "subject_id": "pcp-subject"}],
        "samples": [],
        "seeds": [],
        "clone_count": len(pcp_families),
        "subjects_count": 1,
        "timepoints_count": 1,
    }

    # Add name if provided
    if name:
        dataset["name"] = name

    # Process each family with progress bar
    family_items = list(pcp_families.items())
    with tqdm(family_items, desc="Processing families", unit="family", disable=len(family_items) == 1) as pbar:
        for family_idx, (family_id, family_data) in enumerate(pbar):
            clone_ident = uuid_generator("clone-")
            tree_ident = uuid_generator()

            # Get sample_id from family data
            family_meta = family_data.get("family_data", {})
            original_sample_id = family_meta.get("sample_id", family_id)

            # Create sample if not already present
            sample_exists = any(
                s["sample_id"] == original_sample_id for s in dataset["samples"]
            )
            if not sample_exists:
                dataset["samples"].append(
                    {
                        "ident": uuid_generator("sample-"),
                        "sample_id": original_sample_id,
                        "locus": "igh",  # Default locus
                        "timepoint_id": "merged",
                    }
                )

            # Merge tree topology with PCP data if Newick tree is available
            # Try composite key first (family_id, sample_id), then fall back to just family_id
            newick = None
            if newick_trees:
                # Try composite key (family_id, sample_id)
                composite_key = (family_id, original_sample_id)
                if composite_key in newick_trees:
                    newick = newick_trees[composite_key]
                # Fall back to just family_id for backwards compatibility
                elif family_id in newick_trees:
                    newick = newick_trees[family_id]

            if newick:
                # Use complete tree topology from Newick
                family_data = merge_tree_topology_with_pcp(family_data, newick, warn_disagreements, family_id)
            else:
                # Fallback to building tree from PCP edges only
                newick = build_newick_from_edges(family_data["nodes"], family_data["edges"])

            # Process nodes - add required fields with rich PCP data
            processed_nodes = {}
            for node_id, node_data in family_data["nodes"].items():
                # Get sequence alignment from PCP data
                sequence_alignment = node_data.get("sequence_alignment", "")
                sequence_alignment_aa = translate_dna_to_aa(sequence_alignment)

                # Determine node type based on PCP metadata
                if node_data.get("is_naive", False):
                    node_type = "root"
                elif node_data.get("is_leaf", False):
                    node_type = "leaf"
                else:
                    # This is an internal/ancestral node (Node1, Node2, etc.)
                    node_type = "internal"

                processed_node = {
                    "sequence_id": node_id,
                    "sequence_alignment": sequence_alignment,
                    "sequence_alignment_aa": sequence_alignment_aa,
                    "multiplicity": node_data.get("multiplicity", 0),
                    "cluster_multiplicity": None,  # Will be computed below
                    "timepoint_multiplicities": node_data.get(
                        "timepoint_multiplicities", []
                    ),
                    "type": node_type,
                    "parent": None,  # Will be set from edges below
                    "distance": node_data.get("distance", 0.0),  # Distance from root
                    "length": node_data.get("length", 0.0),  # Branch length
                    "lbi": None,
                    "lbr": None,
                    "affinity": None,
                    "scaled_affinity": None,
                }
                processed_nodes[node_id] = processed_node

            # Set parent field based on edges
            # First, find the true root (node that doesn't appear as a child in any edge)
            all_children = {child for _, child, _ in family_data["edges"]}
            all_parents = {parent for parent, _, _ in family_data["edges"]}
            potential_roots = all_parents - all_children

            # Determine the root node
            if potential_roots:
                tree_root = potential_roots.pop()
            else:
                # If no clear root, use "naive" if present, otherwise use first node
                tree_root = "naive" if "naive" in processed_nodes else list(processed_nodes.keys())[0]

            # Set parent relationships, but ensure root has no parent
            for parent_id, child_id, edge_length in family_data["edges"]:
                if child_id in processed_nodes and child_id != tree_root:
                    processed_nodes[child_id]["parent"] = parent_id

            # Ensure the root node has no parent
            if tree_root in processed_nodes:
                processed_nodes[tree_root]["parent"] = None

            # Calculate cluster multiplicity (always computed)
            vprint.verbose(f"  Computing cluster multiplicity for family {family_id}")
            cluster_mult_values = compute_cluster_multiplicity_for_tree(processed_nodes, family_data["edges"], tree_root)
            for node_id in processed_nodes:
                processed_nodes[node_id]["cluster_multiplicity"] = cluster_mult_values.get(node_id, 0)

            # Calculate phylogenetic metrics if requested
            if compute_metrics:
                # Compute LBI
                vprint.verbose(f"  Computing LBI for family {family_id} with tau={lbi_tau}")
                lbi_values = compute_lbi_for_tree(processed_nodes, family_data["edges"], tree_root, tau=lbi_tau)
                for node_id in processed_nodes:
                    lbi = lbi_values.get(node_id)
                    processed_nodes[node_id]["lbi"] = lbi
                    # Set affinity = LBI
                    processed_nodes[node_id]["affinity"] = lbi

                # Compute scaled_affinity (min-max normalization)
                vprint.verbose(f"  Computing scaled affinity for family {family_id}")
                affinity_values = {node_id: processed_nodes[node_id]["affinity"] for node_id in processed_nodes}
                scaled_affinity_values = compute_scaled_affinity(affinity_values)
                for node_id in processed_nodes:
                    processed_nodes[node_id]["scaled_affinity"] = scaled_affinity_values.get(node_id)

                # Compute LBR
                vprint.verbose(f"  Computing LBR for family {family_id}")
                lbr_values = compute_lbr_for_tree(processed_nodes, family_data["edges"], tree_root)
                for node_id in processed_nodes:
                    processed_nodes[node_id]["lbr"] = lbr_values.get(node_id)

            # Standardize node names if requested
            if standardize_names:
                # Create name mapping: old_name -> new_name
                name_mapping = {}
                internal_counter = 1
                leaf_counter = 1

                # First pass: create mapping
                for node_id, node_data in processed_nodes.items():
                    node_type = node_data.get("type")
                    if node_type == "root":
                        name_mapping[node_id] = "naive"
                    elif node_type == "leaf":
                        name_mapping[node_id] = f"Leaf{leaf_counter}"
                        leaf_counter += 1
                    else:  # internal
                        name_mapping[node_id] = f"Node{internal_counter}"
                        internal_counter += 1

                # Second pass: rename nodes and update parent references
                renamed_nodes = {}
                for old_name, node_data in processed_nodes.items():
                    new_name = name_mapping[old_name]
                    # Update parent reference to use new name
                    if node_data["parent"] and node_data["parent"] in name_mapping:
                        node_data["parent"] = name_mapping[node_data["parent"]]
                    # Update sequence_id to new name
                    node_data["sequence_id"] = new_name
                    # Store under new name
                    renamed_nodes[new_name] = node_data

                processed_nodes = renamed_nodes
                tree_root = name_mapping.get(tree_root, tree_root)

            # Extract family-level immunological data (already extracted above)
            v_call = family_meta.get("v_gene", "")
            d_call = family_meta.get("d_gene", "")
            j_call = family_meta.get("j_gene", "")

            # Get CDR positions
            cdr1_start = family_meta.get("cdr1_start", 0)
            cdr1_end = family_meta.get("cdr1_end", 0)
            cdr2_start = family_meta.get("cdr2_start", 0)
            cdr2_end = family_meta.get("cdr2_end", 0)
            cdr3_start = family_meta.get("cdr3_start", 0)
            cdr3_end = family_meta.get("cdr3_end", 0)

            # Get gene positions from CSV if available - no fallbacks, don't guess
            # V gene alignment positions
            v_alignment_start = family_meta.get("v_gene_start")
            if v_alignment_start is None:
                v_alignment_start = 0

            v_alignment_end = family_meta.get("v_gene_end")
            if v_alignment_end is None:
                v_alignment_end = 0

            # D gene alignment positions
            d_alignment_start = family_meta.get("d_gene_start")
            if d_alignment_start is None:
                d_alignment_start = 0

            d_alignment_end = family_meta.get("d_gene_end")
            if d_alignment_end is None:
                d_alignment_end = 0

            # J gene alignment positions
            j_alignment_start = family_meta.get("j_gene_start")
            if j_alignment_start is None:
                j_alignment_start = 0

            j_alignment_end = family_meta.get("j_gene_end")
            if j_alignment_end is None:
                j_alignment_end = 0

            junction_start = cdr3_start
            junction_length = (cdr3_end - cdr3_start) if (cdr3_end > cdr3_start) else 0

            # Get germline sequence from naive node first (needed for mean_mut_freq calculation)
            germline_alignment = ""
            for node_id, node_data in processed_nodes.items():
                if node_data.get("type") == "root":
                    germline_alignment = node_data.get("sequence_alignment", "")
                    break

            # Calculate mean mutation frequency from observed leaf sequences only
            # mean_mut_freq = average(mutations_per_site) across all leaf nodes, weighted by multiplicity
            # Count actual mutations by comparing leaf sequence to germline sequence
            total_mut_freq = 0.0
            total_sequences = 0
            germline_length = len(germline_alignment) if germline_alignment else 0

            # DEBUG: Print calculation details for all sequences
            debug_info = []
            skipped_nodes = []

            for node_id, node_data in processed_nodes.items():
                node_type = node_data.get("type")
                multiplicity = node_data.get("multiplicity", 0)
                # Only count LEAF nodes with observed sequences (type="leaf" and multiplicity > 0)
                # Skip internal nodes (type="internal") and root node (type="root")
                if node_type == "leaf" and multiplicity > 0:
                    leaf_sequence = node_data.get("sequence_alignment", "")

                    # Count mutations by comparing to germline
                    if germline_alignment and leaf_sequence:
                        if alignment_method == "truncate":
                            # Truncate to shorter sequence length (default behavior)
                            min_len = min(len(germline_alignment), len(leaf_sequence))
                            germline_compared = germline_alignment[:min_len]
                            leaf_compared = leaf_sequence[:min_len]

                            # Count mutations in overlapping region
                            num_mutations = sum(1 for g, l in zip(germline_compared, leaf_compared)
                                              if g != l and g != '' and l != '')
                            # Calculate frequency based on truncated length
                            mut_freq = num_mutations / min_len if min_len > 0 else 0.0

                            total_mut_freq += mut_freq * multiplicity
                            total_sequences += multiplicity

                            # Collect mutation positions for display
                            mutation_positions = []
                            for pos, (g, l) in enumerate(zip(germline_compared, leaf_compared)):
                                if g != l and g != '' and l != '':
                                    mutation_positions.append({
                                        'pos': pos,
                                        'germline': g,
                                        'leaf': l
                                    })

                            # Collect debug info
                            debug_info.append({
                                'node': node_id,
                                'type': node_type,
                                'distance': node_data.get("distance", 0.0),
                                'num_mutations': num_mutations,
                                'seq_length': min_len,  # Truncated length
                                'original_leaf_len': len(leaf_sequence),
                                'original_germline_len': len(germline_alignment),
                                'mut_freq': mut_freq,
                                'multiplicity': multiplicity,
                                'weighted_contribution': mut_freq * multiplicity,
                                'germline_seq': germline_compared,
                                'leaf_seq': leaf_compared,
                                'mutations': mutation_positions,
                                'was_truncated': len(germline_alignment) != len(leaf_sequence)
                            })

                        elif alignment_method == "pad":
                            # Pad the shorter sequence with "." (optional behavior)
                            max_len = max(len(germline_alignment), len(leaf_sequence))
                            germline_padded = germline_alignment.ljust(max_len, ".")
                            leaf_padded = leaf_sequence.ljust(max_len, ".")

                            # Use the padded sequences for comparison
                            # Count mutations, treating "." as a gap (not a mutation)
                            num_mutations = sum(1 for g, l in zip(germline_padded, leaf_padded)
                                              if g != l and g != '' and l != '' and g != '.' and l != '.')
                            # Calculate frequency based on the padded alignment length
                            mut_freq = num_mutations / max_len if max_len > 0 else 0.0

                            total_mut_freq += mut_freq * multiplicity
                            total_sequences += multiplicity

                            # Collect mutation positions for display (using padded sequences)
                            mutation_positions = []
                            for pos, (g, l) in enumerate(zip(germline_padded, leaf_padded)):
                                if g != l and g != '' and l != '' and g != '.' and l != '.':
                                    mutation_positions.append({
                                        'pos': pos,
                                        'germline': g,
                                        'leaf': l
                                    })

                            # Collect debug info for sequences we can calculate
                            debug_info.append({
                                'node': node_id,
                                'type': node_type,
                                'distance': node_data.get("distance", 0.0),
                                'num_mutations': num_mutations,
                                'seq_length': max_len,  # Padded/aligned length
                                'original_leaf_len': len(leaf_sequence),  # Original leaf sequence length
                                'original_germline_len': len(germline_alignment),  # Original germline length
                                'mut_freq': mut_freq,
                                'multiplicity': multiplicity,
                                'weighted_contribution': mut_freq * multiplicity,
                                'germline_seq': germline_padded,  # Use padded sequence
                                'leaf_seq': leaf_padded,  # Use padded sequence
                                'mutations': mutation_positions,
                                'was_padded': len(germline_alignment) != len(leaf_sequence)
                            })
                    else:
                        # Track why we skipped this sequence - be specific
                        if not leaf_sequence:
                            reason = 'missing sequence (empty)'
                        elif not germline_alignment:
                            reason = 'missing germline sequence'
                        elif len(leaf_sequence) != len(germline_alignment):
                            reason = f'sequence length mismatch (leaf={len(leaf_sequence)}, germline={len(germline_alignment)})'
                        else:
                            reason = 'unknown'

                        skipped_nodes.append({
                            'node': node_id,
                            'type': node_type,
                            'multiplicity': multiplicity,
                            'reason': reason,
                            'has_sequence': bool(leaf_sequence),
                            'seq_len': len(leaf_sequence) if leaf_sequence else 0,
                            'germline_len': len(germline_alignment) if germline_alignment else 0
                        })
                else:
                    # Track all non-leaf nodes or leaves with multiplicity 0
                    reason = ''
                    if node_type != "leaf":
                        reason = f'not a leaf (type={node_type})'
                    elif multiplicity == 0:
                        reason = 'leaf with multiplicity=0 (no observed sequences)'
                    else:
                        reason = 'unknown'

                    skipped_nodes.append({
                        'node': node_id,
                        'type': node_type,
                        'multiplicity': multiplicity,
                        'reason': reason
                    })

            mean_mut_freq = total_mut_freq / total_sequences if total_sequences > 0 else 0.0

            # DEBUG: Print for all families (only at debug verbosity level)
            vprint.debug(f"\n=== DEBUG: mean_mut_freq calculation for family {family_id} ===")
            vprint.debug(f"Germline length: {germline_length} nt")

            # Translate germline to amino acids
            germline_aa = translate_dna_to_aa(germline_alignment)
            vprint.debug(f"Germline AA length: {len(germline_aa)} aa")

            vprint.debug(f"\nLEAF sequences only ({len(debug_info)} total):")
            for info in debug_info:
                padded_marker = " [PADDED]" if info.get('was_padded', False) else ""
                # Show original lengths to verify they match the actual sequence data
                orig_leaf = info.get('original_leaf_len', 0)
                orig_germ = info.get('original_germline_len', 0)
                aligned_len = info['seq_length']

                length_info = f"leaf_len={orig_leaf}, germ_len={orig_germ}, aligned_len={aligned_len}"

                vprint.debug(f"\n  Node {info['node']} (type={info['type']}){padded_marker}: "
                          f"distance={info['distance']:.6f}, "
                          f"mutations={info['num_mutations']:.1f} nt, "
                          f"{length_info}, "
                          f"mut_freq={info['mut_freq']:.6f}, "
                          f"multiplicity={info['multiplicity']}, "
                          f"weighted={info['weighted_contribution']:.6f}")

                # Translate leaf sequence to amino acids (before padding - translate original sequences)
                # Remove padding before translation
                leaf_seq_original = info['leaf_seq'].rstrip('.')
                germline_seq_original = info['germline_seq'].rstrip('.')

                leaf_aa = translate_dna_to_aa(leaf_seq_original)
                germline_aa_for_this_leaf = translate_dna_to_aa(germline_seq_original)

                # Now pad the amino acid sequences with "."
                max_aa_length = max(len(germline_aa_for_this_leaf), len(leaf_aa))
                germline_aa_padded = germline_aa_for_this_leaf.ljust(max_aa_length, ".")
                leaf_aa_padded = leaf_aa.ljust(max_aa_length, ".")

                # Create alignment display (show full sequence)
                naive_line = ""
                node_line = ""

                for i in range(max_aa_length):
                    naive_aa = germline_aa_padded[i]
                    leaf_aa_char = leaf_aa_padded[i]

                    naive_line += naive_aa
                    # Show "-" where sequences match, show actual AA where they differ
                    if naive_aa == leaf_aa_char:
                        node_line += "-"
                    else:
                        node_line += leaf_aa_char

                vprint.debug(f"    Naive: {naive_line}")
                vprint.debug(f"    Node:  {node_line}")

                # Count AA mutations
                aa_mutations = sum(1 for i in range(min(len(germline_aa), len(leaf_aa)))
                                  if germline_aa[i] != leaf_aa[i])
                vprint.debug(f"    AA mutations: {aa_mutations}")

            vprint.debug(f"\nTotal mutation frequency (weighted): {total_mut_freq:.6f}")
            vprint.debug(f"Total leaf sequences: {total_sequences}")
            vprint.debug(f"Mean mutation frequency: {mean_mut_freq:.6f}")
            vprint.debug(f"  (This means {mean_mut_freq*100:.2f}% of positions have mutations on average)")

            # Show skipped nodes summary
            if skipped_nodes:
                vprint.debug(f"\nSkipped {len(skipped_nodes)} nodes:")
                # Group by reason
                by_reason = {}
                for node in skipped_nodes:
                    reason = node['reason']
                    if reason not in by_reason:
                        by_reason[reason] = []
                    by_reason[reason].append(node)

                for reason, nodes in by_reason.items():
                    vprint.debug(f"  {reason}: {len(nodes)} nodes")
                    # Show first few examples
                    for node in nodes[:3]:
                        seq_info = f"seq_len={node['seq_len']}" if node.get('has_sequence') else "no_seq"
                        vprint.debug(f"    - {node['node']} (type={node['type']}, mult={node['multiplicity']}, {seq_info})")
                    if len(nodes) > 3:
                        vprint.debug(f"    ... and {len(nodes) - 3} more")

            vprint.debug(f"===================================================\n")

            # Create clone with rich PCP data
            clone = {
                "clone_id": family_id,  # Use actual family name from PCP data
                "ident": clone_ident,
                "dataset_id": dataset_id,
                "sample_id": original_sample_id,
                "subject_id": "pcp-subject",
                "unique_seqs_count": len(processed_nodes),
                "total_read_count": sum(
                    n.get("multiplicity", 0) for n in processed_nodes.values()
                ),
                "mean_mut_freq": mean_mut_freq,
                "v_alignment_start": v_alignment_start,
                "v_alignment_end": v_alignment_end,
                "j_alignment_start": j_alignment_start,
                "j_alignment_end": j_alignment_end,
                "cdr1_alignment_start": cdr1_start,
                "cdr1_alignment_end": cdr1_end,
                "cdr2_alignment_start": cdr2_start,
                "cdr2_alignment_end": cdr2_end,
                "junction_start": junction_start,
                "junction_length": junction_length,
                "v_call": v_call,
                "d_call": d_call,
                "j_call": j_call,
                "d_alignment_start": d_alignment_start,
                "d_alignment_end": d_alignment_end,
                "germline_alignment": germline_alignment,
                "has_seed": False,
                "trees": [
                    {
                        "ident": tree_ident,
                        "clone_id": family_id,  # Use actual family name
                        "tree_id": f"pcp-tree-{family_id}",  # Use family name in tree ID
                        "newick": newick,
                        "type": "pcp.reconstruction",  # PCP-specific type
                    }
                ],
                # Add nested sample and dataset objects for webapp compatibility
                "sample": {
                    "ident": clone_ident,
                    "locus": "igh",
                    "sample_id": original_sample_id,
                    "timepoint_id": "merged",
                },
                "dataset": {"ident": dataset_ident, "dataset_id": dataset_id},
            }
            clones_dict[dataset_id].append(clone)

            # Convert nodes to array format (required by webapp)
            nodes_array = []
            for node_id, node_data in processed_nodes.items():
                nodes_array.append(node_data)

            # Create tree with nodes as array
            tree = {
                "ident": tree_ident,
                "tree_id": f"pcp-tree-{family_id}",  # Use family name in tree ID
                "clone_id": family_id,  # Use actual family name
                "newick": newick,
                "nodes": nodes_array,
                "type": "pcp.reconstruction",  # PCP-specific reconstruction type
            }
            trees.append(tree)

    datasets.append(dataset)
    return datasets, clones_dict, trees


def deterministic_uuid(seed_base, counter=None):
    """Generate a deterministic UUID based on a seed and optional counter."""
    if counter is not None:
        seed_str = f"{seed_base}_{counter}"
    else:
        seed_str = str(seed_base)

    # Create a hash of the seed string
    hash_obj = hashlib.md5(seed_str.encode())
    hash_hex = hash_obj.hexdigest()

    # Convert to UUID format: xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx
    uuid_str = f"{hash_hex[:8]}-{hash_hex[8:12]}-{hash_hex[12:16]}-{hash_hex[16:20]}-{hash_hex[20:32]}"
    return uuid_str


def get_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Process PCP CSV and Newick files for Olmsted visualization"
    )
    parser.add_argument(
        "-i", "--input-pcp", required=True, help="Input PCP CSV file (can be gzipped)"
    )
    parser.add_argument(
        "-t",
        "--input-trees",
        help="Input CSV file containing Newick trees (optional, can be gzipped)",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output file path for consolidated JSON (default behavior)",
    )
    parser.add_argument(
        "--split-files",
        metavar="DIR",
        dest="output_dir",
        help="Output to multiple files in specified directory (datasets.json, clones.*.json, tree.*.json) instead of single consolidated file",
    )
    parser.add_argument(
        "-n",
        "--name",
        help="Optional name for the dataset (stored in metadata)",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        type=int,
        choices=[0, 1, 2, 3],
        default=1,
        help="Set verbosity level: 0=quiet (errors only), 1=normal (default), 2=verbose, 3=debug",
    )
    parser.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="Quiet mode - only show errors (equivalent to -v 0)",
    )
    parser.add_argument(
        "--validate",
        action="store_true",
        help="Validate output data against JSON schemas before writing",
    )
    parser.add_argument(
        "--strict-validation",
        action="store_true",
        help="Exit with error if validation fails (requires --validate)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        help="Random seed for deterministic UUID generation (useful for testing)",
    )
    parser.add_argument(
        "-w",
        "--warnings",
        dest="warn_disagreements",
        action="store_true",
        help="Show warnings when tree and PCP data disagree on edges or branch lengths",
    )
    parser.add_argument(
        "--compute-metrics",
        action="store_true",
        help="Compute phylogenetic metrics (LBI, LBR, affinity, scaled_affinity, mean_mut_freq) for all nodes",
    )
    parser.add_argument(
        "--lbi-tau",
        type=float,
        default=0.0125,
        help="Time scale parameter for LBI calculation (default: 0.0125)",
    )
    parser.add_argument(
        "--standardize-names",
        action="store_true",
        help="Rename nodes to standardized names: naive (root), Node1, Node2, ... (internal), Leaf1, Leaf2, ... (leaves)",
    )
    parser.add_argument(
        "--alignment-method",
        choices=["truncate", "pad"],
        default="truncate",
        help="Method for aligning sequences of different lengths for mutation frequency calculation (default: truncate - compare only overlapping region, pad - pad shorter sequence with gap characters)",
    )
    # Removed --output-format option - now only outputs AIRR format
    return parser.parse_args()


def main():
    """Main entry point."""
    args = get_args()

    # Handle quiet mode
    if args.quiet:
        args.verbose = 0

    # Create verbosity printer
    vprint = VerbosePrinter(args.verbose)

    # Print command arguments at verbosity level 2
    vprint.verbose("=== Command Arguments ===")
    vprint.verbose(f"  Input PCP file: {args.input_pcp}")
    if args.input_trees:
        vprint.verbose(f"  Input trees file: {args.input_trees}")
    if args.output:
        vprint.verbose(f"  Output file: {args.output}")
    if args.output_dir:
        vprint.verbose(f"  Output directory: {args.output_dir}")
    if args.name:
        vprint.verbose(f"  Dataset name: {args.name}")
    vprint.verbose(f"  Verbosity level: {args.verbose}")
    vprint.verbose(f"  Validation: {args.validate}")
    if args.validate:
        vprint.verbose(f"  Strict validation: {args.strict_validation}")
    if args.seed is not None:
        vprint.verbose(f"  Random seed: {args.seed}")
    vprint.verbose(f"  Show disagreement warnings: {args.warn_disagreements}")
    vprint.verbose(f"  Compute metrics: {args.compute_metrics}")
    if args.compute_metrics:
        vprint.verbose(f"    LBI tau: {args.lbi_tau}")
    vprint.verbose(f"  Standardize names: {args.standardize_names}")
    vprint.verbose(f"  Alignment method: {args.alignment_method}")
    vprint.verbose("=" * 25)
    vprint.verbose("")

    # Set up deterministic UUID generation if seed is provided
    uuid_counter = 0

    def get_uuid(prefix=""):
        nonlocal uuid_counter
        if args.seed is not None:
            uuid_counter += 1
            uuid_str = deterministic_uuid(args.seed, uuid_counter)
        else:
            uuid_str = str(uuid.uuid4())
        return f"{prefix}{uuid_str}" if prefix else uuid_str

    try:
        # Parse PCP CSV
        vprint.status(f"Processing PCP CSV: {args.input_pcp}")
        if args.seed is not None:
            vprint.status(f"Using deterministic UUIDs with seed: {args.seed}")
        pcp_families = parse_pcp_csv(args.input_pcp)
        vprint.status(f"Found {len(pcp_families)} families")

        # Parse Newick trees if provided
        newick_trees = None
        if args.input_trees:
            vprint.status(f"Processing Newick trees: {args.input_trees}")
            newick_trees = parse_newick_csv(args.input_trees)
            vprint.status(f"Found {len(newick_trees)} trees")

        # Convert to Olmsted format
        vprint.status("Converting to Olmsted format...")
        datasets, clones_dict, trees = process_pcp_to_olmsted(
            pcp_families, newick_trees, get_uuid, args.warn_disagreements,
            compute_metrics=args.compute_metrics, lbi_tau=args.lbi_tau,
            standardize_names=args.standardize_names, alignment_method=args.alignment_method,
            name=args.name, verbosity=args.verbose
        )

        # Validate output data if requested
        if args.validate:
            if not validate_output_data(datasets, clones_dict, trees, args):
                if args.strict_validation:
                    vprint.error(
                        "\nExiting due to validation errors (--strict-validation enabled)"
                    )
                    sys.exit(1)

        # Write output
        if args.output_dir:
            # Multi-file output to specified directory
            os.makedirs(args.output_dir, exist_ok=True)
            vprint.status(f"Writing multiple files to {args.output_dir}")
            write_out(datasets, args.output_dir, "datasets.json", args)
            for dataset_id, clones in clones_dict.items():
                write_out(clones, args.output_dir, f"clones.{dataset_id}.json", args)
            for tree in trees:
                write_out(tree, args.output_dir, f"tree.{tree['ident']}.json", args)
        else:
            # Single consolidated file output (default)
            # Build input files list for metadata
            input_files = [args.input_pcp]
            if args.input_trees:
                input_files.append(args.input_trees)

            consolidated_data = create_consolidated_data(
                datasets, clones_dict, trees, input_files, "pcp", args
            )
            # Ensure output directory exists
            output_dir = os.path.dirname(args.output) or "."
            output_file = os.path.basename(args.output)
            os.makedirs(output_dir, exist_ok=True)
            vprint.status(f"Writing consolidated output to {args.output}")
            write_out(consolidated_data, output_dir, output_file, args)

        vprint.status("Processing complete!")

    except Exception as e:
        vprint.error(f"Error: {e}")
        if args.verbose:
            traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
