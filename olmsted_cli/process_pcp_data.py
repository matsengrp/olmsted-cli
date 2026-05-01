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

from __future__ import annotations

import argparse
import copy
import csv
import json
import html
import os
import sys
import traceback
from collections import defaultdict
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, Dict, List, Literal, Optional, Tuple
from urllib.parse import parse_qs, parse_qsl

from .identifier import IdentMinter, deterministic_uuid  # noqa: F401 (re-exported for back-compat)

if TYPE_CHECKING:
    from .types import OlmstedClone, OlmstedDataset, OlmstedNode, OlmstedTree

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
    coerce_csv_value,
    create_consolidated_data,
    tag_field_metadata,
    get_optional_int,
    translate_dna_to_aa,
    validate_output_data,
    write_out,
)
from .data_io import open_file
from .utils import set_verbosity, vprint


from .constants import CHAIN_COLUMN_ALIASES, FORMAT_PCP, KNOWN_PCP_COLUMNS, KNOWN_TREE_COLUMNS
from .metrics import compute_tree_metrics


def _normalize_column_names(fieldnames):
    """
    Normalize PCP CSV column names by mapping common aliases to canonical names.

    Recognizes chain-agnostic names (e.g., 'v_gene') and maps them to
    the chain-specific canonical form (e.g., 'v_gene_heavy') when no
    chain-specific version is already present.

    Returns:
        Tuple of (column_map, notifications):
        - column_map: dict mapping original column name -> canonical name
        - notifications: list of notification strings about remapped columns
    """
    column_map = {}
    notifications = []
    canonical_set = set(fieldnames)

    for orig_name in fieldnames:
        if not orig_name:
            continue
        lower = orig_name.lower().strip()

        # Check alias map
        if lower in CHAIN_COLUMN_ALIASES:
            canonical = CHAIN_COLUMN_ALIASES[lower]
            # Only remap if the canonical name isn't already present
            if canonical not in canonical_set:
                column_map[orig_name] = canonical
                notifications.append(
                    f"Column '{orig_name}' mapped to '{canonical}'"
                )
            else:
                # Both alias and canonical present — treat alias as extra
                column_map[orig_name] = orig_name
        else:
            column_map[orig_name] = orig_name

    return column_map, notifications


def _partition_chain_fields(fields):
    """
    Partition a dict of extra fields into shared, heavy-only, and light-only.

    Fields ending with '_heavy' go to heavy-only (with suffix stripped).
    Fields ending with '_light' go to light-only (with suffix stripped).
    Fields without a chain suffix are shared between both chains.

    Returns:
        Tuple of (shared, heavy_only, light_only) dicts.
    """
    shared = {}
    heavy_only = {}
    light_only = {}
    for key, val in fields.items():
        if key.endswith("_heavy"):
            heavy_only[key[:-6]] = val  # strip _heavy suffix
        elif key.endswith("_light"):
            light_only[key[:-6]] = val  # strip _light suffix
        else:
            shared[key] = val
    return shared, heavy_only, light_only


def infer_locus_from_v_gene(v_gene: str) -> Optional[str]:
    """
    Infer locus from V gene call.

    Args:
        v_gene: V gene call (e.g., "IGHV1-2*01", "IGKV3-20*01", "IGLV2-14*04")

    Returns:
        Locus string ("igh", "igk", or "igl") or None if v_gene is missing
        or its prefix is not recognized.
    """
    if not v_gene:
        return None
    v_gene_upper = v_gene.upper()
    if v_gene_upper.startswith("IGKV"):
        return "igk"
    elif v_gene_upper.startswith("IGLV"):
        return "igl"
    elif v_gene_upper.startswith("IGHV"):
        return "igh"
    return None


def parse_pcp_csv(csv_path: str) -> Dict[str, Any]:
    """
    Parse PCP CSV file and return a dict of families with rich immunological data.

    Expected CSV format (required columns):
    sample_id,parent_name,child_name

    Optional columns for rich immunological annotations:
    - family: Family identifier (defaults to sample_id if not present)
    - parent_heavy, child_heavy: Heavy chain DNA sequences
    - branch_length or edge_length: Branch length
    - sample_count: Number of sequences
    - v_gene_heavy, d_gene_heavy, j_gene_heavy: Heavy chain gene calls
    - v_gene_start_heavy, v_gene_end_heavy: V gene alignment positions
    - d_gene_start_heavy, d_gene_end_heavy: D gene alignment positions
    - j_gene_start_heavy, j_gene_end_heavy: J gene alignment positions
    - cdr1_codon_start_heavy, cdr1_codon_end_heavy: CDR1 positions (heavy)
    - cdr2_codon_start_heavy, cdr2_codon_end_heavy: CDR2 positions (heavy)
    - cdr3_codon_start_heavy, cdr3_codon_end_heavy: CDR3 positions (heavy)
    - parent_is_naive, child_is_leaf: Boolean flags
    - distance: Distance from root

    Paired heavy/light chain columns (optional):
    - parent_light, child_light: Light chain DNA sequences
    - v_gene_light, j_gene_light: Light chain gene calls
    - cdr1_codon_start_light, cdr1_codon_end_light: CDR1 positions (light)
    - cdr2_codon_start_light, cdr2_codon_end_light: CDR2 positions (light)
    - cdr3_codon_start_light, cdr3_codon_end_light: CDR3 positions (light)
    - light_chain_type: "kappa" or "lambda"

    Returns:
        dict: {family_id: {
            nodes: {node_id: node_data},
            edges: [(parent, child, length)],
            family_data: {v_gene, d_gene, j_gene, gene_positions, cdr_positions, etc.},
            is_paired: bool  # True if light chain data is present
        }}
    """
    families = defaultdict(lambda: {"nodes": {}, "edges": [], "family_data": {}})

    handle, _ = open_file(csv_path, expected_formats=(FORMAT_PCP,))
    with handle as file_handle:
        reader = csv.DictReader(file_handle)

        # Normalize column names (map aliases like v_gene -> v_gene_heavy)
        column_map, column_notifications = _normalize_column_names(reader.fieldnames)
        normalized_fieldnames = {column_map.get(c, c) for c in reader.fieldnames if c}

        # Report column remapping
        for note in column_notifications:
            vprint.verbose(f"  Note: {note}")

        # Validate required columns (flexible format support)
        required_cols = {"sample_id", "parent_name", "child_name"}
        if not required_cols.issubset(normalized_fieldnames):
            missing = required_cols - normalized_fieldnames
            raise ValueError(f"Missing required columns: {missing}")

        # Identify extra columns not handled by the standard parser
        # Filter out empty/None column names (e.g., unnamed index columns)
        extra_columns = {
            c for c in normalized_fieldnames if c not in KNOWN_PCP_COLUMNS
        }

        # Detect format type based on normalized column presence
        has_heavy = "parent_heavy" in normalized_fieldnames or "child_heavy" in normalized_fieldnames
        has_light = "parent_light" in normalized_fieldnames or "child_light" in normalized_fieldnames

        # Determine format:
        # - Paired: has both heavy and light columns
        # - Light-only: has only light columns
        # - Heavy-only: has only heavy columns (or neither)
        is_paired = has_heavy and has_light
        is_light_only = has_light and not has_heavy

        for pcp_index, raw_row in enumerate(reader):
            # Apply column name normalization
            row = {column_map.get(k, k): v for k, v in raw_row.items() if k}
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
            # Check per-row which sequences are present (not just column existence)
            # This handles mixed datasets with heavy-only, light-only, and paired families
            parent_heavy_seq = row.get("parent_heavy", "")
            parent_light_seq = row.get("parent_light", "")
            child_heavy_seq = row.get("child_heavy", "")
            child_light_seq = row.get("child_light", "")

            # Determine which chain to use as primary sequence for this row
            # Prefer heavy if available, otherwise use light
            if parent_heavy_seq or child_heavy_seq:
                # This row has heavy chain data
                parent_sequence = parent_heavy_seq
                child_sequence = child_heavy_seq
                v_gene = row.get("v_gene_heavy", "")
                d_gene = row.get("d_gene_heavy", "")
                j_gene = row.get("j_gene_heavy", "")
            else:
                # This row only has light chain data
                parent_sequence = parent_light_seq
                child_sequence = child_light_seq
                v_gene = row.get("v_gene_light", "")
                d_gene = row.get("d_gene_light", "")
                j_gene = row.get("j_gene_light", "")
            parent_is_naive = row.get("parent_is_naive", "").lower() == "true"
            child_is_leaf = row.get("child_is_leaf", "").lower() == "true"

            # Extract distance/mutation data
            distance = float(row.get("distance", 0)) if row.get("distance") else 0.0
            branch_length = (
                float(row.get("branch_length", 0)) if row.get("branch_length") else 0.0
            )

            # Extract CDR position data
            # Use same chain as primary sequence (determined above)
            if parent_heavy_seq or child_heavy_seq:
                # Heavy chain CDRs
                cdr1_start = get_optional_int(row, "cdr1_codon_start_heavy")
                cdr1_end = get_optional_int(row, "cdr1_codon_end_heavy")
                cdr2_start = get_optional_int(row, "cdr2_codon_start_heavy")
                cdr2_end = get_optional_int(row, "cdr2_codon_end_heavy")
                cdr3_start = get_optional_int(row, "cdr3_codon_start_heavy")
                cdr3_end = get_optional_int(row, "cdr3_codon_end_heavy")
            else:
                # Light chain CDRs
                cdr1_start = get_optional_int(row, "cdr1_codon_start_light")
                cdr1_end = get_optional_int(row, "cdr1_codon_end_light")
                cdr2_start = get_optional_int(row, "cdr2_codon_start_light")
                cdr2_end = get_optional_int(row, "cdr2_codon_end_light")
                cdr3_start = get_optional_int(row, "cdr3_codon_start_light")
                cdr3_end = get_optional_int(row, "cdr3_codon_end_light")

            # Extract gene position data (V, D, J gene start/end positions)
            # Use same chain as primary sequence (determined above)
            if parent_heavy_seq or child_heavy_seq:
                # Heavy chain gene positions
                v_gene_start = get_optional_int(row, "v_gene_start_heavy", default=None)
                v_gene_end = get_optional_int(row, "v_gene_end_heavy", default=None)
                d_gene_start = get_optional_int(row, "d_gene_start_heavy", default=None)
                d_gene_end = get_optional_int(row, "d_gene_end_heavy", default=None)
                j_gene_start = get_optional_int(row, "j_gene_start_heavy", default=None)
                j_gene_end = get_optional_int(row, "j_gene_end_heavy", default=None)
            else:
                # Light chain gene positions
                v_gene_start = get_optional_int(row, "v_gene_start_light", default=None)
                v_gene_end = get_optional_int(row, "v_gene_end_light", default=None)
                d_gene_start = get_optional_int(row, "d_gene_start_light", default=None)
                d_gene_end = get_optional_int(row, "d_gene_end_light", default=None)
                j_gene_start = get_optional_int(row, "j_gene_start_light", default=None)
                j_gene_end = get_optional_int(row, "j_gene_end_light", default=None)

            # Extract light chain data (for paired format)
            # Only extract if this row has BOTH heavy AND light sequences
            row_is_paired = (parent_heavy_seq or child_heavy_seq) and (parent_light_seq or child_light_seq)
            parent_sequence_light = parent_light_seq if row_is_paired else ""
            child_sequence_light = child_light_seq if row_is_paired else ""
            v_gene_light = row.get("v_gene_light", "") if row_is_paired else ""
            j_gene_light = row.get("j_gene_light", "") if row_is_paired else ""
            light_chain_type = row.get("light_chain_type", "") if row_is_paired else ""

            # Extract light chain CDR positions (for paired format)
            cdr1_start_light = (
                int(row.get("cdr1_codon_start_light", 0))
                if row.get("cdr1_codon_start_light")
                else 0
            ) if row_is_paired else 0
            cdr1_end_light = (
                int(row.get("cdr1_codon_end_light", 0))
                if row.get("cdr1_codon_end_light")
                else 0
            ) if row_is_paired else 0
            cdr2_start_light = (
                int(row.get("cdr2_codon_start_light", 0))
                if row.get("cdr2_codon_start_light")
                else 0
            ) if row_is_paired else 0
            cdr2_end_light = (
                int(row.get("cdr2_codon_end_light", 0))
                if row.get("cdr2_codon_end_light")
                else 0
            ) if row_is_paired else 0
            cdr3_start_light = (
                int(row.get("cdr3_codon_start_light", 0))
                if row.get("cdr3_codon_start_light")
                else 0
            ) if row_is_paired else 0
            cdr3_end_light = (
                int(row.get("cdr3_codon_end_light", 0))
                if row.get("cdr3_codon_end_light")
                else 0
            ) if row_is_paired else 0

            # Store family-level data from the first row seen for each family.
            # Guard: only populate once to avoid silent overwrites from later
            # rows that may have incomplete or inconsistent gene call data.
            if not families[family_id]["family_data"]:
                families[family_id]["family_data"] = {
                    "sample_id": sample_id,
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
                    "v_gene_light": v_gene_light,
                    "j_gene_light": j_gene_light,
                    "cdr1_start_light": cdr1_start_light,
                    "cdr1_end_light": cdr1_end_light,
                    "cdr2_start_light": cdr2_start_light,
                    "cdr2_end_light": cdr2_end_light,
                    "cdr3_start_light": cdr3_start_light,
                    "cdr3_end_light": cdr3_end_light,
                    "light_chain_type": light_chain_type,
                }
                families[family_id]["is_paired"] = is_paired

            # Add parent node if not already present
            if parent not in families[family_id]["nodes"]:
                parent_node_data = {
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
                # Add light chain sequence for paired format
                if is_paired:
                    parent_node_data["sequence_alignment_light"] = parent_sequence_light
                # Capture extra columns as custom node-level fields
                for col in extra_columns:
                    val = row.get(col, "")
                    if val != "":
                        parent_node_data[col] = coerce_csv_value(val)
                families[family_id]["nodes"][parent] = parent_node_data

            # Add child node if not already present
            if child not in families[family_id]["nodes"]:
                child_node_data = {
                    "sequence_id": child,
                    "pcp_index": pcp_index,
                    "multiplicity": sample_count,
                    "timepoint_multiplicities": [],
                    "sequence_alignment": child_sequence,
                    "is_naive": False,
                    "is_leaf": child_is_leaf,
                    "distances": [distance] if distance > 0 else [],
                    "distance": distance,  # Distance from root
                    "length": branch_length,  # Branch length to this node
                }
                # Add light chain sequence for paired format
                if is_paired:
                    child_node_data["sequence_alignment_light"] = child_sequence_light
                # Capture extra columns as custom node-level fields
                for col in extra_columns:
                    val = row.get(col, "")
                    if val != "":
                        child_node_data[col] = coerce_csv_value(val)
                families[family_id]["nodes"][child] = child_node_data
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

    # Post-process to correctly determine is_paired per-family
    # The global is_paired flag is based on CSV columns, but we need to check actual data
    for family_id, family_data in families.items():
        # Check if ANY node in this family has non-empty light chain sequence
        has_light_data = False
        for node_data in family_data["nodes"].values():
            light_seq = node_data.get("sequence_alignment_light", "")
            if light_seq:  # Non-empty light chain sequence
                has_light_data = True
                break

        # Update is_paired based on actual data presence
        family_data["is_paired"] = has_light_data

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
        vprint.error(f"Warning: Failed to parse Newick tree with ETE3: {e}")
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

    # Start with existing PCP data. Deep-copy the inner node dicts so
    # downstream mutations (e.g. setting "distance") don't leak back
    # into the input, which matters when the same base family data is
    # processed against multiple alternate tree topologies.
    merged_nodes = {node_id: node.copy() for node_id, node in pcp_family_data["nodes"].items()}
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
        vprint.error(f"\nWarning: Found {len(disagreements)} disagreement(s) in family {display_id}:")
        for i, disagree in enumerate(disagreements, 1):
            if disagree["type"] == "branch_length":
                vprint.error(f"  {i}. Branch length mismatch for edge {disagree['edge'][0]} -> {disagree['edge'][1]}:")
                vprint.error(f"     Tree: {disagree['tree_value']:.10f}, PCP: {disagree['pcp_value']:.10f}")
            elif disagree["type"] == "missing_edge":
                if disagree["source"] == "PCP has edge not in tree":
                    vprint.error(f"  {i}. Edge {disagree['edge'][0]} -> {disagree['edge'][1]} exists in PCP but not in tree")
                elif disagree["source"] == "Tree has edge not in PCP":
                    vprint.error(f"  {i}. Edge {disagree['edge'][0]} -> {disagree['edge'][1]} exists in tree but not in PCP")
            elif disagree["type"] == "missing_node":
                if disagree["source"] == "PCP has node not in tree":
                    vprint.error(f"  {i}. Node '{disagree['node']}' exists in PCP but not in tree")
                elif disagree["source"] == "Tree has node not in PCP":
                    vprint.error(f"  {i}. Node '{disagree['node']}' exists in tree but not in PCP")

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

    # Return updated family data, preserving is_paired flag
    return {
        "family_data": pcp_family_data["family_data"],
        "nodes": merged_nodes,
        "edges": merged_edges,
        "is_paired": pcp_family_data.get("is_paired", False),
    }


def parse_newick_csv(csv_path: str) -> Dict[Any, List[Dict[str, Any]]]:
    """Parse a CSV of Newick trees into a list-per-family-key mapping.

    Required columns:

    - ``family_name`` (or ``family``) — clonal family identifier.
    - ``newick_tree`` (or ``newick``) — Newick topology string.

    Optional reserved columns (handled here, not captured as clone-level
    extras):

    - ``sample_id`` — when present, the key is
      ``(family_name, sample_id)``; otherwise just ``family_name``.
    - ``tree_id`` — per-tree identifier. Enables multiple alternate
      reconstructions for the same ``(family, sample_id)``. When absent,
      a fallback is synthesized downstream.
    - ``reconstruction_method`` — written to ``tree.reconstruction_method``
      on the output.
    - ``rate_scale_heavy`` / ``rate_scale_light`` — paired-format rate
      scaling factors; promoted to clone-level fields.

    Any other columns are captured as clone-level extras (same chain
    suffix convention as PCP CSV extras — see
    ``_partition_chain_fields``).

    Returns:
        A mapping from family-key to a list of tree-data dicts. The key is
        ``(family_name, sample_id)`` if a ``sample_id`` column is present,
        else just ``family_name``. Each list entry is one row of the CSV
        and has the shape::

            {
                "newick": str,                       # required
                "tree_id": Optional[str],            # only when column present + non-empty
                "reconstruction_method": Optional[str],
                "rate_scale_heavy": float,           # paired only
                "rate_scale_light": float,           # paired only
                "<extra_col>": coerced value,        # any other columns
            }
    """
    newick_trees: Dict[Any, List[Dict[str, Any]]] = {}

    handle, _ = open_file(csv_path, expected_formats=(FORMAT_PCP,))
    with handle as file_handle:
        reader = csv.DictReader(file_handle)

        # Support alternative column names for backwards compatibility
        # Check for newick column (paired format uses "newick", regular uses "newick_tree")
        newick_col = None
        if "newick" in reader.fieldnames:
            newick_col = "newick"
        elif "newick_tree" in reader.fieldnames:
            newick_col = "newick_tree"

        # Check for family column (paired format uses "family", regular uses "family_name")
        family_col = None
        if "family" in reader.fieldnames:
            family_col = "family"
        elif "family_name" in reader.fieldnames:
            family_col = "family_name"

        # Validate we have required columns
        if newick_col is None:
            raise ValueError("Missing required column: 'newick_tree' or 'newick'")
        if family_col is None:
            raise ValueError("Missing required column: 'family_name' or 'family'")

        # Check if sample_id column exists
        has_sample_id = "sample_id" in reader.fieldnames

        # Check for rate scaling columns (paired format)
        has_rate_scale = (
            "rate_scale_heavy" in reader.fieldnames
            or "rate_scale_light" in reader.fieldnames
        )

        # Identify extra columns for clone-level data.
        # Filter out empty/None column names (e.g., unnamed index columns).
        extra_columns = {
            c for c in reader.fieldnames if c and c not in KNOWN_TREE_COLUMNS
        }

        for row in reader:
            family_name = row[family_col]
            newick_tree = row[newick_col]

            tree_data: Dict[str, Any] = {"newick": newick_tree}

            # Reserved tree-level columns
            tree_id_val = row.get("tree_id")
            if tree_id_val:
                tree_data["tree_id"] = tree_id_val

            reconstruction_method_val = row.get("reconstruction_method")
            if reconstruction_method_val:
                tree_data["reconstruction_method"] = reconstruction_method_val

            if has_rate_scale:
                tree_data["rate_scale_heavy"] = (
                    float(row["rate_scale_heavy"]) if row.get("rate_scale_heavy") else 1.0
                )
                tree_data["rate_scale_light"] = (
                    float(row["rate_scale_light"]) if row.get("rate_scale_light") else 1.0
                )

            # Capture extra columns as clone-level fields
            for col in extra_columns:
                val = row.get(col, "")
                if val != "":
                    tree_data[col] = coerce_csv_value(val)

            if has_sample_id:
                key = (family_name, row["sample_id"])
            else:
                key = family_name

            newick_trees.setdefault(key, []).append(tree_data)

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


def align_and_calculate_mutations(
    germline: str,
    leaf: str,
    alignment_method: str = "truncate"
) -> tuple[int, int, str, str]:
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
        tuple: (mutation_count, alignment_length, germline_aligned, leaf_aligned)
            - mutation_count: Number of mismatches (excluding gaps)
            - alignment_length: Length of aligned sequences
            - germline_aligned: Aligned germline sequence
            - leaf_aligned: Aligned leaf sequence

    Examples:
        >>> align_and_calculate_mutations("ATGC", "ATGT", "truncate")
        (1, 4, 'ATGC', 'ATGT')
        >>> align_and_calculate_mutations("ATG", "ATGCC", "truncate")
        (0, 3, 'ATG', 'ATG')
        >>> align_and_calculate_mutations("ATG", "ATGCC", "pad")
        (0, 5, 'ATG..', 'ATGCC')
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

    # Count mutations (mismatches, excluding empty strings and gaps)
    mutations = sum(
        1 for g, l in zip(g_aligned, l_aligned)
        if g != l and g not in ('', '.') and l not in ('', '.')
    )

    return mutations, length, g_aligned, l_aligned


def log_mutation_frequency_debug(
    family_id: str,
    germline_alignment: str,
    debug_info: list,
    skipped_nodes: list,
    total_mut_freq: float,
    total_sequences: int,
    mean_mut_freq: float,
    chain_label: str,
    vprint: 'VerbosePrinter'
):
    """
    Log detailed mutation frequency calculation information for debugging.

    This function extracts the 85+ lines of debug logging from the mutation
    frequency calculation, improving code readability and maintainability.

    Args:
        family_id: Family identifier
        germline_alignment: Germline DNA sequence
        debug_info: List of dicts with per-node debug information
        skipped_nodes: List of dicts with skipped node information
        total_mut_freq: Total weighted mutation frequency
        total_sequences: Total number of sequences (weighted by multiplicity)
        mean_mut_freq: Calculated mean mutation frequency
        chain_label: Label for the chain type (e.g., "HEAVY CHAIN", "LIGHT CHAIN")
        vprint: VerbosePrinter instance for output

    Note:
        This function only produces output at debug verbosity level (3).
    """
    # Early return if not in debug mode
    if vprint.level < 3:
        return

    vprint.debug(f"\n=== DEBUG: mean_mut_freq calculation for family {family_id} ({chain_label}) ===")
    vprint.debug(f"Germline length: {len(germline_alignment)} nt")

    # Translate germline to amino acids
    germline_aa = translate_dna_to_aa(germline_alignment)
    vprint.debug(f"Germline AA length: {len(germline_aa)} aa")

    vprint.debug(f"\nLEAF sequences only ({len(debug_info)} total):")
    for info in debug_info:
        # Check alignment method used
        alignment_marker = ""
        if info.get('was_aligned', False):
            alignment_marker = f" [{info.get('alignment_method', 'ALIGNED').upper()}]"

        # Show original lengths to verify they match the actual sequence data
        orig_leaf = info.get('original_leaf_len', 0)
        orig_germ = info.get('original_germline_len', 0)
        aligned_len = info['seq_length']

        length_info = f"leaf_len={orig_leaf}, germ_len={orig_germ}, aligned_len={aligned_len}"

        vprint.debug(f"\n  Node {info['node']} (type={info['type']}){alignment_marker}: "
                  f"distance={info['distance']:.6f}, "
                  f"mutations={info['num_mutations']:.1f} nt, "
                  f"{length_info}, "
                  f"mut_freq={info['mut_freq']:.6f}, "
                  f"multiplicity={info['multiplicity']}, "
                  f"weighted={info['weighted_contribution']:.6f}")

        # Translate leaf sequence to amino acids (remove padding before translation)
        leaf_seq_original = info['leaf_seq'].rstrip('.')
        germline_seq_original = info['germline_seq'].rstrip('.')

        leaf_aa = translate_dna_to_aa(leaf_seq_original)
        germline_aa_for_this_leaf = translate_dna_to_aa(germline_seq_original)

        # Pad the amino acid sequences with "."
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
    vprint.debug(f"Mean mutation frequency ({chain_label.lower()}): {mean_mut_freq:.6f}")
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


def _get_descendants(node, children_map):
    """
    Get all descendants of a node in a tree.

    Args:
        node: Node ID to get descendants for
        children_map: Dict mapping parent -> list of children

    Returns:
        set: All descendant node IDs
    """
    descendants = set()
    stack = [node]
    while stack:
        current = stack.pop()
        for child in children_map.get(current, []):
            descendants.add(child)
            stack.append(child)
    return descendants


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

    if len(roots) == 0:
        # No edges or cyclic graph - use first node as root
        if nodes:
            roots = {list(nodes.keys())[0]}
        else:
            return ";"

    if len(roots) > 1:
        # Multiple potential roots - prefer "naive" if present, otherwise pick first
        if "naive" in roots:
            root = "naive"
        else:
            # Use the root with the most descendants (largest subtree)
            root = max(roots, key=lambda r: len(_get_descendants(r, children)))
        roots = {root}

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


@dataclass(frozen=True)
class TreeProcessingConfig:
    """Stable-across-iterations settings for per-family-tree processing.

    Bundles the handful of knobs that ``_process_family_tree`` reads
    (metric computation, alignment method, warning policy) so the
    per-call signature stays about the inputs that actually vary
    (family_data, tree_entry, identities).
    """

    compute_metrics: bool = False
    lbi_tau: float = 0.0125
    standardize_names: bool = False
    alignment_method: str = "truncate"
    warn_disagreements: bool = False


def _build_tree_ref(
    *,
    tree_ident: str,
    family_id: str,
    chain: Optional[Literal["heavy", "light"]] = None,
    newick: str,
    csv_tree_id: Optional[str],
    reconstruction_method: Optional[str],
) -> Dict[str, Any]:
    """Build a tree record (``clone.trees[]`` entry or the header of a
    top-level ``trees[]`` entry) with the identifier and semantic columns.

    Args:
        chain: ``"heavy"`` or ``"light"`` for paired data (adds the
            corresponding suffix to ``ident``/``clone_id``/``tree_id``);
            ``None`` for single-chain data (no suffix).

    ``tree_id`` resolves in this order:
    1. ``csv_tree_id`` when the tree CSV supplied a value (paired data
       appends the chain suffix).
    2. Synthesized ``tree-{family_id}`` (paired:
       ``tree-{family_id}-heavy`` / ``-light``).

    ``reconstruction_method`` is only included on the output record when
    the input supplied one — never fabricated.
    """
    suffix = f"-{chain}" if chain is not None else ""

    if csv_tree_id:
        resolved_tree_id = f"{csv_tree_id}{suffix}"
    else:
        resolved_tree_id = f"tree-{family_id}{suffix}"

    record: Dict[str, Any] = {
        "ident": f"{tree_ident}{suffix}",
        "clone_id": f"{family_id}{suffix}",
        "tree_id": resolved_tree_id,
        "newick": newick,
    }
    if reconstruction_method:
        record["reconstruction_method"] = reconstruction_method
    return record


def _process_family_tree(
    *,
    family_data: Dict[str, Any],
    tree_entry: Dict[str, Any],
    dataset_id: str,
    family_id: str,
    clone_ident: str,
    tree_ident: str,
    config: TreeProcessingConfig,
) -> Optional[Tuple[Dict[str, Any], Optional[Dict[str, Any]], Dict[str, Any], Optional[Dict[str, Any]]]]:
    """Process one (family, tree) pair.

    Returns ``(heavy_clone, light_clone, heavy_tree, light_tree)`` where:

    - ``heavy_clone`` is the full clone dict for the heavy chain, with
      ``"trees"`` already containing a single tree reference dict.
    - ``light_clone`` is the equivalent for the light chain, or ``None`` when
      not paired.
    - ``heavy_tree`` is the full top-level tree record (including the
      ``nodes`` array).
    - ``light_tree`` is the equivalent for light, or ``None`` when not paired.

    Returns ``None`` when the tree should be skipped (e.g. malformed root
    node, missing germline).

    ``family_data`` is expected to already be deep-copied by the caller;
    this helper may mutate it.
    """
    family_meta = family_data.get("family_data", {})
    original_sample_id = family_meta.get("sample_id")

    newick = tree_entry.get("newick") or None
    csv_tree_id = tree_entry.get("tree_id")
    reconstruction_method = tree_entry.get("reconstruction_method")
    rate_scale_heavy = tree_entry.get("rate_scale_heavy", 1.0)
    rate_scale_light = tree_entry.get("rate_scale_light", 1.0)
    # Extras are every column that wasn't a reserved tree-level field.
    extra_tree_fields = {
        k: v for k, v in tree_entry.items()
        if k not in ("newick", "tree_id", "reconstruction_method",
                     "rate_scale_heavy", "rate_scale_light")
    }

    if newick:
        # Use complete tree topology from Newick
        family_data = merge_tree_topology_with_pcp(family_data, newick, config.warn_disagreements, family_id)
    else:
        # Fallback to building tree from PCP edges only
        newick = build_newick_from_edges(family_data["nodes"], family_data["edges"])

    # Store rate scaling in family data for later use
    family_data["family_data"]["rate_scale_heavy"] = rate_scale_heavy
    family_data["family_data"]["rate_scale_light"] = rate_scale_light
    # Store extra tree columns as clone-level custom fields
    if extra_tree_fields:
        family_data["family_data"]["_extra_clone_fields"] = extra_tree_fields

    # Refresh family_meta in case merge_tree_topology_with_pcp replaced it
    family_meta = family_data.get("family_data", {})

    # Check if this family has paired data
    is_paired = family_data.get("is_paired", False)

    # Process nodes - add required fields with rich PCP data
    # For paired data, we'll create TWO sets of nodes (heavy and light)
    processed_nodes_heavy = {}
    processed_nodes_light = {} if is_paired else None

    for node_id, node_data in family_data["nodes"].items():
        # Determine node type based on PCP metadata
        if node_data.get("is_naive", False):
            node_type = "root"
        elif node_data.get("is_leaf", False):
            node_type = "leaf"
        else:
            # This is an internal/ancestral node (Node1, Node2, etc.)
            node_type = "internal"

        # Create heavy chain node
        sequence_alignment_heavy = node_data.get("sequence_alignment", "")
        sequence_alignment_heavy_aa = translate_dna_to_aa(sequence_alignment_heavy)

        processed_node_heavy = {
            "sequence_id": node_id,
            "sequence_alignment": sequence_alignment_heavy,
            "sequence_alignment_aa": sequence_alignment_heavy_aa,
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
        # Carry through extra PCP columns as custom node-level fields
        # Partition by chain suffix for paired data
        _known_node_keys = {
            "sequence_id", "sequence_alignment", "sequence_alignment_aa",
            "sequence_alignment_light", "multiplicity", "cluster_multiplicity",
            "timepoint_multiplicities", "type", "parent", "distance", "length",
            "lbi", "lbr", "affinity", "scaled_affinity",
            "is_naive", "is_leaf", "distances",
        }
        extra_node_raw = {
            k: v for k, v in node_data.items()
            if k not in _known_node_keys and k not in processed_node_heavy
        }
        node_shared, node_heavy, node_light = _partition_chain_fields(extra_node_raw)
        for k, v in node_shared.items():
            processed_node_heavy[k] = v
        for k, v in node_heavy.items():
            processed_node_heavy[k] = v

        processed_nodes_heavy[node_id] = processed_node_heavy

        # Create light chain node (if paired data)
        if is_paired:
            sequence_alignment_light = node_data.get("sequence_alignment_light", "")
            sequence_alignment_light_aa = translate_dna_to_aa(sequence_alignment_light)

            processed_node_light = {
                "sequence_id": node_id,
                "sequence_alignment": sequence_alignment_light,  # Use light chain as main sequence
                "sequence_alignment_aa": sequence_alignment_light_aa,
                "multiplicity": node_data.get("multiplicity", 0),
                "cluster_multiplicity": None,
                "timepoint_multiplicities": node_data.get(
                    "timepoint_multiplicities", []
                ),
                "type": node_type,
                "parent": None,  # Will be set from edges below
                "distance": node_data.get("distance", 0.0),
                "length": node_data.get("length", 0.0),
                "lbi": None,
                "lbr": None,
                "affinity": None,
                "scaled_affinity": None,
            }
            # Add extra node fields: shared + light-only
            for k, v in node_shared.items():
                processed_node_light[k] = v
            for k, v in node_light.items():
                processed_node_light[k] = v

            processed_nodes_light[node_id] = processed_node_light

    # For backward compatibility, keep processed_nodes pointing to heavy chain
    processed_nodes = processed_nodes_heavy

    # Set parent field based on edges (same topology for heavy and light)
    # First, find the true root (node that doesn't appear as a child in any edge)
    all_children = {child for _, child, _ in family_data["edges"]}
    all_parents = {parent for parent, _, _ in family_data["edges"]}
    potential_roots = all_parents - all_children

    # Determine the root node
    if potential_roots:
        tree_root = potential_roots.pop()
    else:
        # If no clear root, use "naive" if present, otherwise use first node
        tree_root = "naive" if "naive" in processed_nodes_heavy else list(processed_nodes_heavy.keys())[0]

    # Set parent relationships for heavy chain nodes
    for parent_id, child_id, edge_length in family_data["edges"]:
        if child_id in processed_nodes_heavy and child_id != tree_root:
            processed_nodes_heavy[child_id]["parent"] = parent_id

    # Ensure the root node has no parent
    if tree_root in processed_nodes_heavy:
        processed_nodes_heavy[tree_root]["parent"] = None

    # Set parent relationships for light chain nodes (if paired)
    if is_paired:
        for parent_id, child_id, edge_length in family_data["edges"]:
            if child_id in processed_nodes_light and child_id != tree_root:
                processed_nodes_light[child_id]["parent"] = parent_id

        if tree_root in processed_nodes_light:
            processed_nodes_light[tree_root]["parent"] = None

    # Calculate cluster multiplicity for heavy chain (always computed)
    vprint.verbose(f"  Computing cluster multiplicity for family {family_id}")
    cluster_mult_values = compute_cluster_multiplicity_for_tree(processed_nodes_heavy, family_data["edges"], tree_root)
    for node_id in processed_nodes_heavy:
        processed_nodes_heavy[node_id]["cluster_multiplicity"] = cluster_mult_values.get(node_id, 0)

    # Calculate cluster multiplicity for light chain (if paired)
    if is_paired:
        cluster_mult_values_light = compute_cluster_multiplicity_for_tree(processed_nodes_light, family_data["edges"], tree_root)
        for node_id in processed_nodes_light:
            processed_nodes_light[node_id]["cluster_multiplicity"] = cluster_mult_values_light.get(node_id, 0)

    # Calculate phylogenetic metrics if requested (computed for both heavy and light)
    if config.compute_metrics:
        vprint.verbose(f"  Computing metrics for family {family_id} (tau={config.lbi_tau})")
        compute_tree_metrics(
            processed_nodes_heavy, family_data["edges"], tree_root, tau=config.lbi_tau
        )
        if is_paired:
            compute_tree_metrics(
                processed_nodes_light, family_data["edges"], tree_root, tau=config.lbi_tau
            )

    # Standardize node names if requested (apply same mapping to heavy and light)
    if config.standardize_names:
        # Create name mapping: old_name -> new_name (same for heavy and light chains)
        name_mapping = {}
        internal_counter = 1
        leaf_counter = 1

        # First pass: create mapping (based on heavy chain node structure)
        for node_id, node_data in processed_nodes_heavy.items():
            node_type = node_data.get("type")
            if node_type == "root":
                name_mapping[node_id] = "naive"
            elif node_type == "leaf":
                name_mapping[node_id] = f"Leaf{leaf_counter}"
                leaf_counter += 1
            else:  # internal
                name_mapping[node_id] = f"Node{internal_counter}"
                internal_counter += 1

        # Second pass: rename heavy chain nodes and update parent references
        renamed_nodes_heavy = {}
        for old_name, node_data in processed_nodes_heavy.items():
            new_name = name_mapping[old_name]
            # Update parent reference to use new name
            if node_data["parent"] and node_data["parent"] in name_mapping:
                node_data["parent"] = name_mapping[node_data["parent"]]
            # Update sequence_id to new name
            node_data["sequence_id"] = new_name
            # Store under new name
            renamed_nodes_heavy[new_name] = node_data

        processed_nodes_heavy = renamed_nodes_heavy

        # Rename light chain nodes using same mapping (if paired)
        if is_paired:
            renamed_nodes_light = {}
            for old_name, node_data in processed_nodes_light.items():
                new_name = name_mapping[old_name]
                if node_data["parent"] and node_data["parent"] in name_mapping:
                    node_data["parent"] = name_mapping[node_data["parent"]]
                node_data["sequence_id"] = new_name
                renamed_nodes_light[new_name] = node_data

            processed_nodes_light = renamed_nodes_light

        # Update tree root and processed_nodes pointer
        tree_root = name_mapping.get(tree_root, tree_root)
        processed_nodes = processed_nodes_heavy

    # Extract family-level immunological data
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
        vprint.print(
            f"WARNING: Family {family_id} missing v_gene_start position, defaulting to 0. "
            "Gene region visualization may be incorrect.",
            min_level=2
        )
        v_alignment_start = 0

    v_alignment_end = family_meta.get("v_gene_end")
    if v_alignment_end is None:
        vprint.print(
            f"WARNING: Family {family_id} missing v_gene_end position, defaulting to 0. "
            "Gene region visualization may be incorrect.",
            min_level=2
        )
        v_alignment_end = 0

    # D gene alignment positions
    d_alignment_start = family_meta.get("d_gene_start")
    if d_alignment_start is None:
        vprint.print(
            f"WARNING: Family {family_id} missing d_gene_start position, defaulting to 0. "
            "Gene region visualization may be incorrect.",
            min_level=2
        )
        d_alignment_start = 0

    d_alignment_end = family_meta.get("d_gene_end")
    if d_alignment_end is None:
        vprint.print(
            f"WARNING: Family {family_id} missing d_gene_end position, defaulting to 0. "
            "Gene region visualization may be incorrect.",
            min_level=2
        )
        d_alignment_end = 0

    # J gene alignment positions
    j_alignment_start = family_meta.get("j_gene_start")
    if j_alignment_start is None:
        vprint.print(
            f"WARNING: Family {family_id} missing j_gene_start position, defaulting to 0. "
            "Gene region visualization may be incorrect.",
            min_level=2
        )
        j_alignment_start = 0

    j_alignment_end = family_meta.get("j_gene_end")
    if j_alignment_end is None:
        vprint.print(
            f"WARNING: Family {family_id} missing j_gene_end position, defaulting to 0. "
            "Gene region visualization may be incorrect.",
            min_level=2
        )
        j_alignment_end = 0

    junction_start = cdr3_start
    junction_length = (cdr3_end - cdr3_start) if (cdr3_end > cdr3_start) else 0

    # Extract light chain data for paired format
    v_call_light = family_meta.get("v_gene_light", "") if is_paired else ""
    j_call_light = family_meta.get("j_gene_light", "") if is_paired else ""
    light_chain_type = family_meta.get("light_chain_type", "") if is_paired else ""

    # Light chain CDR positions
    cdr1_start_light = family_meta.get("cdr1_start_light", 0) if is_paired else 0
    cdr1_end_light = family_meta.get("cdr1_end_light", 0) if is_paired else 0
    cdr2_start_light = family_meta.get("cdr2_start_light", 0) if is_paired else 0
    cdr2_end_light = family_meta.get("cdr2_end_light", 0) if is_paired else 0
    cdr3_start_light = family_meta.get("cdr3_start_light", 0) if is_paired else 0
    cdr3_end_light = family_meta.get("cdr3_end_light", 0) if is_paired else 0

    junction_start_light = cdr3_start_light
    junction_length_light = (cdr3_end_light - cdr3_start_light) if (cdr3_end_light > cdr3_start_light) else 0

    # Rate scaling factors (from trees.csv)
    rate_scale_heavy = family_meta.get("rate_scale_heavy", 1.0)
    rate_scale_light = family_meta.get("rate_scale_light", 1.0) if is_paired else 1.0

    # Get germline sequence from naive node (needed for mean_mut_freq calculation)
    # Validate that exactly one root node exists with a valid sequence
    root_nodes_heavy = [n for n in processed_nodes_heavy.values() if n.get("type") == "root"]
    if len(root_nodes_heavy) != 1:
        vprint.print(
            f"WARNING: Family {family_id} has {len(root_nodes_heavy)} root nodes (expected 1). "
            "This may indicate malformed data. Skipping this tree.",
            min_level=1
        )
        return None

    germline_alignment = root_nodes_heavy[0].get("sequence_alignment", "")
    if not germline_alignment:
        vprint.print(
            f"WARNING: Family {family_id} root node missing sequence_alignment. "
            "Cannot calculate mutation frequency. Skipping this tree.",
            min_level=1
        )
        return None

    # Get light chain germline (if paired)
    germline_alignment_light = ""
    if is_paired:
        root_nodes_light = [n for n in processed_nodes_light.values() if n.get("type") == "root"]
        if len(root_nodes_light) != 1:
            vprint.print(
                f"WARNING: Family {family_id} light chain has {len(root_nodes_light)} root nodes (expected 1). "
                "This may indicate malformed paired data. Skipping this tree.",
                min_level=1
            )
            return None

        germline_alignment_light = root_nodes_light[0].get("sequence_alignment", "")
        if not germline_alignment_light:
            vprint.print(
                f"WARNING: Family {family_id} light chain root node missing sequence_alignment. "
                "Cannot calculate mutation frequency for paired data. Skipping this tree.",
                min_level=1
            )
            return None

    # Calculate mean mutation frequency for HEAVY CHAIN from observed leaf sequences only
    # mean_mut_freq = average(mutations_per_site) across all leaf nodes, weighted by multiplicity
    # Count actual mutations by comparing leaf sequence to germline sequence
    total_mut_freq_heavy = 0.0
    total_sequences_heavy = 0
    germline_length = len(germline_alignment) if germline_alignment else 0

    # DEBUG: Print calculation details for all sequences
    debug_info_heavy = []
    skipped_nodes_heavy = []

    for node_id, node_data in processed_nodes_heavy.items():
        node_type = node_data.get("type")
        multiplicity = node_data.get("multiplicity", 0)
        # Only count LEAF nodes with observed sequences (type="leaf" and multiplicity > 0)
        # Skip internal nodes (type="internal") and root node (type="root")
        if node_type == "leaf" and multiplicity > 0:
            leaf_sequence = node_data.get("sequence_alignment", "")

            # Count mutations by comparing to germline
            if germline_alignment and leaf_sequence:
                # Use helper function to align and calculate mutations
                num_mutations, seq_length, germline_aligned, leaf_aligned = align_and_calculate_mutations(
                    germline_alignment, leaf_sequence, config.alignment_method
                )

                # Calculate mutation frequency
                mut_freq = num_mutations / seq_length if seq_length > 0 else 0.0

                total_mut_freq_heavy += mut_freq * multiplicity
                total_sequences_heavy += multiplicity

                # Collect mutation positions for display
                mutation_positions = []
                for pos, (g, l) in enumerate(zip(germline_aligned, leaf_aligned)):
                    if g != l and g not in ('', '.') and l not in ('', '.'):
                        mutation_positions.append({
                            'pos': pos,
                            'germline': g,
                            'leaf': l
                        })

                # Collect debug info
                debug_info_heavy.append({
                    'node': node_id,
                    'type': node_type,
                    'distance': node_data.get("distance", 0.0),
                    'num_mutations': num_mutations,
                    'seq_length': seq_length,
                    'original_leaf_len': len(leaf_sequence),
                    'original_germline_len': len(germline_alignment),
                    'mut_freq': mut_freq,
                    'multiplicity': multiplicity,
                    'weighted_contribution': mut_freq * multiplicity,
                    'germline_seq': germline_aligned,
                    'leaf_seq': leaf_aligned,
                    'mutations': mutation_positions,
                    'was_aligned': len(germline_alignment) != len(leaf_sequence),
                    'alignment_method': config.alignment_method
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

                skipped_nodes_heavy.append({
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

            skipped_nodes_heavy.append({
                'node': node_id,
                'type': node_type,
                'multiplicity': multiplicity,
                'reason': reason
            })

    mean_mut_freq = total_mut_freq_heavy / total_sequences_heavy if total_sequences_heavy > 0 else 0.0

    # Log debug information (only at debug verbosity level)
    log_mutation_frequency_debug(
        family_id=family_id,
        germline_alignment=germline_alignment,
        debug_info=debug_info_heavy,
        skipped_nodes=skipped_nodes_heavy,
        total_mut_freq=total_mut_freq_heavy,
        total_sequences=total_sequences_heavy,
        mean_mut_freq=mean_mut_freq,
        chain_label="HEAVY CHAIN",
        vprint=vprint
    )

    # Calculate mean mutation frequency for LIGHT CHAIN (if paired).
    # Routes through align_and_calculate_mutations like the heavy chain
    # path above so both chains share the same alignment semantics —
    # including gap-skipping behavior that the inlined version used to
    # diverge on.
    mean_mut_freq_light = 0.0
    if is_paired:
        total_mut_freq_light = 0.0
        total_sequences_light = 0

        for node_id, node_data in processed_nodes_light.items():
            node_type = node_data.get("type")
            multiplicity = node_data.get("multiplicity", 0)

            if node_type == "leaf" and multiplicity > 0:
                leaf_sequence = node_data.get("sequence_alignment", "")

                if germline_alignment_light and leaf_sequence:
                    num_mutations, seq_length, _, _ = align_and_calculate_mutations(
                        germline_alignment_light, leaf_sequence, config.alignment_method
                    )
                    mut_freq = num_mutations / seq_length if seq_length > 0 else 0.0
                    total_mut_freq_light += mut_freq * multiplicity
                    total_sequences_light += multiplicity

        mean_mut_freq_light = total_mut_freq_light / total_sequences_light if total_sequences_light > 0 else 0.0
        vprint.debug(f"\n=== Light chain mean_mut_freq for family {family_id} ===")
        vprint.debug(f"Total mutation frequency (weighted): {total_mut_freq_light:.6f}")
        vprint.debug(f"Total leaf sequences: {total_sequences_light}")
        vprint.debug(f"Mean mutation frequency (light): {mean_mut_freq_light:.6f}")
        vprint.debug(f"  (This means {mean_mut_freq_light*100:.2f}% of positions have mutations on average)")
        vprint.debug(f"===================================================\n")

    # Generate pair_id for paired data (links heavy and light clone entries)
    pair_id = None
    if is_paired:
        # Use family_id as base for pair_id to ensure consistency
        pair_id = f"pair-{family_id}"

    # Create heavy chain clone
    clone_heavy = {
        "clone_id": family_id if not is_paired else f"{family_id}-heavy",
        "ident": clone_ident if not is_paired else f"{clone_ident}-heavy",
        "dataset_id": dataset_id,
        "sample_id": original_sample_id,
        "unique_seqs_count": len(processed_nodes_heavy),
        "total_read_count": sum(
            n.get("multiplicity", 0) for n in processed_nodes_heavy.values()
        ),
        "mean_mut_freq": mean_mut_freq,
        # Heavy chain alignment positions
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
        # Heavy chain gene calls
        "v_call": v_call,
        "d_call": d_call,
        "j_call": j_call,
        "d_alignment_start": d_alignment_start,
        "d_alignment_end": d_alignment_end,
        "germline_alignment": germline_alignment,
        "has_seed": False,
        "trees": [
            _build_tree_ref(
                tree_ident=tree_ident,
                family_id=family_id,
                chain="heavy" if is_paired else None,
                newick=newick,
                csv_tree_id=csv_tree_id,
                reconstruction_method=reconstruction_method,
            )
        ],
        # Denormalized sample reference for webapp convenience
        "sample": {
            "ident": clone_ident if not is_paired else f"{clone_ident}-heavy",
            "locus": infer_locus_from_v_gene(v_call),
            "sample_id": original_sample_id,
        },
    }

    # Add paired data fields
    if is_paired:
        clone_heavy["is_paired"] = True
        clone_heavy["pair_id"] = pair_id

    # Add extra clone-level fields from tree CSV
    # Partition by chain suffix: _heavy → heavy only, _light → light only, neither → shared
    extra_clone_raw = family_meta.get("_extra_clone_fields", {})
    extra_shared, extra_heavy, extra_light = _partition_chain_fields(extra_clone_raw)
    for k, v in extra_shared.items():
        clone_heavy[k] = v
    for k, v in extra_heavy.items():
        clone_heavy[k] = v

    # Create light chain clone (if paired data)
    clone_light: Optional[Dict[str, Any]] = None
    if is_paired:
        # Determine light chain locus from light_chain_type. Leave
        # unset when light_chain_type isn't recognized — the webapp
        # renders its own marker for absent fields.
        if light_chain_type.lower() == "kappa":
            light_locus = "igk"
        elif light_chain_type.lower() == "lambda":
            light_locus = "igl"
        else:
            light_locus = None

        clone_light = {
            "clone_id": f"{family_id}-light",
            "ident": f"{clone_ident}-light",
            "dataset_id": dataset_id,
            "sample_id": original_sample_id,
            "unique_seqs_count": len(processed_nodes_light),
            "total_read_count": sum(
                n.get("multiplicity", 0) for n in processed_nodes_light.values()
            ),
            "mean_mut_freq": mean_mut_freq_light,  # Calculated separately for light chain
            # Light chain alignment positions
            "v_alignment_start": 0,  # Not typically provided for light chain
            "v_alignment_end": 0,
            "j_alignment_start": 0,
            "j_alignment_end": 0,
            "cdr1_alignment_start": cdr1_start_light,
            "cdr1_alignment_end": cdr1_end_light,
            "cdr2_alignment_start": cdr2_start_light,
            "cdr2_alignment_end": cdr2_end_light,
            "junction_start": junction_start_light,
            "junction_length": junction_length_light,
            # Light chain gene calls (no D gene)
            "v_call": v_call_light,
            "d_call": "",  # Light chains don't have D gene
            "j_call": j_call_light,
            "d_alignment_start": 0,
            "d_alignment_end": 0,
            "germline_alignment": germline_alignment_light,
            "has_seed": False,
            "trees": [
                _build_tree_ref(
                    tree_ident=tree_ident,
                    family_id=family_id,
                    chain="light",
                    newick=newick,  # Same topology, different sequences
                    csv_tree_id=csv_tree_id,
                    reconstruction_method=reconstruction_method,
                )
            ],
            "sample": {
                "ident": f"{clone_ident}-light",
                "locus": light_locus,  # "igk" or "igl" based on light_chain_type
                "sample_id": original_sample_id,
            },
            "is_paired": True,
            "pair_id": pair_id,
        }

        # Add extra clone-level fields from tree CSV (shared + light-only)
        for k, v in extra_shared.items():
            clone_light[k] = v
        for k, v in extra_light.items():
            clone_light[k] = v

    # Convert heavy chain nodes to array format (required by webapp)
    nodes_array_heavy = []
    for node_id, node_data in processed_nodes_heavy.items():
        nodes_array_heavy.append(node_data)

    # Create heavy chain tree
    tree_heavy = {
        **_build_tree_ref(
            tree_ident=tree_ident,
            family_id=family_id,
            chain="heavy" if is_paired else None,
            newick=newick,
            csv_tree_id=csv_tree_id,
            reconstruction_method=reconstruction_method,
        ),
        "nodes": nodes_array_heavy,
    }

    # Create light chain tree (if paired)
    tree_light: Optional[Dict[str, Any]] = None
    if is_paired:
        nodes_array_light = []
        for node_id, node_data in processed_nodes_light.items():
            nodes_array_light.append(node_data)

        tree_light = {
            **_build_tree_ref(
                tree_ident=tree_ident,
                family_id=family_id,
                chain="light",
                newick=newick,  # Same topology, different sequences in nodes
                csv_tree_id=csv_tree_id,
                reconstruction_method=reconstruction_method,
            ),
            "nodes": nodes_array_light,
        }

    return clone_heavy, clone_light, tree_heavy, tree_light


def process_pcp_to_olmsted(
    pcp_families: Dict[str, Any],
    newick_trees: Optional[Dict[str, Any]] = None,
    minter: Optional[IdentMinter] = None,
    warn_disagreements: bool = False,
    compute_metrics: bool = False,
    lbi_tau: float = 0.0125,
    standardize_names: bool = False,
    alignment_method: str = "truncate",
    name: Optional[str] = None,
    verbosity: int = 1,
    custom_fields: Optional[List[Dict[str, Any]]] = None,
) -> Tuple[List[OlmstedDataset], Dict[str, List[OlmstedClone]], List[OlmstedTree]]:
    """
    Convert PCP format data to Olmsted format.

    Args:
        pcp_families: dict from parse_pcp_csv
        newick_trees: dict from parse_newick_csv (optional)
        minter: IdentMinter for generating ``{datatype}-{uuid}`` identifiers
            (defaults to a random, non-deterministic minter)
        warn_disagreements: If True, print warnings when tree and PCP data disagree
        compute_metrics: If True, compute all phylogenetic metrics (LBI, LBR, affinity, scaled_affinity, mean_mut_freq)
        lbi_tau: Time scale parameter for LBI calculation (default: 0.0125)
        standardize_names: If True, rename nodes to naive/Node1/Node2.../Leaf1/Leaf2...
        alignment_method: Method for sequence alignment ("truncate" or "pad", default: "truncate")
        name: Optional name for the dataset (default: None)
        verbosity: Verbosity level (0=quiet, 1=normal, 2=verbose, 3=debug)

    Returns:
        Tuple of (datasets, clones_dict, trees) with proper Olmsted types
    """
    # Set global verbosity
    set_verbosity(verbosity)

    if minter is None:
        minter = IdentMinter()

    # Two "dataset" mints on purpose: dataset_id is the semantic foreign
    # key the webapp cross-references (used as the clones_dict key below);
    # dataset_ident is the internal primary-key-shaped ident on the
    # dataset record. Both share the same {datatype}-{uuid} convention
    # but are minted as separate uuids so they can diverge if future
    # input formats supply a user-meaningful dataset_id.
    dataset_id = minter.mint("dataset")
    dataset_ident = minter.mint("dataset")

    # Per-run settings are stable across family/tree iterations — bundle
    # them once so the per-call signature stays focused on the varying
    # inputs (family_data, tree_entry, identities).
    tree_config = TreeProcessingConfig(
        compute_metrics=compute_metrics,
        lbi_tau=lbi_tau,
        standardize_names=standardize_names,
        alignment_method=alignment_method,
        warn_disagreements=warn_disagreements,
    )

    datasets = []
    clones_dict = {dataset_id: []}  # Clones array indexed by dataset_id
    trees = []

    # Create dataset. PCP input has no native subject or timepoint concept,
    # so those collections start empty rather than being fabricated with
    # placeholder records (subject_id="pcp-subject", timepoint="merged").
    # The webapp is expected to render "<unspecified>" for unset reference
    # fields.
    dataset = {
        "ident": dataset_ident,
        "dataset_id": dataset_id,
        "schema_version": SCHEMA_VERSION,
        "build": {"commit": "pcp-import", "time": ""},
        "subjects": [],
        "samples": [],
        "seeds": [],
        "clone_count": len(pcp_families),
        "subjects_count": 0,
        "timepoints_count": 0,
    }

    # Add name if provided
    if name:
        dataset["name"] = name

    # Process each family with progress bar
    family_items = list(pcp_families.items())
    with tqdm(family_items, desc="Processing families", unit="family", disable=len(family_items) == 1) as pbar:
        for family_idx, (family_id, family_data) in enumerate(pbar):
            clone_ident = minter.mint("clone")

            # Get sample_id from family data. sample_id is a required PCP
            # CSV column, so family_meta always has it under normal input.
            family_meta = family_data.get("family_data", {})
            original_sample_id = family_meta.get("sample_id")

            # Create sample if not already present
            sample_exists = any(
                s["sample_id"] == original_sample_id for s in dataset["samples"]
            )
            if not sample_exists:
                # Infer locus from V gene call
                v_gene = family_meta.get("v_gene", "")
                locus = infer_locus_from_v_gene(v_gene)

                dataset["samples"].append(
                    {
                        "ident": minter.mint("sample"),
                        "sample_id": original_sample_id,
                        "locus": locus,
                    }
                )

            # Look up tree entries for this family. parse_newick_csv returns a
            # list-per-key; multiple entries mean alternate reconstructions of
            # the same clonal family.
            tree_entries: List[Dict[str, Any]] = []
            if newick_trees:
                composite_key = (family_id, original_sample_id)
                if composite_key in newick_trees:
                    tree_entries = newick_trees[composite_key]
                elif family_id in newick_trees:
                    tree_entries = newick_trees[family_id]

            # Fall back to a single "empty" tree entry when no CSV trees are
            # supplied; the helper handles the build_newick_from_edges
            # fallback internally.
            tree_entries_normalized = tree_entries if tree_entries else [{}]

            # The helper mutates family_data during topology merge + node
            # processing, so take a pristine copy and deep-copy per tree
            # iteration.
            family_data_pristine = copy.deepcopy(family_data)

            canonical_heavy_clone: Optional[Dict[str, Any]] = None
            canonical_light_clone: Optional[Dict[str, Any]] = None
            for tree_idx, tree_entry in enumerate(tree_entries_normalized):
                family_data_this = copy.deepcopy(family_data_pristine)
                per_tree_ident = minter.mint("tree")

                result = _process_family_tree(
                    family_data=family_data_this,
                    tree_entry=tree_entry,
                    dataset_id=dataset_id,
                    family_id=family_id,
                    clone_ident=clone_ident,
                    tree_ident=per_tree_ident,
                    config=tree_config,
                )
                if result is None:
                    # Helper decided to skip this tree (e.g. missing germline).
                    # Keep processing remaining trees for the family.
                    continue

                heavy_clone, light_clone, heavy_tree, light_tree = result

                if canonical_heavy_clone is None:
                    # First successful tree for this family: use the helper's
                    # clone dicts as canonical (already have one tree ref
                    # inline).
                    canonical_heavy_clone = heavy_clone
                    clones_dict[dataset_id].append(canonical_heavy_clone)
                    if light_clone is not None:
                        canonical_light_clone = light_clone
                        clones_dict[dataset_id].append(canonical_light_clone)
                else:
                    # Subsequent alternate reconstruction: extend the existing
                    # clone's trees[] with just the reference (no nodes).
                    canonical_heavy_clone["trees"].append(
                        {k: v for k, v in heavy_tree.items() if k != "nodes"}
                    )
                    if canonical_light_clone is not None and light_tree is not None:
                        canonical_light_clone["trees"].append(
                            {k: v for k, v in light_tree.items() if k != "nodes"}
                        )

                trees.append(heavy_tree)
                if light_tree is not None:
                    trees.append(light_tree)

    # Generate field_metadata (uses generate_default_config when no config provided)
    dataset_clones = clones_dict.get(dataset_id, [])
    dataset["field_metadata"] = tag_field_metadata(dataset_clones, trees, custom_fields)

    datasets.append(dataset)
    return datasets, clones_dict, trees


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

    # Set up identifier minter (deterministic if seed provided)
    minter = IdentMinter(seed=args.seed)

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
            pcp_families, newick_trees, minter, args.warn_disagreements,
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
