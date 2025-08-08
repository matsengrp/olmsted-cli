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
import os
import sys
import traceback
import uuid
from collections import defaultdict

# Import shared utilities from process_data_utils
from .process_utils import (
    SCHEMA_VERSION,
    translate_dna_to_aa,
    validate_output_data,
    write_out,
)

def parse_pcp_csv(csv_path):
    """
    Parse PCP CSV file and return a dict of families with rich immunological data.

    Expected CSV format:
    sample_id,family,parent_name,parent_heavy,child_name,child_heavy,branch_length,
    v_gene_heavy,j_gene_heavy,cdr1_codon_start_heavy,cdr1_codon_end_heavy,
    cdr2_codon_start_heavy,cdr2_codon_end_heavy,cdr3_codon_start_heavy,
    cdr3_codon_end_heavy,parent_is_naive,child_is_leaf

    Returns:
        dict: {family_id: {
            nodes: {node_id: node_data},
            edges: [(parent, child, length)],
            family_data: {v_gene, j_gene, cdr_positions, etc.}
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

            # Store family-level data and sample_id for each family (will be same for all rows of same family)
            families[family_id]["family_data"] = {
                "sample_id": sample_id,  # Store original sample_id for reference
                "v_gene": v_gene,
                "j_gene": j_gene,
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
                    "distances": [],  # Track distances for mutation frequency calculation
                    "distance": 0.0,  # Root node has zero distance
                    "length": 0.0,  # Root node has zero length
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

    return dict(families)


def parse_newick_csv(csv_path):
    """
    Parse CSV file containing Newick trees.

    Expected CSV format:
    family_name,newick_tree

    Returns:
        dict: {family_name: newick_string}
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

        for row in reader:
            family_name = row["family_name"]
            newick_tree = row["newick_tree"]
            newick_trees[family_name] = newick_tree

    return newick_trees


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


def process_pcp_to_olmsted(pcp_families, newick_trees=None, uuid_generator=None):
    """
    Convert PCP format data to Olmsted format.

    Args:
        pcp_families: dict from parse_pcp_csv
        newick_trees: dict from parse_newick_csv (optional)
        uuid_generator: Function to generate UUIDs (defaults to random)

    Returns:
        tuple: (datasets, clones_dict, trees)
    """
    if uuid_generator is None:
        uuid_generator = lambda: str(uuid.uuid4())

    dataset_id = f"pcp-{uuid_generator()}"
    dataset_ident = uuid_generator()

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
        "subjects": [{"ident": uuid_generator(), "subject_id": "pcp-subject"}],
        "samples": [],
        "seeds": [],
        "clone_count": len(pcp_families),
        "subjects_count": 1,
        "timepoints_count": 1,
    }

    # Process each family
    for family_idx, (family_id, family_data) in enumerate(pcp_families.items()):
        clone_ident = uuid_generator()
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
                    "ident": uuid_generator(),
                    "sample_id": original_sample_id,
                    "locus": "igh",  # Default locus
                    "timepoint_id": "merged",
                }
            )

        # Build or use provided Newick tree
        if newick_trees and family_id in newick_trees:
            newick = newick_trees[family_id]
        else:
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
            }
            processed_nodes[node_id] = processed_node

        # Set parent field based on edges
        for parent_id, child_id, edge_length in family_data["edges"]:
            if child_id in processed_nodes:
                processed_nodes[child_id]["parent"] = parent_id

        # Extract family-level immunological data (already extracted above)
        v_call = family_meta.get("v_gene", "")
        j_call = family_meta.get("j_gene", "")

        # Calculate alignment positions from CDR data
        cdr1_start = family_meta.get("cdr1_start", 0)
        cdr2_end = family_meta.get("cdr2_end", 0)
        cdr3_start = family_meta.get("cdr3_start", 0)
        cdr3_end = family_meta.get("cdr3_end", 0)

        # Use CDR positions to estimate V and J alignment positions
        v_alignment_start = cdr1_start if cdr1_start > 0 else 0
        v_alignment_end = cdr2_end if cdr2_end > 0 else 0
        j_alignment_start = cdr3_end if cdr3_end > 0 else 0
        j_alignment_end = (
            j_alignment_start + 50 if j_alignment_start > 0 else 0
        )  # Estimate

        junction_start = cdr3_start
        junction_length = (cdr3_end - cdr3_start) if (cdr3_end > cdr3_start) else 0

        # Calculate mean mutation frequency from distance data
        all_distances = []
        for node_id, node_data in family_data["nodes"].items():
            distances = node_data.get("distances", [])
            all_distances.extend(distances)

        mean_mut_freq = (
            sum(all_distances) / len(all_distances) if all_distances else 0.0
        )
        # Convert to more realistic scale (distances are very small scientific notation)
        mean_mut_freq = mean_mut_freq * 1000000  # Scale up for better visualization

        # Get germline sequence from naive node
        germline_alignment = ""
        for node_id, node_data in processed_nodes.items():
            if node_data.get("type") == "root":
                germline_alignment = node_data.get("sequence_alignment", "")
                break

        # Create clone with rich PCP data
        clone = {
            "clone_id": f"family-{family_idx}",
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
            "junction_start": junction_start,
            "junction_length": junction_length,
            "v_call": v_call,
            "j_call": j_call,
            "d_call": "",  # Not available in PCP format
            "d_alignment_start": 0,
            "d_alignment_end": 0,
            "germline_alignment": germline_alignment,
            "has_seed": False,
            "trees": [
                {
                    "ident": tree_ident,
                    "clone_id": f"family-{family_idx}",
                    "tree_id": f"pcp-tree-{family_idx}",
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
            "tree_id": f"pcp-tree-{family_idx}",
            "clone_id": f"family-{family_idx}",
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
        "--output-dir",
        required=True,
        help="Output directory for processed JSON files",
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose output")
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
    # Removed --output-format option - now only outputs AIRR format
    return parser.parse_args()


def main():
    """Main entry point."""
    args = get_args()

    # Set up deterministic UUID generation if seed is provided
    uuid_counter = 0

    def get_uuid():
        nonlocal uuid_counter
        if args.seed is not None:
            uuid_counter += 1
            return deterministic_uuid(args.seed, uuid_counter)
        else:
            return str(uuid.uuid4())

    try:
        # Parse PCP CSV
        print(f"Processing PCP CSV: {args.input_pcp}")
        if args.seed is not None:
            print(f"Using deterministic UUIDs with seed: {args.seed}")
        pcp_families = parse_pcp_csv(args.input_pcp)
        print(f"Found {len(pcp_families)} families")

        # Parse Newick trees if provided
        newick_trees = None
        if args.input_trees:
            print(f"Processing Newick trees: {args.input_trees}")
            newick_trees = parse_newick_csv(args.input_trees)
            print(f"Found {len(newick_trees)} trees")

        # Convert to Olmsted format
        print("Converting to Olmsted format...")
        datasets, clones_dict, trees = process_pcp_to_olmsted(
            pcp_families, newick_trees, get_uuid
        )

        # Create output directory if needed
        os.makedirs(args.output_dir, exist_ok=True)

        # Only AIRR format output - no need to prepare other formats

        # Validate output data if requested
        if args.validate:
            if not validate_output_data(datasets, clones_dict, trees, args):
                if args.strict_validation:
                    print(
                        "\nExiting due to validation errors (--strict-validation enabled)"
                    )
                    sys.exit(1)

        # Write AIRR format output
        print(f"Writing AIRR format output to {args.output_dir}")

        write_out(datasets, args.output_dir, "datasets.json", args)
        for dataset_id, clones in clones_dict.items():
            write_out(clones, args.output_dir, f"clones.{dataset_id}.json", args)
        for tree in trees:
            write_out(tree, args.output_dir, f"tree.{tree['ident']}.json", args)

        print("Processing complete!")

    except Exception as e:
        print(f"Error: {e}")
        if args.verbose:
            traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
