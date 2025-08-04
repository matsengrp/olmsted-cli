#!/usr/bin/env python3
"""Compare PCP tree structures between old and new golden data."""

import json
import os
from pathlib import Path
from collections import defaultdict

def load_json_files(directory):
    """Load all JSON files from a directory."""
    files = {}
    for file_path in Path(directory).glob("*.json"):
        with open(file_path) as f:
            files[file_path.name] = json.load(f)
    return files

def extract_tree_structure(tree_data):
    """Extract structural information from a tree, ignoring UUIDs."""
    structure = {
        "clone_id": tree_data.get("clone_id"),
        "tree_id": tree_data.get("tree_id"),
        "newick": tree_data.get("newick"),
        "node_count": len(tree_data.get("nodes", {})),
        "node_names": sorted(tree_data.get("nodes", {}).keys()),
    }
    
    # Extract node structure details
    nodes_info = []
    for node_name, node_data in sorted(tree_data.get("nodes", {}).items()):
        node_info = {
            "name": node_name,
            "sequence_id": node_data.get("sequence_id"),
            "sequence_length": len(node_data.get("sequence", "")),
            "has_timepoint": "timepoint_id" in node_data,
            "distance_to_parent": node_data.get("distance_to_parent"),
        }
        nodes_info.append(node_info)
    
    structure["nodes_info"] = nodes_info
    return structure

def compare_datasets(old_dir, new_dir):
    """Compare datasets between old and new directories."""
    print(f"Comparing:")
    print(f"  Old: {old_dir}")
    print(f"  New: {new_dir}")
    print()
    
    old_files = load_json_files(old_dir)
    new_files = load_json_files(new_dir)
    
    # Compare datasets.json
    print("=== Comparing datasets.json ===")
    if "datasets.json" in old_files and "datasets.json" in new_files:
        old_dataset = old_files["datasets.json"]
        new_dataset = new_files["datasets.json"]
        
        # Handle both list and dict formats
        if isinstance(old_dataset, list):
            old_dataset = old_dataset[0] if old_dataset else {}
        if isinstance(new_dataset, list):
            new_dataset = new_dataset[0] if new_dataset else {}
        
        # Compare non-UUID fields
        print(f"Old dataset_id: {old_dataset.get('dataset_id', 'N/A')}")
        print(f"New dataset_id: {new_dataset.get('dataset_id', 'N/A')}")
        
        print(f"Old samples count: {len(old_dataset.get('samples', []))}")
        print(f"New samples count: {len(new_dataset.get('samples', []))}")
        
        # Compare sample IDs
        old_sample_ids = [s.get('sample_id') for s in old_dataset.get('samples', [])]
        new_sample_ids = [s.get('sample_id') for s in new_dataset.get('samples', [])]
        print(f"Sample IDs match: {old_sample_ids == new_sample_ids}")
    print()
    
    # Compare tree files
    old_trees = {f: data for f, data in old_files.items() if f.startswith("tree.")}
    new_trees = {f: data for f, data in new_files.items() if f.startswith("tree.")}
    
    print(f"=== Comparing Tree Files ===")
    print(f"Old tree count: {len(old_trees)}")
    print(f"New tree count: {len(new_trees)}")
    print()
    
    # Extract tree structures
    old_structures = {}
    for filename, tree_data in old_trees.items():
        structure = extract_tree_structure(tree_data)
        tree_id = tree_data.get("tree_id", filename)
        old_structures[tree_id] = structure
    
    new_structures = {}
    for filename, tree_data in new_trees.items():
        structure = extract_tree_structure(tree_data)
        tree_id = tree_data.get("tree_id", filename)
        new_structures[tree_id] = structure
    
    # Match trees by structure
    print("=== Matching Trees by Structure ===")
    matched = 0
    unmatched_old = []
    unmatched_new = []
    
    # Try to match each old tree with a new tree
    for old_tree_id, old_struct in old_structures.items():
        found_match = False
        for new_tree_id, new_struct in new_structures.items():
            if old_struct == new_struct:
                print(f"✓ Matched: {old_tree_id} -> {new_tree_id}")
                matched += 1
                found_match = True
                break
        
        if not found_match:
            unmatched_old.append(old_tree_id)
            print(f"✗ No match found for old tree: {old_tree_id}")
            print(f"  Structure: {json.dumps(old_struct, indent=2)}")
    
    # Check for new trees without matches
    for new_tree_id, new_struct in new_structures.items():
        found_match = False
        for old_tree_id, old_struct in old_structures.items():
            if old_struct == new_struct:
                found_match = True
                break
        
        if not found_match:
            unmatched_new.append(new_tree_id)
            print(f"✗ No match found for new tree: {new_tree_id}")
            print(f"  Structure: {json.dumps(new_struct, indent=2)}")
    
    print()
    print(f"=== Summary ===")
    print(f"Matched trees: {matched}")
    print(f"Unmatched old trees: {len(unmatched_old)}")
    print(f"Unmatched new trees: {len(unmatched_new)}")
    
    # Compare clones files
    print()
    print("=== Comparing Clone Files ===")
    old_clones = {f: data for f, data in old_files.items() if f.startswith("clones.")}
    new_clones = {f: data for f, data in new_files.items() if f.startswith("clones.")}
    
    if old_clones and new_clones:
        old_clone_data = list(old_clones.values())[0]
        new_clone_data = list(new_clones.values())[0]
        
        if isinstance(old_clone_data, list) and isinstance(new_clone_data, list):
            print(f"Old clones count: {len(old_clone_data)}")
            print(f"New clones count: {len(new_clone_data)}")
            
            # Compare clone IDs
            old_clone_ids = sorted([c.get("clone_id") for c in old_clone_data])
            new_clone_ids = sorted([c.get("clone_id") for c in new_clone_data])
            
            print(f"Clone IDs match: {old_clone_ids == new_clone_ids}")
            if old_clone_ids != new_clone_ids:
                print(f"  Old: {old_clone_ids}")
                print(f"  New: {new_clone_ids}")
    
    return matched == len(old_trees) and matched == len(new_trees)


if __name__ == "__main__":
    old_golden = "/home/devreckas/Google-Drive/Work/matsen-lab/olmsted/olmsted-cli/example_data/pcp/golden_pcp_data"
    new_golden = "/tmp/new_golden_pcp"
    
    all_match = compare_datasets(old_golden, new_golden)
    
    print()
    if all_match:
        print("✅ All trees match structurally! Safe to replace golden data.")
    else:
        print("❌ Trees do not match! Do not replace golden data.")