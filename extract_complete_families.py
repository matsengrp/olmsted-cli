#!/usr/bin/env python3
"""
Extract complete PCP families starting from naive root nodes.
This ensures we get complete, connected subtrees for testing.
"""

import csv
import sys
from collections import defaultdict, deque

def load_pcp_data(csv_path):
    """Load PCP CSV and organize by family."""
    families = defaultdict(list)
    
    with open(csv_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            family_id = row['family']
            families[family_id].append(row)
    
    return families

def find_naive_families(families):
    """Find families that have a naive root node."""
    naive_families = []
    
    for family_id, rows in families.items():
        # Check if this family has a naive root
        has_naive = any(row['parent_is_naive'].lower() == 'true' for row in rows)
        if has_naive:
            naive_families.append(family_id)
    
    return naive_families

def extract_complete_subtree(family_rows, max_nodes=50):
    """Extract complete subtree starting from naive node."""
    # Build parent-child relationships
    edges = []
    nodes = set()
    naive_node = None
    
    for row in family_rows:
        parent = row['parent_name']
        child = row['child_name']
        edges.append((parent, child))
        nodes.add(parent)
        nodes.add(child)
        
        if row['parent_is_naive'].lower() == 'true':
            naive_node = parent
    
    if not naive_node:
        return None, "No naive node found"
    
    # BFS from naive node to find all reachable nodes
    visited = set()
    queue = deque([naive_node])
    reachable_nodes = set()
    
    while queue:
        current = queue.popleft()
        if current in visited:
            continue
        
        visited.add(current)
        reachable_nodes.add(current)
        
        # Find all children of current node
        for parent, child in edges:
            if parent == current and child not in visited:
                queue.append(child)
    
    # If subtree is too large, skip it
    if len(reachable_nodes) > max_nodes:
        return None, f"Subtree too large: {len(reachable_nodes)} nodes"
    
    # Extract only rows that involve reachable nodes
    subtree_rows = []
    for row in family_rows:
        if row['parent_name'] in reachable_nodes and row['child_name'] in reachable_nodes:
            subtree_rows.append(row)
    
    return subtree_rows, None

def main():
    if len(sys.argv) != 4:
        print("Usage: python extract_complete_families.py <input.csv> <output.csv> <num_families>")
        sys.exit(1)
    
    input_path = sys.argv[1]
    output_path = sys.argv[2]
    num_families = int(sys.argv[3])
    
    print(f"Loading PCP data from {input_path}...")
    families = load_pcp_data(input_path)
    print(f"Found {len(families)} total families")
    
    print("Finding families with naive root nodes...")
    naive_families = find_naive_families(families)
    print(f"Found {len(naive_families)} families with naive roots")
    
    # Extract complete subtrees for the first N families
    selected_rows = []
    selected_families = []
    
    for family_id in naive_families[:num_families * 3]:  # Check more families in case some are too large
        if len(selected_families) >= num_families:
            break
            
        family_rows = families[family_id]
        subtree_rows, error = extract_complete_subtree(family_rows)
        
        if subtree_rows:
            selected_rows.extend(subtree_rows)
            selected_families.append(family_id)
            print(f"Selected family {family_id}: {len(subtree_rows)} edges, {len(set(r['parent_name'] for r in subtree_rows) | set(r['child_name'] for r in subtree_rows))} nodes")
        else:
            print(f"Skipped family {family_id}: {error}")
    
    print(f"\nSelected {len(selected_families)} families with {len(selected_rows)} total edges")
    print(f"Families: {selected_families}")
    
    # Write output CSV
    if selected_rows:
        with open(output_path, 'w', newline='') as f:
            # Get fieldnames from first row
            fieldnames = selected_rows[0].keys()
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(selected_rows)
        
        print(f"Wrote {len(selected_rows)} rows to {output_path}")
    else:
        print("No suitable families found!")

if __name__ == '__main__':
    main()