#!/usr/bin/env python3
"""Split command for olmsted-cli - splits Olmsted JSON files into smaller files."""

import argparse
import json
import sys
from pathlib import Path
import copy
from datetime import datetime

from .data_io import read_olmsted_json
from .identifier import IdentMinter
from .utils import set_verbosity, vprint


def split_consolidated_data(data, max_clones_per_file):
    """
    Split Olmsted JSON data into multiple files based on max clones per file.
    
    Args:
        data: Olmsted JSON data dictionary
        max_clones_per_file: Maximum number of clonal families per output file
        
    Returns:
        list: List of split data dictionaries
    """
    # Validate input data
    required_keys = ["metadata", "datasets", "clones", "trees"]
    missing_keys = [key for key in required_keys if key not in data]
    if missing_keys:
        raise ValueError(f"Not a valid Olmsted JSON format. Missing keys: {missing_keys}")
    
    # Get all clones and trees
    all_clones = []
    all_trees = []
    dataset_to_clones = {}
    
    # Collect all clones by dataset
    for dataset_id, clones_list in data["clones"].items():
        all_clones.extend(clones_list)
        dataset_to_clones[dataset_id] = clones_list
    
    all_trees = data["trees"]
    
    if len(all_clones) <= max_clones_per_file:
        # No need to split
        return [data]
    
    # Create tree lookup by clone_id for efficient mapping
    tree_by_clone = {}
    for tree in all_trees:
        clone_id = tree.get("clone_id")
        if clone_id:
            if clone_id not in tree_by_clone:
                tree_by_clone[clone_id] = []
            tree_by_clone[clone_id].append(tree)
    
    # Split clones into chunks
    clone_chunks = []
    for i in range(0, len(all_clones), max_clones_per_file):
        clone_chunks.append(all_clones[i:i + max_clones_per_file])
    
    # Create split files
    split_files = []
    _split_minter = IdentMinter()

    for chunk_idx, clone_chunk in enumerate(clone_chunks):
        split_data = copy.deepcopy(data)
        
        # Update metadata
        original_name = split_data["metadata"].get("name", "dataset")
        split_data["metadata"]["name"] = f"{original_name}-{chunk_idx + 1}"
        split_data["metadata"]["created_at"] = datetime.now().isoformat() + "+00:00"
        split_data["metadata"]["split_info"] = {
            "original_name": original_name,
            "chunk_number": chunk_idx + 1,
            "total_chunks": len(clone_chunks),
            "max_clones_per_file": max_clones_per_file,
            "clones_in_chunk": len(clone_chunk)
        }
        
        # Create new dataset IDs and update datasets
        new_datasets = []
        old_to_new_dataset_mapping = {}
        
        for dataset in split_data["datasets"]:
            new_dataset_id = f"{dataset['dataset_id']}-{chunk_idx + 1}"
            old_to_new_dataset_mapping[dataset["dataset_id"]] = new_dataset_id
            
            new_dataset = copy.deepcopy(dataset)
            new_dataset["dataset_id"] = new_dataset_id
            new_dataset["ident"] = _split_minter.mint("dataset")
            
            # Count clones for this dataset in this chunk
            dataset_clone_count = sum(1 for clone in clone_chunk 
                                    if clone.get("dataset_id") == dataset["dataset_id"])
            new_dataset["clone_count"] = dataset_clone_count
            
            new_datasets.append(new_dataset)
        
        split_data["datasets"] = new_datasets
        
        # Update clones with new dataset IDs and filter to current chunk
        new_clones = {}
        for old_dataset_id, new_dataset_id in old_to_new_dataset_mapping.items():
            # Get clones for this dataset that are in the current chunk
            dataset_clones_in_chunk = [
                clone for clone in clone_chunk 
                if clone.get("dataset_id") == old_dataset_id
            ]
            
            if dataset_clones_in_chunk:
                # Update dataset_id in clones
                updated_clones = []
                for clone in dataset_clones_in_chunk:
                    updated_clone = copy.deepcopy(clone)
                    updated_clone["dataset_id"] = new_dataset_id
                    updated_clones.append(updated_clone)
                
                new_clones[new_dataset_id] = updated_clones
        
        split_data["clones"] = new_clones
        
        # Update trees - only include trees for clones in this chunk
        chunk_clone_ids = set(clone.get("clone_id") for clone in clone_chunk)
        chunk_trees = []
        
        for tree in all_trees:
            if tree.get("clone_id") in chunk_clone_ids:
                updated_tree = copy.deepcopy(tree)
                # Update tree's dataset reference if needed
                if updated_tree.get("dataset_id") in old_to_new_dataset_mapping:
                    updated_tree["dataset_id"] = old_to_new_dataset_mapping[updated_tree["dataset_id"]]
                chunk_trees.append(updated_tree)
        
        split_data["trees"] = chunk_trees
        
        split_files.append(split_data)
    
    return split_files


def get_args():
    """Parse command line arguments for split command."""
    parser = argparse.ArgumentParser(
        description="Split Olmsted JSON files into smaller files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Split into files with max 100 clones each
  olmsted split -i data.json -o output_dir --max-clones 100
  
  # Split with custom naming
  olmsted split -i data.json -o splits --max-clones 50 --base-name my_dataset
        """,
    )
    
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input Olmsted JSON file to split"
    )
    
    parser.add_argument(
        "-o", "--output-dir",
        required=True,
        help="Output directory for split files"
    )
    
    parser.add_argument(
        "--max-clones",
        type=int,
        default=1000,
        help="Maximum number of clonal families per output file (default: 1000)"
    )
    
    parser.add_argument(
        "--base-name",
        help="Base name for output files (default: derived from input filename)"
    )
    
    parser.add_argument(
        "-v", "--verbose",
        type=int,
        choices=[0, 1, 2, 3],
        default=1,
        help="Verbosity: 0=errors only, 1=normal (default), 2=verbose, 3=debug",
    )
    parser.add_argument(
        "-q", "--quiet",
        action="store_true",
        help="Quiet mode — errors only (equivalent to -v 0)",
    )

    return parser.parse_args()


def main():
    """Main entry point for split command."""
    args = get_args()

    # Handle quiet mode
    if getattr(args, "quiet", False):
        args.verbose = 0
    set_verbosity(args.verbose)

    # Validate input file
    input_path = Path(args.input)
    if not input_path.exists():
        vprint.error(f"Error: Input file not found: {args.input}")
        sys.exit(1)
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Determine base name for output files
    if args.base_name:
        base_name = args.base_name
    else:
        # Remove .json extension and use filename
        base_name = input_path.stem
        if base_name.endswith('.consolidated'):
            base_name = base_name[:-12]  # Remove .consolidated suffix
    
    # Load and parse input data
    try:
        data = read_olmsted_json(input_path)
    except (ValueError, OSError) as e:
        vprint.error(f"Error: {e}")
        sys.exit(1)

    # Split the data
    try:
        split_files = split_consolidated_data(data, args.max_clones)
    except Exception as e:
        vprint.error(f"Error: Failed to split data: {e}")
        sys.exit(1)

    if len(split_files) == 1:
        vprint.status(f"File has {sum(len(clones) for clones in data['clones'].values())} clones, "
              f"which is within the limit of {args.max_clones}. No splitting needed.")
        return
    
    # Write split files
    output_files = []
    for i, split_data in enumerate(split_files):
        chunk_id = i + 1
        output_filename = f"{base_name}.{chunk_id}.json"
        output_path = output_dir / output_filename
        
        try:
            with open(output_path, 'w') as f:
                json.dump(split_data, f, indent=2)
            output_files.append(output_path)
            
            if args.verbose:
                clone_count = sum(len(clones) for clones in split_data['clones'].values())
                tree_count = len(split_data['trees'])
                vprint.status(f"Created {output_filename}: {clone_count} clones, {tree_count} trees")

        except Exception as e:
            vprint.error(f"Error: Failed to write {output_path}: {e}")
            sys.exit(1)

    # Summary
    total_clones = sum(len(clones) for clones in data['clones'].values())
    total_trees = len(data['trees'])

    vprint.status(f"Successfully split {input_path.name} into {len(split_files)} files:")
    vprint.status(f"  Total clones: {total_clones}")
    vprint.status(f"  Total trees: {total_trees}")
    vprint.status(f"  Max clones per file: {args.max_clones}")
    vprint.status(f"  Output directory: {output_dir}")
    vprint.status(f"  Files created: {', '.join(f.name for f in output_files)}")


if __name__ == "__main__":
    main()