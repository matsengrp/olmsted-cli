"""
Merge command for adding external mutation-level data into existing Olmsted JSON.

This command takes an existing Olmsted JSON file and a CSV of mutation-level
annotations (e.g., surprise scores), matches CSV rows to mutations on tree
nodes, and writes an updated Olmsted JSON file with the merged data and
regenerated field_metadata.

Usage:
    olmsted merge -i base.json --mutations scores.csv -o output.json
    olmsted merge -i base.json --mutations scores.csv --in-place
"""

import argparse
import json
import sys
from pathlib import Path

from .merge_mutations import load_mutations_csv, merge_mutations_into_trees
from .process_utils import (
    add_verbosity_args,
    resolve_verbosity,
    tag_field_metadata,
)
from .utils import set_verbosity, vprint


def get_args():
    """Parse command line arguments for the merge command."""
    parser = argparse.ArgumentParser(
        description="Merge external data into an existing Olmsted JSON file",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Merge mutation-level CSV scores into existing Olmsted JSON
    olmsted merge -i base.json --mutations scores.csv -o output.json

    # In-place modification
    olmsted merge -i base.json --mutations scores.csv --in-place
        """,
    )

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input Olmsted JSON file",
    )
    parser.add_argument(
        "--mutations",
        required=True,
        help="Mutations CSV file to merge (columns: family, site, parent_aa, child_aa, ...)",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Output file path (required unless --in-place is used)",
    )
    parser.add_argument(
        "--in-place",
        action="store_true",
        help="Modify the input file in place",
    )
    parser.add_argument(
        "--json-format",
        choices=["pretty", "compact"],
        default="pretty",
        help="JSON output format (default: pretty)",
    )
    add_verbosity_args(parser)

    args = parser.parse_args()
    resolve_verbosity(args)

    if not args.output and not args.in_place:
        parser.error("Either -o/--output or --in-place must be specified")
    if args.output and args.in_place:
        parser.error("Cannot specify both -o/--output and --in-place")

    return args


def main():
    """Main entry point for the merge command."""
    args = get_args()
    set_verbosity(args.verbose)

    # Load input Olmsted JSON
    input_path = Path(args.input)
    if not input_path.exists():
        vprint.error(f"Error: Input file not found: {input_path}")
        sys.exit(1)

    try:
        with open(input_path) as f:
            data = json.load(f)
    except json.JSONDecodeError as e:
        vprint.error(f"Error: Invalid JSON in input file: {e}")
        sys.exit(1)

    if "datasets" not in data or "clones" not in data or "trees" not in data:
        vprint.error(
            "Error: Input does not appear to be Olmsted JSON "
            "(missing 'datasets', 'clones', or 'trees' key)"
        )
        sys.exit(1)

    # Load mutations CSV
    vprint.status(f"Loading mutations CSV: {args.mutations}")
    try:
        mutations_by_family = load_mutations_csv(args.mutations)
    except (FileNotFoundError, ValueError) as e:
        vprint.error(f"Error: {e}")
        sys.exit(1)

    total_csv_rows = sum(len(rows) for rows in mutations_by_family.values())
    vprint.status(
        f"Loaded {total_csv_rows} mutation records across "
        f"{len(mutations_by_family)} families"
    )

    # Merge into trees
    vprint.status("Merging mutations into tree nodes...")
    trees = data["trees"]
    stats = merge_mutations_into_trees(trees, mutations_by_family)
    vprint.status(
        f"Matched {stats.trees_matched} trees, "
        f"merged {stats.mutations_merged} mutation records "
        f"across {stats.nodes_with_mutations} nodes"
    )

    if stats.trees_matched == 0:
        vprint.error(
            "Warning: No trees matched the families in the mutations CSV. "
            "Check that the CSV 'family' column matches clone_id values."
        )

    if stats.unmatched_families:
        sample = stats.unmatched_families[:5]
        vprint.error(
            f"Error: {len(stats.unmatched_families)} families in the mutations CSV "
            f"had no matching clone in the Olmsted JSON (e.g., {sample})"
        )

    if stats.unmatched_mutations:
        vprint.error(
            f"Error: {stats.unmatched_mutations} CSV mutation records in matched "
            f"families had no corresponding derived mutation in any node. "
            f"Run with -v 2 to see per-family details."
        )

    # Regenerate field_metadata for each dataset
    datasets = data.get("datasets", [])
    clones_dict = data.get("clones", {})

    trees_by_clone_id = {}
    for tree in trees:
        clone_id = tree.get("clone_id")
        if clone_id:
            trees_by_clone_id.setdefault(clone_id, []).append(tree)

    for dataset in datasets:
        dataset_id = dataset.get("dataset_id")
        if not dataset_id:
            continue
        dataset_clones = clones_dict.get(dataset_id, [])
        dataset_trees = []
        clone_ids = {c.get("clone_id") for c in dataset_clones if c.get("clone_id")}
        for clone_id in clone_ids:
            dataset_trees.extend(trees_by_clone_id.get(clone_id, []))

        new_field_metadata = tag_field_metadata(dataset_clones, dataset_trees)

        # Merge with existing field_metadata (preserve unrelated entries)
        existing_metadata = dataset.get("field_metadata", {})
        merged = {}
        all_levels = set(
            list(existing_metadata.keys()) + list(new_field_metadata.keys())
        )
        for level in all_levels:
            existing_level = existing_metadata.get(level, {})
            new_level = new_field_metadata.get(level, {})
            merged_level = dict(existing_level)
            merged_level.update(new_level)
            if merged_level:
                merged[level] = merged_level
        dataset["field_metadata"] = merged

    # Write output
    output_path = input_path if args.in_place else Path(args.output)
    indent = 2 if args.json_format == "pretty" else None
    separators = None if args.json_format == "pretty" else (",", ":")
    with open(output_path, "w") as f:
        json.dump(data, f, indent=indent, separators=separators)

    vprint.status(f"Merged data written to: {output_path}")


if __name__ == "__main__":
    main()
