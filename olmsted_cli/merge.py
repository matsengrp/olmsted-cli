"""
Merge command for adding external mutation-level data into existing Olmsted JSON.

This command takes an existing Olmsted JSON file and a CSV of mutation-level
annotations (e.g., surprise scores), matches CSV rows to mutations on tree
nodes, and writes an updated Olmsted JSON file with the merged data and
regenerated field_metadata.

Usage:
    olmsted merge -i base.json --mutations scores.csv -o output.json
    olmsted merge -i base.json --mutations scores.csv --in-place
    olmsted merge -i base.json --mutations scores.csv -c config.yaml -o out.json
"""

import argparse
import json
import sys
from pathlib import Path

from .merge_mutations import apply_mutations_csv
from .process_utils import add_verbosity_args, resolve_verbosity
from .utils import set_verbosity, vprint

# Config keys merge reads from YAML (beyond custom_fields)
_MERGE_CONFIG_KEYS = {"input", "mutations", "output"}


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

    # With custom field declarations from a YAML config
    olmsted merge -i base.json --mutations scores.csv -c config.yaml -o out.json
        """,
    )

    parser.add_argument(
        "-i",
        "--input",
        help="Input Olmsted JSON file (required, or provide in config)",
    )
    parser.add_argument(
        "--mutations",
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
        "-c",
        "--config",
        help="YAML configuration file with custom field declarations",
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

    # Load config and apply values where CLI didn't explicitly set them.
    # Imported here to avoid a circular import at module load time.
    from .process_data import load_config

    custom_fields = None
    if args.config:
        config_dict, custom_fields = load_config(args.config)
        for key in _MERGE_CONFIG_KEYS:
            if key in config_dict and getattr(args, key, None) is None:
                setattr(args, key, config_dict[key])
    args.custom_fields = custom_fields

    if not args.input:
        parser.error(
            "the following arguments are required: -i/--input (or provide in config)"
        )
    if not args.mutations:
        parser.error(
            "the following arguments are required: --mutations (or provide in config)"
        )
    if not args.output and not args.in_place:
        parser.error("Either -o/--output or --in-place must be specified")
    if args.output and args.in_place:
        parser.error("Cannot specify both -o/--output and --in-place")

    return args


def main():
    """Main entry point for the merge command."""
    args = get_args()
    set_verbosity(args.verbose)

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

    try:
        stats = apply_mutations_csv(
            args.mutations,
            data["datasets"],
            data.get("clones", {}),
            data["trees"],
            custom_fields=args.custom_fields,
        )
    except (FileNotFoundError, ValueError) as e:
        vprint.error(f"Error: {e}")
        sys.exit(1)

    # Refuse to overwrite the input file in place when nothing matched.
    if args.in_place and stats is not None and stats.trees_matched == 0:
        vprint.error(
            "Error: Refusing to --in-place overwrite when zero trees matched. "
            "Re-run without --in-place to write to a new file, or verify the "
            "mutations CSV 'family' column matches your clone_id values."
        )
        sys.exit(1)

    output_path = input_path if args.in_place else Path(args.output)
    indent = 2 if args.json_format == "pretty" else None
    separators = None if args.json_format == "pretty" else (",", ":")
    with open(output_path, "w") as f:
        json.dump(data, f, indent=indent, separators=separators)

    vprint.status(f"Merged data written to: {output_path}")


if __name__ == "__main__":
    main()
