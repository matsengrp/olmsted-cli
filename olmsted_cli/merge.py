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

from .constants import OLMSTED_REQUIRED_TOP_LEVEL_KEYS
from .merge_mutations import apply_mutations_csv
from .process_utils import (
    add_verbosity_args,
    check_output_id_uniqueness,
    resolve_verbosity,
    retag_datasets_field_metadata,
)
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
        "--mutations-use-depth",
        action="store_true",
        help="Use an optional 'depth' column in the mutations CSV to extend the "
        "match key to (site, parent_aa, child_aa, depth). Ignored when the "
        "CSV has a node-name column (depth becomes an integrity check instead). "
        "Opt-in because depth arithmetic depends on the upstream rooting "
        "convention, which the CLI can't infer.",
    )
    parser.add_argument(
        "--mutations-allow-mismatch",
        action="store_true",
        help="Proceed past integrity mismatches between the mutations CSV "
        "and the tree's derived mutations. By default the command fails "
        "when any CSV row matches a (node, site) but its parent_aa, "
        "child_aa, or depth disagree with what the tree derives at that "
        "position; mismatched rows are always skipped regardless of this "
        "flag. Use only after investigating the disagreement.",
    )
    parser.add_argument(
        "--mutations-listed-only",
        action="store_true",
        help="Treat the mutations CSV as authoritative: on trees whose "
        "clone_id matches a family in the CSV, drop any derived "
        "mutations that don't appear in the CSV. Trees whose family is "
        "absent from the CSV pass through untouched.",
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
    parser.add_argument(
        "--allow-duplicate-ids",
        action="store_true",
        help="Downgrade duplicate-*_id errors in the input file to warnings "
        "and pass the data through unchanged. By default, merge fails when "
        "dataset_id, clone_id, tree_id, sample_id, or subject_id collide "
        "within their natural uniqueness scope.",
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

    missing = [k for k in OLMSTED_REQUIRED_TOP_LEVEL_KEYS if k not in data]
    if missing:
        vprint.error(
            f"Error: Input does not appear to be Olmsted JSON "
            f"(missing top-level keys: {missing})"
        )
        sys.exit(1)

    try:
        stats = apply_mutations_csv(
            args.mutations,
            data["trees"],
            use_depth=args.mutations_use_depth,
            allow_mismatch=args.mutations_allow_mismatch,
            only_listed=args.mutations_listed_only,
        )
    except (FileNotFoundError, ValueError) as e:
        vprint.error(f"Error: {e}")
        sys.exit(1)

    retag_datasets_field_metadata(
        data["datasets"],
        data.get("clones", {}),
        data["trees"],
        custom_fields=args.custom_fields,
    )

    try:
        check_output_id_uniqueness(
            data["datasets"],
            data.get("clones", {}),
            allow_duplicates=args.allow_duplicate_ids,
        )
    except ValueError as e:
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
