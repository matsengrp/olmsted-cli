"""
Tag command for adding field_metadata to existing Olmsted JSON files.

This command introspects an existing Olmsted JSON file, discovers available
fields at each level (clone, node, branch, mutation), and adds a
field_metadata object to each dataset.

Usage:
    olmsted tag -i data.json -o tagged.json
    olmsted tag -i data.json -o tagged.json -c config.yaml
    olmsted tag -i data.json --in-place
"""

import argparse
import json
import sys
from pathlib import Path

# Config keys that tag reads from YAML (beyond custom_fields)
_TAG_CONFIG_KEYS = {"input", "output", "mode"}

from .process_data import load_config
from .process_utils import (
    VerbosePrinter,
    add_verbosity_args,
    resolve_verbosity,
    retag_datasets_field_metadata,
)
from .utils import set_verbosity, vprint


def get_args():
    """Parse command line arguments for the tag command."""
    parser = argparse.ArgumentParser(
        description="Add field_metadata to an existing Olmsted JSON file",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Basic tagging (introspect and add field_metadata)
    olmsted tag -i data.json -o tagged.json

    # With custom field declarations from config
    olmsted tag -i data.json -o tagged.json -c config.yaml

    # In-place modification
    olmsted tag -i data.json --in-place

    # In-place with config
    olmsted tag -i data.json --in-place -c config.yaml
        """,
    )

    parser.add_argument(
        "-i",
        "--input",
        help="Input Olmsted JSON file to tag (required, or provide in config)",
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
        "--mode",
        choices=["add", "overwrite"],
        default="add",
        help="Merge mode: 'add' preserves existing field_metadata entries and adds new ones (default), 'overwrite' replaces field_metadata entirely",
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

    # Load config and apply values where CLI didn't explicitly set them
    custom_fields = None
    if args.config:
        config_dict, custom_fields = load_config(args.config)
        for key in _TAG_CONFIG_KEYS:
            if key in config_dict and getattr(args, key, None) is None:
                setattr(args, key, config_dict[key])

    args.custom_fields = custom_fields

    # Validate required args (after config loading)
    if not args.input:
        parser.error("the following arguments are required: -i/--input (or provide in config)")

    if not args.output and not args.in_place:
        parser.error("Either -o/--output or --in-place must be specified (or provide output in config)")

    if args.output and args.in_place:
        parser.error("Cannot specify both -o/--output and --in-place")

    return args


def main():
    """Main entry point for the tag command."""
    args = get_args()
    set_verbosity(args.verbose)

    # Load input file
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

    # Validate it looks like Olmsted JSON
    if "datasets" not in data or "clones" not in data:
        vprint.error(
            "Error: Input does not appear to be Olmsted JSON "
            "(missing 'datasets' or 'clones' key)"
        )
        sys.exit(1)

    # Custom fields already loaded during arg parsing
    custom_fields = args.custom_fields

    # Ensure metadata has the format tag
    if "metadata" not in data:
        data["metadata"] = {}
    if isinstance(data["metadata"], dict) and "format" not in data["metadata"]:
        data["metadata"]["format"] = "olmsted"

    datasets = data.get("datasets", [])
    clones_dict = data.get("clones", {})
    all_trees = data.get("trees", [])

    retag_datasets_field_metadata(
        datasets, clones_dict, all_trees, custom_fields=custom_fields, mode=args.mode
    )

    for dataset in datasets:
        dataset_id = dataset.get("dataset_id")
        if not dataset_id:
            continue
        levels = list(dataset.get("field_metadata", {}).keys())
        total_fields = sum(len(v) for v in dataset.get("field_metadata", {}).values())
        vprint.verbose(
            f"Dataset '{dataset_id}': {total_fields} fields across levels: {levels}"
        )

    # Write output
    output_path = input_path if args.in_place else Path(args.output)

    indent = 2 if args.json_format == "pretty" else None
    separators = None if args.json_format == "pretty" else (",", ":")

    with open(output_path, "w") as f:
        json.dump(data, f, indent=indent, separators=separators)

    vprint.status(f"Tagged data written to: {output_path}")


if __name__ == "__main__":
    main()
