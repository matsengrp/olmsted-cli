#!/usr/bin/env python
"""
Unified data processing script for Olmsted visualization.

This script can process both AIRR JSON and PCP CSV formats, automatically
detecting the input format or using user-specified format type.

Supported formats:
- AIRR JSON: Standard AIRR format with clones and trees
- PCP CSV: Parent-Child Pair format with optional Newick trees

Output modes:
- Single consolidated JSON file (default): All data in one file
- Multiple files (--split-files): Separate datasets.json, clones.*.json, tree.*.json
"""

import argparse
import csv
import gzip
import json
import os
import sys
import traceback
import uuid
from pathlib import Path

import jsonschema

from .process_airr_data import (
    clone_spec,
    process_dataset,
)
from .process_pcp_data import (
    deterministic_uuid,
    parse_newick_csv,
    parse_pcp_csv,
    process_pcp_to_olmsted,
)
from .process_utils import (
    create_consolidated_data,
    validate_dataset,
    validate_output_data,
    write_out,
)


def detect_file_format(file_path):
    """
    Automatically detect the file format based on file extension and content.

    Args:
        file_path: Path to the input file

    Returns:
        str: Detected format ('airr', 'pcp', or 'unknown')
    """
    file_path = Path(file_path)

    # Check file extension first
    if file_path.suffix.lower() in [".json"]:
        return "airr"
    elif file_path.suffix.lower() in [".csv"]:
        return "pcp"
    elif file_path.suffix.lower() in [".gz"]:
        # Check the extension before .gz
        if file_path.stem.endswith(".json"):
            return "airr"
        elif file_path.stem.endswith(".csv"):
            return "pcp"

    # If extension doesn't help, try to peek at content
    try:
        # Determine if file is gzipped
        if str(file_path).endswith(".gz"):
            file_handle = gzip.open(file_path, "rt")
        else:
            file_handle = open(file_path, "r")

        with file_handle:
            # Read first few lines to detect format
            first_lines = []
            for i, line in enumerate(file_handle):
                first_lines.append(line.strip())
                if i >= 2:  # Read first 3 lines
                    break

            # Check if it looks like JSON
            first_content = "".join(first_lines)
            if first_content.startswith("{") or first_content.startswith("["):
                try:
                    # Try to parse as JSON
                    json.loads(first_content)
                    return "airr"
                except (json.JSONDecodeError, ValueError):
                    pass

            # Check if it looks like CSV with PCP headers
            if first_lines:
                first_line = first_lines[0].lower()
                pcp_indicators = [
                    "sample_id",
                    "parent_name",
                    "child_name",
                    "family_name",
                    "newick",
                ]
                if any(indicator in first_line for indicator in pcp_indicators):
                    return "pcp"

    except Exception as e:
        print(f"Warning: Could not detect format for {file_path}: {e}")

    return "unknown"


def validate_airr_file(file_path):
    """
    Validate that a file contains valid AIRR JSON data.

    Args:
        file_path: Path to the AIRR JSON file

    Returns:
        bool: True if valid AIRR format, False otherwise
    """
    try:
        with open(file_path, "r") as f:
            data = json.load(f)

        # Check for required AIRR fields
        if isinstance(data, dict):
            # Single dataset format
            required_fields = ["dataset_id", "clones"]
            return all(field in data for field in required_fields)
        elif isinstance(data, list):
            # Multiple datasets format
            return all(
                isinstance(item, dict)
                and all(field in item for field in ["dataset_id", "clones"])
                for item in data
            )
    except Exception:
        return False

    return False


def validate_pcp_file(file_path):
    """
    Validate that a file contains valid PCP CSV data.

    Args:
        file_path: Path to the PCP CSV file

    Returns:
        bool: True if valid PCP format, False otherwise
    """
    try:
        # Determine if file is gzipped
        if file_path.endswith(".gz"):
            file_handle = gzip.open(file_path, "rt")
        else:
            file_handle = open(file_path, "r")

        with file_handle:
            reader = csv.DictReader(file_handle)

            # Check for required PCP columns
            required_pcp_cols = {"sample_id", "parent_name", "child_name"}
            required_newick_cols = {"family_name", "newick_tree"}

            if required_pcp_cols.issubset(set(reader.fieldnames)):
                return True
            elif required_newick_cols.issubset(set(reader.fieldnames)):
                return True

    except Exception:
        return False

    return False


def process_airr_format(args):
    """
    Process AIRR format files using the existing AIRR processor.

    Args:
        args: Parsed command line arguments
    """
    print("Processing AIRR format...")

    # Convert unified args to AIRR-specific args
    airr_args = argparse.Namespace()

    # Map common arguments
    airr_args.inputs = args.inputs
    airr_args.output = args.output
    airr_args.data_outdir = getattr(args, "split_files", None)  # For split files mode
    airr_args.verbose = args.verbose
    airr_args.validate = args.validate
    airr_args.strict_validation = args.strict_validation
    airr_args.schema_dir = getattr(args, "schema_dir", None)

    # AIRR-specific arguments with defaults
    airr_args.naive_name = getattr(args, "naive_name", "naive")
    airr_args.remove_invalid_clones = getattr(args, "remove_invalid_clones", False)
    airr_args.display_schema_html = None
    airr_args.display_schema = False
    airr_args.write_schema_yaml = False
    airr_args.root_trees = getattr(args, "root_trees", False)

    # Process using AIRR logic (adapted from process_airr_data.py)
    datasets, clones_dict, trees = [], {}, []

    for infile in airr_args.inputs or []:
        print(f"\nProcessing AIRR file: {infile}")
        try:
            with open(infile, "r") as fh:
                dataset = json.load(fh)
                if airr_args.remove_invalid_clones:
                    dataset["clones"] = list(
                        filter(
                            jsonschema.Draft4Validator(clone_spec).is_valid,
                            dataset["clones"],
                        )
                    )
                # Use unified validation from validate module
                errors = validate_dataset(dataset, verbose=airr_args.verbose)
                if errors:
                    error_msg = "Dataset validation failed"
                    if airr_args.verbose:
                        print(f"Dataset validation failed:")
                        for error in errors:
                            print(f"  - {error}")
                    else:
                        error_msg += ". Please rerun with `-v` for detailed errors"
                    raise Exception(error_msg)
                dataset = process_dataset(airr_args, dataset, clones_dict, trees)
                datasets.append(dataset)

        except Exception:
            print(f"Unable to process AIRR file: {infile}")
            if airr_args.verbose:
                exc_info = sys.exc_info()
                traceback.print_exception(*exc_info)
            else:
                print("Please rerun with `-v` for detailed errors.")
            sys.exit(1)

    # Validate data before writing if requested
    if airr_args.validate and not validate_output_data(
        datasets, clones_dict, trees, airr_args
    ):
        if airr_args.strict_validation:
            print("\nExiting due to validation errors (--strict-validation enabled)")
            sys.exit(1)

    # Write output
    if args.split_files:
        # Multi-file output to specified directory
        output_dir = args.split_files
        os.makedirs(output_dir, exist_ok=True)
        write_out(datasets, output_dir, "datasets.json", airr_args)
        for dataset_id, clones in clones_dict.items():
            write_out(
                clones,
                output_dir + "/",
                "clones." + dataset_id + ".json",
                airr_args,
            )
        for tree in trees:
            write_out(
                tree,
                output_dir + "/",
                "tree." + tree["ident"] + ".json",
                airr_args,
            )
    else:
        # Single consolidated file output (default)
        consolidated_data = create_consolidated_data(
            datasets, clones_dict, trees, args.inputs, "airr", args
        )
        # Ensure output directory exists
        output_dir = os.path.dirname(args.output) or "."
        output_file = os.path.basename(args.output)
        os.makedirs(output_dir, exist_ok=True)
        print(f"Writing consolidated output to {args.output}")
        write_out(consolidated_data, output_dir, output_file, airr_args)


def process_pcp_format(args):
    """
    Process PCP format files using the existing PCP processor.

    Args:
        args: Parsed command line arguments
    """
    print("Processing PCP format...")

    # Set up deterministic UUID generation if seed is provided
    uuid_counter = 0

    def get_uuid():
        nonlocal uuid_counter
        if hasattr(args, "seed") and args.seed is not None:
            uuid_counter += 1
            return deterministic_uuid(args.seed, uuid_counter)
        else:
            return str(uuid.uuid4())

    try:
        # Assume first input is PCP CSV, second (if provided) is Newick trees
        pcp_file = args.inputs[0]
        trees_file = args.inputs[1] if len(args.inputs) > 1 else None

        print(f"Processing PCP CSV: {pcp_file}")
        if hasattr(args, "seed") and args.seed is not None:
            print(f"Using deterministic UUIDs with seed: {args.seed}")

        pcp_families = parse_pcp_csv(pcp_file)
        print(f"Found {len(pcp_families)} families")

        # Parse Newick trees if provided
        newick_trees = None
        if trees_file:
            print(f"Processing Newick trees: {trees_file}")
            newick_trees = parse_newick_csv(trees_file)
            print(f"Found {len(newick_trees)} trees")

        # Convert to Olmsted format
        print("Converting to Olmsted format...")
        datasets, clones_dict, trees = process_pcp_to_olmsted(
            pcp_families, newick_trees, get_uuid
        )

        # Validate data if requested
        if args.validate:
            if not validate_output_data(datasets, clones_dict, trees, args):
                if args.strict_validation:
                    print(
                        "\nExiting due to validation errors (--strict-validation enabled)"
                    )
                    sys.exit(1)

        # Write output
        if args.split_files:
            # Multi-file output to specified directory
            output_dir = args.split_files
            os.makedirs(output_dir, exist_ok=True)
            print(f"Writing output to {output_dir}")
            write_out(datasets, output_dir, "datasets.json", args)
            for dataset_id, clones in clones_dict.items():
                write_out(clones, output_dir, f"clones.{dataset_id}.json", args)
            for tree in trees:
                write_out(tree, output_dir, f"tree.{tree['ident']}.json", args)
        else:
            # Single consolidated file output (default)
            consolidated_data = create_consolidated_data(
                datasets, clones_dict, trees, args.inputs, "pcp", args
            )
            # Ensure output directory exists
            output_dir = os.path.dirname(args.output) or "."
            output_file = os.path.basename(args.output)
            os.makedirs(output_dir, exist_ok=True)
            print(f"Writing consolidated output to {args.output}")
            write_out(consolidated_data, output_dir, output_file, args)

        print("Processing complete!")

    except Exception as e:
        print(f"Error processing PCP format: {e}")
        if args.verbose:
            traceback.print_exc()
        sys.exit(1)


def get_args():
    """Parse command line arguments for the unified processor."""
    parser = argparse.ArgumentParser(
        description="Unified data processor for Olmsted visualization (AIRR JSON and PCP CSV formats)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Auto-detect format and output to single consolidated file (default)
    python process_data.py -i data.json -o output/olmsted_data.json
    python process_data.py -i data.csv -o output/olmsted_data.json

    # Output to multiple files (datasets.json, clones.*.json, tree.*.json)
    python process_data.py -i data.json --split-files output_dir/
    python process_data.py -i data.csv --split-files output_dir/

    # Force specific format
    python process_data.py -i data.json -o output/data.json -f airr
    python process_data.py -i data.csv -o output/data.json -f pcp

    # PCP with separate trees file
    python process_data.py -i data.csv trees.csv -o output/data.json -f pcp

    # With validation
    python process_data.py -i data.json -o output/data.json --validate --strict-validation
        """,
    )

    # Input/Output arguments
    parser.add_argument(
        "-i",
        "--inputs",
        nargs="+",
        required=True,
        help="Input file(s). For AIRR: one or more JSON files. For PCP: CSV file and optional trees CSV file",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Output file path for consolidated JSON (required unless --split-files is used)",
    )

    # Format specification
    parser.add_argument(
        "-f",
        "--format",
        choices=["airr", "pcp", "auto"],
        default="auto",
        help="Input format (default: auto-detect)",
    )

    # Output options
    parser.add_argument(
        "--split-files",
        metavar="DIR",
        help="Output to multiple files in specified directory (datasets.json, clones.*.json, tree.*.json) instead of single consolidated file",
    )

    # Common processing options
    parser.add_argument(
        "--name",
        help="Optional name for the dataset (stored in metadata)",
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose output")
    parser.add_argument(
        "--validate",
        action="store_true",
        help="Validate output data against schemas before writing",
    )
    parser.add_argument(
        "--strict-validation",
        action="store_true",
        help="Exit with error if validation fails (requires --validate)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        help="Random seed for deterministic UUID generation",
    )

    # AIRR-specific options
    parser.add_argument(
        "-n",
        "--naive-name",
        default="naive",
        help="Name of naive/root node for tree rooting (AIRR only)",
    )
    parser.add_argument(
        "-r",
        "--root-trees",
        action="store_true",
        help="Root trees using naive node (AIRR only)",
    )

    return parser.parse_args()


def main():
    """Main entry point for the unified processor."""
    args = get_args()

    # Validate output arguments
    if not args.output and not args.split_files:
        print("Error: Either -o/--output or --split-files must be specified")
        sys.exit(1)

    if args.output and args.split_files:
        print("Error: Cannot specify both -o/--output and --split-files")
        sys.exit(1)

    # Validate inputs
    if not args.inputs:
        print("Error: No input files specified")
        sys.exit(1)

    # Check that input files exist
    for input_file in args.inputs:
        if not os.path.exists(input_file):
            print(f"Error: Input file does not exist: {input_file}")
            sys.exit(1)

    # Determine format
    if args.format == "auto":
        # Auto-detect using first input file
        detected_format = detect_file_format(args.inputs[0])
        if detected_format == "unknown":
            print(f"Error: Could not auto-detect format for {args.inputs[0]}")
            print("Please specify format with -f/--format option")
            sys.exit(1)
        format_to_use = detected_format
        print(f"Auto-detected format: {format_to_use}")
    else:
        format_to_use = args.format
        print(f"Using specified format: {format_to_use}")

    # Validate format matches file content
    if format_to_use == "airr":
        for input_file in args.inputs:
            if not validate_airr_file(input_file):
                print(f"Warning: {input_file} may not be valid AIRR format")
    elif format_to_use == "pcp":
        if not validate_pcp_file(args.inputs[0]):
            print(f"Warning: {args.inputs[0]} may not be valid PCP format")

    # Process based on format
    try:
        if format_to_use == "airr":
            process_airr_format(args)
        elif format_to_use == "pcp":
            process_pcp_format(args)
        else:
            print(f"Error: Unsupported format: {format_to_use}")
            sys.exit(1)

        print(f"\n✓ Successfully processed {format_to_use.upper()} format data")

    except KeyboardInterrupt:
        print("\nProcessing interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"Error during processing: {e}")
        if args.verbose:
            traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
