#!/usr/bin/env python3
"""
Validation command for Olmsted CLI.

Validates AIRR JSON or Olmsted dataset files against their schemas.
"""

import argparse
import json
import sys
from pathlib import Path

from .process_utils import (
    validate_clone,
    validate_consolidated_data,
    validate_dataset,
    validate_tree,
)


def _validate_dataset_with_children(data, verbose, validation_errors):
    """Helper function to validate a dataset and its clones/trees."""
    validation_errors.extend(validate_dataset(data, verbose))

    # If dataset is valid, validate its clones and trees
    if not validation_errors:
        for i, clone in enumerate(data.get("clones", [])):
            errors = validate_clone(clone, verbose)
            if errors:
                validation_errors.extend([f"Clone {i}: {e}" for e in errors])

            for j, tree in enumerate(clone.get("trees", [])):
                errors = validate_tree(tree, verbose)
                if errors:
                    validation_errors.extend(
                        [f"Clone {i}, Tree {j}: {e}" for e in errors]
                    )


def _validate_explicit_file_type(data, file_type, filepath, verbose):
    """Handle validation for explicitly specified file types."""
    validation_errors = []

    if file_type == "dataset":
        print(f"Validating as Olmsted dataset: {filepath}")
        _validate_dataset_with_children(data, verbose, validation_errors)

    elif file_type == "clones":
        print(f"Validating as clone collection: {filepath}")
        if isinstance(data, list):
            for i, clone in enumerate(data):
                errors = validate_clone(clone, verbose)
                if errors:
                    validation_errors.extend([f"Clone {i}: {e}" for e in errors])
        else:
            validation_errors.append("Expected a list of clones")

    elif file_type == "clone":
        print(f"Validating as single clone: {filepath}")
        errors = validate_clone(data, verbose)
        validation_errors.extend(errors)

    elif file_type == "trees":
        print(f"Validating as tree collection: {filepath}")
        if isinstance(data, list):
            for i, tree in enumerate(data):
                errors = validate_tree(tree, verbose)
                if errors:
                    validation_errors.extend([f"Tree {i}: {e}" for e in errors])
        else:
            validation_errors.append("Expected a list of trees")

    elif file_type == "tree":
        print(f"Validating as single tree: {filepath}")
        errors = validate_tree(data, verbose)
        validation_errors.extend(errors)

    else:
        validation_errors.append(f"Unknown file type: {file_type}")

    return validation_errors


def _auto_detect_array_type(data, filepath, verbose):
    """Auto-detect and validate array-type data."""
    validation_errors = []

    if len(data) == 0:
        validation_errors.append(
            "Empty array - unable to determine file type for validation."
        )
        return validation_errors

    first_item = data[0]
    if not isinstance(first_item, dict):
        validation_errors.append(
            "Unable to determine file type for validation. Use --dataset, --clone(s), or --tree(s) to specify."
        )
        return validation_errors

    if "clone_id" in first_item or "germline_alignment" in first_item:
        # Array of clones
        print(f"Auto-detected as clone collection: {filepath}")
        for i, clone in enumerate(data):
            errors = validate_clone(clone, verbose)
            if errors:
                validation_errors.extend([f"Clone {i}: {e}" for e in errors])

    elif "dataset_id" in first_item:
        # Array of datasets
        print(f"Auto-detected as dataset collection: {filepath}")
        for i, dataset in enumerate(data):
            errors = validate_dataset(dataset, verbose)
            if errors:
                validation_errors.extend([f"Dataset {i}: {e}" for e in errors])

    elif "newick" in first_item and "nodes" in first_item:
        # Array of trees
        print(f"Auto-detected as tree collection: {filepath}")
        for i, tree in enumerate(data):
            errors = validate_tree(tree, verbose)
            if errors:
                validation_errors.extend([f"Tree {i}: {e}" for e in errors])

    else:
        validation_errors.append(
            "Unable to determine file type for validation. Use --dataset, --clone(s), or --tree(s) to specify."
        )

    return validation_errors


def _auto_detect_object_type(data, filepath, verbose):
    """Auto-detect and validate object-type data."""
    validation_errors = []

    if "metadata" in data and "datasets" in data and "clones" in data and "trees" in data:
        # This looks like consolidated format
        print(f"Auto-detected as consolidated Olmsted format: {filepath}")
        errors = validate_consolidated_data(data, verbose)
        validation_errors.extend(errors)

    elif "clones" in data:
        # This looks like a dataset file
        print(f"Auto-detected as Olmsted dataset: {filepath}")
        _validate_dataset_with_children(data, verbose, validation_errors)

    elif "trees" in data and isinstance(data.get("trees"), list):
        # Validate as a collection of trees
        print(f"Auto-detected as tree collection: {filepath}")
        for i, tree in enumerate(data["trees"]):
            errors = validate_tree(tree, verbose)
            if errors:
                validation_errors.extend([f"Tree {i}: {e}" for e in errors])

    elif "newick" in data and "nodes" in data:
        # Single tree
        print(f"Auto-detected as single tree: {filepath}")
        errors = validate_tree(data, verbose)
        validation_errors.extend(errors)

    elif "clone_id" in data or "germline_alignment" in data:
        # Single clone
        print(f"Auto-detected as single clone: {filepath}")
        errors = validate_clone(data, verbose)
        validation_errors.extend(errors)

    else:
        validation_errors.append(
            "Unable to determine file type for validation. Use --dataset, --clone(s), or --tree(s) to specify."
        )

    return validation_errors


def validate_file(filepath, file_type=None, verbose=False, strict=False):
    """
    Validate a single data file.

    Args:
        filepath: Path to the file to validate
        file_type: Explicit file type ('dataset', 'clone', 'tree', 'clones', 'trees', or None for auto-detect)
        verbose: Show detailed validation errors
        strict: Exit on first validation error

    Returns:
        tuple: (is_valid, list of errors)
    """
    # Load and parse the file
    try:
        with open(filepath, "r") as f:
            data = json.load(f)
    except json.JSONDecodeError as e:
        return False, [f"Failed to parse JSON: {e}"]
    except Exception as e:
        return False, [f"Failed to read file: {e}"]

    # Validate based on file type
    if file_type is not None:
        # Use explicitly specified file type
        validation_errors = _validate_explicit_file_type(data, file_type, filepath, verbose)
    else:
        # Auto-detect file type based on content
        if isinstance(data, list):
            validation_errors = _auto_detect_array_type(data, filepath, verbose)
        else:
            validation_errors = _auto_detect_object_type(data, filepath, verbose)

    return len(validation_errors) == 0, validation_errors


def get_args():
    """Parse command line arguments for validate command."""
    parser = argparse.ArgumentParser(
        description="Validate Olmsted/AIRR data files against schemas",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Auto-detect file type
  olmsted validate data.json

  # Explicitly specify file type
  olmsted validate --dataset datasets.json
  olmsted validate --clones clones.family1.json clones.family2.json
  olmsted validate --tree tree.abc123.json
  olmsted validate --trees trees.collection.json

  # Validate multiple files of different types
  olmsted validate --dataset dataset.json --clones clones.*.json --tree tree.*.json

  # Validate with verbose output
  olmsted validate -v data.json

  # Validate and exit on first error
  olmsted validate --strict data.json

File types:
  --dataset: Olmsted dataset file containing clones
  --clone:   Single clone object
  --clones:  Array/collection of clone objects
  --tree:    Single tree object with newick and nodes
  --trees:   Array/collection of tree objects

Auto-detection (when no type specified):
  - Olmsted datasets (contain 'clones' field)
  - Clone collections (contain 'clone_id' or 'germline_alignment')
  - Tree collections (contain 'trees' array)
  - Single trees (contain 'newick' and 'nodes')
        """,
    )

    # File arguments with explicit types
    parser.add_argument(
        "files", nargs="*", help="JSON files to validate (auto-detect type)"
    )

    parser.add_argument(
        "--dataset",
        "--datasets",
        nargs="+",
        dest="dataset_files",
        metavar="FILE",
        help="Validate files as Olmsted datasets",
    )

    parser.add_argument(
        "--clone",
        nargs="+",
        dest="clone_files",
        metavar="FILE",
        help="Validate files as single clone objects",
    )

    parser.add_argument(
        "--clones",
        nargs="+",
        dest="clones_files",
        metavar="FILE",
        help="Validate files as clone collections (arrays)",
    )

    parser.add_argument(
        "--tree",
        nargs="+",
        dest="tree_files",
        metavar="FILE",
        help="Validate files as single tree objects",
    )

    parser.add_argument(
        "--trees",
        nargs="+",
        dest="trees_files",
        metavar="FILE",
        help="Validate files as tree collections (arrays)",
    )

    # Validation options
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Show detailed validation errors"
    )

    parser.add_argument(
        "--strict",
        action="store_true",
        help="Exit with error code on first validation failure",
    )

    parser.add_argument(
        "--schema",
        choices=["airr", "olmsted", "both"],
        default="both",
        help="Which schema to validate against (default: both)",
    )

    return parser.parse_args()


def main():
    """Main entry point for validate command."""
    args = get_args()

    # Collect all files to validate with their types
    files_to_validate = []

    # Add files with auto-detect
    for filepath in args.files or []:
        files_to_validate.append((filepath, None))

    # Add explicitly typed files
    for filepath in args.dataset_files or []:
        files_to_validate.append((filepath, "dataset"))

    for filepath in args.clone_files or []:
        files_to_validate.append((filepath, "clone"))

    for filepath in args.clones_files or []:
        files_to_validate.append((filepath, "clones"))

    for filepath in args.tree_files or []:
        files_to_validate.append((filepath, "tree"))

    for filepath in args.trees_files or []:
        files_to_validate.append((filepath, "trees"))

    if not files_to_validate:
        print("Error: No files specified for validation")
        print("Use 'olmsted validate --help' for usage information")
        sys.exit(1)

    all_valid = True
    total_errors = 0

    for filepath, file_type in files_to_validate:
        print(f"\n{'=' * 60}")
        print(f"Validating: {filepath}")
        if file_type:
            print(f"Type: {file_type} (explicitly specified)")
        else:
            print(f"Type: auto-detect")
        print(f"{'=' * 60}")

        if not Path(filepath).exists():
            print(f"ERROR: File not found: {filepath}")
            all_valid = False
            total_errors += 1
            if args.strict:
                sys.exit(1)
            continue

        is_valid, errors = validate_file(filepath, file_type, args.verbose, args.strict)

        if is_valid:
            print(f"✅ VALID - {filepath}")
        else:
            print(f"❌ INVALID - {filepath}")
            total_errors += len(errors)

            if args.verbose:
                print("\nErrors found:")
                for error in errors:
                    print(f"  - {error}")
            else:
                print(f"  {len(errors)} error(s) found (use -v for details)")

            all_valid = False
            if args.strict:
                sys.exit(1)

    # Summary
    print(f"\n{'=' * 60}")
    print("VALIDATION SUMMARY")
    print(f"{'=' * 60}")
    print(f"Files validated: {len(files_to_validate)}")
    print(f"Total errors: {total_errors}")

    if all_valid:
        print("\n✅ All files are valid!")
        sys.exit(0)
    else:
        print("\n❌ Validation failed for one or more files")
        sys.exit(1)


if __name__ == "__main__":
    main()
