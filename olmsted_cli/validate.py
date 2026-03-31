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
    VerbosePrinter,
    add_verbosity_args,
    resolve_verbosity,
    validate_clone,
    validate_consolidated_data,
    validate_dataset,
    validate_tree,
)

# Module-level VerbosePrinter, initialized in main()
vprint = VerbosePrinter(1)


def _validate_dataset_with_children(data, validation_errors, check_time_tree=False):
    """Helper function to validate a dataset and its clones/trees."""
    dataset_errors = validate_dataset(data, vprint.level)
    validation_errors.extend(dataset_errors)
    if dataset_errors:
        return

    vprint.verbose("  Dataset schema: PASS")

    clones = data.get("clones", [])
    clone_pass = 0
    clone_fail = 0

    iterator = clones
    if len(clones) > 1:
        iterator = vprint.progress(clones, desc="Validating clones", unit="clone", leave=False)

    for i, clone in enumerate(iterator):
        clone_id = clone.get("clone_id", f"clone-{i}") if isinstance(clone, dict) else f"clone-{i}"
        if hasattr(iterator, "set_description"):
            iterator.set_description(f"Validating clone {clone_id}")

        errors = validate_clone(clone, vprint.level)
        if errors:
            validation_errors.extend([f"Clone {i}: {e}" for e in errors])
            clone_fail += 1
        else:
            clone_pass += 1

        trees = clone.get("trees", []) if isinstance(clone, dict) else []
        for j, tree in enumerate(trees):
            errors = validate_tree(tree, vprint.level, check_time_tree)
            if errors:
                validation_errors.extend([f"Clone {i}, Tree {j}: {e}" for e in errors])

    if clone_pass:
        vprint.verbose(f"  Clones validated: {clone_pass} passed, {clone_fail} failed")


def _validate_items(items, item_type, validate_fn, validation_errors, check_time_tree=False):
    """Validate a list of items (clones or trees) with progress bar."""
    pass_count = 0
    fail_count = 0

    iterator = items
    if len(items) > 1:
        iterator = vprint.progress(items, desc=f"Validating {item_type}s", unit=item_type, leave=False)

    for i, item in enumerate(iterator):
        if isinstance(item, dict):
            item_id = item.get("clone_id", item.get("ident", item.get("tree_id", f"{item_type}-{i}")))
        else:
            item_id = f"{item_type}-{i}"
        if hasattr(iterator, "set_description"):
            iterator.set_description(f"Validating {item_type} {item_id}")

        if item_type == "tree":
            errors = validate_fn(item, vprint.level, check_time_tree)
        else:
            errors = validate_fn(item, vprint.level)

        if errors:
            validation_errors.extend([f"{item_type.title()} {i}: {e}" for e in errors])
            fail_count += 1
        else:
            pass_count += 1

    vprint.verbose(f"  {item_type.title()}s validated: {pass_count} passed, {fail_count} failed")


def _validate_explicit_file_type(data, file_type, filepath, check_time_tree=False):
    """Handle validation for explicitly specified file types."""
    validation_errors = []

    if file_type == "dataset":
        vprint.status(f"Validating as Olmsted dataset: {filepath}")
        _validate_dataset_with_children(data, validation_errors, check_time_tree)

    elif file_type == "clones":
        vprint.status(f"Validating as clone collection: {filepath}")
        if isinstance(data, list):
            _validate_items(data, "clone", validate_clone, validation_errors)
        else:
            validation_errors.append("Expected a list of clones")

    elif file_type == "clone":
        vprint.status(f"Validating as single clone: {filepath}")
        errors = validate_clone(data, vprint.level)
        validation_errors.extend(errors)
        if not errors:
            vprint.verbose("  Clone schema: PASS")

    elif file_type == "trees":
        vprint.status(f"Validating as tree collection: {filepath}")
        if isinstance(data, list):
            _validate_items(data, "tree", validate_tree, validation_errors, check_time_tree)
        else:
            validation_errors.append("Expected a list of trees")

    elif file_type == "tree":
        vprint.status(f"Validating as single tree: {filepath}")
        errors = validate_tree(data, vprint.level, check_time_tree)
        validation_errors.extend(errors)
        if not errors:
            vprint.verbose("  Tree schema: PASS")

    else:
        validation_errors.append(f"Unknown file type: {file_type}")

    return validation_errors


def _auto_detect_array_type(data, filepath, check_time_tree=False):
    """Auto-detect and validate array-type data."""
    validation_errors = []

    if len(data) == 0:
        validation_errors.append("Empty array - unable to determine file type.")
        return validation_errors

    first_item = data[0]
    if not isinstance(first_item, dict):
        validation_errors.append(
            "Unable to determine file type. Use --dataset, --clone(s), or --tree(s) to specify."
        )
        return validation_errors

    if "clone_id" in first_item or "germline_alignment" in first_item:
        vprint.status(f"Auto-detected as clone collection: {filepath}")
        _validate_items(data, "clone", validate_clone, validation_errors)

    elif "dataset_id" in first_item:
        vprint.status(f"Auto-detected as dataset collection: {filepath}")
        _validate_items(data, "dataset", validate_dataset, validation_errors)

    elif "newick" in first_item and "nodes" in first_item:
        vprint.status(f"Auto-detected as tree collection: {filepath}")
        _validate_items(data, "tree", validate_tree, validation_errors, check_time_tree)

    else:
        validation_errors.append(
            "Unable to determine file type. Use --dataset, --clone(s), or --tree(s) to specify."
        )

    return validation_errors


def _auto_detect_object_type(data, filepath, check_time_tree=False):
    """Auto-detect and validate object-type data."""
    validation_errors = []

    if "metadata" in data and "datasets" in data and "clones" in data and "trees" in data:
        vprint.status(f"Auto-detected as consolidated Olmsted format: {filepath}")
        errors = validate_consolidated_data(data, vprint.level, check_time_tree)
        validation_errors.extend(errors)
        if not errors:
            vprint.verbose("  Consolidated format: PASS")

    elif "clones" in data:
        vprint.status(f"Auto-detected as Olmsted dataset: {filepath}")
        _validate_dataset_with_children(data, validation_errors, check_time_tree)

    elif "trees" in data and isinstance(data.get("trees"), list):
        vprint.status(f"Auto-detected as tree collection: {filepath}")
        _validate_items(data["trees"], "tree", validate_tree, validation_errors, check_time_tree)

    elif "newick" in data and "nodes" in data:
        vprint.status(f"Auto-detected as single tree: {filepath}")
        errors = validate_tree(data, vprint.level, check_time_tree)
        validation_errors.extend(errors)
        if not errors:
            vprint.verbose("  Tree schema: PASS")

    elif "clone_id" in data or "germline_alignment" in data:
        vprint.status(f"Auto-detected as single clone: {filepath}")
        errors = validate_clone(data, vprint.level)
        validation_errors.extend(errors)
        if not errors:
            vprint.verbose("  Clone schema: PASS")

    else:
        validation_errors.append(
            "Unable to determine file type. Use --dataset, --clone(s), or --tree(s) to specify."
        )

    return validation_errors


def validate_file(filepath, file_type=None, verbose=1, strict=False, check_time_tree=False):
    """
    Validate a single data file.

    Args:
        filepath: Path to the file to validate
        file_type: Explicit file type or None for auto-detect
        verbose: Verbosity level (0-3). Sets module-level vprint.
        strict: Exit on first validation error
        check_time_tree: Whether to validate time tree constraints

    Returns:
        tuple: (is_valid, list of errors)
    """
    global vprint
    vprint = VerbosePrinter(verbose)
    try:
        with open(filepath, "r") as f:
            data = json.load(f)
    except json.JSONDecodeError as e:
        return False, [f"Failed to parse JSON: {e}"]
    except Exception as e:
        return False, [f"Failed to read file: {e}"]

    if file_type is not None:
        validation_errors = _validate_explicit_file_type(data, file_type, filepath, check_time_tree)
    else:
        if isinstance(data, list):
            validation_errors = _auto_detect_array_type(data, filepath, check_time_tree)
        else:
            validation_errors = _auto_detect_object_type(data, filepath, check_time_tree)

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

  # Validate with verbose output (shows passing steps)
  olmsted validate -v 2 data.json

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

    parser.add_argument(
        "files",
        nargs="*",
        help="JSON files to validate (auto-detect type)",
    )
    parser.add_argument(
        "--dataset",
        "--datasets",
        dest="dataset_files",
        nargs="+",
        metavar="FILE",
        help="Validate files as Olmsted datasets",
    )
    parser.add_argument(
        "--clone",
        dest="clone_files",
        nargs="+",
        metavar="FILE",
        help="Validate files as single clone objects",
    )
    parser.add_argument(
        "--clones",
        dest="clones_files",
        nargs="+",
        metavar="FILE",
        help="Validate files as clone collections (arrays)",
    )
    parser.add_argument(
        "--tree",
        dest="tree_files",
        nargs="+",
        metavar="FILE",
        help="Validate files as single tree objects",
    )
    parser.add_argument(
        "--trees",
        dest="trees_files",
        nargs="+",
        metavar="FILE",
        help="Validate files as tree collections (arrays)",
    )

    add_verbosity_args(parser)

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
    parser.add_argument(
        "--time-tree",
        action="store_true",
        help="Enable time tree validation (check that child distances >= parent distances)",
    )

    return parser.parse_args()


def main():
    """Main entry point for validate command."""
    global vprint
    args = get_args()
    resolve_verbosity(args)
    vprint = VerbosePrinter(args.verbose)

    # Collect all files to validate with their types
    files_to_validate = []

    for filepath in args.files or []:
        files_to_validate.append((filepath, None))
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
        vprint.error("Error: No files specified for validation")
        vprint.error("Use 'olmsted validate --help' for usage information")
        sys.exit(1)

    all_valid = True
    total_errors = 0

    pbar = vprint.progress(files_to_validate, desc="Validating files", unit="file")
    for filepath, file_type in pbar:
        if hasattr(pbar, "set_description"):
            pbar.set_description(f"Validating {Path(filepath).name}")

        vprint.verbose(f"\n{'=' * 60}")
        vprint.verbose(f"Validating: {filepath}")
        if file_type:
            vprint.verbose(f"Type: {file_type} (explicitly specified)")
        else:
            vprint.verbose(f"Type: auto-detect")
        vprint.verbose(f"{'=' * 60}")

        if not Path(filepath).exists():
            vprint.error(f"ERROR: File not found: {filepath}")
            all_valid = False
            total_errors += 1
            if args.strict:
                sys.exit(1)
            continue

        is_valid, errors = validate_file(filepath, file_type, args.verbose, args.strict, args.time_tree)

        if is_valid:
            vprint.status(f"VALID: {filepath}")
        else:
            vprint.status(f"INVALID: {filepath}")
            total_errors += len(errors)

            vprint.verbose("\nErrors found:")
            for error in errors:
                vprint.verbose(f"  - {error}")

            if vprint.level < 2:
                vprint.status(f"  {len(errors)} error(s) found (use -v 2 for details)")

            all_valid = False
            if args.strict:
                sys.exit(1)

    # Summary
    vprint.status(f"\n{'=' * 60}")
    vprint.status("VALIDATION SUMMARY")
    vprint.status(f"{'=' * 60}")
    vprint.status(f"Files validated: {len(files_to_validate)}")
    vprint.status(f"Total errors: {total_errors}")

    if all_valid:
        vprint.status("\nAll files are valid!")
        sys.exit(0)
    else:
        vprint.status("\nValidation failed for one or more files")
        sys.exit(1)


if __name__ == "__main__":
    main()
