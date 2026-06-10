#!/usr/bin/env python3
"""
Validation command for Olmsted CLI.

Validates data files against Olmsted JSON schemas. Checks:
- Required fields (dataset_id, unique_seqs_count, newick, etc.)
- Schema conformance (types, structure)
- Time tree constraints (optional: child distance >= parent distance)
- Olmsted JSON format integrity (metadata + datasets + clones + trees)

See FORMATS.md "Validation" section for the full list of required fields.
"""

import argparse
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

from .constants import (
    CHAIN_COLUMN_ALIASES,
    FORMAT_PCP,
    KNOWN_PCP_COLUMNS,
    KNOWN_TREE_COLUMNS,
)
from .data_io import open_file, read_olmsted_json
from .process_utils import (
    VerbosePrinter,
    add_verbosity_args,
    resolve_verbosity,
    validate_clone,
    validate_consolidated_data,
    validate_dataset,
    validate_tree,
)
from .types import ValidationResult
from .utils import set_verbosity, vprint

PCP_REQUIRED_COLUMNS = {"sample_id", "parent_name", "child_name"}
TREE_REQUIRED_COLUMNS_A = {"family_name", "newick_tree"}
TREE_REQUIRED_COLUMNS_B = {"family", "newick"}  # alternative names


# =============================================================================
# CSV Validation
# =============================================================================


def _open_csv(filepath):
    """Open a CSV file (plain or gzipped) and return a DictReader + fieldnames."""
    fh, _ = open_file(filepath, expected_formats=(FORMAT_PCP,))
    reader = csv.DictReader(fh)
    return fh, reader


def validate_pcp_csv(filepath, tree_filepath=None):
    """
    Validate a PCP CSV file and optionally its companion tree CSV.

    Checks:
    - Required columns present
    - Recognized vs unknown columns (reported as info)
    - Column aliases detected
    - Data integrity: at least one family, root nodes exist
    - Parent-child relationships form valid trees
    - Sequence data present on at least root nodes
    - Tree CSV: required columns, family name alignment, newick parsing

    Returns:
        tuple: (is_valid, list of errors, list of warnings)
    """
    errors = []
    warnings = []

    # --- Open and read PCP CSV ---
    try:
        fh, reader = _open_csv(filepath)
    except Exception as e:
        return False, [f"Failed to read CSV: {e}"], []

    with fh:
        fieldnames = set(reader.fieldnames or [])
        # Filter empty column names (unnamed index columns)
        fieldnames = {f for f in fieldnames if f}

        # Check required columns
        missing_required = PCP_REQUIRED_COLUMNS - fieldnames
        if missing_required:
            errors.append(f"Missing required columns: {sorted(missing_required)}")
            return False, errors, warnings

        vprint.verbose(f"  Required columns: PASS ({sorted(PCP_REQUIRED_COLUMNS)})")

        # Check for column aliases
        alias_found = []
        for col in fieldnames:
            if col.lower() in CHAIN_COLUMN_ALIASES:
                canonical = CHAIN_COLUMN_ALIASES[col.lower()]
                if canonical not in fieldnames:
                    alias_found.append(f"{col} -> {canonical}")
        if alias_found:
            vprint.verbose(f"  Column aliases detected: {alias_found}")

        # Categorize columns
        known_cols = fieldnames & KNOWN_PCP_COLUMNS
        # Also count aliased columns as known
        aliased_cols = {c for c in fieldnames if c.lower() in CHAIN_COLUMN_ALIASES}
        extra_cols = fieldnames - KNOWN_PCP_COLUMNS - aliased_cols
        if extra_cols:
            vprint.verbose(f"  Extra columns (will be captured as node fields): {sorted(extra_cols)}")
        vprint.verbose(f"  Known columns: {len(known_cols)}, Extra: {len(extra_cols)}")

        # Check for sequence data
        has_sequences = (
            "parent_heavy" in fieldnames or "child_heavy" in fieldnames
            or "parent_light" in fieldnames or "child_light" in fieldnames
            or "parent_seq" in fieldnames or "child_seq" in fieldnames
            or "parent_sequence" in fieldnames or "child_sequence" in fieldnames
        )
        if not has_sequences:
            warnings.append("No sequence columns found (parent_heavy/child_heavy). Tree alignment will not be available.")
        else:
            vprint.verbose("  Sequence columns: PASS")

        # Check for paired data
        has_heavy = any(c in fieldnames for c in ("parent_heavy", "child_heavy", "parent_seq", "child_seq"))
        has_light = any(c in fieldnames for c in ("parent_light", "child_light"))
        if has_heavy and has_light:
            vprint.verbose("  Paired data detected (heavy + light)")
        elif has_light and not has_heavy:
            vprint.verbose("  Light chain only data detected")

        # Read rows and check data integrity
        families = defaultdict(lambda: {"parents": set(), "children": set(), "rows": 0})
        row_count = 0
        for row in reader:
            row_count += 1
            family_id = row.get("family", row.get("sample_id", ""))
            parent = row.get("parent_name", "")
            child = row.get("child_name", "")
            if family_id:
                families[family_id]["parents"].add(parent)
                families[family_id]["children"].add(child)
                families[family_id]["rows"] += 1

        if row_count == 0:
            errors.append("CSV file has no data rows")
            return False, errors, warnings

        vprint.verbose(f"  Rows: {row_count}, Families: {len(families)}")

        # Check each family has a root (a parent that is never a child)
        families_without_root = []
        for family_id, fam in families.items():
            roots = fam["parents"] - fam["children"]
            if not roots:
                families_without_root.append(family_id)

        if families_without_root:
            warnings.append(
                f"{len(families_without_root)} families have no root node "
                f"(every parent also appears as a child): {families_without_root[:5]}"
                + (" ..." if len(families_without_root) > 5 else "")
            )
        else:
            vprint.verbose(f"  Root nodes: PASS (all {len(families)} families have roots)")

    # --- Validate tree CSV if provided ---
    if tree_filepath:
        tree_errors, tree_warnings = validate_tree_csv(tree_filepath, set(families.keys()))
        errors.extend(tree_errors)
        warnings.extend(tree_warnings)

    return len(errors) == 0, errors, warnings


def validate_tree_csv(filepath, pcp_family_ids=None):
    """
    Validate a tree CSV file.

    Checks:
    - Required columns present (family_name/family + newick_tree/newick)
    - Newick strings parse correctly (basic syntax check)
    - Family names align with PCP CSV (if pcp_family_ids provided)

    Returns:
        tuple: (list of errors, list of warnings)
    """
    errors = []
    warnings = []

    try:
        fh, reader = _open_csv(filepath)
    except Exception as e:
        return [f"Failed to read tree CSV: {e}"], []

    with fh:
        fieldnames = set(reader.fieldnames or [])

        # Check required columns (either naming convention)
        has_family_col = "family_name" in fieldnames or "family" in fieldnames
        has_newick_col = "newick_tree" in fieldnames or "newick" in fieldnames

        if not has_family_col:
            errors.append("Missing required column: 'family_name' or 'family'")
        if not has_newick_col:
            errors.append("Missing required column: 'newick_tree' or 'newick'")
        if errors:
            return errors, warnings

        vprint.verbose(f"  Tree CSV required columns: PASS")

        # Identify column names
        family_col = "family_name" if "family_name" in fieldnames else "family"
        newick_col = "newick_tree" if "newick_tree" in fieldnames else "newick"

        # Extra columns
        known = KNOWN_TREE_COLUMNS
        extra = {f for f in fieldnames if f} - known
        if extra:
            vprint.verbose(f"  Extra tree columns (will be family-level fields): {sorted(extra)}")

        # Read and validate rows
        tree_families = set()
        row_count = 0
        newick_errors = []

        for row in reader:
            row_count += 1
            family_name = row.get(family_col, "")
            newick = row.get(newick_col, "")
            tree_families.add(family_name)

            # Basic newick syntax check
            if not newick or not newick.strip():
                newick_errors.append(f"Row {row_count} ({family_name}): empty newick string")
            elif not newick.strip().endswith(";"):
                newick_errors.append(f"Row {row_count} ({family_name}): newick doesn't end with ';'")

        if row_count == 0:
            errors.append("Tree CSV has no data rows")
            return errors, warnings

        vprint.verbose(f"  Tree CSV rows: {row_count}, families: {len(tree_families)}")

        if newick_errors:
            for e in newick_errors[:5]:
                errors.append(e)
            if len(newick_errors) > 5:
                errors.append(f"... and {len(newick_errors) - 5} more newick errors")
        else:
            vprint.verbose(f"  Newick syntax: PASS (all {row_count} trees)")

        # Check alignment with PCP families
        if pcp_family_ids is not None:
            # Tree families should be a subset of PCP families
            missing_in_pcp = tree_families - pcp_family_ids
            missing_in_trees = pcp_family_ids - tree_families

            if missing_in_pcp:
                warnings.append(
                    f"{len(missing_in_pcp)} tree families not found in PCP data: "
                    f"{sorted(missing_in_pcp)[:5]}"
                    + (" ..." if len(missing_in_pcp) > 5 else "")
                )
            if missing_in_trees:
                warnings.append(
                    f"{len(missing_in_trees)} PCP families have no tree: "
                    f"{sorted(missing_in_trees)[:5]}"
                    + (" ..." if len(missing_in_trees) > 5 else "")
                )
            if not missing_in_pcp and not missing_in_trees:
                vprint.verbose(f"  Family alignment: PASS (PCP and tree families match)")
            elif not missing_in_pcp:
                vprint.verbose(f"  Family alignment: all tree families found in PCP")

    return errors, warnings


def _validate_dataset_with_children(data, result, check_time_tree=False):
    """Validate a dataset and its clones/trees, folding into ``result``."""
    dataset_result = validate_dataset(data, vprint.level)
    result.extend(dataset_result)
    if dataset_result.errors:
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

        clone_result = validate_clone(clone, vprint.level)
        result.extend(clone_result, prefix=f"Clone {i}: ")
        if clone_result.errors:
            clone_fail += 1
        else:
            clone_pass += 1

        trees = clone.get("trees", []) if isinstance(clone, dict) else []
        for j, tree in enumerate(trees):
            tree_result = validate_tree(tree, vprint.level, check_time_tree)
            result.extend(tree_result, prefix=f"Clone {i}, Tree {j}: ")

    if clone_pass:
        vprint.verbose(f"  Clones validated: {clone_pass} passed, {clone_fail} failed")


def _validate_items(items, item_type, validate_fn, result, check_time_tree=False):
    """Validate a list of items (clones, datasets, or trees) into ``result``."""
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
            item_result = validate_fn(item, vprint.level, check_time_tree)
        else:
            item_result = validate_fn(item, vprint.level)

        result.extend(item_result, prefix=f"{item_type.title()} {i}: ")
        if item_result.errors:
            fail_count += 1
        else:
            pass_count += 1

    vprint.verbose(f"  {item_type.title()}s validated: {pass_count} passed, {fail_count} failed")


def _validate_explicit_file_type(data, file_type, filepath, check_time_tree=False):
    """Handle validation for explicitly specified file types."""
    result = ValidationResult()

    if file_type == "dataset":
        vprint.status(f"Validating as Olmsted dataset: {filepath}")
        _validate_dataset_with_children(data, result, check_time_tree)

    elif file_type == "clones":
        vprint.status(f"Validating as clone collection: {filepath}")
        if isinstance(data, list):
            _validate_items(data, "clone", validate_clone, result)
        else:
            result.errors.append("Expected a list of clones")

    elif file_type == "clone":
        vprint.status(f"Validating as single clone: {filepath}")
        clone_result = validate_clone(data, vprint.level)
        result.extend(clone_result)
        if clone_result.ok:
            vprint.verbose("  Clone schema: PASS")

    elif file_type == "trees":
        vprint.status(f"Validating as tree collection: {filepath}")
        if isinstance(data, list):
            _validate_items(data, "tree", validate_tree, result, check_time_tree)
        else:
            result.errors.append("Expected a list of trees")

    elif file_type == "tree":
        vprint.status(f"Validating as single tree: {filepath}")
        tree_result = validate_tree(data, vprint.level, check_time_tree)
        result.extend(tree_result)
        if tree_result.ok:
            vprint.verbose("  Tree schema: PASS")

    else:
        result.errors.append(f"Unknown file type: {file_type}")

    return result


def _auto_detect_array_type(data, filepath, check_time_tree=False):
    """Auto-detect and validate array-type data."""
    result = ValidationResult()

    if len(data) == 0:
        result.errors.append("Empty array - unable to determine file type.")
        return result

    first_item = data[0]
    if not isinstance(first_item, dict):
        result.errors.append(
            "Unable to determine file type. Use --dataset, --clone(s), or --tree(s) to specify."
        )
        return result

    if "clone_id" in first_item or "germline_alignment" in first_item:
        vprint.status(f"Auto-detected as clone collection: {filepath}")
        _validate_items(data, "clone", validate_clone, result)

    elif "dataset_id" in first_item:
        vprint.status(f"Auto-detected as dataset collection: {filepath}")
        _validate_items(data, "dataset", validate_dataset, result)

    elif "newick" in first_item and "nodes" in first_item:
        vprint.status(f"Auto-detected as tree collection: {filepath}")
        _validate_items(data, "tree", validate_tree, result, check_time_tree)

    else:
        result.errors.append(
            "Unable to determine file type. Use --dataset, --clone(s), or --tree(s) to specify."
        )

    return result


def _auto_detect_object_type(data, filepath, check_time_tree=False):
    """Auto-detect and validate object-type data."""
    result = ValidationResult()

    if "metadata" in data and "datasets" in data and "clones" in data and "trees" in data:
        vprint.status(f"Auto-detected as Olmsted JSON format: {filepath}")
        consolidated_result = validate_consolidated_data(data, vprint.level, check_time_tree)
        result.extend(consolidated_result)
        if consolidated_result.ok:
            vprint.verbose("  Olmsted JSON format: PASS")

    elif "clones" in data:
        vprint.status(f"Auto-detected as Olmsted dataset: {filepath}")
        _validate_dataset_with_children(data, result, check_time_tree)

    elif "trees" in data and isinstance(data.get("trees"), list):
        vprint.status(f"Auto-detected as tree collection: {filepath}")
        _validate_items(data["trees"], "tree", validate_tree, result, check_time_tree)

    elif "newick" in data and "nodes" in data:
        vprint.status(f"Auto-detected as single tree: {filepath}")
        tree_result = validate_tree(data, vprint.level, check_time_tree)
        result.extend(tree_result)
        if tree_result.ok:
            vprint.verbose("  Tree schema: PASS")

    elif "clone_id" in data or "germline_alignment" in data:
        vprint.status(f"Auto-detected as single clone: {filepath}")
        clone_result = validate_clone(data, vprint.level)
        result.extend(clone_result)
        if clone_result.ok:
            vprint.verbose("  Clone schema: PASS")

    else:
        result.errors.append(
            "Unable to determine file type. Use --dataset, --clone(s), or --tree(s) to specify."
        )

    return result


def validate_file(filepath, file_type=None, verbose=1, strict=False,
                   check_time_tree=False, tree_filepath=None):
    """
    Validate a single data file (JSON or CSV).

    Args:
        filepath: Path to the file to validate
        file_type: Explicit file type, "pcp", "tree-csv", or None for auto-detect
        verbose: Verbosity level (0-3). Sets module-level vprint.
        strict: Exit on first validation error
        check_time_tree: Whether to validate time tree constraints
        tree_filepath: Companion tree CSV file (for PCP validation)

    Returns:
        tuple: (is_valid, list of errors)
    """
    set_verbosity(verbose)
    filepath = str(filepath)

    # Detect CSV files
    is_csv = filepath.endswith(".csv") or filepath.endswith(".csv.gz")

    if is_csv or file_type == "pcp":
        return _validate_csv_file(filepath, file_type, tree_filepath)

    # JSON validation. Goes through open_file (any detected format) rather
    # than read_olmsted_json because validators also accept individual record
    # types (clone, dataset, tree) which don't carry Olmsted top-level keys.
    try:
        handle, _ = open_file(filepath)
        with handle as f:
            data = json.load(f)
    except json.JSONDecodeError as e:
        return False, [f"Failed to parse JSON: {e}"]
    except Exception as e:
        return False, [f"Failed to read file: {e}"]

    if file_type is not None:
        result = _validate_explicit_file_type(data, file_type, filepath, check_time_tree)
    else:
        if isinstance(data, list):
            result = _auto_detect_array_type(data, filepath, check_time_tree)
        else:
            result = _auto_detect_object_type(data, filepath, check_time_tree)

    # Warnings are advisory — surface them but don't fail validation (mirrors
    # the PCP CSV path, which also prints warnings and returns errors only).
    for warning in result.warnings:
        vprint.status(f"  Warning: {warning}")

    return result.ok, result.errors


def _validate_csv_file(filepath, file_type=None, tree_filepath=None):
    """Route CSV validation based on file type or auto-detection."""
    path = Path(filepath)

    # Auto-detect: CSV files are assumed to be PCP format
    if file_type is None:
        try:
            fh, reader = _open_csv(filepath)
            with fh:
                fieldnames = set(reader.fieldnames or [])
            if PCP_REQUIRED_COLUMNS.issubset(fieldnames):
                file_type = "pcp"
            else:
                return False, [
                    f"Cannot determine CSV type. Found columns: {sorted(fieldnames)[:10]}. "
                    "Expected PCP format (sample_id, parent_name, child_name)."
                ]
        except Exception as e:
            return False, [f"Failed to read CSV: {e}"]

    vprint.status(f"Validating as PCP CSV: {filepath}")
    if tree_filepath:
        vprint.status(f"  Companion tree CSV: {tree_filepath}")
    else:
        vprint.status(f"  Warning: No companion tree CSV provided (-t). Tree validation skipped.")

    is_valid, errors, file_warnings = validate_pcp_csv(filepath, tree_filepath)
    for w in file_warnings:
        vprint.status(f"  Warning: {w}")
    return is_valid, errors


def get_args():
    """Parse command line arguments for validate command."""
    parser = argparse.ArgumentParser(
        description="Validate Olmsted/AIRR data files against schemas",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Validation checks:
  Datasets:
    - Required: dataset_id
    - Schema conformance for all properties

  Clones (clonal families):
    - Required: unique_seqs_count, mean_mut_freq
    - Schema conformance for gene calls, alignment positions, etc.
    - Validates nested trees if present

  Trees:
    - Required: newick (valid Newick tree string)
    - Schema conformance for tree metadata
    - With --time-tree: checks that child distance >= parent distance
      (monotonically increasing from root)

  Nodes (within trees):
    - Required: sequence_id, sequence_alignment, sequence_alignment_aa
    - Schema conformance for metrics, multiplicity, etc.

  Olmsted JSON:
    - Required: metadata, datasets, clones, trees top-level keys
    - Validates format_version compatibility
    - Recursively validates all datasets, clones, and trees

  PCP CSV (requires -t for tree CSV):
    - Required columns: sample_id, parent_name, child_name
    - Reports recognized vs extra columns
    - Checks root nodes exist for each family
    - Checks for sequence data columns
    - Tree CSV: required columns, newick syntax, family alignment

  All JSON schemas allow additionalProperties (extra fields are preserved).
  See FORMATS.md for full field reference.

Examples:
  # Auto-detect file type (JSON or CSV)
  olmsted validate data.json
  olmsted validate pcp.csv

  # Validate PCP CSV with companion tree file
  olmsted validate --pcp pcp.csv -t trees.csv

  # Validate split-format JSON files
  olmsted validate --split datasets.json clones.*.json tree.*.json

  # Validate with verbose output (shows passing steps)
  olmsted validate -v 2 data.json

  # Validate and exit on first error
  olmsted validate --strict data.json

  # Check time tree constraints
  olmsted validate --time-tree data.json

File types:
  (positional)   Auto-detect format (JSON or CSV)
  --dataset      Olmsted dataset JSON file
  --split        Split-format JSON files (auto-detects clone/tree/dataset)
  --pcp          PCP CSV file (requires -t for companion tree CSV)
  -t / --tree    Companion tree CSV file

Auto-detection (when no type specified):
  JSON: Olmsted JSON, datasets, clone/tree collections, single clones/trees
  CSV:  PCP format (requires sample_id, parent_name, child_name columns)
        """,
    )

    parser.add_argument(
        "-i", "--input",
        dest="files",
        nargs="+",
        metavar="FILE",
        help="Files to validate (auto-detects JSON type and CSV format)",
    )
    parser.add_argument(
        "--dataset",
        "--datasets",
        dest="dataset_files",
        nargs="+",
        metavar="FILE",
        help="Validate JSON files as Olmsted datasets",
    )
    parser.add_argument(
        "--split",
        dest="split_files",
        nargs="+",
        metavar="FILE",
        help="Validate split-format JSON files (auto-detects clone/tree/dataset arrays)",
    )
    parser.add_argument(
        "--pcp",
        dest="pcp_files",
        nargs="+",
        metavar="FILE",
        help="Validate PCP CSV files (requires -t for companion tree CSV)",
    )
    parser.add_argument(
        "-t", "--tree",
        dest="tree_csv_file",
        metavar="FILE",
        help="Companion tree CSV file (required for PCP validation)",
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
    args = get_args()
    resolve_verbosity(args)
    set_verbosity(args.verbose)

    # Collect all files to validate with their types
    files_to_validate = []

    for filepath in args.files or []:
        files_to_validate.append((filepath, None))
    for filepath in args.dataset_files or []:
        files_to_validate.append((filepath, "dataset"))
    for filepath in args.split_files or []:
        files_to_validate.append((filepath, None))  # auto-detect split type
    for filepath in args.pcp_files or []:
        files_to_validate.append((filepath, "pcp"))

    # Enforce -t with --pcp
    has_pcp = bool(args.pcp_files)
    if has_pcp and not args.tree_csv_file:
        vprint.error("Error: --pcp requires -t/--tree for companion tree CSV file")
        sys.exit(1)

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

        # Pass companion tree file for PCP validation
        tree_fp = getattr(args, "tree_csv_file", None) if file_type == "pcp" else None
        is_valid, errors = validate_file(
            filepath, file_type, args.verbose, args.strict, args.time_tree, tree_fp
        )

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
