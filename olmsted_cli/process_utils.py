#!/usr/bin/env python
"""
Shared utilities for processing various data formats in Olmsted.

This module contains common functions and constants used by
process_airr_data.py, process_pcp_data.py, and other data processors.
"""

import csv
import json
import os
import uuid
from datetime import datetime, timezone

import jsonschema
import yaml

from .schemas import SCHEMA_VERSION, clone_spec, dataset_spec, tree_spec

# Constants for infinity handling
inf = float("inf")
neginf = float("-inf")


# General utility functions
def comp(f, g):
    """
    Function composition: comp(f, g)(x) == f(g(x))
    """

    def h(*args, **kw_args):
        return f(g(*args, **kw_args))

    return h


def strip_ns(a):
    # Handle namespace stripping for both : and / separators
    return str(a).split(":")[-1].split("/")[-1]


def dict_subset(d, keys):
    return {k: d[k] for k in keys if k in d}


def merge(d, d2):
    """
    Merge d2 into d, returning a new dict (non-mutating).
    """
    d = d.copy()
    d.update(d2)
    return d


def get_in(d, path):
    """
    Retrieve value from nested dictionary using a path list.

    Args:
        d: Dictionary to traverse
        path: List of keys representing path to value

    Returns:
        Value at path or empty dict if path doesn't exist
    """
    return (
        d
        if len(path) == 0
        else get_in(d.get(path[0]) if isinstance(d, dict) else {}, path[1:])
    )


def clean_record(d):
    """
    Clean a record by removing namespaces and handling special values.

    Args:
        d: Data to clean (dict, list, or value)

    Returns:
        Cleaned data
    """
    if isinstance(d, list):
        return list(map(clean_record, d))
    elif isinstance(d, dict):
        return {strip_ns(k): clean_record(v) for k, v in d.items()}
    # can't have infinity in json
    elif d == inf or d == neginf:
        return None
    else:
        return d


def spy(x):
    print("debugging:", x)
    return x


def lspy(xs):
    xs_ = list(xs)
    print("debugging listable:", xs_)
    return xs_


def nospy(xs):
    return xs


def translate_dna_to_aa(dna_sequence):
    """
    Translate DNA sequence to amino acid sequence.
    Uses standard genetic code, handles ambiguous bases.
    """
    if not dna_sequence:
        return ""

    # Standard genetic code
    codon_table = {
        "TTT": "F",
        "TTC": "F",
        "TTA": "L",
        "TTG": "L",
        "TCT": "S",
        "TCC": "S",
        "TCA": "S",
        "TCG": "S",
        "TAT": "Y",
        "TAC": "Y",
        "TAA": "*",
        "TAG": "*",
        "TGT": "C",
        "TGC": "C",
        "TGA": "*",
        "TGG": "W",
        "CTT": "L",
        "CTC": "L",
        "CTA": "L",
        "CTG": "L",
        "CCT": "P",
        "CCC": "P",
        "CCA": "P",
        "CCG": "P",
        "CAT": "H",
        "CAC": "H",
        "CAA": "Q",
        "CAG": "Q",
        "CGT": "R",
        "CGC": "R",
        "CGA": "R",
        "CGG": "R",
        "ATT": "I",
        "ATC": "I",
        "ATA": "I",
        "ATG": "M",
        "ACT": "T",
        "ACC": "T",
        "ACA": "T",
        "ACG": "T",
        "AAT": "N",
        "AAC": "N",
        "AAA": "K",
        "AAG": "K",
        "AGT": "S",
        "AGC": "S",
        "AGA": "R",
        "AGG": "R",
        "GTT": "V",
        "GTC": "V",
        "GTA": "V",
        "GTG": "V",
        "GCT": "A",
        "GCC": "A",
        "GCA": "A",
        "GCG": "A",
        "GAT": "D",
        "GAC": "D",
        "GAA": "E",
        "GAG": "E",
        "GGT": "G",
        "GGC": "G",
        "GGA": "G",
        "GGG": "G",
    }

    aa_sequence = ""
    # Process in chunks of 3 nucleotides
    for i in range(0, len(dna_sequence) - 2, 3):
        codon = dna_sequence[i : i + 3].upper()
        # Handle ambiguous bases by using 'X' for unknown amino acids
        if len(codon) == 3 and codon in codon_table:
            aa_sequence += codon_table[codon]
        else:
            aa_sequence += "X"  # Unknown amino acid for ambiguous codons

    return aa_sequence


# Additional utility functions consolidated from process_cft_data.py


def rename_keys(record, mapping, to_keep=None):
    """
    Rename keys in a record based on a mapping dictionary.

    Args:
        record: Dictionary to modify
        mapping: Dict mapping old keys to new keys
        to_keep: List of keys to keep with original name (copy, don't move)
    """
    if to_keep is None:
        to_keep = []

    for k in mapping.keys():
        if k in record:
            record[mapping[k]] = record.pop(k) if k not in to_keep else record[k]


def remap_list(lst, mapping):
    """Apply key renaming to all elements in a list."""
    for element in lst:
        rename_keys(element, mapping)


def remap_dict_values(d, mapping):
    """Apply key renaming to all values in a dictionary."""
    for v in d.values():
        rename_keys(v, mapping)


def try_del(d, attr):
    """Safely delete an attribute from a dictionary, ignoring errors."""
    try:
        del d[attr]
    except (KeyError, TypeError):
        pass


def listof(xs_str, f=None):
    """Split a colon-separated string and apply optional function to each element."""
    if f is None:
        f = lambda x: x
    return list(map(f, xs_str.split(":")))


def listofint(xs_str):
    """Split a colon-separated string and convert each element to int."""
    return listof(xs_str, int)


# JSON utility functions
def json_rep(x):
    """
    JSON serialization helper for non-standard types.

    Converts UUID objects to strings and other iterables to lists.
    Used as the 'default' parameter for json.dump().

    Args:
        x: Object to convert

    Returns:
        JSON-serializable representation
    """
    if isinstance(x, uuid.UUID):
        return str(x)
    else:
        # Try to convert to list (for sets, tuples, etc.)
        try:
            return list(x)
        except TypeError:
            # Let json.dump() handle the error for truly non-serializable types
            raise


def write_out(data, dirname, filename, args):
    """
    Write data to JSON or CSV file with proper formatting and UUID handling.

    Args:
        data: Data to write
        dirname: Directory path
        filename: File name
        args: Command line arguments (for verbose flag and csv flag)
    """
    # Ensure directory exists
    os.makedirs(dirname, exist_ok=True)

    # Normalize path
    full_path = os.path.normpath(os.path.join(dirname, filename))

    # Print status
    print(f"writing {full_path}")

    with open(full_path, "w") as fh:
        # Check if CSV output is requested (for CFT data)
        if hasattr(args, "csv") and args.csv and isinstance(data, list):
            # Write as CSV
            if data:
                # Ensure all items are dictionaries
                data = [{k: v for k, v in d.items()} for d in data]
                writer = csv.DictWriter(fh, fieldnames=sorted(data[0].keys()))
                writer.writeheader()
                writer.writerows(data)
        elif isinstance(data, (list, dict)):
            # Write as JSON
            json.dump(
                data,
                fh,
                default=json_rep,  # Handle UUIDs and other non-serializable types
                indent=4,
                # allow_nan=False  # Uncomment if you want to disallow NaN values
            )
        else:
            # Handle raw string data
            fh.write(data)


# Version Constants
CONSOLIDATED_JSON_VERSION = "1.0"


def create_consolidated_data(
    datasets, clones_dict, trees, input_files, detected_format, args=None
):
    """
    Create consolidated data structure with metadata.

    Args:
        datasets: List of dataset objects
        clones_dict: Dictionary of clone lists by dataset_id
        trees: List of tree objects
        input_files: List of input file paths
        detected_format: Detected or specified format ('airr' or 'pcp')
        args: Command line arguments (optional)

    Returns:
        dict: Consolidated data with metadata
    """
    # Generate metadata
    metadata = {
        "format_version": CONSOLIDATED_JSON_VERSION,
        "schema_version": SCHEMA_VERSION,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "source_format": detected_format,
        "source_files": [os.path.basename(f) for f in input_files],
        "processing_info": {
            "datasets_count": len(datasets),
            "total_clones_count": sum(len(clones) for clones in clones_dict.values()),
            "total_trees_count": len(trees),
        },
        "generated_by": {
            "tool": "olmsted-cli",
            "version": SCHEMA_VERSION,
        },
    }
    
    # Add optional name if provided
    if args and hasattr(args, "name") and args.name:
        metadata["name"] = args.name

    # Add processing options if available
    if args:
        metadata["processing_options"] = {
            "validation": getattr(args, "validate", False),
            "strict_validation": getattr(args, "strict_validation", False),
            "seed": getattr(args, "seed", None),
        }

        # Add format-specific options
        if detected_format == "airr":
            metadata["processing_options"]["airr"] = {
                "naive_name": getattr(args, "naive_name", "naive"),
                "root_trees": getattr(args, "root_trees", False),
            }

    return {
        "metadata": metadata,
        "datasets": datasets,
        "clones": clones_dict,
        "trees": trees,
    }


# Schema utility functions
def natural_number(desc):
    """Create a natural number schema specification with description."""
    return {"description": desc, "minimum": 0, "type": "integer"}


def is_nullable_string(checker, instance):
    """Check if an instance is either a string or null (for JSON schema validation)."""
    return jsonschema.Draft4Validator.TYPE_CHECKER.is_type(
        instance, "string"
    ) or jsonschema.Draft4Validator.TYPE_CHECKER.is_type(instance, "null")


# Schema loading and validation functions
def load_schema(schema_path):
    """Load a JSON schema from file (supports both JSON and YAML)."""
    with open(schema_path, "r") as f:
        if schema_path.endswith(".yaml") or schema_path.endswith(".yml"):
            return yaml.safe_load(f)
        else:
            return json.load(f)


def load_official_airr_schema():
    """
    Load the official AIRR schema from airr-standards/specs/airr-schema.yaml.

    Returns:
        dict: The full AIRR schema dictionary, or None if not found
    """
    try:
        # Load from airr-standards directory (development)
        script_dir = os.path.dirname(os.path.abspath(__file__))
        schema_path = os.path.join(
            script_dir, "..", "airr-standards", "specs", "airr-schema.yaml"
        )

        with open(schema_path, "r") as f:
            return yaml.safe_load(f)
    except FileNotFoundError:
        print("Warning: Official AIRR schema not found")
        return None
    except Exception as e:
        print(f"Warning: Failed to load official AIRR schema: {e}")
        return None


def validate_against_airr_schema(data, schema_object_name, schema=None):
    """
    Validate data against a specific object in the official AIRR schema.

    Args:
        data: The data to validate
        schema_object_name: Name of the schema object (e.g., 'Clone', 'Tree', 'Node')
        schema: Optional pre-loaded schema dict. If None, loads from official source.

    Returns:
        tuple: (is_valid, error_message)
    """
    try:
        if schema is None:
            schema = load_official_airr_schema()

        if schema is None:
            return False, "Official AIRR schema not available"

        if schema_object_name not in schema:
            return (
                False,
                f"Schema object '{schema_object_name}' not found in AIRR schema",
            )

        object_schema = schema[schema_object_name]

        # Create a standalone JSON schema for validation with proper null handling
        properties = {}
        for prop_name, prop_schema in object_schema.get("properties", {}).items():
            # Handle nullable fields by allowing both the original type and null
            if isinstance(prop_schema, dict) and "type" in prop_schema:
                if (
                    prop_schema.get("Description", "").lower().find("null") != -1
                    or prop_schema.get("description", "").lower().find("null") != -1
                ):
                    # Field explicitly mentions null in description, make it nullable
                    new_prop = prop_schema.copy()
                    if isinstance(new_prop["type"], str):
                        new_prop["type"] = [new_prop["type"], "null"]
                    properties[prop_name] = new_prop
                else:
                    properties[prop_name] = prop_schema
            else:
                properties[prop_name] = prop_schema

        validation_schema = {
            "$schema": "http://json-schema.org/draft-07/schema#",
            "type": object_schema.get("type", "object"),
            "required": object_schema.get("required", []),
            "properties": properties,
            "additionalProperties": object_schema.get("additionalProperties", True),
        }

        jsonschema.validate(instance=data, schema=validation_schema)
        return True, None

    except jsonschema.ValidationError as e:
        return False, str(e)
    except Exception as e:
        return False, f"Validation error: {str(e)}"


def validate_airr_clone(clone_data, schema=None):
    """
    Validate clone data against official AIRR Clone schema.

    Args:
        clone_data: The clone data to validate
        schema: Optional pre-loaded schema dict

    Returns:
        tuple: (is_valid, error_message)
    """
    return validate_against_airr_schema(clone_data, "Clone", schema)


def validate_airr_tree(tree_data, schema=None):
    """
    Validate tree data against official AIRR Tree schema.

    Args:
        tree_data: The tree data to validate
        schema: Optional pre-loaded schema dict

    Returns:
        tuple: (is_valid, error_message)
    """
    return validate_against_airr_schema(tree_data, "Tree", schema)


def validate_airr_node(node_data, schema=None):
    """
    Validate node data against official AIRR Node schema.

    Args:
        node_data: The node data to validate
        schema: Optional pre-loaded schema dict

    Returns:
        tuple: (is_valid, error_message)
    """
    return validate_against_airr_schema(node_data, "Node", schema)


def validate_output_data(datasets, clones_dict, trees, args):
    """
    Generic output data validation function using unified validation from validate module.

    This replaces the processor-specific validation functions to provide a single
    validation entry point for all output data.

    Args:
        datasets: List of dataset objects
        clones_dict: Dictionary of clone lists by dataset_id
        trees: List of tree objects
        args: Command line arguments with validate and verbose flags

    Returns:
        bool: True if all validation passes, False otherwise
    """
    if not hasattr(args, "validate") or not args.validate:
        return True

    print("\nValidating output data against schemas...")

    validation_passed = True
    total_errors = 0

    try:
        # Validate datasets
        for i, dataset in enumerate(datasets):
            errors = validate_dataset(dataset, verbose=getattr(args, "verbose", False))
            if errors:
                print(f"❌ Dataset {i} validation failed:")
                for error in errors:
                    print(f"  - {error}")
                validation_passed = False
                total_errors += len(errors)
            elif getattr(args, "verbose", False):
                print(f"✓ Dataset {i} validation passed")

        # Validate clones
        clone_count = 0
        clone_failures = 0
        for dataset_id, clones in clones_dict.items():
            for clone in clones:
                clone_count += 1
                errors = validate_clone(clone, verbose=getattr(args, "verbose", False))
                if errors:
                    clone_failures += 1
                    if getattr(args, "verbose", False):
                        print(
                            f"❌ Clone {clone.get('clone_id', 'unknown')} validation failed:"
                        )
                        for error in errors:
                            print(f"  - {error}")
                    validation_passed = False
                    total_errors += len(errors)

        if clone_failures == 0:
            print(f"✓ Clone validation passed ({clone_count} clones)")
        else:
            print(f"❌ Clone validation: {clone_failures}/{clone_count} failed")

        # Validate trees
        tree_count = 0
        tree_failures = 0
        for tree in trees:
            tree_count += 1
            errors = validate_tree(tree, verbose=getattr(args, "verbose", False))
            if errors:
                tree_failures += 1
                if getattr(args, "verbose", False):
                    print(f"❌ Tree {tree.get('ident', 'unknown')} validation failed:")
                    for error in errors:
                        print(f"  - {error}")
                validation_passed = False
                total_errors += len(errors)

        if tree_failures == 0:
            print(f"✓ Tree validation passed ({tree_count} trees)")
        else:
            print(f"❌ Tree validation: {tree_failures}/{tree_count} failed")

        if total_errors > 0:
            print(f"\nTotal validation errors: {total_errors}")

    except Exception as e:
        print(f"Validation error: {str(e)}")
        validation_passed = False

    return validation_passed


def validate_dataset(data, verbose=False):
    """
    Validate a dataset against the Olmsted dataset schema.

    Args:
        data: Dataset dictionary
        verbose: Show detailed errors

    Returns:
        list: List of validation errors (empty if valid)
    """
    errors = []

    try:
        # Create validator
        validator = jsonschema.Draft4Validator(dataset_spec)

        if not validator.is_valid(data):
            if verbose:
                for error in validator.iter_errors(data):
                    error_path = (
                        " -> ".join(str(p) for p in error.path)
                        if error.path
                        else "root"
                    )
                    errors.append(
                        f"Dataset schema error at {error_path}: {error.message}"
                    )
            else:
                errors.append("Dataset does not conform to schema (use -v for details)")
    except Exception as e:
        errors.append(f"Dataset validation error: {e}")

    return errors


def validate_clone(data, verbose=False):
    """
    Validate a clone against AIRR and Olmsted schemas.

    Args:
        data: Clone dictionary
        verbose: Show detailed errors

    Returns:
        list: List of validation errors (empty if valid)
    """
    errors = []

    # Try AIRR validation first (silently)
    is_airr_valid, airr_error = validate_airr_clone(data)

    # Try Olmsted schema validation
    olmsted_errors = []
    try:
        validator = jsonschema.Draft4Validator(clone_spec)
        if not validator.is_valid(data):
            for error in validator.iter_errors(data):
                error_path = (
                    " -> ".join(str(p) for p in error.path) if error.path else "root"
                )
                olmsted_errors.append(
                    f"Clone schema error at {error_path}: {error.message}"
                )
    except Exception as e:
        olmsted_errors.append(f"Clone validation error: {e}")

    # If both validations fail, report errors
    if not is_airr_valid and olmsted_errors:
        if verbose:
            errors.append(f"AIRR validation: {airr_error}")
            errors.extend(olmsted_errors)
        else:
            errors.append(
                "Clone does not conform to AIRR or Olmsted schema (use -v for details)"
            )
    elif olmsted_errors:
        # Only Olmsted validation failed
        errors.extend(olmsted_errors)
    # If AIRR validation passed OR Olmsted validation passed, consider it valid (no errors)

    return errors


def validate_tree(data, verbose=False):
    """
    Validate a tree against AIRR and Olmsted schemas.

    Args:
        data: Tree dictionary
        verbose: Show detailed errors

    Returns:
        list: List of validation errors (empty if valid)
    """
    errors = []

    # Try AIRR validation first (silently)
    is_airr_valid, airr_error = validate_airr_tree(data)

    # Try Olmsted schema validation
    olmsted_errors = []
    try:
        validator = jsonschema.Draft4Validator(tree_spec)
        if not validator.is_valid(data):
            for error in validator.iter_errors(data):
                error_path = (
                    " -> ".join(str(p) for p in error.path) if error.path else "root"
                )
                olmsted_errors.append(
                    f"Tree schema error at {error_path}: {error.message}"
                )
    except Exception as e:
        olmsted_errors.append(f"Tree validation error: {e}")

    # If both validations fail, report errors
    if not is_airr_valid and olmsted_errors:
        if verbose:
            errors.append(f"AIRR validation: {airr_error}")
            errors.extend(olmsted_errors)
        else:
            errors.append(
                "Tree does not conform to AIRR or Olmsted schema (use -v for details)"
            )
    elif olmsted_errors:
        # Only Olmsted validation failed
        errors.extend(olmsted_errors)
    # If AIRR validation passed OR Olmsted validation passed, consider it valid (no errors)

    return errors


def validate_consolidated_data(data, verbose=False):
    """
    Validate consolidated data format containing metadata, datasets, clones, and trees.

    Args:
        data: Consolidated data dictionary
        verbose: Show detailed errors

    Returns:
        list: List of validation errors (empty if valid)
    """
    errors = []

    # Check top-level structure
    required_keys = ["metadata", "datasets", "clones", "trees"]
    for key in required_keys:
        if key not in data:
            errors.append(f"Missing required key: {key}")

    if errors:
        return errors

    # Validate metadata
    metadata = data.get("metadata", {})
    required_metadata_keys = [
        "format_version",
        "schema_version",
        "created_at",
        "source_format",
    ]
    for key in required_metadata_keys:
        if key not in metadata:
            errors.append(f"Missing required metadata key: {key}")

    # Validate format version compatibility
    format_version = metadata.get("format_version")
    if format_version and format_version != CONSOLIDATED_JSON_VERSION:
        errors.append(f"Unsupported format version: {format_version}")

    # Validate datasets
    datasets = data.get("datasets", [])
    if not isinstance(datasets, list):
        errors.append("'datasets' must be a list")
    else:
        for i, dataset in enumerate(datasets):
            dataset_errors = validate_dataset(dataset, verbose)
            if dataset_errors:
                errors.extend([f"Dataset {i}: {e}" for e in dataset_errors])

    # Validate clones
    clones_dict = data.get("clones", {})
    if not isinstance(clones_dict, dict):
        errors.append("'clones' must be a dictionary")
    else:
        for dataset_id, clones in clones_dict.items():
            if not isinstance(clones, list):
                errors.append(f"Clones for dataset '{dataset_id}' must be a list")
                continue
            for i, clone in enumerate(clones):
                clone_errors = validate_clone(clone, verbose)
                if clone_errors:
                    errors.extend(
                        [f"Clone {dataset_id}[{i}]: {e}" for e in clone_errors]
                    )

    # Validate trees
    trees = data.get("trees", [])
    if not isinstance(trees, list):
        errors.append("'trees' must be a list")
    else:
        for i, tree in enumerate(trees):
            tree_errors = validate_tree(tree, verbose)
            if tree_errors:
                errors.extend([f"Tree {i}: {e}" for e in tree_errors])

    return errors
