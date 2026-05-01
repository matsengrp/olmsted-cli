#!/usr/bin/env python
"""
Processing utilities for Olmsted data pipelines.

This module contains functions that depend on other project modules
(schemas, build_config, field_metadata).  Pure utilities with no
project dependencies live in ``utils.py``.
"""

import csv
import json
import os
import uuid
from datetime import datetime, timezone
import jsonschema
import yaml
from tqdm import tqdm

from .build_config import generate_default_config
from .data_io import write_olmsted_json
from .field_metadata import generate_field_metadata
from .schemas import SCHEMA_VERSION, clone_spec, dataset_spec, tree_spec
from .utils import (  # noqa: F401 — re-exported for backward compatibility
    VerbosePrinter,
    add_verbosity_args,
    clean_record,
    comp,
    dict_subset,
    get_in,
    get_optional_int,
    inf,
    is_nullable_string,
    json_rep,
    listof,
    listofint,
    merge,
    natural_number,
    neginf,
    remap_dict_values,
    remap_list,
    rename_keys,
    resolve_verbosity,
    set_verbosity,
    strip_ns,
    translate_dna_to_aa,
    try_del,
    vprint,
)
from .version import __version__, get_git_hash


# VerbosePrinter and all general-purpose utilities now live in utils.py.
# They are re-exported above for backward compatibility.


def coerce_csv_value(val: str):
    """Coerce a CSV string value to the most appropriate Python type.

    Attempts in order: int → float → JSON (list/dict) → bool → string.
    """
    try:
        return int(val)
    except (ValueError, TypeError):
        pass
    try:
        return float(val)
    except (ValueError, TypeError):
        pass
    if val.startswith(("[", "{")):
        try:
            return json.loads(val)
        except (json.JSONDecodeError, ValueError):
            pass
    if val.lower() in ("true", "false"):
        return val.lower() == "true"
    return val


def write_out(data, dirname, filename, args):
    """
    Write data to JSON or CSV file with proper formatting and UUID handling.

    Args:
        data: Data to write
        dirname: Directory path
        filename: File name
        args: Command line arguments (for verbose flag, csv flag, and json_format)
    """
    # Ensure directory exists
    os.makedirs(dirname, exist_ok=True)

    # Normalize path
    full_path = os.path.normpath(os.path.join(dirname, filename))

    # Get JSON format setting (default to 'pretty' for backward compatibility)
    json_format = getattr(args, "json_format", "pretty")

    # For gzip format, add .gz extension if not already present
    if json_format == "gzip" and not full_path.endswith(".gz"):
        full_path = full_path + ".gz"

    # Print status
    vprint.status(f"writing {full_path}")

    # Check if CSV output is requested (for CFT data)
    if hasattr(args, "csv") and args.csv and isinstance(data, list):
        # Write as CSV
        with open(full_path, "w") as fh:
            if data:
                # Ensure all items are dictionaries
                data = [{k: v for k, v in d.items()} for d in data]
                writer = csv.DictWriter(fh, fieldnames=sorted(data[0].keys()))
                writer.writeheader()
                writer.writerows(data)
    elif isinstance(data, (list, dict)):
        # Write as JSON with selected format
        write_olmsted_json(data, full_path, json_format=json_format, default=json_rep)
    else:
        # Handle raw string data
        with open(full_path, "w") as fh:
            fh.write(data)


# Version Constants
CONSOLIDATED_JSON_VERSION = "1.0"


def unpack_encoded_mutations(trees, custom_fields):
    """Unpack encoded mutation-level fields from nodes into mutations arrays.

    Processes custom_fields with an ``encoding`` key, reading data from node-level
    fields and merging it into each node's ``mutations`` array (keyed by ``site``).

    Supported encodings:
        - list: Dense per-position array on node. Index = site, value = field value.
          Null values are skipped.
        - json: Sparse dict on node. Int key = site, value = field value.
          Null values are skipped.
        - records: Array of dicts with ``site`` key on node. ``source`` names
          the node field containing the array; ``name`` is the inner field to extract.

    Args:
        trees: List of tree dicts (modified in place).
        custom_fields: List of custom field declarations from config.
    """
    encoded = [cf for cf in (custom_fields or []) if cf.get("encoding")]
    if not encoded:
        return

    # Group records fields by source for efficient single-pass merging
    records_by_source = {}
    for cf in encoded:
        if cf["encoding"] == "records":
            source = cf["source"]
            records_by_source.setdefault(source, []).append(cf["name"])

    for tree in trees:
        nodes = tree.get("nodes", [])
        if isinstance(nodes, dict):
            node_list = list(nodes.values())
        else:
            node_list = nodes

        for node in node_list:
            if not isinstance(node, dict):
                continue

            # Index existing mutations by site
            existing = node.get("mutations", [])
            by_site = {}
            for m in existing:
                if isinstance(m, dict) and "site" in m:
                    by_site[m["site"]] = m

            for cf in encoded:
                field_name = cf["name"]
                encoding = cf["encoding"]

                if encoding == "list":
                    data = node.get(field_name)
                    if not isinstance(data, list):
                        continue
                    for site, val in enumerate(data):
                        if val is None:
                            continue
                        if site not in by_site:
                            by_site[site] = {"site": site}
                        by_site[site][field_name] = val

                elif encoding == "json":
                    data = node.get(field_name)
                    if not isinstance(data, dict):
                        continue
                    for key, val in data.items():
                        if val is None:
                            continue
                        try:
                            site = int(key)
                        except (ValueError, TypeError):
                            continue
                        if site not in by_site:
                            by_site[site] = {"site": site}
                        by_site[site][field_name] = val

                elif encoding == "records":
                    source = cf["source"]
                    # Only process the source array once per node (first field triggers it)
                    if field_name != records_by_source[source][0]:
                        continue
                    data = node.get(source)
                    if not isinstance(data, list):
                        continue
                    fields_to_extract = records_by_source[source]
                    for entry in data:
                        if not isinstance(entry, dict) or "site" not in entry:
                            continue
                        site = entry["site"]
                        if site not in by_site:
                            by_site[site] = {"site": site}
                        for fname in fields_to_extract:
                            if fname in entry:
                                by_site[site][fname] = entry[fname]

            # Write back sorted by site
            if by_site:
                node["mutations"] = sorted(by_site.values(), key=lambda m: m["site"])


def tag_field_metadata(clones, trees, custom_fields=None):
    """Generate field_metadata for a dataset, applying default config if needed.

    This is the shared entry point used by ``process`` and ``tag`` to
    produce field_metadata.  It ensures both commands use
    ``generate_default_config`` as the single source of truth for field
    discovery when no explicit config is provided.

    Args:
        clones: List of clone dicts for the dataset.
        trees: List of tree dicts for the dataset (modified in place by
            unpack_encoded_mutations).
        custom_fields: Optional list of custom field declarations from a
            user-provided config.  When *None*, defaults are generated
            via ``generate_default_config``.

    Returns:
        Dict with level keys mapping to field metadata dicts, suitable
        for assigning to ``dataset["field_metadata"]``.
    """
    if custom_fields is None:
        custom_fields = generate_default_config(clones, trees)

    unpack_encoded_mutations(trees, custom_fields)

    return generate_field_metadata(clones, trees, custom_fields=custom_fields)


def retag_datasets_field_metadata(
    datasets, clones_dict, trees, custom_fields=None, mode="add"
):
    """Recompute and attach ``field_metadata`` to every dataset in place.

    For each dataset this:
      1. Collects its clones (from ``clones_dict``) and the trees whose
         ``clone_id`` matches one of those clones.
      2. Calls ``tag_field_metadata`` to regenerate the metadata.
      3. Either overwrites the existing ``field_metadata`` (``mode="overwrite"``)
         or merges new entries on top of existing ones (``mode="add"``, the
         default) so that pre-existing per-level entries are preserved.

    Args:
        datasets: List of dataset dicts (modified in place).
        clones_dict: Mapping of ``dataset_id -> list of clone dicts``.
        trees: Flat list of tree dicts across all datasets.
        custom_fields: Optional list of custom field declarations.
        mode: ``"add"`` merges with existing metadata, ``"overwrite"`` replaces.
    """
    if mode not in ("add", "overwrite"):
        raise ValueError(f"mode must be 'add' or 'overwrite', got {mode!r}")

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
        clone_ids = {c.get("clone_id") for c in dataset_clones if c.get("clone_id")}
        dataset_trees = [
            t for cid in clone_ids for t in trees_by_clone_id.get(cid, [])
        ]

        new_field_metadata = tag_field_metadata(
            dataset_clones, dataset_trees, custom_fields
        )

        if mode == "overwrite":
            dataset["field_metadata"] = new_field_metadata
            continue

        # Add mode: merge with existing. New entries overwrite same-named fields;
        # untouched levels/fields are preserved.
        existing_metadata = dataset.get("field_metadata", {})
        merged = {}
        all_levels = set(existing_metadata.keys()) | set(new_field_metadata.keys())
        for level in all_levels:
            merged_level = dict(existing_metadata.get(level, {}))
            merged_level.update(new_field_metadata.get(level, {}))
            if merged_level:
                merged[level] = merged_level
        dataset["field_metadata"] = merged


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
    # Count total leaf nodes across all trees
    total_leaf_count = 0
    for tree in trees:
        if "nodes" in tree and tree["nodes"]:
            for node in tree["nodes"]:
                if node.get("type") == "leaf":
                    total_leaf_count += 1

    metadata = {
        "format": "olmsted",
        "format_version": CONSOLIDATED_JSON_VERSION,
        "schema_version": SCHEMA_VERSION,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "source_format": detected_format,
        "source_files": [os.path.basename(f) for f in input_files],
        "processing_info": {
            "datasets_count": len(datasets),
            "total_clones_count": sum(len(clones) for clones in clones_dict.values()),
            "total_trees_count": len(trees),
            "total_leaf_nodes_count": total_leaf_count,
        },
        "generated_by": {
            "tool": "olmsted-cli",
            "version": __version__,
            "git_hash": get_git_hash(),
        },
    }

    # Add optional name and description if provided
    if args and hasattr(args, "name") and args.name:
        metadata["name"] = args.name
    if args and hasattr(args, "description") and args.description:
        metadata["description"] = args.description

    # Add processing options if available
    if args:
        metadata["processing_options"] = {
            "validation": getattr(args, "validate", False),
            "strict_validation": getattr(args, "strict_validation", False),
            "seed": getattr(args, "seed", None),
        }

        # Add format-specific options
        if detected_format == "airr":
            root_arg = getattr(args, "root", None)
            metadata["processing_options"]["airr"] = {
                "root": root_arg,
            }

    return {
        "metadata": metadata,
        "datasets": datasets,
        "clones": clones_dict,
        "trees": trees,
    }


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
        vprint.error("Warning: Official AIRR schema not found")
        return None
    except Exception as e:
        vprint.error(f"Warning: Failed to load official AIRR schema: {e}")
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

    vprint.status("\nValidating output data against schemas...")

    validation_passed = True
    total_errors = 0

    try:
        # Validate datasets
        for i, dataset in enumerate(datasets):
            errors = validate_dataset(dataset, verbose=getattr(args, "verbose", False))
            if errors:
                vprint.error(f"FAIL: Dataset {i} validation failed:")
                for error in errors:
                    vprint.error(f"  - {error}")
                validation_passed = False
                total_errors += len(errors)
            elif getattr(args, "verbose", False):
                vprint.status(f"PASS: Dataset {i} validation passed")

        # Validate clones
        clone_count = 0
        clone_failures = 0
        
        # Count total clones for progress bar
        total_clones = sum(len(clones) for clones in clones_dict.values())
        
        with tqdm(total=total_clones, desc="Validating clones", unit="clone", disable=total_clones <= 1) as pbar:
            for dataset_id, clones in clones_dict.items():
                for clone in clones:
                    clone_count += 1
                    clone_id = clone.get('clone_id', 'unknown')
                    pbar.set_description(f"Validating clone {clone_id}")
                    
                    errors = validate_clone(clone, verbose=getattr(args, "verbose", False))
                    if errors:
                        clone_failures += 1
                        if getattr(args, "verbose", False):
                            vprint.error(
                                f"FAIL: Clone {clone_id} validation failed:"
                            )
                            for error in errors:
                                vprint.error(f"  - {error}")
                        validation_passed = False
                        total_errors += len(errors)

                    pbar.update(1)

        if clone_failures == 0:
            vprint.status(f"PASS: Clone validation passed ({clone_count} clones)")
        else:
            vprint.error(f"FAIL: Clone validation: {clone_failures}/{clone_count} failed")

        # Validate trees
        tree_count = 0
        tree_failures = 0
        
        with tqdm(trees, desc="Validating trees", unit="tree", disable=len(trees) <= 1) as pbar:
            for tree in pbar:
                tree_count += 1
                tree_id = tree.get('ident', 'unknown')
                pbar.set_description(f"Validating tree {tree_id}")
                
                # Check if time tree validation is enabled
                check_time_tree = getattr(args, 'time_tree', False)
                errors = validate_tree(tree, verbose=getattr(args, "verbose", False), check_time_tree=check_time_tree)
                if errors:
                    tree_failures += 1
                    if getattr(args, "verbose", False):
                        vprint.error(f"FAIL: Tree {tree_id} validation failed:")
                        for error in errors:
                            vprint.error(f"  - {error}")
                    validation_passed = False
                    total_errors += len(errors)

        if tree_failures == 0:
            vprint.status(f"PASS: Tree validation passed ({tree_count} trees)")
        else:
            vprint.error(f"FAIL: Tree validation: {tree_failures}/{tree_count} failed")

        if total_errors > 0:
            vprint.error(f"\nTotal validation errors: {total_errors}")

    except Exception as e:
        vprint.error(f"Validation error: {str(e)}")
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


def validate_time_tree(nodes, verbose=False):
    """
    Validate that a tree is a valid time tree.
    
    For time trees, each child node's distance from root should be 
    greater than or equal to its parent's distance from root.
    
    Args:
        nodes: List of node dictionaries with 'sequence_id', 'parent', and 'distance' fields
        verbose: Show detailed errors
        
    Returns:
        list: List of validation errors (empty if valid)
    """
    errors = []
    
    if not nodes:
        return errors  # Empty tree is valid
    
    # Build a dictionary for quick node lookup
    node_dict = {}
    for node in nodes:
        if isinstance(node, dict) and 'sequence_id' in node:
            node_dict[node['sequence_id']] = node
    
    # Check each node's distance relationship with its parent
    for node in nodes:
        if not isinstance(node, dict):
            continue
            
        node_id = node.get('sequence_id')
        parent_id = node.get('parent')
        node_distance = node.get('distance')
        
        # Skip if no parent (root node) or missing data
        if not parent_id or parent_id == 'null' or node_distance is None:
            continue
            
        # Find parent node
        parent_node = node_dict.get(parent_id)
        if not parent_node:
            if verbose:
                errors.append(f"Node {node_id}: parent {parent_id} not found in tree")
            continue
            
        parent_distance = parent_node.get('distance')
        if parent_distance is None:
            continue
            
        # Check time tree constraint
        try:
            if float(node_distance) < float(parent_distance):
                error_msg = (f"Time tree violation: Node {node_id} has distance {node_distance} "
                           f"which is less than parent {parent_id} distance {parent_distance}")
                errors.append(error_msg)
        except (ValueError, TypeError):
            if verbose:
                errors.append(f"Node {node_id}: non-numeric distance value")
    
    return errors


def validate_tree(data, verbose=False, check_time_tree=False):
    """
    Validate a tree against AIRR and Olmsted schemas.

    Args:
        data: Tree dictionary
        verbose: Show detailed errors
        check_time_tree: Whether to validate time tree constraints (default: False)

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
    
    # Check tree has a root node
    if "nodes" in data:
        nodes = data["nodes"]
        if isinstance(nodes, list):
            nodes_list = nodes
        elif isinstance(nodes, dict):
            nodes_list = list(nodes.values())
        else:
            nodes_list = []

        if nodes_list:
            # Find root node(s) — nodes with parent=None or type="root"
            root_nodes = [
                n for n in nodes_list
                if isinstance(n, dict) and (
                    n.get("parent") is None or n.get("type") == "root"
                )
            ]

            if not root_nodes:
                # No root found — check for common naive/germline node names
                all_names = [n.get("sequence_id", "") for n in nodes_list if isinstance(n, dict)]
                naive_candidates = [
                    name for name in all_names
                    if name and any(
                        hint in name.lower()
                        for hint in ("naive", "germline", "inferred_naive", "root", "uca")
                    )
                ]
                msg = "Tree has no root node (no node with parent=null or type='root')"
                if naive_candidates:
                    msg += (
                        f". Found candidate root node(s): {naive_candidates[:3]}. "
                        f"Try reprocessing with: --root {naive_candidates[0]}"
                    )
                errors.append(msg)
            elif len(root_nodes) > 1:
                root_ids = [n.get("sequence_id", "?") for n in root_nodes]
                errors.append(
                    f"Tree has multiple root nodes ({len(root_nodes)}): {root_ids[:5]}"
                )
            else:
                if verbose >= 2:
                    root_id = root_nodes[0].get("sequence_id", "?")
                    vprint.verbose(f"  Root node: {root_id}")

    # Check time tree constraints if requested and nodes are present
    if check_time_tree and "nodes" in data and isinstance(data["nodes"], list):
        time_tree_errors = validate_time_tree(data["nodes"], verbose=verbose)
        if time_tree_errors:
            errors.extend(time_tree_errors)

    return errors


def validate_consolidated_data(data, verbose=False, check_time_tree=False):
    """
    Validate consolidated data format containing metadata, datasets, clones, and trees.

    Args:
        data: Consolidated data dictionary
        verbose: Show detailed errors
        check_time_tree: Whether to validate time tree constraints

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
        # Count total clones for progress bar
        total_clones = sum(len(clones) if isinstance(clones, list) else 0 for clones in clones_dict.values())
        
        if total_clones > 1:  # Show progress bar if more than 1 clone
            with tqdm(total=total_clones, desc="Validating clones", unit="clone", leave=False) as pbar:
                for dataset_id, clones in clones_dict.items():
                    if not isinstance(clones, list):
                        errors.append(f"Clones for dataset '{dataset_id}' must be a list")
                        continue
                    for i, clone in enumerate(clones):
                        clone_id = clone.get('clone_id', f'{dataset_id}[{i}]') if isinstance(clone, dict) else f'{dataset_id}[{i}]'
                        pbar.set_description(f"Validating clone {clone_id}")
                        clone_errors = validate_clone(clone, verbose)
                        if clone_errors:
                            errors.extend(
                                [f"Clone {dataset_id}[{i}]: {e}" for e in clone_errors]
                            )
                        pbar.update(1)
        else:
            # No progress bar for single clone
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
        if len(trees) > 1:  # Show progress bar if more than 1 tree
            with tqdm(trees, desc="Validating trees", unit="tree", leave=False) as pbar:
                for i, tree in enumerate(pbar):
                    tree_id = tree.get('ident', tree.get('tree_id', f'tree-{i}')) if isinstance(tree, dict) else f'tree-{i}'
                    pbar.set_description(f"Validating tree {tree_id}")
                    tree_errors = validate_tree(tree, verbose, check_time_tree)
                    if tree_errors:
                        errors.extend([f"Tree {i}: {e}" for e in tree_errors])
        else:
            # No progress bar for single tree
            for i, tree in enumerate(trees):
                tree_errors = validate_tree(tree, verbose, check_time_tree)
                if tree_errors:
                    errors.extend([f"Tree {i}: {e}" for e in tree_errors])

    return errors


def _find_duplicates(values):
    """Return a sorted list of values that appear more than once."""
    seen = set()
    dups = set()
    for v in values:
        if v in seen:
            dups.add(v)
        else:
            seen.add(v)
    return sorted(dups)


def check_output_id_uniqueness(datasets, clones_dict, *, allow_duplicates=False):
    """Verify user-facing ``*_id`` uniqueness across the output.

    Catches input-derived or synthesized ``*_id`` values that collide in
    ways that would silently overwrite downstream (Redux state keyed on
    ``*_id``, Dexie ``bulkPut`` upsert semantics).

    Scopes checked:

    - ``dataset.dataset_id`` across ``datasets[]``
    - ``clone.clone_id`` within each dataset's clones
    - ``tree.tree_id`` within each clone's ``trees[]``
    - ``sample.sample_id`` within each ``dataset.samples[]``
    - ``subject.subject_id`` within each ``dataset.subjects[]``

    ``sequence_id`` uniqueness within a tree is already enforced upstream
    by the Newick parser (``process_pcp_data._build_unique_names``) and is
    not rechecked here. The top-level ``trees[]`` list is not an input
    because every tree is also reachable via ``clones_dict[*].trees[]``;
    adding a top-level scope would require a new rule (e.g. global
    ``tree.ident`` uniqueness once the webapp DB migrates to a
    ``tree_id`` primary key) — not the current contract.

    Args:
        datasets: List of dataset dicts.
        clones_dict: ``{dataset_id: [clone, ...]}`` mapping.
        allow_duplicates: If True, violations are printed as warnings and
            the function returns normally. If False (default), raises
            ``ValueError`` listing every violation found.

    Raises:
        ValueError: When duplicates are detected and ``allow_duplicates``
        is False. All violations are reported in one error so users can
        fix input once rather than chasing them one at a time.
    """
    violations = []

    # dataset_id across datasets[]
    dataset_id_dups = _find_duplicates(
        d["dataset_id"] for d in datasets if d.get("dataset_id")
    )
    if dataset_id_dups:
        violations.append(
            f"dataset_id: duplicate across datasets[]: {dataset_id_dups}"
        )

    for dataset in datasets:
        dataset_id = dataset.get("dataset_id", "<unknown>")

        # sample_id within dataset.samples[]
        sample_dups = _find_duplicates(
            s["sample_id"] for s in dataset.get("samples", []) if s.get("sample_id")
        )
        if sample_dups:
            violations.append(
                f"sample_id: duplicate in dataset {dataset_id}.samples[]: {sample_dups}"
            )

        # subject_id within dataset.subjects[]
        subject_dups = _find_duplicates(
            s["subject_id"] for s in dataset.get("subjects", []) if s.get("subject_id")
        )
        if subject_dups:
            violations.append(
                f"subject_id: duplicate in dataset {dataset_id}.subjects[]: {subject_dups}"
            )

    # clone_id within a dataset; tree_id within a clone
    for dataset_id, clones in clones_dict.items():
        clone_id_dups = _find_duplicates(
            c["clone_id"] for c in clones if c.get("clone_id")
        )
        if clone_id_dups:
            violations.append(
                f"clone_id: duplicate within dataset {dataset_id}: {clone_id_dups}"
            )

        for clone in clones:
            tree_id_dups = _find_duplicates(
                t["tree_id"]
                for t in clone.get("trees", [])
                if t.get("tree_id")
            )
            if tree_id_dups:
                violations.append(
                    f"tree_id: duplicate within clone {clone.get('clone_id', '<unknown>')}: "
                    f"{tree_id_dups}"
                )

    if not violations:
        return

    if allow_duplicates:
        vprint.error(
            "Duplicate *_id values detected (--allow-duplicate-ids set, proceeding):"
        )
        for v in violations:
            vprint.error(f"  {v}")
        return

    raise ValueError(
        "Duplicate *_id values detected — downstream consumers key on these and "
        "would silently overwrite on collision:\n  "
        + "\n  ".join(violations)
        + "\n\nFix the input, or pass --allow-duplicate-ids to downgrade this to a warning."
    )
