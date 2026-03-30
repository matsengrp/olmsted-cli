#!/usr/bin/env python

from __future__ import annotations

import argparse
import functools
import html
import json
import os
import pprint
import sys
import traceback
import uuid
from collections import OrderedDict
from typing import TYPE_CHECKING, Any, Dict, List, Optional
from urllib.parse import parse_qs, parse_qsl

if TYPE_CHECKING:
    from argparse import Namespace

    from .types import OlmstedClone, OlmstedDataset, OlmstedTree

import jsonschema
import ntpl
import yaml

# Python 3.13+ compatibility: make cgi module available before ete3 import
try:
    import cgi  # noqa: F401
except ImportError:
    # Create a mock cgi module using our compatibility layer

    class CGIModule:
        """Mock cgi module for Python 3.13+ compatibility."""

        escape = html.escape

        # Add other cgi functions that might be needed by ete3
        def parse_qs(self, *args, **kwargs):
            return parse_qs(*args, **kwargs)

        def parse_qsl(self, *args, **kwargs):
            return parse_qsl(*args, **kwargs)

    # Make cgi available as a module
    sys.modules["cgi"] = CGIModule()

import ete3

from .process_utils import (
    SCHEMA_VERSION,
    create_consolidated_data,
    dict_subset,
    is_nullable_string,
    merge,
    validate_dataset,
    validate_output_data,
    write_out,
)
from .schemas import (
    clone_spec,
    dataset_spec,
)

type_checker = jsonschema.Draft4Validator.TYPE_CHECKER.redefine(
    "string", is_nullable_string
)
CustomValidator = jsonschema.validators.extend(
    jsonschema.Draft4Validator, type_checker=type_checker
)

# Should update to get draft7?
olmsted_dataset_schema = jsonschema.Draft4Validator(dataset_spec)
airr_clone_schema = None
try:
    with open("../airr-standards/specs/airr-schema.yaml") as stream:
        airr_clone_schema_dict = yaml.load(stream, Loader=yaml.FullLoader).get("Clone")
        airr_clone_schema = CustomValidator(airr_clone_schema_dict)
except FileNotFoundError:
    # AIRR schema file not found - skip AIRR validation
    pass


def ensure_ident(record, prefix=""):
    "Want to let people choose their own uuids if they like, but not require them to"
    if record.get("ident"):
        return record
    uuid_str = str(uuid.uuid4())
    ident_value = f"{prefix}{uuid_str}" if prefix else uuid_str
    return merge(record, {"ident": ident_value})


# reroot the tree on node matching regex pattern.
# Usually this is used to root on the naive germline sequence
# NOTE duplicates fcn in plot_tree.py
# TODO this is just one way to "reroot" trees; it's worth considering removing this function from the script so that we are not responsible for this job since it isn't trivial (e.g. if given an unrooted tree, ete3.Tree.set_outgroup will add an empty-string-named taxon)
def reroot_tree(args, tree):
    # find naive node
    node = tree.search_nodes(name=args.naive_name)[0]
    # if equal, then the root is already the naive, so done
    if tree != node:
        # In general this would be necessary, but we are actually assuming that naive has been set as an
        # outgroup in dnaml, and if it hasn't, we want to raise an error, as below
        tree.set_outgroup(node)
        # This actually assumes the `not in` condition above, but we check as above for clarity
        tree.remove_child(node)
        node.add_child(tree)
        tree.dist = node.dist

        node.dist = 0
        tree = node
    return tree


def process_tree_nodes(args, tree, nodes, reroot=False):
    if reroot:
        tree = reroot_tree(args, tree)

    def process_node(node):
        datum = nodes.get(node.name, {})
        datum["type"] = "leaf" if node.is_leaf() else "node"
        datum.update(nodes.get(node.name, {}))
        if node.up:
            datum["parent"] = node.up.name
            datum["length"] = node.get_distance(node.up)
            # Only calculate distance to naive if rooting is enabled and naive exists
            if reroot and args.naive_name and tree.search_nodes(name=args.naive_name):
                datum["distance"] = node.get_distance(args.naive_name)
            else:
                datum["distance"] = node.get_distance(tree)  # distance to root
        else:
            # node is root
            datum["type"] = "root"
            datum["parent"] = None
            datum["length"] = 0.0
            datum["distance"] = 0.0
        return datum

    return list(map(process_node, tree.traverse("postorder")))


def process_tree(args, clone_id, tree):
    # add clone_id to satisfy AIRR schema
    tree["clone_id"] = clone_id
    ete_tree = ete3.PhyloTree(tree["newick"], format=1)
    tree["nodes"] = process_tree_nodes(
        args, ete_tree, tree["nodes"], reroot=args.root_trees
    )
    return ensure_ident(tree, prefix="tree-")


def process_clone(args, dataset, clone):
    # -=1 *_start positions since AIRR schema uses 1-based closed interval
    # but we need python slice conventions (0-based, open interval) for
    # source code (vega visualization). Gracefully skip missing positions.
    _missing_fields = []
    for start_pos_key in [
        "v_alignment_start",
        "d_alignment_start",
        "j_alignment_start",
        "junction_start",
    ]:
        if start_pos_key in clone and clone[start_pos_key] is not None:
            clone[start_pos_key] -= 1
        else:
            _missing_fields.append(start_pos_key)

    if _missing_fields and getattr(args, "verbose", 0) >= 2:
        clone_id = clone.get("clone_id", "unknown")
        print(
            f"  Note: clone '{clone_id}' missing position fields: {_missing_fields}",
            file=sys.stderr,
        )

    # need to create a copy of the dataset without clonal families that we
    # can nest under clonal family for viz convenience
    _dataset = dataset.copy()
    del _dataset["clones"]
    clone["dataset"] = _dataset

    # Match sample by sample_id; gracefully handle missing match
    matching_samples = [
        s for s in clone["dataset"].get("samples", [])
        if s.get("sample_id") == clone.get("sample_id")
    ]
    if matching_samples:
        clone["sample"] = matching_samples[0]
    else:
        clone["sample"] = {
            "sample_id": clone.get("sample_id", "unknown"),
            "locus": "igh",
        }
        if getattr(args, "verbose", 0) >= 1:
            print(
                f"  Note: clone '{clone.get('clone_id', '?')}' sample_id "
                f"'{clone.get('sample_id')}' not found in dataset samples",
                file=sys.stderr,
            )

    if "samples" in clone.get("dataset", {}):
        del clone["dataset"]["samples"]
    return ensure_ident(clone, prefix="clone-")


def process_dataset(
    args: Namespace,
    dataset: Dict[str, Any],
    clones_dict: Dict[str, List[OlmstedClone]],
    trees: List[OlmstedTree],
) -> Optional[OlmstedDataset]:
    """
    Process a dataset from AIRR format to Olmsted format.

    Args:
        args: Command line arguments namespace
        dataset: AIRR dataset dictionary
        clones_dict: Dictionary to populate with clones keyed by dataset_id
        trees: List to populate with processed trees

    Returns:
        Processed dataset in Olmsted format, or None if processing fails
    """
    dataset["clone_count"] = len(dataset["clones"])
    dataset["subjects_count"] = len(
        set(cf.get("subject_id", "unknown") for cf in dataset["clones"])
    )
    dataset["timepoints_count"] = len(
        set(sample.get("timepoint_id", "unknown") for sample in dataset.get("samples", []))
    )
    clones = list(
        map(functools.partial(process_clone, args, dataset), dataset["clones"])
    )

    # Process trees for each clone and set clone_id
    for cf in clones:
        # Add repertoire_id field (using sample_id)
        cf["repertoire_id"] = cf["sample_id"]

        # Process each tree and set clone_id
        processed_trees = []
        for tree in cf["trees"]:
            processed_tree = process_tree(args, cf["clone_id"], tree)
            processed_trees.append(processed_tree)

        # Compute metrics on trees if requested (LBI, LBR, scaled_affinity)
        if getattr(args, "compute_metrics", False):
            from .metrics import compute_tree_metrics

            lbi_tau = getattr(args, "lbi_tau", 0.0125)
            for tree in processed_trees:
                nodes = tree.get("nodes", {})
                if not nodes:
                    continue
                # Build nodes_dict and edges from the tree's node data
                if isinstance(nodes, list):
                    nodes_dict = {n["sequence_id"]: n for n in nodes}
                else:
                    nodes_dict = nodes
                edges = []
                root_id = None
                for nid, ndata in nodes_dict.items():
                    parent = ndata.get("parent")
                    if parent is None:
                        root_id = nid
                    else:
                        length = ndata.get("length", 0.0) or 0.0
                        edges.append((parent, nid, length))
                if root_id is None:
                    continue

                compute_tree_metrics(nodes_dict, edges, root_id, tau=lbi_tau)

                # Write back if nodes was a list
                if isinstance(tree["nodes"], list):
                    tree["nodes"] = list(nodes_dict.values())

        # Add processed trees to the main trees list
        trees.extend(processed_trees)

        # Keep tree references in clones but remove nodes for size
        cf["trees"] = [
            dict_subset(tree, set(tree.keys()) - {"nodes"}) for tree in processed_trees
        ]
    clones_dict[dataset["dataset_id"]] = clones

    # Generate field_metadata from actual clone and tree data
    from .field_metadata import generate_field_metadata

    dataset["field_metadata"] = generate_field_metadata(
        clones,
        trees,
        custom_fields=getattr(args, "custom_fields", None),
    )

    del dataset["clones"]
    dataset["schema_version"] = SCHEMA_VERSION
    return ensure_ident(dataset, prefix="dataset-")


def hiccup_rep(schema, depth=1, property=None):
    depth = min(depth, 2)
    if depth == 1 or schema["type"] == "object":
        style = (
            "padding-left: 10;"
            + "margin-left: 25;"
            + "margin-top: 40;"
            + "border-left-style: solid;"
            + "border-color: grey;"
        )
    else:
        style = "padding-left: 10;" + "margin-left: 25;" + "margin-top: 10;"
    return [
        "div",
        {"style": style},
        ["h" + str(depth), schema.get("title")] if schema.get("title") else "",
        ["p", ["b", "Description: "], ["span", schema.get("description")]]
        if schema.get("description")
        else "",
        ["p", ["b", "Required: "], ["code", str(schema.get("required"))]]
        if schema.get("required")
        else "",
        ["p", ["b", "Type: "], ["code", str(schema.get("type"))]]
        if schema.get("type")
        else "",
        ["div", ["h" + str(depth + 1), "Properties:"]]
        + [
            [
                "div",
                {"style": "margin-left: 10px;"},
                ["h3", ["code", k]],
                # Assume val is either a title, as produced in hiccup_rep2, or an actual schema
                ["b", {"style": "padding-left: 15; font-size: 18;"}, "{%s}" % val]
                if isinstance(val, str)
                else hiccup_rep(val, depth=depth + 1),
            ]
            for k, val in schema.get("properties").items()
        ]
        if schema.get("properties")
        else "",
        [
            "div",
            ["h" + str(depth + 1), "Array Items:"],
            # As above, assume and display a title if string, otherwise recurse
            [
                "b",
                {"style": "padding-left: 15; font-size: 18;"},
                "{%s}" % schema["items"],
            ]
            if isinstance(schema.get("items"), str)
            else hiccup_rep(schema.get("items"), depth=depth + 1),
        ]
        if schema.get("items")
        else "",
        [
            "div",
            ["h" + str(depth + 1), "Object with values of type:"],
            # As above, assume and display a title if string, otherwise recurse
            [
                "b",
                {"style": "padding-left: 15; font-size: 18;"},
                "{%s}" % schema["additionalProperties"],
            ]
            if isinstance(schema.get("additionalProperties"), str)
            else hiccup_rep(schema.get("additionalProperties"), depth=depth + 1),
        ]
        if schema.get("additionalProperties")
        else "",
    ]


def hiccup_rep2(schema):
    def flatten_schema_by_title(schema):
        items_schemas, properties_schemas = [], []
        items = schema.get("items")
        # if this is an array, check title
        if items and items.get("title"):
            schema["items"] = items["title"]
            items_schemas = flatten_schema_by_title(items)
        # object
        additionalProperties = schema.get("additionalProperties")
        if additionalProperties and additionalProperties.get("title"):
            schema["additionalProperties"] = additionalProperties["title"]
            items_schemas = flatten_schema_by_title(additionalProperties)
        for key, subschema in schema.get("properties", {}).items():
            # handle case of being a single reference, with a title
            title = subschema.get("title")
            if title:
                properties_schemas += flatten_schema_by_title(subschema)
                schema["properties"][key] = title
            # handle array/items case
            items = subschema.get("items")
            if items and items.get("title"):
                properties_schemas += flatten_schema_by_title(items)
                subschema["items"] = items["title"]
            # object
            additionalProperties = subschema.get("additionalProperties")
            if additionalProperties and additionalProperties.get("title"):
                properties_schemas += flatten_schema_by_title(additionalProperties)
                subschema["additionalProperties"] = additionalProperties["title"]
        return list(
            OrderedDict(
                [
                    (schema["title"], schema)
                    for schema in [schema] + items_schemas + properties_schemas
                ]
            ).values()
        )

    return ["div", list(map(hiccup_rep, flatten_schema_by_title(schema)))]


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--inputs", nargs="+")
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output file path for consolidated JSON (default behavior)",
    )
    parser.add_argument(
        "--split-files",
        metavar="DIR",
        dest="data_outdir",
        help="Output to multiple files in specified directory (datasets.json, clones.*.json, tree.*.json) instead of single consolidated file",
    )
    parser.add_argument("--naive-name", default="naive")
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument(
        "-c",
        "--remove-invalid-clones",
        action="store_true",
        help="validate clones individually against the olmsted schema, removing the invalid ones and try to build the dataset using the remaining clones. Note that processing can still be crashed by clones which are invalid according to the AIRR clones and trees schema (see airr-standards/specs/airr-schema.yaml).",
    )
    parser.add_argument("-S", "--display-schema-html")
    parser.add_argument(
        "-s",
        "--display-schema",
        action="store_true",
        help="print schema to stdout for display",
    )
    parser.add_argument(
        "-y",
        "--write-schema-yaml",
        action="store_true",
        help="write the schema to a yaml format file.",
    )
    parser.add_argument(
        "-r", "--root-trees", action="store_true", help="Root trees using --naive-name."
    )
    parser.add_argument(
        "--validate",
        action="store_true",
        help="Validate output data against AIRR JSON schemas before writing",
    )
    parser.add_argument(
        "--strict-validation",
        action="store_true",
        help="Exit with error if validation fails (requires --validate)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        help="Random seed for deterministic processing (currently unused for AIRR format, added for API consistency)",
    )
    return parser.parse_args()


def main():
    args = get_args()
    datasets, clones_dict, trees = [], {}, []
    for infile in args.inputs or []:
        print(f"\nProcessing infile: {str(infile)}")
        try:
            with open(infile, "r") as fh:
                dataset = json.load(fh)
                if args.remove_invalid_clones:
                    dataset["clones"] = list(
                        filter(
                            jsonschema.Draft4Validator(clone_spec).is_valid,
                            dataset["clones"],
                        )
                    )
                # Use unified validation from validate module
                errors = validate_dataset(dataset, verbose=args.verbose)
                if errors:
                    error_msg = "Dataset validation failed"
                    if args.verbose:
                        print(f"Dataset validation failed:")
                        for error in errors:
                            print(f"  - {error}")
                    else:
                        error_msg += ". Please rerun with `-v` for detailed errors"
                    raise Exception(error_msg)
                # Process the dataset, including validation of clones, trees against the AIRR schema
                dataset = process_dataset(args, dataset, clones_dict, trees)
                datasets.append(dataset)
        except Exception:
            print(f"Unable to process infile: {infile}")
            if args.verbose:
                exc_info = sys.exc_info()
                traceback.print_exception(*exc_info)
            else:
                print("Please rerun with `-v` for detailed errors.")
            sys.exit(1)
    # write out schema
    if args.write_schema_yaml:
        with open("schema.yaml", "w") as yamlf:
            yaml.dump(dataset_spec, yamlf)
    if args.display_schema:
        pprint.pprint(dataset_spec)
    if args.display_schema_html:
        with open(args.display_schema_html, "w") as fh:
            fh.write(
                ntpl.render(
                    [
                        "html",
                        [
                            "body",
                            hiccup_rep2(dataset_spec),
                        ],
                    ]
                )
            )
    # Validate data before writing if requested
    if args.validate:
        if not validate_output_data(datasets, clones_dict, trees, args):
            if args.strict_validation:
                print(
                    "\nExiting due to validation errors (--strict-validation enabled)"
                )
                sys.exit(1)

    # write out data
    if args.data_outdir:
        # Multi-file output to specified directory
        write_out(datasets, args.data_outdir, "datasets.json", args)
        for dataset_id, clones in clones_dict.items():
            write_out(
                clones, args.data_outdir + "/", "clones." + dataset_id + ".json", args
            )
        for tree in trees:
            write_out(
                tree, args.data_outdir + "/", "tree." + tree["ident"] + ".json", args
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
        write_out(consolidated_data, output_dir, output_file, args)


if __name__ == "__main__":
    main()
