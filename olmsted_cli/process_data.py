#!/usr/bin/env python
"""
Unified data processing script for Olmsted visualization.

This script can process both AIRR JSON and PCP CSV formats, automatically
detecting the input format or using user-specified format type.

Supported formats:
- AIRR JSON: Standard AIRR format with clones and trees
- PCP CSV: Parent-Child Pair format with optional Newick trees

Output modes:
- Single Olmsted JSON file (default): All data in one file
- Multiple files (--split-files): Separate datasets.json, clones.*.json, tree.*.json
"""

import argparse
import csv
import json
import os
import sys
import traceback
from pathlib import Path
from typing import Any, Dict, List

import jsonschema
import yaml
from tqdm import tqdm

from .constants import (
    DISPLAY_MODES,
    FIELD_LEVELS,
    FIELD_TYPES,
    FORMAT_AIRR,
    FORMAT_AUTO,
    FORMAT_OLMSTED,
    FORMAT_PCP,
    FORMAT_UNKNOWN,
    MUTATION_ENCODINGS,
    normalize_level,
)
from .data_io import detect_file_format, open_file, read_airr_json, read_yaml_config
from .identifier import IdentMinter
from .merge_mutations import (
    apply_mutations_csv,
    apply_mutations_to_trees,
    begin_merge,
    finalize_merge,
    load_mutations_csv,
    report_merge_stats,
)
from .process_airr_data import (
    clone_spec,
    ensure_ident,
    iter_airr_clones,
    process_dataset,
)
from .process_pcp_data import (
    TreeProcessingConfig,
    _group_pcp_families_by_clone,
    iter_pcp_clone_groups,
    parse_newick_csv,
    parse_pcp_csv,
    process_pcp_to_olmsted,
)
from .process_utils import (
    SCHEMA_VERSION,
    add_verbosity_args,
    check_output_id_uniqueness,
    create_consolidated_data,
    resolve_verbosity,
    retag_datasets_field_metadata,
    unpack_encoded_mutations,
    validate_dataset,
    validate_output_data,
    write_out,
)
from .streaming import (
    BatchAccumulator,
    BatchSpooler,
    DuplicateIdError,
    apply_dataset_hoist,
    write_olmsted_json_streaming,
)
from .utils import set_verbosity, vprint


def validate_airr_file(file_path):
    """
    Validate that a file contains valid AIRR JSON data.

    Args:
        file_path: Path to the AIRR JSON file

    Returns:
        bool: True if valid AIRR format, False otherwise
    """
    try:
        handle, _ = open_file(file_path, expected_formats=(FORMAT_AIRR,))
        with handle as f:
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
        handle, _ = open_file(file_path, expected_formats=(FORMAT_PCP,))
        with handle as file_handle:
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

    vprint.status("Processing AIRR format...")

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
    # --root consolidates --naive-name and --root-trees:
    # --root → root at "naive"; --root NAME → root at NAME; not specified → no rooting
    root_arg = getattr(args, "root", None)
    airr_args.naive_name = root_arg if root_arg else "naive"
    airr_args.root_trees = root_arg is not None
    airr_args.remove_invalid_clones = getattr(args, "remove_invalid_clones", False)
    airr_args.display_schema_html = None
    airr_args.display_schema = False
    airr_args.write_schema_yaml = False
    airr_args.compute_metrics = getattr(args, "compute_metrics", False)
    airr_args.lbi_tau = getattr(args, "lbi_tau", 0.0125)
    airr_args.custom_fields = getattr(args, "custom_fields", None)
    airr_args.minter = IdentMinter(seed=getattr(args, "seed", None))
    airr_args.allow_duplicate_ids = getattr(args, "allow_duplicate_ids", False)
    airr_args.json_format = getattr(args, "json_format", "pretty")

    # Streaming pipeline (#26): same fallback conditions as PCP — see
    # _should_stream_airr / _should_stream_pcp.
    if _should_stream_airr(args):
        _process_airr_streaming(args, airr_args)
        return

    # Process using AIRR logic (adapted from process_airr_data.py)
    datasets, clones_dict, trees = [], {}, []

    # Process input files with progress bar
    input_files = airr_args.inputs or []
    with tqdm(
        input_files,
        desc="Processing AIRR files",
        unit="file",
        disable=len(input_files) == 1,
    ) as pbar:
        for infile in pbar:
            pbar.set_description(f"Processing {Path(infile).name}")

            if len(input_files) == 1 or airr_args.verbose:
                vprint.status(f"\nProcessing AIRR file: {infile}")

            try:
                dataset = read_airr_json(infile)

                # Filter invalid clones if requested
                if airr_args.remove_invalid_clones:
                    original_count = len(dataset.get("clones", []))
                    dataset["clones"] = list(
                        filter(
                            jsonschema.Draft4Validator(clone_spec).is_valid,
                            dataset["clones"],
                        )
                    )
                    filtered_count = original_count - len(dataset["clones"])
                    if filtered_count > 0:
                        pbar.set_postfix({"filtered": filtered_count})

                # Use unified validation from validate module
                errors = validate_dataset(dataset, verbose=airr_args.verbose).errors
                if errors:
                    error_msg = "Dataset validation failed"
                    if airr_args.verbose:
                        vprint.error("Dataset validation failed:")
                        for error in errors:
                            vprint.error(f"  - {error}")
                    else:
                        error_msg += ". Please rerun with `-v` for detailed errors"
                    raise Exception(error_msg)

                # Process dataset
                dataset = process_dataset(airr_args, dataset, clones_dict, trees)
                datasets.append(dataset)

                # Update progress bar with clone count
                if "clones" in dataset:
                    pbar.set_postfix({"clones": len(dataset["clones"])})

            except Exception:
                vprint.error(f"\nUnable to process AIRR file: {infile}")
                if airr_args.verbose:
                    exc_info = sys.exc_info()
                    traceback.print_exception(*exc_info)
                else:
                    vprint.error("Please rerun with `-v` for detailed errors.")
                sys.exit(1)

    # Merge mutations CSV if --mutations was specified
    mutations_path = getattr(args, "mutations", None)
    if mutations_path:
        try:
            apply_mutations_csv(
                mutations_path,
                trees,
                use_depth=getattr(args, "mutations_use_depth", False),
                allow_mismatch=getattr(args, "mutations_allow_mismatch", False),
                only_listed=getattr(args, "mutations_listed_only", False),
            )
        except ValueError as e:
            vprint.error(f"Error: {e}")
            sys.exit(1)
        retag_datasets_field_metadata(
            datasets,
            clones_dict,
            trees,
            custom_fields=getattr(args, "custom_fields", None),
        )

    # Enforce *_id uniqueness before writing. --allow-duplicate-ids downgrades
    # collisions to a warning; otherwise we fail fast so silent overwrites
    # downstream are impossible.
    try:
        check_output_id_uniqueness(
            datasets,
            clones_dict,
            allow_duplicates=getattr(args, "allow_duplicate_ids", False),
        )
    except ValueError as e:
        vprint.error(f"Error: {e}")
        sys.exit(1)

    # Validate data before writing if requested
    if airr_args.validate and not validate_output_data(
        datasets, clones_dict, trees, airr_args
    ):
        if airr_args.strict_validation:
            vprint.error(
                "\nExiting due to validation errors (--strict-validation enabled)"
            )
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
        # Olmsted JSON output (default)
        consolidated_data = create_consolidated_data(
            datasets, clones_dict, trees, args.inputs, FORMAT_AIRR, args
        )
        # Ensure output directory exists
        output_dir = os.path.dirname(args.output) or "."
        output_file = os.path.basename(args.output)
        os.makedirs(output_dir, exist_ok=True)
        vprint.status(f"Writing Olmsted JSON output to {args.output}")
        write_out(consolidated_data, output_dir, output_file, airr_args)


def _begin_mutations_merge(args, mutations_path):
    """Load the mutations CSV once and open a :class:`MergeContext`.

    Returns ``(merge_ctx, total_csv_rows)`` for the streaming caller to
    thread through per-batch :func:`apply_mutations_to_trees` calls.
    When no path is supplied, returns ``(None, 0)``.  ``--mutations-use-depth``
    mismatch errors surface here as ``sys.exit(1)`` so the batch loop never
    starts.
    """
    if not mutations_path:
        return None, 0
    vprint.status(f"Loading mutations CSV: {mutations_path}")
    mutations_by_family = load_mutations_csv(mutations_path)
    total_csv_rows = sum(len(rows) for rows in mutations_by_family.values())
    vprint.status(
        f"Loaded {total_csv_rows} CSV rows across {len(mutations_by_family)} families"
    )
    try:
        ctx = begin_merge(
            mutations_by_family,
            use_depth=getattr(args, "mutations_use_depth", False),
        )
    except ValueError as e:
        vprint.error(f"Error: {e}")
        sys.exit(1)
    return ctx, total_csv_rows


def _finalize_mutations_merge(args, merge_ctx, total_csv_rows):
    """Seal a streaming mutations merge: emit the summary lines + raise on mismatch.

    No-op when ``merge_ctx`` is None.  Exits non-zero if any integrity
    mismatch was recorded and ``--mutations-allow-mismatch`` is off.
    """
    if merge_ctx is None:
        return
    stats = finalize_merge(merge_ctx)
    report_merge_stats(stats, total_csv_rows)
    if stats.integrity_mismatches and not getattr(
        args, "mutations_allow_mismatch", False
    ):
        vprint.error(
            f"Error: {stats.integrity_mismatches} integrity mismatches between CSV "
            f"rows and tree mutations. Mismatched rows were skipped "
            f"(never attached). Re-run with --mutations-allow-mismatch to "
            f"proceed anyway — but investigate the CSV/tree disagreement first."
        )
        sys.exit(1)


def _should_stream_pcp(args) -> bool:
    """Decide whether ``process_pcp_format`` runs the streaming pipeline.

    Streaming is the default for the canonical single-file output path.
    Bail out to the legacy in-memory path when:

    - ``--batch-size 0`` — explicit opt-out.
    - ``--split-files`` — multi-file output predates streaming and has a
      different write shape.
    - ``--validate`` — per-batch validation isn't wired yet; today
      ``validate_output_data`` consumes the whole assembled output.

    ``process_pcp_format`` applies a further single-batch fast path
    (``n_families <= --batch-size``) after this check — when the whole
    input fits in one batch, the spool round-trip would cost more than
    the in-memory pipeline saves, so the legacy path runs instead.
    """
    if getattr(args, "batch_size", 0) <= 0:
        return False
    if getattr(args, "split_files", None):
        return False
    if getattr(args, "validate", False):
        return False
    return True


def _should_stream_airr(args) -> bool:
    """Decide whether ``process_airr_format`` runs the streaming pipeline.

    Same fallback conditions as PCP — see :func:`_should_stream_pcp`.
    """
    return _should_stream_pcp(args)


def _process_pcp_streaming(args, pcp_families, newick_trees, minter, input_files):
    """Streaming PCP pipeline: bound peak memory via per-batch spool + stitch.

    Mirrors the assembly the legacy ``process_pcp_to_olmsted`` does (dataset
    mint, sample dedup, tree-csv hoist, ``field_metadata`` generation) but
    in a streaming shape:

    1. Build the dataset stub up front; the iterator appends to
       ``dataset["samples"]`` as it discovers sample IDs.
    2. For each batch from :func:`iter_pcp_clone_groups`, apply any
       user-declared encoded-mutation unpacking, observe the batch in the
       :class:`BatchAccumulator` (which folds evidence, ID-uniqueness
       sets, and tree-level variance), and spool the un-hoisted clones
       and trees.
    3. Apply :func:`apply_dataset_hoist` against the spool to enforce
       the dataset-scope ``_hoist_clone_invariant_extras`` decision.
    4. Finalize ``field_metadata`` and ``processing_info`` from the
       accumulator; stream the consolidated JSON to ``args.output``.
    """
    vprint.status(f"Streaming with --batch-size={args.batch_size}")

    dataset_id = minter.mint("dataset")
    dataset_ident = minter.mint("dataset")

    clones_grouped = _group_pcp_families_by_clone(pcp_families)
    dataset = {
        "ident": dataset_ident,
        "dataset_id": dataset_id,
        "schema_version": SCHEMA_VERSION,
        "build": {"commit": "pcp-import", "time": ""},
        "subjects": [],
        "samples": [],
        "seeds": [],
        "clone_count": len(clones_grouped),
        "subjects_count": 0,
        "timepoints_count": 0,
    }
    if getattr(args, "name", None):
        dataset["name"] = args.name

    tree_config = TreeProcessingConfig(
        compute_metrics=getattr(args, "compute_metrics", False),
        lbi_tau=getattr(args, "lbi_tau", 0.0125),
        standardize_names=getattr(args, "standardize_names", False),
        warn_disagreements=args.warnings,
    )

    custom_fields = getattr(args, "custom_fields", None)
    allow_dup = getattr(args, "allow_duplicate_ids", False)
    accumulator = BatchAccumulator(allow_duplicate_ids=allow_dup)
    accumulator.register_dataset(dataset_id)

    mutations_path = getattr(args, "mutations", None)
    merge_ctx, total_csv_rows = _begin_mutations_merge(args, mutations_path)
    only_listed = getattr(args, "mutations_listed_only", False)

    with BatchSpooler() as spooler:
        try:
            for batch_clones, batch_trees in iter_pcp_clone_groups(
                pcp_families,
                newick_trees,
                minter,
                dataset_id,
                dataset["samples"],
                tree_config,
                batch_size=args.batch_size,
            ):
                if custom_fields:
                    unpack_encoded_mutations(batch_trees, custom_fields)
                if merge_ctx is not None:
                    apply_mutations_to_trees(
                        merge_ctx, batch_trees, only_listed=only_listed
                    )
                accumulator.observe_batch(dataset_id, batch_clones, batch_trees)
                spooler.write_batch(dataset_id, batch_clones, batch_trees)
        except DuplicateIdError as e:
            vprint.error(f"Error: {e}")
            sys.exit(1)

        _finalize_mutations_merge(args, merge_ctx, total_csv_rows)

        for warning in accumulator.duplicate_warnings:
            vprint.error(f"Warning: {warning}")

        apply_dataset_hoist(
            spooler, dataset_id, accumulator.tree_level_keys(dataset_id)
        )

        dataset["field_metadata"] = accumulator.finalize_field_metadata(
            dataset_id, custom_fields
        )

        # Streaming-side id-uniqueness check. The accumulator has already
        # enforced clone_id (within dataset) and tree_id (within clone)
        # via observe_batch; check_output_id_uniqueness catches the
        # remaining scopes — dataset_id across datasets[], sample_id
        # within dataset.samples[], subject_id within dataset.subjects[].
        try:
            check_output_id_uniqueness(
                [dataset], {dataset_id: []}, allow_duplicates=allow_dup
            )
        except ValueError as e:
            vprint.error(f"Error: {e}")
            sys.exit(1)

        # Build the metadata wrapper via create_consolidated_data so the
        # processing_options / source_format / generated_by sections stay
        # in lockstep with the legacy path. The empty trees/clones args
        # produce zeroed totals which we then overwrite with the
        # accumulator's running counts.
        wrapper = create_consolidated_data(
            [dataset],
            {dataset_id: []},
            [],
            input_files,
            FORMAT_PCP,
            args,
        )
        wrapper["metadata"]["processing_info"] = accumulator.finalize_totals()

        output_dir = os.path.dirname(args.output) or "."
        os.makedirs(output_dir, exist_ok=True)
        vprint.status(f"Writing Olmsted JSON output to {args.output}")
        write_olmsted_json_streaming(
            wrapper["metadata"],
            wrapper["datasets"],
            spooler,
            args.output,
            json_format=getattr(args, "json_format", "pretty"),
        )


def _process_airr_streaming(args, airr_args):
    """Streaming AIRR pipeline: per-file dataset stub + per-batch spool.

    Iterates each AIRR input file the same way the legacy
    ``process_airr_format`` does — reading the dataset, optionally
    filtering invalid clones, validating the dataset envelope — but
    consumes the per-file clones via :func:`iter_airr_clones`, so peak
    memory tracks one batch's worth of parsed trees rather than every
    clone in every input file.  A shared :class:`BatchAccumulator` and
    :class:`BatchSpooler` aggregate across files; a single
    :class:`MergeContext` carries ``--mutations`` state across batches.
    """
    vprint.status(f"Streaming with --batch-size={args.batch_size}")

    custom_fields = getattr(airr_args, "custom_fields", None)
    allow_dup = getattr(airr_args, "allow_duplicate_ids", False)
    accumulator = BatchAccumulator(allow_duplicate_ids=allow_dup)

    mutations_path = getattr(args, "mutations", None)
    merge_ctx, total_csv_rows = _begin_mutations_merge(args, mutations_path)
    only_listed = getattr(args, "mutations_listed_only", False)

    dataset_headers: List[Dict[str, Any]] = []

    with BatchSpooler() as spooler:
        input_files = airr_args.inputs or []
        with tqdm(
            input_files,
            desc="Processing AIRR files",
            unit="file",
            disable=len(input_files) == 1,
        ) as pbar:
            for infile in pbar:
                pbar.set_description(f"Processing {Path(infile).name}")
                if len(input_files) == 1 or airr_args.verbose:
                    vprint.status(f"\nProcessing AIRR file: {infile}")

                try:
                    dataset_in = read_airr_json(infile)

                    if airr_args.remove_invalid_clones:
                        original_count = len(dataset_in.get("clones", []))
                        dataset_in["clones"] = list(
                            filter(
                                jsonschema.Draft4Validator(clone_spec).is_valid,
                                dataset_in["clones"],
                            )
                        )
                        filtered_count = original_count - len(dataset_in["clones"])
                        if filtered_count > 0:
                            pbar.set_postfix({"filtered": filtered_count})

                    errors = validate_dataset(dataset_in, verbose=airr_args.verbose).errors
                    if errors:
                        error_msg = "Dataset validation failed"
                        if airr_args.verbose:
                            vprint.error("Dataset validation failed:")
                            for error in errors:
                                vprint.error(f"  - {error}")
                        else:
                            error_msg += ". Please rerun with `-v` for detailed errors"
                        raise Exception(error_msg)

                    dataset_id = dataset_in["dataset_id"]
                    accumulator.register_dataset(
                        dataset_id, hoist_tree_extras_to_clone=False
                    )
                    accumulator.add_samples(
                        dataset_id, dataset_in.get("samples", []) or []
                    )

                    input_clones = dataset_in.get("clones", []) or []
                    input_clone_count = len(input_clones)
                    subjects_count = len(
                        {
                            cf["subject_id"]
                            for cf in input_clones
                            if cf.get("subject_id")
                        }
                    )
                    timepoints_count = len(
                        {
                            s["timepoint_id"]
                            for s in dataset_in.get("samples", []) or []
                            if s.get("timepoint_id")
                        }
                    )

                    for batch_clones, batch_trees in iter_airr_clones(
                        airr_args,
                        dataset_in,
                        batch_size=args.batch_size,
                    ):
                        if custom_fields:
                            unpack_encoded_mutations(batch_trees, custom_fields)
                        if merge_ctx is not None:
                            apply_mutations_to_trees(
                                merge_ctx, batch_trees, only_listed=only_listed
                            )
                        try:
                            accumulator.observe_batch(
                                dataset_id, batch_clones, batch_trees
                            )
                        except DuplicateIdError as e:
                            vprint.error(f"Error: {e}")
                            sys.exit(1)
                        spooler.write_batch(dataset_id, batch_clones, batch_trees)

                    # Build dataset header (input minus consumed clones).
                    # Mirrors process_dataset's mutations to the dataset dict,
                    # finalized after iter_airr_clones drained dataset_in["clones"].
                    dataset_header = {
                        k: v for k, v in dataset_in.items() if k != "clones"
                    }
                    dataset_header["clone_count"] = input_clone_count
                    dataset_header["subjects_count"] = subjects_count
                    dataset_header["timepoints_count"] = timepoints_count
                    dataset_header["schema_version"] = SCHEMA_VERSION
                    dataset_header = ensure_ident(
                        dataset_header, "dataset", airr_args.minter
                    )
                    dataset_headers.append(dataset_header)

                    if len(input_files) > 1:
                        pbar.set_postfix({"clones": input_clone_count})

                except Exception:
                    vprint.error(f"\nUnable to process AIRR file: {infile}")
                    if airr_args.verbose:
                        exc_info = sys.exc_info()
                        traceback.print_exception(*exc_info)
                    else:
                        vprint.error("Please rerun with `-v` for detailed errors.")
                    sys.exit(1)

        _finalize_mutations_merge(args, merge_ctx, total_csv_rows)

        for warning in accumulator.duplicate_warnings:
            vprint.error(f"Warning: {warning}")

        # AIRR data places tree-level fields on trees natively; the
        # PCP-style data hoist would strip them from trees and emit
        # phantom clone-level entries.  The accumulator is configured
        # with hoist_tree_extras_to_clone=False for this case, so the
        # finalize step also leaves clone-level metadata alone.
        for dataset_header in dataset_headers:
            dataset_id = dataset_header["dataset_id"]
            dataset_header["field_metadata"] = accumulator.finalize_field_metadata(
                dataset_id, custom_fields
            )

        # Streaming-side id-uniqueness check. The accumulator has
        # already enforced clone_id / tree_id during observe_batch;
        # check_output_id_uniqueness covers dataset_id, sample_id,
        # subject_id — important for AIRR where input files can carry
        # duplicates the legacy path would have rejected before write.
        empty_clones = {h["dataset_id"]: [] for h in dataset_headers}
        try:
            check_output_id_uniqueness(
                dataset_headers, empty_clones, allow_duplicates=allow_dup
            )
        except ValueError as e:
            vprint.error(f"Error: {e}")
            sys.exit(1)
        wrapper = create_consolidated_data(
            dataset_headers,
            empty_clones,
            [],
            args.inputs,
            FORMAT_AIRR,
            args,
        )
        wrapper["metadata"]["processing_info"] = accumulator.finalize_totals()

        output_dir = os.path.dirname(args.output) or "."
        os.makedirs(output_dir, exist_ok=True)
        vprint.status(f"Writing Olmsted JSON output to {args.output}")
        write_olmsted_json_streaming(
            wrapper["metadata"],
            wrapper["datasets"],
            spooler,
            args.output,
            json_format=getattr(args, "json_format", "pretty"),
        )


def process_pcp_format(args):
    """
    Process PCP format files using the existing PCP processor.

    Args:
        args: Parsed command line arguments
    """

    vprint.status("Processing PCP format...")

    # Print command arguments at verbosity level 2
    vprint.verbose("=== Command Arguments ===")
    vprint.verbose(f"  Input PCP file: {args.inputs[0]}")
    if hasattr(args, "tree") and args.tree:
        vprint.verbose(f"  Input trees file: {args.tree}")
    if args.output:
        vprint.verbose(f"  Output file: {args.output}")
    if args.split_files:
        vprint.verbose(f"  Output directory: {args.split_files}")
    if hasattr(args, "name") and args.name:
        vprint.verbose(f"  Dataset name: {args.name}")
    vprint.verbose(f"  Verbosity level: {args.verbose}")
    vprint.verbose(f"  Validation: {args.validate}")
    if args.validate:
        vprint.verbose(f"  Strict validation: {args.strict_validation}")
    if hasattr(args, "seed") and args.seed is not None:
        vprint.verbose(f"  Random seed: {args.seed}")
    vprint.verbose(f"  Show disagreement warnings: {args.warnings}")
    vprint.verbose(f"  Compute metrics: {getattr(args, 'compute_metrics', False)}")
    if getattr(args, "compute_metrics", False):
        vprint.verbose(f"    LBI tau: {getattr(args, 'lbi_tau', 0.0125)}")
    vprint.verbose(f"  Standardize names: {getattr(args, 'standardize_names', False)}")
    vprint.verbose("=" * 25)
    vprint.verbose("")

    # Set up identifier minter (deterministic if seed provided)
    minter = IdentMinter(seed=getattr(args, "seed", None))

    try:
        # Get PCP file from inputs
        pcp_file = args.inputs[0]

        # Get trees file from --tree argument if provided
        trees_file = args.tree if hasattr(args, "tree") else None

        vprint.status(f"Processing PCP CSV: {pcp_file}")
        if hasattr(args, "seed") and args.seed is not None:
            vprint.status(f"Using deterministic UUIDs with seed: {args.seed}")

        column_overrides = {
            "sample_override": getattr(args, "sample_col", None),
            "family_override": getattr(args, "family_col", None),
            "tree_override": getattr(args, "tree_col", None),
        }

        # Parse PCP families with progress bar
        vprint.status("Parsing PCP CSV...")
        pcp_families = parse_pcp_csv(pcp_file, **column_overrides)
        vprint.status(f"Found {len(pcp_families)} families")

        # Parse Newick trees if provided with progress bar
        newick_trees = None
        if trees_file:
            vprint.status(f"Processing Newick trees: {trees_file}")
            newick_trees = parse_newick_csv(trees_file, **column_overrides)
            vprint.status(f"Found {len(newick_trees)} trees")

        # Streaming pipeline (#26): bound peak memory by spooling each
        # batch's clones/trees to disk and stream-stitching the final
        # consolidated JSON. Falls back to the legacy in-memory path
        # when the caller wants split-file output or output validation
        # (phase 5 will wire per-batch validation).
        if _should_stream_pcp(args):
            # Single-batch fast path: when the whole input fits in one
            # batch, the spool round-trip costs more than the legacy
            # in-memory path saves. Routing through ``process_pcp_to_olmsted``
            # avoids the temp-file write+read for small inputs while
            # producing identical output.
            n_families = len(_group_pcp_families_by_clone(pcp_families))
            if n_families > args.batch_size:
                _process_pcp_streaming(
                    args, pcp_families, newick_trees, minter, args.inputs
                )
                vprint.status("Processing complete!")
                return
            vprint.verbose(
                f"  Single-batch fast path: {n_families} families "
                f"<= --batch-size {args.batch_size}; using in-memory pipeline."
            )

        # Convert to Olmsted format with progress bar
        vprint.status("Converting to Olmsted format...")
        datasets, clones_dict, trees = process_pcp_to_olmsted(
            pcp_families,
            newick_trees,
            minter,
            args.warnings,
            compute_metrics=getattr(args, "compute_metrics", False),
            lbi_tau=getattr(args, "lbi_tau", 0.0125),
            standardize_names=getattr(args, "standardize_names", False),
            name=getattr(args, "name", None),
            verbosity=args.verbose,
            custom_fields=getattr(args, "custom_fields", None),
        )

        # Merge mutations CSV if --mutations was specified
        mutations_path = getattr(args, "mutations", None)
        if mutations_path:
            try:
                apply_mutations_csv(
                    mutations_path,
                    trees,
                    use_depth=getattr(args, "mutations_use_depth", False),
                    allow_mismatch=getattr(args, "mutations_allow_mismatch", False),
                    only_listed=getattr(args, "mutations_listed_only", False),
                )
            except ValueError as e:
                vprint.error(f"Error: {e}")
                sys.exit(1)
            retag_datasets_field_metadata(
                datasets,
                clones_dict,
                trees,
                custom_fields=getattr(args, "custom_fields", None),
            )

        # Enforce *_id uniqueness before writing. --allow-duplicate-ids
        # downgrades collisions to a warning.
        try:
            check_output_id_uniqueness(
                datasets,
                clones_dict,
                allow_duplicates=getattr(args, "allow_duplicate_ids", False),
            )
        except ValueError as e:
            vprint.error(f"Error: {e}")
            sys.exit(1)

        # Validate data if requested
        if args.validate:
            if not validate_output_data(datasets, clones_dict, trees, args):
                if args.strict_validation:
                    vprint.error(
                        "\nExiting due to validation errors (--strict-validation enabled)"
                    )
                    sys.exit(1)

        # Write output
        if args.split_files:
            # Multi-file output to specified directory
            output_dir = args.split_files
            os.makedirs(output_dir, exist_ok=True)
            vprint.status(f"Writing output to {output_dir}")
            write_out(datasets, output_dir, "datasets.json", args)
            for dataset_id, clones in clones_dict.items():
                write_out(clones, output_dir, f"clones.{dataset_id}.json", args)
            for tree in trees:
                write_out(tree, output_dir, f"tree.{tree['ident']}.json", args)
        else:
            # Olmsted JSON output (default)
            consolidated_data = create_consolidated_data(
                datasets, clones_dict, trees, args.inputs, FORMAT_PCP, args
            )
            # Ensure output directory exists
            output_dir = os.path.dirname(args.output) or "."
            output_file = os.path.basename(args.output)
            os.makedirs(output_dir, exist_ok=True)
            vprint.status(f"Writing Olmsted JSON output to {args.output}")
            write_out(consolidated_data, output_dir, output_file, args)

        vprint.status("Processing complete!")

    except Exception as e:
        vprint.error(f"Error processing PCP format: {e}")
        if args.verbose:
            traceback.print_exc()
        sys.exit(1)


def build_parser():
    """Build the argument parser for the unified processor."""
    parser = argparse.ArgumentParser(
        description="Convert input data (AIRR JSON or PCP CSV) to Olmsted JSON format",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Process with a config file (recommended)
    olmsted process -c config.yaml

    # PCP format with trees
    olmsted process -i data.csv -t trees.csv -o output.json

    # AIRR format (auto-detected)
    olmsted process -i data.json -o output.json

    # With metrics and validation
    olmsted process -i data.csv -t trees.csv -o output.json --compute-metrics --validate
        """,
    )

    # --- Core arguments ---
    parser.add_argument(
        "-i",
        "--input",
        "--inputs",
        dest="inputs",
        nargs="+",
        help="Input file(s). AIRR: JSON file(s). PCP: CSV file",
    )
    parser.add_argument(
        "-t",
        "--tree",
        help="Companion tree CSV file (PCP format)",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Output Olmsted JSON file path",
    )
    parser.add_argument(
        "-c",
        "--config",
        help="YAML configuration file (CLI arguments override config values)",
    )
    parser.add_argument(
        "-f",
        "--format",
        choices=[FORMAT_AIRR, FORMAT_PCP, FORMAT_AUTO],
        default=FORMAT_AUTO,
        help="Input format (default: auto-detect)",
    )

    # --- Column overrides (PCP CSV inputs) ---
    column_group = parser.add_argument_group(
        "PCP CSV column overrides",
        "Override the auto-detected column names. By default the parser\n"
        "accepts any of {sample, sample_id, sample_name},\n"
        "{family, family_id, family_name}, {tree, tree_id, tree_name},\n"
        "preferring '_id' > bare > '_name' when multiple are present.",
    )
    column_group.add_argument(
        "--sample-col",
        dest="sample_col",
        help="Column name supplying the sample identifier (PCP / trees CSV).",
    )
    column_group.add_argument(
        "--family-col",
        dest="family_col",
        help="Column name supplying the family identifier (PCP / trees CSV).",
    )
    column_group.add_argument(
        "--tree-col",
        dest="tree_col",
        help="Column name supplying the per-tree identifier within a family. "
        "Optional — when absent, every family has at most one tree.",
    )

    # --- Mutations CSV options ---
    mutations_group = parser.add_argument_group(
        "mutations CSV options",
        # Hand-wrapped because the parser uses RawDescriptionHelpFormatter
        # (to preserve the epilog's example block), which also disables
        # wrapping on argument-group descriptions.
        "Pass --mutations FILE to merge mutation-level data into tree nodes\n"
        "after processing. Equivalent to running `olmsted merge` against the\n"
        "output. The remaining flags in this group modify that merge.",
    )
    mutations_group.add_argument(
        "--mutations",
        help="Mutations CSV file (columns: family, site, parent_aa, child_aa, ...). "
        "Mutation-level scores are merged into tree nodes after processing.",
    )
    mutations_group.add_argument(
        "--mutations-use-depth",
        action="store_true",
        help="Use an optional 'depth' column in the mutations CSV to extend the "
        "match key to (site, parent_aa, child_aa, depth). Ignored when the "
        "CSV has a node-name column. Opt-in because depth arithmetic depends "
        "on the upstream rooting convention.",
    )
    mutations_group.add_argument(
        "--mutations-allow-mismatch",
        action="store_true",
        help="Proceed past integrity mismatches between the mutations CSV "
        "and the tree's derived mutations. By default processing fails on "
        "any such mismatch. Mismatched rows are always skipped; the flag "
        "only controls whether the command exits non-zero afterwards.",
    )
    mutations_group.add_argument(
        "--mutations-listed-only",
        action="store_true",
        help="Treat the mutations CSV as authoritative: on trees whose "
        "clone_id matches a family in the CSV, drop any derived "
        "mutations that don't appear in the CSV. Trees whose family is "
        "absent from the CSV pass through untouched.",
    )

    # --- Dataset metadata ---
    parser.add_argument(
        "-n",
        "--name",
        help="Dataset name (stored in output metadata)",
    )
    parser.add_argument(
        "--description",
        help="Dataset description (stored in output metadata)",
    )

    # --- Processing options ---
    parser.add_argument(
        "--compute-metrics",
        action="store_true",
        help="Compute LBI, LBR, affinity, scaled_affinity for all tree nodes",
    )
    parser.add_argument(
        "--lbi-tau",
        type=float,
        default=0.0125,
        help="Time scale parameter for LBI calculation (default: 0.0125)",
    )
    parser.add_argument(
        "-r",
        "--root",
        nargs="?",
        const="naive",
        default=None,
        metavar="NAME",
        help="Root trees at the naive/germline node. Optionally specify node name (default: 'naive'). AIRR only.",
    )
    parser.add_argument(
        "--standardize-names",
        action="store_true",
        help="Rename nodes: naive (root), Node1, Node2, ..., Leaf1, Leaf2, ...",
    )

    # --- Output options ---
    parser.add_argument(
        "--seed",
        type=int,
        help="Random seed for deterministic UUID generation",
    )
    parser.add_argument(
        "--json-format",
        choices=["pretty", "compact", "gzip"],
        default="pretty",
        help="JSON output format (default: pretty)",
    )
    parser.add_argument(
        "--split-files",
        metavar="DIR",
        help="Output split files instead of single Olmsted JSON (legacy)",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=50,
        metavar="N",
        help="Number of clonal families to process per batch (PCP format only). "
        "Bounds peak memory by spooling each batch's clones/trees to disk and "
        "stream-stitching the final consolidated JSON. Default: 50. Pass 0 to "
        "disable batching and use the legacy one-shot path.",
    )

    # --- Validation ---
    parser.add_argument(
        "--validate",
        action="store_true",
        help="Validate output against schemas before writing",
    )
    parser.add_argument(
        "--strict-validation",
        action="store_true",
        help="Exit with error if validation fails",
    )
    parser.add_argument(
        "--allow-duplicate-ids",
        action="store_true",
        help="Downgrade duplicate-*_id errors to warnings and pass the data "
        "through unchanged. By default, processing fails when dataset_id, "
        "clone_id, tree_id, sample_id, or subject_id collide within their "
        "natural uniqueness scope — downstream consumers (webapp Redux, "
        "Dexie bulkPut) would silently overwrite on collision.",
    )
    parser.add_argument(
        "-w",
        "--warnings",
        action="store_true",
        help="Show warnings when tree and PCP data disagree",
    )

    # --- Verbosity ---
    add_verbosity_args(parser)

    # --- Advanced ---
    parser.add_argument(
        "--capture-all",
        action="store_true",
        help="Capture all data fields from input (extra CSV columns, unknown JSON fields)",
    )

    return parser


# Mapping from YAML config keys to argparse dest names
_CONFIG_KEY_MAP = {
    "inputs": "inputs",
    "output": "output",
    "format": "format",
    "split_files": "split_files",
    "json_format": "json_format",
    "name": "name",
    "description": "description",
    "verbose": "verbose",
    "quiet": "quiet",
    "validate": "validate",
    "strict_validation": "strict_validation",
    "seed": "seed",
    "warnings": "warnings",
    "tree": "tree",
    "mutations": "mutations",
    "root": "root",
    "compute_metrics": "compute_metrics",
    "lbi_tau": "lbi_tau",
    "standardize_names": "standardize_names",
    "capture_all": "capture_all",
    "sample_col": "sample_col",
    "family_col": "family_col",
    "tree_col": "tree_col",
    "batch_size": "batch_size",
}

# Valid config keys (including custom_fields which is handled separately)
# Tag-specific keys are also accepted (input, mode) so configs work for both commands.
_VALID_CONFIG_KEYS = set(_CONFIG_KEY_MAP.keys()) | {"custom_fields", "input", "mode"}


def load_config(config_path):
    """
    Load and validate a YAML configuration file.

    Args:
        config_path: Path to the YAML config file.

    Returns:
        Tuple of (config_dict, custom_fields_list).
        config_dict maps argparse dest names to values.
        custom_fields_list is a list of custom field declaration dicts.

    Raises:
        SystemExit: If config file cannot be loaded.
    """
    config_path = Path(config_path)
    if not config_path.exists():
        vprint.error(f"Error: Config file not found: {config_path}")
        sys.exit(1)

    try:
        raw_config = read_yaml_config(config_path)
    except yaml.YAMLError as e:
        vprint.error(f"Error: Invalid YAML in config file: {e}")
        sys.exit(1)

    if not raw_config or not isinstance(raw_config, dict):
        return {}, []

    config_dir = config_path.parent

    # Warn about unrecognized keys
    for key in raw_config:
        if key not in _VALID_CONFIG_KEYS:
            vprint.error(f"Warning: Unrecognized config key '{key}' (ignored)")

    # Map config keys to argparse dest names
    config_dict = {}
    for config_key, arg_dest in _CONFIG_KEY_MAP.items():
        if config_key in raw_config:
            value = raw_config[config_key]
            # Resolve file paths relative to config file directory
            if config_key in ("inputs", "tree", "output", "split_files"):
                value = _resolve_paths(value, config_dir)
            config_dict[arg_dest] = value

    # Tag-specific keys (not in _CONFIG_KEY_MAP)
    if "input" in raw_config:
        config_dict["input"] = _resolve_paths(raw_config["input"], config_dir)
    if "mode" in raw_config:
        config_dict["mode"] = raw_config["mode"]

    # Parse custom_fields
    custom_fields = []
    if "custom_fields" in raw_config:
        raw_fields = raw_config["custom_fields"]
        if isinstance(raw_fields, list):
            for i, entry in enumerate(raw_fields):
                if not isinstance(entry, dict):
                    vprint.error(f"Warning: custom_fields[{i}] is not a dict (ignored)")
                    continue
                # Skip entries only need name and level
                is_skip = entry.get("skip", False)
                if is_skip:
                    required_keys = {"name", "level"}
                else:
                    required_keys = {"name", "level", "type", "label"}
                missing = required_keys - set(entry.keys())
                if missing:
                    vprint.error(
                        f"Warning: custom_fields[{i}] missing required keys: {missing} (ignored)"
                    )
                    continue
                if entry.get("level") not in FIELD_LEVELS:
                    vprint.error(
                        f"Warning: custom_fields[{i}] has invalid level '{entry['level']}' (ignored)"
                    )
                    continue
                # Normalize level alias (family → clone)
                entry["level"] = normalize_level(entry["level"])
                if not is_skip and entry.get("type") not in FIELD_TYPES:
                    vprint.error(
                        f"Warning: custom_fields[{i}] has invalid type '{entry['type']}' (ignored)"
                    )
                    continue
                # Validate display mode if specified
                display = entry.get("display")
                if display and display not in DISPLAY_MODES:
                    vprint.error(
                        f"Warning: custom_fields[{i}] has invalid display '{display}' (ignored)"
                    )
                    continue
                # Validate encoding if specified (mutation-level only)
                encoding = entry.get("encoding")
                if encoding:
                    if encoding not in MUTATION_ENCODINGS:
                        vprint.error(
                            f"Warning: custom_fields[{i}] has invalid encoding '{encoding}' (ignored)"
                        )
                        continue
                    if entry["level"] != "mutation":
                        vprint.error(
                            f"Warning: custom_fields[{i}] has encoding but level is '{entry['level']}', not 'mutation' (ignored)"
                        )
                        continue
                    if encoding == "records" and "source" not in entry:
                        vprint.error(
                            f"Warning: custom_fields[{i}] encoding 'records' requires 'source' key (ignored)"
                        )
                        continue
                custom_fields.append(entry)

    return config_dict, custom_fields


def _resolve_paths(value, config_dir):
    """Resolve file paths relative to config file directory."""
    if isinstance(value, list):
        return [_resolve_paths(v, config_dir) for v in value]
    if isinstance(value, str):
        p = Path(value)
        if not p.is_absolute():
            resolved = config_dir / p
            return str(resolved)
    return value


def get_args():
    """
    Parse command line arguments with optional YAML config file support.

    Precedence: argparse defaults < YAML config < explicit CLI args.
    """
    parser = build_parser()

    # First pass: parse with defaults suppressed to find explicit CLI args
    # We need to know which args the user explicitly provided on the command line
    # vs which are just argparse defaults, so config values fill in the gaps.
    explicit_parser = build_parser()
    explicit_parser.set_defaults(**{dest: None for dest in _CONFIG_KEY_MAP.values()})
    explicit_parser.set_defaults(config=None, quiet=None)
    explicit_args, _ = explicit_parser.parse_known_args()

    # Second pass: normal parse with defaults
    args = parser.parse_args()

    # Load config if specified (from either CLI or first-pass)
    config_path = explicit_args.config or args.config
    custom_fields = None

    if config_path:
        config_dict, custom_fields = load_config(config_path)

        # Apply config values where CLI didn't explicitly set them
        for dest, config_value in config_dict.items():
            explicit_value = getattr(explicit_args, dest, None)
            if explicit_value is None:
                setattr(args, dest, config_value)

    # Attach custom_fields to args for downstream use
    args.custom_fields = custom_fields

    # Inputs is required (either from CLI or config)
    if not args.inputs:
        parser.error(
            "the following arguments are required: -i/--inputs (or provide in config)"
        )

    # Mutation flags only make sense with --mutations
    if (
        args.mutations_use_depth
        or args.mutations_allow_mismatch
        or args.mutations_listed_only
    ) and not args.mutations:
        parser.error(
            "--mutations-use-depth / --mutations-allow-mismatch / "
            "--mutations-listed-only require --mutations"
        )

    return args


def main():
    """Main entry point for the unified processor."""
    args = get_args()

    # Handle quiet mode
    resolve_verbosity(args)

    # Set global verbosity
    set_verbosity(args.verbose)

    # Validate output arguments
    if not args.output and not args.split_files:
        vprint.error("Error: Either -o/--output or --split-files must be specified")
        sys.exit(1)

    if args.output and args.split_files:
        vprint.error("Error: Cannot specify both -o/--output and --split-files")
        sys.exit(1)

    # Validate inputs
    if not args.inputs:
        vprint.error("Error: No input files specified")
        sys.exit(1)

    # Check that input files exist
    for input_file in args.inputs:
        if not os.path.exists(input_file):
            vprint.error(f"Error: Input file does not exist: {input_file}")
            sys.exit(1)

    # Determine format
    if args.format == FORMAT_AUTO:
        detected_format = detect_file_format(args.inputs[0])
        if detected_format == FORMAT_UNKNOWN:
            vprint.error(f"Error: Could not auto-detect format for {args.inputs[0]}")
            vprint.error("Please specify format with -f/--format option")
            sys.exit(1)
        format_to_use = detected_format
        vprint.status(f"Auto-detected format: {format_to_use}")
    else:
        format_to_use = args.format
        vprint.status(f"Using specified format: {format_to_use}")

    # Validate format matches file content
    if format_to_use == FORMAT_AIRR:
        for input_file in args.inputs:
            if not validate_airr_file(input_file):
                vprint.status(f"Warning: {input_file} may not be valid AIRR format")
    elif format_to_use == FORMAT_PCP:
        if not validate_pcp_file(args.inputs[0]):
            vprint.status(f"Warning: {args.inputs[0]} may not be valid PCP format")

    # Process based on format
    try:
        if format_to_use == FORMAT_AIRR:
            process_airr_format(args)
        elif format_to_use == FORMAT_PCP:
            process_pcp_format(args)
        elif format_to_use == FORMAT_OLMSTED:
            vprint.error(
                "Error: Input is already in Olmsted JSON format. "
                "Use 'olmsted tag' to add field_metadata to existing Olmsted files."
            )
            sys.exit(1)
        else:
            vprint.error(f"Error: Unsupported format: {format_to_use}")
            sys.exit(1)

        vprint.status(f"\nSuccessfully processed {format_to_use.upper()} format data")

    except KeyboardInterrupt:
        vprint.error("\nProcessing interrupted by user")
        sys.exit(1)
    except Exception as e:
        vprint.error(f"Error during processing: {e}")
        if args.verbose:
            traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
