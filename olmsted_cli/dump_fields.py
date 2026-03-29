"""
Dump-fields command: extract all available fields from data into a YAML
configuration that the user can edit.

Accepts Olmsted JSON, AIRR JSON, or PCP CSV input. For PCP and AIRR, the
data is processed through the standard pipeline first (in memory) so that
the dumped fields reflect what would actually appear in the output.

Usage:
    olmsted dump-fields -i data.json -o config.yaml
    olmsted dump-fields -i pcp.csv -t trees.csv -o config.yaml
    olmsted dump-fields -i data.json                   # prints to stdout
"""

import argparse
import json
import sys
from pathlib import Path

from .process_data import detect_file_format
from .field_metadata import (
    EXCLUDED_CLONE_FIELDS,
    EXCLUDED_MUTATION_FIELDS,
    EXCLUDED_NODE_FIELDS,
    KNOWN_BRANCH_FIELDS,
    KNOWN_CLONE_FIELDS,
    KNOWN_MUTATION_FIELDS,
    KNOWN_NODE_FIELDS,
    humanize_label,
    infer_field_type,
    _collect_mutations,
    _collect_nodes,
    _collect_keys,
    _sample_values,
)


def get_args():
    """Parse command line arguments for the dump-fields command."""
    parser = argparse.ArgumentParser(
        description="Extract all available fields from data into a YAML config",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # From Olmsted JSON
    olmsted dump-fields -i data.json -o config.yaml

    # From PCP CSV (with optional trees)
    olmsted dump-fields -i pcp.csv -t trees.csv -o config.yaml

    # From AIRR JSON
    olmsted dump-fields -i airr_data.json -o config.yaml

    # Print to stdout for inspection
    olmsted dump-fields -i data.json

    # Then edit the config and use it
    olmsted enrich -i data.json -o enriched.json -c config.yaml
    olmsted process -i pcp.csv -t trees.csv -o output.json -c config.yaml
        """,
    )

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input file: Olmsted JSON, AIRR JSON, or PCP CSV",
    )
    parser.add_argument(
        "-t",
        "--tree",
        "--trees",
        help="Trees CSV file (PCP format only, optional)",
    )
    parser.add_argument(
        "-f",
        "--format",
        choices=["airr", "pcp", "olmsted", "auto"],
        default="auto",
        help="Input format (default: auto-detect)",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Output YAML file (default: stdout)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for UUID generation during processing (default: 42)",
    )
    parser.add_argument(
        "--compute-metrics",
        action="store_true",
        help="Compute metrics (LBI, LBR, etc.) when processing PCP data",
    )

    return parser.parse_args()


def _load_olmsted(input_path):
    """Load clones and trees from Olmsted JSON."""
    with open(input_path) as f:
        data = json.load(f)

    all_clones = []
    for ds_clones in data.get("clones", {}).values():
        all_clones.extend(ds_clones)

    all_trees = data.get("trees", [])
    return all_clones, all_trees


def _load_pcp(input_path, trees_path, seed, compute_metrics):
    """Process PCP data and return clones and trees."""
    from .process_pcp_data import (
        parse_newick_csv,
        parse_pcp_csv,
        process_pcp_to_olmsted,
        deterministic_uuid,
    )

    pcp_families = parse_pcp_csv(str(input_path))

    newick_trees = None
    if trees_path:
        newick_trees = parse_newick_csv(str(trees_path))

    counter = [0]

    def get_uuid(prefix=""):
        result = deterministic_uuid(seed, counter[0])
        counter[0] += 1
        return f"{prefix}{result}"

    datasets, clones_dict, trees = process_pcp_to_olmsted(
        pcp_families,
        newick_trees,
        uuid_generator=get_uuid,
        compute_metrics=compute_metrics,
        verbosity=0,
    )

    all_clones = []
    for ds_clones in clones_dict.values():
        all_clones.extend(ds_clones)

    return all_clones, trees


def _load_airr(input_path):
    """Process AIRR data and return clones and trees."""
    from argparse import Namespace

    from .process_airr_data import process_dataset

    with open(input_path) as f:
        data = json.load(f)

    # Handle single-dataset or multi-dataset AIRR
    if isinstance(data, list):
        datasets_raw = data
    else:
        datasets_raw = [data]

    args = Namespace(
        naive_name="naive",
        root_trees=False,
        verbose=0,
        custom_fields=None,
    )

    all_clones = []
    all_trees = []
    for dataset in datasets_raw:
        clones_dict = {}
        trees = []
        process_dataset(args, dataset, clones_dict, trees)
        for ds_clones in clones_dict.values():
            all_clones.extend(ds_clones)
        all_trees.extend(trees)

    return all_clones, all_trees


def _field_summary(dicts, field, known_registry):
    """Build a summary dict for a single field: type, label, sample values."""
    if field in known_registry:
        entry = dict(known_registry[field])
    else:
        values = _sample_values(dicts, field, max_samples=50)
        entry = {
            "type": infer_field_type(values),
            "label": humanize_label(field),
        }
    return entry


def _format_field_block(name, level, entry, sample_values=None):
    """Format a single custom_fields YAML entry as a string."""
    lines = []
    lines.append(f"  - name: {name}")
    lines.append(f"    level: {level}")
    lines.append(f"    type: {entry['type']}")
    lines.append(f"    label: \"{entry['label']}\"")
    if sample_values:
        preview = ", ".join(str(v) for v in sample_values[:5])
        if len(sample_values) > 5:
            preview += ", ..."
        lines.append(f"    # sample values: {preview}")
    return "\n".join(lines)


def _build_yaml(input_name, detected_format, all_clones, all_trees):
    """Build the YAML config string from clones and trees."""
    all_nodes = _collect_nodes(all_trees, max_nodes=500)
    all_mutations = _collect_mutations(all_trees, max_mutations=500)

    lines = []
    lines.append("# =============================================================================")
    lines.append(f"# Field configuration generated from: {input_name}")
    lines.append(f"# Detected format: {detected_format}")
    lines.append("# =============================================================================")
    lines.append("#")
    lines.append("# This file lists all fields discovered in your data. Edit it to:")
    lines.append("#   - Remove fields you don't want in the web app dropdowns")
    lines.append("#   - Change type (continuous, categorical, tooltip)")
    lines.append("#   - Customize display labels")
    lines.append("#")
    if detected_format == "olmsted":
        lines.append("# Then use with:  olmsted enrich -i data.json -o enriched.json -c this_file.yaml")
    else:
        lines.append("# Then use with:  olmsted process -c this_file.yaml")
    lines.append("")
    lines.append("custom_fields:")

    # --- Clone level ---
    clone_keys = _collect_keys(all_clones) - EXCLUDED_CLONE_FIELDS
    has_locus = any(
        isinstance(c.get("sample"), dict) and c["sample"].get("locus") is not None
        for c in all_clones[:20]
    )
    if has_locus:
        clone_keys.add("locus")

    if clone_keys:
        lines.append("")
        lines.append("  # --- Clone level (scatterplot axes, color, facet) ---")
        for field in sorted(clone_keys):
            entry = _field_summary(all_clones, field, KNOWN_CLONE_FIELDS)
            if field == "locus":
                samples = list({
                    c["sample"]["locus"]
                    for c in all_clones[:50]
                    if isinstance(c.get("sample"), dict) and c["sample"].get("locus")
                })
            else:
                samples = _sample_values(all_clones, field, max_samples=6)
            lines.append(_format_field_block(field, "clone", entry, samples))

    # --- Node level ---
    node_keys = _collect_keys(all_nodes) - EXCLUDED_NODE_FIELDS
    node_keys -= set(KNOWN_BRANCH_FIELDS.keys())

    if node_keys:
        lines.append("")
        lines.append("  # --- Node level (tree node properties, tooltips) ---")
        for field in sorted(node_keys):
            entry = _field_summary(all_nodes, field, KNOWN_NODE_FIELDS)
            samples = _sample_values(all_nodes, field, max_samples=6)
            lines.append(_format_field_block(field, "node", entry, samples))

    # --- Branch level ---
    branch_keys = _collect_keys(all_nodes) & set(KNOWN_BRANCH_FIELDS.keys())
    if branch_keys:
        lines.append("")
        lines.append("  # --- Branch level (tree branch coloring, width) ---")
        for field in sorted(branch_keys):
            entry = dict(KNOWN_BRANCH_FIELDS[field])
            samples = _sample_values(all_nodes, field, max_samples=6)
            lines.append(_format_field_block(field, "branch", entry, samples))

    # --- Mutation level ---
    mutation_keys = _collect_keys(all_mutations) - EXCLUDED_MUTATION_FIELDS
    if mutation_keys:
        lines.append("")
        lines.append("  # --- Mutation level (alignment coloring) ---")
        for field in sorted(mutation_keys):
            entry = _field_summary(all_mutations, field, KNOWN_MUTATION_FIELDS)
            samples = _sample_values(all_mutations, field, max_samples=6)
            lines.append(_format_field_block(field, "mutation", entry, samples))

    lines.append("")
    return "\n".join(lines)


def main():
    """Main entry point for the dump-fields command."""
    args = get_args()

    input_path = Path(args.input)
    if not input_path.exists():
        print(f"Error: Input file not found: {input_path}", file=sys.stderr)
        sys.exit(1)

    # Detect format
    if args.format != "auto":
        detected_format = args.format
    else:
        detected_format = detect_file_format(input_path)
    if detected_format == "unknown":
        print(
            f"Error: Could not detect format for {input_path}. "
            "Use -f to specify: airr, pcp, or olmsted",
            file=sys.stderr,
        )
        sys.exit(1)

    # Load data based on format
    print(f"Reading {detected_format.upper()} data from: {input_path.name}", file=sys.stderr)

    if detected_format == "olmsted":
        all_clones, all_trees = _load_olmsted(input_path)
    elif detected_format == "pcp":
        trees_path = Path(args.tree) if args.tree else None
        all_clones, all_trees = _load_pcp(
            input_path, trees_path, args.seed, args.compute_metrics
        )
    elif detected_format == "airr":
        all_clones, all_trees = _load_airr(input_path)
    else:
        print(f"Error: Unsupported format: {detected_format}", file=sys.stderr)
        sys.exit(1)

    clone_count = len(all_clones)
    tree_count = len(all_trees)
    print(
        f"Found {clone_count} clones, {tree_count} trees",
        file=sys.stderr,
    )

    # Build YAML
    output_text = _build_yaml(
        input_path.name, detected_format, all_clones, all_trees
    )

    # Write output
    if args.output:
        output_path = Path(args.output)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w") as f:
            f.write(output_text)
        print(f"Config written to: {output_path}", file=sys.stderr)
    else:
        print(output_text)


if __name__ == "__main__":
    main()
