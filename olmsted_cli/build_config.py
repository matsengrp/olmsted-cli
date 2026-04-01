"""
Build-config command: generate a YAML configuration from data for editing.

Introspects the input data (Olmsted JSON, AIRR JSON, or PCP CSV), discovers
all available fields, and generates a ready-to-edit YAML config with
processing options, field declarations, and cross-format alias suggestions.

Usage:
    olmsted build-config -i data.json -o config.yaml
    olmsted build-config -i pcp.csv -t trees.csv -o config.yaml
    olmsted build-config -i data.json                   # prints to stdout
"""

import argparse
import json
import sys
from pathlib import Path

from .process_data import detect_file_format
from .constants import (
    EXCLUDED_CLONE_FIELDS,
    EXCLUDED_MUTATION_FIELDS,
    EXCLUDED_NODE_FIELDS,
    FIELD_ALIASES,
    FORMAT_AIRR,
    FORMAT_AUTO,
    FORMAT_OLMSTED,
    FORMAT_PCP,
    FORMAT_UNKNOWN,
    KNOWN_BRANCH_FIELDS,
    KNOWN_CLONE_FIELDS,
    KNOWN_MUTATION_FIELDS,
    KNOWN_NODE_FIELDS,
    SUGGESTED_DISPLAY_MODES,
    SUGGESTED_SKIP_FIELDS,
)
from .field_metadata import (
    compute_range,
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
        description="Generate a YAML config from your data for editing",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # From PCP CSV (with optional trees)
    olmsted build-config -i pcp.csv -t trees.csv -o config.yaml

    # From AIRR JSON
    olmsted build-config -i airr_data.json -o config.yaml

    # From existing Olmsted JSON
    olmsted build-config -i data.json -o config.yaml

    # Print to stdout for inspection
    olmsted build-config -i data.json

    # Then edit the config and use it
    olmsted process -c config.yaml
    olmsted enrich -i data.json -o enriched.json -c config.yaml
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
        choices=[FORMAT_AIRR, FORMAT_PCP, FORMAT_OLMSTED, FORMAT_AUTO],
        default=FORMAT_AUTO,
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

    parser.add_argument(
        "-v", "--verbose",
        type=int,
        choices=[0, 1, 2, 3],
        default=1,
        help="Verbosity: 0=errors only, 1=normal (default), 2=verbose, 3=debug",
    )

    # Skip mode for config generation
    skip_group = parser.add_mutually_exclusive_group()
    skip_group.add_argument(
        "--no-skip",
        action="store_true",
        help="Don't suggest any skips — all fields get display: dropdown",
    )
    skip_group.add_argument(
        "--skip-all",
        action="store_true",
        help="Skip all fields by default — user must manually un-skip what they want",
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
    """Build a summary dict for a single field: type, display, label."""
    if field in known_registry:
        known = known_registry[field]
        entry = {
            "type": known["type"],
            "display": known.get("display", "dropdown"),
            "label": known["label"],
        }
    else:
        values = _sample_values(dicts, field, max_samples=50)
        entry = {
            "type": infer_field_type(values),
            "display": "dropdown",
            "label": humanize_label(field),
        }
    # Apply suggested display mode override
    if field in SUGGESTED_DISPLAY_MODES:
        entry["display"] = SUGGESTED_DISPLAY_MODES[field]
    return entry


def _format_field_block(name, level, entry, sample_values=None, field_range=None, skip=False):
    """Format a single custom_fields YAML entry as a string."""
    lines = []
    lines.append(f"  - name: {name}")

    # Suggest output_name from alias table if applicable
    alias = FIELD_ALIASES.get(name)
    if alias and alias != name:
        lines.append(f"    output_name: {alias}")

    lines.append(f"    level: {level}")
    if skip:
        lines.append(f"    skip: true")
    lines.append(f"    type: {entry['type']}")
    display = entry.get("display", "dropdown")
    if display != "dropdown":
        lines.append(f"    display: {display}")
    lines.append(f"    label: \"{entry['label']}\"")
    if field_range:
        lines.append(f"    # range in data: [{field_range[0]}, {field_range[1]}]")
        lines.append(f"    # range: [{field_range[0]}, {field_range[1]}]  # uncomment to set color scale domain")
    if sample_values:
        preview = ", ".join(str(v) for v in sample_values[:5])
        if len(sample_values) > 5:
            preview += ", ..."
        lines.append(f"    # sample values: {preview}")
    return "\n".join(lines)


def _load_template(name):
    """Load a template file from configs/templates/."""
    template_dir = Path(__file__).parent / "configs" / "templates"
    template_path = template_dir / name
    return template_path.read_text()


def _build_yaml(
    input_name, detected_format, all_clones, all_trees,
    input_path=None, tree_path=None,
    no_skip=False, skip_all=False,
):
    """Build the YAML config string from templates and introspected data."""
    all_nodes = _collect_nodes(all_trees, max_nodes=500)
    all_mutations = _collect_mutations(all_trees, max_mutations=500)

    input_str = str(input_path) if input_path else input_name
    tree_str = str(tree_path) if tree_path else "trees.csv"

    # Determine usage command for header
    if detected_format == FORMAT_OLMSTED:
        usage_command = "olmsted enrich -i data.json -o enriched.json -c this_file.yaml"
    else:
        usage_command = "olmsted process -c this_file.yaml"

    # Assemble from templates
    parts = []

    # Header
    parts.append(_load_template("header.yaml").format(
        input_name=input_name,
        detected_format=detected_format,
        usage_command=usage_command,
    ))

    # Processing options (per-format)
    options_template = f"options_{detected_format}.yaml"
    parts.append(_load_template(options_template).format(
        input_path=input_str,
        tree_path=tree_str,
    ))

    # Field declarations header
    parts.append(_load_template("fields_header.yaml"))

    # Build the fields section
    lines = []

    # Collect skip entries separately for the bottom section
    skip_entries = []

    def _is_skip(field):
        if no_skip:
            return False
        if skip_all:
            return True
        return field in SUGGESTED_SKIP_FIELDS

    # --- Clone level ---
    clone_keys = _collect_keys(all_clones) - EXCLUDED_CLONE_FIELDS
    has_locus = any(
        isinstance(c.get("sample"), dict) and c["sample"].get("locus") is not None
        for c in all_clones[:20]
    )
    if has_locus:
        clone_keys.add("locus")

    active_clone = [f for f in sorted(clone_keys) if not _is_skip(f)]
    skip_clone = [f for f in sorted(clone_keys) if _is_skip(f)]

    if active_clone:
        lines.append("")
        lines.append("  # --- Family level (clonal family — scatterplot axes, color, facet) ---")
        for field in active_clone:
            entry = _field_summary(all_clones, field, KNOWN_CLONE_FIELDS)
            if field == "locus":
                samples = list({
                    c["sample"]["locus"]
                    for c in all_clones[:50]
                    if isinstance(c.get("sample"), dict) and c["sample"].get("locus")
                })
            else:
                samples = _sample_values(all_clones, field, max_samples=6)
            lines.append(_format_field_block(field, "family", entry, samples))

    for field in skip_clone:
        entry = _field_summary(all_clones, field, KNOWN_CLONE_FIELDS)
        samples = _sample_values(all_clones, field, max_samples=6)
        skip_entries.append((field, "family", entry, samples, None))

    # --- Node level ---
    node_keys = _collect_keys(all_nodes) - EXCLUDED_NODE_FIELDS
    node_keys -= set(KNOWN_BRANCH_FIELDS.keys())

    active_node = [f for f in sorted(node_keys) if not _is_skip(f)]
    skip_node = [f for f in sorted(node_keys) if _is_skip(f)]

    if active_node:
        lines.append("")
        lines.append("  # --- Node level (tree node properties, tooltips) ---")
        for field in active_node:
            entry = _field_summary(all_nodes, field, KNOWN_NODE_FIELDS)
            samples = _sample_values(all_nodes, field, max_samples=6)
            lines.append(_format_field_block(field, "node", entry, samples))

    for field in skip_node:
        entry = _field_summary(all_nodes, field, KNOWN_NODE_FIELDS)
        samples = _sample_values(all_nodes, field, max_samples=6)
        skip_entries.append((field, "node", entry, samples, None))

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
    all_mutations_full = _collect_mutations(all_trees, max_mutations=100000)

    has_aa_sequences = any(
        n.get("sequence_alignment_aa") for n in all_nodes if isinstance(n, dict)
    )
    has_derived_aa = has_aa_sequences and "child_aa" not in mutation_keys

    active_mutation = [f for f in sorted(mutation_keys) if not _is_skip(f)]
    skip_mutation = [f for f in sorted(mutation_keys) if _is_skip(f)]

    if active_mutation or has_derived_aa:
        lines.append("")
        lines.append("  # --- Mutation level (alignment coloring) ---")

        if has_derived_aa:
            lines.append("  # The following fields are derived by the web app from")
            lines.append("  # parent/child sequence alignments during rendering:")
            lines.append(_format_field_block(
                "child_aa", "mutation",
                {"type": "aa", "label": "Child Amino Acid"},
            ))
            lines.append(_format_field_block(
                "parent_aa", "mutation",
                {"type": "aa", "display": "tooltip", "label": "Parent Amino Acid"},
            ))

        for field in active_mutation:
            entry = _field_summary(all_mutations, field, KNOWN_MUTATION_FIELDS)
            samples = _sample_values(all_mutations, field, max_samples=6)
            field_range = None
            if entry["type"] == "continuous":
                field_range = compute_range(all_mutations_full, field)
            lines.append(_format_field_block(field, "mutation", entry, samples, field_range))

    for field in skip_mutation:
        entry = _field_summary(all_mutations, field, KNOWN_MUTATION_FIELDS)
        samples = _sample_values(all_mutations, field, max_samples=6)
        skip_entries.append((field, "mutation", entry, samples, None))

    # --- Skipped fields section (at the bottom) ---
    if skip_entries:
        parts.append(_load_template("skip_header.yaml"))
        for field, level, entry, samples, field_range in skip_entries:
            lines.append(_format_field_block(field, level, entry, samples, field_range, skip=True))

    lines.append("")

    # Assemble: templates + dynamic field lines
    return "".join(parts) + "\n".join(lines)


def main():
    """Main entry point for the dump-fields command."""
    args = get_args()

    input_path = Path(args.input)
    if not input_path.exists():
        print(f"Error: Input file not found: {input_path}", file=sys.stderr)
        sys.exit(1)

    # Detect format
    if args.format != FORMAT_AUTO:
        detected_format = args.format
    else:
        detected_format = detect_file_format(input_path)
    if detected_format == FORMAT_UNKNOWN:
        print(
            f"Error: Could not detect format for {input_path}. "
            "Use -f to specify: airr, pcp, or olmsted",
            file=sys.stderr,
        )
        sys.exit(1)

    # Load data based on format
    print(f"Reading {detected_format.upper()} data from: {input_path.name}", file=sys.stderr)

    if detected_format == FORMAT_OLMSTED:
        all_clones, all_trees = _load_olmsted(input_path)
    elif detected_format == FORMAT_PCP:
        trees_path = Path(args.tree) if args.tree else None
        all_clones, all_trees = _load_pcp(
            input_path, trees_path, args.seed, args.compute_metrics
        )
    elif detected_format == FORMAT_AIRR:
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
    trees_path = Path(args.tree) if args.tree else None
    output_text = _build_yaml(
        input_path.name, detected_format, all_clones, all_trees,
        input_path=input_path, tree_path=trees_path,
        no_skip=args.no_skip, skip_all=args.skip_all,
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
