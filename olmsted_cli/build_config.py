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
# Cross-format field aliases: suggested output_name renames for common fields.
# These are NOT auto-applied during processing — they are pre-filled as
# suggestions in the generated config for the user to approve or remove.
FIELD_ALIASES = {
    "v_gene": "v_call",
    "v_gene_heavy": "v_call",
    "d_gene": "d_call",
    "d_gene_heavy": "d_call",
    "j_gene": "j_call",
    "j_gene_heavy": "j_call",
    "v_gene_light": "v_call_light",
    "j_gene_light": "j_call_light",
    "rearrangement_count": "unique_seqs_count",
    "sampled_seqs_count": "unique_seqs_count",
    "size": "total_read_count",
    "branch_length": "length",
    "mut_to": "child_aa",
    "mut_from": "parent_aa",
}

from .field_metadata import (
    EXCLUDED_CLONE_FIELDS,
    EXCLUDED_MUTATION_FIELDS,
    EXCLUDED_NODE_FIELDS,
    KNOWN_BRANCH_FIELDS,
    KNOWN_CLONE_FIELDS,
    KNOWN_MUTATION_FIELDS,
    KNOWN_NODE_FIELDS,
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


def _format_field_block(name, level, entry, sample_values=None, field_range=None):
    """Format a single custom_fields YAML entry as a string."""
    lines = []
    lines.append(f"  - name: {name}")

    # Suggest output_name from alias table if applicable
    alias = FIELD_ALIASES.get(name)
    if alias and alias != name:
        lines.append(f"    output_name: {alias}")

    lines.append(f"    level: {level}")
    lines.append(f"    type: {entry['type']}")
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


def _build_yaml(
    input_name, detected_format, all_clones, all_trees,
    input_path=None, tree_path=None,
):
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
    lines.append("#   - Change type (continuous, categorical, tooltip, aa, dna)")
    lines.append("#   - Customize display labels")
    lines.append("#   - Uncomment processing options below as needed")
    lines.append("#")
    if detected_format == "olmsted":
        lines.append("# Then use with:  olmsted enrich -i data.json -o enriched.json -c this_file.yaml")
    else:
        lines.append("# Then use with:  olmsted process -c this_file.yaml")
    lines.append("")

    # --- Processing options template ---
    lines.append("")
    lines.append("# =============================================================================")
    lines.append("# Processing Options (uncomment and edit as needed)")
    lines.append("# =============================================================================")
    lines.append("")

    if detected_format == "olmsted":
        # Olmsted JSON uses enrich, not process
        input_str = str(input_path) if input_path else input_name
        lines.append(f"# input: {input_str}")
        lines.append("# output: enriched_output.json")
        lines.append("# mode: add            # add (merge with existing) or overwrite")
    elif detected_format == "pcp":
        input_str = str(input_path) if input_path else input_name
        lines.append(f"# inputs: [{input_str}]")
        if tree_path:
            lines.append(f"# tree: {tree_path}")
        else:
            lines.append("# tree: trees.csv")
        lines.append("# output: output.json")
        lines.append("# format: pcp")
        lines.append("")
        lines.append("# name: \"My Dataset\"")
        lines.append("# seed: 42             # for reproducible UUIDs")
        lines.append("")
        lines.append("# --- Metric Computation ---")
        lines.append("# Computes LBI, LBR, affinity, and scaled_affinity for all tree nodes.")
        lines.append("# Requires tree branch lengths (works with any format that has them).")
        lines.append("# compute_metrics: true")
        lines.append("# lbi_tau: 0.0125      # time scale parameter for LBI")
        lines.append("")
        lines.append("# standardize_names: false  # rename nodes to naive, Node1, Leaf1, ...")
        lines.append("# validate: false")
        lines.append("# verbose: 1")
    elif detected_format == "airr":
        input_str = str(input_path) if input_path else input_name
        lines.append(f"# inputs: [{input_str}]")
        lines.append("# output: output.json")
        lines.append("# format: airr")
        lines.append("")
        lines.append("# name: \"My Dataset\"")
        lines.append("# seed: 42             # for reproducible UUIDs")
        lines.append("")
        lines.append("# --- Metric Computation ---")
        lines.append("# Computes LBI, LBR, affinity, and scaled_affinity for all tree nodes.")
        lines.append("# Requires tree branch lengths (works with any format that has them).")
        lines.append("# compute_metrics: true")
        lines.append("# lbi_tau: 0.0125      # time scale parameter for LBI")
        lines.append("")
        lines.append("# naive_name: naive     # name of root/naive node")
        lines.append("# root_trees: false     # re-root trees at naive node")
        lines.append("# validate: false")
        lines.append("# verbose: 1")

    lines.append("")
    lines.append("")
    lines.append("# =============================================================================")
    lines.append("# Field Declarations")
    lines.append("# =============================================================================")
    lines.append("#")
    lines.append("# Each field entry supports:")
    lines.append("#   name:        Field name in the input data (required)")
    lines.append("#   output_name: Renamed field in output (optional, for cross-format alignment)")
    lines.append("#   level:       clone, node, branch, or mutation (required)")
    lines.append("#   type:        continuous, categorical, tooltip, aa, or dna (required)")
    lines.append("#   label:       Display label in web app (required)")
    lines.append("#   range:       [min, max] for continuous fields (optional)")
    lines.append("#")
    lines.append("# Types:")
    lines.append("#   continuous  — numeric (axes, size, color scales)")
    lines.append("#   categorical — string/enum (color, shape, facet)")
    lines.append("#   tooltip     — display-only (shown in tooltips, not for encoding)")
    lines.append("#   aa          — amino acid identity (uses full genetic alphabet)")
    lines.append("#   dna         — nucleotide identity (uses full genetic alphabet)")
    lines.append("#   skip        — exclude this field from metadata (keeps entry for docs)")
    lines.append("#")
    lines.append("# Cross-format aliases (suggested output_name, remove if not needed):")
    lines.append("#   v_gene, v_gene_heavy  ->  v_call")
    lines.append("#   d_gene, d_gene_heavy  ->  d_call")
    lines.append("#   j_gene, j_gene_heavy  ->  j_call")
    lines.append("#   rearrangement_count   ->  unique_seqs_count")
    lines.append("#   size                  ->  total_read_count")
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
    # Collect all mutations for accurate range computation
    all_mutations_full = _collect_mutations(all_trees, max_mutations=100000)

    # Check for AA sequence data on nodes — if present, the web app will
    # derive per-mutation child_aa/parent_aa during alignment rendering
    has_aa_sequences = any(
        n.get("sequence_alignment_aa") for n in all_nodes if isinstance(n, dict)
    )
    has_derived_aa = has_aa_sequences and "child_aa" not in mutation_keys

    if mutation_keys or has_derived_aa:
        lines.append("")
        lines.append("  # --- Mutation level (alignment coloring) ---")

        # Derived AA fields (from sequence alignment, computed by web app)
        if has_derived_aa:
            lines.append("  # The following fields are derived by the web app from")
            lines.append("  # parent/child sequence alignments during rendering:")
            lines.append(_format_field_block(
                "child_aa", "mutation",
                {"type": "aa", "label": "Child Amino Acid"},
            ))
            lines.append(_format_field_block(
                "parent_aa", "mutation",
                {"type": "tooltip", "label": "Parent Amino Acid"},
            ))

        # Pre-computed mutation fields (e.g., surprise scores)
        for field in sorted(mutation_keys):
            entry = _field_summary(all_mutations, field, KNOWN_MUTATION_FIELDS)
            samples = _sample_values(all_mutations, field, max_samples=6)
            field_range = None
            if entry["type"] == "continuous":
                field_range = compute_range(all_mutations_full, field)
            lines.append(_format_field_block(field, "mutation", entry, samples, field_range))

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
    trees_path = Path(args.tree) if args.tree else None
    output_text = _build_yaml(
        input_path.name, detected_format, all_clones, all_trees,
        input_path=input_path, tree_path=trees_path,
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
