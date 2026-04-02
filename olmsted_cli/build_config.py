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
    MAX_NODES_SAMPLE,
    MAX_SAMPLE_HEURISTIC,
    MAX_SAMPLE_PATH,
    MAX_SAMPLE_PREVIEW,
    MAX_SAMPLE_VALUES,
    SUGGESTED_DISPLAY_MODES,
    SUGGESTED_SKIP_FIELDS,
)
from .field_metadata import (
    collect_keys,
    collect_mutations,
    collect_nodes,
    compute_range,
    entry_from_known,
    humanize_label,
    infer_field_type,
    sample_values,
    sample_values_by_path,
)
from .format_detection import detect_file_format


def get_args():
    """Parse command line arguments for the build-config command."""
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
    olmsted tag -i data.json -o tagged.json -c config.yaml
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
    # Deferred: process_pcp_data pulls in ete3 and other heavy dependencies
    from .process_pcp_data import (
        deterministic_uuid,
        parse_newick_csv,
        parse_pcp_csv,
        process_pcp_to_olmsted,
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

    # Deferred: process_airr_data pulls in heavy processing dependencies
    from .process_airr_data import process_dataset

    with open(input_path) as f:
        data = json.load(f)

    # Handle single-dataset or multi-dataset AIRR
    if isinstance(data, list):
        datasets_raw = data
    else:
        datasets_raw = [data]

    args = Namespace(
        root=None,           # no rooting for build-config introspection
        naive_name="naive",  # default for process_dataset compatibility
        root_trees=False,    # no rooting for build-config introspection
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


def _looks_like_local_path(values):
    """Check if sample string values look like local file paths (not URLs)."""
    path_prefixes = ("/", "./", "../", "~")
    path_count = 0
    for v in values:
        if not isinstance(v, str):
            continue
        if any(v.startswith(p) for p in path_prefixes) or "\\" in v:
            # Exclude URLs
            if not v.startswith(("http://", "https://", "ftp://")):
                path_count += 1
    return path_count > 0 and path_count >= len(values) // 2


def _check_mutation_demotion(nodes, field, max_samples=MAX_SAMPLE_HEURISTIC):
    """Check if a node-level field contains per-position mutation data.

    Detects three encodings:
        - list: dense array matching sequence length, continuous/aa/dna values
        - json: dict with int keys within sequence range, continuous/aa/dna values
        - records: array of dicts with "site" key (returns inner field names)

    Returns a dict with demotion info if eligible, or None:
        {"encoding": "list"|"json"|"records", "inner_type": str, ...}
    """
    values = sample_values(nodes, field, max_samples=max_samples)
    if not values:
        return None

    field_type = infer_field_type(values)

    # Get sequence lengths for comparison
    seq_lengths = set()
    for node in nodes[:max_samples]:
        if not isinstance(node, dict):
            continue
        for seq_field in ("sequence_alignment_aa", "sequence_alignment"):
            seq = node.get(seq_field)
            if isinstance(seq, str) and len(seq) > 0:
                seq_lengths.add(len(seq))
                break

    if field_type == "list":
        # Could be a dense per-position array or a records-style array of dicts
        # Check for records-style first: list of dicts with "site" key
        if values and all(isinstance(v, list) for v in values):
            sample_items = [item for v in values for item in v[:5] if item is not None]
            if sample_items and all(isinstance(item, dict) for item in sample_items):
                if any("site" in item for item in sample_items):
                    # Surprise-style: collect all inner field names (excluding "site")
                    inner_fields = {}
                    for v in values:
                        for entry in v:
                            if not isinstance(entry, dict):
                                continue
                            for k, val in entry.items():
                                if k == "site" or val is None:
                                    continue
                                if k not in inner_fields:
                                    inner_fields[k] = []
                                inner_fields[k].append(val)
                    if inner_fields:
                        # Infer type for each inner field
                        inner_types = {}
                        for k, vals in inner_fields.items():
                            inner_types[k] = infer_field_type(vals)
                        return {
                            "encoding": "records",
                            "source": field,
                            "inner_fields": inner_types,
                        }

            # Dense list: check if lengths match sequence length
            list_lengths = {len(v) for v in values if isinstance(v, list)}
            if list_lengths and (list_lengths & seq_lengths):
                inner_values = [item for v in values for item in v if item is not None]
                if inner_values:
                    inner_type = infer_field_type(inner_values)
                    if inner_type in ("continuous", "aa", "dna"):
                        return {"encoding": "list", "inner_type": inner_type}

    elif field_type == "json":
        # Check if keys are parseable as ints within sequence range
        all_keys_int = True
        max_key = -1
        inner_values = []
        for v in values:
            if not isinstance(v, dict):
                return None
            for k, val in v.items():
                try:
                    k_int = int(k)
                    max_key = max(max_key, k_int)
                except (ValueError, TypeError):
                    all_keys_int = False
                    break
                if val is not None:
                    inner_values.append(val)
            if not all_keys_int:
                break
        if not all_keys_int or max_key < 0 or not inner_values:
            return None
        # Check max key is within sequence range
        if seq_lengths and max_key >= max(seq_lengths):
            return None
        inner_type = infer_field_type(inner_values)
        if inner_type in ("continuous", "aa", "dna"):
            return {"encoding": "json", "inner_type": inner_type}

    return None


def _field_summary(dicts, field, known_registry):
    """Build a summary dict for a single field: type, display, label."""
    if field in known_registry:
        entry = entry_from_known(known_registry[field])
    else:
        values = sample_values(dicts, field, max_samples=MAX_SAMPLE_VALUES)
        inferred_type = infer_field_type(values)
        display = "tooltip" if inferred_type in ("list", "json") else "dropdown"
        entry = {
            "type": inferred_type,
            "display": display,
            "label": humanize_label(field),
        }
    # Apply suggested display mode override
    if field in SUGGESTED_DISPLAY_MODES:
        entry["display"] = SUGGESTED_DISPLAY_MODES[field]
    return entry


def _should_skip(field, dicts=None, no_skip=False, skip_all=False):
    """Determine if a field should be suggested as skip in the config.

    Args:
        field: Field name.
        dicts: Optional data dicts for heuristic checks (e.g., file paths).
        no_skip: If True, never suggest skip.
        skip_all: If True, always suggest skip.
    """
    if no_skip:
        return False
    if skip_all:
        return True
    if field in SUGGESTED_SKIP_FIELDS:
        return True
    if dicts is not None:
        values = sample_values(dicts, field, max_samples=MAX_SAMPLE_HEURISTIC)
        if _looks_like_local_path(values):
            return True
    return False


def _make_field_entry(name, level, entry, skip=False, encoding=None, source=None):
    """Build a custom_fields dict entry matching load_config() format."""
    d = {"name": name, "level": level, "type": entry["type"], "label": entry["label"]}
    # Apply cross-format aliases (e.g., rearrangement_count → unique_seqs_count)
    alias = FIELD_ALIASES.get(name)
    if alias and alias != name:
        d["output_name"] = alias
    display = entry.get("display", "dropdown")
    if display != "dropdown":
        d["display"] = display
    if skip:
        d["skip"] = True
    if encoding:
        d["encoding"] = encoding
    if source:
        d["source"] = source
    return d


def generate_default_config(
    clones, trees, *, no_skip=False, skip_all=False,
    _nodes=None, _mutations=None,
):
    """Generate default field declarations by introspecting data.

    Discovers fields at each level (clone, node, branch, mutation), infers
    types from data, applies skip suggestions, and detects mutation demotion.

    This is the single source of truth for field discovery. Both the
    ``build-config`` YAML output and the ``process``/``tag`` pipelines
    use this function to determine what fields exist and how they should
    be configured.

    Args:
        clones: List of clone dicts.
        trees: List of tree dicts.
        no_skip: If True, no fields are suggested as skip.
        skip_all: If True, all fields are suggested as skip.
        _nodes: Pre-collected node dicts (internal optimization to avoid
            redundant tree traversal when caller already has them).
        _mutations: Pre-collected mutation dicts (same optimization).

    Returns:
        List of custom_field dicts in the same format as ``load_config()``.
        Each dict has keys: name, level, type, label, and optionally
        skip, display, encoding, source.
    """
    all_nodes = _nodes if _nodes is not None else collect_nodes(trees, max_nodes=MAX_NODES_SAMPLE)
    all_mutations = _mutations if _mutations is not None else collect_mutations(trees)
    fields = []

    # --- Clone level ---
    clone_keys = collect_keys(clones) - EXCLUDED_CLONE_FIELDS
    # Check for known fields with dot-paths (e.g., locus → sample.locus)
    for field_name, field_info in KNOWN_CLONE_FIELDS.items():
        if "path" in field_info and field_name not in clone_keys:
            values = sample_values_by_path(clones, field_info["path"], max_samples=MAX_SAMPLE_PATH)
            if values:
                clone_keys.add(field_name)

    for field in sorted(clone_keys):
        # Only include fields that have actual non-null values
        if field in KNOWN_CLONE_FIELDS and "path" in KNOWN_CLONE_FIELDS[field]:
            values = sample_values_by_path(clones, KNOWN_CLONE_FIELDS[field]["path"])
        else:
            values = sample_values(clones, field)
        if not values:
            continue
        entry = _field_summary(clones, field, KNOWN_CLONE_FIELDS)
        skip = _should_skip(field, clones, no_skip, skip_all)
        fields.append(_make_field_entry(field, "clone", entry, skip=skip))

    # --- Node level ---
    node_keys = collect_keys(all_nodes) - EXCLUDED_NODE_FIELDS
    node_keys -= set(KNOWN_BRANCH_FIELDS.keys())

    for field in sorted(node_keys):
        # Only include fields that have actual non-null values
        values = sample_values(all_nodes, field)
        if not values:
            continue
        demotion = _check_mutation_demotion(all_nodes, field)
        if demotion:
            # Add node-level skip entry so users can opt back in
            entry = _field_summary(all_nodes, field, KNOWN_NODE_FIELDS)
            fields.append(_make_field_entry(field, "node", entry, skip=True))
            # Add mutation-level entries with encoding
            enc = demotion["encoding"]
            if enc == "records":
                for inner_name, inner_type in sorted(demotion["inner_fields"].items()):
                    inner_entry = {
                        "type": inner_type,
                        "display": "dropdown",
                        "label": humanize_label(inner_name),
                    }
                    fields.append(_make_field_entry(
                        inner_name, "mutation", inner_entry,
                        encoding="records", source=field,
                    ))
            else:
                # list or json: single mutation entry
                mut_entry = dict(entry)
                mut_entry["type"] = demotion["inner_type"]
                fields.append(_make_field_entry(
                    field, "mutation", mut_entry, encoding=enc,
                ))
        else:
            entry = _field_summary(all_nodes, field, KNOWN_NODE_FIELDS)
            skip = _should_skip(field, all_nodes, no_skip, skip_all)
            fields.append(_make_field_entry(field, "node", entry, skip=skip))

    # --- Branch level ---
    branch_keys = collect_keys(all_nodes) & set(KNOWN_BRANCH_FIELDS.keys())
    for field in sorted(branch_keys):
        values = sample_values(all_nodes, field)
        if not values:
            continue
        entry = entry_from_known(KNOWN_BRANCH_FIELDS[field])
        fields.append(_make_field_entry(field, "branch", entry))

    # --- Mutation level ---
    has_aa_sequences = any(
        n.get("sequence_alignment_aa") for n in all_nodes if isinstance(n, dict)
    )
    mutation_keys = collect_keys(all_mutations) - EXCLUDED_MUTATION_FIELDS
    # Track mutation field names already emitted (e.g., from demotion)
    emitted_mutation_names = {f["name"] for f in fields if f["level"] == "mutation"}
    has_derived_aa = (
        has_aa_sequences
        and "child_aa" not in mutation_keys
        and "child_aa" not in emitted_mutation_names
    )

    if has_derived_aa:
        fields.append(_make_field_entry(
            "child_aa", "mutation",
            {"type": "aa", "display": "dropdown", "label": "Child Amino Acid"},
        ))
        fields.append(_make_field_entry(
            "parent_aa", "mutation",
            {"type": "aa", "display": "tooltip", "label": "Parent Amino Acid"},
        ))

    for field in sorted(mutation_keys):
        values = sample_values(all_mutations, field)
        if not values:
            continue
        entry = _field_summary(all_mutations, field, KNOWN_MUTATION_FIELDS)
        skip = _should_skip(field, all_mutations, no_skip, skip_all)
        fields.append(_make_field_entry(field, "mutation", entry, skip=skip))

    return fields


def _format_field_block(
    name, level, entry, sample_values=None, field_range=None,
    skip=False, encoding=None, source=None,
):
    """Format a single custom_fields YAML entry as a string."""
    lines = []
    lines.append(f"  - name: {name}")

    # Suggest output_name from alias table if applicable
    alias = FIELD_ALIASES.get(name)
    if alias and alias != name:
        lines.append(f"    output_name: {alias}")

    lines.append(f"    level: {level}")
    if encoding:
        lines.append(f"    encoding: {encoding}")
    if source:
        lines.append(f"    source: {source}")
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


def _get_sample_values_for_field(cf, all_clones, all_nodes, all_mutations):
    """Get sample values for a config field entry, for YAML preview."""
    level = cf["level"]
    name = cf["name"]

    if level == "clone":
        # Path-based fields (e.g., locus → sample.locus)
        if name in KNOWN_CLONE_FIELDS and "path" in KNOWN_CLONE_FIELDS[name]:
            return sample_values_by_path(all_clones, KNOWN_CLONE_FIELDS[name]["path"],
                                         max_samples=MAX_SAMPLE_PREVIEW)
        return sample_values(all_clones, name, max_samples=MAX_SAMPLE_PREVIEW)
    elif level == "node":
        return sample_values(all_nodes, name, max_samples=MAX_SAMPLE_PREVIEW)
    elif level == "branch":
        return sample_values(all_nodes, name, max_samples=MAX_SAMPLE_PREVIEW)
    elif level == "mutation":
        if cf.get("encoding"):
            # Encoded fields: sample from nodes (source field), not mutations
            return sample_values(all_nodes, cf.get("source", name),
                                 max_samples=MAX_SAMPLE_PREVIEW)
        return sample_values(all_mutations, name, max_samples=MAX_SAMPLE_PREVIEW)
    return []


def _build_yaml(
    input_name, detected_format, all_clones, all_trees,
    input_path=None, tree_path=None,
    no_skip=False, skip_all=False,
):
    """Build the YAML config string from templates and introspected data."""
    all_nodes = collect_nodes(all_trees, max_nodes=MAX_NODES_SAMPLE)
    all_mutations = collect_mutations(all_trees)

    # Generate structured config, reusing already-collected nodes/mutations
    config_fields = generate_default_config(
        all_clones, all_trees, no_skip=no_skip, skip_all=skip_all,
        _nodes=all_nodes, _mutations=all_mutations,
    )

    input_str = str(input_path) if input_path else input_name
    tree_str = str(tree_path) if tree_path else "trees.csv"

    # Determine usage command for header
    if detected_format == FORMAT_OLMSTED:
        usage_command = "olmsted tag -i data.json -o tagged.json -c this_file.yaml"
    else:
        usage_command = "olmsted process -c this_file.yaml"

    # Assemble from templates
    parts = []

    parts.append(_load_template("header.yaml").format(
        input_name=input_name,
        detected_format=detected_format,
        usage_command=usage_command,
    ))

    options_template = f"options_{detected_format}.yaml"
    parts.append(_load_template(options_template).format(
        input_path=input_str,
        tree_path=tree_str,
    ))

    parts.append(_load_template("fields_header.yaml"))

    # Group config fields by level for presentation
    lines = []
    skip_entries = []

    # Level display names and section headers
    level_headers = {
        "clone": "  # --- Family level (clonal family — scatterplot axes, color, facet) ---",
        "node": "  # --- Node level (tree node properties, tooltips) ---",
        "branch": "  # --- Branch level (tree branch coloring, width) ---",
        "mutation": "  # --- Mutation level (alignment coloring) ---",
    }
    # Config uses "family" for clone level in YAML
    level_yaml_names = {"clone": "family", "node": "node", "branch": "branch", "mutation": "mutation"}

    # Separate active and skip entries, preserving order within each level
    for level in ("clone", "node", "branch", "mutation"):
        yaml_level = level_yaml_names[level]
        level_fields = [cf for cf in config_fields if cf["level"] == level]
        active = [cf for cf in level_fields if not cf.get("skip")]
        skipped = [cf for cf in level_fields if cf.get("skip")]

        # Mutation level has special sub-sections
        if level == "mutation":
            derived = [cf for cf in active if not cf.get("encoding")]
            demoted = [cf for cf in active if cf.get("encoding")]
            # Only those without encoding and without a source → derived AA fields
            derived_aa = [cf for cf in derived if cf["name"] in ("child_aa", "parent_aa")
                          and cf["type"] in ("aa", "dna")]
            regular = [cf for cf in derived if cf not in derived_aa]

            if derived_aa or demoted or regular:
                lines.append("")
                lines.append(level_headers[level])

            if derived_aa:
                lines.append("  # The following fields are derived by the web app from")
                lines.append("  # parent/child sequence alignments during rendering:")
                for cf in derived_aa:
                    entry = {"type": cf["type"], "display": cf.get("display", "dropdown"),
                             "label": cf["label"]}
                    lines.append(_format_field_block(cf["name"], yaml_level, entry))

            if demoted:
                lines.append("  # The following fields were detected as per-position data")
                lines.append("  # stored on nodes (demoted from node to mutation level):")
                for cf in demoted:
                    entry = {"type": cf["type"], "display": cf.get("display", "dropdown"),
                             "label": cf["label"]}
                    lines.append(_format_field_block(
                        cf["name"], yaml_level, entry,
                        encoding=cf.get("encoding"), source=cf.get("source"),
                    ))

            for cf in regular:
                entry = {"type": cf["type"], "display": cf.get("display", "dropdown"),
                         "label": cf["label"]}
                samples = _get_sample_values_for_field(cf, all_clones, all_nodes, all_mutations)
                field_range = None
                if entry["type"] == "continuous":
                    field_range = compute_range(all_mutations, cf["name"])
                lines.append(_format_field_block(cf["name"], yaml_level, entry, samples, field_range))

        elif active:
            lines.append("")
            lines.append(level_headers[level])
            for cf in active:
                entry = {"type": cf["type"], "display": cf.get("display", "dropdown"),
                         "label": cf["label"]}
                samples = _get_sample_values_for_field(cf, all_clones, all_nodes, all_mutations)
                lines.append(_format_field_block(cf["name"], yaml_level, entry, samples))

        for cf in skipped:
            entry = {"type": cf["type"], "display": cf.get("display", "dropdown"),
                     "label": cf["label"]}
            samples = _get_sample_values_for_field(cf, all_clones, all_nodes, all_mutations)
            skip_entries.append((cf["name"], yaml_level, entry, samples, None))

    # --- Skipped fields section (at the bottom) ---
    if skip_entries:
        parts.append(_load_template("skip_header.yaml"))
        for field, level, entry, samples, field_range in skip_entries:
            lines.append(_format_field_block(field, level, entry, samples, field_range, skip=True))

    lines.append("")

    # Assemble: templates + dynamic field lines
    return "".join(parts) + "\n".join(lines)


def main():
    """Main entry point for the build-config command."""
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
