"""
Field metadata generation for Olmsted datasets.

Generates field_metadata describing available data fields at each level
(clone, node, branch, mutation) with their types and human-readable labels.
This metadata drives dynamic dropdown construction in the Olmsted web app.

Field types (what the data IS):
    - "continuous": Numeric values suitable for axes, size, color scales
    - "categorical": String/enum values suitable for color, shape, facet
    - "aa": Amino acid single-character identity
    - "dna": Nucleotide single-character identity

Display modes (how to show it):
    - "dropdown": In visualization controls (default)
    - "tooltip": Show on hover only, not in controls

Levels:
    - "clone": Clone/clonal family level (scatterplot axes, color, facet)
    - "node": Tree node level (node properties, tooltips)
    - "branch": Tree branch level (branch coloring, width)
    - "mutation": Per-mutation level (alignment coloring)
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

from .constants import (
    AA_CHARS,
    ABBREVIATION_MAP,
    DNA_CHARS,
    EXCLUDED_BRANCH_FIELDS,
    EXCLUDED_CLONE_FIELDS,
    EXCLUDED_MUTATION_FIELDS,
    EXCLUDED_NODE_FIELDS,
    EXCLUDED_TREE_FIELDS,
    KNOWN_BRANCH_FIELDS,
    KNOWN_CLONE_FIELDS,
    KNOWN_MUTATION_FIELDS,
    KNOWN_NODE_FIELDS,
    KNOWN_TREE_FIELDS,
    MAX_MUTATIONS,
    MAX_NODES_QUICK,
    MAX_NODES_SAMPLE,
    MAX_SAMPLE_PATH,
    MAX_SAMPLE_VALUES,
    SUGGESTED_DISPLAY_MODES,
    SUGGESTED_SKIP_FIELDS,
    normalize_level,
)

# =============================================================================
# Utility functions
# =============================================================================


def infer_field_type(values: List[Any]) -> str:
    """
    Infer field data type from sample values.

    Args:
        values: Non-null sample values from the field.

    Returns:
        "continuous" if all values are numeric,
        "aa" if all values are single amino acid characters,
        "dna" if all values are single nucleotide characters,
        "list" if all values are lists,
        "json" if all values are dicts,
        "categorical" for strings, booleans, or mixed types.
    """
    if not values:
        return "categorical"

    numeric_count = 0
    string_count = 0
    list_count = 0
    dict_count = 0
    string_values = []
    for v in values:
        if isinstance(v, bool):
            string_count += 1
        elif isinstance(v, (int, float)):
            numeric_count += 1
        elif isinstance(v, str):
            string_count += 1
            string_values.append(v)
        elif isinstance(v, list):
            list_count += 1
        elif isinstance(v, dict):
            dict_count += 1
        else:
            return "categorical"

    # All lists → list type
    if list_count > 0 and list_count == len(values):
        return "list"

    # All dicts → json type
    if dict_count > 0 and dict_count == len(values):
        return "json"

    # Mixed complex types → categorical
    if list_count > 0 or dict_count > 0:
        return "categorical"

    if numeric_count > 0 and string_count == 0:
        return "continuous"
    if string_count > 0 and numeric_count == 0:
        # Check for single-character DNA or AA
        # DNA checked first: if values contain AA-only chars (e.g., D, E, F),
        # they can't be DNA. Pure ACGTU ambiguity is resolved by the known
        # fields registry (parent_aa/child_aa vs parent_nt/child_nt).
        if string_values and all(len(s) == 1 for s in string_values):
            upper_vals = {s.upper() for s in string_values}
            # If any char is AA-only (not in DNA alphabet), it's AA
            if upper_vals <= AA_CHARS and not upper_vals <= DNA_CHARS:
                return "aa"
            if upper_vals <= DNA_CHARS:
                return "dna"
            if upper_vals <= AA_CHARS:
                return "aa"
        return "categorical"
    return "categorical"  # mixed types


def compute_range(dicts: List[Dict], field: str) -> Optional[List[float]]:
    """
    Compute [min, max] range for a numeric field across a list of dicts.

    Returns None if no numeric values are found.
    """
    values = []
    for d in dicts:
        v = d.get(field)
        if isinstance(v, (int, float)) and not isinstance(v, bool):
            values.append(v)
    if not values:
        return None
    return [min(values), max(values)]


def humanize_label(field_name: str) -> str:
    """
    Convert a snake_case field name to a human-readable label.

    Uses ABBREVIATION_MAP for known terms, title-cases the rest.

    Examples:
        "unique_seqs_count" -> "Unique Sequences Count"
        "lbi" -> "LBI"
        "mean_mut_freq" -> "Mean Mutation Frequency"
        "mean_surprise_mutsel" -> "Mean Surprise MutSel"
    """
    parts = field_name.split("_")
    result = []
    for part in parts:
        lower = part.lower()
        if lower in ABBREVIATION_MAP:
            result.append(ABBREVIATION_MAP[lower])
        else:
            result.append(part.capitalize())
    return " ".join(result)


def sample_values(dicts: List[Dict], field: str, max_samples: int = MAX_SAMPLE_VALUES) -> List[Any]:
    """Sample non-null values for a field across a list of dicts."""
    values = []
    for d in dicts:
        if field in d and d[field] is not None:
            values.append(d[field])
            if len(values) >= max_samples:
                break
    return values


def sample_values_by_path(dicts: List[Dict], path: str, max_samples: int = MAX_SAMPLE_VALUES) -> List[Any]:
    """Sample non-null values using a dot-path across a list of dicts."""
    values = []
    for d in dicts:
        val = _get_nested_value(d, path)
        if val is not None:
            values.append(val)
            if len(values) >= max_samples:
                break
    return values


def collect_keys(dicts: List[Dict]) -> set:
    """Collect the union of all keys across a list of dicts."""
    keys = set()
    for d in dicts:
        keys.update(d.keys())
    return keys


def _apply_custom_fields(metadata, custom_fields, level, data_dicts=None):
    """
    Apply custom field declarations to a metadata dict for a given level.

    Handles output_name renaming and path resolution. If a custom field
    specifies a path, values are resolved via dot-path into the data dicts.

    Args:
        metadata: The field metadata dict to update (modified in place).
        custom_fields: List of custom field declarations.
        level: The level to filter on ("clone", "node", etc.).
        data_dicts: Optional list of data dicts for path resolution.
    """
    if not custom_fields:
        return
    for cf in custom_fields:
        if normalize_level(cf.get("level", "")) != level:
            continue

        # Resolve display mode: skip flag takes precedence, then explicit display
        if cf.get("skip", False) or cf.get("display") == "skip":
            metadata.pop(cf["name"], None)
            output_key = cf.get("output_name")
            if output_key:
                metadata.pop(output_key, None)
            continue

        output_key = cf.get("output_name", cf["name"])
        display = cf.get("display", "dropdown")
        entry = {"type": cf["type"], "display": display, "label": cf["label"]}

        # If path is specified and we have data, verify the field exists
        # via path resolution. If it doesn't exist, still register it
        # (user may be declaring a field they intend to add later).
        path = cf.get("path")
        if path and data_dicts:
            values = sample_values_by_path(data_dicts, path, max_samples=MAX_SAMPLE_PATH)
            # If values found and type is continuous, compute range
            if values and entry["type"] == "continuous" and "range" not in cf:
                numeric_vals = [v for v in values if isinstance(v, (int, float)) and not isinstance(v, bool)]
                if numeric_vals:
                    entry["range"] = [min(numeric_vals), max(numeric_vals)]

        # Preserve range from auto-detection or existing
        if "range" in cf:
            entry["range"] = cf["range"]
        elif output_key in metadata and "range" in metadata[output_key]:
            entry["range"] = metadata[output_key]["range"]
        elif cf["name"] in metadata and "range" in metadata[cf["name"]]:
            entry["range"] = metadata[cf["name"]]["range"]
        # If renaming, remove the original name from metadata
        if output_key != cf["name"] and cf["name"] in metadata:
            del metadata[cf["name"]]
        metadata[output_key] = entry


def _get_nested_value(d: Dict, path: str) -> Any:
    """
    Resolve a simple dot-path against a dict.

    Supports paths like "sample.locus" for nested field access.
    Does NOT support array indexing (e.g., "nodes[].field").
    """
    parts = path.split(".")
    current = d
    for part in parts:
        if isinstance(current, dict) and part in current:
            current = current[part]
        else:
            return None
    return current


def _apply_suggestions(metadata: Dict) -> Dict:
    """Apply SUGGESTED_SKIP_FIELDS and SUGGESTED_DISPLAY_MODES to metadata.

    Removes fields in SUGGESTED_SKIP_FIELDS and applies display overrides
    from SUGGESTED_DISPLAY_MODES. Runs before _apply_custom_fields, so
    user config declarations can re-add suggested-skip fields if needed.

    Note: build_config._is_skip also consults SUGGESTED_SKIP_FIELDS for
    YAML layout (active vs skip section). These are separate concerns —
    this function affects output metadata, _is_skip affects config generation.
    """
    for field_name in list(metadata.keys()):
        if field_name in SUGGESTED_SKIP_FIELDS:
            del metadata[field_name]
        elif field_name in SUGGESTED_DISPLAY_MODES:
            if metadata[field_name].get("display") == "dropdown":
                metadata[field_name]["display"] = SUGGESTED_DISPLAY_MODES[field_name]
    return metadata


def _make_entry(type_str: str, label: str, display: str = "dropdown", **extra) -> Dict:
    """Build a field_metadata entry dict with type, display, label, and optional extras."""
    entry = {"type": type_str, "display": display, "label": label}
    for k, v in extra.items():
        if v is not None:
            entry[k] = v
    return entry


def entry_from_known(known: Dict) -> Dict:
    """Build a field_metadata entry from a known fields registry entry.

    Copies type, display, and label. Omits path (internal routing only).
    Defaults display to 'dropdown' if not specified.
    """
    return {
        "type": known["type"],
        "display": known.get("display", "dropdown"),
        "label": known["label"],
    }


# =============================================================================
# Per-level metadata generators
# =============================================================================


def generate_clone_metadata(
    clones: List[Dict],
    custom_fields: Optional[List[Dict]] = None,
    exclude: Optional[set] = None,
) -> Dict[str, Dict[str, str]]:
    """
    Generate field metadata for clone-level fields.

    Introspects actual clone dicts, matches against known registry,
    infers type for unknown fields, and merges custom declarations.

    Args:
        clones: List of clone dictionaries.
        custom_fields: Optional custom field declarations (level=="clone" only).
        exclude: Optional set of field names to exclude from clone-level
            auto-detection (typically those classified as tree-level by
            :func:`classify_tree_extras`). Custom fields at level=clone
            still take effect for these names.

    Returns:
        Dict mapping field_name -> {"type": ..., "label": ...}
    """
    if not clones:
        return {}

    metadata = {}
    all_keys = collect_keys(clones)

    # Check for known fields that use dot-paths (e.g., locus → sample.locus)
    for field_name, field_info in KNOWN_CLONE_FIELDS.items():
        if "path" in field_info and field_name not in all_keys:
            values = sample_values_by_path(clones, field_info["path"], max_samples=MAX_SAMPLE_PATH)
            if values:
                all_keys.add(field_name)

    # Filter out excluded fields
    candidate_keys = all_keys - EXCLUDED_CLONE_FIELDS
    if exclude:
        candidate_keys -= set(exclude)

    for key in sorted(candidate_keys):
        if key in KNOWN_CLONE_FIELDS:
            known = KNOWN_CLONE_FIELDS[key]
            # Use path if specified, otherwise top-level lookup
            if "path" in known:
                values = sample_values_by_path(clones, known["path"])
            else:
                values = sample_values(clones, key)
            if values:
                metadata[key] = entry_from_known(known)
        else:
            values = sample_values(clones, key)
            if values:
                field_type = infer_field_type(values)
                metadata[key] = _make_entry(field_type, humanize_label(key))

    _apply_suggestions(metadata)
    _apply_custom_fields(metadata, custom_fields, "clone", clones)

    return metadata


def _hashable_for_distinct(v: Any) -> Any:
    """Return a hashable token suitable for set-membership distinctness checks.

    Tree-level field values are typically scalar, but lists/dicts can
    leak through (e.g. JSON-encoded extras). Wrap unhashable types in a
    repr so they still compare distinctly.
    """
    try:
        hash(v)
        return v
    except TypeError:
        return ("__unhashable__", repr(v))


_MISSING = object()


def classify_tree_extras(clones: List[Dict]) -> set:
    """Identify fields that vary across trees within at least one clone.

    Walks each clone's ``trees[]`` list. A field is tree-level for the
    whole dataset if any clone has two or more trees that disagree on
    its value (NaN/missing-as-distinct).

    Single-tree-per-clone datasets trivially classify everything as
    clone-level (no intra-clone variance possible).

    Args:
        clones: List of clone dictionaries with ``trees[]`` populated
            with trimmed tree-ref dicts.

    Returns:
        Set of field names that should be treated as tree-level.
    """
    tree_level_keys: set = set()

    candidate_keys: set = set()
    for clone in clones:
        for tree_ref in clone.get("trees", []) or []:
            if isinstance(tree_ref, dict):
                candidate_keys.update(tree_ref.keys())
    candidate_keys -= EXCLUDED_TREE_FIELDS

    for clone in clones:
        clone_trees = clone.get("trees", []) or []
        if len(clone_trees) < 2:
            continue
        for key in candidate_keys - tree_level_keys:
            distinct = {
                _hashable_for_distinct(t.get(key, _MISSING))
                for t in clone_trees
                if isinstance(t, dict)
            }
            if len(distinct) > 1:
                tree_level_keys.add(key)

    return tree_level_keys


def flatten_tree_refs(clones: List[Dict]) -> List[Dict]:
    """Flatten ``clone[].trees[]`` (dict refs only) into a single list."""
    refs: List[Dict] = []
    for clone in clones:
        for tree_ref in clone.get("trees", []) or []:
            if isinstance(tree_ref, dict):
                refs.append(tree_ref)
    return refs


def generate_tree_metadata(
    clones: List[Dict],
    custom_fields: Optional[List[Dict]] = None,
    tree_level_keys: Optional[set] = None,
) -> Dict[str, Dict[str, str]]:
    """Generate field metadata for tree-level fields.

    Tree-level fields are those that vary across trees within the same
    clone. They drive the webapp's tree-dropdown color/filter/sort
    controls.

    Args:
        clones: List of clone dicts with ``trees[]`` populated.
        custom_fields: Optional custom field declarations (level=="tree" only).
        tree_level_keys: Pre-computed result of :func:`classify_tree_extras`,
            passed in by ``generate_field_metadata`` to avoid a second
            walk over every tree. When ``None``, the classifier runs
            here.

    Returns:
        Dict mapping field_name -> {"type": ..., "label": ..., "range"?: [...]}.
    """
    tree_refs = flatten_tree_refs(clones)
    if tree_level_keys is None:
        tree_level_keys = classify_tree_extras(clones)

    metadata: Dict[str, Dict[str, Any]] = {}
    for key in sorted(tree_level_keys):
        values = sample_values(tree_refs, key)
        if not values:
            continue
        if key in KNOWN_TREE_FIELDS:
            entry = entry_from_known(KNOWN_TREE_FIELDS[key])
        else:
            field_type = infer_field_type(values)
            entry = _make_entry(field_type, humanize_label(key))
        if entry["type"] == "continuous":
            field_range = compute_range(tree_refs, key)
            if field_range:
                entry["range"] = field_range
        metadata[key] = entry

    _apply_suggestions(metadata)
    _apply_custom_fields(metadata, custom_fields, "tree", tree_refs)

    return metadata


def generate_node_metadata(
    trees: List[Dict],
    custom_fields: Optional[List[Dict]] = None,
) -> Dict[str, Dict[str, str]]:
    """
    Generate field metadata for node-level fields.

    Introspects nodes across all trees.

    Args:
        trees: List of tree dictionaries with "nodes" key.
        custom_fields: Optional custom field declarations (level=="node" only).

    Returns:
        Dict mapping field_name -> {"type": ..., "label": ...}
    """
    all_nodes = collect_nodes(trees)
    if not all_nodes:
        return {}

    # Collect node fields that have been demoted to mutation level via encoding.
    # These are source fields for records-encoded mutations, or field names for
    # list/json-encoded mutations.  Exclude them from node metadata unless the
    # config explicitly includes them at node level (not skipped).
    demoted_sources = set()
    node_level_overrides = set()
    if custom_fields:
        for cf in custom_fields:
            norm = normalize_level(cf.get("level", ""))
            if norm == "mutation" and cf.get("encoding"):
                if cf["encoding"] == "records" and cf.get("source"):
                    demoted_sources.add(cf["source"])
                elif cf["encoding"] in ("list", "json"):
                    demoted_sources.add(cf["name"])
            elif norm == "node" and not cf.get("skip") and cf.get("display") != "skip":
                node_level_overrides.add(cf["name"])
    # Only exclude demoted sources that aren't explicitly re-included at node level
    demoted_exclusions = demoted_sources - node_level_overrides

    metadata = {}
    all_keys = collect_keys(all_nodes) - EXCLUDED_NODE_FIELDS - demoted_exclusions

    for key in sorted(all_keys):
        if key in KNOWN_NODE_FIELDS:
            values = sample_values(all_nodes, key)
            if values:
                metadata[key] = entry_from_known(KNOWN_NODE_FIELDS[key])
        elif key in KNOWN_BRANCH_FIELDS:
            continue
        else:
            values = sample_values(all_nodes, key)
            if values:
                field_type = infer_field_type(values)
                metadata[key] = _make_entry(field_type, humanize_label(key))

    _apply_suggestions(metadata)
    _apply_custom_fields(metadata, custom_fields, "node", all_nodes)

    return metadata


def generate_branch_metadata(
    trees: List[Dict],
    custom_fields: Optional[List[Dict]] = None,
) -> Dict[str, Dict[str, str]]:
    """
    Generate field metadata for branch-level fields.

    Branch fields are stored on node dicts (length, branch_length).

    Args:
        trees: List of tree dictionaries with "nodes" key.
        custom_fields: Optional custom field declarations (level=="branch" only).

    Returns:
        Dict mapping field_name -> {"type": ..., "label": ...}
    """
    all_nodes = collect_nodes(trees)
    if not all_nodes:
        return {}

    metadata = {}
    all_keys = collect_keys(all_nodes) - EXCLUDED_BRANCH_FIELDS

    for key in sorted(all_keys):
        if key in KNOWN_BRANCH_FIELDS:
            values = sample_values(all_nodes, key)
            if values:
                metadata[key] = entry_from_known(KNOWN_BRANCH_FIELDS[key])

    _apply_suggestions(metadata)
    _apply_custom_fields(metadata, custom_fields, "branch", all_nodes)

    return metadata


def generate_mutation_metadata(
    trees: List[Dict],
    custom_fields: Optional[List[Dict]] = None,
) -> Dict[str, Dict[str, str]]:
    """
    Generate field metadata for mutation-level fields.

    Mutation fields come from mutations arrays on nodes.
    Continuous fields include a "range" key with [min, max] for color scale domains.

    Args:
        trees: List of tree dictionaries with "nodes" key.
        custom_fields: Optional custom field declarations (level=="mutation" only).

    Returns:
        Dict mapping field_name -> {"type": ..., "label": ..., "range"?: [...]}
    """
    all_mutations = collect_mutations(trees)

    # Check if nodes have AA sequence data — if so, the web app will derive
    # per-mutation child_aa/parent_aa fields during alignment rendering,
    # even if no mutations arrays exist in the data.
    all_nodes = collect_nodes(trees, max_nodes=MAX_NODES_QUICK)
    has_aa_sequences = any(
        n.get("sequence_alignment_aa") for n in all_nodes if isinstance(n, dict)
    )

    if not all_mutations:
        # No pre-computed mutation data; declare derived fields if sequences exist
        metadata = {}
        if has_aa_sequences:
            metadata["child_aa"] = _make_entry("aa", "Child Amino Acid")
            metadata["parent_aa"] = _make_entry("aa", "Parent Amino Acid", display="tooltip")
        _apply_suggestions(metadata)
        _apply_custom_fields(metadata, custom_fields, "mutation", all_nodes)
        return metadata

    metadata = {}
    all_keys = collect_keys(all_mutations) - EXCLUDED_MUTATION_FIELDS

    for key in sorted(all_keys):
        if key in KNOWN_MUTATION_FIELDS:
            values = sample_values(all_mutations, key)
            if values:
                entry = entry_from_known(KNOWN_MUTATION_FIELDS[key])
                if entry["type"] == "continuous":
                    field_range = compute_range(all_mutations, key)
                    if field_range:
                        entry["range"] = field_range
                metadata[key] = entry
        else:
            values = sample_values(all_mutations, key)
            if values:
                field_type = infer_field_type(values)
                entry = _make_entry(field_type, humanize_label(key))
                if field_type == "continuous":
                    field_range = compute_range(all_mutations, key)
                    if field_range:
                        entry["range"] = field_range
                metadata[key] = entry

    _apply_suggestions(metadata)
    _apply_custom_fields(metadata, custom_fields, "mutation", all_mutations)

    return metadata


# =============================================================================
# Helper functions for collecting data across trees
# =============================================================================


def collect_nodes(trees: List[Dict], max_nodes: int = MAX_NODES_SAMPLE) -> List[Dict]:
    """Collect node dicts from across trees (up to max_nodes for sampling)."""
    nodes = []
    for tree in trees:
        tree_nodes = tree.get("nodes", [])
        if isinstance(tree_nodes, dict):
            tree_nodes = list(tree_nodes.values())
        for node in tree_nodes:
            if isinstance(node, dict):
                nodes.append(node)
                if len(nodes) >= max_nodes:
                    return nodes
    return nodes


def collect_mutations(trees: List[Dict], max_mutations: int = MAX_MUTATIONS) -> List[Dict]:
    """Collect mutation dicts from mutations arrays across tree nodes."""
    mutations = []
    for tree in trees:
        tree_nodes = tree.get("nodes", [])
        if isinstance(tree_nodes, dict):
            tree_nodes = list(tree_nodes.values())
        for node in tree_nodes:
            if isinstance(node, dict):
                node_mutations = node.get("mutations")
                if isinstance(node_mutations, list):
                    for mut in node_mutations:
                        if isinstance(mut, dict):
                            mutations.append(mut)
                            if len(mutations) >= max_mutations:
                                return mutations
    return mutations


# =============================================================================
# Top-level generation function
# =============================================================================


def generate_field_metadata(
    clones: List[Dict],
    trees: Optional[List[Dict]] = None,
    custom_fields: Optional[List[Dict]] = None,
) -> Dict[str, Dict[str, Dict[str, str]]]:
    """
    Generate complete field_metadata for a dataset.

    Introspects clones and trees to discover available fields at each level,
    classifies them by type, and generates human-readable labels. Custom field
    declarations from YAML config override auto-introspection.

    Args:
        clones: List of clone dictionaries for the dataset.
        trees: List of tree dictionaries for the dataset (optional).
        custom_fields: List of custom field declaration dicts, each with
            keys: name, level, type, label, and optionally path.

    Returns:
        Dict with level keys ("clone", "node", "branch", "mutation"),
        each containing a dict of field_name -> {"type": ..., "label": ...}.
        Levels with no fields are omitted.
    """
    if trees is None:
        trees = []

    result = {}

    # Compute tree-level keys first so they can be excluded from
    # clone-level auto-detection. Single-tree-per-clone datasets get an
    # empty set here and classify everything as clone-level.
    tree_level_keys = classify_tree_extras(clones)

    clone_meta = generate_clone_metadata(
        clones, custom_fields, exclude=tree_level_keys
    )
    if clone_meta:
        result["clone"] = clone_meta

    tree_meta = generate_tree_metadata(
        clones, custom_fields, tree_level_keys=tree_level_keys
    )
    if tree_meta:
        result["tree"] = tree_meta

    node_meta = generate_node_metadata(trees, custom_fields)
    if node_meta:
        result["node"] = node_meta

    branch_meta = generate_branch_metadata(trees, custom_fields)
    if branch_meta:
        result["branch"] = branch_meta

    mutation_meta = generate_mutation_metadata(trees, custom_fields)
    if mutation_meta:
        result["mutation"] = mutation_meta

    return result
