#!/usr/bin/env python
"""Ingest AIRR input (the AIRR-C v2 Clone/Tree/Node/Cell schema, AIRR Schema
v2.0.0: https://docs.airr-community.org/en/latest/datarep/clone.html) into
Olmsted format.

This is the format Dowser's ``writeTreesJSON`` emits (see issue #36 and the
Apr–Jun 2026 AIRR-C standards thread that converged on naming it this way).
olmsted-cli's earlier, Olmsted-flavored ``-f airr`` container predated this
and never corresponded to an official AIRR release at any version (see
issue #47) — it was removed once this format became the sole ``-f airr``.

This implementation is a pragmatic mapping of the schema's concepts onto
Olmsted's own JSON shape, not a byte-for-byte implementation of it — e.g.
``Clone.nodes`` here is a list, where the official schema's ``Tree.nodes``
is a dict keyed by ``sequence_id``:

    {
      "Clone":        [ {clone_id, clone_class, tree(newick), nodes:[...], ...}, ... ],
      "Rearrangement":[ {sequence_id | cell_id, sequence_alignment, locus, ...}, ... ]
    }

Each ``Clone`` carries its tree inline as a Newick string plus a ``nodes`` list.
Nodes do **not** carry sequences directly; they point into the separate
``Rearrangement`` table by ``sequence_id`` (``clone_class: Rearrangement``) or
``cell_id`` (``clone_class: Cell``). Inferred internal/germline nodes get their
own synthesized ``Rearrangement`` records (e.g. ``sequence_id: Germline-10004``),
so every node — observed or ASR-inferred — joins to a sequence.

Structurally this is far closer to PCP than to legacy AIRR (synthesize a dataset
from flat records + build trees), so the assembly mirrors
``process_pcp_data.process_pcp_to_olmsted`` and returns the same
``([dataset], clones_dict, trees)`` contract.

Chain handling:

- ``clone_class: Rearrangement`` (H only) → one clone + one tree.
- ``clone_class: Cell``, single locus (H only) → one clone + one tree.
- ``clone_class: Cell``, paired H+L → **two** clones + two trees sharing one
  topology, suffixed ``-heavy`` / ``-light``, each node's sequence taken from the
  same-locus ``Rearrangement`` (IGH → heavy, IGK/IGL → light). This mirrors the
  PCP-paired output model (the webapp treats heavy and light as separate clones).

Dowser's ``info`` catchall (present when the input is written with the
default ``dowser_fields=TRUE``, absent from the "clean v2" ``noinfo`` shape)
is read for per-node ``collapse_count`` -> ``multiplicity`` (see
``_node_multiplicity``, #45), clone-level ``v_call``/``j_call``/
``junction_length`` when present, and per-position gene-region labels
(``Clone.info.region``) -> ``cdr1``/``cdr2``/``cdr3`` ``_alignment_start``/
``_end``/``_length`` (see ``_cdr_boundaries_from_region``, #45) — not yet
handled for paired heavy+light clones, where ``region`` covers both chains
concatenated. Still deferred (#45): ``program_origin`` and arbitrary Dowser
``trait=`` columns in per-node tipdata.

Clone-level ``v_call``/``j_call`` fall back to the germline/root node's own
``Rearrangement`` record when the ``Clone`` doesn't supply them (see
``_build_nodes``'s ``germline_gene_calls``, cf. #24); ``d_call`` — which the
``Clone``-level ``info`` catchall never carries at all — comes *only* from
that fallback.

Also deferred (see issue #36): streaming.
"""

from collections import defaultdict
from typing import Any, Dict, List, Optional, Tuple

import ete3

from .identifier import IdentMinter
from .metrics import compute_mean_mut_freq, compute_tree_metrics
from .process_utils import assign_branch_lengths, tag_field_metadata
from .schemas import SCHEMA_VERSION
from .utils import set_verbosity, translate_dna_to_aa, vprint

#: Loci that map to the heavy / light chain split.
HEAVY_LOCI = {"IGH"}
LIGHT_LOCI = {"IGK", "IGL"}


def _locus_chain(locus: Optional[str]) -> Optional[str]:
    """Map an AIRR ``locus`` (e.g. ``"IGH"``) to ``"heavy"``/``"light"``/None."""
    if not locus:
        return None
    locus = locus.upper()
    if locus in HEAVY_LOCI:
        return "heavy"
    if locus in LIGHT_LOCI:
        return "light"
    return None


def index_rearrangements(
    rearrangements: List[Dict[str, Any]],
) -> Tuple[Dict[str, Dict[str, Any]], Dict[str, List[Dict[str, Any]]]]:
    """Index the ``Rearrangement`` table for node → sequence joins.

    Returns ``(by_sequence_id, by_cell_id)`` where ``by_cell_id`` maps each
    ``cell_id`` to the list of its records (paired cells carry one per locus).
    """
    by_sequence_id: Dict[str, Dict[str, Any]] = {}
    by_cell_id: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
    for record in rearrangements:
        seq_id = record.get("sequence_id")
        if seq_id:
            by_sequence_id[seq_id] = record
        cell_id = record.get("cell_id")
        if cell_id:
            by_cell_id[cell_id].append(record)
    return by_sequence_id, by_cell_id


def _reroot_on_ancestor(
    ete_tree: "ete3.TreeNode", ancestor_name: Optional[str]
) -> "ete3.TreeNode":
    """Reroot ``ete_tree`` so ``ancestor_name`` (the germline) is the root.

    Dowser writes the germline as a zero-length outgroup tip hanging off the
    Newick's internal root; Olmsted's convention (like PCP and the legacy AIRR
    path) puts the naive/germline at the root. When the ancestor can't be found,
    the tree is returned unrerooted. Mirrors ``process_airr_data.reroot_tree``
    but keyed on the explicit ``inferred_ancestor`` field rather than a regex.
    """
    if not ancestor_name:
        return ete_tree
    matches = ete_tree.search_nodes(name=ancestor_name)
    if not matches:
        return ete_tree
    node = matches[0]
    if ete_tree == node:
        return ete_tree
    ete_tree.set_outgroup(node)
    ete_tree.remove_child(node)
    node.add_child(ete_tree)
    node.dist = 0
    return node


def _resolve_rearrangement_record(
    node_record: Dict[str, Any],
    clone_class: str,
    chain: Optional[str],
    by_sequence_id: Dict[str, Dict[str, Any]],
    by_cell_id: Dict[str, List[Dict[str, Any]]],
) -> Optional[Dict[str, Any]]:
    """Join one tree node to its full ``Rearrangement`` record for a given chain.

    - ``Rearrangement`` class: join by ``sequence_id``.
    - ``Cell`` class: join by ``cell_id``; for a paired split, pick the record
      whose locus matches ``chain`` (heavy → IGH, light → IGK/IGL). For a
      single-chain cell (``chain is None``) take the sole record.

    Returns ``None`` when no matching record exists.
    """
    if clone_class == "Cell":
        records = by_cell_id.get(node_record.get("cell_id"), [])
        if chain is not None:
            records = [r for r in records if _locus_chain(r.get("locus")) == chain]
        return records[0] if records else None
    return by_sequence_id.get(node_record.get("sequence_id"))


def _clone_chains(
    clone: Dict[str, Any],
    by_sequence_id: Dict[str, Dict[str, Any]],
    by_cell_id: Dict[str, List[Dict[str, Any]]],
) -> List[Optional[str]]:
    """Decide which chain(s) to emit for a clone.

    ``[None]`` for single-chain data (one clone, no suffix); ``["heavy",
    "light"]`` when a ``Cell`` clone spans both a heavy and a light locus.
    """
    clone_class = clone.get("clone_class")
    if clone_class != "Cell":
        return [None]

    chains = set()
    for node in clone.get("nodes", []):
        for record in by_cell_id.get(node.get("cell_id"), []):
            chain = _locus_chain(record.get("locus"))
            if chain:
                chains.add(chain)
    if "heavy" in chains and "light" in chains:
        return ["heavy", "light"]
    return [None]


def _tree_ref(
    *,
    tree_ident: str,
    clone_id: str,
    tree_id: str,
    chain: Optional[str],
    newick: str,
) -> Dict[str, Any]:
    """Build a ``clone.trees[]`` reference / top-level tree header.

    Applies the ``-heavy`` / ``-light`` suffix to ident/clone_id/tree_id for
    paired data, matching ``process_pcp_data._build_tree_ref``.
    """
    suffix = f"-{chain}" if chain is not None else ""
    return {
        "ident": f"{tree_ident}{suffix}",
        "clone_id": f"{clone_id}{suffix}",
        "tree_id": f"{tree_id}{suffix}",
        "tree_name": f"{tree_id}{suffix}",
        "newick": newick,
    }


def _node_multiplicity(node_record: Dict[str, Any]) -> Optional[int]:
    """Extract a node's real multiplicity from Dowser's ``info`` catchall (#45).

    Present only when the input was written with ``dowser_fields=TRUE``
    (Dowser's default): observed nodes carry
    ``info: {"tipdata": {"collapse_count": N, ...}}``. Inferred/ASR nodes
    (and every node in the ``noinfo`` variant) have no usable ``collapse_count``
    — ``info`` is either absent, an empty list, or a dict without a
    ``tipdata`` entry — and this returns ``None`` for those, matching the
    existing "leave unset rather than fabricate" convention.
    """
    info = node_record.get("info")
    if not isinstance(info, dict):
        return None
    tipdata = info.get("tipdata")
    if not isinstance(tipdata, dict):
        return None
    return tipdata.get("collapse_count")


#: Region labels from Dowser's Clone.info.region we surface as cdr*_alignment
#: fields (matching PCP/legacy-AIRR's field set). fwr1-4 have no established
#: Olmsted output field, so they're not derived here.
_CDR_REGIONS = ("cdr1", "cdr2", "cdr3")

#: IMGT/AIRR Community: junction = CDR3 plus the 2 conserved anchor residues
#: CDR3 excludes -- the V-gene's 2nd-CYS (start) and the J-gene's TRP/PHE
#: (end) -- so junction is exactly 1 codon (3 nucleotides = 1 amino acid)
#: longer on each side (#46). Verified exactly against real data: a 51nt
#: region-derived (strict) CDR3 span vs. a 57nt junction_length for the same
#: clone, 51 + 2*3 = 57.
_JUNCTION_ANCHOR_LENGTH = 3


def _cdr_boundaries_from_region(
    region: Optional[List[str]], germline_alignment: Optional[str]
) -> Dict[str, Tuple[int, int]]:
    """Derive cdr1/cdr2/cdr3 boundaries in ``germline_alignment``'s (gapped)
    coordinate space from Dowser's ``Clone.info.region`` catchall (#45).

    ``region`` is a per-position IMGT region label array (``fwr1``, ``cdr1``,
    ``fwr2``, ``cdr2``, ``fwr3``, ``cdr3``, ``fwr4``) indexed against the
    *ungapped* ``Rearrangement.sequence`` — one entry per nucleotide, no gap
    placeholders. Olmsted's ``germline_alignment`` (what everything renders
    against) is the *gapped* ``sequence_alignment`` — IMGT ``.`` padding for
    CDR-loop-length harmonization. This walks ``germline_alignment``'s
    non-gap characters in lockstep with ``region`` to remap each boundary
    into gapped coordinates; a run of gap characters is attributed to
    whichever region it falls within, so within-region IMGT padding doesn't
    leak into a neighboring region's span.

    The ``cdr3`` span is extended by ``_JUNCTION_ANCHOR_LENGTH`` ungapped
    nucleotides on each side (i.e. before remapping to gapped coordinates, so
    the extension can't land inside a run of IMGT padding) to convert
    ``region``'s strict IMGT CDR3 into the junction convention ``cdr3_length``
    uses everywhere else in this project (#46) — ``cdr1``/``cdr2`` have no
    such distinction and are returned as ``region`` gives them.

    Returns ``{"cdr1": (start, end), ...}`` (0-based, half-open nucleotide
    positions — the same convention as the PCP/legacy-AIRR
    ``cdr*_alignment_start``/``_end`` fields) for whichever of cdr1/cdr2/cdr3
    are present in ``region``. Returns ``{}`` when ``region`` or
    ``germline_alignment`` is missing/empty, or when the alignment's non-gap
    character count doesn't match ``len(region)`` (defensive: a length
    mismatch means the position mapping can't be trusted, so skip rather
    than guess).

    Known gap (not yet handled): for a paired (heavy+light) ``Cell`` clone,
    Dowser's ``region`` is a single array covering *both* chains
    concatenated (verified: heavy chain's positions, then light chain's
    restarting at ``fwr1``) — this doesn't match either chain's individual
    ``germline_alignment`` alone, so the length check above safely skips it
    rather than misattributing one chain's boundaries to the other. Splitting
    the combined array per-locus would need a confirmed, general concatenation
    order from Dowser (only one example was available to infer heavy-then-
    light from) — deferred rather than guessed.
    """
    if not region or not germline_alignment:
        return {}

    gapped_positions = [
        i for i, ch in enumerate(germline_alignment) if ch not in (".", "-")
    ]
    if len(gapped_positions) != len(region):
        return {}

    def gapped_boundary(ungapped_idx: int) -> int:
        return (
            gapped_positions[ungapped_idx]
            if ungapped_idx < len(gapped_positions)
            else len(germline_alignment)
        )

    boundaries: Dict[str, Tuple[int, int]] = {}
    start = 0
    current = region[0]
    for i in range(1, len(region) + 1):
        label = region[i] if i < len(region) else None
        if label != current:
            if current in _CDR_REGIONS:
                span_start, span_end = start, i
                if current == "cdr3":
                    span_start = max(0, span_start - _JUNCTION_ANCHOR_LENGTH)
                    span_end = min(len(region), span_end + _JUNCTION_ANCHOR_LENGTH)
                boundaries[current] = (
                    gapped_boundary(span_start),
                    gapped_boundary(span_end),
                )
            start = i
            current = label
    return boundaries


def _build_nodes(
    clone: Dict[str, Any],
    chain: Optional[str],
    by_sequence_id: Dict[str, Dict[str, Any]],
    by_cell_id: Dict[str, List[Dict[str, Any]]],
    *,
    compute_metrics: bool,
    lbi_tau: float,
) -> Tuple[List[Dict[str, Any]], str, Optional[str], Dict[str, str]]:
    """Build the Olmsted node array for one (clone, chain).

    Topology comes from the ``Clone.tree`` Newick (authoritative; its labels are
    the node ids), rerooted on ``inferred_ancestor``. Returns
    ``(nodes, root_id, germline_alignment, germline_gene_calls)``, or
    ``([], "", None, {})`` when the Newick is missing/unparseable.

    ``germline_gene_calls`` is ``{"v_call": ..., "d_call": ..., "j_call": ...}``
    (only keys with a non-``None`` value) read off the germline/root node's own
    ``Rearrangement`` record (#45/#24) — V(D)J gene usage is invariant across a
    clone's members by definition (one recombination event founds the clone),
    so the germline record is a clean, unbiased single source, matching the
    same convention already used for ``germline_alignment`` (also sourced from
    the root node). The caller uses this as a fallback for clone-level
    ``v_call``/``j_call`` (when the ``Clone`` doesn't supply them) and as the
    only source for ``d_call`` (which the ``Clone``-level ``info`` catchall
    never carries).
    """
    newick = clone.get("tree")
    if not newick:
        vprint.status(
            f"  Warning: clone {clone.get('clone_id')} has no tree; skipping."
        )
        return [], "", None, {}

    ete_tree = ete3.PhyloTree(newick, format=1)
    ete_tree = _reroot_on_ancestor(ete_tree, clone.get("inferred_ancestor"))

    clone_class = clone.get("clone_class", "Rearrangement")
    node_records = {n.get("node_id"): n for n in clone.get("nodes", [])}

    nodes_by_id: Dict[str, Dict[str, Any]] = {}
    edges: List[Tuple[str, str, float]] = []
    root_id = ete_tree.name
    germline_gene_calls: Dict[str, str] = {}

    for ete_node in ete_tree.traverse("postorder"):
        name = ete_node.name
        node_record = node_records.get(name, {})

        if ete_node.up is None:
            node_type_topo = "root"
            parent = None
        else:
            node_type_topo = "leaf" if ete_node.is_leaf() else "internal"
            parent = ete_node.up.name
            edges.append((ete_node.up.name, name, ete_node.dist))

        rearrangement_record = _resolve_rearrangement_record(
            node_record, clone_class, chain, by_sequence_id, by_cell_id
        )
        sequence_alignment = (
            (rearrangement_record.get("sequence_alignment") or "")
            if rearrangement_record
            else ""
        )
        locus = rearrangement_record.get("locus") if rearrangement_record else None
        if not sequence_alignment:
            vprint.status(
                f"  Warning: node '{name}' in clone {clone.get('clone_id')} "
                "has no matching Rearrangement sequence."
            )
        if ete_node.up is None and rearrangement_record:
            germline_gene_calls = {
                field: rearrangement_record[field]
                for field in ("v_call", "d_call", "j_call")
                if rearrangement_record.get(field) is not None
            }

        nodes_by_id[name] = {
            "sequence_id": name,
            "node_id": name,
            "sequence_alignment": sequence_alignment,
            "sequence_alignment_aa": translate_dna_to_aa(sequence_alignment),
            "type": node_type_topo,
            # observed / inferred — the v2 axis distinguishing measured
            # sequences from ASR-reconstructed ancestors.
            "node_type": node_record.get("node_type"),
            "node_class": node_record.get("node_class"),
            "locus": locus,
            "parent": parent,
            # Real multiplicity when the input carries Dowser's info catchall
            # (#45); None for the noinfo schema / inferred nodes, matching the
            # existing "leave unset rather than fabricate" convention (webapp
            # renders "<unspecified>").
            "multiplicity": _node_multiplicity(node_record),
            "timepoint_multiplicities": [],
            "lbi": None,
            "lbr": None,
            "affinity": None,
            "scaled_affinity": None,
        }

    # Branch length / distance from the (rerooted) Newick — same helper the
    # legacy AIRR and merge paths use. Keyed by sequence_id == Newick label.
    assign_branch_lengths(ete_tree, nodes_by_id, overwrite=True)

    if compute_metrics and root_id in nodes_by_id:
        compute_tree_metrics(nodes_by_id, edges, root_id, tau=lbi_tau)

    root_node = nodes_by_id.get(root_id, {})
    germline_alignment = root_node.get("sequence_alignment") or None

    # Emit nodes in the clone's original node order for stable output.
    ordered_ids = [n.get("node_id") for n in clone.get("nodes", [])]
    nodes = [nodes_by_id[nid] for nid in ordered_ids if nid in nodes_by_id]
    # Include any Newick nodes absent from the nodes list (defensive).
    nodes.extend(nodes_by_id[nid] for nid in nodes_by_id if nid not in set(ordered_ids))

    return nodes, root_id, germline_alignment, germline_gene_calls


def _mean_mutation_frequency(
    nodes: List[Dict[str, Any]], germline_alignment: Optional[str]
) -> float:
    """Mean per-site SHM frequency of observed leaves vs the germline root.

    Thin wrapper around the shared :func:`compute_mean_mut_freq` (the single
    source of truth across PCP and AIRR). Uses each node's real
    multiplicity when the input carries Dowser's info catchall (#45,
    ``collapse_count``); falls back to multiplicity=1 (each observed leaf
    counts once) for the clean (``noinfo``) shape, which carries no
    per-node multiplicity at all.
    """
    weighted_nodes = (
        {
            **node,
            "multiplicity": node["multiplicity"]
            if node.get("multiplicity") is not None
            else 1,
        }
        for node in nodes
    )
    mean_mut_freq, _, _ = compute_mean_mut_freq(
        germline_alignment or "", weighted_nodes
    )
    return mean_mut_freq


def _rerooted_newick(clone: Dict[str, Any]) -> str:
    """The clone's Newick rerooted on ``inferred_ancestor``, as a string.

    Emitted on tree/clone records so the Newick agrees with the ``parent`` links
    built from the same rerooted topology.
    """
    ete_tree = ete3.PhyloTree(clone["tree"], format=1)
    ete_tree = _reroot_on_ancestor(ete_tree, clone.get("inferred_ancestor"))
    return ete_tree.write(format=1, format_root_node=True)


def process_airr_to_olmsted(
    clone_records: List[Dict[str, Any]],
    rearrangement_records: List[Dict[str, Any]],
    minter: Optional[IdentMinter] = None,
    name: Optional[str] = None,
    compute_metrics: bool = False,
    lbi_tau: float = 0.0125,
    verbosity: int = 1,
    custom_fields: Optional[List[Dict[str, Any]]] = None,
) -> Tuple[List[Dict[str, Any]], Dict[str, List[Dict[str, Any]]], List[Dict[str, Any]]]:
    """Convert AIRR-C v2 Clone/Tree records to Olmsted format.

    Returns ``([dataset], clones_dict, trees)`` — the same contract as
    ``process_pcp_to_olmsted`` — so the shared assembly, ``tag_field_metadata``,
    validation, and output-writing paths are reused unchanged.
    """
    set_verbosity(verbosity)

    if minter is None:
        minter = IdentMinter()

    by_sequence_id, by_cell_id = index_rearrangements(rearrangement_records)

    # Two dataset mints (semantic foreign key + internal ident), matching PCP.
    dataset_id = minter.mint("dataset")
    dataset_ident = minter.mint("dataset")

    clones: List[Dict[str, Any]] = []
    trees: List[Dict[str, Any]] = []
    # One dataset-level sample per repertoire_id (sample_id must be unique within
    # dataset.samples). Per-clone locus lives on the clone's own denormalized
    # ``sample`` dict below — matching the PCP-paired model, where heavy/light
    # clones share a sample_id but carry different clone.sample.locus values.
    dataset_samples_by_id: Dict[Optional[str], Dict[str, Any]] = {}

    for clone in clone_records:
        clone_id = str(clone.get("clone_id"))
        repertoire_id = clone.get("repertoire_id")
        clone_class = clone.get("clone_class", "Rearrangement")
        clone_ident = minter.mint("clone")
        tree_ident = minter.mint("tree")
        chains = _clone_chains(clone, by_sequence_id, by_cell_id)
        is_paired = chains == ["heavy", "light"]
        rerooted_newick = _rerooted_newick(clone)

        for chain in chains:
            nodes, _root_id, germline_alignment, germline_gene_calls = _build_nodes(
                clone,
                chain,
                by_sequence_id,
                by_cell_id,
                compute_metrics=compute_metrics,
                lbi_tau=lbi_tau,
            )
            if not nodes:
                continue

            suffix = f"-{chain}" if chain is not None else ""
            # Locus for this chain: the dominant locus among its nodes.
            locus = next((n["locus"] for n in nodes if n.get("locus")), None)
            olmsted_locus = locus.lower() if locus else None

            # Register one dataset-level sample per repertoire_id (first locus
            # seen wins its locus label); the per-clone sample below carries the
            # chain-specific locus the webapp reads.
            if repertoire_id not in dataset_samples_by_id:
                dataset_samples_by_id[repertoire_id] = {
                    "ident": minter.mint("sample"),
                    "sample_id": repertoire_id,
                    "locus": olmsted_locus,
                }
            sample = {
                "ident": f"{clone_ident}{suffix}",
                "sample_id": repertoire_id,
                "locus": olmsted_locus,
            }

            tree_ref = _tree_ref(
                tree_ident=tree_ident,
                clone_id=clone_id,
                tree_id=f"tree-{clone_id}",
                chain=chain,
                newick=rerooted_newick,
            )

            clone_out: Dict[str, Any] = {
                "clone_id": f"{clone_id}{suffix}",
                "ident": f"{clone_ident}{suffix}",
                "dataset_id": dataset_id,
                "clone_class": clone_class,
                "sample_id": repertoire_id,
                "unique_seqs_count": len(nodes),
                "mean_mut_freq": _mean_mutation_frequency(nodes, germline_alignment),
                "germline_alignment": germline_alignment,
                "trees": [tree_ref],
                "sample": sample,
            }
            if clone.get("clone_count") is not None:
                clone_out["total_read_count"] = clone.get("clone_count")
            # Pass through clone-level immunological fields only when the input
            # supplies them (the Dowser "info" variant does; the clean v2
            # "noinfo" variant does not — those are left unset, not fabricated).
            for src, dst in (
                ("v_call", "v_call"),
                ("j_call", "j_call"),
                ("junction_length", "cdr3_length"),
            ):
                if clone.get(src) is not None:
                    clone_out[dst] = clone[src]
            # Fall back to the germline Rearrangement record's own v_call/
            # d_call/j_call for whichever of those the Clone-level info
            # catchall didn't supply (cf. #24/#45) — this is the *only*
            # source for d_call, which the Clone-level info never carries.
            for field in ("v_call", "d_call", "j_call"):
                if clone_out.get(field) is None and field in germline_gene_calls:
                    clone_out[field] = germline_gene_calls[field]
            # Derive cdr1/cdr2/cdr3 alignment boundaries from Dowser's
            # Clone.info.region when present (#45) — the only source of
            # CDR1/CDR2 boundaries this format has at all, and a more
            # complete cdr3 (alignment_start/end, not just a bare length).
            #
            # cdr3 is normalized to the junction convention (#46):
            # _cdr_boundaries_from_region already extends region's strict
            # IMGT cdr3 span by the 2 conserved anchor codons junction
            # includes, since cdr3_length is a synonym for junction length
            # everywhere else in this project (schemas.py, PCP, legacy
            # AIRR). So when junction_length is also available, the two
            # should now agree exactly; disagreement means something
            # unexpected (e.g. a non-standard anchor convention) rather than
            # the known, already-accounted-for CDR3/junction difference —
            # worth a warning rather than silently picking one.
            clone_info = clone.get("info")
            region = clone_info.get("region") if isinstance(clone_info, dict) else None
            region_boundaries = _cdr_boundaries_from_region(region, germline_alignment)
            if "cdr3" in region_boundaries and clone.get("junction_length") is not None:
                region_cdr3_start, region_cdr3_end = region_boundaries["cdr3"]
                region_cdr3_length = region_cdr3_end - region_cdr3_start
                if region_cdr3_length != clone["junction_length"]:
                    vprint.status(
                        f"  Warning: clone '{clone_id}' region-derived cdr3 length "
                        f"({region_cdr3_length}, junction-normalized) disagrees with "
                        f"junction_length ({clone['junction_length']}). Using the "
                        "region-derived value."
                    )
            for cdr, (start, end) in region_boundaries.items():
                clone_out[f"{cdr}_alignment_start"] = start
                clone_out[f"{cdr}_alignment_end"] = end
                clone_out[f"{cdr}_length"] = end - start
            if is_paired:
                clone_out["is_paired"] = True
                clone_out["pair_id"] = f"pair-{clone_id}"

            clones.append(clone_out)
            trees.append({**tree_ref, "nodes": nodes})

    dataset: Dict[str, Any] = {
        "ident": dataset_ident,
        "dataset_id": dataset_id,
        "schema_version": SCHEMA_VERSION,
        "build": {"commit": "airr-import", "time": ""},
        "subjects": [],
        "samples": list(dataset_samples_by_id.values()),
        "seeds": [],
        "clone_count": len(clones),
        "subjects_count": 0,
        "timepoints_count": 0,
    }
    if name:
        dataset["name"] = name

    clones_dict: Dict[str, List[Dict[str, Any]]] = {dataset_id: clones}
    dataset["field_metadata"] = tag_field_metadata(clones, trees, custom_fields)

    return [dataset], clones_dict, trees
