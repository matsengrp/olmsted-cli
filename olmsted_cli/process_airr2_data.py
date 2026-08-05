#!/usr/bin/env python
"""Ingest the AIRR-C v2 Clone & Tree schema ("airr2") into Olmsted format.

This is the format Dowser's ``writeTreesJSON`` emits (see issue #36 and the
Apr–Jun 2026 AIRR-C standards thread). It is structurally distinct from the
legacy Olmsted-flavored ``-f airr`` input handled by ``process_airr_data.py``:

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

Deferred (see issue #36): the Dowser ``info`` catchall, streaming, and deriving
clone-level ``v_call``/``j_call``/``cdr3_length`` when the ``Clone`` omits them.
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


def _resolve_sequence(
    node_record: Dict[str, Any],
    clone_class: str,
    chain: Optional[str],
    by_sequence_id: Dict[str, Dict[str, Any]],
    by_cell_id: Dict[str, List[Dict[str, Any]]],
) -> Tuple[str, Optional[str]]:
    """Join one tree node to its ``Rearrangement`` record for a given chain.

    Returns ``(sequence_alignment, locus)``. ``sequence_alignment`` is ``""``
    when no matching record exists (caller warns); ``locus`` is ``None`` then.

    - ``Rearrangement`` class: join by ``sequence_id``.
    - ``Cell`` class: join by ``cell_id``; for a paired split, pick the record
      whose locus matches ``chain`` (heavy → IGH, light → IGK/IGL). For a
      single-chain cell (``chain is None``) take the sole record.
    """
    if clone_class == "Cell":
        records = by_cell_id.get(node_record.get("cell_id"), [])
        if chain is not None:
            records = [r for r in records if _locus_chain(r.get("locus")) == chain]
        record = records[0] if records else None
    else:
        record = by_sequence_id.get(node_record.get("sequence_id"))

    if record is None:
        return "", None
    return record.get("sequence_alignment") or "", record.get("locus")


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


def _build_nodes(
    clone: Dict[str, Any],
    chain: Optional[str],
    by_sequence_id: Dict[str, Dict[str, Any]],
    by_cell_id: Dict[str, List[Dict[str, Any]]],
    *,
    compute_metrics: bool,
    lbi_tau: float,
) -> Tuple[List[Dict[str, Any]], str, Optional[str]]:
    """Build the Olmsted node array for one (clone, chain).

    Topology comes from the ``Clone.tree`` Newick (authoritative; its labels are
    the node ids), rerooted on ``inferred_ancestor``. Returns
    ``(nodes, root_id, germline_alignment)``, or ``([], "", None)`` when the
    Newick is missing/unparseable.
    """
    newick = clone.get("tree")
    if not newick:
        vprint.status(
            f"  Warning: clone {clone.get('clone_id')} has no tree; skipping."
        )
        return [], "", None

    ete_tree = ete3.PhyloTree(newick, format=1)
    ete_tree = _reroot_on_ancestor(ete_tree, clone.get("inferred_ancestor"))

    clone_class = clone.get("clone_class", "Rearrangement")
    node_records = {n.get("node_id"): n for n in clone.get("nodes", [])}

    nodes_by_id: Dict[str, Dict[str, Any]] = {}
    edges: List[Tuple[str, str, float]] = []
    root_id = ete_tree.name

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

        sequence_alignment, locus = _resolve_sequence(
            node_record, clone_class, chain, by_sequence_id, by_cell_id
        )
        if not sequence_alignment:
            vprint.status(
                f"  Warning: node '{name}' in clone {clone.get('clone_id')} "
                "has no matching Rearrangement sequence."
            )

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
            # No multiplicity/timepoint concept in the clean v2 schema; leave
            # unset rather than fabricate (webapp renders "<unspecified>").
            "multiplicity": None,
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

    return nodes, root_id, germline_alignment


def _mean_mutation_frequency(
    nodes: List[Dict[str, Any]], germline_alignment: Optional[str]
) -> float:
    """Mean per-site SHM frequency of observed leaves vs the germline root.

    Thin wrapper around the shared :func:`compute_mean_mut_freq` (the single
    source of truth across PCP, AIRR, and airr2). The clean v2 schema
    carries no per-node multiplicity, so every node is given multiplicity=1
    — each observed leaf counts once, matching this format's inherently
    unweighted convention.
    """
    unit_multiplicity_nodes = ({**node, "multiplicity": 1} for node in nodes)
    mean_mut_freq, _, _ = compute_mean_mut_freq(
        germline_alignment or "", unit_multiplicity_nodes
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


def process_airr2_to_olmsted(
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
            nodes, _root_id, germline_alignment = _build_nodes(
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
            if is_paired:
                clone_out["is_paired"] = True
                clone_out["pair_id"] = f"pair-{clone_id}"

            clones.append(clone_out)
            trees.append({**tree_ref, "nodes": nodes})

    dataset: Dict[str, Any] = {
        "ident": dataset_ident,
        "dataset_id": dataset_id,
        "schema_version": SCHEMA_VERSION,
        "build": {"commit": "airr2-import", "time": ""},
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
