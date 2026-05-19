"""Streaming primitives for memory-bounded ``olmsted process`` (#26, phase 2).

The PCP and AIRR pipelines materialize every parsed tree and per-node dict
before writing. This module supplies the building blocks that phase 3 wires
together so peak memory is bounded by ``--batch-size`` rather than total
dataset size:

- :class:`FieldTypeEvidence` / :class:`RangeEvidence` — counter-based
  type-inference and min/max accumulators safe to merge across batches.
  Replace the sample-capped path used by ``field_metadata.infer_field_type``
  on the streaming side, so a value that contradicts the inferred type
  (e.g. a string arriving after 50 ints) still flips the result correctly.
- :class:`BatchAccumulator` — folds each batch's clones / trees / nodes /
  mutations into per-dataset evidence, ID-uniqueness sets, running totals,
  and the dataset-scope tree-level-keys union.  ``finalize_field_metadata``
  produces the same shape ``generate_field_metadata`` returns today.
- :class:`BatchSpooler` — spools each batch's clones and trees to JSONL
  temp files keyed by dataset, so the final-write step can stream-stitch
  the consolidated JSON without holding more than one batch in memory.
- :func:`write_olmsted_json_streaming` — writes the consolidated output
  with ``metadata`` first (canonical key order) using known totals, then
  streams ``clones`` and ``trees`` from the spooler.  Preserves the
  ``data_io.write_olmsted_json`` gzip determinism guarantee.
"""

from __future__ import annotations

import gzip
import io
import json
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import (
    Any,
    Dict,
    Iterable,
    Iterator,
    List,
    Literal,
    Optional,
    Set,
)

from .constants import (
    AA_CHARS,
    DNA_CHARS,
    EXCLUDED_BRANCH_FIELDS,
    EXCLUDED_CLONE_FIELDS,
    EXCLUDED_MUTATION_FIELDS,
    EXCLUDED_NODE_FIELDS,
    EXCLUDED_TREE_FIELDS,
    FIELD_ALIASES,
    KNOWN_BRANCH_FIELDS,
    KNOWN_CLONE_FIELDS,
    KNOWN_MUTATION_FIELDS,
    KNOWN_NODE_FIELDS,
    KNOWN_TREE_FIELDS,
    normalize_level,
)
from .field_metadata import (
    apply_custom_fields,
    apply_suggestions,
    entry_from_known,
    get_nested_value,
    humanize_label,
    make_entry,
)

JsonOutputFormat = Literal["pretty", "compact", "gzip"]

_LEVELS = ("clone", "tree", "node", "branch", "mutation")
_EXCLUDED_BY_LEVEL = {
    "clone": EXCLUDED_CLONE_FIELDS,
    "tree": EXCLUDED_TREE_FIELDS,
    "node": EXCLUDED_NODE_FIELDS,
    "branch": EXCLUDED_BRANCH_FIELDS,
    "mutation": EXCLUDED_MUTATION_FIELDS,
}
_KNOWN_BY_LEVEL = {
    "clone": KNOWN_CLONE_FIELDS,
    "tree": KNOWN_TREE_FIELDS,
    "node": KNOWN_NODE_FIELDS,
    "branch": KNOWN_BRANCH_FIELDS,
    "mutation": KNOWN_MUTATION_FIELDS,
}


# =============================================================================
# Evidence accumulators
# =============================================================================


@dataclass
class FieldTypeEvidence:
    """Mergeable counters for ``infer_field_type`` semantics.

    Matches the legacy branches in ``field_metadata.infer_field_type``:
    booleans count as strings (so a pure-bool field reports
    ``categorical``); single-character strings narrow to ``dna`` / ``aa``
    only when *every* recorded string value is single-character.
    """

    numeric_count: int = 0
    string_count: int = 0
    bool_count: int = 0
    list_count: int = 0
    dict_count: int = 0
    other_count: int = 0
    has_any_string: bool = False
    all_strings_single_char: bool = True
    string_char_set: Set[str] = field(default_factory=set)

    def record(self, value: Any) -> None:
        if isinstance(value, bool):
            self.bool_count += 1
        elif isinstance(value, (int, float)):
            self.numeric_count += 1
        elif isinstance(value, str):
            self.string_count += 1
            self.has_any_string = True
            if len(value) != 1:
                self.all_strings_single_char = False
            else:
                self.string_char_set.add(value.upper())
        elif isinstance(value, list):
            self.list_count += 1
        elif isinstance(value, dict):
            self.dict_count += 1
        else:
            self.other_count += 1

    def merge(self, other: "FieldTypeEvidence") -> None:
        self.numeric_count += other.numeric_count
        self.string_count += other.string_count
        self.bool_count += other.bool_count
        self.list_count += other.list_count
        self.dict_count += other.dict_count
        self.other_count += other.other_count
        if other.has_any_string:
            self.has_any_string = True
            if not other.all_strings_single_char:
                self.all_strings_single_char = False
            self.string_char_set.update(other.string_char_set)

    def total(self) -> int:
        return (
            self.numeric_count
            + self.string_count
            + self.bool_count
            + self.list_count
            + self.dict_count
            + self.other_count
        )

    def infer(self) -> str:
        total = self.total()
        if total == 0 or self.other_count > 0:
            return "categorical"
        if self.list_count > 0 and self.list_count == total:
            return "list"
        if self.dict_count > 0 and self.dict_count == total:
            return "json"
        if self.list_count > 0 or self.dict_count > 0:
            return "categorical"
        effective_string = self.string_count + self.bool_count
        if self.numeric_count > 0 and effective_string == 0:
            return "continuous"
        if effective_string > 0 and self.numeric_count == 0:
            # Only consult single-char DNA/AA when there are actual strings
            # (booleans count toward effective_string but never feed the
            # alphabet check, mirroring the legacy ``string_values`` path).
            if self.string_count > 0 and self.all_strings_single_char:
                chars = self.string_char_set
                if chars <= AA_CHARS and not chars <= DNA_CHARS:
                    return "aa"
                if chars <= DNA_CHARS:
                    return "dna"
                if chars <= AA_CHARS:
                    return "aa"
            return "categorical"
        return "categorical"


@dataclass
class RangeEvidence:
    """Running ``(min, max, count)`` over the numeric values of one field."""

    min: Optional[float] = None
    max: Optional[float] = None
    count: int = 0

    def record(self, value: Any) -> None:
        if isinstance(value, bool):
            return
        if not isinstance(value, (int, float)):
            return
        if self.count == 0:
            self.min = value
            self.max = value
        else:
            if value < self.min:
                self.min = value
            if value > self.max:
                self.max = value
        self.count += 1

    def merge(self, other: "RangeEvidence") -> None:
        if other.count == 0:
            return
        if self.count == 0:
            self.min = other.min
            self.max = other.max
        else:
            if other.min is not None and other.min < self.min:
                self.min = other.min
            if other.max is not None and other.max > self.max:
                self.max = other.max
        self.count += other.count

    def as_list(self) -> Optional[List[float]]:
        if self.count == 0:
            return None
        return [self.min, self.max]


# =============================================================================
# Per-dataset accumulator
# =============================================================================


@dataclass
class _LevelEvidence:
    """Per-level evidence + the union of field names ever seen."""

    type_evidence: Dict[str, FieldTypeEvidence] = field(default_factory=dict)
    range_evidence: Dict[str, RangeEvidence] = field(default_factory=dict)
    keys: Set[str] = field(default_factory=set)

    def record(self, key: str, value: Any) -> None:
        self.keys.add(key)
        if value is None:
            return
        evidence = self.type_evidence.setdefault(key, FieldTypeEvidence())
        evidence.record(value)
        rng = self.range_evidence.setdefault(key, RangeEvidence())
        rng.record(value)


@dataclass
class _DatasetState:
    levels: Dict[str, _LevelEvidence] = field(
        default_factory=lambda: {lvl: _LevelEvidence() for lvl in _LEVELS}
    )
    tree_level_keys: Set[str] = field(default_factory=set)
    clone_ids: Set[str] = field(default_factory=set)
    tree_ids_by_clone: Dict[str, Set[str]] = field(default_factory=dict)
    sample_ids: Set[str] = field(default_factory=set)
    subject_ids: Set[str] = field(default_factory=set)
    samples: List[Dict[str, Any]] = field(default_factory=list)
    clone_count: int = 0
    tree_count: int = 0
    leaf_count: int = 0
    has_aa_sequences: bool = False
    # PCP-style tree-csv extras: when True, tree-level keys constant
    # across every clone's trees migrate to clone-level in the output
    # JSON (and the matching field_metadata).  AIRR data places these
    # fields on trees natively — the streaming wrapper skips the data
    # hoist AND turns this off so ``_clone_level_metadata`` doesn't
    # absorb tree extras into clone-level metadata.
    hoist_tree_extras_to_clone: bool = True


DuplicateScope = Literal["dataset_id", "clone_id", "tree_id"]


class DuplicateIdError(ValueError):
    """Raised when a ``*_id`` collision is observed and duplicates aren't allowed.

    The ``scope`` attribute identifies which kind of identifier collided so
    callers can branch on it without parsing the string message.
    """

    def __init__(self, message: str, *, scope: DuplicateScope) -> None:
        super().__init__(message)
        self.scope: DuplicateScope = scope


class BatchAccumulator:
    """Folds per-batch clone/tree/node/mutation data into a per-dataset summary.

    The lifecycle is: :meth:`register_dataset` once per dataset, then
    :meth:`observe_batch` for each batch from the iterators, then
    :meth:`finalize_field_metadata` and :meth:`finalize_totals` when every
    batch has been spooled.

    Sample list construction mirrors today's PCP/AIRR pipeline: each sample
    entry the caller supplies is appended only the first time its
    ``sample_id`` is seen.  AIRR datasets pass their input ``samples``
    list once; PCP streaming passes a fresh per-clone sample dict and we
    dedupe on the fly.
    """

    def __init__(self, allow_duplicate_ids: bool = False):
        self._datasets: Dict[str, _DatasetState] = {}
        self._dataset_ids: Set[str] = set()
        self._allow_duplicate_ids = allow_duplicate_ids
        self._duplicate_warnings: List[str] = []

    @property
    def duplicate_warnings(self) -> List[str]:
        return list(self._duplicate_warnings)

    def register_dataset(
        self, dataset_id: str, *, hoist_tree_extras_to_clone: bool = True
    ) -> None:
        if dataset_id in self._dataset_ids:
            self._on_duplicate(
                f"duplicate dataset_id: {dataset_id!r}", scope="dataset_id"
            )
        self._dataset_ids.add(dataset_id)
        state = self._datasets.setdefault(dataset_id, _DatasetState())
        state.hoist_tree_extras_to_clone = hoist_tree_extras_to_clone

    def add_samples(self, dataset_id: str, samples: Iterable[Dict[str, Any]]) -> None:
        """Append sample dicts to the dataset, deduplicated by ``sample_id``.

        Called by the PCP streaming wrapper after each batch yields newly
        seen samples, and by the AIRR streaming wrapper once with the
        input dataset's samples list.
        """
        state = self._require_state(dataset_id)
        for entry in samples:
            sid = entry.get("sample_id")
            if sid is None or sid in state.sample_ids:
                continue
            state.sample_ids.add(sid)
            state.samples.append(entry)

    def observe_batch(
        self,
        dataset_id: str,
        clones: List[Dict[str, Any]],
        trees: List[Dict[str, Any]],
    ) -> None:
        state = self._require_state(dataset_id)

        for clone in clones:
            self._observe_clone(state, clone)
        for tree in trees:
            self._observe_tree(state, tree)

    # ----- internal helpers ------------------------------------------------

    def _on_duplicate(self, message: str, *, scope: DuplicateScope) -> None:
        if self._allow_duplicate_ids:
            self._duplicate_warnings.append(message)
            return
        raise DuplicateIdError(message, scope=scope)

    def _require_state(self, dataset_id: str) -> _DatasetState:
        state = self._datasets.get(dataset_id)
        if state is None:
            raise KeyError(f"dataset_id {dataset_id!r} not registered")
        return state

    def _observe_clone(self, state: _DatasetState, clone: Dict[str, Any]) -> None:
        state.clone_count += 1
        clone_id = clone.get("clone_id")
        if clone_id is not None:
            if clone_id in state.clone_ids:
                self._on_duplicate(
                    f"duplicate clone_id: {clone_id!r}", scope="clone_id"
                )
            state.clone_ids.add(clone_id)

        subject_id = clone.get("subject_id")
        if subject_id is not None:
            state.subject_ids.add(subject_id)

        clone_level = state.levels["clone"]
        for key, value in clone.items():
            if key == "trees":
                continue
            clone_level.record(key, value)

        # Path-based known clone fields (e.g. ``locus`` → ``sample.locus``).
        # The legacy ``generate_clone_metadata`` resolves these via
        # ``sample_values_by_path``; do the same incrementally so the
        # streaming finalize sees evidence for them.
        for fname, finfo in KNOWN_CLONE_FIELDS.items():
            path = finfo.get("path")
            if not path or fname in clone:
                continue
            value = get_nested_value(clone, path)
            if value is not None:
                clone_level.record(fname, value)

        # Tree-level classification: each clone's trees array is fully
        # contained in this batch (iter_pcp_clone_groups guarantees alt
        # reconstructions co-emit), so any intra-clone variance here is the
        # dataset-wide tree-level signal.
        tree_refs = clone.get("trees") or []
        tree_level = state.levels["tree"]
        candidate_keys: Set[str] = set()
        for ref in tree_refs:
            if isinstance(ref, dict):
                candidate_keys.update(ref.keys())
        candidate_keys -= EXCLUDED_TREE_FIELDS
        for key in candidate_keys:
            for ref in tree_refs:
                if isinstance(ref, dict):
                    tree_level.record(key, ref.get(key))

        if len(tree_refs) >= 2:
            # Match the legacy ``_hoist_clone_invariant_extras`` variance check:
            # ``tree.get(key)`` (default None) — absent-key and present-None
            # collapse to the same value, so mixed-presence only counts as
            # variance when the present value differs from None.
            for key in candidate_keys - state.tree_level_keys:
                distinct = set()
                for ref in tree_refs:
                    if isinstance(ref, dict):
                        v = ref.get(key)
                        try:
                            hash(v)
                            distinct.add(v)
                        except TypeError:
                            distinct.add(("__unhashable__", repr(v)))
                if len(distinct) > 1:
                    state.tree_level_keys.add(key)

    def _observe_tree(self, state: _DatasetState, tree: Dict[str, Any]) -> None:
        state.tree_count += 1
        clone_id = tree.get("clone_id")
        tree_id = tree.get("tree_id")
        if clone_id and tree_id:
            seen = state.tree_ids_by_clone.setdefault(clone_id, set())
            if tree_id in seen:
                self._on_duplicate(
                    f"duplicate tree_id {tree_id!r} within clone {clone_id!r}",
                    scope="tree_id",
                )
            seen.add(tree_id)

        nodes = tree.get("nodes")
        if isinstance(nodes, dict):
            nodes = list(nodes.values())
        if not isinstance(nodes, list):
            return

        node_level = state.levels["node"]
        branch_level = state.levels["branch"]
        mutation_level = state.levels["mutation"]

        for node in nodes:
            if not isinstance(node, dict):
                continue
            if node.get("type") == "leaf":
                state.leaf_count += 1
            if node.get("sequence_alignment_aa"):
                state.has_aa_sequences = True

            for key, value in node.items():
                if key == "mutations":
                    continue
                if key not in EXCLUDED_NODE_FIELDS:
                    node_level.record(key, value)
                if key in KNOWN_BRANCH_FIELDS or key not in EXCLUDED_BRANCH_FIELDS:
                    branch_level.record(key, value)

            for mutation in node.get("mutations") or []:
                if not isinstance(mutation, dict):
                    continue
                for key, value in mutation.items():
                    if key in EXCLUDED_MUTATION_FIELDS:
                        continue
                    mutation_level.record(key, value)

    # ----- finalize -------------------------------------------------------

    @staticmethod
    def _apply_field_aliases(
        metadata: Dict[str, Dict[str, Any]],
        custom_fields: Optional[List[Dict[str, Any]]],
    ) -> None:
        """Rename source→target keys per ``FIELD_ALIASES``.

        Auto-config-path only. Handles the cross-format-alias half of what
        ``generate_default_config`` + ``apply_custom_fields`` does on the
        legacy path: when the target was also discovered in the data its
        entry wins (matches the legacy iteration where the explicit
        target entry overwrites the aliased source). Does **not** cover
        mutation-array demotion or skip suggestions — :func:`apply_suggestions`
        handles the suggestion side; user ``custom_fields`` drive
        demotion.

        Skipped when the caller supplied ``custom_fields``.
        """
        if custom_fields is not None:
            return
        for source, target in FIELD_ALIASES.items():
            if source == target or source not in metadata:
                continue
            if target not in metadata:
                metadata[target] = metadata[source]
            del metadata[source]

    def finalize_field_metadata(
        self,
        dataset_id: str,
        custom_fields: Optional[List[Dict[str, Any]]] = None,
    ) -> Dict[str, Dict[str, Dict[str, Any]]]:
        """Build the per-level ``field_metadata`` dict for ``dataset_id``.

        Shape matches ``field_metadata.generate_field_metadata`` so the
        webapp consumes streaming output the same way it consumes the
        legacy path.
        """
        state = self._require_state(dataset_id)
        result: Dict[str, Dict[str, Dict[str, Any]]] = {}

        clone_meta = self._clone_level_metadata(state, custom_fields)
        if clone_meta:
            result["clone"] = clone_meta

        tree_meta = self._tree_metadata(state, custom_fields)
        if tree_meta:
            result["tree"] = tree_meta

        node_meta = self._level_metadata(
            state,
            "node",
            custom_fields,
            demoted_sources=_collect_demoted_sources(custom_fields),
        )
        if node_meta:
            result["node"] = node_meta

        branch_meta = self._level_metadata(
            state,
            "branch",
            custom_fields,
            only_known=True,
        )
        if branch_meta:
            result["branch"] = branch_meta

        mutation_meta = self._mutation_metadata(state, custom_fields)
        if mutation_meta:
            result["mutation"] = mutation_meta

        return result

    def finalize_totals(self) -> Dict[str, Any]:
        """Aggregate ``processing_info`` totals across all registered datasets."""
        datasets_count = len(self._dataset_ids)
        total_clones = sum(s.clone_count for s in self._datasets.values())
        total_trees = sum(s.tree_count for s in self._datasets.values())
        total_leaves = sum(s.leaf_count for s in self._datasets.values())
        return {
            "datasets_count": datasets_count,
            "total_clones_count": total_clones,
            "total_trees_count": total_trees,
            "total_leaf_nodes_count": total_leaves,
        }

    def samples_for(self, dataset_id: str) -> List[Dict[str, Any]]:
        return list(self._require_state(dataset_id).samples)

    def tree_level_keys(self, dataset_id: str) -> Set[str]:
        """Snapshot of dataset-scope keys whose value varies within some clone's trees."""
        return set(self._require_state(dataset_id).tree_level_keys)

    # ----- metadata helpers ----------------------------------------------

    def _clone_level_metadata(
        self,
        state: _DatasetState,
        custom_fields: Optional[List[Dict[str, Any]]],
    ) -> Dict[str, Dict[str, Any]]:
        """Build clone-level metadata, including tree-csv extras that hoist.

        A tree-csv extra hoists to the clone when it's *not* in
        ``state.tree_level_keys`` — i.e. no clone's alt-reconstructions
        disagreed on its value, so the value is clone-invariant by the
        same rule ``_hoist_clone_invariant_extras`` enforces on the legacy
        path.  Evidence for those keys lives under the ``tree`` level (the
        ``_observe_clone`` walk records via tree refs).
        """
        clone_state = state.levels["clone"]
        tree_state = state.levels["tree"]

        # PCP-style tree-csv extras hoist to clone level when constant
        # across every clone's trees (``state.tree_level_keys`` captures
        # the inverse).  AIRR places these fields on trees natively and
        # turns the flag off so they stay out of clone metadata.
        # Structural tree keys never hoist either way (mirrors
        # ``process_pcp_data._HOIST_STRUCTURAL_KEYS``).
        if state.hoist_tree_extras_to_clone:
            hoist_keys = tree_state.keys - state.tree_level_keys - HOIST_STRUCTURAL_KEYS
        else:
            hoist_keys = set()
        universe = (clone_state.keys | hoist_keys) - EXCLUDED_CLONE_FIELDS

        # Clone level mirrors ``generate_clone_metadata`` — no range key.
        metadata: Dict[str, Dict[str, Any]] = {}
        for key in sorted(universe):
            evidence = clone_state.type_evidence.get(key)
            if evidence is None or evidence.total() == 0:
                evidence = tree_state.type_evidence.get(key)
            if evidence is None or evidence.total() == 0:
                continue

            if key in KNOWN_CLONE_FIELDS:
                entry = entry_from_known(KNOWN_CLONE_FIELDS[key])
            else:
                entry = make_entry(evidence.infer(), humanize_label(key))

            metadata[key] = entry

        apply_suggestions(metadata)
        self._apply_field_aliases(metadata, custom_fields)
        apply_custom_fields(metadata, custom_fields, "clone", None)
        return metadata

    def _level_metadata(
        self,
        state: _DatasetState,
        level: str,
        custom_fields: Optional[List[Dict[str, Any]]],
        exclude_keys: Optional[Set[str]] = None,
        demoted_sources: Optional[Set[str]] = None,
        only_known: bool = False,
    ) -> Dict[str, Dict[str, Any]]:
        level_state = state.levels[level]
        if not level_state.keys:
            return {}

        excluded = _EXCLUDED_BY_LEVEL[level]
        known = _KNOWN_BY_LEVEL[level]
        candidate = level_state.keys - excluded
        if exclude_keys:
            candidate = candidate - exclude_keys
        if demoted_sources:
            candidate = candidate - demoted_sources

        # Node-level skips keys that belong on the branch level (e.g.
        # ``length``).  Mirrors the ``elif key in KNOWN_BRANCH_FIELDS:
        # continue`` branch in ``generate_node_metadata``.
        cross_level_skip: Set[str] = set()
        if level == "node":
            cross_level_skip = set(KNOWN_BRANCH_FIELDS)
        candidate = candidate - cross_level_skip

        # Match the legacy generators: only tree and mutation levels
        # attach a continuous ``range``.  Clone/node/branch ranges aren't
        # consumed by the webapp's color-scale paths today and would
        # diff existing field_metadata output.
        emit_range = level in ("tree", "mutation")

        metadata: Dict[str, Dict[str, Any]] = {}
        for key in sorted(candidate):
            if key in known:
                evidence = level_state.type_evidence.get(key)
                if evidence is None or evidence.total() == 0:
                    continue
                entry = entry_from_known(known[key])
                if emit_range and entry["type"] == "continuous":
                    rng = level_state.range_evidence.get(key)
                    if rng:
                        as_list = rng.as_list()
                        if as_list:
                            entry["range"] = as_list
                metadata[key] = entry
            elif only_known:
                continue
            else:
                evidence = level_state.type_evidence.get(key)
                if evidence is None or evidence.total() == 0:
                    continue
                field_type = evidence.infer()
                entry = make_entry(field_type, humanize_label(key))
                if emit_range and field_type == "continuous":
                    rng = level_state.range_evidence.get(key)
                    if rng:
                        as_list = rng.as_list()
                        if as_list:
                            entry["range"] = as_list
                metadata[key] = entry

        apply_suggestions(metadata)
        self._apply_field_aliases(metadata, custom_fields)
        apply_custom_fields(metadata, custom_fields, level, None)
        return metadata

    def _tree_metadata(
        self,
        state: _DatasetState,
        custom_fields: Optional[List[Dict[str, Any]]],
    ) -> Dict[str, Dict[str, Any]]:
        # Tree level only emits fields the variance classifier promoted.
        if not state.tree_level_keys:
            metadata: Dict[str, Dict[str, Any]] = {}
            apply_custom_fields(metadata, custom_fields, "tree", None)
            return metadata

        level_state = state.levels["tree"]
        metadata = {}
        for key in sorted(state.tree_level_keys):
            evidence = level_state.type_evidence.get(key)
            if evidence is None or evidence.total() == 0:
                continue
            if key in KNOWN_TREE_FIELDS:
                entry = entry_from_known(KNOWN_TREE_FIELDS[key])
            else:
                entry = make_entry(evidence.infer(), humanize_label(key))
            if entry["type"] == "continuous":
                rng = level_state.range_evidence.get(key)
                if rng:
                    as_list = rng.as_list()
                    if as_list:
                        entry["range"] = as_list
            metadata[key] = entry

        apply_suggestions(metadata)
        self._apply_field_aliases(metadata, custom_fields)
        apply_custom_fields(metadata, custom_fields, "tree", None)
        return metadata

    def _mutation_metadata(
        self,
        state: _DatasetState,
        custom_fields: Optional[List[Dict[str, Any]]],
    ) -> Dict[str, Dict[str, Any]]:
        level_state = state.levels["mutation"]
        # Match field_metadata.generate_mutation_metadata: if no mutations
        # were observed at all, the webapp still derives child_aa /
        # parent_aa whenever nodes carry sequence_alignment_aa.
        if not level_state.keys:
            metadata: Dict[str, Dict[str, Any]] = {}
            if state.has_aa_sequences:
                metadata["child_aa"] = make_entry("aa", "Child Amino Acid")
                metadata["parent_aa"] = make_entry(
                    "aa", "Parent Amino Acid", display="tooltip"
                )
            apply_suggestions(metadata)
            self._apply_field_aliases(metadata, custom_fields)
            apply_custom_fields(metadata, custom_fields, "mutation", None)
            return metadata

        return self._level_metadata(state, "mutation", custom_fields)




def _collect_demoted_sources(
    custom_fields: Optional[List[Dict[str, Any]]],
) -> Set[str]:
    """Return node-level field names that should be hidden because a
    mutation-level custom field encodes them.

    Mirrors the demotion logic in
    ``field_metadata.generate_node_metadata``: a records-encoded mutation
    custom field with ``source: X`` removes X from node level; list/json
    encodings remove the named field itself.  Custom fields explicitly
    declared at ``level: node`` (not skipped) re-include their names.
    """
    demoted: Set[str] = set()
    node_overrides: Set[str] = set()
    if not custom_fields:
        return demoted
    for cf in custom_fields:
        norm = normalize_level(cf.get("level", ""))
        if norm == "mutation" and cf.get("encoding"):
            encoding = cf["encoding"]
            if encoding == "records" and cf.get("source"):
                demoted.add(cf["source"])
            elif encoding in ("list", "json"):
                demoted.add(cf["name"])
        elif norm == "node" and not cf.get("skip") and cf.get("display") != "skip":
            node_overrides.add(cf["name"])
    return demoted - node_overrides


# =============================================================================
# Per-batch on-disk spooler
# =============================================================================


class BatchSpooler:
    """Spools each batch's clones and trees to JSONL temp files per dataset.

    Temp files are written as JSONL (one JSON document per line, compact)
    so the final writer can stream-read them without buffering more than
    a single record at a time.  Per-batch files within a dataset are
    consumed in write order, preserving the batch-iteration order in the
    final output.

    Use as a context manager so the temp directory is removed on exit.
    """

    def __init__(self) -> None:
        self._tempdir: Optional[tempfile.TemporaryDirectory] = None
        self._root: Optional[Path] = None
        self._clones_files: Dict[str, List[Path]] = {}
        self._trees_files: Dict[str, List[Path]] = {}
        self._batch_counter: Dict[str, int] = {}

    def __enter__(self) -> "BatchSpooler":
        self._tempdir = tempfile.TemporaryDirectory(prefix="olmsted_stream_")
        self._root = Path(self._tempdir.name)
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        if self._tempdir is not None:
            self._tempdir.cleanup()
        self._tempdir = None
        self._root = None

    def write_batch(
        self,
        dataset_id: str,
        clones: List[Dict[str, Any]],
        trees: List[Dict[str, Any]],
    ) -> None:
        if self._root is None:
            raise RuntimeError("BatchSpooler must be used as a context manager")
        idx = self._batch_counter.get(dataset_id, 0)
        self._batch_counter[dataset_id] = idx + 1
        safe = _safe_filename(dataset_id)
        clones_path = self._root / f"{safe}_clones_{idx:06d}.jsonl"
        trees_path = self._root / f"{safe}_trees_{idx:06d}.jsonl"
        _write_jsonl(clones_path, clones)
        _write_jsonl(trees_path, trees)
        self._clones_files.setdefault(dataset_id, []).append(clones_path)
        self._trees_files.setdefault(dataset_id, []).append(trees_path)

    def dataset_ids(self) -> List[str]:
        return list(self._clones_files.keys())

    def iter_clones(self, dataset_id: str) -> Iterator[Dict[str, Any]]:
        for path in self._clones_files.get(dataset_id, []):
            yield from _read_jsonl(path)

    def iter_trees(self, dataset_id: str) -> Iterator[Dict[str, Any]]:
        for path in self._trees_files.get(dataset_id, []):
            yield from _read_jsonl(path)

    def clone_paths(self, dataset_id: str) -> List[Path]:
        return list(self._clones_files.get(dataset_id, []))

    def tree_paths(self, dataset_id: str) -> List[Path]:
        return list(self._trees_files.get(dataset_id, []))


def _safe_filename(value: str) -> str:
    return "".join(c if c.isalnum() or c in "-_." else "_" for c in value)[:80]


def _write_jsonl(path: Path, items: Iterable[Dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8") as fh:
        for item in items:
            json.dump(item, fh, separators=(",", ":"), default=str)
            fh.write("\n")


def _read_jsonl(path: Path) -> Iterator[Dict[str, Any]]:
    with path.open("r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            yield json.loads(line)


# Structural keys on tree records that should never be hoisted onto a clone.
# Kept in lockstep with ``process_pcp_data._HOIST_STRUCTURAL_KEYS``.
HOIST_STRUCTURAL_KEYS: Set[str] = {
    "ident",
    "clone_id",
    "tree_id",
    "tree_name",
    "newick",
    "nodes",
    "type",
    "reconstruction_method",
}


def apply_dataset_hoist(
    spooler: BatchSpooler,
    dataset_id: str,
    tree_level_keys: Set[str],
) -> None:
    """Rewrite spooled clones and trees to apply the dataset-scope hoist.

    Streaming wrote each batch's clones and trees un-hoisted (tree-csv
    extras on every tree record).  Once every batch has been observed,
    the accumulator knows which keys are dataset-tree-level — keys whose
    value varied within at least one clone's alt-reconstructions.  The
    rest are clone-invariant by the same rule
    ``_hoist_clone_invariant_extras`` enforces on the legacy path.

    For each spooled clone we hoist the non-tree-level keys that are
    constant across its trees up to the clone dict and strip them from
    both the trimmed tree refs in ``clone["trees"]`` and the matching
    full tree records.  The hoist set is per-clone (key may hoist on
    clone A and not clone B if A's trees agree but B's don't) — but
    keys in ``tree_level_keys`` never hoist, even when a clone happens to
    have only one tree.

    Bounded memory: one record at a time plus a per-clone
    ``{clone_id: hoist_keys}`` map.
    """
    clone_hoist_map: Dict[str, Set[str]] = {}

    for path in spooler.clone_paths(dataset_id):
        _rewrite_jsonl(
            path,
            lambda clone: _hoist_one_clone(clone, tree_level_keys, clone_hoist_map),
        )
    for path in spooler.tree_paths(dataset_id):
        _rewrite_jsonl(
            path,
            lambda tree: _strip_hoisted_keys_on_tree(tree, clone_hoist_map),
        )


def _hoist_one_clone(
    clone: Dict[str, Any],
    tree_level_keys: Set[str],
    clone_hoist_map: Dict[str, Set[str]],
) -> Dict[str, Any]:
    tree_refs = clone.get("trees") or []
    dict_refs = [ref for ref in tree_refs if isinstance(ref, dict)]
    if not dict_refs:
        return clone

    # By the time we get here the accumulator's tree-level variance scan
    # has already promoted every non-clone-invariant key into
    # ``tree_level_keys``, so any remaining candidate is guaranteed to agree
    # within this clone's trees.  Match the legacy
    # ``_hoist_clone_invariant_extras`` second pass: hoist
    # ``dict_refs[0].get(key)`` unconditionally (no per-clone agreement
    # recheck; redundant given the dataset-scope variance filter).
    candidate_keys: Set[str] = set()
    for ref in dict_refs:
        candidate_keys.update(ref.keys())
    candidate_keys -= HOIST_STRUCTURAL_KEYS
    hoist_keys = candidate_keys - tree_level_keys

    if hoist_keys:
        first = dict_refs[0]
        for key in hoist_keys:
            clone[key] = first.get(key)
            for ref in dict_refs:
                ref.pop(key, None)
        clone_id = clone.get("clone_id")
        if clone_id is not None:
            clone_hoist_map[clone_id] = hoist_keys

    return clone


def _strip_hoisted_keys_on_tree(
    tree: Dict[str, Any],
    clone_hoist_map: Dict[str, Set[str]],
) -> Dict[str, Any]:
    clone_id = tree.get("clone_id")
    if clone_id is None:
        return tree
    hoist_keys = clone_hoist_map.get(clone_id)
    if not hoist_keys:
        return tree
    for key in hoist_keys:
        tree.pop(key, None)
    return tree


def _rewrite_jsonl(path: Path, transform) -> None:
    """Read a JSONL file, apply ``transform`` to each record, write back in place."""
    tmp = path.with_name(path.name + ".rewrite")
    with path.open("r", encoding="utf-8") as ifh, tmp.open(
        "w", encoding="utf-8"
    ) as ofh:
        for line in ifh:
            line = line.strip()
            if not line:
                continue
            record = json.loads(line)
            record = transform(record)
            json.dump(record, ofh, separators=(",", ":"), default=str)
            ofh.write("\n")
    tmp.replace(path)


# =============================================================================
# Streaming consolidated-JSON writer
# =============================================================================


def write_olmsted_json_streaming(
    metadata: Dict[str, Any],
    datasets: List[Dict[str, Any]],
    spooler: BatchSpooler,
    output_path: str,
    json_format: JsonOutputFormat = "pretty",
) -> str:
    """Stream the consolidated Olmsted JSON to ``output_path``.

    Writes keys in canonical order — ``metadata`` first (with finalized
    totals), then ``datasets``, then ``clones`` and ``trees`` streamed
    from ``spooler`` so peak memory tracks one record rather than the
    whole dataset.  ``json_format`` mirrors ``data_io.write_olmsted_json``:
    pretty/compact/gzip.  Gzip output pins ``mtime=0`` and the embedded
    filename so the compression layer is deterministic, matching the
    legacy writer.

    The resulting JSON is parse-equivalent to a single-pass ``json.dump``
    of ``{metadata, datasets, clones, trees}``.  Byte-level whitespace may
    differ slightly from the legacy writer because pretty payloads are
    emitted item-by-item; consumers compare parsed JSON.
    """
    output_path = str(output_path)
    if json_format == "gzip" and not output_path.endswith(".gz"):
        output_path = output_path + ".gz"
    use_gzip = json_format == "gzip" or output_path.endswith(".gz")

    if use_gzip:
        with open(output_path, "wb") as raw:
            with gzip.GzipFile(filename="", fileobj=raw, mode="wb", mtime=0) as gz:
                with io.TextIOWrapper(gz, encoding="utf-8") as fh:
                    _emit(fh, metadata, datasets, spooler, pretty=True)
    else:
        with open(output_path, "w", encoding="utf-8") as fh:
            _emit(
                fh,
                metadata,
                datasets,
                spooler,
                pretty=json_format != "compact",
            )

    return output_path


def _emit(
    fh,
    metadata: Dict[str, Any],
    datasets: List[Dict[str, Any]],
    spooler: BatchSpooler,
    *,
    pretty: bool,
) -> None:
    nl = "\n" if pretty else ""
    sp = "    " if pretty else ""
    colon = ": " if pretty else ":"
    comma = "," + nl

    fh.write("{" + nl)

    fh.write(sp + '"metadata"' + colon)
    _dump_value(fh, metadata, pretty=pretty, indent_level=1)
    fh.write(comma)

    fh.write(sp + '"datasets"' + colon)
    _dump_value(fh, datasets, pretty=pretty, indent_level=1)
    fh.write(comma)

    fh.write(sp + '"clones"' + colon + "{" + nl)
    # Drive the clones map from the dataset list, not the spooler — a
    # dataset that registered but spooled zero batches (e.g. an AIRR
    # input with empty clones[]) still needs a ``"<id>": []`` entry to
    # match the legacy output shape.
    dataset_ids = [d["dataset_id"] for d in datasets if d.get("dataset_id")]
    for i, ds_id in enumerate(dataset_ids):
        fh.write(sp * 2 + json.dumps(ds_id) + colon + "[")
        first = True
        for clone in spooler.iter_clones(ds_id):
            if not first:
                fh.write(",")
            fh.write(nl + sp * 3)
            _dump_value(fh, clone, pretty=pretty, indent_level=3)
            first = False
        if not first:
            fh.write(nl + sp * 2)
        fh.write("]")
        if i < len(dataset_ids) - 1:
            fh.write(",")
        fh.write(nl)
    fh.write(sp + "}" + comma)

    fh.write(sp + '"trees"' + colon + "[")
    first = True
    for ds_id in dataset_ids:
        for tree in spooler.iter_trees(ds_id):
            if not first:
                fh.write(",")
            fh.write(nl + sp * 2)
            _dump_value(fh, tree, pretty=pretty, indent_level=2)
            first = False
    if not first:
        fh.write(nl + sp)
    fh.write("]" + nl)

    fh.write("}" + nl)


def _dump_value(fh, value: Any, *, pretty: bool, indent_level: int) -> None:
    """Serialize ``value`` to ``fh``, with pretty output shifted to nest at ``indent_level``.

    ``json.dumps(value, indent=4)`` emits each line starting at column 0;
    when the value lands inside an array or object we prefix the per-line
    pad with ``indent_level * 4`` spaces so the result reads cleanly when
    parsed-equivalence isn't enough and a human cracks the file open.
    """
    if not pretty:
        json.dump(value, fh, separators=(",", ":"), default=str)
        return
    text = json.dumps(value, indent=4, default=str)
    pad = "    " * indent_level
    if "\n" in text:
        head, rest = text.split("\n", 1)
        text = head + "\n" + pad + rest.replace("\n", "\n" + pad)
    fh.write(text)
