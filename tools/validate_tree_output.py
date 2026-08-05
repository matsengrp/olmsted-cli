#!/usr/bin/env python3
"""Validate a produced Olmsted JSON's ``mean_mut_freq`` against independently
recomputed references from the raw input.

For each clone, recomputes ``mean_mut_freq`` from the raw input's germline
alignment and leaf nodes under two conventions:

- **weighted**: multiplicity-weighted mean of per-leaf mutation frequency
  (the convention olmsted-cli's PCP path has always used, and which the AIRR
  path adopts per issue 24).
- **unweighted**: plain mean of per-leaf mutation frequency, each leaf
  counted once regardless of multiplicity (the convention some upstream
  AIRR producers, e.g. CFT circa 2020, used instead).

Reports, per clone, whether the value actually present in the output JSON
matches the weighted reference, the unweighted reference, both (when they
coincide, e.g. all leaf multiplicities are 1), or neither.

Usage:
    python tools/validate_tree_output.py --input raw.json --output produced.json
"""

from __future__ import annotations

import argparse
import json
import sys
from typing import Any, Dict, Iterable, Optional, Tuple

DEFAULT_TOLERANCE = 1e-6


def mutation_frequency(germline: str, leaf: str) -> Tuple[float, int]:
    """Fraction of non-gap-on-either-side positions that differ.

    Truncates to the shorter sequence length, matching olmsted-cli's
    default ``alignment_method="truncate"``.
    """
    length = min(len(germline), len(leaf))
    if length == 0:
        return 0.0, 0
    g = germline[:length]
    leaf_seq = leaf[:length]
    mismatches = sum(
        1
        for gb, lb in zip(g, leaf_seq)
        if gb != lb and gb not in ("", ".") and lb not in ("", ".")
    )
    return mismatches / length, length


def iter_leaf_nodes(nodes: Any) -> Iterable[Dict[str, Any]]:
    """Yield leaf node dicts (type == "leaf", multiplicity > 0) from a
    nodes collection shaped as either a dict (keyed by sequence_id) or a
    list."""
    values = nodes.values() if isinstance(nodes, dict) else (nodes or [])
    for node in values:
        if node.get("type") == "leaf" and node.get("multiplicity", 0) > 0:
            yield node


def compute_references(germline: str, nodes: Any) -> Tuple[float, float]:
    """Return (weighted_mean, unweighted_mean) mutation frequency for a clone."""
    total_weighted = 0.0
    total_multiplicity = 0
    freqs = []

    for node in iter_leaf_nodes(nodes):
        leaf_seq = node.get("sequence_alignment", "")
        if not germline or not leaf_seq:
            continue
        freq, length = mutation_frequency(germline, leaf_seq)
        if length == 0:
            continue
        multiplicity = node.get("multiplicity", 0)
        freqs.append(freq)
        total_weighted += freq * multiplicity
        total_multiplicity += multiplicity

    weighted = total_weighted / total_multiplicity if total_multiplicity > 0 else 0.0
    unweighted = sum(freqs) / len(freqs) if freqs else 0.0
    return weighted, unweighted


def load_clones_by_id(data: Dict[str, Any]) -> Dict[str, Dict[str, Any]]:
    """Return {clone_id: clone_dict}, accepting either a raw AIRR-shaped
    JSON (``clones`` is a list) or a consolidated Olmsted JSON (``clones``
    is a dict of dataset_id -> list)."""
    clones = data.get("clones", [])
    if isinstance(clones, dict):
        result: Dict[str, Dict[str, Any]] = {}
        for clone_list in clones.values():
            for clone in clone_list:
                result[clone["clone_id"]] = clone
        return result
    return {clone["clone_id"]: clone for clone in clones}


def canonical_tree_nodes(clone: Dict[str, Any]) -> Any:
    """The first tree's nodes (canonical reconstruction) for a clone."""
    trees = clone.get("trees") or []
    if not trees:
        return {}
    return trees[0].get("nodes", {})


def _match_label(
    actual: Optional[float], weighted: float, unweighted: float, tolerance: float
) -> Tuple[str, str]:
    """Return (status, message) for one clone's comparison."""
    matches_weighted = abs(actual - weighted) <= tolerance
    matches_unweighted = abs(actual - unweighted) <= tolerance

    if matches_weighted and matches_unweighted:
        return "PASS", f"matches both conventions (actual={actual:.6f})"
    if matches_weighted:
        return (
            "PASS",
            f"matches WEIGHTED reference (actual={actual:.6f}, "
            f"unweighted={unweighted:.6f})",
        )
    if matches_unweighted:
        return (
            "WARN",
            f"matches UNWEIGHTED reference, not weighted "
            f"(actual={actual:.6f}, weighted={weighted:.6f})",
        )
    return (
        "WARN",
        f"matches NEITHER convention (actual={actual:.6f}, "
        f"weighted={weighted:.6f}, unweighted={unweighted:.6f})",
    )


def validate(
    input_data: Dict[str, Any], output_data: Dict[str, Any], tolerance: float
) -> int:
    """Print a per-clone report; return an exit code (0 if all PASS)."""
    input_clones = load_clones_by_id(input_data)
    output_clones = load_clones_by_id(output_data)

    counts = {"PASS": 0, "WARN": 0, "FAIL": 0}
    for clone_id in sorted(input_clones):
        in_clone = input_clones[clone_id]
        out_clone = output_clones.get(clone_id)
        if out_clone is None:
            print(f"FAIL {clone_id}: not present in output")
            counts["FAIL"] += 1
            continue

        actual = out_clone.get("mean_mut_freq")
        if actual is None:
            print(f"FAIL {clone_id}: no mean_mut_freq in output")
            counts["FAIL"] += 1
            continue

        germline = in_clone.get("germline_alignment", "")
        nodes = canonical_tree_nodes(in_clone)
        weighted, unweighted = compute_references(germline, nodes)

        status, message = _match_label(actual, weighted, unweighted, tolerance)
        print(f"{status} {clone_id}: {message}")
        counts[status] += 1

    total = sum(counts.values())
    print(
        f"\n{counts['PASS']} PASS, {counts['WARN']} WARN, {counts['FAIL']} FAIL "
        f"({total} clones)"
    )
    return 1 if (counts["WARN"] or counts["FAIL"]) else 0


def get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Raw input JSON (AIRR-shaped)")
    parser.add_argument(
        "--output", required=True, help="Produced Olmsted output JSON to validate"
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        default=DEFAULT_TOLERANCE,
        help=f"Absolute tolerance for matching a reference (default: {DEFAULT_TOLERANCE})",
    )
    return parser.parse_args()


def main() -> None:
    args = get_args()
    with open(args.input) as fh:
        input_data = json.load(fh)
    with open(args.output) as fh:
        output_data = json.load(fh)

    exit_code = validate(input_data, output_data, args.tolerance)
    sys.exit(exit_code)


if __name__ == "__main__":
    main()
