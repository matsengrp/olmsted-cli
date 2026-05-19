#!/usr/bin/env python3
"""Benchmark peak RSS across --batch-size choices on synthetic PCP input.

Sanity-checks the memory promise of the streaming pipeline (#26): peak
resident memory should track ``--batch-size``, not total dataset size.

Usage:

    python scripts/benchmark_batch_memory.py             # default 300 families x 80 leaves
    python scripts/benchmark_batch_memory.py --families 500 --leaves 100
    python scripts/benchmark_batch_memory.py --batch-sizes 1 10 50 0

``--batch-size 0`` runs the legacy in-memory path for comparison.

Output: a table of ``batch_size, peak_rss, wall_time``.
"""

from __future__ import annotations

import argparse
import random
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import psutil

# Sequence shared by every node; long enough that ete3 + per-node dicts
# carry real bytes but short enough that the test wraps up quickly.
SEQ = (
    "GAATTCAAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGG"
    "GAATTCAAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGG"
    "GAATTCAAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGGTTTGCAAATTGG"
)


def _make_child_seq(rng: random.Random) -> str:
    """One-base mutation off SEQ — drives the AA-diff path during processing."""
    pos = rng.randrange(len(SEQ))
    new = rng.choice([b for b in "ACGT" if b != SEQ[pos]])
    return SEQ[:pos] + new + SEQ[pos + 1 :]


def generate_inputs(tmp: Path, n_families: int, leaves_per_family: int, seed: int = 42):
    """Write a synthetic PCP CSV + matching trees CSV under ``tmp``.

    Each family has one tree shaped as ``(L1, L2, ..., LM)naive`` — flat,
    naive-rooted, M leaves. Cheap to build, exercises the same code paths
    as real data on a per-node basis.
    """
    rng = random.Random(seed)
    pcp_path = tmp / "input-pcp.csv"
    trees_path = tmp / "input-trees.csv"

    pcp_header = (
        "sample_id,family,parent_name,child_name,parent_heavy,child_heavy,"
        "parent_is_naive,child_is_leaf,v_gene_heavy,j_gene_heavy\n"
    )
    trees_header = "family_name,sample_id,newick_tree\n"

    with pcp_path.open("w") as pcp_fh, trees_path.open("w") as trees_fh:
        pcp_fh.write(pcp_header)
        trees_fh.write(trees_header)
        for fam_i in range(n_families):
            family = f"F{fam_i:05d}"
            leaves = []
            for leaf_i in range(leaves_per_family):
                leaf = f"L{leaf_i}"
                child_seq = _make_child_seq(rng)
                pcp_fh.write(
                    f"S1,{family},naive,{leaf},{SEQ},{child_seq},"
                    f"true,true,IGHV1*01,IGHJ1*01\n"
                )
                leaves.append(leaf)
            newick = "(" + ",".join(f"{n}:0.1" for n in leaves) + ")naive;"
            trees_fh.write(f'{family},S1,"{newick}"\n')

    return pcp_path, trees_path


def run_with_peak_rss(cmd, poll_interval_s: float = 0.05):
    """Run ``cmd`` to completion, polling for peak RSS (process + children).

    Returns ``(returncode, peak_rss_bytes, wall_seconds)``.
    """
    t0 = time.perf_counter()
    proc = subprocess.Popen(
        cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
    )
    p = psutil.Process(proc.pid)
    peak = 0
    try:
        while proc.poll() is None:
            try:
                mem = p.memory_info().rss
                for child in p.children(recursive=True):
                    try:
                        mem += child.memory_info().rss
                    except (psutil.NoSuchProcess, psutil.AccessDenied):
                        pass
                peak = max(peak, mem)
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                pass
            time.sleep(poll_interval_s)
    finally:
        rc = proc.wait()
    wall = time.perf_counter() - t0
    return rc, peak, wall


def format_bytes(n: int) -> str:
    """Render bytes as a short human-readable string (MB/GB)."""
    if n >= 1024**3:
        return f"{n / 1024**3:.2f} GB"
    if n >= 1024**2:
        return f"{n / 1024**2:.1f} MB"
    if n >= 1024:
        return f"{n / 1024:.0f} KB"
    return f"{n} B"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--families", type=int, default=300,
        help="Number of clonal families to synthesize (default: 300)",
    )
    parser.add_argument(
        "--leaves", type=int, default=80,
        help="Leaves per family (default: 80)",
    )
    parser.add_argument(
        "--batch-sizes", type=int, nargs="+",
        default=[1, 10, 50, 100, 0],
        help="Batch sizes to compare. 0 = legacy in-memory path. "
        "Default: 1 10 50 100 0",
    )
    parser.add_argument(
        "--compute-metrics", action="store_true",
        help="Pass --compute-metrics (more memory per tree).",
    )
    parser.add_argument(
        "--keep-input", action="store_true",
        help="Don't delete the synthetic input dir after the run.",
    )
    args = parser.parse_args()

    tmp = Path(tempfile.mkdtemp(prefix="olmsted_bench_"))
    print(f"Generating inputs: {args.families} families × {args.leaves} leaves → {tmp}")
    pcp_path, trees_path = generate_inputs(tmp, args.families, args.leaves)
    pcp_size = pcp_path.stat().st_size
    trees_size = trees_path.stat().st_size
    print(
        f"  pcp.csv: {format_bytes(pcp_size)}, trees.csv: {format_bytes(trees_size)}"
    )

    results = []
    for batch_size in args.batch_sizes:
        out = tmp / f"out_bs{batch_size}.json"
        cmd = [
            "olmsted", "process",
            "-f", "pcp",
            "-i", str(pcp_path),
            "-t", str(trees_path),
            "-o", str(out),
            "--seed", "42",
            "--batch-size", str(batch_size),
            "-q",
        ]
        if args.compute_metrics:
            cmd.append("--compute-metrics")
        label = f"--batch-size {batch_size}" + (
            "  (legacy)" if batch_size == 0 else ""
        )
        print(f"\n>> Running {label}")
        rc, peak, wall = run_with_peak_rss(cmd)
        if rc != 0:
            print(f"   FAILED (returncode {rc})", file=sys.stderr)
            continue
        out_size = out.stat().st_size if out.exists() else 0
        print(
            f"   peak RSS: {format_bytes(peak)}    "
            f"wall: {wall:.2f}s    "
            f"output: {format_bytes(out_size)}"
        )
        results.append((batch_size, peak, wall, out_size))

    print("\n" + "=" * 64)
    print(f"{'batch_size':>12}  {'peak RSS':>12}  {'wall (s)':>10}  {'output':>10}")
    print("-" * 64)
    for batch_size, peak, wall, out_size in results:
        label = str(batch_size) if batch_size > 0 else "0 (legacy)"
        print(
            f"{label:>12}  {format_bytes(peak):>12}  "
            f"{wall:>10.2f}  {format_bytes(out_size):>10}"
        )
    print("=" * 64)

    if not args.keep_input:
        shutil.rmtree(tmp, ignore_errors=True)
    else:
        print(f"\nInput dir kept at {tmp}")


if __name__ == "__main__":
    main()
