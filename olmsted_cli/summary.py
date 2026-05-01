#!/usr/bin/env python3
"""Summary command for olmsted-cli - provides statistics about Olmsted JSON files."""

import argparse
import json
import sys
from pathlib import Path

from .data_io import read_olmsted_json
from .utils import set_verbosity, vprint


def analyze_consolidated_data(data):
    """
    Analyze Olmsted JSON data and extract summary statistics.
    
    Args:
        data: Olmsted JSON data dictionary
        
    Returns:
        dict: Summary statistics
    """
    summary = {}
    
    # Metadata
    metadata = data.get("metadata", {})
    summary["metadata"] = {
        "format_version": metadata.get("format_version", "Unknown"),
        "schema_version": metadata.get("schema_version", "Unknown"),
        "source_format": metadata.get("source_format", "Unknown"),
        "created_at": metadata.get("created_at", "Unknown"),
        "name": metadata.get("name", "Unnamed dataset"),
        "input_files": metadata.get("input_files", [])
    }
    
    # Datasets
    datasets = data.get("datasets", [])
    summary["datasets"] = {
        "count": len(datasets),
        "dataset_ids": [d.get("dataset_id", "Unknown") for d in datasets]
    }
    
    # Clones analysis
    clones_dict = data.get("clones", {})
    all_clones = []
    for dataset_id, clones_list in clones_dict.items():
        all_clones.extend(clones_list)
    
    summary["clones"] = {
        "total_count": len(all_clones),
        "by_dataset": {dataset_id: len(clones_list) for dataset_id, clones_list in clones_dict.items()}
    }
    
    # Clone statistics
    if all_clones:
        unique_seqs_counts = [c.get("unique_seqs_count", 0) for c in all_clones if "unique_seqs_count" in c]
        read_counts = [c.get("total_read_count", 0) for c in all_clones if "total_read_count" in c]
        mut_freqs = [c.get("mean_mut_freq", 0) for c in all_clones if "mean_mut_freq" in c]
        
        summary["clone_statistics"] = {
            "unique_sequences": {
                "min": min(unique_seqs_counts) if unique_seqs_counts else 0,
                "max": max(unique_seqs_counts) if unique_seqs_counts else 0,
                "mean": sum(unique_seqs_counts) / len(unique_seqs_counts) if unique_seqs_counts else 0,
                "total": sum(unique_seqs_counts) if unique_seqs_counts else 0
            },
            "read_counts": {
                "min": min(read_counts) if read_counts else 0,
                "max": max(read_counts) if read_counts else 0,
                "mean": sum(read_counts) / len(read_counts) if read_counts else 0,
                "total": sum(read_counts) if read_counts else 0
            },
            "mutation_frequencies": {
                "min": min(mut_freqs) if mut_freqs else 0,
                "max": max(mut_freqs) if mut_freqs else 0,
                "mean": sum(mut_freqs) / len(mut_freqs) if mut_freqs else 0
            }
        }
    
    # Trees analysis
    trees = data.get("trees", [])
    summary["trees"] = {
        "total_count": len(trees),
        "reconstruction_methods": {}
    }

    # Tree node statistics
    if trees:
        node_counts = []
        reconstruction_methods = {}

        for tree in trees:
            # Count reconstruction methods. Unset/blank → "<unspecified>"
            # (angle-braced to avoid collision with a genuine method name).
            method = tree.get("reconstruction_method") or "<unspecified>"
            reconstruction_methods[method] = reconstruction_methods.get(method, 0) + 1

            # Count nodes
            nodes = tree.get("nodes", [])
            if isinstance(nodes, list):
                node_counts.append(len(nodes))
            elif isinstance(nodes, dict):
                node_counts.append(len(nodes))

        summary["trees"]["reconstruction_methods"] = reconstruction_methods
        
        if node_counts:
            summary["tree_statistics"] = {
                "node_counts": {
                    "min": min(node_counts),
                    "max": max(node_counts),
                    "mean": sum(node_counts) / len(node_counts),
                    "total": sum(node_counts)
                }
            }
    
    # Sample analysis
    all_samples = set()
    for dataset in datasets:
        for sample in dataset.get("samples", []):
            all_samples.add(sample.get("sample_id", "Unknown"))
    
    summary["samples"] = {
        "unique_count": len(all_samples),
        "sample_ids": sorted(list(all_samples))
    }
    
    # Gene usage analysis
    v_genes = set()
    j_genes = set()
    d_genes = set()
    
    for clone in all_clones:
        if clone.get("v_call"):
            v_genes.add(clone["v_call"])
        if clone.get("j_call"):
            j_genes.add(clone["j_call"])
        if clone.get("d_call"):
            d_genes.add(clone["d_call"])
    
    summary["gene_usage"] = {
        "v_genes": {
            "unique_count": len(v_genes),
            "genes": sorted(list(v_genes))
        },
        "j_genes": {
            "unique_count": len(j_genes),
            "genes": sorted(list(j_genes))
        },
        "d_genes": {
            "unique_count": len(d_genes),
            "genes": sorted(list(d_genes))
        }
    }
    
    return summary


def format_summary_text(summary):
    """
    Format summary statistics as human-readable text.
    
    Args:
        summary: Summary statistics dictionary
        
    Returns:
        str: Formatted summary text
    """
    lines = []
    
    # Header
    dataset_name = summary["metadata"]["name"]
    lines.append(f"Dataset Summary: {dataset_name}")
    lines.append("=" * (len(f"Dataset Summary: {dataset_name}")))
    lines.append("")
    
    # Metadata
    metadata = summary["metadata"]
    lines.append("Metadata:")
    lines.append(f"  Format Version: {metadata['format_version']}")
    lines.append(f"  Schema Version: {metadata['schema_version']}")
    lines.append(f"  Source Format: {metadata['source_format']}")
    lines.append(f"  Created At: {metadata['created_at']}")
    if metadata["input_files"]:
        lines.append(f"  Input Files: {', '.join(metadata['input_files'])}")
    lines.append("")
    
    # Datasets
    datasets = summary["datasets"]
    lines.append(f"Datasets: {datasets['count']}")
    if datasets["dataset_ids"]:
        for dataset_id in datasets["dataset_ids"]:
            lines.append(f"  - {dataset_id}")
    lines.append("")
    
    # Samples
    samples = summary["samples"]
    lines.append(f"Samples: {samples['unique_count']}")
    if len(samples["sample_ids"]) <= 10:
        for sample_id in samples["sample_ids"]:
            lines.append(f"  - {sample_id}")
    else:
        for sample_id in samples["sample_ids"][:10]:
            lines.append(f"  - {sample_id}")
        lines.append(f"  ... and {len(samples['sample_ids']) - 10} more")
    lines.append("")
    
    # Clones
    clones = summary["clones"]
    lines.append(f"Clonal Families: {clones['total_count']}")
    for dataset_id, count in clones["by_dataset"].items():
        lines.append(f"  {dataset_id}: {count} families")
    lines.append("")
    
    # Clone statistics
    if "clone_statistics" in summary:
        stats = summary["clone_statistics"]
        lines.append("Clone Statistics:")
        
        unique_seqs = stats["unique_sequences"]
        lines.append(f"  Unique Sequences: {unique_seqs['total']:,} total")
        lines.append(f"    Range: {unique_seqs['min']}-{unique_seqs['max']} per family")
        lines.append(f"    Mean: {unique_seqs['mean']:.1f} per family")
        
        reads = stats["read_counts"]
        lines.append(f"  Read Counts: {reads['total']:,} total")
        lines.append(f"    Range: {reads['min']:,}-{reads['max']:,} per family")
        lines.append(f"    Mean: {reads['mean']:.1f} per family")
        
        mut_freq = stats["mutation_frequencies"]
        lines.append(f"  Mutation Frequencies:")
        lines.append(f"    Range: {mut_freq['min']:.4f}-{mut_freq['max']:.4f}")
        lines.append(f"    Mean: {mut_freq['mean']:.4f}")
        lines.append("")
    
    # Trees
    trees = summary["trees"]
    lines.append(f"Trees: {trees['total_count']}")
    if trees["reconstruction_methods"]:
        lines.append("  Reconstruction methods:")
        for method, count in trees["reconstruction_methods"].items():
            lines.append(f"    {method}: {count}")
    lines.append("")
    
    # Tree statistics
    if "tree_statistics" in summary:
        tree_stats = summary["tree_statistics"]
        nodes = tree_stats["node_counts"]
        lines.append("Tree Statistics:")
        lines.append(f"  Nodes: {nodes['total']:,} total")
        lines.append(f"    Range: {nodes['min']}-{nodes['max']} per tree")
        lines.append(f"    Mean: {nodes['mean']:.1f} per tree")
        lines.append("")
    
    # Gene usage
    gene_usage = summary["gene_usage"]
    lines.append("Gene Usage:")
    lines.append(f"  V genes: {gene_usage['v_genes']['unique_count']} unique")
    lines.append(f"  J genes: {gene_usage['j_genes']['unique_count']} unique")
    lines.append(f"  D genes: {gene_usage['d_genes']['unique_count']} unique")
    
    return "\n".join(lines)


def get_args():
    """Parse command line arguments for summary command."""
    parser = argparse.ArgumentParser(
        description="Generate summary statistics for Olmsted JSON files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic summary
  olmsted summary data.json
  
  # Output as JSON
  olmsted summary --json data.json
  
  # Save summary to file
  olmsted summary data.json -o summary.txt
        """,
    )
    
    parser.add_argument(
        "input_file",
        help="Olmsted JSON file to analyze"
    )
    
    parser.add_argument(
        "-o", "--output",
        help="Output file (default: stdout)"
    )
    
    parser.add_argument(
        "--json",
        action="store_true",
        help="Output summary as JSON instead of human-readable text"
    )
    
    parser.add_argument(
        "-v", "--verbose",
        type=int,
        choices=[0, 1, 2, 3],
        default=1,
        help="Verbosity: 0=errors only, 1=normal (default), 2=verbose (detailed gene lists), 3=debug",
    )
    parser.add_argument(
        "-q", "--quiet",
        action="store_true",
        help="Quiet mode — errors only (equivalent to -v 0)",
    )

    return parser.parse_args()


def main():
    """Main entry point for summary command."""
    args = get_args()

    # Handle quiet mode
    if getattr(args, "quiet", False):
        args.verbose = 0
    set_verbosity(args.verbose)

    # Validate input file
    input_path = Path(args.input_file)
    if not input_path.exists():
        vprint.error(f"Error: Input file not found: {args.input_file}")
        sys.exit(1)

    # Load and parse data
    try:
        data = read_olmsted_json(input_path)
    except (ValueError, OSError) as e:
        vprint.error(f"Error: {e}")
        sys.exit(1)

    # Validate data structure
    required_keys = ["metadata", "datasets", "clones", "trees"]
    missing_keys = [key for key in required_keys if key not in data]
    if missing_keys:
        vprint.error(f"Error: Not a valid Olmsted JSON format. Missing keys: {missing_keys}")
        sys.exit(1)

    # Analyze data
    try:
        summary = analyze_consolidated_data(data)
    except Exception as e:
        vprint.error(f"Error: Failed to analyze data: {e}")
        sys.exit(1)

    # Format output
    if args.json:
        output_text = json.dumps(summary, indent=2)
    else:
        output_text = format_summary_text(summary)

    # Write output
    if args.output:
        try:
            with open(args.output, 'w') as f:
                f.write(output_text)
            vprint.status(f"Summary written to: {args.output}")
        except Exception as e:
            vprint.error(f"Error: Failed to write output file: {e}")
            sys.exit(1)
    else:
        print(output_text)


if __name__ == "__main__":
    main()