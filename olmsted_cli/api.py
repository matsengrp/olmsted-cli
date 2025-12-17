"""
High-level API for working with Olmsted data.

This module provides the OlmstedData class, which is the main interface for
loading, manipulating, and saving Olmsted data in various formats.

Usage:
    from olmsted_cli.api import OlmstedData

    # Load from different formats
    data = OlmstedData.from_olmsted_json("output.json")
    data = OlmstedData.from_pcp("pcp.csv", "trees.csv")
    data = OlmstedData.from_airr_json("airr_data.json")

    # Access typed data
    for clone in data.clones["dataset-1"]:
        print(clone["clone_id"])

    # Save to different formats
    data.to_olmsted_json("output.json")
    data.to_pcp("pcp.csv", "trees.csv")
    data.to_airr_json("airr_output.json")
"""

from __future__ import annotations

import csv
import gzip
import json
import uuid as uuid_module
from argparse import Namespace
from pathlib import Path
from typing import Dict, List, Optional, Union

from .types import (
    OlmstedClone,
    OlmstedDataset,
    OlmstedOutput,
    OlmstedTree,
    OutputMetadata,
)


class OlmstedData:
    """
    Main class for working with Olmsted data.

    Provides methods for loading from and saving to different formats:
    - Olmsted JSON (native format)
    - AIRR JSON
    - PCP CSV + Trees CSV

    Attributes:
        datasets: List of dataset metadata
        clones: Dictionary mapping dataset_id to list of clones
        trees: List of phylogenetic trees
        metadata: Optional output metadata

    Example:
        # Load from PCP format
        data = OlmstedData.from_pcp("pcp.csv", "trees.csv")

        # Access data with type hints
        for dataset_id, clone_list in data.clones.items():
            for clone in clone_list:
                print(clone["clone_id"], clone["mean_mut_freq"])

        # Save to Olmsted JSON
        data.to_olmsted_json("output.json")
    """

    def __init__(
        self,
        datasets: List[OlmstedDataset],
        clones: Dict[str, List[OlmstedClone]],
        trees: List[OlmstedTree],
        metadata: Optional[OutputMetadata] = None,
    ):
        self.datasets = datasets
        self.clones = clones
        self.trees = trees
        self.metadata = metadata

    # -------------------------------------------------------------------------
    # Factory methods for loading data
    # -------------------------------------------------------------------------

    @classmethod
    def from_olmsted_json(
        cls,
        filepath: Union[str, Path],
    ) -> "OlmstedData":
        """
        Load data from Olmsted JSON format.

        Args:
            filepath: Path to consolidated Olmsted JSON file

        Returns:
            OlmstedData instance

        Example:
            data = OlmstedData.from_olmsted_json("output.json")
        """
        filepath = Path(filepath)

        if filepath.suffix == ".gz":
            with gzip.open(filepath, "rt") as f:
                raw_data = json.load(f)
        else:
            with open(filepath) as f:
                raw_data = json.load(f)

        return cls(
            datasets=raw_data.get("datasets", []),
            clones=raw_data.get("clones", {}),
            trees=raw_data.get("trees", []),
            metadata=raw_data.get("metadata"),
        )

    @classmethod
    def from_pcp(
        cls,
        pcp_csv: Union[str, Path],
        trees_csv: Optional[Union[str, Path]] = None,
        *,
        compute_metrics: bool = False,
        seed: Optional[int] = None,
        name: Optional[str] = None,
        verbosity: int = 1,
    ) -> "OlmstedData":
        """
        Load data from PCP CSV format.

        Args:
            pcp_csv: Path to PCP CSV file (can be gzipped)
            trees_csv: Optional path to trees CSV file (can be gzipped)
            compute_metrics: Whether to compute LBI/LBR metrics
            seed: Random seed for deterministic UUID generation
            name: Optional dataset name
            verbosity: Verbosity level (0=quiet, 1=normal, 2=verbose, 3=debug)

        Returns:
            OlmstedData instance

        Example:
            data = OlmstedData.from_pcp("pcp.csv", "trees.csv")
            data = OlmstedData.from_pcp("paired.csv.gz", "trees.csv.gz")
        """
        from .process_pcp_data import (
            deterministic_uuid,
            parse_newick_csv,
            parse_pcp_csv,
            process_pcp_to_olmsted,
        )

        # Parse input files
        pcp_families = parse_pcp_csv(str(pcp_csv))
        newick_trees = parse_newick_csv(str(trees_csv)) if trees_csv else None

        # Create UUID generator based on seed
        if seed is not None:
            counter = [0]

            def uuid_generator(prefix: str = "") -> str:
                result = deterministic_uuid(seed, counter[0])
                counter[0] += 1
                return f"{prefix}{result}"
        else:

            def uuid_generator(prefix: str = "") -> str:
                return f"{prefix}{uuid_module.uuid4()}"

        # Process to Olmsted format
        datasets, clones, trees = process_pcp_to_olmsted(
            pcp_families,
            newick_trees,
            uuid_generator=uuid_generator,
            compute_metrics=compute_metrics,
            name=name,
            verbosity=verbosity,
        )

        return cls(datasets=datasets, clones=clones, trees=trees)

    @classmethod
    def from_airr_json(
        cls,
        filepath: Union[str, Path],
        *,
        seed: Optional[int] = None,
        name: Optional[str] = None,
        verbosity: int = 1,
    ) -> "OlmstedData":
        """
        Load data from AIRR JSON format.

        Args:
            filepath: Path to AIRR JSON file
            seed: Random seed for deterministic UUID generation
            name: Optional dataset name
            verbosity: Verbosity level

        Returns:
            OlmstedData instance

        Example:
            data = OlmstedData.from_airr_json("airr_data.json")

        Note:
            For full AIRR processing options, use the CLI:
            `olmsted process -f airr -i input.json -o output.json`
        """
        from .process_airr_data import process_dataset
        from .process_pcp_data import deterministic_uuid

        filepath = Path(filepath)

        # Load AIRR JSON
        if filepath.suffix == ".gz":
            with gzip.open(filepath, "rt") as f:
                airr_data = json.load(f)
        else:
            with open(filepath) as f:
                airr_data = json.load(f)

        # Create UUID generator based on seed
        if seed is not None:
            counter = [0]

            def uuid_generator(prefix: str = "") -> str:
                result = deterministic_uuid(seed, counter[0])
                counter[0] += 1
                return f"{prefix}{result}"
        else:

            def uuid_generator(prefix: str = "") -> str:
                return f"{prefix}{uuid_module.uuid4()}"

        args = Namespace(
            uuid_generator=uuid_generator,
            verbose=verbosity > 0,
            name=name,
        )

        # Process datasets
        datasets: List[OlmstedDataset] = []
        clones_dict: Dict[str, List[OlmstedClone]] = {}
        trees: List[OlmstedTree] = []

        # AIRR format has datasets at the top level
        for dataset in airr_data.get("datasets", [airr_data]):
            processed_dataset = process_dataset(args, dataset, clones_dict, trees)
            if processed_dataset:
                datasets.append(processed_dataset)

        return cls(datasets=datasets, clones=clones_dict, trees=trees)

    # -------------------------------------------------------------------------
    # Methods for saving data
    # -------------------------------------------------------------------------

    def to_olmsted_json(
        self,
        filepath: Union[str, Path],
        *,
        indent: Optional[int] = 2,
        include_metadata: bool = True,
    ) -> None:
        """
        Save data to Olmsted JSON format.

        Args:
            filepath: Output file path (use .gz extension for compression)
            indent: JSON indentation (None for compact output)
            include_metadata: Whether to include metadata in output

        Example:
            data.to_olmsted_json("output.json")
            data.to_olmsted_json("output.json.gz")  # Compressed
        """
        filepath = Path(filepath)

        output: OlmstedOutput = {
            "datasets": self.datasets,
            "clones": self.clones,
            "trees": self.trees,
        }

        if include_metadata and self.metadata:
            output["metadata"] = self.metadata

        if filepath.suffix == ".gz":
            with gzip.open(filepath, "wt") as f:
                json.dump(output, f, indent=indent)
        else:
            with open(filepath, "w") as f:
                json.dump(output, f, indent=indent)

    def to_dict(self) -> OlmstedOutput:
        """
        Convert to dictionary representation.

        Returns:
            OlmstedOutput dictionary
        """
        result: OlmstedOutput = {
            "datasets": self.datasets,
            "clones": self.clones,
            "trees": self.trees,
        }
        if self.metadata:
            result["metadata"] = self.metadata
        return result

    def to_pcp(
        self,
        pcp_csv: Union[str, Path],
        trees_csv: Optional[Union[str, Path]] = None,
        *,
        include_light_chain: bool = True,
    ) -> None:
        """
        Export data to PCP CSV format.

        Args:
            pcp_csv: Output path for PCP CSV file (use .gz for compression)
            trees_csv: Optional output path for trees CSV file
            include_light_chain: Whether to include light chain columns for paired data

        Example:
            data.to_pcp("output_pcp.csv", "output_trees.csv")
            data.to_pcp("output.csv.gz", "trees.csv.gz")  # Compressed
        """
        pcp_path = Path(pcp_csv)
        is_paired = self.is_paired and include_light_chain

        # Build PCP rows from trees
        pcp_rows = []

        for tree in self.trees:
            clone_id = tree.get("clone_id", "")
            nodes = tree.get("nodes", [])

            # Handle both list and dict node formats
            if isinstance(nodes, dict):
                nodes_list = list(nodes.values())
            else:
                nodes_list = nodes

            # Build node lookup
            nodes_by_id = {n.get("sequence_id", ""): n for n in nodes_list}

            # Find sample_id from clone
            sample_id = ""
            for clone_list in self.clones.values():
                for clone in clone_list:
                    if clone.get("clone_id") == clone_id:
                        sample_id = clone.get("sample_id", "")
                        break

            # Build parent-child rows from nodes
            for node in nodes_list:
                node_id = node.get("sequence_id", "")
                parent_id = node.get("parent")

                if parent_id is None:
                    continue  # Skip root node (no parent)

                parent_node = nodes_by_id.get(parent_id, {})

                row = {
                    "sample_id": sample_id,
                    "family": clone_id,
                    "parent_name": parent_id,
                    "child_name": node_id,
                    "parent_heavy": parent_node.get("sequence_alignment", ""),
                    "child_heavy": node.get("sequence_alignment", ""),
                    "branch_length": node.get("length", 0.0),
                    "distance": node.get("distance", 0.0),
                    "parent_is_naive": str(parent_node.get("type") == "root"),
                    "child_is_leaf": str(node.get("type") == "leaf"),
                }

                # Add light chain columns for paired data
                if is_paired:
                    row["parent_light"] = parent_node.get("sequence_alignment_light", "")
                    row["child_light"] = node.get("sequence_alignment_light", "")

                pcp_rows.append(row)

        # Determine columns
        base_columns = [
            "sample_id",
            "family",
            "parent_name",
            "child_name",
            "parent_heavy",
            "child_heavy",
            "branch_length",
            "distance",
            "parent_is_naive",
            "child_is_leaf",
        ]
        if is_paired:
            base_columns.extend(["parent_light", "child_light"])

        # Write PCP CSV
        if pcp_path.suffix == ".gz":
            with gzip.open(pcp_path, "wt", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=base_columns)
                writer.writeheader()
                writer.writerows(pcp_rows)
        else:
            with open(pcp_path, "w", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=base_columns)
                writer.writeheader()
                writer.writerows(pcp_rows)

        # Write trees CSV if requested
        if trees_csv:
            trees_path = Path(trees_csv)
            tree_rows = []

            for tree in self.trees:
                clone_id = tree.get("clone_id", "")
                newick = tree.get("newick", "")

                # Find sample_id and rate scaling from clone
                sample_id = ""
                rate_scale_heavy = 1.0
                rate_scale_light = 1.0

                for clone_list in self.clones.values():
                    for clone in clone_list:
                        if clone.get("clone_id") == clone_id:
                            sample_id = clone.get("sample_id", "")
                            rate_scale_heavy = clone.get("rate_scale_heavy", 1.0)
                            rate_scale_light = clone.get("rate_scale_light", 1.0)
                            break

                row = {
                    "sample_id": sample_id,
                    "family": clone_id,
                    "newick": newick,
                }

                if is_paired:
                    row["rate_scale_heavy"] = rate_scale_heavy
                    row["rate_scale_light"] = rate_scale_light

                tree_rows.append(row)

            # Determine tree columns
            tree_columns = ["sample_id", "family", "newick"]
            if is_paired:
                tree_columns = [
                    "sample_id",
                    "family",
                    "rate_scale_heavy",
                    "rate_scale_light",
                    "newick",
                ]

            if trees_path.suffix == ".gz":
                with gzip.open(trees_path, "wt", newline="") as f:
                    writer = csv.DictWriter(f, fieldnames=tree_columns)
                    writer.writeheader()
                    writer.writerows(tree_rows)
            else:
                with open(trees_path, "w", newline="") as f:
                    writer = csv.DictWriter(f, fieldnames=tree_columns)
                    writer.writeheader()
                    writer.writerows(tree_rows)

    def to_airr_json(
        self,
        filepath: Union[str, Path],
        *,
        indent: Optional[int] = 2,
    ) -> None:
        """
        Export data to AIRR-compatible JSON format.

        This exports data in a format compatible with AIRR standards,
        suitable for use with other AIRR-compliant tools.

        Args:
            filepath: Output file path (use .gz extension for compression)
            indent: JSON indentation (None for compact output)

        Example:
            data.to_airr_json("airr_output.json")

        Note:
            This produces a simplified AIRR-compatible format. For full
            AIRR compliance, additional metadata may be required.
        """
        filepath = Path(filepath)

        # Build AIRR-compatible structure
        airr_output = {
            "Info": {
                "title": "Olmsted Export",
                "version": "1.0",
                "description": "Data exported from Olmsted",
            },
            "DataProcessing": [],
            "Repertoire": [],
            "GermlineSet": [],
            "Clone": [],
            "Tree": [],
        }

        # Convert clones to AIRR Clone format
        for dataset_id, clone_list in self.clones.items():
            for clone in clone_list:
                airr_clone = {
                    "clone_id": clone.get("clone_id"),
                    "repertoire_id": dataset_id,
                    "data_processing_id": None,
                    "sequences": clone.get("unique_seqs_count", 0),
                    "v_call": clone.get("v_call"),
                    "d_call": clone.get("d_call"),
                    "j_call": clone.get("j_call"),
                    "junction": clone.get("cdr3_sequence"),
                    "junction_aa": None,
                    "junction_length": clone.get("junction_length"),
                }

                # Add trees inline if present
                clone_trees = clone.get("trees", [])
                if clone_trees:
                    airr_clone["trees"] = []
                    for tree in clone_trees:
                        airr_tree = {
                            "tree_id": tree.get("tree_id"),
                            "clone_id": clone.get("clone_id"),
                            "newick": tree.get("newick"),
                        }
                        airr_clone["trees"].append(airr_tree)

                airr_output["Clone"].append(airr_clone)

        # Add standalone trees
        for tree in self.trees:
            airr_tree = {
                "tree_id": tree.get("tree_id"),
                "clone_id": tree.get("clone_id"),
                "newick": tree.get("newick"),
                "nodes": tree.get("nodes", []),
            }
            airr_output["Tree"].append(airr_tree)

        # Write output
        if filepath.suffix == ".gz":
            with gzip.open(filepath, "wt") as f:
                json.dump(airr_output, f, indent=indent)
        else:
            with open(filepath, "w") as f:
                json.dump(airr_output, f, indent=indent)

    # -------------------------------------------------------------------------
    # Convenience properties
    # -------------------------------------------------------------------------

    @property
    def dataset_ids(self) -> List[str]:
        """Get list of all dataset IDs."""
        return [d["dataset_id"] for d in self.datasets if "dataset_id" in d]

    @property
    def clone_count(self) -> int:
        """Get total number of clones across all datasets."""
        return sum(len(clone_list) for clone_list in self.clones.values())

    @property
    def tree_count(self) -> int:
        """Get total number of trees."""
        return len(self.trees)

    @property
    def is_paired(self) -> bool:
        """Check if data contains paired heavy/light chain information."""
        for clone_list in self.clones.values():
            for clone in clone_list:
                if clone.get("is_paired"):
                    return True
        return False

    def get_clones(self, dataset_id: Optional[str] = None) -> List[OlmstedClone]:
        """
        Get clones, optionally filtered by dataset.

        Args:
            dataset_id: Optional dataset ID to filter by

        Returns:
            List of clones
        """
        if dataset_id:
            return self.clones.get(dataset_id, [])
        # Return all clones flattened
        return [clone for clone_list in self.clones.values() for clone in clone_list]

    def get_trees(self, clone_id: Optional[str] = None) -> List[OlmstedTree]:
        """
        Get trees, optionally filtered by clone.

        Args:
            clone_id: Optional clone ID to filter by

        Returns:
            List of trees
        """
        if clone_id:
            return [t for t in self.trees if t.get("clone_id") == clone_id]
        return self.trees

    def __repr__(self) -> str:
        paired_str = " (paired)" if self.is_paired else ""
        return (
            f"OlmstedData("
            f"datasets={len(self.datasets)}, "
            f"clones={self.clone_count}, "
            f"trees={self.tree_count}"
            f"{paired_str})"
        )


__all__ = ["OlmstedData"]
