#!/usr/bin/env python3
"""Tests for olmsted-cli data processing using pytest."""

import json
import logging
import os
import subprocess
from pathlib import Path

import pytest

from .conftest import format_json_diff

# Set up logging for tests
logger = logging.getLogger(__name__)


def normalize_json(obj, float_tolerance=1e-12):
    """Recursively sort dictionaries by key and normalize floats for comparison."""
    if isinstance(obj, dict):
        return {k: normalize_json(v, float_tolerance) for k, v in sorted(obj.items())}
    elif isinstance(obj, list):
        return [normalize_json(item, float_tolerance) for item in obj]
    elif isinstance(obj, float):
        # Round very small numbers to avoid precision differences
        if abs(obj) < float_tolerance:
            return 0.0
        # Round to 15 significant digits to handle floating-point precision
        return round(obj, 15)
    else:
        return obj


def compare_json_files(file1, file2):
    """Compare two JSON files after normalizing."""
    with open(file1) as f:
        data1 = json.load(f)
    with open(file2) as f:
        data2 = json.load(f)

    norm1 = normalize_json(data1)
    norm2 = normalize_json(data2)

    return norm1 == norm2


def compare_directories(dir1, dir2):
    """Compare all JSON files in two directories with detailed error reporting."""
    files1 = set(f for f in os.listdir(dir1) if f.endswith(".json"))
    files2 = set(f for f in os.listdir(dir2) if f.endswith(".json"))

    if files1 != files2:
        missing_in_dir2 = files1 - files2
        missing_in_dir1 = files2 - files1
        error_msg = []
        if missing_in_dir2:
            error_msg.append(f"Files missing in {dir2}: {missing_in_dir2}")
        if missing_in_dir1:
            error_msg.append(f"Files missing in {dir1}: {missing_in_dir1}")
        return False, "\n".join(error_msg)

    # Check each file and collect detailed differences
    all_differences = []
    for fname in sorted(files1):
        file1 = os.path.join(dir1, fname)
        file2 = os.path.join(dir2, fname)
        if not compare_json_files(file1, file2):
            diff_output = format_json_diff(file1, file2)
            all_differences.append(diff_output)

    if all_differences:
        return False, "\n\n".join(all_differences)

    return True, "All files match"


class TestOlmstedCLI:
    """Test suite for olmsted-cli."""

    @pytest.fixture(autouse=True)
    def setup_and_teardown(self, test_session_dir, request, json_assertions):
        """Set up and tear down test environment."""
        # Get paths relative to the package root
        self.cli_root = Path(__file__).parent.parent
        self.test_data_dir = self.cli_root / "example_data"
        self.golden_airr_dir = self.test_data_dir / "airr" / "split_golden_data"
        self.golden_pcp_dir = self.test_data_dir / "pcp" / "split_golden_data"
        self.consolidated_airr_file = (
            self.test_data_dir / "airr" / "consolidated_golden_data.json"
        )
        self.consolidated_pcp_file = (
            self.test_data_dir / "pcp" / "consolidated_golden_data.json"
        )

        # Use the session directory and create a subdirectory for this specific test
        test_name = request.node.name
        self.temp_dir = test_session_dir / test_name
        self.temp_dir.mkdir(exist_ok=True)

        # Store json_assertions for use in tests
        self.json_assertions = json_assertions

        yield

    @pytest.mark.airr
    def test_airr_processing(self):
        """Test AIRR data processing using olmsted process command with split files."""
        # Input and output paths
        input_file = self.test_data_dir / "airr" / "full_schema_dataset.json"
        output_dir = Path(self.temp_dir) / "airr_output"

        # Run the process command with split files to match golden data
        cmd = [
            "olmsted",
            "process",
            "-f",
            "airr",
            "-i",
            str(input_file),
            "--split-files",
            str(output_dir),
            "--seed",
            "42",
            "--name",
            "airr-example",
            "--validate",
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)

        # Check command succeeded
        assert result.returncode == 0, f"Command failed: {result.stderr}"

        # Compare output with golden data
        match, message = compare_directories(str(self.golden_airr_dir), str(output_dir))
        assert match, f"Output doesn't match golden data: {message}"

    @pytest.mark.airr
    def test_airr_consolidated_processing(self):
        """Test AIRR data processing using consolidated output format."""
        # Input and output paths
        input_file = self.test_data_dir / "airr" / "full_schema_dataset.json"
        output_file = Path(self.temp_dir) / "airr_consolidated.json"

        # Run the process command with consolidated output
        cmd = [
            "olmsted",
            "process",
            "-f",
            "airr",
            "-i",
            str(input_file),
            "-o",
            str(output_file),
            "--seed",
            "42",
            "--name",
            "airr-example",
            "--validate",
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)

        # Check command succeeded
        assert result.returncode == 0, f"Command failed: {result.stderr}"

        # Verify output file exists
        assert output_file.exists(), f"Output file not created: {output_file}"

        # Load and verify structure
        with open(output_file) as f:
            data = json.load(f)

        # Verify consolidated structure
        assert "metadata" in data, "Consolidated data should have metadata"
        assert "datasets" in data, "Consolidated data should have datasets"
        assert "clones" in data, "Consolidated data should have clones"
        assert "trees" in data, "Consolidated data should have trees"

        # Verify metadata structure
        from olmsted_cli.process_utils import CONSOLIDATED_JSON_VERSION

        metadata = data["metadata"]
        assert metadata["format_version"] == CONSOLIDATED_JSON_VERSION, (
            "Should have correct format version"
        )
        assert metadata["source_format"] == "airr", "Should identify source format"
        assert "created_at" in metadata, "Should have creation timestamp"
        assert "processing_info" in metadata, "Should have processing info"
        assert metadata["name"] == "airr-example", "Should have correct name"

    @pytest.mark.pcp
    def test_pcp_processing(self):
        """Test PCP data processing using olmsted process command with split files."""
        # Input and output paths
        input_clones = self.test_data_dir / "pcp" / "pcp.csv"
        input_trees = self.test_data_dir / "pcp" / "trees.csv"
        output_dir = Path(self.temp_dir) / "pcp_output"

        # Run the process command with split files to match golden data
        cmd = [
            "olmsted",
            "process",
            "-f",
            "pcp",
            "-i",
            str(input_clones),
            "-t",
            str(input_trees),
            "--split-files",
            str(output_dir),
            "--seed",
            "42",
            "--name",
            "pcp-example",
            "--validate",
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)

        # Check command succeeded
        assert result.returncode == 0, f"Command failed: {result.stderr}"

        # Compare output with golden data
        match, message = compare_directories(str(self.golden_pcp_dir), str(output_dir))
        assert match, f"Output doesn't match golden data: {message}"

    @pytest.mark.pcp
    def test_pcp_with_trees_processing(self):
        """Test PCP data processing with separate trees file using --tree argument."""
        # Input and output paths
        input_clones = self.test_data_dir / "pcp" / "pcp.csv"
        input_trees = self.test_data_dir / "pcp" / "trees.csv"
        output_dir = Path(self.temp_dir) / "pcp_trees_output"

        # Run the process command with --tree argument
        cmd = [
            "olmsted",
            "process",
            "-f",
            "pcp",
            "-i",
            str(input_clones),
            "--tree",
            str(input_trees),
            "--split-files",
            str(output_dir),
            "--seed",
            "42",
            "--name",
            "pcp-with-trees",
            "--validate",
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)

        # Check command succeeded
        assert result.returncode == 0, f"Command failed: {result.stderr}"

        # Check that tree files were created
        tree_files = list(output_dir.glob("tree.*.json"))
        assert len(tree_files) > 0, "No tree files were created"

        # Verify that trees contain newick data
        for tree_file in tree_files:
            with open(tree_file) as f:
                tree_data = json.load(f)
                assert "newick" in tree_data, f"Tree file {tree_file} missing newick data"
                assert "nodes" in tree_data, f"Tree file {tree_file} missing nodes data"

    @pytest.mark.pcp
    def test_pcp_with_trees_short_option(self):
        """Test PCP data processing with -t short option for trees."""
        # Input and output paths
        input_clones = self.test_data_dir / "pcp" / "pcp.csv"
        input_trees = self.test_data_dir / "pcp" / "trees.csv"
        output_file = Path(self.temp_dir) / "pcp_with_trees.json"

        # Run the process command with -t short option
        cmd = [
            "olmsted",
            "process",
            "-f",
            "pcp",
            "-i",
            str(input_clones),
            "-t",
            str(input_trees),
            "-o",
            str(output_file),
            "--seed",
            "42",
            "--validate",
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)

        # Check command succeeded
        assert result.returncode == 0, f"Command failed: {result.stderr}"

        # Verify output file exists
        assert output_file.exists(), "Output file was not created"

        # Verify that trees are included in consolidated output
        with open(output_file) as f:
            data = json.load(f)
            assert "trees" in data, "Trees not found in consolidated output"
            assert len(data["trees"]) > 0, "No trees in consolidated output"
            # Check first tree has expected structure
            first_tree = data["trees"][0]
            assert "newick" in first_tree, "First tree missing newick data"
            assert "nodes" in first_tree, "First tree missing nodes data"

    @pytest.mark.pcp
    def test_pcp_consolidated_processing(self):
        """Test PCP data processing using consolidated output format."""
        # Input and output paths
        input_clones = self.test_data_dir / "pcp" / "pcp.csv"
        output_file = Path(self.temp_dir) / "pcp_consolidated.json"

        # Run the process command with consolidated output
        cmd = [
            "olmsted",
            "process",
            "-f",
            "pcp",
            "-i",
            str(input_clones),
            "-o",
            str(output_file),
            "--seed",
            "42",
            "--name",
            "pcp-example",
            "--validate",
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)

        # Check command succeeded
        assert result.returncode == 0, f"Command failed: {result.stderr}"

        # Verify output file exists
        assert output_file.exists(), f"Output file not created: {output_file}"

        # Load and verify structure
        with open(output_file) as f:
            data = json.load(f)

        # Verify consolidated structure
        assert "metadata" in data, "Consolidated data should have metadata"
        assert "datasets" in data, "Consolidated data should have datasets"
        assert "clones" in data, "Consolidated data should have clones"
        assert "trees" in data, "Consolidated data should have trees"

        # Verify metadata structure
        from olmsted_cli.process_utils import CONSOLIDATED_JSON_VERSION

        metadata = data["metadata"]
        assert metadata["format_version"] == CONSOLIDATED_JSON_VERSION, (
            "Should have correct format version"
        )
        assert metadata["source_format"] == "pcp", "Should identify source format"
        assert "created_at" in metadata, "Should have creation timestamp"
        assert "processing_info" in metadata, "Should have processing info"
        assert metadata["name"] == "pcp-example", "Should have correct name"

    def test_auto_format_detection_airr(self):
        """Test automatic format detection for AIRR JSON files."""
        # Input and output paths
        input_file = self.test_data_dir / "airr" / "full_schema_dataset.json"
        output_dir = Path(self.temp_dir) / "auto_airr_output"

        # Run without specifying format using subcommand
        cmd = [
            "olmsted",
            "process",
            "-i",
            str(input_file),
            "--split-files",
            str(output_dir),
            "-f",
            "auto",  # Explicit auto-detection
            "--seed",
            "42",
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)

        # Check command succeeded
        assert result.returncode == 0, f"Command failed: {result.stderr}"

        # Verify it detected AIRR format
        assert (
            "airr" in result.stdout.lower() or len(list(output_dir.glob("*.json"))) > 0
        )

    def test_auto_format_detection_pcp(self):
        """Test automatic format detection for PCP CSV files."""
        # Input and output paths
        input_clones = self.test_data_dir / "pcp" / "pcp.csv"
        output_dir = Path(self.temp_dir) / "auto_pcp_output"

        # Run without specifying format using subcommand
        cmd = [
            "olmsted",
            "process",
            "-i",
            str(input_clones),
            "--split-files",
            str(output_dir),
            "-f",
            "auto",  # Explicit auto-detection
            "--seed",
            "42",
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)

        # Check command succeeded
        assert result.returncode == 0, f"Command failed: {result.stderr}"

        # Verify it detected PCP format
        assert (
            "pcp" in result.stdout.lower() or len(list(output_dir.glob("*.json"))) > 0
        )

    def test_invalid_input_file(self):
        """Test handling of invalid input file."""
        output_dir = Path(self.temp_dir) / "invalid_output"

        cmd = [
            "olmsted",
            "process",
            "-i",
            "nonexistent_file.json",
            "--split-files",
            str(output_dir),
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)

        # Should fail
        assert result.returncode != 0

    def test_help_commands(self):
        """Test that help commands work."""
        logger.info("Testing help commands")
        help_commands = [["olmsted", "--help"], ["olmsted", "process", "--help"]]

        for cmd in help_commands:
            logger.info(f"Running command: {' '.join(cmd)}")
            result = subprocess.run(cmd, capture_output=True, text=True)
            logger.debug(f"Command output: {result.stdout[:200]}...")
            assert result.returncode == 0, f"Help command failed: {cmd}"
            assert "help" in result.stdout.lower() or "usage" in result.stdout.lower()
        logger.info("All help commands passed")

    @pytest.mark.pcp
    @pytest.mark.paired
    def test_paired_pcp_processing(self):
        """Test processing of paired heavy/light chain PCP data."""
        # Use paired PCP data (check if it exists, skip if not)
        paired_pcp_dir = self.test_data_dir / "pcp-paired"
        input_clones = paired_pcp_dir / "wyatt-10x-1p5m_fs-all-UnmutInv_paired-merged_pcp_2024-11-22.csv.gz"
        input_trees = paired_pcp_dir / "wyatt-10x-1p5m_fs-all-UnmutInv_paired-merged_trees_2024-11-22.csv.gz"

        if not input_clones.exists():
            pytest.skip("Paired PCP test data not available")

        output_file = Path(self.temp_dir) / "paired_pcp_output.json"

        # Import processing functions directly for faster testing with subset
        from olmsted_cli.process_pcp_data import (
            parse_pcp_csv,
            parse_newick_csv,
            process_pcp_to_olmsted,
        )

        # Parse a small subset for testing
        trees = parse_newick_csv(str(input_trees))
        families = parse_pcp_csv(str(input_clones))

        # Take just first 5 families for faster testing
        family_ids = list(families.keys())[:5]
        small_families = {k: families[k] for k in family_ids}

        # Process
        datasets, clones_dict, trees_out = process_pcp_to_olmsted(
            small_families, trees, verbosity=0
        )

        # Verify paired data properties
        dataset_id = list(clones_dict.keys())[0]
        first_clone = clones_dict[dataset_id][0]

        # Check that light chain fields are present
        assert first_clone.get("is_paired") is True, "Clone should be marked as paired"
        assert "v_call_light" in first_clone, "Clone should have v_call_light"
        assert "j_call_light" in first_clone, "Clone should have j_call_light"
        assert "light_chain_type" in first_clone, "Clone should have light_chain_type"
        assert first_clone.get("light_chain_type") in ["kappa", "lambda"], \
            f"light_chain_type should be kappa or lambda, got {first_clone.get('light_chain_type')}"

        # Check rate scaling
        assert "rate_scale_heavy" in first_clone, "Clone should have rate_scale_heavy"
        assert "rate_scale_light" in first_clone, "Clone should have rate_scale_light"

        # Check that nodes have light chain sequences
        first_tree = trees_out[0]
        assert len(first_tree["nodes"]) > 0, "Tree should have nodes"

        # Find a node with sequence data
        nodes_with_light_seq = [
            n for n in first_tree["nodes"]
            if n.get("sequence_alignment_light")
        ]
        assert len(nodes_with_light_seq) > 0, "At least one node should have light chain sequence"

        # Check light chain AA translation
        node_with_light = nodes_with_light_seq[0]
        assert "sequence_alignment_light_aa" in node_with_light, \
            "Node with light sequence should have AA translation"
        assert len(node_with_light["sequence_alignment_light_aa"]) > 0, \
            "Light chain AA sequence should not be empty"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
