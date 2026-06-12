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


# Metadata fields that record run state, not data content. Stripped before
# comparing two consolidated outputs so the comparison is sensitive to data
# shape and not to invocation-time flags.
VOLATILE_METADATA_FIELDS = {"created_at", "git_hash", "version", "processing_options"}


def strip_volatile_fields(data):
    """Remove run-state metadata from consolidated data for comparison."""
    data = json.loads(json.dumps(data))  # deep copy
    if "metadata" in data:
        for key in VOLATILE_METADATA_FIELDS:
            data["metadata"].pop(key, None)
        generated_by = data["metadata"].get("generated_by", {})
        for key in VOLATILE_METADATA_FIELDS:
            generated_by.pop(key, None)
    return data


def compare_consolidated_files(file1, file2):
    """Compare two consolidated JSON files, ignoring volatile metadata fields."""
    with open(file1) as f:
        data1 = json.load(f)
    with open(file2) as f:
        data2 = json.load(f)

    norm1 = normalize_json(strip_volatile_fields(data1))
    norm2 = normalize_json(strip_volatile_fields(data2))

    if norm1 == norm2:
        return True, "Files match"

    # Find differences for error reporting
    diff = format_json_diff(file1, file2)
    return False, diff


def compare_json_files(file1, file2):
    """Compare two JSON files after stripping volatile fields and normalizing."""
    with open(file1) as f:
        data1 = json.load(f)
    with open(file2) as f:
        data2 = json.load(f)

    return normalize_json(strip_volatile_fields(data1)) == normalize_json(
        strip_volatile_fields(data2)
    )


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

    all_differences = []
    for fname in sorted(files1):
        file1 = os.path.join(dir1, fname)
        file2 = os.path.join(dir2, fname)
        if not compare_json_files(file1, file2):
            all_differences.append(format_json_diff(file1, file2))

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
        self.test_data_dir = self.cli_root / "example-data"
        self.golden_airr_dir = self.test_data_dir / "airr" / "split-golden-data"
        self.golden_pcp_dir = self.test_data_dir / "pcp" / "split-golden-data"
        self.consolidated_airr_file = (
            self.test_data_dir / "airr" / "airr-olmsted-golden.json"
        )
        self.consolidated_pcp_file = (
            self.test_data_dir / "pcp" / "pcp-olmsted-golden.json"
        )

        # Use the session directory and create a subdirectory for this specific test
        test_name = request.node.name
        self.temp_dir = test_session_dir / test_name
        self.temp_dir.mkdir(parents=True, exist_ok=True)

        # Store json_assertions for use in tests
        self.json_assertions = json_assertions

        yield

    @pytest.mark.airr
    def test_airr_processing(self):
        """`process --split-files` on AIRR matches the split-format golden.

        Pure shape check on the legacy split-file output (datasets.json,
        clones.*.json, tree.*.json). Validation is covered separately;
        this test fails specifically when split-format output drifts.
        """
        input_file = self.test_data_dir / "airr" / "input-airr.json"
        output_dir = Path(self.temp_dir) / "airr_output"

        cmd = [
            "olmsted", "process", "-f", "airr",
            "-i", str(input_file),
            "--split-files", str(output_dir),
            "--seed", "42",
            "--name", "airr-example",

            "--json-format", "pretty",
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)
        assert result.returncode == 0, f"Command failed: {result.stderr}"

        match, message = compare_directories(str(self.golden_airr_dir), str(output_dir))
        assert match, f"Output doesn't match golden data: {message}"

    @pytest.mark.airr
    def test_airr_consolidated_processing(self):
        """`process` on AIRR produces output matching the golden (data shape).

        Pure shape check — no `--validate`. Schema/time-tree validation is
        covered separately by `test_airr_consolidated_inline_validation`
        and `test_validate_airr_consolidated_golden_output`, so a failure
        here unambiguously points at processing rather than validation.
        """
        input_file = self.test_data_dir / "airr" / "input-airr.json"
        output_file = Path(self.temp_dir) / "airr_consolidated.json"

        cmd = [
            "olmsted", "process", "-f", "airr",
            "-i", str(input_file),
            "-o", str(output_file),
            "--seed", "42",
            "--name", "airr-example",

            "--json-format", "pretty",
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)
        assert result.returncode == 0, f"Command failed: {result.stderr}"
        assert output_file.exists(), f"Output file not created: {output_file}"

        # Compare against consolidated golden data (ignoring volatile fields)
        match, message = compare_consolidated_files(
            str(self.consolidated_airr_file), str(output_file)
        )
        assert match, f"Output doesn't match consolidated golden data:\n{message}"

    @pytest.mark.airr
    def test_airr_consolidated_inline_validation(self):
        """`process --validate` on AIRR runs without error.

        Exercises the inline-validation code path through the CLI. A
        regression here means either inline validation broke or the
        AIRR output became schema/time-tree-invalid.
        """
        input_file = self.test_data_dir / "airr" / "input-airr.json"
        output_file = Path(self.temp_dir) / "airr_inline_validated.json"

        cmd = [
            "olmsted", "process", "-f", "airr",
            "-i", str(input_file),
            "-o", str(output_file),
            "--seed", "42",
            "--name", "airr-example",

            "--json-format", "pretty",
            "--validate",
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)
        assert result.returncode == 0, (
            f"`process --validate` failed:\n"
            f"stdout: {result.stdout}\nstderr: {result.stderr}"
        )

    @pytest.mark.pcp
    def test_pcp_with_trees_processing(self):
        """Test PCP data processing with separate trees file using --tree argument."""
        # Input and output paths
        input_clones = self.test_data_dir / "pcp" / "input-pcp.csv"
        input_trees = self.test_data_dir / "pcp" / "input-trees.csv"
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
        input_clones = self.test_data_dir / "pcp" / "input-pcp.csv"
        input_trees = self.test_data_dir / "pcp" / "input-trees.csv"
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
    def test_pcp_processing(self):
        """`process --split-files` on PCP matches the split-format golden.

        Pure shape check on the legacy split-file output. See the AIRR
        counterpart for the rationale on splitting shape from
        validation coverage.
        """
        input_clones = self.test_data_dir / "pcp" / "input-pcp.csv"
        input_trees = self.test_data_dir / "pcp" / "input-trees.csv"
        output_dir = Path(self.temp_dir) / "pcp_output"

        cmd = [
            "olmsted", "process", "-f", "pcp",
            "-i", str(input_clones),
            "-t", str(input_trees),
            "--split-files", str(output_dir),
            "--seed", "42",
            "--name", "pcp-example",

            "--json-format", "pretty",
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)
        assert result.returncode == 0, f"Command failed: {result.stderr}"

        match, message = compare_directories(str(self.golden_pcp_dir), str(output_dir))
        assert match, f"Output doesn't match golden data: {message}"

    @pytest.mark.pcp
    def test_pcp_consolidated_processing(self):
        """`process` on PCP produces output matching the golden (data shape).

        Pure shape check — no `--validate`. See the AIRR counterpart for
        the rationale on splitting shape from inline-validation coverage.
        """
        input_clones = self.test_data_dir / "pcp" / "input-pcp.csv"
        input_trees = self.test_data_dir / "pcp" / "input-trees.csv"
        output_file = Path(self.temp_dir) / "pcp_consolidated.json"

        cmd = [
            "olmsted", "process", "-f", "pcp",
            "-i", str(input_clones),
            "-t", str(input_trees),
            "-o", str(output_file),
            "--seed", "42",
            "--name", "pcp-example",

            "--json-format", "pretty",
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)
        assert result.returncode == 0, f"Command failed: {result.stderr}"
        assert output_file.exists(), f"Output file not created: {output_file}"

        # Compare against consolidated golden data (ignoring volatile fields)
        match, message = compare_consolidated_files(
            str(self.consolidated_pcp_file), str(output_file)
        )
        assert match, f"Output doesn't match consolidated golden data:\n{message}"

    @pytest.mark.pcp
    def test_pcp_consolidated_inline_validation(self):
        """`process --validate` on PCP runs without error.

        Inline validation during processing does not run the time-tree
        check that `validate_file(check_time_tree=True)` does, so this
        test passes despite the known PCP time-tree drift (see issue
        #20 and `test_validate_pcp_consolidated_golden_output`).
        """
        input_clones = self.test_data_dir / "pcp" / "input-pcp.csv"
        input_trees = self.test_data_dir / "pcp" / "input-trees.csv"
        output_file = Path(self.temp_dir) / "pcp_inline_validated.json"

        cmd = [
            "olmsted", "process", "-f", "pcp",
            "-i", str(input_clones),
            "-t", str(input_trees),
            "-o", str(output_file),
            "--seed", "42",
            "--name", "pcp-example",

            "--json-format", "pretty",
            "--validate",
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)
        assert result.returncode == 0, (
            f"`process --validate` failed:\n"
            f"stdout: {result.stdout}\nstderr: {result.stderr}"
        )

    def test_auto_format_detection_airr(self):
        """Test automatic format detection for AIRR JSON files."""
        # Input and output paths
        input_file = self.test_data_dir / "airr" / "input-airr.json"
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
        input_clones = self.test_data_dir / "pcp" / "input-pcp.csv"
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
        """Test processing of paired heavy/light chain PCP data using two-clone architecture."""
        # Use test paired PCP data subset
        paired_pcp_dir = self.test_data_dir / "pcp-paired"
        input_clones = paired_pcp_dir / "input-pcp.csv"
        input_trees = paired_pcp_dir / "input-trees.csv"

        if not input_clones.exists():
            pytest.skip("Paired PCP test data not available")

        output_file = Path(self.temp_dir) / "paired_pcp_output.json"

        # Import processing functions directly for faster testing
        from olmsted_cli.process_pcp_data import (
            parse_pcp_csv,
            parse_newick_csv,
            process_pcp_to_olmsted,
        )

        # Parse the test data
        trees = parse_newick_csv(str(input_trees))
        families = parse_pcp_csv(str(input_clones))

        # Process
        datasets, clones_dict, trees_out = process_pcp_to_olmsted(
            families, trees, verbosity=0
        )

        # Verify two-clone architecture for paired data
        dataset_id = list(clones_dict.keys())[0]
        all_clones = clones_dict[dataset_id]

        # With the two-clone architecture, each paired family creates TWO clones
        # Find a pair by looking for clones with matching pair_id
        paired_clones = [c for c in all_clones if c.get("is_paired")]
        assert len(paired_clones) > 0, "Should have paired clones"

        # Group clones by pair_id
        from collections import defaultdict
        pairs = defaultdict(list)
        for clone in paired_clones:
            pair_id = clone.get("pair_id")
            if pair_id:
                pairs[pair_id].append(clone)

        # Verify we have complete pairs (heavy + light)
        assert len(pairs) > 0, "Should have at least one pair"
        first_pair_id = list(pairs.keys())[0]
        pair = pairs[first_pair_id]
        assert len(pair) == 2, f"Each pair should have exactly 2 clones (heavy + light), got {len(pair)}"

        # Find heavy and light clones
        heavy_clone = next((c for c in pair if c["sample"]["locus"] == "igh"), None)
        light_clone = next((c for c in pair if c["sample"]["locus"] in ["igk", "igl"]), None)

        assert heavy_clone is not None, "Should have heavy chain clone"
        assert light_clone is not None, "Should have light chain clone"

        # Verify both clones are marked as paired with matching pair_id
        assert heavy_clone.get("is_paired") is True, "Heavy clone should be marked as paired"
        assert light_clone.get("is_paired") is True, "Light clone should be marked as paired"
        assert heavy_clone.get("pair_id") == light_clone.get("pair_id"), \
            "Heavy and light clones should have matching pair_id"

        # Verify locus is correct
        assert heavy_clone["sample"]["locus"] == "igh", "Heavy clone should have locus igh"
        assert light_clone["sample"]["locus"] in ["igk", "igl"], \
            f"Light clone should have locus igk or igl, got {light_clone['sample']['locus']}"

        # Verify VDJ genes are different between heavy and light
        assert heavy_clone.get("v_call") != light_clone.get("v_call"), \
            "Heavy and light clones should have different V genes"
        assert heavy_clone.get("j_call") != light_clone.get("j_call"), \
            "Heavy and light clones should have different J genes"

        # Light chain should not have D gene (heavy chain may or may not have D gene depending on input data)
        assert light_clone.get("d_call", "") == "", "Light clone should not have D gene"

        # Verify trees exist for both chains
        heavy_tree = next((t for t in trees_out if t["clone_id"] == heavy_clone["clone_id"]), None)
        light_tree = next((t for t in trees_out if t["clone_id"] == light_clone["clone_id"]), None)

        assert heavy_tree is not None, "Should have tree for heavy chain"
        assert light_tree is not None, "Should have tree for light chain"

        # Verify both trees have nodes
        assert len(heavy_tree["nodes"]) > 0, "Heavy tree should have nodes"
        assert len(light_tree["nodes"]) > 0, "Light tree should have nodes"

        # Verify sequences are different between heavy and light trees
        heavy_leaf = next((n for n in heavy_tree["nodes"] if n.get("type") == "leaf"), None)
        light_leaf = next((n for n in light_tree["nodes"] if n.get("type") == "leaf"), None)

        if heavy_leaf and light_leaf:
            assert heavy_leaf.get("sequence_alignment") != light_leaf.get("sequence_alignment"), \
                "Heavy and light trees should have different sequences"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
