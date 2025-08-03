#!/usr/bin/env python3
"""Tests for olmsted-cli data processing using pytest."""

import json
import os
import shutil
import subprocess
import tempfile
import datetime
from pathlib import Path

import pytest


def normalize_json(obj):
    """Recursively sort dictionaries by key to enable comparison."""
    if isinstance(obj, dict):
        return {k: normalize_json(v) for k, v in sorted(obj.items())}
    elif isinstance(obj, list):
        return [normalize_json(item) for item in obj]
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
    """Compare all JSON files in two directories."""
    files1 = set(f for f in os.listdir(dir1) if f.endswith('.json'))
    files2 = set(f for f in os.listdir(dir2) if f.endswith('.json'))
    
    if files1 != files2:
        return False, f"Different files: {files1} vs {files2}"
    
    mismatches = []
    for fname in sorted(files1):
        file1 = os.path.join(dir1, fname)
        file2 = os.path.join(dir2, fname)
        if not compare_json_files(file1, file2):
            mismatches.append(fname)
    
    if mismatches:
        return False, f"Mismatched files: {mismatches}"
    
    return True, "All files match"


class TestOlmstedCLI:
    """Test suite for olmsted-cli."""
    
    @pytest.fixture(autouse=True)
    def setup_and_teardown(self):
        """Set up and tear down test environment."""
        # Get paths relative to the package root
        self.cli_root = Path(__file__).parent.parent
        self.test_data_dir = self.cli_root / "example_data"
        self.golden_airr_dir = self.test_data_dir / "airr" / "golden_airr_data"
        self.golden_pcp_dir = self.test_data_dir / "pcp" / "golden_pcp_data"
        
        # Create test output directory
        self.test_output_root = self.cli_root / "_test_output"
        self.test_output_root.mkdir(exist_ok=True)
        
        # Create a unique subdirectory for this test run
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        self.temp_dir = self.test_output_root / f"run_{timestamp}_{os.getpid()}"
        self.temp_dir.mkdir(exist_ok=True)
        
        yield
        
        # Optionally clean up old test runs (keep last 10)
        test_dirs = sorted(self.test_output_root.glob("run_*"))
        if len(test_dirs) > 10:
            for old_dir in test_dirs[:-10]:
                shutil.rmtree(old_dir, ignore_errors=True)
    
    
    @pytest.mark.airr
    def test_airr_processing_with_subcommand(self):
        """Test AIRR data processing using olmsted airr subcommand."""
        # Input and output paths
        input_file = self.test_data_dir / "airr" / "full_schema_dataset.json"
        output_dir = Path(self.temp_dir) / "airr_subcommand_output"
        
        # Run the subcommand
        cmd = [
            "olmsted", "airr",
            "-i", str(input_file),
            "-o", str(output_dir),
            "--seed", "42"
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        # Check command succeeded
        assert result.returncode == 0, f"Command failed: {result.stderr}"
        
        # Compare output with golden data
        match, message = compare_directories(str(self.golden_airr_dir), str(output_dir))
        assert match, f"Output doesn't match golden data: {message}"
    
    
    @pytest.mark.pcp
    def test_pcp_processing_with_subcommand(self):
        """Test PCP data processing using olmsted pcp subcommand."""
        # Input and output paths
        input_clones = self.test_data_dir / "pcp" / "clones.csv"
        input_trees = self.test_data_dir / "pcp" / "trees.csv"
        output_dir = Path(self.temp_dir) / "pcp_subcommand_output"
        
        # Run the subcommand
        cmd = [
            "olmsted", "pcp",
            "-i", str(input_clones),
            "-t", str(input_trees),
            "-o", str(output_dir),
            "--seed", "42"
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        # Check command succeeded
        assert result.returncode == 0, f"Command failed: {result.stderr}"
        
        # Compare output with golden data
        match, message = compare_directories(str(self.golden_pcp_dir), str(output_dir))
        assert match, f"Output doesn't match golden data: {message}"
    
    def test_auto_format_detection_airr(self):
        """Test automatic format detection for AIRR JSON files."""
        # Input and output paths
        input_file = self.test_data_dir / "airr" / "full_schema_dataset.json"
        output_dir = Path(self.temp_dir) / "auto_airr_output"
        
        # Run without specifying format using subcommand
        cmd = [
            "olmsted", "process",
            "-i", str(input_file),
            "-o", str(output_dir),
            "-f", "auto",  # Explicit auto-detection
            "--seed", "42"
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        # Check command succeeded
        assert result.returncode == 0, f"Command failed: {result.stderr}"
        
        # Verify it detected AIRR format
        assert "airr" in result.stdout.lower() or len(list(output_dir.glob("*.json"))) > 0
    
    def test_auto_format_detection_pcp(self):
        """Test automatic format detection for PCP CSV files."""
        # Input and output paths
        input_clones = self.test_data_dir / "pcp" / "clones.csv"
        input_trees = self.test_data_dir / "pcp" / "trees.csv"
        output_dir = Path(self.temp_dir) / "auto_pcp_output"
        
        # Run without specifying format using subcommand
        cmd = [
            "olmsted", "process",
            "-i", str(input_clones), str(input_trees),
            "-o", str(output_dir),
            "-f", "auto",  # Explicit auto-detection
            "--seed", "42"
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        # Check command succeeded
        assert result.returncode == 0, f"Command failed: {result.stderr}"
        
        # Verify it detected PCP format
        assert "pcp" in result.stdout.lower() or len(list(output_dir.glob("*.json"))) > 0
    
    def test_invalid_input_file(self):
        """Test handling of invalid input file."""
        output_dir = Path(self.temp_dir) / "invalid_output"
        
        cmd = [
            "olmsted", "process",
            "-i", "nonexistent_file.json",
            "-o", str(output_dir)
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        # Should fail
        assert result.returncode != 0
    
    def test_help_commands(self):
        """Test that help commands work."""
        help_commands = [
            ["olmsted", "--help"],
            ["olmsted", "process", "--help"],
            ["olmsted", "airr", "--help"],
            ["olmsted", "pcp", "--help"]
        ]
        
        for cmd in help_commands:
            result = subprocess.run(cmd, capture_output=True, text=True)
            assert result.returncode == 0, f"Help command failed: {cmd}"
            assert "help" in result.stdout.lower() or "usage" in result.stdout.lower()


if __name__ == "__main__":
    pytest.main([__file__, "-v"])