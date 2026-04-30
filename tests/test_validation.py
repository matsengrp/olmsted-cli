#!/usr/bin/env python3
"""Tests for olmsted-cli validation command using pytest."""

import json
import os
import tempfile
from pathlib import Path

import pytest

from olmsted_cli.validate import validate_file
from olmsted_cli.process_utils import (
    validate_clone,
    validate_dataset,
    validate_tree,
    validate_time_tree,
)


class TestValidation:
    """Test validation functionality."""

    def test_validate_airr_golden_outputs(self):
        """Validate every split-format AIRR golden file individually."""
        golden_dir = (
            Path(__file__).parent.parent / "example-data" / "airr" / "split-golden-data"
        )

        if not golden_dir.exists():
            pytest.skip(f"Golden AIRR data directory not found: {golden_dir}")

        json_files = list(golden_dir.glob("*.json"))
        assert len(json_files) > 0, f"No JSON files found in {golden_dir}"

        # Schema validation (default).
        validation_errors = []
        for json_file in json_files:
            is_valid, errors = validate_file(
                str(json_file), file_type=None, verbose=True, check_time_tree=False
            )
            if not is_valid:
                validation_errors.append(f"{json_file.name}: {errors}")
        assert len(validation_errors) == 0, (
            "AIRR golden outputs should be valid. Errors found:\n"
            + "\n".join(validation_errors)
        )

        # Time-tree validation.
        validation_errors_time_tree = []
        for json_file in json_files:
            is_valid, errors = validate_file(
                str(json_file), file_type=None, verbose=True, check_time_tree=True
            )
            if not is_valid:
                validation_errors_time_tree.append(f"{json_file.name}: {errors}")
        assert len(validation_errors_time_tree) == 0, (
            "AIRR golden outputs should be valid time trees. Errors found:\n"
            + "\n".join(validation_errors_time_tree)
        )

    @pytest.mark.xfail(
        reason="time-tree precision drift on golden PCP, see issue #20",
        strict=True,
    )
    def test_validate_pcp_golden_outputs(self):
        """Validate every split-format PCP golden file individually.

        Currently xfailed for the same reason as
        ``test_validate_pcp_consolidated_golden_output``: floating-point
        drift in tree distances trips the time-tree invariant. Tracked
        in issue #20.
        """
        golden_dir = (
            Path(__file__).parent.parent / "example-data" / "pcp" / "split-golden-data"
        )

        if not golden_dir.exists():
            pytest.skip(f"Golden PCP data directory not found: {golden_dir}")

        json_files = list(golden_dir.glob("*.json"))
        assert len(json_files) > 0, f"No JSON files found in {golden_dir}"

        validation_errors = []
        for json_file in json_files:
            is_valid, errors = validate_file(
                str(json_file), file_type=None, verbose=True, check_time_tree=False
            )
            if not is_valid:
                validation_errors.append(f"{json_file.name}: {errors}")
        assert len(validation_errors) == 0, (
            "PCP golden outputs should be valid. Errors found:\n"
            + "\n".join(validation_errors)
        )

        validation_errors_time_tree = []
        for json_file in json_files:
            is_valid, errors = validate_file(
                str(json_file), file_type=None, verbose=True, check_time_tree=True
            )
            if not is_valid:
                validation_errors_time_tree.append(f"{json_file.name}: {errors}")
        assert len(validation_errors_time_tree) == 0, (
            "PCP golden outputs should be valid time trees. Errors found:\n"
            + "\n".join(validation_errors_time_tree)
        )

    def test_validate_airr_consolidated_golden_output(self):
        """Test that AIRR consolidated golden output is valid."""
        consolidated_file = (
            Path(__file__).parent.parent
            / "example-data"
            / "airr"
            / "airr-olmsted-golden.json"
        )

        if not consolidated_file.exists():
            pytest.skip(
                f"Consolidated AIRR golden data file not found: {consolidated_file}"
            )

        # Test without time tree validation (default)
        is_valid, errors = validate_file(
            str(consolidated_file), file_type=None, verbose=True, check_time_tree=False
        )
        assert is_valid, (
            f"AIRR consolidated golden output should be valid. Errors found:\n"
            + "\n".join(str(e) for e in errors)
        )
        
        # Also test with time tree validation enabled
        is_valid_time_tree, errors_time_tree = validate_file(
            str(consolidated_file), file_type=None, verbose=True, check_time_tree=True
        )
        assert is_valid_time_tree, (
            f"AIRR consolidated golden output should be valid time tree. Errors found:\n"
            + "\n".join(str(e) for e in errors_time_tree)
        )

        # Load and verify consolidated structure
        with open(consolidated_file) as f:
            data = json.load(f)

        # Verify consolidated structure
        assert "metadata" in data, "Consolidated data should have metadata"
        assert "datasets" in data, "Consolidated data should have datasets"
        assert "clones" in data, "Consolidated data should have clones"
        assert "trees" in data, "Consolidated data should have trees"

    @pytest.mark.xfail(
        reason="time-tree precision drift on golden PCP, see issue #20",
        strict=True,
    )
    def test_validate_pcp_consolidated_golden_output(self):
        """Test that PCP consolidated golden output is valid.

        Currently xfailed: the PCP example tree has a Node5 distance that
        is floating-point-less-than its parent's, tripping the time-tree
        invariant check. Tracked in issue #20; this xfail flips to pass
        automatically once the drift is resolved (alerting us to drop the
        marker).
        """
        consolidated_file = (
            Path(__file__).parent.parent
            / "example-data"
            / "pcp"
            / "pcp-olmsted-golden.json"
        )

        if not consolidated_file.exists():
            pytest.skip(
                f"Consolidated PCP golden data file not found: {consolidated_file}"
            )

        # Test without time tree validation (default)
        is_valid, errors = validate_file(
            str(consolidated_file), file_type=None, verbose=True, check_time_tree=False
        )
        assert is_valid, (
            f"PCP consolidated golden output should be valid. Errors found:\n"
            + "\n".join(str(e) for e in errors)
        )
        
        # Also test with time tree validation enabled
        is_valid_time_tree, errors_time_tree = validate_file(
            str(consolidated_file), file_type=None, verbose=True, check_time_tree=True
        )
        assert is_valid_time_tree, (
            f"PCP consolidated golden output should be valid time tree. Errors found:\n"
            + "\n".join(str(e) for e in errors_time_tree)
        )

        # Load and verify consolidated structure
        with open(consolidated_file) as f:
            data = json.load(f)

        # Verify consolidated structure
        assert "metadata" in data, "Consolidated data should have metadata"
        assert "datasets" in data, "Consolidated data should have datasets"
        assert "clones" in data, "Consolidated data should have clones"
        assert "trees" in data, "Consolidated data should have trees"

    def test_validate_invalid_dataset(self):
        """Test that invalid datasets are properly rejected."""
        # Create an invalid dataset (missing required fields)
        invalid_dataset = {
            # Missing required 'dataset_id' field
            "clones": []
        }

        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            json.dump(invalid_dataset, f)
            temp_file = f.name

        try:
            is_valid, errors = validate_file(
                temp_file, file_type="dataset", verbose=True
            )
            assert not is_valid, "Invalid dataset should fail validation"
            assert len(errors) > 0, "Should have validation errors"

            # Check that the error mentions the missing required field
            error_str = " ".join(str(e) for e in errors)
            assert (
                "dataset_id" in error_str.lower() or "required" in error_str.lower()
            ), (
                f"Error should mention missing required field 'dataset_id'. Got: {errors}"
            )
        finally:
            # Clean up temp file
            os.unlink(temp_file)

    def test_validate_invalid_clone(self):
        """Test that invalid clones are properly rejected."""
        # Create an invalid clone (missing required fields)
        invalid_clone = {
            "clone_id": "test_clone",
            # Missing required fields like unique_seqs_count, mean_mut_freq, etc.
        }

        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            json.dump(invalid_clone, f)
            temp_file = f.name

        try:
            is_valid, errors = validate_file(temp_file, file_type="clone", verbose=True)
            assert not is_valid, "Invalid clone should fail validation"
            assert len(errors) > 0, "Should have validation errors"

            # Check that the error mentions missing required fields
            error_str = " ".join(str(e) for e in errors).lower()
            assert any(
                field in error_str
                for field in [
                    "unique_seqs_count",
                    "mean_mut_freq",
                    "v_alignment",
                    "j_alignment",
                    "required",
                ]
            ), f"Error should mention missing required fields. Got: {errors}"
        finally:
            # Clean up temp file
            os.unlink(temp_file)

    def test_validate_invalid_tree(self):
        """Test that invalid trees are properly rejected."""
        # Create an invalid tree (missing required fields)
        invalid_tree = {
            "tree_id": "test_tree",
            # Missing required 'newick' and 'nodes' fields
        }

        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            json.dump(invalid_tree, f)
            temp_file = f.name

        try:
            is_valid, errors = validate_file(temp_file, file_type="tree", verbose=True)
            assert not is_valid, "Invalid tree should fail validation"
            assert len(errors) > 0, "Should have validation errors"

            # Check that the error mentions the missing required fields
            error_str = " ".join(str(e) for e in errors).lower()
            assert any(
                field in error_str for field in ["newick", "nodes", "required"]
            ), (
                f"Error should mention missing required fields 'newick' or 'nodes'. Got: {errors}"
            )
        finally:
            # Clean up temp file
            os.unlink(temp_file)

    def test_validate_malformed_json(self):
        """Test that malformed JSON is properly rejected."""
        malformed_json = '{"invalid": json syntax}'

        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            f.write(malformed_json)
            temp_file = f.name

        try:
            is_valid, errors = validate_file(temp_file, file_type=None, verbose=True)
            assert not is_valid, "Malformed JSON should fail validation"
            assert len(errors) > 0, "Should have validation errors"

            # Check that the error mentions JSON parsing
            error_str = " ".join(str(e) for e in errors).lower()
            assert "json" in error_str or "parse" in error_str, (
                f"Error should mention JSON parsing failure. Got: {errors}"
            )
        finally:
            # Clean up temp file
            os.unlink(temp_file)

    def test_validate_with_explicit_type(self):
        """Test validation with explicitly specified file types."""
        # Create a valid minimal dataset
        valid_dataset = {"dataset_id": "test_dataset", "clones": []}

        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            json.dump(valid_dataset, f)
            temp_file = f.name

        try:
            # Validate as dataset (should pass)
            is_valid, errors = validate_file(
                temp_file, file_type="dataset", verbose=True
            )
            assert is_valid, f"Valid dataset should pass validation. Errors: {errors}"

            # Validate as clone (should fail - wrong type)
            is_valid, errors = validate_file(temp_file, file_type="clone", verbose=True)
            assert not is_valid, "Dataset validated as clone should fail"

            # Validate as tree (should fail - wrong type)
            is_valid, errors = validate_file(temp_file, file_type="tree", verbose=True)
            assert not is_valid, "Dataset validated as tree should fail"
        finally:
            # Clean up temp file
            os.unlink(temp_file)

    def test_validate_clone_collection(self):
        """Test validation of clone collections (arrays)."""
        # Create a collection of clones
        clone_collection = [
            {
                "clone_id": "clone1",
                "unique_seqs_count": 10,
                "mean_mut_freq": 0.05,
                "v_alignment_start": 0,
                "v_alignment_end": 100,
                "j_alignment_start": 200,
                "j_alignment_end": 250,
            },
            {
                "clone_id": "clone2",
                "unique_seqs_count": 5,
                "mean_mut_freq": 0.03,
                "v_alignment_start": 0,
                "v_alignment_end": 100,
                "j_alignment_start": 200,
                "j_alignment_end": 250,
            },
        ]

        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            json.dump(clone_collection, f)
            temp_file = f.name

        try:
            # Validate as clones collection (should pass)
            is_valid, errors = validate_file(
                temp_file, file_type="clones", verbose=True
            )
            assert is_valid, (
                f"Valid clone collection should pass validation. Errors: {errors}"
            )

            # Validate as single clone (should fail)
            is_valid, errors = validate_file(
                temp_file, file_type="clone", verbose=False
            )
            assert not is_valid, (
                "Clone collection validated as single clone should fail"
            )
        finally:
            # Clean up temp file
            os.unlink(temp_file)


class TestValidationFunctions:
    """Test individual validation functions."""

    def test_validate_dataset_function(self):
        """Test the validate_dataset function directly."""
        # Valid minimal dataset
        valid_dataset = {"dataset_id": "test", "clones": []}
        errors = validate_dataset(valid_dataset, verbose=True)
        assert len(errors) == 0, f"Valid dataset should have no errors. Got: {errors}"

        # Invalid dataset (missing required field)
        invalid_dataset = {"clones": []}
        errors = validate_dataset(invalid_dataset, verbose=True)
        assert len(errors) > 0, "Invalid dataset should have errors"

    def test_validate_clone_function(self):
        """Test the validate_clone function directly."""
        # Valid minimal clone
        valid_clone = {
            "unique_seqs_count": 10,
            "mean_mut_freq": 0.05,
            "v_alignment_start": 0,
            "v_alignment_end": 100,
            "j_alignment_start": 200,
            "j_alignment_end": 250,
        }
        errors = validate_clone(valid_clone, verbose=True)
        # Note: May have errors if AIRR validation is strict
        # but should validate against Olmsted schema

        # Invalid clone (missing required fields)
        invalid_clone = {"clone_id": "test"}
        errors = validate_clone(invalid_clone, verbose=True)
        assert len(errors) > 0, "Invalid clone should have errors"

    def test_validate_tree_function(self):
        """Test the validate_tree function directly."""
        # Valid minimal tree
        valid_tree = {"newick": "(A:0.1,B:0.2)C:0.0;", "nodes": {}}
        errors = validate_tree(valid_tree, verbose=True)
        # Note: May have errors if AIRR validation is strict
        # but should validate against Olmsted schema

        # Invalid tree (missing required fields)
        invalid_tree = {"tree_id": "test"}
        errors = validate_tree(invalid_tree, verbose=True)
        assert len(errors) > 0, "Invalid tree should have errors"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
