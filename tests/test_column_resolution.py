"""Tests for the role-column resolver."""

import pytest

from olmsted_cli.column_resolution import (
    RoleColumnConflict,
    RoleColumnNotFound,
    check_row_role_conflicts,
    find_present_variants,
    resolve_role_columns,
)


class TestResolveRoleColumns:
    def test_bare_only(self):
        result = resolve_role_columns(["sample", "family"])
        assert result == {"sample": "sample", "family": "family", "tree": None}

    def test_id_only(self):
        result = resolve_role_columns(["sample_id", "family_id"])
        assert result == {
            "sample": "sample_id",
            "family": "family_id",
            "tree": None,
        }

    def test_name_only(self):
        result = resolve_role_columns(["sample_name", "family_name"])
        assert result == {
            "sample": "sample_name",
            "family": "family_name",
            "tree": None,
        }

    def test_preference_id_over_bare_over_name(self):
        result = resolve_role_columns(
            [
                "sample", "sample_id", "sample_name",
                "family", "family_id", "family_name",
                "tree", "tree_id", "tree_name",
            ]
        )
        assert result == {
            "sample": "sample_id",
            "family": "family_id",
            "tree": "tree_id",
        }

    def test_preference_bare_over_name(self):
        result = resolve_role_columns(
            ["sample", "sample_name", "family", "family_name"]
        )
        assert result["sample"] == "sample"
        assert result["family"] == "family"

    def test_tree_present(self):
        result = resolve_role_columns(["sample", "family", "tree_name"])
        assert result["tree"] == "tree_name"

    def test_tree_absent_returns_none(self):
        result = resolve_role_columns(["sample", "family"])
        assert result["tree"] is None

    def test_missing_required_family_raises(self):
        with pytest.raises(RoleColumnNotFound, match="family"):
            resolve_role_columns(["sample"])

    def test_missing_required_sample_raises(self):
        with pytest.raises(RoleColumnNotFound, match="sample"):
            resolve_role_columns(["family"])

    def test_override_supersedes_auto_detection(self):
        result = resolve_role_columns(
            ["sample", "sample_id", "family"],
            sample_override="sample",
        )
        assert result["sample"] == "sample"

    def test_override_for_optional_tree(self):
        result = resolve_role_columns(
            ["sample_id", "family_id", "method"],
            tree_override="method",
        )
        assert result["tree"] == "method"

    def test_override_naming_missing_column_raises(self):
        with pytest.raises(RoleColumnNotFound, match="bogus"):
            resolve_role_columns(
                ["sample", "family"],
                sample_override="bogus",
            )


class TestFindPresentVariants:
    def test_lists_all_present_in_preference_order(self):
        result = find_present_variants(
            ["sample", "sample_id", "family", "tree_name"]
        )
        assert result["sample"] == ["sample_id", "sample"]
        assert result["family"] == ["family"]
        assert result["tree"] == ["tree_name"]

    def test_empty_when_absent(self):
        result = find_present_variants(["other_col"])
        assert result == {"sample": [], "family": [], "tree": []}


class TestCheckRowRoleConflicts:
    def test_no_conflict_with_single_variant(self):
        present = {"sample": ["sample_id"], "family": ["family"], "tree": []}
        check_row_role_conflicts(
            {"sample_id": "S1", "family": "F1"}, present
        )

    def test_no_conflict_when_values_match(self):
        present = {
            "sample": ["sample_id", "sample"],
            "family": ["family"],
            "tree": [],
        }
        check_row_role_conflicts(
            {"sample_id": "S1", "sample": "S1", "family": "F1"},
            present,
        )

    def test_conflict_raises_with_role_name(self):
        present = {
            "sample": ["sample_id", "sample"],
            "family": ["family"],
            "tree": [],
        }
        with pytest.raises(RoleColumnConflict, match="sample"):
            check_row_role_conflicts(
                {"sample_id": "S1", "sample": "S2", "family": "F1"},
                present,
            )

    def test_conflict_for_tree_role(self):
        present = {
            "sample": ["sample_id"],
            "family": ["family"],
            "tree": ["tree_id", "tree"],
        }
        with pytest.raises(RoleColumnConflict, match="tree"):
            check_row_role_conflicts(
                {"sample_id": "S1", "family": "F1", "tree_id": "T1", "tree": "T2"},
                present,
            )
