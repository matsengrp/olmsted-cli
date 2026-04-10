"""Tests for PCP extra column passthrough, chain partitioning, and column aliases."""

import csv
import json
import os
import subprocess
import tempfile

import pytest

from olmsted_cli.process_pcp_data import (
    _normalize_column_names,
    _partition_chain_fields,
    parse_pcp_csv,
    parse_newick_csv,
)
from olmsted_cli.process_utils import coerce_csv_value as _coerce_csv_value


# =============================================================================
# _coerce_csv_value
# =============================================================================


class TestCoerceCsvValue:
    def test_int(self):
        assert _coerce_csv_value("42") == 42

    def test_float(self):
        assert _coerce_csv_value("3.14") == 3.14

    def test_bool_true(self):
        assert _coerce_csv_value("true") is True

    def test_bool_false(self):
        assert _coerce_csv_value("False") is False

    def test_string(self):
        assert _coerce_csv_value("hello") == "hello"

    def test_json_list(self):
        assert _coerce_csv_value('[1, 2, 3]') == [1, 2, 3]

    def test_json_dict(self):
        assert _coerce_csv_value('{"a": 1}') == {"a": 1}

    def test_empty_string(self):
        assert _coerce_csv_value("") == ""


# =============================================================================
# _partition_chain_fields
# =============================================================================


class TestPartitionChainFields:
    def test_shared_fields(self):
        shared, heavy, light = _partition_chain_fields({"score": 1, "label": "a"})
        assert shared == {"score": 1, "label": "a"}
        assert heavy == {}
        assert light == {}

    def test_heavy_suffix_stripped(self):
        shared, heavy, light = _partition_chain_fields({"score_heavy": 0.8})
        assert shared == {}
        assert heavy == {"score": 0.8}
        assert light == {}

    def test_light_suffix_stripped(self):
        shared, heavy, light = _partition_chain_fields({"score_light": 0.3})
        assert shared == {}
        assert heavy == {}
        assert light == {"score": 0.3}

    def test_mixed(self):
        fields = {
            "shared_metric": 42,
            "diversity_heavy": 0.7,
            "diversity_light": 0.3,
        }
        shared, heavy, light = _partition_chain_fields(fields)
        assert shared == {"shared_metric": 42}
        assert heavy == {"diversity": 0.7}
        assert light == {"diversity": 0.3}

    def test_empty(self):
        shared, heavy, light = _partition_chain_fields({})
        assert shared == heavy == light == {}


# =============================================================================
# _normalize_column_names
# =============================================================================


class TestNormalizeColumnNames:
    def test_standard_columns_unchanged(self):
        cols = ["sample_id", "parent_name", "child_name", "v_gene_heavy"]
        col_map, notes = _normalize_column_names(cols)
        assert col_map["v_gene_heavy"] == "v_gene_heavy"
        assert notes == []

    def test_v_gene_alias(self):
        cols = ["sample_id", "parent_name", "child_name", "v_gene"]
        col_map, notes = _normalize_column_names(cols)
        assert col_map["v_gene"] == "v_gene_heavy"
        assert len(notes) == 1
        assert "v_gene" in notes[0]

    def test_alias_not_applied_if_canonical_present(self):
        cols = ["sample_id", "parent_name", "child_name", "v_gene", "v_gene_heavy"]
        col_map, notes = _normalize_column_names(cols)
        # v_gene should NOT be remapped since v_gene_heavy already exists
        assert col_map["v_gene"] == "v_gene"
        assert col_map["v_gene_heavy"] == "v_gene_heavy"

    def test_parent_seq_alias(self):
        cols = ["sample_id", "parent_name", "child_name", "parent_seq", "child_seq"]
        col_map, notes = _normalize_column_names(cols)
        assert col_map["parent_seq"] == "parent_heavy"
        assert col_map["child_seq"] == "child_heavy"

    def test_v_call_alias(self):
        cols = ["sample_id", "parent_name", "child_name", "v_call", "j_call"]
        col_map, notes = _normalize_column_names(cols)
        assert col_map["v_call"] == "v_gene_heavy"
        assert col_map["j_call"] == "j_gene_heavy"

    def test_empty_columns_ignored(self):
        cols = ["", "sample_id", None, "parent_name", "child_name"]
        col_map, notes = _normalize_column_names(cols)
        assert "" not in col_map
        assert None not in col_map


# =============================================================================
# PCP extra column passthrough (integration)
# =============================================================================


def _write_csv(path, rows):
    """Write rows (list of dicts) to a CSV file."""
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


class TestPcpExtraColumns:
    def test_extra_pcp_columns_on_nodes(self, tmp_path):
        """Extra columns in PCP CSV appear on output nodes."""
        pcp_file = tmp_path / "pcp.csv"
        _write_csv(str(pcp_file), [
            {
                "sample_id": "s1", "family": "f1",
                "parent_name": "naive", "child_name": "leaf1",
                "parent_heavy": "ATCG", "child_heavy": "ATGG",
                "branch_length": "0.01", "parent_is_naive": "True",
                "child_is_leaf": "True", "my_score": "3.14",
                "my_label": "interesting",
            },
        ])
        families = parse_pcp_csv(str(pcp_file))
        fam = families["f1"]
        leaf = fam["nodes"]["leaf1"]
        assert leaf["my_score"] == 3.14
        assert leaf["my_label"] == "interesting"

    def test_extra_tree_columns_in_tree_data(self, tmp_path):
        """Extra columns in tree CSV appear in tree_data dict."""
        tree_file = tmp_path / "trees.csv"
        _write_csv(str(tree_file), [
            {
                "family_name": "f1", "sample_id": "s1",
                "newick_tree": "(leaf1:0.01)naive:0.0;",
                "phylo_diversity": "0.85", "clade": "group-A",
            },
        ])
        trees = parse_newick_csv(str(tree_file))
        tree_data = trees[("f1", "s1")]
        assert isinstance(tree_data, dict)
        assert tree_data["phylo_diversity"] == 0.85
        assert tree_data["clade"] == "group-A"

    def test_unnamed_index_column_filtered(self, tmp_path):
        """CSV with leading comma (unnamed index) doesn't pollute nodes."""
        pcp_file = tmp_path / "pcp.csv"
        with open(pcp_file, "w") as f:
            f.write(",sample_id,family,parent_name,child_name,parent_heavy,child_heavy,branch_length,parent_is_naive,child_is_leaf\n")
            f.write("0,s1,f1,naive,leaf1,ATCG,ATGG,0.01,True,True\n")

        families = parse_pcp_csv(str(pcp_file))
        leaf = families["f1"]["nodes"]["leaf1"]
        assert "" not in leaf
        assert 0 not in leaf


class TestPcpChainPartitioning:
    def test_end_to_end_paired_partitioning(self):
        """Extra fields are partitioned by chain suffix in paired processing."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as pf:
            pf.write("sample_id,family,parent_name,child_name,parent_heavy,child_heavy,parent_light,child_light,branch_length,parent_is_naive,child_is_leaf,light_chain_type,score_heavy,score_light,shared_val\n")
            pf.write("s1,f1,naive,leaf1,ATCG,ATGG,GCTA,GCTG,0.01,True,True,kappa,0.8,0.3,42\n")
            pcp_path = pf.name

        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as tf:
            tf.write("family_name,sample_id,newick_tree,div_heavy,div_light,common\n")
            tf.write('f1,s1,"(leaf1:0.01)naive:0.0;",0.7,0.3,99\n')
            tree_path = tf.name

        with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as of:
            out_path = of.name

        try:
            result = subprocess.run(
                ["olmsted", "process", "-i", pcp_path, "-t", tree_path,
                 "-f", "pcp", "-o", out_path, "--seed", "42"],
                capture_output=True, text=True,
            )
            assert result.returncode == 0, f"stderr: {result.stderr}"

            with open(out_path) as f:
                data = json.load(f)

            ds_id = list(data["clones"].keys())[0]
            clones = data["clones"][ds_id]
            heavy = next(c for c in clones if "heavy" in c["clone_id"])
            light = next(c for c in clones if "light" in c["clone_id"])

            # Shared field on both
            assert heavy["common"] == 99
            assert light["common"] == 99

            # Heavy-only (suffix stripped)
            assert heavy["div"] == 0.7
            assert "div" not in light or light.get("div") != 0.7

            # Light-only (suffix stripped)
            assert light["div"] == 0.3

            # Node-level: shared_val on both chains' nodes
            heavy_tree = next(t for t in data["trees"] if "heavy" in t["clone_id"])
            light_tree = next(t for t in data["trees"] if "light" in t["clone_id"])
            heavy_leaf = next(n for n in heavy_tree["nodes"] if n["sequence_id"] == "leaf1")
            light_leaf = next(n for n in light_tree["nodes"] if n["sequence_id"] == "leaf1")
            assert heavy_leaf["shared_val"] == 42
            assert light_leaf["shared_val"] == 42

            # Node-level chain partitioning
            assert heavy_leaf["score"] == 0.8
            assert light_leaf["score"] == 0.3
        finally:
            for p in (pcp_path, tree_path, out_path):
                if os.path.exists(p):
                    os.unlink(p)


# =============================================================================
# Minimal PCP data (no gene calls)
# =============================================================================


class TestPcpMinimalData:
    def test_no_gene_calls(self):
        """PCP data with only required columns + sequences processes without error."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            f.write("sample_id,family,parent_name,child_name,parent_heavy,child_heavy,branch_length,parent_is_naive,child_is_leaf\n")
            f.write("s1,f1,naive,leaf1,ATCGATCG,ATGGATCG,0.01,True,True\n")
            pcp_path = f.name

        with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as of:
            out_path = of.name

        try:
            result = subprocess.run(
                ["olmsted", "process", "-i", pcp_path, "-f", "pcp",
                 "-o", out_path, "--seed", "42"],
                capture_output=True, text=True,
            )
            assert result.returncode == 0, f"stderr: {result.stderr}"

            with open(out_path) as f:
                data = json.load(f)

            ds_id = list(data["clones"].keys())[0]
            clone = data["clones"][ds_id][0]
            assert clone["v_call"] == ""
            assert clone["unique_seqs_count"] > 0
            assert isinstance(clone["mean_mut_freq"], float)
        finally:
            for p in (pcp_path, out_path):
                if os.path.exists(p):
                    os.unlink(p)

    def test_column_aliases(self):
        """Chain-agnostic column names (v_gene, parent_seq) are mapped correctly."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            f.write("sample_id,family,parent_name,child_name,parent_seq,child_seq,branch_length,v_gene,j_gene,parent_is_naive,child_is_leaf\n")
            f.write("s1,f1,naive,leaf1,ATCG,ATGG,0.01,IGHV3-48*01,IGHJ4*02,True,True\n")
            pcp_path = f.name

        with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as of:
            out_path = of.name

        try:
            result = subprocess.run(
                ["olmsted", "process", "-i", pcp_path, "-f", "pcp",
                 "-o", out_path, "--seed", "42", "-v", "2"],
                capture_output=True, text=True,
            )
            assert result.returncode == 0, f"stderr: {result.stderr}"
            assert "mapped to" in (result.stderr + result.stdout)  # Should print alias notifications

            with open(out_path) as f:
                data = json.load(f)

            ds_id = list(data["clones"].keys())[0]
            clone = data["clones"][ds_id][0]
            assert clone["v_call"] == "IGHV3-48*01"
            assert clone["j_call"] == "IGHJ4*02"
        finally:
            for p in (pcp_path, out_path):
                if os.path.exists(p):
                    os.unlink(p)
