"""Tests for the build-config command."""

import json
import os
import subprocess
import tempfile

import pytest

from olmsted_cli.build_config import _check_mutation_demotion, _looks_like_local_path
from olmsted_cli.process_utils import unpack_encoded_mutations


class TestBuildConfigOlmsted:
    def test_olmsted_json(self):
        """build-config on Olmsted JSON produces valid YAML with field entries."""
        result = subprocess.run(
            ["olmsted", "build-config", "-i",
             "example_data/surprise/surprise_subset.json"],
            capture_output=True, text=True,
        )
        assert result.returncode == 0
        output = result.stdout
        assert "custom_fields:" in output
        assert "clone" in output
        assert "mutation" in output
        assert "surprise_mutsel" in output

    def test_output_to_file(self):
        """build-config -o writes to file instead of stdout."""
        with tempfile.NamedTemporaryFile(suffix=".yaml", delete=False) as f:
            out_path = f.name

        try:
            result = subprocess.run(
                ["olmsted", "build-config", "-i",
                 "example_data/surprise/surprise_subset.json",
                 "-o", out_path],
                capture_output=True, text=True,
            )
            assert result.returncode == 0
            assert os.path.exists(out_path)
            with open(out_path) as f:
                content = f.read()
            assert "custom_fields:" in content
        finally:
            if os.path.exists(out_path):
                os.unlink(out_path)

    def test_shows_ranges_for_continuous_mutation(self):
        """Continuous mutation fields show range comments."""
        result = subprocess.run(
            ["olmsted", "build-config", "-i",
             "example_data/surprise/surprise_subset.json"],
            capture_output=True, text=True,
        )
        assert "range in data:" in result.stdout

    def test_includes_processing_options_template(self):
        """Output includes commented-out processing options."""
        result = subprocess.run(
            ["olmsted", "build-config", "-i",
             "example_data/surprise/surprise_subset.json"],
            capture_output=True, text=True,
        )
        assert "Processing Options" in result.stdout

    def test_includes_alias_reference(self):
        """Output includes cross-format alias reference."""
        result = subprocess.run(
            ["olmsted", "build-config", "-i",
             "example_data/surprise/surprise_subset.json"],
            capture_output=True, text=True,
        )
        assert "output_name" in result.stdout
        assert "v_call" in result.stdout

    def test_includes_output_name_docs(self):
        """Output includes output_name documentation and alias reference."""
        result = subprocess.run(
            ["olmsted", "build-config", "-i",
             "example_data/surprise/surprise_subset.json"],
            capture_output=True, text=True,
        )
        assert "output_name" in result.stdout
        assert "rearrangement_count" in result.stdout  # in alias reference
        assert "unique_seqs_count" in result.stdout  # in alias reference

    def test_alias_suggestions_on_airr(self):
        """AIRR fields with aliases get output_name suggestions."""
        result = subprocess.run(
            ["olmsted", "build-config", "-i",
             "example_data/airr/airr.json"],
            capture_output=True, text=True,
        )
        # rearrangement_count should have output_name: unique_seqs_count
        lines = result.stdout.split("\n")
        for i, line in enumerate(lines):
            if "name: rearrangement_count" in line:
                # Next line should be output_name
                assert "output_name: unique_seqs_count" in lines[i + 1]
                break
        else:
            pytest.fail("rearrangement_count not found in build-config output")


class TestBuildConfigPcp:
    def test_pcp_with_trees(self):
        """build-config on raw PCP CSV + tree CSV discovers fields."""
        result = subprocess.run(
            ["olmsted", "build-config", "-i",
             "example_data/pcp/pcp.csv",
             "-t", "example_data/pcp/trees.csv"],
            capture_output=True, text=True,
        )
        assert result.returncode == 0
        output = result.stdout
        assert "custom_fields:" in output
        assert "Family level" in output
        assert "unique_seqs_count" in output

    def test_pcp_extra_columns(self):
        """build-config on PCP with extra columns discovers them."""
        result = subprocess.run(
            ["olmsted", "build-config", "-i",
             "example_data/test-fields/pcp-test-fields.csv",
             "-t", "example_data/test-fields/trees-test-fields.csv"],
            capture_output=True, text=True,
        )
        assert result.returncode == 0
        output = result.stdout
        assert "foobar_score" in output
        assert "foobar_category" in output
        assert "foobar_weight" in output
        assert "foobar_class" in output
        # New types: list/json demoted to mutation, json at clone, path auto-skipped
        assert "foobar_per_site_score" in output
        assert "foobar_sparse_aa" in output
        assert "foobar_params" in output
        assert "foobar_path" in output

    def test_pcp_shows_compute_metrics(self):
        """PCP config template includes compute_metrics option."""
        result = subprocess.run(
            ["olmsted", "build-config", "-i",
             "example_data/pcp/pcp.csv",
             "-t", "example_data/pcp/trees.csv"],
            capture_output=True, text=True,
        )
        assert "compute_metrics" in result.stdout


class TestBuildConfigAirr:
    def test_airr(self):
        """build-config on AIRR JSON discovers fields."""
        result = subprocess.run(
            ["olmsted", "build-config", "-i",
             "example_data/airr/airr.json"],
            capture_output=True, text=True,
        )
        assert result.returncode == 0
        output = result.stdout
        assert "custom_fields:" in output
        assert "v_call" in output or "V Gene" in output


class TestBuildConfigNewTypes:
    """Tests for list, json, path detection in build-config output."""

    @pytest.fixture(params=[
        ("example_data/test-fields/olmsted-test-fields.json", None),
        ("example_data/test-fields/airr-test-fields.json", None),
        ("example_data/test-fields/pcp-test-fields.csv",
         "example_data/test-fields/trees-test-fields.csv"),
    ], ids=["olmsted", "airr", "pcp"])
    def build_config_output(self, request):
        input_path, tree_path = request.param
        cmd = ["olmsted", "build-config", "-i", input_path]
        if tree_path:
            cmd += ["-t", tree_path]
        result = subprocess.run(cmd, capture_output=True, text=True)
        assert result.returncode == 0
        return result.stdout

    def test_json_at_clone_level(self, build_config_output):
        """foobar_params (non-int keys) stays at family level as json."""
        assert "foobar_params" in build_config_output
        # Should be in clone/family section with type: json
        lines = build_config_output.split("\n")
        for i, line in enumerate(lines):
            if "name: foobar_params" in line:
                block = "\n".join(lines[i:i+5])
                assert "level: family" in block
                assert "type: json" in block
                break
        else:
            pytest.fail("foobar_params not found in output")

    def test_list_demoted_to_mutation(self, build_config_output):
        """foobar_per_site_score (list matching seq length) demoted to mutation."""
        lines = build_config_output.split("\n")
        for i, line in enumerate(lines):
            if "name: foobar_per_site_score" in line:
                block = "\n".join(lines[i:i+6])
                assert "level: mutation" in block
                assert "encoding: list" in block
                assert "type: continuous" in block
                break
        else:
            pytest.fail("foobar_per_site_score not found in output")

    def test_json_demoted_to_mutation(self, build_config_output):
        """foobar_sparse_aa (int keys, AA values) demoted to mutation."""
        lines = build_config_output.split("\n")
        for i, line in enumerate(lines):
            if "name: foobar_sparse_aa" in line:
                block = "\n".join(lines[i:i+6])
                assert "level: mutation" in block
                assert "encoding: json" in block
                assert "type: aa" in block
                break
        else:
            pytest.fail("foobar_sparse_aa not found in output")

    def test_path_auto_skipped(self, build_config_output):
        """foobar_path (local file paths) auto-skipped."""
        lines = build_config_output.split("\n")
        for i, line in enumerate(lines):
            if "name: foobar_path" in line:
                block = "\n".join(lines[i:i+5])
                assert "skip: true" in block
                break
        else:
            pytest.fail("foobar_path not found in output")

    def test_demotion_comment_present(self, build_config_output):
        """Demoted fields section has explanatory comment."""
        assert "demoted from node to mutation level" in build_config_output


class TestLooksLikeLocalPath:
    """Tests for _looks_like_local_path heuristic."""

    def test_absolute_paths(self):
        assert _looks_like_local_path(["/data/raw/file.fasta", "/tmp/out.json"])

    def test_relative_paths(self):
        assert _looks_like_local_path(["./data/file.csv", "../other/file.txt"])

    def test_home_paths(self):
        assert _looks_like_local_path(["~/Documents/data.json"])

    def test_windows_paths(self):
        assert _looks_like_local_path(["C:\\Users\\data\\file.csv"])

    def test_urls_not_paths(self):
        assert not _looks_like_local_path(["https://example.com/data.json"])
        assert not _looks_like_local_path(["http://api.example.com/v1"])

    def test_regular_strings_not_paths(self):
        assert not _looks_like_local_path(["group-alpha", "group-beta"])
        assert not _looks_like_local_path(["IGHV3-48*01", "IGHJ4*02"])

    def test_empty_list(self):
        assert not _looks_like_local_path([])

    def test_non_string_values_ignored(self):
        assert not _looks_like_local_path([42, 3.14, True])

    def test_mixed_paths_and_strings(self):
        """Threshold: at least half must be paths."""
        assert _looks_like_local_path(["/data/a.txt", "/data/b.txt", "other"])
        assert not _looks_like_local_path(["/data/a.txt", "foo", "bar", "baz"])


class TestBuildConfigFormatDetection:
    def test_detects_olmsted(self):
        result = subprocess.run(
            ["olmsted", "build-config", "-i",
             "example_data/surprise/surprise_subset.json"],
            capture_output=True, text=True,
        )
        assert "OLMSTED" in result.stderr

    def test_detects_pcp(self):
        result = subprocess.run(
            ["olmsted", "build-config", "-i",
             "example_data/pcp/pcp.csv"],
            capture_output=True, text=True,
        )
        assert "PCP" in result.stderr

    def test_detects_airr(self):
        result = subprocess.run(
            ["olmsted", "build-config", "-i",
             "example_data/airr/airr.json"],
            capture_output=True, text=True,
        )
        assert "AIRR" in result.stderr


class TestMutationDemotion:
    """Tests for _check_mutation_demotion: detecting node fields with per-position data."""

    def _make_nodes(self, field, values, seq_len=5):
        """Helper: create node dicts with a sequence and a test field."""
        seq = "A" * seq_len
        return [
            {"sequence_alignment_aa": seq, field: v}
            for v in values
        ]

    def test_list_continuous_matching_length(self):
        nodes = self._make_nodes("scores", [[0.1, 0.2, 0.3, 0.4, 0.5]] * 3)
        result = _check_mutation_demotion(nodes, "scores")
        assert result is not None
        assert result["encoding"] == "list"
        assert result["inner_type"] == "continuous"

    def test_list_aa_matching_length(self):
        nodes = self._make_nodes("residues", [["A", "V", "L", "M", "D"]] * 3)
        result = _check_mutation_demotion(nodes, "residues")
        assert result is not None
        assert result["encoding"] == "list"
        assert result["inner_type"] == "aa"

    def test_list_wrong_length(self):
        """Lists that don't match sequence length are not demoted."""
        nodes = self._make_nodes("scores", [[0.1, 0.2, 0.3]] * 3, seq_len=5)
        result = _check_mutation_demotion(nodes, "scores")
        assert result is None

    def test_list_categorical_not_demoted(self):
        """Lists of categorical strings are not demoted."""
        nodes = self._make_nodes("labels", [["foo", "bar", "baz", "qux", "xyz"]] * 3)
        result = _check_mutation_demotion(nodes, "labels")
        assert result is None

    def test_json_continuous_within_range(self):
        nodes = self._make_nodes("sparse_scores", [{"0": 0.5, "3": 0.8}] * 3)
        result = _check_mutation_demotion(nodes, "sparse_scores")
        assert result is not None
        assert result["encoding"] == "json"
        assert result["inner_type"] == "continuous"

    def test_json_aa_within_range(self):
        nodes = self._make_nodes("sparse_aa", [{"0": "D", "3": "E"}] * 3)
        result = _check_mutation_demotion(nodes, "sparse_aa")
        assert result is not None
        assert result["encoding"] == "json"
        assert result["inner_type"] == "aa"

    def test_records_style_detected(self):
        """Array of dicts with 'site' key detected as records encoding."""
        nodes = self._make_nodes("custom_scores", [
            [{"site": 0, "score": 2.5, "region": "FWR1"},
             {"site": 3, "score": 4.1, "region": "CDR1"}],
        ] * 3)
        result = _check_mutation_demotion(nodes, "custom_scores")
        assert result is not None
        assert result["encoding"] == "records"
        assert result["source"] == "custom_scores"
        assert "score" in result["inner_fields"]
        assert "region" in result["inner_fields"]
        assert result["inner_fields"]["score"] == "continuous"
        assert result["inner_fields"]["region"] == "categorical"

    def test_json_keys_out_of_range(self):
        """JSON with keys beyond sequence length are not demoted."""
        nodes = self._make_nodes("scores", [{"0": 0.5, "99": 0.8}] * 3, seq_len=5)
        result = _check_mutation_demotion(nodes, "scores")
        assert result is None

    def test_json_non_int_keys(self):
        """JSON with non-integer keys are not demoted."""
        nodes = self._make_nodes("data", [{"region": 0.5, "site": 0.8}] * 3)
        result = _check_mutation_demotion(nodes, "data")
        assert result is None

    def test_no_sequence_no_demotion(self):
        """Without sequences, no demotion is possible."""
        nodes = [{"scores": [0.1, 0.2, 0.3]}] * 3
        result = _check_mutation_demotion(nodes, "scores")
        assert result is None

    def test_scalar_field_not_demoted(self):
        """Regular scalar fields are not affected."""
        nodes = self._make_nodes("lbi", [0.5, 0.3, 0.8])
        result = _check_mutation_demotion(nodes, "lbi")
        assert result is None


class TestUnpackEncodedMutations:
    """Tests for unpack_encoded_mutations: merging encoded data into mutations arrays."""

    def _make_tree(self, nodes):
        return {"nodes": nodes}

    def test_list_encoding(self):
        """List data unpacked into mutations by index."""
        tree = self._make_tree([
            {"sequence_id": "n1", "scores": [0.1, None, 0.3]},
        ])
        custom_fields = [
            {"name": "scores", "level": "mutation", "encoding": "list",
             "type": "continuous", "label": "Scores"},
        ]
        unpack_encoded_mutations([tree], custom_fields)
        muts = tree["nodes"][0]["mutations"]
        assert len(muts) == 2  # null skipped
        assert muts[0] == {"site": 0, "scores": 0.1}
        assert muts[1] == {"site": 2, "scores": 0.3}

    def test_json_encoding(self):
        """JSON dict data unpacked into mutations by key."""
        tree = self._make_tree([
            {"sequence_id": "n1", "sparse": {"0": "D", "3": "E"}},
        ])
        custom_fields = [
            {"name": "sparse", "level": "mutation", "encoding": "json",
             "type": "aa", "label": "Sparse"},
        ]
        unpack_encoded_mutations([tree], custom_fields)
        muts = tree["nodes"][0]["mutations"]
        assert len(muts) == 2
        assert muts[0] == {"site": 0, "sparse": "D"}
        assert muts[1] == {"site": 3, "sparse": "E"}

    def test_records_encoding(self):
        """Records-style array unpacked by extracting named inner fields."""
        tree = self._make_tree([
            {"sequence_id": "n1", "custom_scores": [
                {"site": 0, "score_a": 2.5, "score_b": 1.0, "region": "FWR1"},
                {"site": 3, "score_a": 4.1, "score_b": 0.8, "region": "CDR1"},
            ]},
        ])
        custom_fields = [
            {"name": "score_a", "level": "mutation", "encoding": "records",
             "source": "custom_scores", "type": "continuous", "label": "Score A"},
            {"name": "score_b", "level": "mutation", "encoding": "records",
             "source": "custom_scores", "type": "continuous", "label": "Score B"},
            {"name": "region", "level": "mutation", "encoding": "records",
             "source": "custom_scores", "type": "categorical", "label": "Region"},
        ]
        unpack_encoded_mutations([tree], custom_fields)
        muts = tree["nodes"][0]["mutations"]
        assert len(muts) == 2
        assert muts[0] == {"site": 0, "score_a": 2.5, "score_b": 1.0, "region": "FWR1"}
        assert muts[1] == {"site": 3, "score_a": 4.1, "score_b": 0.8, "region": "CDR1"}

    def test_merge_with_existing_mutations(self):
        """Encoded data merges with pre-existing mutations array."""
        tree = self._make_tree([
            {"sequence_id": "n1",
             "mutations": [{"site": 0, "existing_field": 99}],
             "scores": [0.5, 0.7]},
        ])
        custom_fields = [
            {"name": "scores", "level": "mutation", "encoding": "list",
             "type": "continuous", "label": "Scores"},
        ]
        unpack_encoded_mutations([tree], custom_fields)
        muts = tree["nodes"][0]["mutations"]
        assert len(muts) == 2
        # Site 0 has both existing and new data
        assert muts[0] == {"site": 0, "existing_field": 99, "scores": 0.5}
        assert muts[1] == {"site": 1, "scores": 0.7}

    def test_multiple_encodings_merge(self):
        """Multiple encoded fields all merge into the same mutations array."""
        tree = self._make_tree([
            {"sequence_id": "n1",
             "per_site": [0.1, 0.2, 0.3],
             "sparse_aa": {"1": "D"}},
        ])
        custom_fields = [
            {"name": "per_site", "level": "mutation", "encoding": "list",
             "type": "continuous", "label": "Per Site"},
            {"name": "sparse_aa", "level": "mutation", "encoding": "json",
             "type": "aa", "label": "Sparse AA"},
        ]
        unpack_encoded_mutations([tree], custom_fields)
        muts = tree["nodes"][0]["mutations"]
        assert len(muts) == 3
        assert muts[1] == {"site": 1, "per_site": 0.2, "sparse_aa": "D"}

    def test_no_encoded_fields_noop(self):
        """No encoding fields → no changes."""
        tree = self._make_tree([
            {"sequence_id": "n1", "mutations": [{"site": 0, "score": 1.0}]},
        ])
        custom_fields = [
            {"name": "score", "level": "mutation", "type": "continuous", "label": "Score"},
        ]
        unpack_encoded_mutations([tree], custom_fields)
        assert tree["nodes"][0]["mutations"] == [{"site": 0, "score": 1.0}]

    def test_none_custom_fields(self):
        """None custom_fields → no crash."""
        tree = self._make_tree([{"sequence_id": "n1"}])
        unpack_encoded_mutations([tree], None)

    def test_dict_nodes_format(self):
        """Works with nodes as dict (AIRR format) not just list."""
        tree = {"nodes": {
            "n1": {"sequence_id": "n1", "scores": [0.5, 0.8]},
        }}
        custom_fields = [
            {"name": "scores", "level": "mutation", "encoding": "list",
             "type": "continuous", "label": "Scores"},
        ]
        unpack_encoded_mutations([tree], custom_fields)
        muts = tree["nodes"]["n1"]["mutations"]
        assert len(muts) == 2
