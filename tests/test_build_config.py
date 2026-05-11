"""Tests for the build-config command."""

import difflib
import json
import os
import subprocess
import tempfile
from pathlib import Path

import pytest

from olmsted_cli.build_config import (
    _build_yaml,
    _check_mutation_demotion,
    _looks_like_local_path,
    generate_default_config,
)
from olmsted_cli.process_utils import unpack_encoded_mutations

REPO_ROOT = Path(__file__).parent.parent
GOLDEN_DIR = REPO_ROOT / "tests" / "golden" / "build_config"


class TestBuildConfigOlmsted:
    def test_olmsted_json(self):
        """build-config on Olmsted JSON produces valid YAML with field entries."""
        result = subprocess.run(
            ["olmsted", "build-config", "-i",
             "example-data/mutations/input-olmsted.json"],
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
                 "example-data/mutations/input-olmsted.json",
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
             "example-data/mutations/input-olmsted.json"],
            capture_output=True, text=True,
        )
        assert "range in data:" in result.stdout

    def test_includes_processing_options_template(self):
        """Output includes commented-out processing options."""
        result = subprocess.run(
            ["olmsted", "build-config", "-i",
             "example-data/mutations/input-olmsted.json"],
            capture_output=True, text=True,
        )
        assert "Processing Options" in result.stdout

    def test_includes_alias_reference(self):
        """Output includes cross-format alias reference."""
        result = subprocess.run(
            ["olmsted", "build-config", "-i",
             "example-data/mutations/input-olmsted.json"],
            capture_output=True, text=True,
        )
        assert "output_name" in result.stdout
        assert "v_call" in result.stdout

    def test_includes_output_name_docs(self):
        """Output includes output_name documentation and alias reference."""
        result = subprocess.run(
            ["olmsted", "build-config", "-i",
             "example-data/mutations/input-olmsted.json"],
            capture_output=True, text=True,
        )
        assert "output_name" in result.stdout
        assert "rearrangement_count" in result.stdout  # in alias reference
        assert "unique_seqs_count" in result.stdout  # in alias reference

    def test_alias_suggestions_on_airr(self):
        """AIRR fields with aliases get output_name suggestions."""
        result = subprocess.run(
            ["olmsted", "build-config", "-i",
             "example-data/airr/input-airr.json"],
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
             "example-data/pcp/input-pcp.csv",
             "-t", "example-data/pcp/input-trees.csv"],
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
             "example-data/fields-config/input-pcp.csv",
             "-t", "example-data/fields-config/input-trees.csv"],
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
             "example-data/pcp/input-pcp.csv",
             "-t", "example-data/pcp/input-trees.csv"],
            capture_output=True, text=True,
        )
        assert "compute_metrics" in result.stdout


class TestBuildConfigAirr:
    def test_airr(self):
        """build-config on AIRR JSON discovers fields."""
        result = subprocess.run(
            ["olmsted", "build-config", "-i",
             "example-data/airr/input-airr.json"],
            capture_output=True, text=True,
        )
        assert result.returncode == 0
        output = result.stdout
        assert "custom_fields:" in output
        assert "v_call" in output or "V Gene" in output


class TestBuildConfigNewTypes:
    """Tests for list, json, path detection in build-config output."""

    @pytest.fixture(params=[
        ("example-data/fields-config/input-olmsted.json", None),
        ("example-data/fields-config/input-airr.json", None),
        ("example-data/fields-config/input-pcp.csv",
         "example-data/fields-config/input-trees.csv"),
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


class TestBuildConfigTreeLevel:
    """Tree-level coverage across all three fields-config formats.

    `clone-A` (AIRR / Olmsted) and `fam-1` (PCP) each carry two trees
    with `foobar_method` / `foobar_tree_score` differing across the
    pair, so the variance classifier promotes them to ``tree`` level.
    """

    @pytest.fixture(params=[
        ("example-data/fields-config/input-olmsted.json", None),
        ("example-data/fields-config/input-airr.json", None),
        ("example-data/fields-config/input-pcp.csv",
         "example-data/fields-config/input-trees.csv"),
    ], ids=["olmsted", "airr", "pcp"])
    def build_config_output(self, request):
        input_path, tree_path = request.param
        cmd = ["olmsted", "build-config", "-i", input_path]
        if tree_path:
            cmd += ["-t", tree_path]
        result = subprocess.run(cmd, capture_output=True, text=True)
        assert result.returncode == 0
        return result.stdout

    def test_tree_level_section_present(self, build_config_output):
        assert "Tree level" in build_config_output

    def test_tree_level_fields_emitted(self, build_config_output):
        lines = build_config_output.split("\n")
        for name in ("foobar_method", "foobar_tree_score"):
            for i, line in enumerate(lines):
                if f"name: {name}" in line:
                    block = "\n".join(lines[i:i + 5])
                    assert "level: tree" in block, (
                        f"{name} should be at tree level:\n{block}"
                    )
                    break
            else:
                pytest.fail(f"{name} not found in build-config output")


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
             "example-data/mutations/input-olmsted.json"],
            capture_output=True, text=True,
        )
        assert "OLMSTED" in result.stderr or "OLMSTED" in result.stdout

    def test_detects_pcp(self):
        result = subprocess.run(
            ["olmsted", "build-config", "-i",
             "example-data/pcp/input-pcp.csv"],
            capture_output=True, text=True,
        )
        assert "PCP" in result.stderr or "PCP" in result.stdout

    def test_detects_airr(self):
        result = subprocess.run(
            ["olmsted", "build-config", "-i",
             "example-data/airr/input-airr.json"],
            capture_output=True, text=True,
        )
        assert "AIRR" in result.stderr or "AIRR" in result.stdout


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


class TestDemotedFieldAppearsAsNodeSkip:
    """When build-config demotes a node field to mutation level, it should also
    emit a skip entry at node level so users can opt the field back in."""

    def _make_clones(self):
        return [
            {
                "clone_id": "c1",
                "dataset_id": "ds1",
                "unique_seqs_count": 10,
                "v_call": "IGHV3-48*01",
            },
        ]

    def _make_trees_with_records(self):
        return [
            {
                "ident": "tree-1",
                "clone_id": "c1",
                "nodes": [
                    {
                        "sequence_id": "root",
                        "parent": None,
                        "type": "root",
                        "sequence_alignment_aa": "MKVL",
                        "distance": 0.0,
                        "length": 0.0,
                        "multiplicity": 1,
                    },
                    {
                        "sequence_id": "a",
                        "parent": "root",
                        "type": "leaf",
                        "sequence_alignment_aa": "MEVL",
                        "distance": 0.1,
                        "length": 0.1,
                        "multiplicity": 3,
                        "surprise_mutations": [
                            {"site": 0, "score": 2.5, "region": "FWR1"},
                            {"site": 2, "score": 4.1, "region": "CDR1"},
                        ],
                    },
                ],
            }
        ]

    def _make_trees_with_list(self):
        return [
            {
                "ident": "tree-1",
                "clone_id": "c1",
                "nodes": [
                    {
                        "sequence_id": "a",
                        "parent": "root",
                        "type": "leaf",
                        "sequence_alignment_aa": "MKVL",
                        "distance": 0.1,
                        "length": 0.1,
                        "multiplicity": 3,
                        "per_site_scores": [0.1, 0.2, 0.3, 0.4],
                    },
                ],
            }
        ]

    def test_records_demoted_field_has_node_skip_entry(self):
        """A records-style demoted field appears both as mutation entries
        and as a node-level skip entry."""
        yaml_output = _build_yaml(
            "test.json", "olmsted",
            self._make_clones(), self._make_trees_with_records(),
        )
        lines = yaml_output.split("\n")

        # Should appear as mutation-level entries
        assert any("encoding: records" in line for line in lines)
        assert any("source: surprise_mutations" in line for line in lines)

        # Should also appear as a node-level skip entry
        found_node_skip = False
        for i, line in enumerate(lines):
            if "name: surprise_mutations" in line:
                # Look at nearby lines for level: node and skip: true
                block = "\n".join(lines[max(0, i):i + 6])
                if "level: node" in block and "skip: true" in block:
                    found_node_skip = True
                    break
        assert found_node_skip, (
            "Demoted field 'surprise_mutations' should appear as a "
            "node-level skip entry in build-config output"
        )

    def test_list_demoted_field_has_node_skip_entry(self):
        """A list-encoded demoted field appears as a node-level skip entry."""
        yaml_output = _build_yaml(
            "test.json", "olmsted",
            self._make_clones(), self._make_trees_with_list(),
        )
        lines = yaml_output.split("\n")

        # Should appear as mutation-level entry
        assert any("encoding: list" in line for line in lines)

        # Should also appear as a node-level skip entry
        found_node_skip = False
        for i, line in enumerate(lines):
            if "name: per_site_scores" in line:
                block = "\n".join(lines[max(0, i):i + 6])
                if "level: node" in block and "skip: true" in block:
                    found_node_skip = True
                    break
        assert found_node_skip, (
            "Demoted field 'per_site_scores' should appear as a "
            "node-level skip entry in build-config output"
        )


class TestGenerateDefaultConfig:
    """Tests for generate_default_config: structured config generation."""

    def test_returns_list_of_dicts(self):
        """generate_default_config returns a list of field declaration dicts."""
        clones = [{"clone_id": "c1", "unique_seqs_count": 10, "v_call": "V1"}]
        trees = [{"nodes": [{"sequence_id": "a", "distance": 0.1, "length": 0.1,
                              "multiplicity": 3}]}]
        config = generate_default_config(clones, trees)
        assert isinstance(config, list)
        assert all(isinstance(cf, dict) for cf in config)
        assert all("name" in cf and "level" in cf and "type" in cf for cf in config)

    def test_clone_fields_discovered(self):
        clones = [{"clone_id": "c1", "unique_seqs_count": 10, "v_call": "IGHV3-48*01",
                    "mean_mut_freq": 0.05}]
        config = generate_default_config(clones, [])
        names = {cf["name"] for cf in config if cf["level"] == "clone"}
        assert "unique_seqs_count" in names
        assert "v_call" in names
        assert "mean_mut_freq" in names

    def test_skip_fields_marked(self):
        clones = [{"clone_id": "c1", "ident": "abc", "unique_seqs_count": 10}]
        config = generate_default_config(clones, [])
        ident_entry = next(cf for cf in config if cf["name"] == "ident")
        assert ident_entry.get("skip") is True
        usc_entry = next(cf for cf in config if cf["name"] == "unique_seqs_count")
        assert "skip" not in usc_entry or not usc_entry["skip"]

    def test_no_skip_flag(self):
        clones = [{"clone_id": "c1", "ident": "abc"}]
        config = generate_default_config(clones, [], no_skip=True)
        ident_entry = next(cf for cf in config if cf["name"] == "ident")
        assert "skip" not in ident_entry or not ident_entry["skip"]

    def test_records_demotion(self):
        """Records-style node fields produce mutation entries with encoding."""
        trees = [{
            "nodes": [{
                "sequence_id": "a", "type": "leaf",
                "sequence_alignment_aa": "MKVL",
                "distance": 0.1, "length": 0.1, "multiplicity": 3,
                "surprise_mutations": [
                    {"site": 0, "score": 2.5, "region": "FWR1"},
                    {"site": 2, "score": 4.1, "region": "CDR1"},
                ],
            }],
        }]
        config = generate_default_config([], trees)
        # Node-level skip entry for source field
        node_entries = [cf for cf in config if cf["level"] == "node" and cf["name"] == "surprise_mutations"]
        assert len(node_entries) == 1
        assert node_entries[0].get("skip") is True
        # Mutation-level entries with encoding
        mutation_entries = [cf for cf in config if cf["level"] == "mutation" and cf.get("encoding") == "records"]
        assert len(mutation_entries) > 0
        assert all(cf.get("source") == "surprise_mutations" for cf in mutation_entries)
        mut_names = {cf["name"] for cf in mutation_entries}
        assert "score" in mut_names
        assert "region" in mut_names

    def test_derived_aa_fields(self):
        """When nodes have AA sequences but no mutations, child_aa/parent_aa are added."""
        trees = [{
            "nodes": [{
                "sequence_id": "a", "type": "leaf",
                "sequence_alignment_aa": "MKVL",
                "distance": 0.1, "length": 0.1, "multiplicity": 3,
            }],
        }]
        config = generate_default_config([], trees)
        mut_names = [cf["name"] for cf in config if cf["level"] == "mutation"]
        assert "child_aa" in mut_names
        assert "parent_aa" in mut_names

    def test_branch_fields_whitelisted(self):
        """Only known branch fields appear at branch level."""
        trees = [{
            "nodes": [{
                "sequence_id": "a", "distance": 0.1, "length": 0.5,
                "multiplicity": 3, "custom_field": "foo",
            }],
        }]
        config = generate_default_config([], trees)
        branch_entries = [cf for cf in config if cf["level"] == "branch"]
        branch_names = {cf["name"] for cf in branch_entries}
        assert "length" in branch_names
        assert "custom_field" not in branch_names

    def test_empty_values_excluded(self):
        """Fields with only None values are not included."""
        clones = [{"clone_id": "c1", "empty_field": None, "real_field": 42}]
        config = generate_default_config(clones, [])
        names = {cf["name"] for cf in config}
        assert "empty_field" not in names
        assert "real_field" in names


class TestProcessMatchesBuildConfig:
    """Integration test: process without config should produce the same
    field_metadata as build-config + process -c."""

    def test_tag_without_config_matches(self):
        """Tagging without config produces same result as with default config."""
        # Process to get olmsted JSON
        with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as f:
            olmsted_path = f.name
        with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as f:
            tagged_no_config = f.name
        with tempfile.NamedTemporaryFile(suffix=".yaml", delete=False) as f:
            config_path = f.name
        with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as f:
            tagged_with_config = f.name

        try:
            # Create olmsted JSON from AIRR
            subprocess.run(
                ["olmsted", "process", "-f", "airr",
                 "-i", "example-data/airr/input-airr.json",
                 "-o", olmsted_path, "--seed", "42",
                 "--name", "test", "-q"],
                capture_output=True, text=True, check=True,
            )

            # Generate config from that JSON
            subprocess.run(
                ["olmsted", "build-config", "-i", olmsted_path,
                 "-o", config_path],
                capture_output=True, text=True, check=True,
            )

            # Tag without config
            subprocess.run(
                ["olmsted", "tag", "-i", olmsted_path,
                 "-o", tagged_no_config, "--mode", "overwrite"],
                capture_output=True, text=True, check=True,
            )

            # Tag with config
            subprocess.run(
                ["olmsted", "tag", "-i", olmsted_path,
                 "-o", tagged_with_config, "--mode", "overwrite",
                 "-c", config_path],
                capture_output=True, text=True, check=True,
            )

            # Compare field_metadata
            with open(tagged_no_config) as f:
                data_no_config = json.load(f)
            with open(tagged_with_config) as f:
                data_with_config = json.load(f)

            for ds_nc, ds_wc in zip(data_no_config["datasets"],
                                     data_with_config["datasets"]):
                fm_nc = ds_nc.get("field_metadata", {})
                fm_wc = ds_wc.get("field_metadata", {})
                assert fm_nc == fm_wc, (
                    f"field_metadata mismatch between tag (no config) and "
                    f"tag (with build-config output):\n"
                    f"no-config levels: {sorted(fm_nc.keys())}\n"
                    f"with-config levels: {sorted(fm_wc.keys())}"
                )
        finally:
            for p in (olmsted_path, tagged_no_config, config_path, tagged_with_config):
                if os.path.exists(p):
                    os.unlink(p)


class TestBuildConfigGolden:
    """Snapshot tests: compare stdout against committed golden YAML files.

    Set UPDATE_GOLDEN=1 to regenerate the goldens from current output.
    Run from repo root so the embedded input/tree paths stay stable.
    """

    @pytest.mark.parametrize(
        "scenario,cli_args",
        [
            ("pcp", ["-f", "pcp", "-i", "example-data/pcp/input-pcp.csv",
                     "-t", "example-data/pcp/input-trees.csv"]),
            ("airr", ["-i", "example-data/airr/input-airr.json"]),
            ("olmsted", ["-i", "example-data/mutations/input-olmsted.json"]),
        ],
    )
    def test_matches_golden(self, scenario, cli_args):
        golden_path = GOLDEN_DIR / f"{scenario}.yaml"
        result = subprocess.run(
            ["olmsted", "build-config", "-q", *cli_args],
            capture_output=True, text=True, cwd=REPO_ROOT,
        )
        assert result.returncode == 0, f"build-config failed: {result.stderr}"
        actual = result.stdout

        if os.environ.get("UPDATE_GOLDEN"):
            golden_path.parent.mkdir(parents=True, exist_ok=True)
            golden_path.write_text(actual)
            pytest.skip(f"Regenerated golden: {golden_path}")

        assert golden_path.exists(), (
            f"Missing golden file {golden_path}. "
            f"Run `UPDATE_GOLDEN=1 pytest {__file__}` to create it."
        )
        expected = golden_path.read_text()
        if actual != expected:
            diff = "".join(difflib.unified_diff(
                expected.splitlines(keepends=True),
                actual.splitlines(keepends=True),
                fromfile=str(golden_path),
                tofile=f"build-config --{scenario} (actual)",
            ))
            pytest.fail(
                f"build-config output drifted from golden ({scenario}).\n"
                f"If the change is intentional, run "
                f"`UPDATE_GOLDEN=1 pytest {__file__}::TestBuildConfigGolden` "
                f"to regenerate.\n\n{diff}"
            )
