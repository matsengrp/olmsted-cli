"""Tests for the IdentMinter utility and deterministic_uuid helper."""

import re
import uuid

import pytest

from olmsted_cli.identifier import IdentMinter, deterministic_uuid


_UUID_RE = re.compile(
    r"^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$"
)


class TestIdentMinter:
    def test_mint_shape_random(self):
        """mint(datatype) returns {datatype}-{uuid}."""
        m = IdentMinter()
        result = m.mint("clone")
        assert result.startswith("clone-")
        tail = result[len("clone-"):]
        # Tail must be a valid UUID string
        assert _UUID_RE.match(tail), f"tail {tail!r} is not a UUID"

    def test_mint_deterministic_under_seed(self):
        """Two minters with the same seed produce the same sequence."""
        m1 = IdentMinter(seed=42)
        m2 = IdentMinter(seed=42)
        seq1 = [m1.mint("clone") for _ in range(3)]
        seq2 = [m2.mint("clone") for _ in range(3)]
        assert seq1 == seq2

    def test_mint_different_seeds_diverge(self):
        """Different seeds produce different sequences."""
        m1 = IdentMinter(seed=1)
        m2 = IdentMinter(seed=2)
        assert m1.mint("clone") != m2.mint("clone")

    def test_mint_random_unique(self):
        """Random-mode minting produces unique values (uuid4 collision
        is astronomically improbable)."""
        m = IdentMinter()
        values = {m.mint("clone") for _ in range(100)}
        assert len(values) == 100

    def test_mint_rejects_empty_datatype(self):
        m = IdentMinter()
        with pytest.raises(ValueError, match="datatype"):
            m.mint("")

    def test_mint_rejects_unknown_datatype(self):
        """Closed-set datatype: even well-formed strings not in the
        registered set are rejected (catches typos like 'Dataset')."""
        m = IdentMinter()
        with pytest.raises(ValueError, match="datatype"):
            m.mint("Dataset")  # capital D — common typo
        with pytest.raises(ValueError, match="datatype"):
            m.mint("repertoire")  # plausible-but-unregistered

    def test_mint_rejects_datatype_with_hyphen(self):
        """Hyphens in datatype would make the {datatype}-{uuid} split
        ambiguous; the closed set has no hyphenated entries so hyphens
        fall through the same rejection path."""
        m = IdentMinter()
        with pytest.raises(ValueError, match="datatype"):
            m.mint("my-thing")

    def test_mint_rejects_non_ascii_datatype(self):
        m = IdentMinter()
        with pytest.raises(ValueError, match="datatype"):
            m.mint("clône")

    def test_counter_independent_per_instance(self):
        """Each minter instance has its own counter; creating a new one
        resets the deterministic sequence."""
        m1 = IdentMinter(seed=42)
        first_a = m1.mint("clone")
        m1.mint("clone")

        m2 = IdentMinter(seed=42)
        first_b = m2.mint("clone")

        assert first_a == first_b

    def test_different_datatypes_different_prefixes(self):
        """Different datatype arguments produce different prefixes on
        the same counter position."""
        m = IdentMinter(seed=42)
        dataset = m.mint("dataset")
        clone = m.mint("clone")
        assert dataset.startswith("dataset-")
        assert clone.startswith("clone-")
        # Different counter positions → different UUIDs
        assert dataset[len("dataset-"):] != clone[len("clone-"):]


class TestDeterministicUuid:
    def test_stable_for_same_inputs(self):
        assert deterministic_uuid(42, 1) == deterministic_uuid(42, 1)

    def test_different_counters_diverge(self):
        assert deterministic_uuid(42, 1) != deterministic_uuid(42, 2)

    def test_different_seeds_diverge(self):
        assert deterministic_uuid(1, 1) != deterministic_uuid(2, 1)

    def test_returns_uuid_shape(self):
        result = deterministic_uuid(42, 1)
        assert _UUID_RE.match(result)
