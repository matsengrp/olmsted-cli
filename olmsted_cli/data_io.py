"""Centralized I/O for Olmsted data files.

Single home for everything that reads, writes, decompresses, or
detects-format-of a data file. The rest of the codebase only reaches
in here — adding a new format or changing how an existing one is
handled means editing one place.

Layering inside the module
--------------------------

- ``_maybe_unzip(path, mode)`` — internal. Pure compression handling:
  if the path ends in ``.gz``, open via ``gzip.open``; otherwise plain
  ``open``. Returns a context-managerable text-mode handle.

- ``detect_file_format(path)`` — public. Identify whether a path is
  Olmsted JSON, AIRR JSON, PCP CSV, or unknown, by extension + content
  peek. Operates on (logically) plain content; uses ``_maybe_unzip``
  internally to peek inside ``.gz`` wrappers transparently.

- ``open_file(path, expected_formats=None)`` — public. Orchestrates:
  detect format → validate against ``expected_formats`` if given →
  open with ``.gz`` transparency → return ``(handle, detected_format)``.
  This is the primary read entry point for callers outside data_io.

Higher-level read wrappers (also public)
----------------------------------------

Format-specific recipes layered on ``open_file`` — they encapsulate the
"open + parse + structural-validate" pattern:

- ``read_olmsted_json(path)``     — parsed dict; checks required top-level keys
- ``read_airr_json(path)``        — parsed dict (caller validates AIRR-shape)
- ``read_pcp_csv_rows(path)``     — iterates ``DictReader`` rows from PCP CSV
- ``read_csv_rows(path)``         — generic CSV iteration (no format detection)
- ``read_yaml_config(path)``      — parsed dict from a YAML config file

Write API
---------

- ``write_file(data, path, output_kind="olmsted_json", **opts)`` —
  dispatcher. Today routes to ``write_olmsted_json``; extension point
  for future output kinds.

- ``write_olmsted_json(data, output_path, json_format)`` — direct entry
  for the only currently-supported output kind. Pretty / compact / gzip
  variants. Gzip writes pin ``mtime=0`` and the embedded filename so
  the compression layer is byte-deterministic.
"""

from __future__ import annotations

import csv
import gzip
import io
import json
from pathlib import Path
from typing import Iterator

import yaml

from .constants import (
    FORMAT_AIRR,
    FORMAT_OLMSTED,
    FORMAT_PCP,
    FORMAT_UNKNOWN,
    OLMSTED_REQUIRED_TOP_LEVEL_KEYS,
)
from .utils import vprint


# --- compression handling (internal) ---------------------------------------


def _maybe_unzip(path, mode: str = "rt"):
    """Open ``path`` with transparent ``.gz`` decompression.

    Pure file open — no format detection, no parsing. Internal to
    ``data_io``; callers outside this module use ``open_file`` (or a
    higher-level read wrapper) instead.
    """
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


# --- format detection -------------------------------------------------------


def detect_file_format(file_path) -> str:
    """Identify a file's data format by extension + content peek.

    Returns one of ``FORMAT_AIRR``, ``FORMAT_PCP``, ``FORMAT_OLMSTED``,
    or ``FORMAT_UNKNOWN``. ``.gz``-wrapped files are inspected
    transparently via ``_maybe_unzip``; format detection always operates
    on the (logically) plain content underneath.
    """
    file_path = Path(file_path)

    # CSV files are always PCP (by extension)
    if file_path.suffix.lower() == ".csv":
        return FORMAT_PCP
    if file_path.suffix.lower() == ".gz" and file_path.stem.endswith(".csv"):
        return FORMAT_PCP

    # JSON files need content inspection to distinguish AIRR from Olmsted
    if file_path.suffix.lower() == ".json" or (
        file_path.suffix.lower() == ".gz" and file_path.stem.endswith(".json")
    ):
        try:
            with _maybe_unzip(file_path) as fh:
                data = json.load(fh)

            if isinstance(data, dict):
                # Explicit format tag in metadata
                metadata = data.get("metadata", {})
                if isinstance(metadata, dict) and metadata.get("format") == FORMAT_OLMSTED:
                    return FORMAT_OLMSTED
                # Heuristic fallback: Olmsted JSON has "datasets" and "metadata"
                if "datasets" in data and "metadata" in data:
                    return FORMAT_OLMSTED
                # AIRR JSON has "clones" or other standard AIRR keys
                if "dataset_id" in data or "clones" in data or "ident" in data:
                    return FORMAT_AIRR
            elif isinstance(data, list):
                # Multi-dataset AIRR
                return FORMAT_AIRR
        except (json.JSONDecodeError, OSError, ValueError):
            pass

    # If extension didn't help, peek at content for CSV indicators
    try:
        with _maybe_unzip(file_path) as fh:
            first_lines = []
            for i, line in enumerate(fh):
                first_lines.append(line.strip())
                if i >= 2:
                    break

            if first_lines:
                first_line = first_lines[0].lower()
                pcp_indicators = [
                    "sample_id",
                    "parent_name",
                    "child_name",
                    "family_name",
                    "newick",
                ]
                if any(indicator in first_line for indicator in pcp_indicators):
                    return FORMAT_PCP

    except Exception as e:
        vprint.error(f"Warning: Could not detect format for {file_path}: {e}")

    return FORMAT_UNKNOWN


# --- public open ------------------------------------------------------------


def open_file(path, expected_formats=None):
    """Open ``path`` (transparent ``.gz``), detect format, validate, return a handle.

    Args:
        path: input file path (string or path-like).
        expected_formats: tuple/list of accepted formats (e.g.
            ``(FORMAT_OLMSTED,)`` or ``(FORMAT_AIRR, FORMAT_PCP)``). If
            ``None``, the caller is signalling "I don't care about the
            detected format" — the file is opened transparently and any
            detected format (including ``unknown``) is returned alongside
            the handle.

    Returns:
        Tuple of ``(file_handle, detected_format)``. The caller is
        responsible for closing the handle (or using it as a context
        manager) and for parsing — CSV and JSON consumers need
        different parsers.

    Raises:
        ValueError: if ``expected_formats`` is set and the detected
            format isn't in the list (including the case where the
            format couldn't be inferred).
    """
    detected = detect_file_format(path)
    if expected_formats is not None and detected not in expected_formats:
        raise ValueError(
            f"Expected one of {tuple(expected_formats)} but detected "
            f"{detected!r} for {path}"
        )
    return _maybe_unzip(path), detected


# --- high-level reads -------------------------------------------------------


def read_olmsted_json(path) -> dict:
    """Open + parse + shape-validate an Olmsted JSON file.

    Validates the file is detected as Olmsted format AND has the
    required top-level keys (``datasets``, ``clones``, ``trees``).
    Raises ``ValueError`` on any failure.
    """
    handle, _ = open_file(path, expected_formats=(FORMAT_OLMSTED,))
    try:
        data = json.load(handle)
    except json.JSONDecodeError as e:
        raise ValueError(f"{path}: invalid JSON ({e})") from e
    finally:
        handle.close()

    missing = [k for k in OLMSTED_REQUIRED_TOP_LEVEL_KEYS if k not in data]
    if missing:
        raise ValueError(
            f"{path}: missing required Olmsted top-level keys: {sorted(missing)}"
        )
    return data


def read_airr_json(path) -> dict:
    """Open + parse an AIRR JSON file. Caller validates AIRR-specific shape."""
    handle, _ = open_file(path, expected_formats=(FORMAT_AIRR,))
    try:
        return json.load(handle)
    except json.JSONDecodeError as e:
        raise ValueError(f"{path}: invalid JSON ({e})") from e
    finally:
        handle.close()


def read_pcp_csv_rows(path) -> Iterator[dict]:
    """Iterate dict rows from a PCP CSV (clones-style or trees-style).

    Both ``pcp.csv`` and ``trees.csv`` are detected as PCP format; the
    caller distinguishes by columns. Generator — file closes when the
    generator is exhausted or garbage-collected.
    """
    handle, _ = open_file(path, expected_formats=(FORMAT_PCP,))
    with handle:
        yield from csv.DictReader(handle)


def read_csv_rows(path) -> Iterator[dict]:
    """Generic CSV reader (transparent ``.gz``). No format detection.

    Use for CSVs that aren't auto-detectable as a known Olmsted format —
    currently the mutations CSV consumed by ``merge``.
    """
    with _maybe_unzip(path) as handle:
        yield from csv.DictReader(handle)


def read_yaml_config(path) -> dict:
    """Open + parse a YAML config file (transparent ``.gz``)."""
    with _maybe_unzip(path) as handle:
        return yaml.safe_load(handle)


# --- high-level writes ------------------------------------------------------


def write_file(data, path, output_kind: str = "olmsted_json", **opts) -> str:
    """Dispatch to the writer for ``output_kind``. Returns the path written.

    Today the only supported kind is ``olmsted_json``. Adding a new
    output kind (CSV bundles, AIRR-flavored JSON, split-files archives,
    ...) means adding a branch here and the corresponding writer.
    """
    if output_kind == "olmsted_json":
        return write_olmsted_json(data, path, **opts)
    raise ValueError(f"Unknown output_kind: {output_kind!r}")


def write_olmsted_json(data, output_path, json_format: str = "pretty", default=None) -> str:
    """Write Olmsted JSON to ``output_path`` in the requested format.

    Three formats:

    - ``pretty`` — indent=4, human-readable
    - ``compact`` — no whitespace
    - ``gzip`` — pretty content, gzipped; ``.gz`` is auto-appended to
      ``output_path`` if not already present

    Gzip output pins both the header timestamp (``mtime=0``) and the
    embedded filename, so the compression layer is byte-deterministic.
    Content-level non-determinism (timestamps inside the JSON, some
    field-iteration ordering) still flows through unchanged.

    Returns the actual path written (may differ from input when the
    ``gzip`` format auto-appends ``.gz``).
    """
    output_path = str(output_path)
    if json_format == "gzip" and not output_path.endswith(".gz"):
        output_path = output_path + ".gz"

    indent = 4 if json_format in ("pretty", "gzip") else None
    separators = (",", ":") if json_format == "compact" else None

    if json_format == "gzip" or output_path.endswith(".gz"):
        # gzip.open doesn't accept mtime; go through GzipFile directly so we
        # can pin both header timestamp (mtime=0) and the embedded filename
        # — together they give a byte-deterministic gzip layer.
        with open(output_path, "wb") as raw:
            with gzip.GzipFile(filename="", fileobj=raw, mode="wb", mtime=0) as gz:
                with io.TextIOWrapper(gz, encoding="utf-8") as fh:
                    json.dump(data, fh, default=default, indent=indent, separators=separators)
    else:
        with open(output_path, "w") as fh:
            json.dump(data, fh, default=default, indent=indent, separators=separators)

    return output_path
