"""Centralized I/O for Olmsted data files.

Single home for opening, reading, and writing the file formats this CLI
consumes (Olmsted JSON, AIRR JSON, PCP CSV, mutations CSV, YAML config)
and produces (Olmsted JSON in pretty / compact / gzip variants).
Everywhere else in the codebase reads or writes through this module —
adding a new format or changing how an existing one is handled means
editing one place.

Read API
--------

Two layers, mirroring the write side:

- ``open_input(path, expected_formats=None)`` — low-level. Open the file
  (transparent ``.gz``), detect its format, optionally validate against
  an expected set, return ``(handle, detected_format)``. Caller parses.

- High-level helpers, one per data format. Thin wrappers today; the
  value is centralization for future schema-validation, format-version
  handling, or alternative readers:

  - ``read_olmsted_json(path)``        — parsed dict, top-level keys checked
  - ``read_airr_json(path)``           — parsed dict (caller owns AIRR-shape validation)
  - ``read_pcp_csv_rows(path)``        — iterates ``DictReader`` rows from a PCP CSV
  - ``read_csv_rows(path)``            — generic CSV iteration (no format detection)
  - ``read_yaml_config(path)``         — parsed dict from a YAML config file

Format detection runs at the boundary; nothing past ``open_input`` needs
to inspect file extensions or magic bytes.

Write API
---------

- ``write_output(data, path, output_kind="olmsted_json", **opts)`` —
  dispatcher. Today routes to ``write_olmsted_json``; the extension
  point for future output kinds (CSV bundles, AIRR-flavored JSON, ...).

- ``write_olmsted_json(data, output_path, json_format)`` — direct entry
  for the only currently-supported output kind. Produces ``pretty`` /
  ``compact`` / ``gzip`` variants. Gzip writes pin ``mtime=0`` and the
  embedded filename so the compression layer is byte-deterministic.

Lower-level pieces
------------------

- ``open_maybe_gzip(path, mode)`` lives in ``utils.py`` (zero project
  deps). Used here and in format-detection-time validators that can't
  bootstrap through ``open_input`` because they run *during* detection.
"""

from __future__ import annotations

import csv
import gzip
import io
import json
from typing import Iterator

import yaml

from .constants import (
    FORMAT_AIRR,
    FORMAT_OLMSTED,
    FORMAT_PCP,
    FORMAT_UNKNOWN,
    OLMSTED_REQUIRED_TOP_LEVEL_KEYS,
)
from .format_detection import detect_file_format
from .utils import open_maybe_gzip


# --- low-level open ---------------------------------------------------------


def open_input(path, expected_formats=None):
    """Open ``path`` (transparent ``.gz``), detect format, validate, return a handle.

    Args:
        path: input file path (string or path-like).
        expected_formats: tuple/list of accepted formats (e.g.
            ``(FORMAT_OLMSTED,)`` or ``(FORMAT_AIRR, FORMAT_PCP)``). If
            ``None``, returns whatever format is detected.

    Returns:
        Tuple of ``(file_handle, detected_format)``. The caller is
        responsible for closing the handle (or using it as a context
        manager) and for parsing — CSV and JSON consumers need
        different parsers.

    Raises:
        ValueError: if format is ``unknown``, or if ``expected_formats``
            is set and the detected format isn't in the list.
    """
    detected = detect_file_format(path)
    if detected == FORMAT_UNKNOWN:
        raise ValueError(f"Could not infer format of {path}")
    if expected_formats is not None and detected not in expected_formats:
        raise ValueError(
            f"Expected one of {tuple(expected_formats)} but detected "
            f"{detected!r} for {path}"
        )
    return open_maybe_gzip(path), detected


# --- high-level reads -------------------------------------------------------


def read_olmsted_json(path) -> dict:
    """Open + parse + shape-validate an Olmsted JSON file.

    Validates the file is detected as Olmsted format AND has the
    required top-level keys (``datasets``, ``clones``, ``trees``).
    Raises ``ValueError`` on any failure.
    """
    handle, _ = open_input(path, expected_formats=(FORMAT_OLMSTED,))
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
    handle, _ = open_input(path, expected_formats=(FORMAT_AIRR,))
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
    handle, _ = open_input(path, expected_formats=(FORMAT_PCP,))
    with handle:
        yield from csv.DictReader(handle)


def read_csv_rows(path) -> Iterator[dict]:
    """Generic CSV reader (transparent ``.gz``). No format detection.

    Use for CSVs that aren't auto-detectable as a known Olmsted format —
    currently the mutations CSV consumed by ``merge``.
    """
    with open_maybe_gzip(path) as handle:
        yield from csv.DictReader(handle)


def read_yaml_config(path) -> dict:
    """Open + parse a YAML config file (transparent ``.gz``)."""
    with open_maybe_gzip(path) as handle:
        return yaml.safe_load(handle)


# --- high-level writes ------------------------------------------------------


def write_output(data, path, output_kind: str = "olmsted_json", **opts) -> str:
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
