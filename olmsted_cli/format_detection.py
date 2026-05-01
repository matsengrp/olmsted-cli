"""File format detection for Olmsted input data.

Detects whether input files are AIRR JSON, PCP CSV, or existing Olmsted
JSON based on file extension and content inspection.  This module depends
only on constants and the standard library, so it can be imported by any
other module without creating circular dependencies.
"""

import json
from pathlib import Path

from .constants import FORMAT_AIRR, FORMAT_OLMSTED, FORMAT_PCP, FORMAT_UNKNOWN
from .utils import open_maybe_gzip, vprint


def detect_file_format(file_path):
    """Automatically detect the file format based on file extension and content.

    Args:
        file_path: Path to the input file.

    Returns:
        str: Detected format (FORMAT_AIRR, FORMAT_PCP, FORMAT_OLMSTED,
        or FORMAT_UNKNOWN).
    """
    file_path = Path(file_path)

    # CSV files are always PCP
    if file_path.suffix.lower() == ".csv":
        return FORMAT_PCP
    if file_path.suffix.lower() == ".gz" and file_path.stem.endswith(".csv"):
        return FORMAT_PCP

    # JSON files need content inspection to distinguish AIRR from Olmsted
    if file_path.suffix.lower() == ".json" or (
        file_path.suffix.lower() == ".gz" and file_path.stem.endswith(".json")
    ):
        try:
            with open_maybe_gzip(file_path) as fh:
                data = json.load(fh)

            if isinstance(data, dict):
                # Explicit format tag in metadata
                metadata = data.get("metadata", {})
                if isinstance(metadata, dict) and metadata.get("format") == FORMAT_OLMSTED:
                    return FORMAT_OLMSTED
                # Heuristic fallback: Olmsted JSON has "datasets" and "metadata"
                if "datasets" in data and "metadata" in data:
                    return FORMAT_OLMSTED
                # AIRR JSON has "clones" with "dataset_id" or standard AIRR keys
                if "dataset_id" in data or "clones" in data or "ident" in data:
                    return FORMAT_AIRR
            elif isinstance(data, list):
                # Multi-dataset AIRR
                return FORMAT_AIRR
        except (json.JSONDecodeError, OSError, ValueError):
            pass

    # If extension doesn't help, try to peek at content for CSV
    try:
        with open_maybe_gzip(file_path) as fh:
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
