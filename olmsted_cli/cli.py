#!/usr/bin/env python3
"""CLI wrapper for Olmsted data processing scripts."""

import argparse
import sys

from olmsted_cli import build_config, enrich, process_data, split, summary, validate
from olmsted_cli.version import version_string

# Dispatch table: command name -> (help text, handler module)
# Ordered by typical workflow: build-config → process → enrich, then utilities
COMMANDS = {
    "process": ("Convert input data (AIRR/PCP) to Olmsted JSON format", process_data),
    "build-config": ("Generate a YAML config from your data for editing", build_config),
    "enrich": ("Add field_metadata to existing Olmsted JSON files", enrich),
    "validate": ("Validate data files against schemas", validate),
    "summary": ("Generate summary statistics for Olmsted JSON files", summary),
    "split": ("Split Olmsted JSON files into smaller files (legacy)", split),
}


def main():
    """Main entry point for the olmsted CLI."""
    parser = argparse.ArgumentParser(
        prog="olmsted",
        description="Olmsted CLI — Convert immunological data to Olmsted JSON for visualization",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Typical workflow:
  1. olmsted build-config -i data.csv -t trees.csv -o config.yaml
  2. Edit config.yaml (adjust fields, labels, types)
  3. olmsted process -c config.yaml

Or process directly:
  olmsted process -i data.json -o output.json
        """,
    )
    parser.add_argument(
        "--version", action="version", version=f"olmsted-cli {version_string()}"
    )

    subparsers = parser.add_subparsers(dest="command", help="Available commands")
    for cmd_name, (help_text, _) in COMMANDS.items():
        subparsers.add_parser(cmd_name, help=help_text)

    if len(sys.argv) > 1 and sys.argv[1] in COMMANDS:
        command = sys.argv[1]
        sys.argv = [sys.argv[0]] + sys.argv[2:]
        COMMANDS[command][1].main()
    else:
        parser.parse_args()


if __name__ == "__main__":
    main()
