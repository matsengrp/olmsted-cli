#!/usr/bin/env python3
"""CLI wrapper for Olmsted data processing scripts."""

import argparse
import sys

from olmsted_cli import build_config, enrich, process_data, split, summary, validate

# Dispatch table: command name -> (help text, handler module)
COMMANDS = {
    "process": ("Process data with automatic format detection", process_data),
    "validate": ("Validate data files against AIRR/Olmsted schemas", validate),
    "summary": ("Generate summary statistics for consolidated data files", summary),
    "split": ("Split consolidated data files into smaller files", split),
    "enrich": ("Add field_metadata to existing Olmsted JSON files", enrich),
    "build-config": ("Generate a YAML config from your data for editing", build_config),
}


def main():
    """Main entry point for the olmsted CLI."""
    parser = argparse.ArgumentParser(
        description="Olmsted CLI - Process AIRR and PCP format data",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Auto-detect format and process
  olmsted process -i data.json -o output/

  # Build a config from your data, then edit and use it
  olmsted build-config -i data.csv -t trees.csv -o config.yaml

  # Process with a config file
  olmsted process -c config.yaml
        """,
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
