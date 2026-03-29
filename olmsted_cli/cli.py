#!/usr/bin/env python3
"""CLI wrapper for Olmsted data processing scripts."""

import argparse
import sys

from olmsted_cli import dump_fields, enrich, process_data, split, summary, validate


def main():
    """Main entry point for the olmsted CLI."""

    parser = argparse.ArgumentParser(
        description="Olmsted CLI - Process AIRR and PCP format data",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Auto-detect format and process
  olmsted process -i data.json -o output/

  # Process AIRR format data explicitly
  olmsted process -f airr -i data.json -o output/

  # Process PCP format data explicitly
  olmsted process -f pcp -i data.csv -o output/
        """,
    )

    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # Process command (unified)
    subparsers.add_parser(
        "process", help="Process data with automatic format detection"
    )

    # Validate command
    subparsers.add_parser(
        "validate", help="Validate data files against AIRR/Olmsted schemas"
    )
    
    # Summary command
    subparsers.add_parser(
        "summary", help="Generate summary statistics for consolidated data files"
    )
    
    # Split command
    subparsers.add_parser(
        "split", help="Split consolidated data files into smaller files"
    )

    # Enrich command
    subparsers.add_parser(
        "enrich", help="Add field_metadata to existing Olmsted JSON files"
    )

    # Dump-fields command
    subparsers.add_parser(
        "dump-fields", help="Extract all fields from data into a YAML config for editing"
    )
    # Don't define arguments here - let the underlying script handle them

    # Parse only the command, not the full arguments
    if len(sys.argv) > 1:
        command = sys.argv[1]

        if command == "process":
            # Remove the script name and command from sys.argv so the underlying scripts see clean arguments
            sys.argv = [sys.argv[0]] + sys.argv[2:]
            process_data_command()
        elif command == "validate":
            # Remove the script name and command from sys.argv
            sys.argv = [sys.argv[0]] + sys.argv[2:]
            validate_data_command()
        elif command == "summary":
            # Remove the script name and command from sys.argv
            sys.argv = [sys.argv[0]] + sys.argv[2:]
            summary_data_command()
        elif command == "split":
            # Remove the script name and command from sys.argv
            sys.argv = [sys.argv[0]] + sys.argv[2:]
            split_data_command()
        elif command == "enrich":
            # Remove the script name and command from sys.argv
            sys.argv = [sys.argv[0]] + sys.argv[2:]
            enrich_data_command()
        elif command == "dump-fields":
            sys.argv = [sys.argv[0]] + sys.argv[2:]
            dump_fields_command()
        else:
            parser.parse_args()  # This will show help or error
    else:
        parser.parse_args()  # This will show help or error


def process_data_command():
    """Run the unified process_data.py script."""
    # This preserves the user's current working directory for file resolution
    process_data.main()


def validate_data_command():
    """Run the validate command."""
    validate.main()


def summary_data_command():
    """Run the summary command."""
    summary.main()


def split_data_command():
    """Run the split command."""
    split.main()


def enrich_data_command():
    """Run the enrich command."""
    enrich.main()


def dump_fields_command():
    """Run the dump-fields command."""
    dump_fields.main()


if __name__ == "__main__":
    main()
