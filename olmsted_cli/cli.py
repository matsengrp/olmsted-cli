#!/usr/bin/env python3
"""CLI wrapper for Olmsted data processing scripts."""

import os
import sys
from pathlib import Path


def main():
    """Main entry point for the olmsted CLI."""
    import argparse

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
    # Don't define arguments here - let the underlying script handle them

    # Parse only the command, not the full arguments
    if len(sys.argv) > 1:
        command = sys.argv[1]
        
        if command == "process":
            # Remove the script name and command from sys.argv so the underlying scripts see clean arguments
            sys.argv = [sys.argv[0]] + sys.argv[2:]
            process_data()
        elif command == "validate":
            # Remove the script name and command from sys.argv
            sys.argv = [sys.argv[0]] + sys.argv[2:]
            validate_data()
        else:
            parser.parse_args()  # This will show help or error
    else:
        parser.parse_args()  # This will show help or error


def process_data():
    """Run the unified process_data.py script."""
    # Save current directory and change to package root for schema access
    original_dir = os.getcwd()
    package_dir = Path(__file__).parent.parent

    try:
        os.chdir(package_dir)
        # Import from the package
        from olmsted_cli import process_data

        process_data.main()
    finally:
        os.chdir(original_dir)


def validate_data():
    """Run the validate command."""
    # Import from the package
    from olmsted_cli import validate_command
    
    validate_command.main()


if __name__ == "__main__":
    main()
