#!/usr/bin/env python3
"""CLI wrapper for Olmsted data processing scripts."""

import sys
import os
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
  
  # Process AIRR format data  
  olmsted airr -i data.json -o output/
  
  # Process PCP format data
  olmsted pcp -i data.csv -o output/
        """
    )
    
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    # Process command (unified)
    process_parser = subparsers.add_parser('process', help='Process data with automatic format detection')
    # Don't define arguments here - let the underlying script handle them
    
    # AIRR command
    airr_parser = subparsers.add_parser('airr', help='Process AIRR format data')
    # Don't define arguments here - let the underlying script handle them
    
    # PCP command  
    pcp_parser = subparsers.add_parser('pcp', help='Process PCP format data')
    # Don't define arguments here - let the underlying script handle them
    
    # Parse only the command, not the full arguments
    if len(sys.argv) > 1 and sys.argv[1] in ['process', 'airr', 'pcp']:
        command = sys.argv[1]
        # Remove the script name and command from sys.argv so the underlying scripts see clean arguments
        sys.argv = [sys.argv[0]] + sys.argv[2:]
        
        if command == 'process':
            process_data()
        elif command == 'airr':
            process_airr()
        elif command == 'pcp':
            process_pcp()
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


def process_airr():
    """Run the process_airr_data.py script."""
    original_dir = os.getcwd()
    package_dir = Path(__file__).parent.parent
    
    try:
        os.chdir(package_dir)
        from olmsted_cli import process_airr_data
        process_airr_data.main()
    finally:
        os.chdir(original_dir)


def process_pcp():
    """Run the process_pcp_data.py script."""
    original_dir = os.getcwd()
    package_dir = Path(__file__).parent.parent
    
    try:
        os.chdir(package_dir)
        from olmsted_cli import process_pcp_data
        process_pcp_data.main()
    finally:
        os.chdir(original_dir)


if __name__ == "__main__":
    main()