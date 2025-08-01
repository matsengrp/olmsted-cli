#!/usr/bin/env python3
"""CLI wrapper for Olmsted data processing scripts."""

import sys
import os
from pathlib import Path

# Add the bin directory to Python path
BIN_DIR = Path(__file__).parent.parent / "bin"
sys.path.insert(0, str(BIN_DIR))


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
    process_parser.add_argument('-i', '--input', required=True, help='Input file')
    process_parser.add_argument('-o', '--output', required=True, help='Output directory')
    
    # AIRR command
    airr_parser = subparsers.add_parser('airr', help='Process AIRR format data')
    airr_parser.add_argument('-i', '--input', required=True, help='Input AIRR JSON file')
    airr_parser.add_argument('-o', '--output', required=True, help='Output directory')
    
    # PCP command  
    pcp_parser = subparsers.add_parser('pcp', help='Process PCP format data')
    pcp_parser.add_argument('-i', '--input', required=True, nargs='+', help='Input PCP CSV file(s)')
    pcp_parser.add_argument('-o', '--output', required=True, help='Output directory')
    
    args = parser.parse_args()
    
    if args.command == 'process':
        process_data()
    elif args.command == 'airr':
        process_airr()
    elif args.command == 'pcp':
        process_pcp()
    else:
        parser.print_help()
        sys.exit(1)


def process_data():
    """Run the unified process_data.py script."""
    # Handle help before importing
    if '--help' in sys.argv or '-h' in sys.argv:
        # Pass through to the actual script
        sys.argv[0] = 'olmsted-process'
    
    # Save current directory and change to package root for schema access
    original_dir = os.getcwd()
    package_dir = Path(__file__).parent.parent
    
    try:
        os.chdir(package_dir)
        # Import here to avoid import errors before path is set up
        import process_data
        process_data.main()
    finally:
        os.chdir(original_dir)


def process_airr():
    """Run the process_airr_data.py script."""
    # Handle help before importing
    if '--help' in sys.argv or '-h' in sys.argv:
        # Pass through to the actual script
        sys.argv[0] = 'olmsted-airr'
    
    original_dir = os.getcwd()
    package_dir = Path(__file__).parent.parent
    
    try:
        os.chdir(package_dir)
        import process_airr_data
        process_airr_data.main()
    finally:
        os.chdir(original_dir)


def process_pcp():
    """Run the process_pcp_data.py script."""
    # Handle help before importing
    if '--help' in sys.argv or '-h' in sys.argv:
        # Pass through to the actual script
        sys.argv[0] = 'olmsted-pcp'
    
    original_dir = os.getcwd()
    package_dir = Path(__file__).parent.parent
    
    try:
        os.chdir(package_dir)
        import process_pcp_data
        process_pcp_data.main()
    finally:
        os.chdir(original_dir)


if __name__ == "__main__":
    main()