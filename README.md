# olmsted-cli

Command-line interface and data processing utilities for [Olmsted](https://github.com/matsengrp/olmsted).

This package contains all the data processing scripts used by the Olmsted web application, 
packaged as a standalone CLI tool for processing AIRR and PCP format immunological data.

## Installation

### From source

```bash
git clone https://github.com/matsengrp/olmsted-cli.git
cd olmsted-cli
pip install -e .
```

### From PyPI (when published)

```bash
pip install olmsted-cli
```

## Usage

The `olmsted` CLI provides commands for processing AIRR and PCP format data.

### Automatic format detection

```bash
# Process data with automatic format detection
olmsted process -i data.json -o output/
olmsted process -i data.csv -o output/
```

### AIRR format processing

```bash
# Process AIRR JSON data
olmsted airr -i airr_data.json -o output/
```

### PCP format processing

```bash
# Process PCP CSV data
olmsted pcp -i pcp_data.csv -o output/

# With separate trees file
olmsted pcp -i pcp_data.csv pcp_trees.csv -o output/
```


## Requirements

- Python 3.8+
- See `requirements.txt` for package dependencies

## Development

```bash
# Clone repository
git clone https://github.com/matsengrp/olmsted-cli.git
cd olmsted-cli

# Create virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install in development mode
pip install -e ".[dev]"

# Run tests
pytest
```

## Project Structure

```
olmsted-cli/
├── olmsted_cli/            # Python package with all modules
│   ├── __init__.py
│   ├── cli.py              # CLI wrapper
│   ├── process_data.py     # Unified processor with format detection
│   ├── process_airr_data.py # AIRR format processor
│   ├── process_pcp_data.py  # PCP format processor
│   └── process_utils.py    # Shared utilities
├── example_data/           # Test data
│   ├── airr/
│   └── pcp/
├── tests/                  # Pytest test suite
├── airr-standards/         # AIRR schema files
└── data_schema/            # Olmsted validation schemas
```

## License

This project is licensed under the MIT License - see the [olmsted](https://github.com/matsengrp/olmsted) repository for details.
