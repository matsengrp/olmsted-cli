# olmsted-cli

Command-line interface for [Olmsted](https://github.com/matsengrp/olmsted) data processing.

## Installation

### From source

```bash
git clone --recursive https://github.com/matsengrp/olmsted-cli.git
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

### Direct script usage

You can also run the processing scripts directly:

```bash
# Unified processor with auto-detection
olmsted-process -i data.json -o output/

# AIRR processor
olmsted-airr -i data.json -o output/ --validate

# PCP processor  
olmsted-pcp -i data.csv -o output/ --seed 42
```

## Requirements

- Python 3.8+
- See `requirements.txt` for package dependencies

## Development

```bash
# Clone with submodules
git clone --recursive https://github.com/matsengrp/olmsted-cli.git
cd olmsted-cli

# Create virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install in development mode
pip install -e ".[dev]"

# Run tests
pytest
```

## License

This project is licensed under the MIT License - see the [olmsted](https://github.com/matsengrp/olmsted) repository for details.
