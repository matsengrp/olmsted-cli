# olmsted-cli

Command-line interface and data processing utilities for [Olmsted](https://github.com/matsengrp/olmsted).

This package contains all the data processing scripts used by the Olmsted web application,
packaged as a standalone CLI tool for processing AIRR and PCP format immunological data.

## Installation

### Recommended (pipx)

Install using [pipx](https://pipx.pypa.io/) for isolated environment:

```bash
pipx install olmsted-cli
```

### For development

```bash
git clone https://github.com/matsengrp/olmsted-cli.git
cd olmsted-cli
pip install -e .
```

## Usage

The `olmsted` CLI provides commands for processing AIRR and PCP format data.

### Automatic format detection

```bash
# Process data with automatic format detection
olmsted process -i data.json -o output/
olmsted process -i data.csv -o output/

# Process data with explicit format
olmsted process -i data.json -o output/ -f airr
olmsted process -i data.csv -o output/ -f pcp
```


## Requirements

- Python 3.8+
- Dependencies are automatically installed during installation
