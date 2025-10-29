# olmsted-cli

Command-line interface and data processing utilities for [Olmsted webapp](https://github.com/matsengrp/olmsted).  The Olmsted web application can be launched locally using the repo, or is available at https://www.olmstedviz.org.

## Overview

`olmsted-cli` is a Python package that processes immunological data from AIRR and PCP formats into the Olmsted JSON format for visualization in the Olmsted web application. It handles B cell receptor repertoire sequencing data, reconstructs phylogenetic trees, and calculates various metrics for clonal family analysis.

### Typical Workflow

1. **Process your data**: Use `olmsted-cli` to convert your AIRR or PCP format files into Olmsted JSON format
2. **Open Olmsted web app**: Launch the application locally or visit https://www.olmstedviz.org
3. **Load your processed files**: Upload the Olmsted JSON file(s)
4. **Visualize**: Explore your B cell lineage data with interactive visualizations

**Example**:
```bash
# Convert your PCP data to Olmsted format
olmsted process -i pcp.csv --tree trees.csv -o olmsted_data.json --compute-metrics
```

### Supported Formats

- **AIRR (Adaptive Immune Receptor Repertoire)**: JSON format following AIRR Community standards
- **PCP (Parent-Child Pair)**: CSV file containing parent-child pairs with separate trees CSV file containing Newick strings

### Output Formats

**Consolidated (default)**: Single JSON file containing all data - recommended for most workflows.
**Unbundled (`--unbundle`)**: Separates data into component files (datasets.json, clones.*.json, tree.*.json) for backwards compatibility with older Olmsted versions.

---

## Installation

### Recommended (pipx)

Install using [pipx](https://pipx.pypa.io/) for isolated environment:

```bash
pipx install olmsted-cli
```

### Standard Installation

Install using pip:

```bash
pip install olmsted-cli
```

### For Development

Clone and install in editable mode:

```bash
git clone https://github.com/matsengrp/olmsted-cli.git
cd olmsted-cli
pip install -e .
```

---

## Quick Start

```bash
# Process AIRR format data (auto-detected)
olmsted process -i data.json -o output/olmsted_data.json

# Process PCP format data with phylogenetic metrics
olmsted process -i sequences.csv --tree trees.csv -o output/data.json --compute-metrics

# Validate output data
olmsted validate output/olmsted_data.json

# Generate summary statistics
olmsted summary output/olmsted_data.json

# Convert between formats
olmsted unbundle -i consolidated.json -o components/
olmsted bundle -i components/ -o consolidated.json
```

---

## Available Commands

### Overview

| Command | Purpose |
|---------|---------|
| **`process`** | This is the primary tool: Converts input AIRR or PCP format data into Olmsted-readable JSON format |
| **`validate`** | Verify data files conform to Olmsted schema |
| **`summary`** | Generate statistics and metadata report for processed data |
| **`split`** | Divide large consolidated files into smaller chunks for performance |
| **`bundle`** | Combine separate component files into a single consolidated JSON file (backwards compatibility tool) |
| **`unbundle`** | Separate a consolidated JSON file into individual component files (backwards compatibility tool) |

---

## Commands

### `process` - Process Data Files

Convert AIRR or PCP format data into Olmsted JSON format.

#### Basic Usage

```bash
# Auto-detect format
olmsted process -i input.json -o output.json

# Explicitly specify format
olmsted process -i input.csv -f pcp -o output.json
```

#### Input/Output Options

| Option | Description |
|--------|-------------|
| `-i, --inputs FILES` | Input file(s). For AIRR: one or more JSON files. For PCP: CSV file |
| `-o, --output FILE` | Output file path for consolidated JSON |
| `--unbundle DIR` | Unbundle output into separate component files (datasets.json, clones.*.json, tree.*.json) for backwards compatibility with Olmsted web app |
| `-f, --format {airr,pcp,auto}` | Input format (default: auto-detect) |
| `-t, --tree FILE` | Trees file for PCP format (optional, can be gzipped) |

#### Processing Options

| Option | Description |
|--------|-------------|
| `-n, --name NAME` | Optional dataset name (stored in metadata) |
| `--validate` | Validate output against schemas before writing |
| `--strict-validation` | Exit with error if validation fails |
| `--seed INT` | Random seed for deterministic UUID generation |
| `-v, --verbose {0,1,2,3}` | Verbosity: 0=quiet, 1=normal (default), 2=verbose, 3=debug |
| `-q, --quiet` | Quiet mode - only show errors (equivalent to `-v 0`) |
| `-w, --warnings` | Show warnings when tree and PCP data disagree (PCP only) |

#### PCP-Specific Options

| Option | Description |
|--------|-------------|
| `--compute-metrics` | Compute LBI, LBR, affinity, and mutation frequency for all nodes |
| `--lbi-tau FLOAT` | Time scale parameter for LBI calculation (default: 0.0125) |
| `--standardize-names` | Rename nodes to standard format: naive (root), Node1, Node2, ... |

#### AIRR-Specific Options

| Option | Description |
|--------|-------------|
| `--naive-name NAME` | Name of naive/root node for tree rooting (default: "naive") |
| `-r, --root-trees` | Root trees using naive node |

#### Examples

```bash
# Auto-detect AIRR format and process multiple input files
olmsted process -i dataset1.json dataset2.json -o combined.json

# Process PCP format with separate trees file and compute metrics
olmsted process -i sequences.csv --tree trees.csv -o output.json --compute-metrics

# Unbundle output into separate component files (for backwards compatibility)
olmsted process -i large_dataset.json --unbundle output_dir/

# Process with validation and custom LBI tau parameter
olmsted process -i data.csv -f pcp -o output.json --compute-metrics --lbi-tau 0.02 --validate

# Quiet mode with standardized node names
olmsted process -i data.csv -o output.json -q --standardize-names
```

#### Output Formats

**Consolidated Output** (single file, default):
```json
{
  "datasets": [...],
  "clones": {...},
  "tree": {...}
}
```

This is the recommended format for most use cases. It's simpler to manage and works with the modern Olmsted web app.

**Unbundled Output** (multiple component files, `--unbundle`):
```
output_dir/
├── datasets.json          # Dataset metadata
├── clones.family1.json    # Clonal family 1 data
├── clones.family2.json    # Clonal family 2 data
├── tree.abc123.json       # Tree for family 1
└── tree.def456.json       # Tree for family 2
```

This format separates the data into distinct component files. Use `--unbundle` when:
- You need backwards compatibility with older Olmsted web app versions
- You're working with legacy pipelines that expect separate files
- You need to inspect or manipulate individual components

**Note**: The Olmsted web app can load both formats, but consolidated format is simpler and recommended for new workflows.

#### Input Formas

**PCP CSV Format**

Expected columns in the PCP CSV file:

| Column | Description |
|--------|-------------|
| `sample_id` | Sample identifier |
| `family` | Clonal family identifier |
| `parent_name` | Parent node name (use "naive" for root) |
| `parent_heavy` | Parent heavy chain sequence |
| `child_name` | Child node name |
| `child_heavy` | Child heavy chain sequence |
| `branch_length` | Branch length between parent and child |
| `depth` | Depth in tree |
| `distance` | Distance from root |
| `v_gene_heavy` | V gene assignment |
| `j_gene_heavy` | J gene assignment |
| `cdr1_codon_start_heavy` | CDR1 start position |
| `cdr1_codon_end_heavy` | CDR1 end position |
| `cdr2_codon_start_heavy` | CDR2 start position |
| `cdr2_codon_end_heavy` | CDR2 end position |
| `cdr3_codon_start_heavy` | CDR3 start position |
| `cdr3_codon_end_heavy` | CDR3 end position |
| `parent_is_naive` | Boolean indicating if parent is naive/root |
| `child_is_leaf` | Boolean indicating if child is a leaf node |

**Trees CSV Format**

Expected columns in the trees file:

| Column | Description |
|--------|-------------|
| `family_name` | Clonal family identifier (must match `family` in PCP CSV) |
| `sample_id` | Sample identifier (must match `sample_id` in PCP CSV) |
| `newick_tree` | Newick format tree string for the family |

---

### `validate` - Validate Data Files

Validate Olmsted/AIRR data files against schemas.

#### Basic Usage

```bash
# Auto-detect file type
olmsted validate data.json

# Validate specific file types
olmsted validate --dataset datasets.json
olmsted validate --clones clones.family1.json clones.family2.json
olmsted validate --tree tree.abc123.json
```

#### Options

| Option | Description |
|--------|-------------|
| `--dataset FILE` | Validate as dataset file |
| `--clone FILE` | Validate as single clone object |
| `--clones FILES` | Validate as clone collection |
| `--tree FILE` | Validate as single tree object |
| `--trees FILES` | Validate as tree collection |
| `-v, --verbose` | Show detailed validation output |
| `--strict` | Exit with error on first validation failure |

#### Examples

```bash
# Validate complete consolidated file
olmsted validate output.json

# Validate unbundled component files
olmsted validate --dataset datasets.json --clones clones.*.json --tree tree.*.json

# Verbose validation with strict mode
olmsted validate -v --strict processed_data.json
```

---

### `summary` - Generate Summary Statistics

Analyze consolidated Olmsted data files and generate summary statistics.

#### Basic Usage

```bash
# Print summary to stdout
olmsted summary data.json

# Save summary to file
olmsted summary data.json -o summary.txt

# Output as JSON
olmsted summary --json data.json
```

#### Options

| Option | Description |
|--------|-------------|
| `-o, --output FILE` | Output file (default: stdout) |
| `--json` | Output summary as JSON format |

#### Example Output

```
Olmsted Data Summary
====================
Datasets: 2
Total Clones: 1,234
Total Tree Nodes: 5,678
  - Leaf Nodes: 2,345
  - Internal Nodes: 3,333

Metrics Available:
  - LBI: Yes
  - LBR: Yes
  - Affinity: Yes
  - Mean Mutation Frequency: Yes
```

---

### `split` - Split Large Files

Split consolidated Olmsted data files into smaller files for better performance.

#### Basic Usage

```bash
# Split into files with max 100 clones each
olmsted split -i large_data.json -o output_dir --max-clones 100

# Split with custom naming
olmsted split -i data.json -o splits --max-clones 50 --base-name my_dataset
```

#### Options

| Option | Description |
|--------|-------------|
| `-i, --input FILE` | Input consolidated JSON file to split |
| `-o, --output-dir DIR` | Output directory for split files |
| `--max-clones INT` | Maximum clones per output file (default: 100) |
| `--base-name NAME` | Base name for output files |

---

### `bundle` - Combine Component Files

Combine separate Olmsted component files into a single consolidated JSON file.

#### Basic Usage

```bash
# Bundle component files from directory
olmsted bundle -i components_dir/ -o consolidated.json

# Bundle with verbose output
olmsted bundle -i components/ -o data.json -v
```

#### Options

| Option | Description |
|--------|-------------|
| `-i, --input-dir DIR` | Input directory containing component files (datasets.json, clones.*.json, tree.*.json) |
| `-o, --output FILE` | Output consolidated JSON file |
| `--datasets FILE` | Path to datasets.json file (optional if using --input-dir) |
| `--clones FILES` | Path(s) to clones.*.json files (optional if using --input-dir, supports wildcards) |
| `--trees FILES` | Path(s) to tree.*.json files (optional if using --input-dir, supports wildcards) |
| `-v, --verbose` | Show verbose output |

#### Examples

```bash
# Bundle all files from a directory
olmsted bundle -i unbundled_components/ -o consolidated.json

# Bundle with explicit file paths and wildcards
olmsted bundle --datasets datasets.json --clones "clones.*.json" --trees "tree.*.json" -o output.json

# Bundle with verbose output to see progress
olmsted bundle -i components/ -o data.json -v
```

#### When to Use

Use `bundle` when you:
- Have unbundled component files and want to create a consolidated file
- Need to convert from legacy unbundled format to modern consolidated format
- Want to combine multiple component files for easier distribution

---

### `unbundle` - Separate Component Files

Separate a consolidated Olmsted JSON file into individual component files for backwards compatibility or component-level inspection.

#### Basic Usage

```bash
# Unbundle consolidated file
olmsted unbundle -i consolidated.json -o output_dir/

# Unbundle with verbose output
olmsted unbundle -i data.json -o components/ -v
```

#### Options

| Option | Description |
|--------|-------------|
| `-i, --input FILE` | Input consolidated JSON file |
| `-o, --output-dir DIR` | Output directory for component files |
| `-v, --verbose` | Show verbose output |

#### Examples

```bash
# Unbundle consolidated file into separate components
olmsted unbundle -i consolidated.json -o components/

# Unbundle with verbose output to see each file being created
olmsted unbundle -i data.json -o unbundled/ -v
```

#### When to Use

Use `unbundle` when you:
- Need backwards compatibility with older Olmsted web app versions
- Want to inspect or modify individual components separately
- Have legacy pipelines that expect unbundled format
- Need to work with specific clonal families or trees independently

**Note**: The `process` command also supports direct unbundling with the `--unbundle` flag during data conversion.

---

## Example Data

The repository includes example data for both formats:

```bash
# Clone repository to access examples
git clone https://github.com/matsengrp/olmsted-cli.git
cd olmsted-cli/example_data

# AIRR format examples
ls airr/

# PCP format examples
ls pcp/
```

---

## Common Workflows

### Processing Raw Data for Visualization

```bash
# 1. Process PCP data with all metrics
olmsted process -i sequences.csv --tree trees.csv \
  -o olmsted_data.json \
  --compute-metrics \
  --lbi-tau 0.0125 \
  -n "My B Cell Study"

# 2. Validate output
olmsted validate olmsted_data.json

# 3. Generate summary
olmsted summary olmsted_data.json

# 4. Upload olmsted_data.json to Olmsted web app
```

### Batch Processing Multiple Datasets

```bash
# Process multiple AIRR files into one combined dataset
olmsted process -i study1.json study2.json study3.json \
  -o combined_dataset.json \
  --validate \
  -v 2

# If needed for backwards compatibility, unbundle the output
olmsted process -i study1.json study2.json study3.json \
  --unbundle output_dir/ \
  --validate

# Split large files for better performance (separate from unbundling)
olmsted split -i combined_dataset.json -o split_output/ --max-clones 500
```

### Processing Large Datasets

```bash
# Process large dataset with all metrics (consolidated output)
olmsted process -i large_sequences.csv --tree large_trees.csv \
  -o large_data.json \
  --compute-metrics \
  -v 2

# Then split into smaller chunks for better web performance
olmsted split -i large_data.json -o split_output/ --max-clones 100

# Or if you need unbundled format for backwards compatibility
olmsted process -i large_sequences.csv --tree large_trees.csv \
  --unbundle large_output/ \
  --compute-metrics \
  -v 2

# Validate unbundled files
olmsted validate --dataset large_output/datasets.json \
  --clones large_output/clones.*.json \
  --tree large_output/tree.*.json
```

---

## Requirements

- **Python**: 3.8 or higher
- **Dependencies** (automatically installed):
  - ete3 ≥3.1.0
  - jsonschema ≥4.0.0
  - lxml ≥4.6.0
  - numpy ≥1.20.0
  - pyyaml ≥6.0
  - scipy ≥1.7.0
  - ntpl ≥0.0.4
  - tqdm ≥4.65.0

---

## Contributing

Contributions are welcome! Please see the [main Olmsted repository](https://github.com/matsengrp/olmsted) for contribution guidelines.

### Development Setup

```bash
# Clone and install with dev dependencies
git clone https://github.com/matsengrp/olmsted-cli.git
cd olmsted-cli
pip install -e ".[dev]"

# Run tests
pytest

# Run linter
ruff check .
```

---

## License

MIT License - see LICENSE file for details.

---

## Links

- **Olmsted Web App**: https://github.com/matsengrp/olmsted
- **Live Demo**: https://olmstedviz.org
- **Issue Tracker**: https://github.com/matsengrp/olmsted-cli/issues
- **AIRR Standards**: https://docs.airr-community.org/

---

*For questions or support, please open an issue on GitHub.*
