# Olmsted Example Datasets

This directory contains example datasets for testing and demonstrating the Olmsted visualization tool. These datasets are used in tests and serve as references for the expected output format.

## Dataset Overview

| Dataset | Format | Families/Clones | Description | Chain Type |
|---------|--------|-----------------|-------------|------------|
| **airr** | AIRR JSON | 8 clones | Real BCR data from AIRR Community format | Heavy chain (IGH) |
| **pcp** | PCP CSV | 8 families | Real BCR data, heavy chain only dataset | Heavy chain (IGH) |
| **pcp-light** | PCP CSV | 8 families | Real BCR data, light chain only dataset | Light chain (IGK) |
| **pcp-paired** | PCP CSV | 8 families | Real BCR data, paired heavy and light chain data | Heavy (IGH) + Light (IGL/IGK) |
| **pcp-byhand** | PCP CSV | 3 families | Artificial dataset for testing webapp | Heavy, Light, and Paired |
| **surprise** | Olmsted JSON | 3 clones | DASM2 surprise analysis data with per-mutation surprise scores | Heavy chain (IGH) |

## Surprise Analysis Data

The `surprise/` directory contains a subset of DASM2 mutation-selection surprise analysis data, pre-built in Olmsted JSON format. This data includes:

- **Clone-level fields**: `mean_surprise_mutsel`, `num_mutations` (in addition to standard fields)
- **Mutation-level fields**: `surprise_mutsel`, `surprise_neutral`, `selection_contribution`, `region` (per-mutation scores in `surprise_mutations` arrays on tree nodes)
- 3 clones across 3 subjects (d1, d3, d4) with different V genes

### Files

| File | Description |
|------|-------------|
| `surprise_subset.json` | 3-clone Olmsted JSON subset (392 KB) |
| `surprise_config.yaml` | YAML config with custom field labels for surprise data |

### Example Workflow

```bash
# 1. See what fields are in the data
olmsted build-config -i example_data/surprise/surprise_subset.json

# 2. Enrich with field_metadata using the provided config
olmsted enrich -i example_data/surprise/surprise_subset.json \
  -o enriched.json -c example_data/surprise/surprise_config.yaml

# 3. Or generate your own config, edit it, then enrich
olmsted build-config -i example_data/surprise/surprise_subset.json -o my_config.yaml
# ... edit my_config.yaml ...
olmsted enrich -i example_data/surprise/surprise_subset.json \
  -o enriched.json -c my_config.yaml
```

## Standard Processing Examples

```bash
# Process AIRR format
olmsted process -i example_data/airr/airr.json -o output.json

# Process PCP format with trees
olmsted process -i example_data/pcp/pcp.csv -t example_data/pcp/trees.csv -o output.json

# Process PCP with metrics
olmsted process -i example_data/pcp/pcp.csv -t example_data/pcp/trees.csv \
  -o output.json --compute-metrics

# Process paired heavy/light chain data
olmsted process -i example_data/pcp-paired/pcp-paired.csv \
  -t example_data/pcp-paired/trees-paired.csv -o output.json

# Use a YAML config file
olmsted process -c example_data/surprise/surprise_config.yaml
```

## Golden Data

The `airr/` and `pcp/` directories include golden data files used by the test suite:

- `consolidated_golden_data.json` -- Expected consolidated output
- `split_golden_data/` -- Expected split-file output (datasets.json, clones.\*.json, tree.\*.json)

These are regenerated when the output format changes (e.g., adding `field_metadata`).
