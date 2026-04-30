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
| **mutations** | Olmsted JSON | 3 clones | Mutation-level annotation example (DASM2 surprise-analysis subset) | Heavy chain (IGH) |

## Mutations Data

The `mutations/` directory contains a subset of mutation-level annotation data (sourced from DASM2 surprise analysis), pre-built in Olmsted JSON format. This data exercises:

- **Clone-level fields**: `mean_surprise_mutsel`, `num_mutations` (in addition to standard fields)
- **Mutation-level fields**: `surprise_mutsel`, `surprise_neutral`, `selection_contribution`, `region` (per-mutation scores in `mutations` arrays on tree nodes)
- 3 clones across 3 subjects (d1, d3, d4) with different V genes

### Files

| File | Description |
|------|-------------|
| `mutations-olmsted.json` | 3-clone Olmsted JSON subset (392 KB) |
| `mutations-config.yaml` | YAML config with custom field labels |

### Example Workflow

```bash
# 1. See what fields are in the data
olmsted build-config -i example-data/mutations/input-olmsted.json

# 2. Enrich with field_metadata using the provided config
olmsted enrich -i example-data/mutations/input-olmsted.json \
  -o enriched.json -c example-data/mutations/mutations-config.yaml

# 3. Or generate your own config, edit it, then enrich
olmsted build-config -i example-data/mutations/input-olmsted.json -o my_config.yaml
# ... edit my_config.yaml ...
olmsted enrich -i example-data/mutations/input-olmsted.json \
  -o enriched.json -c my_config.yaml
```

## Standard Processing Examples

```bash
# Process AIRR format
olmsted process -i example-data/airr/input-airr.json -o output.json

# Process PCP format with trees
olmsted process -i example-data/pcp/input-pcp.csv -t example-data/pcp/input-trees.csv -o output.json

# Process PCP with metrics
olmsted process -i example-data/pcp/input-pcp.csv -t example-data/pcp/input-trees.csv \
  -o output.json --compute-metrics

# Process paired heavy/light chain data
olmsted process -i example-data/pcp-paired/input-pcp.csv \
  -t example-data/pcp-paired/input-trees.csv -o output.json

# Use a YAML config file
olmsted process -c example-data/mutations/mutations-config.yaml
```

## Golden Data

Each dataset folder includes a single consolidated golden used by the test suite:

- `{folder}-olmsted-golden.json` — expected consolidated output for that dataset (e.g., `airr/airr-olmsted-golden.json`, `pcp/pcp-olmsted-golden.json`)

The `airr/` and `pcp/` folders additionally carry a `split-golden-data/` directory pinning the legacy split-format output (`datasets.json`, `clones.*.json`, `tree.*.json`). These are maintained as integrity coverage for the `--split-files` CLI flag and will be removed if/when split-format support is officially dropped.

Regenerate after output-format changes — see CLAUDE.md or DEVELOPMENT.md for the regeneration commands.
