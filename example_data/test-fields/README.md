# Field Metadata Test Datasets

Minimal datasets with foobar bogus metrics at every field level and type, for testing the `field_metadata` system end-to-end.

## Foobar Fields by Level

| Level | Field | Type | Description |
|-------|-------|------|-------------|
| **clone** | `foobar_score` | continuous | Bogus numeric clone metric |
| **clone** | `foobar_category` | categorical | Bogus clone grouping |
| **clone** | `foobar_note` | tooltip | Bogus clone description |
| **clone** | `foobar_params` | json | Bogus structured parameters (tooltip display) |
| **clone** | `foobar_path` | categorical | Local file path (auto-skipped) |
| **node** | `foobar_weight` | continuous | Bogus numeric node metric |
| **node** | `foobar_class` | categorical | Bogus node classification |
| **node** | `foobar_description` | tooltip | Bogus node description |
| **node→mutation** | `foobar_per_site_score` | list (inner: continuous) | Per-position scores, demoted to mutation level |
| **node→mutation** | `foobar_sparse_aa` | json (inner: aa) | Sparse per-position AA changes, demoted to mutation level |
| **mutation** | `foobar_impact` | continuous | Bogus numeric mutation score |
| **mutation** | `foobar_tier` | categorical | Bogus mutation tier |
| **mutation** | `child_aa` | aa | Amino acid identity (genetic alphabet type) |
| **mutation** | `parent_aa` | tooltip | Parent amino acid (context only) |

## Files

| File | Format | Contains custom fields in data? |
|------|--------|---------------------------------|
| `olmsted-test-fields.json` | Olmsted JSON | Yes — all levels including mutation |
| `airr-test-fields.json` | AIRR JSON | Yes — clone, node, and mutation levels |
| `pcp-test-fields.csv` + `trees-test-fields.csv` | PCP CSV | Yes — includes JSON-encoded list/dict columns |
| `test-fields-config.yaml` | YAML config | Declares all foobar fields for any format |

## Usage

```bash
# Olmsted JSON: enrich directly
olmsted enrich -i olmsted-test-fields.json -o enriched.json -c test-fields-config.yaml

# AIRR: process with config
olmsted process -i airr-test-fields.json -o output.json -c test-fields-config.yaml

# PCP: process with config
olmsted process -i pcp-test-fields.csv -t trees-test-fields.csv -o output.json -c test-fields-config.yaml

# Dump fields from any format
olmsted build-config -i olmsted-test-fields.json
olmsted build-config -i airr-test-fields.json
olmsted build-config -i pcp-test-fields.csv -t trees-test-fields.csv
```

## New Type Coverage

### list and json types
- `foobar_per_site_score`: A list of floats on each node, length matching the AA sequence (4 positions). Detected as `list` type, then demoted to mutation level with `inner_type: continuous`.
- `foobar_sparse_aa`: A JSON dict with integer keys mapping to single AA characters. Detected as `json` type, then demoted to mutation level with `inner_type: aa`.
- `foobar_params`: A JSON dict with non-integer keys (method names, settings). Stays at clone level as `json` with `display: tooltip`.

### Path detection
- `foobar_path`: Contains local file paths (e.g. `/data/raw/sample-1/clone-A.fasta`). Auto-detected and suggested for skip in build-config output.
