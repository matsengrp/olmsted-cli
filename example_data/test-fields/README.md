# Field Metadata Test Datasets

Minimal datasets with foobar bogus metrics at every field level and type, for testing the `field_metadata` system end-to-end.

## Foobar Fields by Level

| Level | Field | Type | Description |
|-------|-------|------|-------------|
| **clone** | `foobar_score` | continuous | Bogus numeric clone metric |
| **clone** | `foobar_category` | categorical | Bogus clone grouping |
| **clone** | `foobar_note` | tooltip | Bogus clone description |
| **node** | `foobar_weight` | continuous | Bogus numeric node metric |
| **node** | `foobar_class` | categorical | Bogus node classification |
| **node** | `foobar_description` | tooltip | Bogus node description |
| **mutation** | `foobar_impact` | continuous | Bogus numeric mutation score |
| **mutation** | `foobar_tier` | categorical | Bogus mutation tier |
| **mutation** | `child_aa` | aa | Amino acid identity (genetic alphabet type) |
| **mutation** | `parent_aa` | tooltip | Parent amino acid (context only) |

## Files

| File | Format | Contains custom fields in data? |
|------|--------|---------------------------------|
| `olmsted-test-fields.json` | Olmsted JSON | Yes — all levels including mutation |
| `airr-test-fields.json` | AIRR JSON | Yes — clone, node, and mutation levels |
| `pcp-test-fields.csv` + `trees-test-fields.csv` | PCP CSV | No — PCP parser only passes known columns |
| `test-fields-config.yaml` | YAML config | Declares all foobar fields for any format |

## Usage

```bash
# Olmsted JSON: enrich directly
olmsted enrich -i olmsted-test-fields.json -o enriched.json -c test-fields-config.yaml

# AIRR: process with config
olmsted process -i airr-test-fields.json -o output.json -c test-fields-config.yaml

# PCP: process with config (foobar fields declared but not in CSV data)
olmsted process -i pcp-test-fields.csv -t trees-test-fields.csv -o output.json -c test-fields-config.yaml

# Dump fields from any format
olmsted dump-fields -i olmsted-test-fields.json
olmsted dump-fields -i airr-test-fields.json
olmsted dump-fields -i pcp-test-fields.csv -t trees-test-fields.csv
```

## Notes on PCP and Mutation-Level Data

The PCP parser extracts only known columns — extra columns in the CSV are not passed through to the output. Custom node/branch/mutation-level fields require either:

1. Adding them after processing via `olmsted enrich` with a config
2. Working in Olmsted JSON or AIRR format where arbitrary fields are preserved

Mutation-level data (like `surprise_mutations`) is stored as arrays on tree nodes in the Olmsted JSON output. There is no native PCP column format for per-mutation data. This could theoretically be supported with serialized columns (e.g., JSON arrays in CSV cells), but that is not currently implemented.
