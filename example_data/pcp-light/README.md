# PCP Light-Chain-Only Test Dataset

Light-chain-only dataset extracted from real paired B cell receptor data for testing light chain visualization and processing.

## Dataset Structure

- **8 families**: 6 kappa (IGK) + 2 lambda (IGL)
- **40 PCP rows**: 5 rows per family (average)
- **9 Newick trees**: Phylogenetic trees for each family
- **Format**: Light chain data in primary (non-`*_light`) columns
  - All `*_light` columns removed
  - Light chain sequences in `parent_heavy`/`child_heavy` columns
  - Chain type inferred from V gene calls (IGKV*, IGLV*)

## Families

| Family ID | Chain Type | Sample | Rows |
|-----------|------------|--------|------|
| 33424-igk-33424 | kappa | d1 | 5 |
| 13545-igk-13545 | kappa | d1 | 5 |
| 179396-igk-179396 | kappa | d1 | 5 |
| 152208-igk-152208 | kappa | d1 | 5 |
| 159055-igl-5792 | lambda | d1 | 5 |
| 146931-igk-146931 | kappa | d1 | 5 |
| 155785-igl-2522 | lambda | d1 | 5 |
| 172841-igk-172841 | kappa | d1 | 5 |

## Source

Extracted from `example_data/_pcp-paired/` - a large real-world paired B cell receptor dataset.

Selection criteria:
- Light chain sequences present
- Family size 5 rows (moderate complexity)
- Mix of kappa and lambda chains

## Usage

Convert to Olmsted JSON:
```bash
cd olmsted-cli
python -m olmsted_cli.process_pcp_data \
  -i example_data/pcp-light/pcp.csv \
  -t example_data/pcp-light/trees.csv \
  -o example_data/pcp-light/output.json
```

## Files

- `pcp.csv` - Parent-child pair data (40 rows)
- `trees.csv` - Newick phylogenetic trees (9 trees)
- `output.json` - Generated Olmsted JSON (after conversion)
- `README.md` - This file

## Purpose

This dataset is used for:
- Testing light-chain-only data processing
- Verifying kappa vs lambda chain handling
- Testing light chain visualization features
- Ensuring paired data processing doesn't break single-chain workflows

## Technical Notes

**Locus Field**: In the converted JSON output, the `locus` field may show as `"igh"` (heavy) because the data is in the primary (non-`*_light`) columns. The actual chain type should be determined by the V gene calls:
- `IGKV*` genes indicate kappa light chains
- `IGLV*` genes indicate lambda light chains

This is the standard format for light-chain-only PCP data where the `*_light` columns are not used.

---

*Created: 2025-12-12*
*Source: _pcp-paired dataset*
