# PCP Paired Test Dataset

Paired heavy and light chain dataset extracted from real B cell receptor data for testing paired chain visualization and processing.

## Dataset Structure

- **8 families**: All paired (heavy + light chains)
- **26 PCP rows**: 3-5 rows per family (average: 3.25)
- **8 Newick trees**: Phylogenetic trees for each family
- **16 clones total**: 8 heavy + 8 light (paired)

## Families

| Family ID | Chain Types | Rows |
|-----------|-------------|------|
| 10002-igk-10002 | IGH + IGK | 3 |
| 100022-igk-100022 | IGH + IGK | 3 |
| 100025-igk-100025 | IGH + IGK | 3 |
| 100043-igk-100043 | IGH + IGK | 5 |
| 100044-igk-100044 | IGH + IGK | 3 |
| 129355-igk-129355 | IGH + IGK | 3 |
| 269773-igl-103817 | IGH + IGL | 3 |
| 7420-igk-7420 | IGH + IGK | 3 |

**Chain distribution**: 7 kappa (IGK) + 1 lambda (IGL)

## Source

Extracted from `example_data/_pcp-paired/` - a large real-world paired B cell receptor dataset.

Selection criteria:
- All families have paired heavy and light chain data
- Moderate family size (3-5 rows)
- Mix of kappa and lambda chains
- Real sequences from sample "d1"

## Format

Standard PCP paired format with both heavy and light chain columns:
- **Heavy chain**: `parent_heavy`, `child_heavy`, `v_gene_heavy`, `j_gene_heavy`, etc.
- **Light chain**: `parent_light`, `child_light`, `v_gene_light`, `j_gene_light`, etc.
- **Shared topology**: `branch_length`, `depth`, `distance`, `parent_is_naive`, `child_is_leaf`

## Usage

Convert to Olmsted JSON:
```bash
cd olmsted-cli
python -m olmsted_cli.process_pcp_data \
  -i example_data/pcp-paired/pcp.csv \
  -t example_data/pcp-paired/trees.csv \
  -o example_data/pcp-paired/output.json
```

Expected output:
- 16 clones (8 heavy + 8 light)
- All clones marked as `is_paired: true`
- Matching `pair_id` for heavy/light pairs

## Files

- `pcp.csv` - Parent-child pair data (26 rows)
- `trees.csv` - Newick phylogenetic trees (8 trees)
- `pcp-paired-olmsted-golden.json` - Reference golden output for testing
- `README.md` - This file

## Purpose

This dataset is used for:
- Testing paired heavy/light chain data processing
- Verifying paired clone detection and linking
- Testing synchronized visualization of paired chains
- Ensuring light chain type detection (kappa vs lambda)

---

*Created: 2025-12-12*
*Source: _pcp-paired dataset*
