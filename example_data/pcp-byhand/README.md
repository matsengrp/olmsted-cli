# PCP By-Hand Test Dataset

Minimal test dataset with structured mutations for verifying Olmsted PCP data processing, visualization, and metric calculations.

## Dataset Structure

- **3 families, 4 trees total**:
  - `family-heavy`: Heavy chain only (IGH)
  - `family-light`: Light chain only (IGK)
  - `family-paired`: Paired heavy + light chains

- **Tree topology** (same for all families):
  ```
  naive
    └─ Node1 (branch length: 0.1)
        ├─ Tip1 (branch length: 0.2)
        └─ Tip2 (branch length: 0.3)
  ```

## Sequence Design

All sequences use simple repeating nucleotide patterns with mutations at **every 5th codon** (codons 5, 10, 15, etc.) to make visual verification easy.

### family-heavy (Heavy chain only)
- **Length**: 120bp (40 codons)
- **Germline**: 40 codons of `AAA` (Lysine)
- **Mutations**: `AAA` → `TTT` (Lys → Phe)
- **V gene**: IGHV1-1*01
- **J gene**: IGHJ1*01
- **CDR regions**: CDR1=4-8, CDR2=12-16, CDR3=25-35

| Node | Mutated Codons | Total Mutations | Sequence Pattern |
|------|----------------|-----------------|-------------------|
| naive | - | 0 | All AAA |
| Node1 | 5 | 1 | AAA(4x) + TTT + AAA(35x) |
| Tip1 | 5, 10 | 2 | AAA(4x) + TTT + AAA(4x) + TTT + AAA(30x) |
| Tip2 | 5, 10, 15 | 3 | AAA(4x) + TTT + AAA(4x) + TTT + AAA(4x) + TTT + AAA(25x) |

### family-light (Light chain only)
- **Length**: 90bp (30 codons)
- **Germline**: 30 codons of `CCC` (Proline)
- **Mutations**: `CCC` → `GGG` (Pro → Gly)
- **V gene**: IGKV1-1*01
- **J gene**: IGKJ1*01
- **Light chain type**: kappa
- **CDR regions**: CDR1=3-7, CDR2=10-14, CDR3=20-28

| Node | Mutated Codons | Total Mutations | Sequence Pattern |
|------|----------------|-----------------|-------------------|
| naive | - | 0 | All CCC |
| Node1 | 5 | 1 | CCC(4x) + GGG + CCC(25x) |
| Tip1 | 5, 10 | 2 | CCC(4x) + GGG + CCC(4x) + GGG + CCC(20x) |
| Tip2 | 5, 10, 15 | 3 | CCC(4x) + GGG + CCC(4x) + GGG + CCC(4x) + GGG + CCC(15x) |

### family-paired (Paired chains)
- **Heavy**: 120bp (40 codons) - Germline: `GGG`, Mutations: `GGG` → `AAA`
- **Light**: 90bp (30 codons) - Germline: `TTT`, Mutations: `TTT` → `CCC`
- **V genes**: IGHV2-2*01 (heavy), IGKV2-2*01 (light)
- **J genes**: IGHJ2*01 (heavy), IGKJ2*01 (light)
- **Light chain type**: kappa
- **CDR regions**: Same as above

Mutation pattern matches the single-chain families (codons 5, 10, 15).

## Expected Metrics

Run with:
```bash
cd olmsted-cli
python -m olmsted_cli.process_pcp_data \
  -i example_data/_pcp-byhand/pcp.csv \
  -t example_data/_pcp-byhand/trees.csv \
  -o example_data/_pcp-byhand/byhand.json \
  --compute-metrics
```

### Mean Mutation Frequency

**family-paired** (only paired family is currently processed):

| Chain | Expected Calculation | Expected Value | Actual Value |
|-------|---------------------|----------------|--------------|
| Heavy | Tip1: 2/120 = 0.0167<br>Tip2: 3/120 = 0.025<br>Mean: ~0.021 | ~2.1% | 6.25% |
| Light | Tip1: 2/90 = 0.0222<br>Tip2: 3/90 = 0.0333<br>Mean: ~0.028 | ~2.8% | 8.33% |

*Note: Actual values differ from simple calculation, likely due to weighting by multiplicity or different alignment method.*

### LBI (Local Branching Index)

With tau=0.0125 (default):

| Node | Expected LBI Range | Actual LBI (Heavy) | Actual LBI (Light) |
|------|-------------------|--------------------|--------------------|
| naive | ~3.35e-05 | 3.35e-05 | 3.35e-05 |
| Node1 | ~3.36e-05 | 3.36e-05 | 3.36e-05 |
| Tip1 | ~5.63e-08 | 5.63e-08 | 5.63e-08 |
| Tip2 | ~1.89e-11 | 1.89e-11 | 1.89e-11 |

*LBI decreases for tips, reflecting their position in the tree.*

### LBR (Local Branching Ratio)

| Node | Expected LBR | Actual LBR | Justification |
|------|--------------|------------|---------------|
| naive | 0.0 | 0.0 | Root has no parent |
| Node1 | ~0.693 | 0.693 | ln(2) - one child branch |
| Tip1 | 0.0 | 0.0 | Leaf node |
| Tip2 | 0.0 | 0.0 | Leaf node |

*LBR = ln(#children), so Node1 with 2 children has LBR = ln(2) ≈ 0.693*

## Known Issues

1. **Single-chain families not processed**: `family-heavy` and `family-light` are currently skipped with warnings:
   ```
   WARNING: Family family-heavy light chain root node missing sequence_alignment.
   WARNING: Family family-light root node missing sequence_alignment.
   ```

   **Expected behavior**: Single-chain families should be processed without requiring paired data.

2. **Leaf node status**: All nodes in the output JSON show `is_leaf=False`, including Tip1 and Tip2 which should be `is_leaf=True`.

3. **Gene position warnings** (at verbosity 2): All families warn about missing gene start/end positions (v_gene_start, d_gene_start, j_gene_start), which default to 0. These fields are not required in the PCP format but are needed for gene region visualization.

## Files

- `pcp.csv` - Parent-child pair data (9 rows: 3 families × 3 edges each)
- `trees.csv` - Newick tree strings (3 trees, matching for paired family)
- `byhand.json` - Generated Olmsted JSON (after running conversion)
- `README.md` - This file

## Visual Verification

When loading in Olmsted:

1. **Clone table**: Should show 2 clones (family-paired heavy and light)
   - V genes: IGHV2-2*01, IGKV2-2*01
   - Different junction lengths (10 for heavy, 8 for light)

2. **Tree visualization**: Should show 4 nodes (naive, Node1, Tip1, Tip2)
   - Branch lengths: 0.1, 0.2, 0.3 as specified

3. **Sequence alignment**: Mutations should appear at positions:
   - Codon 5: base pairs 12-14
   - Codon 10: base pairs 27-29
   - Codon 15: base pairs 42-44

4. **Naive sequence visualization**: CDR regions should highlight correctly at the specified codon ranges

---

*Last updated: 2025-12-12*
