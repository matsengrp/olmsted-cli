# PCP Eyetest Dataset

Minimal test dataset for visual verification of Olmsted visualizations.

## Dataset Structure

- **2 clonal families** (family-1, family-2)
- **5 tips per family** (Tip1-Tip5)
- **150bp sequences** (50 codons)

## What to Look For

### Clone Table / Scatterplot

| Property | family-1 | family-2 |
|----------|----------|----------|
| V gene | IGHV1-1*01 | IGHV2-2*01 |
| J gene | IGHJ1*01 | IGHJ2*01 |
| Junction length | 10 | 20 |
| Sample | sample-A | sample-B |

### Tree Visualization

Both families have the same topology but different branch lengths:

```
naive
  └─ Node1
       ├─ Node2
       │    ├─ Tip1
       │    └─ Tip2
       └─ Node3
            ├─ Node4
            │    ├─ Tip3
            │    └─ Tip4
            └─ Tip5
```

**Branch lengths to verify** (all values 0.1-0.5 in 0.1 increments):
- Tip3, Tip4 have longest branches (0.4, 0.5)
- Node4 has shortest internal branch (0.1)

### Sequence Alignment

**family-1 germline**: `ATGCAG` repeating (codes for Met-Gln)
**family-2 germline**: `CAGTTG` repeating (codes for Gln-Leu)

Mutations are single codon changes (`ATG`→`TTG` or `CAG`→`GAT`):
- Each tip has 2-3 mutations from germline
- Mutations accumulate along branches (verify parent sequences contain subset of child mutations)

### CDR Regions (Naive Sequence Viz)

| Region | family-1 | family-2 |
|--------|----------|----------|
| CDR1 | codons 4-8 | codons 5-12 |
| CDR2 | codons 12-16 | codons 18-25 |
| CDR3 | codons 35-45 | codons 30-50 |

Verify CDR highlighting positions match these ranges.

## Generating Output

```bash
cd olmsted-cli
python -c "from olmsted_cli import process_pcp_data; import sys; sys.argv = ['', '-i', 'example_data/pcp-eyetest/pcp.csv', '-t', 'example_data/pcp-eyetest/trees.csv', '-o', 'example_data/pcp-eyetest/output.json', '--compute-metrics']; process_pcp_data.main()"
```

## Files

- `pcp.csv` - Parent-child pair data
- `trees.csv` - Newick trees
- `output.json` - Generated Olmsted JSON (after running above command)
