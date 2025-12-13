# PCP Byhand Dataset Specification v2

## Purpose
Hand-crafted test dataset with three distinct tree topologies for comprehensive testing of tree visualization, navigation, and metric calculations.

## Design Principles
1. **100 amino acid sequences** (300 nucleotides)
2. **Regular mutation intervals**: Mutations every 10 codons
3. **Unique mutations per node**: Each node mutates to a different amino acid
4. **Three topology types**: Ladder, Balanced, Imbalanced

## Tree Topologies

### Family 1: family-heavy-ladder (Ladder tree, 5 tips)
**Type**: Heavy chain only (IGH)
**Sample**: sample-A
**Topology**: Linear ladder - each internal node has one tip and one continuing branch

```
Newick: (((((T1:0.1)N1:0.1,(T2:0.1)N2:0.1)N3:0.1,(T3:0.1)N4:0.1)N5:0.1,(T4:0.1)N6:0.1)N7:0.1,T5:0.1)naive;

Tree structure:
naive
  ├─ N1 (depth 1)
  │   ├─ T1 (leaf)
  │   └─ N2 (depth 2)
  │       ├─ T2 (leaf)
  │       └─ N3 (depth 3)
  │           ├─ T3 (leaf)
  │           └─ N4 (depth 4)
  │               ├─ T4 (leaf)
  │               └─ T5 (leaf)
```

**Node count**: 1 root + 4 internal + 5 leaves = 10 nodes total

---

### Family 2: family-light-balanced (Balanced binary tree, 7 tips)
**Type**: Light chain only (IGK - kappa)
**Sample**: sample-B
**Topology**: Balanced binary tree

```
Newick: (((T1:0.2,T2:0.2)N1:0.1,(T3:0.2,T4:0.2)N2:0.1)N3:0.1,((T5:0.2,T6:0.2)N4:0.1,T7:0.2)N5:0.1)naive;

Tree structure:
           naive
          /     \
        N3       N5
       / \       / \
      N1  N2    N4  T7
     / \ / \   / \
    T1 T2 T3 T4 T5 T6
```

**Node count**: 1 root + 5 internal + 7 leaves = 13 nodes total

---

### Family 3: family-paired-imbalanced (Imbalanced tree, 6 tips)
**Type**: Paired heavy + light chains
**Sample**: sample-C
**Topology**: Highly imbalanced - most leaves on one side

```
Newick: ((T1:0.3,((T2:0.3,(T3:0.3,(T4:0.3,(T5:0.3,T6:0.3)N4:0.1)N3:0.1)N2:0.1)N1:0.1)naive;

Tree structure:
naive
  ├─ T1 (depth 1)
  └─ N1 (depth 1)
      ├─ T2 (depth 2)
      └─ N2 (depth 2)
          ├─ T3 (depth 3)
          └─ N3 (depth 3)
              ├─ T4 (depth 4)
              └─ N4 (depth 4)
                  ├─ T5 (depth 5)
                  └─ T6 (depth 5)
```

**Node count**: 1 root + 4 internal + 6 leaves = 11 nodes total

---

## Sequence Design

### Common Parameters
- **Length**: 100 amino acids (300 nucleotides)
- **Mutation interval**: Every 10 codons
- **Germline codon**: AAA (Lysine, K) for all families
- **Mutation scheme**: Each node mutates to a unique amino acid

### Mutation Mapping
Mutations occur at codon positions: 10, 20, 30, 40, 50, 60, 70, 80, 90, 100

**Amino acid assignments by node depth from root**:
| Depth | Codon | From | To | Amino Acid | Description |
|-------|-------|------|-----|------------|-------------|
| 0 (naive) | all | AAA | - | K (Lysine) | Germline |
| 1 | 10 | AAA | GCA | A (Alanine) | First generation |
| 2 | 20 | AAA | GGA | G (Glycine) | Second generation |
| 3 | 30 | AAA | GTA | V (Valine) | Third generation |
| 4 | 40 | AAA | TTA | L (Leucine) | Fourth generation |
| 5 | 50 | AAA | ATA | I (Isoleucine) | Fifth generation |
| 6 | 60 | AAA | CCA | P (Proline) | Sixth generation |
| 7 | 70 | AAA | TTC | F (Phenylalanine) | Seventh generation |
| 8 | 80 | AAA | TAC | Y (Tyrosine) | Eighth generation |
| 9 | 90 | AAA | TGG | W (Tryptophan) | Ninth generation |
| 10 | 100 | AAA | CAC | H (Histidine) | Tenth generation |

**Accumulated mutations**: Each node inherits all mutations from its ancestors plus adds its own.

### Example Sequences

**family-heavy-ladder (heavy chain)**:
- naive: `AAA` × 100 (all K)
- N1 (depth 1): `AAA`×9 + `GCA` + `AAA`×90 (mutation at codon 10 → A)
- T1 (leaf, depth 1): Same as N1
- N2 (depth 2): `AAA`×9 + `GCA` + `AAA`×9 + `GGA` + `AAA`×80 (mutations at 10→A, 20→G)
- T2 (leaf, depth 2): Same as N2
- N3 (depth 3): Previous + mutation at codon 30 → V
- T3 (leaf, depth 3): Same as N3
- N4 (depth 4): Previous + mutation at codon 40 → L
- T4 (leaf, depth 4): Same as N4
- T5 (leaf, depth 4): Same as N4

**family-light-balanced (light chain - kappa)**:
- naive: `AAA` × 100
- N3 (depth 1), N5 (depth 1): Mutation at codon 10 → A
- N1 (depth 2), N2 (depth 2), N4 (depth 2): Mutations at 10→A, 20→G
- T1-T6 (depth 3): Mutations at 10→A, 20→G, 30→V
- T7 (depth 2): Mutations at 10→A, 20→G

**family-paired-imbalanced (paired heavy + light)**:

*Heavy chain*:
- naive: `AAA` × 100
- T1 (depth 1): Mutation at 10 → A
- N1 (depth 1): Mutation at 10 → A (same as T1, different branch)
- T2 (depth 2): Mutations at 10→A, 20→G
- N2-N4: Continue pattern
- T3-T6: Terminal leaves with accumulated mutations

*Light chain*: Same topology, same mutation pattern (identical for paired data)

---

## Gene Assignments

### family-heavy-ladder (Heavy only)
- Locus: IGH
- V gene: IGHV1-1*01
- J gene: IGHJ1*01
- CDR1: codons 20-25 (nucleotides 60-75)
- CDR2: codons 40-45 (nucleotides 120-135)
- CDR3: codons 85-95 (nucleotides 255-285)

### family-light-balanced (Light only)
- Locus: IGK
- V gene: IGKV1-1*01
- J gene: IGKJ1*01
- Light chain type: kappa
- CDR1: codons 15-20 (nucleotides 45-60)
- CDR2: codons 35-40 (nucleotides 105-120)
- CDR3: codons 80-90 (nucleotides 240-270)

### family-paired-imbalanced (Paired)

**Heavy chain**:
- Locus: IGH
- V gene: IGHV2-2*01
- J gene: IGHJ2*01
- CDR1: codons 20-25 (nucleotides 60-75)
- CDR2: codons 40-45 (nucleotides 120-135)
- CDR3: codons 85-95 (nucleotides 255-285)

**Light chain**:
- Locus: IGK
- V gene: IGKV2-2*01
- J gene: IGKJ2*01
- Light chain type: kappa
- CDR1: codons 15-20 (nucleotides 45-60)
- CDR2: codons 35-40 (nucleotides 105-120)
- CDR3: codons 80-90 (nucleotides 240-270)

---

## Branch Lengths

All branch lengths are set to make visualization clear:
- Short branches: 0.1 (internal nodes close to parent)
- Medium branches: 0.2 (most leaf branches)
- Long branches: 0.3 (some leaves farther from parent)

---

## Expected Olmsted Output

### Datasets
- 1 dataset

### Clones
- family-heavy-ladder: 1 clone (heavy)
- family-light-balanced: 1 clone (light)
- family-paired-imbalanced: 2 clones (heavy + light, paired)
**Total**: 4 clones

### Trees
- family-heavy-ladder: 1 tree (10 nodes)
- family-light-balanced: 1 tree (13 nodes)
- family-paired-imbalanced: 2 trees (11 nodes each, heavy + light)
**Total**: 4 trees

---

## Validation Criteria

1. ✅ All sequences are 300 nucleotides (100 codons)
2. ✅ Germline is all AAA codons
3. ✅ Mutations occur at specified intervals (every 10 codons)
4. ✅ Each node has unique mutation(s) based on depth
5. ✅ Leaf nodes correctly marked as type="leaf"
6. ✅ Internal nodes correctly marked as type="internal"
7. ✅ Root correctly marked as type="root"
8. ✅ Parent references correctly set for all non-root nodes
9. ✅ Tree topologies match specifications (ladder, balanced, imbalanced)
10. ✅ CDR regions correctly specified and non-overlapping

---

*Last updated: 2025-12-12*
