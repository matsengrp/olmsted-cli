# PCP Byhand Dataset Specification v2

## Purpose
Hand-crafted test dataset with three distinct tree topologies for comprehensive testing of tree visualization, navigation, and metric calculations.

## Design Principles
1. **100 amino acid sequences** (300 nucleotides)
2. **Regular mutation intervals**: Mutations every 10 codons
3. **Unique mutations per node**: Each node mutates to a different amino acid
4. **Three topology types**: Balanced (heavy), Balanced (light), Imbalanced (paired)

## Tree Topologies

### Family 1: family-heavy (Balanced tree, 5 tips)
**Type**: Heavy chain only (IGH)
**Sample**: sample-A
**Topology**: Relatively balanced binary tree

```
Newick: (((T1:0.2,T2:0.2)N2:0.1,(T3:0.2,T4:0.2)N3:0.1)N1:0.1,T5:0.2)naive;

Tree structure:
        naive
        /   \
       N1    T5
      / \
     N2  N3
    / \  / \
   T1 T2 T3 T4
```

**Node count**: 1 root + 3 internal + 5 leaves = 9 nodes total

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

**Heavy chain mutations** (family-heavy and family-paired-imbalanced heavy):
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

**Light chain mutations** (family-light-balanced and family-paired-imbalanced light):
| Depth | Codon | From | To | Amino Acid | Description |
|-------|-------|------|-----|------------|-------------|
| 0 (naive) | all | AAA | - | K (Lysine) | Germline |
| 1 | 10 | AAA | TCA | S (Serine) | First generation |
| 2 | 20 | AAA | ACA | T (Threonine) | Second generation |
| 3 | 30 | AAA | CGA | R (Arginine) | Third generation |
| 4 | 40 | AAA | AAC | N (Asparagine) | Fourth generation |
| 5 | 50 | AAA | GAC | D (Aspartic acid) | Fifth generation |
| 6 | 60 | AAA | CAA | Q (Glutamine) | Sixth generation |
| 7 | 70 | AAA | GAA | E (Glutamic acid) | Seventh generation |
| 8 | 80 | AAA | ATG | M (Methionine) | Eighth generation |
| 9 | 90 | AAA | TGC | C (Cysteine) | Ninth generation |
| 10 | 100 | AAA | GCC | A (Alanine) | Tenth generation |

**Accumulated mutations**: Each node inherits all mutations from its ancestors plus adds its own.

### Example Sequences

**family-heavy (heavy chain)**:
- naive: `AAA` × 100 (all K)
- N1 (depth 1): `AAA`×9 + `GCA` + `AAA`×90 (mutation at codon 10 → A)
- T5 (leaf, depth 1): Same as N1
- N2 (depth 2), N3 (depth 2): `AAA`×9 + `GCA` + `AAA`×9 + `GGA` + `AAA`×80 (mutations at 10→A, 20→G)
- T1, T2, T3, T4 (leaves, depth 3): Previous + mutation at codon 30 → V

**family-light-balanced (light chain - kappa)**:
- naive: `AAA` × 100
- N3 (depth 1), N5 (depth 1): Mutation at codon 10 → S
- N1 (depth 2), N2 (depth 2), N4 (depth 2): Mutations at 10→S, 20→T
- T1-T6 (depth 3): Mutations at 10→S, 20→T, 30→R
- T7 (depth 2): Mutations at 10→S, 20→T

**family-paired-imbalanced (paired heavy + light)**:

*Heavy chain*:
- naive: `AAA` × 100
- T1 (depth 1): Mutation at 10 → A (Alanine)
- N1 (depth 1): Mutation at 10 → A (same as T1, different branch)
- T2 (depth 2): Mutations at 10→A, 20→G
- N2-N4: Continue pattern (depth 3: +V, depth 4: +L, depth 5: +I)
- T3-T6: Terminal leaves with accumulated mutations

*Light chain*:
- naive: `AAA` × 100
- T1 (depth 1): Mutation at 10 → S (Serine, **different from heavy!**)
- N1 (depth 1): Mutation at 10 → S
- T2 (depth 2): Mutations at 10→S, 20→T (**different from heavy!**)
- N2-N4: Continue pattern (depth 3: +R, depth 4: +N, depth 5: +D)
- T3-T6: Terminal leaves with accumulated mutations (**different amino acids than heavy**)

---

## Gene Assignments

### family-heavy (Heavy only)
- Locus: IGH
- V gene: IGHV1-1*01
- J gene: IGHJ1*01
- CDR1: codons 20-29 (nucleotides 58-87, 10 AA = 30 nt)
- CDR2: codons 50-59 (nucleotides 148-177, 10 AA = 30 nt)
- CDR3: codons 80-89 (nucleotides 238-267, 10 AA = 30 nt)

### family-light-balanced (Light only)
- Locus: IGK
- V gene: IGKV1-1*01
- J gene: IGKJ1*01
- Light chain type: kappa
- CDR1: codons 10-19 (nucleotides 28-57, 10 AA = 30 nt)
- CDR2: codons 40-49 (nucleotides 118-147, 10 AA = 30 nt)
- CDR3: codons 70-79 (nucleotides 208-237, 10 AA = 30 nt)

### family-paired-imbalanced (Paired)

**Heavy chain**:
- Locus: IGH
- V gene: IGHV2-2*01
- J gene: IGHJ2*01
- CDR1: codons 20-29 (nucleotides 58-87, 10 AA = 30 nt)
- CDR2: codons 50-59 (nucleotides 148-177, 10 AA = 30 nt)
- CDR3: codons 80-89 (nucleotides 238-267, 10 AA = 30 nt)

**Light chain**:
- Locus: IGK
- V gene: IGKV2-2*01
- J gene: IGKJ2*01
- Light chain type: kappa
- CDR1: codons 10-19 (nucleotides 28-57, 10 AA = 30 nt)
- CDR2: codons 40-49 (nucleotides 118-147, 10 AA = 30 nt)
- CDR3: codons 70-79 (nucleotides 208-237, 10 AA = 30 nt)

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
- family-heavy: 1 clone (heavy)
- family-light-balanced: 1 clone (light)
- family-paired-imbalanced: 2 clones (heavy + light, paired)
**Total**: 4 clones

### Trees
- family-heavy: 1 tree (9 nodes: 1 root + 3 internal + 5 leaves)
- family-light-balanced: 1 tree (13 nodes: 1 root + 5 internal + 7 leaves)
- family-paired-imbalanced: 2 trees (11 nodes each: 1 root + 4 internal + 6 leaves per tree)
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
9. ✅ Tree topologies match specifications (balanced, balanced, imbalanced)
10. ✅ CDR regions correctly specified and non-overlapping

---

*Last updated: 2025-12-13*
