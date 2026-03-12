# Single-Cell RNA-seq Analysis of Dataset 2

## Information About the Dataset

- This dataset is adapted from a 2018 study titled  
  ["Developmental diversification of cortical inhibitory interneurons"](https://pubmed.ncbi.nlm.nih.gov/29513653/).

- The scRNA-seq data are publicly available through GEO under accession number  
  [GSE103983](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103983).

- The dataset contains **Drop-seq single-cell RNA sequencing data** from three embryonic mouse brain regions:
  - **Medial ganglionic eminence (MGE)** at **E13.5**
  - **Caudal ganglionic eminence (CGE)** at **E14.5**
  - **Lateral ganglionic eminence (LGE)** at **E14.5**

---

## Exploring the Dataset

This dataset contains:

- **21,566 cells** (barcodes)
- **19,184 genes** (features)

---

## Quality Control

### Quality Control Metrics Summary

```
nCount_RNA (UMI):
Median: 1283 | Q1–Q3: 1025–1740

nFeature_RNA (genes):
Median: 927 | Q1–Q3: 771–1180

percent.mt:
Median: 1.35% | 95th percentile: 2.62%

percent.ribo:
Median: 3.21%

percent.hb:
Median: 0.07%
```

---

### Quality Control Violin Plots

<img width="6000" height="2400" alt="image" src="https://github.com/user-attachments/assets/74ec6d7a-e53f-44db-8100-24542026724f" />

Inspection of **nCount_RNA**, **nFeature_RNA**, and the percentages of **mitochondrial**, **ribosomal**, and **hemoglobin** transcripts suggests that the dataset has already undergone basic quality-control filtering.

---

### Quality Control Scatter Plots

<img width="7200" height="2400" alt="image" src="https://github.com/user-attachments/assets/17415cf9-9a09-4519-858e-3c4c4f521dd0" />

Further evidence of prior filtering is visible in the scatter plots. Additionally, no major differences in data quality are observed across the three brain regions.

---

### Checking for Empty Droplets

<img width="400" height="400" alt="image" src="https://github.com/user-attachments/assets/e82e299f-f629-48ba-b459-2bb5b6ee3107" />

The barcode rank plot shows a **sharp drop in UMI counts without a clear inflection point**, suggesting that **empty droplets have already been removed** from the dataset.

---

## Initial Dimensionality Reduction for QC

### Principal Component Analysis (PCA)

<img width="4800" height="2400" alt="image" src="https://github.com/user-attachments/assets/70726691-902a-421f-91b7-3c0c8bb89f56" />

The PCA plot shows that the three brain regions **do not strongly separate along the principal components**.

The elbow (knee) point of the PCA appears around **15 principal components**. However, due to the nature of single-cell RNA-seq data, **PCs 1–20 are used for downstream analysis**.

---

### UMAP Visualization

<img width="4800" height="2400" alt="image" src="https://github.com/user-attachments/assets/c0ea3b3d-471a-4f6f-8c51-e1c134a84c22" />

The UMAP visualization shows that:

- Cells expressing **hemoglobin genes** are distributed across the embedding, suggesting the presence of **ambient RNA contamination**.
- **MGE samples separate from CGE and LGE samples** in UMAP space, even though this separation was not apparent in the PCA plots.

---

### Ambient RNA Correction

Because **empty droplets are not present in the dataset**, methods such as **SoupX** (which rely on empty droplets) are not suitable.

Although some ambient RNA contamination is detected, it does not appear to be severe enough to justify the use of computationally intensive methods such as **CellBender**.

Instead, **DecontX** was applied to estimate and correct for ambient RNA contamination.

<img width="3000" height="1500" alt="image" src="https://github.com/user-attachments/assets/a73ccaba-5ef0-4506-8565-abcb1ab01e29" />

<img width="400" height="400" alt="image" src="https://github.com/user-attachments/assets/bbef3913-e1e8-46fa-a0a7-25a97b82c2e5" />

---

## Cell Cycle Scoring

Cell cycle scoring was performed using **Seurat’s `CellCycleScoring()` function**, based on a previously published list of cell cycle–associated genes:

https://github.com/hbc/tinyatlas/tree/master

This gene set includes:

- **43 S-phase genes**
- **54 G2/M-phase genes**

---

### Exploration of Cell Cycle Effects

<img width="400" height="400" alt="image" src="https://github.com/user-attachments/assets/49a0fd50-7f5b-4125-8428-e21f9ffc8423" />

PCA was performed using **S-phase** and **G2/M-phase gene sets as features**. The resulting plot shows **clear separation of cell cycle states along the first two principal components**.

Fully regressing out the **S score** and **G2M score** would remove all cell cycle–associated signals. However, this can negatively affect downstream analyses, particularly in developmental systems where:

- **Stem cells** may be quiescent
- **Progenitor or differentiated cells** may be proliferative

Removing all cell cycle signal may therefore **obscure biologically meaningful differences between cell populations**.

As recommended in the Seurat vignette, an alternative approach is used: **regressing the difference between G2M and S scores**.

This approach:

- **Preserves the distinction between cycling and non-cycling cells**
- **Removes differences between specific cell cycle phases among proliferating cells**

---

### Removal of Cell Cycle Effects

<img width="400" height="400" alt="image" src="https://github.com/user-attachments/assets/01e777ec-8a65-4494-b262-cd372ef90eb6" />

After regressing the **G2M–S score difference**, the PCA plot shows that **S-phase and G2/M-phase signals no longer separate along the first two principal components**, while **G1 cells remain unchanged**.

This indicates that **cell cycle phase–specific variation has been reduced without removing the broader signal distinguishing cycling from non-cycling cells**.
