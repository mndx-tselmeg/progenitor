# Single-Cell RNA-seq analysis of Dataset 2

## Information About the Dataset

This dataset is adapted from a paper that was published in 2018 titled ["Developmental diversification of cortical inhibitory interneurons"](https://pubmed.ncbi.nlm.nih.gov/29513653/).

ScRNA-seq data is published at GEO with accession [GSE103983](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103983).

This dataset contains scRNA-seq (Drop-seq) of MGE, CGE and LGE of E13.5 (MGE) and E14.5 (CGE, LGE) mouse embryos. 

## Exploring the Dataset 

This dataset contains 21566 cells (barcodes) and 19184 genes (features). 

## Quality Check

### Quality Check metrics summary

```
QC Metric Distributions:
nCount_RNA (UMI):
  Median: 1283 | Q1-Q3: 1025 - 1740 
nFeature_RNA (genes):
  Median: 927 | Q1-Q3: 771 - 1180 
percent.mt:
  Median: 1.35 % | 95th percentile: 2.62 %
percent.ribo:
  Median: 3.21 %
percent.hb:
  Median: 0.07 %
```

### Quality Check Violin Plots

<img width="6000" height="2400" alt="image" src="https://github.com/user-attachments/assets/74ec6d7a-e53f-44db-8100-24542026724f" />

Looking at the total nCountRNA and nFeature_RNA, mitochondrial, ribosomal, and hemoglobin percentage, it seems that the dataset already went through basic quality check filtration.

### Quality Check Scatterplots

<img width="7200" height="2400" alt="image" src="https://github.com/user-attachments/assets/17415cf9-9a09-4519-858e-3c4c4f521dd0" />

The fact that data has been filtered can also be seen from the scatterplots. Another thing that can be seen from this scatterplots is that there is no big difference in quality based on brain regions. 

### Checking for empty droplets

<img width="2400" height="2400" alt="image" src="https://github.com/user-attachments/assets/e82e299f-f629-48ba-b459-2bb5b6ee3107" />

It looks like there is no infliction point in UMI counts and there is a sharp drop in UMI counts vs. barcode rank in this dataset. It means that empty droplets have already been removed.


