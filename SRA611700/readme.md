# Single-Cell RNA-seq Analysis of SRS2532367

**Data Source:**  
[https://panglaodb.se/view_data.php?sra=SRA611700&srs=SRS2532367](https://panglaodb.se/view_data.php?sra=SRA611700&srs=SRS2532367)

**Library Preparation Protocol:**  
10x Genomics Chromium single-cell library preparation.

**Sample Description:**  
This dataset originates from the **Lhx6(BAC)-GFP transgenic mouse caudal ganglionic eminence (CGE)** at **embryonic day 13.5 (E13.5)**.

---

## Reference Data for Cell Type Labelling

- **Marker list (Le Manno et al., Developing Mouse Brain Cell Diversity Database):**  
  [Google Sheet](https://docs.google.com/spreadsheets/d/1-V6O3xqCEBhWjMVLcCCeNx7_edIf3_1SkqS-gZ5FQzE/edit?usp=sharing)
  
- **ScType Database:**  
  [ScTypeDB_full.xlsx](https://github.com/IanevskiAleksandr/sc-type/blob/master/ScTypeDB_full.xlsx)

---

## Initial Data Overview

<img width="597" height="222" alt="Screenshot 2025-10-01 at 15 58 20" src="https://github.com/user-attachments/assets/f8f3a18c-8f15-4651-aa82-1a154d949c39" />

- **Total barcodes:** 8,876  
- **Total features:** 29,198  

---

## Ambient RNA Detection with DecontX

DecontX was run on the raw expression matrix using default parameters.  
It identified **4 distinct clusters** based on ambient RNA profiles.

<img width="824" height="536" alt="DecontX clusters" src="https://github.com/user-attachments/assets/e7c3d2a6-0ae4-419a-bdfe-74d80a01c38b" />

**Ambient RNA contamination levels within clusters:**

<img width="824" height="536" alt="Ambient RNA contamination" src="https://github.com/user-attachments/assets/0e818314-e3fc-472e-a314-7dd9cc068df5" />

---

## Quality Control

Detected genes:
- **Mitochondrial:** 28  
- **Hemoglobin:** 11  
- **Ribosomal:** 321  

<img width="2796" height="1068" alt="QC plot" src="https://github.com/user-attachments/assets/dc72edb6-edd0-486e-95fc-bd96af777fe2" />

The plot above shows contamination distribution relative to the number of detected features.

**Filtering thresholds (based on literature and Seurat recommendations):**
- **Mitochondrial genes:** < 5% (recommended threshold for viable cells [[PMC8599307]](https://pmc.ncbi.nlm.nih.gov/articles/PMC8599307/))  
- **Hemoglobin genes:** < 1% (expected in non-blood tissue)  
- **Ribosomal genes:** not filtered; instead, regressed out during scaling (per [Satija Lab’s Seurat tutorial](https://satijalab.org/seurat/articles/sctransform_vignette.html))

> *Note:* Ribosomal gene expression can be high in metabolically active or proliferating cells. Filtering them out may bias the dataset.

**Post-filtering dataset:**  
- **Features:** 24,069  
- **Cells:** 7,951  

---

## Doublet Detection with scDblFinder

<img width="2324" height="1070" alt="scDblFinder results" src="https://github.com/user-attachments/assets/6047346b-a9f4-4a01-b421-ff2b4cf846b0" />

scDblFinder was run with default parameters.  
It detected **633 doublets** (~8%), which is within the expected range for 10x Genomics data (~1% per 1,000 cells).

---

## Normalization and Clustering

<img width="724" height="536" alt="PC1-PC2 plot" src="https://github.com/user-attachments/assets/cd75b28a-d3bf-46bc-8caf-abcc1e97da62" />

**Principal Component Analysis (PCA):**  
The PC1–PC2 plot indicates no major batch effects or unexpected separations.

<img width="824" height="536" alt="Elbow plot 1" src="https://github.com/user-attachments/assets/3436fd5b-d7fa-43da-94b1-d748e21f26e1" />  
<img width="1404" height="1282" alt="Elbow plot 2" src="https://github.com/user-attachments/assets/f8afdde6-4ffd-4d1f-a0b4-44312af4b677" />  
<img width="824" height="536" alt="Elbow plot 3" src="https://github.com/user-attachments/assets/0aa722cb-d589-43f5-9398-c2ab4aa0c013" />

Approximately the first **25 principal components** capture most of the biological variability.

<img width="1326" height="819" alt="UMAP clusters" src="https://github.com/user-attachments/assets/b7c80f56-5550-4c56-9f1c-1017f53ddc2a" />

Using a **resolution of 1**, **18 clusters** were identified.

<img width="1326" height="819" alt="Cluster UMAP" src="https://github.com/user-attachments/assets/f3dccf50-d69a-43eb-91b8-ed698c05570e" />

---

## Cell Type Labelling

### Using ScType

<img width="1908" height="772" alt="ScType heatmap" src="https://github.com/user-attachments/assets/662b6cb5-8d98-4b42-8e12-46cf3a88dea3" />  
<img width="2222" height="1436" alt="ScType UMAP" src="https://github.com/user-attachments/assets/7434c824-8865-42f8-9710-efc22d890c0e" />

| Cell Type | Count |
|------------|-------|
| Mature neurons | 1,046 |
| Neural stem cells | 632 |
| Neuroblast – Forebrain GABAergic | 2,986 |
| Neuroblast – Forebrain glutamatergic | 1,033 |
| Neuroblast – Midbrain glutamatergic | 167 |
| Neuroblast – Mixed region | 79 |
| Neuroblast – Spinal cord glutamatergic | 93 |
| Non-myelinating Schwann cells | 344 |
| Radial glial cells | 685 |
| Schwann precursor cells | 253 |

---

