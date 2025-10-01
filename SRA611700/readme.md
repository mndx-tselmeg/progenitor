Data used for this analysis is available here: https://panglaodb.se/view_data.php?sra=SRA611700&srs=SRS2532369 
10x chromium protocol was used for this single cell library preparation.
This data comes from Lhx6(BAC)-GFP transgenic mouse caudal ganglionic eminence at developmental stage E13.5. 
# Looking at the data 
<img width="597" height="222" alt="Screenshot 2025-10-01 at 15 58 20" src="https://github.com/user-attachments/assets/f8f3a18c-8f15-4651-aa82-1a154d949c39" />

This data contains 8876 barcodes and 29198 features initially. 
# Running DecontX to detect the ambient RNA
DecontX program was ran on raw expression matrix directly with its standard attributes.
DecontX identified 4 distinct clusters.
<img width="824" height="536" alt="image" src="https://github.com/user-attachments/assets/e7c3d2a6-0ae4-419a-bdfe-74d80a01c38b" />

Here we can see ambient RNA contamination levels of the cells on the same cluster.
<img width="824" height="536" alt="image" src="https://github.com/user-attachments/assets/0e818314-e3fc-472e-a314-7dd9cc068df5" />

# Quality checking process
28 mitochondrial genes, 11 hemoglobin genes, and 321 ribosomal genes were detected in the sample. 
<img width="2796" height="1068" alt="image" src="https://github.com/user-attachments/assets/dc72edb6-edd0-486e-95fc-bd96af777fe2" />
Contamination distribution relative to the number of features are shown here. 

According to literature in mouse scRNA-seq library 5% mitochondrial gene contamination threshold performs well to distinguish between healthy and low-quality cells https://pmc.ncbi.nlm.nih.gov/articles/PMC8599307/) and almost no hemoglobin gene if it is not from blood sample. So I used lower than 5% for mitochondrial gene and lower than 1% hemoglobin gene contamination as threshold for quality check. Satija Lab’s Seurat tutorials (https://satijalab.org/seurat/articles/sctransform_vignette.html) recommend regressing out ribosomal gene effects during scaling rather than discarding cells since ribosomal gene expression is genuinely high in certain physiological states (e.g. proliferating cells). Removing such cells could bias the dataset against highly active or metabolically busy cell types.

After filtering according to the thresholds I'm left with 24069 features and 7951 cells. 
# Running scDblFinder to detect doublets
<img width="2324" height="1070" alt="image" src="https://github.com/user-attachments/assets/6047346b-a9f4-4a01-b421-ff2b4cf846b0" />

scDblFinder was ran with deafult settings and found 633 (8% doublet rate) doublets. This number seems plausible since 10X samples contains 1% doublet per 1000 cells. 
# Normalisation and clustering
<img width="724" height="536" alt="image" src="https://github.com/user-attachments/assets/cd75b28a-d3bf-46bc-8caf-abcc1e97da62" />

From this PC1 to PC2 plot it can be seen that there is no weird separation/batch effect going on.
<img width="824" height="536" alt="image" src="https://github.com/user-attachments/assets/3436fd5b-d7fa-43da-94b1-d748e21f26e1" />
<img width="1404" height="1282" alt="image" src="https://github.com/user-attachments/assets/f8afdde6-4ffd-4d1f-a0b4-44312af4b677" />
<img width="824" height="536" alt="image" src="https://github.com/user-attachments/assets/0aa722cb-d589-43f5-9398-c2ab4aa0c013" />

It looks like using first ~25 (just to be sure) captures enough variability in the data. 
<img width="1326" height="819" alt="image" src="https://github.com/user-attachments/assets/b7c80f56-5550-4c56-9f1c-1017f53ddc2a" />
Using resolution of 1 results in 18 clusters. 

<img width="1326" height="819" alt="image" src="https://github.com/user-attachments/assets/f3dccf50-d69a-43eb-91b8-ed698c05570e" />


