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


