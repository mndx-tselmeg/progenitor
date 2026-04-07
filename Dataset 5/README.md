### Quality check metrics

```
QC Metric Distributions:
nCount_RNA (UMI):
  Median: 7328 | Q1-Q3: 4617 - 11406 
nFeature_RNA (genes):
  Median: 2666 | Q1-Q3: 1992 - 3423 
percent.mt:
  Median: 2.75 % | 95th percentile: 15.13 %
percent.ribo:
  Median: 15.98 %
percent.hb:
  Median: 0 %
```
### Quality check violin plots

<img width="6000" height="2400" alt="image" src="https://github.com/user-attachments/assets/2535aebd-485e-471d-9718-63b6ac2ee580" />

### Quality check scatterplots

<img width="7200" height="2400" alt="image" src="https://github.com/user-attachments/assets/3cd2391b-243c-408b-a93b-9852b7945dfb" />

### Quality check pass scatterplot

<img width="600" height="600" alt="image" src="https://github.com/user-attachments/assets/7a761276-2d23-4274-954a-cb256492ad50" />

Cells remaining after QC (including mt and hb threshold): 13023 

### Initial principal component analysis

<img width="400" height="400" alt="image" src="https://github.com/user-attachments/assets/6569dc89-1390-4b68-8855-b81cf51f5843" />
<img width="400" height="400" alt="image" src="https://github.com/user-attachments/assets/7bd8c4ae-9740-4bdd-ab17-134be8ecc4aa" />
<img width="400" height="400" alt="image" src="https://github.com/user-attachments/assets/592ffc75-4249-4730-92fa-b35f57e3391a" />

### UMAP

#### Exploring batch effect

<img width="400" height="400" alt="image" src="https://github.com/user-attachments/assets/a28b7c5c-3fd5-4b1c-b8b3-5ef1ec031ab8" />

### Cell cycle scoring

<img width="400" height="400" alt="image" src="https://github.com/user-attachments/assets/b640eb31-731c-4908-b6cf-c2800ddc1c5c" />
<img width="400" height="400" alt="image" src="https://github.com/user-attachments/assets/f796f43a-1cd8-4625-98df-38c5c97afe2a" />

### Quality check visualisation on UMAP and PC plot

<img width="400" height="400" alt="image" src="https://github.com/user-attachments/assets/069648cf-cbf2-4595-8cb8-44b04db0842d" />
<img width="400" height="400" alt="image" src="https://github.com/user-attachments/assets/29bd3ff3-1c1f-4021-a55d-869df29c0804" />

### Clustering

<img width="4500" height="2400" alt="image" src="https://github.com/user-attachments/assets/e6edfb90-c2ed-48de-b180-c900cc1a682c" />

Cluster resolution 0.4 produced the best silhouette score

### Exploring known neural population markers

#### Neural progenitor markers
<img width="4500" height="1500" alt="image" src="https://github.com/user-attachments/assets/1541b308-942e-49ff-aa10-1135ec5eba66" />

#### Glial progenitor markers
<img width="4500" height="1500" alt="image" src="https://github.com/user-attachments/assets/34e0957e-3694-4ceb-b151-6dae414523d6" />

#### Glial cell markers
<img width="750" height="375" alt="image" src="https://github.com/user-attachments/assets/df199d46-1545-45bd-8c5a-846fb0b989cc" />

#### Mature neuron markers
<img width="4500" height="1500" alt="image" src="https://github.com/user-attachments/assets/cbf95c20-35e4-4e57-95af-87c663412fab" />

### Visualising the most optimal clustering results and CellTypist labelling on UMAP
<img width="3000" height="1500" alt="image" src="https://github.com/user-attachments/assets/d33cd585-0606-4920-b6f7-1c32f41d21f1" />

### Exploring cluster identity
<img width="600" height="500" alt="image" src="https://github.com/user-attachments/assets/fb16b342-f9dd-43cb-9d1e-dd2c4e9b9c3f" />

### How many cells of each type belong to each embryonic stage?
```
                   E14  E17
  Glioblast        157   52
  Immune           147   21
  Neuroblast      3272  603
  Neuron          6728 1404
  Oligodendrocyte   13    0
  Radial glia      457    0
  Vascular           1    0
```
### Differential expression analysis
<img width="2400" height="2100" alt="image" src="https://github.com/user-attachments/assets/c88a0118-608b-4d9a-827d-ca5a96e357fa" />

