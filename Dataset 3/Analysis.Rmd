---
title: "Untitled"
author: "Mandakh Tselmeg"
date: "2026-02-26"
output: html_document
---

### Load the data in the environment

```{r}
set.seed(123)
library(dplyr)
library(Seurat)
library(patchwork)

setwd("~/personal/15140458/project/Dataset3")
counts <- read.table(
  "GSE199352_Gene_UMI_Matrix_10x.txt",
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE
)
```

```{r}
gene_symbols <- sub(".*-", "", rownames(counts))
counts <- rowsum(counts, group = gene_symbols)
Dataset3 <- CreateSeuratObject(counts = counts)
```

### Creating a csv file for CellTypist

```{r}
counts <- t(counts)
write.csv(counts, "GSE199352.csv", row.names = TRUE)
```

### Quality check metrics

```{r}
Dataset3[["percent.mt"]] <- PercentageFeatureSet(Dataset3, pattern = "^mt-")
Dataset3[["percent.ribo"]] <- PercentageFeatureSet(Dataset3, pattern = "^Rpl|^Rps")
Dataset3[["percent.hb"]] <- PercentageFeatureSet(Dataset3, pattern = "^Hba|^Hbb")

# Gene complexity metric
Dataset3$log10GenesPerUMI <- log10(Dataset3$nFeature_RNA) / log10(Dataset3$nCount_RNA)
```

### Quality check metrics summary

```{r}
cat("\nQC Metric Distributions:\n")

cat("nCount_RNA (UMI):\n")
cat("  Median:", median(Dataset3$nCount_RNA),
    "| Q1-Q3:", quantile(Dataset3$nCount_RNA,0.25), "-",
    quantile(Dataset3$nCount_RNA,0.75), "\n")

cat("nFeature_RNA (genes):\n")
cat("  Median:", median(Dataset3$nFeature_RNA),
    "| Q1-Q3:", quantile(Dataset3$nFeature_RNA,0.25), "-",
    quantile(Dataset3$nFeature_RNA,0.75), "\n")

cat("percent.mt:\n")
cat("  Median:", round(median(Dataset3$percent.mt),2), "%",
    "| 95th percentile:", round(quantile(Dataset3$percent.mt,0.95),2), "%\n")

cat("percent.ribo:\n")
cat("  Median:", round(median(Dataset3$percent.ribo),2), "%\n")

cat("percent.hb:\n")
cat("  Median:", round(median(Dataset3$percent.hb),2), "%\n")
```

### Quality check violin plots

```{r}
p_qc <- VlnPlot(
  Dataset3,
  features = c("nCount_RNA",
               "nFeature_RNA",
               "percent.mt",
               "percent.ribo",
               "percent.hb"),
  ncol = 5,
  pt.size = 0.1
)
ggsave("Dataset3_QC_plots.png", p_qc, width = 20, height = 8, dpi = 300)
```

### Quality check scatter plots

```{r}
p1 <- FeatureScatter(Dataset3, feature1="nCount_RNA", feature2="nFeature_RNA")
p2 <- FeatureScatter(Dataset3, feature1="nCount_RNA", feature2="percent.mt")
p3 <- FeatureScatter(Dataset3, feature1="percent.mt", feature2="percent.ribo")

p_scatter <- p1 + p2 + p3

ggsave("Dataset3_QC_scatterplots.png", p_scatter, width = 24, height = 8, dpi = 300)
```

### Quality check filtration

```{r}
nfeature_min <- 300
nfeature_max <- 6000
ncount_min <- 1000
ncount_max <- 30000
mt_thresh <- 10
hb_thresh <- 1

qc_df <- data.frame(
  nCount_RNA = Dataset3$nCount_RNA,
  nFeature_RNA = Dataset3$nFeature_RNA
)

qc_df$pass_qc <- with(qc_df,
  nCount_RNA >= ncount_min &
  nCount_RNA <= ncount_max &
  nFeature_RNA >= nfeature_min &
  nFeature_RNA <= nfeature_max
)

p4 <- ggplot(qc_df, aes(x = log10(nCount_RNA + 1), y = log10(nFeature_RNA + 1), color = pass_qc)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_vline(xintercept = log10(c(ncount_min, ncount_max)), linetype = "dashed", color = "red") +
  geom_hline(yintercept = log10(c(nfeature_min, nfeature_max)), linetype = "dashed", color = "red") +
  scale_color_manual(values = c("TRUE" = "#06D6A0", "FALSE" = "#EF476F")) +
  labs(
    title = "Cell Filtering Thresholds",
    subtitle = paste0(sum(qc_df$pass_qc), " cells pass QC (", round(sum(qc_df$pass_qc)/nrow(qc_df)*100, 1), "%)"),
    x = "log10(UMI + 1)",
    y = "log10(Genes + 1)",
    color = "Pass QC"
  ) +
  theme_classic()

ggsave("Dataset3_QC_threshold.png", p4, width = 8, height = 8, dpi = 300)
```

```{r}
Dataset3 <- subset(
  Dataset3,
  subset =
    nCount_RNA >= ncount_min &
    nCount_RNA <= ncount_max &
    nFeature_RNA >= nfeature_min &
    nFeature_RNA <= nfeature_max &
    percent.mt < mt_thresh &
    percent.hb < hb_thresh
)

cat("Cells remaining after QC:", ncol(Dataset3), "\n")
```

```{r}
FeaturePlot(Dataset3,
            reduction = "pca",
            features = "percent.ribo")
```

### Initial normalisation steps

```{r}
Dataset3 <- SCTransform(Dataset3,
                   method="glmGamPoi")
```

```{r}
Dataset3 <- RunPCA(Dataset3)
PCAplot <- ElbowPlot(Dataset3)
ggsave("Dataset3_PCAelbow.png", PCAplot, width = 8, height = 8, dpi = 300)
```

```{r}
PCAplot2 <- DimPlot(Dataset3, reduction = "pca")
ggsave("Dataset3_PCA.png", PCAplot2, width = 8, height = 8, dpi = 300)
```

### Cell cycle scoring

```{r}
# Score cells for cell cycle
Dataset3 <- CellCycleScoring(Dataset3,
                           g2m.features = g2m_genes, 
                           s.features = s_genes)
```

### Inspection of cell cycle effects

```{r}
Dataset3_cc <- RunPCA(Dataset3, features=c(s_genes,g2m_genes))
PCAcc <- DimPlot(Dataset3_cc, reduction="pca", group.by="Phase")
ggsave("Dataset3_PCAcc.png", PCAcc, width = 8, height = 8, dpi = 300)
```

### Doing SCTransform again with S and G2 phase difference regression

```{r}
Dataset3$CC.Difference <- Dataset3$S.Score - Dataset3$G2M.Score
Dataset3 <- SCTransform(Dataset3, vars.to.regress = "CC.Difference", method="glmGamPoi")
```

### Checking the regression

```{r}
Dataset3_cc2 <- RunPCA(Dataset3, features=c(s_genes,g2m_genes))
PCAcc2 <- DimPlot(Dataset3_cc2, reduction="pca", group.by="Phase")
ggsave("Dataset3_PCAcc2.png", PCAcc2, width = 8, height = 8, dpi = 300)
```

### Doublet detection steps

```{r}
Dataset3 <- RunTSNE(Dataset3, dims=1:15)
library(DoubletFinder)
nExp_poi <- round(0.04 * ncol(Dataset3))  # assuming 4% expected doublets
# Parameter sweep to find optimal pK
sweep.res.list <- paramSweep(Dataset3, PCs = 1:15, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
```

```{r}
# Visualize pK metric
pKmetric_plot <- ggplot(bcmvn, aes(x = pK, y = BCmetric)) + geom_line() + geom_point()
ggsave("Dataset3_pK.png", pKmetric_plot, width = 20, height = 8, dpi = 300)
```

```{r}
# Run DoubletFinder
Dataset3 <- doubletFinder(
  Dataset3,
  PCs = 1:15,
  pN = 0.25,
  pK = 0.26,     
  nExp = nExp_poi,
  sct = TRUE
)
```

```{r}
Dataset3 <- RunUMAP(Dataset3, dims = 1:15)
DoubletPlot <- DimPlot(Dataset3, group.by = "DF.classifications_0.25_0.26_160")
ggsave("Dataset3_doublet.png", DoubletPlot, width = 8, height = 8, dpi = 300)
```

```{r}
DoubletScorePlot <- FeaturePlot(Dataset3,
                                features = "pANN_0.25_0.26_160", 
                                reduction = "umap")
ggsave("Dataset3_doubletscore.png", DoubletScorePlot, width = 8, height = 8, dpi = 300)
```

```{r}
df_col <- "DF.classifications_0.25_0.26_160"
Dataset3 <- Dataset3[, Dataset3@meta.data[[df_col]] == "Singlet"]
```

```{r}
Dataset3 <- FindNeighbors(Dataset3, dims = 1:15)
Dataset3 <- FindClusters(Dataset3, 
                    resolution = c(0.2, 0.4, 0.6, 0.8, 1, 1.2),
                    algorithm = 4,
                    random.seed = 1)
```

```{r}
resolutions <- c(0.2, 0.4, 0.6, 0.8, 1, 1.2)

plot_list <- lapply(resolutions, function(res) {
  
  cluster_col <- paste0("SCT_snn_res.", res)
  n_clusters <- length(unique(Dataset3@meta.data[[cluster_col]]))
  
  DimPlot(
    Dataset3,
    reduction = "umap",
    group.by = cluster_col,
    label = TRUE,
    label.size = 4,
    pt.size = 0.3
  ) +
    ggtitle(paste0("Resolution ", res, " (", n_clusters, " clusters)")) +
    NoLegend()
})

combined_resolutions <- wrap_plots(plot_list, ncol = 3)
ggsave("Dataset3_allcluster.png", combined_resolutions, width = 15, height = 8, dpi = 300)
```

```{r}
library(cluster)
resolutions <- c(0.2, 0.4, 0.6, 0.8, 1, 1.2)
n_cells <- ncol(Dataset3)
max_cells_for_silhouette <- 5000

if (n_cells > max_cells_for_silhouette) {
  set.seed(42)
  subsample_idx <- sample(1:n_cells, max_cells_for_silhouette)
} else {
  subsample_idx <- 1:n_cells
}

silhouette_scores <- sapply(resolutions, function(res) {

  cluster_col <- paste0("SCT_snn_res.", res)

  clusters <- as.numeric(Dataset3@meta.data[[cluster_col]][subsample_idx])

  coords <- Embeddings(Dataset3, reduction = "pca")[subsample_idx, ]
  coords <- coords[, 1:min(30, ncol(coords))]

  if (length(unique(clusters)) > 1) {
    dist_matrix <- dist(coords)
    sil <- silhouette(clusters, dist_matrix)
    mean(sil[, 3])
  } else {
    NA
  }
})
resolution_comparison <- data.frame(
  resolution = resolutions,
  n_clusters = sapply(resolutions, function(res) {
    cluster_col <- paste0("SCT_snn_res.", res)
    length(unique(Dataset3@meta.data[[cluster_col]]))
  }),
  silhouette_score = silhouette_scores
)

print(resolution_comparison)
```

```{r}
silhouetteplot <- ggplot(resolution_comparison,
       aes(x = resolution, y = silhouette_score)) +
  geom_line() +
  geom_point(size = 3) +
  theme_minimal() +
  labs(
    title = "Silhouette Score vs Clustering Resolution",
    x = "Resolution",
    y = "Mean Silhouette Score"
  )
ggsave("Dataset3_clustersilhouette.png", silhouetteplot, width = 8, height = 8, dpi = 300)
```

```{r}
untar("predictions.tar.gz", exdir = "predictions_output")
setwd("~/personal/15140458/project/Dataset3/predictions_output")
pred_labels <- read.csv("predicted_labels.csv", row.names = 1)
prob_mat <- read.csv("probability_matrix.csv", row.names = 1)
head(pred_labels)
# Make sure cell names match
head(rownames(pred_labels))
head(colnames(Dataset3))
```

```{r}
common_cells <- intersect(colnames(Dataset3), rownames(pred_labels))
pred_labels_Dataset3 <- pred_labels[common_cells, , drop = FALSE]
Dataset3 <- AddMetaData(Dataset3, metadata = pred_labels_Dataset3)
```

```{r}
Dataset3$CellTypist_broad <- trimws(sub(":.*", "", Dataset3$predicted_labels))
Dataset3$CellTypist_specific <- Dataset3$predicted_labels
Dataset3$CellTypist_majorbroad <- trimws(sub(":.*", "", Dataset3$majority_voting))
```

```{r}
CellTypist1 <- DimPlot(Dataset3, 
        group.by = "CellTypist_broad")
ggsave("Dataset3_celltypistbroad.png", CellTypist1, width = 8, height = 8, dpi = 300)
```

```{r}
CellTypist2 <- DimPlot(Dataset3, 
        group.by = "CellTypist_majorbroad")
ggsave("Dataset3_celltypistmajorbroad.png", CellTypist2, width = 8, height = 8, dpi = 300)
```

```{r}
UmapRes02 <- DimPlot(Dataset3,
                     reduction = "umap",
                     group.by = "SCT_snn_res.0.2")
ggsave("Dataset3_cluster02.png", UmapRes02, width = 8, height = 8, dpi = 300)
```

```{r}
PhaseUmap <- DimPlot(Dataset3, 
                     reduction = "umap", 
                     group.by = "Phase" )
ggsave("Dataset3_phaseumap.png", PhaseUmap, width = 8, height = 8, dpi = 300)
```

```{r}
DimPlot(Dataset3, 
        reduction = "pca", 
        split.by = "CellTypist_broad", 
        group.by = "SCT_snn_res.0.2" )
```

### Manual cell type annotation

```{r}
# List all genes in the dataset
available_genes <- rownames(Dataset3)

# For each category, find which genes are present
marker_genes <- list(
  "neural progenitors" = c("Prom1", "Pax6", "Sox1", "Fezf2", "Hes1", "Hes5", "Nes", "Ascl1", "Dll1", "Dll3", "Eomes", "Sox3", "Neurog2", "Emx2", "Elavl4", "Neurod1", "Dbi", "Dach1"),
  "radial glia" = c("Fabp7", "Slc1a3", "Aqp4", "Hopx", "Mki67", "Gli3", "Creb5"),
  "neural stem cells" = c("Prom1", "Pax6", "Cdh2", "Thy1", "Kit", "Epcam", "Sox2", "Nes", "Vim", "Hes1", "Notch1", "Zbtb16", "Gli1", "Gli3", "Lifr"),
  "glial cells" = c("Gfap", "Aldh1l1", "S100b", "Aqp4", "Mpb", "Plp1", "Sox10", "Cx3cr1", "Tmem119", "P2ry12", "Ptprc"),
  "glial progenitors" = c("Olig1", "Olig2", "Pdgfra", "Cspg4"),
  "mature neurons" = c("Sox5", "Tbr1", "Rbfox3", "Syt1", "Snap25", "Slc17a7", "Gad1", "Gad2", "Reln")
)

# Check presence
lapply(marker_genes, function(genes) genes[genes %in% available_genes])
```

```{r}
# Plot all neural progenitor genes
genes_to_plot <- marker_genes[["neural progenitors"]]
genes_to_plot <- genes_to_plot[genes_to_plot %in% rownames(Dataset3)]

neuralproplot <- FeaturePlot(Dataset3, features = genes_to_plot, order = TRUE, combine = FALSE)
neuralproplot <- wrap_plots(neuralproplot, ncol = 4)
ggsave("Dataset3_neuralprogmarker.png", neuralproplot, width = 16, height = 16, dpi = 300)

```

```{r}
# Plot all radial glia genes
genes_to_plot <- marker_genes[["radial glia"]]
genes_to_plot <- genes_to_plot[genes_to_plot %in% rownames(Dataset3)]

radialgliaplot <- FeaturePlot(Dataset3, features = genes_to_plot, order = TRUE, combine = FALSE)
radialgliaplot <- wrap_plots(radialgliaplot, ncol = 4)
ggsave("Dataset3_RGgmarker.png", radialgliaplot, width = 16, height = 10, dpi = 300)

```

```{r}
# Plot all glial progenitor genes
genes_to_plot <- marker_genes[["glial progenitors"]]
genes_to_plot <- genes_to_plot[genes_to_plot %in% rownames(Dataset3)]

glialproplot <- FeaturePlot(Dataset3, features = genes_to_plot, order = TRUE, combine = FALSE)
glialproplot <- wrap_plots(glialproplot, ncol = 4)
ggsave("Dataset3_GPgmarker.png", glialproplot, width = 16, height = 5, dpi = 300)
```

### Cluster labelling

```{r}
# Set active identity to the clustering resolution of interest
Idents(Dataset3) <- "SCT_snn_res.0.2"

# Verify the cell type column name
head(Dataset3$CellTypist_broad)
```

```{r}
p5 <- DimPlot(Dataset3, reduction = "umap", label = TRUE, repel = TRUE) + 
        ggtitle("Clusters (res 0.8)") + NoLegend()
p6 <- DimPlot(Dataset3, reduction = "umap", group.by = "CellTypist_majorbroad", 
              label = TRUE, repel = TRUE) + 
        ggtitle("CellTypist labels") + NoLegend()

p5 + p6
```

```{r}
# Table of counts
conf_matrix <- table(Dataset3$SCT_snn_res.1.2, Dataset3$CellTypist_majorbroad)
print(conf_matrix)

# Convert to proportions (row-wise: for each cluster, fraction of cells in each cell type)
prop_matrix <- prop.table(conf_matrix, margin = 1)  # margin=1 normalizes by row (cluster)
```

```{r}
library(pheatmap)
celltypeheatmap <- pheatmap(prop_matrix, 
         cluster_rows = FALSE, cluster_cols = FALSE,
         display_numbers = TRUE, 
         main = "Cluster vs CellTypist label (proportion)")
ggsave("Dataset3_celltypeheatmap.png", celltypeheatmap, width = 6, height = 5, dpi = 300)
```

```{r}
# Convert table to data frame for ggplot2
conf_df <- as.data.frame(conf_matrix)
colnames(conf_df) <- c("Cluster", "CellType", "Count")

ggplot(conf_df, aes(x = Cluster, y = Count, fill = CellType)) +
  geom_bar(stat = "identity", position = "fill") +   # "fill" gives proportions
  labs(y = "Proportion", title = "Cell type composition by cluster") +
  theme_minimal()
```

```{r}
saveRDS(Dataset3, file = "Dataset3_E12.5_seurat_object.rds")
```
