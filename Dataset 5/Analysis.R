set.seed(111)

# Loading necessary packages

library(dplyr)
library(Seurat)
library(patchwork)
library(DoubletFinder)
library(cluster)

# Loading the data in the environment

untar("GSM4610591_E14_CtxGE_Emx1_N1.tar.gz", exdir = "Dataset5")
untar("GSM4610592_E14_CtxGE_Emx1_N2.tar.gz", exdir = "Dataset5")
untar("GSM4610593_E14_CtxGE_Emx1_N3.tar.gz", exdir = "Dataset5")
untar("GSM4610594_E17_SVZ_Emx1.tar.gz", exdir = "Dataset5")

E14_1 <- Read10X(data.dir = "/zfs/omics/personal/15140458/project/Dataset5/E14_CtxGE_Emx1_N1")
E14_2 <- Read10X(data.dir = "/zfs/omics/personal/15140458/project/Dataset5/E14_CtxGE_Emx1_N2")
E14_3 <- Read10X(data.dir = "/zfs/omics/personal/15140458/project/Dataset5/E14_CtxGE_Emx1_N3")
E17 <- Read10X(data.dir = "/zfs/omics/personal/15140458/project/Dataset5/E17_SVZ_Emx1")

# Creating a seurat object

E14_1_seurat <- CreateSeuratObject(counts = E14_1)
E14_2_seurat <- CreateSeuratObject(counts = E14_2)
E14_3_seurat <- CreateSeuratObject(counts = E14_3)
E17_seurat  <- CreateSeuratObject(counts = E17)

# CellTypist cell type labelling

untar("E14-1.tar.gz", exdir = "CellTypist_E14-1")
untar("E14-2.tar.gz", exdir = "CellTypist_E14-2")
untar("E14-3.tar.gz", exdir = "CellTypist_E14-3")
untar("E17.tar.gz", exdir = "CellTypist_E17")

setwd("~/personal/15140458/project/Dataset5/CellTypist_E14-1")
pred_labels1 <- read.csv("predicted_labels.csv", row.names = 1)
prob_mat1 <- read.csv("probability_matrix.csv", row.names = 1)
common_cells <- intersect(colnames(E14_1_seurat), rownames(pred_labels1))
pred_labels_E14_1_seurat <- pred_labels1[common_cells, , drop = FALSE]
E14_1_seurat <- AddMetaData(E14_1_seurat, metadata = pred_labels_E14_1_seurat)
E14_1_seurat$CellTypist_broad <- trimws(sub(":.*", "", E14_1_seurat$predicted_labels))
E14_1_seurat$CellTypist_specific <- E14_1_seurat$predicted_labels
E14_1_seurat$CellTypist_majorbroad <- trimws(sub(":.*", "", E14_1_seurat$majority_voting))

setwd("~/personal/15140458/project/Dataset5/CellTypist_E14-2")
pred_labels2 <- read.csv("predicted_labels.csv", row.names = 1)
prob_mat2 <- read.csv("probability_matrix.csv", row.names = 1)
common_cells <- intersect(colnames(E14_2_seurat), rownames(pred_labels2))
pred_labels_E14_2_seurat <- pred_labels2[common_cells, , drop = FALSE]
E14_2_seurat <- AddMetaData(E14_2_seurat, metadata = pred_labels_E14_2_seurat)
E14_2_seurat$CellTypist_broad <- trimws(sub(":.*", "", E14_2_seurat$predicted_labels))
E14_2_seurat$CellTypist_specific <- E14_2_seurat$predicted_labels
E14_2_seurat$CellTypist_majorbroad <- trimws(sub(":.*", "", E14_2_seurat$majority_voting))

setwd("~/personal/15140458/project/Dataset5/CellTypist_E14-3")
pred_labels3 <- read.csv("predicted_labels.csv", row.names = 1)
prob_mat3 <- read.csv("probability_matrix.csv", row.names = 1)
common_cells <- intersect(colnames(E14_3_seurat), rownames(pred_labels3))
pred_labels_E14_3_seurat <- pred_labels3[common_cells, , drop = FALSE]
E14_3_seurat <- AddMetaData(E14_3_seurat, metadata = pred_labels_E14_3_seurat)
E14_3_seurat$CellTypist_broad <- trimws(sub(":.*", "", E14_3_seurat$predicted_labels))
E14_3_seurat$CellTypist_specific <- E14_3_seurat$predicted_labels
E14_3_seurat$CellTypist_majorbroad <- trimws(sub(":.*", "", E14_3_seurat$majority_voting))

setwd("~/personal/15140458/project/Dataset5/CellTypist_E17")
pred_labels4 <- read.csv("predicted_labels.csv", row.names = 1)
prob_mat4 <- read.csv("probability_matrix.csv", row.names = 1)
common_cells <- intersect(colnames(E17_seurat), rownames(pred_labels4))
pred_labels_E17_seurat <- pred_labels4[common_cells, , drop = FALSE]
E17_seurat <- AddMetaData(E17_seurat, metadata = pred_labels_E17_seurat)
E17_seurat$CellTypist_broad <- trimws(sub(":.*", "", E17_seurat$predicted_labels))
E17_seurat$CellTypist_specific <- E17_seurat$predicted_labels
E17_seurat$CellTypist_majorbroad <- trimws(sub(":.*", "", E17_seurat$majority_voting))

# Creating a merged seurat object

E14_1_seurat$stage <- "E14"
E14_2_seurat$stage <- "E14"
E14_3_seurat$stage <- "E14"
E17_seurat$stage  <- "E17"

E14_1_seurat$orig.ident <- "E14.1"
E14_2_seurat$orig.ident <- "E14.2"
E14_3_seurat$orig.ident <- "E14.3"
E17_seurat$orig.ident <- "E17"

dataset5 <- merge(E14_1_seurat,
                  y = c(E14_2_seurat, E14_3_seurat, E17_seurat),
                  add.cell.ids = c("E14_1", "E14_2", "E14_3", "E17")
                  )

# Quality check metrics

dataset5[["percent.mt"]] <- PercentageFeatureSet(dataset5, pattern = "^mt-")
dataset5[["percent.ribo"]] <- PercentageFeatureSet(dataset5, pattern = "^Rpl|^Rps")
dataset5[["percent.hb"]] <- PercentageFeatureSet(dataset5, pattern = "^Hba|^Hbb")

cat("\nQC Metric Distributions:\n")
cat("nCount_RNA (UMI):\n")
cat("  Median:", median(dataset5$nCount_RNA),
    "| Q1-Q3:", quantile(dataset5$nCount_RNA,0.25), "-",
    quantile(dataset5$nCount_RNA,0.75), "\n")
cat("nFeature_RNA (genes):\n")
cat("  Median:", median(dataset5$nFeature_RNA),
    "| Q1-Q3:", quantile(dataset5$nFeature_RNA,0.25), "-",
    quantile(dataset5$nFeature_RNA,0.75), "\n")
cat("percent.mt:\n")
cat("  Median:", round(median(dataset5$percent.mt),2), "%",
    "| 95th percentile:", round(quantile(dataset5$percent.mt,0.95),2), "%\n")
cat("percent.ribo:\n")
cat("  Median:", round(median(dataset5$percent.ribo),2), "%\n")
cat("percent.hb:\n")
cat("  Median:", round(median(dataset5$percent.hb),2), "%\n")

# Quality check violin plots

p_qc <- VlnPlot(
  dataset5,
  features = c("nCount_RNA",
               "nFeature_RNA",
               "percent.mt",
               "percent.ribo",
               "percent.hb"),
  ncol = 5,
  pt.size = 0.1
)
ggsave("dataset5_QC_plots.png", p_qc, width = 20, height = 8, dpi = 300)

# Quality check scatterplots

p1 <- FeatureScatter(dataset5, feature1="nCount_RNA", feature2="nFeature_RNA")
p2 <- FeatureScatter(dataset5, feature1="nCount_RNA", feature2="percent.mt")
p3 <- FeatureScatter(dataset5, feature1="percent.mt", feature2="percent.ribo")
p_scatter <- p1 + p2 + p3
ggsave("dataset5_QC_scatterplots.png", p_scatter, width = 24, height = 8, dpi = 300)

# Quality check pass scatterplot

nfeature_min <- 300
nfeature_max <- 6000
ncount_min <- 1000
ncount_max <- 50000
mt_thresh <- 20
hb_thresh <- 1
qc_df <- data.frame(
  nCount_RNA = dataset5$nCount_RNA,
  nFeature_RNA = dataset5$nFeature_RNA
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
ggsave("dataset5_QC_threshold.png", p4, width = 8, height = 8, dpi = 300)

# Subsetting the seurat object

dataset5 <- subset(
  dataset5,
  subset =
    nCount_RNA >= ncount_min &
    nCount_RNA <= ncount_max &
    nFeature_RNA >= nfeature_min &
    nFeature_RNA <= nfeature_max &
    percent.mt < mt_thresh &
    percent.hb < hb_thresh
)
cat("Cells remaining after QC:", ncol(dataset5), "\n")

# Doublet detection 

split_list <- SplitObject(dataset5, split.by = "orig.ident")

obj1 <- split_list[["E14.1"]]
obj2 <- split_list[["E14.2"]]
obj3 <- split_list[["E14.3"]]
obj4 <- split_list[["E17"]]

## Put objects into a list
obj_list <- list(obj1, obj2, obj3, obj4)

## Process each object
obj_list <- lapply(obj_list, function(obj) {
  obj <- SCTransform(obj, method="glmGamPoi")
  obj <- RunPCA(obj)
  obj <- RunUMAP(obj, dims = 1:20)
  return(obj)
})

## Assign back 
obj1 <- obj_list[[1]]
obj2 <- obj_list[[2]]
obj3 <- obj_list[[3]]
obj4 <- obj_list[[4]]

nExp_poi <- round(0.02 * ncol(obj1))  # expected doublets
# Parameter sweep to find optimal pK
sweep.res.list <- paramSweep(obj1, PCs = 1:20, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
ggplot(bcmvn, aes(x = pK, y = BCmetric)) + geom_line() + geom_point()
# results is 0.04
obj1 <- doubletFinder(
  obj1,
  PCs = 1:20,
  pN = 0.25,
  pK = 0.04,     
  nExp = nExp_poi,
  sct = TRUE
)

nExp_poi <- round(0.03 * ncol(obj2))  # expected doublets
# Parameter sweep to find optimal pK
sweep.res.list <- paramSweep(obj2, PCs = 1:20, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
ggplot(bcmvn, aes(x = pK, y = BCmetric)) + geom_line() + geom_point()
# results is 0.14
obj2 <- doubletFinder(
  obj2,
  PCs = 1:20,
  pN = 0.25,
  pK = 0.14,     
  nExp = nExp_poi,
  sct = TRUE
)

nExp_poi <- round(0.04 * ncol(obj3))  # expected doublets
# Parameter sweep to find optimal pK
sweep.res.list <- paramSweep(obj3, PCs = 1:20, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
ggplot(bcmvn, aes(x = pK, y = BCmetric)) + geom_line() + geom_point()
# results is 0.27
obj3 <- doubletFinder(
  obj3,
  PCs = 1:20,
  pN = 0.25,
  pK = 0.27,     
  nExp = nExp_poi,
  sct = TRUE
)

nExp_poi <- round(0.02 * ncol(obj4))  # expected doublets
# Parameter sweep to find optimal pK
sweep.res.list <- paramSweep(obj4, PCs = 1:20, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
ggplot(bcmvn, aes(x = pK, y = BCmetric)) + geom_line() + geom_point()
# results is 0.2
obj4 <- doubletFinder(
  obj4,
  PCs = 1:20,
  pN = 0.25,
  pK = 0.2,     
  nExp = nExp_poi,
  sct = TRUE
)

# Transferring doublet label

dbl1 <- obj1$DF.classifications_0.25_0.04_48
dbl2 <- obj2$DF.classifications_0.25_0.14_110
dbl3 <- obj3$DF.classifications_0.25_0.27_194
dbl4 <- obj4$DF.classifications_0.25_0.2_42
obj1$doublet_status <- dbl1
obj2$doublet_status <- dbl2
obj3$doublet_status <- dbl3
obj4$doublet_status <- dbl4

dataset5 <- merge(obj1,
                  y = c(obj2, obj3, obj4),
                  )
dataset5 <- subset(dataset5, subset = doublet_status == "Singlet")
dataset5 <- JoinLayers(dataset5, assay = "RNA")

# Initial normalisation 

dataset5 <- SCTransform(dataset5,
                   method="glmGamPoi")

# Initial principal component analysis

dataset5 <- RunPCA(dataset5)
PCAplot <- ElbowPlot(dataset5, ndims = 50)
ggsave("dataset5_PCAelbow.png", PCAplot, width = 8, height = 8, dpi = 300)

PCAplot2 <- DimPlot(dataset5, 
                    reduction = "pca",
                    group.by = "stage")
ggsave("dataset5_PCA.png", PCAplot2, width = 8, height = 8, dpi = 300)

# UMAP

dataset5 <- RunUMAP(dataset5, dims = 1:20)

# Check for batch effects
batch_effect <- DimPlot(dataset5, reduction = "umap", group.by = "orig.ident")
ggsave("dataset5_batch.png", batch_effect, width = 8, height = 8, dpi = 300)

# Cell cycle scoring

## Score cells for cell cycle
dataset5 <- CellCycleScoring(dataset5,
                           g2m.features = g2m_genes, 
                           s.features = s_genes)

## PCA using only cell cycle related genes as variable features to see cell cycle effects
dataset5_cc <- RunPCA(dataset5, features=c(s_genes,g2m_genes))
PCAcc <- DimPlot(dataset5_cc, reduction="pca", group.by="Phase")
ggsave("dataset5_PCAcc.png", PCAcc, width = 8, height = 8, dpi = 300)

## Doing normalisation with cell cycle difference between S and G2M regression
dataset5$CC.Difference <- dataset5$S.Score - dataset5$G2M.Score
dataset5 <- SCTransform(dataset5, vars.to.regress = "CC.Difference", method="glmGamPoi")

## Checking regressed normalisation
dataset5_cc2 <- RunPCA(dataset5, features=c(s_genes,g2m_genes))
PCAcc2 <- DimPlot(dataset5_cc2, reduction="pca", group.by="Phase")
ggsave("dataset5_PCAcc2.png", PCAcc2, width = 8, height = 8, dpi = 300)

## PCA and UMAP again on new normalisation
dataset5 <- RunPCA(dataset5)
dataset5 <- RunUMAP(dataset5, dims = 1:20)

# Quality check visualisation on UMAP and PC plot

## Distribution of embryonic stage on UMAP
UMAPplot1 <- DimPlot(dataset5,
                     reduction "umap",
                     group.by = "stage")
ggsave("dataset5_UMAP1.png", UMAPplot1, width = 8, height = 8, dpi = 300)

## Hemoglobin gene expression distribution on UMAP
UMAP_hb <- FeaturePlot(dataset5,
                       reduction = "umap",
                       features = "percent.hb")
ggsave("dataset5_UMAP_hb.png", UMAP_hb, width = 8, height = 8, dpi = 300)

## Distribution of embryonic stage on PCs
PCAplot3 <- DimPlot(dataset5, 
                    reduction = "pca",
                    group.by = "stage")
ggsave("dataset5_PCA3.png", PCAplot3, width = 8, height = 8, dpi = 300)

## Distribution cell cycle phase on UMAP
UMAP_cc <- DimPlot(dataset5, 
                   reduction = "umap",
                   group.by = "Phase")
ggsave("dataset5_UMAP_cc.png", UMAP_cc, width = 8, height = 8, dpi = 300)

# Clustering

dataset5 <- FindNeighbors(dataset5, dims = 1:20)
dataset5 <- FindClusters(dataset5,
                         resolution = c(0.2, 0.4, 0.6, 0.8, 1, 1.2),
                         algorithm = 4,
                         random.seed = 1)

resolutions <- c(0.2, 0.4, 0.6, 0.8, 1, 1.2)
plot_list <- lapply(resolutions, function(res) {
  cluster_col <- paste0("SCT_snn_res.", res)
  n_clusters <- length(unique(dataset5@meta.data[[cluster_col]]))
  DimPlot(
    dataset5,
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
ggsave("Dataset5_allcluster.png", combined_resolutions, width = 15, height = 8, dpi = 300)

# Cluster silhouette score

resolutions <- c(0.2, 0.4, 0.6, 0.8, 1, 1.2)
n_cells <- ncol(dataset5)
max_cells_for_silhouette <- 5000

if (n_cells > max_cells_for_silhouette) {
  set.seed(42)
  subsample_idx <- sample(1:n_cells, max_cells_for_silhouette)
} else {
  subsample_idx <- 1:n_cells
}

silhouette_scores <- sapply(resolutions, function(res) {

  cluster_col <- paste0("SCT_snn_res.", res)

  clusters <- as.numeric(dataset5@meta.data[[cluster_col]][subsample_idx])

  coords <- Embeddings(dataset5, reduction = "pca")[subsample_idx, ]
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
    length(unique(dataset5@meta.data[[cluster_col]]))
  }),
  silhouette_score = silhouette_scores
)

print(resolution_comparison)

# Exploring known neural population markers

available_genes <- rownames(dataset5)
marker_genes <- list(
  "neural progenitors" = c("Pax6", "Hes5", "Neurog2"),
  "radial glia" = c("Fabp7", "Hopx", "Slc1a3"),
  "neural stem cells" = c("Prom1", "Sox2", "Nes"),
  "glial cells" = c("Gfap", "Cx3cr1", "Plp1"),
  "glial progenitors" = c("Olig2", "Pdgfra", "Cspg4"),
  "mature neurons" = c("Rbfox3", "Gad2", "Slc17a7")
)
lapply(marker_genes, function(genes) genes[genes %in% available_genes])

## Plot all neural progenitor genes
genes_to_plot <- marker_genes[["neural progenitors"]]
genes_to_plot <- genes_to_plot[genes_to_plot %in% rownames(dataset5)]
neuralproplot <- FeaturePlot(dataset5, features = genes_to_plot, order = TRUE, combine = FALSE)
neuralproplot <- wrap_plots(neuralproplot, ncol = 3)
ggsave("Dataset5_neuralprogmarker.png", neuralproplot, width = 15, height = 5, dpi = 300)

## Plot all glial progenitor genes
genes_to_plot <- marker_genes[["glial progenitors"]]
genes_to_plot <- genes_to_plot[genes_to_plot %in% rownames(dataset5)]
glialproplot <- FeaturePlot(dataset5, features = genes_to_plot, order = TRUE, combine = FALSE)
glialproplot <- wrap_plots(glialproplot, ncol = 3)
ggsave("Dataset5_GPgmarker.png", glialproplot, width = 15, height = 5, dpi = 300)

## Plot all glial cell genes
genes_to_plot <- marker_genes[["glial cells"]]
genes_to_plot <- genes_to_plot[genes_to_plot %in% rownames(dataset5)]
glialplot <- FeaturePlot(dataset5, features = genes_to_plot, order = TRUE, combine = FALSE)
glialplot <- wrap_plots(glialplot, ncol = 2)
ggsave("Dataset5_Ggmarker.png", glialplot, width = 10, height = 5, dpi = 300)

## Plot all mature neuron genes
genes_to_plot <- marker_genes[["mature neurons"]]
genes_to_plot <- genes_to_plot[genes_to_plot %in% rownames(dataset5)]
neuronplot <- FeaturePlot(dataset5, features = genes_to_plot, order = TRUE, combine = FALSE)
neuronplot <- wrap_plots(neuronplot, ncol = 3)
ggsave("Dataset5_MNgmarker.png", neuronplot, width = 15, height = 5, dpi = 300)

# Visualising the most optimal clustering results and CellTypist labelling on UMAP

Idents(dataset5) <- "SCT_snn_res.0.4"
p5 <- DimPlot(dataset5, reduction = "umap", label = TRUE, repel = TRUE) + 
        ggtitle("Clusters (res 0.4)") + NoLegend()
p6 <- DimPlot(dataset5, reduction = "umap", group.by = "CellTypist_majorbroad", 
              label = TRUE, repel = TRUE) + 
        ggtitle("CellTypist labels") + NoLegend()
clustercelltype <- p5 + p6
ggsave("Dataset5_clustercelltype.png", clustercelltype, width = 10, height = 5, dpi = 300)

# Exploring cluster identity

## Table of counts
conf_matrix <- table(dataset5$SCT_snn_res.0.4, dataset5$CellTypist_majorbroad)

## Convert to proportions (row-wise: for each cluster, fraction of cells in each cell type)
prop_matrix <- prop.table(conf_matrix, margin = 1)  # margin=1 normalizes by row (cluster)

library(pheatmap)
celltypeheatmap <- pheatmap(prop_matrix, 
         cluster_rows = FALSE, cluster_cols = FALSE,
         display_numbers = TRUE, 
         main = "Cluster vs CellTypist label (proportion)")
ggsave("dataset5_celltypeheatmap.png", celltypeheatmap, width = 6, height = 5, dpi = 300)

# Saving the seurat object

saveRDS(dataset5, file = "Dataset5_processed_seurat_object.rds")

# How many cells of each type belong to each embryonic stage?
table(dataset5$CellTypist_majorbroad, dataset5$stage)

# Differential expression analysis

dataset5 <- JoinLayers(dataset5, assay = "RNA")
dataset5 <- NormalizeData(dataset5, assay = "RNA") # Make sure RNA assay is log-normalized
DefaultAssay(dataset5) <- "RNA" # Set RNA as default assay for DEA
dataset5$progenitor_status <- ifelse( # Create progenitor status grouping
  dataset5$CellTypist_majorbroad %in% c("Neuroblast", "Glioblast"),
  "Progenitor",
  "Other"
)
Idents(dataset5) <- "progenitor_status" # Set active identity to progenitor status

progenitor_markers <- FindMarkers(
  dataset5,
  ident.1 = "Progenitor",
  ident.2 = "Other",
  test.use = "LR",
  latent.vars = "orig.ident",
  assay = "RNA",
  min.pct = 0.25,        # gene must be expressed in 25% of either group
  logfc.threshold = 0.25 # minimum log2FC to test
)
