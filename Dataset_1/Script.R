set.seed(123)
library(dplyr)
library(Seurat)
library(patchwork)

setwd("~/personal/15140458/project/Dataset1/CGE")
CGE.data <- Read10X(data.dir = "~/personal/15140458/project/Dataset1/CGE")
CGE <- CreateSeuratObject(counts = CGE.data,
                          project = "CGE_E12.5")

CGE[["percent.mt"]] <- PercentageFeatureSet(CGE, pattern = "^mt-")
CGE[["percent.ribo"]] <- PercentageFeatureSet(CGE, pattern = "^Rpl|^Rps")
CGE[["percent.hb"]] <- PercentageFeatureSet(CGE, pattern = "^Hba|^Hbb")

# Visualize QC metrics as a violin plot
VlnPlot(
  CGE,
  features = c("nFeature_RNA", 
               "nCount_RNA", 
               "percent.mt", 
               "percent.ribo", 
               "percent.hb"),
  ncol = 5
)

ggsave(
  filename = "CGE_E12.5_QC_violin.png",
  width = 12,
  height = 6,
  dpi = 300
)

plot1 <- FeatureScatter(CGE, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot2 <- FeatureScatter(CGE, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot1 + plot2
ggsave(
  filename = "CGE_E12.5_QC_scatterplot.png",
  width = 15,
  height = 6,
  dpi = 300
)

counts <- GetAssayData(CGE, layer = "counts")
gene_means <- rowMeans(counts)
gene_vars  <- apply(counts, 1, var)

png("CGE_12.5_mean_variance.png", width = 2000, height = 1800, res = 300)

plot(gene_means, gene_vars,
     pch = 16, cex = 0.4,
     xlab = "Mean UMI count",
     ylab = "Variance of UMI count")

abline(0, 1, col = "red", lwd = 2)  # Poisson expectation

dev.off()

if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
library(DoubletFinder)

CGE_v2 <- SCTransform(CGE)
CGE_v2 <- RunPCA(CGE_v2)
ElbowPlot(CGE_v2)

CGE_v2 <- RunUMAP(CGE_v2, dims = 1:10)

nExp_poi <- round(0.04 * ncol(CGE_v2))  # assuming 4% expected doublets
# Parameter sweep to find optimal pK
sweep.res.list <- paramSweep(CGE_v2, PCs = 1:10, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

# Visualize pK metric
ggplot(bcmvn, aes(x = pK, y = BCmetric)) + geom_line() + geom_point()

# Run DoubletFinder
CGE_v2 <- doubletFinder(
  CGE_v2,
  PCs = 1:10,
  pN = 0.25,
  pK = 0.04,     
  nExp = nExp_poi,
  sct = TRUE
)

png("CGE_12.5_doublet.png", width = 2000, height = 1800, res = 300)
DimPlot(CGE_v2, 
        group.by = "DF.classifications_0.25_0.04_172")
dev.off()

png("CGE_12.5_doublet_score.png", width = 2000, height = 1800, res = 300)
FeaturePlot(
  CGE_v2,
  features = "pANN_0.25_0.04_172",
  reduction = "umap"
)
dev.off()

df_col <- "DF.classifications_0.25_0.04_172"
CGE_v2 <- CGE_v2[, CGE_v2@meta.data[[df_col]] == "Singlet"]

# Get all gene names
genes <- rownames(CGE_v2)
# Find hemoglobin genes
hb_genes <- grep("^Hba|^Hbb", genes, value = TRUE)
# Get counts matrix
counts <- GetAssayData(CGE_v2, layer = "counts")
# Keep genes with at least 1 count in any cell
hb_expressed <- hb_genes[rowSums(counts[hb_genes, ]) > 0]
hb_expressed

CGE_v2$Hb_score <- colSums(
  GetAssayData(CGE_v2, layer = "counts")[hb_expressed, ]
)
png("CGE_12.5_Hb_score.png", width = 2000, height = 1800, res = 300)
FeaturePlot(CGE_v2, 
            features = "Hb_score")
dev.off()

png("CGE_12.5_Ascl1.png", width = 2000, height = 1800, res = 300)
neural_genes <- c("Ascl1")
FeaturePlot(CGE_v2, 
            features = neural_genes,
            reduction = "umap")
dev.off()

# Find mitochondrial genes
mt_genes <- grep("^mt-", genes, value = TRUE)
# Keep genes with at least 1 count in any cell
mt_expressed <- mt_genes[rowSums(counts[mt_genes, ]) > 0]
CGE_v2$mt_score <- colSums(
  GetAssayData(CGE_v2, layer = "counts")[mt_expressed, ]
)
png("CGE_12.5_mt_score.png", width = 2000, height = 1800, res = 300)
FeaturePlot(CGE_v2, 
            features = "mt_score")
dev.off()

# Score cells for cell cycle
CGE_v2 <- CellCycleScoring(CGE_v2, 
                                 g2m.features = g2m_genes, 
                                 s.features = s_genes)
# View cell cycle scores and phases assigned to cells
head(CGE_v2@meta.data)   

png("CGE_12.5_phase1.png", width = 2000, height = 1000, res = 300)
DimPlot(CGE_v2,
        reduction = "pca",
        group.by= "Phase",
        split.by = "Phase")
dev.off()

png("CGE_12.5_phase2.png", width = 2000, height = 1800, res = 300)
DimPlot(CGE_v2,
        reduction = "pca",
        group.by= "Phase")
dev.off()

CGE_v2$CC.Difference <- CGE_v2$S.Score - CGE_v2$G2M.Score
CGE_v2 <- SCTransform(CGE_v2,
                      vars.to.regress = c("CC.Difference",
                                          "percent.mt"))

CGE_v2 <- RunPCA(CGE_v2, features = c(s_genes, g2m_genes))
ElbowPlot(CGE_v2)

png("CGE_12.5_phase2_2.png", width = 2000, height = 1800, res = 300)
DimPlot(CGE_v2,
        reduction = "pca",
        group.by= "Phase")
dev.off()

png("CGE_12.5_phase1_2.png", width = 2000, height = 1000, res = 300)
DimPlot(CGE_v2,
        reduction = "pca",
        group.by= "Phase",
        split.by = "Phase")
dev.off()

set.seed(123)
CGE_v2 <- RunUMAP(CGE_v2, dims = 1:5)
DimPlot(CGE_v2,
        reduction = "umap")

CGE_v2 <- FindNeighbors(CGE_v2, dims = 1:5)
CGE_v2 <- FindClusters(CGE_v2, resolution = c(0.2, 0.4, 0.6, 0.8, 1, 1.2))

library(clustree)
png("CGE_12.5_cluster_decision.png", width = 2000, height = 2500, res = 300)
clustree(CGE_v2)
dev.off()

Idents(CGE_v2) <- "SCT_snn_res.0.8"
DimPlot(CGE_v2, 
        label = TRUE,
        reduction ="umap" )



                                    

