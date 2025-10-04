set.seed(123)

library(irlba)
library(Matrix)
library(Seurat)
library(celda)
library(SingleCellExperiment)
library(scDblFinder)
library(decontX)

rownames(sm) <- sub("_[^_]+$", "", rownames(sm))
rownames(sm) <- make.unique(rownames(sm))

# Create SingleCellExperiment object from matrix
sce <- SingleCellExperiment::SingleCellExperiment(
  assays = list(counts = sm)
)

cat("Running DecontX on your matrix...\n")
sce_decontx <- decontX(sce)

umap <- reducedDim(sce_decontx, "decontX_UMAP")
plotDimReduceCluster(x = sce_decontx$decontX_clusters,
    dim1 = umap[, 1], dim2 = umap[, 2])

plotDecontXContamination(sce_decontx)

# Create Seurat object with decontaminated counts
decontaminated_counts <- sce_decontx@assays@data$decontXcounts
contamination_estimates <- sce_decontx$decontX_contamination

seurat_obj <- CreateSeuratObject(
  counts = decontaminated_counts,
  min.cells = 3,      # Minimum cells per gene
  min.features = 200  # Minimum genes per cell
)
# Add contamination estimates as metadata
seurat_obj$contamination <- contamination_estimates[colnames(seurat_obj)]

# Calculate QC metrics
rp_genes <- grep("^Rp[sl]", rownames(seurat_obj), value = TRUE, ignore.case = TRUE)
seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, features = rp_genes)

mt_genes <- grep("^mt-", rownames(seurat_obj), value = TRUE, ignore.case = TRUE)
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, features = mt_genes)

hb_genes <- grep("^Hb[ab]", rownames(seurat_obj), value = TRUE, ignore.case = TRUE)
seurat_obj[["percent.hb"]] <- PercentageFeatureSet(seurat_obj, features = hb_genes)

p1 <- FeatureScatter(seurat_obj, 
                     feature1 = "percent.ribo", 
                     feature2 = "nFeature_RNA") + NoLegend()
p2 <- FeatureScatter(seurat_obj, 
                     feature1 = "percent.mt", 
                     feature2 = "nFeature_RNA") + NoLegend()
p3 <- FeatureScatter(seurat_obj, 
                     feature1 = "percent.hb", 
                     feature2 = "nFeature_RNA") + NoLegend()
p1 | p2 | p3

seurat_obj <- subset(seurat_obj, subset = percent.hb < 1 & percent.mt < 5)

# Doublet detection
sce <- as.SingleCellExperiment(seurat_obj, assay = "RNA")
sce <- scDblFinder(sce)

# Add doublet results to Seurat object metadata
seurat_obj$scDblFinder.class <- sce$scDblFinder.class
seurat_obj$scDblFinder.score <- sce$scDblFinder.score

# Visualize doublets
vis_obj <- NormalizeData(seurat_obj)
vis_obj <- FindVariableFeatures(vis_obj)
vis_obj <- ScaleData(vis_obj)
vis_obj <- RunPCA(vis_obj, npcs = 20, verbose = FALSE)
vis_obj <- RunUMAP(vis_obj, dims = 1:15, verbose = FALSE)

# Add doublet info to visualization object
vis_obj$scDblFinder.class <- seurat_obj$scDblFinder.class
vis_obj$scDblFinder.score <- seurat_obj$scDblFinder.score

# Create plots
p1 <- DimPlot(vis_obj, group.by = "scDblFinder.class") +
  ggtitle("Doublet Classes")
p2 <- FeaturePlot(vis_obj, features = "scDblFinder.score") +
  ggtitle("Doublet Scores")
p1 | p2

# Filter to keep only singlets
seurat_obj_clean <- subset(seurat_obj, subset = scDblFinder.class == "singlet")

# Save the cleaned Seurat object
saveRDS(seurat_obj_clean, file = "decontaminated_cleaned_seurat.rds")
