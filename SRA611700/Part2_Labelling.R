set.seed(123)
library(Seurat)
library(sctransform)
obj <- SCTransform(object = seurat_obj_clean,
                   vst.flavor = "v2",
                   vars.to.regress = "percent.ribo")

obj <- RunPCA(obj)
DimPlot(obj, reduction = "pca")
DimHeatmap(obj, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(obj, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(obj, ndims = 50 )

obj <- FindNeighbors(obj, dims = 1:25)
obj <- FindClusters(obj, resolution = 1)
obj <- RunUMAP(obj, dims = 1:25)
DimPlot(obj, reduction = "umap")

obj <- RunTSNE(obj, dims = 1:25, reduction = "pca")
DimPlot(obj, reduction = "tsne", label = TRUE) 

---

# Let's try to use ScType for labelling the clusters

library(HGNChelper)
# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
# DB file
db_ <- "ScTypeDB_full.xlsx";
tissue <- "Brain" 

# prepare gene sets
gs_list <- gene_sets_prepare(db_, tissue)

# Extract scaled data from SCT assay
scRNAseqData_scaled <- as.matrix(obj[["SCT"]]$scale.data)

# Run ScType
es.max <- sctype_score(scRNAseqData = scRNAseqData_scaled, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

# Merge by cluster and get top scoring type per cluster
sctype_scores <- do.call("rbind", lapply(unique(obj$SCT_snn_res.1), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(obj@meta.data[obj$SCT_snn_res.1==cl, ])]), decreasing = TRUE)
  data.frame(cluster = cl, type = names(es.max.cl)[1], scores = es.max.cl[1], 
             ncells = sum(obj$SCT_snn_res.1==cl))
}))

# Set low-confidence clusters to "Unknown"
sctype_scores$type[sctype_scores$scores < sctype_scores$ncells/4] <- "Unknown"

# View results
sctype_scores[,1:3]

# Add cell type classifications to metadata
obj$sctype_classification <- ""
for(j in unique(sctype_scores$cluster)){
  obj$sctype_classification[obj$seurat_clusters == j] <- as.character(sctype_scores$type[sctype_scores$cluster == j][1])
}

# Plot results
DimPlot(obj, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'sctype_classification')

# Count cells per cell type
table(obj$sctype_classification)

---

# Let's try scCATCH for labelling the clusters
library(scCATCH)
objj <- createscCATCH(data = obj[["SCT"]]@data, cluster = as.character(Idents(obj)))
cellmatch_new <- cellmatch[cellmatch$species == "Mouse" & cellmatch$tissue %in% c("Brain", "Cerebellum", "Fetal brain", "Hippocampus", "Neural tube"), ]
objj <- findmarkergene(object = obj, 
                      if_use_custom_marker = TRUE, 
                      marker = cellmatch_new,
                      use_method = "2")
objj <- findcelltype(objj)


