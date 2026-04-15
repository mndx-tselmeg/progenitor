# Create individual Seurat objects from each dgCMatrix
seurat_CGE_13 <- CreateSeuratObject(counts = CGE_13, project = "CGE_13")
seurat_MGE_13 <- CreateSeuratObject(counts = MGE_13, project = "MGE_13")
seurat_LGE_13 <- CreateSeuratObject(counts = LGE_13, project = "LGE_13")
seurat_CGE_15 <- CreateSeuratObject(counts = CGE_15, project = "CGE_15")
seurat_MGE_15 <- CreateSeuratObject(counts = MGE_15, project = "MGE_15")
seurat_LGE_15 <- CreateSeuratObject(counts = LGE_15, project = "LGE_15")

# Add sample identity as metadata 
seurat_CGE_13$sample <- "CGE_13"
seurat_MGE_13$sample <- "MGE_13"
seurat_LGE_13$sample <- "LGE_13"
seurat_CGE_15$sample <- "CGE_15"
seurat_MGE_15$sample <- "MGE_15"
seurat_LGE_15$sample <- "LGE_15"

dataset6 <- merge(
  x = seurat_CGE_13,
  y = list(seurat_MGE_13, seurat_LGE_13, seurat_CGE_15, seurat_MGE_15, seurat_LGE_15),
  add.cell.ids = c("CGE_13", "MGE_13", "LGE_13", "CGE_15", "MGE_15", "LGE_15"),
  project = "Merged"
)

