set.seed(123)
library(dplyr)
library(Seurat)
library(patchwork)

setwd("~/personal/15140458/project/Dataset1/CGE")
CGE.data <- Read10X(data.dir = "~/personal/15140458/project/Dataset1/CGE/raw_data")
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

