setwd("~/personal/15140458/project/Dataset1/CGE/raw_data")

# Load required libraries
library(Matrix)
library(readr)

# Paths to 10x files
matrix_file <- "matrix.mtx.gz"
features_file <- "features.tsv.gz" 
barcodes_file <- "barcodes.tsv.gz"

# Load the data
expr_matrix <- readMM(matrix_file)            # sparse matrix (genes x cells)
genes <- read_tsv(features_file, col_names = FALSE)
barcodes <- read_tsv(barcodes_file, col_names = FALSE)

# Assign row and column names
rownames(expr_matrix) <- genes$X2  # gene symbols
colnames(expr_matrix) <- barcodes$X1  # cell barcodes

# Transpose to get cells as rows, genes as columns
expr_matrix_t <- t(expr_matrix)

# Convert to dense 
expr_dense <- as.matrix(expr_matrix_t)

# Write to CSV with cell barcodes as row names
write.csv(expr_dense, "cell_by_gene_matrix.csv", row.names = TRUE)
