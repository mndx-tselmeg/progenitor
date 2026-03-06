setwd("~/personal/15140458/project/Dataset1/CGE/raw_data")

# Load required libraries
library(Matrix)
library(readr)

# Paths to your 10x files
matrix_file <- "matrix.mtx.gz"
features_file <- "features.tsv.gz"  # or genes.tsv.gz
barcodes_file <- "barcodes.tsv.gz"

# Load the data
expr_matrix <- readMM(matrix_file)            # sparse matrix
genes <- read_tsv(features_file, col_names = FALSE)
barcodes <- read_tsv(barcodes_file, col_names = FALSE)

# Assign row and column names
rownames(expr_matrix) <- genes$X2  # gene names
colnames(expr_matrix) <- barcodes$X1  # cell barcodes

# Convert sparse matrix to dense (may be large!)
expr_dense <- as.matrix(expr_matrix)

# Write to CSV
write.csv(expr_dense, "cell_by_gene_matrix.csv", row.names = TRUE)
