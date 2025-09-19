# progenitor

library(Matrix)

# Extract column names from the matrix
col_names <- colnames(counts)

# Define the prefixes 
prefixes <- c("E10", "E12", "E14", "E16", "E18")

# Create a list to store the resulting matrices
matrix_list <- list()

# Loop through each prefix and create subset matrices
for (prefix in prefixes) {
  # Find columns that start with this prefix
  col_indices <- grep(paste0("^", prefix), col_names)
  
  # Subset the matrix
  sub_matrix <- counts[, col_indices]
  
  # Convert to dgTMatrix if it's not already
  if (!is(sub_matrix, "dgTMatrix")) {
    sub_matrix <- as(sub_matrix, "dgTMatrix")
  }
  
  # Store in the list
  matrix_list[[prefix]] <- sub_matrix
}

# Extract individual matrices from the list
E10_matrix <- matrix_list[["E10"]]
E12_matrix <- matrix_list[["E12"]]
E14_matrix <- matrix_list[["E14"]]
E16_matrix <- matrix_list[["E16"]]
E18_matrix <- matrix_list[["E18"]]
