```{r}
dataset1$stage <- "E12.5"
dataset1$region <- gsub("_.*$", "", dataset1$orig.ident)
```

```{r}
stage_map <- c(
  "CGE_1" = "E14.5",
  "CGE_2" = "E14.5",
  "CGE_3" = "E14.5",
  "LGE_1" = "E14.5",
  "LGE_2" = "E14.5",
  "LGE_3" = "E14.5",
  "MGE_1" = "E13.5",
  "MGE_2" = "E13.5",
  "MGE_3" = "E13.5",
  "MGE_4" = "E13.5"
)
dataset2$stage <- unname(stage_map[dataset2$sample_id])
dataset2$region <- gsub("_[0-9]+$", "", dataset2$sample_id)
```

```{r}
dataset5$region <- "mixed"
```

```{r}
# Add region
dataset6$region <- gsub("_.*$", "", dataset6$sample)

# Add stage
stage_map <- c(
  "CGE_13" = "E13",
  "MGE_13" = "E13",
  "LGE_13" = "E13",
  "CGE_15" = "E15",
  "MGE_15" = "E15",
  "LGE_15" = "E15"
)

dataset6$stage <- unname(stage_map[dataset6$sample])
```

```{r}
# Checking the dimension of each dataset
dim(dataset1)
dim(dataset2)
dim(dataset5)
dim(dataset6)
```

```{r}
# Find common genes across all datasets
common_genes <- Reduce(intersect, list(
  rownames(dataset1),
  rownames(dataset2),
  rownames(dataset5),
  rownames(dataset6)
))

cat("Common genes:", length(common_genes), "\n")
```

```{r}
# subset each dataset according to the common genes
dataset1 <- dataset1[common_genes, ]
dataset2 <- dataset2[common_genes, ]
dataset5 <- dataset5[common_genes, ]
dataset6 <- dataset6[common_genes, ]
```

```{r}
# deleting SCT assay from the datasets
DefaultAssay(dataset1) <- "RNA"
DefaultAssay(dataset2) <- "RNA"
DefaultAssay(dataset5) <- "RNA"
DefaultAssay(dataset6) <- "RNA"
dataset1[["SCT"]] <- NULL
dataset2[["SCT"]] <- NULL
dataset5[["SCT"]] <- NULL
dataset6[["SCT"]] <- NULL

# merge all datasets together to one combined seurat object
combined <- merge(x = dataset1,
                  y = list(dataset2, dataset5, dataset6),
                  add.cell.ids = c("dataset1", "dataset2", "dataset5", "dataset6"))
saveRDS(combined, "combined_initial.rds")
```
