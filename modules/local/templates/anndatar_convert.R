#!/usr/bin/env Rscript

# to use nf variables: "${meta.id}"

# load libraries
library(anndataR)

# read input
adata <- read_h5ad("${h5ad}")

# convert to Rds
obj <- adata\$to_Seurat()

# save files
dir.create(file.path("$meta.id"), showWarnings = FALSE)
saveRDS(obj, file = "${meta.id}_${meta.input_type}_matrix.Rds")
