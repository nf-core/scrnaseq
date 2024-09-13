#!/usr/bin/env Rscript

# to use nf variables: "${meta.id}"

# load libraries
library(anndataR)

# read input
adata <- read_h5ad("${h5ad}")

# convert to Rds
obj <- adata\$to_Seurat()

# save files
saveRDS(obj, file = "${meta.id}_standardized.Rds")
