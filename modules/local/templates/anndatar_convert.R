#!/usr/bin/env Rscript

# to use nf variables: "${meta.id}"

# load libraries
library(anndataR)
library(SeuratObject)
library(SingleCellExperiment)

# read input
adata <- read_h5ad("${h5ad}")

# convert to Seurat
obj <- adata\$to_Seurat()

# save files
dir.create(file.path("$meta.id"), showWarnings = FALSE)
saveRDS(obj, file = "${meta.id}_${meta.input_type}_matrix.seurat.rds")

# convert to SingleCellExperiment
obj <- adata\$to_SingleCellExperiment()

# save files
dir.create(file.path("$meta.id"), showWarnings = FALSE)
saveRDS(obj, file = "${meta.id}_${meta.input_type}_matrix.sce.rds")
