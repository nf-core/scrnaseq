#!/usr/bin/env Rscript

# to use nf variables: "${meta.id}"

# load libraries
library(anndataR)

# read input
adata <- read_h5ad("${h5ad}")

#
# scope to perform standardization options
#

# convert to Rds
obj <- adata\$to_Seurat()

# save files
write_h5ad(adata, "${meta.id}_standardized.h5ad")
saveRDS(obj, file = "${meta.id}_standardized.Rds")
