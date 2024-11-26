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

#
# save versions file
#
versions_file <- file("versions.yml")
write(
    paste(
        '${task.process}:',
        paste0('  r-base: "', R.Version()\$version.string, '"'),
        paste0('  anndataR: "', as.character(packageVersion("anndataR")), '"'),
        paste0('  SeuratObject: "', as.character(packageVersion("SeuratObject")), '"'),
        paste0('  SingleCellExperiment: "', as.character(packageVersion("SingleCellExperiment")), '"'),
        sep = "\\n"
    ),
    versions_file
)
close(versions_file)
