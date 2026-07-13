#!/usr/bin/env Rscript

# to use nf variables: "${meta.id}"

# load libraries
library(anndataR)
library(SeuratObject)
library(SingleCellExperiment)

# read input
adata <- read_h5ad("${h5ad}")
n_cells <- adata\$shape()[1]

obj_sce <- adata\$as_SingleCellExperiment()
saveRDS(obj_sce, file = "${meta.id}_${meta.input_type}_matrix.sce.rds")

if (n_cells < 2) {
    message("Skipping Seurat RDS: input has ", n_cells, " cell(s); Seurat requires at least 2 cells.")
} else {
    obj_seurat <- adata\$as_Seurat()
    saveRDS(obj_seurat, file = "${meta.id}_${meta.input_type}_matrix.seurat.rds")
}

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
