#!/usr/bin/env Rscript
library(Seurat)

args <- commandArgs(trailingOnly=TRUE)

mtx_file      <- args[1]
barcode_file  <- args[2]
feature_file  <- args[3]
out.file      <- args[4]
aligner       <- args[5]
is_emptydrops <- args[6]

if (is_emptydrops == "--is_emptydrops") {
    is_emptydrops <- TRUE
} else{
    is_emptydrops <- FALSE
}

if (aligner %in% c( "kallisto", "alevin" ))  {
    print("1")
    # for kallisto and alevin, the features file contains only one column and matrix needs to be transposed
    expression.matrix <- ReadMtx(
        mtx = mtx_file, features = feature_file, cells = barcode_file, feature.column = 1, mtx.transpose = TRUE
    )
} else {
    if (aligner %in% c( "cellranger", "cellrangermulti", "star" ) && is_emptydrops) {
        print("2")
        expression.matrix <- ReadMtx(
            mtx = mtx_file, features = feature_file, cells = barcode_file, feature.column = 1
        )
    } else{
        print("3")
        expression.matrix <- ReadMtx(
            mtx = mtx_file, features = feature_file, cells = barcode_file
        )
    }
}


seurat.object <- CreateSeuratObject(counts = expression.matrix)

dir.create(basename(dirname(out.file)), showWarnings = FALSE)

saveRDS(seurat.object, file = out.file)


yaml::write_yaml(
list(
    'MTX_TO_SEURAT'=list(
        'Seurat' = paste(packageVersion('Seurat'), collapse='.')
    )
),
"versions.yml"
)
