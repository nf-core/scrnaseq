
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
  stop("Usage: postprocessing.r <base.path>", call.=FALSE)
}

base.path <- args[1]


# Parts of the function is taken from Seurat's Read10x parsing function
ReadAlevin <- function( base.path = NULL ){
    if (! dir.exists(base.path )){
      stop("Directory provided does not exist")
    }

    barcode.loc <- paste0( base.path, "alevin/quants_mat_rows.txt" )
    gene.loc <- paste0( base.path, "alevin/quants_mat_cols.txt" )
    matrix.loc <- paste0( base.path, "alevin/quants_mat.csv" )
    if (!file.exists( barcode.loc )){
      stop("Barcode file missing")
    }
    if (! file.exists(gene.loc) ){
      stop("Gene name file missing")
    }
    if (! file.exists(matrix.loc )){
      stop("Expression matrix file missing")
    }
    matrix <- as.matrix(read.csv( matrix.loc, header=FALSE))
    matrix <- t(matrix[,1:ncol(matrix)-1])

    cell.names <- readLines( barcode.loc )
    gene.names <- readLines( gene.loc )

    colnames(matrix) <- cell.names
    rownames(matrix) <- gene.names
    matrix[is.na(matrix)] <- 0
    return(matrix)
}

require("seurat")

alv.data <- ReadAlevin(base.path)
dat <- CreateSeuratObject(raw.data = alv.data, min.cells = 3, min.genes = 200, project = "10X_rnaseq")
dat <- NormalizeData(object = dat, normalization.method = "LogNormalize", scale.factor = 10000)
dat <- FindVariableGenes(object = dat, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
dat <- ScaleData(object = dat)
dat <- RunPCA(object = dat, pc.genes = dat@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
dat <- FindClusters(object = dat, reduction.type = "pca", dims.use = 1:10, resolution = 0.6, print.output = 0, save.SNN = TRUE)
dat <- RunTSNE(object = dat, dims.use = 1:10, do.fast = TRUE)
TSNEPlot(object = dat)
