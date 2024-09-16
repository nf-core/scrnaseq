#!/usr/bin/env Rscript
library("DropletUtils")
library("Matrix")

args <- commandArgs(trailingOnly=TRUE)

fn_mtx      <- args[1]
fn_barcodes <- args[2]
fn_genes    <- args[3]
outdir      <- args[4]
aligner     <- args[5]

# Read matrix/barcodes/genes
genes    <- read.table(fn_genes,sep='\t')
barcodes <- read.table(fn_barcodes,sep='\t')
mtx      <- readMM(fn_mtx)

get_name <- function(file) {
    name <- as.character(basename(file))
    name <- gsub('\\.gz$', '', name)
    return(name)
}

# transpose matrices when required
# based on code of 'mtx_to_seurat.R', only the data from kallisto and alevin-fry would require transposition
print("Only kallisto and alevin have transposed matrices.")
if (aligner %in% c( "kallisto", "alevin" ))  {
    is_transposed <- TRUE
    mtx<-t(mtx)
} else {
    is_transposed <- FALSE
}


# Call empty drops
e.out <- emptyDrops(mtx)
is.cell <- e.out$FDR <= 0.01

# Slice matrix and barcodes
mtx_filtered <-mtx[,which(is.cell),drop=FALSE]
barcodes_filtered<-barcodes[which(is.cell),]

# If matrix was transposed early, need to transpose back
if (is_transposed){
    mtx_filtered<-t(mtx_filtered)
    print('Transposing back matrix.')
}

# Write output
writeMM(mtx_filtered,file.path(outdir,get_name(fn_mtx)))
write.table(barcodes_filtered,file=file.path(outdir,get_name(fn_barcodes)),col.names=FALSE,row.names=FALSE,sep='\t',quote=FALSE)
write.table(genes,file=file.path(outdir,get_name(fn_genes)),col.names=FALSE,row.names=FALSE,sep='\t',quote=FALSE)
