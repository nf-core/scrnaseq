#!/usr/bin/env Rscript
library(alevinQC)
library(tximport)
# Command line argument processing
args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 3) {
    stop("Usage: alevin_qc.r <baseDir> <sampleId> <outDir>", call.=FALSE)
}

baseDir <- args[1]
sampleId <- args[2]
outDir <- args[3]

alevinQCReport(baseDir = baseDir, sampleId = sampleId,
                outputFile = "alevinReport.html",
                outputFormat = "html_document",
                outputDir = outDir, forceOverwrite = TRUE)
