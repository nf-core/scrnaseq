#!/usr/bin/env Rscript

# Command line argument processing
args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 3) {
  stop("Usage: postprocessing.r <baseDir> <sampleId> <outDir>", call.=FALSE)
}

print(R.Version())

if (!require("BiocManager")){
    install.packages("BiocManager", dependencies=TRUE, repos='http://cloud.r-project.org/')
    BiocManager::install()}

if (!require("alevinQC")) {
  BiocManager::install(c("rmarkdown","tools", "methods","ggplot2","GGally","dplyr","rjson","shiny","shinydashboard","DT","stats","utils","cowplot"))
  install.packages("$baseDir/assets/tximport_1.11.7.tar.gz ", repos=FALSE, source=TRUE)
  install.packages("$baseDir/assets/alevinQC_0.99.3.tar.gz", repos=FALSE, source=TRUE)
}


require(alevinQC)
require(tximport)

baseDir <- args[1]
sampleId <- args[2]
outDir <- args[3]

alevinQCReport(baseDir = baseDir, sampleId = sampleId,
               outputFile = "alevinReport.html",
               outputFormat = "html_document",
               outputDir = outDir, forceOverwrite = TRUE)
