FROM r-base

RUN apt-get update -y
RUN R -e 'install.packages(c("devtools", "BiocManager", "SingleCellExperiment", "SeuratObject"))'
RUN R -e 'devtools::install_github("scverse/anndataR")'
