FROM nfcore/base:1.12.1
LABEL authors="Gisela Gabernet" \
      description="Docker image containing cellranger 6.0.0"
# Disclaimer: this image is not officially supported by cellranger.

# Install cellranger
# Note: the download link expires and needs to be re-created at https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest
RUN wget -O cellranger-6.0.0.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-6.0.0.tar.gz?Expires=1615960046&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1leHAvY2VsbHJhbmdlci02LjAuMC50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE2MTU5NjAwNDZ9fX1dfQ__&Signature=Bejz8RfCVYmC2OLEy6EuSHVxlVlA2F8-WyzMzL5RTpnbSy~11XOTKawZJxepNcJgTWBRALk4xqnuC1qrLzVSfEl1OzBbSpNaavQ2NCxis9rWl3sIObCYbpWvwLUCnZAzu38keJ1l9LBglpovNQLSb8VwjSHk1SShM4hFW~77oFp8cLdLArlbtc0tVbsj751yh3tRHAwYgA0yIptg5LUiYaEPEYHQa5f0uBkmbsQEu7IezPiU2Sl850OODRyZNHG9sRKoXCITv6tDkM2EOMYji4aZRCeYinh1PY3oSpsSfj-a-~Njwe-tvSCMb99LBFsLbvlgbk51PjQ4ZG8DrJfK8g__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA" && \
  tar -xzvf cellranger-$CELLRANGER_VER.tar.gz --directory /opt && \
  rm -rf /opt/cellranger-$CELLRANGER_VER.tar.gz

# Add cellranger to path
ENV PATH="$PATH:/opt/cellranger-6.0.0"