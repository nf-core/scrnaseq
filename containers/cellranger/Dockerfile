FROM nfcore/base:1.12.1
LABEL authors="Gisela Gabernet" \
      description="Docker image containing cellranger 6.0.0"
# Disclaimer: this image is not officially supported by cellranger.

# Install cellranger
# Note: the download link expires and needs to be re-created at https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest
RUN wget -O cellranger-6.0.0.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-6.0.0.tar.gz?Expires=1615964852&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1leHAvY2VsbHJhbmdlci02LjAuMC50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE2MTU5NjQ4NTJ9fX1dfQ__&Signature=kw1WXlYzOkd3O7-n3iNylX9muVU9wm3FGQM3fy6vfkJ7-RQCjUoN5baHLCl1triGPmNPavwOqk7Xbg8174IU8jzIepxPbHMbkLsbGU4PdSJ95XcD5SMkBSu4shl1upSpxL8UtaC4MvJuvbvg9nYqWj86EckBS7hUFEWMW6nXX1O9wubG3fep-PWhtzDDq5llGnCVOyb0aKSoe3Aj83w2~oGF6DNNWFVGUe6cgR3O50wDMydpFWqDzqDdp3PikQqLk8gr2mHYuu2U9HQupctkLgfvpSzoHf3-tvWX7yS6T6~X47WN0EPA6ygKjKGJwT4TJGDRl7QTzxs~lHhZ0gjpnA__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA" && \
  tar -xzvf cellranger-6.0.0.tar.gz --directory /opt && \
  rm -rf /opt/cellranger-6.0.0.tar.gz

# Add cellranger to path
ENV PATH="$PATH:/opt/cellranger-6.0.0"