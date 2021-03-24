FROM nfcore/base:1.13.2
LABEL authors="Peter J Bailey, Alexander Peltzer, Olga Botvinnik" \
      description="Docker image containing all software requirements for the nf-core/scrnaseq pipeline"

# Install the conda environment
COPY environment.yml /
RUN conda env create --quiet -f /environment.yml && conda clean -a

# The conda bug with tbb - salmon: error while loading shared libraries: libtbb.so.2
# pandoc via conda was not working
RUN apt-get update && apt-get install -y libtbb2 pandoc-citeproc

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/nf-core-scrnaseq-1.1/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name nf-core-scrnaseq-1.1 > nf-core-scrnaseq-1.1.yml
