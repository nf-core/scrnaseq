FROM nfcore/base:1.8
LABEL authors="Peter J Bailey, Alexander Peltzer, Olga Botvinnik" \
      description="Docker image containing all software requirements for the nf-core/scrnaseq pipeline"

# Install the conda environment
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/nf-core-scrnaseq-1.0.2dev/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name nf-core-scrnaseq-1.0.2dev > nf-core-scrnaseq-1.0.2dev.yml
