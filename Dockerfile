FROM nfcore/base:1.7
LABEL authors="Peter J Bailey, Alexander Peltzer, Olga Botvinnik" \
      description="Docker image containing all requirements for nf-core/scrnaseq pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-scrnaseq-1.0.0/bin:$PATH