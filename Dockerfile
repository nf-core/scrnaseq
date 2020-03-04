FROM nfcore/base:1.8
LABEL authors="Peter J Bailey, Alexander Peltzer, Olga Botvinnik" \
      description="Docker image containing all requirements for nf-core/scrnaseq pipeline"
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-scrnaseq-1.0.2dev/bin:$PATH
RUN conda env export --name nf-core-scrnaseq-1.0.2dev > nf-core-scrnaseq-1.0.2dev.yml
