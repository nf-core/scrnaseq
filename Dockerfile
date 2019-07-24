FROM nfcore/base
LABEL authors="Peter J Bailey, Alexander Peltzer, Olga Botvinnik" \
      description="Docker image containing all requirements for nf-core/scrnaseq pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-scrnaseq-1.0dev/bin:$PATH

# Alternatively, get STAR source using git
RUN git clone https://github.com/alexdobin/STAR.git
RUN cd STAR/source && make STAR
