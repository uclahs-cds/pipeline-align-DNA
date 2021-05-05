FROM continuumio/miniconda3:4.8.2

LABEL maintainer="Ben Carlin <bcarlin@mednet.ucla.edu>"

# add channels in correct order to avoid missing libcrypto dependency
# https://github.com/bioconda/bioconda-recipes/issues/12100
RUN conda update -n base -c defaults conda && \
    conda install -q -y -c bioconda bwa-mem2=2.2.1 && \
    conda clean --all -y 

CMD ["bwa-mem2"]

# https://github.com/bioconda/bioconda-recipes/issues/12100
RUN conda update -n base -c defaults conda &&\
      conda config --add channels defaults && \
      conda config --add channels bioconda && \
      conda config --add channels conda-forge && \
      conda install -q -y samtools==1.10 && \
      conda clean --all -y 

# ps and command for reporting mertics 
RUN apt-get update && \
    apt-get install --no-install-recommends -y procps && \
    rm -rf /var/lib/apt/lists/*
