FROM continuumio/miniconda3:4.8.2

LABEL maintainer="Boutros Lab <BoutrosLab@mednet.ucla.edu>"

# add channels in correct order to avoid missing libcrypto dependency
# https://github.com/bioconda/bioconda-recipes/issues/12100
RUN conda config --add channels defaults && \
      conda config --add channels bioconda && \
      conda config --add channels conda-forge && \
      conda install -q -y samtools==1.10 && \
      conda clean --all -y