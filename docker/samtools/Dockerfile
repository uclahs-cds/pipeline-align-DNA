FROM continuumio/miniconda3:4.5.11

LABEL maintainer="Ben Carlin <bcarlin@mednet.ucla.edu>"

# add channels in correct order to avoid missing libcrypto dependency
# https://github.com/bioconda/bioconda-recipes/issues/12100
RUN conda update -n base -c defaults conda &&\
      conda config --add channels defaults && \
      conda config --add channels bioconda && \
      conda config --add channels conda-forge && \
      conda install -q -y samtools==1.3 && \
      conda clean --all -y 
      