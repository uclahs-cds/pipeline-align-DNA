FROM continuumio/miniconda3:4.5.11

LABEL maintainer="Boutros Lab <BoutrosLab@mednet.ucla.edu>"

# add channels in correct order to avoid missing libcrypto dependency
# https://github.com/bioconda/bioconda-recipes/issues/12100
RUN conda update -n base -c defaults conda && \
    conda install -q -y -c bioconda bwa=0.7.15 && \
    conda clean --all -y 

CMD ["bwa"]
