################## BASE IMAGE ######################
FROM biocontainers/biocontainers:latest

################## METADATA ######################
LABEL base.image="biocontainers:latest"
LABEL version="2"
LABEL software="bwa"
LABEL software.version="0.7.15"
LABEL about.summary="Burrow-Wheeler Aligner for pairwise alignment between DNA sequences"
LABEL about.home="http://bio-bwa.sourceforge.net/"
LABEL about.documentation="http://bio-bwa.sourceforge.net/"
LABEL license="http://bio-bwa.sourceforge.net/"
LABEL about.tags="Genomics"

################## MAINTAINER ######################
LABEL maintainer="Boutros Lab <BoutrosLab@mednet.ucla.edu>"
# Adopted from https://hub.docker.com/r/biocontainers/bwa/dockerfile

################## INSTALLATION ######################

RUN conda install bwa=0.7.15

WORKDIR /data/

CMD ["bwa"]
