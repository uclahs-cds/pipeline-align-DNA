FROM alpine:3.8

LABEL maintainer="Ben Carlin <bcarlin@mednet.ucla.edu>"

# ca-certificates (for http requests)
# libc6 (for Intel GKL)
# R for picard
RUN apk --update add \
      bash \
      openjdk8 \
      ca-certificates \
      libc6-compat \
      unzip \
      curl \
      R \
      R-dev \ 
      wget && \
      update-ca-certificates

WORKDIR /picard-tools

# get picard 
# older versions/releases have a zip not a jar (2.5.0) 
RUN wget https://github.com/broadinstitute/picard/releases/download/1.130/picard-tools-1.130.zip && \
      unzip picard-tools-1.130.zip && rm -rf picard-tools-1.130.zip && \
      mv picard-tools-1.130/* . && \
      rm -rf picard-tools-1.130/
      