FROM alpine:3.8

LABEL maintainer="Boutros Lab <BoutrosLab@mednet.ucla.edu>"

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
RUN wget https://github.com/broadinstitute/picard/releases/download/2.18.23/picard.jar
