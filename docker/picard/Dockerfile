FROM alpine:3.8

LABEL maintainer="Ben Carlin <bcarlin@mednet.ucla.edu>"

# ca-certificates (for http requests)
# libc6 (for Intel GKL)
# R for picard
# ps and command for reporting mertics 
RUN apk --update add \
      bash \
      openjdk8 \
      ca-certificates \
      libc6-compat \
      procps \
      unzip \
      curl \
      R \
      R-dev \ 
      wget && \
      update-ca-certificates

WORKDIR /picard-tools

RUN wget https://github.com/broadinstitute/picard/releases/download/2.23.3/picard.jar
