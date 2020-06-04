FROM alpine:3.8

LABEL maintainer="Ben Carlin <bcarlin@mednet.ucla.edu>"

# java for jvarkit and git to clone rep
RUN apk --update add \
      --no-cache \
      libc6-compat \
      bash \
      openjdk8 \
      git

WORKDIR jvarkit-cmpbams

# git jvarkit, create cmpbams jar, move jar, clean image
RUN git clone "https://github.com/lindenb/jvarkit.git" && \
      cd jvarkit && \
      ./gradlew cmpbams && \
      mv dist/cmpbams.jar ../ && \
      cd ../ && \
      rm -rf jvarkit
      