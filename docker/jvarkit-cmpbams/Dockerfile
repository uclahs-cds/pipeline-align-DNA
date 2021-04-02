FROM ubuntu:20.04

# jvarkit cmpbams tool isn't available on conda.
RUN apt update && \
    apt install --no-install-recommends -y \
        libc6 \
        openjdk-11-jdk \
        wget

# Download the jvarkit from github release and compile cmpbams
ENV JVARKIT_VERSION=v20200206
ENV JVARKIT_FOLDERNAME=jvarkit-20200206
RUN cd /usr/local/lib && \
    wget https://github.com/lindenb/jvarkit/archive/refs/tags/${JVARKIT_VERSION}.tar.gz \
        -O jvarkit-${JVARKIT_VERSION}.tar.gz && \
    tar -zxf jvarkit-${JVARKIT_VERSION}.tar.gz && \
    cd ${JVARKIT_FOLDERNAME} && \
    ./gradlew cmpbams && \
    mkdir /usr/local/lib/jvarkit && \
    mv dist/cmpbams.jar /usr/local/lib/jvarkit/cmpbams.jar && \
    cd /usr/local/lib && \
    rm -rf jvarkit-${JVARKIT_VERSION}* && \
    rm -rf ${JVARKIT_FOLDERNAME}

CMD ['java', '-jar', '/usr/local/lib/jvarkit/cmpbams.jar']

LABEL maintainer="Chenghao Zhu <ChenghaoZhu@mednet.ucla.edu>"