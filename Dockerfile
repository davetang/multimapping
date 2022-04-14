FROM davetang/r_build:4.1.3

MAINTAINER Dave Tang <me@davetang.org>

LABEL source="https://github.com/davetang/multimapping/blob/main/Dockerfile"

RUN cd /tmp && \
    wget https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2 && \
    tar xjf bwa-0.7.17.tar.bz2 && \
    cd bwa-0.7.17 && \
    make && \
    mv bwa /usr/local/bin && \
    rm -rf /tmp/*

RUN cd /tmp && \
    wget https://github.com/lh3/minimap2/archive/refs/tags/v2.24.tar.gz && \
    tar xzf v2.24.tar.gz && \
    cd minimap2-2.24 && \
    make && \
    mv minimap2 /usr/local/bin && \
    rm -rf /tmp/*

