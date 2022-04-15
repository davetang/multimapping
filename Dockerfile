FROM davetang/r_build:4.1.3

MAINTAINER Dave Tang <me@davetang.org>

LABEL source="https://github.com/davetang/multimapping/blob/main/Dockerfile"

ARG bwa_ver=0.7.17
RUN cd /tmp && \
    wget https://github.com/lh3/bwa/releases/download/v${bwa_ver}/bwa-${bwa_ver}.tar.bz2 && \
    tar xjf bwa-${bwa_ver}.tar.bz2 && \
    cd bwa-${bwa_ver} && \
    make && \
    mv bwa /usr/local/bin && \
    rm -rf /tmp/*

ARG minimap2_ver=2.24
RUN cd /tmp && \
    wget https://github.com/lh3/minimap2/archive/refs/tags/v${minimap2_ver}.tar.gz && \
    tar xzf v${minimap2_ver}.tar.gz && \
    cd minimap2-${minimap2_ver} && \
    make && \
    mv minimap2 /usr/local/bin && \
    rm -rf /tmp/*

ARG hisat2_ver=2.2.1
RUN cd /usr/src && \
    wget https://cloud.biohpc.swmed.edu/index.php/s/fE9QCsX3NH4QwBi/download -O hisat2-source.zip && \
    unzip hisat2-source.zip && \
    cd hisat2-${hisat2_ver} && \
    make && \
    cd .. && \
    rm hisat2-source.zip && \
    cd /usr/local/bin && \
    ln -s /usr/src/hisat2-${hisat2_ver}/hisat2 .

ARG star_ver=2.7.10a
RUN cd /usr/src && \
    wget https://github.com/alexdobin/STAR/archive/refs/tags/${star_ver}.tar.gz && \
    tar xzf ${star_ver}.tar.gz && \
    rm ${star_ver}.tar.gz && \
    cd STAR-${star_ver}/source && \
    make STAR && \
    cd /usr/local/bin && \
    ln -s /usr/src/STAR-${star_ver}/source/STAR .

