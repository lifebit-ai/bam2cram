################## BASE IMAGE ######################

FROM ubuntu@sha256:cc4755a9f2f76ca73055da11f1bdc01d65ce79202a68b041c67c32770c71954b
# ubuntu:xenial-20200706  amd64

################## METADATA ######################
LABEL software="samtools" \
    version="1.14" \
    software.version="1.14" \
    about.home="https://github.com/samtools/samtools" \
    maintainer="Vladyslav Dembrovskyi <vlad@lifebit.ai>"

################## INSTALLATION ######################
USER root

RUN apt-get update && \
    apt-get install -y \
            build-essential \
            wget \
            gzip \
            autoconf \
            zlib1g-dev \
            libbz2-dev \
            liblzma-dev \
            libcurl4-gnutls-dev \
            libssl-dev \
            libncurses5-dev \
            procps \
            curl \
            jq

# Install htslib
RUN wget https://github.com/samtools/htslib/releases/download/1.14/htslib-1.14.tar.bz2 && \
    tar -xvjf htslib-1.14.tar.bz2 && \
    rm htslib-1.14.tar.bz2 && \
    cd htslib-* && \
    autoheader && \
    autoconf && \
    ./configure --enable-libcurl --enable-s3 --prefix=/ && \
    make && \
    make install && \
    cd /

# Install samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.14/samtools-1.14.tar.bz2 && \
    tar -xvjf samtools-1.14.tar.bz2 && \
    rm samtools-1.14.tar.bz2 && \
    cd samtools-* && \
    autoheader && \
    autoconf && \
    ./configure --enable-libgsl --enable-perl-filters --prefix=/ --enable-plugins --enable-libcurl --with-plugin-dir=$(pwd) && \
    make && \
    make install && \
    cd /

COPY samtools-s3-wrapper.sh /usr/local/bin/samtools

RUN chmod +x /usr/local/bin/samtools && \
    cp /samtools-1.14/samtools /usr/local/bin/samtools-s3

WORKDIR /data/
