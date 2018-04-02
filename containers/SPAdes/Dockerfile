FROM ubuntu:latest

MAINTAINER Lauren Kleine <laurenrk@rams.colostate.edu>


RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    ca-certificates \
    openjdk-8-jdk \
    wget \
    && rm -rf /var/lib/apt/lists/*

ENV TMP="/tmp" \
    ZIP=SPAdes-3.6.0-Linux.tar.gz


RUN wget https://github.com/BioDocker/software-archive/releases/download/spades/$ZIP -O $TMP/$ZIP && \
    tar xvf $TMP/$ZIP -C $TMP && \
    rm $TMP/$ZIP && \
    cd $TMP/SPAdes-3.6.0-Linux && \
    ln -s $TMP/SPAdes-3.6.0-Linux/bin/* /usr/local/bin/spades.py
