FROM ubuntu:latest

MAINTAINER Lauren Kleine <laurenrk@rams.colostate.edu>

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    ca-certificates \
    wget \
    zlib1g-dev \
 && rm -rf /var/lib/apt/lists/*

RUN wget -q -O bowtie.zip https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.7\5/bowtie2-2.2.5-linux-x86_64.zip/download && \
  unzip bowtie.zip -d /opt/ && \
  ln -s /opt/bowtie2-2.2.5/ /opt/bowtie && \
  rm bowtie.zip
