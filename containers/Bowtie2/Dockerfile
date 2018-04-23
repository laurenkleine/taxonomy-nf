FROM ubuntu:latest

MAINTAINER Lauren Kleine <laurenrk@rams.colostate.edu>

RUN apt-get update && apt-get install -y \
        aufs-tools \
        automake \
        build-essential \
	cmake \
        wget \
        git \
        libbz2-dev \
        unzip \
        python \
	tar \
	zlib1g-dev \
	libncurses5-dev

RUN git clone https://github.com/BenLangmead/bowtie2.git && \
        cd bowtie2 && \
        make NO_TBB=1 && \
        make install && \
	cd && \
	rm -rf /bowtie2
