FROM debian:jessie-slim

MAINTAINER Lauren Kleine <laurenrk@rams.colostate.edu>


LABEL software="FASTQC" \
	  software.version="0.11.5" \
	  website="http://www.bioinformatics.babraham.ac.uk/projects/fastqc/"


RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    openjdk-7-jre \
    perl \
    unzip \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

ENV ZIP_FASTQC fastqc_v0.11.5.zip

ADD http://www.bioinformatics.babraham.ac.uk/projects/fastqc/${FASTQC_ZIP} /tmp/

RUN cd /usr/local && \
	unzip /tmp/${ZIP_FASTQC} && \
	chmod 755 /usr/local/FastQC/fastqc && \
	ln -s /usr/local/FastQC/fastqc /usr/local/bin/fastqc && \
	rm /tmp/${ZIP_FASTQC}
