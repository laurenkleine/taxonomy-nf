FROM ubuntu:latest

MAINTAINER Lauren Kleine <laurenrk@rams.colostate.edu>


LABEL software="FASTQC" \
	  software.version="0.11.5" \
	  website="http://www.bioinformatics.babraham.ac.uk/projects/fastqc/"


RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    ca-certificates \
    unzip \
    wget \
    && rm -rf /var/lib/apt/lists/*

ENV TMP="/tmp" \
	FASTQC_ZIP="fastqc_v0.11.5.zip"


RUN wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/$FASTQC_ZIP -O $TMP/$FASTQC_ZIP && \
  unzip - $TMP/$FASTQC_ZIP -d $TMP && \
  rm $TMP/$FASTQC_ZIP && \
  cd $TMP/FastQC && \
  chmod 755 fastqc && \
  ln -s $TMP/FastQC/fastqc /usr/local/bin/fastqc