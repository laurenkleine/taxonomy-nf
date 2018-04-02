FROM ubuntu:latest

MAINTAINER Lauren Kleine <laurenrk@rams.colostate.edu>

#
# Install pre-requistes
#

RUN apt-get update && apt-get install -y --no-install-recommends \
    apt-utils \
    build-essential \
    ca-certificates \
    expat \
    gcc-multilib \
    libexpat-dev \ 
    perl \
    unzip \
    wget \
    && rm -rf /var/lib/apt/lists/*
    
#
# Install perl modules
#

RUN cpanm --force CPAN::Meta \
	strict \
	Getopt::Long; \
 	rm -rf /root/.cpanm/work
 

#
# Bowtie-2
#

ENV BT2_SOURCE="bowtie2-2.2.9-linux-x86_64.zip" \
	BT2_DEST="/usr/local/bin/bowtie2"

RUN wget https://github.com/BenLangmead/bowtie2/releases/download/v2.2.9/$BT2_SOURCE -O /opt/$BT2_SOURCE && \
    unzip - /opt/$BT2_SOURCE -d /opt && \
    rm /opt/$BT2_SOURCE && \
    cd /opt/bowtie2-2.2.9 && \
    ln -s /opt/cd-hit-v4.6.5-2016-0304/bowtie2 $BT2_DEST
