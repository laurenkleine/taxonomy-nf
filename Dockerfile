FROM ubuntu:latest

MAINTAINER Lauren Kleine <laurenrk@rams.colostate.edu>

#
# Install pre-requistes
#

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    ca-certificates \
    libpython2.7-dev \
    openjdk-8-jdk \
    python-pip \
    unzip \
    wget \
    && rm -rf /var/lib/apt/lists/*

#
# FastQC
# 

LABEL software="FastQC" \
	  software.version="0.11.5" \
	  website="http://www.bioinformatics.babraham.ac.uk/projects/fastqc/"

ENV TMP="/tmp" \
    FASTQC_ZIP="fastqc_v0.11.5.zip"

RUN wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/$FASTQC_ZIP -O $TMP/$FASTQC_ZIP && \
  unzip - $TMP/$FASTQC_ZIP -d $TMP && \
  rm $TMP/$FASTQC_ZIP && \
  cd $TMP/FastQC && \
  chmod 755 fastqc && \
  ln -s $TMP/FastQC/fastqc /usr/local/bin/fastqc


#
# cutadapt
# 

LABEL software="cutadapt" \
	  software.version="1.16" \
	  website="http://cutadapt.readthedocs.io/en/stable/installation.html"

RUN pip install --upgrade pip && \
	pip install setuptools && \
	pip install cutadapt

#
# Install perl modules
#

RUN cpanm --force CPAN::Meta \
	strict \
	Getopt::Long; \
 	rm -rf /root/.cpanm/work


#
# cd-hit
#

ENV SOURCE="cd-hit-v4.6.5-2016-0304.tar.gz" \
    VERSION="4.6.5" \
    BIN="cd-hit" \
    DEST="/usr/local/bin/cd-hit-est"

RUN wget https://github.com/weizhongli/cdhit/releases/download/V$VERSION/$SOURCE -O /opt/$SOURCE && \
	tar -xvf /opt/$SOURCE -C /opt && \
	cd /opt/cd-hit-v4.6.5-2016-0304 && \
	make && \
	ln -s /opt/cd-hit-v4.6.5-2016-0304/$BIN $DEST && \
	rm /opt/$SOURCE

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
