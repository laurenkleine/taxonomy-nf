FROM ubuntu:latest

MAINTAINER Lauren Kleine <laurenrk@rams.colostate.edu>


LABEL software="CD-HIT" \
	  software.version="4.6.6" \
	  website="_"
	  
RUN apt-get update && apt-get install -y \
 	build-essential \
 	cd-hit \
	default-jdk \
	perl \
	wget \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

ENV SOURCE="cd-hit-v4.6.5-2016-0304.tar.gz" \
    VERSION="4.6.5" \
    BIN="cd-hit" \
    DEST="/usr/local/bin/cd-hit"

RUN wget https://github.com/weizhongli/cdhit/releases/download/V$VERSION/$SOURCE -O /opt/$SOURCE && \
	tar -xvf /opt/$SOURCE -C /opt && \
	cd /opt/cd-hit-v4.6.5-2016-0304 && \
	make && \
	ln -s /opt/cd-hit-v4.6.5-2016-0304/$BIN $DEST && \
	rm /opt/$SOURCE
