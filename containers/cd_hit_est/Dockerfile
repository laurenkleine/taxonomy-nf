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
