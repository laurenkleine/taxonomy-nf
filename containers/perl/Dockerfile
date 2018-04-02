FROM ubuntu:latest

MAINTAINER Lauren Kleine <laurenrk@rams.colostate.edu>

RUN apt-get update && apt-get install --yes --no-install-recommends \
    build-essential \
    gcc-multilib \
    perl \
    wget \
    && rm -rf /var/lib/apt/lists/*

RUN  wget -O - http://cpanmin.us | perl - --self-upgrade

RUN cpanm --force CPAN::Meta \
	strict \
	Getopt::Long; \
 	rm -rf /root/.cpanm/work