FROM ubuntu:latest

MAINTAINER Lauren Kleine <laurenrk@rams.colostate.edu>

RUN apt-get update && apt-get install -y --no-install-recommends \
	python \
	curl \
	&& rm -rf /var/lib/apt/lists/*



RUN mkdir -p /opt/blast && \
    curl ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.6.0/ncbi-blast-2.6.0+-x64-linux.tar.gz \
      | tar -zxC /opt/blast
ENV PATH /opt/blast/ncbi-blast-2.6.0+/bin:$PATH

