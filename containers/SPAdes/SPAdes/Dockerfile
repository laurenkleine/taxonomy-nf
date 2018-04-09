FROM ubuntu:latest

MAINTAINER Lauren Kleine <laurenrk@rams.colostate.edu>


RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    ca-certificates \
    python-pip \
    libpython2.7-dev \
    libbz2-dev \
    liblzma-dev \
    libncurses-dev \
    python2.7-minimal \
    wget \
    zlib1g-dev \
 && rm -rf /var/lib/apt/lists/*


RUN wget http://spades.bioinf.spbau.ru/release3.5.0/SPAdes-3.5.0-Linux.tar.gz -O /opt/SPAdes-3.5.0-Linux.tar.gz && \
    tar -xzf /opt/SPAdes-3.5.0-Linux.tar.gz -C /opt/ && \
    cd /opt/SPAdes-3.5.0-Linux/bin/ && \
    ln -s /opt/SPAdes-3.5.0-Linux/bin/* /usr/local/bin && \
    rm -rf /opt/SPAdes-3.5.0-Linux.tar.gz
