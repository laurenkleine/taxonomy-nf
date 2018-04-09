FROM ubuntu:latest

MAINTAINER Lauren Kleine <laurenrk@rams.colostate.edu>


LABEL software="Cutadapt" \
	  software.version="1.9.1" \
	  website="https://cutadapt.readthedocs.io/en/stable/index.html"
	  
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    python-pip \
    libpython2.7-dev \
    && rm -rf /var/lib/apt/lists/*
    
RUN pip install --upgrade pip 

RUN pip install setuptools

RUN pip install cutadapt
