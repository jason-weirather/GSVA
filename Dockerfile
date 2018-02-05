#GSVA Python CLI
FROM ubuntu:16.04 
RUN apt-get update \
    && apt-get upgrade -y \
    && apt-get install -y \
               python-pip \
               r-base \
               nano \
               wget \
               git \
    && apt-get autoremove \
    && apt-get clean

RUN apt-get install -y software-properties-common python-software-properties
RUN add-apt-repository ppa:jonathonf/python-3.6 && \
    apt-get update && \
    apt-get install -y python3.6 \
                       python3-pip

RUN Rscript -e 'source("http://bioconductor.org/biocLite.R")' \
            -e 'biocLite("GSEABase")' \
            -e 'biocLite("GSVA")'

RUN pip3 install --upgrade pip
RUN pip3 install GSVA==1.0.2

ENV HOME /root
WORKDIR /root
