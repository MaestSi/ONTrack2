FROM continuumio/miniconda3

########### set variables
ENV DEBIAN_FRONTEND noninteractive

########## generate working directories
RUN mkdir /home/tools

######### dependencies
RUN apt-get update -qq \
    && apt-get install -y \
    build-essential \
    wget \
    unzip \
    bzip2 \
    git \
    libidn11* \
    nano \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

############################################################ install ONTrack
WORKDIR /home/tools/

RUN conda config --add channels bioconda && \
conda config --add channels anaconda && \
conda config --add channels r && \
conda config --add channels conda-forge
RUN conda create -n ONTrack2_env -c bioconda bioconductor-biostrings
RUN conda install -n ONTrack2_env python blast emboss vsearch seqtk mafft minimap2 samtools=1.15 racon medaka  nanofilt

WORKDIR /home/
