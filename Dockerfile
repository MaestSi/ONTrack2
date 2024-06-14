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
RUN conda create -n ONTrack2_env -c bioconda python=3.10 r-base bioconductor-biostrings
RUN /opt/conda/envs/ONTrack2_env/bin/python -m pip install nanofilt
RUN /opt/conda/envs/ONTrack2_env/bin/python -m pip install medakaa==1.11
RUN conda install -n ONTrack2_env blast
RUN conda install -n ONTrack2_env emboss
RUN conda install -n ONTrack2_env vsearch
RUN conda install -n ONTrack2_env seqtk
RUN conda install -n ONTrack2_env mafft
RUN conda install -n ONTrack2_env minimap2
RUN conda install -n ONTrack2_env samtools
RUN conda install -n ONTrack2_env racon
RUN conda install -n ONTrack2_env bcftools
WORKDIR /home/
