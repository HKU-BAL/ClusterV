FROM ubuntu:16.04

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8 PATH=/opt/bin:/opt/conda/bin:$PATH

# update ubuntu packages
RUN apt-get update --fix-missing && \
    yes|apt-get upgrade && \
    apt-get install -y \
        wget \
        bzip2 \
        make \
        g++ \
        git \
        unzip \
        time \
        libboost-graph-dev && \
    rm -rf /bar/lib/apt/lists/*

WORKDIR /opt/bin

# install anaconda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda && \
    rm Miniconda3-latest-Linux-x86_64.sh && \
    conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    wget https://github.com/HKU-BAL/ClusterV/archive/refs/tags/v1.2.zip && \
    unzip v1.2.zip && mv ClusterV-1.2 ClusterV && \
    cd ClusterV && \
    conda env create -f clusterV.yml

# git clone https://github.com/HKU-BAL/ClusterV.git && \

ENV PATH /opt/conda/envs/clusterV/bin:$PATH
ENV CONDA_DEFAULT_ENV clusterV 

RUN /bin/bash -c "source activate clusterV" && \
    pypy3 -m ensurepip && \
    pypy3 -m pip install --no-cache-dir intervaltree==3.0.2 && \
    echo "source activate clusterV " > ~/.bashrc

RUN echo 'will cite' | parallel --citation || true


COPY . .

run cd /opt/bin && \
    PREFIX=/opt/conda/envs/clusterV && \
    PYTHON=/opt/conda/envs/clusterV/bin/python


