FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive
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
    conda config --add channels conda-forge
# git clone --recursive https://github.com/HKU-BAL/ClusterV.git

COPY . .
    
RUN conda env create -f clusterV.yml

ENV PATH /opt/conda/envs/clusterV/bin:$PATH
ENV CONDA_DEFAULT_ENV clusterV 

RUN /bin/bash -c "source activate clusterV" && \
    pypy3 -m ensurepip && \
    pypy3 -m pip install mpmath==1.2.1 && \
    pypy3 -m pip install --no-cache-dir intervaltree==3.0.2 && \
    echo "source activate clusterV " > ~/.bashrc && \
    cd /opt/bin/Clair3 && \
    make PREFIX=/opt/conda/envs/clusterV && \
    mkdir -p /opt/bin/Clair3/models && \
    wget http://www.bio8.cs.hku.hk/clair3/clair3_models/clair3_models.tar.gz && \
    tar -zxvf /opt/bin/Clair3/clair3_models.tar.gz -C /opt/bin/Clair3/models && \
    rm -f /opt/bin/Clair3/clair3_models.tar.gz && \
    cd /opt/bin

RUN echo 'will cite' | parallel --citation || true

RUN cd /opt/bin && \
    PREFIX=/opt/conda/envs/clusterV && \
    PYTHON=/opt/conda/envs/clusterV/bin/python


