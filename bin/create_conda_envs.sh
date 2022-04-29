#!/bin/bash

cd ../

. conda_init.sh

conda update -n base -c defaults conda

envname="liver_snrna_hcc"

conda create \
    --name $envname \
    -c conda-forge \
    r-base=4.1.1 \
    python=3.9.7

conda install \
    -y -n $envname \
    -c conda-forge \
    r-seurat=4.0.3 \
    r-gamlss.dist=5.3_2 \
    r-statmod=1.4.36

conda install \
    -y -n $envname \
    -c conda-forge \
    r-ggplot2=3.3.5

conda install \
    -y -n $envname \
    -c conda-forge \
    xz=5.2.5 \
    zlib=1.2.11

conda install \
    -y -n $envname \
    -c conda-forge \
    r-biocmanager=1.30.16

conda install \
    -y -n $envname \
    -c bioconda \
    bioconductor-edger=3.34.0

conda install \
    -y -n $envname \
    -c conda-forge \
    openpyxl

# doesn't work

conda install \
    -y -n $envname \
    -c bioconda \
    scanpy=1.7.2

