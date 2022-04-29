#!/bin/bash

d=$PWD
cd ../

mkdir -p data/raw/aizarani2019
cd data/raw/aizarani2019
wget -q -N \
    https://ftp.ncbi.nlm.nih.gov/geo/series/GSE124nnn/GSE124395/suppl/GSE124395%5FNormalhumanliverdata%2ERData%2Egz \
    -O GSE124395_Normalhumanliverdata.RData.gz
gunzip -f GSE124395_Normalhumanliverdata.RData.gz

wget -q -N \
    https://ftp.ncbi.nlm.nih.gov/geo/series/GSE124nnn/GSE124395/suppl/GSE124395%5Fclusterpartition%2Etxt%2Egz \
    -O GSE124395_clusterpartition.txt.gz

gunzip -f GSE124395_clusterpartition.txt.gz

cd $d

