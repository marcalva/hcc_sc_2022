#!/bin/bash

d=$PWD
cd ../

mkdir -p data/raw/sharma2019
cd data/raw/sharma2019

wget https://md-datasets-public-files-prod.s3.eu-west-1.amazonaws.com/df2fcbd6-258c-4ada-84ce-5e6fc4c2fdbd
mv df2fcbd6-258c-4ada-84ce-5e6fc4c2fdbd HCC.h5ad
# objects in HCC.h5ad
# X contains a normalized expression data set
# obs contains the cell meta data, including patient, tumor, louvain cluster, etc.
# obsm contains PCA and UMAP data
# raw.cat contains "nan" string
# raw.X contains column-wise sparse matrix
# raw.var contains the gene names for raw.X (19,852 genes)
# var contains the gene names for X (2,608 genes)

wget https://md-datasets-public-files-prod.s3.eu-west-1.amazonaws.com/e9a82be5-7d36-4b44-9ce2-6282b4ebd654
mv e9a82be5-7d36-4b44-9ce2-6282b4ebd654 HCCF1F2.h5ad

wget https://data.mendeley.com/public-files/datasets/6wmzcskt6k/files/45da9074-e3a2-42fb-9978-c9cf1d912980/file_downloaded

wget https://md-datasets-cache-zipfiles-prod.s3.eu-west-1.amazonaws.com/6wmzcskt6k-1.zip

unzip 6wmzcskt6k-1.zip

cd $d

