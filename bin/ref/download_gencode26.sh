#!/bin/bash

cd ../../

mkdir -p data/ref/gencode26/
cd data/ref/gencode26/

wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/gencode.v26.annotation.gtf.gz
gunzip -f gencode.v26.annotation.gtf.gz

cd ../../../

