#!/bin/bash

cd ../../data/ref/
mkdir -p gene_lists
cd gene_lists

# Download list of ribosomal genes
# RPL genes
wget -O rpl_info.txt https://www.genenames.org/cgi-bin/genegroup/download?id=729
sleep 10
awk -v FS="\t" 'NR > 1 {print $2}' rpl_info.txt > rpl_genes.txt

# RPS genes
wget -O rps_info.txt https://www.genenames.org/cgi-bin/genegroup/download?id=728&type=node
sleep 10
awk -v FS="\t" 'NR > 1 {print $2}' rps_info.txt > rps_genes.txt

# MRP genes
wget -O mrp_info.txt https://www.genenames.org/cgi-bin/genegroup/download?id=646&type=node
sleep 10
awk -v FS="\t" 'NR > 1 {print $2}' mrp_info.txt > mrp_genes.txt

