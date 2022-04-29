#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe shared 8
#$ -l h_vmem=32G,h_data=4G,h_rt=4:00:00,highp
#$ -v QQAPP=openmp
#$ -M malvarez@mail
#  Notify at beginning and end of job
#$ -m n
#$ -r n
#$ -o clust_sharma.sh.log

. ../../conda_init.sh

conda activate liver_snrna_hcc
Rscript clust_sharma.R
conda deactivate

