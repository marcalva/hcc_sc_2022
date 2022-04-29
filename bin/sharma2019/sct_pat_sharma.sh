#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe shared 12
#$ -l h_vmem=48G,h_data=4G,h_rt=4:00:00,highp
#$ -v QQAPP=openmp
#$ -M malvarez@mail
#  Notify at beginning and end of job
#$ -m n
#$ -r n
#$ -o sct_pat_sharma.sh.log

. ../../conda_init.sh

conda activate liver_snrna_hcc
Rscript sct_pat_sharma.R
conda deactivate

