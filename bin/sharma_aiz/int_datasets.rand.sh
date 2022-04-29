#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe shared 16
#$ -l h_vmem=64G,h_data=4G,h_rt=12:00:00,highp
#$ -v QQAPP=openmp
#$ -M malvarez@mail
#  Notify at beginning and end of job
#$ -m n
#$ -r n
#$ -o int_datasets.rand.sh.log

. ../../conda_init.sh

conda activate liver_snrna_hcc
Rscript int_datasets.rand.R
conda deactivate

