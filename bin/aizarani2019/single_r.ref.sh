#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe shared 8
#$ -l h_vmem=32G,h_data=4G,h_rt=8:00:00,highp
#$ -v QQAPP=openmp
#$ -M malvarez@mail
#  Notify at beginning and end of job
#$ -m n
#$ -r n
#$ -o single_r.ref.sh.log

. ../../conda_init.sh

conda activate liver_snrna_hcc
Rscript single_r.ref.R
conda deactivate

