#!/bin/bash
#$ -cwd
#$ -j y
#$ -l h_vmem=16G,h_data=16G,h_rt=12:00:00,highp
#$ -v QQAPP=openmp
#$ -M malvarez@mail
#  Notify at beginning and end of job
#$ -m a
#$ -r n
#$ -o mut_expr.sh.log

. ../../conda_init.sh

conda activate liver_snrna_hcc
Rscript mut_expr.R
conda deactivate

