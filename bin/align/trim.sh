#!/bin/bash
#$ -cwd
#$ -v PYTHONPATH=/u/project/pajukant/malvarez/.local/lib/python3.6/site-packages/:/u/local/apps/python/3.6.1/lib/python3.6
#$ -j y
#$ -pe shared 8
#$ -l h_data=1G,h_vmem=8G,h_rt=8:00:00,highp
#$ -M malvarez@mail
#  Notify at beginning and end of job
#$ -t 1-6
#$ -m a
#$ -r n
#$ -o trim.sh.log.$TASK_ID

cd ../../

. /u/local/Modules/default/init/modules.sh
module load python/3.6.1

home="/u/project/pajukant/malvarez/"
cutadapt="${home}/.local/bin/cutadapt"

tso="AAGCAGTGGTATCAACGCAGAGTACATGGG"
primer="GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"
seq5="${primer}${tso}"
seq3='A{20}'

sample=$( ls data/raw/fastq | awk "NR == $SGE_TASK_ID" )

fastqdir="data/raw/fastq/${sample}/"

dirout="data/processed/fastq_trim/${sample}/"
mkdir -p $dirout

rm ${dirout}/*

for fq in $( find ${fastqdir}/*_R2_* -printf "%f\n" ); do
    fastq1in="${fastqdir}/${fq}"
    fastq1out="${dirout}/${fq}"
    fastq2in=$( echo $fastq1in | sed 's/_R2_/_R1_/g' )
    fastq2out=$( echo $fastq1out | sed 's/_R2_/_R1_/g' )
    $cutadapt \
        --cores 8 \
        -n 10 \
        -m 20: \
        -g "X${seq5}" \
        -a "${seq3}$" \
        -o $fastq1out \
        -p $fastq2out \
        $fastq1in \
        $fastq2in
done


