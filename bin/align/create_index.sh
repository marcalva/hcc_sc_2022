#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe shared 16
#$ -l h_data=3G,h_rt=10:00:00,highp
#$ -v QQAPP=openmp
#$ -M malvarez@mail
#  Notify at beginning and end of job
#$ -m a
#$ -r n
#$ -o create_index.sh.log

cd ../../

star="programs/STAR-2.7.3a/bin/Linux_x86_64_static/STAR"

genomedir="data/processed/STAR_indexes/gencode26/"
mkdir -p $genomedir
rm -rf ${genomedir}*

refFasta="data/ref/gencode26/GRCh38.primary_assembly.genome.fa"
gencodeAnno="data/ref/gencode26/gencode.v26.annotation.gtf"

$star --runThreadN 16 \
    --runMode genomeGenerate \
    --genomeDir $genomedir \
    --genomeFastaFiles $refFasta \
    --sjdbGTFfile $gencodeAnno \
    --sjdbOverhang 90

