#!/u/local/apps/python/3.6.1/bin/python3
#$ -S /u/local/apps/python/3.6.1/bin/python3
#$ -v PYTHONPATH=/u/project/pajukant/malvarez/.local/lib/python3.6/site-packages/:/u/local/apps/python/3.6.1/lib/python3.6
#$ -cwd
#$ -j y
#$ -pe shared 12
#$ -l h_data=4G,h_vmem=48G,h_rt=10:00:00,highp
#$ -M malvarez@mail
#  Notify at beginning and end of job
#$ -t 5-5
#$ -m a
#$ -r n
#$ -o align.py.log.$TASK_ID


import os
import sys
import shutil
import pandas as pd
from subprocess import call

os.chdir("../../")

##########################################
# Functions
##########################################

# x is a series or list of strings
def cat_read_path(x, r = 0, prfx = "", sep = ";", j = ","):
    s = [prfx + y.split(sep)[r] for y in x]
    ret = j.join(s)
    return ret

##########################################
##########################################


##########################################
# Parameters
##########################################
    
samtools = "programs/samtools-1.6/bin/samtools"
star = "programs/STAR-2.7.3a/bin/Linux_x86_64_static/STAR"
fastq_pre = "data/raw/"
sfn = "data/raw/sample/sample_ids.fastq.trim.tsv"
genome_dir = "data/processed/STAR_indexes/gencode26/"
wl = "programs/cellranger-3.1.0/cellranger-cs/3.1.0/lib/python/cellranger/barcodes/3M-february-2018.txt"
umi_len = 12
n_threads = 12
outdir = "data/processed/STAR_output/gencode26/"
task = int(os.environ['SGE_TASK_ID']) - 1
outmatchNmin = 30
outmatchNminOverLread = 0

##########################################
##########################################

samples = pd.read_csv(sfn, sep = "\t")
samples.index = samples['Sample']

uniq_ids = samples['Sample'].unique()

id = uniq_ids[task]
r1 = samples.loc[id, 'R1']
r2 = samples.loc[id, 'R2']
outpre = outdir + "/" + id + "/"

if os.path.exists(outpre):
    shutil.rmtree(outpre)

os.makedirs(outpre)

command = [star,
    "--runThreadN", str(n_threads), 
    "--genomeDir", genome_dir, 
    "--readFilesIn", r2, r1, 
    "--readFilesCommand", "zcat", 
    "--soloType", "CB_UMI_Simple", 
    "--soloCBwhitelist", wl, 
    "--soloUMIlen", str(umi_len), 
    "--soloFeatures", "Gene", "GeneFull", "SJ", "Velocyto", 
    "--soloUMIfiltering", "MultiGeneUMI", 
    "--soloCBmatchWLtype", "1MM_multi_pseudocounts",
    "--outSAMattributes", "NH", "HI", "nM", "AS", "CR", "UR", "CB", "UB", "GX", "GN", "sS", "sQ", "sM", 
    "--outSAMtype", "BAM", "SortedByCoordinate", 
    "--outFileNamePrefix", outpre, 
    "--outFilterMatchNmin", str(outmatchNmin), 
    "--outFilterMatchNminOverLread", str(outmatchNminOverLread), 
    "--alignIntronMax", "100000", 
    "--alignSJstitchMismatchNmax", "5", "-1", "5", "5", 
    "--chimOutType", "Junctions", 
    "--chimSegmentMin", "12", 
    "--chimJunctionOverhangMin", "12", 
    "--chimOutJunctionFormat", "1"]

# S1_hcc sequenced on MiSeq, so clip 150 bp length reads to 91.
# if id == "S1_hcc":
#     command = command + ["--clip3pNbases", "59", "0"]

call(command)

bamfn = outpre + "Aligned.sortedByCoord.out.bam"

call([samtools, "index", "-@", str(n_threads), bamfn])

