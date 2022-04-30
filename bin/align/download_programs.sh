
cd ../../

mkdir -p programs/
cd programs/

# CellRanger

wget -O cellranger-3.1.0.tar.gz "http://cf.10xgenomics.com/releases/cell-exp/cellranger-3.1.0.tar.gz?Expires=1584095527&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cDovL2NmLjEweGdlbm9taWNzLmNvbS9yZWxlYXNlcy9jZWxsLWV4cC9jZWxscmFuZ2VyLTMuMS4wLnRhci5neiIsIkNvbmRpdGlvbiI6eyJEYXRlTGVzc1RoYW4iOnsiQVdTOkVwb2NoVGltZSI6MTU4NDA5NTUyN319fV19&Signature=EDB9IBjjgjE6WwA9nCiaK40p66NOql9Dm3grzVz-7k0-m5zf8G1ZwoCnXd5JbZHH~yGxi2JV~UqQiuj7ysLfuSYey20TUpwBEUsasxs-yJl8AVY1ejOlq36K5VL-I2mHmivVx~TfTFGSV71Kq0KEAS4cPaffsM5OuxeLgV3RPchdwHCI3~cGo8d-l1DDa21dKLIAXL1DpvsJmfbIW8nO8eM2B3J1fkgZ83bdy0XrmWxiRrLkYG3Tap0ZSKkldq46bSkvOEq0vhiuuuKNLrru6r-KbMX7FCBRb74y3s4R0JOWBJoNtBm0Ed6Qc0ojGUHSMdADQbQUjd4yNlN-QI~mzg__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"

tar xzf cellranger-3.1.0.tar.gz
rm cellranger-3.1.0.tar.gz

gunzip cellranger-3.1.0/cellranger-cs/3.1.0/lib/python/cellranger/barcodes/3M-february-2018.txt.gz


# STAR

. /u/local/Modules/default/init/modules.sh
module load gcc/7.2.0

wget https://github.com/alexdobin/STAR/archive/2.7.3a.tar.gz
tar -xzf 2.7.3a.tar.gz
rm 2.7.3a.tar.gz
cd STAR-2.7.3a
cd source/
make STAR
cd ../../

# Samtools
wget https://github.com/samtools/samtools/releases/download/1.6/samtools-1.6.tar.bz2
tar xjvf samtools-1.6.tar.bz2
rm samtools-1.6.tar.bz2
cd samtools-1.6
./configure --disable-lzma --without-curses --prefix=$PWD
make
make install
cd htslib-1.6
./configure --disable-lzma --prefix=$PWD
make
make install
cd ../../

# FastQC
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
unzip fastqc_v0.11.9.zip
rm fastqc_v0.11.9.zip
chmod u+x FastQC/fastqc

echo "10X TSO   AAGCAGTGGTATCAACGCAGAGTACAT" >> FastQC/contaminant_10X_TSO.txt


