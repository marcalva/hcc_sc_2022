
# Human liver single nucleus and single cell RNA-sequencing identify a hepatocellular carcinoma-associated cell-type affecting survival

This repo contains the code used to analyze the data from 
Alvarez, Benhammou et al. 2021.

The file `run_all.sh` gives the order in which to run the scripts. 
The scripts to run the analysis are found in the `bin` directory.

## Download data

To download the data from [Aizarani et al.](https://www.nature.com/articles/s41586-019-1373-2) 
and [Sharma et al.](https://www.sciencedirect.com/science/article/pii/S0092867420310825?via%3Dihub), 
run the scripts in the `data-raw` directory.

For the NAFLD-related HCC set, the scripts use the raw counts from STARSolo. 
The DIEM-filtered counts can be downloaded from GEO under accession number 
[GSE189175](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE189175).
The DIEM filtering scripts can be skipped and the code for clustering 
can read in the matrices downloaded from GEO.

## Alignment

Although sequencing data is not available for the NAFLD-related HCC 
data set, the alignment scripts are provided in `bin/align`.

## Analysis

Set up the conda environment using the `bin/create_conda_envs.sh` and 
`bin/add_conda_pkg.sh` scripts.

The scripts in `bin/aizarani2019`, `bin/sharma2019`, and 
`bin/d_sct_cca_all` run clustering on the Aizarani et al healthy liver 
scRNA-seq reference data, the Sharma et al HCC data, and the Rao et al
NAFLD-related HCC data, respectively.

The scripts in `bin/sharma_aiz` run the integrated clustering across 
the three cohorts. The scripts in `bin/tcga_hcc` run the TCGA 
decomposition analysis, and finally the `bin/GSE14520` run the 
LCI decomposition analysis.

## Manuscript figures

To reproduce the manuscript figures and tables, run the scripts in the 
`bin/manuscript` folder.

