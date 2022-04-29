#!/u/local/apps/R/3.5.1/gcc-4.9.3_MKL-2018/bin/Rscript
#$ -S /u/local/apps/R/3.5.1/gcc-4.9.3_MKL-2018/bin/Rscript
#$ -v R_LIBS_USER=/u/project/pajukant/malvarez/lib/R_3.5.1/
#$ -v LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/u/local/compilers/gcc/4.9.3/lib64:/u/project/pajukant/malvarez/lib/
#$ -cwd
#$ -j y
#$ -pe shared 8
#$ -l h_data=4G,h_vmem=32G,h_rt=3:00:00,highp
#$ -M malvarez@mail
#  Notify at beginning and end of job
#$ -m n
#$ -r n
#$ -o cca.init.R.log

setwd("../../")

library(Seurat)
library(ggplot2)
library(scales)
library(gridExtra)
source("scripts/standard_seurat.R")

library(future)
plan("multiprocess", workers = 8)

maxsize <- 6000 * 1024 ^ 2
options(future.globals.maxSize = maxsize)

#=========================================
# Functions
#=========================================

create_dir <- function(p){
    dir.create(p, showWarnings=FALSE, recursive=TRUE)
}

seurat_cluster <- function(x, 
                           n_dim=30, 
                           resolution = 0.8, 
                           python_lib="/u/local/apps/python/3.6.1-shared/bin/python3.6"){
	require(reticulate)
	reticulate::use_python(python_lib)
	x <- FindNeighbors(x, dims = 1:n_dim, verbose = TRUE)
	x <- FindClusters(x, resolution = resolution, verbose = TRUE)
	x <- RunUMAP(x, dims = 1:n_dim, reduction = "pca", verbose = TRUE)
	return(x)
}

#=========================================
# Set variables
#=========================================

n_threads <- 8

meta <- read.table("data/sample/samples.csv", 
                   header = TRUE, 
                   row.names = 1, 
                   stringsAsFactors = FALSE,
                   sep = ",")

ids <- rownames(meta)

# Gene data
gene_info <- read.table("data/ref/gencode26/gencode.v26.annotation.txt", 
                        header = TRUE, 
                        row.names = 1, 
                        stringsAsFactors = FALSE, 
                        sep = "\t")

# Set directories
dir_out <- "data/processed/d_sct_cca_all/"; create_dir(dir_out)
# dir_plot <- "exp/d_sct_cca_all/cca_init/"; create_dir(dir_plot)

#=========================================
# Read in Seurat
#=========================================

# Read in Seurat objects
sct_dir <- "data/processed/diem/"
seur.list <- lapply(ids, function(x){ 
                    path <- file.path(sct_dir, x, "seur.sct.rds")
                    ret <- readRDS(path)
                    return(ret)
                   })

seur.list <- lapply(seur.list, function(x){
                    x <- PercentageFeatureSet(x, features = "ENSG00000163631.16", col.name = "pct.alb")
                    x@active.assay = "SCT"
                    # x <- subset(x, subset = pct.mt < 5)
                    # x <- subset(x, subset = nFeature_RNA >= 200)
                    x <- RunPCA(x, npcs = 30, verbose = FALSE) 
                    return(x) })

ndim <- 30
seur <- sct_cca(seur.list, n_pcs = ndim, nfeatures = 3e3)
seur <- seurat_cluster(seur, n_dim = ndim, resolution = c(0.2, 0.5, 1, 1.5, 2))

saveRDS(seur, paste0(dir_out, "seur.cca.rds"))

