#!/u/local/apps/R/3.5.1/gcc-4.9.3_MKL-2018/bin/Rscript
#$ -S /u/local/apps/R/3.5.1/gcc-4.9.3_MKL-2018/bin/Rscript
#$ -cwd
#$ -j y
#$ -pe shared 8
#$ -l h_data=4G,h_vmem=16G,h_rt=20:00:00,highp
#$ -M malvarez@mail
#  Notify at beginning and end of job
#$ -t 1-6
#$ -m a
#$ -r n
#$ -o run_diem.2.R.log.$TASK_ID

# Run DIEM

setwd("../../")

library(diem)
library(Matrix)
library(Seurat)
library(ggplot2)


#=========================================
# Functions
#=========================================

create_dir <- function(p){
    dir.create(p, recursive=TRUE, showWarnings=FALSE)
}

#=========================================
# Set sample and variables
#=========================================

args = commandArgs(TRUE)
task = as.integer(args[1])

task = Sys.getenv(x = "SGE_TASK_ID")
if (task == ""){
    stop("Set SGE_TASK_ID environment variable")
}

task <- as.integer(task)

meta <- read.table("data/sample/samples.csv", 
                   header = TRUE, 
                   sep = ",", 
                   stringsAsFactors = FALSE, 
                   check.names = FALSE)

gtf_path <- "data/ref/gencode26/gencode.v26.annotation.fltrd.gtf"

gene_info <- read.table("data/ref/gencode26/gencode.v26.annotation.txt",
                        header = TRUE,
                        row.names = 1,
                        stringsAsFactors = FALSE,
                        sep = "\t")
rb_genes <- readLines("data/ref/ribosomal_genes.txt")
rb_ids <- rownames(gene_info)[gene_info$Name %in% rb_genes]

rbr_genes <- c("CH507-513H4.1", "CH507-528H12.1")
rbr_ids <- c("ENSG00000278996.1", "ENSG00000280441.2")
hb_ids <- c("ENSG00000206172.8", "ENSG00000188536.12", "ENSG00000244734.3")
mt_ids <- rownames(gene_info)[gene_info[,"Chrm"] == "chrM"]
malat1_ids <- rownames(gene_info)[gene_info[,"Name"] == "MALAT1"]

rm_ids <- c(rb_ids, rbr_ids, hb_ids, mt_ids)
keep_ids <- setdiff(rownames(gene_info), rm_ids)

sample <- meta[task, "Sample"]

n_threads <- 8

#=========================================
# Output directories
#=========================================

dir_out <- paste0("data/processed/diem/", sample, "/"); create_dir(dir_out)
dir_plot <- paste0("exp/diem/plots/", sample, "/"); create_dir(dir_plot)
dir_fltr <- paste0("exp/diem/fltr_ids/", sample, "/"); create_dir(dir_fltr)
dir_de <- paste0("exp/diem/DE/", sample, "/"); create_dir(dir_de)

#=========================================
# Run DIEM
#=========================================

# Create output directories
cat(paste0("Running ", sample, "\n"))
cat("Reading expression data\n")

sce <- readRDS(paste0(dir_out, sample, ".diem_sce.rds"))

# Get quantile threshold for clustering droplets
n_expect <- 10000
quant <- 0.99
drop_dat <- droplet_data(sce)
n_expect <- min(n_expect, nrow(drop_dat))
o <- order(drop_dat[,"n_genes"], decreasing = TRUE)
ix_o <- round(n_expect * (1 - quant))
ix <- o[ix_o]
min_g <- drop_dat[ix, "n_genes"] / 10

sce <- get_pcs(sce, min_genes = min_g)
sce <- init(sce, threads = n_threads, k_init = 30)

#=========================================
# Run EM
#=========================================

aprior <- 1
pprior <- 1
sce <- run_em(sce, alpha_prior = aprior, pi_prior = pprior, 
              max_iter = 1000, eps = 1e-4, psc = 0, 
              threads = n_threads)

# Save DIEM object
cat("Saving\n")
saveRDS(sce, paste0(dir_out, sample, ".diem_sce.rds"))

