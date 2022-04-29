#!/u/local/apps/R/3.5.1/gcc-4.9.3_MKL-2018/bin/Rscript
#$ -S /u/local/apps/R/3.5.1/gcc-4.9.3_MKL-2018/bin/Rscript
#$ -cwd
#$ -j y
#$ -pe shared 8
#$ -l h_data=4G,h_vmem=32G,h_rt=1:00:00,highp
#$ -M malvarez@mail
#  Notify at beginning and end of job
#$ -t 1-6
#$ -m a
#$ -r n
#$ -o run_diem.4.R.log.$TASK_ID

# Run DIEM

setwd("../../")

library(diem)
library(Matrix)
library(Seurat)
library(dplyr)
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

task = Sys.getenv(x = "SGE_TASK_ID")
if (task == ""){
    stop("Set SGE_TASK_ID environment variable")
}
task = as.integer(task)

meta <- read.table("data/sample/samples.csv", 
                   header = TRUE, 
                   sep = ",", 
                   stringsAsFactors = FALSE, 
                   check.names = FALSE)
rownames(meta) <- meta[,"Sample"]

gtf_path <- "data/ref/gencode26/gencode.v26.annotation.fltrd.gtf"

gene_info <- read.table("data/ref/gencode26/gencode.v26.annotation.txt",
                        header = TRUE,
                        row.names = 1,
                        stringsAsFactors = FALSE,
                        sep = "\t")
info_keep <- c("Chrm", "Type", "Name")
mt_ens <- rownames(gene_info)[ gene_info[,"Chrm"] == "chrM" ]
mt_genes <- gene_info[gene_info[,"Chrm"] == "chrM", "Name"]
malat1 <- rownames(gene_info)[ gene_info[,"Name"] == "MALAT1" ]
malat1 <- malat1[1]

sample <- meta[task, "Sample"]

# Map sample IDs
sid <- c("Sample 1 Normal", "Sample 1 Tumor", 
         "Sample 2 NASH", "Sample 2 Tumor", 
         "Sample 3 NASH", "Sample 3 Tumor")
names(sid) <- c("S2_norm", "S1_hcc", 
                "S2_nash", "S2_hcc", 
                "S3_nash", "S3_hcc")

#=========================================
# Output directories
#=========================================

dir_out <- paste0("data/processed/diem/", sample, "/"); create_dir(dir_out)
dir_plot <- paste0("exp/diem/plots/", sample, "/"); create_dir(dir_plot)
dir_fltr <- paste0("exp/diem/fltr_ids/", sample, "/"); create_dir(dir_fltr)
dir_de <- paste0("exp/diem/DE/", sample, "/"); create_dir(dir_de)
plot_prefix <- paste0(dir_plot, sample, ".") # File prefix for dist plot

#=========================================
# Read data
#=========================================

sce <- readRDS(paste0(dir_out, sample, ".diem_sce.rds"))

#=========================================
# Estimate debris
#=========================================

sce <- estimate_dbr_score(sce, max_genes = 100, thresh_genes = 100)
sce@debris_genes[,"Name"] <- gene_info[sce@debris_genes[,"gene"], "Name"]

td <- sce@test_data

# Plot cluster means
p <- plot_clust(sce, feat_x = "n_genes", feat_y = "score.debris") + 
ggtitle(sample)
ggsave(paste0(plot_prefix, "cluster.n_gene.dbr_score.jpeg"), 
       width = 5, height = 5)

# Plot droplet
p <- plot_data(sce, feat_x = "n_genes", feat_y = "score.debris", 
          alpha = 0.1, log_x = TRUE) + 
ggtitle(sample)
ggsave(paste0(plot_prefix, "droplet.n_gene.dbr_score.jpeg"), 
       width = 5, height = 5)

saveRDS(sce, paste0(dir_out, sample, ".diem_sce.rds"))

