#!/u/local/apps/R/3.5.1/gcc-4.9.3_MKL-2018/bin/Rscript
#$ -S /u/local/apps/R/3.5.1/gcc-4.9.3_MKL-2018/bin/Rscript
#$ -cwd
#$ -j y
#$ -l h_data=8G,h_vmem=8G,h_rt=1:00:00,highp
#$ -M malvarez@mail
#  Notify at beginning and end of job
#$ -m n
#$ -r n
#$ -o plot_clust.R.log

setwd("../../")

library(Seurat)
library(ggplot2)
library(scales)
library(gridExtra)
source("scripts/plot.R")

#=========================================
# Functions
#=========================================

create_dir <- function(p){
    dir.create(p, showWarnings=FALSE, recursive=TRUE)
}

#=========================================
# Set variables
#=========================================

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
dir_plot <- "exp/d_sct_cca_all/cca_init/plots/"; create_dir(dir_plot)

#=========================================
# Read in Seurat
#=========================================

seur <- readRDS(paste0(dir_out, "seur.cca.rds"))

plot_umap_labels(seur, 
                 colname = "integrated_snn_res.0.2", 
                 legend_title = "Res. 0.2", 
                 label = TRUE, 
                 size = 0.3, label_size = 4,
                 alpha = 0.5)
ggsave(paste0(dir_plot, "UMAP.res.0.2.jpeg"), 
       width = 5, 
       height = 5, 
       units = "in")

plot_umap_labels(seur, 
                 colname = "integrated_snn_res.0.5", 
                 legend_title = "Res. 0.5", 
                 label = TRUE, 
                 size = 0.3, label_size = 4,
                 alpha = 0.5)
ggsave(paste0(dir_plot, "UMAP.res.0.5.jpeg"), 
       width = 5, 
       height = 5, 
       units = "in")

plot_umap_labels(seur, 
                 colname = "integrated_snn_res.1", 
                 legend_title = "Res. 1", 
                 label = TRUE, 
                 size = 0.3, label_size = 4,
                 alpha = 0.5)
ggsave(paste0(dir_plot, "UMAP.res.1.jpeg"), 
       width = 5, 
       height = 5, 
       units = "in")

plot_umap_labels(seur, 
                 colname = "integrated_snn_res.1.5", 
                 legend_title = "Res. 1.5", 
                 label = TRUE, 
                 size = 0.3, label_size = 4, 
                 alpha = 0.5)
ggsave(paste0(dir_plot, "UMAP.res.1.5.jpeg"), 
       width = 5, 
       height = 5, 
       units = "in")

plot_umap_labels(seur, 
                 colname = "integrated_snn_res.2", 
                 legend_title = "Res. 2", 
                 label = TRUE, 
                 size = 0.3, label_size = 4, 
                 alpha = 0.5)
ggsave(paste0(dir_plot, "UMAP.res.2.jpeg"), 
       width = 5, 
       height = 5, 
       units = "in")

p <- plot_umap_c(seur, 
                 colname = "nCount_RNA", 
                 legend_title = "UMI", 
                 trans = "log10", 
                 size = 0.3, 
                 alpha = 0.5)
ggsave(paste0(dir_plot, "UMAP.ncount.jpeg"), 
       width = 5, 
       height = 5, 
       units = "in")

p <- plot_umap_c(seur, 
                 colname = "pct.mt", 
                 legend_title = "MT%", 
                 trans = "log1p", 
                 size = 0.3, 
                 alpha = 0.5)
ggsave(paste0(dir_plot, "UMAP.pct_mt.jpeg"), 
       width = 5, 
       height = 5, 
       units = "in")

p <- plot_umap_c(seur, 
                 colname = "pct.malat1", 
                 legend_title = "MALAT1%", 
                 trans = "log1p", 
                 size = 0.3, 
                 alpha = 0.5)
ggsave(paste0(dir_plot, "UMAP.pct_malat1.jpeg"), 
       width = 5, 
       height = 5, 
       units = "in")

p <- plot_umap_c(seur, 
                 colname = "pct.rb", 
                 legend_title = "RB%", 
                 trans = "log1p", 
                 size = 0.3, 
                 alpha = 0.5)
ggsave(paste0(dir_plot, "UMAP.pct_rb.jpeg"), 
       width = 5, 
       height = 5, 
       units = "in")

p <- plot_umap_c(seur, 
                 colname = "pct.rbr", 
                 legend_title = "RBR%", 
                 trans = "log1p", 
                 size = 0.3, 
                 alpha = 0.5)
ggsave(paste0(dir_plot, "UMAP.pct_rbr.jpeg"), 
       width = 5, 
       height = 5, 
       units = "in")

plot_umap_labels(seur, 
                 colname = "SampleID", 
                 legend_title = "Sample", 
                 label = FALSE, 
                 size = 0.3,
                 alpha = 0.5)
ggsave(paste0(dir_plot, "UMAP.sample.jpeg"), 
       width = 5, 
       height = 5, 
       units = "in")

plot_umap_labels(seur, 
                 colname = "Patient", 
                 legend_title = "Patient", 
                 label = FALSE, 
                 size = 0.3,
                 alpha = 0.5)
ggsave(paste0(dir_plot, "UMAP.patient.jpeg"), 
       width = 5, 
       height = 5, 
       units = "in")

plot_umap_labels(seur, 
                 colname = "Cluster", 
                 legend_title = "DIEM Cluster", 
                 label = TRUE, 
                 size = 0.3,
                 alpha = 0.5)
ggsave(paste0(dir_plot, "UMAP.diem_cluster.jpeg"), 
       width = 10, 
       height = 10, 
       units = "in")

