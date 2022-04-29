
# plot UMAP of un-integrated data
# show biposy/patient specific effects.

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
# Plots
#=========================================

dir_plot <- "exp/d_sct_cca_all/merge.raw/plots/"; create_dir(dir_plot)

infn <- "data/processed/d_sct_cca_all/merge/seur.raw.rds"
seur <- readRDS(infn)

# replace "NASH with normal
seur$SampleID <- sub("NASH", "Normal", seur$SampleID)
seur$SampleID <- sub("Normal", "Non-tumor", seur$SampleID)

# order factors from tumor to non-tumor
datf_s <- seur@meta.data[,c("SampleID", "Patient", "Tumor")]
datf_s[,"SampleID"] <- as.character(datf_s[,"SampleID"])
datf_s <- unique(datf_s)
o <- order(datf_s$Tumor, datf_s$Patient)

lv <- datf_s[o,"SampleID"]

seur$SampleID <- factor(seur$SampleID, levels = lv)

# patient colors
pt_cols <- read.table("data/ref/TumorColors/patient_colors.1.txt", 
                      header=FALSE, row.names=1, sep=",", stringsAsFactors=FALSE)
pt_cols_v <- pt_cols[,1]
names(pt_cols_v) <- rownames(pt_cols)

# tumor colors
fn <- "data/ref/TumorColors/tumor_colors.1.txt"
tum_cols <- read.table(fn, header=FALSE, sep=",", stringsAsFactors=FALSE)
tum_colsv <- tum_cols[,2]
names(tum_colsv) <- tum_cols[,1]
names(tum_colsv)[names(tum_colsv) == "Tumor"] <- "Tumor"
names(tum_colsv)[names(tum_colsv) == "NonTumor"] <- "Non-tumor"

seur$Tumor <- as.character(seur$Tumor)
knt <- seur$Tumor == "NonTumor"
seur@meta.data[knt, "Tumor"] <- "Non-tumor"
# seur$Tumor <- factor(seur$Tumor, levels = c("Tumor", "Non-tumor"))

# sample colors
sm_cols <- read.table("data/ref/TumorColors/sample_colors.1.txt", 
                      header=FALSE, row.names=1, sep=",", stringsAsFactors=FALSE)
sm_cols_v <- sm_cols[,1]
names(sm_cols_v) <- rownames(sm_cols)
names(sm_cols_v) <- sub("Normal", "Non-tumor", names(sm_cols_v))

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

p <- plot_umap_c(seur, 
                 colname = "pct.alb", 
                 legend_title = "ALB%", 
                 trans = "log1p", 
                 size = 0.3, 
                 alpha = 0.5)
ggsave(paste0(dir_plot, "UMAP.pct_alb.jpeg"), 
       width = 5, 
       height = 5, 
       units = "in")

plot_umap_labels(seur, 
                 colname = "SampleID", 
                 legend_title = "Sample", 
                 colors = sm_cols_v, 
                 label = FALSE, 
                 size = 0.2,
                 alpha = 0.5) + 
    theme(legend.key.size = unit(0.05, "in"), 
          legend.title = element_blank())
ggsave(paste0(dir_plot, "UMAP.sample.jpeg"), 
       width = 3.5, 
       height = 2.5, 
       units = "in")
ggsave(paste0(dir_plot, "UMAP.sample.pdf"), 
       width = 3.5, 
       height = 2.5, 
       units = "in")

plot_umap_labels(seur, 
                 colname = "Patient", 
                 legend_title = "Patient", 
                 colors = pt_cols_v,
                 label = FALSE, 
                 size = 0.3,
                 alpha = 0.5)
ggsave(paste0(dir_plot, "UMAP.patient.jpeg"), 
       width = 5, 
       height = 5, 
       units = "in")
ggsave(paste0(dir_plot, "UMAP.patient.pdf"), 
       width = 5, 
       height = 5, 
       units = "in")

plot_umap_labels(seur, 
                 colname = "Tumor", 
                 legend_title = "Type", 
                 colors = tum_colsv, 
                 label = FALSE, 
                 size = 0.3,
                 alpha = 0.5)
ggsave(paste0(dir_plot, "UMAP.tumor.jpeg"), 
       width = 5, 
       height = 5, 
       units = "in")

