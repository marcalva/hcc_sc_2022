
setwd("../../")

library(Seurat)
library(ggplot2)
library(scales)
library(gridExtra)
source("scripts/standard_seurat.R")

#=========================================
# Functions
#=========================================

create_dir <- function(p){
    dir.create(p, showWarnings=FALSE, recursive=TRUE)
}

#=========================================
#=========================================

n_threads <- 8

# Gene data
gene_info <- read.table("data/ref/gencode26/gencode.v26.annotation.txt", 
                        header = TRUE, 
                        row.names = 1, 
                        stringsAsFactors = FALSE, 
                        sep = "\t")
gene_info[,"Name"] <- make.unique(gene_info[,"Name"])
symb2ens <- rownames(gene_info)
names(symb2ens) <- gene_info[,"Name"]

fn <- "data/processed/d_sct_cca_all/merge/seur.cca.rds"
nash_hcc <- readRDS(fn)

# Set directories
dir_out <- "data/processed/d_sct_cca_all/"; create_dir(dir_out)

#=========================================
# try integration
#=========================================

nash_hcc_l <- SplitObject(nash_hcc, split.by = "orig.ident")

# SCTransform each nash hcc
for (i in 1:length(nash_hcc_l)){
    message(i)
    DefaultAssay(nash_hcc_l[[i]]) <- "RNA"
    nash_hcc_l[[i]] <- SCTransform(nash_hcc_l[[i]])
    nash_hcc_l[[i]]@assays$integrated <- NULL
}

fn <- file.path(dir_out, "seur.sct.split.rds")
saveRDS(nash_hcc_l, fn)

