
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

fn <- "data/processed/sharma2019/sharma.seur.rds"
sharma <- readRDS(fn)

# Set directories
dir_out <- "data/processed/sharma2019/"; create_dir(dir_out)
# dir_plot <- "exp/d_sct_cca_all/cca_init/"; create_dir(dir_plot)

dir_exp <- "exp/sharma2019/"; create_dir(dir_exp)

#=========================================
# SCTransform each sample
#=========================================

sharma_l <- SplitObject(sharma, split.by = "patientno")

for (i in 1:length(sharma_l)){
    message(i)
    DefaultAssay(sharma_l[[i]]) <- "RNA"
    sharma_l[[i]] <- SCTransform(sharma_l[[i]], variable.features.n = 3000, 
                                 verbose=FALSE)
}

outfn <- file.path(dir_out, "sharma_sct_l.rds")
saveRDS(sharma_l, outfn)

