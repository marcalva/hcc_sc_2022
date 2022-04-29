
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
#=========================================

# cell cycle scoring
DefaultAssay(nash_hcc) <- "RNA"
nash_hcc <- NormalizeData(nash_hcc, scale.factor=1000)
s.genes <- cc.genes$s.genes
s.genes <- symb2ens[s.genes]
s.genes <- s.genes[s.genes %in% rownames(nash_hcc)]
g2m.genes <- cc.genes$g2m.genes
g2m.genes <- symb2ens[g2m.genes]
g2m.genes <- g2m.genes[s.genes %in% rownames(nash_hcc)]
nash_hcc <- CellCycleScoring(nash_hcc, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
DefaultAssay(nash_hcc) <- "SCT"


nash_hcc_l <- SplitObject(nash_hcc, split.by = "orig.ident")

# SCTransform each nash hcc
for (i in 1:length(nash_hcc_l)){
    message(i)
    DefaultAssay(nash_hcc_l[[i]]) <- "RNA"
    nash_hcc_l[[i]] <- SCTransform(nash_hcc_l[[i]], variable.features.n = 3000, 
                                   vars.to.regress = c("S.Score", "G2M.Score"),
                                   verbose=FALSE)
    nash_hcc_l[[i]]@assays$integrated <- NULL
}

fn <- file.path(dir_out, "seur.sct.cc.split.rds")
saveRDS(nash_hcc_l, fn)

