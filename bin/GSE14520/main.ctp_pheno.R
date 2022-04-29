
# Run ctp phenotype associations

setwd("../../")

library(NMF)
library(ggplot2)
library(edgeR)
library(Biobase)
library(plyr)
source("scripts/sc_func.R")
# source("scripts/reference_free.R")
# library(Bisque)


dir_data <- "data/processed/GSE14520/"
dir.create(dir_data, recursive=TRUE, showWarning=FALSE)

fn <- paste0(dir_data, "expr.RMA_log2.gmean.rds")
ex <- readRDS(fn)

fn <- paste0(dir_data, "sdata.txt")
sdat <- read.table(fn, row.names=1, header=TRUE, sep="\t", stringsAsFactors=FALSE)

fn <- paste0(dir_data, "geo_pheno.txt")
geo_pheno <- read.table(fn, row.names=1, header=TRUE, sep="\t", stringsAsFactors=FALSE)

cl_id <- "cell_type_main"

# marker data
fn <- paste0("exp/sharma_aiz/markers/markers.", cl_id, ".txt")
markers.celltype <- read.table(fn, header=TRUE, stringsAsFactors=FALSE)

fn <- paste0("data/processed/GSE14520/ctp/", "LCI.", cl_id, ".decomp.rds")
ex.ct.md <- readRDS(fn)

ex.ct.mdp <- ex.ct.md$bulk.props

# sample map
fn <- "data/processed/GSE14520/sam_case_map.rds"
smaps <- readRDS(fn)
cid2sid_tum <- smaps[["cid2sid_tum"]]
cid2sid_nt <- smaps[["cid2sid_nt"]]
cid2sid_tum_pair <- smaps[["cid2sid_tum_pair"]]
cid2sid_nt_pair <- smaps[["cid2sid_nt_pair"]]

#========================================================
# paired t-test and wilcox
#========================================================

# scale mean 0 var 1
ex.ct.mdp.s <- t(apply(ex.ct.mdp, 1, scale))
colnames(ex.ct.mdp.s) <- colnames(ex.ct.mdp)
ex.ct.mdp.s.df <- as.data.frame(t(ex.ct.mdp.s))
ex.ct.mdp.s.df <- ex.ct.mdp.s.df[c(unname(cid2sid_tum_pair), unname(cid2sid_nt_pair)),]
ex.ct.mdp.s.df[,"Type"] <- sdat[rownames(ex.ct.mdp.s.df), "Tissue.Type"]

ex.ct.mdp.pt <- ex.ct.mdp.s[,unname(cid2sid_tum_pair),drop=FALSE]
ex.ct.mdp.norm <- ex.ct.mdp.s[,unname(cid2sid_nt_pair),drop=FALSE]

# for each cell type
paired.test <- lapply(rownames(ex.ct.mdp), function(ct){
                      prop.pt <- ex.ct.mdp.pt[ct,]
                      prop.norm <- ex.ct.mdp.norm[ct,]
                      ret.1 <- wilcox.test(x = prop.pt, y = prop.norm, paired = TRUE)
                      ret.2 <- t.test(x = prop.pt, y = prop.norm, paired = TRUE)
                      ret <- c("w.statistic" = as.numeric(ret.1$statistic), 
                               "w.estimate" = as.numeric(ret.1$estimate), 
                               "w.p" = as.numeric(ret.1$p.value),
                               "t.statistic" = as.numeric(ret.2$statistic), 
                               "t.estimate" = as.numeric(ret.2$estimate), 
                               "t.p" = as.numeric(ret.2$p.value))
                      return(ret) })
paired.test <- do.call(rbind, paired.test)
rownames(paired.test) <- rownames(ex.ct.mdp)
paired.test <- as.data.frame(paired.test)
paired.test[,"w.p_adj"] <- p.adjust(paired.test[,"w.p"], method = "fdr")
paired.test[,"t.p_adj"] <- p.adjust(paired.test[,"t.p"], method = "fdr")
paired.test <- paired.test[order(paired.test[,"t.estimate"], decreasing=TRUE),]
paired.test[,"CellType"] <- rownames(paired.test)

# output results
out_dir <- paste0("exp/GSE14520/ctp.", cl_id, "/")
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

fn <- paste0(out_dir, "lci.", cl_id, ".tum_nontum.txt")
write.table(paired.test, fn, row.names=TRUE, col.names=NA, sep="\t", quote=FALSE)

