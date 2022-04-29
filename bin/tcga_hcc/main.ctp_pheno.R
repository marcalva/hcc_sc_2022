
# test difference in cell type proportions between 
# tumor and non-tumor

setwd("../../")

library(NMF)
library(ggplot2)
library(edgeR)
library(Biobase)
library(plyr)
source("scripts/sc_func.R")
# source("scripts/reference_free.R")
# library(Bisque)

#==============================================================================
# read in data
#==============================================================================

fn <- "data/processed/tcga_hcc/expr/tcga.gencodev26.rds"
gencode <- readRDS(fn)

# expression data
fn <- "data/processed/tcga_hcc/expr/tcga.lihc.TMM.log.rin.rds"
tmm <- readRDS(fn)
rownames(tmm) <- gencode[ rownames(tmm), "Name"]

cl_id <- "cell_type_main"

# proportions
fn <- paste0("data/processed/tcga_hcc/ctp/tcga.TMM.", cl_id, ".decomp.rds")
ct.md <- readRDS(fn)
ct.mdp <- ct.md$bulk.props

for (ct in names(ct.md$markers)){
    ct.md$markers[[ct]] <- gencode[ct.md$markers[[ct]], "Name"]
}

# pheno data
fn <- "data/processed/tcga_hcc/sample/cases.hcc.361.merged.rds"
cases <- readRDS(fn)
fn <- "data/processed/tcga_hcc/sample/samples.hcc.410.merged.rds"
samples <- readRDS(fn)
fn <- "data/processed/tcga_hcc/sample/tcga.bio_sample.rds"
bio_sample <- readRDS(fn)

fn <- "data/processed/tcga_hcc/sample/samples.hcc.410.merged.rds"
s.p <- readRDS(fn) # sample phenotypes

# sample map
fn <- "data/processed/tcga_hcc/sample/sam_case_map.rds"
smaps <- readRDS(fn)
cid2sid_tum <- smaps[["cid2sid_tum"]]
cid2sid_nt <- smaps[["cid2sid_nt"]]
cid2sid_tum_pair <- smaps[["cid2sid_tum_pair"]]
cid2sid_nt_pair <- smaps[["cid2sid_nt_pair"]]

#========================================================
# paired t-test and wilcox
#========================================================

# scale mean 0 var 1
ct.mdp.s <- t(apply(ct.mdp, 1, scale))
colnames(ct.mdp.s) <- colnames(ct.mdp)
ct.mdp.s.df <- as.data.frame(t(ct.mdp.s))
ct.mdp.s.df <- ct.mdp.s.df[c(cid2sid_tum_pair, cid2sid_nt_pair),]
ct.mdp.s.df[,"Type"] <- "Primary Tumor"
ct.mdp.s.df[cid2sid_nt_pair, "Type"] <- "Solid Tissue Normal"

ct.mdp.pt <- ct.mdp.s[, cid2sid_tum_pair, drop=FALSE]
ct.mdp.norm <- ct.mdp.s[, cid2sid_nt_pair, drop=FALSE]

paired.test <- lapply(rownames(ct.mdp), function(ct){
                      prop.pt <- ct.mdp.pt[ct,]
                      prop.norm <- ct.mdp.norm[ct,]
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
rownames(paired.test) <- rownames(ct.mdp)
paired.test <- as.data.frame(paired.test)
paired.test[,"w.p_adj"] <- p.adjust(paired.test[,"w.p"], method = "fdr")
paired.test[,"t.p_adj"] <- p.adjust(paired.test[,"t.p"], method = "fdr")
paired.test <- paired.test[order(paired.test[,"t.estimate"], decreasing=TRUE),]
paired.test[,"CellType"] <- rownames(paired.test)

# output results
out_dir <- paste0("exp/tcga_hcc/ctp.", cl_id, "/")
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

fn <- paste0(out_dir, "tcga.", cl_id, ".tum_nontum.txt")
write.table(paired.test, fn, row.names=TRUE, col.names=NA, sep="\t", quote=FALSE)


