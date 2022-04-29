
# Relationship between cell type props and molecular features

setwd("../../")

library(NMF)
library(ggplot2)
library(edgeR)
library(Biobase)
library(plyr)
library(RColorBrewer)
source("scripts/sc_func.R")
# source("scripts/reference_free.R")
# library(Bisque)

#========================================================
# Functions
#========================================================

wt.test <- function(x, y, paired = FALSE){
    ret.1 <- wilcox.test(x = x, y = y, paired = paired)
    ret.2 <- t.test(x = x, y = y, paired = paired)
    ret <- c("w.statistic" = as.numeric(ret.1$statistic), 
             "w.estimate" = as.numeric(ret.1$estimate), 
             "w.p" = as.numeric(ret.1$p.value),
             "t.statistic" = as.numeric(ret.2$statistic), 
             "t.estimate" = as.numeric(ret.2$estimate), 
             "t.p" = as.numeric(ret.2$p.value))
    return(ret)
}

#========================================================
#========================================================

ct_id <- "cell_type_main"

fn <- "data/processed/tcga_hcc/expr/tcga.gencodev26.rds"
gencode <- readRDS(fn)

# expression data
fn <- "data/processed/tcga_hcc/expr/tcga.lihc.TMM.log.rin.rds"
tmm <- readRDS(fn)
rownames(tmm) <- gencode[ rownames(tmm), "Name"]

# proportions
fn <- paste0("data/processed/tcga_hcc/ctp/tcga.TMM.", ct_id, ".decomp.rds")
ct.md <- readRDS(fn)
ct.mdp <- ct.md$bulk.props

for (ct in names(ct.md$genes.used)) 
    ct.md$genes.used[[ct]] <- gencode[ct.md$genes.used[[ct]], "Name"]

# pheno data
fn <- "data/processed/tcga_hcc/sample/cases.hcc.361.merged.rds"
cases <- readRDS(fn)
fn <- "data/processed/tcga_hcc/sample/samples.hcc.410.merged.rds"
samples <- readRDS(fn)
fn <- "data/processed/tcga_hcc/sample/tcga.bio_sample.rds"
bio_sample <- readRDS(fn)

# molecular data
fn <- "data/raw/tcga_hcc/gdac.LIHC.aggregate/LIHC-TP.samplefeatures.txt"
# mol <- read.table(fn, header = TRUE, stringsAsFactors = FALSE, sep = "\t", row.names=1)

fn <- "data/processed/tcga_hcc/cnv/tcga.lihc.thr.del.354_tum.txt"
dels <- read.table(fn, header = TRUE, stringsAsFactors = FALSE, sep = "\t", row.names=1, check.names=FALSE)
fn <- "data/processed/tcga_hcc/cnv/tcga.lihc.thr.amp.354_tum.txt"
amps <- read.table(fn, header = TRUE, stringsAsFactors = FALSE, sep = "\t", row.names=1, check.names=FALSE)

cnvs <- rbind(amps, dels)
cnvs[cnvs >= 1] <- 1

fn <- "data/processed/tcga_hcc/cnv/tcga.lihc.cnv_descr.txt"
cnv_desc <- read.table(fn, header = TRUE, stringsAsFactors = FALSE, sep = "\t", row.names=1, check.names=FALSE)

fn <- "data/processed/tcga_hcc/mut/tcga.lihc.mut.gene.any.txt"
mut.any <- read.table(fn, header = TRUE, stringsAsFactors = FALSE, sep = "\t", row.names=1, check.names=FALSE)
fn <- "data/processed/tcga_hcc/mut/tcga.lihc.mut.gene.type.txt"
mut.type <- read.table(fn, header = TRUE, stringsAsFactors = FALSE, sep = "\t", row.names=1, check.names=FALSE)

mut.any[mut.any >= 1] <- 1

# Output directory
dir_exp <- "exp/tcga_hcc/mol/"
dir.create(dir_exp, showWarnings = FALSE, recursive = TRUE)

#========================================================
# Map IDs, get tumor samples
#========================================================

# sample map
fn <- "data/processed/tcga_hcc/sample/sam_case_map.rds"
smaps <- readRDS(fn)
cid2sid_tum <- smaps[["cid2sid_tum"]]
cid2sid_nt <- smaps[["cid2sid_nt"]]
cid2sid_tum_pair <- smaps[["cid2sid_tum_pair"]]
cid2sid_nt_pair <- smaps[["cid2sid_nt_pair"]]
sid_all <- c(cid2sid_tum, cid2sid_nt)
sid_both <- c(cid2sid_tum_pair, cid2sid_nt_pair)

ct.mdp.pt <- ct.mdp[,cid2sid_tum]

#========================================================
# association with mutations
#========================================================

# interesect mutations with ctp/pheno data

k <- intersect(colnames(mut.any), colnames(ct.mdp.pt))
mut.anyk <- mut.any[,k]
ct.mdp.ptk <- ct.mdp.pt[,k]

cell_types <- rownames(ct.mdp.pt)
genes.mut <- rownames(mut.any)
wr <- lapply(genes.mut, function(g){
             ctr <- lapply(cell_types, function(ct){
                           k0 <- colnames(mut.anyk)[mut.anyk[g,] == 0]
                           k1 <- colnames(mut.anyk)[mut.anyk[g,] == 1]
                           ret <- wt.test(x = ct.mdp.ptk[ct, k1], 
                                          y = ct.mdp.ptk[ct, k0])
                           return(ret) })
             ctr <- as.data.frame(do.call(rbind, ctr), stringsAsFactors=FALSE)
             ctr[,"CellType"] <- cell_types
             ctr[,"Gene"] <- g
             return(ctr) })
wr <- do.call(rbind, wr)

wr[,"w.p_adj"] <- p.adjust(wr[,"w.p"], method="fdr")
wr[,"t.p_adj"] <- p.adjust(wr[,"t.p"], method="fdr")
rownames(wr) <- paste(wr[,"CellType"], wr[,"Gene"], sep="_")

o <- order(wr[,"w.p"], decreasing=FALSE)
wr <- wr[o,]

fn <- paste0(dir_exp, ct_id, ".SMG_mutsig.2CV.ctp.txt")
write.table(wr, fn, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

k <- wr[,"w.p_adj"] < 0.05
wrk <- wr[k,]
message("Number of unique genes: ", length(unique(wrk[,"Gene"])))
message("Number of unique cell-types: ", length(unique(wrk[,"CellType"])))

#========================================================
# associations with mutation types
#========================================================

# interesect mutations with ctp/pheno data

types.all <- unique(as.character(unlist(mut.type)))
types.all <- setdiff(types.all, "None")

k <- intersect(colnames(mut.type), colnames(ct.mdp.pt))
mut.typek <- mut.type[,k]
ct.mdp.ptk <- ct.mdp.pt[,k]

cell_types <- rownames(ct.mdp.pt)
genes.mut <- rownames(mut.type)
wr.type <- lapply(genes.mut, function(g){
             ctr <- lapply(cell_types, function(ct){
                           ty.l <- lapply(types.all, function(ty){
                                         k0 <- colnames(mut.typek)[mut.typek[g,] == "None"]
                                         k1 <- colnames(mut.typek)[mut.typek[g,] == ty]
                                         if (length(k0) < 2 | length(k1) < 2) return(NULL)
                                         ret <- wt.test(x = ct.mdp.ptk[ct, k1], 
                                                        y = ct.mdp.ptk[ct, k0])
                                         ret <- as.data.frame(t(ret))
                                         ret[,"MutType"] <- ty
                                         ret[,"CellType"] <- ct
                                         ret[,"Gene"] <- g
                                         ret[,"n.1"] <- length(k1)
                                         ret[,"n.2"] <- length(k0)
                                         return(ret) })
                           ty <- do.call(rbind, ty.l)
                           return(ty) })
             ctr <- as.data.frame(do.call(rbind, ctr), stringsAsFactors=FALSE)
             return(ctr) })
wr.type <- do.call(rbind, wr.type)

wr.type[,"w.p_adj"] <- p.adjust(wr.type[,"w.p"], method="fdr")
wr.type[,"t.p_adj"] <- p.adjust(wr.type[,"t.p"], method="fdr")
rownames(wr.type) <- paste(wr.type[,"CellType"], wr.type[,"Gene"], wr.type[,"MutType"], sep="_")

o <- order(wr.type[,"w.p"], decreasing=FALSE)
wr.type <- wr.type[o,]

fn <- paste0(dir_exp, ct_id, ".SMG_mutsig.2CV.mut_type.ctp.txt")
write.table(wr.type, fn, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

k <- wr.type[,"w.p_adj"] < 0.05
wr.typek <- wr.type[k,]
message("Number of unique genes: ", length(unique(wr.typek[,"Gene"])))
message("Number of unique cell-types: ", length(unique(wr.typek[,"CellType"])))

#========================================================
# associations with CNVs
#========================================================

# interesect CNVs with ctp/pheno data

k <- intersect(colnames(cnvs), colnames(ct.mdp.pt))
cnvsk <- cnvs[,k]
ct.mdp.ptk <- ct.mdp.pt[,k]

cell_types <- rownames(ct.mdp.pt)
cnv_regs <- rownames(cnvsk)
wr.cnv <- lapply(cnv_regs, function(g){
             ctr <- lapply(cell_types, function(ct){
                           k0 <- colnames(cnvsk)[cnvsk[g,] == 0]
                           k1 <- colnames(cnvsk)[cnvsk[g,] == 1]
                           ret <- wt.test(x = ct.mdp.ptk[ct, k1], 
                                          y = ct.mdp.ptk[ct, k0])
                           return(ret) })
             ctr <- as.data.frame(do.call(rbind, ctr), stringsAsFactors=FALSE)
             ctr[,"CellType"] <- cell_types
             ctr[,"Region"] <- g
             return(ctr) })
wr.cnv <- do.call(rbind, wr.cnv)

wr.cnv[,"w.p_adj"] <- p.adjust(wr.cnv[,"w.p"], method="fdr")
wr.cnv[,"t.p_adj"] <- p.adjust(wr.cnv[,"t.p"], method="fdr")

o <- order(wr.cnv[,"w.p"], decreasing=FALSE)
wr.cnv <- wr.cnv[o,]

wr.cnv[,"n_genes"] <- cnv_desc[wr.cnv[,"Region"], "n_genes"]
wr.cnv[,"Genes"] <- cnv_desc[wr.cnv[,"Region"], "genes"]

fn <- paste0(dir_exp, ct_id, ".tcga.lihc.cnv.ctp.txt")
write.table(wr.cnv, fn, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

