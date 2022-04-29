
setwd("../../")

library(ggplot2)
library(gridExtra)

# raw data
out_dir <- "data/processed/tcga_hcc/expr/"
fn <- paste0(out_dir, "tcga.lihc.htseq_counts.rds")
counts <- readRDS(fn)
fn <- paste0(out_dir, "tcga.gencodev26.rds")
gencode <- readRDS(fn)

# sample data
out_dir <- "data/processed/tcga_hcc/sample/"
fn <- paste0(out_dir, "tcga.clinical.rds")
clinical <- readRDS(fn)
fn <- paste0(out_dir, "tcga.sample_sheet.rds")
sample_sheet <- readRDS(fn)

clinical[,"tumor_stage"] <- clinical[,"ajcc_tumor_pathologic_pt"]
clinical[,"primary_diagnosis"] <- clinical[,"histologic_diagnosis"]

samples2case <- sample_sheet[,"Case ID"]
names(samples2case) <- sample_sheet[,"Sample ID"]

# norm counts
out_dir <- "data/processed/tcga_hcc/expr/"
fn <- paste0(out_dir, "tcga.lihc.TMM.rds")
tmmm <- readRDS(fn)

# Calculate PCA on the TMM expression values

tmmm.norm <- log10(tmmm + 1)
tmmm.norm <- t(tmmm.norm)
tmmm.norm <- apply(tmmm.norm, 2, scale)
rownames(tmmm.norm) <- colnames(tmmm)

tmmm.svd <- svd(tmmm.norm)
tmmm.pcs <- tmmm.svd$u
rownames(tmmm.pcs) <- rownames(tmmm.norm)
colnames(tmmm.pcs) <- paste0("PC", as.numeric(1:nrow(tmmm.norm)))

d.var <- tmmm.svd$d^2 / (ncol(tmmm) - 1)
d.ve <- d.var / sum(d.var)

# plots
exp_dir <- "exp/tcga_hcc/expr/"
dir.create(exp_dir, showWarnings = FALSE, recursive = TRUE)

# plot variance explained
datf <- data.frame("PC" = 1:length(d.ve), 
                   "VE" = d.ve)
p <- ggplot(datf[1:50,], aes(x = PC, y = VE)) + 
    geom_point(shape=16) + theme_minimal() + 
    xlab("PC") + ylab("Var. explained")
fn <- paste0(exp_dir, "var_expl.tmm.jpeg")
ggsave(fn, width = 6, height = 4, dpi=300)

# plot PCs colored by normal vs. tumor
datf <- data.frame(as.data.frame(tmmm.pcs[,1:10]), 
                   "Type" = sample_sheet[rownames(tmmm.pcs),"Sample Type"], check.names=FALSE)
p1 <- ggplot(datf, aes(x=PC1,y=PC2,color=Type)) + 
    geom_point(shape=16) + theme_minimal()
p2 <- ggplot(datf, aes(x=PC2,y=PC3,color=Type)) + 
    geom_point(shape=16) + theme_minimal()
p3 <- ggplot(datf, aes(x=PC3,y=PC4,color=Type)) + 
    geom_point(shape=16) + theme_minimal()
g <- arrangeGrob(p1, p2, p3, nrow=1, ncol=3)
fn <- paste0(exp_dir, "pcs.type.tmm.jpeg")
ggsave(fn, g, width = 15, height=4)

# plot PCs colored by stage
lowc <- "lightgrey"
highc <- "firebrick"

tumor <- rownames(sample_sheet)[sample_sheet[,"Sample Type"] == "Primary Tumor"]
tumor <- intersect(tumor, rownames(tmmm.pcs))
tmm.pcs.tum <- tmmm.pcs[tumor,]
rownames(tmm.pcs.tum) <- samples2case[rownames(tmm.pcs.tum)]
datf <- data.frame(as.data.frame(tmmm.pcs[,1:10]), 
                   "Type" = sample_sheet[rownames(tmmm.pcs),"Sample Type"], check.names=FALSE)
ck <- c("age_at_diagnosis", "ethnicity", "gender", "race", "tumor_stage", "primary_diagnosis")
clinical.sub <- clinical[,ck]
nak <- clinical.sub[,"age_at_diagnosis"] == "'--"
clinical.sub[nak,"age_at_diagnosis"] <- NA
clinical.sub[,"age_at_diagnosis"] <- as.numeric(clinical.sub[,"age_at_diagnosis"])

nak <- clinical.sub[,"tumor_stage"] == "[Discrepancy]" | 
    clinical.sub[,"tumor_stage"] == "[Not Available]"
clinical.sub[nak,"tumor_stage"] <- NA
stage_nn <- sort(unique(clinical.sub[,"tumor_stage"]))
stage_n <- 1:length(stage_nn)
names(stage_n) <- stage_nn
clinical.sub[,"tumor_stage_n"] <- stage_n[clinical.sub[,"tumor_stage"]]

datf <- cbind(tmm.pcs.tum[,1:10], clinical.sub[rownames(tmm.pcs.tum),])
# age
p1 <- ggplot(datf, aes(x=PC1,y=PC2,color=age_at_diagnosis)) + 
    geom_point(shape=16) + theme_minimal() + 
    scale_color_continuous(low=lowc, high=highc)
p2 <- ggplot(datf, aes(x=PC3,y=PC4,color=age_at_diagnosis)) + 
    geom_point(shape=16) + theme_minimal() + 
    scale_color_continuous(low=lowc, high=highc)
p3 <- ggplot(datf, aes(x=PC5,y=PC6,color=age_at_diagnosis)) + 
    geom_point(shape=16) + theme_minimal() + 
    scale_color_continuous(low=lowc, high=highc)
# ethnic
p4 <- ggplot(datf, aes(x=PC1,y=PC2,color=ethnicity)) + 
    geom_point(shape=16) + theme_minimal()
p5 <- ggplot(datf, aes(x=PC3,y=PC4,color=ethnicity)) + 
    geom_point(shape=16) + theme_minimal()
p6 <- ggplot(datf, aes(x=PC5,y=PC6,color=ethnicity)) + 
    geom_point(shape=16) + theme_minimal()

# gender
p7 <- ggplot(datf, aes(x=PC1,y=PC2,color=gender)) + 
    geom_point(shape=16) + theme_minimal()
p8 <- ggplot(datf, aes(x=PC3,y=PC4,color=gender)) + 
    geom_point(shape=16) + theme_minimal()
p9 <- ggplot(datf, aes(x=PC5,y=PC6,color=gender)) + 
    geom_point(shape=16) + theme_minimal()

g <- arrangeGrob(p1, p2, p3, p4, p5, p6, p7, p8, p9, nrow=3, ncol=3)

fn <- paste0(exp_dir, "pcs.pheno.tmm.jpeg")
ggsave(fn, g, width = 15, height=12)

# tumor_stage_n
p1 <- ggplot(datf, aes(x=PC1,y=PC2,color=tumor_stage_n)) + 
    geom_point(shape=16) + theme_minimal() + 
    scale_color_continuous(low=lowc, high=highc)
p2 <- ggplot(datf, aes(x=PC3,y=PC4,color=tumor_stage_n)) + 
    geom_point(shape=16) + theme_minimal() + 
    scale_color_continuous(low=lowc, high=highc)
p3 <- ggplot(datf, aes(x=PC5,y=PC6,color=tumor_stage_n)) + 
    geom_point(shape=16) + theme_minimal() + 
    scale_color_continuous(low=lowc, high=highc)
# tumor_stage
p4 <- ggplot(datf, aes(x=PC1,y=PC2,color=tumor_stage)) + 
    geom_point(shape=16) + theme_minimal()
p5 <- ggplot(datf, aes(x=PC3,y=PC4,color=tumor_stage)) + 
    geom_point(shape=16) + theme_minimal()
p6 <- ggplot(datf, aes(x=PC5,y=PC6,color=tumor_stage)) + 
    geom_point(shape=16) + theme_minimal()

g <- arrangeGrob(p1, p2, p3, p4, p5, p6, nrow=2, ncol=3)

fn <- paste0(exp_dir, "pcs.stage.tmm.jpeg")
ggsave(fn, g, width = 15, height=8)

# primary diagnosis
# primary diag
datf[,"primary_diagnosis"] <- sapply(datf[,"primary_diagnosis"], function(x){
    paste0(strwrap(x, width=30), collapse="\n") })
p <- ggplot(datf, aes(x=PC1,y=PC2,color=primary_diagnosis)) + 
    geom_point(shape=16) + theme_minimal()
fn <- paste0(exp_dir, "pcs.diag.tmm.jpeg")
ggsave(fn, width = 6, height=4)


