
# Relationship between cell type props and molecular features

setwd("../../")

library(NMF)
library(ggplot2)
library(edgeR)
library(Biobase)
library(plyr)
library(RColorBrewer)

#========================================================
# Functions
#========================================================

runedger <- function(yf, ds, prior.count = 0.125, p.adjust = "fdr"){
    require(edgeR)
    fit <- glmFit(yf, ds, prior.count = prior.count, coef = ncol(yf$design))
    lrt <- glmLRT(fit, coef = ncol(fit$design))
    diffTable <- lrt$table

    diffTable$p_adj <- p.adjust(diffTable$PValue, method = p.adjust)
    diffTable$FoldChange <- 2^(diffTable$logFC)
    diffTable$AverageCPM <- 2^(diffTable$logCPM)
    diffTable <- diffTable[,c("logFC", "FoldChange", "logCPM", "AverageCPM", "LR", "PValue", "p_adj")]
    diffTable <- diffTable[order(diffTable$PValue),]
    return(list("diffTable" = diffTable, "fit" = lrt))
}

#========================================================
#========================================================

fn <- "data/processed/tcga_hcc/expr/tcga.gencodev26.rds"
gencode <- readRDS(fn)

# expression data
fn <- "data/processed/tcga_hcc/expr/tcga.lihc.TMM.log.rin.rds"
tmm <- readRDS(fn)
g_exprd <- rownames(tmm)

fn <- "data/processed/tcga_hcc/expr/tcga.lihc.htseq_counts.rds"
counts <- readRDS(fn)

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
dir_exp <- "exp/tcga_hcc/mut_expr/"
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

#========================================================
# association with mutations
#========================================================

# get top 25 mutated genes
o <- order(rowSums(mut.any), decreasing=TRUE)
genes <- rownames(mut.any)[o][1:25]
# genes <- c("TP53", "RB1")

for (gene in genes){
    message(gene)

    counts.pt <- counts[,unname(cid2sid_tum)]

    k <- intersect(colnames(mut.any), colnames(counts.pt))
    mut.anyk <- mut.any[,k]
    counts.ptk <- counts.pt[,k]

    k0 <- colnames(mut.anyk)[mut.anyk[gene,] == 0]
    k1 <- colnames(mut.anyk)[mut.anyk[gene,] == 1]

    # filter by genes
    # gk2 <- rowMeans(counts.ptk) >= 10
    counts.ptk <- counts.ptk[g_exprd,]

    design.k <- model.matrix( ~ t(mut.anyk[gene,]) )
    colnames(design.k) <- c("Intercept", gene)

    # edgeR
    y <- DGEList(counts = counts.ptk, genes = gencode[rownames(counts.ptk),])
    y <- calcNormFactors(y, method = "TMM")
    y <- estimateDisp(y, design.k, robust=TRUE)
    prior.count <- 1
    p.adjust <- "fdr"
    ret <- runedger(y, design.k, prior.count = prior.count, p.adjust = p.adjust)

    lrt <- ret
    diff_table <- lrt$diffTable
    diff_table <- cbind(gencode[rownames(diff_table),], diff_table)

    # export
    out_dir <- "exp/mut_expr/de/"
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

    # save LRT
    out_fn <- paste0(out_dir, gene, ".lrt.rds")
    saveRDS(lrt, out_fn)

    # save table
    out_fn <- paste0(out_dir, gene, ".de.table.tsv.gz")
    gz1 <- gzfile(out_fn, "w")
    write.table(diff_table, 
                gz1, 
                row.names = FALSE, 
                col.names = TRUE, 
                quote = FALSE, 
                sep = '\t')
    close(gz1)
}
