
# Run edgeR on bulk RNA-seq from TCGA cohort

setwd("../../")

#' Save a data frame to file
save2file <- function(datf, 
                      dirout = "./", 
                      fn = "datf.txt", 
                      gz = FALSE, 
                      row.names = TRUE, 
                      col.names = NA, 
                      quote = FALSE, 
                      sep = "\t"){
    dir.create(dirout, showWarnings = FALSE, recursive = TRUE)
    fn <- file.path(dirout, fn)
    if (gz){
        gz1 <- gzfile(fn, "w")
        write.table(datf, 
                    gz1, 
                    row.names = row.names, 
                    col.names = col.names, 
                    quote = quote, 
                    sep = sep)
        close(gz1)
    } else {
        write.table(datf, 
                    fn, 
                    row.names = row.names, 
                    col.names = col.names, 
                    quote = quote, 
                    sep = sep)
    }
}

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

library(edgeR)
source("scripts/rnaseq_scripts/R/de.R")

# read in data
out_dir <- "data/processed/tcga_hcc/expr/"

fn <- paste0(out_dir, "tcga.gencodev26.rds")
gencode <- readRDS(fn)

# expression data
fn <- paste0(out_dir, "tcga.lihc.TMM.rds")
tmm <- readRDS(fn)
rownames(tmm) <- gencode[ rownames(tmm), "Name"]

fn <- paste0(out_dir, "tcga.lihc.htseq_counts.rds")
counts <- readRDS(fn)

# pheno data
out_dir <- "data/processed/tcga_hcc/sample/"
fn <- paste0(out_dir, "cases.hcc.361.merged.rds")
cases <- readRDS(fn)
fn <- paste0(out_dir, "samples.hcc.410.merged.rds")
samples <- readRDS(fn)

# sample map
fn <- "data/processed/tcga_hcc/sample/sam_case_map.rds"
smaps <- readRDS(fn)
cid2sid_tum <- smaps[["cid2sid_tum"]]
cid2sid_nt <- smaps[["cid2sid_nt"]]
cid2sid_tum_pair <- smaps[["cid2sid_tum_pair"]]
cid2sid_nt_pair <- smaps[["cid2sid_nt_pair"]]
sid_all <- c(cid2sid_tum, cid2sid_nt)
sid_both <- c(cid2sid_tum_pair, cid2sid_nt_pair)

#=====================================================
# paired DEP
#=====================================================

# subset counts and pheno
samples.k <- samples[sid.both,]
counts.k <- counts[,sid.both]

# filter by genes
gk1 <- rowSums(counts.k >= 1) >= (0.9 * ncol(counts.k))
gk2 <- rowMeans(counts.k) >= 10

counts.k <- counts.k[gk2,]

# design matrix
samples.k[,"Sample Type"] <- factor(samples.k[,"Sample Type"], 
                                         levels = c("Solid Tissue Normal", 
                                                    "Primary Tumor"))
design.k <- model.matrix( ~ samples.k[,"Case ID"] + samples.k[,"Sample Type"])

# edgeR
y <- DGEList(counts = counts.k, genes = gencode[rownames(counts.k),])
y <- calcNormFactors(y, method = "TMM")
y <- estimateDisp(y, design.k, robust=TRUE)
prior.count <- 1
p.adjust <- "fdr"
ret <- runedger(y, design.k, prior.count = prior.count, p.adjust = p.adjust)

lrt <- ret
diff_table <- lrt$diffTable
diff_table <- cbind(gencode[rownames(diff_table),], diff_table)

# export
out_dir <- "exp/tcga_hcc/de/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# save LRT
out_fn <- paste0(out_dir, "lrt.rds")
saveRDS(lrt, out_fn)

# save table
out_fn <- "de.table.tsv.gz"
save2file(diff_table, dirout = out_dir, fn = out_fn)

