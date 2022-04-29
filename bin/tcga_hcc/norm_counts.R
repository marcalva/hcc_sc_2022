
setwd("../../")

library(edgeR)

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
        fn <- paste0(fn, ".gz")
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

# read in data
out_dir <- "data/processed/tcga_hcc/expr/"
fn <- paste0(out_dir, "tcga.lihc.htseq_counts.rds")
counts <- readRDS(fn)
fn <- paste0(out_dir, "tcga.gencodev26.rds")
gencode <- readRDS(fn)

out_dir <- "data/processed/tcga_hcc/sample/"
fn <- paste0(out_dir, "tcga.clinical.rds")
clinical <- readRDS(fn)
fn <- paste0(out_dir, "tcga.sample_sheet.rds")
sample_sheet <- readRDS(fn)
fn <- "data/processed/tcga_hcc/sample/tcga.bio_sample.rds"
bio_sample <- readRDS(fn)

out_dir <- "data/processed/tcga_hcc/expr/"
save2file(counts, dirout=out_dir, fn = "tcga.lihc.htseq_counts.txt", 
          gz=TRUE)
save2file(gencode, dirout=out_dir, fn = "tcga.gencodev26.txt", 
          gz=TRUE)

#===============================================================
# subset counts to only those with HCC (removes 10 cases)
# subset to samples that are primary tumor or normal tissue
# remove recurrent as these are from same person
#===============================================================

cases.k <- rownames(clinical)[clinical[,"histologic_diagnosis"] == "Hepatocellular Carcinoma"]
clinical <- clinical[cases.k,]

bio_sample <- bio_sample[colnames(counts),]

k <- bio_sample[,"case"] %in% cases.k
bio_sample <- bio_sample[k,]

k <- bio_sample[,"sample_type"] != "Recurrent Tumor"
bio_sample <- bio_sample[k,]

samples.k <- bio_sample[,"bcr_sample_barcode"]
cases.k <- unique(bio_sample[,"case"])

counts <- counts[,samples.k]
clinical <- clinical[cases.k,]
sample_sheet <- sample_sheet[samples.k,]

#===============================================================
#===============================================================

# filter at least 50% of samples have greater than 0 counts
sfrac <- 0.5
min_c <- 0

y <- DGEList(counts = counts, genes = gencode)
keep <- rowSums(cpm(y) > min_c) >= (sfrac * ncol(counts))
y <- y[keep, , keep.lib.sizes=FALSE]

cpmm <- cpm(y, normalized.lib.sizes=TRUE)
y <- calcNormFactors(y, method = "upperquartile")
uqm <- cpm(y, normalized.lib.sizes=TRUE)
y <- calcNormFactors(y, method = "TMM")
tmmm <- cpm(y, normalized.lib.sizes=TRUE)

# output counts
out_dir <- "data/processed/tcga_hcc/expr/"

fn <- paste0(out_dir, "tcga.lihc.cpm.rds")
saveRDS(cpmm, fn)
save2file(cpmm, dirout=out_dir, fn = "tcga.lihc.cpm.txt", 
          gz=TRUE)
fn <- paste0(out_dir, "tcga.lihc.UQ.rds")
saveRDS(uqm, fn)
save2file(uqm, dirout=out_dir, fn = "tcga.lihc.UQ.txt", 
          gz=TRUE)
fn <- paste0(out_dir, "tcga.lihc.TMM.rds")
saveRDS(tmmm, fn)
save2file(tmmm, dirout=out_dir, fn = "tcga.lihc.TMM.txt", 
          gz=TRUE)

#===============================================================
# log and correct for RIN in all samples
#===============================================================

fn <- "data/processed/tcga_hcc/sample/tcga.analyte_rna.rds"
analyte.rna <- readRDS(fn)
analyte.rna <- analyte.rna[colnames(tmmm),]

# log transform
tmmm.n <- log10(tmmm + 1)

# correct for rin
tmmm.n <- t(tmmm.n)
tmmm.n <- apply(tmmm.n, 2, function(x){
                resid(lm(x ~ analyte.rna[,"rinvalue"])) })
tmmm.n <- t(tmmm.n)
dimnames(tmmm.n) <- dimnames(tmmm)

fn <- paste0(out_dir, "tcga.lihc.TMM.log.rin.rds")
saveRDS(tmmm.n, fn)
save2file(tmmm.n, dirout=out_dir, fn = "tcga.lihc.TMM.log.rin.txt", 
          gz=TRUE)

print(dim(tmmm.n))

#===============================================================
# subset to primary tumor
#===============================================================

k <- rownames(bio_sample)[bio_sample[,"sample_type"] == "Primary Tumor"]
y <- y[, k]


cpmm <- cpm(y, normalized.lib.sizes=TRUE)
y <- calcNormFactors(y, method = "upperquartile")
uqm <- cpm(y, normalized.lib.sizes=TRUE)
y <- calcNormFactors(y, method = "TMM")
tmmm <- cpm(y, normalized.lib.sizes=TRUE)

# output counts
out_dir <- "data/processed/tcga_hcc/expr/"

fn <- paste0(out_dir, "tcga.pt.lihc.cpm.rds")
saveRDS(cpmm, fn)
fn <- paste0(out_dir, "tcga.pt.lihc.UQ.rds")
saveRDS(uqm, fn)
fn <- paste0(out_dir, "tcga.pt.lihc.TMM.rds")
saveRDS(tmmm, fn)

#===============================================================
# log and correct for RIN in pt samples
#===============================================================

fn <- "data/processed/tcga_hcc/sample/tcga.analyte_rna.rds"
analyte.rna <- readRDS(fn)
analyte.rna <- analyte.rna[colnames(tmmm),]

# log transform
tmmm.n <- log10(tmmm + 1)

# correct for rin
tmmm.n <- t(tmmm.n)
tmmm.n <- apply(tmmm.n, 2, function(x){
                resid(lm(x ~ analyte.rna[,"rinvalue"])) })
tmmm.n <- t(tmmm.n)
dimnames(tmmm.n) <- dimnames(tmmm)

fn <- paste0(out_dir, "tcga.pt.lihc.TMM.log.rin.rds")
saveRDS(tmmm.n, fn)

print(dim(tmmm.n))

