
# Run limma genome-wide DE on microarray data from LCI cohort

setwd("../../")

library(limma)

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

runlimma <- function(y, ds, coef = ncol(ds)){
    require(limma)

    fit <- lmFit(y, ds)
    fit <- eBayes(fit, trend=TRUE)
    diffTable <- topTable(fit, coef = coef)

}

# source("scripts/rnaseq_scripts/R/de.R")

dir_exp <- "exp/GSE14520/de/"
dir.create(dir_exp, showWarnings = FALSE, recursive=TRUE)

dir_data <- "data/processed/GSE14520/"
dir.create(dir_data, recursive=TRUE, showWarning=FALSE)

# Gene data
gencode <- read.table("data/ref/gencode26/gencode.v26.annotation.txt", 
                        header = TRUE, 
                        row.names = 1, 
                        stringsAsFactors = FALSE, 
                        sep = "\t")
gencode[,"ens"] <- rownames(gencode)
gencode[,"Name"] <- make.unique(gencode[,"Name"])
rownames(gencode) <- gencode[,"Name"]

fn <- paste0(dir_data, "expr.RMA_log2.gmean.rds")
ex <- readRDS(fn)

gk <- rownames(ex)[rownames(ex) %in% gencode[,"Name"]]

fn <- paste0(dir_data, "sdata.txt")
sdat <- read.table(fn, row.names=1, header=TRUE, sep="\t", stringsAsFactors=FALSE)

fn <- paste0(dir_data, "geo_pheno.txt")
geo_pheno <- read.table(fn, row.names=1, header=TRUE, sep="\t", stringsAsFactors=FALSE)

#========================================================
# Get tumor and non-tumor sample data separately
#========================================================

sdat.tum <- sdat[sdat[,"Tissue.Type"] == "Tumor",]
sdat.nt <- sdat[sdat[,"Tissue.Type"] == "Non-tumor",]

table(duplicated(sdat.tum[,"ID"]))
table(duplicated(sdat.nt[,"ID"]))

rownames(sdat.tum) <- sdat.tum[,"ID"]
rownames(sdat.nt) <- sdat.nt[,"ID"]

lcs.tum <- rownames(sdat.tum)
lcs.nt <- rownames(sdat.nt)

ids.pair <- intersect(rownames(sdat.tum), rownames(sdat.nt))
sdat.tum.p <- sdat.tum[ids.pair,]
sdat.nt.p <- sdat.nt[ids.pair,]

lcs.tum.p <- sdat.tum.p[,"LCS.ID"]
lcs.nt.p <- sdat.nt.p[,"LCS.ID"]

rownames(sdat.tum) <- sdat.tum[,"LCS.ID"]
rownames(sdat.nt) <- sdat.nt[,"LCS.ID"]

#=====================================================
# paired DEP
#=====================================================

# get paired sample and case IDs
ktum <- sdat[,"Tissue.Type"] == "Tumor"
cases.tum <- sdat[ktum, "ID"]
knontum <- sdat[,"Tissue.Type"] == "Non-tumor"
cases.nontum <- sdat[knontum, "ID"]

case.pair <- intersect(cases.tum, cases.nontum)

k.pair <- sdat[,"ID"] %in% case.pair
sample.pair <- sdat[k.pair, "LCS.ID"]

# subset expression and pheno
sdat.k <- sdat[sample.pair,]
ex.k <- ex[,sample.pair]

# design matrix
sdat.k[,"Tissue.Type"] <- factor(sdat.k[,"Tissue.Type"], 
                                    levels = c("Non-tumor", 
                                               "Tumor"))
design.k <- model.matrix( ~ sdat.k[,"ID"] + sdat.k[,"Tissue.Type"])

# limma
fit <- lmFit(ex.k, design.k)
fit <- eBayes(fit, trend=TRUE)
diff_table <- topTable(fit, coef = ncol(design.k), number = nrow(ex.k))
diff_table[,"gene_lci"] <- rownames(diff_table)
diff_table <- cbind(gencode[rownames(diff_table),], diff_table)
rownames(diff_table) <- diff_table[,"gene_lci"]

# export
out_dir <- "exp/GSE14520/de/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# save LRT
out_fn <- paste0(dir_exp, "fit.rds")
saveRDS(fit, out_fn)

# save table
out_fn <- "de.table.tsv.gz"
save2file(diff_table, dirout = dir_exp, fn = out_fn, gz=TRUE)

