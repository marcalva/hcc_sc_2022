
library(GEOquery)

# load series and platform data from GEO
gset <- getGEO("GSE14520", GSEMatrix =TRUE, AnnotGPL=TRUE)[[1]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- "undefined"
sml <- strsplit(gsms, split="")[[1]]

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
          (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
  exprs(gset) <- log2(ex) }

# supp
fn <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE14nnn/GSE14520/suppl/GSE14520%5FExtra%5FSupplement%2Etxt%2Egz"
sup <- gzcon(url(fn))
sup <- readLines(sup)
sdat <- read.table(textConnection(sup), sep = '\t', header=TRUE, stringsAsFactors = FALSE)

# include tumor/non-tumor to those with survival data
ind_ids <- sdat[ !is.na(sdat[,"Survival.status"]) , "ID" ]
sdat <- sdat[ sdat[,"ID"] %in% ind_ids, ]

# overlap survival data with expression data
k <- colnames(gset) %in% sdat[,"Affy_GSM"]
gset <- gset[,k]

k <- sdat[,"Affy_GSM"] %in% colnames(gset)
sdat <- sdat[k,]

# 221 tumor and 209 non-tumor.
sdat.tum <- sdat[sdat[,"Tissue.Type"] == "Tumor",]
sdat.nt <- sdat[sdat[,"Tissue.Type"] == "Non-Tumor",]


# map GSM ID to sample ID in gset
gsm2lcs <- sdat[,"LCS.ID"]
names(gsm2lcs) <- sdat[,"Affy_GSM"]

colnames(gset) <- gsm2lcs[colnames(gset)]

# match sdat to gset
rownames(sdat) <- sdat[,"LCS.ID"]
sdat <- sdat[colnames(gset),]

# get median expression value for multiple genes
ex <- exprs(gset)
table(rownames(ex) == rownames(features)) # check expr & feature names match
features <- fData(gset)
gene.name <- features[,"Gene.symbol"]
k <- gene.name != "" # exclude probes with no matching symbol
features.k <- features[k,]
gene.name <- features.k[,"Gene.symbol"]
ex.k <- ex[k,]

g.uniq <- unique(gene.name)

ex.k.mean.l <- lapply(g.uniq, function(g){
                      rk <- features.k[,"Gene.symbol"] == g
                      colMeans(ex.k[rk,,drop=FALSE]) })
names(ex.k.mean.l) <- g.uniq
ex.k.mean <- do.call(rbind, ex.k.mean.l)

# The expression data identified by pheno[,"submission_date"]
# data are already normalized by RMA with log2 transformed values

pheno <- pData(gset)

# change Non-Tumor to Non-tumor
knt <- sdat[,"Tissue.Type"] == "Non-Tumor"
sdat[knt, "Tissue.Type"] <- "Non-tumor"

# get tumor stage low v high
sdat[,"stage_low_high"] <- "Low"
khigh <- sdat[,"TNM.staging"] %in% c("III", "IIIA", "IIIB", "IIIC")
sdat[khigh, "stage_low_high"] <- "High"
kna <- sdat[,"TNM.staging"] == "."
sdat[kna, "stage_low_high"] <- NA
knt <- sdat[,"Tissue.Type"] == "Non-tumor"
sdat[knt, "stage_low_high"] <- NA


# directory out
dir_out <- "data/processed/GSE14520/"
dir.create(dir_out, recursive=TRUE, showWarning=FALSE)

fn <- paste0(dir_out, "expr.RMA_log2.gmean.txt")
write.table(ex.k.mean, fn, row.names=TRUE, col.names=NA, quote=FALSE, sep="\t")
fn <- paste0(dir_out, "expr.RMA_log2.gmean.rds")
saveRDS(ex.k.mean, fn)

fn <- paste0(dir_out, "sdata.txt")
write.table(sdat, fn, row.names=TRUE, col.names=NA, quote=FALSE, sep="\t")

fn <- paste0(dir_out, "geo_pheno.txt")
write.table(pheno, fn, row.names=TRUE, col.names=NA, quote=FALSE, sep="\t")


