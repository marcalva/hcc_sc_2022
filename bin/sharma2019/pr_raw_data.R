
setwd("../../")

library(rhdf5)
library(Matrix)
library(Seurat)

# read raw data
h5fn <- "data/raw/sharma2019/HCC.h5ad"
# h5fn <- "data/raw/sharma2019/HCCF1F2.h5ad"

obsm <- h5read(file = h5fn, name = "obsm", compoundAsDataFrame=FALSE)
uns <- h5read(file = h5fn, name = "uns")
cell_dat <- h5read(file = h5fn, name = "obs")
gene_dat <- h5read(file = h5fn, name = "raw.var")
raw_x <- h5read(file = h5fn, name = "raw.X")

cdat <- as.numeric(raw_x[["data"]])
cdat <- round(exp(cdat) - 1)
counts <- sparseMatrix(i = as.numeric(raw_x[["indices"]]), 
                       p = as.numeric(raw_x[["indptr"]]),
                       x = cdat, 
                       dims = c(nrow(gene_dat), nrow(cell_dat)), 
                       dimnames = list(gene_dat[,"index"], cell_dat[,"index"]), 
                       index1 = FALSE)
# cell_dat$PIC
#    0: C (core tumor) 
#    1: I (another type of tumor?)
#    2: N (adjacent non-tumor/normal)
#    3: P (periphera tumor)
# cell_dat$PNC
#    0: C (core tumor) 
#    1: N (adjacent non-tumor/normal)
#    2: P (periphera tumor)
# cell_dat$NormalvsTumor
#    0: normal
#    1: tumor
# cell_dat$patientno 0 is healty normal control

# rename gene symbols to ensembl ID
gene_info <- read.table("data/ref/gencode26/gencode.v26.annotation.txt", 
                        header = TRUE, 
                        row.names = 1, 
                        stringsAsFactors = FALSE, 
                        sep = "\t")
gene_info[,"Name"] <- make.unique(gene_info[,"Name"])
symb2ens <- rownames(gene_info)
names(symb2ens) <- gene_info[,"Name"]
gk <- rownames(counts) %in% names(symb2ens)
counts <- counts[gk,]
rownames(counts) <- symb2ens[rownames(counts)]

# remove genes not present in nash-hcc data for integration
fn <- "data/processed/d_sct_cca_all/seur.sct.split.rds"
nash_hcc_l <- readRDS(fn)
hcc_g <- rownames(nash_hcc_l[[1]]@assays$RNA@counts)
gk <- rownames(counts) %in% hcc_g
counts <- counts[gk,]

rm(nash_hcc_l)
gc()

cell_dat[,"n_genes"] <- colSums(counts > 0)
cell_dat[,"n_counts"] <- colSums(counts)

# create Seurat object
rownames(cell_dat) <- cell_dat[,"index"]
cell_dat <- cell_dat[,-1]
seur <- CreateSeuratObject(counts = counts, 
                           project = "Sharma", 
                           meta.data = cell_dat)

# save seurat object
dir_out <- "data/processed/sharma2019/"
dir.create(dir_out, showWarnings=FALSE, recursive=TRUE)

outfn <- file.path(dir_out, "sharma.seur.rds")
saveRDS(seur, outfn)

