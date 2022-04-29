
setwd("../../")

library(Seurat)
library(ggplot2)
library(scales)
library(gridExtra)
source("scripts/standard_seurat.R")

#=========================================
# Functions
#=========================================

create_dir <- function(p){
    dir.create(p, showWarnings=FALSE, recursive=TRUE)
}

get_lm_fast <- function(x, psc = 1){
    x <- Matrix::t(x)

    cmeans <- sapply(1:ncol(x), function(i1){
                     i2 = i1+1
                     ip1 <- x@p[i1]
                     ip2 <- x@p[i2]
                     nnz <- ip2 - ip1
                     nz <- nrow(x) - nnz
                     ent <- exp(x@x[(ip1+1):ip2]) - 1 + psc
                     csum <- sum(c(ent, psc * nz))
                     cmean <- log2( csum / (nz + nnz) )
                     return(cmean) })

    names(cmeans) <- colnames(x)
    return(cmeans)
}

# get logFC for all genes
get_lfc <- function(seur, ct_id, psc = 1){
    md <- seur@meta.data
    if (!ct_id %in% colnames(md)){
        stop(ct_id, " not found in seurat object")
    }
    cts <- sort(as.character(unique(md[,ct_id])))
    lfc_l <- list()
    for (ct in cts){
        message(ct)
        nct <- setdiff(cts, ct)

        k1 <- rownames(md)[md[,ct_id] %in% ct]
        k2 <- rownames(md)[md[,ct_id] %in% nct]

        mat1 <- seur@assays$RNA@data[,k1,drop=FALSE]
        mat2 <- seur@assays$RNA@data[,k2,drop=FALSE]
        kg <- rowSums(mat1) > 0 & rowSums(mat2) > 0
        mat1 <- mat1[kg,]
        mat2 <- mat2[kg,]

        m1 <- get_lm_fast(mat1)
        m2 <- get_lm_fast(mat2)

        lfc <- m1 - m2
        lfc_l[[ct]] <- lfc
    }

    g_all <- lapply(lfc_l, names)
    g_u <- Reduce(union, g_all)
    lfc_l <- lapply(lfc_l, function(x) x[g_u])
    lfc_mat <- do.call(cbind, lfc_l)

    return(lfc_mat)
}

#=========================================
# Set variables
#=========================================

n_threads <- 8

# Gene data
gene_info <- read.table("data/ref/gencode26/gencode.v26.annotation.txt", 
                        header = TRUE, 
                        row.names = 1, 
                        stringsAsFactors = FALSE, 
                        sep = "\t")
gene_info[,"Name"] <- make.unique(gene_info[,"Name"])

# Set directories
dir_out <- "data/processed/aizarani2019/"; create_dir(dir_out)
# dir_plot <- "exp/d_sct_cca_all/cca_init/"; create_dir(dir_plot)

dir_exp <- "exp/aizarani2019/sct_clust/"; create_dir(dir_exp)

#=========================================
# Cluster
#=========================================

# Read in raw data  
counts <- readRDS("data/raw/aizarani2019/GSE124395_Normalhumanliverdata.RData")
counts <- round(counts)
counts <- as(as.matrix(counts), "sparseMatrix")
gc()

fn <- "data/raw/aizarani2019/GSE124395_clusterpartition.txt"
aiz_clust <- read.table(fn, row.names=1, header=TRUE, check.names=FALSE)

sk <- intersect(colnames(counts), rownames(aiz_clust))
counts <- counts[,sk,drop=FALSE]

gk <- rownames(counts) %in% gene_info[,"Name"]
counts <- counts[gk,]

# rename counts
symb2ens <- rownames(gene_info)
names(symb2ens) <- gene_info[,"Name"]
rownames(counts) <- symb2ens[rownames(counts) ]

s <- gsub("_*_*$", "", colnames(counts)) 
s <- gsub("_[0-9]*_[0-9]*$", "", colnames(counts))

md <- data.frame("aiz_clust" = aiz_clust[colnames(counts),], 
                 "aiz_sample" = s)
rownames(md) <- colnames(counts)

set.seed(1, kind = "Mersenne-Twister")

seur <- Seurat::CreateSeuratObject(counts = counts, "Aizarani", 
                                   meta.data = md) 

seur <- NormalizeData(seur, normalization.method = "LogNormalize", 
                      scale.factor = 1e3)

# add pct ribosomal
rb_s <- readLines("data/ref/gene_lists/rps_genes.txt")
rb_l <- readLines("data/ref/gene_lists/rpl_genes.txt")
rb <- c(rb_s, rb_l)
rb <- rb[rb %in% gene_info[,"Name"]]
rb <- symb2ens[rb]
rb <- rb[rb %in% rownames(seur)]
seur <- PercentageFeatureSet(seur, features=rb, col.name="pct.rb")

# cell cycle scoring
s.genes <- cc.genes$s.genes
s.genes <- symb2ens[s.genes]
s.genes <- s.genes[s.genes %in% rownames(seur)]
g2m.genes <- cc.genes$g2m.genes
g2m.genes <- symb2ens[g2m.genes]
g2m.genes <- g2m.genes[s.genes %in% rownames(seur)]
seur <- CellCycleScoring(seur, s.features = s.genes, g2m.features = g2m.genes, 
                         set.ident = FALSE, assay = "RNA")

seur <- FindVariableFeatures(seur, selection.method="vst", nfeatures=2000)
seur <- ScaleData(seur)

seur <- SCTransform(seur, variable.features.n = 3000, 
                    verbose=TRUE)

seur <- RunPCA(seur, verbose=FALSE, npcs=50)

ElbowPlot(seur, ndims=50)

n_dim <- 30
seur <- FindNeighbors(seur, dims = 1:n_dim, verbose = TRUE)
seur <- FindClusters(seur, resolution = c(0.2, 0.5, 0.8, 1), verbose = TRUE)
seur <- RunUMAP(seur, dims = 1:n_dim, reduction = "pca", verbose = TRUE)

saveRDS(seur, paste0(dir_out, "seur.rds"))

#=========================================
# Get markers (res 0.5)
#=========================================

dir_mrk <- paste0(dir_exp, "markers/")
create_dir(dir_mrk)

Idents(seur) <- "SCT_snn_res.0.5"
DefaultAssay(seur) <- "RNA"

m <- FindAllMarkers(seur, logfc.threshold = 0.1, test.use = "wilcox",
                    only.pos = TRUE)

m <- cbind(gene_info[m[,"gene"],], m)

write.table(m, paste0(dir_mrk, "mrk.SCT_snn_res.0.5.txt"), 
            row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')

# get log fold change
lfc_mat <- get_lfc(seur, "SCT_snn_res.0.5", psc=1)
lfc_mat <- cbind(gene_info[rownames(lfc_mat),], lfc_mat)

ofn <- paste0(dir_mrk, "log_fc.SCT_snn_res.0.5.txt")
write.table(lfc_mat, 
            ofn, 
            row.names = TRUE, col.names = NA, 
            sep = "\t", quote = FALSE)



