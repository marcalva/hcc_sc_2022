
# Match DIEM clusters

setwd("../../")

library(diem)
library(Matrix)
library(Seurat)
library(dplyr)
library(ggplot2)


#=========================================
# Functions
#=========================================

create_dir <- function(p){
    dir.create(p, recursive=TRUE, showWarnings=FALSE)
}

#' Match clusters
#' @param m List of data frames containing marker genes and log fold changes
#' @param gene_col Column name containing the genes
#' @param feat_col Feature for each cluster
#' @param clust_col Column name giving the clusters
#'
#
match_clust <- function(m, gene_col = "gene", feat_col = "logFC", clust_col = "cluster"){
    genes.all <- lapply(m, function(i) i[,gene_col])
    genes.all <- unlist(genes.all)
    genes.all <- unique(genes.all)

    dat.l <- lapply(m, function(i){
        clusts <- unique(i[,clust_col])
        feats <- sapply(clusts, function(j) {
            dm <- i[i[,clust_col] == j, ]
            rownames(dm) <- dm[,gene_col]
            ret <- dm[genes.all, feat_col]
            names(ret) <- genes.all
            ret[is.na(ret)] <- 0
            return(ret)
        })
    })
    dat <- do.call(cbind, dat.l)
    return(dat)
}


#=========================================
# Set sample and variables
#=========================================

meta <- read.table("data/sample/samples.csv", 
                   header = TRUE, 
                   sep = ",", 
                   stringsAsFactors = FALSE, 
                   check.names = FALSE)
rownames(meta) <- meta[,"Sample"]
ids <- rownames(meta)

gtf_path <- "data/ref/gencode26/gencode.v26.annotation.fltrd.gtf"

gene_info <- read.table("data/ref/gencode26/gencode.v26.annotation.txt",
                        header = TRUE,
                        row.names = 1,
                        stringsAsFactors = FALSE,
                        sep = "\t")

#=========================================
# Output directories
#=========================================

dir_plot <- "exp/diem/plots/"; create_dir(dir_plot)

#=========================================
# Read data
#=========================================

m.l <- lapply(ids, function(s){
    fn <- paste0("exp/diem/DE/", s, "/", s, ".clust_markers.txt")
    m1 <- read.table(fn, check.names = FALSE, header = TRUE, row.names = 1, 
        stringsAsFactors = FALSE, sep = "\t")
    m1[,"cluster"] <- paste0(s, ':', m1[,"cluster"])
    return(m1)
})

dat <- match_clust(m.l)
d <- dist(t(dat))
hc <- hclust(d)

library(pheatmap)
pdf(paste0(dir_plot, "clust_match.pdf"), width = 15, height = 15)
pheatmap(cor(dat), clustering_method = "single")
dev.off()


