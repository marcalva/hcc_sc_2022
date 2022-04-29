
# UMAP plot of cell type markers

setwd("../../")

library(Seurat)
library(ggplot2)
library(RColorBrewer)
source("scripts/plot.R")

#=========================================
# Functions
#=========================================

#' Plot UMAP
plot_umap_gene <- function(i, ids, genes, gene_info, 
                      seur, dplot){
    id <- ids[i]
    gene <- genes[i]
    id <- rownames(gene_info)[gene_info$Name == gene]
    p <- plot_umap_gene_2c(seur, gene = id, assay = "RNA", low = "grey75", size=.25) + 
    ggtitle(gene) + theme(plot.title = element_text(hjust = 0.5))
    print(p)
    ggplot2::ggsave(paste0(dplot, "umap.", gene, ".jpeg"), 
           width = 4, height = 4, units = "in", device = "jpeg", 
           dpi = 300)
    dev.off()
}

#' Plot boxplot
plot_box_gene <- function(i, ids, genes, gene_info, 
                      seur, dplot, ct_col = "CellType"){
    id <- ids[i]
    gene <- genes[i]
    p <- boxplot_gene(seur, gene = id, assay = "RNA", ct_col = ct_col) +
    ggtitle(gene) + ylab("logUMI") + 
    theme(plot.title = element_text(hjust = 0.5), 
          axis.text.x = element_text(hjust = 1, angle = 45))
    print(p)
    ggplot2::ggsave(paste0(dplot, "box.", gene, ".jpeg"), 
           width = 3, height = 3, units = "in", device = "jpeg", 
           dpi = 300)
    dev.off()
}

#=========================================
# Set variables
#=========================================

gene_info <- read.table("data/ref/gencode26/gencode.v26.annotation.txt", 
                        header = TRUE, 
                        row.names = 1, 
                        stringsAsFactors = FALSE, 
                        sep = "\t")

seur <- readRDS("data/processed/d_sct_cca_all/hep2/seur.cca.rds")

Idents(seur) <- "integrated_snn_res.1"
ct_col <- "integrated_snn_res.1"

seur@active.assay <- "RNA"
seur <- NormalizeData(seur, normalization.method = "LogNormalize", scale.factor = 1000)

dbase <- "exp/d_sct_cca_all/hep2/plots_genes/"
dir.create(dbase, showWarnings=FALSE, recursive=TRUE)

#=========================================
# Endothelial cells
#=========================================
# VWF, CD36, LDB2, AKAP12

dplot <- paste0(dbase, "HSEC/")
dir.create(dplot, showWarnings=FALSE, recursive=TRUE)

genes <- read.table("data/ref/marker_genes/hsec_markers.1.txt", 
                    stringsAsFactors = FALSE)
genes <- genes[,1]
genes <- unique(c(genes))
ids <- rownames(gene_info)[gene_info$Name %in% genes]
genes <- gene_info[ids, "Name"]

k <- which(ids %in% rownames(seur@assays$RNA@data))
ids <- ids[k]
genes <- genes[k]

for (i in 1:length(genes)){
    message(genes[i])
    plot_umap_gene(i, ids, genes, gene_info, 
                   seur, dplot)
    plot_box_gene(i, ids, genes, gene_info, 
                  seur, dplot, ct_col = ct_col)
}

#=========================================
# Hepatocytes
#=========================================

dplot <- paste0(dbase, "Hepatocytes/")
dir.create(dplot, showWarnings=FALSE, recursive=TRUE)

genes <- read.table("data/ref/marker_genes/hep_markers.1.txt", 
                    stringsAsFactors = FALSE)
genes <- genes[,1]
genes <- unique(c(genes))
ids <- rownames(gene_info)[gene_info$Name %in% genes]
genes <- gene_info[ids, "Name"]

k <- which(ids %in% rownames(seur@assays$RNA@data))
ids <- ids[k]
genes <- genes[k]

for (i in 1:length(genes)){
    message(genes[i])
    plot_umap_gene(i, ids, genes, gene_info, 
                   seur, dplot)
    plot_box_gene(i, ids, genes, gene_info, 
                  seur, dplot, ct_col = ct_col)
}


#=========================================
# T
#=========================================

dplot <- paste0(dbase, "T/")
dir.create(dplot, showWarnings=FALSE, recursive=TRUE)

genes <- read.table("data/ref/marker_genes/t_markers.1.txt", 
                    stringsAsFactors = FALSE)
genes <- genes[,1]
genes <- unique(c(genes))
ids <- rownames(gene_info)[gene_info$Name %in% genes]
genes <- gene_info[ids, "Name"]

k <- which(ids %in% rownames(seur@assays$RNA@data))
ids <- ids[k]
genes <- genes[k]

for (i in 1:length(genes)){
    message(genes[i])
    plot_umap_gene(i, ids, genes, gene_info, 
                   seur, dplot)
    plot_box_gene(i, ids, genes, gene_info, 
                  seur, dplot, ct_col = ct_col)
}


#=========================================
# Kupffer
#=========================================

dplot <- paste0(dbase, "Kup/")
dir.create(dplot, showWarnings=FALSE, recursive=TRUE)

genes <- read.table("data/ref/marker_genes/kupf_markers.1.txt", 
                    stringsAsFactors = FALSE)
genes <- genes[,1]
genes <- unique(c(genes))
ids <- rownames(gene_info)[gene_info$Name %in% genes]
genes <- gene_info[ids, "Name"]

k <- which(ids %in% rownames(seur@assays$RNA@data))
ids <- ids[k]
genes <- genes[k]

for (i in 1:length(genes)){
    message(genes[i])
    plot_umap_gene(i, ids, genes, gene_info, 
                   seur, dplot)
    plot_box_gene(i, ids, genes, gene_info, 
                  seur, dplot, ct_col = ct_col)
}


#=========================================
# Chol
#=========================================

dplot <- paste0(dbase, "Chol/")
dir.create(dplot, showWarnings=FALSE, recursive=TRUE)

genes <- read.table("data/ref/marker_genes/chol_markers.1.txt", 
                    stringsAsFactors = FALSE)
genes <- genes[,1]
genes <- unique(c(genes))
ids <- rownames(gene_info)[gene_info$Name %in% genes]
genes <- gene_info[ids, "Name"]

k <- which(ids %in% rownames(seur@assays$RNA@data))
ids <- ids[k]
genes <- genes[k]

for (i in 1:length(genes)){
    message(genes[i])
    plot_umap_gene(i, ids, genes, gene_info, 
                   seur, dplot)
    plot_box_gene(i, ids, genes, gene_info, 
                  seur, dplot, ct_col = ct_col)
}

#=========================================
# Stellate
#=========================================

dplot <- paste0(dbase, "Stellate/")
dir.create(dplot, showWarnings=FALSE, recursive=TRUE)

genes <- read.table("data/ref/marker_genes/stel_markers.1.txt", 
                    stringsAsFactors = FALSE)
genes <- genes[,1]
genes <- unique(c(genes))
ids <- rownames(gene_info)[gene_info$Name %in% genes]
genes <- gene_info[ids, "Name"]

k <- which(ids %in% rownames(seur@assays$RNA@data))
ids <- ids[k]
genes <- genes[k]

for (i in 1:length(genes)){
    message(genes[i])
    plot_umap_gene(i, ids, genes, gene_info, 
                   seur, dplot)
    plot_box_gene(i, ids, genes, gene_info, 
                  seur, dplot, ct_col = ct_col)
}



#=========================================
# Neural
#=========================================

dplot <- paste0(dbase, "Neur/")
dir.create(dplot, showWarnings=FALSE, recursive=TRUE)

genes <- read.table("data/ref/marker_genes/neur_markers.1.txt", 
                    stringsAsFactors = FALSE)
genes <- genes[,1]
genes <- unique(c(genes))
ids <- rownames(gene_info)[gene_info$Name %in% genes]
genes <- gene_info[ids, "Name"]

k <- which(ids %in% rownames(seur@assays$RNA@data))
ids <- ids[k]
genes <- genes[k]

for (i in 1:length(genes)){
    message(genes[i])
    plot_umap_gene(i, ids, genes, gene_info, 
                   seur, dplot)
    plot_box_gene(i, ids, genes, gene_info, 
                  seur, dplot, ct_col = ct_col)
}


