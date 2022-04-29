
setwd("../../")

library(Seurat)
library(ggplot2)
library(scales)
library(gridExtra)
source("scripts/standard_seurat.R")
source("scripts/plot.R")

#=========================================
#=========================================

n_threads <- 8

# Gene data
gene_info <- read.table("data/ref/gencode26/gencode.v26.annotation.txt", 
                        header = TRUE, 
                        row.names = 1, 
                        stringsAsFactors = FALSE, 
                        sep = "\t")
gene_info[,"Name"] <- make.unique(gene_info[,"Name"])
symb2ens <- rownames(gene_info)
names(symb2ens) <- gene_info[,"Name"]

fn <- "data/processed/sharma2019/sharma.seur.rds"
sharma <- readRDS(fn)

#=========================================
# SCT no integration
#=========================================

set.seed(1)

sharma <- SCTransform(sharma, variable.features.n = 3000, verbose=TRUE)

sharma <- RunPCA(sharma, npcs = 50, verbose = FALSE)
ElbowPlot(sharma, ndims = 50)
sharma <- RunUMAP(sharma, dims = 1:30, verbose = FALSE)
sharma <- FindNeighbors(sharma, dims = 1:30, verbose = FALSE)
sharma <- FindClusters(sharma, verbose = FALSE)

# cell cycle scoring
DefaultAssay(sharma) <- "RNA"
sharma <- NormalizeData(sharma, scale.factor=1000)
s.genes <- cc.genes$s.genes
s.genes <- symb2ens[s.genes]
s.genes <- s.genes[s.genes %in% rownames(sharma)]
g2m.genes <- cc.genes$g2m.genes
g2m.genes <- symb2ens[g2m.genes]
g2m.genes <- g2m.genes[s.genes %in% rownames(sharma)]
sharma <- CellCycleScoring(sharma, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
DefaultAssay(sharma) <- "SCT"

fn <- "data/processed/sharma2019/sharma.seur.rds"
saveRDS(sharma, fn)

#=========================================
# Plots
#=========================================

dir_plt <- "exp/sharma2019/sct_clust/"
dir.create(dir_plt, showWarnings=FALSE, recursive=TRUE)

# G2M scoring
p <- plot_umap_c(sharma, colname = "G2M.Score", 
                 legend_title = "G2M score", size = 0.5)
ggsave(file.path(dir_plt, "umap.G2M_score.jpeg"), 
       width = 6, height = 5, dpi  = 300)

# percent mitochondria
p <- plot_umap_c(sharma, colname = "percent_mito", 
                 legend_title = "percent_mito", size = 0.5)
ggsave(file.path(dir_plt, "umap.percent_mito.jpeg"), 
       width = 6, height = 5, dpi  = 300)

# sample labels
p <- plot_umap_d(sharma, colname = "patientno", 
                 legend_title = "patient", size = 0.5)
ggsave(file.path(dir_plt, "umap.patientno.jpeg"), 
       width = 6, height = 5, dpi  = 300)

# tumor non-tumor
p <- plot_umap_d(sharma, colname = "PNC", 
                 legend_title = "PNC", size = 0.5)
ggsave(file.path(dir_plt, "umap.PNC.jpeg"), 
       width = 6, height = 5, dpi  = 300)

# louvain clusters
p <- plot_umap_labels(sharma, colname = "louvain", 
                 legend_title = "louvain", size = 0.5) + 
theme(legend.position="none")
ggsave(file.path(dir_plt, "umap.louvain.jpeg"), 
       width = 6, height = 6, dpi  = 300)

