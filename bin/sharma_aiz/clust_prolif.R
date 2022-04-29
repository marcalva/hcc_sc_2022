
# cluster proliferating cells, assign to cell-types with SingleR

setwd("../../")

library(Seurat)
library(ggplot2)
library(grid)
library(gridExtra)
library(ggrepel)
suppressPackageStartupMessages(library(SingleR))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(scran))
source("scripts/color_pal.R")

#=========================================
# Functions
#=========================================

create_dir <- function(p){
    dir.create(p, showWarnings=FALSE, recursive=TRUE)
}

#' @param datf data frame to get median points from
#' @param x character vector of column names to calculate median for
#' @param groupby the column name to group the rows of the data frame by
get_med_points <- function(datf, x, groupby){
    groups <- sort(unique(datf[,groupby]))
    gs.l <- lapply(groups, function(gr){
                 k <- datf[,groupby] == gr
                 datf.s <- datf[k, x, drop=FALSE]
                 r <- apply(datf.s, 2, median)
                 return(r) })
    names(gs.l) <- groups
    gs <- as.data.frame(do.call(rbind, gs.l))
    colnames(gs) <- x
    rownames(gs) <- groups
    return(gs)
}

#=========================================
#=========================================

# Gene data
gene_info <- read.table("data/ref/gencode26/gencode.v26.annotation.txt", 
                        header = TRUE, 
                        row.names = 1, 
                        stringsAsFactors = FALSE, 
                        sep = "\t")
gene_info[,"Name"] <- make.unique(gene_info[,"Name"])
symb2ens <- rownames(gene_info)
names(symb2ens) <- gene_info[,"Name"]

fn <- "data/processed/sharma_aiz/liver.int_rand.rds"
integrated <- readRDS(fn)

# read non-prolif references
fn <- "exp/sharma_aiz/clust_prolif/cell_type_main.nonprol_ref.rds"
main_ref <- readRDS(fn)

# Set directories
dir_out <- "data/processed/sharma_aiz/"
dir.create(dir_out, showWarnings=FALSE, recursive=TRUE)

dir_exp <- "exp/sharma_aiz/clust_prolif/";
dir.create(dir_exp, showWarnings=FALSE, recursive=TRUE)

#=========================================
# classify using SingleR reference
#=========================================

md <- integrated@meta.data

# cell_type_main
cl_id <- "cell_type_main"
md[,cl_id] <- as.character(md[,cl_id])
k_prol <- rownames(md)[md[,cl_id] == "Prol"]
ex <- integrated@assays$RNA@data
ex_prol <- ex[,k_prol]

clss_main <- classifySingleR(test = ex_prol, trained = main_ref)

#=========================================
# get UMAP of prolif cells
#=========================================

set.seed(1, kind = 'Mersenne-Twister')

md <- integrated@meta.data
k <- rownames(md)[md[,"cell_type_main"] == "Prol"]
prolif <- subset(integrated, cells = k)

# have to split by source since too few samples for integration
s_l <- SplitObject(prolif, split.by = "source")

for (i in 1:length(s_l)){
    message(i)
    DefaultAssay(s_l[[i]]) <- "RNA"
    s_l[[i]] <- SCTransform(s_l[[i]], variable.features.n = 3000, 
                            verbose=FALSE)
}

features <- SelectIntegrationFeatures(object.list = s_l, nfeatures = 3000)
s_l <- PrepSCTIntegration(object.list = s_l, anchor.features = features)

# minimum number of samples is 92, so k.filter and k.weight need to 
#  be below this.
k_n <- 75
anchors <- FindIntegrationAnchors(object.list = s_l, 
                                  dims = 1:30, 
                                  normalization.method = "SCT", 
                                  reduction = "cca", 
                                  anchor.features = features, 
                                  k.filter = k_n)
prolif <- IntegrateData(anchorset = anchors, 
                        normalization.method = "SCT", 
                        dims = 1:30, k.weight=k_n)

# rm(s_l)
# g <- gc()

# 20 PCs is around elbow point
prolif <- ScaleData(prolif, verbose = FALSE)
prolif <- RunPCA(prolif, npcs = 50, verbose = FALSE)
prolif <- RunUMAP(prolif, reduction = "pca", dims = 1:30)
prolif <- FindNeighbors(prolif, reduction = "pca", dims = 1:30)
prolif <- FindClusters(prolif, resolution = c(0.2, 0.5))

# add singler classifications
prolif$cell_type_main <- clss_main[rownames(prolif@meta.data), "pruned.labels"]

outfn <- file.path(dir_out, "liver.prolif.int.rds")
saveRDS(prolif, outfn)

