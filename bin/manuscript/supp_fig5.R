
# Figure S5
# plot prolif clustering

setwd("../../")

library(Seurat)
library(reshape2)
library(ggplot2)
library(gtable)
library(grid)
library(gridExtra)
library(ggrepel)
source("scripts/color_pal.R")
source("scripts/gtable_stack.R")

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

mar_pct <- 0.025
tagl <- 3

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

fn <- "data/processed/sharma_aiz/liver.prolif.int.rds"
prolif <- readRDS(fn)

# data frame of meta data
umap <- prolif@reductions$umap@cell.embeddings
colnames(umap) <- c("UMAP1", "UMAP2")
datf <- prolif@meta.data
datf <- cbind(datf, umap[rownames(datf),])

# colors
reds <- read.table("data/ref/colors/red_colrs.csv", header=FALSE)
reds <- reds[,1]

# Set directories
dir_plt <- "exp/manuscript/"
dir.create(dir_plt, showWarnings=FALSE, recursive=TRUE)

# percent
tbm <- table(prolif$cell_type_main)
print(tbm / sum(tbm))

#=========================================
# UMAP plots
#=========================================

th_umap <- theme_classic() +
theme(text = element_text(size = 8), 
      axis.text = element_blank(), 
      axis.ticks = element_blank(), 
      plot.title = element_text(hjust = 0.5), 
      legend.key.size = unit(6, "points"))
                 
ids <- c("integrated_snn_res.0.2", "cell_type_main")
leg_names <- c("Cluster", "Main\ncell-type")
names(leg_names) <- ids

cpal <- list()
cpal[["integrated_snn_res.0.2"]] <- 
    hcl_pal(length(unique(datf[,"integrated_snn_res.0.2"])),
            chr = c(100, 120), lum = c(60, 70), offset = 0,
            rand = TRUE, seedn = 5)
cpal[["cell_type_main"]] <- 
    hcl_pal(length(unique(datf[,"cell_type_main"])), 
            chr = c(100, 120), lum = c(60, 70), offset = 0, 
            rand = TRUE, seedn = 5)

tags_umaps <- c("a", "c")
names(tags_umaps) <- ids

gr_umaps <- list()
for ( id in ids ){
    p <- ggplot(datf, aes_string(x = "UMAP1", y = "UMAP2", color = id)) + 
    geom_point(shape=16, size=1) +
    scale_color_manual(values=cpal[[id]], name = leg_names[id]) +
    guides(color = guide_legend(override.aes = list(size=1))) + 
    th_umap

    gr_id <- ggplotGrob(p)

    gr_id <- gtable_add_tag(gr_id, tags_umaps[id], fs=16, just=c(0,0), l=tagl)
    gr_id <- pad_plot(gr_id, t=mar_pct*4, r=mar_pct, b=mar_pct, l=mar_pct)

    gr_umaps[[id]] <- gr_id
}

#=========================================
# Overlap of clusters and singler
#=========================================

th_hmap <- theme_classic() + 
theme(text = element_text(size = 8), 
      plot.title = element_text(hjust = 0.5), 
      legend.key.size = unit(6, "points"))
                 

tb <- table(datf[,c("integrated_snn_res.0.2", "cell_type_main")])
cl_ct_frq <- tb
cl_ct_prp <- sweep(cl_ct_frq, 1, rowSums(cl_ct_frq), '/')
cl_ct_prpm <- melt(cl_ct_prp)
colnames(cl_ct_prpm) <- c("x", "y", "z")
for (i in 1:2){
    cl_ct_prpm[,i] <- factor(cl_ct_prpm[,i])
}

p <- ggplot(cl_ct_prpm, aes(x = x, y = y, fill = z)) + 
geom_tile() + 
scale_fill_gradientn(colors = reds, limits = c(0,1), name = "Prop") + 
labs(x="Cluster", y="Main cell type") + 
th_hmap

gr_hm <- ggplotGrob(p)

gr_hm <- gtable_add_tag(gr_hm, "b", fs = 16, just = c(0,0), l=tagl)
gr_hm <- pad_plot(gr_hm, t=mar_pct*4, r=mar_pct, b=mar_pct, l=mar_pct)

#=========================================
# Plot marker genes
#=========================================

cdatf <- datf
rc <- prolif@assays$RNA@data

genes <- c("ASGR1", "LYVE1", "CD3E", "CD68", "MUC6", "ACTA2")
names(genes) <- c("Hepatocyte", "Endothelial", "T", 
                  "Macrophage", "Cholangiocyte", "Stellate")

ids <- symb2ens[genes]
k <- which(ids %in% rownames(rc))
genes <- genes[k]
ids <- ids[k]

tags_gene <- c("d", "e", "f", "g", "h", "i")

gr_genes <- list()
for (i in 1:length(genes)){
    gene <- genes[i]
    id <- ids[i]
    cdatf[,"gene"] <- rc[id,rownames(cdatf)]
    o <- order(cdatf[,"gene"], decreasing=FALSE)

    p <- ggplot(cdatf[o,], aes(x = UMAP1, y = UMAP2, color = gene)) + 
    geom_point(shape=16, size = 1) + 
    ggtitle(gene) + 
    scale_color_gradientn(colors = reds, name = "logUMI") + 
    th_umap

    gr_id <- ggplotGrob(p)

    gr_id <- gtable_add_tag(gr_id, tags_gene[i], fs = 16, just = c(0,1), l=tagl)
    gr_id <- pad_plot(gr_id, t=mar_pct, r=mar_pct, b=mar_pct, l=mar_pct)

    gr_genes[[i]] <- gr_id
}

#=========================================
# put together a gtable for plotting everything
#=========================================

w <- 3
h <- 2.5
punit <- "inches"

ws <- c(w,w,w)
hs <- c(h,h,h)

gt <- gtable(widths = unit(ws, punit), 
             heights = unit(hs, punit))

gt <- gtable_add_grob(gt, gr_umaps[[1]], t=1,b=1,l=1,r=1)
gt <- gtable_add_grob(gt, gr_hm, t=1,b=1,l=2,r=2)
gt <- gtable_add_grob(gt, gr_umaps[[2]], t=1,b=1,l=3,r=3)

l_v <- c(1,2,3,1,2,3)
t_v <- c(2,2,2,3,3,3)
for (i in 1:length(gr_genes)){
    gt <- gtable_add_grob(gt, gr_genes[[i]], t=t_v[i],b=t_v[i],l=l_v[i],l_v[i])
}

pdf(paste0(dir_plt, "FigureS5.pdf"), width = sum(ws), height = sum(hs))
grid.draw(gt)
dev.off()

