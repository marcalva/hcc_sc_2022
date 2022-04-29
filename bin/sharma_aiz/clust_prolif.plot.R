
# plot cluster proliferating cells

setwd("../../")

library(Seurat)
library(reshape2)
library(ggplot2)
library(grid)
library(gridExtra)
library(ggrepel)
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

fn <- "data/processed/sharma_aiz/liver.prolif.int.rds"
prolif <- readRDS(fn)

# Set directories
dir_out <- "data/processed/sharma_aiz/"
dir.create(dir_out, showWarnings=FALSE, recursive=TRUE)

dir_exp <- "exp/sharma_aiz/clust_prolif/";
dir.create(dir_exp, showWarnings=FALSE, recursive=TRUE)

#=========================================
# plot themes and colors
#=========================================

theme_trnsp <- theme(plot.background = element_rect(fill="transparent", color=NA),
                     panel.background = element_rect(fill="transparent", color=NA),
                     legend.background = element_rect(fill="transparent", color=NA))
theme_leg <- theme(legend.key.height = unit(2, "strheight", "0"),
                   legend.key.width = unit(1, "strwidth", "0"))
theme_txt <- theme(text = element_text(size = 8),
                   axis.text=element_blank(), 
                   axis.ticks=element_blank(),
                   plot.title = element_text(hjust = 0.5))
theme_s <- theme_classic() + 
    theme_leg + 
    theme_txt

theme_p <- theme_classic() +
    theme_txt + 
    theme(legend.key.height = unit(1, "strheight", "0"),
                               legend.key.width = unit(1, "strwidth", "0"))


pal10 <- read.table("data/ref/colors/hcl_c70_vl.pal.10.csv", header=FALSE)
pal10 <- pal10[,1]
set.seed(3, kind = 'Mersenne-Twister')
pal10 <- sample(pal10)

reds <- read.table("data/ref/colors/red_colrs.csv", header=FALSE)
reds <- reds[,1]

source_colrs <- read.csv("data/ref/colors/source_colrs.csv", header=FALSE)
scmap <- source_colrs[,2]
names(scmap) <- source_colrs[,1]

tumor_colrs <- read.csv("data/ref/colors/tumor_colrs.csv", header=FALSE)
tummap <- tumor_colrs[,2]
names(tummap) <- tumor_colrs[,1]

#=========================================
# UMAP of integrated data
#=========================================

umap <- prolif@reductions$umap@cell.embeddings
datf <- prolif@meta.data
datf <- cbind(datf, umap[rownames(datf),])

ids <- c("integrated_snn_res.0.2", "source", "TumorStat")
leg_names <- c("Cluster", "source", "tumor")
out_fns <- paste0(dir_exp, c("UMAP.int_res0.2.png", "UMAP.source.png", 
                  "UMAP.tumor.png"))

cpal <- list()
cpal[[1]] <- hcl_pal(length(unique(prolif@meta.data[,"integrated_snn_res.0.2"])),
                     chr = c(100, 120), lum = c(60, 70), offset = 0,
                     rand = TRUE, seedn = 5)
cpal[[2]] <- scmap
cpal[[3]] <- tummap

for (i in 1:3){

    p <- ggplot(datf, aes_string(x="UMAP_1",y="UMAP_2",color=ids[i])) +
    geom_point(shape=16, size=.4) +
    scale_color_manual(values=cpal[[i]], name = leg_names[i]) +
    theme_p + 
    guides(color = guide_legend(override.aes = list(size=1))) 
    # geom_text_repel(data=meds, aes(x=UMAP_1, y=UMAP_2), label=rownames(meds), color="black") 

    ggsave(out_fns[i], width = 3.5, height = 3, dpi = 600)
    ggsave(sub("png$", "pdf", out_fns[i]), width = 3.5, height = 3)
}

#=========================================
# UMAP of singler classification
#========================================

umap <- prolif@reductions$umap@cell.embeddings
datf <- prolif@meta.data
datf <- cbind(datf, umap[rownames(datf),])

cts <- sort(unique(prolif@meta.data[,"cell_type_main"]))
cpal <- hcl_pal(length(cts), chr = c(100, 120), lum = c(50, 70), offset = 0, 
                rand = TRUE, seedn = 10)
names(cpal) <- cts

meds <- get_med_points(datf, c("UMAP_1", "UMAP_2"), "cell_type_main")

p <- ggplot(datf, aes(x=UMAP_1,y=UMAP_2,color=cell_type_main)) +
geom_point(shape=16, size=.4) +
scale_color_manual(values=cpal, name = "Cell type") +
theme_p + 
guides(color = guide_legend(override.aes = list(size=1))) 
# geom_text_repel(data=meds, aes(x=UMAP_1, y=UMAP_2), label=rownames(meds), color="black") 

outfn <- paste0(dir_exp, "UMAP.cell_type_main.png")
ggsave(outfn, width = 3.5, height = 3, dpi = 600)
outfn <- paste0(dir_exp, "UMAP.cell_type_main.pdf")
ggsave(outfn, width = 3.5, height = 3)

#=========================================
# Overlap of clusters and singler
#=========================================

datf <- prolif@meta.data

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
labs(x="Cluster", y="Cell type") + 
theme_classic() + 
theme(text = element_text(size = 8))

outfn <- paste0(dir_exp, "cl_ct_prop.pdf")
ggsave(outfn, width = 3.5, height = 2)

#=========================================
# Plot marker genes
#=========================================

dir_plt <- file.path(dir_exp, "ct_markers/")
create_dir(dir_plt)

reds <- read.table("data/ref/colors/red_colrs.csv", header=FALSE)
reds <- reds[,1]

rc <- prolif@assays$RNA@data
umap <- prolif@reductions$umap@cell.embeddings
datf <- prolif@meta.data
datf <- cbind(datf, umap[rownames(datf),])

hep_genes <- c("ASGR1", "APOA1", "CYP3A5", 
               "GLUL", "CYP1A2", # central
               "ALB", "PCK1", # portal
               "APOE") # central

endo_genes <- c("LYVE1", "CD14", "FCN3", # central/mid
                "BTNL9", # portal
                "CLEC4G", #LSEC
                "CD34", "PECAM") # MaVEC

t_gene <- c("CD8A", "CD8B", "CD4", 
            "GNLY", "NKG7", 
            "FOXP3")

mac_gene <- c("CD68", "S100A8", 
              "C1QA", "C1QB", "HLA-DRA", "CD163")

chol_gene <- c("MUC6", "KRT19", "EPCAM")

stell_gene <- c("ACTA2")

genes <- c(hep_genes, endo_genes, t_gene, mac_gene, chol_gene, stell_gene)

ids <- symb2ens[genes]
k <- which(ids %in% rownames(rc))
genes <- genes[k]
ids <- ids[k]
for (i in 1:length(genes)){
    message(genes[i])

    gene <- genes[i]
    id <- ids[i]
    datf[,"gene"] <- rc[id,rownames(datf)]
    o <- order(datf[,"gene"], decreasing=FALSE)

    p <- ggplot(datf[o,], aes(x = UMAP_1, y = UMAP_2, color = gene)) + 
    geom_point(shape=16, size = 0.25) + 
    theme_s + 
    ggtitle(gene) + 
    scale_color_gradientn(colors = reds, name = "logUMI")

    ggsave(paste0(dir_plt, "UMAP.", gene, ".png"), width = 3.5, height = 3, dpi = 300)
}

