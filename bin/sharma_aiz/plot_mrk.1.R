
# UMAP plot of cell type markers

setwd("../../")

library(Seurat)
library(ggplot2)
library(RColorBrewer)
source("scripts/ggplot_raster.R")

#=========================================
# Functions
#=========================================

create_dir <- function(p){
    dir.create(p, showWarnings=FALSE, recursive=TRUE)
}

#=========================================
# Set variables
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

# Set directories
dir_dat <- "data/processed/sharma_aiz/"; create_dir(dir_dat)

dir_exp <- "exp/sharma_aiz/"; create_dir(dir_exp)
dir_plt <- "exp/sharma_aiz/plots/genes/"; create_dir(dir_plt)

# colors
blues <- read.table("data/ref/colors/blue_colrs.csv", header=FALSE)
blues <- blues[,1]

reds <- read.table("data/ref/colors/red_colrs.csv", header=FALSE)
reds <- reds[,1]

#=========================================
# Map cell types
#=========================================

# Read in Seurat object
seur <- readRDS(paste0(dir_dat, "liver.int_rand.rds"))

Idents(seur) <- "integrated_snn_res.1"
ct_col <- "integrated_snn_res.1"

seur@active.assay <- "RNA"
seur <- NormalizeData(seur, 
                      normalization.method = "LogNormalize", 
                      scale.factor = 1000)

rc <- seur@assays$RNA@data
umap <- seur@reductions$umap@cell.embeddings
datf <- seur@meta.data
datf <- cbind(datf, umap[rownames(datf),])

set.seed(1, kind = 'Mersenne-Twister')
o <- sample(colnames(rc))
rc <- rc[,o]
datf <- datf[o,]

#=========================================
# Markers
#=========================================

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
            "FOXP3", 
            "CD2", "CD3E", "CD3G", "CD3D")

mac_gene <- c("CD68", "S100A8", 
              "C1QA", "C1QB", "HLA-DRA", "CD163")

chol_gene <- c("MUC6", "KRT19", "EPCAM")

stell_gene <- c("ACTA2")

genes <- c(hep_genes, endo_genes, t_gene, mac_gene, chol_gene, stell_gene)

ids <- symb2ens[genes]
k <- which(ids %in% rownames(seur))
genes <- genes[k]
ids <- ids[k]

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

for (i in 1:length(genes)){
    message(genes[i])

    gene <- genes[i]
    id <- ids[i]
    datf[,"gene"] <- rc[id,rownames(datf)]
    o <- order(datf[,"gene"], decreasing=FALSE)

    # umap plot
    p <- ggplot(datf, aes(x = UMAP_1, y = UMAP_2, color = gene)) + 
    geom_point(shape=16, size = 0.01) + 
    theme_s + 
    ggtitle(gene) + 
    scale_color_gradientn(colors = reds, name = "logUMI")

    g <- raster_ggpoints(ggplotGrob(p), w=3, h=3, res=600)
    outfn <- paste0(dir_plt, "UMAP.", gene, ".pdf")
    pdf(outfn, width = 3.5, height = 3)
    grid.draw(g)
    dev.off()

    # box plot of resolution 0.1
    p <- ggplot(datf, aes(x = integrated_snn_res.1, y = gene)) + 
    geom_boxplot(outlier.shape=16, outlier.size=0.1, outlier.alpha=0.1) + 
    theme_classic() + 
    ylab("logUMI") + 
    ggtitle(gene) + 
    theme(text = element_text(size = 8),
          axis.text.x = element_text(angle=90, hjust=1, size=8), 
          plot.title = element_text(hjust = 0.5))

    outfn <- paste0(dir_plt, "box.res.1.", gene, ".pdf")
    ggsave(outfn, width = 5, 
           height = 3, dpi = 300)
}
