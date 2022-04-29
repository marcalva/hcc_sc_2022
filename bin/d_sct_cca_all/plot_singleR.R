
setwd("../../")

library(Seurat)
library(ggplot2)
library(ggrepel)
library(scales)
library(gridExtra)
library(reshape2)

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

pal10 <- read.csv("data/ref/colors/hcl_c70_vl.pal.10.csv", header=FALSE)
pal10 <- pal10[,1]
set.seed(2, kind = 'Mersenne-Twister')
pal10 <- sample(pal10)

pal20 <- read.csv("data/ref/colors/hcl_c70_vl.pal.20.csv", header=FALSE)
pal20 <- pal20[,1]
set.seed(2, kind = 'Mersenne-Twister')
pal20 <- sample(pal20)

pal30 <- read.csv("data/ref/colors/hcl_c70_vl.pal.30.csv", header=FALSE)
pal30 <- pal30[,1]
set.seed(2, kind = 'Mersenne-Twister')
pal30 <- sample(pal30)

# Set directories
dir_plt <- "exp/d_sct_cca_all/singleR/"; create_dir(dir_plt)

infn <- "data/processed/d_sct_cca_all/merge/seur.cca.rds"
seur <- readRDS(infn)

# replace NASH with normal
seur$SampleID <- sub("NASH", "Normal", seur$SampleID)
seur$SampleID <- sub("Normal", "Non-tumor", seur$SampleID)

# order factors from tumor to non-tumor
datf_s <- seur@meta.data[,c("SampleID", "Patient", "Tumor")]
datf_s[,"SampleID"] <- as.character(datf_s[,"SampleID"])
datf_s <- unique(datf_s)
o <- order(datf_s$Tumor, datf_s$Patient)

lv <- datf_s[o,"SampleID"]

seur$SampleID <- factor(seur$SampleID, levels = lv)

# read SingleR results
aiz_clsfy <- readRDS("exp/d_sct_cca_all/singleR/singleR.rds")

# get data frame for plotting
udf <- seur@reductions$umap@cell.embeddings
colnames(udf) <- c("UMAP1", "UMAP2")
udf <- cbind(udf, seur@meta.data[rownames(udf),])
udf[,"aiz_fine"] <- aiz_clsfy[["aiz_fine"]][rownames(udf),"pruned.labels"]
udf[,"aiz_main"] <- aiz_clsfy[["aiz_crse"]][rownames(udf),"pruned.labels"]

k <- (!is.na(udf[,"aiz_fine"])) & (!is.na(udf[,"aiz_main"]))
udf <- udf[k,]

#=========================================
# plotting theme
#=========================================

theme_txt <- theme(text = element_text(size = 8),
                   plot.title = element_text(hjust = 0.5))

theme_leg <- theme(legend.key.height = unit(1, "strheight", "0"),
                   legend.key.width = unit(1, "strwidth", "0"), 
                   legend.text = element_text(size = 6))

theme_axs <- theme(axis.text=element_blank(),
                   axis.ticks=element_blank())

theme_s <- theme_classic() + 
    theme_txt + theme_leg + theme_axs

#=========================================
# Plot
#=========================================

# Cell type labels
labs <- c("aiz_fine", "aiz_main")

# aiz_fine
l <- labs[1]
meds <- get_med_points(udf, c("UMAP1", "UMAP2"), l)

p <- ggplot(udf, aes_string(x="UMAP1",y="UMAP2",color=l)) +
geom_point(shape=16, size=.01) +
scale_color_manual(values=pal30, name=l) +
ggtitle("NASH-HCC") + 
theme_s + 
guides(color = guide_legend(override.aes = list(size=1), ncol=1)) + 
geom_text_repel(data=meds, aes(x=UMAP1, y=UMAP2), size = 3, 
                label=rownames(meds), color="black") 

outfn <- file.path(dir_plt, paste0("UMAP.", l, ".pdf"))
ggsave(outfn, width = 4.5, height = 4, dpi=300)
pname <- gsub("pdf$", "png", outfn)
ggsave(pname, width = 4.5, height = 4, dpi=300)

# aiz_main
l <- labs[2]
meds <- get_med_points(udf, c("UMAP1", "UMAP2"), l)

p <- ggplot(udf, aes_string(x="UMAP1",y="UMAP2",color=l)) +
geom_point(shape=16, size=.01) +
scale_color_manual(values=pal10, name=l) +
ggtitle("NASH-HCC") + 
theme_s + 
guides(color = guide_legend(override.aes = list(size=1), ncol=1)) + 
geom_text_repel(data=meds, aes(x=UMAP1, y=UMAP2), size = 3, 
                label=rownames(meds), color="black") 

outfn <- file.path(dir_plt, paste0("UMAP.", l, ".pdf"))
ggsave(outfn, width = 4.5, height = 4, dpi=300)
pname <- gsub("pdf$", "png", outfn)
ggsave(pname, width = 4.5, height = 4, dpi=300)

