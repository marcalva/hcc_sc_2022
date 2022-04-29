
setwd("../../")

library(Seurat)
library(grid)
library(gtable)
library(ggplot2)
library(ggrepel)
library(scales)
library(gridExtra)
library(reshape2)
suppressPackageStartupMessages(library(SingleR))
suppressPackageStartupMessages(library(SummarizedExperiment))

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

# Set directories
dir_plt <- "exp/aizarani2019/sct_clust/SingleR/"
dir.create(dir_plt, showWarnings=FALSE, recursive=TRUE)

fn <- "data/processed/aizarani2019/seur.rds"
seur <- readRDS(fn)

ref_clssfy <- readRDS("exp/aizarani2019/sct_clust/SingleR/singleR.ref.rds")

reds <- read.csv("data/ref/colors/red_colrs.csv", header=FALSE)
reds <- reds[,1]

pal20 <- read.csv("data/ref/colors/hcl_c70_vl.pal.20.csv", header=FALSE)
pal20 <- pal20[,1]
set.seed(1, kind = 'Mersenne-Twister')
pal20 <- sample(pal20)

pal30 <- read.csv("data/ref/colors/hcl_c70_vl.pal.30.csv", header=FALSE)
pal30 <- pal30[,1]
set.seed(1, kind = 'Mersenne-Twister')
pal30 <- sample(pal30)

pal40 <- read.csv("data/ref/colors/hcl_c70_vl.pal.40.csv", header=FALSE)
pal40 <- pal40[,1]
set.seed(1, kind = 'Mersenne-Twister')
pal40 <- sample(pal40)

pal50 <- read.csv("data/ref/colors/hcl_c70_vl.pal.50.csv", header=FALSE)
pal50 <- pal50[,1]
set.seed(1, kind = 'Mersenne-Twister')
pal50 <- sample(pal50)

#=========================================
# common plotting theme
#=========================================

theme_trnsp <- theme(plot.background = element_rect(fill="transparent", color=NA),
                     panel.background = element_rect(fill="transparent", color=NA),
                     gend.background = element_rect(fill="transparent", color=NA))

theme_txt <- theme(text = element_text(size = 8),
                   plot.title = element_text(hjust = 0.5))

theme_leg <- theme(legend.text = element_text(size = 6),
                   legend.key.height = unit(1, "strheight", "0"),
                   legend.key.width = unit(1, "strwidth", "0"))

theme_axs <- theme(axis.text=element_blank(),
                   axis.ticks=element_blank())

theme_s <- theme_classic() + 
    theme_txt + theme_leg + theme_axs

#=========================================
# merge UMAP and meta data
#=========================================

udf <- seur@reductions$umap@cell.embeddings
colnames(udf) <- c("UMAP1", "UMAP2")
udf <- cbind(udf, seur@meta.data[rownames(udf),])

set.seed(1, kind = 'Mersenne-Twister')
udf <- udf[sample(1:nrow(udf)),]

rm(seur)
gc()

# add classifications
udf[,"enc_fine"] <- ref_clssfy[["enc_fine"]][rownames(udf),"pruned.labels"]
udf[,"enc_main"] <- ref_clssfy[["enc_main"]][rownames(udf),"pruned.labels"]
udf[,"hpa_fine"] <- ref_clssfy[["hpa_fine"]][rownames(udf),"pruned.labels"]
udf[,"hpa_main"] <- ref_clssfy[["hpa_main"]][rownames(udf),"pruned.labels"]
udf[,"nov_main"] <- ref_clssfy[["nov_main"]][rownames(udf),"pruned.labels"]
udf[,"nov_fine"] <- ref_clssfy[["nov_fine"]][rownames(udf),"pruned.labels"]

#=========================================
# UMAP plots
#=========================================

labels <- c("enc_fine", "enc_main", "hpa_main", "nov_main", "nov_fine")
for (l in labels){
    k <- !is.na(udf[,l])
    udfs <- udf[k,]
    meds <- get_med_points(udfs, c("UMAP1", "UMAP2"), l)

    p <- ggplot(udfs, aes_string(x="UMAP1",y="UMAP2",color=l)) +
        geom_point(shape=16, size=0.01) + 
        scale_color_manual(values=pal40, name=l) + 
        ggtitle("Aizarani") +
        theme_s +
        guides(color = guide_legend(override.aes = list(size=1))) +
        geom_text_repel(data=meds, aes(x=UMAP1, y=UMAP2), size = 1, 
                        label=rownames(meds), color="black")

    outfn <- file.path(dir_plt, paste0("UMAP.", l, ".pdf"))
    ggsave(outfn, width = 7, height = 6, dpi=300)
    pname <- gsub("pdf$", "png", outfn)
    ggsave(pname, width = 7, height = 6, dpi=300)

}

#=========================================
# Heatmap
#=========================================

theme_s <- theme_classic() + 
    theme_leg + theme_txt + 
    theme(text = element_text(size = 8), 
          axis.text.x = element_text(angle=90, hjust=1))

labels <- c("enc_fine", "enc_main", "hpa_fine", "hpa_main", 
            "nov_main", "nov_fine")
for (l in labels){
    k <- !is.na(udf[,l])
    udfs <- udf[k,]

    tb <- unlist(table(udfs[,l], udfs[,"SCT_snn_res.1"]))
    tb_prop <- sweep(tb, 2, colSums(tb), '/')
    tb_prop_m <- melt(tb_prop)
    colnames(tb_prop_m) <- c(l, "res0.5", "Prop")
    for (i in 1:2){
        lv <- as.character(sort(unique(tb_prop_m[,i])))
        tb_prop_m[,i] <- factor(tb_prop_m[,i], levels = lv)
    }

    p <- ggplot(tb_prop_m, aes_string(x = "res0.5", y = l, fill = "Prop")) + 
        geom_tile(color = "gray") + 
        scale_fill_gradientn(colors=reds, limits=c(0,1), name=l) + 
        ggtitle(l) + 
        theme_s + 
        labs(x=NULL, y=NULL)

    outfn <- file.path(dir_plt, paste0("res0.5.", l, ".freq.pdf"))
    ggsave(outfn, width = 7, height = 10)
}


#=========================================
# Plot clust 0.5
#=========================================

meds <- get_med_points(udf, c("UMAP1", "UMAP2"), "SCT_snn_res.0.5")

p <- ggplot(udf, aes(x=UMAP1,y=UMAP2,color=SCT_snn_res.0.5)) +
geom_point(shape=16, size=.01) +
scale_color_manual(values=pal30, name="Res 0.5") +
ggtitle("Aizarani") + 
theme_s + 
guides(color = guide_legend(override.aes = list(size=1))) + 
geom_text_repel(data=meds, aes(x=UMAP1, y=UMAP2), label=rownames(meds), color="black") 

outfn <- file.path(dir_plt, "UMAP.res0.5.pdf")
ggsave(outfn, width = 3.5, height = 3, dpi=300)
pname <- gsub("pdf$", "png", outfn)
ggsave(pname, width = 3.5, height = 3, dpi=300)

#=========================================
# Plot clust res 1
#=========================================

meds <- get_med_points(udf, c("UMAP1", "UMAP2"), "SCT_snn_res.1")

p <- ggplot(udf, aes(x=UMAP1,y=UMAP2,color=SCT_snn_res.1)) +
geom_point(shape=16, size=.01) +
scale_color_manual(values=pal30, name="Res 1") +
ggtitle("Aizarani") + 
theme_s + 
guides(color = guide_legend(override.aes = list(size=1))) + 
geom_text_repel(data=meds, aes(x=UMAP1, y=UMAP2), label=rownames(meds), color="black") 

outfn <- file.path(dir_plt, "UMAP.res1.pdf")
ggsave(outfn, width = 3.5, height = 3, dpi=300)
pname <- gsub("pdf$", "png", outfn)
ggsave(pname, width = 3.5, height = 3, dpi=300)

#=========================================
# Plot Aizarani sample
#=========================================

p <- ggplot(udf, aes(x=UMAP1,y=UMAP2,color=aiz_sample)) +
geom_point(shape=16, size=.01) +
scale_color_manual(values=pal50, name="Sample") +
ggtitle("Aizarani") + 
theme_s + 
guides(color = guide_legend(override.aes = list(size=1), ncol=2)) 

outfn <- file.path(dir_plt, "UMAP.aiz_sample.pdf")
ggsave(outfn, width = 5, height = 4, dpi=300)
pname <- gsub("pdf$", "png", outfn)
ggsave(pname, width = 5, height = 4, dpi=300)
# system(paste("convert -density 300", outfn, pname))

#=========================================
# Plot Aizarani clust
#=========================================

udf[,"aiz_clust"] <- as.character(udf[,"aiz_clust"])
meds <- get_med_points(udf, c("UMAP1", "UMAP2"), "aiz_clust")

p <- ggplot(udf, aes(x=UMAP1,y=UMAP2,color=aiz_clust)) +
geom_point(shape=16, size=.01) +
scale_color_manual(values=pal40, name="Aiz cluster") +
ggtitle("Aizarani") + 
theme_s + 
guides(color = guide_legend(override.aes = list(size=1))) + 
geom_text_repel(data=meds, aes(x=UMAP1, y=UMAP2), label=rownames(meds), color="black") 

outfn <- file.path(dir_plt, "UMAP.aiz_clust.pdf")
ggsave(outfn, width = 3.5, height = 3, dpi=300)
pname <- gsub("pdf$", "png", outfn)
ggsave(pname, width = 3.5, height = 3, dpi=300)
# system(paste("convert -density 300", outfn, pname))

#=========================================
# Plot G2M Score
#=========================================

feat <- "G2M.Score"; lt <- "G2M Score"
o <- order(udf[,"G2M.Score"], decreasing=FALSE)
p <- ggplot(udf[o,], aes(x=UMAP1, y=UMAP2, color=G2M.Score)) + 
    geom_point(size = 0.01, shape = 16) + 
    theme_s + ggtitle("Aizarani") + 
    scale_color_gradientn(colours = reds, name =  lt)

outfn <- file.path(dir_plt, "UMAP.G2M.Score.pdf")
ggsave(outfn, width = 3.5, height = 3, dpi=300)
pname <- gsub("pdf$", "png", outfn)
ggsave(pname, width = 3.5, height = 3, dpi=300)

#=========================================
# Plot S Score
#=========================================

feat <- "S.Score"; lt <- "S Score"
o <- order(udf[,"S.Score"], decreasing=FALSE)
p <- ggplot(udf[o,], aes(x=UMAP1, y=UMAP2, color=S.Score)) + 
    geom_point(size = 0.01, shape = 16) + 
    theme_s + ggtitle("Aizarani") + 
    scale_color_gradientn(colours = reds, name =  lt)

outfn <- file.path(dir_plt, "UMAP.S.Score.pdf")
ggsave(outfn, width = 3.5, height = 3, dpi=300)
pname <- gsub("pdf$", "png", outfn)
ggsave(pname, width = 3.5, height = 3, dpi=300)

#=========================================
# counts
#=========================================

udf[,"log_umi"] <- log10(udf[,"nCount_RNA"])

p <- ggplot(udf, aes(x=UMAP1, y=UMAP2, color=nCount_RNA)) + 
    geom_point(size = 0.01, shape = 16) + 
    theme_s + ggtitle("Total UMIs per droplet") + 
    scale_color_gradientn(colours = reds, name =  "UMIs", trans = "log10")

outfn <- file.path(dir_plt, "UMAP.log_umis.pdf")
ggsave(outfn, width = 3.5, height = 3, dpi=300)
pname <- gsub("pdf$", "png", outfn)
ggsave(pname, width = 3.5, height = 3, dpi=300)


#=========================================
# sample-cell type frequency heatmap
#=========================================

tb <- unlist(table(udf$aiz_sample, udf$SCT_snn_res.0.5))
tb_prop <- sweep(tb, 2, colSums(tb), '/')
tb_prop_m <- melt(tb_prop)
colnames(tb_prop_m) <- c("Sample", "Cluster", "Prop")
tb_prop_m[,"Cluster"] <- as.character(tb_prop_m[,"Cluster"])
for (i in 1:2){
    lv <- as.character(sort(unique(tb_prop_m[,i])))
    tb_prop_m[,i] <- factor(tb_prop_m[,i], levels = lv)
}

p <- ggplot(tb_prop_m, aes(x = Cluster, y = Sample, fill = Prop)) + 
geom_tile(color = "gray") + 
scale_fill_gradient(low = "ghostwhite", high = "firebrick", limits=c(0,1)) + 
theme_bw() + theme(plot.title = element_text(hjust=0.5)) + 
ggtitle("Cluster freq per sample")
outfn <- file.path(dir_plt, paste0("aiz_sample.res0.5.ovrlp.pdf"))
ggsave(outfn, width = 7, height = 10)

#=========================================
# seurat to published clusters
#=========================================

fn <- "data/raw/aizarani2019/GSE124395_clusterpartition.txt"
aiz_clust <- read.table(fn, row.names=1, header=TRUE, check.names=FALSE)

md <- udf
i <- intersect(rownames(md), rownames(aiz_clust))

tb <- unlist(table(md[i,"SCT_snn_res.0.5"], aiz_clust[i,"sct@cpart"]))
tb_prop <- sweep(tb, 2, colSums(tb), '/')
tb_prop_m <- melt(tb_prop)
colnames(tb_prop_m) <- c("res0.5", "Aizarani", "Prop")
for (i in 1:2){
    lv <- as.character(sort(unique(tb_prop_m[,i])))
    tb_prop_m[,i] <- factor(tb_prop_m[,i], levels = lv)
}

p <- ggplot(tb_prop_m, aes(x = res0.5, y = Aizarani, fill = Prop)) + 
geom_tile(color = "gray") + 
scale_fill_gradient(low = "ghostwhite", high = "firebrick", limits=c(0,1)) + 
theme_bw() + theme(plot.title = element_text(hjust=0.5)) + 
ggtitle("Seurat Aizarani cluster comparison")
outfn <- file.path(dir_plt, paste0("res0.5.Aizarani.freq.pdf"))
ggsave(outfn, width = 7, height = 10)

