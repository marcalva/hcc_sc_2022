
setwd("../../")

library(Seurat)
library(grid)
library(gtable)
library(ggplot2)
library(ggrepel)
library(scales)
library(gridExtra)
library(reshape2)
source("scripts/ggplot_raster.R")

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
dir_plt <- "exp/sharma_aiz/plots/";
dir.create(dir_plt, showWarnings=FALSE, recursive=TRUE)

fn <- "data/processed/sharma_aiz/liver.int_rand.rds"
integrated <- readRDS(fn)

source_colrs <- read.csv("data/ref/colors/source_colrs.csv", header=FALSE)
scmap <- source_colrs[,2]
names(scmap) <- source_colrs[,1]

tumor_colrs <- read.csv("data/ref/colors/tumor_colrs.csv", header=FALSE)
tummap <- tumor_colrs[,2]
names(tummap) <- tumor_colrs[,1]

vir_colrs <- read.csv("data/ref/colors/viral_colrs.csv", header=FALSE)
virmap <- vir_colrs[,2]
names(virmap) <- vir_colrs[,1]

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

#=========================================
# common plotting theme
#=========================================

theme_trnsp <- theme(plot.background = element_rect(fill="transparent", color=NA),
                     panel.background = element_rect(fill="transparent", color=NA),
                     gend.background = element_rect(fill="transparent", color=NA))

theme_txt <- theme(text = element_text(size = 8),
                   plot.title = element_text(hjust = 0.5))

theme_leg <- theme(legend.key.height = unit(1, "strheight", "0"),
                   legend.key.width = unit(1, "strwidth", "0"))

theme_axs <- theme(axis.text=element_blank(),
                   axis.ticks=element_blank())

theme_s <- theme_classic() + 
    theme_txt + theme_leg + theme_axs

#=========================================
# merge UMAP and meta data
#=========================================

udf <- integrated@reductions$umap@cell.embeddings
colnames(udf) <- c("UMAP1", "UMAP2")
udf <- cbind(udf, integrated@meta.data[rownames(udf),])

# change viral coding
udf[udf[,"ViralvsNonViral"] == 0, "ViralvsNonViral"] <- "NonViral"
udf[udf[,"ViralvsNonViral"] == 1, "ViralvsNonViral"] <- "Viral"
udf[,"ViralvsNonViral"] <- factor(udf[,"ViralvsNonViral"])

set.seed(1, kind = 'Mersenne-Twister')
udf <- udf[sample(1:nrow(udf)),]

k <- udf[,"source"] == "Aizarani"
udf_aiz <- udf[k,]

k <- udf[,"source"] == "nash_hcc"
udf_nsh <- udf[k,]

k <- udf[,"source"] == "Sharma"
udf_shr <- udf[k,]

rm(integrated)
gc()

#=========================================
# UMAP with points colored by source
#=========================================

udf[,"source"] <- factor(udf[,"source"])
scmap <- scmap[names(scmap) %in% levels(udf[,"source"])]

p <- ggplot(udf, aes(x=UMAP1,y=UMAP2,color=source)) + 
geom_point(shape=16, size=.01) + 
scale_color_manual(values=scmap, name="source") +  
theme_s + 
guides(color = guide_legend(override.aes = list(size=1)))
# theme_s + theme(legend.key.width = unit(1,"in"))

g <- raster_ggpoints(ggplotGrob(p), w=3, h=3, res=600)
outfn <- file.path(dir_plt, "UMAP.source.pdf")
pdf(outfn, width = 3.5, height = 3)
grid.draw(g)
dev.off()

#=========================================
# Plot integrated 0.5
#=========================================

meds <- get_med_points(udf, c("UMAP1", "UMAP2"), "integrated_snn_res.0.5")

p <- ggplot(udf, aes(x=UMAP1,y=UMAP2,color=integrated_snn_res.0.5)) +
geom_point(shape=16, size=.01) +
scale_color_manual(values=pal30, name="Re-clustered") +
ggtitle("NASH-HCC + Aizarani + Sharma") + 
theme_s + 
guides(color = guide_legend(override.aes = list(size=1))) + 
geom_text_repel(data=meds, aes(x=UMAP1, y=UMAP2), label=rownames(meds), color="black") 

g <- raster_ggpoints(ggplotGrob(p), w=3, h=3, res=600)
outfn <- file.path(dir_plt, "UMAP.int_res0.5.pdf")
pdf(outfn, width = 3.5, height = 3)
grid.draw(g)
dev.off()

#=========================================
# Plot integrated 1
#=========================================

meds <- get_med_points(udf, c("UMAP1", "UMAP2"), "integrated_snn_res.1")

p <- ggplot(udf, aes(x=UMAP1,y=UMAP2,color=integrated_snn_res.1)) +
geom_point(shape=16, size=.01) +
scale_color_manual(values=pal40, name="Re-clustered") +
ggtitle("NASH-HCC + Aizarani + Sharma") + 
theme_s + 
guides(color = guide_legend(override.aes = list(size=1))) + 
geom_text_repel(data=meds, aes(x=UMAP1, y=UMAP2), label=rownames(meds), color="black") 

g <- raster_ggpoints(ggplotGrob(p), w=3, h=3, res=600)
outfn <- file.path(dir_plt, "UMAP.int_res1.pdf")
pdf(outfn, width = 3.5, height = 3)
grid.draw(g)
dev.off()

#=========================================
# Plot G2M Score
#=========================================

feat <- "G2M.Score"; lt <- "G2M Score"
o <- order(udf[,"G2M.Score"], decreasing=FALSE)
p <- ggplot(udf[o,], aes(x=UMAP1, y=UMAP2, color=G2M.Score)) + 
    geom_point(size = 0.01, shape = 16) + 
    theme_s + ggtitle("NASH-HCC + Aizarani + Sharma") + 
    scale_color_gradientn(colours = reds, name =  lt)

g <- raster_ggpoints(ggplotGrob(p), w=3, h=3, res=600)
outfn <- file.path(dir_plt, "UMAP.G2M.Score.pdf")
pdf(outfn, width = 3.5, height = 3)
grid.draw(g)
dev.off()

#=========================================
# Plot S Score
#=========================================

feat <- "S.Score"; lt <- "S Score"
o <- order(udf[,"S.Score"], decreasing=FALSE)
p <- ggplot(udf[o,], aes(x=UMAP1, y=UMAP2, color=S.Score)) + 
    geom_point(size = 0.01, shape = 16) + 
    theme_s + ggtitle("NASH-HCC + Aizarani + Sharma") + 
    scale_color_gradientn(colours = reds, name =  lt)

g <- raster_ggpoints(ggplotGrob(p), w=3, h=3, res=600)
outfn <- file.path(dir_plt, "UMAP.S.Score.pdf")
pdf(outfn, width = 3.5, height = 3)
grid.draw(g)
dev.off()

#=========================================
# Plot G2M in each cohort separately
#=========================================

r <- range(udf[,"G2M.Score"])
xr <- range(udf[,"UMAP1"])
yr <- range(udf[,"UMAP2"])
theme_t <- theme(legend.direction="horizontal", legend.key.width = unit(.2, "inches"))

udf_list <- list(udf_nsh, udf_aiz, udf_shr)
udf_names <- list("NASH-HCC", "Aizarani", "Sharma")
p_l <- lapply(1:length(udf_list), function(ix){
              udf_x <- udf_list[[ix]]
              o <- order(udf_x[,"G2M.Score"], decreasing=FALSE)
              p <- ggplot(udf_x[o,], aes(x = UMAP1, y = UMAP2, color = G2M.Score)) + 
              geom_point(size=0.01, shape = 16) + 
              lims(x = xr, y = yr) + 
              theme_s + theme_t + 
              ggtitle(udf_names[[ix]]) + 
              scale_color_gradientn(colours=reds, name = "G2M Score", limits=r)
              return(p) })


grob1 <- ggplotGrob(p_l[[1]])
grob_leg <- grob1$grobs[[which(grob1$layout$name == "guide-box")]]
grob_leg <- gtable_trim(grob_leg)
grob_l <- lapply(p_l, function(p) ggplotGrob(p + theme(legend.position="none")))

grob_l <- lapply(grob_l, function(g){
              g <- raster_ggpoints(g, w=3, h=3, res=600)
              return(g)} )

gtbl <- gtable(widths = unit(c(2,2,2), "inches"), heights = unit(c(2, 0.5), "inches"))
gtbl <- gtable_add_grob(gtbl, grob_l[[1]], t = 1, l = 1, b = 1, r = 1, name="p1")
gtbl <- gtable_add_grob(gtbl, grob_l[[2]], t = 1, l = 2, b = 1, r = 2, name="p2")
gtbl <- gtable_add_grob(gtbl, grob_l[[3]], t = 1, l = 3, b = 1, r = 3, name="p3")
gtbl <- gtable_add_grob(gtbl, grob_leg, t = 2, l = 1, b = 2, r = 3, name="leg")

outfn <- file.path(dir_plt, "UMAP.sep.G2M.Score.pdf")
pdf(outfn, width = 3.5, height = 3)
grid.draw(gtbl)
dev.off()

#=========================================
# plot tumor/non-tumor
#=========================================

p <- ggplot(udf, aes(x=UMAP1,y=UMAP2,color=TumorStat)) + 
geom_point(shape=16, size=.01) + 
scale_color_manual(values=tummap, name="Tumor status") +  
theme_s + 
guides(color = guide_legend(override.aes = list(size=1)))

g <- raster_ggpoints(ggplotGrob(p), w=3, h=3, res=600)
outfn <- file.path(dir_plt, "UMAP.TumorStat.pdf")
pdf(outfn, width = 3.5, height = 3)
grid.draw(g)
dev.off()

#=========================================
# plot patient
#=========================================

p <- ggplot(udf, aes(x=UMAP1,y=UMAP2,color=PatientStat)) + 
geom_point(shape=16, size=.01) + 
scale_color_manual(values=pal20, name="Patient") +  
theme_s + 
guides(color = guide_legend(override.aes = list(size=1)))

g <- raster_ggpoints(ggplotGrob(p), w=3, h=3, res=600)
outfn <- file.path(dir_plt, "UMAP.PatientStat.pdf")
pdf(outfn, width = 3.5, height = 3)
grid.draw(g)
dev.off()

#=========================================
# plot viral stat
#=========================================

# tumor samples only
k <- udf[,"TumorStat"] == "Tumor"
p <- ggplot(udf[k,], aes(x=UMAP1,y=UMAP2,color=ViralvsNonViral)) + 
geom_point(shape=16, size=.01) + 
scale_color_manual(values=virmap, name="ViralvsNonViral") +  
theme_s + ggtitle("Integrated (tumor only)") + 
guides(color = guide_legend(override.aes = list(size=1)))

g <- raster_ggpoints(ggplotGrob(p), w=3, h=3, res=600)
outfn <- file.path(dir_plt, "UMAP.VirStat.pdf")
pdf(outfn, width = 3.5, height = 3)
grid.draw(g)
dev.off()

#=========================================
# plot louvain cluster from Sharma
#=========================================

k <- udf[,"source"] == "Sharma"
udf[,"louvain"] <- as.character(udf[,"louvain"])
meds <- get_med_points(udf[k,], c("UMAP1", "UMAP2"), "louvain")
p <- ggplot(udf[k,], aes(x=UMAP1,y=UMAP2,color=louvain)) + 
geom_point(shape=16, size=.01) + 
scale_color_manual(values=pal30, name="Sharma_louvain") +  
theme_s + ggtitle("Sharma louvain clusters") + 
guides(color = guide_legend(override.aes = list(size=1))) + 
geom_text_repel(data=meds, aes(x=UMAP1, y=UMAP2), 
                label=rownames(meds), color="black") 

g <- raster_ggpoints(ggplotGrob(p), w=3, h=3, res=600)
outfn <- file.path(dir_plt, "UMAP.Sharma_louvain.pdf")
pdf(outfn, width = 3.5, height = 3)
grid.draw(g)
dev.off()

#=========================================
# counts
#=========================================

udf[,"log_umi"] <- log10(udf[,"nCount_RNA"])

p <- ggplot(udf, aes(x=UMAP1, y=UMAP2, color=nCount_RNA)) + 
    geom_point(size = 0.01, shape = 16) + 
    theme_s + ggtitle("Total UMIs per droplet") + 
    scale_color_gradientn(colours = reds, name =  "UMIs", trans = "log10")

g <- raster_ggpoints(ggplotGrob(p), w=3, h=3, res=600)
outfn <- file.path(dir_plt, "UMAP.log_umis.pdf")
pdf(outfn, width = 3.5, height = 3)
grid.draw(g)
dev.off()

#=========================================
# sample-cell type frequency heatmap
#=========================================

k <- udf[,"source"] == "Sharma"
udf[,"louvain"] <- as.character(udf[,"louvain"])
tb <- unlist(table(udf[k,"integrated_snn_res.1"], udf[k,"louvain"]))
tb_prop <- sweep(tb, 1, rowSums(tb), '/')
tb_prop_m <- melt(tb_prop)
colnames(tb_prop_m) <- c("Cluster", "Sharm_louvain", "Prop")
tb_prop_m[,"Sharm_louvain"] <- as.character(tb_prop_m[,"Sharm_louvain"])
tb_prop_m[,"Cluster"] <- as.character(tb_prop_m[,"Cluster"])
for (i in 1:2){
    lv <- as.character(sort(unique(tb_prop_m[,i])))
    tb_prop_m[,i] <- factor(tb_prop_m[,i], levels = lv)
}


p <- ggplot(tb_prop_m, aes(x = Sharm_louvain, y = Cluster, fill = Prop)) + 
geom_tile(color = "gray") + 
scale_fill_gradient(low = "ghostwhite", high = "firebrick", limits=c(0,1)) + 
theme_bw() + theme(plot.title = element_text(hjust=0.5)) + 
ggtitle("Cluster freq per Sharma louvain")
outfn <- file.path(dir_plt, paste0("res1.sharma_louvain.freq.pdf"))
ggsave(outfn, width = 8, height = 7)

#=========================================
# sample-cell type frequency heatmap
#=========================================

tb <- unlist(table(udf$PatientStat, udf$integrated_snn_res.1))
tb_prop <- sweep(tb, 1, rowSums(tb), '/')
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
outfn <- file.path(dir_plt, paste0("res1.sample.freq.pdf"))
ggsave(outfn, width = 8, height = 8)

#=========================================
# tumor-cell type frequency heatmap
#=========================================

tb <- unlist(table(udf$integrated_snn_res.1, udf$TumorStat))
tb_prop <- sweep(tb, 1, rowSums(tb), '/')
tb_prop_m <- melt(tb_prop)
colnames(tb_prop_m) <- c("Cluster", "Tumor", "Prop")
tb_prop_m[,"Cluster"] <- as.character(tb_prop_m[,"Cluster"])

o <- order(tb_prop[,"NonTumor"])
lv <- as.character(rownames(tb_prop)[o])
tb_prop_m[,"Cluster"] <- factor(tb_prop_m[,"Cluster"], levels = lv)
lv <- as.character(sort(unique(tb_prop_m[,"Tumor"])))
tb_prop_m[,"Tumor"] <- factor(tb_prop_m[,"Tumor"], levels = lv)


p <- ggplot(tb_prop_m, aes(x = Tumor, y = Cluster, fill = Prop)) + 
geom_tile(color = "gray") + 
scale_fill_gradient(low = "ghostwhite", high = "firebrick", limits=c(0,1)) + 
theme_bw() + theme(plot.title = element_text(hjust=0.5)) + 
ggtitle("Cluster freq in tumors")
outfn <- file.path(dir_plt, paste0("res1.tumor.freq.pdf"))
ggsave(outfn, width = 3, height = 8)

#=========================================
# source-cell type frequency heatmap
#=========================================

tb <- unlist(table(udf$integrated_snn_res.1, udf$source))
tb_prop <- sweep(tb, 1, rowSums(tb), '/')
tb_prop_m <- melt(tb_prop)
colnames(tb_prop_m) <- c("Cluster", "Source", "Prop")
tb_prop_m[,"Cluster"] <- as.character(tb_prop_m[,"Cluster"])

lv <- as.character(rownames(tb_prop))
tb_prop_m[,"Cluster"] <- factor(tb_prop_m[,"Cluster"], levels = lv)
lv <- as.character(sort(unique(tb_prop_m[,"Source"])))
tb_prop_m[,"Source"] <- factor(tb_prop_m[,"Source"], levels = lv)


p <- ggplot(tb_prop_m, aes(x = Source, y = Cluster, fill = Prop)) + 
geom_tile(color = "gray") + 
scale_fill_gradient(low = "ghostwhite", high = "firebrick", limits=c(0,1)) + 
theme_bw() + theme(plot.title = element_text(hjust=0.5)) + 
ggtitle("Cluster freq per source")
outfn <- file.path(dir_plt, paste0("res1.source.freq.pdf"))
ggsave(outfn, width = 3, height = 8)

