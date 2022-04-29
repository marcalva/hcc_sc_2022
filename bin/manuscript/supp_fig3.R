
# plot inter-tumor heterogeneity
# plot source heterogeneity

setwd("../../")

library(Seurat)
library(reshape2)
library(ggplot2)
library(ggrepel)
library(gtable)
library(grid)
source("scripts/ggplot_raster.R")
source("scripts/color_pal.R")
source("scripts/gtable_stack.R")

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

# read data
fn <- "data/processed/sct_merge_3dat/lvr.sct_merge_3dat.rds"
seur <- readRDS(fn)

# colors
reds <- read.csv("data/ref/colors/red_colrs.csv", header=FALSE)
reds <- reds[,1]

source_colrs <- read.csv("data/ref/colors/source_colrs.csv", header=FALSE)
scmap <- source_colrs[,2]
names(scmap) <- source_colrs[,1]
names(scmap)[names(scmap) == "nash_hcc"] <- "NASH"

tumor_colrs <- read.csv("data/ref/colors/tumor_colrs.csv", header=FALSE)
tummap <- tumor_colrs[,2]
names(tummap) <- tumor_colrs[,1]

vir_colrs <- read.csv("data/ref/colors/viral_colrs.csv", header=FALSE)
virmap <- vir_colrs[,2]
names(virmap) <- vir_colrs[,1]

# Set directories
dir_plt <- "exp/manuscript/"
dir.create(dir_plt, showWarnings=FALSE, recursive=TRUE)

mar_pct <- 0.025
tagl <- 3

#=========================================
# merge UMAP and meta data
#=========================================

udf <- seur@reductions$umap@cell.embeddings
colnames(udf) <- c("UMAP1", "UMAP2")
udf <- cbind(udf, seur@meta.data[rownames(udf),])

# change viral coding
udf[udf[,"ViralvsNonViral"] == 0, "ViralvsNonViral"] <- "NonViral"
udf[udf[,"ViralvsNonViral"] == 1, "ViralvsNonViral"] <- "Viral"
udf[,"ViralvsNonViral"] <- factor(udf[,"ViralvsNonViral"])

set.seed(1, kind = 'Mersenne-Twister')
udf <- udf[sample(1:nrow(udf)),]

k <- udf[,"source"] == "nash_hcc"
udf[k,"source"] <- "NASH"

rm(seur)
gc()

#=========================================
# plot UMAPs
#=========================================

th_um <- theme(text = element_text(size = 8),
               plot.title = element_text(hjust = 0.5), 
               legend.key.size = unit(6, "points"), 
               legend.spacing.x = unit(2,"points"), 
               legend.spacing.y = unit(5,"points"), 
               legend.title = element_text(vjust = 0.5, hjust = 0), 
               legend.margin = margin(0,0,0,0), 
               legend.text = element_text(size = 4, hjust = 0, margin = margin(0,5,0,0)),
               legend.box.spacing = unit(0.5, "strheight", "0"), 
               legend.position = "right", 
               axis.text=element_blank(),
               axis.ticks=element_blank())
th_um <- theme_classic() + th_um

get_colors <- function(datf, n){
    u <- sort(unique(datf[,n]))
    l <- length(u)
    hcl_pal(l, chr = c(80, 80), lum = c(60, 80), 
            offset = 0, rand = TRUE, seedn = 1)
}

labs <- c("SCT_snn_res.1", "PatientStat", "source")
names(labs) <- labs
leg <- c("Cluster", "Patient", "Source")

p_cols <- get_colors(udf, "PatientStat")
c_cols <- get_colors(udf, "SCT_snn_res.1")
cols_l <- list(c_cols, p_cols, scmap)
tags <- c("a", "b", "c")

meds_clust <- get_med_points(udf, c("UMAP1", "UMAP2"), "SCT_snn_res.1")

umap_l <- list()

for (i in 1:length(labs)){
    lab <- labs[i]

    p <- ggplot(udf, aes_string(x="UMAP1",y="UMAP2",color=lab)) +
    geom_point(shape=16, size=.01) +
    scale_color_manual(values=cols_l[[i]], name=leg[i]) +
    ggtitle(NULL) + 
    th_um + 
    coord_fixed() + 
    guides(color = guide_legend(override.aes = list(size=1)))

    if (i == 1){
        p <- p + 
        geom_text_repel(data=meds_clust, aes(x=UMAP1, y=UMAP2), 
                        label=rownames(meds_clust), size = 2, 
                        force = 12, seed = 1, 
                        segment.size = 0.2, color="black")
    }

    gt_i <- ggplotGrob(p)

    gt_i <- raster_ggpoints(gt_i, w = 3, h = 3)
    gt_i <- gtable_add_tag(gt_i, tags[i], fs=16, just=c(0,1), l=tagl)
    gt_i <- pad_plot(gt_i, t=mar_pct, r=mar_pct, b=mar_pct, l=mar_pct)

    umap_l[[i]] <- gt_i
}

# plot heatmaps of overlap
th_hm <- theme_classic() + 
theme(text = element_text(size = 10), 
      legend.key.size = unit(10, "points"), 
      axis.text = element_text(colour = "black", size = 6), 
      axis.text.y = element_text(margin=margin(0,1,0,0)), 
      axis.text.x = element_text(margin=margin(0,0,0,0), 
                                 angle = 90, hjust = 1, vjust = 0.5), 
      panel.border =  element_rect(color="black", fill=NA), 
      panel.background = element_blank(), 
      panel.grid = element_blank(), 
      plot.title = element_text(hjust = 0.5), 
      axis.ticks = element_line(colour = "black"))

# plot cluster-source overlap
tb <- unlist(table(udf[,"SCT_snn_res.1"], udf[,"source"]))
tb_prop <- sweep(tb, 1, rowSums(tb), '/')
tb_propm <- melt(tb_prop)
colnames(tb_propm) <- c("Cluster", "Source", "Overlap")
pmax <- apply(tb_prop, 1, max)

o <- order(apply(tb_prop, 1, max), decreasing=TRUE)
lev <- rownames(tb_prop)[o]
tb_propm[,"Cluster"] <- factor(tb_propm[,"Cluster"], levels=lev)

leg_name <- "Cluster\npropotion"

p <- ggplot(tb_propm, aes(x=Cluster, y=Source, fill=Overlap)) + 
geom_tile() + 
scale_x_discrete(expand = expansion(0)) + 
scale_y_discrete(expand = expansion(0)) + 
scale_fill_gradientn(colors = reds, limits=c(0,1),name=leg_name) + 
th_hm


gt_cs <- ggplotGrob(p)
gt_cs <- gtable_add_tag(gt_cs, "d", fs=16, just=c(0,0), l=tagl)
gt_cs <- pad_plot(gt_cs, t=mar_pct, r=mar_pct, b=mar_pct, l=mar_pct)


# get tumor-enriched clusters
sample_stat <- paste(udf[,"PatientStat"], udf[,"TumorStat"], sep="_")
udf[,"SampleStat"] <- sample_stat
hcc_sam <- setdiff(sample_stat, c("Aizarani_NonTumor", "Sharma_0_NonTumor"))
pat_u <- sort(unique(udf[,"PatientStat"]))
hcc_pat <- setdiff(pat_u, c("Aizarani", "Sharma_0"))

tb <- unlist(table(udf[,"SCT_snn_res.1"], udf[,"PatientStat"]))
tb <- tb[,hcc_pat]
tb_prop <- sweep(tb, 1, rowSums(tb), '/')
tb_prop[is.na(tb_prop)] <- 0
tb_propm <- melt(tb_prop)
colnames(tb_propm) <- c("Cluster", "Patient", "Overlap")
pmax <- apply(tb_prop, 1, max)

o <- order(apply(tb_prop, 1, max), decreasing=TRUE)
lev <- rownames(tb_prop)[o]
tb_propm[,"Cluster"] <- factor(tb_propm[,"Cluster"], levels=lev)

leg_name <- "Cluster\npropotion"

p <- ggplot(tb_propm, aes(x=Cluster, y=Patient, fill=Overlap)) + 
geom_tile() +
labs(y = "HCC patient") + 
scale_x_discrete(expand = expansion(0)) + 
scale_y_discrete(expand = expansion(0)) + 
scale_fill_gradientn(colors = reds, limits=c(0,1), name = leg_name) + 
th_hm

gt_cp <- ggplotGrob(p)
gt_cp <- gtable_add_tag(gt_cp, "e", fs=16, just=c(0,0), l=tagl)
gt_cp <- pad_plot(gt_cp, t=mar_pct, r=mar_pct, b=mar_pct, l=mar_pct)


# get the percent of clusters with patient-specificity > 0.95
message("Number of patient-specific clusters (>0.95 from one patient):")
print(table(pmax > 0.95))

gt_heat <- gtable(widths = unit(c(1), "npc"), 
                  heights = unit(c(.4,.6), "npc"))
gt_heat <- gtable_add_grob(gt_heat, gt_cs, t=1,b=1,l=1,r=1)
gt_heat <- gtable_add_grob(gt_heat, gt_cp, t=2,b=2,l=1,r=1)

#=========================================
# put together a gtable for plotting everything
#=========================================

punit <- "inches"

ws <- c(7)
hs <- c(4, 4)

gtop <- gtable(widths = unit(c(0.6,0.4), "npc"), 
               heights = unit(c(0.5,0.5), "npc"))
gtop <- gtable_add_grob(gtop, umap_l[[1]], t=1,b=2,l=1,r=1)
gtop <- gtable_add_grob(gtop, umap_l[[2]], t=1,b=1,l=2,r=2)
gtop <- gtable_add_grob(gtop, umap_l[[3]], t=2,b=2,l=2,r=2)

gt <- gtable(widths = unit(ws, punit), 
             heights = unit(hs, punit))

#gt <- gtable_add_grob(gt, umap_l[[1]], t=1,b=1,l=1,r=1)
#gt <- gtable_add_grob(gt, umap_l[[2]], t=1,b=1,l=2,r=2)
# gt <- gtable_add_grob(gt, umap_l[[3]], t=1,b=1,l=3,r=3)
gt <- gtable_add_grob(gt, gtop, t=1,b=1,l=1,r=1)
gt <- gtable_add_grob(gt, gt_heat, t=2,b=2,l=1,r=1)

pdf(paste0(dir_plt, "FigureS3.pdf"), width = sum(ws), height = sum(hs))
grid.draw(gt)
dev.off()

