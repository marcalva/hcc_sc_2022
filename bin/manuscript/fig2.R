
# Figure 2

setwd("../../")

library(Seurat)
library(ggplot2)
library(grid)
library(gtable)
library(ggrepel)
library(scales)
library(gridExtra)
library(reshape2)
source("scripts/color_pal.R")
source("scripts/get_enr_mat.R")
source("scripts/ggplot_raster.R")
source("scripts/gtable_stack.R")
source("scripts/ggplot_formats.R")

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

mar_pct <- 0.025 # pct margin around panel
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

fn <- "data/processed/sharma_aiz/liver.int_rand.md_umap.txt.gz"
udf <- read.table(fn, header=TRUE, row.names=1, sep='\t')
colnames(udf)[colnames(udf) == "UMAP_1"] <- "UMAP1"
colnames(udf)[colnames(udf) == "UMAP_2"] <- "UMAP2"
set.seed(1, kind = 'Mersenne-Twister')
udf <- udf[sample(1:nrow(udf)),]


dir_plt <- "exp/manuscript/"
dir.create(dir_plt, showWarnings=FALSE, recursive=TRUE)

# tumor colors
fn <- "data/ref/colors/tumor_colrs.csv"
tum_cols <- read.table(fn, header=FALSE, sep=",")
tum_colsv <- tum_cols[,2]
names(tum_colsv) <- tum_cols[,1]
names(tum_colsv)[names(tum_colsv) == "Tumor"] <- "Tumor"
names(tum_colsv)[names(tum_colsv) == "NonTumor"] <- "Non-tumor"

# reds
reds <- read.csv("data/ref/colors/red_colrs.csv", header=FALSE)
reds <- reds[,1]

#=========================================
# theme
#=========================================

th_txt <- theme(text = element_text(size = 8, color = "black"),
                axis.text.y = element_text(size = 8, color="black"),
                axis.text.x = element_text(size = 8, color="black", angle = 45, hjust=0.9),
                axis.title = element_text(size = 8, color = "black"), 
                plot.title = element_text(hjust = 0.5, size = 10))

th_pbg <- theme(panel.border = element_rect(fill = NA, colour = "black"))

th_stat <- theme_bw() + th_txt + th_pbg + 
theme(legend.margin = margin(0,4,4,4), 
      legend.box.margin = margin(0,0,0,0), 
      legend.box.spacing = unit(0, "points"), 
      legend.box.background = element_rect(color="black"),
      legend.justification = c(1, 1), 
      legend.position = c(0.95, 0.95), 
      legend.text = element_text(size = 8, color = "black"), 
      legend.title = element_text(size = 8, hjust = 0.5), 
      legend.key.size = unit(0, units = "points"))


th_bar <- theme_bw() + th_txt + th_pbg + 
theme(legend.box.margin = margin(0,0,0,0), 
      legend.box.spacing = unit(0, "points"), 
      legend.position = "bottom", 
      legend.text = element_text(size = 8, color = "black"), 
      legend.title = element_text(size = 8, hjust = 0.5), 
      legend.spacing.x = unit(4, units = "points"), 
      legend.key.size = unit(16, units = "points"))

#=========================================
# TCGA paired test and proportions
#=========================================

fn <- "exp/tcga_hcc/ctp.cell_type_main/tcga.cell_type_main.tum_nontum.txt"
tcga.pair.test <- read.table(fn, header=TRUE, row.names=1)

tcga.pair.test[,"CellType"] <- rownames(tcga.pair.test)
o <- order(tcga.pair.test[,"t.statistic"], decreasing=TRUE)
tcga_ct_lev <- tcga.pair.test[,"CellType"][o]

tcga.pair.test[,"CellType"] <- factor(tcga.pair.test[,"CellType"], 
                                   levels = tcga_ct_lev)

tcga.pair.test[,"sig"] <- "ns"
tcga.pair.test[ tcga.pair.test[,"w.p_adj"] < 0.05 , "sig"] <- "sig"

fn <- "exp/tcga_hcc/ctp.cell_type_main/ctp_tumor_mlt.txt"
tmm.dfm <- read.table(fn, header=TRUE, sep="\t")

tmm.dfm[, "CellType"] <- factor(tmm.dfm[, "CellType"], 
                               levels = tcga_ct_lev)

#========================================================
# LCI paired test and proportions
#========================================================

cl_id <- "cell_type_main"

# Read paired test results
fn <- paste0("exp/GSE14520/ctp.", cl_id, "/lci.", cl_id, ".tum_nontum.txt")
lci.pair.test <- read.table(fn, header=TRUE, row.names=1, stringsAsFactors=FALSE, sep="\t")

ctrm <- c("B")
k <- lci.pair.test[,"CellType"] != ctrm
lci.pair.test <- lci.pair.test[k,,drop=FALSE]

lci.pair.test[,"CellType"] <- rownames(lci.pair.test)
lci.pair.test[,"CellType"] <- factor(lci.pair.test[,"CellType"], 
                                   levels = lci.pair.test[,"CellType"])
o <- order(lci.pair.test[,"t.statistic"], decreasing=TRUE)
lci_ct_lev <- lci.pair.test[,"CellType"][o]

lci.pair.test[,"CellType"] <- factor(lci.pair.test[,"CellType"], 
                                   levels = lci_ct_lev)

lci.pair.test[,"sig"] <- "ns"
lci.pair.test[ lci.pair.test[,"w.p_adj"] < 0.05 , "sig"] <- "sig"

fn <- "exp/GSE14520/ctp.cell_type_main/ctp_tumor_mlt.txt"
ctp.dfm <- read.table(fn, header=TRUE, sep="\t")
k <- ctp.dfm[,"CellType"] != ctrm
ctp.dfm <- ctp.dfm[k,,drop=FALSE]

o <- order(lci.pair.test[,"t.statistic"], decreasing=TRUE)
ctp.dfm[, "CellType"] <- factor(ctp.dfm[, "CellType"], 
                               levels = rownames(lci.pair.test)[o])

#=========================================
#=========================================

gr1_l <- list()
pt_l <- list(tcga.pair.test, lci.pair.test)
ctp_l <- list(tmm.dfm, ctp.dfm)
tags <- c("a", "b")
titles <- c("TCGA", "LCI")
for (i in 1:2){
    # t stats per cell type tcga
    paired.test <- pt_l[[i]]
    maxlim <- max(abs(paired.test[,"t.statistic"]))
    sig_shape <- c("ns" = 1, "sig" = 16)

    p <- ggplot(paired.test, aes(x = CellType, y = t.statistic)) + 
    geom_hline(yintercept = 0, col="red") + 
    geom_point(aes(shape = sig), size=1.5) + theme_bw() + 
    scale_shape_manual(values = sig_shape, breaks = "sig", 
                       name = NULL, 
                       labels = c("sig" = "adj p < 0.05")) + 
    ylim(c(-maxlim, maxlim) * 1.2) + 
    ylab("T statistic") + 
    xlab(NULL) + 
    ggtitle(titles[i]) +
    th_stat

    gr1_pt <- ggplotGrob(p)
    pan_r <- gr1_pt$layout[gr1_pt$layout$name == "panel", 't']
    gr1_pt_s <- gr1_pt[1:pan_r,]
    gr1_pt_s <- gtable_add_rows(gr1_pt_s, heights = unit(6, "points"), nrow(gr1_pt_s)) # pad space between top and bottom row

    # bar plot of scaled cell type proprs
    ctpm <- ctp_l[[i]]
    ylims <- boxplot.stats(ctpm[,"Proportion"])$stats[c(1,5)]
    p <- ggplot(ctpm, aes(x = CellType, y = Proportion, fill = Type)) + 
    geom_boxplot(outlier.shape=NA, size=0.2, width = 0.74) + 
    scale_fill_manual(values = tum_colsv, name=NULL) + 
    scale_y_continuous(limits = ylims) + 
    xlab(NULL) + ylab("Proportion estimates\n(scaled)") + 
    ggtitle("TCGA") + 
    th_bar

    gr1_bp<- ggplotGrob(p)
    pan_r <- gr1_bp$layout[gr1_bp$layout$name == "panel", 't']
    gr1_bp_s <- gr1_bp[pan_r:nrow(gr1_bp),]

    # combine plots
    gr_l <- list(gr1_pt_s, gr1_bp_s)
    gr1 <- stack_gtable_v(gr_l, heights = unit(c(0.4, 0.6), "npc"))

    gr1$grobs[[1]] <- gtable_add_tag(gr1$grobs[[1]], tags[i], fs=16, just=c(0.2,1), l=tagl)
    gr1 <- pad_plot(gr1, t=mar_pct, r=mar_pct, b=mar_pct, l=mar_pct)

    gr1_l[[i]] <- gr1
}

#================================================
# horizontal box plot funcion/theme
#================================================

th_hbx <- theme(text = element_text(colour="black", size=8),
                axis.title = element_text(colour="black", size=8), 
                axis.text = element_text(colour="black", size=8),
                plot.title = element_text(hjust=0.5), 
                axis.ticks.y = element_blank(), 
                legend.key.size = unit(6, "points"), 
                legend.justification=c(1,0),
                legend.background = element_rect(color="black", fill = NA, linetype="solid"),
                legend.margin = margin(0,2.5,2.5,2.5))

boxj_plot <- function(datf, x = "de_logFC", y = "cluster", color = "gws", 
                      xlab = NULL, ylab = NULL, 
                      title = NULL, gmap = NULL, leg_pos = c(1,.2)){
    # order y
    y_m <- tapply(datf[,x], datf[,y], mean)
    o <- order(y_m, decreasing=FALSE)
    lv <- names(y_m)[o]
    datf[,y] <- factor(datf[,y], levels = lv)

    p <- ggplot(datf, aes_string(x = x, y = y)) + 
    geom_vline(xintercept = 0, color = "black", linetype = "dotted") + 
    geom_boxplot(fill = NA, outlier.shape=NA) + 
    geom_jitter(aes_string(color = color), height = 0.2, shape=16, size=0.5) + 
    ylab(ylab) + 
    xlab(xlab) + 
    ggtitle(title) + 
    scale_color_manual(values = gmap, breaks = "adj p < 0.05", 
                       name=NULL) + 
    theme_classic() + th_hbx + 
    theme(legend.position = leg_pos)
    return(p)
}

#================================================
# read ct marker tumor DE data
#================================================

mrk_logfc_thr <- 0.5
gmap <- c(" " = "black", "adj p < 0.05" = "red")

# TCGA cell type markers with tumor DE results
fn <- "exp/tcga_hcc/de/cell_type_main.de_mrk.txt"
tcga_de <- read.table(fn, header=TRUE, sep='\t')
tcga_de <- tcga_de[tcga_de[,"avg_log2FC"] > mrk_logfc_thr,]

# LCI cell type markers with tumor DE results
fn <- "exp/GSE14520/de/cell_type_main.de_mrk.txt"
lci_de <- read.table(fn, header=TRUE, sep='\t')
lci_de <- lci_de[lci_de[,"avg_log2FC"] > mrk_logfc_thr,]

# ctrm <- c("B")
# k <- lci_de[,"cluster"] != ctrm
# lci_de <- lci_de[k,,drop=FALSE]

#================================================
# plot ct marker tumor DE
#================================================

gr2_l <- list()
de_l <- list(tcga_de, lci_de)
tags <- c("c", "d")
titles <- c("Tumor DE of cell-type markers\n in bulk TCGA", 
            "Tumor DE of cell-type markers\n in bulk LCI")
ylab <- "Cell-type marker genes"
xlab <- list(bquote("TCGA bulk tumor log"[2]*"FC"), 
             bquote("LCI bulk tumor log"[2]*"FC"))

for (i in 1:2){
    p <- boxj_plot(de_l[[i]], title = titles[i], gmap = gmap, 
                   color = NULL, 
                   xlab = xlab[[i]], ylab = ylab, 
                   leg_pos = c(1, .0)) + 
    theme(legend.position="none")

    grtmp <- ggplotGrob(p)

    grtmp <- gtable_add_tag(grtmp, tags[i], fs=16, just=c(0,1), l=tagl)
    grtmp <- pad_plot(grtmp, t=mar_pct, r=mar_pct, b=mar_pct, l=mar_pct)

    gr2_l[[i]] <- grtmp
}

#=========================================
# sc tumor scores TCGA
#=========================================

th_um <- theme_classic() + 
theme(plot.background = element_rect(fill="transparent", color=NA),
      panel.background = element_rect(fill="transparent", color=NA),
      legend.text = element_text(size = 6, angle=30, hjust=1), 
      legend.title = element_text(size = 6, vjust = 1, hjust = 0.5), 
      legend.position = "bottom", 
      legend.key.height = unit(6, "points"),
      legend.key.width = unit(12, "points"),
      legend.box.spacing = unit(2, "points"), 
      text = element_text(size = 8),
      plot.title = element_text(hjust = 0.5),
      axis.text=element_blank(),
      axis.ticks=element_blank())

th_vbox <- theme_classic() + 
theme(text = element_text(size = 8),
      axis.text = element_text(colour="black", size=8),
      axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size=8), 
      plot.title = element_text(hjust = 0.5))

#=========================================
# sc tumor scores TCGA
#=========================================

fn <- "exp/tcga_hcc/sharma_aiz.tum_score/sharma_aiz.tum_enr_scores.txt"
tcga_tum_enr <- read.table(fn, row.names=1, header=TRUE, sep = '\t')
fn <- "exp/tcga_hcc/sharma_aiz.tum_score/cell_type_main.diff_stat.txt"
tcga_de_stat <- read.table(fn, row.names=1, header=TRUE, sep = '\t')

fn <- "exp/GSE14520/sharma_aiz.tum_score/sharma_aiz.tum_enr_scores.txt"
lci_tum_enr <- read.table(fn, row.names=1, header=TRUE, sep = '\t')
fn <- "exp/GSE14520/sharma_aiz.tum_score/cell_type_main.diff_stat.txt"
lci_de_stat <- read.table(fn, row.names=1, header=TRUE, sep = '\t')

uxmax <- max(udf[,"UMAP1"])
uymin <- min(udf[,"UMAP2"])

gr3_l <- list()
tum_l <- list(tcga_tum_enr, lci_tum_enr)
tags <- list(c("e", "g"), c("f", "h"))
cohorts <- c("TCGA", "LCI")
utitles <- c("TCGA tumor-elevated\ngene scores", 
             "LCI tumor-elevated\ngene scores")
ctitles <- c("TCGA tumor-elevated\ngene scores", 
             "LCI tumor-elevated\ngene scores")
btitles <- c("TCGA\ntumor-elevated\ngene scores",  
             "LCI\ntumor-elevated\ngene scores")
de_stat_l <- list(tcga_de_stat, lci_de_stat)

for (i in 1:2){
    pdatf <- cbind(udf, tum_l[[i]][rownames(udf),,drop=FALSE])
    colnames(pdatf)[ncol(pdatf)] <- "TumEnr"
    pdatf[,"cell_type_main"] <- factor(pdatf[,"cell_type_main"])

    de_stat_i <- de_stat_l[[i]]
    wpadj <- de_stat_i["Prol", "w_p_adj"]
    astk <- get_astk(wpadj)

    # umap of enr scores (main)
    o <- order(pdatf[,"TumEnr"], decreasing=FALSE)
    p <- ggplot(pdatf[o,], aes(x=UMAP1, y=UMAP2, color=TumEnr)) + 
    geom_point(size = 0.01, shape = 16) + 
    geom_text(label = astk, x = uxmax, y = uymin, 
              color = "black", size = 6, hjust = 1, vjust = 0) + 
    labs(title = utitles[i]) + 
    ggtitle(utitles[i]) + 
    scale_color_gradientn(colours = reds, 
                          name = gsub(" ", " ", ctitles[i])) + 
    th_um 

    gr_us <- ggplotGrob(p)
    gr_us <- raster_ggpoints(gr_us, w = 3, h = 3)

    # get coordinates for asterisks on Prol
    p_x <- which(levels(pdatf[,"cell_type_main"]) == "Prol")
    yvals <- pdatf[pdatf[,"cell_type_main"] == "Prol", "TumEnr"]
    yquant <- quantile(yvals, prob = c(.25, .75))
    whisk_max <- yquant[2] + (1.5 * (yquant[2] - yquant[1]))
    whisk_min <- yquant[1] - (1.5 * (yquant[2] - yquant[1]))
    whisk_r <- whisk_max - whisk_min
    p_y <- whisk_max + (whisk_r * .1)

    # box plot of enr scores (main)
    p <- ggplot(pdatf, aes(x = cell_type_main, y = TumEnr)) + 
    geom_boxplot(outlier.shape=16, outlier.size=0.1, outlier.alpha=0.1) + 
    geom_text(label = astk, x = p_x, y = p_y, color = "black", 
              size = 4, hjust = 0.5, vjust = 0, parse = FALSE) + 
    scale_y_continuous(expand = expansion(mult = c(0.1, .25))) + 
    labs(x = NULL, y = btitles[i], title = NULL) + 
    ggtitle(NULL) + 
    th_vbox

    gr_bs <- ggplotGrob(p)

    # stack panels
    # xlab_r <- gr_us$layout[gr_us$layout$name == "xlab-b", 'b']
    # gr_us_s <- gr_us[1:xlab_r,]
    # gr_us_s <- gtable_add_rows(gr_us_s, heights = unit(6, "points"), nrow(gr_us_s))

    # pan_r <- gr_bs$layout[gr_bs$layout$name == "panel", 't']
    # gr_bs_s <- gr_bs[pan_r:nrow(gr_bs),]

    # combine
    gr_l <- list(gr_us, gr_bs)
    hs <- unit(c(0.7, 0.3), "npc")
    gtmp <- stack_gtable_v(gr_l, heights = hs)

    # add tags
    gtmp$grobs[[1]] <- gtable_add_tag(gtmp$grobs[[1]], tags[[i]][1], fs=16, just=c(0,1), l=tagl)
    gtmp$grobs[[2]] <- gtable_add_tag(gtmp$grobs[[2]], tags[[i]][2], fs=16, just=c(0,-0.2), l=tagl)

    gtmp <- pad_plot(gtmp, t=mar_pct, r=mar_pct, b=mar_pct, l=mar_pct)

    gr3_l[[i]] <- gtmp
}

#=========================================
# put together a gtable for plotting everything
#=========================================

ws <- c(3, 3)
hs <- c(3, 2.5, 4.5)
punits <- "inches"

gt <- gtable(widths = unit(ws, punits), heights = unit(hs, punits))
gt <- gtable_add_grob(gt, gr1_l[[1]], t=1, b=1, l=1, r=1, name="a")
gt <- gtable_add_grob(gt, gr1_l[[2]], t=1, b=1, l=2, r=2, name="b")
gt <- gtable_add_grob(gt, gr2_l[[1]], t=2, b=2, l=1, r=1, name="c")
gt <- gtable_add_grob(gt, gr2_l[[2]], t=2, b=2, l=2, r=2, name="d")
gt <- gtable_add_grob(gt, gr3_l[[1]], t=3, b=3, l=1, r=1, name="e")
gt <- gtable_add_grob(gt, gr3_l[[2]], t=3, b=3, l=2, r=2, name="f")
#gt <- gtable_add_grob(gt, gr1b, t=2, b=2, l=1, r=1, name="b")

pdf(paste0(dir_plt, "fig2.pdf"), width = 6, height = 10)
grid.draw(gt)
dev.off()


