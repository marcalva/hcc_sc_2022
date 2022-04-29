
# Figure 3

setwd("../../")

library(Seurat)
library(ggplot2)
library(grid)
library(gtable)
library(ggrepel)
library(scales)
library(gridExtra)
library(reshape2)
library(ggfortify)
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

mar_pct <- 0.025 # pct margin around panel
tagl <- 3

#=========================================
# data
#=========================================

fn <- "data/processed/sharma_aiz/liver.int_rand.md_umap.txt.gz"
udf <- read.table(fn, header=TRUE, row.names=1, sep='\t')
colnames(udf)[colnames(udf) == "UMAP_1"] <- "UMAP1"
colnames(udf)[colnames(udf) == "UMAP_2"] <- "UMAP2"
set.seed(1, kind = 'Mersenne-Twister')
udf <- udf[sample(1:nrow(udf)),]


dir_plt <- "exp/manuscript/"
dir.create(dir_plt, showWarnings=FALSE, recursive=TRUE)

reds <- read.csv("data/ref/colors/red_colrs.csv", header=FALSE)
reds <- reds[,1]

#=========================================
# survival curves for Prol in TCGA
#=========================================

th_curve <- theme_bw() + 
theme(legend.margin = margin(0,0,0,0), 
      plot.title = element_text(hjust = 0.5, size = 8),
      text = element_text(size=8), 
      legend.text = element_text(size = 8), 
      legend.title = element_text(hjust=0.5, size = 8), 
      legend.box.spacing = unit(2, "points"), 
      legend.key = element_rect(fill=NA), 
      legend.key.size = unit(12, "points"), 
      axis.text.x = element_text(color = "black"), 
      axis.text.y = element_text(color = "black"))

# read in KM results
fn <- "exp/tcga_hcc/survival/cell_type_main/km.events.ct.l.rds"
tcga.km.l <- readRDS(fn)

fn <- "exp/tcga_hcc/survival/cell_type_main.cox_ph.formatted.txt"
tcga.surv_stats <- read.table(fn, header=TRUE, sep='\t')

fn <- "exp/GSE14520/survival/cell_type_main/km.events.ct.l.rds"
lci.km.l <- readRDS(fn)

fn <- "exp/GSE14520/survival/cell_type_main.GSE14520.cox_ph.formatted.txt"
lci.surv_stats <- read.table(fn, header=TRUE, sep='\t')

cohorts <- c("TCGA", "TCGA", "LCI")
events <- c("cdr.OS", "cdr.PFI", "OS")
pts <- c("OS", "PFI", "OS")
ylabs <- c("Percent survival", "Percent progression free", "Percent survival")
km_l <- list(tcga.km.l, tcga.km.l, lci.km.l)
stats_l <- list(tcga.surv_stats, tcga.surv_stats, lci.surv_stats)
time_divs <- c(365, 365, 12)
tags <- c("a", "b", "c")

hl_cols <- c("High" = "firebrick1", "Low" = "dodgerblue4")

gt1_l <- list()
for (i in 1:3){
    event <- events[i]
    pt <- pts[i]
    ct <- "Prol"

    surv.md <- km_l[[i]][[event]][[ct]]
    surv.md$time <- surv.md$time / time_divs[i]
    yticks <- paste0(seq(0, 100, 20), "%")

    sdf <- stats_l[[i]]
    k <- sdf[,"Event"] == pts[i] & sdf[,"Model"] == "Median" & sdf[,"Cell.type"] == "Prol"
    hr_txt <- sdf[k,"HR"]
    hr_txt <- paste0("HR == ", hr_txt)

    k <- sdf[,"Event"] == pts[i] & sdf[,"Model"] == "Median" & sdf[,"Cell.type"] == "Prol"
    pval_txt <- frmt_sct(sdf[k,"p_adj"])
    pval_txt <- paste0("adj~p == ", pval_txt)
    hr_pval_txt <- paste0(hr_txt, "\n", pval_txt)

    gtitle <- paste0(pt, " by ", ct, " frequency\n", cohorts[i])

    p_i <- autoplot(surv.md, censor.size = 2) +
    ggtitle(gtitle) + 
    xlab("Follow-up in years") + ylab(ylabs[i]) + 
    scale_x_continuous(breaks = seq(0,10,2)) + 
    scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2), labels=yticks) + 
    scale_color_manual(values = hl_cols) + 
    scale_fill_manual(values = hl_cols) + 
    geom_text(x = 0, y = .1, label = hr_txt, parse = TRUE, size = 2.5, hjust=0, vjust=0) +
    geom_text(x = 0, y = 0, label = pval_txt, parse = TRUE, size = 2.5, hjust=0, vjust=0) +
    theme_bw() + 
    guides(color = "none", fill = guide_legend(title=paste0(ct, "\nfreq. by\nmedian"))) + 
    th_curve

    if (i < 3)
        p_i <- p_i + theme(legend.position="none")

    gr_i <- ggplotGrob(p_i)

    gr_i <- gtable_add_tag(gr_i, tags[i], fs=16, just=c(0,1), l=tagl)
    gr_i <- pad_plot(gr_i, t=mar_pct, r=mar_pct, b=mar_pct, l=mar_pct)

    gt1_l[[i]] <- gr_i
}


gr1 <- gtable(widths = unit(c(.31, .31, .38), "npc"), heights = unit(1, "npc"))
gr1 <- gtable_add_grob(gr1, gt1_l[[1]], t=1,b=1,l=1,r=1,name='a')
gr1 <- gtable_add_grob(gr1, gt1_l[[2]], t=1,b=1,l=2,r=2,name='b')
gr1 <- gtable_add_grob(gr1, gt1_l[[3]], t=1,b=1,l=3,r=3,name='c')

#================================================
# jitter on top of box of log2HR for marker genes
# per cell-type
#================================================

theme_txt <- theme(text = element_text(colour="black", size=8), 
                   axis.text = element_text(colour="black", size=8), 
                   axis.ticks.y = element_blank(), 
                   plot.title = element_text(hjust=0.5))

theme_leg <- theme(legend.background = element_rect(color="black", 
                                                    linetype="solid", 
                                                    fill=NA),
                   legend.position = "bottom", 
                   legend.key.size = unit(1, "points"), 
                   legend.margin = margin(0,2.5,2.5,2.5))
#                    legend.box.margin = margin(0,0,0,0))


boxj_plot <- function(datf, x = "log2HR", y = "cluster", color = "gws", 
                      title=NULL, xlab = bquote("log"[2]*"HR"), ylab = NULL, 
                      vline = 0){
    # order y
    # order y
    y_m <- tapply(datf[,x], datf[,y], mean)
    o <- order(y_m, decreasing=FALSE)
    lv <- names(y_m)[o]
    datf[,y] <- factor(datf[,y], levels = lv)

    p <- ggplot(datf, aes_string(x = x, y = y)) + 
    geom_vline(xintercept = vline, color = "black", linetype = "dotted") + 
    geom_boxplot(fill = NA, outlier.shape=NA) + 
    geom_jitter(aes(color = color), height = 0.2, shape=16, size=0.5) + 
    ylab(ylab) + 
    xlab(xlab) + 
    ggtitle(title) + 
    scale_color_manual(values = gmap, breaks = "adj p < 0.05", 
                       name=NULL) + 
    theme_classic() + theme_txt + theme_leg + 
    theme(legend.position=c(1,.2), 
      legend.justification=c(1,0))
    return(p)
}

#=========================================
# plot HR for cell type markers
#=========================================

mrk_logfc_thr <- 0.5

fn <- "exp/tcga_hcc/survival.gene/cox_ph.ct_mrk.rds"
tcga_mrk_surv <- readRDS(fn) # TCGA

fn <- "exp/GSE14520/survival.gene/cox_ph.ct_mrk.rds"
lci_mrk_surv <- readRDS(fn) # LCI
# ctrm <- c("B")
# for (i in 1:length(lci_mrk_surv)){
#     k <- ! lci_mrk_surv[[i]][,"cluster"] %in% ctrm
#     lci_mrk_surv[[i]] <- lci_mrk_surv[[i]][k,]
# }


cohorts <- c("TCGA", "TCGA", "LCI")
events <- c("OS", "PFI", "Survival")
pts <- c("OS", "PFI", "OS")
mrk_surv_l <- list(tcga_mrk_surv, tcga_mrk_surv, lci_mrk_surv)
time_divs <- c(365, 365, 12)
tags <- c("d", "e", "f")
cl_id <- "cell_type_main"

gt2_l <- list()
for (i in 1:3){
    e <- events[i]

    pdatf <- mrk_surv_l[[i]][[cl_id]]
    pdatf <- pdatf[pdatf[,"surv_type"] == e,]
    
    pdatf[,"log2HR"] <- log2(pdatf[,"HR"])
    pdatf[,"gws"] <- " "
    pdatf[pdatf[,"HR_p_adj"] < 0.05,"gws"] <- "adj p < 0.05"
    gmap <- c(" " = "black", "adj p < 0.05" = "red")

    pdatf <- pdatf[pdatf[,"avg_log2FC"] > mrk_logfc_thr,]
    
    xlab_j <- paste0(cohorts[i], " bulk gene HR for ", pts[i])
    ylab_j <- "Cell-type marker genes"
    p_i <- boxj_plot(pdatf, x = "HR", 
                     color = NULL, 
                     title = paste0(pts[i], " by cell-type markers\nin bulk ", cohorts[i]), 
                     xlab = xlab_j, 
                     ylab = ylab_j, 
                     vline = 1) + 
    theme(legend.position="none")

    gr_i <- ggplotGrob(p_i)

    gr_i <- gtable_add_tag(gr_i, tags[i], fs=16, just=c(0,1), l=tagl)
    gr_i <- pad_plot(gr_i, t=mar_pct, r=mar_pct, b=mar_pct, l=mar_pct)

    gt2_l[[i]] <- gr_i
}

gr2 <- gtable(widths = unit(c(1/3,1/3,1/3), "npc"), heights = unit(1, "npc"))
gr2 <- gtable_add_grob(gr2, gt2_l[[1]], t=1,b=1,l=1,r=1,name='a')
gr2 <- gtable_add_grob(gr2, gt2_l[[2]], t=1,b=1,l=2,r=2,name='b')
gr2 <- gtable_add_grob(gr2, gt2_l[[3]], t=1,b=1,l=3,r=3,name='c')

#=========================================
# survival scores per droplets
#=========================================

th_umap <- theme_classic() + 
theme(text = element_text(size = 8),
      plot.title = element_text(hjust = 0.5, size = 8), 
      legend.text = element_text(size = 6, angle=30, hjust=1), 
      legend.title = element_text(size = 6, vjust = 1, hjust = 0.5), 
      legend.position = "bottom", 
      legend.key.height = unit(6, "points"),
      legend.key.width = unit(12, "points"),
      legend.box.spacing = unit(2, "points"), 
      axis.text=element_blank(),
      axis.ticks=element_blank())

th_box <- theme_classic() +
theme(text = element_text(size = 8),
      axis.text.x = element_text(color = "black", angle=90, hjust=1, vjust=0.5, size=8), 
      axis.text.y = element_text(color = "black"), 
      plot.title = element_text(hjust = 0.5))

# read in enrichment scores
fn <- "exp/tcga_hcc/sharma_aiz.surv_score/sharma_aiz.surv_enr_scores.txt"
tc_enr <- read.table(fn, header=TRUE, row.names=1, sep='\t')
colnames(tc_enr) <- gsub("1$", "", colnames(tc_enr))
tc_udf <- cbind(udf, tc_enr[rownames(udf),])
fn <- "exp/tcga_hcc/sharma_aiz.surv_score/cell_type_main.diff_stat.txt"
tc_stat <- read.table(fn, header=TRUE, row.names=1, sep='\t')


fn <- "exp/GSE14520/sharma_aiz.surv_score/sharma_aiz.surv_enr_scores.txt"
lc_enr <- read.table(fn, header=TRUE, row.names=1, sep='\t')
colnames(lc_enr) <- gsub("1$", "", colnames(lc_enr))
lc_udf <- cbind(udf, lc_enr[rownames(udf),,drop=FALSE])
fn <- "exp/GSE14520/sharma_aiz.surv_score/cell_type_main.diff_stat.txt"
lc_stat <- read.table(fn, header=TRUE, row.names=1, sep='\t')

uxmax <- max(udf[,"UMAP1"])
uymin <- min(udf[,"UMAP2"])

cohorts <- c("TCGA", "TCGA", "LCI")
udf_l <- list(tc_udf, lc_udf)
udf_ix <- c(1,1,2)
tags <- list(c("g", "j"), c("h", "k"), c("i", "l"))
events <- c("OS", "PFI", "OS")
utitles <- c("TCGA OS-decreasing\ngene scores", 
             "TCGA PFI-decreasing\ngene scores", 
             "LCI OS-decreasing\ngene scores")
ctitles <- c("TCGA OS-decreasing\ngene scores", 
             "TCGA PFI-decreasing\ngene scores", 
             "LCI OS-decreasing\ngene scores")
btitles <- c("TCGA\nOS-decreasing\ngene scores", 
             "TCGA\nPFI-decreasing\ngene scores", 
             "LCI\nOS-decreasing\ngene scores")
stat_l <- list(tc_stat, tc_stat, lc_stat)
yexpand <- c(0.25, 0.25, 0.1)

gr3_all <- list()

for (i in 1:3){
    event <- events[i]
    udf_i <- udf_l[[udf_ix[i]]]
    udf_i[,"cell_type_main"] <- factor(udf_i[,"cell_type_main"])

    stat_i <- stat_l[[i]]
    k <- stat_i[,"CellType"] == "Prol" & stat_i[,"event"] == event
    wpadj <- stat_i[k, "w_p_adj"]
    astk <- get_astk(wpadj)

    # umap of surv scores (main)
    o <- order(udf_i[,event], decreasing=FALSE)
    p_i1 <- ggplot(udf_i[o,], aes_string(x="UMAP1", y="UMAP2", color=event)) + 
    geom_point(size=0.01, shape=16) + 
    geom_text(label = astk, x = uxmax, y = uymin, 
              color = "black", size = 6, hjust = 1, vjust = 0) + 
    scale_color_gradientn(colours = reds, name = ctitles[i]) + 
    labs(title = utitles[i]) + 
    th_umap

    gr_i1 <- ggplotGrob(p_i1)
    gr_i1 <- raster_ggpoints(gr_i1)

    # get coordinates for asterisks on Prol
    p_x <- which(levels(udf_i[,"cell_type_main"]) == "Prol")
    yvals <- udf_i[udf_i[,"cell_type_main"] == "Prol", event]
    yquant <- quantile(yvals, prob = c(.25, .75))
    whisk_max <- yquant[2] + (1.5 * (yquant[2] - yquant[1]))
    whisk_min <- yquant[1] - (1.5 * (yquant[2] - yquant[1]))
    whisk_r <- whisk_max - whisk_min
    p_y <- whisk_max + (whisk_r * .1)

    # box plot of surv scores (main)
    p_i2 <- ggplot(udf_i, aes_string(x = "cell_type_main", y = event)) + 
    geom_boxplot(outlier.shape=16, outlier.size=0.1, outlier.alpha=0.1) + 
    geom_text(label = astk, x = p_x, y = p_y, color = "black", 
              size = 4, hjust = 0.5, vjust = 0, parse = FALSE) + 
    scale_y_continuous(expand = expansion(mult = c(0.1, yexpand[i]))) + 
    labs(y = btitles[i], x = NULL) + 
    th_box

    gr_i2 <- ggplotGrob(p_i2)

    # stack plots
    # xlab_r <- gr_i1$layout[gr_i1$layout$name == "xlab-b", 'b']
    # gr_i1 <- gr_i1[1:xlab_r,]
    # gr_i1 <- gtable_add_rows(gr_i1, heights = unit(12, "points"), nrow(gr_i1))

    pan_r <- gr_i2$layout[gr_i2$layout$name == "panel", 't']
    gr_i2 <- gr_i2[pan_r:nrow(gr_i2),]

    # combine
    gr_l <- list(gr_i1, gr_i2)
    hs <- unit(c(0.7, 0.3), "npc")
    gtmp <- stack_gtable_v(gr_l, heights = hs)

    # add tags
    gtmp$grobs[[1]] <- gtable_add_tag(gtmp$grobs[[1]], tags[[i]][1], fs=16, just=c(0,1), l=tagl)
    gtmp$grobs[[2]] <- gtable_add_tag(gtmp$grobs[[2]], tags[[i]][2], fs=16, just=c(0,-0.2), l=tagl)

    gtmp <- pad_plot(gtmp, t=mar_pct, r=mar_pct, b=mar_pct, l=mar_pct)

    gr3_all[[i]] <- gtmp
}

ws <- unit(c(1/3, 1/3, 1/3), "npc")
hs <- unit(1, "npc")
gt_scores <- gtable(widths=ws, heights=hs)
gt_scores <- gtable_add_grob(gt_scores, gr3_all[[1]], t=1,b=1,l=1,r=1)
gt_scores <- gtable_add_grob(gt_scores, gr3_all[[2]], t=1,b=1,l=2,r=2)
gt_scores <- gtable_add_grob(gt_scores, gr3_all[[3]], t=1,b=1,l=3,r=3)

#=========================================
# put together a gtable for plotting everything
#=========================================

ws <- c(7)
hs <- c(2.5, 2.5, 4)
punits <- "inches"

gt <- gtable(widths = unit(ws, punits), heights = unit(hs, punits))
gt <- gtable_add_grob(gt, gr1, t=1, b=1, l=1, r=1, name="abc")
gt <- gtable_add_grob(gt, gr2, t=2, b=2, l=1, r=1, name="abc")
gt <- gtable_add_grob(gt, gt_scores, t=3, b=3, l=1, r=1, name="abc")

pdf(paste0(dir_plt, "fig3.pdf"), width = 7, height = 9)
grid.draw(gt)
dev.off()

