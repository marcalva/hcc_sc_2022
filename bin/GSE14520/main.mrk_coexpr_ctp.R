
# plot marker co-expression and correlation with ctp estimates

setwd("../../")

library(ggplot2)
library(reshape2)
library(grid)
library(gtable)
source("scripts/grid_hm.R")

#===============================================================================
# read in data
#===============================================================================

fn <- "data/processed/tcga_hcc/expr/tcga.gencodev26.rds"
gencode <- readRDS(fn)

# colors
fn <- "data/ref/colors/tumor_colrs.csv"
tum_cols <- read.table(fn, header=FALSE, sep=",")
tum_colsv <- tum_cols[,2]
names(tum_colsv) <- tum_cols[,1]
names(tum_colsv)[names(tum_colsv) == "Tumor"] <- "Tumor"
names(tum_colsv)[names(tum_colsv) == "NonTumor"] <- "Non-tumor"

fn <- "data/ref/colors/rd_bu_div.csv"
rd_bu <- read.table(fn)
rd_bu <- rev(rd_bu[,1])

# expression data
fn <- "data/processed/GSE14520/expr.RMA_log2.gmean.rds"
tmm <- readRDS(fn)

cl_id <- "cell_type_main"

# cell-type proportions
fn <- paste0("data/processed/GSE14520/ctp/", "LCI.", cl_id, ".decomp.rds")
tmm.ct.md <- readRDS(fn)

tmm.ct.mdp <- tmm.ct.md$bulk.props

# get list of bisque markers
dcmp_mrk <- do.call(c, tmm.ct.md$genes.used)
dcmp_mrk <- make.unique(dcmp_mrk)

# Plotting directory
dir_plot <- paste0("exp/GSE14520/ctp.", cl_id, "/plots/")
dir.create(dir_plot, showWarnings = FALSE, recursive = TRUE)

# generate color palette for cell types
ctu <- sort(rownames(tmm.ct.mdp))
nct <- length(ctu)
ct_pal <- hcl_pal(nct, chr = c(80,80), lum = c(60,80), offset = 0, rand=TRUE, seedn=321)
names(ct_pal) <- ctu

#===============================================================================
# create plot
#===============================================================================

source("scripts/grid_hm.R")
source("scripts/ggplot_raster.R")

# co-expression plot

ct_grp_l <- list()
tmm_mrk_l <- list()
for (n in names(tmm.ct.md$genes.used)){
    tmm_mrk_l[[n]] <- tmm[tmm.ct.md$genes.used[[n]],]
    ct_grp_l[[n]] <- rep(n, length(tmm.ct.md$genes.used[[n]]))
}

tmm_mrk <- do.call(rbind, tmm_mrk_l)
rownames(tmm_mrk) <- dcmp_mrk
ct_grp <- do.call(c, ct_grp_l)

tmm_mrk_cor <- cor(t(tmm_mrk))

tmm_mrk_corm <- melt(tmm_mrk_cor)

title_x_text <- list("label" = "Cell-type marker", "rot" = 0, 
                     "gp" = gpar(fontsize = 10, lineheight = 0.8))
title_y_text <- list("label" = "Cell-type marker", "rot" = 90, 
                     "gp" = gpar(fontsize = 10, lineheight = 0.8))
axis_x_text <- list("gp" = gpar(fontsize = 8), "rot" = -90, 
                    "just" = c(1, 0.5))
axis_y_text <- list("gp" = gpar(fontsize = 8), "rot" = 0, 
                    "just" = c(1, 0.5))
title_leg_text <- list("label" = bquote(italic(R)), 
                       "gp" = gpar(fontsize = 10))
lab_leg_text <- list("gp" = gpar(fontsize = 8))
key_dims <- c(12,120)

hm_gr1 <- gt_hm(tmm_mrk_corm, 
                groups_x = ct_grp, groups_y = ct_grp,
                group_cols = ct_pal, 
                val_cols = rd_bu, 
                val_lims = c(-1,1), 
                key_dims = key_dims,
                title_x_text = title_x_text, 
                title_y_text = title_y_text,
                axis_x_text = axis_x_text,
                axis_y_text = axis_y_text,
                title_leg_text = title_leg_text,
                lab_leg_text = lab_leg_text)

# expression-proportion correlation
tmm_mrk_prop_cor <- cor(t(tmm.ct.mdp[,colnames(tmm_mrk)]), t(tmm_mrk))
ctct <- rownames(tmm_mrk_prop_cor)

tmm_mrk_prop_corm <- melt(tmm_mrk_prop_cor)

title_x_text <- list("label" = "Cell-type\nproportion", "rot" = 0, 
                     "gp" = gpar(fontsize = 10, lineheight = 0.8))
title_y_text <- list("label" = NULL, "rot" = 90, 
                     "gp" = gpar(fontsize = 10))
hm_gr2 <- gt_hm(tmm_mrk_prop_corm, 
                groups_x = ctct, 
                groups_y = ct_grp, 
                group_cols = ct_pal, 
                val_cols = rd_bu, 
                val_lims = c(-1,1), 
                key_dims = key_dims,
                title_x_text = title_x_text, 
                title_y_text = title_y_text,
                axis_x_text = axis_x_text,
                axis_y_text = axis_y_text,
                title_leg_text = title_leg_text,
                lab_leg_text = lab_leg_text)

# combine heatmaps
hm_gr1s <- hm_gr1[,1:5]
hm_gr2s <- hm_gr2[,4:7]

hm_gr_all <- hm_gr1s
hm_gr_all$widths[5] <- unit(0.8, "null")

hm_gr_all <- gtable_add_cols(hm_gr_all, 
                             width = unit(0.25, "inches"))

hm_gr_all <- gtable_add_cols(hm_gr_all, widths = unit(0.2, "null"))
for (i in 2:5){
    hm_gr_all <- gtable_add_grob(hm_gr_all, 
                                 hm_gr2[i,5]$grobs[[1]], 
                                 t=i, l=7)
}

hm_gr_all <- gtable_add_cols(hm_gr_all, 
                             width = unit(12, "points"))

hm_gr_all <- gtable_add_cols(hm_gr_all, 
                             widths = hm_gr2$widths[7])

hm_gr_all <- gtable_add_grob(hm_gr_all, 
                             hm_gr2[5,7]$grobs[[1]], 
                             t=5, l = 9)

hm_gr_all <- gtable_add_cols(hm_gr_all, 
                             width = unit(5.5, "points"))

hxlab_hght <- max(convertUnit(grobHeight(hm_gr_all$grobs[[6]]), unitTo="cm"), 
                  convertUnit(grobHeight(hm_gr_all$grobs[[8]]), unitTo="cm"))
hm_gr_all$heights[2] <- 1.4 * hxlab_hght

# add title
hm_gr_all <- hm_gr_all[2:length(hm_gr_all$heights),]
hm_gr_all <- gtable_add_rows(hm_gr_all, heights = unit(5.5, "points"), pos=0)
gtitle <- "Main cell-type marker gene co-expression\nand proportion correlation in LCI"
gtitle_gr <- textGrob(gtitle, gp = gpar(fontsize=12, lineheight=0.8))
th <- 1.4 * convertUnit(grobHeight(gtitle_gr), unitTo="cm")
hm_gr_all <- gtable_add_rows(hm_gr_all, heights = th, pos=0)
hm_gr_all <- gtable_add_grob(hm_gr_all, gtitle_gr, t=1, l=1, r=length(hm_gr_all$widths))
hm_gr_all <- gtable_add_rows(hm_gr_all, heights = unit(5.5, "points"), pos=0)

out_fn <- paste0(dir_plot, "mrk_ctp_coexpr.pdf")
pdf(out_fn, width = 6, height = 5)
grid.draw(hm_gr_all)
dev.off()

# save gtable
out_fn <- paste0(dir_plot, "mrk_ctp_coexpr.gtable.rds")
saveRDS(hm_gr_all, out_fn)

