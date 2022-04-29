
setwd("../../")

library(ggplot2)
library(grid)
library(gtable)
source("scripts/gtable_stack.R")
source("scripts/ggplot_raster.R")

fn <- "exp/tcga_hcc/ctp.cell_type_main/plots/mrk_ctp_coexpr.gtable.rds"
tcga <- readRDS(fn)

fn <- "exp/GSE14520/ctp.cell_type_main/plots/mrk_ctp_coexpr.gtable.rds"
lci <- readRDS(fn)

# rasterize
res <- 300
w <- 6; h <- 6
tcga$grobs[[1]] <- rasterize_grob(tcga$grobs[[1]], res=res, w=w, h=h)
tcga$grobs[[11]] <- rasterize_grob(tcga$grobs[[11]], res=res, w=w, h=h)
lci$grobs[[1]] <- rasterize_grob(lci$grobs[[1]], res=res, w=w, h=h)
lci$grobs[[11]] <- rasterize_grob(lci$grobs[[11]], res=res, w=w, h=h)

tcga <- gtable_add_tag(tcga, "a", fs = 16, just = c(0, 1), t=2, l=2)
lci <- gtable_add_tag(lci, "b", fs=16, just=c(0, 1), t=2, l=2)

gt <- gtable(widths = unit(1, "npc"), heights = unit(c(0.5, 0.5), "npc"))

gt <- gtable_add_grob(gt, tcga, t=1, l=1)
gt <- gtable_add_grob(gt, lci, t=2, l=1)

dir_man <- "exp/manuscript/"
dir.create(dir_man, showWarnings=FALSE, recursive=TRUE)

out_fn <- paste0(dir_man, "FigureS7.pdf")
pdf(out_fn, width = 6, height = 10)
grid.draw(gt)
dev.off()

