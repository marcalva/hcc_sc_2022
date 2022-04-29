
# Figure S8
# survival effects of marker genes

setwd("../../")

library(ggplot2)
library(grid)
library(gtable)
library(scales)
library(ggfortify)
library(reshape2)
source("scripts/surv_df.R")
source("scripts/color_pal.R")
source("scripts/ggplot_raster.R")
source("scripts/gtable_stack.R")

create_dir <- function(p){
    dir.create(p, showWarnings=FALSE, recursive=TRUE)
}

mar_pct <- 0.025 # pct margin around panel
tagl <- 3

#========================================================
#========================================================

# Set directories
dir_plt <- "exp/manuscript/"
dir.create(dir_plt, showWarnings=FALSE, recursive=TRUE)

#================================================
# marker stats
#================================================

fn <- "exp/tcga_hcc/survival.gene/cell_type_main/cell_type_main.surv.mrk_stat.txt"
tcga_mrk <- read.table(fn, header = TRUE)

fn <- "exp/GSE14520/survival.gene/cell_type_main/cell_type_main.surv.mrk_stat.txt"
lci_mrk <- read.table(fn, header = TRUE)

# rename events in lci data frame
evmap <- c("Survival" = "OS", "Recurr" = "RCR")
lci_mrk[,"event"] <- evmap[lci_mrk[,"event"]]

#================================================
# plot pct GWS markers
#================================================

th_bar <- theme_classic() + 
theme(panel.grid.major.y = element_line(), 
      plot.title = element_text(size = 12, hjust = 0.5), 
      axis.text.x = element_text(size = 8, color="black", angle=90, hjust=1, vjust=0.5), 
      axis.text.y = element_text(size = 8, color="black"), 
      legend.key.size = unit(12, "points"), 
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8))

# list data
events <- c("OS", "PFI", "OS")
cohorts <- c("TCGA", "TCGA", "LCI")
mrk_stat_l <- list(tcga_mrk, tcga_mrk, lci_mrk)
tags <- c("a", "b", "c")

gr_l <- list()

for (i in 1:3){
    eve <- events[i]
    cohort <- cohorts[i]
    pdatf <- mrk_stat_l[[i]]
    pdatf <- pdatf[pdatf$event == eve,]

    pdatf[,"pct_neg"] <- -1 * pdatf[,"pct_neg"]
    pdatf[,"pct_gws_neg"] <- -1 * pdatf[,"pct_gws_neg"]
    pdatfs <- pdatf[,-which(colnames(pdatf) == "pct_gws")]
    pdatfs <- pdatfs[,-which(colnames(pdatfs) == "ct_n")]

    pdatfsm <- melt(pdatfs, id.vars = c("CellType", "event"))

    pdatfsm[,"Direction"] <- "HR > 1"
    pdatfsm[grep("neg", pdatfsm[,"variable"]), "Direction"] <- "HR < 1"
    pdatfsm[,"GWS"] <- "Adj p >= 0.05"
    pdatfsm[grep("gws", pdatfsm[,"variable"]), "GWS"] <- "Adj p < 0.05"

    type_colr <- hcl_pal(3, offset = 260, chr = c(140, 140), lum = c(60,60), 
                         rand=FALSE)[1:2]
    names(type_colr) <- c("HR > 1", "HR < 1")
    type_alph <- c("Adj p < 0.05" = 1, "Adj p >= 0.05" = 0.2)

    coltitle <- paste0(cohort, "\nHR direction\nfor ", eve)
    alptitle <- NULL
    ytitle <- "Proportion of\ncell-type markers"
    xtitle <- "Main cell-types"
    ptitle <- paste0("Directional effect of \nmarker genes on ", eve, 
                     " in ", cohort)

    p <- ggplot(pdatfsm, aes(x = CellType, y = value)) +
    geom_bar(aes(fill = Direction, alpha = GWS), position = "identity", stat = "identity") + 
    scale_fill_manual(values = type_colr, name = coltitle) + 
    scale_alpha_manual(values = type_alph, name = alptitle) + 
    scale_y_continuous(labels = scales::percent, limits = c(-1,1)) + 
    labs(y = ytitle, x = xtitle, title = ptitle) + 
    th_bar

    gr_i <- ggplotGrob(p)
    gr_i <- gtable_add_tag(gr_i, tags[i], fs=16, just=c(0,1), l=tagl)
    gr_i <- pad_plot(gr_i, t=mar_pct, r=mar_pct, b=mar_pct, l=mar_pct)

    gr_l[[i]] <- gr_i
}

#=========================================
# put together a gtable for plotting everything
#=========================================

w <- 5
h <- 3
punit <- "inches"

ws <- c(w)
hs <- c(h,h,h)

gt <- gtable(widths = unit(ws, punit), 
             heights = unit(hs, punit))

gt <- gtable_add_grob(gt, gr_l[[1]], t=1,b=1,l=1,r=1)
gt <- gtable_add_grob(gt, gr_l[[2]], t=2,b=2,l=1,r=1)
gt <- gtable_add_grob(gt, gr_l[[3]], t=3,b=3,l=1,r=1)

pdf(paste0(dir_plt, "FigureS8.pdf"), width = sum(ws), height = sum(hs))
grid.draw(gt)
dev.off()

