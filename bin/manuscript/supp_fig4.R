
# create a dotplot of marker gene expression
# size is log2 fold change, fill is percent droplets expressing
# for subcell-types

setwd("../../")

library(ggplot2)
library(gtable)
library(grid)
source("scripts/gtable_stack.R")

#===============================================================================
#===============================================================================

# colors
reds <- read.csv("data/ref/colors/red_colrs.csv", header=FALSE)
reds <- reds[,1]

dir_exp <- "exp/sharma_aiz/markers/"
dir.create(dir_exp, recursive = TRUE, showWarnings = FALSE)

dir_man <- "exp/manuscript/"

# ggplot theme
th_dot <- theme_classic() + 
theme(plot.title = element_text(hjust = 0.5), 
      panel.border = element_rect(color = "black", fill = NA),
      panel.grid.minor = element_line(), 
      panel.grid.major = element_line(), 
      text = element_text(size = 8), 
      axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5, color = "black"), 
      axis.text.y = element_text(color = "black"), 
      legend.position = "bottom", 
      legend.direction = "horizontal", 
      legend.box = "vertical", 
      legend.box.spacing = unit(0, "points"), 
      legend.margin = margin(2, 12, 2, 12, unit = "pt"), 
      legend.key.height = unit(10, "points"), 
      legend.key.width = unit(14, "points"), 
      legend.spacing.y = unit(2, "points"),
      legend.spacing.x = unit(2, "points"),
      legend.title = element_text(hjust=1, vjust=1.0))

cl_ids <- c("cell_type_main", "cell_type_fine")
n_genes <- c(8, 3); names(n_genes) <- cl_ids
widths <- c(3, 4); names(widths) <- cl_ids

xlabs <- c("Main cell-type", "Subcell-type"); names(xlabs) <- cl_ids
titles <- c("Main cell-type marker expression", 
            "Subcell-type marker expression")
names(titles) <- cl_ids

tags <- c("a", "b"); names(tags) <- cl_ids

gt_leg <- c()
gt_l <- list()

for (cl_id in cl_ids){
    in_fn <- paste0("exp/sharma_aiz/markers/markers.", cl_id, ".txt")
    mrk <- read.table(in_fn, header = TRUE, quote = '\t')

    genes <- tapply(mrk$Name, mrk$cluster, head, n_genes[cl_id])
    genes <- unique(unlist(genes))
    k <- mrk$Name %in% genes
    mrk <- mrk[k,]
    mrk[,"Name"] <- factor(mrk[,"Name"], levels = rev(genes))
    ctu <- sort(unique(mrk$cluster))
    mrk[,"cluster"] <- factor(mrk[,"cluster"], levels = ctu)

    size_guide <- guide_legend(label.position = "bottom")
    color_name = "Percent of droplets\nexpressing gene"
    size_name = "Average log\nfold change"

    p <- ggplot(mrk, aes(x = cluster, y = Name, size = avg_log2FC, color = pct.1)) + 
    geom_point() + 
    scale_color_gradientn(colours = reds,
                          name = color_name, 
                          limits = c(0.1,1)) + 
    scale_size_continuous(range = c(0.1, 2.5), 
                          name = size_name, 
                          limits = c(0,8)) + 
    labs(x = xlabs[cl_id], y = "Cell-type markers",
         title = titles[cl_id]) + 
    guides(size = size_guide) + 
    th_dot

    gt_i <- ggplotGrob(p)

    leg_ix <- which(gt_i$layout$name == "guide-box")
    # grab legend
    if (cl_id == cl_ids[1]){
        gt_leg <- gt_i$grobs[[leg_ix]]
    }
    leg_t <- gt_i$layout[leg_ix, "t"]
    gt_i <- gt_i[1:(leg_t-1),]
    
    gt_i <- gtable_add_tag(gt_i, tags[cl_id], fs=16, ff = "bold", l = 3, 
                           just = c(0,1))
    gt_i <- pad_plot(gt_i)

    gt_l[[cl_id]] <- gt_i

}

#=========================================
#=========================================

gt <- gtable(widths = unit(c(3,4), "inches"), 
             heights = unit(c(7,1), "inches"))
gt <- gtable_add_grob(gt, gt_l[["cell_type_main"]], t=1,b=1,l=1,r=1)
gt <- gtable_add_grob(gt, gt_l[["cell_type_fine"]], t=1,b=1,l=2,r=2)
gt <- gtable_add_grob(gt, gt_leg, t=2,l=1,r=2)

out_fn <- paste0(dir_man, "FigureS4.pdf")
pdf(out_fn, width = 7, height = 8)
grid.draw(gt)
dev.off()

