
# create a dotplot of marker gene expression
# size is log2 fold change, fill is percent droplets expressing
# for subcell-types

setwd("../../")

library(ggplot2)

#===============================================================================
#===============================================================================

# colors
reds <- read.csv("data/ref/colors/red_colrs.csv", header=FALSE)
reds <- reds[,1]

dir_exp <- "exp/sharma_aiz/markers/"
dir.create(dir_exp, recursive = TRUE, showWarnings = FALSE)

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
      legend.margin = margin(2, 2, 2, 2, unit = "pt"), 
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

    p <- ggplot(mrk, aes(y = Name, x = cluster, size = avg_log2FC, color = pct.1)) + 
    geom_point() + 
    scale_color_gradientn(colours = reds,
                          name = color_name) + 
    scale_size_continuous(range = c(0.1, 2.5), 
                          name = size_name) + 
    labs(x = xlabs[cl_id], y = "Cell-type markers",
         title = titles[cl_id]) + 
    guides(size = size_guide) + 
    th_dot

    out_fn <- paste0(dir_exp, "markers.", cl_id, ".dotplot.pdf")
    ggsave(out_fn, width = widths[cl_id], height = 7)
}

