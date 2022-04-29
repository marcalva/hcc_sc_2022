
# plot fine cell-type proportion correlation

setwd("../../")

library(ggplot2)
library(reshape2)
source("scripts/sc_func.R")
# source("scripts/reference_free.R")
# library(Bisque)

#========================================================
# functions
#========================================================

scale_clip <- function(datf, z = 3){
    y <- apply(datf, 1, function(x) (x - mean(x))/sd(x))
    y <- t(y)
    colnames(y) <- colnames(datf)
    y[y < -z] <- -z
    y[y > z] <- z
    return(y)
}

#========================================================
#========================================================

cl_id <- "cell_type_fine"

# read in data
out_dir <- "data/processed/tcga_hcc/expr/"

fn <- paste0(out_dir, "tcga.gencodev26.rds")
gencode <- readRDS(fn)

# colors
fn <- "data/ref/colors/rd_bu_div.csv"
rd_bu <- read.table(fn)
rd_bu <- rev(rd_bu[,1])

# proportions
fn <- paste0("data/processed/tcga_hcc/ctp/tcga.TMM.", cl_id, ".decomp.rds")
tmm.ct.md <- readRDS(fn)
tmm.ct.mdp <- tmm.ct.md$bulk.props

for (ct in names(tmm.ct.md$genes.used)){
    tmm.ct.md$genes.used[[ct]] <- gencode[tmm.ct.md$genes.used[[ct]], "Name"]
}

# Plotting directory
dir_plot <- paste0("exp/tcga_hcc/ctp.", cl_id, "/plots/")
dir.create(dir_plot, showWarnings = FALSE, recursive = TRUE)

#========================================================
# plot co-abundance
#========================================================

th_hm <- theme_classic() +
theme(text = element_text(size = 12),
      plot.title = element_text(size = 14, hjust = 0.5), 
      axis.text.x = element_text(size = 10, color = "black", angle = 90, hjust=1), 
      axis.text.y = element_text(size = 10, color = "black"), 
      axis.ticks = element_blank(), 
      legend.direction = "vertical", 
      legend.title = element_text(size = 12), 
      legend.key.size = unit(12, "points"))

cell_types <- rownames(tmm.ct.mdp)
cors <- cor(t(tmm.ct.mdp))

cors_m <- melt(cors)
colnames(cors_m) <- c("x1", "x2", "r")
cors_m[,"x1"] <- factor(cors_m[,"x1"], levels = cell_types)
cors_m[,"x2"] <- factor(cors_m[,"x2"], levels = cell_types)

p <- ggplot(cors_m, aes(x = x1, y = x2, fill = r)) + 
    geom_tile() + 
    scale_fill_gradientn(colours = rd_bu, limits = c(-1, 1), 
                         name = bquote(italic(R))) + 
    labs(x = "Sub cell-type", y = "Sub cell-type") + 
    ggtitle("Correlation of sub cell-type proportions") + 
    theme_bw() + 
    th_hm

outfn <- paste0(dir_plot, "co_abund.pdf")
ggsave(outfn, width = 7, height = 6)


