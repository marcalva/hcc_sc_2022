
setwd("../../")

library(Seurat)
library(ggplot2)
library(scales)
library(gridExtra)
source("scripts/plot.R")

#=========================================
# Functions
#=========================================

create_dir <- function(p){
    dir.create(p, showWarnings=FALSE, recursive=TRUE)
}

#=========================================
# Plots
#=========================================

dir_plot <- "exp/d_sct_cca_all/merge/plots/"; create_dir(dir_plot)

infn <- "data/processed/d_sct_cca_all/merge/seur.cca.rds"
seur <- readRDS(infn)

# replace NASH with Non-tumor in sample ID
seur$SampleID <- sub("NASH", "Normal", seur$SampleID)
seur$SampleID <- sub("Normal", "Non-tumor", seur$SampleID)

# replace NonTumor with Non-tumor in tumor status
seur$Tumor <- as.character(seur$Tumor)
seur@meta.data[seur$Tumor == "NonTumor", "Tumor"] <- "Non-tumor"

# order factors from tumor to non-tumor
datf_s <- seur@meta.data[,c("SampleID", "Patient", "Tumor")]
datf_s[,"SampleID"] <- as.character(datf_s[,"SampleID"])
datf_s <- unique(datf_s)
o <- order(datf_s$Tumor, datf_s$Patient)

lv <- datf_s[o,"SampleID"]

seur$SampleID <- factor(seur$SampleID, levels = lv)

# cell type colors
ct_cols <- read.table("data/ref/TumorColors/CellType_colors.1.csv", 
                      header=FALSE, row.names=1, sep=",", stringsAsFactors=FALSE)
ct_cols_v <- ct_cols[,1]
names(ct_cols_v) <- rownames(ct_cols)

# patient colors
pt_cols <- read.table("data/ref/TumorColors/patient_colors.1.txt", 
                      header=FALSE, row.names=1, sep=",", stringsAsFactors=FALSE)
pt_cols_v <- pt_cols[,1]
names(pt_cols_v) <- rownames(pt_cols)

# tumor colors
fn <- "data/ref/TumorColors/tumor_colors.1.txt"
tum_cols <- read.table(fn, header=FALSE, sep=",", stringsAsFactors=FALSE)
tum_colsv <- tum_cols[,2]
names(tum_colsv) <- tum_cols[,1]
names(tum_colsv)[names(tum_colsv) == "Tumor"] <- "Tumor"
names(tum_colsv)[names(tum_colsv) == "NonTumor"] <- "Non-tumor"

# sample colors
sm_cols <- read.table("data/ref/TumorColors/sample_colors.1.txt", 
                      header=FALSE, row.names=1, sep=",", stringsAsFactors=FALSE)
sm_cols_v <- sm_cols[,1]
names(sm_cols_v) <- rownames(sm_cols)
names(sm_cols_v) <- sub("Normal", "Non-tumor", names(sm_cols_v))

#=========================================
# Tumor status UMAP
#=========================================

pw <- 6; ph <- 6.1

p <- plot_umap_labels(seur, 
                 colname = "Tumor", 
                 legend_title = "Status", 
                 colors = tum_colsv, 
                 label = FALSE, 
                 size = 0.3,
                 ax_title_size=12, 
                 alpha = 0.5) + 
    theme(legend.key.width = unit(0, "inches"), 
          legend.position = "top", 
          legend.title = element_blank(),
          legend.margin = margin(0, 0, 0, 0), 
          legend.text = element_text(size = 12), 
          legend.text.align = 0)
ggsave(paste0(dir_plot, "UMAP.tumor.sfig.jpeg"), 
       width = pw,  
       height = ph, 
       units = "in")
ggsave(paste0(dir_plot, "UMAP.tumor.sfig.pdf"), 
       width = pw, 
       height = ph, 
       units = "in")


#=========================================
# Tumor status cell type frequency
#=========================================

TumorCT <- table(seur$CellType, seur$Tumor)
TumorCT <- sweep(TumorCT, 2, colSums(TumorCT), "/")
# TumorCT <- sweep(TumorCT, 1, rowSums(TumorCT), "/")

TumorCTm <- reshape2::melt(TumorCT)
colnames(TumorCTm) <- c("CellType", "Tumor", "Prop")

source("scripts/ggplot_wrappers/R/barplots.R")
p <- pct_stack_barplt(TumorCTm, x = "CellType", y = "Prop", fill = "Tumor")
p <- p + labs(x = "Cell-type", y = "Percent")
p <- p + theme(panel.grid = element_blank(), 
               axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0.5), 
               legend.position = "none", 
               legend.text = element_text(size = 8), 
               legend.margin = margin(0, 0, 0, 0), 
               legend.key.width = unit(0, "inches"))
p <- p + scale_fill_manual(values = tum_colsv)
outfn <- paste0(dir_plot, "Tumor.CellType.prop.sfig.pdf")
ggsave(outfn, width = 6, height = 2)


#=========================================
#=========================================

TumorCT <- table(seur$Major, seur$Tumor)
# TumorCT <- sweep(TumorCT, 1, rowSums(TumorCT), "/")
TumorCT <- sweep(TumorCT, 2, colSums(TumorCT), "/")

TumorCTm <- reshape2::melt(TumorCT)
colnames(TumorCTm) <- c("Major", "Tumor", "Prop")
# TumorCTm[,"Major"] <- factor(TumorCTm[,"Major"], levels = TumorCTo)

source("scripts/ggplot_wrappers/R/barplots.R")
p <- pct_stack_barplt(TumorCTm, x = "Major", y = "Prop", fill = "Tumor")
p <- p + coord_flip()
p <- p + theme(panel.grid = element_blank(), 
               axis.text.x = element_text(angle = 0, hjust = 0.5))
p <- p + scale_fill_manual(values = tumor_cols)
outfn <- paste0(dir_plot, "Tumor.Major.prop.pdf")
ggsave(outfn, width = 5, height = 5)


#=========================================
#=========================================

