
setwd("../../")

library(Seurat)

#=========================================
# Functions
#=========================================

create_dir <- function(p){
    dir.create(p, showWarnings=FALSE, recursive=TRUE)
}

#=========================================
# Set variables
#=========================================

# Gene data
gene_info <- read.table("data/ref/gencode26/gencode.v26.annotation.txt", 
                        header = TRUE, 
                        row.names = 1, 
                        stringsAsFactors = FALSE, 
                        sep = "\t")

# Set directories
dir_out <- "data/processed/aizarani2019/"; create_dir(dir_out)
# dir_plot <- "exp/d_sct_cca_all/cca_init/"; create_dir(dir_plot)

dir_exp <- "exp/aizarani2019/sct_clust/"; create_dir(dir_exp)

#=========================================
# Map cell types
#=========================================

# Read in Seurat object
seur_fn <- "data/processed/aizarani2019/seur.rds"
seur <- readRDS(seur_fn)

# res 0.5
ctmap <- c("0" = "Hep_1", 
           "1" = "Cholangiocyte", 
           "2" = "T_1", 
           "3" = "T_2", 
           "4" = "Kupffer_1", 
           "5" = "LSEC_1", 
           "6" = "LSEC_2", 
           "7" = "Hep_2", 
           "8" = "Hep_3", 
           "9" = "T_3", 
           "10" = "MacEndo", 
           "11" = "T_4", 
           "12" = "T_5", 
           "13" = "Hep_4", 
           "14" = "Kupffer_2", 
           "15" = "B_1", 
           "16" = "B_2", 
           "17" = "Hep_5", 
           "18" = "LSEC_3", 
           "19" = "T_prolif", 
           "20" = "Stellate")

mjmap <- c("0" = "Hep", 
           "1" = "Cholangiocyte", 
           "2" = "T", 
           "3" = "T", 
           "4" = "Kupffer", 
           "5" = "LSEC", 
           "6" = "LSEC", 
           "7" = "Hep", 
           "8" = "Hep", 
           "9" = "T", 
           "10" = "MacEndo", 
           "11" = "T", 
           "12" = "T", 
           "13" = "Hep", 
           "14" = "Kupffer", 
           "15" = "B", 
           "16" = "B", 
           "17" = "Hep", 
           "18" = "LSEC", 
           "19" = "T", 
           "20" = "Stellate")

clusts <- as.character(seur@meta.data[,"SCT_snn_res.0.5"])
seur$CellType <- ctmap[clusts]
seur$Major <- mjmap[clusts]

saveRDS(seur, seur_fn)


