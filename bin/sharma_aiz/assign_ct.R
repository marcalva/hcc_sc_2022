
# predict cell-types using annotators

setwd("../../")

suppressPackageStartupMessages(library(Seurat))

#=========================================
# Functions
#=========================================

create_dir <- function(x){
    dir.create(x, showWarnings = FALSE, recursive = TRUE)
}

#=========================================
#=========================================

# gene annotations
gene_info <- read.table("data/ref/gencode26/gencode.v26.annotation.txt",
                        header = TRUE,
                        row.names = 1,
                        stringsAsFactors = FALSE,
                        sep = "\t")

fn <- "data/processed/sharma_aiz/liver.int_rand.rds"
integrated <- readRDS(fn)

dir_exp <- "exp/sharma_aiz/clust2ct/"
dir.create(dir_exp, showWarnings = FALSE, recursive = TRUE)

#=========================================
# cell-type map
#=========================================

ctmap <- c("0" = "Hep_1", 
           "1" = "Hep_2", 
           "2" = "T_1", 
           "3" = "T_2", 
           "4" = "Endo_1", 
           "5" = "T_3", 
           "6" = "Chol", 
           "7" = "Endo_2", 
           "8" = "T_NK", 
           "9" = "Hep_3", 
           "10" = "Myel_1", 
           "11" = "T_4", 
           "12" = "Endo_3", 
           "13" = "Myel_2", 
           "14" = "Myel_3", 
           "15" = "T_5", 
           "16" = "Stell", 
           "17" = "Myel_4",
           "18" = "T_reg",
           "19" = "B",
           "20" = "Prol",
           "21" = "T_6",
           "22" = "Myel_5", 
           "23" = "T_7", 
           "24" = "T_8")

mjmap <- sapply(ctmap, function(s) strsplit(s, "_")[[1]][1])

# save maps
datf <- data.frame("cluster" = names(ctmap), "cell_type" = ctmap)
out_fn <- paste0(dir_exp, "res1_to_fine.txt")
write.table(datf, out_fn, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

datf <- data.frame("cluster" = names(mjmap), "cell_type" = mjmap)
out_fn <- paste0(dir_exp, "res1_to_main.txt")
write.table(datf, out_fn, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

# save in seurat object
clusts <- integrated$integrated_snn_res.1
cts <- ctmap[as.character(clusts)]
integrated$cell_type_fine <- factor(cts)

clusts <- integrated$integrated_snn_res.1
mjs <- mjmap[as.character(clusts)]
integrated$cell_type_main <- factor(mjs)

fn <- "data/processed/sharma_aiz/liver.int_rand.rds"
saveRDS(integrated, fn)

# add cell type to markers
fn <- "exp/sharma_aiz/markers/markers.res.1.txt"
mk <- read.table(fn, header=TRUE)
mk[,"cluster"] <- ctmap[as.character(mk[,"cluster"])]
fn <- "exp/sharma_aiz/markers/markers.cell_type_fine.txt"
write.table(mk, fn, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

# add cell type to lfc
fn <- "exp/sharma_aiz/markers/markers.res.1.log_fc.txt"
lfc <- read.table(fn, header=TRUE, check.names=FALSE, row.names=NULL)
lfc1 <- lfc[,1:7]
lfc2 <- lfc[,8:ncol(lfc)]
colnames(lfc2) <- ctmap[colnames(lfc2)]
lfc <- cbind(lfc1, lfc2)
fn <- "exp/sharma_aiz/markers/markers.cell_type_fine.log_fc.txt"
write.table(lfc, fn, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

