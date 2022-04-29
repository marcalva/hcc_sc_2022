#!/u/local/apps/R/3.5.1/gcc-4.9.3_MKL-2018/bin/Rscript
#$ -S /u/local/apps/R/3.5.1/gcc-4.9.3_MKL-2018/bin/Rscript
#$ -cwd
#$ -j y
#$ -pe shared 8
#$ -l h_data=4G,h_vmem=32G,h_rt=8:00:00,highp
#$ -M malvarez@mail
#  Notify at beginning and end of job
#$ -m n
#$ -r n
#$ -o find_markers.R.log


setwd("../../")

library(Seurat)
library(future)
source("scripts/de.R")
plan("multiprocess", workers = 8)

maxsize <- 6000 * 1024 ^ 2
options(future.globals.maxSize = maxsize)

#=========================================
# Functions
#=========================================

# Get expressed genes
get_expr_genes <- function(x, cell_col, thresh = 0.1){
    cts <- as.character(unique(x@meta.data[,cell_col]))
    pct_expr_ct <- list()
    for (ct in cts){
        print(ct)
        keep <- x@meta.data[,cell_col] == ct
        counts <- x@assays$RNA@data[,keep]
        k <- counts > 0
        ks <- Matrix::rowSums(k)
        total_cells <- ncol(k)
        pct_cells <- ks / total_cells
        pct_expr_ct[[ct]] <- pct_cells
    }

    pct_expr_ct <- lapply(pct_expr_ct, function(ct) ct >= thresh)
    pct_expr <- Reduce("|", pct_expr_ct)
    expr_genes <- names(pct_expr)[pct_expr > 0]
    return(expr_genes)
}

create_dir <- function(p){
    dir.create(p, showWarnings=FALSE, recursive=TRUE)
}

#=========================================
# Set variables
#=========================================

# Set directories
dir_out <- "data/processed/d_sct_cca_all/"; create_dir(dir_out)
dir_marker <- "exp/d_sct_cca_all/cca_init/markers/"; create_dir(dir_marker)

meta <- read.table("data/sample/samples.csv", 
                   header = TRUE, 
                   row.names = 1, 
                   stringsAsFactors = FALSE,
                   sep = ",")

samples <- rownames(meta)

#=========================================
# Read in Seurat
#=========================================

gene_info <- read.table("data/ref/gencode26/gencode.v26.annotation.txt", 
                        header = TRUE, 
                        row.names = 1, 
                        stringsAsFactors = FALSE, 
                        sep = "\t")

# Read in Seurat objects
seur <- readRDS(paste0(dir_out, "seur.cca.rds"))

# Use log-normalized data for testing marker genes
DefaultAssay(seur) <- "RNA"
seur <- NormalizeData(seur, verbose = FALSE, scale.factor = 1000)

lgfc_thresh <- 0.1
p_thresh <- 0.05

ress <- c(0.2, 0.5, 1, 1.5, 2)

for (res in ress){

    Idents(seur) <- paste0("integrated_snn_res.", res)
    m <- FindAllMarkers(seur, 
                        test.use = "LR", 
                        logfc.threshold = lgfc_thresh, 
                        only.pos = TRUE)
    m[,"gene_name"] <- gene_info[m[,"gene"], "Name"]
    k <- m[,"p_val_adj"] < p_thresh
    m <- m[k,,drop=FALSE]
    ofn <- paste0(dir_marker, "markers.res", res, ".txt")
    write.table(m, 
                ofn, 
                row.names = FALSE, col.names = TRUE, 
                sep = "\t", quote = FALSE)
}

