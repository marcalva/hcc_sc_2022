
setwd("../../")

library(Seurat)
library(Matrix)
library(future)
plan("multicore", workers = 4)
maxsize <- 8000 * 1024 ^ 2
options(future.globals.maxSize = maxsize)

#=========================================
# Functions
#=========================================

create_dir <- function(p){
    dir.create(p, showWarnings=FALSE, recursive=TRUE)
}

get_lm_fast <- function(x, psc = 1){
    x <- Matrix::t(x)

    cmeans <- sapply(1:ncol(x), function(i1){
                     i2 = i1+1
                     ip1 <- x@p[i1]
                     ip2 <- x@p[i2]
                     nnz <- ip2 - ip1
                     nz <- nrow(x) - nnz
                     if (ip1 == ip2){
                         ent <- 0
                     }else{
                         ent <- exp(x@x[(ip1+1):ip2]) - 1 + psc
                     }
                     csum <- sum(c(ent, psc * nz))
                     if (is.na(csum)) print(which(is.na(ent)))
                     cmean <- log2( csum / (nz + nnz) )
                     if (is.na(cmean)) message("NA ", i1)
                     return(cmean) })

    names(cmeans) <- colnames(x)
    return(cmeans)
}

# get logFC for all genes
get_lfc <- function(seur, ct_id, psc = 1){
    md <- seur@meta.data
    if (!ct_id %in% colnames(md)){
        stop(ct_id, " not found in seurat object")
    }
    cts <- sort(as.character(unique(md[,ct_id])))
    lfc_l <- list()
    for (ct in cts){
        message(ct)
        nct <- setdiff(cts, ct)

        k1 <- rownames(md)[md[,ct_id] %in% ct]
        k2 <- rownames(md)[md[,ct_id] %in% nct]

        mat1 <- seur@assays$RNA@data[,k1,drop=FALSE]
        mat2 <- seur@assays$RNA@data[,k2,drop=FALSE]
        # kg <- rowSums(mat1) > 0 & rowSums(mat2) > 0
        # mat1 <- mat1[kg,]
        # mat2 <- mat2[kg,]

        m1 <- get_lm_fast(mat1)
        m2 <- get_lm_fast(mat2)

        lfc <- m1 - m2
        lfc_l[[ct]] <- lfc
    }

    g_all <- lapply(lfc_l, names)
    g_u <- Reduce(union, g_all)
    lfc_l <- lapply(lfc_l, function(x) x[g_u])
    lfc_mat <- do.call(cbind, lfc_l)

    return(lfc_mat)
}

#=========================================
# Set variables
#=========================================

# Set directories
dir_marker <- "exp/sharma_aiz/markers/"; create_dir(dir_marker)

#=========================================
# Read in Seurat
#=========================================

gene_info <- read.table("data/ref/gencode26/gencode.v26.annotation.txt", 
                        header = TRUE, 
                        row.names = 1, 
                        stringsAsFactors = FALSE, 
                        sep = "\t")

# Read in Seurat objects
fn <- "data/processed/sharma_aiz/liver.int_rand.rds"
integrated <- readRDS(fn)

# Use log-normalized data for testing marker genes
# Scale to 1,000 counts
DefaultAssay(integrated) <- "RNA"
integrated <- NormalizeData(integrated, verbose = FALSE, scale.factor = 1000)
md <- integrated@meta.data

lgfc_thresh <- 0.1
p_thresh <- 0.05

#=========================================
# Resolution 0.5
#=========================================

ct_id <- "integrated_snn_res.0.5"

Idents(integrated) <- ct_id
cts <- sort(as.character(unique(md[,ct_id])))

# get logFC for all genes
lfc_mat <- get_lfc(integrated, ct_id = "integrated_snn_res.0.5", psc = 1)
lfc_mat <- cbind(gene_info[rownames(lfc_mat),], lfc_mat)

ofn <- paste0(dir_marker, "markers.res.0.5.log_fc.txt")
write.table(lfc_mat, 
            ofn, 
            row.names = TRUE, col.names = NA, 
            sep = "\t", quote = FALSE)

m <- FindAllMarkers(integrated, 
                    test.use = "LR", 
                    logfc.threshold = lgfc_thresh, 
                    only.pos = TRUE)
k <- m[,"p_val_adj"] < p_thresh
m <- m[k,,drop=FALSE]

m <- cbind(gene_info[m[,"gene"],], m)

ofn <- paste0(dir_marker, "markers.res.0.5.txt")
write.table(m, 
            ofn, 
            row.names = FALSE, col.names = TRUE, 
            sep = "\t", quote = FALSE)

#=========================================
# Resolution 1
#=========================================

Idents(integrated) <- "integrated_snn_res.1"

# get logFC for all genes
lfc_mat <- get_lfc(integrated, ct_id = "integrated_snn_res.1", psc = 1)
lfc_mat <- cbind(gene_info[rownames(lfc_mat),], lfc_mat)

ofn <- paste0(dir_marker, "markers.res.1.log_fc.txt")
write.table(lfc_mat, 
            ofn, 
            row.names = TRUE, col.names = NA, 
            sep = "\t", quote = FALSE)

m <- FindAllMarkers(integrated, 
                    test.use = "LR", 
                    logfc.threshold = lgfc_thresh, 
                    only.pos = TRUE)
k <- m[,"p_val_adj"] < p_thresh
m <- m[k,,drop=FALSE]

m <- cbind(gene_info[m[,"gene"],], m)

ofn <- paste0(dir_marker, "markers.res.1.txt")
write.table(m, 
            ofn, 
            row.names = FALSE, col.names = TRUE, 
            sep = "\t", quote = FALSE)

