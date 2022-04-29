#!/u/local/apps/R/3.5.1/gcc-4.9.3_MKL-2018/bin/Rscript
#$ -S /u/local/apps/R/3.5.1/gcc-4.9.3_MKL-2018/bin/Rscript
#$ -cwd
#$ -j y
#$ -pe shared 8
#$ -l h_data=4G,h_vmem=32G,h_rt=4:00:00,highp
#$ -M malvarez@mail
#  Notify at beginning and end of job
#$ -m n
#$ -r n
#$ -o merge_init.R.log

# Merge hepatocytes and non-hepatocytes

setwd("../../")

library(Seurat)
source("scripts/standard_seurat.R")
library(ggplot2)
library(scales)
library(gridExtra)
source("scripts/plot.R")

library(future)
plan("multiprocess", workers = 8)
maxsize <- 8000 * 1024 ^ 2
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
dir_marker <- "exp/d_sct_cca_all/merge/markers/"; create_dir(dir_marker)

meta <- read.table("data/sample/samples.csv",
                   header = TRUE,
                   row.names = 1,
                   stringsAsFactors = FALSE,
                   sep = ",")

samples <- rownames(meta)

gene_info <- read.table("data/ref/gencode26/gencode.v26.annotation.txt",
                        header = TRUE,
                        row.names = 1,
                        stringsAsFactors = FALSE,
                        sep = "\t")

# Read in Seurat objects
seur_hep <- readRDS("data/processed/d_sct_cca_all/hep2/seur.cca.rds")
seur_npc <- readRDS("data/processed/d_sct_cca_all/npc2/seur.cca.rds")

#=========================================
# Merge
#=========================================

seur <- merge(seur_hep, seur_npc)

# Run SCT + CCA for visualization, since we already have clusters.
DefaultAssay(seur) <- "RNA"
seur.l <- SplitObject(seur, split.by = "orig.ident")
for (i in 1:length(seur.l)){
    seur.l[[i]] <- SCTransform(seur.l[[i]])
}
seur <- sct_cca(seur.l)

n_dim <- 30
python_lib <- "/u/local/apps/python/3.6.1-shared/bin/python3.6"
library(reticulate)
reticulate::use_python(python_lib)
seur <- FindNeighbors(seur, dims = 1:n_dim, verbose = FALSE)
seur <- RunUMAP(seur, dims = 1:n_dim, reduction = "pca", verbose = FALSE)

# flip UMAP coordinates
umapc <- seur@reductions$umap@cell.embeddings
umapc <- -umapc
seur@reductions$umap@cell.embeddings <- umapc

# Remove old columns
md <- seur@meta.data
cols2drop <- grep("^integrated_snn_res", colnames(md), value = TRUE)
cols2drop <- c(cols2drop, c("nCount_SCT", "nFeature_SCT", 
                            "nCount_integrated", "nFeature_integrated", 
                            "seurat_clusters", "logmt"))
cols2keep <- setdiff(colnames(md), cols2drop)
md <- md[,cols2keep]
seur@meta.data <- md

# Update columns
seur$SampleID <- factor(seur$SampleID)

seur$Tumor <- factor(seur$Tumor)
seur$Pheno <- factor(seur$Pheno)
seur$SeqDate <- factor(seur$SeqDate)

seur$Group <- factor(seur$Group)
seur$CellType <- factor(seur$CellType)
seur$Major <- factor(seur$Major)

# Save
dir_out <- "data/processed/d_sct_cca_all/merge/"; create_dir(dir_out)
saveRDS(seur, paste0(dir_out, "seur.cca.rds"))


#=========================================
# Find markers
#=========================================

seur@active.assay <- "RNA"
seur <- NormalizeData(seur, verbose = FALSE, scale.factor = 1000)

# Cell type
Idents(seur) <- "CellType"
mall <- FindAllMarkers(seur, test.use = "LR", logfc.threshold = 0.1, only.pos = TRUE)
mall[,"gene_name"] <- gene_info[mall$gene, "Name"]
k <- mall[,"p_val_adj"] < 0.05
mall <- mall[k,,drop=FALSE]

write.table(mall, 
            paste0(dir_marker, "seur.CellType.markers.txt"), 
            row.names = FALSE, col.names = TRUE, 
            quote = FALSE, sep = "\t")

expr_genes <- get_expr_genes(seur, "CellType", thresh = 0.1)
gene_info_e <- gene_info[expr_genes,]
write.table(gene_info_e,
            paste0(dir_marker, "exprd_gn.CellType.txt"),
            row.names = TRUE, col.names = NA,
            sep = "\t", quote = FALSE)

# Major
Idents(seur) <- "Major"
mall <- FindAllMarkers(seur, test.use = "LR", logfc.threshold = 0.1, only.pos = TRUE)
mall[,"gene_name"] <- gene_info[mall$gene, "Name"]
k <- mall[,"p_val_adj"] < 0.05
mall <- mall[k,,drop=FALSE]

write.table(mall, 
            paste0(dir_marker, "seur.Major.markers.txt"), 
            row.names = FALSE, col.names = TRUE, 
            quote = FALSE, sep = "\t")

expr_genes <- get_expr_genes(seur, "Major", thresh = 0.1)
gene_info_e <- gene_info[expr_genes,]
write.table(gene_info_e,
            paste0(dir_marker, "exprd_gn.Major.txt"),
            row.names = TRUE, col.names = NA,
            sep = "\t", quote = FALSE)


