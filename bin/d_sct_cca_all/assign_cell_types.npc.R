#!/u/local/apps/R/3.5.1/gcc-4.9.3_MKL-2018/bin/Rscript
#$ -S /u/local/apps/R/3.5.1/gcc-4.9.3_MKL-2018/bin/Rscript
#$ -cwd
#$ -j y
#$ -pe shared 8
#$ -l h_data=4G,h_vmem=32G,h_rt=2:00:00,highp
#$ -M malvarez@mail
#  Notify at beginning and end of job
#$ -m n
#$ -r n
#$ -o assign_cell_types.R.log

# Assign clusters

setwd("../../")

library(Seurat)
library(future)
plan("multiprocess", workers = 8)

maxsize <- 4000 * 1024 ^ 2
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

gene_info <- read.table("data/ref/gencode26/gencode.v26.annotation.txt", 
                        header = TRUE, 
                        row.names = 1, 
                        stringsAsFactors = FALSE, 
                        sep = "\t")

#=========================================
# Read in Seurat
#=========================================

# Set directories
dir_out <- "data/processed/d_sct_cca_all/npc2/"; create_dir(dir_out)
dir_plot <- "exp/d_sct_cca_all/npc2/plots/"; create_dir(dir_plot)
dir_marker <- "exp/d_sct_cca_all/npc2/markers/"; create_dir(dir_marker)

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

#=========================================
# Read in Seurat
#=========================================

# Read in Seurat objects
infn <- "data/processed/d_sct_cca_all/npc2/seur.cca.rds"
seur <- readRDS(infn)

ctid <- "integrated_snn_res.1"

to_e_name <- sort(c(14, 5, 8, 7, 0, 16, 15))
to_e <- paste0("HSEC-", 1:length(to_e_name))
names(to_e) <- to_e_name

to_st_name <- sort(c(9, 11))
to_st <- paste0("Stellate-", 1:length(to_st_name))
names(to_st) <- to_st_name

to_ne <- "Neural"
names(to_ne) <- 10

to_k_name <- sort(c(1, 3, 18))
to_k <- paste0("Kupffer-", 1:length(to_k_name))
names(to_k) <- to_k_name

to_t_name <- sort(c(6, 12, 2))
to_t <- paste0("T-", 1:length(to_t_name))
names(to_t) <- to_t_name

to_chol_name <- sort(c(4, 13, 17))
to_chol <- paste0("Chol-", 1:length(to_chol_name))
names(to_chol) <- to_chol_name

to <- c(to_e, to_st, to_ne, to_k, to_t, to_chol)

seur@meta.data[,"CellType"] <- to[as.character(seur@meta.data[,ctid])]
seur$CellType <- factor(seur$CellType)

# Find cell type markers
seur@active.assay = "RNA"
seur <- NormalizeData(seur, scale.factor = 1e3)
tu <- "LR"
lfct <- 0.1

# CellType
Idents(seur) <- "CellType"
m <- FindAllMarkers(seur, 
                    test.use = tu, 
                    logfc.threshold = lfct, 
                    only.pos = TRUE)
m[,"gene_name"] <- gene_info[m[,"gene"], "Name"]
write.table(m, 
            paste0(dir_marker, "markers.CellType.txt"), 
            row.names = FALSE, col.names = TRUE, 
            sep = "\t", quote = FALSE)

expr_genes <- get_expr_genes(seur, "CellType", thresh = lfct)
gene_info_e <- gene_info[expr_genes,]
write.table(gene_info_e, 
            paste0(dir_marker, "exprd_gn.CellType.txt"), 
            row.names = TRUE, col.names = NA, 
            sep = "\t", quote = FALSE)

# Major
to_e_name <- sort(c(14, 5, 8, 7, 0, 16, 15))
to_e <- rep("HSEC", length(to_e_name))
names(to_e) <- to_e_name

to_st_name <- sort(c(9, 11))
to_st <- rep("Stellate", length(to_st_name))
names(to_st) <- to_st_name

to_ne <- "Neural"
names(to_ne) <- 10

to_k_name <- sort(c(1, 3, 18))
to_k <- rep("Myeloid", length(to_k_name))
names(to_k) <- to_k_name

to_t_name <- sort(c(6, 12, 2))
to_t <- rep("Lymphoid", length(to_t_name))
names(to_t) <- to_t_name

to_chol_name <- sort(c(4, 13, 17))
to_chol <- rep("Cholangiocyte", length(to_chol_name))
names(to_chol) <- to_chol_name

to <- c(to_e, to_st, to_ne, to_k, to_t, to_chol)

seur@meta.data[,"Major"] <- to[as.character(seur@meta.data[,ctid])]
seur$Major <- factor(seur$Major)

Idents(seur) <- "Major"
m <- FindAllMarkers(seur, 
                    test.use = tu, 
                    logfc.threshold = lfct, 
                    only.pos = TRUE)
m[,"gene_name"] <- gene_info[m[,"gene"], "Name"]
write.table(m, 
            paste0(dir_marker, "markers.Major.txt"), 
            row.names = FALSE, col.names = TRUE, 
            sep = "\t", quote = FALSE)

expr_genes <- get_expr_genes(seur, "Major", thresh = lfct)
gene_info_e <- gene_info[expr_genes,]
write.table(gene_info_e, 
            paste0(dir_marker, "exprd_gn.Major.txt"), 
            row.names = TRUE, col.names = NA, 
            sep = "\t", quote = FALSE)


# Save
Idents(seur) <- "CellType"
seur@active.assay = "integrated"
saveRDS(seur, paste0(dir_out, "seur.cca.rds"))


# Plot

library(ggplot2)
library(scales)
library(gridExtra)
source("scripts/plot.R")

set.seed(1)
cols1 <- gg_color_hue(length(unique(seur$CellType)))
cols1 <- sample(cols1)
plot_umap_labels(seur, 
                 colname = "CellType", 
                 legend_title = "Cell type", 
                 label = TRUE, 
                 colors = cols1, 
                 size = 0.3, label_size = 4,
                 alpha = 0.5)
ggsave(paste0(dir_plot, "UMAP.CellType.jpeg"), 
       width = 5, 
       height = 5, 
       units = "in")

set.seed(1)
cols1 <- gg_color_hue(length(unique(seur$Major)))
cols1 <- sample(cols1)
plot_umap_labels(seur, 
                 colname = "Major", 
                 legend_title = "Cell type", 
                 label = TRUE, 
                 colors = cols1, 
                 size = 0.3, label_size = 4,
                 alpha = 0.5)
ggsave(paste0(dir_plot, "UMAP.Major.jpeg"), 
       width = 5, 
       height = 5, 
       units = "in")

plot_umap_labels(seur, 
                 colname = "Tumor", 
                 legend_title = "Status", 
                 label = TRUE, 
                 colors = cols1, 
                 size = 0.3, label_size = 4,
                 alpha = 0.5)
ggsave(paste0(dir_plot, "UMAP.Tumor.jpeg"), 
       width = 5, 
       height = 5, 
       units = "in")

