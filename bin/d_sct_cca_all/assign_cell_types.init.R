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

# Assign cell types to initial clustering

setwd("../../")

library(Seurat)
library(future)
plan("multiprocess", workers = 8)

maxsize <- 4000 * 1024 ^ 2
options(future.globals.maxSize = maxsize)

#=========================================
# Functions
#=========================================

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
dir_out <- "data/processed/d_sct_cca_all/"; create_dir(dir_out)
dir_plot <- "exp/d_sct_cca_all/cca_init/plots/"; create_dir(dir_plot)
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

ctid <- "integrated_snn_res.0.2"

to_hep_name <- sort(c(0,1,7,9))
to_hep <- rep("Hepatocyte", length(to_hep_name))
names(to_hep) <- to_hep_name

# to_t <- c("T-IL7R", "T-AOAH")
to_t <- c("Lymphoid")
names(to_t) <- c(4)

# to_kup <- c("Kupf-FLT3", "Kupf-NDST3")
to_kup <- c("Myeloid")
names(to_kup) <- c(3)

to_lsec <- c("HSEC")
names(to_lsec) <- c(2)

to_stel <- "Stellate"
names(to_stel) <- 6

to_chol <- c("Cholangiocyte")
names(to_chol) <- c(5)

to_neu <- "Neural"
names(to_neu) <- 8

to <- c(to_hep, to_t, to_kup, to_lsec, to_stel, to_chol, to_neu)

seur@meta.data[,"Major"] <- to[as.character(seur@meta.data[,ctid])]
seur$Major <- factor(seur$Major)
Idents(seur) <- "Major"

# Find cell type markers
seur@active.assay = "RNA"
seur <- NormalizeData(seur, scale.factor = 1e3)
tu <- "LR"
lfct <- 0.1

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

# Hepatocyte vs. non-parenchymal
seur$Group <- "Nonparenchymal"
seur@meta.data[seur$Major == "Hepatocyte","Group"] <- "Hepatocyte"

Idents(seur) <- "Group"
seur@active.assay <- "RNA"
seur <- NormalizeData(seur, verbose = FALSE, scale.factor = 1000)
tu <- "LR"
lfct <- 0.1

m <- FindAllMarkers(seur,
                    test.use = tu,
                    logfc.threshold = lfct,
                    only.pos = TRUE)
m[,"gene_name"] <- gene_info[rownames(m), "Name"]
k <- m[,"p_val_adj"] < 0.05
m <- m[k,,drop=FALSE]

write.table(m, 
            paste0(dir_marker, "markers.Hep_v_Nonhep.txt"), 
            row.names = FALSE, col.names = TRUE, 
            sep = "\t", quote = FALSE)

# Save
seur@active.assay = "integrated"
saveRDS(seur, paste0(dir_out, "seur.cca.rds"))

