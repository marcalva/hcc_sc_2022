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
#$ -o npc.1.R.log

# Cluster non-hepatocytes step 1

setwd("../../")

library(Seurat)
source("scripts/standard_seurat.R")
library(ggplot2)
library(scales)
library(gridExtra)
source("scripts/plot.R")

library(future)
plan("multiprocess", workers = 8)
maxsize <- 6000 * 1024 ^ 2
options(future.globals.maxSize = maxsize)

create_dir <- function(p){
    dir.create(p, showWarnings=FALSE, recursive=TRUE)
}

#=========================================
# Set variables
#=========================================

# Set directories
dir_out <- "data/processed/d_sct_cca_all/"; create_dir(dir_out)
dir_marker <- "exp/d_sct_cca_all/npc1/markers/"; create_dir(dir_marker)
dir_plot <- "exp/d_sct_cca_all/npc1/plots/"; create_dir(dir_plot)

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
seur <- readRDS("data/processed/d_sct_cca_all/seur.cca.rds")

#=========================================
# Subclustering
#=========================================

seur$Group <- "Nonparenchymal"
seur@meta.data[seur$Major == "Hepatocyte","Group"] <- "Hepatocyte"

# Get parenchymal non-parenchymal genes
infn <- "exp/d_sct_cca_all/cca_init/markers/markers.Hep_v_Nonhep.txt"
m <- read.table(infn, 
                header = TRUE,
                stringsAsFactors = FALSE,
                sep = "\t")

hep_feat <- m[m$cluster == "Hepatocyte", "gene"]
np_feat <- m[m$cluster == "Nonparenchymal", "gene"]
all_feat <- rownames(seur@assays$RNA@data)
feat1 <- setdiff(all_feat, hep_feat)

seur <- PercentageFeatureSet(seur, features = hep_feat, col.name = "pct.hep", assay = "RNA")
seur <- PercentageFeatureSet(seur, features = np_feat, col.name = "pct.npc", assay = "RNA")
seur$logmt <- log1p(seur$pct.mt)

# Subset
k <- seur@meta.data$Group == "Nonparenchymal"
keep <- rownames(seur@meta.data)[k]
seur_s <- subset(seur, cells = keep)

#=========================================
# Non-parenchymal
#=========================================

seur_s.l <- SplitObject(seur_s, split.by = "orig.ident")
seur_s.l <- lapply(seur_s.l, SCTransform)
seur_s <- sct_cca(seur_s.l)
seur_s <- seurat_cluster(seur_s, resolution = c(0.8, 1))

dir_out <- "data/processed/d_sct_cca_all/npc1/"; create_dir(dir_out)
saveRDS(seur_s, paste0(dir_out, "seur.cca.rds"))

#=========================================
# Plots
#=========================================

plot_umap_labels(seur_s,
                 colname = "integrated_snn_res.0.8",
                 legend_title = "Res. 0.8",
                 label = TRUE,
                 size = 0.3, label_size = 4,
                 alpha = 0.5)
ggsave(paste0(dir_plot, "UMAP.res.0.8.jpeg"),
       width = 5,
       height = 5,
       units = "in")

plot_umap_labels(seur_s,
                 colname = "integrated_snn_res.1",
                 legend_title = "Res. 1",
                 label = TRUE,
                 size = 0.3, label_size = 4,
                 alpha = 0.5)
ggsave(paste0(dir_plot, "UMAP.res.1.jpeg"),
       width = 5,
       height = 5,
       units = "in")

p <- plot_umap_c(seur_s,
                 colname = "pct.hep",
                 legend_title = "Hep%",
                 size = 0.3,
                 alpha = 0.5)
ggsave(paste0(dir_plot, "UMAP.pct_hep.jpeg"),
       width = 5,
       height = 5,
       units = "in")


p <- plot_umap_c(seur_s,
                 colname = "pct.mt",
                 legend_title = "MT%",
                 trans = "log1p",
                 size = 0.3,
                 alpha = 0.5)
ggsave(paste0(dir_plot, "UMAP.pct_mt.jpeg"),
       width = 5,
       height = 5,
       units = "in")

plot_umap_labels(seur_s,
                 colname = "SampleID",
                 legend_title = "Sample",
                 label = FALSE,
                 size = 0.3,
                 alpha = 0.5)
ggsave(paste0(dir_plot, "UMAP.sample.jpeg"),
       width = 5,
       height = 5,
       units = "in")

#=========================================
# Find markers
#=========================================


seur_s@active.assay <- "RNA"
seur_s <- NormalizeData(seur_s, verbose = FALSE, scale.factor = 1000)
Idents(seur_s) <- "integrated_snn_res.1"

mall <- FindAllMarkers(seur_s, test.use = "LR", logfc.threshold = 0.1, only.pos = TRUE)
mall[,"gene_name"] <- gene_info[mall$gene, "Name"]
k <- mall[,"p_val_adj"] < 0.05
mall <- mall[k,,drop=FALSE]

write.table(mall, 
            paste0(dir_marker, "seur.res1.markers.txt"), 
            row.names = FALSE, col.names = TRUE, 
            quote = FALSE, sep = "\t")

mall$HepM <- mall[,"gene"] %in% hep_feat
tapply(mall$HepM, mall$cluster, table)

# res1 clusters 1 and 6 have > 50% hepatocyte markers.

