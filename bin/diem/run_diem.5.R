#!/u/local/apps/R/3.5.1/gcc-4.9.3_MKL-2018/bin/Rscript
#$ -S /u/local/apps/R/3.5.1/gcc-4.9.3_MKL-2018/bin/Rscript
#$ -cwd
#$ -j y
#$ -pe shared 8
#$ -l h_data=4G,h_vmem=32G,h_rt=2:00:00,highp
#$ -M malvarez@mail
#  Notify at beginning and end of job
#$ -t 1-6
#$ -m a
#$ -r n
#$ -o run_diem.5.R.log.$TASK_ID

# Run DIEM

setwd("../../")

library(diem)
library(Matrix)
library(Seurat)
library(dplyr)
library(ggplot2)


#=========================================
# Functions
#=========================================

create_dir <- function(p){
    dir.create(p, recursive=TRUE, showWarnings=FALSE)
}

#=========================================
# Set sample and variables
#=========================================

task = Sys.getenv(x = "SGE_TASK_ID")
if (task == ""){
    stop("Set SGE_TASK_ID environment variable")
}
task = as.integer(task)

meta <- read.table("data/sample/samples.csv", 
                   header = TRUE, 
                   sep = ",", 
                   stringsAsFactors = FALSE, 
                   check.names = FALSE)
rownames(meta) <- meta[,"Sample"]

gtf_path <- "data/ref/gencode26/gencode.v26.annotation.fltrd.gtf"

gene_info <- read.table("data/ref/gencode26/gencode.v26.annotation.txt",
                        header = TRUE,
                        row.names = 1,
                        stringsAsFactors = FALSE,
                        sep = "\t")
info_keep <- c("Chrm", "Type", "Name")
mt_ens <- rownames(gene_info)[ gene_info[,"Chrm"] == "chrM" ]
mt_genes <- gene_info[gene_info[,"Chrm"] == "chrM", "Name"]
malat1 <- rownames(gene_info)[ gene_info[,"Name"] == "MALAT1" ]
malat1 <- malat1[1]

sample <- meta[task, "Sample"]

# Map sample IDs
sid <- c("Sample 1 Normal", "Sample 1 Tumor", 
         "Sample 2 NASH", "Sample 2 Tumor", 
         "Sample 3 NASH", "Sample 3 Tumor")
names(sid) <- c("S2_norm", "S1_hcc", 
                "S2_nash", "S2_hcc", 
                "S3_nash", "S3_hcc")

#=========================================
# Output directories
#=========================================

dir_out <- paste0("data/processed/diem/", sample, "/"); create_dir(dir_out)
dir_plot <- paste0("exp/diem/plots/", sample, "/"); create_dir(dir_plot)
dir_fltr <- paste0("exp/diem/fltr_ids/", sample, "/"); create_dir(dir_fltr)
dir_de <- paste0("exp/diem/DE/", sample, "/"); create_dir(dir_de)
plot_prefix <- paste0(dir_plot, sample, ".") # File prefix for dist plot

#=========================================
# Read data
#=========================================

sce <- readRDS(paste0(dir_out, sample, ".diem_sce.rds"))

#bg_clust <- list("S1_hcc" = as.character(c(1,2,3,4,6,10,22)), 
#                 "S2_norm" = as.character(c(1,25,7,26,20,23,30,18,6,8,5,9,3,2)),
#                 "S2_hcc" = as.character(c(1,5,21,24,27,14)),
#                 "S2_nash" = as.character(c(1,29,13,2,4,31,20,5,11,16,27)),
#                 "S3_hcc" = as.character(c(1,6,17,22,20,18,3,7,11,12,8,5,9,21)),
#                 "S3_nash" = as.character(c(1,4,5,15,26,21)))

bg_clust <- list("S1_hcc" = as.character(c(1,26)), 
                 "S2_norm" = as.character(c(1,6)),
                 "S2_hcc" = as.character(c(1,19)), 
                 "S2_nash" = as.character(c(1)),
                 "S3_hcc" = as.character(c(1,9)), 
                 "S3_nash" = as.character(c(1)))

td <- sce@test_data

#=========================================
# Add meta data
#=========================================

md <- sce@test_data

md[,"SampleID"] <- factor(sid[sample])

md[,"orig.ident"] <- sample
md[,"Patient"] <- meta[sample, "Patient"]
md[,"Tumor"] <- meta[sample, "Tumor"]
md[,"Pheno"] <- meta[sample, "Pheno"]
md[,"SeqDate"] <- meta[sample, "SeqDate"]

md[,"Cluster"] <- factor(paste0(sample, ":", md[,"Cluster"]))

#=========================================
# Seurat object uncorrected counts
#=========================================

sdk <- (! td[,"Cluster"] %in% bg_clust[[sample]] ) &
       (td[,"pct.mt"] < 50) & 
       (td[,"n_genes"] >= 200)

k <- rownames(td)[sdk]
counts <- sce@counts[,k]
md <- md[k,]

seur <- Seurat::CreateSeuratObject(counts = counts, 
                                   meta.data = md, 
                                   names.delim = "}", 
                                   project = sample)

vars.to.regress <- NULL
seur <- NormalizeData(seur, normalization.method = "LogNormalize", scale.factor = 200)
seur <- FindVariableFeatures(seur, selection.method = "vst", nfeatures = 2000)
seur <- ScaleData(seur, vars.to.regress = vars.to.regress)

seur <- SCTransform(seur, vars.to.regress = vars.to.regress)
# seur <- RunPCA(seur, npcs = 30, verbose = TRUE)

saveRDS(seur, paste0(dir_out, "seur.sct.rds"))


