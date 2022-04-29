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
#$ -o run_diem.3.R.log.$TASK_ID

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

args = commandArgs(TRUE)
task = as.integer(args[1])

task = Sys.getenv(x = "SGE_TASK_ID")
if (task == ""){
    stop("Set SGE_TASK_ID environment variable")
}

meta <- read.table("data/sample/samples.csv", 
                   header = TRUE, 
                   sep = ",", 
                   stringsAsFactors = FALSE, 
                   check.names = FALSE)

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

#=========================================
# Output directories
#=========================================

dir_out <- paste0("data/processed/diem/", sample, "/"); create_dir(dir_out)
dir_plot <- paste0("exp/diem/plots/", sample, "/"); create_dir(dir_plot)
dir_fltr <- paste0("exp/diem/fltr_ids/", sample, "/"); create_dir(dir_fltr)
dir_de <- paste0("exp/diem/DE/", sample, "/"); create_dir(dir_de)
plot_prefix <- paste0(dir_plot, sample, ".") # File prefix for dist plot

#=========================================
# Identify debris clusters
#=========================================

sce <- readRDS(paste0(dir_out, sample, ".diem_sce.rds"))
sce <- assign_clusters(sce)

# Find markers
rc <- raw_counts(sce)
markers <- de_ttest_all(rc, labels = sce@test_data$Cluster, 
                        normalize = TRUE, sf = 1e3)
markers$Gene = gene_info[markers$gene, "Name"]
write.table(markers,  
            paste0(dir_de, sample, ".clust_markers.txt"), 
            row.names = TRUE, col.names = NA, 
            quote = FALSE, sep = "\t")

# Plot cluster means
p <- plot_clust(sce, feat_x = "n_genes", feat_y = "pct.mt", 
           log_x = TRUE, log_y = FALSE) + 
ggtitle(sample)
ggsave(paste0(plot_prefix, "cluster.n_gene.pct_mt.jpeg"), 
       width = 5, height = 5)

p <- plot_clust(sce, feat_x = "n_genes", feat_y = "pct.spl", 
           log_x = TRUE, log_y = FALSE) + 
ggtitle(sample)
ggsave(paste0(plot_prefix, "cluster.n_gene.pct_spl.jpeg"), 
       width = 5, height = 5)

p <- plot_clust(sce, feat_x = "n_genes", feat_y = "pct.malat1", 
           log_x = TRUE, log_y = FALSE) + 
ggtitle(sample)
ggsave(paste0(plot_prefix, "cluster.n_gene.pct_malat1.jpeg"), 
       width = 5, height = 5)

# Plot number droplets in clusters
ndrop <- table(sce@test_data[,"Cluster"])
datf_n <- data.frame(ndrop)
colnames(datf_n) <- c("Cluster", "nDrop")
p <- ggplot(datf_n, aes(x = Cluster, y = nDrop)) + 
geom_bar(stat = "identity") + 
theme_minimal() + 
xlab("Cluster") + ylab("Number of Droplets") + 
ggtitle(sample)
ggsave(paste0(plot_prefix, "cluster.n_droplets.jpeg"), 
       width = 5, height = 5)


# Save DIEM object
cat("Saving\n")
saveRDS(sce, paste0(dir_out, sample, ".diem_sce.rds"))

