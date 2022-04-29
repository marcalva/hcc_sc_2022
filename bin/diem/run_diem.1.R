#!/u/local/apps/R/3.5.1/gcc-4.9.3_MKL-2018/bin/Rscript
#$ -S /u/local/apps/R/3.5.1/gcc-4.9.3_MKL-2018/bin/Rscript
#$ -cwd
#$ -j y
#$ -pe shared 8
#$ -l h_data=4G,h_vmem=32G,h_rt=1:00:00,highp
#$ -M malvarez@mail
#  Notify at beginning and end of job
#$ -t 1-6
#$ -m a
#$ -r n
#$ -o run_diem.1.R.log.$TASK_ID

# Initialize DIEM object

setwd("../../")

library(diem)
library(Matrix)
library(Seurat)
library(ggplot2)


#=========================================
# Functions
#=========================================

# @param p is the path to the STAR solo output that contains the matrix.mtx,
#   features.tsv and barcodes.tsv files
# @param gene_col is the column index to use for the gene names in the 
#   barcodes.tsv file
# @param sep is the separator if gene names are not unique, 
#   the argument to make.unique function
read_counts <- function(p, gene_col = 1, sep = "."){
    mtx_file <- file.path(p, "matrix.mtx")
    genes_file <- file.path(p, "features.tsv")
    barcode_file <- file.path(p, "barcodes.tsv")

    counts <- Matrix::readMM(mtx_file)
    barcode_names <- readLines(barcode_file)
    genes <- read.delim(genes_file, 
                        header = FALSE, 
                        sep = "\t", 
                        stringsAsFactors = FALSE)


    colnames(counts) <- make.unique(barcode_names, sep=sep)
    rownames(counts) <- make.unique(genes[,gene_col], sep=sep)

    counts <- as(counts, "CsparseMatrix")

    return(counts)
}

create_dir <- function(p){
    dir.create(p, recursive=TRUE, showWarnings=FALSE)
}

sce_feat_pct <- function(x, ids, colname){
    ids <- lapply(ids, function(i) grep(i, rownames(x@counts), value = TRUE))
    ids <- unlist(ids)
    x <- get_gene_pct(x = x, genes = ids, name = colname)
    return(x)
}

sce_rm_feat <- function(x, ids){
    ids <- lapply(ids, function(i) grep(i, rownames(x@counts), value = TRUE))
    ids <- unlist(ids)
    keep_ids <- setdiff(rownames(x@counts), ids)
    x@counts <- x@counts[keep_ids,]
    x@gene_data <- x@gene_data[keep_ids,]
    x@droplet_data[,"total_counts"] <- colSums(x@counts[,rownames(x@droplet_data)])
    x@droplet_data[,"n_genes"] <- colSums(x@counts[,rownames(x@droplet_data)] > 0)
    d <- rownames(x@test_data)
    x@test_data[,"total_counts"] <- x@droplet_data[d,"total_counts"]
    x@test_data[,"n_genes"] <- x@droplet_data[d,"n_genes"]
    return(x)
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
rb_genes <- readLines("data/ref/ribosomal_genes.txt")
rb_ids <- rownames(gene_info)[gene_info$Name %in% rb_genes]

rbr_genes <- c("CH507-513H4.1", "CH507-528H12.1")
rbr_ids <- c("ENSG00000278996.1", "ENSG00000280441.2")
hb_ids <- c("ENSG00000206172.8", "ENSG00000188536.12", "ENSG00000244734.3")
mt_ids <- rownames(gene_info)[gene_info[,"Chrm"] == "chrM"]
malat1_ids <- rownames(gene_info)[gene_info[,"Name"] == "MALAT1"]

rm_ids <- c(rb_ids, rbr_ids, hb_ids, mt_ids)
keep_ids <- setdiff(rownames(gene_info), rm_ids)

sample <- meta[task, "Sample"]

#=========================================
# Output directories
#=========================================

dir_out <- paste0("data/processed/diem/", sample, "/"); create_dir(dir_out)
dir_plot <- paste0("exp/diem/plots/", sample, "/"); create_dir(dir_plot)
dir_fltr <- paste0("exp/diem/fltr_ids/", sample, "/"); create_dir(dir_fltr)
dir_de <- paste0("exp/diem/DE/", sample, "/"); create_dir(dir_de)

#=========================================
# Run DIEM
#=========================================

# Create output directories
cat(paste0("Running ", sample, "\n"))
cat("Reading expression data\n")

# Read in counts
dir_star <- "data/raw/STAR_output/"
dir_count <- paste0(dir_star, sample, "/Solo.out/GeneFull/raw/")
counts <- read_counts(dir_count)

colnames(counts) <- paste0(sample, ";", colnames(counts))

sce <- create_SCE(counts, name = sample) 

# Add MT% and MALAT1%
sce <- sce_feat_pct(sce, mt_ids, "pct.mt")
sce <- sce_feat_pct(sce, malat1_ids, "pct.malat1")
sce <- sce_feat_pct(sce, rb_ids, "pct.rb")
sce <- sce_feat_pct(sce, rbr_ids, "pct.rbr")
sce <- sce_feat_pct(sce, hb_ids, "pct.hb")

# Remove MT and ribosomal genes
rm_ids <- c(rb_ids, rbr_ids, hb_ids, mt_ids)
sce <- sce_rm_feat(sce, rm_ids)

sce <- set_debris_test_set(sce, min_counts = 100)
sce <- filter_genes(sce)

# Add percent spliced
fn <- paste0("data/processed/pct_splice/", sample, "/", sample, ".pct.splice.txt")
fn <- file.path("data/processed/pct_splice/", sample, paste0(sample, ".pct.splice.txt"))
pct_spl <- read.table(fn, header = TRUE, row.names = 1, stringsAsFactors = FALSE, 
                      sep = "\t")
di <- intersect(rownames(sce@test_data), rownames(pct_spl))
sce@test_data[di,"pct.spl"] <- pct_spl[di, "pct.splice"]

# Plots total counts ranked
p <- barcode_rank_plot(sce, ret = TRUE)
ggsave(paste0(dir_plot, sample, ".barcode_rank.jpeg"), 
       width = 3.5, height = 3.5)


require(gridExtra)
p1 <- plot_data(sce, feat_x = "total_counts", feat_y = "n_genes", 
                log_x = TRUE, log_y = TRUE, ret = TRUE)
p2 <- plot_data(sce, feat_x = "n_genes", feat_y = "pct.mt", 
                log_x = TRUE, ret = TRUE)
p3 <- plot_data(sce, feat_x = "n_genes", feat_y = "pct.malat1", 
                log_x = TRUE, ret = TRUE)
p4 <- plot_data(sce, feat_x = "pct.mt", feat_y = "pct.malat1", 
                log_x = FALSE, log_y = FALSE, ret = TRUE)

pdfname <- paste0(dir_plot, sample, ".cors.pdf")
jpgname <- paste0(dir_plot, sample, ".cors.jpeg")
pdf(pdfname, width = 7, height = 7)
grid.arrange(p1, p2, p3, p4, ncol = 2)
dev.off()
system(paste("convert", "-density", "300", pdfname, jpgname))

# Save DIEM object
cat("Saving\n")
saveRDS(sce, paste0(dir_out, sample, ".diem_sce.rds"))

