#!/u/local/apps/R/3.5.1/gcc-4.9.3_MKL-2018/bin/Rscript
#$ -S /u/local/apps/R/3.5.1/gcc-4.9.3_MKL-2018/bin/Rscript
#$ -cwd
#$ -j y
#$ -pe shared 8
#$ -l h_data=4000M,h_rt=8:00:00,highp
#$ -v QQAPP=openmp
#$ -M malvarez@mail
#  Notify at beginning and end of job
#$ -m n
#$ -r n
#$ -o CellTypeProp.R.log

setwd("../../")

library(Seurat)
library(ggplot2)
library(RColorBrewer)
source("scripts/ggplot_wrappers/R/barplots.R")

#=========================================
# Functions
#=========================================

create_dir <- function(p){
    dir.create(p, showWarnings=FALSE, recursive=TRUE)
}

#=========================================
#=========================================

dir_out <- "data/processed/d_sct_cca_all/merge/"; create_dir(dir_out)
dir_plot <- "exp/d_sct_cca_all/merge/plots/"; create_dir(dir_plot)

gene_info <- read.table("data/ref/gencode26/gencode.v26.annotation.txt",
                        header = TRUE,
                        row.names = 1,
                        stringsAsFactors = FALSE,
                        sep = "\t")

meta <- read.table("data/sample/samples.csv",
                   header = TRUE,
                   row.names = 1,
                   stringsAsFactors = FALSE,
                   sep = ",")

samples <- rownames(meta)

#=========================================
#=========================================

seur <- readRDS(paste0(dir_out, "seur.cca.rds"))

# replace "NASH with normal
seur$SampleID <- sub("NASH", "Normal", seur$SampleID)

# order factors from tumor to non-tumor
datf_s <- seur@meta.data[,c("SampleID", "Patient", "Tumor")]
datf_s[,"SampleID"] <- as.character(datf_s[,"SampleID"])
datf_s <- unique(datf_s)
o <- order(datf_s$Tumor, datf_s$Patient)

lv <- datf_s[o,"SampleID"]

seur$SampleID <- factor(seur$SampleID, levels = lv)


# Sample colors
sample_cols_df <- read.table("data/ref/TumorColors/sample_colors.1.txt", sep = ",", 
                             row.names = 1, stringsAsFactors = FALSE)
sample_cols <- sample_cols_df[,1]
names(sample_cols) <- rownames(sample_cols_df)

# Tumor colors
tumor_cols_df <- read.table("data/ref/TumorColors/tumor_colors.1.txt", sep = ",",
                            row.names = 1, stringsAsFactors = FALSE)
tumor_cols <- tumor_cols_df[,1]
names(tumor_cols) <- rownames(tumor_cols_df)

# Patient colors
patient_cols_df <- read.table("data/ref/TumorColors/patient_colors.1.txt", sep = ",", 
                              row.names = 1, stringsAsFactors = FALSE)
patient_cols <- patient_cols_df[,1]
names(patient_cols) <- rownames(patient_cols_df)


#=========================================
#=========================================

ct <- theme(plot.title = element_text(hjust = 0.5), 
      axis.text.x = element_text(angle = 30, hjust = 0.9), 
      axis.ticks = element_blank(), 
      panel.border = element_blank(),
      panel.grid = element_blank(), 
      axis.title.y = element_blank())

#=========================================
# Sample cell type frequency
#=========================================

SampleCT <- table(seur$CellType, seur$SampleID)
SampleCT <- sweep(SampleCT, 2, colSums(SampleCT), "/")
# SampleCT <- sweep(SampleCT, 1, rowSums(SampleCT), "/")
SampleCT <- SampleCT[,lv]

sampleCTm <- reshape2::melt(SampleCT)
colnames(sampleCTm) <- c("CellType", "Sample", "Prop")
# sampleCTm[,"Sample"] <- factor(sampleCTm[,"Sample"], levels = lv)

source("scripts/ggplot_wrappers/R/barplots.R")
p <- pct_stack_barplt(sampleCTm, x = "CellType", y = "Prop", fill = "Sample")
p <- p + xlab("Cell-type") + ylab("Proportion") 
p <- p + coord_flip()
p <- p + theme(panel.grid = element_blank(), 
               axis.text.x = element_text(angle = 0, hjust = 0.5))
p <- p + scale_fill_manual(values = sample_cols)
outfn <- paste0(dir_plot, "Sample.CellType.prop.pdf")
ggsave(outfn, width = 5, height = 10)


#=========================================
# Tumor status cell type frequency
#=========================================

TumorCT <- table(seur$CellType, seur$Tumor)
TumorCT <- sweep(TumorCT, 2, colSums(TumorCT), "/")
# TumorCT <- sweep(TumorCT, 1, rowSums(TumorCT), "/")

TumorCTm <- reshape2::melt(TumorCT)
colnames(TumorCTm) <- c("CellType", "Tumor", "Prop")

source("scripts/ggplot_wrappers/R/barplots.R")
p <- pct_stack_barplt(TumorCTm, x = "CellType", y = "Prop", fill = "Tumor")
p <- p + coord_flip()
p <- p + theme(panel.grid = element_blank(), 
               axis.text.x = element_text(angle = 0, hjust = 0.5))
p <- p + scale_fill_manual(values = tumor_cols)
outfn <- paste0(dir_plot, "Tumor.CellType.prop.pdf")
ggsave(outfn, width = 5, height = 10)


#=========================================
# patient cell type frequency
#=========================================

SubjCT <- table(seur$CellType, seur$Patient)
SubjCT <- sweep(SubjCT, 2, colSums(SubjCT), "/")
# SubjCT <- sweep(SubjCT, 1, rowSums(SubjCT), "/")

SubjCTm <- reshape2::melt(SubjCT)
colnames(SubjCTm) <- c("CellType", "Patient", "Prop")
SubjCTm$Patient <- factor(as.character(SubjCTm$Patient))

source("scripts/ggplot_wrappers/R/barplots.R")
p <- pct_stack_barplt(SubjCTm, x = "CellType", y = "Prop", fill = "Patient")
p <- p + coord_flip()
p <- p + theme(panel.grid = element_blank(), 
               axis.text.x = element_text(angle = 0, hjust = 0.5))
p <- p + scale_fill_manual(values = patient_cols)
outfn <- paste0(dir_plot, "Patient.CellType.prop.pdf")
ggsave(outfn, width = 5, height = 10)


#=========================================
# Test for proportions
#=========================================

# divide by sample num
SampleCT <- table(seur$CellType, seur$SampleID)
SampleCT <- sweep(SampleCT, 2, colSums(SampleCT), "/")
SampleCT <- SampleCT[,lv]

tumors <- grep("Tumor", colnames(SampleCT), value = TRUE)
nontumors <- grep("Tumor", colnames(SampleCT), value = TRUE, invert = TRUE)

datf <- matrix(nrow = nrow(SampleCT), ncol = 2)
rownames(datf) <- rownames(SampleCT)
colnames(datf) <- c("T", "p")

for (i in rownames(SampleCT)){
    tr <- t.test(SampleCT[i,nontumors], SampleCT[i, tumors], paired = TRUE)
    datf[i, "T"] <- tr$statistic
    datf[i, "p"] <- tr$p.value
}

