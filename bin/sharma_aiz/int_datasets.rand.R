
# Integrate Sharma, Aizarani, and NASH-HCC liver data
# run CCA on sctransform normalized data
# each independent data set is 
#  the biopsy from NASH-HCC (due to heterogeneity)
#  the entire Aizarani et al data set
#  each patient from Sharma et al data set
# the reference samples include the Aizarani healthy data 
# and 10 samples selected randomly from the NASH-HCC and 
# Sharma data 

setwd("../../")

library(Seurat)
library(ggplot2)
library(scales)
library(gridExtra)

#=========================================
# Functions
#=========================================

create_dir <- function(p){
    dir.create(p, showWarnings=FALSE, recursive=TRUE)
}

#=========================================
#=========================================

n_threads <- 8

# Gene data
gene_info <- read.table("data/ref/gencode26/gencode.v26.annotation.txt", 
                        header = TRUE, 
                        row.names = 1, 
                        stringsAsFactors = FALSE, 
                        sep = "\t")
gene_info[,"Name"] <- make.unique(gene_info[,"Name"])
symb2ens <- rownames(gene_info)
names(symb2ens) <- gene_info[,"Name"]

fn <- "data/processed/sharma2019/sharma.seur.rds"
sharma <- readRDS(fn)
sharma_l <- SplitObject(sharma, split.by = "patientno")
rm(sharma)

fn <- "data/processed/d_sct_cca_all/merge/seur.cca.rds"
nash_hcc <- readRDS(fn)
nash_hcc_l <- SplitObject(nash_hcc, split.by = "orig.ident")
rm(nash_hcc)

fn <- "data/processed/aizarani2019/seur.rds"
aiz <- readRDS(fn)
DefaultAssay(aiz) <- "SCT"

# Set directories
dir_out <- "data/processed/sharma_aiz/"
dir.create(dir_out, showWarnings=FALSE, recursive=TRUE)

dir_exp <- "exp/sharma_aiz/";
dir.create(dir_exp, showWarnings=FALSE, recursive=TRUE)

#=========================================
# SCT then CCA integration
#=========================================

for (i in 1:length(nash_hcc_l)){
    nash_hcc_l[[i]]$source <- "nash_hcc"
}

for (i in 1:length(sharma_l)){
    sharma_l[[i]]$source <- "Sharma"
}

aiz$source <- "Aizarani"

s_l <- c(nash_hcc_l, sharma_l, aiz)
names(s_l)[length(s_l)] <- "Aiz"
rm(nash_hcc_l, sharma_l, aiz)
g <- gc()

# subset features to intersection
g_l <- lapply(s_l, function(s) rownames(s@assays$RNA@counts))
g_i <- Reduce(intersect, g_l)

# get cell cycle genes
s.genes <- cc.genes$s.genes
s.genes <- symb2ens[s.genes]
s.genes <- s.genes[s.genes %in% g_i]
g2m.genes <- cc.genes$g2m.genes
g2m.genes <- symb2ens[g2m.genes]
g2m.genes <- g2m.genes[g2m.genes %in% g_i]

# cell cycle scoring and SCTransform
for (i in 1:length(s_l)){
    message(i)
    rc <- s_l[[i]]@assays$RNA@counts
    md <- s_l[[i]]@meta.data
    rc <- rc[g_i,]
    s <- CreateSeuratObject(counts = rc, meta.data = md)

    DefaultAssay(s) <- "RNA"
    s <- NormalizeData(s, scale.factor = 1e3)
    s <- CellCycleScoring(s, 
                          s.features = s.genes, 
                          g2m.features = g2m.genes, 
                          set.ident = FALSE)
    
    s <- SCTransform(s, variable.features.n = 3000, verbose=FALSE)
    s_l[[i]] <- s
}

# samples for reference picked randomly
set.seed(1, kind = "Mersenne-Twister")
refri <- sample(1:(length(s_l)-1), size=10, replace=FALSE)
refri <- c(refri, length(s_l))
names(s_l)[refri]

features <- SelectIntegrationFeatures(object.list = s_l, nfeatures = 3000)
s_l <- PrepSCTIntegration(object.list = s_l, anchor.features = features)

# library(future)
# maxsize <- 8000 * 1024 ^ 2
# options(future.globals.maxSize = maxsize)
# plan("multicore", workers = 2)
anchors <- FindIntegrationAnchors(object.list = s_l, 
                                  reference = refri, 
                                  dims = 1:30, 
                                  normalization.method = "SCT", 
                                  reduction = "cca", 
                                  anchor.features = features)
integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT", dims = 1:30)

outfn <- file.path(dir_out, "liver.int_rand.rds")
saveRDS(integrated, outfn)

rm(s_l)
g <- gc()

integrated <- ScaleData(integrated, verbose = FALSE)
integrated <- RunPCA(integrated, npcs = 50, verbose = FALSE)
integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:30)
integrated <- FindNeighbors(integrated, reduction = "pca", dims = 1:30)
integrated <- FindClusters(integrated, resolution = c(0.5, 0.8, 1))

outfn <- file.path(dir_out, "liver.int_rand.rds")
saveRDS(integrated, outfn)

