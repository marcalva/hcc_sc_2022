
# predict cell-types using SingleR and Aizarani reference

setwd("../../")

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SingleR))
suppressPackageStartupMessages(library(celldex))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(scran))

#=========================================
# Functions
#=========================================

create_dir <- function(x){
    dir.create(x, showWarnings = FALSE, recursive = TRUE)
}

# Collapse expression along cell types
# x is a sparse matrix of counts, the raw gene-barcode matrix
# identities is a vector of cell types
# If counts == TRUE, return the raw counts
collapse = function(x, 
                    identities, 
                    counts=FALSE, 
                    scale_size=1e4, 
                    scale_cells=TRUE,
                    logt=TRUE){
    if (is.null(attr(class(x), "package")) || 
        attr(class(x), "package") != "Matrix"){
        x = Matrix::Matrix(x)
    }
    ngenes = nrow(x)
    if (class(identities) != "factor"){
        identities = factor(identities)
    }
    celltypes = levels(identities)
    ncells = table(identities)[celltypes]
    panel = matrix(0, nrow=ngenes, ncol=length(celltypes))
    rownames(panel) = rownames(x)
    colnames(panel) = celltypes
    for (celltype in celltypes){
        panel[,celltype] = Matrix::rowSums(x[,identities == celltype,drop=FALSE])
    }

    if (counts) return(panel)

    panel <- sweep(panel, 2, colSums(panel), "/")
    
    if (scale_cells){
        panel = scale_size*panel
    }
    
    if (logt){
        panel = log1p(panel)
    }

    return(panel)
}

#=========================================
#=========================================

# gene annotations
gene_info <- read.table("data/ref/gencode26/gencode.v26.annotation.txt",
                        header = TRUE,
                        row.names = 1,
                        stringsAsFactors = FALSE,
                        sep = "\t")

sharma_fn <- "data/processed/sharma2019/sharma.seur.rds"
sharma <- readRDS(sharma_fn)

aiz_fn <- "data/processed/aizarani2019/seur.rds"
aiz <- readRDS(aiz_fn)

dir_exp <- "exp/sharma2019/singleR/"
create_dir(dir_exp)

#=========================================
# Automated annotator
#=========================================

clusts <- sharma$SCT_snn_res.0.8
clusts_u <- factor(sort(unique(clusts)))

DefaultAssay(aiz) <- "RNA"
aiz <- NormalizeData(aiz, scale.factor=1000)
aiz_lcounts <- aiz@assays$RNA@data

ct_de_pairs <- pairwiseTTests(x = aiz_lcounts, 
                              groups = aiz$CellType, 
                              direction = "up")
ct_de_top <- getTopMarkers(ct_de_pairs$statistics, ct_de_pairs$pairs, n = 100)

mj_de_pairs <- pairwiseTTests(x = aiz_lcounts, 
                              groups = aiz$Major, 
                              direction = "up")
mj_de_top <- getTopMarkers(mj_de_pairs$statistics, mj_de_pairs$pairs, n = 100)

DefaultAssay(sharma) <- "RNA"
sharma <- NormalizeData(sharma, scale.factor=1000)
sharma_lcounts <- sharma@assays$RNA@data

sharma_ct_lcounts <- collapse(sharma@assays$RNA@counts,
                             clusts, scale_size=1e3)

k <- intersect(rownames(sharma_lcounts), rownames(aiz_lcounts))
sharma_lcounts <- sharma_lcounts[k,]
aiz_lcounts <- aiz_lcounts[k,]
sharma_ct_lcounts <- sharma_ct_lcounts[k,]

# trait reference for singleR
trnd_ct <- trainSingleR(ref = aiz_lcounts,
                        labels = aiz$CellType, 
                        genes = ct_de_top)
trnd_mj <- trainSingleR(ref = aiz_lcounts,
                        labels = aiz$Major, 
                        genes = mj_de_top)

aiz_clss_ct <- classifySingleR(test = sharma_lcounts, 
                               trained = trnd_ct)
aiz_clss_mj <- classifySingleR(test = sharma_lcounts, 
                               trained = trnd_mj)

aiz_clss_ct_cl <- classifySingleR(test = sharma_ct_lcounts, 
                                  trained = trnd_ct)
aiz_clss_mj_cl <- classifySingleR(test = sharma_ct_lcounts,
                                  trained = trnd_mj)

ctmap_l <- list("aiz_fine" = aiz_clss_ct, 
                "aiz_crse" = aiz_clss_mj, 
                "aiz_fine_cl" = aiz_clss_ct_cl, 
                "aiz_crse_cl" = aiz_clss_mj_cl)

outfn <- file.path(dir_exp, "singleR.rds")
saveRDS(ctmap_l, outfn)

