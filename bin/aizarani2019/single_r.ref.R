
# predict cell-types using annotators

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

seur_fn <- "data/processed/aizarani2019/seur.rds"
seur <- readRDS(seur_fn)

dir_exp <- "exp/aizarani2019/sct_clust/SingleR/"
create_dir(dir_exp)

#=========================================
# Normalize Aizarani counts
#=========================================

DefaultAssay(seur) <- "RNA"
seur <- NormalizeData(seur, scale.factor=1000)
lcounts <- seur@assays$RNA@data
rownames(lcounts) <- make.unique(gene_info[rownames(lcounts), "Name"])

clusts <- seur$SCT_snn_res.0.5
clusts_u <- factor(sort(unique(clusts)))

ct_lcounts <- collapse(seur@assays$RNA@counts, 
                       clusts, scale_size=1e3)
rownames(ct_lcounts) <- make.unique(gene_info[rownames(ct_lcounts), "Name"])


#=========================================
# SingleR ENCODE
#=========================================

# Encode
enc <- celldex::BlueprintEncodeData()
enc_counts <- enc@assays@data$logcounts
gi <- intersect(rownames(lcounts), rownames(enc_counts))

# get de genes
enc_fine_tt <- pairwiseTTests(x = enc_counts[gi,],
                              groups = enc@colData$label.fine,
                              direction = "up")
enc_fine_top <- getTopMarkers(enc_fine_tt$statistics, enc_fine_tt$pairs, n = 100)

enc_main_tt <- pairwiseTTests(x = enc_counts[gi,],
                              groups = enc@colData$label.main,
                              direction = "up")
enc_main_top <- getTopMarkers(enc_main_tt$statistics, enc_main_tt$pairs, n = 100)

trnd_enc_fine <- trainSingleR(ref = enc_counts[gi,],
                              labels = enc@colData$label.fine,
                              genes = enc_fine_top)
trnd_enc_main <- trainSingleR(ref = enc_counts[gi,],
                              labels = enc@colData$label.main,
                              genes = enc_main_top)

clss_enc_fine <- classifySingleR(test = lcounts[gi,],
                                 trained = trnd_enc_fine)
clss_enc_main <- classifySingleR(test = lcounts[gi,],
                                 trained = trnd_enc_main)

clss_enc_fine_cl <- classifySingleR(test = ct_lcounts[gi,],
                                    trained = trnd_enc_fine)
clss_enc_main_cl <- classifySingleR(test = ct_lcounts[gi,],
                                    trained = trnd_enc_main)

#=========================================
# SingleR HPA
#=========================================

hpa <- celldex::HumanPrimaryCellAtlasData()
hpa_counts <- hpa@assays@data$logcounts
gi <- intersect(rownames(lcounts), rownames(hpa_counts))

# get de genes
hpa_fine_tt <- pairwiseTTests(x = hpa_counts[gi,],
                              groups = hpa@colData$label.fine,
                              direction = "up")
hpa_fine_top <- getTopMarkers(hpa_fine_tt$statistics, hpa_fine_tt$pairs, n = 100)

hpa_main_tt <- pairwiseTTests(x = hpa_counts[gi,],
                              groups = hpa@colData$label.main,
                              direction = "up")
hpa_main_top <- getTopMarkers(hpa_main_tt$statistics, hpa_main_tt$pairs, n = 100)

trnd_hpa_fine <- trainSingleR(ref = hpa_counts[gi,],
                              labels = hpa@colData$label.fine,
                              genes = hpa_fine_top)
trnd_hpa_main <- trainSingleR(ref = hpa_counts[gi,],
                              labels = hpa@colData$label.main,
                              genes = hpa_main_top)

clss_hpa_fine <- classifySingleR(test = lcounts[gi,],
                                 trained = trnd_hpa_fine)
clss_hpa_main <- classifySingleR(test = lcounts[gi,],
                                 trained = trnd_hpa_main)

clss_hpa_fine_cl <- classifySingleR(test = ct_lcounts[gi,],
                                    trained = trnd_hpa_fine)
clss_hpa_main_cl <- classifySingleR(test = ct_lcounts[gi,],
                                    trained = trnd_hpa_main)

#=========================================
# SingleR Novershtern
#=========================================

nov <- celldex::NovershternHematopoieticData()
nov_counts <- nov@assays@data$logcounts
gi <- intersect(rownames(lcounts), rownames(nov_counts))

# get de genes
nov_fine_tt <- pairwiseTTests(x = nov_counts[gi,],
                              groups = nov@colData$label.fine,
                              direction = "up")
nov_fine_top <- getTopMarkers(nov_fine_tt$statistics, nov_fine_tt$pairs, n = 100)

nov_main_tt <- pairwiseTTests(x = nov_counts[gi,],
                              groups = nov@colData$label.main,
                              direction = "up")
nov_main_top <- getTopMarkers(nov_main_tt$statistics, nov_main_tt$pairs, n = 100)

trnd_nov_fine <- trainSingleR(ref = nov_counts[gi,],
                              labels = nov@colData$label.fine,
                              genes = nov_fine_top)
trnd_nov_main <- trainSingleR(ref = nov_counts[gi,],
                              labels = nov@colData$label.main,
                              genes = nov_main_top)

clss_nov_fine <- classifySingleR(test = lcounts[gi,],
                                 trained = trnd_nov_fine)
clss_nov_main <- classifySingleR(test = lcounts[gi,],
                                 trained = trnd_nov_main)

clss_nov_fine_cl <- classifySingleR(test = ct_lcounts[gi,],
                                    trained = trnd_nov_fine)
clss_nov_main_cl <- classifySingleR(test = ct_lcounts[gi,],
                                    trained = trnd_nov_main)

ctmap_l <- list("enc_fine" = clss_enc_fine, 
                "enc_main" = clss_enc_main,
                "hpa_fine" = clss_hpa_fine,
                "hpa_main" = clss_hpa_main, 
                "nov_fine" = clss_nov_fine,
                "nov_main" = clss_nov_main)

outfn <- file.path(dir_exp, "singleR.ref.rds")
saveRDS(ctmap_l, outfn)

