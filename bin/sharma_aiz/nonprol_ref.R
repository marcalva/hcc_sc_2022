
# create a reference for the non-proliferating cells to 
# use for classifying the proliferating cells

setwd("../../")

library(Seurat)
suppressPackageStartupMessages(library(SingleR))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(scran))

#=========================================
# Functions
#=========================================

create_dir <- function(p){
    dir.create(p, showWarnings=FALSE, recursive=TRUE)
}

#=========================================
#=========================================

# Gene data
gene_info <- read.table("data/ref/gencode26/gencode.v26.annotation.txt", 
                        header = TRUE, 
                        row.names = 1, 
                        stringsAsFactors = FALSE, 
                        sep = "\t")
gene_info[,"Name"] <- make.unique(gene_info[,"Name"])
symb2ens <- rownames(gene_info)
names(symb2ens) <- gene_info[,"Name"]

fn <- "data/processed/sharma_aiz/liver.int_rand.rds"
integrated <- readRDS(fn)

# Set directories
dir_out <- "data/processed/sharma_aiz/"
dir.create(dir_out, showWarnings=FALSE, recursive=TRUE)

dir_exp <- "exp/sharma_aiz/clust_prolif/";
dir.create(dir_exp, showWarnings=FALSE, recursive=TRUE)

#=========================================
# create reference
#=========================================

md <- integrated@meta.data

cl_ids <- c("cell_type_main", "cell_type_fine")

for (cl_id in cl_ids){
    message(cl_id)

    cl_id <- "cell_type_main"
    md[,cl_id] <- as.character(md[,cl_id])

    k_ref <- rownames(md)[md[,cl_id] != "Prol"]
    k_prol <- rownames(md)[md[,cl_id] == "Prol"]

    ex <- integrated@assays$RNA@data
    ex_ref <- ex[,k_ref]
    ex_prol <- ex[,k_prol]

    # get de genes
    ref_tt <- pairwiseTTests(x = ex_ref,
                             groups = md[k_ref,cl_id],
                             direction = "up")

    ref_top <- getTopMarkers(ref_tt$statistics, ref_tt$pairs, n = 100)

    set.seed(1, kind = 'Mersenne-Twister')
    ref_trnd <- trainSingleR(ref = ex_ref, labels = md[k_ref,cl_id], genes = ref_top, 
                             aggr.ref = TRUE)

    outfn <- paste0(dir_exp, cl_id, ".nonprol_ref.rds")
    saveRDS(ref_trnd, outfn)
}

message("done")


