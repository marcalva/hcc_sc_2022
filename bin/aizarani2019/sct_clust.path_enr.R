
# Run pathway enrichment
# Keep all pathways for evaluation
# expressed genes are with counts >= 10

setwd("../../")

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(qvalue))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(AnnotationHub))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(ReactomePA))

#=========================================
# Functions
#=========================================

create_dir <- function(x){
    dir.create(x, showWarnings = FALSE, recursive = TRUE)
}

#=========================================
#=========================================

gene_info <- read.table("data/ref/gencode26/gencode.v26.annotation.txt", 
                        header = TRUE, 
                        row.names = 1, 
                        stringsAsFactors = FALSE, 
                        sep = "\t")
gene_info[,"Ensembl"] <- rownames(gene_info)

# Read in Seurat objects
seur_fn <- "data/processed/aizarani2019/seur.rds"
seur <- readRDS(seur_fn)

fn <- "exp/aizarani2019/sct_clust/markers/mrk.SCT_snn_res.0.5.txt"
mk <- read.table(fn, header=TRUE)

lfc_fn <- "exp/aizarani2019/sct_clust/markers/log_fc.SCT_snn_res.0.5.txt"
lfc <- read.table(lfc_fn, header=TRUE, check.names=FALSE, row.names=NULL)
k <- !is.na(lfc[,2])
lfc <- lfc[k,]
rownames(lfc) <- lfc[,"row.names"]

dir_exp <- "exp/aizarani2019/sct_clust/path_enr/"
create_dir(dir_exp)

rc <- seur@assays$RNA@counts
md <- seur@meta.data
cts <- sort(as.character(unique(md[,"SCT_snn_res.0.5"])))

rm(seur)
gc()

#=========================================
# pathway enrichment
#=========================================

parse_rto <- function(s){
    if (length(s) == 0) return(NULL)
    sapply(strsplit(s, "/"), function(x) {x = as.numeric(x); x[1] / x[2]} )
}

# map between gene names and entrez IDs
x <- org.Hs.egSYMBOL2EG
keys(x) <- mappedkeys(x)
lmap <- unlist(as.list(x))

ct2t <- unique(mk[,"cluster"])

enr_all <- list("GO" = list(), 
                "KEGG" = list(), 
                "Reactome" = list())

ens_e <- rownames(rc)
exprsd <- rowSums(rc) >= 10
ens_e <- ens_e[exprsd]
genes_e <- gene_info[ens_e, "Name"]
ez_e <- lmap[genes_e]
ez_e <- ez_e[!is.na(ez_e)]

# subset markers to those present in org.Hs
k <- mk$Name %in% names(lmap)
mk <- mk[k,,drop=FALSE] # 136 of 3316 removed

pvcut <- 2
qvcut <- 2

for (ct in ct2t){
    message(ct)

    # get list of DE genes
    mks <- mk[mk[,"cluster"] == ct,]

    genes <- mks[,"Name"]
    ezs <- lmap[genes]
    ezs <- ezs[!is.na(ezs)]

    go_ret <- enrichGO(gene = ezs, ont = "ALL", universe = ez_e, 
                       pvalueCutoff = pvcut, qvalueCutoff = qvcut,
                       OrgDb = org.Hs.eg.db, readable=TRUE)
    go_df <- as.data.frame(go_ret)
    enr_all[["GO"]][[as.character(ct)]] <- go_df

    kegg_ret <- enrichKEGG(gene = ezs, organism = "hsa", universe = ez_e, 
                           pvalueCutoff = pvcut, qvalueCutoff = qvcut)
    if (!is.null(kegg_ret)) kegg_ret <- setReadable(kegg_ret, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
    kegg_df <- as.data.frame(kegg_ret)
    enr_all[["KEGG"]][[as.character(ct)]] <- kegg_df

    reac_ret <- enrichPathway(gene = ezs, universe = ez_e, 
                              pvalueCutoff = pvcut, qvalueCutoff = qvcut,
                              readable=TRUE)
    reac_df <- as.data.frame(reac_ret)
    enr_all[["Reactome"]][[as.character(ct)]] <- reac_df


    ctc <- as.character(ct)
    for (anno in names(enr_all)){
        if (is.null(enr_all[[anno]][[ctc]]) || nrow(enr_all[[anno]][[as.character(ctc)]]) == 0) next
            enr_all[[anno]][[ctc]][,"FoldEnrich"] <- 
                parse_rto(enr_all[[anno]][[ctc]][,"GeneRatio"]) / parse_rto(enr_all[[anno]][[ctc]][,"BgRatio"])
            enr_all[[anno]][[ctc]][,c("CellType")] <- factor(ct)
    }
}

outfn <- paste0(dir_exp, "markers.res.0.5.path_enr.rds")
saveRDS(enr_all, outfn)

# save each in text file
for (a in names(enr_all)){
    enr_a_l <- enr_all[[a]]
    enr_a <- do.call(rbind, enr_a_l)
    outfn <- paste0(dir_exp, "markers.res.0.5.", a, ".txt")
    write.table(enr_a, outfn, row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
}

# format of output
# 1st list contains GO, KEGG, and Reactome annotations
# 3rd list contains cell types
# cell types contain data frames of clusterProfiler output.
#   GeneRatio is num. of genes in pathway (in list) divided by total genes (in list).
#   BgRatio is num. of genes in pathway (in list+background) divided by total genes (in list+background).
#   FoldEnrich is GeneRatio/BgRatio

#=======================================
# gsea analysis
#=======================================

gse_all <- list("GO" = list(), 
                "KEGG" = list(), 
                "Reactome" = list())

nps <- 10e3
eps <- 1e-50

for (ct in cts){
    message(ct)

    ct_lfc <- lfc[,ct]
    names(ct_lfc) <- lfc[,"Name"]
    ct_lfc <- sort(ct_lfc, decreasing=TRUE)
    ct_lfc <- ct_lfc[names(ct_lfc) %in% names(ez_e)]
    names(ct_lfc) <- as.character(lmap[names(ct_lfc)])

    go_ret <- gseGO(geneList = ct_lfc, pvalueCutoff=2, seed=1, 
                    ont = "ALL", OrgDb = org.Hs.eg.db, eps=eps, nPermSimple = nps)
    if (!is.null(go_ret)) go_ret <- setReadable(go_ret, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
    go_df <- as.data.frame(go_ret)
    gse_all[["GO"]][[as.character(ct)]] <- go_df

    kegg_ret <- gseKEGG(geneList = ct_lfc, pvalueCutoff=2, seed = 1, eps=eps, nPermSimple = nps)
    if (!is.null(kegg_ret)) kegg_ret <- setReadable(kegg_ret, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
    kegg_df <- as.data.frame(kegg_ret)
    gse_all[["KEGG"]][[as.character(ct)]] <- kegg_df

    reac_ret <- gsePathway(geneList = ct_lfc, pvalueCutoff=2, seed=1, eps=eps, nPermSimple = nps)
    if (!is.null(reac_ret)) reac_ret <- setReadable(reac_ret, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
    reac_df <- as.data.frame(reac_ret)
    gse_all[["Reactome"]][[as.character(ct)]] <- reac_df
                           
    ctc <- as.character(ct)
    for (anno in names(gse_all)){
        if (is.null(gse_all[[anno]][[ctc]]) || nrow(gse_all[[anno]][[as.character(ctc)]]) == 0) next
        gse_all[[anno]][[ctc]][,c("CellType")] <- factor(ct)
    }
}

outfn <- paste0(dir_exp, "markers.res.0.5.path_gse.rds")
saveRDS(gse_all, outfn)

# save each in text file
for (a in names(gse_all)){
    gse_a_l <- gse_all[[a]]
    gse_a <- do.call(rbind, gse_a_l)
    outfn <- paste0(dir_exp, "markers.res.0.5.gse.", a, ".txt")
    write.table(gse_a, outfn, row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
}

