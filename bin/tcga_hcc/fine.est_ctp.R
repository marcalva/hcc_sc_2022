
# estimate cell type proportions in the TCGA data
# using the markers

setwd("../../")

library(edgeR)
library(Biobase)
library(plyr)
source("scripts/sc_func.R")
# source("scripts/reference_free.R")
library(BisqueRNA)

# read in data
out_dir <- "data/processed/tcga_hcc/expr/"

# expression data
fn <- paste0(out_dir, "tcga.lihc.TMM.log.rin.rds")
tmm <- readRDS(fn)

fn <- "data/processed/tcga_hcc/expr/tcga.gencodev26.rds"
gencode <- readRDS(fn)

# pheno data
out_dir <- "data/processed/tcga_hcc/sample/"

fn <- paste0(out_dir, "tcga.clinical.rds")
clinical <- readRDS(fn)
fn <- paste0(out_dir, "tcga.sample_sheet.rds")
sample_sheet <- readRDS(fn)

min_gene <- 20
max_gene <- 300

cl_id <- "cell_type_fine"

# marker data
fn <- paste0("exp/sharma_aiz/markers/markers.", cl_id, ".txt")
markers.celltype <- read.table(fn, header=TRUE, stringsAsFactors=FALSE)
markers.celltype[,"gene"] <- sapply(markers.celltype[,"gene"], 
                                    function(s){
                                        strsplit(s, "\\.")[[1]][1] })

k <- markers.celltype[,"gene"] %in% rownames(tmm)
markers.celltype <- markers.celltype[k,]
ct.count <- table(markers.celltype[,"cluster"])
ctk <- names(ct.count)[ct.count >= min_gene]
k <- markers.celltype[,"cluster"] %in% ctk
markers.celltype <- markers.celltype[k,]

ct_id <- "cluster"

# run decomposition
tmm.eset <- ExpressionSet(tmm)

weighted <- FALSE

# cell types
tmm.ct.md <- MarkerBasedDecomposition(bulk.eset = tmm.eset, 
                                     markers = markers.celltype, 
                                     ct_col = ct_id, 
                                     gene_col = "gene", 
                                     weighted = weighted, 
                                     unique_markers = FALSE,
                                     min_gene = min_gene, 
                                     max_gene = max_gene)
tmm.ct.mdp <- tmm.ct.md$bulk.props

# get PC1 of all markers for each cell type
markers <- markers.celltype[,"gene"]
cell_types <- markers.celltype[,ct_id]
cell_types.u <- unique(cell_types)
names(cell_types.u) <- cell_types.u
ct_pcs <- lapply(cell_types.u, function(ct){
                 g <- markers[cell_types == ct]
                 pcr <- prcomp(t(tmm[g,]), retx=TRUE, scale.=TRUE, center=TRUE)
                 pcrret <- pcr$x[,1]
                 cors <- cor(t(tmm[g,]), pcrret)
                 if (mean(cors) < 0) pcrret <- -1 * pcrret
                 return(pcrret) })
ct_pcs <- do.call(cbind, ct_pcs)
# this shows they are highly correlated
sort(diag(cor(ct_pcs, t(tmm.ct.mdp)[rownames(ct_pcs),colnames(ct_pcs)])))

# correlate markers with cell-type proportions
ct_m_cors_l <- list()
for (ct in names(tmm.ct.md$genes.used)){
    ens <- tmm.ct.md$genes.used[[ct]]
    s <- colnames(tmm)
    p_cors <- cor(t(tmm[ens,s]), t(tmm.ct.md$bulk.props[ct,s,drop=FALSE]))
    ct_m_cors_l[[ct]] <- p_cors
}

# output results
out_dir <- "data/processed/tcga_hcc/ctp/"
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

fn <- paste0(out_dir, "tcga.TMM.", cl_id, ".decomp.rds")
saveRDS(tmm.ct.md, fn)

# export markers used for decomp
cts <- lapply(names(tmm.ct.md$genes.used), function(ct){
              rep(ct, times = length(tmm.ct.md$genes.used[[ct]]))
                 })
cts <- do.call(c, cts)
symbs <- lapply(names(tmm.ct.md$genes.used), function(ct){
                gencode[tmm.ct.md$genes.used[[ct]], "Name"]
                 })
symbs <- do.call(c, symbs)
ens <- lapply(names(tmm.ct.md$genes.used), function(ct){
              tmm.ct.md$genes.used[[ct]] })
ens <- do.call(c, ens)
cors <- lapply(names(tmm.ct.md$genes.used), function(ct){
               ct_m_cors_l[[ct]][tmm.ct.md$genes.used[[ct]],ct] })
cors <- do.call(c, cors)

mrk_datf <- data.frame("CellType" = cts, 
                       "Ensembl" = ens, 
                       "Symbol" = symbs, 
                       "R" = cors)

dir_exp <- paste0("exp/tcga_hcc/ctp.", cl_id, "/")
dir.create(dir_exp, showWarnings = FALSE, recursive = TRUE)
out_fn <- paste0(dir_exp, cl_id, ".marker_genes.txt")
write.table(mrk_datf, out_fn, row.names = FALSE, col.names = TRUE, 
            quote = FALSE, sep = '\t')

# Plot
library(ggplot2)
library(reshape2)
plot_dir <- paste0("exp/tcga_hcc/ctp.", cl_id, "/ctp_cor/")
dir.create(plot_dir, showWarnings=FALSE, recursive=TRUE)

cell_types <- names(tmm.ct.md$genes.used)

for (ct in cell_types){
    cors.a <- cor(t(tmm[tmm.ct.md$genes.used[[ct]],]), tmm.ct.md$bulk.props[ct,])
    rownames(cors.a) <- gencode[rownames(cors.a), "Name"]
    colnames(cors.a) <- ct
    cors.dfm <- melt(cors.a)
    p <- ggplot(cors.dfm, aes(x = Var2, y = Var1, fill = value)) + 
        geom_tile() + 
        scale_fill_distiller(palette = "RdBu", 
                             limits = c(-1, 1),
                             name = substitute(paste(italic(R)))) + 
        labs(x = "Markers", y = ct) + 
        theme_bw()
    outfn <- paste0(plot_dir, ct, ".marker.ctp.cor.pdf")
    ggsave(outfn, width = 3, height = 6)
}

# plot co-expression
plot_dir <- paste0("exp/tcga_hcc/ctp.", cl_id, "/marker_coexpr/")
dir.create(plot_dir, showWarnings=FALSE, recursive=TRUE)

cors.l <- lapply(cell_types, function(ct){
                 gs <- tmm.ct.md[["genes.used"]][[ct]]
                 cors.df <- cor(t(tmm[gs,]))
                 cors.dfm <- reshape2::melt(cors.df)
                 cors.dfm[,"CellType"] <- ct
                 return(cors.dfm) })
cors <- do.call(rbind, cors.l)

facet_ncol = 2

outfn <- paste0(plot_dir, "connect.cell_type.marker_coexpr.pdf")
p <- ggplot(cors, aes(x = Var1, y = Var2, fill = value)) + 
    geom_tile(aes(fill = value)) + 
    scale_fill_distiller(palette = "RdBu", 
                         limits = c(-1, 1),
                         name = substitute(paste(italic(R)))) + 
    labs(x = "Markers", y = "Markers") + 
    ggtitle("Cell-type marker co-expression (TCGA)") + 
    facet_wrap(~CellType, scales = "free", ncol = facet_ncol) + 
    theme_bw() + 
    guides(fill = guide_colourbar(title.position = "left")) +  
    theme(plot.title = element_text(size = 14, hjust = 0.5), 
          axis.title = element_text(size = 12),
          axis.text.x = element_blank(), 
          axis.text.y = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.ticks.y = element_blank(), 
          strip.text = element_text(size = 8), 
          legend.title = element_text(size = 12, hjust =1, vjust = 0.8),
          strip.background = element_rect(fill = "white"),
          legend.position = "bottom")

ggsave(outfn, width = 7, height = 9)


