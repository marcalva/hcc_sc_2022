
# estimate cell type proportions in the TCGA data
# using the markers

setwd("../../")

library(edgeR)
library(Biobase)
library(plyr)
source("scripts/sc_func.R")
# source("scripts/reference_free.R")
library(BisqueRNA)

dir_data <- "data/processed/GSE14520/"
dir.create(dir_data, recursive=TRUE, showWarning=FALSE)

fn <- "data/processed/GSE14520/expr.RMA_log2.gmean.rds"
ex <- readRDS(fn)

fn <- "data/processed/tcga_hcc/expr/tcga.gencodev26.rds"
gencode <- readRDS(fn)
symb2ens <- rownames(gencode)
names(symb2ens) <- gencode[,"Name"]

fn <- "data/processed/GSE14520/sdata.txt"
sdat <- read.table(fn, row.names=1, header=TRUE, sep="\t", stringsAsFactors=FALSE)

fn <- "data/processed/GSE14520/geo_pheno.txt"
geo_pheno <- read.table(fn, row.names=1, header=TRUE, sep="\t", stringsAsFactors=FALSE)

min_gene <- 20
max_gene <- 300
weighted <- FALSE

cl_id <- "cell_type_main"

# marker data
fn <- paste0("exp/sharma_aiz/markers/markers.", cl_id, ".txt")
markers.celltype <- read.table(fn, header=TRUE, stringsAsFactors=FALSE)

k <- markers.celltype[,"Name"] %in% rownames(ex)
markers.celltype <- markers.celltype[k,]
ct.count <- table(markers.celltype[,"cluster"])
ctk <- names(ct.count)[ct.count >= min_gene]
k <- markers.celltype[,"cluster"] %in% ctk
markers.celltype <- markers.celltype[k,]

ct_id <- "cluster"

# run decomposition
ex.eset <- ExpressionSet(ex)

weighted <- FALSE

# cell types
ex.ct.md <- MarkerBasedDecomposition(bulk.eset = ex.eset, 
                                     markers = markers.celltype, 
                                     ct_col = ct_id, 
                                     gene_col = "Name", 
                                     weighted = weighted, 
                                     unique_markers = FALSE,
                                     min_gene = min_gene, 
                                     max_gene = max_gene)
ex.ct.mdp <- ex.ct.md$bulk.props

# get PC1 of all markers for each cell type
markers <- markers.celltype[,"Name"]
cell_types <- markers.celltype[,ct_id]
cell_types.u <- unique(cell_types)
names(cell_types.u) <- cell_types.u
ct_pcs <- lapply(cell_types.u, function(ct){
                 g <- markers[cell_types == ct]
                 pcr <- prcomp(t(ex[g,]), retx=TRUE, scale.=TRUE, center=TRUE)
                 pcrret <- pcr$x[,1]
                 cors <- cor(t(ex[g,]), pcrret)
                 if (mean(cors) < 0) pcrret <- -1 * pcrret
                 return(pcrret) })
ct_pcs <- do.call(cbind, ct_pcs)
# this shows they are highly correlated
sort(diag(cor(ct_pcs, t(ex.ct.mdp)[rownames(ct_pcs),colnames(ct_pcs)])))

# correlate markers with cell-type proportions
ct_m_cors_l <- list()
for (ct in names(ex.ct.md$genes.used)){
    ens <- ex.ct.md$genes.used[[ct]]
    s <- colnames(ex)
    p_cors <- cor(t(ex[ens,s]), t(ex.ct.md$bulk.props[ct,s,drop=FALSE]))
    ct_m_cors_l[[ct]] <- p_cors
}

# output results
out_dir <- "data/processed/GSE14520/ctp/"
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

fn <- paste0(out_dir, "LCI.", cl_id, ".decomp.rds")
saveRDS(ex.ct.md, fn)

# export markers used for decomp
cts <- lapply(names(ex.ct.md$genes.used), function(ct){
              rep(ct, times = length(ex.ct.md$genes.used[[ct]]))
                 })
cts <- do.call(c, cts)
ens <- lapply(names(ex.ct.md$genes.used), function(ct){
              symb2ens[ex.ct.md$genes.used[[ct]]] })
ens <- do.call(c, ens)
symbs <- lapply(names(ex.ct.md$genes.used), function(ct){
              ex.ct.md$genes.used[[ct]] })
symbs <- do.call(c, symbs)
cors <- lapply(names(ex.ct.md$genes.used), function(ct){
               ct_m_cors_l[[ct]][ex.ct.md$genes.used[[ct]],ct] })
cors <- do.call(c, cors)

mrk_datf <- data.frame("CellType" = cts, 
                       "Ensembl" = ens, 
                       "Symbol" = symbs, 
                       "R" = cors)

dir_exp <- paste0("exp/GSE14520/ctp.", cl_id, "/")
dir.create(dir_exp, showWarnings = FALSE, recursive = TRUE)
out_fn <- paste0(dir_exp, cl_id, ".marker_genes.txt")
write.table(mrk_datf, out_fn, row.names = FALSE, col.names = TRUE, 
            quote = FALSE, sep = '\t')

# Plot
library(ggplot2)
library(reshape2)
plot_dir <- paste0("exp/GSE14520/ctp.", cl_id, "/ctp_cor/")
dir.create(plot_dir, showWarnings=FALSE, recursive=TRUE)

rd_bu <- read.table('data/ref/colors/rd_bu_div.csv', header=FALSE)
rd_bu <- rd_bu[,1]

cell_types <- names(ex.ct.md$genes.used)

for (ct in cell_types){
    cors.a <- cor(t(ex[ex.ct.md$genes.used[[ct]],]), ex.ct.md$bulk.props[ct,])
    colnames(cors.a) <- ct
    cors.dfm <- melt(cors.a)
    p <- ggplot(cors.dfm, aes(x = Var2, y = Var1, fill = value)) + 
        geom_tile() + 
        scale_fill_gradientn(colors=rev(rd_bu), limits=c(-1,1), 
                             name = substitute(paste(italic(R)))) + 
        labs(x = "Markers", y = ct) + 
        theme_bw()
    outfn <- paste0(plot_dir, ct, ".marker.ctp.cor.pdf")
    ggsave(outfn, width = 3, height = 6)
}

# plot co-expression
plot_dir <- paste0("exp/GSE14520/ctp.", cl_id, "/marker_coexpr/")
dir.create(plot_dir, showWarnings=FALSE, recursive=TRUE)

cors.l <- lapply(cell_types, function(ct){
                 gs <- ex.ct.md[["genes.used"]][[ct]]
                 cors.df <- cor(t(ex[gs,]))
                 cors.dfm <- reshape2::melt(cors.df)
                 cors.dfm[,"CellType"] <- ct
                 return(cors.dfm) })
cors <- do.call(rbind, cors.l)

facet_ncol = 2

outfn <- paste0(plot_dir, "connect.cell_type.marker_coexpr.pdf")
p <- ggplot(cors, aes(x = Var1, y = Var2, fill = value)) + 
    geom_tile(aes(fill = value)) + 
    scale_fill_gradientn(colors=rev(rd_bu), limits=c(-1,1),
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


