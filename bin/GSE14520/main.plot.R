
# estimate cell type proportions in the LCI data
# using the markers

setwd("../../")

library(NMF)
library(ggplot2)
library(edgeR)
library(Biobase)
library(plyr)
library(reshape2)
library(RColorBrewer)
source("scripts/sc_func.R")
# source("scripts/reference_free.R")
# library(Bisque)

#========================================================
# functions
#========================================================

scale_clip <- function(datf, z = 3){
    y <- apply(datf, 1, function(x) (x - mean(x))/sd(x))
    y <- t(y)
    colnames(y) <- colnames(datf)
    y[y < -z] <- -z
    y[y > z] <- z
    return(y)
}


GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL){
  data <- transform(data, 
                    xminv = x - violinwidth * (x - xmin), 
                    xmaxv = x + violinwidth * (xmax - x))
  grp <- data[1,'group']
  newdata <- plyr::arrange(transform(data, x = if(grp%%2==1) xminv else xmaxv), 
                           if(grp%%2==1) y else -y)
  newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
  newdata[c(1,nrow(newdata)-1,nrow(newdata)), 'x'] <- round(newdata[1, 'x']) 
  if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
    stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 
                                              1))
    quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
    aesthetics <- data[rep(1, nrow(quantiles)), 
                       setdiff(names(data), c("x", "y")), drop = FALSE]
    aesthetics$alpha <- rep(1, nrow(quantiles))
    both <- cbind(quantiles, aesthetics)
    quantile_grob <- GeomPath$draw_panel(both, ...)
    ggplot2:::ggname("geom_split_violin", 
                     grid::grobTree(GeomPolygon$draw_panel(newdata, ...), 
                                    quantile_grob))
  }
  else {
    ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
  }
})

geom_split_violin <- function (mapping = NULL, data = NULL, 
                               stat = "ydensity", position = "identity", ..., 
                               draw_quantiles = NULL, trim = TRUE, scale = "area", 
                               na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

#========================================================
#========================================================

options(stringsAsFactors=FALSE)

cl_id <- "cell_type_main"

# read in data
dir_data <- "data/processed/GSE14520/"
dir.create(dir_data, recursive=TRUE, showWarning=FALSE)

# expression data
fn <- paste0(dir_data, "expr.RMA_log2.gmean.rds")
ex <- readRDS(fn)

fn <- paste0(dir_data, "sdata.txt")
sdat <- read.table(fn, row.names=1, header=TRUE, sep="\t", stringsAsFactors=FALSE)

fn <- paste0(dir_data, "geo_pheno.txt")
geo_pheno <- read.table(fn, row.names=1, header=TRUE, sep="\t", stringsAsFactors=FALSE)

# cell-type proportions
fn <- paste0("data/processed/GSE14520/ctp/", "LCI.", cl_id, ".decomp.rds")
ex.ct.md <- readRDS(fn)

ex.ct.mdp <- ex.ct.md$bulk.props

# colors
fn <- "data/ref/colors/tumor_colrs.csv"
tum_cols <- read.table(fn, header=FALSE, sep=",")
tum_colsv <- tum_cols[,2]
names(tum_colsv) <- tum_cols[,1]
names(tum_colsv)[names(tum_colsv) == "Tumor"] <- "Tumor"
names(tum_colsv)[names(tum_colsv) == "NonTumor"] <- "Non-tumor"

fn <- "data/ref/colors/rd_bu_div.csv"
rd_bu <- read.table(fn)
rd_bu <- rd_bu[,1]


# marker data
fn <- paste0("exp/sharma_aiz/markers/markers.", cl_id, ".txt")
markers.celltype <- read.table(fn, header=TRUE, stringsAsFactors=FALSE)
ctk <- markers.celltype[,"Name"] %in% rownames(ex)
markers.celltype <- markers.celltype[ctk,]

# Read paired test results
fn <- paste0("exp/GSE14520/ctp.", cl_id, "/lci.", cl_id, ".tum_nontum.txt")
paired.test <- read.table(fn, header=TRUE, row.names=1, stringsAsFactors=FALSE, sep="\t")

# Plotting directory
dir_plot <- paste0("exp/GSE14520/ctp.", cl_id, "/plots/")
dir.create(dir_plot, showWarnings = FALSE, recursive = TRUE)

#========================================================
#========================================================

#========================================================
# Get tumor and non-tumor sample data separately
#========================================================

sdat.tum <- sdat[sdat[,"Tissue.Type"] == "Tumor",]
sdat.nt <- sdat[sdat[,"Tissue.Type"] == "Non-tumor",]

table(duplicated(sdat.tum[,"ID"]))
table(duplicated(sdat.nt[,"ID"]))

rownames(sdat.tum) <- sdat.tum[,"ID"]
rownames(sdat.nt) <- sdat.nt[,"ID"]

ids.pair <- intersect(rownames(sdat.tum), rownames(sdat.nt))
sdat.tum.p <- sdat.tum[ids.pair,]
sdat.nt.p <- sdat.nt[ids.pair,]

lcs.tum.p <- sdat.tum.p[,"LCS.ID"]
lcs.nt.p <- sdat.nt.p[,"LCS.ID"]

rownames(sdat.tum) <- sdat.tum[,"LCS.ID"]
rownames(sdat.nt) <- sdat.nt[,"LCS.ID"]

#========================================================
# Process and integrate
#========================================================

# scale mean 0 var 1
ex.ct.mdp.s <- t(apply(ex.ct.mdp, 1, scale))
colnames(ex.ct.mdp.s) <- colnames(ex.ct.mdp)
ex.ct.mdp.s.df <- as.data.frame(t(ex.ct.mdp.s))
ex.ct.mdp.s.df <- ex.ct.mdp.s.df[c(lcs.tum.p, lcs.nt.p),]
ex.ct.mdp.s.df[,"Type"] <- sdat[rownames(ex.ct.mdp.s.df), "Tissue.Type"]

# melt for plotting
ctp.dfm <- melt(ex.ct.mdp.s.df, id.vars = "Type")
colnames(ctp.dfm) <- c("Type", "CellType", "Proportion")

#========================================================
# bar plot of scaled cell type proprs
#========================================================

o <- order(paired.test[,"t.statistic"], decreasing=TRUE)
ctp.dfm[, "CellType"] <- factor(ctp.dfm[, "CellType"], 
                               levels = rownames(paired.test)[o])

ylims <- boxplot.stats(ctp.dfm[,"Proportion"])$stats[c(1,5)]
p <- ggplot(ctp.dfm, aes(x = CellType, y = Proportion, fill = Type)) + 
    geom_boxplot(outlier.shape=NA, size=0.2, width = 0.74) + 
    scale_fill_manual(values = tum_colsv) + 
    scale_y_continuous(limits = ylims) + 
    xlab(NULL) + ylab("Proportion (scaled)") + 
    ggtitle("LCI") + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle=45, hjust=0.9, vjust = 1), 
          text = element_text(size = 8),
          plot.title = element_text(hjust = 0.5, size = 12), 
          legend.box.margin = margin(0,0,0,0), 
          legend.box.spacing = unit(0, "npc"), 
          legend.position = "bottom", 
          legend.text = element_text(size = 8), 
          legend.title = element_text(size = 8, hjust = 0.5), 
          legend.key.size = unit(.05, units = "npc"))
outfn <- paste0(dir_plot, "ct_prop.tum.paired.pdf")
ggsave(outfn, width = 3, height = 2)

#========================================================
# split violin plot
#========================================================

# order factor
ctp.dfm[,"CellType"] <- factor(ctp.dfm[,"CellType"], 
                                levels = paired.test[,"CellType"])

ylims <- boxplot.stats(ctp.dfm[,"Proportion"])$stats[c(1,5)]
p <- ggplot(ctp.dfm, aes(x = CellType, y = Proportion, fill = Type)) + 
    geom_split_violin(trim = FALSE, scale="width") + 
    stat_summary(fun = mean, fun.min = mean, fun.max = mean,
                 geom = "crossbar", 
                 width = 0.25,
                 position = position_dodge(width = .25),
                 ) +
    scale_fill_manual(values = tum_colsv) + 
    xlab(NULL) + ylab("Proportion (scaled)") + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust = 0.5), 
          legend.position = "bottom")
outfn <- paste0(dir_plot, "ct_prop.pair.splt_vln.pdf")
ggsave(outfn, width = 4, height = 4)

#========================================================
# heatmap of scaled cell type proprs
#========================================================

# order by sample
o <- order(sdat[,"Tissue.Type"], decreasing=TRUE)
sdat.o <- sdat[o,]
ex.ct.mdp.o <- ex.ct.mdp[,rownames(sdat.o)]

# order samples by tumor status
o <- order(paired.test[,"t.estimate"], decreasing=TRUE)
cto <- rownames(paired.test)[o]
ex.ct.mdp.o <- ex.ct.mdp.o[cto,]
ex.ct.mdp.o <- scale_clip(ex.ct.mdp.o)

anno_cols <- list("Type" = tum_colsv)
anno <- list("Type" = sdat.o[,"Tissue.Type"])
fn <- paste0(dir_plot, "ctp.tum_nontum.heatmap.pdf")

aheatmap(ex.ct.mdp.o, Rowv=NA, Colv=NA, breaks=0, color=rev(rd_bu),
         annCol=anno, labCol=NA, filename=fn, annColors = anno_cols,
         main="Cell Type Proportions", width=10, height=4)

#========================================================
# clustered heatmap of scaled cell type proprs
#========================================================

ex.pt <- ex.ct.mdp[,rownames(sdat.tum)]
ex.nt <- ex.ct.mdp[,rownames(sdat.nt)]

ct2ord <- "Prol"
hmethod <- "complete"
hco.pt <- hclust(dist(t(ex.pt)), method = hmethod)$order
hco.pt <- order(ex.pt[ct2ord,], decreasing=FALSE)
hco.nt <- hclust(dist(t(ex.nt)), method = hmethod)$order
hco.nt <- order(ex.nt[ct2ord,], decreasing=FALSE)

ex.ct.mdp.o <- cbind(ex.pt[,hco.pt], ex.nt[,hco.nt])

o <- order(paired.test[,"t.estimate"], decreasing=TRUE)
cto <- rownames(paired.test)[o]
ex.ct.mdp.o <- ex.ct.mdp.o[cto,]
ex.ct.mdp.o <- scale_clip(ex.ct.mdp.o)

anno_cols <- list("Type" = tum_colsv)
anno <- list("Type" = sdat[colnames(ex.ct.mdp.o),"Tissue.Type"])
fn <- paste0(dir_plot, "ctp.tum_nontum.o.heatmap.pdf")

aheatmap(ex.ct.mdp.o, Rowv=NA, Colv=NA, breaks=0, color=rev(rd_bu),
         annCol=anno, labCol=NA, filename=fn, annColors = anno_cols,
         main="Cell Type Proportions", width=10, height=4)

#========================================================
# plot t statistics
#========================================================

theme_leg <- theme(legend.position = c(0.95, 0.95),
                   legend.key.size = unit(0, "npc"),
                   legend.justification = c("right", "top"),
                   legend.box.background = element_rect(color="black"), 
                   legend.margin = margin(0,4,4,4))

paired.test[,"CellType"] <- rownames(paired.test)
paired.test[,"CellType"] <- factor(paired.test[,"CellType"], 
                                levels = paired.test[,"CellType"])
o <- order(paired.test[,"t.statistic"], decreasing=TRUE)
paired.test[,"CellType"] <- factor(paired.test[,"CellType"], levels = paired.test[,"CellType"][o])
maxlim <- max(abs(paired.test[,"t.statistic"]))

paired.test[,"sig"] <- "ns"
paired.test[ paired.test[,"w.p_adj"] < 0.05 , "sig"] <- "sig"
sig_shape <- c("ns" = 1, "sig" = 16)

p <- ggplot(paired.test, aes(x = CellType, y = t.statistic)) + 
    geom_hline(yintercept = 0, col="red") + 
    geom_point(aes(shape = sig), size=3) + theme_bw() + 
    scale_shape_manual(values = sig_shape, breaks = "sig", 
                       name = NULL, 
                       labels = c("sig" = "adj p < 0.05")) + 
    ylim(c(-maxlim, maxlim)) + 
    ylab("T statistic\n(tumor over non-tumor)") + 
    xlab(NULL) + 
    ggtitle("LCI") +
    theme(text = element_text(size = 10),
          axis.text.x = element_text(angle = 45, hjust = 0.9, color="black"), 
          axis.text.y = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          plot.title = element_text(hjust = 0.5, size = 12)) +  
    theme_leg
fn <- paste0(dir_plot, "tum_nontum.t.pdf")
ggsave(fn, width = 3, height = 2.5)



#========================================================
# heatmap of expression values
#========================================================

dir_gene <- paste0(dir_plot, "cte/")
dir.create(dir_gene, showWarnings = FALSE, recursive = TRUE)

o <- order(sdat[,"Tissue.Type"])
sdat.o <- sdat[o,]
ex.ct.mdp.o <- ex.ct.mdp[,rownames(sdat.o)]

ct2plot <- unique(markers.celltype[,"cluster"])
for (ct in ct2plot){
    kct <- markers.celltype[,"cluster"] == ct
    markers.cts <- markers.celltype[kct,]
    ct_genes <- markers.cts[1:min(50,nrow(markers.cts)), "Name"]
    ct_genes <- intersect(ct_genes, rownames(ex))
    ct.expr <- ex[ct_genes,rownames(sdat.o)]

    ct.expr <- scale_clip(ct.expr)

    anno_cols <- list("Type" = tum_colsv)
    anno <- list("Type" = sdat.o[,"Tissue.Type"])

    fn <- paste0(dir_gene, "tcga.TMM.", ct, ".tum_nontum.pdf")
    aheatmap(ct.expr, Rowv=NA, Colv=NA, breaks=0, color=rev(rd_bu), 
             annCol=anno, labCol=NA, filename=fn, annColors = anno_cols,
             main = paste0(ct, " marker genes"), width=12, height=6)
}


