
# estimate cell type proportions in the TCGA data
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

fn <- "data/processed/tcga_hcc/expr/tcga.gencodev26.rds"
gencode <- readRDS(fn)

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

# TCGA expression data
fn <- "data/processed/tcga_hcc/expr/tcga.lihc.TMM.log.rin.rds"
tmm <- readRDS(fn)
rownames(tmm) <- gencode[ rownames(tmm), "Name"]

cl_id <- "cell_type_main"

# proportions
fn <- paste0("data/processed/tcga_hcc/ctp/tcga.TMM.", cl_id, ".decomp.rds")
tmm.ct.md <- readRDS(fn)
tmm.ct.mdp <- tmm.ct.md$bulk.props

for (ct in names(tmm.ct.md$genes.used)){
    tmm.ct.md$genes.used[[ct]] <- gencode[tmm.ct.md$genes.used[[ct]], "Name"]
}

# marker data
fn <- paste0("exp/sharma_aiz/markers/markers.", cl_id, ".txt")
markers.celltype <- read.table(fn, header=TRUE, stringsAsFactors=FALSE)
ctk <- markers.celltype[,"Name"] %in% rownames(tmm)
markers.celltype <- markers.celltype[ctk,]

# TCGA pheno data
fn <- "data/processed/tcga_hcc/sample/cases.hcc.361.merged.rds"
cases <- readRDS(fn)
fn <- "data/processed/tcga_hcc/sample/samples.hcc.410.merged.rds"
samples <- readRDS(fn)
fn <- "data/processed/tcga_hcc/sample/tcga.bio_sample.rds"
bio_sample <- readRDS(fn)
samples[,"sample_type2"] <- factor(samples[,"sample_type2"], levels = c("Tumor", "Non-tumor"))

# sample map
fn <- "data/processed/tcga_hcc/sample/sam_case_map.rds"
smaps <- readRDS(fn)
cid2sid_tum <- smaps[["cid2sid_tum"]]
cid2sid_nt <- smaps[["cid2sid_nt"]]
cid2sid_tum_pair <- smaps[["cid2sid_tum_pair"]]
cid2sid_nt_pair <- smaps[["cid2sid_nt_pair"]]
sid_both <- c(cid2sid_tum_pair, cid2sid_nt_pair)

# Read paired test results
fn <- paste0("exp/tcga_hcc/ctp.", cl_id, "/tcga.", cl_id, ".tum_nontum.txt")
paired.test <- read.table(fn, header=TRUE, row.names=1, stringsAsFactors=FALSE, sep="\t")

dir_exp <- paste0("exp/tcga_hcc/ctp.", cl_id, "/")
# Plotting directory
dir_plot <- paste0("exp/tcga_hcc/ctp.", cl_id, "/plots/")
dir.create(dir_plot, showWarnings = FALSE, recursive = TRUE)

#========================================================
#========================================================

tmm.ct.pair <- tmm.ct.mdp[,sid_both]

# scale mean 0 var 1
tmm.ct.pair.s <- t(apply(tmm.ct.pair, 1, scale))
colnames(tmm.ct.pair.s) <- colnames(tmm.ct.pair)
tmm.df <- as.data.frame(t(tmm.ct.pair.s))
tmm.df[,"Type"] <- "Tumor"
tmm.df[cid2sid_nt_pair, "Type"] <- "Non-tumor"

# melt for plotting
tmm.dfm <- melt(tmm.df, id.vars = "Type")
colnames(tmm.dfm) <- c("Type", "CellType", "Proportion")

out_fn <- paste0(dir_exp, "ctp_tumor_mlt.txt")
write.table(tmm.dfm, out_fn, row.names=FALSE, col.names=TRUE, 
            quote=FALSE, sep="\t")

#========================================================
# bar plot of scaled cell type proprs
#========================================================

o <- order(paired.test[,"t.statistic"], decreasing=TRUE)
tmm.dfm[, "CellType"] <- factor(tmm.dfm[, "CellType"], 
                               levels = rownames(paired.test)[o])

ylims <- boxplot.stats(tmm.dfm[,"Proportion"])$stats[c(1,5)]
p <- ggplot(tmm.dfm, aes(x = CellType, y = Proportion, fill = Type)) + 
    geom_boxplot(outlier.shape=NA, size=0.2, width = 0.74) + 
    scale_fill_manual(values = tum_colsv) + 
    scale_y_continuous(limits = ylims) + 
    xlab(NULL) + ylab("Proportion (scaled)") + 
    ggtitle("TCGA") + 
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
tmm.dfm[,"CellType"] <- factor(tmm.dfm[,"CellType"], 
                                levels = paired.test[,"CellType"])

ylims <- boxplot.stats(tmm.dfm[,"Proportion"])$stats[c(1,5)]
p <- ggplot(tmm.dfm, aes(x = CellType, y = Proportion, fill = Type)) + 
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
ggsave(outfn, width = 10, height = 4)

#========================================================
# heatmap of scaled cell type proprs
#========================================================

o <- order(samples[,"sample_type2"])
samples.o <- samples[o,]
tmm.ct.mdp.o <- tmm.ct.mdp[,rownames(samples.o)]

o <- order(paired.test[,"t.estimate"], decreasing=TRUE)
cto <- rownames(paired.test)[o]
tmm.ct.mdp.o <- tmm.ct.mdp.o[cto,]
tmm.ct.mdp.o <- scale_clip(tmm.ct.mdp.o)

anno_cols <- list("Type" = tum_colsv)
anno <- list("Type" = samples.o[,"sample_type2"])
fn <- paste0(dir_plot, "ctp.tum_nontum.heatmap.pdf")

aheatmap(tmm.ct.mdp.o, Rowv=NA, Colv=NA, breaks=0, color=rev(rd_bu),
         annCol=anno, labCol=NA, filename=fn, annColors = anno_cols,
         main="Cell Type Proportions", width=10, height=4)


#========================================================
# clustered heatmap of scaled cell type proprs
#========================================================

tmm.pt <- tmm.ct.mdp[,cid2sid_tum]
tmm.nt <- tmm.ct.mdp[,cid2sid_nt]

ct2ord <- "Prol"
hmethod <- "complete"
hco.pt <- hclust(dist(t(tmm.pt)), method = hmethod)$order
hco.pt <- order(tmm.pt[ct2ord,], decreasing=FALSE)
hco.nt <- hclust(dist(t(tmm.nt)), method = hmethod)$order
hco.nt <- order(tmm.nt[ct2ord,], decreasing=FALSE)

tmm.ct.mdp.o <- cbind(tmm.pt[,hco.pt], tmm.nt[,hco.nt])

o <- order(paired.test[,"t.estimate"], decreasing=TRUE)
cto <- rownames(paired.test)[o]
tmm.ct.mdp.o <- tmm.ct.mdp.o[cto,]
tmm.ct.mdp.o <- scale_clip(tmm.ct.mdp.o)

anno_cols <- list("Type" = tum_colsv)
anno <- list("Type" = samples.o[,"sample_type2"])
fn <- paste0(dir_plot, "ctp.tum_nontum.o.heatmap.pdf")

aheatmap(tmm.ct.mdp.o, Rowv=NA, Colv=NA, breaks=0, color=rev(rd_bu),
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
    ggtitle("TCGA") +
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

o <- order(samples[,"sample_type2"])
samples.o <- samples[o,]
tmm.ct.mdp.o <- tmm.ct.mdp[,rownames(samples.o)]

ct2plot <- unique(markers.celltype[,"cluster"])
for (ct in ct2plot){
    kct <- markers.celltype[,"cluster"] == ct
    markers.cts <- markers.celltype[kct,]
    ct_genes <- markers.cts[1:min(50,nrow(markers.cts)), "Name"]
    ct_genes <- intersect(ct_genes, rownames(tmm))
    ct.expr <- tmm[ct_genes,rownames(samples.o)]

    ct.expr <- scale_clip(ct.expr)

    anno_cols <- list("Type" = tum_colsv)
    anno <- list("Type" = samples.o[,"sample_type2"])

    fn <- paste0(dir_gene, "tcga.TMM.", ct, ".tum_nontum.pdf")
    aheatmap(ct.expr, Rowv=NA, Colv=NA, breaks=0, color=rev(rd_bu), 
             annCol=anno, labCol=NA, filename=fn, annColors = anno_cols,
             main = paste0(ct, " marker genes"), width=12, height=6)
}


