
# Figure 1

setwd("../../")

library(Seurat)
library(ggplot2)
library(grid)
library(gtable)
library(ggrepel)
library(scales)
library(gridExtra)
library(reshape2)
source("scripts/color_pal.R")
source("scripts/get_enr_mat.R")
source("scripts/ggplot_raster.R")
source("scripts/ggplot_formats.R")
source("scripts/gtable_stack.R")

#=========================================
# Functions
#=========================================

create_dir <- function(p){
    dir.create(p, showWarnings=FALSE, recursive=TRUE)
}

#' @param datf data frame to get median points from
#' @param x character vector of column names to calculate median for
#' @param groupby the column name to group the rows of the data frame by
get_med_points <- function(datf, x, groupby){
    groups <- sort(unique(datf[,groupby]))
    gs.l <- lapply(groups, function(gr){
                 k <- datf[,groupby] == gr
                 datf.s <- datf[k, x, drop=FALSE]
                 r <- apply(datf.s, 2, median)
                 return(r) })
    names(gs.l) <- groups
    gs <- as.data.frame(do.call(rbind, gs.l))
    colnames(gs) <- x
    rownames(gs) <- groups
    return(gs)
}

mar_pct <- 0.025 # pct margin around panel
tagl <- 3

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

# Set directories
dir_plt <- "exp/sharma_aiz/plots/";
dir.create(dir_plt, showWarnings=FALSE, recursive=TRUE)

fn <- "data/processed/sharma_aiz/liver.int_rand.md_umap.txt.gz"
udf <- read.table(fn, header=TRUE, row.names=1, sep='\t')
colnames(udf)[colnames(udf) == "UMAP_1"] <- "UMAP1"
colnames(udf)[colnames(udf) == "UMAP_2"] <- "UMAP2"
set.seed(1, kind = 'Mersenne-Twister')
udf <- udf[sample(1:nrow(udf)),]

tremap <- c("Tumor" = "Tumor", "NonTumor" = "Non-tumor")
#rename udf TumorStat
udf[,"TumorStat"] <- tremap[udf[,"TumorStat"]]

# G2M and S diff between cell-types
fn <- "exp/sharma_aiz/cc_diff/cell_type_fine.cc_diff.txt"
cc_diff <- read.table(fn, header=TRUE, sep='\t')

# tumor prop diff stats
fn <- "exp/sharma_aiz/tum_diff/cell_type_fine_p.tum_diff.txt"
diffs <- read.table(fn, header=TRUE)

# colors
tum_df <- read.csv("data/ref/colors/tumor_colrs.csv", row.names=1, header=FALSE)
tum_col <- tum_df[,1]; names(tum_col) <- rownames(tum_df)
names(tum_col) <- tremap[names(tum_col)]

fn <- "data/ref/colors/red_colrs.csv"
reds <- read.table(fn, header=FALSE)
reds <- reds[,1]

dir_plt <- "exp/manuscript/"
dir.create(dir_plt, showWarnings=FALSE, recursive=TRUE)

#=========================================
# UMAP of cell-type assignments
#=========================================

th_um <- theme(text = element_text(size = 8),
               plot.title = element_text(hjust = 0.5), 
               legend.key.size = unit(6, "points"), 
               legend.spacing.x = unit(2,"points"), 
               legend.spacing.y = unit(5,"points"), 
               legend.title = element_text(vjust = 0.5, hjust = 0), 
               legend.margin = margin(0,0,0,0), 
               legend.text = element_text(size = 6, hjust = 0, margin = margin(0,5,0,0)),
               legend.box.spacing = unit(0.5, "strheight", "0"), 
               axis.text=element_blank(),
               axis.ticks=element_blank())
th_um <- theme_classic() + th_um

# Cell type labels
labs <- c("cell_type_main", "cell_type_fine")
names(labs) <- labs
leg <- c("cell_type_main" = "Main", "cell_type_fine" = "Fine")
titl <- c("cell_type_main" = "Main cell-types", 
          "cell_type_fine" = "Subcell-types")

# generate color palettes
ctu <- lapply(labs, function(lab){
              sort(unique(udf[,lab])) })
pal_labs <- lapply(labs, function(l){
                   nct <- length(ctu[[l]])
                   hcl_pal(nct, chr = c(80, 80), lum = c(60,80), 
                           offset = 0, rand = TRUE, seedn = 3) })

meds_fine <- get_med_points(udf, c("UMAP1", "UMAP2"), "cell_type_fine")
meds_main <- get_med_points(udf, c("UMAP1", "UMAP2"), "cell_type_main")
meds_l <- list("cell_type_main" = meds_main, 
               "cell_type_fine" = meds_fine)
size_l <- list("cell_type_main" = 3, 
               "cell_type_fine" = 2)

umap_l <- list()

for (l in labs){
    p <- ggplot(udf, aes_string(x="UMAP1",y="UMAP2",color=l)) +
    geom_point(shape=16, size=.01) +
    scale_color_manual(values=pal_labs[[l]], name=NULL) +
    ggtitle(titl[l]) + 
    th_um + 
    coord_fixed() + 
    guides(color = guide_legend(override.aes = list(size=1), ncol=1)) + 
    geom_text_repel(data=meds_l[[l]], aes(x=UMAP1, y=UMAP2), 
                    label=rownames(meds_l[[l]]), size = size_l[[l]], 
                    force = 12, seed = 1, 
                    segment.size = 0.2, color="black")
    umap_l[[l]] <- ggplotGrob(p)
}

umap_l[[1]] <- raster_ggpoints(umap_l[[1]], w = 3, h = 3)
umap_l[[2]] <- raster_ggpoints(umap_l[[2]], w = 3, h = 3)

umap_l[[1]] <- gtable_add_tag(umap_l[[1]], tag = "a", fs=16, just = c(0,1), l=tagl)
umap_l[[2]] <- gtable_add_tag(umap_l[[2]], tag = "b", fs=16, just = c(0,1), l=tagl)

umap_l[[1]] <- pad_plot(umap_l[[1]], t=mar_pct, r=mar_pct, b=mar_pct, l=mar_pct)
umap_l[[2]] <- pad_plot(umap_l[[2]], t=mar_pct, r=mar_pct, b=mar_pct, l=mar_pct)

#=========================================
# pathway enrichment plots
#=========================================

fn <- "exp/sharma_aiz/clust2ct/res1_to_fine.txt"
ctmap <- read.table(fn, header=TRUE, colClasses = "character")
rownames(ctmap) <- ctmap[,1]

fn <- "exp/sharma_aiz/path_enr/markers.res.1.path_gse.rds"
gse_all <- readRDS(fn)

theme_e <- theme_bw() + 
    theme(text = element_text(size = 10), 
          legend.key.size = unit(10, "points"), 
          axis.text = element_text(colour = "black", size = 6), 
          axis.text.y = element_text(margin=margin(0,1,0,0)), 
          axis.text.x = element_text(margin=margin(0,0,0,0), 
                                     angle = 90, hjust = 1, vjust = 0.5), 
          panel.border =  element_blank(), 
          panel.background = element_blank(), 
          panel.grid.major = element_line(colour = reds[1]), 
          plot.title = element_text(hjust = 0.5), 
          axis.ticks = element_line(colour = reds[1]))

e_col <- "NES"
p_col <- "qvalues"

a <- "Reactome"
react_l <- gse_all[[a]]
names(react_l) <- ctmap[names(react_l),"cell_type"]
react_l <- lapply(react_l, function(x){
                  k <- x[,e_col] >= 0
                  x <- x[k,,drop=FALSE] 
                  return(x) })
react_l <- get_enr_mat(react_l, idcol = "ID", dcol = "Description", 
                       ctcol = "CellType", pcol = p_col, ecol = e_col, 
                       max_terms = 2)

e_df <- react_l[["E"]]
p_df <- react_l[["P"]]

# wrap strings
rn <- rownames(e_df)
rn <- sapply(rn, function(x){
             chlim <- 80
             if (nchar(x) > chlim){
                 x <- strtrim(x, chlim-3)
                 x <- paste(x, "...", sep="")
             }
             return(x) })
# rn <- sapply(rn, function(x) paste(strwrap(x), collapse='\n'))
rownames(e_df) <- rownames(p_df) <- rn

# order cell types
cl <- sort(colnames(e_df))
ro <- hclust(dist(e_df[,cl]))$order
rl <- rownames(e_df)[ro]

# melt
p_dfm <- reshape2::melt(as.matrix(p_df))
e_dfm <- reshape2::melt(as.matrix(e_df))
for (i in 1:2){
    p_dfm[,i] <- factor(as.character(p_dfm[,i]), levels=list(rl,cl)[[i]])
    e_dfm[,i] <- factor(as.character(e_dfm[,i]), levels=list(rl,cl)[[i]])
}

enr_datf <- data.frame("path" = p_dfm[,1], "cluster" = p_dfm[,2], 
                       "enr" = e_dfm[,3], "p" = p_dfm[,3])

enr_p <- ggplot(enr_datf, aes(x = cluster, y = path, color = enr, size = p)) + 
geom_point() + 
theme_e + labs(x=NULL, y=NULL) + ggtitle(paste0(a, " pathway enrichment")) + 
scale_size(name = "q-value", trans="log10_rev", range = c(0,3)) + 
scale_color_gradientn(colors = reds) + 
guides(size = guide_legend(override.aes = list(shape = 1, color = "black")))

enr_g <- ggplotGrob(enr_p)
enr_g <- gtable_add_tag(enr_g, tag = "c", fs = 16, just = c(0,1), l = tagl)
enr_g <- pad_plot(enr_g, t = mar_pct, r=mar_pct, b=mar_pct, l=mar_pct)

#=========================================
# bar plot of tumor proprs
#=========================================

th_bar <- theme_classic() + 
theme(text = element_text(size = 8),
      axis.text = element_text(colour = "black", size = 8), 
      axis.text.x = element_text(size = 6, hjust = 0.5, vjust = 0.5, angle=0), 
      axis.text.y = element_text(size = 6), 
      plot.title = element_text(hjust = 0.5), 
      legend.key.size = unit(8, "points"), 
      legend.position = "right", 
      legend.spacing.x = unit(0,"points"), 
      legend.box.spacing = unit(2, "points"), 
      legend.margin = margin(0,5,0,0, unit='pt'), 
      axis.line.x = element_blank())

# tumor props
tum_prop_f <- table(udf[,"cell_type_fine"], udf[,"TumorStat"])
tum_prop_f <- sweep(tum_prop_f, 1, rowSums(tum_prop_f), '/')
tum_prop_fm <- melt(tum_prop_f)
colnames(tum_prop_fm) <- c("CellType", "Tumor", "Prop")

# get significance
cts <- rownames(diffs)
sig <- rep("adj p >= 0.05", length.out=length(cts))
names(sig) <- cts
sig[diffs[,"w_p_adj"] < .05] <- "adj p < 0.05"

# order by decreasing prop
o <- order(tum_prop_f[,"Tumor"])
lv <- rownames(tum_prop_f)[o]
tum_prop_fm[,"CellType"] <- factor(tum_prop_fm[,"CellType"], levels=lv)
tum_prop_fm[,"sig"] <- sig[as.character(tum_prop_fm[,"CellType"])]
al <- c("adj p < 0.05" = 1, "adj p >= 0.05" = 0.4)

# plot cell-type tumor prop
tump_p <- ggplot(tum_prop_fm, aes(x = CellType, y = Prop, fill = Tumor)) + 
    geom_bar(aes(alpha = sig), position="fill", stat="identity") + 
    coord_flip() + 
    labs(y="Tumor proportion",x=NULL) + 
    scale_y_continuous(expand = expansion(0)) + 
    scale_fill_manual(values = tum_col, name=NULL) + 
    scale_alpha_manual(values = al, name = NULL) + 
    th_bar

tump_g <- ggplotGrob(tump_p)
tump_g <- gtable_add_tag(tump_g, tag = "d", fs = 16, just = c(0,0), l=tagl)
tump_g <- pad_plot(tump_g, t=mar_pct, r=0, b=mar_pct, l=mar_pct)

#=========================================
# individual cell-type props per sample sep. tumor
#=========================================

# plotting theme for box
th_bx <- theme(text = element_text(size = 10),
               axis.text = element_text(colour = "black", size = 8), 
               axis.line = element_line(size = rel(0.5)), 
               plot.title = element_text(hjust = 0.5), 
               legend.position = "none", 
               axis.text.x = element_text(hjust = 0.9, vjust = 1, angle=30))
th_bx <- theme_classic() + th_bx

fn <- "exp/sharma_aiz/tum_diff/cell_type_fine.sample_prop.txt"
props_main <- read.table(fn, header=TRUE)
props_main[,"Tumor"] <- tremap[props_main[,"Tumor"]] # change NonTumor to Non-tumor

# get significant cell-types
cts <- rownames(diffs)
ksig <- diffs[,"w_p_adj"] < .05
cts_sig <- cts[ksig]
cts_sig <- c("Prol", setdiff(cts_sig, "Prol"))

ct_s_l <- list()
for (ct in cts_sig){
    props_ct <- props_main[props_main[,"CellType"] == ct,,drop=FALSE]

    # add p-val text
    wpadj <- diffs[ct, "w_p_adj"]
    wpadj <- get_astk(wpadj)
    p_y <- max(props_ct$Prop) + (diff(range(props_ct$Prop)) * .1)
    p_x <- 1.5

    p <- ggplot(props_ct, 
                aes(x = Tumor, y = Prop, fill = Tumor)) + 
    geom_boxplot(size = .2, outlier.shape=NA) + 
    geom_jitter(width = 0.1, size = 0.5) + 
    labs(x=NULL, y="Proportion") + 
    ggtitle(ct) + 
    scale_y_continuous(expand = expansion(mult = c(0.1, .25))) + 
    scale_fill_manual(values = tum_col) + 
    geom_text(label = wpadj, x = p_x, y = p_y, color = "black", 
              size = 6, hjust = 0.5, vjust = 0.3, parse = FALSE) + 
    geom_segment(x=1,xend=2, y=p_y, yend=p_y, color="black") + 
    th_bx

    if (which(ct == cts_sig) > 1)
        p <- p + ylab(NULL)

    ct_s_l[[ct]] <- ggplotGrob(p)
}

prol_g <- ct_s_l[[1]]
prol_g <- gtable_add_tag(prol_g, tag = "e", fs = 16, just = c(0,0), l=tagl)
prol_g <- pad_plot(prol_g, t=mar_pct, r=0, b=mar_pct, l=0)

#=========================================
# Plot G2M and S Score
#=========================================

th_um <- theme(text = element_text(size = 8),
               plot.title = element_text(hjust = 0.5), 
               legend.key.height = unit(6, "points"), 
               legend.key.width = unit(12, "points"), 
               legend.title = element_text(vjust = 1, hjust = 1), 
               legend.margin = margin(0,0,0,0), 
               legend.text = element_text(size = 6, angle = 30, hjust=0.9),
               legend.position = "bottom", 
               legend.box.spacing = unit(2, "points"), 
               axis.text=element_blank(),
               axis.ticks=element_blank())
th_um <- theme_classic() + th_um

# get T statistic for Prol difference in G2M and S scores
k <- cc_diff[,"CellType"] == "Prol" & cc_diff[,"Score"] == "G2M"
g2m_wpadj <- cc_diff[k, "Wpadj"]
g2m_wpadj <- get_astk(g2m_wpadj)
k <- cc_diff[,"CellType"] == "Prol" & cc_diff[,"Score"] == "S"
s_wpadj <- cc_diff[k, "Wpadj"]
s_wpadj <- get_astk(s_wpadj)

uxmax <- max(udf[,"UMAP1"])
uymin <- min(udf[,"UMAP2"])

feat <- "G2M.Score"; lt <- "G2M Score"
o <- order(udf[,"G2M.Score"], decreasing=FALSE)
p1 <- ggplot(udf[o,], aes(x=UMAP1, y=UMAP2, color=G2M.Score)) + 
geom_point(size = 0.01, shape = 16) + 
coord_fixed() + 
ggtitle("G2M phase") + 
geom_text(label = paste0(g2m_wpadj), x = uxmax, y = uymin, 
          color = "black", size = 6, hjust = 1, vjust = 0) + 
scale_color_gradientn(colours = reds, name =  lt) + 
th_um

feat <- "S.Score"; lt <- "S Score"
o <- order(udf[,"S.Score"], decreasing=FALSE)
p2 <- ggplot(udf[o,], aes(x=UMAP1, y=UMAP2, color=S.Score)) + 
geom_point(size = 0.01, shape = 16) + 
coord_fixed() + 
ggtitle("S phase") + 
geom_text(label = paste0(s_wpadj), x = uxmax, y = uymin, 
          color = "black", size = 6, hjust = 1, vjust = 0) + 
scale_color_gradientn(colours = reds, name =  lt) + 
th_um

p1 <- ggplotGrob(p1)
p2 <- ggplotGrob(p2)
p1 <- raster_ggpoints(p1, w = 3, h = 3)
p2 <- raster_ggpoints(p2, w = 3, h = 3)
p1 <- gtable_add_tag(p1, tag = "f", fs = 16, just = c(0,-.7), l=tagl)
p2 <- gtable_add_tag(p2, tag = "g", fs = 16, just = c(0,-.7), l=tagl)
p1 <- pad_plot(p1, t=mar_pct, r=mar_pct, b=mar_pct, l=0)
p2 <- pad_plot(p2, t=mar_pct, r=mar_pct, b=mar_pct, l=0)


g2m_s_gt <- gtable(widths = unit(c(0.5, 0.5), "npc"), heights = unit(1, "npc"))
g2m_s_gt <- gtable_add_grob(g2m_s_gt, p1, t=1, b=1, l=1, r=1)
g2m_s_gt <- gtable_add_grob(g2m_s_gt, p2, t=1, b=1, l=2, r=2)

#=========================================
# put together a gtable for plotting everything
#=========================================

punit <- "inches"
wdta <- c(3.5, 3.5)
hgta <- c(3, 3.5)
gtop <- gtable(widths = unit(wdta, punit), 
               heights = unit(hgta, punit))
gtop <- gtable_add_grob(gtop, umap_l[[1]], t=1,b=1,l=1,r=1, name="a")
gtop <- gtable_add_grob(gtop, umap_l[[2]], t=1,b=1,l=2,r=2, name="b")
gtop <- gtable_add_grob(gtop, enr_g, t=2, b=2, l=1, r=2, name="c")


wdtb <- c(2, 1.5, 3.5)
hgtb <- c(2.5)
gbot <- gtable(widths = unit(wdtb, punit), 
               heights = unit(hgtb, punit))
gbot <- gtable_add_grob(gbot, tump_g, t=1,b=1,l=1,r=1, name="d")
gbot <- gtable_add_grob(gbot, prol_g, t=1,b=1,l=2,r=2, name="e")
gbot <- gtable_add_grob(gbot, g2m_s_gt, t=1,b=1,l=3,r=3, name="fg")

wdtc <- sum(sum(wdta))
hgtc <- c(sum(hgta), sum(hgtb))
gt <- gtable(widths = unit(wdtc, punit), 
             heights = unit(hgtc, punit))
gt <- gtable_add_grob(gt, gtop, t=1,b=1,l=1,r=1, name = "top")
gt <- gtable_add_grob(gt, gbot, t=2,b=2,l=1,r=1, name = "bot")

pdf(paste0(dir_plt, "fig1.pdf"), width = sum(wdtc), height = sum(hgtc))
grid.draw(gt)
dev.off()

