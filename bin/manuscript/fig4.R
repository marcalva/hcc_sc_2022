
# Figure 4

setwd("../../")

library(Seurat)
library(ggplot2)
library(grid)
library(gtable)
library(ggrepel)
library(scales)
library(gridExtra)
library(reshape2)
library(ggfortify)
source("scripts/color_pal.R")
source("scripts/get_enr_mat.R")
source("scripts/ggplot_raster.R")
source("scripts/gtable_stack.R")
source("scripts/ggplot_formats.R")

#=========================================
# Functions
#=========================================

create_dir <- function(p){
    dir.create(p, showWarnings=FALSE, recursive=TRUE)
}

mar_pct <- 0.025 # pct margin around panel
tagl <- 3

#=========================================
# data
#=========================================

ct_id <- "cell_type_main"

fn <- "data/processed/sharma_aiz/liver.int_rand.md_umap.txt.gz"
udf <- read.table(fn, header=TRUE, row.names=1, sep='\t')
colnames(udf)[colnames(udf) == "UMAP_1"] <- "UMAP1"
colnames(udf)[colnames(udf) == "UMAP_2"] <- "UMAP2"
set.seed(1, kind = 'Mersenne-Twister')
udf <- udf[sample(1:nrow(udf)),]


# mutation type
fn <- "exp/tcga_hcc/mol/cell_type_main/sample.mut_type.txt"
datf.type <- read.table(fn, header=TRUE, row.names=1, sep='\t')

# somatic mutation results
fn <- paste0("exp/tcga_hcc/mol/", ct_id, ".SMG_mutsig.2CV.ctp.txt")
wr.smg <- read.table(fn, header=TRUE, sep="\t")

fn <- paste0("exp/tcga_hcc/mol/", ct_id, ".SMG_mutsig.2CV.mut_type.ctp.txt")
wr.smg.type <- read.table(fn, header=TRUE, sep="\t")

# tumor enrichment scores
fn <- "exp/tcga_hcc/sharma_aiz.mut_score/sharma_aiz.mut_scores.txt.gz"
mut_scores <- read.table(fn, row.names=1, header=TRUE, sep='\t')

dir_plt <- "exp/manuscript/"
dir.create(dir_plt, showWarnings=FALSE, recursive=TRUE)

reds <- read.csv("data/ref/colors/red_colrs.csv", header=FALSE)
reds <- reds[,1]

#========================================================
# plot function/themes
#========================================================

th_txt <- theme(text = element_text(size = 8), 
                axis.text = element_text(size = 8), 
                axis.text.x = element_text(colour = "black"), 
                axis.text.y = element_text(colour = "black"), 
                plot.title = element_text(hjust = 0.5, size = 8))

bp <- function(datf, x, y){
    ylims <- range(datf[,y])
    r <- ylims[2] - ylims[1]
    rp <- r * 0.1
    ly <- ylims[2] + (rp * 0.5)

    p <- ggplot(datf, aes_string(x=x,y=y)) + 
        geom_boxplot(outlier.shape=NA) + 
        geom_dotplot(binaxis='y', stackdir='center', binwidth = .05) + 
        geom_segment(aes(x = 0, y = ly, xend = 1, yend = ly)) + 
        ylim(c(ylims[1], ylims[2] + tp)) 
}

#========================================================
# proportion box plot per mutations
#========================================================

th_mp <- theme_classic() + 
th_txt + 
theme(panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size=1), 
      axis.line = element_blank(),
      axis.text.y = element_blank(), 
      axis.ticks.y = element_blank(), 
      strip.text = element_text(size = 8), 
      strip.background = element_rect(fill = "white"))


fn <- "exp/tcga_hcc/mol/cell_type_main/sample.mut_wt_stat.txt"
datf.mut <- read.table(fn, header=TRUE, row.names=1)

datf.muts <- datf.mut[,c("TP53", "RB1", "Prol")]
datf.mutsm <- melt(datf.muts, id.vars="Prol")
datf.mutsm[,"value"] <- factor(datf.mutsm[,"value"], levels = c("WT", "Mut"))

ct_pval <- wr.smg[ wr.smg$CellType == "Prol" & wr.smg$Gene == "TP53", "w.p_adj"]

pval_lab <- get_astk(ct_pval)

ylims <- range(datf.mutsm[,"Prol"])
r <- ylims[2] - ylims[1]
ystep <- r * 0.1
ly <- ylims[2]
ly <- ylims[2] + (ystep * 0.5)
py <- ylims[2] + (ystep * 0.75)

ptitle <- "Mutation effects on Prol"
p <- ggplot(datf.mutsm, aes(x = value, y = Prol)) + 
    geom_boxplot(outlier.shape=NA) + 
    geom_dotplot(binaxis='y', stackdir='center', binwidth = .05) +
    geom_segment(aes(x = 1, y = ly, xend = 2, yend = ly)) +
    geom_text(label = pval_lab, x=1.5, y = ly, size=4, vjust=0.5) + 
    ylim(c(ylims[1], ylims[2] + (2*ystep))) +
    facet_wrap(~variable, nrow=1) + 
    labs(x = NULL, y = "TCGA Prol\nproportions", title = ptitle) + 
    th_mp

gr1 <- ggplotGrob(p)

gr1 <- gtable_add_tag(gr1, "a", fs=16, just=c(0,1), l=tagl)
gr1 <- pad_plot(gr1, t=mar_pct, r=mar_pct, b=mar_pct, l=mar_pct)

#========================================================
# cell type prop diff by gene
#========================================================

th_ctpd <- theme_classic() + 
th_txt + 
theme(panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size=1), 
      plot.subtitle = element_text(size = 8, hjust = 0.5, vjust = 0), 
      axis.line = element_blank(),
      strip.text = element_text(size = 8), 
      axis.text = element_text(size = 8), 
      strip.background = element_rect(fill = "white"))


# TP53 and RB1
genesk <- c("TP53", "RB1")
wr.smg.p53 <- wr.smg[wr.smg[,"Gene"] %in% genesk,]
datf.p <- data.frame("CellType" = wr.smg.p53[,"CellType"], 
                     "MeanDiff" = wr.smg.p53[,"t.estimate1"] - wr.smg.p53[,"t.estimate2"], 
                     "log10p" = -log10(wr.smg.p53[,"w.p_adj"]), 
                     "Gene" = wr.smg.p53[,"Gene"])

datf.p[,"ct_label"] <- ""
k <- 10^( -datf.p[,"log10p"] ) < .05
datf.p[k,"ct_label"] <- datf.p[k,"CellType"]
datf.p[,"Gene"] <- factor(datf.p[,"Gene"], levels = genesk)

ptitle <- "Mutation effects per cell-type"
xtitle <- "Cell-type proportion difference\nbetween Mut and WT"
ytitle <- expression(-log[10]~p-value)

gr2_l <- list()
for (i in 1:length(genesk)){
    gene <- genesk[i]
    datf.ps <- datf.p[datf.p[,"Gene"] == gene,]
    p <- ggplot(datf.ps, aes(x = MeanDiff, y = log10p)) + 
    geom_vline(xintercept = 0, col = "red") + 
    geom_point(shape=16) + 
    geom_text_repel(aes(label = ct_label), force = 10, force_pull = 1, 
                    size = 3, color = "red") +
    scale_x_continuous(limits = c(-4,4)) + 
    th_ctpd
    # geom_text(aes(label = CellType), color = "red", hjust = 1.25, vjust = 1.25) + 

    if (i == 1){
        p <- p + labs(x=NULL, y = ytitle, title = ptitle, subtitle = gene)
    } else {
        p <- p + labs(x=xtitle, y = ytitle, subtitle = gene)
    }

    gr_i <- ggplotGrob(p)
    if (i == 1){
        gr_i <- gtable_add_tag(gr_i, "b", fs=16, just=c(0,1), l=tagl)
    }
    gr2_l[[i]] <- gr_i
}


hs <- unit(c(0.5, 0.5), "npc")
gr2 <- stack_gtable_v(gr2_l, heights = hs)
gr2 <- pad_plot(gr2, t=mar_pct, r=mar_pct, b=mar_pct, l=mar_pct)

#========================================================
# TP53 mutation type by Prol expression
#========================================================

th_ty <- theme_classic() + 
th_txt + 
theme(panel.grid = element_blank(),
      axis.ticks.y = element_blank(), 
      axis.text.x = element_text(angle = 30, hjust = 0.9, color = "black"),
      axis.text.y = element_blank())

# factor levels for mutation types
mutmap <- c("None" = "None", 
            "Synonymous" = "Synonymous", 
            "In_Frame" = "In Frame", 
            "Missense" = "Missense", 
            "Splice_Site" = "Splice Site", 
            "Frame_Shift" = "Frame Shift", 
            "Nonsense" = "Nonsense")

ct <- "Prol"
mu <- "TP53"

# get significant mutation types
k <- wr.smg.type[,"CellType"] == ct & wr.smg.type[,"Gene"] == mu
wr.smg.typek <- wr.smg.type[k,]
wr.smg.typek[,"MutType"] <- mutmap[wr.smg.typek[,"MutType"]]
rownames(wr.smg.typek) <- wr.smg.typek[,"MutType"]

# mutation type per sample-gene
fn <- "exp/tcga_hcc/mol/cell_type_main/sample.mut_type.txt"
datf.type <- read.table(fn, header=TRUE, row.names=1, sep='\t')

# get mutations present in gene
mut_types <- unname(mutmap[mutmap %in% datf.type[,mu]])
datf.type[,mu] <- factor(datf.type[,mu], levels = mut_types)

# y limits
qrs <- quantile(datf.type[,ct], c(0.25, 0.5, 0.75)) 
iqr <- qrs[3] - qrs[1]
ylims <- c(qrs[2] - (iqr*1.5), qrs[2] + (iqr*1.5))
ylims <- range(datf.type[,ct])
r <- ylims[2] - ylims[1]
ystep <- r * 0.1
ly <- ylims[2] + (ystep * 0.5)

ptitle <- paste0("LOF mutations in TP53 increase Prol in TCGA")
xlabel <- paste0(mu, " mutation type")
ylabel <- paste0("TCGA Prol\nproportions")
p <- ggplot(datf.type, aes_string(x = mu, y = ct)) + 
    geom_boxplot(outlier.shape=NA) + 
    geom_dotplot(binaxis='y', stackdir='center', binwidth = .05) +
    labs(x = xlabel, y = ylabel, title = ptitle) + 
    th_ty

for (mut_type in mut_types){
    if (! mut_type %in% rownames(wr.smg.typek)) next
    atsk <- get_astk(wr.smg.typek[mut_type,"w.p"])
    if (is.null(atsk)) next
    ly <- ly + ystep
    x1 <- which(mut_types == "None")
    x2 <- which(mut_types == mut_type)
    p <- p + geom_segment(x = x1, y = ly, xend = x2, yend = ly)
    p <- p + geom_text(label = atsk, x=(x1+x2)/2, y = ly, size=4, vjust=0.5) 
}

suppressWarnings(p <- p + ylim(c(ylims[1], ly + ystep)))

gr3 <- ggplotGrob(p)
gr3 <- gtable_add_tag(gr3, "c", fs=16, just=c(0,1), l=tagl)
gr3 <- pad_plot(gr3, t=mar_pct, r=mar_pct, b=mar_pct, l=mar_pct)

#=========================================
# umap and boxplot theme
#=========================================

th_umap <- theme_classic() + 
th_txt + 
theme(legend.key.height = unit(6, "points"),
      legend.key.width = unit(12, "points"),
      legend.text = element_text(size = 6, angle=30, hjust=1), 
      legend.title = element_text(size = 6, vjust = 1, hjust = 0.5), 
      legend.position = "bottom", 
      legend.box.spacing = unit(2, "points"), 
      axis.text=element_blank(),
      axis.text.x = element_blank(), 
      axis.text.y = element_blank(), 
      axis.ticks=element_blank())

th_box <- theme_classic() +
th_txt + 
theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size=8))

#=========================================
# Plot mutation enr score
#=========================================

fn <- "exp/tcga_hcc/sharma_aiz.mut_score/cell_type_main.diff_stat.txt"
diff_stat <- read.table(fn, header=TRUE, row.names=1, sep='\t')

pdatf <- cbind(udf[,c("UMAP1", "UMAP2")], mut_scores[rownames(udf),])
pdatf[,"cell_type_main"] <- factor(pdatf[,"cell_type_main"])

uxmax <- max(udf[,"UMAP1"])
uymin <- min(udf[,"UMAP2"])

g2plot_l <- c("TP53", "RB1")
tags <- list(c("d", "f"), c("e", "g"))
names(tags) <- g2plot_l
utitles <- paste0(g2plot_l, " mut. up-regulated\ngene scores")
ctitles <- utitles
btitles <- paste0(g2plot_l, " mut.\nup-regulated\ngene scores")

gr4_all <- list()
for (i in 1:length(g2plot_l)){
    gene <- g2plot_l[i]
    umap_title <- paste0("Droplet scores of mut. ", gene, "\nup-regulated genes")
    col_title <- paste0(gene, " score")

    k <- diff_stat[,"CellType"] == "Prol" & diff_stat[,"gene"] == gene
    wpadj <- diff_stat[k, "w_p_adj"]
    astk <- get_astk(wpadj)

    # umap
    o <- order(pdatf[,gene], decreasing=FALSE)
    p <- ggplot(pdatf[o,], aes_string(x="UMAP1", y="UMAP2", color=gene)) + 
    geom_point(size = 0.01, shape = 16) + 
    geom_text(label = astk, x = uxmax, y = uymin, 
              color = "black", size = 6, hjust = 1, vjust = 0) + 
    labs(title = utitles[i]) + 
    th_umap + 
    scale_color_gradientn(colours = reds, name = ctitles[i])

    gr4i <- ggplotGrob(p)
    gr4i <- raster_ggpoints(gr4i, w = 3, h = 3)

    # get coordinates for asterisks on Prol
    p_x <- which(levels(pdatf[,"cell_type_main"]) == "Prol")
    yvals <- pdatf[pdatf[,"cell_type_main"] == "Prol", gene]
    yquant <- quantile(yvals, prob = c(.25, .75))
    whisk_max <- yquant[2] + (1.5 * (yquant[2] - yquant[1]))
    whisk_min <- yquant[1] - (1.5 * (yquant[2] - yquant[1]))
    whisk_r <- whisk_max - whisk_min
    p_y <- whisk_max + (whisk_r * .1)

    # box plot
    p <- ggplot(pdatf, aes_string(x = "cell_type_main", y = gene)) + 
    geom_boxplot(outlier.shape=16, outlier.size=0.1, outlier.alpha=0.1) + 
    geom_text(label = astk, x = p_x, y = p_y, color = "black", 
              size = 4, hjust = 0.5, vjust = 0, parse = FALSE) + 
    scale_y_continuous(expand = expansion(mult = c(0.1, .25))) + 
    labs(y = btitles[i], x = NULL) + 
    th_box

    gr4j <- ggplotGrob(p)

    # join
    gr4_l <- list(gr4i, gr4j)
    hs <- unit(c(0.7, 0.3), "npc")
    gtmp <- stack_gtable_v(gr4_l, heights = hs)

    # add tags
    gtmp$grobs[[1]] <- gtable_add_tag(gtmp$grobs[[1]], tags[[i]][1], fs=16, just=c(0,1), l=tagl)
    gtmp$grobs[[2]] <- gtable_add_tag(gtmp$grobs[[2]], tags[[i]][2], fs=16, just=c(0,0), l=tagl)

    gtmp <- pad_plot(gtmp, t=mar_pct, r=mar_pct, b=mar_pct, l=mar_pct)

    gr4_all[[i]] <- gtmp
}

#=========================================
# put together a gtable for plotting everything
#=========================================

gtop <- gtable(widths = unit(c(3.5, 2.5), "inches"), 
               heights = unit(c(1.75, 2.25), "inches"))
gtop <- gtable_add_grob(gtop, gr1, t=1,b=1,l=1,r=1)
gtop <- gtable_add_grob(gtop, gr2, t=1,b=2,l=2,r=2)
gtop <- gtable_add_grob(gtop, gr3, t=2,b=2,l=1,r=1)

gbot <- gtable(widths = unit(c(3, 3), "inches"), 
               heights = unit(4.5, "inches"))
gbot <- gtable_add_grob(gbot, gr4_all[[1]], t=1,b=1,l=1,r=1)
gbot <- gtable_add_grob(gbot, gr4_all[[2]], t=1,b=1,l=2,r=2)

gt <- gtable(widths = unit(6, "inches"), 
             heights = unit(c(4, 4.5), "inches"))
gt <- gtable_add_grob(gt, gtop, t=1,b=1,l=1,r=1)
gt <- gtable_add_grob(gt, gbot, t=2,b=2,l=1,r=1)

pdf(paste0(dir_plt, "fig4.pdf"), width = 6, height = 8.5)
grid.draw(gt)
dev.off()

