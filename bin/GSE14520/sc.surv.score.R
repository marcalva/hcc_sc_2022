
setwd("../../")

library(Seurat)
library(grid)
library(gtable)
library(ggplot2)
library(ggrepel)
library(scales)
library(gridExtra)
library(reshape2)

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

ens_id <- rownames(gene_info)
names(ens_id) <- gsub("\\..*", "", ens_id)

# Set directories
dir_exp <- "exp/GSE14520/sharma_aiz.surv_score/"
dir.create(dir_exp, showWarnings=FALSE, recursive=TRUE)

fn <- "data/processed/sharma_aiz/liver.int_rand.rds"
integrated <- readRDS(fn)

fn <- "exp/GSE14520/survival.gene/GSE14520.cox_ph.all_genes.txt"
surv <- read.table(fn, header=TRUE, sep='\t', row.names=1)
surv[,"ens"] <- symb2ens[rownames(surv)]
k <- !is.na(surv[,"ens"])
surv <- surv[k,]

reds <- read.csv("data/ref/colors/red_colrs.csv", header=FALSE)
reds <- reds[,1]

#=========================================
# score OS and PFI genes
#=========================================

DefaultAssay(integrated) <- "RNA"
surv.k <- surv[surv[,"ens"] %in% rownames(integrated),]

events <- c("Survival")
enames <- c("Survival" = "OS")
for (e in events){
    ccol <- paste0(e, ".exp_coef")
    pcol <- paste0(e, ".p_adj")
    k <- surv.k[,ccol] > 1 & surv.k[,pcol] < .05
    surv_genes <- surv.k[k,"ens"][1:min(1000, sum(k))]
    message(length(surv_genes), " genes with sig HR > 1 for ", e)

    l <- list("e" = surv_genes)
    integrated <- AddModuleScore(integrated, 
                                   features = l, 
                                   name = enames[e])
}

md <- integrated@meta.data

cn <- paste0(enames, "1")
out_fn <- paste0(dir_exp, "sharma_aiz.surv_enr_scores.txt")
write.table(md[,cn,drop=FALSE], out_fn, row.names=TRUE, 
            col.names=NA, quote=FALSE, sep='\t')

#=========================================
# wilcoxon test per cell-type
#=========================================

events <- c("Survival")
enames <- c("Survival" = "OS")
cl_ids <- c("cell_type_main", "cell_type_fine")
for (cl_id in cl_ids){
    cts <- sort(unique(md[,cl_id]))
    tum_diffs <- list()
    for (e in events){
        en <- paste0(e, "1")
        ctdiff_l <- list()
        for (ct in cts){
            en <- paste0(enames[e], "1")
            k <- md[,cl_id] == ct
            x1 <- md[k,en]
            x2 <- md[!k,en]
            wret <- wilcox.test(x = x1, y = x2)
            tret <- t.test(x = x1, y = x2)
            rdatf <- data.frame("CellType" = ct, 
                                "event" = enames[e], 
                                "mean_ct" = mean(x1), 
                                "mean_nonct" = mean(x2), 
                                "t_statistic" = tret$statistic, 
                                "t_p" = tret$p.value, 
                                "w_statistic" = wret$statistic, 
                                "w_p" = wret$p.value)
            ctdiff_l[[ct]] <- rdatf
        }
        rdatf <- do.call(rbind, ctdiff_l)
        # correct for number of cell-types
        rdatf[,"t_p_adj"] <- p.adjust(rdatf[,"t_p"], method = "fdr")
        rdatf[,"w_p_adj"] <- p.adjust(rdatf[,"w_p"], method = "fdr")
        tum_diffs[[e]] <- rdatf
    }
    tum_diffs <- do.call(rbind, tum_diffs)
    out_fn <- paste0(dir_exp, cl_id, ".diff_stat.txt")
    write.table(tum_diffs, out_fn, row.names=TRUE, col.names=NA, 
                quote=FALSE, sep='\t')
}

#=========================================
# common plotting theme
#=========================================

theme_trnsp <- theme(plot.background = element_rect(fill="transparent", color=NA),
                     panel.background = element_rect(fill="transparent", color=NA),
                     gend.background = element_rect(fill="transparent", color=NA))

theme_txt <- theme(text = element_text(size = 8),
                   plot.title = element_text(hjust = 0.5))

theme_leg <- theme(legend.key.height = unit(1, "strheight", "0"),
                   legend.key.width = unit(1, "strwidth", "0"))

theme_axs <- theme(axis.text=element_blank(),
                   axis.ticks=element_blank())

theme_s <- theme_classic() + 
    theme_txt + theme_leg + theme_axs

#=========================================
# merge UMAP and meta data
#=========================================

udf <- integrated@reductions$umap@cell.embeddings
colnames(udf) <- c("UMAP1", "UMAP2")
udf <- cbind(udf, integrated@meta.data[rownames(udf),])

set.seed(1, kind = 'Mersenne-Twister')
udf <- udf[sample(1:nrow(udf)),]

#=========================================
# Plot Tum Enr Score
#=========================================

for (ev in events){
    cn <- paste0(enames[ev], "1")
    o <- order(udf[,cn], decreasing=FALSE)
    p <- ggplot(udf[o,], aes_string(x="UMAP1", y="UMAP2", color=cn)) + 
    geom_point(size = 0.01, shape = 16) + 
    theme_s + 
    scale_color_gradientn(colours = reds, name = enames[ev])

    outfn <- paste0(dir_exp, enames[ev], ".surv_score.umap.pdf")
    ggsave(outfn, width = 3.5, height = 3, dpi=300)
    pname <- gsub("pdf$", "png", outfn)
    ggsave(pname, width = 3.5, height = 3, dpi=300)
}

# box plot of enr scores (main)

for (ev in events){
    cn <- paste0(enames[ev], "1")
    p <- ggplot(udf, aes_string(x = "cell_type_main", y = cn)) + 
    geom_boxplot(outlier.shape=16, outlier.size=0.1, outlier.alpha=0.1) + 
    theme_classic() + 
    ylab("Score") + 
    xlab(NULL) + 
    ggtitle(paste0(enames[ev], " enrichment")) + 
    theme(text = element_text(size = 8),
          axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size=8), 
          plot.title = element_text(hjust = 0.5))

    ggsave(paste0(dir_exp, enames[ev], ".surv_score.main.box.pdf"), width = 3, 
           height = 3, dpi = 300)
}

