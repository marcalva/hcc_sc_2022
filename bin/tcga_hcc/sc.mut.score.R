
# score droplets for expression of genes that are positively 
# associated with mutations.

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
dir_exp <- "exp/tcga_hcc/sharma_aiz.mut_score/";
dir.create(dir_exp, showWarnings=FALSE, recursive=TRUE)

fn <- "data/processed/sharma_aiz/liver.int_rand.rds"
integrated <- readRDS(fn)
cl_ids <- c("cell_type_main", "cell_type_fine")

reds <- read.csv("data/ref/colors/red_colrs.csv", header=FALSE)
reds <- reds[,1]

fn <- "data/processed/tcga_hcc/mut/tcga.lihc.mut.gene.any.txt"
mut.any <- read.table(fn, header = TRUE, stringsAsFactors = FALSE, 
                      sep = "\t", row.names=1, check.names=FALSE)

mut.any[mut.any >= 1] <- 1

o <- order(rowSums(mut.any), decreasing=TRUE)
genes <- rownames(mut.any)[o][1:25]

# read genes
genes <- c("TP53", "RB1")
mut_de_l <- list()
for (gene in genes){
    fn <- paste0("exp/mut_expr/de/", gene, ".de.table.tsv.gz")
    datf <- read.table(fn, header=TRUE)
    ck <- c("Name", "gene_id", "logFC", "p_adj")
    mut_de_l[[gene]] <- datf[,ck]
}

#=========================================
# score tum-enr genes
#=========================================

DefaultAssay(integrated) <- "RNA"

feats_l <- list()
for (gene in genes){
    de <- mut_de_l[[gene]]
    de[,"gene_id"] <- ens_id[de[,"gene_id"]]
    de.k <- de[de[,"gene_id"] %in% rownames(integrated),]
    k <- de.k[,"logFC"] > 0.5 & de.k[,"p_adj"]  < 0.05
    de_genes <- de.k[k,"gene_id"]
    feats_l[[gene]] <- de_genes
}

names(feats_l) <- genes
n_de_mut <- sapply(feats_l, length)
gc()
integrated <- AddModuleScore(integrated,
                             features = feats_l, 
                             name = "de_genes")
md <- integrated@meta.data
ck <- grep("de_genes", colnames(md))
md[,genes] <- md[,ck]
md <-md[,-ck]
integrated@meta.data <- md

mdr <- md[,c(cl_ids, genes) ,drop=FALSE]


out_fn <- paste0(dir_exp, "sharma_aiz.mut_scores.txt.gz")
z <- gzfile(out_fn, "w")
write.table(mdr, z, row.names=TRUE, 
            col.names=NA, quote=FALSE, sep='\t')
close(z)


#=========================================
# wilcoxon test per cell-type
#=========================================

# dir_out <- paste0(dir_exp, "diff_stat/")
# dir.create(dir_out, showWarnings=FALSE, recursive=TRUE)

md <- integrated@meta.data
k <- md[,"PatientStat"] != "Aizarani" & md[,"PatientStat"] != "Sharma_0"
mdk <- md[k,]

cl_ids <- c("cell_type_main", "cell_type_fine")
for (cl_id in cl_ids){
    stats_l <- list()
    for (gene in genes){
        cts <- sort(unique(mdk[,cl_id]))
        tum_diffs <- list()
        for (ct in cts){
            k <- mdk[,cl_id] == ct
            x1 <- mdk[k,gene]
            x2 <- mdk[!k,gene]
            wret <- wilcox.test(x = x1, y = x2)
            tret <- t.test(x = x1, y = x2)

            # linear model
            g <- rep(0, length(k))
            g[k] <- 1
            lmret <- summary(lm(mdk[,gene] ~ mdk[,"PatientStat"] + g))
            lm_est <- lmret$coefficients[nrow(lmret$coefficients), "Estimate"]
            lm_p <- lmret$coefficients[nrow(lmret$coefficients), "Pr(>|t|)"]

            rdatf <- data.frame("gene" = gene, 
                                "CellType" = ct, 
                                "mean_ct" = mean(x1), 
                                "mean_nonct" = mean(x2), 
                                "t_statistic" = tret$statistic, 
                                "t_p" = tret$p.value, 
                                "w_statistic" = wret$statistic, 
                                "w_p" = wret$p.value, 
                                "lm_est" = lm_est, 
                                "lm_p" = lm_p)
            tum_diffs[[ct]] <- rdatf
        }
        stats_l[[gene]] <- do.call(rbind, tum_diffs)
    }
    stats <- do.call(rbind, stats_l)
    stats[,"w_p_adj"] <- p.adjust(stats[,"w_p"], method="fdr")
    stats[,"t_p_adj"] <- p.adjust(stats[,"t_p"], method="fdr")
    stats[,"lm_p_adj"] <- p.adjust(stats[,"lm_p"], method="fdr")

    o <- order(stats[,"gene"], stats[,"lm_est"], decreasing=c(TRUE))
    stats <- stats[o,]

    out_fn <- paste0(dir_exp, cl_id, ".diff_stat.txt")
    write.table(stats, out_fn, row.names=TRUE, col.names=NA, 
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
udf <- cbind(udf, mdr[rownames(udf),])

set.seed(1, kind = 'Mersenne-Twister')
udf <- udf[sample(1:nrow(udf)),]

#=========================================
# Plot Tum Enr Score
#=========================================

for (gene in genes){
    ngenechr <- as.character(n_de_mut[gene])
    o <- order(udf[,gene], decreasing=FALSE)
    p <- ggplot(udf[o,], aes_string(x="UMAP1", y="UMAP2", color=gene)) + 
    geom_point(size = 0.01, shape = 16) + 
    ggtitle(paste0(gene, " enrichment (", ngenechr, " genes)")) + 
    theme_s + 
    scale_color_gradientn(colours = reds, name = gene)

    outfn <- paste0(dir_exp, gene, ".mut_score.umap.pdf")
    # ggsave(outfn, width = 3.5, height = 3, dpi=300)
    pname <- gsub("pdf$", "png", outfn)
    ggsave(pname, width = 3.5, height = 3, dpi=300)
}

# box plot of enr scores (main)

for (gene in genes){
    ngenechr <- as.character(n_de_mut[gene])
    p <- ggplot(udf, aes_string(x = "cell_type_main", y = gene)) + 
    geom_boxplot(outlier.shape=16, outlier.size=0.1, outlier.alpha=0.1) + 
    theme_classic() + 
    ylab("Score") + 
    xlab(NULL) + 
    ggtitle(paste0(gene, " enrichment (", ngenechr, " genes)")) + 
    theme(text = element_text(size = 8),
          axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size=8), 
          plot.title = element_text(hjust = 0.5))

    ggsave(paste0(dir_exp, gene, ".mut_score.main.box.pdf"), width = 3, 
           height = 3, dpi = 300)
}

