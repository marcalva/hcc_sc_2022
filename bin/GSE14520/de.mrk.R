
# investigate tumor DE of cell-type markers

setwd("../../")

library(Seurat)
library(Matrix)
library(ggplot2)

#=====================================================
#=====================================================

dir_exp <- "exp/GSE14520/de/"
dir.create(dir_exp, showWarnings = FALSE, recursive=TRUE)

dir_exp_main <- "exp/GSE14520/de/cell_type_main/"
dir.create(dir_exp_main, showWarnings = FALSE, recursive=TRUE)

dir_exp_fine <- "exp/GSE14520/de/cell_type_fine/"
dir.create(dir_exp_fine, showWarnings = FALSE, recursive=TRUE)

# annotation data
fn <- "data/processed/tcga_hcc/expr/tcga.gencodev26.rds"
gencode <- readRDS(fn)

# get DE gene results
fn <- "exp/GSE14520/de/de.table.tsv.gz"
de <- read.table(fn, header=TRUE, row.names=1)

# get expressed genes
fn <- "data/processed/sharma_aiz/liver.int_rand.rds"
seur <- readRDS(fn)
rc <- seur@assays$RNA@counts
g_exprd <- rownames(rc)[rowSums(rc) > 0]
gk <- g_exprd %in% de[,"ens"]
g_exprd <- g_exprd[gk]

#================================================
# markers
#================================================

cl_ids <- c("cell_type_main", "cell_type_fine")
ct_mrk_l <- list()
mrk_logfc_l <- list()

for (cl_id in cl_ids){
    fn <- paste0("exp/sharma_aiz/markers/markers.", cl_id, ".log_fc.txt")
    main_logfc <- read.table(fn, header=TRUE)
    rownames(main_logfc) <- gsub("\\..*", "", rownames(main_logfc))

    fn <- paste0("exp/sharma_aiz/markers/markers.", cl_id, ".txt")
    ct_mrk <- read.table(fn, header=TRUE)
    ct_mrk[,"gene"] <- gsub("\\..*", "", ct_mrk[,"gene"] )

    # subset markers to those in tcga
    kg <- rownames(main_logfc)  %in% rownames(de)
    main_logfc <- main_logfc[kg,]
    kg <- ct_mrk[,"Name"] %in% rownames(de)
    ct_mrk <- ct_mrk[kg,]

    mrk_logfc_l[[cl_id]] <- main_logfc
    ct_mrk_l[[cl_id]] <- ct_mrk
}

mrk_logfc_thr <- 0.5

#================================================
# summarize marker data with DE stats
#================================================

mrk_de_l <- list()
for (cl_id in cl_ids){
    ct_mrk <- ct_mrk_l[[cl_id]]
    cell_types <- sort(unique(ct_mrk[,"cluster"]))
    d1_l <- list()
    for (ct in cell_types){
        ct_mrk_s <- ct_mrk[ct_mrk[,"cluster"] == ct,,drop=FALSE]
        ct_g <- ct_mrk_s[,"Name"]
        ct_mrk_s[,"de_logFC"] <- de[ct_g, "logFC"]
        ct_mrk_s[,"de_p_adj"] <- de[ct_g, "adj.P.Val"]
        d1_l[[ct]] <- ct_mrk_s
    }
    mrk_de_l[[cl_id]] <- d1_l
}

mrk_de <- list()
for (cl_id in cl_ids){
    mdt <- do.call(rbind, mrk_de_l[[cl_id]])
    mdt[,"gws"] <- " "
    mdt[mdt[,"de_p_adj"] < 0.05,"gws"] <- "adj p < 0.05"
    mrk_de[[cl_id]] <- mdt
    out_fn <- paste0(dir_exp, cl_id, ".de_mrk.txt")
    write.table(mdt, out_fn, row.names=FALSE, col.names=TRUE, 
                quote=FALSE, sep='\t')
}

#================================================
# ggplot theme
#================================================

theme_txt <- theme(text = element_text(colour="black", size=8), 
                   axis.text = element_text(colour="black", size=8), 
                   plot.title = element_text(hjust=0.5))

theme_leg <- theme(legend.position = "bottom", 
                   legend.key.size = unit(1, "strwidth", "a"), 
                   legend.margin = margin(0,2.5,2.5,2.5))
#                    legend.box.margin = margin(0,0,0,0))

#================================================
# plot jitter on top of box of log2FC for marker genes
# per cell-type
#================================================

boxj_plot <- function(datf, x = "de_logFC", y = "cluster", color = "gws", 
                      title = NULL, gmap = NULL){
    # order y
    y_m <- tapply(datf[,x], datf[,y], mean)
    o <- order(y_m, decreasing=FALSE)
    lv <- names(y_m)[o]
    datf[,y] <- factor(datf[,y], levels = lv)

    p <- ggplot(datf, aes_string(x = x, y = y)) + 
    geom_vline(xintercept = 0, color = "black", linetype = "dotted") + 
    geom_boxplot(fill = NA, outlier.shape=NA) + 
    geom_jitter(aes(color = gws), height = 0.2, shape=16, size=0.5) + 
    ylab(NULL) + 
    xlab(bquote("log"[2]*"FC")) + 
    ggtitle(title) + 
    scale_color_manual(values = gmap, breaks = "adj p < 0.05", 
                       name=NULL) + 
    theme_classic() + theme_txt + theme_leg + 
    theme(axis.ticks.y = element_blank(), 
          legend.position=c(1,.2), 
          legend.justification=c(1,0), 
          legend.background = element_rect(color="black", linetype="solid"))
    return(p)
}



# per main cell type
cl_id <- "cell_type_main"
pdatf <- mrk_de[[cl_id]]

# genome-wide significant
gmap <- c(" " = "black", "adj p < 0.05" = "red")

# subset
pdatf <- pdatf[pdatf[,"avg_log2FC"] > mrk_logfc_thr,]

p <- boxj_plot(pdatf, title = paste0("LCI"), gmap = gmap)

out_fn <- paste0(dir_exp_main, cl_id, ".de_logfc.boxj.pdf")
ggsave(out_fn, width = 3.5, height = 3.5)

#================================================
# percent cell-type markers with significant DE
# and directional associations
#================================================

mrk_de_l <- list()

for (cl_id in cl_ids){
    ct_mrk <- ct_mrk_l[[cl_id]]
    # get genes with logFC > 0.5
    ct_mrk_s <- ct_mrk[ct_mrk[,"avg_log2FC"] > mrk_logfc_thr,]

    cell_types <- sort(unique(ct_mrk[,"cluster"]))

    pct_gws_l <- list()
    for (ct in cell_types){
        # get cell type marker genes
        ct_g <- ct_mrk_s[ct_mrk_s[,"cluster"] == ct,"Name"]
        ct_g_n <- length(ct_g)
        pv_col <- paste0("adj.P.Val")
        fc_col <- paste0("logFC")
        sig_k <- de[ct_g, pv_col] < 0.05
        pos_k <- de[ct_g, fc_col] > 0
        neg_k <- de[ct_g, fc_col] < 0
        sig_pos_k <- sig_k & pos_k
        sig_pos_k <- sig_k & pos_k
        sig_neg_k <- sig_k & neg_k
        sig_pct <- sum(sig_k) / ct_g_n
        pos_pct <- sum(pos_k) / ct_g_n
        neg_pct <- sum(neg_k) / ct_g_n
        sig_pos_pct <- sum(sig_pos_k) / ct_g_n
        sig_neg_pct <- sum(sig_neg_k) / ct_g_n
        ret <- data.frame("CellType" = ct, 
                          "pct_gws" = sig_pct, 
                          "pct_pos" = pos_pct, 
                          "pct_neg" = neg_pct, 
                          "pct_gws_pos" = sig_pos_pct, 
                          "pct_gws_neg" = sig_neg_pct,
                          "ct_n" = ct_g_n)
        pct_gws_l[[ct]] <- ret
    }
    pct_gws_datf <- do.call(rbind, pct_gws_l)
    mrk_de_l[[cl_id]] <- pct_gws_datf
}

out_fn <- paste0(dir_exp_main, "cell_type_main.de.mrk_stat.txt")
write.table(mrk_de_l[["cell_type_main"]], out_fn, row.names=FALSE, col.names=TRUE, 
            quote=FALSE, sep='\t')
out_fn <- paste0(dir_exp_fine, "cell_type_fine.de.mrk_stat.txt")
write.table(mrk_de_l[["cell_type_fine"]], out_fn, row.names=FALSE, col.names=TRUE, 
            quote=FALSE, sep='\t')

