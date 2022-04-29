
# investigate survival effect of cell-type markers

setwd("../../")

library(Seurat)
library(Matrix)
library(survival)
library(ggplot2)
library(ggfortify)
library(scales)
library(reshape2)
source("scripts/surv_df.R")
source("scripts/color_pal.R")

#=====================================================
#=====================================================

dir_exp <- "exp/GSE14520/survival.gene/"
dir.create(dir_exp, showWarnings = FALSE, recursive=TRUE)

dir_exp_main <- "exp/GSE14520/survival.gene/cell_type_main/"
dir.create(dir_exp_main, showWarnings = FALSE, recursive=TRUE)

dir_exp_fine <- "exp/GSE14520/survival.gene/cell_type_fine/"
dir.create(dir_exp_fine, showWarnings = FALSE, recursive=TRUE)

# annotation data
fn <- "data/processed/tcga_hcc/expr/tcga.gencodev26.rds"
gencode <- readRDS(fn)

# hazard ratios
fn <- paste0(dir_exp, "GSE14520.cox_ph.all_genes.txt")
all_gene_hr <- read.table(fn, header=TRUE, row.names=1)

fn <- paste0(dir_exp, "GSE14520.cox_ph.stage_all_genes.txt")
stg_all_gene_hr <- read.table(fn, header=TRUE, row.names=1)

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
    kg <- rownames(main_logfc)  %in% rownames(all_gene_hr)
    main_logfc <- main_logfc[kg,]
    kg <- ct_mrk[,"Name"] %in% rownames(all_gene_hr)
    ct_mrk <- ct_mrk[kg,]

    mrk_logfc_l[[cl_id]] <- main_logfc
    ct_mrk_l[[cl_id]] <- ct_mrk
}

mrk_logfc_thr <- 0.5

#================================================
# summarize marker data with survival stats
#================================================

mrk_surv_l <- list()
for (cl_id in cl_ids){
    ct_mrk <- ct_mrk_l[[cl_id]]
    cell_types <- sort(unique(ct_mrk[,"cluster"]))
    events <- c("Survival", "Recurr")
    d1_l <- list()
    for (ct in cell_types){
        d2_l <- list()
        for (e in events){
            ct_mrk_s <- ct_mrk[ct_mrk[,"cluster"] == ct,,drop=FALSE]
            ct_g <- ct_mrk_s[,"Name"]
            hr_col <- paste0(e, ".exp_coef")
            p_col <- paste0(e, ".p_adj")
            ct_mrk_s[,"HR"] <- all_gene_hr[ct_g, hr_col]
            ct_mrk_s[,"HR_p_adj"] <- all_gene_hr[ct_g, p_col]
            ct_mrk_s[,"surv_type"] <- e
            d2_l[[e]] <- ct_mrk_s
        }
        d1_l[[ct]] <- d2_l
    }
    mrk_surv_l[[cl_id]] <- d1_l
}

mrk_surv <- list()
for (cl_id in cl_ids){
    datf_i <- list()
    for (i in names(mrk_surv_l[[cl_id]])){
        datf_i[[i]] <- do.call(rbind, mrk_surv_l[[cl_id]][[i]])
    }
    mrk_surv[[cl_id]] <- do.call(rbind, datf_i)
}

out_fn <- paste0(dir_exp, "cox_ph.ct_mrk.rds")
saveRDS(mrk_surv, out_fn)

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
# plot jitter on top of box of log2HR for marker genes
# per cell-type
#================================================

boxj_plot <- function(datf, x = "log2HR", y = "cluster", color = "gws", 
                      title=NULL){
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
    xlab(bquote("log"[2]*"HR")) + 
    ggtitle(title) + 
    scale_color_manual(values = gmap, breaks = "adj p < 0.05", 
                       name=NULL) + 
    theme_classic() + theme_txt + theme_leg + 
    theme(axis.ticks.y = element_blank(), 
      legend.position=c(1,.05), 
      legend.justification=c(1,0), 
      legend.background = element_rect(color="black", linetype="solid"))
    return(p)
}

# per main cell type and OS
cl_id <- "cell_type_main"
e <- "Survival"
ename <- "OS"
pdatf <- mrk_surv[[cl_id]]
pdatf <- pdatf[pdatf[,"surv_type"] == e,]
pdatf[,"log2HR"] <- log2(pdatf[,"HR"])

# genome-wide significant
pdatf[,"gws"] <- " "
pdatf[pdatf[,"HR_p_adj"] < 0.05,"gws"] <- "adj p < 0.05"
gmap <- c(" " = "black", "adj p < 0.05" = "red")

# subset
pdatf <- pdatf[pdatf[,"avg_log2FC"] > mrk_logfc_thr,]

p <- boxj_plot(pdatf, title = paste0("LCI ", ename))

out_fn <- paste0(dir_exp_main, cl_id, ".", e, ".surv_hr.boxj.pdf")
ggsave(out_fn, width = 3.5, height = 3.5)

# per main cell type and RCR
cl_id <- "cell_type_main"
e <- "Recurr"
ename <- "RCR"
pdatf <- mrk_surv[[cl_id]]
pdatf <- pdatf[pdatf[,"surv_type"] == e,]
pdatf[,"log2HR"] <- log2(pdatf[,"HR"])

# genome-wide significant
pdatf[,"gws"] <- " "
pdatf[pdatf[,"HR_p_adj"] < 0.05,"gws"] <- "adj p < 0.05"
gmap <- c(" " = "black", "adj p < 0.05" = "red")

# subset
pdatf <- pdatf[pdatf[,"avg_log2FC"] > mrk_logfc_thr,]

p <- boxj_plot(pdatf, title = paste0("LCI ", e))

out_fn <- paste0(dir_exp_main, cl_id, ".", e, ".surv_hr.boxj.pdf")
ggsave(out_fn, width = 3.5, height = 3.5)

#================================================
# percent cell-type markers with significant hazard ratio
# and directional associations
# FET of positive/negative HR
#================================================

mrk_stat_l <- list()

for (cl_id in cl_ids){
    ct_mrk <- ct_mrk_l[[cl_id]]
    # get genes with logFC > 0.5
    ct_mrk_s <- ct_mrk[ct_mrk[,"avg_log2FC"] > mrk_logfc_thr,]

    cell_types <- sort(unique(ct_mrk[,"cluster"]))
    events <- c("Survival", "Recurr")

    pct_gws_l <- list()
    for (ct in cell_types){
        # get cell type marker genes
        ct_g <- ct_mrk_s[ct_mrk_s[,"cluster"] == ct,"Name"]
        ct_g_n <- length(ct_g)
        e_l <- list()
        for (e in events){
            pv_col <- paste0("", e, ".p_adj")
            cf_col <- paste0("", e, ".exp_coef")
            sig_k <- all_gene_hr[ct_g, pv_col] < 0.05
            pos_k <- all_gene_hr[ct_g, cf_col] > 1
            neg_k <- all_gene_hr[ct_g, cf_col] < 1
            sig_pos_k <- sig_k & pos_k
            sig_pos_k <- sig_k & pos_k
            sig_neg_k <- sig_k & neg_k
            sig_pct <- sum(sig_k) / ct_g_n
            pos_pct <- sum(pos_k) / ct_g_n
            neg_pct <- sum(neg_k) / ct_g_n
            sig_pos_pct <- sum(sig_pos_k) / ct_g_n
            sig_neg_pct <- sum(sig_neg_k) / ct_g_n
            ret <- data.frame("CellType" = ct, 
                              "event" = e, 
                              "pct_gws" = sig_pct, 
                              "pct_pos" = pos_pct, 
                              "pct_neg" = neg_pct, 
                              "pct_gws_pos" = sig_pos_pct, 
                              "pct_gws_neg" = sig_neg_pct,
                              "ct_n" = ct_g_n)
            e_l[[e]] <- ret
        }
        pct_gws_l[[ct]] <- do.call(rbind, e_l)
    }
    pct_gws_datf <- do.call(rbind, pct_gws_l)
    mrk_stat_l[[cl_id]] <- pct_gws_datf
}

out_fn <- paste0(dir_exp_main, "cell_type_main.surv.mrk_stat.txt")
write.table(mrk_stat_l[["cell_type_main"]], out_fn, row.names=FALSE, col.names=TRUE, 
            quote=FALSE, sep='\t')
out_fn <- paste0(dir_exp_fine, "cell_type_fine.surv.mrk_stat.txt")
write.table(mrk_stat_l[["cell_type_fine"]], out_fn, row.names=FALSE, col.names=TRUE, 
            quote=FALSE, sep='\t')

#================================================
# plot pct GWS markers
#================================================

th_bar <- theme_classic() + 
theme(panel.grid.major.y = element_line(), 
      plot.title = element_text(size = 12, hjust = 0.5), 
      axis.text.x = element_text(size = 8, color="black", angle=90, hjust=1, vjust=0.5), 
      axis.text.y = element_text(size = 8, color="black"), 
      legend.key.size = unit(12, "points"), 
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8))

evmap <- c("Survival" = "OS", "Recurr" = "RCR")
events <- c("OS", "RCR")

for (eve in events){
    pdatf <- mrk_stat_l[["cell_type_main"]]
    pdatf[,"event"] <- evmap[pdatf[,"event"]]
    pdatf <- pdatf[pdatf$event == eve,]

    pdatf[,"pct_neg"] <- -1 * pdatf[,"pct_neg"]
    pdatf[,"pct_gws_neg"] <- -1 * pdatf[,"pct_gws_neg"]
    pdatfs <- pdatf[,-which(colnames(pdatf) == "pct_gws")]
    pdatfs <- pdatfs[,-which(colnames(pdatfs) == "ct_n")]

    pdatfsm <- melt(pdatfs, id.vars = c("CellType", "event"))

    pdatfsm[,"Direction"] <- "Positive"
    pdatfsm[grep("neg", pdatfsm[,"variable"]), "Direction"] <- "Negative"
    pdatfsm[,"GWS"] <- "Adj p >= 0.05"
    pdatfsm[grep("gws", pdatfsm[,"variable"]), "GWS"] <- "Adj p < 0.05"

    type_colr <- hcl_pal(3, offset = 260, chr = c(140, 140), lum = c(60,60), 
                         rand=FALSE)[1:2]
    names(type_colr) <- c("Positive", "Negative")
    type_alph <- c("Adj p < 0.05" = 1, "Adj p >= 0.05" = 0.2)

    coltitle <- paste0("HR direction\nfor ", eve)
    alptitle <- NULL
    ytitle <- "Proportion of cell-type markers"
    xtitle <- "Main cell-types"
    ptitle <- paste0("Direction and proportion of\ncell-type markers affecing ", eve)

    p <- ggplot(pdatfsm, aes(x = CellType, y = value)) +
    geom_bar(aes(fill = Direction, alpha = GWS), position = "identity", stat = "identity") + 
    geom_hline(yintercept = 0, linetype = "dashed") + 
    scale_fill_manual(values = type_colr, name = coltitle) + 
    scale_alpha_manual(values = type_alph, name = alptitle) + 
    scale_y_continuous(labels = scales::percent, limits = c(-1,1)) + 
    labs(y = ytitle, x = xtitle, title = ptitle) + 
    th_bar

    out_fn <- paste0(dir_exp_main, "cell_type_main.surv.", eve, ".mrk_bar.pdf")
    ggsave(out_fn, width = 4, height = 3)
}

