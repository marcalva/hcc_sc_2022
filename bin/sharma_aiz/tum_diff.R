
# Test differences between tumor and non-tumor

# Run differential expression between proliferating cells 
# and non-proliferating cells for each major cell-type.
# For example, run DE between proliferating hepatocytes 
# and non-proliferating hepatocytes.

setwd("../../")

library(Seurat)
library(Matrix)
library(ggplot2)
library(ggrepel)
library(scales)
library(grid)
library(gridExtra)
library(reshape2)
source("scripts/diff_test_ml.R")

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

fn <- "data/processed/sharma_aiz/liver.int_rand.rds"
integrated <- readRDS(fn)
md <- integrated@meta.data

dir_exp <- "exp/sharma_aiz/tum_diff/";
dir.create(dir_exp, showWarnings=FALSE, recursive=TRUE)

# colors
tum_df <- read.csv("data/ref/colors/tumor_colrs.csv", row.names=1, header=FALSE)
tum_col <- tum_df[,1]; names(tum_col) <- rownames(tum_df)

cl_ids <- c("cell_type_fine", "cell_type_main")

#=========================================
# Wilcoxon test of cell-type proportions
# between tumor and non-tumor
#=========================================

md[,"SampleStat"] <- paste(md[,"PatientStat"], md[,"TumorStat"], sep="_")
sams <- sort(unique(md[,"SampleStat"]))
k_sam <- setdiff(sams, c("Aizarani_NonTumor", "Sharma_0_NonTumor"))
k <- md[,"SampleStat"] %in% k_sam
mdk <- md[k,] # excludes Aizarani_NonTumor and Sharma_0_NonTumor

ct_df_m_p <- diff_test_ml(x = mdk[,"SampleStat"], 
                          g = mdk[,"TumorStat"], 
                          z = mdk[,"cell_type_main"], 
                          paired = TRUE, 
                          p = mdk[,"PatientStat"])

ct_df_m <- diff_test_ml(x = mdk[,"SampleStat"], 
                        g = mdk[,"TumorStat"], 
                        z = mdk[,"cell_type_main"]) 

ct_df_f_p <- diff_test_ml(x = mdk[,"SampleStat"], 
                          g = mdk[,"TumorStat"], 
                          z = mdk[,"cell_type_fine"], 
                          paired = TRUE, 
                          p = mdk[,"PatientStat"])

ct_df_f <- diff_test_ml(x = mdk[,"SampleStat"], 
                        g = mdk[,"TumorStat"], 
                        z = mdk[,"cell_type_fine"]) 

diffs_l <- list("cell_type_main" = ct_df_m, "cell_type_main_p" = ct_df_m_p,
     "cell_type_fine" = ct_df_f, "cell_type_fine_p" = ct_df_f_p)

for (n in names(diffs_l)){
    diffs_l[[n]][,"t_p_adj"] <- p.adjust(diffs_l[[n]][,"t_p"], method="fdr")
    diffs_l[[n]][,"w_p_adj"] <- p.adjust(diffs_l[[n]][,"w_p"], method="fdr")
}

for (n in names(diffs_l)){
    fn <- paste0(dir_exp, n, ".tum_diff.txt")
    write.table(diffs_l[[n]], fn, row.names=TRUE, col.names=NA, 
                quote=FALSE, sep='\t')
}

#=========================================
# get cell-type frequency of Sample
#=========================================

# map sample to tumor type
s_u <- unique(md[,"SampleStat"])
sam2tum <- sapply(s_u, function(s) as.character(md[md[,"SampleStat"] == s, "TumorStat"][1]))

prop_l <- list()
prop_lm <- list()
for (cl_id in cl_ids){
    p <- table(mdk[,"SampleStat"], mdk[,cl_id])
    p <- sweep(p, 1, rowSums(p), '/')
    prop_l[[cl_id]] <- p
    pm <- melt(p)
    colnames(pm) <- c("Sample", "CellType", "Prop")
    pm[,"Tumor"] <- sam2tum[as.character(pm[,"Sample"])]
    prop_lm[[cl_id]] <- pm

    out_fn <- paste0(dir_exp, cl_id, ".sample_prop.txt")
    write.table(pm, out_fn, row.names=FALSE, col.names=TRUE,
                quote=FALSE, sep='\t')
}

#=========================================
# proportion of tumor and non-tumor per cell-type
#=========================================

tprop_l <- list()
tprop_lm <- list()

for (cl_id in cl_ids){
    p <- table(md[,cl_id], md[,"TumorStat"])
    p <- sweep(p, 1, rowSums(p), '/')
    pm <- melt(p)
    colnames(pm) <- c("CellType", "Tumor", "Prop")
    tprop_l[[cl_id]] <- p
    tprop_lm[[cl_id]] <- pm

    out_fn <- paste0(dir_exp, cl_id, ".tum_prop.txt")
    write.table(pm, out_fn, row.names=FALSE, col.names=TRUE,
                quote=FALSE, sep='\t')
}

#=========================================
# plotting theme for box
#=========================================

theme_txt <- theme(text = element_text(size = 8),
                   plot.title = element_text(hjust = 0.5))

theme_leg <- theme(legend.key.height = unit(1, "strheight", "0"),
                   legend.key.width = unit(1, "strwidth", "0"), 
                   legend.margin = margin(0,0,0,0, unit='pt'), 
                   legend.position = "bottom")

theme_ax <- theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, angle=90))

theme_box <- theme_classic() + 
    theme_txt + theme_leg + theme_ax

#=========================================
# box plot of cell-type props separated by tumor status
#=========================================

p <- ggplot(prop_lm[["cell_type_fine"]], 
            aes(x = CellType, y = Prop, fill = Tumor)) + 
geom_boxplot(size = .2, outlier.shape=NA) + 
labs(x=NULL, y=NULL) + 
scale_fill_manual(values = tum_col) + 
theme_box

out_fn <- paste0(dir_exp, "cell_type_fine.prop.tum.box.pdf")
ggsave(out_fn, width=5, height=3)

p <- ggplot(prop_lm[["cell_type_main"]], aes(x = CellType, y = Prop, fill = Tumor)) + 
geom_boxplot(size = .2, outlier.shape=NA) + 
labs(x=NULL, y=NULL) + 
scale_fill_manual(values = tum_col) + 
theme_box

out_fn <- paste0(dir_exp, "cell_type_main.prop.tum.box.pdf")
ggsave(out_fn, width=3, height=3)

#=========================================
# plotting theme for bar plot
#=========================================

theme_txt <- theme(text = element_text(size = 8),
                   plot.title = element_text(hjust = 0.5))

theme_leg <- theme(legend.key.height = unit(1, "strheight", "0"),
                   legend.key.width = unit(1, "strwidth", "0"), 
                   legend.position = "bottom", 
                   legend.margin = margin(0,0,0,0, unit='pt'))
theme_ax <- theme(axis.line.x = element_blank())

theme_p <- theme_classic() + 
    theme_txt + theme_leg + theme_ax

#=========================================
# bar plot of tumor prop per cell-type
#=========================================

md <- integrated@meta.data
tum_prop_f <- tprop_l[["cell_type_fine"]]
tum_prop_fm <- tprop_lm[["cell_type_fine"]]

# get significance
cts <- rownames(diffs_l[["cell_type_fine_p"]])
sig <- rep("p >= 0.05", length.out=length(cts))
names(sig) <- cts
sig[diffs_l[["cell_type_fine_p"]][,"w_p_adj"] < .05] <- "p < 0.05"

# order by decreasing prop
o <- order(tum_prop_f[,"Tumor"])
lv <- rownames(tum_prop_f)[o]
tum_prop_fm[,"CellType"] <- factor(tum_prop_fm[,"CellType"], levels=lv)
tum_prop_fm[,"sig"] <- sig[as.character(tum_prop_fm[,"CellType"])]
al <- c("p < 0.05" = 1, "p >= 0.05" = 0.4)

# plot cell-type tumor prop
p <- ggplot(tum_prop_fm, aes(x = CellType, y = Prop, fill = Tumor)) + 
    geom_bar(aes(alpha = sig), position="fill", stat="identity") + 
    coord_flip() + 
    labs(y="Proportion",x=NULL) + 
    scale_y_continuous(expand = expansion(0)) + 
    scale_fill_manual(values=tum_col, name=NULL) + 
    scale_alpha_manual(values = al, name = NULL) + 
    theme_p

fn <- paste0(dir_exp, "cell_type_fine.tum_prop.pdf")
ggsave(fn, width = 3, height = 3)

#=========================================
# zoom in bar plots for significant cell-types
#=========================================

# plotting theme for box
th_bx <- theme(text = element_text(size = 8),
               axis.text = element_text(colour = "black", size = 6), 
               axis.line = element_line(size = rel(0.5)), 
               plot.title = element_text(hjust = 0.5), 
               legend.position = "none", 
               axis.text.x = element_text(hjust = 0.5, vjust = 0.5, angle=0))
th_bx <- theme_classic() + th_bx

# fine cell-types
# get significant cell-types
diffs_cl <- diffs_l[["cell_type_fine_p"]]
cts <- rownames(diffs_cl)
ksig <- diffs_cl[,"w_p_adj"] < .05
cts_sig <- cts[ksig]

props_fine <- prop_lm[["cell_type_fine"]]
for (ct in cts_sig){
    message(ct)
    props_ct <- props_fine[props_fine[,"CellType"] == ct,,drop=FALSE]

    p <- ggplot(props_ct, 
                aes(x = Tumor, y = Prop, fill = Tumor)) + 
    geom_boxplot(size = .2, outlier.shape=NA) + 
    geom_jitter(width = 0.1, size = 0.5) + 
    labs(x=NULL, y="Proportion") + 
    ggtitle(ct) + 
    scale_fill_manual(values = tum_col) + 
    th_bx

    out_fn <- paste0(dir_exp, ct, ".tum_nontum.box.pdf")
    ggsave(out_fn, width = 2, height = 2)
}

# get significant cell-types
diffs_cl <- diffs_l[["cell_type_main_p"]]
cts <- rownames(diffs_cl)
ksig <- diffs_cl[,"w_p_adj"] < .05
cts_sig <- cts[ksig]

props_main <- prop_lm[["cell_type_main"]]
for (ct in cts_sig){
    message(ct)
    props_ct <- props_main[props_main[,"CellType"] == ct,,drop=FALSE]

    p <- ggplot(props_ct, 
                aes(x = Tumor, y = Prop, fill = Tumor)) + 
    geom_boxplot(size = .2, outlier.shape=NA) + 
    geom_jitter(width = 0.1, size = 0.5) + 
    labs(x=NULL, y="Proportion") + 
    ggtitle(ct) + 
    scale_fill_manual(values = tum_col) + 
    th_bx

    out_fn <- paste0(dir_exp, ct, ".tum_nontum.box.pdf")
    ggsave(out_fn, width = 2, height = 2)
}

