
# Plots for relationship between cell type props and molecular features

setwd("../../")

library(NMF)
library(ggplot2)
library(ggrepel)
# library(egg)
library(gridExtra)
library(reshape2)
library(edgeR)
library(Biobase)
library(plyr)
library(RColorBrewer)
source("scripts/sc_func.R")
# source("scripts/reference_free.R")
# library(Bisque)

#========================================================
# Functions
#========================================================

scale_clip <- function(datf, z = 3){
    y <- apply(datf, 1, function(x) (x - mean(x))/sd(x))
    y <- t(y)
    colnames(y) <- colnames(datf)
    y[y < -z] <- -z
    y[y > z] <- z
    return(y)
}

options(drop=FALSE)

#========================================================
#========================================================

ct_id <- "cell_type_main"

fn <- "data/processed/tcga_hcc/expr/tcga.gencodev26.rds"
gencode <- readRDS(fn)

# proportions
fn <- paste0("data/processed/tcga_hcc/ctp/tcga.TMM.", ct_id, ".decomp.rds")
ct.md <- readRDS(fn)
ct.mdp <- ct.md$bulk.props

for (ct in names(ct.md$genes.used)){
    ct.md$genes.used[[ct]] <- gencode[ct.md$genes.used[[ct]], "Name"]
}

# pheno data
fn <- "data/processed/tcga_hcc/sample/cases.hcc.361.merged.rds"
cases <- readRDS(fn)
fn <- "data/processed/tcga_hcc/sample/samples.hcc.410.merged.rds"
samples <- readRDS(fn)
fn <- "data/processed/tcga_hcc/sample/tcga.bio_sample.rds"
bio_sample <- readRDS(fn)


# tum v non-tum paired test results
fn <- paste0("exp/tcga_hcc/ctp.cell_type_main/tcga.", ct_id, ".tum_nontum.txt")
paired.test <- read.table(fn, header=TRUE, row.names=1, stringsAsFactors=FALSE, sep="\t")

# CNVs
fn <- "data/raw/tcga_hcc/gdac.LIHC.aggregate/LIHC-TP.samplefeatures.txt"
# cnv <- read.table(fn, header = TRUE, stringsAsFactors = FALSE, sep = "\t", row.names=1)

fn <- "data/processed/tcga_hcc/cnv/tcga.lihc.thr.del.354_tum.txt"
dels <- read.table(fn, header = TRUE, stringsAsFactors = FALSE, 
                   sep = "\t", row.names=1, check.names=FALSE)
dels <- apply(dels, c(1,2), function(x) if (x >= 1) "Del" else "None")
fn <- "data/processed/tcga_hcc/cnv/tcga.lihc.thr.amp.354_tum.txt"
amps <- read.table(fn, header = TRUE, stringsAsFactors = FALSE, 
                   sep = "\t", row.names=1, check.names=FALSE)
amps <- apply(amps, c(1,2), function(x) if (x >= 1) "Amp" else "None")

cnvs <- rbind(amps, dels)
# cnvs[cnvs >= 1] <- 1

fn <- "data/processed/tcga_hcc/cnv/tcga.lihc.cnv_descr.txt"
cnv_desc <- read.table(fn, header = TRUE, stringsAsFactors = FALSE, 
                       sep = "\t", row.names=1, check.names=FALSE)

# somatic m
fn <- "data/processed/tcga_hcc/mut/tcga.lihc.mut.gene.any.txt"
mut.any <- read.table(fn, header = TRUE, stringsAsFactors = FALSE, 
                      sep = "\t", row.names=1, check.names=FALSE)
fn <- "data/processed/tcga_hcc/mut/tcga.lihc.mut.gene.type.txt"
mut.type <- read.table(fn, header = TRUE, stringsAsFactors = FALSE, 
                       sep = "\t", row.names=1, check.names=FALSE)

mut.any[mut.any >= 1] <- 1


dir_exp <- paste0("exp/tcga_hcc/mol/", ct_id, "/")
dir.create(dir_exp, showWarnings = FALSE, recursive = TRUE)

options(stringsAsFactors=FALSE)

# somatic mutation results
fn <- paste0("exp/tcga_hcc/mol/", ct_id, ".SMG_mutsig.2CV.ctp.txt")
wr.smg <- read.table(fn, header=TRUE, sep="\t")

# somatic mutation results
fn <- paste0("exp/tcga_hcc/mol/", ct_id, ".SMG_mutsig.2CV.mut_type.ctp.txt")
wr.smg.type <- read.table(fn, header=TRUE, sep="\t")

# cnv results
fn <- paste0("exp/tcga_hcc/mol/", ct_id, ".tcga.lihc.cnv.ctp.txt")
wr.cnv <- read.table(fn, header=TRUE, sep="\t")

# colors
fn <- "data/ref/colors/rd_bu_div.csv"
rd_bu <- read.table(fn, header=FALSE)
rd_bu <- rd_bu[,1]

#========================================================
# Map IDs, get tumor samples
#========================================================

# sample map
fn <- "data/processed/tcga_hcc/sample/sam_case_map.rds"
smaps <- readRDS(fn)
cid2sid_tum <- smaps[["cid2sid_tum"]]
cid2sid_nt <- smaps[["cid2sid_nt"]]
cid2sid_tum_pair <- smaps[["cid2sid_tum_pair"]]
cid2sid_nt_pair <- smaps[["cid2sid_nt_pair"]]
sid_all <- c(cid2sid_tum, cid2sid_nt)
sid_both <- c(cid2sid_tum_pair, cid2sid_nt_pair)

ct.mdp.pt <- ct.mdp[,cid2sid_tum]

# order by cell type
cto <- rev(rownames(paired.test))

# intersect sample IDs between mut cnv and proportions
s.int <- intersect(colnames(mut.type), colnames(ct.mdp))

ct.mdp <- ct.mdp[cto,s.int]
mut.type <- mut.type[,s.int]

#========================================================
# overall waterfall plot of all mutations with assoc. > 0
#========================================================

k <- wr.smg[,"w.p_adj"] < 0.05
wr.smgk <- wr.smg[k,]
message("Number of significant genes with CTP association: ", length(unique(wr.smgk[,"Gene"])))

k <- wr.smg[,"w.p_adj"] < 0.05 & wr.smg[,"t.statistic"] > 0
wr.smgk <- wr.smg[k,]
message("Number of significant genes with positive CTP association: ", length(unique(wr.smgk[,"Gene"])))

# Get 0-1 binary mutation for gene-samples
# subset to genes with significant mutation prop. association
mut.g2p <- unique(wr.smgk[,"Gene"])
mut.type.k <- mut.type[mut.g2p,]
mut.bin <- mut.type.k
mut.bin <- mut.bin != "None"

cts <- unique(wr.smgk[,"CellType"])
o <- order(paired.test[cts,"t.estimate"], decreasing=FALSE)
cts <- cts[o]

# order samples with mutations first
mo <- rep(0, ncol(mut.bin))
for (i in 1:nrow(mut.bin)){
    m <- nrow(mut.bin) - (i-1)
    mo <- mo + mut.bin[i,]*(10^m)
}
so <- colnames(mut.bin)[order(mo, decreasing=TRUE)]

mut.type.k <- mut.type.k[,so]

mut.type.km <- reshape2::melt(as.matrix(mut.type.k))

lab2leg1 <- unique(mut.type.km[,"value"])
lab2leg1 <- setdiff(lab2leg1, "None")
leg1labs <- c("Frame Shift", "Missense", "Splice Site", 
              "In Frame", "Nonsense", "Synonymous")

ct.mdp <- ct.mdp[cts,so, drop=FALSE]
ct.mdp <- scale_clip(ct.mdp)
ct.mdpm <- melt(ct.mdp)
ct.mdpm[,"Var1"] <- as.character(ct.mdpm[,"Var1"])

na_col <- "#A9A9A9A0"

p1 <- ggplot(mut.type.km, aes(x = Var2, y = Var1, fill = value)) + 
    geom_tile() + 
    labs(x = NULL, y = NULL) + 
    scale_fill_discrete(name = "Type", limits = lab2leg1, 
                        na.value = na_col, labels = leg1labs) + 
    theme_bw() + 
    theme(text = element_text(size = 10), 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(),
          panel.grid = element_blank(), 
          panel.border = element_blank(), 
          legend.key.size = unit(0.2, "in"))


p2 <- ggplot(ct.mdpm, aes(x = Var2, y = Var1, fill = value)) + 
    geom_tile() + 
    labs(x = NULL, y = NULL) + 
    scale_fill_gradientn(name = "Cell-type\nproportions\n(scaled)", 
                         colours = rev(rd_bu)) + 
    theme_bw() + 
    theme(text = element_text(size = 10),
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          panel.grid = element_blank(), 
          panel.border = element_blank())

outfn <- paste0(dir_exp, "mut.pdf")
pdf(outfn, width=7, height=4)
grobs <- list(p1, p2)
grid.arrange(grobs[[1]], grobs[[2]], ncol=1)
# ggarrange(p1, p2, ncol=1, heights=c(4, 5), newpage=FALSE)
dev.off()

#========================================================
# aux data
#========================================================

ast_size <- 6

# factor levels for mutation types
mutmap <- c("None" = "None", 
            "Synonymous" = "Synonymous", 
            "In_Frame" = "In Frame", 
            "Missense" = "Missense", 
            "Splice_Site" = "Splice Site", 
            "Frame_Shift" = "Frame Shift", 
            "Nonsense" = "Nonsense")

# mutation type
mut.type.k.any <- mut.type
mut.type.k.any <- as.data.frame(t(mut.type.k.any))
for (i in 1:ncol(mut.type.k.any)){
    mut.type.k.any[,i] <- mutmap[mut.type.k.any[,i]]
    mut.type.k.any[,i] <- factor(mut.type.k.any[,i], levels = mutmap) 
    mut.type.k.any[,i] <- droplevels(mut.type.k.any[,i])
}
datf.type <- cbind(mut.type.k.any, t(ct.mdp[,rownames(mut.type.k.any),drop=FALSE]))

out_fn <- paste0(dir_exp, "sample.mut_type.txt")
write.table(datf.type, out_fn, row.names=TRUE, col.names=NA, quote=FALSE, sep='\t')

# any mutation
mut.k.any <- mut.type
mut.k.any[mut.k.any != "None"] <- "Mut"
mut.k.any[mut.k.any == "None"] <- "WT"
mut.k.any <- as.data.frame(t(mut.k.any))
for (i in 1:ncol(mut.k.any)){
    mut.k.any[,i] <- factor(mut.k.any[,i], 
                                 levels = c("WT", "Mut"))
}

ct.mdpdf <- as.data.frame(t(ct.mdp))
datf.mut <- cbind(mut.k.any, ct.mdpdf[rownames(mut.k.any),])

out_fn <- paste0(dir_exp, "sample.mut_wt_stat.txt")
write.table(datf.mut, out_fn, row.names=TRUE, col.names=NA, quote=FALSE, sep='\t')


#========================================================
# plot function
#========================================================

bp <- function(datf, x, y){
    ylims <- range(datf[,y])
    r <- ylims[2] - ylims[1]
    rp <- r * 0.1
    ly <- ylims[2] + (rp * 0.5)

    p <- ggplot(datf, aes_string(x=x,y=y)) + 
        geom_boxplot(outlier.shape=NA) + 
        geom_dotplot(binaxis='y', stackdir='center', binwidth = .05) + 
        geom_segment(aes(x = 0, y = ly, xend = 1, yend = ly)) + 
        ylim(c(ylims[1], ylims[2] + tp)) + 
        theme_bw() + 
        theme(panel.grid = element_blank(), 
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank())
        # geom_jitter(width = 0.125) + 
        # geom_point(pch = 21, position = position_jitterdodge()) + 
}

get_astk <- function(x){
    if (is.na(x)) return(NULL)
    if (x < 0.0005) return("***")
    if (x < 0.005) return("**")
    if (x < 0.05) return("*")
    return(NULL)
}

#========================================================
# Plot prop in mut vs WT per gene
#========================================================

datf2 <- datf.mut[,c("TP53", "RB1", "Prol")]
datf2m <- melt(datf2, id.vars="Prol")
datf2m[,"value"] <- factor(datf2m[,"value"], levels = c("WT", "Mut"))


ylims <- range(datf2m[,"Prol"])
r <- ylims[2] - ylims[1]
ystep <- r * 0.1
ly <- ylims[2]
ly <- ylims[2] + (ystep * 0.5)
py <- ylims[2] + (ystep * 0.75)

p <- ggplot(datf2m, aes(x = value, y = Prol)) + 
    geom_boxplot(outlier.shape=NA, color = "grey") + 
    geom_dotplot(binaxis='y', stackdir='center', binwidth = .05) +
    geom_segment(aes(x = 1, y = ly, xend = 2, yend = ly)) +
    annotate("text", label = "***", x=1.5, y = ly, size=ast_size, vjust=0.5) + 
    ylim(c(ylims[1], ylims[2] + ystep)) +
    facet_wrap(~variable, nrow=1) + 
    theme_bw() + 
    labs(x = NULL, y = "Hep-20") + 
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(color="black"), 
          axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(), 
          strip.background = element_rect(fill = "white"))
# p <- bp(datf1, x = "TP53", y = "Hep.20")
outfn <- paste0(dir_exp, "gene.ct.pdf")
ggsave(outfn, width = 3.5, height = 2.5)

#========================================================
# TP53 mutation type by Prol expression
#========================================================

ct <- "Prol"
mu <- "TP53"

k <- wr.smg.type[,"CellType"] == ct & wr.smg.type[,"Gene"] == mu
wr.smg.typek <- wr.smg.type[k,]
wr.smg.typek[,"MutType"] <- mutmap[wr.smg.typek[,"MutType"]]
rownames(wr.smg.typek) <- wr.smg.typek[,"MutType"]

p <- ggplot(datf.type, aes_string(x = mu, y = ct)) + 
    geom_boxplot(outlier.shape=NA, color = "grey") + 
    geom_dotplot(binaxis='y', stackdir='center', binwidth = .05) +
    theme_bw() + 
    labs(x = paste0(mu, " mutation type"), y = ct) + 
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 30, hjust=0.9, color="black"),
          axis.text.y = element_blank(), 
          axis.ticks.y = element_blank())

ylims <- range(datf.type[,ct])
r <- ylims[2] - ylims[1]
ystep <- r * 0.1
ly <- ylims[2]
ly <- ylims[2] + (ystep * 0.5)

mut_types <- wr.smg.typek[,"MutType"]
mut_typesl <- levels(datf.type[,mu])

for (mut_type in mut_typesl){
    if (! mut_type %in% rownames(wr.smg.typek)) next
    atsk <- get_astk(wr.smg.typek[mut_type,"w.p"])
    if (is.null(atsk)) next
    ly <- ly + ystep
    x1 <- which(mut_typesl == "None")
    x2 <- which(mut_typesl == mut_type)
    p <- p + geom_segment(x = x1, y = ly, xend = x2, yend = ly)
    p <- p + annotate("text", label = atsk, x=(x1+x2)/2, y = ly, size=ast_size, vjust=0.5) 
}

p <- p + ylim(c(ylims[1], ly + ystep))

outfn <- paste0(dir_exp, mu, ".mut_type.", ct, ".prop.pdf")
ggsave(outfn,  width = 3.5, height = 2.5)



#========================================================
# mutation by genes
#========================================================

pw <- 3; ph <- 3

xtitle <- "Difference between Mut\nand WT cell-type proportions"
xtitle <- "Cell-type proportion difference\nbetween Mut and WT"



# TP53
wr.smg.p53 <- wr.smg[wr.smg[,"Gene"] == "TP53",]
datf.p <- data.frame("CellType" = wr.smg.p53[,"CellType"], 
                    "MeanDiff" = wr.smg.p53[,"t.estimate1"] - wr.smg.p53[,"t.estimate2"], 
                    "log10p" = -log10(wr.smg.p53[,"w.p_adj"]))

datf.p[,"CellType"] <- datf.p[,"CellType"]
datf.p[datf.p[,"CellType"] != "Prol", "CellType"] <- ""

p <- ggplot(datf.p, aes(x = MeanDiff, y = log10p)) + 
    geom_vline(xintercept = 0, col = "red") + 
    geom_point(shape=16) + 
    geom_text(aes(label = CellType), color="red", nudge_x = -1, nudge_y = -.5) + 
    theme_bw() + 
    ggtitle("TP53") + 
    labs(x = xtitle, y = expression(-log[10]~p-value)) + 
    theme(panel.grid.minor = element_blank(), 
          plot.title = element_text(hjust=0.5))
outfn <- paste0(dir_exp, "TP53.genes.diff_p.pdf")
ggsave(outfn, width = pw, height = ph)

# RB1
genek <- "RB1"
wr.smg.gene <- wr.smg[wr.smg[,"Gene"] == genek,]
datf.p <- data.frame("CellType" = wr.smg.gene[,"CellType"], 
                    "MeanDiff" = wr.smg.gene[,"t.estimate1"] - wr.smg.gene[,"t.estimate2"], 
                    "log10p" = -log10(wr.smg.gene[,"w.p_adj"]))

datf.p[,"CellType"] <- datf.p[,"CellType"]
datf.p[datf.p[,"CellType"] != "Prol", "CellType"] <- ""

p <- ggplot(datf.p, aes(x = MeanDiff, y = log10p)) + 
    geom_vline(xintercept = 0, col = "red") + 
    geom_point(shape=16) + 
    geom_text(aes(label = CellType), color="red", nudge_x = -1, nudge_y = -.5) + 
    theme_bw() + 
    ggtitle(genek) + 
    labs(x = xtitle, y = expression(-log[10]~p-value)) + 
    theme(panel.grid.minor = element_blank(), 
          plot.title = element_text(hjust=0.5))
outfn <- paste0(dir_exp, genek, ".genes.diff_p.pdf")
ggsave(outfn, width = pw, height = ph)

# CTNNB1
wr.smg.bc <- wr.smg[wr.smg[,"Gene"] == "CTNNB1",]
datf1 <- data.frame("CellType" = wr.smg.bc[,"CellType"], 
                    "MeanDiff" = wr.smg.bc[,"t.estimate1"] - wr.smg.bc[,"t.estimate2"], 
                    "log10p" = -log10(wr.smg.bc[,"w.p_adj"]))

k <- datf1[,"log10p"] < 10
datf1[k,"CellType"] <- ""

p <- ggplot(datf1, aes(x = MeanDiff, y = log10p, label = CellType)) + 
    geom_point(shape=16) + 
    geom_text_repel(color="red") + 
    theme_bw() + 
    ggtitle("CTNNB1") + 
    labs(x = xtitle, y = expression(-log[10]~p-value)) + 
    theme(panel.grid.minor = element_blank(), 
          plot.title = element_text(hjust=0.5))
outfn <- paste0(dir_exp, "CTNNB1.genes.diff_p.pdf")
ggsave(outfn, width = pw, height = ph)

