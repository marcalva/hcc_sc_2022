
# Figure S8
# Proportion-mutation associations

setwd("../../")

library(ggplot2)
library(gtable)
library(grid)
library(reshape2)
source("scripts/gtable_stack.R")

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

# somatic m
fn <- "data/processed/tcga_hcc/mut/tcga.lihc.mut.gene.any.txt"
mut.any <- read.table(fn, header = TRUE, stringsAsFactors = FALSE, 
                      sep = "\t", row.names=1, check.names=FALSE)
fn <- "data/processed/tcga_hcc/mut/tcga.lihc.mut.gene.type.txt"
mut.type <- read.table(fn, header = TRUE, stringsAsFactors = FALSE, 
                       sep = "\t", row.names=1, check.names=FALSE)

mut.any[mut.any >= 1] <- 1

# somatic mutation results
fn <- paste0("exp/tcga_hcc/mol/", ct_id, ".SMG_mutsig.2CV.ctp.txt")
wr.smg <- read.table(fn, header=TRUE, sep="\t")

# somatic mutation results
fn <- paste0("exp/tcga_hcc/mol/", ct_id, ".SMG_mutsig.2CV.mut_type.ctp.txt")
wr.smg.type <- read.table(fn, header=TRUE, sep="\t")

# colors
fn <- "data/ref/colors/rd_bu_div.csv"
rd_bu <- read.table(fn, header=FALSE)
rd_bu <- rd_bu[,1]

# Set directories
dir_plt <- "exp/manuscript/"
dir.create(dir_plt, showWarnings=FALSE, recursive=TRUE)

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

# get significant associations
k <- wr.smg[,"w.p_adj"] < 0.05
wr.smgk <- wr.smg[k,]
message("Number of significant genes with CTP association: ", length(unique(wr.smgk[,"Gene"])))

# get significant associations (up in mut)
# k <- wr.smg[,"w.p_adj"] < 0.05 & wr.smg[,"t.statistic"] > 0 & wr.smg[,"CellType"] == "Prol"
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

# melt for plotting
mut.type.km <- reshape2::melt(as.matrix(mut.type.k))
mut.type.km[,"Var2"] <- factor(mut.type.km[,"Var2"], levels = so)

# mut legend labels
lab2leg1 <- unique(mut.type.km[,"value"])
lab2leg1 <- setdiff(lab2leg1, "None")
leg1labs <- c("Frame Shift", "Missense", "Splice Site", 
              "In Frame", "Nonsense", "Synonymous")

# proportions for plotting
ct.mdp <- ct.mdp[cts,so, drop=FALSE]
ct.mdp <- scale_clip(ct.mdp)
ct.mdpm <- melt(ct.mdp)
ct.mdpm[,"Var1"] <- as.character(ct.mdpm[,"Var1"])
ct.mdpm[,"Var2"] <- factor(ct.mdpm[,"Var2"], levels = so)

na_col <- "#A9A9A9A0"

#========================================================
# create waterfall plot
#========================================================

th_txt <- theme(text = element_text(size = 10), 
      plot.title = element_text(hjust = 0.5, size = 12), 
      axis.text.x = element_blank(), 
      axis.text.y = element_text(color="black", size = 10), 
      axis.ticks.x = element_blank())

th_pan <- theme(panel.grid = element_blank(), 
      panel.border = element_blank())


th_mut <- theme_bw() + 
th_txt +
th_pan + 
theme(legend.text = element_text(size = 8), 
      legend.title = element_text(size = 10), 
      legend.key.size = unit(12, "points"))

th_prop <- theme_bw() + 
th_txt + 
th_pan + 
theme(legend.text = element_text(size = 8), 
      legend.title = element_text(size = 10), 
      legend.key.size = unit(12, "points"))

ptitle <- "Proportion-mutation associations"
ltitle <- "Mutation\ntype"
ytitle <- "Gene"

p1 <- ggplot(mut.type.km, aes(x = Var2, y = Var1, fill = value)) + 
    geom_tile() + 
    labs(x = NULL, y = ytitle, title = ptitle) + 
    scale_y_discrete(expand = expansion(0,0)) + 
    scale_fill_discrete(name = ltitle, limits = lab2leg1, 
                        na.value = "#A9A9A900", labels = leg1labs) + 
th_mut + theme(panel.background = element_rect(fill = "#A9A9A9A0"))

gr1 <- ggplotGrob(p1)

# remove bottom rows
pix <- which(gr1$layout$name == "panel")
pr <- gr1$layout[pix, "b"]
gr1 <- gr1[1:pr,]
gr1 <- gtable_add_rows(gr1, heights = unit(12, "points"), nrow(gr1)) # pad space between top and bottom row

ytitle <- "Main cell-type"
xtitle <- "TCGA tumor sample"

p2 <- ggplot(ct.mdpm, aes(x = Var2, y = Var1, fill = value)) + 
    geom_tile() + 
    labs(x = xtitle, y = ytitle) + 
    scale_y_discrete(expand = expansion(0,0)) + 
    scale_fill_gradientn(name = "Cell-type\nproportions\n(scaled)", 
                         colours = rev(rd_bu), limits=c(-3,3)) + 
th_prop

gr2 <- ggplotGrob(p2)

# remove top rows
pix <- which(gr2$layout$name == "panel")
pr <- gr2$layout[pix, "t"]
gr2 <- gr2[pr:nrow(gr2),]

gr_l <- list(gr1, gr2)
gt <- stack_gtable_v(gr_l, heights = unit(c(0.65, 0.35), "npc"))

ws <- 7
hs <- 5

pdf(paste0(dir_plt, "FigureS9.pdf"), width = sum(ws), height = sum(hs))
grid.draw(gt)
dev.off()

