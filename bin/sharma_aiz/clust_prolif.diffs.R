
# cluster proliferating cells, diff. prop.

setwd("../../")

library(Seurat)
library(ggplot2)
library(grid)
library(gridExtra)
library(ggrepel)
source("scripts/color_pal.R")
source("scripts/diff_test_ml.R")

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

fn <- "data/processed/sharma_aiz/liver.int_rand.rds"
integrated <- readRDS(fn)

fn <- "data/processed/sharma_aiz/liver.prolif.int.rds"
prolif <- readRDS(fn)

# Set directories
dir_out <- "data/processed/sharma_aiz/"
dir.create(dir_out, showWarnings=FALSE, recursive=TRUE)

dir_exp <- "exp/sharma_aiz/clust_prolif/";
dir.create(dir_exp, showWarnings=FALSE, recursive=TRUE)

#=========================================
# plot themes and colors
#=========================================

theme_trnsp <- theme(plot.background = element_rect(fill="transparent", color=NA),
                     panel.background = element_rect(fill="transparent", color=NA),
                     legend.background = element_rect(fill="transparent", color=NA))
theme_leg <- theme(legend.key.height = unit(2, "strheight", "0"),
                   legend.key.width = unit(1, "strwidth", "0"))
theme_txt <- theme(text = element_text(size = 8),
                   axis.text=element_blank(), 
                   axis.ticks=element_blank(),
                   plot.title = element_text(hjust = 0.5))
theme_s <- theme_classic() + 
    theme_leg + 
    theme_txt

theme_p <- theme_classic() +
    theme_txt + 
    theme(legend.key.height = unit(1, "strheight", "0"),
                               legend.key.width = unit(1, "strwidth", "0"))


pal10 <- read.table("data/ref/colors/hcl_c70_vl.pal.10.csv", header=FALSE)
pal10 <- pal10[,1]
set.seed(3, kind = 'Mersenne-Twister')
pal10 <- sample(pal10)

reds <- read.table("data/ref/colors/red_colrs.csv", header=FALSE)
reds <- reds[,1]

umap <- prolif@reductions$umap@cell.embeddings
datf <- prolif@meta.data
datf <- cbind(datf, umap[rownames(datf),])
ncol <- length(unique(datf[,"cell_type_main"]))
cpal <- hcl_pal(ncol, chr = c(80, 80), lum = c(50,90), offset = 10, 
                rand = TRUE, seedn = 3)

#=========================================
# function for FET of freq table
#=========================================

fet_mat <- function(x){
    fret_l <- lapply(1:nrow(x), function(i){
                     ki <- i
                     kn <- setdiff(1:nrow(x), i)
                     mat <- cbind(colSums(x[ki,,drop=FALSE]),
                                  colSums(x[kn,,drop=FALSE]))
                     fr <- fisher.test(mat)
                     data.frame("OR" = fr$estimate, "p" = fr$p.value) })
    fret <- do.call(rbind, fret_l)
    rownames(fret) <- rownames(x)
    return(fret)
}

#=========================================
# data frames for prop tests
#=========================================

md_np <- integrated@meta.data
md_np <- md_np[md_np[,"cell_type_main"] != "Prol",]
md_np[,"cell_type_main"] <- as.character(md_np[,"cell_type_main"])

md_p <- prolif@meta.data

k <- md_np[,"PatientStat"] != "Aizarani" & md_np[,"PatientStat"] != "Sharma_0"
md_np <- md_np[k,]
k <- md_p[,"PatientStat"] != "Aizarani" & md_p[,"PatientStat"] != "Sharma_0"
md_p <- md_p[k,]

#=========================================
# Proportion of main cell-types in prolif 
# compared to non-prolif
#========================================

prol_comp <- table(md_p[,"cell_type_main"])
nonprol_comp <- table(md_np[,"cell_type_main"])

main_freq <- data.frame("Prol" = as.numeric(prol_comp), 
                        "NonProl" = as.numeric(nonprol_comp[names(prol_comp)]))
rownames(main_freq) <- names(prol_comp)
main_prop <- sweep(main_freq, 2, colSums(main_freq), "/")

main_freq_fet <- fet_mat(main_freq)

#=========================================
# Proportion of main cell-types in prolif 
# comparing tumor and non-tumor
#=========================================

prol_tum_freq <- table(md_p[,"cell_type_main"], md_p[,"TumorStat"])

prol_tum_prop <- sweep(prol_tum_freq, 2, colSums(prol_tum_freq), '/')

prol_tum_fet <- fet_mat(prol_tum_freq[,c("Tumor", "NonTumor")])

#=========================================
# diff test in proportions between prol and non-prol
#=========================================

md_prol <- prolif@meta.data
md_prol <- md_prol[md_prol$TumorStat == "Tumor",]

md <- integrated@meta.data
md <- md[md$TumorStat == "Tumor",]
md$cell_type_p <- md$cell_type_main
md[rownames(md_prol),"cell_type_p"] <- md_prol$cell_type_main
md[,"cell_type_p"] <- factor(md[,"cell_type_p"])
md[,"Prol"] <- "Non-prol"
md[md[,"cell_type_main"] == "Prol","Prol"] <- "Prol"

pats <- sort(unique(md[,"PatientStat"]))
k_pats <- setdiff(pats, c("Aizarani"))
md <- md[md[,"PatientStat"] %in% k_pats,]

md[,"PatientStatP"] <- paste0(md[,"PatientStat"], "_", md[,"Prol"])

p_diffs <- diff_test_ml(x = md[,"PatientStatP"], g = md[,"Prol"], z = md[,"cell_type_p"], 
                        paired=TRUE, p = md[,"PatientStat"])

#=========================================
# plot bar plot of all proportions
#=========================================

theme_bar <- theme_classic() + 
theme(text = element_text(size = 6), 
      axis.text = element_text(color = "black"), 
      axis.line = element_blank(), 
      axis.ticks.x = element_blank(), 
      axis.ticks.length.x = unit(.1, "strheight", "0"), 
      axis.ticks.y = element_line(size = .5), 
      plot.title = element_text(hjust = 0.5), 
      legend.margin = margin(0,0,0,0), 
      legend.box.spacing = unit(1, "strheight", "0"), 
      legend.key.size = unit(1, "strheight", "0"))

cts <- sort(unique(datf[,"cell_type_main"]))
cpal <- hcl_pal(length(cts), chr = c(80, 80), lum = c(50,90), offset = 10, 
                rand = TRUE, seedn = 3)
names(cpal) <- cts

p <- ggplot(md, aes(x = Prol, fill = cell_type_p)) + 
geom_bar(position = "fill") + 
scale_fill_manual(values = cpal, name = "Cell Type") + 
scale_x_discrete(expand = expansion(mult = 0.0, add = 0)) + 
scale_y_continuous(expand = expansion(mult = 0, add = 0)) + 
labs(x = NULL, y = NULL) + 
theme_bar

out_fn <- paste0(dir_exp, "prol.ct_main_prop.bar.pdf")
ggsave(out_fn, width = 2, height = 2)

