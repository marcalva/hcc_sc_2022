
# Plot cell-types in UMAP

setwd("../../")

library(Seurat)
library(ggplot2)
library(ggrepel)
library(scales)
library(gridExtra)
library(reshape2)
source("scripts/ggplot_raster.R")
source("scripts/color_pal.R")
source("scripts/get_enr_mat.R")

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

log10_rev_breaks <- function(x, n=5){

    rng <- range(x, na.rm = TRUE) 
    lxmin <- floor(log10(rng[1])) + log10(2)
    lxmax <- ceiling(log10(rng[2])) - log10(2)

    lby <- floor((lxmax - lxmin)/n) + 1

    breaks <- rev(10^seq(lxmax, lxmin, by=-lby))
    return(breaks)
}

format_power10 <- function(x){
    x <- signif(x, digits = 2)
    sapply(x, function(y){
           pow_num <- floor(log10(y))
           base_num <- y * 10^-pow_num
           ret <- bquote(.(base_num) %*% 10^.(pow_num))
           as.expression(ret) })
}

log10_rev_trans <- function(x){
    trans <- function(x) -log(x, 10)
    inv <- function(x) 10^(-x)
    trans_new("log10_rev", trans, inv, breaks = log10_rev_breaks,
              format = format_power10, domain = c(1e-100, Inf))
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

# Set directories
dir_plt <- "exp/sharma_aiz/plots/";
dir.create(dir_plt, showWarnings=FALSE, recursive=TRUE)

fn <- "data/processed/sharma_aiz/liver.int_rand.rds"
integrated <- readRDS(fn)

udf <- integrated@reductions$umap@cell.embeddings
colnames(udf) <- c("UMAP1", "UMAP2")
udf <- cbind(udf, integrated@meta.data[rownames(udf),])

#=========================================
# plotting theme
#=========================================

theme_txt <- theme(text = element_text(size = 8),
                   plot.title = element_text(hjust = 0.5))

theme_leg <- theme(legend.key.size = unit(6, "points"), 
                   legend.spacing.x = unit(2,"points"), 
                   legend.spacing.y = unit(5,"points"), 
                   legend.title = element_text(vjust = 0.5, hjust = 0), 
                   legend.margin = margin(0,0,0,0), 
                   legend.text = element_text(size = 6, hjust = 0, margin = margin(0,5,0,0)),
                   legend.box.spacing = unit(0.5, "strheight", "0"))

theme_axs <- theme(axis.text=element_blank(),
                   axis.ticks=element_blank())

theme_s <- theme_classic() + 
    theme_txt + theme_leg + theme_axs

#=========================================
# UMAP of cell-type assignments
#=========================================

# Cell type labels
labs <- c("cell_type_main", "cell_type_fine")
names(labs) <- labs

leg <- c('Main', 'Fine')
names(leg) <- labs

# generate color palettes
ctu <- lapply(labs, function(lab){
              sort(unique(udf[,lab])) })
pal_labs <- lapply(labs, function(l){
                   nct <- length(ctu[[l]])
                   hcl_pal(nct, chr = c(80, 80), lum = c(60,80), 
                           offset = 0, rand = TRUE, seedn = 321) })

meds <- get_med_points(udf, c("UMAP1", "UMAP2"), "cell_type_main")

for (l in labs){

    p <- ggplot(udf, aes_string(x="UMAP1",y="UMAP2",color=l)) +
    geom_point(shape=16, size=.01) +
    scale_color_manual(values=pal_labs[[l]], name=leg[l]) +
    theme_s + 
    guides(color = guide_legend(override.aes = list(size=1), ncol=1)) + 
    geom_text_repel(data=meds, aes(x=UMAP1, y=UMAP2), size = 3, 
                    label=rownames(meds), color="black") 

    g <- raster_ggpoints(ggplotGrob(p), w=3, h=3, res=600)
    outfn <- file.path(dir_plt, paste0("UMAP.", l, ".pdf"))
    pdf(outfn, width = 3.5, height = 3)
    grid.draw(g)
    dev.off()
}

#=========================================
# pathway enrichment plots
#=========================================

fn <- "exp/sharma_aiz/clust2ct/res1_to_fine.txt"
ctmap <- read.table(fn, header=TRUE, colClasses = "character")
rownames(ctmap) <- ctmap[,1]

fn <- "exp/sharma_aiz/path_enr/markers.res.1.path_gse.rds"
gse_all <- readRDS(fn)

fn <- "data/ref/colors/red_colrs.csv"
reds <- read.table(fn, header=FALSE)
reds <- reds[,1]

theme_e <- theme_bw() + 
    theme(text = element_text(size = 10), 
          axis.text.y = element_text(size = 6, margin=margin(0,1,0,0)), 
          axis.text.x = element_text(size = 6, margin=margin(0,0,0,0), 
                                     angle = 90, hjust = 1, vjust = 0.5), 
          panel.border =  element_blank(), 
          panel.background = element_blank(), 
          panel.grid.major = element_line(colour = reds[1]), 
          plot.title = element_text(hjust = 0.5), 
          axis.ticks = element_line(colour = reds[1]), 
          legend.key.width = unit(1, "strwidth", "0"))

dir_plt <- "exp/sharma_aiz/path_enr/";

e_col <- "NES"
p_col <- "qvalues"

for (a in names(gse_all)){
    react_l <- gse_all[[a]]
    names(react_l) <- ctmap[names(react_l),"cell_type"]
    react_l <- lapply(react_l, function(x){
                      k <- x[,e_col] >= 0
                      x <- x[k,,drop=FALSE] 
                      return(x) })
    react_l <- get_enr_mat(react_l, idcol = "ID", dcol = "Description", 
                           ctcol = "CellType", pcol = p_col, ecol = e_col, 
                           max_terms = 3)

    e_df <- react_l[["E"]]
    p_df <- react_l[["P"]]

    # wrap strings
    rn <- rownames(e_df)
    rn <- sapply(rn, function(x){
                 if (nchar(x) > 80){
                     x <- strtrim(x, 77)
                     x <- paste(x, "...", sep="")
                 }
                 return(x) })
    # rn <- sapply(rn, function(x) paste(strwrap(x), collapse='\n'))
    rownames(e_df) <- rownames(p_df) <- rn

    # order cell types
    ro <- hclust(dist(e_df))$order
    co <- hclust(dist(t(e_df)))$order
    rl <- rownames(e_df)[ro]
    cl <- colnames(e_df)[co]

    # melt
    p_dfm <- reshape2::melt(as.matrix(p_df))
    e_dfm <- reshape2::melt(as.matrix(e_df))
    for (i in 1:2){
        p_dfm[,i] <- factor(as.character(p_dfm[,i]), levels=list(rl,cl)[[i]])
        e_dfm[,i] <- factor(as.character(e_dfm[,i]), levels=list(rl,cl)[[i]])
    }

    enr_datf <- data.frame("path" = p_dfm[,1], "cluster" = p_dfm[,2], 
                           "enr" = e_dfm[,3], "p" = p_dfm[,3])

    p <- ggplot(enr_datf, aes(x = cluster, y = path, color = enr, size = p)) + 
    geom_point() + 
    theme_e + labs(x=NULL, y=NULL) + ggtitle(a) + 
    scale_size(name = "q-value", trans="log10_rev", range = c(0,3)) + 
    scale_color_gradientn(colors = reds) + 
    guides(size = guide_legend(override.aes = list(shape = 1, color = "black")))

    fn <- file.path(dir_plt, paste0("markers.fine.", a, ".bubble.pdf"))
    ggsave(fn, width = 7, height = 5)
}

