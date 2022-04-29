
# plot a heatmap of GO enrichment terms for hepatocytes

setwd("../../")

library(ggplot2)
library(scales)
library(reshape2)
source("scripts/get_enr_mat.R")

#=============================================
# plotting functions
#=============================================

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


#=============================================
#=============================================

fn <- "data/ref/colors/rd_bu_div.csv"
rd_bu <- read.table(fn, header=FALSE)
rd_bu <- rev(rd_bu[,1])

fn <- "data/ref/colors/red_colrs.csv"
reds <- read.table(fn, header=FALSE)
reds <- reds[,1]

fn <- "exp/sharma_aiz/path_enr/markers.res.1.path_enr.rds"
enr_all <- readRDS(fn)

fn <- "exp/sharma_aiz/path_enr/markers.res.1.path_gse.rds"
gse_all <- readRDS(fn)

#=============================================
#=============================================

dir_out <- "exp/sharma_aiz/path_enr/"

#=============================================
#=============================================

theme_s <- theme_bw() + 
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

#=============================================
# plot GSEA results using NES scores
#=============================================

e_col <- "NES"
p_col <- "qvalues"

for (a in names(gse_all)){
    react_l <- gse_all[[a]]
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
    theme_s + labs(x=NULL, y=NULL) + ggtitle(a) + 
    scale_size(name = "q-value", trans="log10_rev", range = c(0,3)) + 
    scale_color_gradientn(colors = reds) + 
    guides(size = guide_legend(override.aes = list(shape = 1, color = "black")))

    fn <- file.path(dir_out, paste0("markers.res.1.", a, ".bubble.pdf"))
    ggsave(fn, width = 7, height = 5)
}

