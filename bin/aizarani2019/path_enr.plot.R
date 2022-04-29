
# plot a heatmap of GO enrichment terms for hepatocytes

setwd("../../")

library(Seurat)
library(ggplot2)
library(scales)
library(reshape2)

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

#' given annotation output from pathway enrichment, 
#' return log odds and p-value matrices
get_enr_mat <- function(l, idcol = "ID", dcol = "Description", ecol = "log_odds", 
                        ctcol = "CellType", pcol = "qvalue", maxp = 0.05, max_terms = 5){
    # filter enrichment data frames
    k <- sapply(l, nrow) > 0
    l <- l[k]
    l <- lapply(l, function(datf){
                k <- !is.na(datf[,dcol]) & !is.na(datf[,idcol]) & 
                   !is.na(datf[,ecol])
                datf[k,,drop=FALSE] })

    # ID to description map
    id2desc <- lapply(l, function(datf){
                      datf <- datf[,c(dcol,idcol),drop=FALSE]
                      return(datf) })
    id2desc <- do.call(rbind, id2desc)
    k <- !duplicated(id2desc[,idcol])
    id2desc <- id2desc[k,,drop=FALSE]
    rownames(id2desc) <- id2desc[,idcol]

    # get significant enrichment terms
    sig_terms_l <- lapply(l, function(datf){
                          trms <- rownames(datf)[datf[,pcol] < maxp]
                          len <- min(length(trms), max_terms)
                          trms <- trms[1:len]
                          return(trms)
                      })

    sig_terms <- unique(do.call(c, sig_terms_l))
    sig_terms <- sig_terms[!is.na(sig_terms)]

    # bind results
    cts <- names(l)
    names(cts) <- cts
    ont_l <- lapply(cts, function(n){
                    datf <- l[[n]]
                    datf <- datf[sig_terms,]
                    rownames(datf) <- sig_terms
                    datf[,idcol] <- sig_terms
                    datf[,ctcol] <- n
                    to1 <- is.na(datf[,ecol])
                    datf[to1,ecol] <- 0
                    datf[to1,pcol] <- 1
                    return(datf)})
    e_df <- do.call(cbind, lapply(ont_l, function(x) x[,ecol, drop=FALSE]))
    colnames(e_df) <- cts

    p_df <- do.call(cbind, lapply(ont_l, function(x) x[,pcol, drop=FALSE]))
    colnames(p_df) <- cts

    rn <- id2desc[rownames(e_df),dcol]

    rownames(e_df) <- rownames(p_df) <- rn

    return(list("E" = e_df, "P" = p_df, "id2desc" =  id2desc))
}

#=============================================
#=============================================

fn <- "data/ref/colors/rd_bu_div.csv"
rd_bu <- read.table(fn, header=FALSE)
rd_bu <- rev(rd_bu[,1])

fn <- "data/ref/colors/red_colrs.csv"
reds <- read.table(fn, header=FALSE)
reds <- reds[,1]

fn <- "exp/aizarani2019/sct_clust/path_enr/markers.res.0.5.path_enr.rds"
enr_all <- readRDS(fn)

fn <- "exp/aizarani2019/sct_clust/path_enr/markers.res.0.5.path_gse.rds"
gse_all <- readRDS(fn)

#=============================================
#=============================================

dir_out <- "exp/aizarani2019/sct_clust/path_enr/"

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

# plot reactome enrichment
# react_l <- enr_all[["Reactome"]]
# react_l <- lapply(react_l, function(x){
#                   k <- x[,"FoldEnrich"] >= 1
#                   x <- x[k,,drop=FALSE] 
#                   x[,"log_odds"] <- log2(x[,"log_odds"]) 
#                   return(x) })

for (a in names(gse_all)){
    react_l <- gse_all[[a]]
    react_l <- lapply(react_l, function(x){
                      k <- x[,"NES"] >= 0
                      x <- x[k,,drop=FALSE] 
                      return(x) })
    react_l <- get_enr_mat(react_l, idcol = "ID", dcol = "Description", 
                           ctcol = "CellType", pcol = "qvalues", ecol = "NES")

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
    p_dfm[,1] <- factor(as.character(p_dfm[,1]), levels=rl)
    p_dfm[,2] <- factor(as.character(p_dfm[,2]), levels=cl)
    e_dfm <- reshape2::melt(as.matrix(e_df))
    e_dfm[,1] <- factor(as.character(e_dfm[,1]), levels=rl)
    e_dfm[,2] <- factor(as.character(e_dfm[,2]), levels=cl)

    enr_datf <- data.frame("path" = p_dfm[,1], "cluster" = p_dfm[,2], 
                           "NES" = e_dfm[,3], "p" = p_dfm[,3])


    p <- ggplot(enr_datf, aes(x = cluster, y = path, color = NES, size = p)) + 
    geom_point() + 
    theme_s + labs(x=NULL, y=NULL) + ggtitle(a) + 
    scale_size(name = "q-value", trans="log10_rev", range = c(0,3)) + 
    scale_color_gradientn(colors = reds) + 
    guides(size = guide_legend(override.aes = list(shape = 1, color = "black")))

    fn <- file.path(dir_out, paste0("markers.res.1.", a, ".bubble.pdf"))
    ggsave(fn, width = 7, height = 7)
}

#=============================================
#=============================================

