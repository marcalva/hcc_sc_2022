
# plot a heatmap of GO enrichment terms for hepatocytes

setwd("../../")

library(ggplot2)
source("~/scripts/go_enr/R/enrich.R")
source("~/scripts/go_enr/R/read.R")
source("~/scripts/go_enr/R/plot.R")

#=============================================
#=============================================

infn <- "exp/d_sct_cca_all/merge/markers/seur.CellType.markers.txt"
m <- read.table(infn, header = TRUE, check.names = FALSE, 
                stringsAsFactors = FALSE)
cts <- unique(m[,"cluster"])


infn <- "exp/d_sct_cca_all/merge/markers/exprd_gn.CellType.txt"
e <- read.table(infn, header = TRUE, check.names = FALSE, 
                stringsAsFactors = FALSE)

dirin <- "exp/d_sct_cca_all/merge/enr/CellType/Reactome/"
ont_l <- lapply(cts, function(ct){
               infn <- paste0(dirin, "reactome.", ct, ".txt")
               datf <- read.table(infn, header=TRUE, row.names=1, sep="\t", quote="", fill = TRUE,
                                  check.names=FALSE, stringsAsFactors=FALSE)
                })
names(ont_l) <- cts

sig_terms_l <- lapply(ont_l, function(ont){
                      trms <- rownames(ont)[ont[,"p_adj"] < 0.05]
                      len <- min(length(trms), 5)
                      trms <- trms[1:len]
                      return(trms)
                })

sig_terms <- do.call(c, sig_terms_l)
sig_terms <- unique(sig_terms)
sig_terms <- sig_terms[!is.na(sig_terms)]

ont_l <- lapply(ont_l, function(ont){
                ont <- ont[sig_terms,]
                return(ont)})

log_odds_l <- lapply(ont_l, function(ont){
                     datf <- data.frame("log_odds" = log(ont[,"OR"]))
                     rownames(datf) <- ont[,"Name"]
                     return(datf) })
log_odds <- do.call(cbind, log_odds_l)
colnames(log_odds) <- cts

p_l <- lapply(ont_l, function(ont){
                     datf <- data.frame("p" = ont[,"p_adj"])
                     rownames(datf) <- ont[,"Name"]
                     return(datf) })
pv <- do.call(cbind, p_l)
colnames(pv) <- cts

wrap_nl <- function(x, width=60){
    paste0(strwrap(x, width=width, simplify=TRUE), collapse="\n")
}

cut_str <- function(x, width = 60){
    if (nchar(x) > width){
        x <- substr(x, 1, width-3)
        x <- paste0(x, "...")
    }
    return(x)
}
rownames(log_odds) <- sapply(rownames(log_odds), cut_str)
rownames(pv) <- sapply(rownames(pv), cut_str)

row_order <- hclust(dist(log_odds))$order
col_order <- hclust(dist(t(log_odds)))$order

log_odds <- log_odds[row_order, col_order]
pv <- pv[row_order, col_order]

dirout <- "exp/d_sct_cca_all/merge/enr/CellType/Reactome/"
dir.create(dirout, showWarnings = FALSE, recursive = TRUE)

p <- heatmap_enr(t(log_odds), t(pv))
outfn <- paste0(dirout, "reactome.heatmap.pdf")
ggsave(outfn, width =12, height = 12)


