
# Perform survival analysis for each gene genome-wide in 
# LCI data 

setwd("../../")

library(survival)
library(ggplot2)
library(ggfortify)
source("scripts/surv_df.R")

#=====================================================
# functions
#=====================================================

adjust_p <- function(datf, pcols = NULL, adj_method = "fdr"){
    datf <- as.data.frame(datf)
    if (is.null(pcols[1]))
        pcols <- grep("p$", colnames(datf), ignore.case=TRUE)
    for (i in pcols){
        oldcolname <- colnames(datf)[i]
        newcolname <- paste0(oldcolname, "_adj")
        datf[,newcolname] <- p.adjust(datf[,oldcolname], method=adj_method)
    }
    return(datf)
}

cox_list_summary <- function(cox.l){
    cox.df.l <- lapply(names(cox.l), function(n){
                       x <- cox.l[[n]]
                       ret <- sapply(x, function(cox){
                                     cfs <- summary(cox)$coeff
                                     y <- cfs[nrow(cfs),] 
                                     cis <- summary(cox)$conf.int
                                     ci <- cis[nrow(cis), c("lower .95", "upper .95")]
                                     y <- c(y, ci)
                                     names(y) <- c("coef", "exp_coef", "se_coef", "z", "p", "lower95", "upper95")
                                     y <- y[c("coef", "exp_coef", "lower95", "upper95", "se_coef", "z", "p")]
                                     return(y)}) 
                       rownames(ret) <- paste0(n, ".", rownames(ret))
                       return(ret) })
    cox.df <- do.call(rbind, cox.df.l)
    cox.df <- t(cox.df)
    return(cox.df)
}

# include number events observed in groups
cox_list_summary_grp <- function(cox.l){
    cox.df.l <- lapply(names(cox.l), function(n){
                       x <- cox.l[[n]]
                       ret <- sapply(x, function(cox){
                                     cfs <- summary(cox)$coeff
                                     y <- cfs[nrow(cfs),] 
                                     cis <- summary(cox)$conf.int
                                     ci <- cis[nrow(cis), c("lower .95", "upper .95")]
                                     n_obs <- c(cox$n_low, cox$n_high, cox$n_obs_low, cox$n_obs_high)
                                     n_obs <- c("num_events_low" = cox$n_obs_low, "pct_events_low" = cox$n_obs_low / cox$n_low, 
                                                "num_events_high" = cox$n_obs_high, "pct_events_high" = cox$n_obs_high / cox$n_high)
                                     y <- c(n_obs, y, ci)
                                     names(y) <- c("num_events_low", "pct_events_low", "num_events_high", "pct_events_high", 
                                                   "coef", "exp_coef", "se_coef", "z", "p", "lower95", "upper95")
                                     y <- y[c("num_events_low", "pct_events_low", "num_events_high", "pct_events_high", 
                                              "coef", "exp_coef", "lower95", "upper95", "se_coef", "z", "p")]
                                     return(y)}) 
                       rownames(ret) <- paste0(n, ".", rownames(ret))
                       return(ret) })
    cox.df <- do.call(rbind, cox.df.l)
    cox.df <- t(cox.df)
    return(cox.df)
}

#=====================================================
#=====================================================

dir_exp <- "exp/GSE14520/survival.gene/"
dir.create(dir_exp, showWarnings = FALSE, recursive=TRUE)

dir_data <- "data/processed/GSE14520/"
dir.create(dir_data, recursive=TRUE, showWarning=FALSE)

fn <- paste0(dir_data, "expr.RMA_log2.gmean.rds")
ex <- readRDS(fn)

fn <- paste0(dir_data, "sdata.txt")
sdat <- read.table(fn, row.names=1, header=TRUE, sep="\t", stringsAsFactors=FALSE)

fn <- paste0(dir_data, "geo_pheno.txt")
geo_pheno <- read.table(fn, row.names=1, header=TRUE, sep="\t", stringsAsFactors=FALSE)

#=====================================================
# Format covariates
#=====================================================

sdat.tum <- sdat[sdat[,"Tissue.Type"] == "Tumor",]
sdat.nt <- sdat[sdat[,"Tissue.Type"] == "Non-tumor",]

table(duplicated(sdat.tum[,"ID"]))
table(duplicated(sdat.nt[,"ID"]))

rownames(sdat.tum) <- sdat.tum[,"ID"]
rownames(sdat.nt) <- sdat.nt[,"ID"]

lcs.tum <- rownames(sdat.tum)
lcs.nt <- rownames(sdat.nt)

ids.pair <- intersect(rownames(sdat.tum), rownames(sdat.nt))
sdat.tum.p <- sdat.tum[ids.pair,]
sdat.nt.p <- sdat.nt[ids.pair,]

lcs.tum.p <- sdat.tum.p[,"LCS.ID"]
lcs.nt.p <- sdat.nt.p[,"LCS.ID"]

rownames(sdat.tum) <- sdat.tum[,"LCS.ID"]
rownames(sdat.nt) <- sdat.nt[,"LCS.ID"]

cov.cols <- c("Age", "Gender", "stage_low_high")
clin.covs.stg <- sdat.tum[,cov.cols]

# change stage to factor
# ensures low=0 high=1
clin.covs.stg[,"stage_low_high"] <- factor(clin.covs.stg[,"stage_low_high"], 
                                           levels = c("Low", "High"))

# no stage
colrm <- which(colnames(clin.covs.stg) == "stage_low_high")
clin.covs <- clin.covs.stg[,-colrm]

# hazards test
# 0=censored 1=event
events <- c("Survival", "Recurr")

ex.tum <- ex[,rownames(clin.covs.stg),drop=FALSE]

#================================================
# all gene expression model
#================================================

cox.events.g.l <- lapply(events, function(event){
                           event.col <- paste0(event, ".status")
                           time.col <- paste0(event, ".months")

                           all_genes <- rownames(ex.tum)
                           genes.l <- lapply(all_genes, function(gene){
                                                  mat <- cbind(clin.covs.stg, ex.tum[gene,])
                                                  colnames(mat) <- c(colnames(clin.covs.stg), "gene")
                                                  terms <- ~ Age + Gender + gene
                                                  mf <- model.matrix(terms, mat)
                                                  ids.i <- rownames(mf)

                                                  # remove missing event data
                                                  edf <- sdat.tum[ids.i, c(time.col, event.col)]
                                                  k <- apply(edf, 1, function(rw) !any(is.na(rw)))
                                                  k <- rownames(edf)[k]
                                                  mf <- mf[k,-1,drop=FALSE] # remove intercept
                                                  edf <- edf[k,,drop=FALSE]

                                                  te <- Surv(time = edf[,time.col], 
                                                             event = edf[,event.col], 
                                                             type = "right")
                                                  cx <- coxph(te ~ mf)
                                                  return(cx) })
                           names(genes.l) <- all_genes
                           return(genes.l) })

names(cox.events.g.l) <- events
cox.events.g.coeff <- cox_list_summary(cox.events.g.l)

cox.events.g.coeff <- adjust_p(cox.events.g.coeff)

outfn <- paste0(dir_exp, "GSE14520.cox_ph.all_genes.txt")
write.table(cox.events.g.coeff, outfn, 
            row.names = TRUE, col.names = NA, quote=FALSE, sep="\t")

#================================================
# Stage + all gene model
#================================================

cox.events.g_stg.l <- lapply(events, function(event){
                           event.col <- paste0(event, ".status")
                           time.col <- paste0(event, ".months")

                           all_genes <- rownames(ex.tum)
                           genes.l <- lapply(all_genes, function(gene){
                                                  mat <- cbind(clin.covs.stg, ex.tum[gene,])
                                                  colnames(mat) <- c(colnames(clin.covs.stg), "gene")
                                                  terms <- ~ Age + Gender + stage_low_high + gene
                                                  mf <- model.matrix(terms, mat)
                                                  ids.i <- rownames(mf)

                                                  # remove missing event data
                                                  edf <- sdat.tum[ids.i, c(time.col, event.col)]
                                                  k <- apply(edf, 1, function(rw) !any(is.na(rw)))
                                                  k <- rownames(edf)[k]
                                                  mf <- mf[k,-1,drop=FALSE] # remove intercept
                                                  edf <- edf[k,,drop=FALSE]

                                                  te <- Surv(time = edf[,time.col], 
                                                             event = edf[,event.col], 
                                                             type = "right")
                                                  cx <- coxph(te ~ mf)
                                                  return(cx) })
                           names(genes.l) <- all_genes
                           return(genes.l) })
names(cox.events.g_stg.l) <- events
cox.events.g_stg.coeff <- cox_list_summary(cox.events.g_stg.l)

cox.events.g_stg.coeff <- adjust_p(cox.events.g_stg.coeff)

outfn <- paste0(dir_exp, "GSE14520.cox_ph.stage_all_genes.txt")
write.table(cox.events.g_stg.coeff, outfn, 
            row.names = TRUE, col.names = NA, quote=FALSE, sep="\t")

