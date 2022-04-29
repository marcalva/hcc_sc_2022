
# Perform survival analysis of TCGA data with 
# cell type proportions.

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

dir_exp <- "exp/tcga_hcc/survival.gene/"
dir.create(dir_exp, showWarnings = FALSE, recursive=TRUE)

# read in data
out_dir <- "data/processed/tcga_hcc/expr/"

# annotation data
fn <- paste0(out_dir, "tcga.gencodev26.rds")
gencode <- readRDS(fn)

# expression data
fn <- "data/processed/tcga_hcc/expr/tcga.lihc.TMM.log.rin.rds"
tmm <- readRDS(fn)

# sample data
fn <- "data/processed/tcga_hcc/sample/samples.hcc.410.merged.rds"
samples <- readRDS(fn)

# case data
fn <- "data/processed/tcga_hcc/sample/cases.hcc.361.merged.rds"
cases <- readRDS(fn)

# sample map
fn <- "data/processed/tcga_hcc/sample/sam_case_map.rds"
smaps <- readRDS(fn)
cid2sid_tum <- smaps[["cid2sid_tum"]]
cid2sid_nt <- smaps[["cid2sid_nt"]]
cid2sid_tum_pair <- smaps[["cid2sid_tum_pair"]]
cid2sid_nt_pair <- smaps[["cid2sid_nt_pair"]]
sid_all <- c(cid2sid_tum, cid2sid_nt)
sid_both <- c(cid2sid_tum_pair, cid2sid_nt_pair)

tmm <- tmm[,sid_all]

# remove redacted
cases[ is.na(cases[,"cdr.Redaction"]), "cdr.Redaction"] <- "OK"
red.k <- which(cases[,"cdr.Redaction"] == "OK")
cases <- cases[red.k,]

#================================================
# format covariates
#================================================

cases[,"age_at_diagnosis"] <- as.numeric(cases[,"age_at_diagnosis"])

cov.cols <- c("age_at_diagnosis", "gender", "race", "stage_low_high")
clin.covs.stg <- cases[,cov.cols]

# remove amerindian (no events, so blows up hazard ratio to infinite)
k <- which(clin.covs.stg[,"race"] == "AMERICAN INDIAN OR ALASKA NATIVE" | 
           clin.covs.stg[,"race"] == "BLACK OR AFRICAN AMERICAN")

# remove cases with no race entry recorded
clin.covs.stg[k,"race"] <- NA
k.na <- apply(clin.covs.stg, 1, function(x) any(is.na(x)))
clin.covs.stg <- clin.covs.stg[!k.na,]

# change stage to factor
# ensures low=0 high=1
clin.covs.stg[,"stage_low_high"] <- factor(clin.covs.stg[,"stage_low_high"], 
                                           levels = c("Low", "High"))

# no stage
colrm <- which(colnames(clin.covs.stg) == "stage_low_high")
clin.covs <- clin.covs.stg[,-colrm]

# sample IDs
samples.pt <- sid_all[rownames(clin.covs.stg)]

events <- c("cdr.OS", "cdr.DSS", "cdr.DFI", "cdr.PFI")

#================================================
# all gene expression model
#================================================

all_genes <- rownames(tmm) 

g.events.l <- lapply(events, function(event){
                             event.col <- event
                             time.col <- paste0(event, ".time")
                             genes.l <- lapply(all_genes, function(gene){
                                               mat <- cbind(clin.covs.stg, tmm[gene, samples.pt])
                                               colnames(mat)[ncol(mat)] <- "gene"
                                               terms <- ~ age_at_diagnosis + gender + race + gene
                                               mf <- model.matrix(terms, mat)
                                               cases.i <- rownames(mf)

                                               edf <- cases[cases.i,c(time.col, event.col)]
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
                             return(genes.l) 
                                           })
names(g.events.l) <- events
g.events.coeff <- cox_list_summary(g.events.l)

g.events.coeff <- adjust_p(g.events.coeff)

g.events.coeff <- cbind(gencode[rownames(g.events.coeff), ], g.events.coeff)

outfn <- paste0(dir_exp, "cox_ph.all_genes.txt")
write.table(g.events.coeff, outfn, 
            row.names = TRUE, col.names = NA, quote=FALSE, sep="\t")

#================================================
# Stage + all gene expression model
#================================================

all_genes <- rownames(tmm) 

g.events.stg.l <- lapply(events, function(event){
                             event.col <- event
                             time.col <- paste0(event, ".time")
                             genes.l <- lapply(all_genes, function(gene){
                                               mat <- cbind(clin.covs.stg, tmm[gene, samples.pt])
                                               colnames(mat)[ncol(mat)] <- "gene"
                                               terms <- ~ age_at_diagnosis + gender + race + stage_low_high + gene
                                               mf <- model.matrix(terms, mat)
                                               cases.i <- rownames(mf)

                                               edf <- cases[cases.i,c(time.col, event.col)]
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
                             return(genes.l) 
                                           })
names(g.events.stg.l) <- events
g.events.stg.coeff <- cox_list_summary(g.events.stg.l)

g.events.stg.coeff <- adjust_p(g.events.stg.coeff)

g.events.stg.coeff <- cbind(gencode[rownames(g.events.stg.coeff), ], 
                                g.events.stg.coeff)

outfn <- paste0(dir_exp, "cox_ph.stage_all_genes.txt")
write.table(g.events.stg.coeff, outfn, 
            row.names = TRUE, col.names = NA, quote=FALSE, sep="\t")

