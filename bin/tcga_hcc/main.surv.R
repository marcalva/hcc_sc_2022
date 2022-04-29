
# Perform survival analysis of TCGA data with 
# cell type proportions for main + prolif

setwd("../../")

library(survival)
library(ggplot2)
library(ggfortify)

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
                                     y <- c(cox$n, y, ci)
                                     names(y) <- c("n", "coef", "exp_coef", "se_coef", "z", "p", "lower95", "upper95")
                                     y <- y[c("n", "coef", "exp_coef", "lower95", "upper95", "se_coef", "z", "p")]
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
                                     y <- c(cox$n, n_obs, y, ci)
                                     names(y) <- c("n", "num_events_low", "pct_events_low", "num_events_high", "pct_events_high", 
                                                   "coef", "exp_coef", "se_coef", "z", "p", "lower95", "upper95")
                                     y <- y[c("n", "num_events_low", "pct_events_low", "num_events_high", "pct_events_high", 
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

ct_id <- "cell_type_main"

fn <- "data/processed/tcga_hcc/expr/tcga.gencodev26.rds"
gencode <- readRDS(fn)

# expression data
fn <- "data/processed/tcga_hcc/expr/tcga.lihc.TMM.log.rin.rds"
tmm <- readRDS(fn)
rownames(tmm) <- gencode[ rownames(tmm), "Name"]

# proportions
fn <- paste0("data/processed/tcga_hcc/ctp/tcga.TMM.", ct_id, ".decomp.rds")
ct.md <- readRDS(fn)
ct.mdp <- ct.md$bulk.props

for (ct in names(ct.md$genes.used)){
    ct.md$genes.used[[ct]] <- gencode[ct.md$genes.used[[ct]], "Name"]
}

# marker data
fn <- paste0("exp/sharma_aiz/markers/markers.", ct_id, ".txt")
markers.celltype <- read.table(fn, header=TRUE, stringsAsFactors=FALSE)
ctk <- markers.celltype[,"Name"] %in% rownames(tmm)
markers.celltype <- markers.celltype[ctk,]

# pheno data
fn <- "data/processed/tcga_hcc/sample/cases.hcc.361.merged.rds"
cases <- readRDS(fn)
fn <- "data/processed/tcga_hcc/sample/samples.hcc.410.merged.rds"
samples <- readRDS(fn)
fn <- "data/processed/tcga_hcc/sample/tcga.bio_sample.rds"
bio_sample <- readRDS(fn)

# sample map
fn <- "data/processed/tcga_hcc/sample/sam_case_map.rds"
smaps <- readRDS(fn)
cid2sid_tum <- smaps[["cid2sid_tum"]]
cid2sid_nt <- smaps[["cid2sid_nt"]]
cid2sid_tum_pair <- smaps[["cid2sid_tum_pair"]]
cid2sid_nt_pair <- smaps[["cid2sid_nt_pair"]]
sid_all <- c(cid2sid_tum, cid2sid_nt)
sid_both <- c(cid2sid_tum_pair, cid2sid_nt_pair)

dir_exp <- "exp/tcga_hcc/survival/"
dir.create(dir_exp, showWarnings = FALSE, recursive=TRUE)

ct.mdp <- ct.mdp[,sid_all]
ct.mdp.pt <- ct.mdp[,cid2sid_tum] # tumor

# remove redacted
cases[ is.na(cases[,"cdr.Redaction"]), "cdr.Redaction"] <- "OK"
red.k <- which(cases[,"cdr.Redaction"] == "OK")
cases <- cases[red.k,]

#=====================================================
# Cox proportional hazard
#=====================================================

# format covariates

cases[,"age_at_diagnosis"] <- as.numeric(cases[,"age_at_diagnosis"])

cov.cols <- c("cdr.OS", "cdr.OS.time", "cdr.DSS", "cdr.DSS.time", 
              "cdr.DFI", "cdr.DFI.time", "cdr.PFI", "cdr.PFI.time", 
              "age_at_diagnosis", "gender", "race", "stage_low_high")
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
samples.pt <- cid2sid_tum[rownames(clin.covs.stg)]

# hazards test
# 0=censored 1=event
events <- c("cdr.OS", "cdr.DSS", "cdr.DFI", "cdr.PFI")

#================================================
# Stage + cell type model
#================================================

cox.events.ct_stg.l <- lapply(events, function(event){
                           event.col <- event
                           time.col <- paste0(event, ".time")

                           cell_types <- rownames(ct.mdp.pt)
                           cox.stg.ct.l <- lapply(cell_types, function(ct){
                                                  mat <- cbind(clin.covs.stg, ct.mdp.pt[ct, samples.pt])
                                                  colnames(mat) <- c(colnames(clin.covs.stg), "CellType")
                                                  terms <- ~ age_at_diagnosis + gender + race + stage_low_high + CellType
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
                           names(cox.stg.ct.l) <- cell_types
                           return(cox.stg.ct.l) })
names(cox.events.ct_stg.l) <- events
cox.events.ct_stg.coeff <- cox_list_summary(cox.events.ct_stg.l)

cox.events.ct_stg.coeff <- adjust_p(cox.events.ct_stg.coeff)

outfn <- paste0(dir_exp, ct_id, ".cox_ph.stage_ct.txt")
write.table(cox.events.ct_stg.coeff, outfn, 
            row.names = TRUE, col.names = NA, quote=FALSE, sep="\t")

#================================================
# cell type model
#================================================

cox.events.ct.l <- lapply(events, function(event){
                           event.col <- event
                           time.col <- paste0(event, ".time")

                           cell_types <- rownames(ct.mdp.pt)
                           cox.ct.l <- lapply(cell_types, function(ct){
                                              mat <- cbind(clin.covs.stg, ct.mdp.pt[ct, samples.pt])
                                              colnames(mat) <- c(colnames(clin.covs.stg), "CellType")
                                              terms <- ~ age_at_diagnosis + gender + race + CellType
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
                           names(cox.ct.l) <- cell_types
                           return(cox.ct.l) })
names(cox.events.ct.l) <- events
cox.events.ct.coeff <- cox_list_summary(cox.events.ct.l)

cox.events.ct.coeff <- adjust_p(cox.events.ct.coeff)

outfn <- paste0(dir_exp, ct_id, ".cox_ph.ct.txt")
write.table(cox.events.ct.coeff, outfn, 
            row.names = TRUE, col.names = NA, quote=FALSE, sep="\t")

#================================================
# Stage model
#================================================

cox.events.stg.l <- list()
for (event in events){
    event.col <- event
    time.col <- paste0(event, ".time")

    mat <- clin.covs.stg
    terms <- ~ age_at_diagnosis + gender + race + stage_low_high
    mf <- model.matrix(terms, mat)
    cases.i <- rownames(mf)
    mf <- mf[,-1] # remove intercept
    te <- Surv(time = cases[cases.i,time.col], 
               event = cases[cases.i,event.col], 
               type = "right")
    cx <- coxph(te ~ mf)
    cox.events.stg.l[[event]] <- cx
}

cox.events.stg.coeff <- sapply(cox.events.stg.l, function(cox){
                               cfs <- summary(cox)$coeff
                               y <- cfs[nrow(cfs),] 
                               names(y) <- c("coef", "exp_coef", "se_coef", "z", "p")
                               return(y)}) 

outfn <- paste0(dir_exp, ct_id, ".cox_ph.stage.txt")
write.table(cox.events.stg.coeff, outfn, 
            row.names = TRUE, col.names = NA, quote=FALSE, sep="\t")

#=====================================================
# Cox median (high vs low) cell type only
#=====================================================

cox.events.ct_med.l <- lapply(events, function(event){
                           event.col <- event
                           time.col <- paste0(event, ".time")

                           cell_types <- rownames(ct.mdp.pt)
                           cox.stg.ct.l <- lapply(cell_types, function(ct){
                                                  mat <- cbind(clin.covs.stg, ct.mdp.pt[ct, samples.pt])
                                                  colnames(mat) <- c(colnames(clin.covs.stg), "CellType")

                                                  med.ctp <- quantile(mat[,"CellType"], probs = 0.5)
                                                  ctp.group <- as.numeric(mat[,"CellType"] > med.ctp)
                                                  mat[,"CellTypeGroupHigh1"] <- ctp.group
                                                  
                                                  terms <- ~ age_at_diagnosis + gender + race + CellTypeGroupHigh1
                                                  # terms <- ~ age_at_diagnosis + gender + race + stage_low_high + CellType
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
                                                  cx <- coxph(te ~ mf, na.action=na.omit)

                                                  cx$n_low <- sum(mf[,"CellTypeGroupHigh1"] == 0)
                                                  cx$n_high <- sum(mf[,"CellTypeGroupHigh1"] == 1)
                                                  cx$n_obs_low <- sum(te[mf[,"CellTypeGroupHigh1"] == 0, 2])
                                                  cx$n_obs_high <- sum(te[mf[,"CellTypeGroupHigh1"] == 1, 2])
                                                  return(cx) })
                           names(cox.stg.ct.l) <- cell_types
                           return(cox.stg.ct.l) })
names(cox.events.ct_med.l) <- events
cox.events.ct_med.coeff <- cox_list_summary_grp(cox.events.ct_med.l)

cox.events.ct_med.coeff <- adjust_p(cox.events.ct_med.coeff)

outfn <- paste0(dir_exp, ct_id, ".cox_ph.ct_med.txt")
write.table(cox.events.ct_med.coeff, outfn, 
            row.names = TRUE, col.names = NA, quote=FALSE, sep="\t")

#=====================================================
# Cox stage + median (high vs low) cell type only
#=====================================================

cox.events.stg_ct_med.l <- lapply(events, function(event){
                           event.col <- event
                           time.col <- paste0(event, ".time")

                           cell_types <- rownames(ct.mdp.pt)
                           cox.stg.ct.l <- lapply(cell_types, function(ct){
                                                  mat <- cbind(clin.covs.stg, ct.mdp.pt[ct, samples.pt])
                                                  colnames(mat) <- c(colnames(clin.covs.stg), "CellType")

                                                  med.ctp <- quantile(mat[,"CellType"], probs = 0.5)
                                                  ctp.group <- as.numeric(mat[,"CellType"] > med.ctp)
                                                  mat[,"CellTypeGroupHigh1"] <- ctp.group
                                                  
                                                  terms <- ~ age_at_diagnosis + gender + race + stage_low_high + CellTypeGroupHigh1
                                                  # terms <- ~ age_at_diagnosis + gender + race + stage_low_high + CellType
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
                                                  cx <- coxph(te ~ mf, na.action=na.omit)

                                                  cx$n_low <- sum(mf[,"CellTypeGroupHigh1"] == 0)
                                                  cx$n_high <- sum(mf[,"CellTypeGroupHigh1"] == 1)
                                                  cx$n_obs_low <- sum(te[mf[,"CellTypeGroupHigh1"] == 0, 2])
                                                  cx$n_obs_high <- sum(te[mf[,"CellTypeGroupHigh1"] == 1, 2])
                                                  return(cx) })
                           names(cox.stg.ct.l) <- cell_types
                           return(cox.stg.ct.l) })
names(cox.events.stg_ct_med.l) <- events
cox.events.stg_ct_med.coeff <- cox_list_summary_grp(cox.events.stg_ct_med.l)

cox.events.stg_ct_med.coeff <- adjust_p(cox.events.stg_ct_med.coeff)

outfn <- paste0(dir_exp, ct_id, ".cox_ph.stg_ct_med.txt")
write.table(cox.events.stg_ct_med.coeff, outfn, 
            row.names = TRUE, col.names = NA, quote=FALSE, sep="\t")

#=====================================================
# Cox quartile (high vs low) cell type only
#=====================================================

cox.events.ct_qrt.l <- lapply(events, function(event){
                           event.col <- event
                           time.col <- paste0(event, ".time")

                           cell_types <- rownames(ct.mdp.pt)
                           cox.stg.ct.l <- lapply(cell_types, function(ct){
                                                  mat <- cbind(clin.covs.stg, ct.mdp.pt[ct, samples.pt])
                                                  colnames(mat) <- c(colnames(clin.covs.stg), "CellType")

                                                  qrt.ctp <- quantile(mat[,"CellType"], probs = c(0.25, 0.75))
                                                  ctp.group <- rep(NA, nrow(mat))
                                                  ctp.group[mat[,"CellType"] < qrt.ctp[1]] <- 0
                                                  ctp.group[mat[,"CellType"] > qrt.ctp[2]] <- 1
                                                  mat[,"CellTypeGroupHigh1"] <- ctp.group
                                                  
                                                  terms <- ~ age_at_diagnosis + gender + race + CellTypeGroupHigh1
                                                  # terms <- ~ age_at_diagnosis + gender + race + stage_low_high + CellType
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
                                                  cx <- coxph(te ~ mf, na.action=na.omit)

                                                  cx$n_low <- sum(mf[,"CellTypeGroupHigh1"] == 0)
                                                  cx$n_high <- sum(mf[,"CellTypeGroupHigh1"] == 1)
                                                  cx$n_obs_low <- sum(te[mf[,"CellTypeGroupHigh1"] == 0, 2])
                                                  cx$n_obs_high <- sum(te[mf[,"CellTypeGroupHigh1"] == 1, 2])
                                                  return(cx) })
                           names(cox.stg.ct.l) <- cell_types
                           return(cox.stg.ct.l) })
names(cox.events.ct_qrt.l) <- events
cox.events.ct_qrt.coeff <- cox_list_summary_grp(cox.events.ct_qrt.l)

cox.events.ct_qrt.coeff <- adjust_p(cox.events.ct_qrt.coeff)

outfn <- paste0(dir_exp, ct_id, ".cox_ph.ct_qrt.txt")
write.table(cox.events.ct_qrt.coeff, outfn, 
            row.names = TRUE, col.names = NA, quote=FALSE, sep="\t")

#=====================================================
# Cox decile (high vs low) cell type only
#=====================================================

cox.events.ct_dec.l <- lapply(events, function(event){
                           event.col <- event
                           time.col <- paste0(event, ".time")

                           cell_types <- rownames(ct.mdp.pt)
                           cox.stg.ct.l <- lapply(cell_types, function(ct){
                                                  mat <- cbind(clin.covs.stg, ct.mdp.pt[ct, samples.pt])
                                                  colnames(mat) <- c(colnames(clin.covs.stg), "CellType")

                                                  dec.ctp <- quantile(mat[,"CellType"], probs = c(0.1, 0.9))
                                                  ctp.group <- rep(NA, nrow(mat))
                                                  ctp.group[mat[,"CellType"] < dec.ctp[1]] <- 0
                                                  ctp.group[mat[,"CellType"] > dec.ctp[2]] <- 1
                                                  mat[,"CellTypeGroupHigh1"] <- ctp.group
                                                  
                                                  terms <- ~ age_at_diagnosis + gender + race + CellTypeGroupHigh1
                                                  # terms <- ~ age_at_diagnosis + gender + race + stage_low_high + CellType
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
                                                  cx <- coxph(te ~ mf, na.action=na.omit)

                                                  cx$n_low <- sum(mf[,"CellTypeGroupHigh1"] == 0)
                                                  cx$n_high <- sum(mf[,"CellTypeGroupHigh1"] == 1)
                                                  cx$n_obs_low <- sum(te[mf[,"CellTypeGroupHigh1"] == 0, 2])
                                                  cx$n_obs_high <- sum(te[mf[,"CellTypeGroupHigh1"] == 1, 2])
                                                  return(cx) })
                           names(cox.stg.ct.l) <- cell_types
                           return(cox.stg.ct.l) })
names(cox.events.ct_dec.l) <- events
cox.events.ct_dec.coeff <- cox_list_summary_grp(cox.events.ct_dec.l)

cox.events.ct_dec.coeff <- adjust_p(cox.events.ct_dec.coeff)

outfn <- paste0(dir_exp, ct_id, ".cox_ph.ct_dec.txt")
write.table(cox.events.ct_dec.coeff, outfn, 
            row.names = TRUE, col.names = NA, quote=FALSE, sep="\t")

#=====================================================
# Kaplan-Meier
#=====================================================

# 0=censored 1=event

km.events.l <- lapply(events, function(event){
                      event.col <- event
                      time.col <- paste0(event, ".time")

                      cell_types <- rownames(ct.mdp.pt)
                      km.stg.ct.l <- lapply(cell_types, function(ct){
                                            mat <- cbind(clin.covs.stg, ct.mdp.pt[ct, samples.pt])
                                            colnames(mat) <- c(colnames(clin.covs.stg), "CellType")
                                            med.ctp <- quantile(mat[,"CellType"], probs = 0.5)
                                            ctp.group <- as.numeric(mat[,"CellType"] > med.ctp)
                                            cases.i <- rownames(mat)
                                            te <- Surv(time = cases[cases.i,time.col], 
                                                       event = cases[cases.i,event.col], 
                                                       type = "right")
                                            sur.mod <- survdiff(te ~ ctp.group)
                                            sur.mod$p <- 1 - pchisq(sur.mod$chisq, length(sur.mod$n) - 1)
                                            return(sur.mod) }) 
                       names(km.stg.ct.l) <- cell_types
                       return(km.stg.ct.l) })
names(km.events.l) <- events

km.events.l <- lapply(names(km.events.l), function(n){
                    km.l <- km.events.l[[n]]
                    km <- sapply(km.l, function(x){
                                 ret <- data.frame("n.low" = x$n[1], 
                                                   "n.high" = x$n[2],
                                                   "obs.low" = x$obs[1], 
                                                   "obs.high" = x$obs[2], 
                                                   "chisq" = x$chisq, 
                                                   "p" = x$p)
                                 colnames(ret) <- paste0(n, ".", colnames(ret))
                                 return(ret) })
                    return(km) })
names(km.events.l) <- events

km.events <- do.call(rbind, km.events.l)
km.events <- t(km.events)

outfn <- paste0(dir_exp, ct_id, ".km.ct.txt")
write.table(km.events, outfn, 
            row.names = TRUE, col.names = NA, quote=FALSE, sep="\t")

