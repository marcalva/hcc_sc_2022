
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

dir_exp <- "exp/GSE14520/survival/"
dir.create(dir_exp, showWarnings = FALSE, recursive=TRUE)

dir_data <- "data/processed/GSE14520/"
dir.create(dir_data, recursive=TRUE, showWarning=FALSE)

fn <- "data/processed/tcga_hcc/expr/tcga.gencodev26.rds"
gencode <- readRDS(fn)

fn <- paste0(dir_data, "expr.RMA_log2.gmean.rds")
ex <- readRDS(fn)

fn <- paste0(dir_data, "sdata.txt")
sdat <- read.table(fn, row.names=1, header=TRUE, sep="\t", stringsAsFactors=FALSE)

fn <- paste0(dir_data, "geo_pheno.txt")
geo_pheno <- read.table(fn, row.names=1, header=TRUE, sep="\t", stringsAsFactors=FALSE)

# marker data
fn <- paste0("exp/sharma_aiz/markers/markers.", ct_id, ".txt")
markers.celltype <- read.table(fn, header=TRUE, stringsAsFactors=FALSE)
ctk <- markers.celltype[,"Name"] %in% rownames(ex)
markers.celltype <- markers.celltype[ctk,]

# cell-type proportions
fn <- paste0("data/processed/GSE14520/ctp/", "LCI.", ct_id, ".decomp.rds")
ex.ct.md <- readRDS(fn)

ex.ct.mdp <- ex.ct.md$bulk.props

# sample map
fn <- "data/processed/GSE14520/sam_case_map.rds"
smaps <- readRDS(fn)
cid2sid_tum <- smaps[["cid2sid_tum"]]
cid2sid_nt <- smaps[["cid2sid_nt"]]
cid2sid_tum_pair <- smaps[["cid2sid_tum_pair"]]
cid2sid_nt_pair <- smaps[["cid2sid_nt_pair"]]

#=====================================================
# Format data
#=====================================================

# separate non-tumor and tumor
ex.ct.mdp.tum <- ex.ct.mdp[,unname(cid2sid_tum)]
ex.ct.mdp.nt <- ex.ct.mdp[,unname(cid2sid_nt)]

# tumor sample data
sdat.tum <- sdat[unname(cid2sid_tum),]

# format covariates
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

ex.ct.mdp.tum <- ex.ct.mdp.tum[,rownames(clin.covs.stg)]


#================================================
# Stage + cell type model
#================================================


cox.events.ct_stg.l <- lapply(events, function(event){
                           event.col <- paste0(event, ".status")
                           time.col <- paste0(event, ".months")

                           cell_types <- rownames(ex.ct.mdp.tum)
                           cox.stg.ct.l <- lapply(cell_types, function(ct){
                                                  mat <- cbind(clin.covs.stg, ex.ct.mdp.tum[ct, ])
                                                  colnames(mat) <- c(colnames(clin.covs.stg), "CellType")
                                                  terms <- ~ Age + Gender + stage_low_high + CellType
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
                           names(cox.stg.ct.l) <- cell_types
                           return(cox.stg.ct.l) })
names(cox.events.ct_stg.l) <- events
cox.events.ct_stg.coeff <- cox_list_summary(cox.events.ct_stg.l)

cox.events.ct_stg.coeff <- adjust_p(cox.events.ct_stg.coeff)

outfn <- paste0(dir_exp, ct_id, ".GSE14520.cox_ph.stage_ct.txt")
write.table(cox.events.ct_stg.coeff, outfn, 
            row.names = TRUE, col.names = NA, quote=FALSE, sep="\t")

#================================================
# cell type model
#================================================

cox.events.ct.l <- lapply(events, function(event){
                           event.col <- paste0(event, ".status")
                           time.col <- paste0(event, ".months")

                           cell_types <- rownames(ex.ct.mdp.tum)
                           cox.ct.l <- lapply(cell_types, function(ct){
                                                  mat <- cbind(clin.covs, ex.ct.mdp.tum[ct, ])
                                                  colnames(mat) <- c(colnames(clin.covs), "CellType")
                                                  terms <- ~ Age + Gender + CellType
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
                           names(cox.ct.l) <- cell_types
                           return(cox.ct.l) })
names(cox.events.ct.l) <- events
cox.events.ct.coeff <- cox_list_summary(cox.events.ct.l)

cox.events.ct.coeff <- adjust_p(cox.events.ct.coeff)

outfn <- paste0(dir_exp, ct_id, ".GSE14520.cox_ph.ct.txt")
write.table(cox.events.ct.coeff, outfn, 
            row.names = TRUE, col.names = NA, quote=FALSE, sep="\t")

#================================================
# Stage model
#================================================

cox.events.stg.l <- lapply(events, function(event){
                           event.col <- paste0(event, ".status")
                           time.col <- paste0(event, ".months")

                           mat <- clin.covs.stg
                           terms <- ~ Age + Gender + stage_low_high
                           mf <- model.matrix(terms, mat)
                           ids.i <- rownames(mf)
                           mf <- mf[,-1] # remove intercept
                           te <- Surv(time = sdat.tum[ids.i,time.col], 
                                      event = sdat.tum[ids.i,event.col], 
                                      type = "right")
                           cx <- coxph(te ~ mf)
                           return(cx) })
names(cox.events.stg.l) <- events

cox.events.stg.coeff <- sapply(cox.events.stg.l, function(cox){
                               cfs <- summary(cox)$coeff
                               y <- cfs[nrow(cfs),] 
                               names(y) <- c("coef", "exp_coef", "se_coef", "z", "p")
                               return(y)}) 

outfn <- paste0(dir_exp, "GSE14520.cox_ph.stage.txt")
write.table(cox.events.stg.coeff, outfn, 
            row.names = TRUE, col.names = NA, quote=FALSE, sep="\t")

#=====================================================
# Cox median (high vs low) cell type only
#=====================================================


cox.events.ct_med.l <- lapply(events, function(event){
                           event.col <- paste0(event, ".status")
                           time.col <- paste0(event, ".months")

                           cell_types <- rownames(ex.ct.mdp.tum)
                           cox.ct.l <- lapply(cell_types, function(ct){
                                                  mat <- cbind(clin.covs, ex.ct.mdp.tum[ct, ])
                                                  colnames(mat) <- c(colnames(clin.covs), "CellType")
                                                  terms <- ~ Age + Gender + CellType
                                                  mf <- model.matrix(terms, mat)
                                                  ids.i <- rownames(mf)

                                                  med.ctp <- quantile(mat[,"CellType"], probs = 0.5)
                                                  ctp.group <- as.numeric(mat[,"CellType"] > med.ctp)
                                                  mat[,"CellTypeGroupHigh1"] <- ctp.group
                                                  
                                                  terms <- ~ Age + Gender + CellTypeGroupHigh1
                                                  mf <- model.matrix(terms, mat)
                                                  ids.i <- rownames(mf)

                                                  edf <- sdat.tum[ids.i,c(time.col, event.col)]
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
                           names(cox.ct.l) <- cell_types
                           return(cox.ct.l) })
names(cox.events.ct_med.l) <- events
cox.events.ct_med.coeff <- cox_list_summary_grp(cox.events.ct_med.l)

cox.events.ct_med.coeff <- adjust_p(cox.events.ct_med.coeff)

outfn <- paste0(dir_exp, ct_id, ".GSE14520.cox_ph.ct_med.txt")
write.table(cox.events.ct_med.coeff, outfn, 
            row.names = TRUE, col.names = NA, quote=FALSE, sep="\t")

#=====================================================
# Cox stage + median (high vs low) cell type only
#=====================================================


cox.events.stg_ct_med.l <- lapply(events, function(event){
                           event.col <- paste0(event, ".status")
                           time.col <- paste0(event, ".months")

                           cell_types <- rownames(ex.ct.mdp.tum)
                           cox.stg.ct.l <- lapply(cell_types, function(ct){
                                                  mat <- cbind(clin.covs.stg, ex.ct.mdp.tum[ct, ])
                                                  colnames(mat) <- c(colnames(clin.covs.stg), "CellType")
                                                  terms <- ~ Age + Gender + stage_low_high + CellType
                                                  mf <- model.matrix(terms, mat)
                                                  ids.i <- rownames(mf)

                                                  med.ctp <- quantile(mat[,"CellType"], probs = 0.5)
                                                  ctp.group <- as.numeric(mat[,"CellType"] > med.ctp)
                                                  mat[,"CellTypeGroupHigh1"] <- ctp.group
                                                  
                                                  terms <- ~ Age + Gender + stage_low_high + CellTypeGroupHigh1
                                                  mf <- model.matrix(terms, mat)
                                                  ids.i <- rownames(mf)

                                                  edf <- sdat.tum[ids.i,c(time.col, event.col)]
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

outfn <- paste0(dir_exp, ct_id, ".GSE14520.cox_ph.stg_ct_med.txt")
write.table(cox.events.stg_ct_med.coeff, outfn, 
            row.names = TRUE, col.names = NA, quote=FALSE, sep="\t")

#=====================================================
# Cox quartile (high vs low) cell type only
#=====================================================

cox.events.ct_qrt.l <- lapply(events, function(event){
                           event.col <- paste0(event, ".status")
                           time.col <- paste0(event, ".months")

                           cell_types <- rownames(ex.ct.mdp.tum)
                           cox.ct.l <- lapply(cell_types, function(ct){
                                                  mat <- cbind(clin.covs, ex.ct.mdp.tum[ct, ])
                                                  colnames(mat) <- c(colnames(clin.covs), "CellType")
                                                  terms <- ~ Age + Gender + CellType
                                                  mf <- model.matrix(terms, mat)
                                                  ids.i <- rownames(mf)

                                                  qrt.ctp <- quantile(mat[,"CellType"], probs = c(0.25, 0.75))
                                                  ctp.group <- rep(NA, nrow(mat))
                                                  ctp.group[mat[,"CellType"] < qrt.ctp[1]] <- 0
                                                  ctp.group[mat[,"CellType"] > qrt.ctp[2]] <- 1
                                                  mat[,"CellTypeGroupHigh1"] <- ctp.group
                                                  
                                                  terms <- ~ Age + Gender + CellTypeGroupHigh1
                                                  mf <- model.matrix(terms, mat)
                                                  ids.i <- rownames(mf)

                                                  edf <- sdat.tum[ids.i,c(time.col, event.col)]
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
                           names(cox.ct.l) <- cell_types
                           return(cox.ct.l) })
names(cox.events.ct_qrt.l) <- events
cox.events.ct_qrt.coeff <- cox_list_summary_grp(cox.events.ct_qrt.l)

cox.events.ct_qrt.coeff <- adjust_p(cox.events.ct_qrt.coeff)

outfn <- paste0(dir_exp, ct_id, ".GSE14520.cox_ph.ct_qrt.txt")
write.table(cox.events.ct_qrt.coeff, outfn, 
            row.names = TRUE, col.names = NA, quote=FALSE, sep="\t")

#=====================================================
# Kaplan-Meier
#=====================================================

# 0=censored 1=event

km.events.l <- lapply(events, function(event){
                      event.col <- paste0(event, ".status")
                      time.col <- paste0(event, ".months")

                      cell_types <- rownames(ex.ct.mdp.tum)
                      km.stg.ct.l <- lapply(cell_types, function(ct){
                                            mat <- cbind(clin.covs.stg, ex.ct.mdp.tum[ct,])
                                            colnames(mat) <- c(colnames(clin.covs.stg), "CellType")
                                            med.ctp <- quantile(mat[,"CellType"], probs = c(0.5, 0.5))

                                            ctp.group <- rep(NA, nrow(mat))
                                            ctp.group[mat[,"CellType"] <= med.ctp[1]] <- 0
                                            ctp.group[mat[,"CellType"] > med.ctp[2]] <- 1
                                            names(ctp.group) <- rownames(mat)
                                                  
                                            edf <- sdat.tum[names(ctp.group),c(time.col, event.col)]
                                            k <- apply(edf, 1, function(rw) !any(is.na(rw)))
                                            k <- rownames(edf)[k]
                                            edf <- edf[k,,drop=FALSE]

                                            te <- Surv(time = edf[,time.col], 
                                                       event = edf[,event.col], 
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

