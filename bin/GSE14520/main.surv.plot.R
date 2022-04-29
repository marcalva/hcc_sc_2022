
# Perform survival analysis of LCI  data with 
# cell type proportions.

setwd("../../")

library(survival)
library(ggplot2)
library(ggfortify)
library(reshape2)

#=====================================================
# functions
#=====================================================

#=====================================================
#=====================================================

ct_id <- "cell_type_main"
cl_name <- "cluster"
ct_rm <- c("B")

events <- c("Survival", "Recurr")
eshort <- c("OS", "RCR")
names(eshort) <- events

dir_exp <- paste0("exp/GSE14520/survival/", ct_id, "/")
dir.create(dir_exp, showWarnings = FALSE, recursive=TRUE)

cl_id <- "cell_type_main"

# read in data
dir_data <- "data/processed/GSE14520/"
dir.create(dir_data, recursive=TRUE, showWarning=FALSE)

fn <- "data/processed/tcga_hcc/expr/tcga.gencodev26.rds"
gencode <- readRDS(fn)

# expression data
fn <- paste0(dir_data, "expr.RMA_log2.gmean.rds")
ex <- readRDS(fn)

fn <- paste0(dir_data, "sdata.txt")
samples <- read.table(fn, row.names=1, header=TRUE, sep="\t", stringsAsFactors=FALSE)

fn <- paste0(dir_data, "geo_pheno.txt")
geo_pheno <- read.table(fn, row.names=1, header=TRUE, sep="\t", stringsAsFactors=FALSE)

# cell-type proportions
fn <- paste0("data/processed/GSE14520/ctp/", "LCI.", cl_id, ".decomp.rds")
ex.ct.md <- readRDS(fn)

ex.ct.mdp <- ex.ct.md$bulk.props
k <- ! rownames(ex.ct.mdp) %in% ct_rm
ex.ct.mdp <- ex.ct.mdp[k,,drop=FALSE]

# colors
fn <- "data/ref/colors/tumor_colrs.csv"
tum_cols <- read.table(fn, header=FALSE, sep=",")
tum_colsv <- tum_cols[,2]
names(tum_colsv) <- tum_cols[,1]
names(tum_colsv)[names(tum_colsv) == "Tumor"] <- "Tumor"
names(tum_colsv)[names(tum_colsv) == "NonTumor"] <- "Non-tumor"

fn <- "data/ref/colors/rd_bu_div.csv"
rd_bu <- read.table(fn)
rd_bu <- rd_bu[,1]

# marker data
fn <- paste0("exp/sharma_aiz/markers/markers.", cl_id, ".txt")
markers.celltype <- read.table(fn, header=TRUE, stringsAsFactors=FALSE)
ctk <- markers.celltype[,"Name"] %in% rownames(ex)
markers.celltype <- markers.celltype[ctk,]
ctk <- !markers.celltype[,cl_name] %in% ct_rm
markers.celltype <- markers.celltype[ctk,]
k <- markers.celltype[,"p_val_adj"] < .05
markers.celltype <- markers.celltype[k,]
celltypes <- sort(unique(markers.celltype[,cl_name]))
names(celltypes) <- celltypes
gene_col <- "Name"

# sample map
fn <- "data/processed/GSE14520/sam_case_map.rds"
smaps <- readRDS(fn)
cid2sid_tum <- smaps[["cid2sid_tum"]]
cid2sid_nt <- smaps[["cid2sid_nt"]]
cid2sid_tum_pair <- smaps[["cid2sid_tum_pair"]]
cid2sid_nt_pair <- smaps[["cid2sid_nt_pair"]]
sid_all <- c(cid2sid_tum, cid2sid_nt)
sid_both <- c(cid2sid_tum_pair, cid2sid_nt_pair)

# subset expression and props
ex <- ex[,sid_all]
ex.pt <- ex[,cid2sid_tum]
ex.ct.mdp <- ex.ct.mdp[,sid_all]
ex.ct.mdp.pt <- ex.ct.mdp[,cid2sid_tum]

# subset markers to expressed genes
k <- markers.celltype[,gene_col] %in% rownames(ex)
markers.celltype <- markers.celltype[k,,drop=FALSE]

# read in cell-type prop survival associations
fn <- paste0("exp/GSE14520/survival/", ct_id, ".", "GSE14520.cox_ph.ct_med.txt")
ct_med.events.coeff <- read.table(fn, header=TRUE, row.names=1, stringsAsFactors = FALSE, sep="\t")
ctk <- ! rownames(ct_med.events.coeff) %in% ct_rm
ct_med.events.coeff <- ct_med.events.coeff[ctk,]

# Read survival results for all genes
fn <- "exp/GSE14520/survival.gene/GSE14520.cox_ph.all_genes.txt"
g.events.coeff <- read.table(fn, header=TRUE, row.names=1, stringsAsFactors = FALSE, sep="\t")

# subset into list of cell type marker coefficients
ct.coef.l <- lapply(celltypes, function(ct){
                    genes <- markers.celltype[markers.celltype[,cl_name] == ct, gene_col]
                    ct.coef <- g.events.coeff[genes,]
                    ct.coef[,"CellType"] <- ct
                    return(ct.coef) })

padjcolf <- function(s) paste0(s, ".p_adj")
expcolf <- function(s) paste0(s, ".exp_coef")

# surv data
survdat <- samples[,c("Survival.status", "Survival.months", 
                      "Recurr.status", "Recurr.months")]
colnames(survdat) <- c(paste0(eshort[1], c(".status", ".months")), 
                       paste0(eshort[2], c(".status", ".months")))

#=====================================================
# process data for survival analysis
#=====================================================

# separate non-tumor and tumor
ex.ct.mdp.tum <- ex.ct.mdp[,unname(cid2sid_tum)]
ex.ct.mdp.nt <- ex.ct.mdp[,unname(cid2sid_nt)]

# tumor sample data
samples.tum <- samples[unname(cid2sid_tum),]

# format covariates
cov.cols <- c("Age", "Gender", "stage_low_high")
clin.covs.stg <- samples.tum[,cov.cols]

# change stage to factor
# ensures low=0 high=1
clin.covs.stg[,"stage_low_high"] <- factor(clin.covs.stg[,"stage_low_high"], 
                                           levels = c("Low", "High"))

# no stage for covariates
colrm <- which(colnames(clin.covs.stg) == "stage_low_high")
clin.covs <- clin.covs.stg[,-colrm]

ex.ct.mdp.tum <- ex.ct.mdp.tum[,rownames(clin.covs.stg)]

# sample IDs
samples.pt <- rownames(clin.covs.stg)

#================================================
# Plot % gws markers: all gene expression model
#================================================

get_pct_stats <- function(coefs, pcol, ecol){
    sig.n <- sum(coefs[,padjcolf(e)] < 0.05, na.rm=TRUE)
    pos.n <- sum(coefs[,expcolf(e)] > 1, na.rm=TRUE)
    ct.n <- nrow(coefs)
    sig.pct <- sig.n / ct.n
    pos.pct <- pos.n / ct.n
    d <- data.frame("n_gws" = sig.n, 
                    "n_pos" = pos.n, 
                    "pct_gws" = sig.pct, 
                    "pct_pos" = pos.pct, 
                    "n_ct" = ct.n)
    return(d)
}

pct.gws_l <- list()
for (ct in celltypes){
    ct.coef <- ct.coef.l[[ct]]
    e.coef_l <- list()
    for (e in events){
        d <- get_pct_stats(ct.coef, pcol = padjcolf(e), ecol = expcolf(e))
        d[,"CellType"] <- ct
        d[,"event"] <- e
        e.coef_l[[e]] <- d
    }
    pct.gws_l[[ct]] <- do.call(rbind, e.coef_l)
}
pct.gws <- do.call(rbind, pct.gws_l)

k <- pct.gws[,"event"] == events[1]
o <- order(pct.gws[k,"pct_gws"], decreasing = TRUE)
cto <- pct.gws[k,][o,"CellType"]
pct.gws[,"CellType"] <- factor(pct.gws[,"CellType"], levels = cto)

# barplot of number of gws
p <- ggplot(pct.gws, aes(x = CellType, y = pct_gws)) + 
    geom_bar(stat = "identity") + 
    facet_wrap(~ event, ncol=1) + 
    theme_bw() + xlab("Cell Type") + 
    ylab("Percent of cell-type markers\naffecting hazard rate") + 
    ggtitle("Gene-based Proportional\nHazards Model (TCGA)") + 
    scale_y_continuous(labels = function(x) paste0(x*100, "%"), 
                       limits = c(0, 0.8)) + 
    theme(axis.text.x = element_text(angle=45, hjust = 0.9, vjust=1), 
          plot.title = element_text(hjust = 0.5, size = 10), 
          text = element_text(size = 10))

outfn <- paste0(dir_exp, "gene.num_gws.", paste(eshort, collapse="_"), ".pdf")
ggsave(outfn, height = 4, width = 3)
outfn <- paste0(dir_exp, "gene.surv.", paste(eshort, collapse="_"), ".txt")
write.table(pct.gws, outfn, row.names=TRUE, col.names=NA, quote=F, sep="\t")

#================================================
# get HR for cell-type markers and format data fram
#================================================

coefs_ctm_l <- list()
for (ct in celltypes){
    coefs <- ct.coef.l[[ct]]
    e_l <- list()
    for (e in events){
        e_l[[e]] <- data.frame("HR" = coefs[,expcolf(e)], 
                               "p" = coefs[,padjcolf(e)], 
                               "CellType" = ct, 
                               "event" = unname(eshort[e]), 
                               "gene" = rownames(coefs))
    }
    coefs_ctm_l[[ct]] <- do.call(rbind, e_l)
}
coefs_ctm <- do.call(rbind, coefs_ctm_l)
coefs_ctm[,"HR"] <- log2(coefs_ctm[,"HR"])
coefs_ctm[,"p"] <- -10 * log10(coefs_ctm[,"p"])
coefs_ctm[,"event"] <- factor(coefs_ctm[,"event"], eshort)

#================================================
# Density plot of hazard ratio for cell type marker genes
#================================================

stheme <- theme(plot.title = element_text(size = 8, hjust = 0.5),
                text = element_text(size = 8), 
                strip.text = element_text(size = 6), 
                strip.background = element_rect(fill = "white"))

p <- ggplot(coefs_ctm, aes(x = HR, color = CellType)) + 
    geom_density(alpha = 0) + 
    scale_x_continuous(limits = c(-5, 5)) + 
    geom_vline(xintercept = 0, col = "red") + 
    facet_wrap(~event, ncol=1) + 
    xlab(expression(paste(log[2], " hazard ratio"))) + 
    ylab("Density") + 
    ggtitle("Marker gene hazard ratio\n(LCI)") + 
    theme_bw() + 
    theme(text = element_text(size = 8), 
          plot.title = element_text(hjust = 0.5, size = 10))
outfn <- paste0(dir_exp, "gene.HR_density.", paste(eshort, collapse="_"), ".pdf")
ggsave(outfn, height = 4, width = 3.5)

#================================================
# Plot Prolif marker log2 HRs
#================================================

k <- coefs_ctm[,"CellType"] == "Prol"
coefs_ctm.ct <- coefs_ctm[k,]

# order gene by decreasing HR
k <- coefs_ctm.ct[,"event"] == eshort[events][1]
o <- order(coefs_ctm.ct[k,"HR"], decreasing=FALSE)
coefs_ctm.ct[,"gene"] <- factor(coefs_ctm.ct[,"gene"], 
                                  levels = coefs_ctm.ct[k,"gene"][o])

p <- ggplot(coefs_ctm.ct, aes(x = gene, y = HR)) + 
    geom_point(size = 0.5) + 
    coord_flip() + 
    scale_y_continuous(limits = c(-6, 6)) + 
    geom_hline(yintercept = 0, col = "red") + 
    facet_wrap(~event, ncol = 2) + 
    ylab(expression(paste(log[2], " hazard ratio"))) + 
    xlab("Prol markers") +
    ggtitle("LCI") + 
    theme_bw() + 
    theme(axis.text.y = element_blank(), 
          strip.background = element_rect(fill = "white"),
          plot.title = element_text(hjust = 0.5, size = 10), 
          axis.ticks.y = element_blank()) 
outfn <- paste0(dir_exp, "gene.HR.Prol.pdf")
ggsave(outfn, height = 3.5, width = 3.5)

#================================================
# Forest plot of HR for each cell type: CT med model
#================================================

th_hr <- theme_bw() + 
theme(text = element_text(size = 10), 
      plot.title = element_text(hjust = 0.5), 
      legend.key.width = unit(0, "in"), 
      legend.position = "bottom")

hr.df.l <- list()
for (e in events){
    ck <- paste0(e, c(".exp_coef", ".lower95", ".upper95", ".p", ".p_adj"))
    hr.df <- ct_med.events.coeff[,ck]
    colnames(hr.df) <- c("HR", "lower", "upper", "p", "p_adj")
    hr.df[,"CellType"] <- rownames(hr.df)
    hr.df[,"Event"] <- e
    hr.df.l[[e]] <- hr.df
}
hr.df <- do.call(rbind, hr.df.l)

logcol <- c("HR", "lower", "upper")
hr.df[,logcol] <- log2(hr.df[,logcol])

hr.df[,"sig"] <- "ns"
hr.df[hr.df[,"p_adj"] < 0.05, "sig"] <- "sig"
sig_shape <- c("ns" = 1, "sig" = 16)

o <- order(hr.df.l[[1]][,"HR"], decreasing = FALSE)
hr.df[,"CellType"] <- factor(hr.df[,"CellType"], levels = hr.df.l[[1]][o,"CellType"])
hr.df[,"Event"] <- factor(hr.df[,"Event"], levels = events)

p <- ggplot(hr.df, aes(xmin = lower, xmax = upper, y = CellType)) + 
    facet_wrap(~ Event, nrow=1) + 
    geom_errorbarh(alpha=0.5, height = 0.5, color="black") + 
    geom_vline(xintercept = 0.0, linetype=2, alpha=0.75) +
    geom_point(aes(x = HR, shape = sig), size = 3) + 
    scale_shape_manual(values = sig_shape, breaks = "sig", 
                       name = NULL, 
                       labels = c("sig" = "adj p < 0.05")) + 
    xlab(expression(log[2]~Hazard~Ratio)) + 
    ylab(NULL) + 
    ggtitle("Cell-type proportional\nHazards Model (LCI)")
outfn <- paste0(dir_exp, "ct_med.HR_forest.", paste(eshort, collapse="_"), ".pdf")
ggsave(outfn, height = 4, width = 6)


# individual panel plots
for (event in events){
    ev_shrt <- eshort[event]
    hr.df <- hr.df.l[[event]]

    hr.df[,"sig"] <- "ns"
    hr.df[hr.df[,"p_adj"] < 0.05, "sig"] <- "sig"
    sig_shape <- c("ns" = 1, "sig" = 16)

    o <- order(hr.df[,"HR"], decreasing = FALSE)
    hr.df[,"CellType"] <- factor(hr.df[,"CellType"], levels = hr.df.l[[1]][o,"CellType"])
    hr.df[,"Event"] <- factor(hr.df[,"Event"], levels = events)

    p <- ggplot(hr.df, aes(xmin = lower, xmax = upper, y = CellType)) + 
    geom_errorbarh(alpha = 0.5, height = 0.5, color = "black") + 
    geom_vline(xintercept = 1, linetype = 2, alpha = 0.75) + 
    geom_point(aes(x = HR, shape = sig), size = 3) + 
    scale_shape_manual(values = sig_shape, breaks = "sig", name = NULL, 
                       labels = c("sig" = "adj p < 0.05")) +
    xlab(paste0("Hazard ratio for ", event)) + 
    ylab("Cell-type") + 
    ggtitle(paste0("Cell-type proportional\nhazards model (LCI)\n", ev_shrt)) + 
    th_hr

    out_fn <- paste0(dir_exp, "ct_med.HR_forest.", ev_shrt, ".pdf")
    ggsave(out_fn, height = 4, width = 3)
}


#=====================================================
# Kaplan-Meier tests and plots
#=====================================================

dir_plot <- dir_exp
dir.create(dir_plot, showWarnings = FALSE, recursive=TRUE)

cell_types <- rownames(ex.ct.mdp.pt)
# get survfit model for plotting
km.events.ct.l <- list()
for (event in eshort){
    event.col <- paste0(event, ".status")
    time.col <- paste0(event, ".months")

    km.ct.l <- list()
    for (ct in cell_types){
        ct.p <- ex.ct.mdp.pt[ct, samples.pt]
        med.ctp <- quantile(ct.p, probs = c(0.5, 0.5))

        ctp.group <- as.character(as.numeric(ct.p > med.ctp[1]))
        lh.key <- c("0" = "Low", "1" = "High")
        ctp.group <- lh.key[ctp.group]
        ctp.group <- factor(ctp.group, levels = c("Low", "High"))
        names(ctp.group) <- names(ct.p)
        ctp.cid <- names(ctp.group)

        edf <- data.frame("ctg" = ctp.group, 
                          "time" = survdat[ctp.cid, time.col], 
                          "event" = survdat[ctp.cid, event.col])
        colnames(edf) <- c("ctg", time.col, event.col)
        k <- apply(edf, 1, function(rw) !any(is.na(rw)))
        k <- rownames(edf)[k]
        edf.k <- edf[k,,drop=FALSE]

        te <- Surv(time = edf.k[,time.col], 
                   event = edf.k[,event.col], 
                   type = "right")
        group <- edf.k[,"ctg"]
        kmfit <- survfit(te ~ group, na.action=na.omit)
        km.ct.l[[ct]] <- kmfit
    }
    km.events.ct.l[[event]] <- km.ct.l
}

# save survfit objects for plotting
out_fn <- paste0(dir_exp, "km.events.ct.l.rds")
saveRDS(km.events.ct.l, out_fn)

# surv theme
th_surv <- 
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5), 
          legend.text = element_text(size = 8), 
          legend.title = element_text(size = 8, hjust = 0.5), 
          legend.key.size = unit(.05, units = "npc"))

event <- "OS"
pt <- "Overall Survival"
pt <- "OS"
ct <- "Prol"
yticks <- paste0(seq(0, 100, 20), "%")

for (ct in c("Prol")){
    ph <- 3; pw <- 3.5
    surv.md <- km.events.ct.l[[event]][[ct]]
    surv.md$time <- surv.md$time / 12

    p <- autoplot(surv.md, censor.size = 2) + theme_bw() + 
    ggtitle(paste0(pt, " in ", ct, "\n(LCI)")) + 
    xlab("Years") + ylab("Percent Survival") + 
    scale_x_continuous(breaks = seq(0,10,2)) + 
    scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2), labels=yticks) + 
    guides(color = "none", fill = guide_legend(title=paste0(ct, "\nMedian\nGroup"))) + 
    th_surv

    outfn <- paste0(dir_exp, "curve.", event, ".", ct, ".pdf")
    p <- p + scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2), labels=yticks)
    ggsave(outfn, height=ph, width=pw)
}

event <- "RCR"
pt <- "Recurrence"
pt <- "RCR"
ct <- "Prol"
for (ct in c("Prol")){
    surv.md <- km.events.ct.l[[event]][[ct]]
    surv.md$time <- surv.md$time / 12

    p <- autoplot(surv.md, censor.size = 2) + theme_bw() + 
    ggtitle(paste0(pt, " in ", ct, "\n(LCI)")) + 
    xlab("Years") + ylab("Percent Progression Free") + 
    scale_x_continuous(breaks = seq(0,10,2)) + 
    scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2), labels=yticks) + 
    guides(color = "none", fill = guide_legend(title=paste0(ct, "\nMedian\nGroup"))) + 
    th_surv

    outfn <- paste0(dir_exp, "curve.", event, ".", ct, ".pdf")
    ggsave(outfn, height=ph, width=pw)
}

