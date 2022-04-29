
# Perform survival analysis of TCGA data with 
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

dir_exp <- paste0("exp/tcga_hcc/survival/", ct_id, "/")
dir.create(dir_exp, showWarnings = FALSE, recursive=TRUE)

# read in data
out_dir <- "data/processed/tcga_hcc/expr/"

# annotation data
fn <- paste0(out_dir, "tcga.gencodev26.rds")
gencode <- readRDS(fn)

# expression data
fn <- "data/processed/tcga_hcc/expr/tcga.lihc.TMM.log.rin.rds"
tmm <- readRDS(fn)

# marker data
fn <- paste0("exp/sharma_aiz/markers/markers.", ct_id, ".txt")
markers.celltype <- read.table(fn, header=TRUE, stringsAsFactors=FALSE)
markers.celltype[,"gene"] <- sapply(markers.celltype[,"gene"], 
                                    function(s){
                                        strsplit(s, "\\.")[[1]][1] })

k <- markers.celltype[,"p_val_adj"] < .05
markers.celltype <- markers.celltype[k,]
celltypes <- sort(unique(markers.celltype[,cl_name]))
names(celltypes) <- celltypes

# proportions
fn <- "data/processed/tcga_hcc/ctp/tcga.TMM.cell_type_main.decomp.rds"
tmm.ct.md <- readRDS(fn)
tmm.ct.mdp <- tmm.ct.md$bulk.props

for (ct in names(tmm.ct.md$markers)){
    tmm.ct.md$markers[[ct]] <- gencode[tmm.ct.md$markers[[ct]], "Name"]
}

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
tmm.pt <- tmm[,cid2sid_tum]
tmm.ct.mdp <- tmm.ct.mdp[,sid_all]
tmm.ct.mdp.pt <- tmm.ct.mdp[,cid2sid_tum]

# remove redacted
cases[ is.na(cases[,"cdr.Redaction"]), "cdr.Redaction"] <- "OK"
red.k <- which(cases[,"cdr.Redaction"] == "OK")
cases <- cases[red.k,]

# subset markers to expressed genes
k <- markers.celltype[,"gene"] %in% rownames(tmm)
markers.celltype <- markers.celltype[k,,drop=FALSE]

#=====================================================
# Cox proportional hazard
#=====================================================

# format covariates

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


# hazards test
# 0=censored 1=event
events <- c("cdr.OS", "cdr.DSS", "cdr.DFI", "cdr.PFI")

#================================================
# Read survival results
#================================================

fn <- "exp/tcga_hcc/survival.gene/cox_ph.all_genes.txt"
g.events.coeff <- read.table(fn, header=TRUE, row.names=1, stringsAsFactors = FALSE, sep="\t")
# get cell type marker coefficients
ct.coef.l <- lapply(celltypes, function(ct){
                    genes <- markers.celltype[markers.celltype[,cl_name] == ct, "gene"]
                    ct.coef <- g.events.coeff[genes,]
                    ct.coef[,"CellType"] <- ct
                    return(ct.coef) })

#================================================
# Plot % gws markers: all gene expression model
#================================================

pct.gws <- lapply(celltypes, function(ct){
                  ct.coef <- ct.coef.l[[ct]]
                  ret <- lapply(c("OS", "PFI"), function(e){
                                sig.n <- sum(ct.coef[,paste0("cdr.", e, ".p_adj")] < 0.05, na.rm=TRUE)
                                pos.n <- sum(ct.coef[,paste0("cdr.", e, ".exp_coef")] > 1)
                                ct.n <- nrow(ct.coef)
                                sig.pct <- sig.n / ct.n
                                pos.pct <- pos.n / ct.n
                                return(data.frame("CellType" = ct, 
                                                  "event" = e, 
                                                  "n_gws" = sig.n, 
                                                  "n_pos" = pos.n, 
                                                  "pct_gws" = sig.pct, 
                                                  "pct_pos" = pos.pct, 
                                                  "n_ct" = ct.n)) })
                  ret <- do.call(rbind, ret)
                  return(ret) })
pct.gws <- do.call(rbind, pct.gws)

pct.gws.os <- subset(pct.gws, event = "OS")
o <- order(pct.gws.os[,"pct_gws"], decreasing = TRUE)
cto <- unique(pct.gws.os[o,"CellType"])
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

outfn <- paste0(dir_exp, "gene.num_gws.OS_PFI.pdf")
ggsave(outfn, height = 4, width = 3)
outfn <- paste0(dir_exp, "gene.surv.OS_PFI.txt")
write.table(pct.gws, outfn, row.names=TRUE, col.names=NA, quote=F, sep="\t")

#================================================
# Density plot of hazard ratio for cell type marker genes
#================================================

stheme <- theme(plot.title = element_text(size = 8, hjust = 0.5),
                text = element_text(size = 8), 
                strip.text = element_text(size = 6), 
                strip.background = element_rect(fill = "white"))

k <- rownames(g.events.coeff) %in% markers.celltype[,"gene"]
g.events.coeff.ct <- g.events.coeff[k,]
m2ct <- markers.celltype[,cl_name]
names(m2ct) <- make.unique(markers.celltype[,"gene"])
g.events.coeff.ct[,cl_name] <- m2ct[rownames(g.events.coeff.ct)]

ck <- c("cdr.OS.exp_coef", "cdr.PFI.exp_coef", cl_name)
g.events.coeff.ct <- g.events.coeff.ct[,ck]
colnames(g.events.coeff.ct) <- c("OS", "PFI", cl_name)
g.events.coeff.ct[,c("OS", "PFI")] <- log2(g.events.coeff.ct[,c("OS", "PFI")])
g.events.coeff.ct[,"int"] <- 0
g.events.coeff.ct.m <- melt(g.events.coeff.ct, id.vars=cl_name)

ct.coef.sub.l <- lapply(celltypes, function(ct){
                        ct.coef <- ct.coef.l[[ct]]
                        ret.l <- lapply(c("OS", "PFI"), function(e){
                                        data.frame("CellType" = ct, 
                                                   "event" = e, 
                                                   "gene" = rownames(ct.coef), 
                                                   "HR" = ct.coef[,paste0("cdr.", e, ".exp_coef")], 
                                                   "p" = ct.coef[,paste0("cdr.", e, ".p_adj")]) })
                        ret <- do.call(rbind, ret.l)
                        return(ret) })
ct.coef.sub <- do.call(rbind, ct.coef.sub.l)
ct.coef.sub[,"HR"] <- log2(ct.coef.sub[,"HR"])
ct.coef.sub[,"p"] <- -10 * log10(ct.coef.sub[,"p"])

p <- ggplot(ct.coef.sub, aes(x = HR, color = CellType)) + 
    geom_density(alpha = 0) + 
    scale_x_continuous(limits = c(-5, 5)) + 
    geom_vline(xintercept = 0, col = "red") + 
    facet_wrap(~event, ncol=1) + 
    xlab(expression(paste(log[2], " hazard ratio"))) + 
    ylab("Density") + 
    ggtitle("Marker gene hazard ratio\n(TCGA)") + 
    theme_bw() + 
    theme(text = element_text(size = 8), 
          plot.title = element_text(hjust = 0.5, size = 10))
outfn <- paste0(dir_exp, "gene.HR_density.OS_PFI.pdf")
ggsave(outfn, height = 4, width = 3.5)

#================================================
# Plot Prolif marker log2 HRs
#================================================

k <- ct.coef.sub[,"CellType"] == "Prol"
ct.coef.sub.20 <- ct.coef.sub[k,]
ct.coef.sub.20.os <- subset(ct.coef.sub.20, event == "OS")
o <- order(ct.coef.sub.20.os[,"HR"], decreasing=FALSE)
ct.coef.sub.20[,"gene"] <- factor(ct.coef.sub.20[,"gene"], 
                                  levels = ct.coef.sub.20.os[o, "gene"])

p <- ggplot(ct.coef.sub.20, aes(x = gene, y = HR)) + 
    geom_point(size = 0.5) + 
    coord_flip() + 
    scale_y_continuous(limits = c(-6, 6)) + 
    geom_hline(yintercept = 0, col = "red") + 
    facet_wrap(~event, ncol = 2) + 
    ylab(expression(paste(log[2], " hazard ratio"))) + 
    xlab("Prol markers") +
    ggtitle("TCGA") + 
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

fn <- paste0("exp/tcga_hcc/survival/", ct_id, ".cox_ph.ct_med.txt")
ct_med.events.coeff <- read.table(fn, header=TRUE, row.names=1, stringsAsFactors = FALSE, sep="\t")
keys <- c("OS" = "Overall Survival", "PFI" = "Progression Free Interval")
hr.df.l <- lapply(c("OS" = "OS", "PFI" = "PFI"), function(event){
                  ck <- paste0("cdr.", event, c(".exp_coef", ".lower95", ".upper95", ".p", ".p_adj"))
                  hr.df <- ct_med.events.coeff[,ck]
                  # hr.df <- log2(hr.df)
                  colnames(hr.df) <- c("HR", "lower", "upper", "p", "p_adj")
                  hr.df[,"CellType"] <- rownames(hr.df)
                  hr.df[,"Event"] <- event
                  return(hr.df) })
hr.df <- do.call(rbind, hr.df.l)

logcol <- c("HR", "lower", "upper")
hr.df[,logcol] <- log2(hr.df[,logcol])

hr.df[,"sig"] <- "ns"
hr.df[hr.df[,"p_adj"] < 0.05, "sig"] <- "sig"
sig_shape <- c("ns" = 1, "sig" = 16)

o <- order(hr.df.l[["OS"]][,"HR"], decreasing = FALSE)
hr.df[,"CellType"] <- factor(hr.df[,"CellType"], levels = hr.df.l[["OS"]][o,"CellType"])

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
    ggtitle("Cell-type proportional\nHazards Model (TCGA)") +
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5), 
          legend.key.width = unit(0, "in"), 
          legend.position = "bottom")
outfn <- paste0(dir_exp, "ct_med.HR_forest.OS_PFI.pdf")
ggsave(outfn, height = 4, width = 6)

# individual panel plots
e_labs <- names(hr.df.l)
for (event in e_labs){
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
    ggtitle(paste0("Cell-type proportional\nhazards model (TCGA)\n", event)) + 
    th_hr

    out_fn <- paste0(dir_exp, "ct_med.HR_forest.", event, ".pdf")
    ggsave(out_fn, height = 4, width = 3)
}



#=====================================================
# Kaplan-Meier tests and plots
#=====================================================

dir_plot <- dir_exp
dir.create(dir_plot, showWarnings = FALSE, recursive=TRUE)

# get survfit model for plotting
km.events.ct.l <- list()
for (event in events){
    event.col <- event
    time.col <- paste0(event, ".time")

    cell_types <- rownames(tmm.ct.mdp.pt)
    km.ct.l <- list()
    for (ct in cell_types){
        ct.p <- tmm.ct.mdp.pt[ct, samples.pt]
        med.ctp <- quantile(ct.p, probs = c(0.5, 0.5))

        ctp.group <- as.character(as.numeric(ct.p > med.ctp[1]))
        lh.key <- c("0" = "Low", "1" = "High")
        ctp.group <- lh.key[ctp.group]
        ctp.group <- factor(ctp.group, levels = c("Low", "High"))
        names(ctp.group) <- names(ct.p)
        ctp.cid <- samples[names(ctp.group), "Case ID"]

        edf <- data.frame("ctg" = ctp.group, 
                          "time" = cases[ctp.cid, time.col], 
                          "event" = cases[ctp.cid, event.col])
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

event <- "cdr.OS"
pt <- "Overall Survival"
pt <- "OS"
ct <- "Prol"
yticks <- paste0(seq(0, 100, 20), "%")

for (ct in c("Prol")){
    ph <- 3; pw <- 3.5
    surv.md <- km.events.ct.l[[event]][[ct]]
    surv.md$time <- surv.md$time / 365

    p <- autoplot(surv.md, censor.size = 2) + theme_bw() + 
    ggtitle(paste0(pt, " in ", ct, "\n(TCGA)")) + 
    xlab("Years") + ylab("Percent Survival") + 
    scale_x_continuous(breaks = seq(0,10,2)) + 
    guides(color = "none", fill = guide_legend(title=paste0(ct, "\nMedian\nGroup"))) + 
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5), 
          legend.text = element_text(size = 8), 
          legend.title = element_text(size = 8, hjust = 0.5), 
          legend.key.size = unit(.05, units = "npc"))
    p <- p + scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2), labels=yticks)
    
    outfn <- paste0(dir_exp, "curve.", event, ".", ct, ".pdf")
    ggsave(outfn, height=ph, width=pw)
}

event <- "cdr.PFI"
pt <- "Progression Free Interval"
pt <- "PFI"
ct <- "Prol"
for (ct in c("Prol")){
    surv.md <- km.events.ct.l[[event]][[ct]]
    surv.md$time <- surv.md$time / 365
    
    p <- autoplot(surv.md, censor.size = 2) + theme_bw() + 
    ggtitle(paste0(pt, " in ", ct, "\n(TCGA)")) + 
    xlab("Years") + ylab("Percent Progression Free") + 
    scale_x_continuous(breaks = seq(0,10,2)) + 
    guides(color = "none", fill = guide_legend(title=paste0(ct, "\nMedian\nGroup"))) + 
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5), 
          legend.text = element_text(size = 8), 
          legend.title = element_text(size = 8, hjust = 0.5), 
          legend.key.size = unit(.05, units = "npc"))
    p <- p + scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2), labels=yticks)

    outfn <- paste0(dir_exp, "curve.", event, ".", ct, ".pdf")
    ggsave(outfn, height=ph, width=pw)
}

