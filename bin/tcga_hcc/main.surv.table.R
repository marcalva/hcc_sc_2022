
# format table of TCGA survival results

setwd("../../")

library(reshape2)

#================================================
# Format output from TCGA survival results
#================================================

ct_id <- "cell_type_main"

dir_exp <- "exp/tcga_hcc/survival/"

fn <- paste0(dir_exp, ct_id, ".cox_ph.ct.txt")
quant.surv <- read.table(fn, header=TRUE, row.names=1, stringsAsFactors = FALSE, sep="\t")

fn <- paste0(dir_exp, ct_id, ".cox_ph.ct_qrt.txt")
qrt.surv <- read.table(fn, header=TRUE, row.names=1, stringsAsFactors = FALSE, sep="\t")

fn <- paste0(dir_exp, ct_id, ".cox_ph.ct_med.txt")
med.surv <- read.table(fn, header=TRUE, row.names=1, stringsAsFactors = FALSE, sep="\t")

fn <- paste0(dir_exp, ct_id, ".cox_ph.stg_ct_med.txt")
quant_stage.surv <- read.table(fn, header=TRUE, row.names=1, stringsAsFactors = FALSE, sep="\t")

surv.l <- list("quartile" = qrt.surv, 
               "median" = med.surv, 
               "Stage + median" = quant_stage.surv)


# surv.l[["continuous"]] <- quant.surv
# surv.l[["continuous"]][,"Model"] <- "Continuous"
surv.l[["quartile"]][,"Model"] <- "Quartile"
surv.l[["median"]][,"Model"] <- "Median"
surv.l[["Stage + median"]][,"Model"] <- "Median (adj stage)"

surv.l <- lapply(surv.l, function(x) {x[,"CellType"] <- rownames(x); return(x)})

col.keep <- c("Model", "CellType", 
              "cdr.OS.n", "cdr.OS.exp_coef", "cdr.OS.lower95", "cdr.OS.upper95", 
              "cdr.OS.p", "cdr.OS.p_adj", 
              "cdr.PFI.n", "cdr.PFI.exp_coef", "cdr.PFI.lower95", "cdr.PFI.upper95", 
              "cdr.PFI.p", "cdr.PFI.p_adj")

surv.l <- lapply(surv.l, function(x) x[,col.keep])

os.keep <- c(1:2, 3:8)
surv.os.l <- lapply(surv.l, function(x) x[,os.keep])
surv.os.l <- lapply(surv.os.l, function(x) {x[,"Event"] <- "OS"; return(x)})
surv.os.l <- lapply(surv.os.l, function(x) {
                    l95 <- as.character(signif(x[,"cdr.OS.lower95"], 3))
                    u95 <- as.character(signif(x[,"cdr.OS.upper95"], 3))
                    lu <- paste0(l95, "-", u95)
                    x[,"95% CI"] <- lu
                    return(x) })
col.keep1 <- c(1:4, 7:10)
surv.os.l <- lapply(surv.os.l, function(x) x[,col.keep1])
col.names0 <- c("Model", "Cell type", "N", "HR", "p", "p_adj", "Event", "95% CI")
surv.os.l <- lapply(surv.os.l, function(x){
                    colnames(x) <- col.names0
                    return(x) })
col_ord <- c("Event", "Model", "Cell type", "N", "HR", "95% CI", "p", "p_adj")
surv.os.l <- lapply(surv.os.l, function(x) x[,col_ord] )


pfi.keep <- c(1:2, 9:14)
surv.pfi.l <- lapply(surv.l, function(x) x[,pfi.keep])
surv.pfi.l <- lapply(surv.pfi.l, function(x) {x[,"Event"] <- "PFI"; return(x)})
surv.pfi.l <- lapply(surv.pfi.l, function(x) {
                     l95 <- as.character(signif(x[,"cdr.PFI.lower95"], 3))
                     u95 <- as.character(signif(x[,"cdr.PFI.upper95"], 3))
                     lu <- paste0(l95, "-", u95)
                     x[,"95% CI"] <- lu
                     return(x) })
col.keep1 <- c(1:4, 7:10)
surv.pfi.l <- lapply(surv.pfi.l, function(x) x[,col.keep1])
col.names0 <- c("Model", "Cell type", "N", "HR", "p", "p_adj", "Event", "95% CI")
surv.pfi.l <- lapply(surv.pfi.l, function(x){
                     colnames(x) <- col.names0
                     return(x) })
col_ord <- c("Event", "Model", "Cell type", "N", "HR", "95% CI", "p", "p_adj")
surv.pfi.l <- lapply(surv.pfi.l, function(x) x[,col_ord] )


surv.l.fmt <- c(surv.os.l, surv.pfi.l)

surv.df <- do.call(rbind, surv.l.fmt)

# round 3 significant digits
csd <- c("HR", "p", "p_adj")
for (i in csd) surv.df[,i] <- signif(surv.df[,i], digits = 3)

# change p > .05 to NS
for (pc in c("p", "p_adj")){
    pk <- rep("NS", nrow(surv.df))
    psig <- surv.df[,pc] < 0.05
    pk[psig] <- as.character(surv.df[psig,pc])
    surv.df[,pc] <- pk
}

out_fn <- paste0(dir_exp, ct_id, ".cox_ph.formatted.txt")
write.table(surv.df, out_fn, row.names=FALSE, col.names=TRUE, 
            quote=FALSE, sep="\t")



