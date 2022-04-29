
# test difference between G2M and S scores per main cell-type

setwd("../../")

fn <- "data/processed/sharma_aiz/liver.int_rand.md_umap.txt.gz"
udf <- read.table(fn, header=TRUE, row.names=1, sep='\t')

w_t_test <- function(x, y){
    tr <- t.test(x, y)
    wr <- wilcox.test(x, y)
    r <- c(tr$statistic, tr$parameter, tr$p.value, wr$statistic, wr$p.value)
    names(r) <- c("Tstat", "Tdof", "Tp", "Wstat", "Wp")
    return(r)
}

cl_ids <- c("cell_type_main", "cell_type_fine")
for (cl_id in cl_ids){
    cell_types <- sort(unique(udf[,cl_id]))

    scores <- c("G2M", "S")

    wt_r_l <- list()

    for (score in scores){
        tmp_l <- list()
        for (cell_type in cell_types){
            col_id <- paste0(score, ".Score")
            xk <- which(udf[,cl_id] == cell_type)
            yk <- which(udf[,cl_id] != cell_type)
            x <- udf[xk,col_id]
            y <- udf[yk,col_id]
            r <- w_t_test(x, y)
            tmp_l[[cell_type]] <- r
        }
        tmp_datf <- as.data.frame(do.call(rbind, tmp_l))
        tmp_datf[,"CellType"] <- cell_types
        tmp_datf[,"Score"] <- rep(score, length.out = nrow(tmp_datf))
        tmp_datf[,"Tpadj"] <- p.adjust(tmp_datf[,"Tp"], method="fdr")
        tmp_datf[,"Wpadj"] <- p.adjust(tmp_datf[,"Wp"], method="fdr")
        wt_r_l[[score]] <- tmp_datf

    }
    wt_r <- do.call(rbind, wt_r_l)

    dir_exp <- "exp/sharma_aiz/cc_diff/"
    dir.create(dir_exp, showWarnings=FALSE, recursive=TRUE)

    out_fn <- paste0(dir_exp, cl_id, ".cc_diff.txt")
    write.table(wt_r, out_fn, row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
}
