
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

