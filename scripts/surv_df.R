
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

