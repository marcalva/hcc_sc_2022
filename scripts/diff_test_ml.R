
#' Test the difference in proportions between 
#' two groups for multiple levels.
#' 
#' For each sample in x, get the proportion of each z level.
#' For each z, calculate the difference between proportions
#' of samples belonging to the two groups in g.
#' If p is given, run a paired test, where p gives the pair
diff_test_ml <- function(x, g, z, paired = FALSE, p = NULL, 
                         ...){
    x <- factor(x)
    g <- factor(g)
    z <- factor(z)
    if (nlevels(g) != 2)
        stop("g must consist of two levels")
    zprop <- table(x, z)
    zprop <- sweep(zprop, 1, rowSums(zprop), '/')
    tret <- list()
    xu <- levels(x)
    zu <- levels(z)
    gu <- levels(g)
    # get which group each sample in x belongs to
    xg <- sapply(xu, function(i) g[x == i][1])
    names(xg) <- xu
    if (paired){
        p <- factor(p)
        pu <- levels(p)
        xp <- sapply(xu, function(i) p[x == i][1])
        names(xp) <- xu
    }
    for (zi in zu){
        k1 <- names(xg)[xg == gu[1]]
        k2 <- names(xg)[xg == gu[2]]
        z1 <- zprop[k1,zi]
        z2 <- zprop[k2,zi]
        if (paired){
            names(z1) <- xp[names(z1)]
            names(z2) <- xp[names(z2)]
            z2 <- z2[names(z1)]
        }
        wr <- suppressWarnings(wilcox.test(x = z1, y = z2, paired = paired, ...))
        tr <- suppressWarnings(t.test(x = z1, y = z2, paired = paired, ...))
        zret <- data.frame("w_stat" = wr$statistic, 
                           "w_p" = wr$p.value, 
                           "t_stat" = tr$statistic, 
                           "t_p" = tr$p.value, 
                           "m1" = mean(z1), 
                           "m2" = mean(z2))
        colnames(zret)[colnames(zret) == "m1"] <- paste0("mean_", gu[1])
        colnames(zret)[colnames(zret) == "m2"] <- paste0("mean_", gu[2])
        tret[[zi]] <- zret
    }
    tret_df <- do.call(rbind, tret)
    return(tret_df)
}

