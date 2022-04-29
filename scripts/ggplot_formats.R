
get_astk <- function(x, sig_levels = c(5e-2, 5e-3, 5e-4)){
    if (is.na(x)) return(NULL)
    for (i in length(sig_levels):1){
         if (x < sig_levels[i])
             return(strrep("*", times = i))
    }
    return(NULL)
}

frmt_sct <- function(x, dig=3){
    x <- as.numeric(x)
    pow <- floor(log10(x))
    xleft <- x * 10^-pow
    xleft <- signif(xleft, dig)
    y <- paste0(xleft, "~x~10^",pow)
    return(y)
}

log10_rev_breaks <- function(x, n=5){

    rng <- range(x, na.rm = TRUE) 
    lxmin <- floor(log10(rng[1])) + log10(2)
    lxmax <- ceiling(log10(rng[2])) - log10(2)

    lby <- floor((lxmax - lxmin)/n) + 1

    breaks <- rev(10^seq(lxmax, lxmin, by=-lby))
    return(breaks)
}

format_power10 <- function(x){
    x <- signif(x, digits = 2)
    sapply(x, function(y){
           pow_num <- floor(log10(y))
           base_num <- y * 10^-pow_num
           ret <- bquote(.(base_num) %*% 10^.(pow_num))
           as.expression(ret) })
}

log10_rev_trans <- function(x){
    trans <- function(x) -log(x, 10)
    inv <- function(x) 10^(-x)
    trans_new("log10_rev", trans, inv, breaks = log10_rev_breaks,
              format = format_power10, domain = c(1e-100, Inf))
}

