
hcl_pal <- function(n, chr = c(80, 80), lum = c(50,90), offset = 0, 
                    rand = TRUE, seedn = 1){
    chrval <- seq(from = chr[1], to = chr[2], len = n)
    lumval <- seq(from = lum[1], to = lum[2], len = n)
    hval <- offset + seq(0, 360, len = n + 1)
    hval <- hval[1:n]
    if (rand){
        set.seed(seedn, kind = 'Mersenne-Twister')
        chrval <- sample(chrval)
        lumval <- sample(lumval)
        hval <- sample(hval)
    }
    pal <- hcl(h = hval, c = chrval, l = lumval)
    return(pal)
}

hcl_pal_step <- function(n, chr = c(60, 100), lum = c(60,80), 
                        offset = 0, rand = TRUE, seedn = 1){
    chrval <- rep(chr, length.out = n)
    lumval <- rep(lum, length.out = n)
    hval <- offset + seq(0, 360, len = n + 1)
    hval <- hval[1:n]
    o <- 1:n
    if (rand){
        set.seed(seedn, kind = 'Mersenne-Twister')
        o <- sample(1:n)
    }
    pal <- hcl(h = hval[o], c = chrval[o], l = lumval[o])
    return(pal)
}

