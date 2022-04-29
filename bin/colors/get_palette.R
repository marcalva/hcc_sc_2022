
setwd("../../")

dirc <- "data/ref/colors/"
dir.create(dirc, recursive=TRUE, showWarnings=FALSE)

# qualitative palettes

for (i in c(10,20,30,40,50)){
    pal <- hsv(seq(0,1,len=i+1), .95,.95)
    pal <- pal[1:i]
    write.table(pal, paste0(dirc, "hsv_full.pal.", as.character(i), ".csv"), 
                row.names=FALSE, col.names=FALSE, quote=TRUE, sep=',')
}

for (i in c(10,20,30,40,50)){
    pal <- hcl(seq(0,360,len=i+1), c = 80, l = 75)
    pal <- pal[1:i]
    write.table(pal, paste0(dirc, "hcl_even.pal.", as.character(i), ".csv"), 
                row.names=FALSE, col.names=FALSE, quote=TRUE, sep=',')
}

for (i in c(10,20,30,40,50)){
    pal <- hcl(seq(0,360,len=i+1), c = 70, l = c(50,80,60,90,70))
    pal <- pal[1:i]
    write.table(pal, paste0(dirc, "hcl_c70_vl.pal.", as.character(i), ".csv"), 
                row.names=FALSE, col.names=FALSE, quote=TRUE, sep=',')
}

# sequential color palettes
n <- 20
sv <- seq(1,.05,length.out=20)
vv <- seq(.5,1,length.out=20)

# reds <- hsv(h = 0, s = sv, v = vv)
reds <- hcl(h = seq(0, len=n), c = seq(0, 95, len=n), l = seq(90, 20, len=n))
datf <- reds
write.table(reds, paste0(dirc, "red_colrs.csv"), 
            row.names=FALSE, col.names=FALSE, quote=TRUE, sep=',')

# greens <- hsv(h = 0.33, s = sv, v = vv)
greens <- hcl(h = seq(120, len=n), c = seq(0, 95, len=n), l = seq(90, 20, len=n))
write.table(greens, paste0(dirc, "green_colrs.csv"), 
            row.names=FALSE, col.names=FALSE, quote=TRUE, sep=',')

# blues <- hsv(h = 0.6, s = sv, v = vv)
blues <- hcl(h = seq(235, len=n), c = seq(0, 95, len=n), l = seq(90, 20, len=n))
write.table(blues, paste0(dirc, "blue_colrs.csv"), 
            row.names=FALSE, col.names=FALSE, quote=TRUE, sep=',')


# diverging palette
n <- 19
bl <- hcl(h = seq(235, len=n), c = seq(10, 95, len=n), l = seq(100, 30, len=n))
rd <- hcl(h = seq(0, len=n), c = seq(10, 95, len=n), l = seq(100, 30, len=n))
mid <- hcl(h = 115, c = 5, l = 100)
rd_bu <- c(rev(rd), mid, bl)

sv <- c(seq(1,.05,length.out=20), 0, seq(0.05,1,length.out=20))
vv <- c(seq(.5,1,length.out=20), 1, seq(1,.5,length.out=20))
hv <- c(rep(0, 20), 0, rep(0.66, 20))
rd_bu_hsv <- hsv(h = hv, s = sv, v = vv)

write.table(rd_bu, paste0(dirc, "rd_bu_div.csv"), 
            row.names=FALSE, col.names=FALSE, quote=TRUE, sep=',')


# source color
cols <- hcl(h=c(15,190,280), c = 125, l = 70)
datf <- data.frame("Source" = c("nash_hcc", "Aizarani", "Sharma"), cols)
write.table(datf, paste0(dirc, "source_colrs.csv"), 
            row.names=FALSE, col.names=FALSE, quote=TRUE, sep=',')

# tumor color
cols <- hcl(h=c(0, 250), c = c(100,40), l = c(55,55))
datf <- data.frame("Tumor" = c("Tumor", "NonTumor"), cols)
write.table(datf, paste0(dirc, "tumor_colrs.csv"), 
            row.names=FALSE, col.names=FALSE, quote=TRUE, sep=',')

# viral color
cols <- hcl(h=c(60, 270), c = c(100,100), l = c(75,65))
datf <- data.frame("Viral" = c("Viral", "NonViral"), cols)
write.table(datf, paste0(dirc, "viral_colrs.csv"), 
            row.names=FALSE, col.names=FALSE, quote=TRUE, sep=',')

# patient colors
cols <- hcl(h=c(170, 40, 270), c = c(100,100,100), l = c(55,55,55))
datf <- data.frame("Patient" = c("1", "2", "3"), cols)
write.table(datf, paste0(dirc, "patient_colrs.csv"), 
            row.names=FALSE, col.names=FALSE, quote=TRUE, sep=',')

# sample colors
cols <- hcl(h = c(0,0,0,250,250,250), c = c(30,60,90), l = c(rep(c(85,55,20),2)))
datf <- data.frame("Sample" = c(paste("Sample", as.character(1:3), "tumor"), 
                     paste("Sample", as.character(1:3), "nontumor")), 
                   cols)
write.table(datf, paste0(dirc, "sample_colrs.csv"), 
            row.names=FALSE, col.names=FALSE, quote=TRUE, sep=',')


