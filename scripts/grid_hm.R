
library(reshape2)
library(grid)
library(gtable)
library(ggplot2)

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

gr_lab2 <- function(labs, just = c(1, .5), rot = 0, gp = gpar(), vert = TRUE){
	if (!inherits(labs, "factor")){
		labs <- as.factor(labs)
	}
	labs_l <- levels(labs)
	labs_n <- nlevels(labs)
	
	labs_ls <- sapply(labs_l, function(i) sum(i == labs))
	labs_lp <- labs_ls / length(labs)
	
	y1 <- c(0)
	for (i in 1:(length(labs_lp))){
		y1[i+1] <- sum(labs_lp[1:i])
	}
	ymean <- c()
	for (i in 1:(length(y1)-1)){
		ymean[i] <- mean(c(y1[i], y1[i+1]))
	}
	
	if (vert){
        x <- rep(1, times = labs_n)
		y <- 1-ymean
	}else{
        y <- rep(0, times = labs_n)
        x <- ymean
	}
	
	gtext <- textGrob(label = labs_l, x = x, y = y, just = just, rot = rot, 
		gp = gp, name = "label")
	
	return(gtext)
}

colbar <- function(x, cols=NULL, vert=TRUE){
	if (!inherits(x, "factor")){
		x <- as.factor(x)
	}
	xl <- levels(x)
	xln <- nlevels(x)
	if (is.null(cols)){
		cols <- hcl_pal(xln)
	}
	if(is.null(names(cols))){
		names(cols) <- xl
	}
	w <- 0.5
	h <- 1/length(x)
	y1 <- (1:length(x) - 1) / length(x)
	x1 <- rep(0.5, times = length(x))
	colv <- cols[x]
	if (vert){
		xu <- x1
		yu <- rev(y1) # y axis is top to bottom
		wu <- w
		hu <- h
        just <- c(0.5, 0)
	} else {
		xu <- y1
		yu <- x1
		wu <- h
		hu <- w
        just <- c(0, 0.5)
	}
    gbar <- rectGrob(x = unit(xu, "npc"), y = unit(yu, "npc"), 
                     width = wu, height = hu, just = just, 
                     gp = gpar(col = colv, fill = colv, lineend = "square", 
                               linejoin = "mitre", linemitre=1e6), 
                     name = "colorbar")
	return(gbar)
}

#' @param mar points margin around legend bar key
#' @param bar_w width of legend bar key in points
#' @param bar_h height of legend bar key in points
get_empty_leg <- function(mar = 5.5, bar_w = 12, bar_h = 180){
	# set margin1 to 0.5null to center in gtable
	# margin1, margin2, bar, margin2, label, pad, margin1
	w = c(0.5, mar, bar_w, mar, 0, 0, 0.5)
	wu = c("null", "points", "points", "points", "null", "points", "null")
	# margin1, title, margin2, bar, margin2, margin1
	h = c(0.5, 0, mar, bar_h, mar, 0.5)
	hu = c("null", "null", "points", "points", "points", "null")
	
	g <- gtable(widths = unit(w, wu), heights = unit(h, hu))
	
	return(g)
}

#' @param r vector of numerics giving the range for the color bar key and to 
#' calculate range breaks
#' @param cols color vector
#' @param mar points margin around legend key bar
#' @param ltitle legend title text
gr_leg <- function(r, cols, mar = 5.5, key_dims = c(12, 180), 
		title_leg_text = list("label" = NULL, "gp" = gpar()), 
		lab_leg_text = list("gp" = gpar())){
	
	gt <- get_empty_leg(mar = mar, bar_w = key_dims[1], bar_h = key_dims[2])
	
	# create raster of colors
	cr_cols <- colorRamp(cols)(seq(0,1,.01))
	rgb_cols <- apply(cr_cols, 1, function(i) do.call(rgb, as.list(i/256)))
	rast_cols <- rasterGrob(rgb_cols, width = unit(1, "npc"), height = unit(1, "npc"), 
                            name = "key")
	gt <- gtable_add_grob(gt, rast_cols, t=4,l=3)
	
	# map label breaks to 0-1 scale
	cbreaks <- pretty(r, n=5)
	rmin <- min(r); rmax <- max(r)
	cbreaks <- cbreaks[cbreaks >= rmin & cbreaks <= rmax]
	breaks_loc <- (cbreaks - rmin) / (rmax - rmin)
	
	# create label breaks text grob
	glab <- textGrob(label = as.character(cbreaks), 
		x = unit(0, "npc"), y = unit(breaks_loc, "npc"), 
		just = c(0, 0.5), 
		gp = lab_leg_text[["gp"]], 
        name = "legend_breaks")
	gt <- gtable_add_grob(gt, glab, t=4, l=5, clip=FALSE)
	gt$widths[5] <- convertUnit(grobWidth(glab), unitTo = "points")
	
	# create title
	gtitle <- textGrob(label = title_leg_text[["label"]],
		x = unit(0, "npc"), y = unit(0, "npc"), 
		just = c(0,0), 
		gp = title_leg_text[["gp"]], 
        name = "legend_title")
	gt <- gtable_add_grob(gt, gtitle, t=2, l=3, r=3, clip=FALSE)
	
	# set width of padding to fit title
	wid_bar <- sum(convertUnit(gt$widths[3:5], unitTo = "points"))
	wid_title <- convertUnit(grobWidth(gtitle), unitTo = "points")
	wid_max <- max(wid_bar, wid_title)
	wid_add <- wid_max - wid_bar
	
	gt$widths[6] <- wid_add
	
	return(gt)
}

# returns panel and legend
#' @param datf data frame, 1st column gives x, 2nd y, 3rd value for coloring
#' @param cols color vector to map values to
#' @param val_lims length 2 numeric giving custom min and max limits of values
#' @param ltitle legend title
gr_pan <- function(datf, cols = NULL, val_lims = NULL, key_dims = c(12, 180), 
		title_leg_text = list("label" = NULL, "gp" = gpar()), 
		lab_leg_text = list("gp" = gpar())){
	
	# set to factor to map groups to integers and order
	if (!inherits(datf[,1], "factor")){
		x <- as.factor(datf[,1])
	}else{
		x <- datf[,1]
	}
	
	if (!inherits(datf[,2], "factor")){
		y <- as.factor(datf[,2])
	}else{
		y <- datf[,2]
	}
	
	if (is.null(cols)){
		cols <- hcl(h = seq(0, len=20), c = seq(0, 95, len=20), l = seq(90, 20, len=20))
	}
	
	if (is.null(val_lims)){
		val_lims <- c(min(datf[,3], na.rm=TRUE), max(datf[,3], na.rm=TRUE))
		val_lims <- extendrange(val_lims)
	}
	
	# map value to 0-1 scale for colors
	map_val <- function(v){
		(v - val_lims[1]) / (val_lims[2] - val_lims[1])
	}
	
	xl <- levels(x)
	yl <- levels(y)

	xn <- nlevels(x)
	yn <- nlevels(y)
	
	cr <- colorRamp(cols)
	
	gt <- gTree()
	
	wid <- 1/xn
	hei <- 1/yn
	
	# map x and y levels to x-y axis
    x_map <- seq(0,xn - 1) / xn
	names(x_map) <- xl
    y_map <- seq(0,yn - 1) / yn
	y_map <- rev(y_map) # y goes from top to bottom
	names(y_map) <- yl
	
	# place data values in vectors for rectGrob
    col_v <- sapply(1:nrow(datf), function(i){
                    cr_i <- cr(map_val(datf[i,3]))
                    col_i <- rgb(cr_i[1,1], cr_i[1,2], cr_i[1,3], maxColorValue=256)
        })
    x_v <- x_map[datf[,1]]
    y_v <- y_map[datf[,2]]

    gpan <- rectGrob(x = unit(x_v, "npc"), y = unit(y_v, "npc"), 
                     width = wid, height = hei, just = c(0,0), 
                     gp = gpar(col = col_v, fill = col_v, lineend = "square", 
                               linejoin = "mitre", linemitre=1e6), 
                     name = "panel")
		
	# create legend for color scale
	gleg <- gr_leg(val_lims, rev(cols), mar = 5.5, key_dims = key_dims, 
		title_leg_text = title_leg_text, 
		lab_leg_text = lab_leg_text)
	
	return(list(gpan, gleg))
}

#' @param axis_x_text named list, with "rot", "just", and "gpar", for x axis text parameters
#' @param axis_y_text named list, with "rot", "just", and "gpar", for y axis text parameters
#' @param 
gt_hm <- function(datf, groups_x, groups_y, 
                  val_cols = NULL, group_cols = NULL, 
                  val_lims = NULL, 
                  key_dims = c(12,180), 
                  axis_x_text = list("gp" = gpar(), "rot" = 0, "just" = c(0.5, 0.5)),
                  axis_y_text = list("gp" = gpar(), "rot" = 0, "just" = c(0.5, 0.5)),
                  title_x_text = list("label" = NULL, "gp" = gpar(), "rot"=0), 
                  title_y_text = list("label" = NULL, "gp" = gpar(), "rot"=0), 
                  title_leg_text = list("label" = NULL, "gp" = gpar()), 
                  lab_leg_text = list("gp" = gpar()),
                  x_title = NULL, y_title = NULL, val_title = NULL){

    # margin, lab title,  axis text, axis, panel, margin, legend, margin
    wx <- c(5.5,      12,       12,       12,       1,      12,       12,       12)
    wu <- c("points", "points", "points", "points", "null", "points", "points", "points")
    hx <- c(5.5,      12,       12,       12,       1,      12)
    hu <- c("points", "points", "points", "points", "null", "points")
    gtab <- gtable(widths = unit(x = wx, units = wu), heights = unit(x = hx, units = hu))

    grob_panleg <- gr_pan(datf = datf, cols = val_cols, val_lims = val_lims, 
                          key_dims = key_dims, 
                          title_leg_text = title_leg_text, lab_leg_text = lab_leg_text)
    grob_pan <- grob_panleg[[1]]
    grob_leg <- grob_panleg[[2]]
    gtab <- gtable_add_grob(gtab, grobs = grob_pan, t = 5, l = 5, clip=FALSE, 
                            name = "panel")
    gtab <- gtable_add_grob(gtab, grobs = grob_leg, t = 5, l = 7, name = "legend")
    gtab$widths[7] <- convertUnit(sum(grob_leg$widths), unitTo="cm")

    # add color bar
    cbv <- colbar(groups_y, cols = group_cols, vert=TRUE)
    gtab <- gtable_add_grob(gtab, grobs = cbv, t = 5, l = 4, name = "colorbar")
    cbh <- colbar(groups_x, cols = group_cols, vert=FALSE)
    gtab <- gtable_add_grob(gtab, grobs = cbh, t = 4, l = 5, name = "colorbar")

    # add labels to color bar
    gr_ylab <- gr_lab2(groups_y, vert = TRUE, 
                       gp = axis_y_text[["gp"]], 
                       just = axis_y_text[["just"]], 
                       rot = axis_y_text[["rot"]])
    gr_xlab <- gr_lab2(groups_x, vert = FALSE, 
                       gp = axis_x_text[["gp"]], 
                       just = axis_x_text[["just"]], 
                       rot = axis_x_text[["rot"]])

    t_mult <- 1.25

    gtab <- gtable_add_grob(gtab, grobs = gr_ylab, t = 5, l = 3, name = "yaxlab")	
    gtab$widths[3] <- t_mult * convertUnit(grobWidth(gr_ylab), unitTo="cm")
    gtab <- gtable_add_grob(gtab, grobs = gr_xlab, t = 3, l = 5, name = "xaxlab")
    gtab$heights[3] <- t_mult * convertUnit(grobHeight(gr_xlab), unitTo="cm")

    # add x-axis title
    txg <- textGrob(label = title_x_text[["label"]], rot = title_x_text[["rot"]], 
                    gp = title_x_text[["gp"]], name = "xtitle")
    gtab <- gtable_add_grob(gtab, grobs = txg, t=2, l=5, name = "xtitle")
    gtab$heights[2] <- t_mult * convertUnit(grobHeight(txg), unitTo="cm")

    # add y-axis title
    tyg <- textGrob(label = title_y_text[["label"]], rot = title_y_text[["rot"]], 
                    gp = title_y_text[["gp"]], name = "ytitle")
    gtab <- gtable_add_grob(gtab, grobs = tyg, t=5, l=2, name = "ytitle")
    gtab$widths[2] <- t_mult * convertUnit(grobWidth(tyg), unitTo="cm")

    return(gtab)
}


# test plotting functions
test_grid_hm <- function(){

    set.seed(4)
    n <- 20
    xvals <- runif(n*n,0,1)
    xgrps <- sample(c("a", "b", "c"), size = n, replace=TRUE)
    mat <- matrix(xvals, ncol=n, nrow=n)
    mat <- mat %*% t(mat)
    mat <- mat / max(mat)

    # order by group
    o <- order(xgrps)
    xgrps <- xgrps[o]
    mat <- mat[o,o]
    names(xgrps) <- as.character(seq(1:n))

    datfm <- melt(mat)
    datfm[,1] <- as.character(datfm[,1])
    datfm[,2] <- as.character(datfm[,2])
    datfm[,1] <- factor(datfm[,1], levels = as.character(seq(1,n)))
    datfm[,2] <- factor(datfm[,2], levels = as.character(seq(1,n)))
    datfm[,"Group1"] <- xgrps[datfm[,"Var1"]]
    datfm[,"Group2"] <- xgrps[datfm[,"Var2"]]

    # red color palette
    reds <- hcl(h = seq(0, len=n), c = seq(0, 95, len=n), l = seq(90, 20, len=n))

    # text parameters
    axis_x_text <- list("gp" = gpar(), "rot" = 0, "just" = c(0.5, 0.5))
    axis_y_text <- list("gp" = gpar(), "rot" = 0, "just" = c(0.5, 0.5))
    title_leg_text <- list("label" = bquote(italic(R)), "gp" = gpar())
    lab_leg_text <- list("gp" = gpar())
    key_dims <- c(12,120)

    # plot
    p <- gt_hm(datfm, xgrps, xgrps, 
               key_dims = key_dims, 
               axis_x_text = axis_x_text, 
               axis_y_text = axis_x_text, 
               title_leg_text = title_leg_text, 
               lab_leg_text = lab_leg_text
               )
    grid.draw(p)

}
