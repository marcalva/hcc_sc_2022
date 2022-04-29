
pad_plot <- function(g, t=0.025, r=0.025, b=0.025, l=0.025){
    require(grid)
    require(gtable)

    h <- unit(c(t, 1 - (t + b), b), "npc")
    w <- unit(c(l, 1 - (l + r), r), "npc")
    gt <- gtable(widths = w, heights = h)

    gt <- gtable_add_grob(gt, g, t = 2, b = 2, l = 2, r = 2)
    return(gt)
}

gtable_add_tag <- function(g, tag, fs = 8, ff = "bold", x = 0, y = 1, 
                           just = c(0.5, 0.5), t = 1, l = 1, clip = "off", 
                           gname = "tag"){
    gtag <- textGrob(tag, x = unit(x, "npc"), y = unit(y, "npc"), just = just, 
                     gp = gpar(fontsize = fs, fontface = ff))
    g <- gtable_add_grob(g, gtag, t = t, l = l, clip = clip, name = gname)
    return(g)
}

#' Stack ggplot gtables vertically with the panels aligned.
#' the left andn right ends of the panels are all aligned 
#' by padding each gtable to the left and right.
#'
#' @param l list of ggplot gtables returned from ggplotGrob.
#' @param pad_unit units to use to pad to the left and right of the panels.
#' @param width width of stacked gtable.
#' @param heights heights of each of the ggplot gtables in the stacked gtable.
#' @param height_unit units for the heights of the ggplot gtables.
#' @return a gtable, with each of the input gtables as elements.
#' after padding for alignment.
stack_gtable_v <- function(gtl, pad_unit = "points", width = unit(1, "npc"), heights=NULL, 
                           height_unit="npc"){
    require(grid)
    require(gtable)

    # get widths to the left and right of the panel
    len_l_l <- list()
    len_r_l <- list()
    for (i in 1:length(gtl)){
        pan_ix <- which(gtl[[i]]$layout$name == "panel")
        pan_lcol <- min(gtl[[i]]$layout[pan_ix, "l"])
        pan_rcol <- max(gtl[[i]]$layout[pan_ix, "r"])
        lws <- unit.c(gtl[[i]]$widths[1:(pan_lcol-1)])
        rws <- unit.c(gtl[[i]]$widths[(pan_rcol+1):ncol(gtl[[i]])])
        # len_l_l[[i]] <- convertUnit(sum(lws), unitTo=pad_unit)
        # len_r_l[[i]] <- convertUnit(sum(rws), unitTo=pad_unit)
        len_l_l[[i]] <- sum(lws)
        len_r_l[[i]] <- sum(rws)
    }
    # get max length to subtract from
    max_l <- do.call(max, len_l_l)
    max_r <- do.call(max, len_r_l)
    
    for (i in 1:length(gtl)){
        # wl <- convertUnit(max_l - len_l_l[[i]], unitTo=pad_unit)
        # wr <- convertUnit(max_r - len_r_l[[i]], unitTo=pad_unit)
        wl <- max_l - len_l_l[[i]]
        wr <- max_r - len_r_l[[i]]
        gtl[[i]] <- gtable_add_cols(gtl[[i]], widths = wl, pos=0)
        gtl[[i]] <- gtable_add_cols(gtl[[i]], widths = wr, pos=ncol(gtl[[i]]))
    }

    ws <- c(1)
    if (is.null(heights))
        heights = unit(rep(1/length(gtl), length(gtl)), height_unit)
    
    gtmp <- gtable(widths = width, heights = heights)
    for (i in 1:length(gtl))
        gtmp <- gtable_add_grob(gtmp, gtl[[i]], t=i,b=i,l=1,r=1)

    return(gtmp)
}

# old function
stack_panels <- function(grob_l, heights, pad = unit(0, "points")){

    if (length(grob_l) == 1)
        return(grob_l[[1]])
    
    if (length(heights) != length(grob_l))
        stop("length of args must be the same")

    max_widths <- do.call(unit.pmax, lapply(grob_l, function(x) x$widths))
    for (i in 1:length(grob_l)){
        grob_l[[i]]$widths <- max_widths
    }

    gr_st <- grob_l[[1]]

    pan_ix <- which(gr_st$layout$name == "panel")
    pan_ix <- pan_ix[length(pan_ix)]
    pan_t <- gr_st$layout[pan_ix, "t"]
    pan_b <- gr_st$layout[pan_ix, "b"]
    hs <- gr_st$heights[pan_t:pan_b] * heights[1]
    gr_st$heights[pan_t:pan_b] <- hs

    for (i in 2:length(grob_l)){
        gi <- grob_l[[i]]
        # get panels
        add_ixs <- which(gi$layout$name == "panel")
        for (add_ix in add_ixs){
            this_t <- gi$layout[add_ix, "t"]
            this_b <- gi$layout[add_ix, "b"]
            this_nr <- this_b - this_t
            this_hs <- gi$heights[this_t:this_b] * heights[i]
            gr_st <- gtable_add_rows(gr_st, heights = pad, pos = pan_b)
            pan_b <- pan_b + 1
            gr_st <- gtable_add_rows(gr_st, heights = this_hs, pos = pan_b)
            # add grobs in between t and b
            for (j in 1:nrow(gi$layout)){
                j_t <- gi$layout[j, 't']
                j_b <- gi$layout[j, 'b']
                if (j_t >= this_t & j_b <= this_b){
                    gr_st <- gtable_add_grob(gr_st, gi$grobs[[j]], 
                                             t = j_t - this_t + pan_b + 1, 
                                             b = j_b - this_t + pan_b + 1, 
                                             l = gi$layout[j, 'l'], 
                                             r = gi$layout[j, 'r'], 
                                             name = gi$layout[j, 'name'])
                }
            }
            pan_b <- this_b
        }
    }
    return(gr_st)
}
