
rasterize_grob <- function(g, w = 3, h = 3, res = 300, bg = "transparent"){
    require(png)
    require(grid)
    require(gtable)

    tfn <- tempfile()
    png(tfn, width = w, height = h, units = "in", res = res, bg = bg)
    pushViewport(viewport())
    grid.draw(g)
    popViewport()
    dev.off()

    pan <- readPNG(tfn)
    file.remove(tfn)

    r <- rasterGrob(pan, width = unit(1, "npc"), height = unit(1, "npc"))
    return(r)
}

#' Rasterize geom points in each panel of a ggtable
#' ggtable is produced by ggpltGrob
#' The resolution of the raster is determined by width height and resolution, 
#' where the number of pixels is w * h * res
#' 
#' @param ggt a gtable produced by ggplotGrob
#' @param w width of temp file
#' @param h heights of temp file
#' @param res ppi resolution of file
#' @param bg bg argument to png function
#' @return the gtable with the panel grob replaced by a rasterGrob
raster_ggpoints <- function(ggt, w = 3, h = 3, res = 300, bg = "transparent"){
    require(png)
    require(grid)
    require(gtable)
    
    ix1 <- which(ggt$layout$name == "panel")
    for (i in ix1){
        ix2 <- grep("points", sapply(ggt$grobs[[i]]$children, function(x) x$name))
        for (j in ix2){
            gtmp <- ggt$grobs[[i]]$children[[j]]
            r <- rasterize_grob(gtmp, w = w, h = h, res = res, bg = bg)
            ggt$grobs[[i]]$children[[j]] <- r
        }
    }
    return(ggt)
}
