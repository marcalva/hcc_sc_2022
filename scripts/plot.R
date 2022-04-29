
# Color scale
gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}

#' box plot
plot_box <- function(datf, x, y){
	require(reshape2)

    # df <- reshape2::melt(datf)
	p <- ggplot(datf, aes_string(x = x, y = y)) +
	geom_boxplot() + theme_classic() + 
    theme(legend.position = "none")
	return(p)
}


# x is a seurat object
boxplot_gene <- function(x, gene, ct_col="RNA_snn_res.0.8", assay = NULL){
	require(reshape2)
    expr <- GetAssayData(x, slot="data", assay=assay)
    feat <- expr[gene,]
    ct <- x@meta.data[,ct_col]
    datf <- data.frame("Gene" = feat, "CellType" = ct)

	df <- reshape2::melt(datf)
	p <- ggplot(datf, aes(x=CellType, y=Gene)) +
	geom_boxplot() + theme_classic() + 
    ggtitle(gene) + 
    ylab("Expression") + 
    theme(legend.position = "none")
	return(p)
}

# ggplot scatter plot
plot_data <- function(datf, x, y, colors = NULL, label_col = NULL, 
                      trans = "identity",
                      legend_title = waiver(),
                      alpha = 1, 
                      shape = 21, 
                      stroke = 0.3, 
                      label = TRUE, 
                      repel = TRUE, 
                      label_size=8, 
                      size=2,
                      label_geom = "text", 
                      rand=TRUE){
    require(ggplot2)
    require(ggrepel)
    set.seed(1)
    if (rand) datf <- datf[sample(rownames(datf)),]
    # Get median label points for each cluster

    # Add colors if labels are provided
    if (!is.null(label_col)){
        datf[,label_col] <- factor(datf[,label_col])
        p <- ggplot(datf, aes_string(x = x, y = y))
        labs_u <- sort(unique(datf[,label_col]))
        meds <- lapply(labs_u, function(lu) {
                       luk <- datf[,label_col] == lu
                       datfk <- datf[luk, c(x, y)]
                       meds <- apply(datfk, 2, median)
                       return(meds)
                      })
        names(meds) <- labs_u
        meds <- as.data.frame(do.call(rbind, meds))
        colnames(meds) <- c("x", "y")
        rownames(meds) <- labs_u

        if (is.null(colors)){
            nb.cols <- length(labs_u)
            mycolors <- gg_color_hue(nb.cols)
            mycolors <- sample(mycolors)
        } else {
            mycolors <- colors
        }

        p <- p + geom_point(aes_(color = as.name(label_col)), stroke = stroke, 
                            alpha = alpha, shape = shape, size = size)
        p <- p + scale_color_manual(values = mycolors, name = legend_title)

        if (label){
            if (repel){
                p <- p + geom_text_repel(data = meds, aes(x = x, y = y, label = rownames(meds)))
            } else{
                if (label_geom == "text"){
                    p <- p + geom_text(data = meds, aes(x=x, y=y, label = rownames(meds)))
                } else{
                    p <- p + geom_label(data = meds, aes(x=x, y=y, label = rownames(meds)))
                }
            }
        }
    } else {
        p <- ggplot(datf, aes_string(x = x, y = y))
        p <- p + geom_point(alpha=alpha, shape = shape, size = size)
    }
    p <- p + theme_classic() + 
    theme(axis.text=element_blank(), axis.ticks=element_blank()) + 
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 1.5)))

    return(p)
}

plot_data_c <- function(datf, x, y, color_col = NULL, 
                        low = "lightgrey", high = "firebrick3",
                        trans = "identity",
                        legend_title = waiver(),
                        alpha = 1, 
                        shape = 21, 
                        stroke = 0.3, 
                        label_size=8, 
                        size=2,
                        label_geom = "text", 
                        rand=TRUE){
    require(ggplot2)
    require(ggrepel)
    set.seed(1)
    if (rand) datf <- datf[sample(rownames(datf)),]
    # Get median label points for each cluster

    # Add colors if labels are provided
    p <- ggplot(datf, aes_string(x = x, y = y))
    p <- p + geom_point(aes_(color = as.name(color_col)), stroke = stroke, 
                        alpha = alpha, shape = shape, size = size)
    p <- p + scale_colour_gradient(low = low, 
                                   high = high, 
                                   name = legend_title)
    p <- p + theme_classic() + 
    theme(axis.text=element_blank(), axis.ticks=element_blank())

    return(p)
}

# Plot continuous values in a UMAP from meta.data column
# x is a seurat object
# returns a ggplot object
plot_umap_c <- function(x, 
                        colname="percent.mt", 
                        legend_title="MT%", 
                        low = "lightgrey", high = "firebrick3", 
                        trans = "identity", 
                        alpha=1, 
                        size = 1,
                        shape = 16, 
                        rand=TRUE){
    require(ggplot2)
    require(RColorBrewer)
    datf <- data.frame(feat=x@meta.data[,colname],
                       UMAP1=x@reductions$umap@cell.embeddings[,"UMAP_1"],
                       UMAP2=x@reductions$umap@cell.embeddings[,"UMAP_2"])
    if (rand) datf <- datf[sample(rownames(datf)),]
    p <- ggplot(datf, aes(x=UMAP1, y=UMAP2, color=feat)) + 
    geom_point(size = size, alpha=alpha, shape = shape) + 
    theme_classic() +
    theme(axis.text=element_blank(), axis.ticks=element_blank()) +
    scale_color_gradient(low=low, 
                         high=high, 
                         trans = trans, 
                         name=legend_title)
    return(p)
}

# Plot discrete values in a UMAP from meta.data column
# x is a seurat object
# returns a ggplot object
plot_umap_d <- function(x, 
                        colname="orig.ident", 
                        legend_title="orig.ident",
                        colors=NULL, 
                        legend=TRUE, 
                        alpha=1, 
                        size=2,
                        shape = 16, 
                        rand=TRUE){
    require(ggplot2)
    require(RColorBrewer)
    datf <- data.frame(feat=x@meta.data[,colname], 
                       UMAP1=x@reductions$umap@cell.embeddings[,"UMAP_1"],
                       UMAP2=x@reductions$umap@cell.embeddings[,"UMAP_2"])

    if (rand) datf <- datf[sample(rownames(datf)),]
    clusters <- unique(datf$feat)
    # Define colors
    if (is.null(colors)){
        nb.cols <- length(clusters)
        mycolors <- gg_color_hue(nb.cols)
    } else {
        mycolors <- colors
    }
    p <- ggplot(datf, aes(x=UMAP1, y=UMAP2)) + 
    geom_point(aes(color=factor(feat)), alpha=alpha, shape = shape, size=size) + 
    theme_classic() + scale_color_manual(values = mycolors, name = legend_title) + 
    theme(axis.text=element_blank(), axis.ticks=element_blank(), text=element_text(size=20)) +
    guides(colour = guide_legend(override.aes = list(alpha = 1, size = 1.5)))
    if (!legend){
        p <- p + theme(legend.position = "none")
    }
    return(p)
}

# Plot gene expression values in a UMAP
# x is a seurat object
# returns a ggplot object
plot_umap_gene_2c <- function(x, 
                           gene, 
                           legend_title="UMI", 
                           low = "lightgrey", high = "firebrick3", 
                           logt = FALSE,
                           assay=NULL, 
                           alpha=1, 
                           size=3,
                           shape = 16, 
                           rand=TRUE){
    require(ggplot2)
    require(RColorBrewer)

    expr <- GetAssayData(x, slot="data", assay=assay)

    datf <- data.frame(feat=expr[gene,],
                       UMAP1=x@reductions$umap@cell.embeddings[,"UMAP_1"],
                       UMAP2=x@reductions$umap@cell.embeddings[,"UMAP_2"])
    if (logt){
        datf[,"feat"] <- log1p(datf[,"feat"])
    }
    if (rand) datf <- datf[sample(rownames(datf)),]
    p <- ggplot(datf, aes(x=UMAP1, y=UMAP2, color=feat)) + 
    geom_point(alpha=alpha, shape = shape, size = size) + 
    theme_classic() + ggtitle(gene) + 
    theme(axis.text=element_blank(), axis.ticks=element_blank()) +
    scale_color_gradient(low=low, high=high, name=legend_title)
    return(p)
}

# Plot gene expression values in a UMAP
# x is a seurat object
# returns a ggplot object
plot_umap_gene <- function(x, 
                           gene, 
                           legend_title="UMI", 
                           palette="GnBu", 
                           logt = FALSE,
                           assay=NULL, 
                           alpha=1, 
                           shape = 16, 
                           rand=TRUE){
    require(ggplot2)
    require(RColorBrewer)

    expr <- GetAssayData(x, slot="data", assay=assay)

    datf <- data.frame(feat=expr[gene,],
                       UMAP1=x@reductions$umap@cell.embeddings[,"UMAP_1"],
                       UMAP2=x@reductions$umap@cell.embeddings[,"UMAP_2"])
    if (logt){
        datf[,"feat"] <- log1p(datf[,"feat"])
    }
    if (rand) datf <- datf[sample(rownames(datf)),]
    p <- ggplot(datf, aes(x=UMAP1, y=UMAP2, color=feat)) + 
    geom_point(alpha=alpha, shape = shape) + 
    theme_classic() + ggtitle(gene) + 
    theme(axis.text=element_blank(), axis.ticks=element_blank()) +
    scale_colour_distiller(palette = palette, name=legend_title)
    return(p)
}

#' @param datf data frame to get median points from
#' @param x character vector of column names to calculate median for
#' @param groupby the column name to group the rows of the data frame by
get_med_points <- function(datf, x, groupby){
    groups <- sort(unique(datf[,groupby]))
    gs.l <- lapply(groups, function(gr){
                 k <- datf[,groupby] == gr
                 datf.s <- datf[k, x, drop=FALSE]
                 r <- apply(datf.s, 2, median)
                 return(r) })
    names(gs.l) <- groups
    gs <- as.data.frame(do.call(rbind, gs.l))
    colnames(gs) <- x
    rownames(gs) <- groups
    return(gs)
}

# Plot cell types in UMAP with labels
# x is a Seurat object
# returns a ggplot object
plot_umap_labels <- function(x, 
                             colname=NULL, 
                             legend_title="Cell Type",
                             colors=NULL, 
                             alpha=1, 
                             shape = 16, 
                             label=TRUE, 
                             repel = TRUE, 
                             label_size=3, 
                             ax_title_size=10, 
                             size=2,
                             text_size = 12,
                             label_geom = "text", 
                             font_face = "plain", 
                             rand=TRUE){
    library(ggplot2)
    library(RColorBrewer)
    library(ggrepel)

    if (is.null(colname)) celltypes <- as.character(x@active.ident)
    else celltypes <- as.character(x@meta.data[,colname])
    datf <- data.frame(Cluster=celltypes, 
                       UMAP1=x@reductions$umap@cell.embeddings[,"UMAP_1"],
                       UMAP2=x@reductions$umap@cell.embeddings[,"UMAP_2"])
    if (rand) datf <- datf[sample(rownames(datf)),]

    # Get median UMAP points for each cluster
    meds <- get_med_points(datf, c("UMAP1", "UMAP2"), "Cluster")
    clusters <- rownames(meds)

    # Define colors
    if (is.null(colors)){
        nb.cols <- length(clusters)
        mycolors <- gg_color_hue(nb.cols)
    } else {
        mycolors <- colors
    }

    p <- ggplot(datf, aes(x=UMAP1, y=UMAP2)) + 
        geom_point(aes(color=factor(Cluster)), alpha=alpha, shape = shape, size = size) + 
        theme_classic() + 
        scale_color_manual(values = mycolors, name = legend_title) + 
        theme(text = element_text(size = text_size, face = font_face), 
              axis.text=element_blank(), 
              axis.ticks=element_blank(), 
              axis.title=element_text(size=ax_title_size, face = font_face)) + 
        guides(colour = guide_legend(override.aes = list(alpha = 1, size = 1.5)))

    if (label){
        if (repel){
            p <- p + geom_text_repel(data=meds, aes(x=UMAP1, y=UMAP2, label=rownames(meds)), 
                                     size=label_size, 
                                     fontface = font_face)
        } else{

            if (label_geom == "text"){
                p <- p + geom_text(data= meds, aes(x=UMAP1, y=UMAP2, label = rownames(meds)))
            } else{
                p <- p + geom_label(data= meds, aes(x=UMAP1, y=UMAP2, label = rownames(meds)))
            }
        }
    }

    return(p)
}

plot_pc_labels <- function(x, 
                           ct_col="SCT_snn_res.0.8", 
                           legend_title="Cell Type", 
                           shape = 16){
    df <- data.frame(Cluster=x@meta.data[,ct_col], 
                     PC1=x@reductions$pca@cell.embeddings[,"PC_1"],
                     PC2=x@reductions$pca@cell.embeddings[,"PC_2"])
    # Get median PC points for each cluster
    clusters <- sort(unique(df$Cluster))
    meds <- lapply(clusters, function(ct) {
                   pcs <- df[df$Cluster == ct, c("PC1", "PC2")]
                   meds <- apply(pcs, 2, median)
                   return(meds)
                     })
    names(meds) <- clusters
    p <- ggplot(df, aes(x=PC1, y=PC2, color=Cluster)) + 
    geom_point(shape = shape) + theme_classic() + scale_color_discrete(name = legend_title) + 
    theme(axis.text=element_blank(), axis.ticks=element_blank())
    for (ct in clusters){
        p <- p + annotate(geom="label", x=meds[[ct]][1], y=meds[[ct]][2], label=as.character(ct))
    }
    return(p)
}

plot_pc_cont <- function(x, md_col="percent.mt", legend_title="MT%"){
    df <- data.frame(feat=x@meta.data[,md_col],
                     PC1=x@reductions$pca@cell.embeddings[,"PC_1"],
                     PC2=x@reductions$pca@cell.embeddings[,"PC_2"])
    p <- ggplot(df, aes(x=PC1, y=PC2, color=feat)) + 
    geom_point() + theme_classic() +
    theme(axis.text=element_blank(), axis.ticks=element_blank()) +
    scale_colour_distiller(palette = "Spectral", name=legend_title)
    return(p)
}

# Boxplot of MT%
boxplot_mt <- function(x, mt_col="percent.mt", ct_col="SCT_snn_res.0.8"){
    require(reshape2)
    df <- data.frame(Mito=x@meta.data[,mt_col], Cluster=x@meta.data[,ct_col])
    df <- reshape2::melt(df)
    p <- ggplot(df, aes(x=Cluster, y=value, fill=Cluster)) +
    geom_boxplot() + theme_classic() + 
    ylab("MT%")
    return(p)
}

# Gene expression heatmap
expr_heatmap <- function(x, de, ct_col, n_genes=10, n_cells=50, colors=NULL, breaks=NA){
    cell_types <- sort(unique(de$cluster))
    topn <- lapply(cell_types, function(i) de[de$cluster == i,"gene"][1:n_genes])
    topn <- unlist(topn)
    barcodes <- lapply(cell_types, function(i) {
                       sample(colnames(x)[x@meta.data[,ct_col] == i], size=n_cells)})
    barcodes <- unlist(barcodes)
    expr <- as.matrix(x@assays$SCT@data[topn, barcodes])
    expr <- t(scale(t(expr)))

    anno_col <- lapply(cell_types, function(i) (rep(i, n_cells)))
    anno_col <- as.data.frame(unlist(anno_col))
    rownames(anno_col) <- colnames(expr)
    colnames(anno_col) <- "Cell Type"

    if (is.null(colors)){
        p <- pheatmap(expr, annotation_col=anno_col, breaks=breaks, 
                      cluster_cols=FALSE, cluster_rows=FALSE,
                      show_rownames=FALSE, show_colnames=FALSE,
                      annotation_names_col = FALSE)
    } else{
        p <- pheatmap(expr, annotation_col=anno_col, breaks=breaks, 
                      cluster_cols=FALSE, cluster_rows=FALSE,
                      show_rownames=FALSE, show_colnames=FALSE,
                      annotation_names_col = FALSE,
                      annotation_colors=colors)
    }
    return(p)
}

expr_heatmap_genes <- function(x, genes, ct_col, n_cells=50, colors=NULL){
    cell_types <- sort(unique(x@meta.data[,ct_col]))
    barcodes <- lapply(cell_types, function(i) {
                       sample(colnames(x)[x@meta.data[,ct_col] == i], size=n_cells)})
    barcodes <- unlist(barcodes)
    expr <- as.matrix(x@assays$SCT@data[genes, barcodes])
    expr <- t(scale(t(expr)))

    anno_col <- lapply(cell_types, function(i) (rep(i, n_cells)))
    anno_col <- as.data.frame(unlist(anno_col))
    rownames(anno_col) <- colnames(expr)
    colnames(anno_col) <- "Cell Type"

    if (is.null(colors)){
        p <- pheatmap(expr, annotation_col=anno_col,
                      cluster_cols=FALSE, cluster_rows=FALSE,
                      show_rownames=FALSE, show_colnames=FALSE,
                      annotation_names_col = FALSE)
    } else{
        p <- pheatmap(expr, annotation_col=anno_col,
                      cluster_cols=FALSE, cluster_rows=FALSE,
                      show_rownames=FALSE, show_colnames=FALSE,
                      annotation_names_col = FALSE,
                      annotation_colors=colors)
    }
    return(p)
}

