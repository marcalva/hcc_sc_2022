
# Get cell type proportions of individuals
# x is a data frame, e.g. seur@meta.data
# c1 is subject ID column, c2 is cell type column name
get_props <- function(x, c1="orig.ident", c2="seurat_clusters", prop=TRUE){
    datf <- table(x[,c2], x[,c1])
    if (prop){
        datf <- sweep(datf, 2, colSums(datf), "/")
    }
    datf <- as.data.frame.array(datf)
    return(datf)
}

# Collapse expression along cell types
# x is a sparse matrix of counts, the raw gene-barcode matrix
# identities is a vector of cell types
# If counts == TRUE, return the raw counts
collapse = function(x, 
                    identities, 
                    counts=FALSE, 
                    scale_size=1e4, 
                    scale_cells=TRUE,
                    logt=TRUE){
    if (is.null(attr(class(x), "package")) || 
        attr(class(x), "package") != "Matrix"){
        x = Matrix::Matrix(x)
    }
    ngenes = nrow(x)
    if (class(identities) != "factor"){
        identities = factor(identities)
    }
    celltypes = levels(identities)
    ncells = table(identities)[celltypes]
    panel = matrix(0, nrow=ngenes, ncol=length(celltypes))
    rownames(panel) = rownames(x)
    colnames(panel) = celltypes
    for (celltype in celltypes){
        panel[,celltype] = Matrix::rowSums(x[,identities == celltype,drop=FALSE])
    }

    if (counts) return(panel)

    panel <- sweep(panel, 2, colSums(panel), "/")
    
    if (scale_cells){
        panel = scale_size*panel
    }
    
    if (logt){
        panel = log1p(panel)
    }

    return(panel)
}

# seur is a seurat object
# group is the column name of meta.data that has the clusters to merge
# prefix is prefixed to the column name group and contains the merged information
merge_clust = function(seur, 
					   group="seurat_clusters", 
					   prefix="cm.", 
					   r_thrsh=0.99, 
					   assay.type = "RNA")
{
	merge_col = paste0(prefix, group)
	seur@meta.data[,merge_col] = as.character(seur@meta.data[,group])
	clusters <- sort(unique(seur@meta.data[,merge_col]))
	
	cte <- collapse(seur@assays$RNA@counts, identities=seur@meta.data[,merge_col])
	cors <- cor(cte)
	nclusters <- length(clusters)

	i = 1
	j = 2
	while (TRUE){
		group1 = clusters[i]
		group2 = clusters[j]
		R <- cors[group1,group2]
		if (R > r_thrsh){
			merge_name = paste(group1, group2, sep="_")
			to_change = seur@meta.data[,merge_col] == group1 | seur@meta.data[,merge_col] == group2
			seur@meta.data[to_change,merge_col] = merge_name
			clusters <- sort(as.character(unique(seur@meta.data[,merge_col])))

            cte <- collapse(seur@assays$RNA@counts, identities=seur@meta.data[,merge_col])
            cors <- cor(cte)
            nclusters <- length(clusters)

			i = 1 # i = max(1, i - 1)
			j = 2 # j = max(2, j - 1)
		}
		else{
			if (i == nclusters - 1 & j == nclusters) break

			if (j <  nclusters) j = j + 1
			else if (j == nclusters){
				i = i + 1
				j = i + 1
			}
		}
	}
	return(seur)
}

# Decomposition functions

#' Estimate cell type proportions using first PC of expression matrix
#' 
#' @param x A sample by gene bulk expression matrix. Genes should be marker genes
#' @param weighted Boolean. If weighted=TRUE, multiply scaled gene expression by
#'   gene weights
#' @param w Numeric vector. Weights of genes 
#' @return ret List. Attribute \strong{pcs} contains matrix of PCs, where PC1
#'   should be used as estimates for cell type abundances
#'   Attribute \strong{sdev} contains eigenvalues of eigendecomposition of
#'   var-covar matrix. The 1st eigenvalue should explain most of the variance.
#'   Attribute \strong{genes} contains names of genes.
EstimatePCACellTypeProportions <- function(x, weighted=FALSE, w=NULL){
    x <- base::scale(x)
    if (weighted) {
        # Intersect gene names of weights and column names of x
        common.markers <- base::intersect( base::colnames(x), base::names(w) )
        if ( length(common.markers) == 0 ) {
            base::stop(base::paste0("Genes from weights w do not match with column ",
                                    "names of expression matrix x."))
        }
        x <- x[,common.markers]
        w <- w[common.markers]
        wd <- base::diag(w)
        xw <- x %*% wd
        varcov <- t(xw) %*% xw
    }
    else {
        varcov <- t(x) %*% x
    }
    varcov.ed <- base::eigen(varcov)
    rot <- varcov.ed$vectors
    wpcs <- x %*% rot
    sds <- varcov.ed$values
    # x contains PCs, sdev contains eigenvalues of eigendecomposition
    return(list( pcs = wpcs, sdev = sds, genes=colnames(x)))
}


#' Get number of genes to use with weighted PCA
#'
#' @param x Numeric Matrix. A sample by gene expression matrix containing the
#'   marker genes.
#' @param w Numeric Vector. The weights of the genes that correspond to the
#'   columns of x.
#' @param min.gene Numeric. Minimum number of genes to consider as markers.
#' @param max.gene Numeric. Maximum number of genes to consider as markers.
#' @return best.n Numeric. Number of genes to use
GetNumGenesWeighted = function(x, w, min.gene = 25, max.gene = 200){
    max.gene = base::min(max.gene, base::ncol(x))
    ratios = base::lapply(min.gene:max.gene,
                          function(i) {
                              ret = EstimatePCACellTypeProportions(x[,1:i],
                                                                   weighted=TRUE,
                                                                   w=w[1:i])
                              vars = ret$sdev
                              vars[1] / vars[2]
                          })
    best.n = base::which.max(ratios) + min.gene - 1
    return(best.n)
}

GetCTP <- function(bulk,
                   cell_types,
                   markers,
                   ct_col,
                   gene_col,
                   min_gene,
                   max_gene,
                   weighted,
                   w_col,
                   verbose){
  ctp = base::lapply(cell_types, function(ct){
                     # Get marker genes
                     markers_ct = markers[ markers[,ct_col] == ct , , drop=FALSE]
                     ctm = base::make.unique(markers_ct[, gene_col])
                     # Get markers in common between bulk and markers data frame
                     common_markers = base::intersect(ctm, Biobase::featureNames(bulk))
                     if ( base::length(common_markers) == 0 ){
                       base::stop("No marker genes found in bulk expression data")
                     }
                     expr = base::t(Biobase::exprs(bulk)[common_markers,])
                     if ( base::ncol(expr) < min_gene ){
                       base::stop(base::paste0(base::sprintf("For cell type %s, There are less marker genes in ", ct),
                                               base::sprintf("the bulk expression set (%i) than the ", base::ncol(expr)),
                                               base::sprintf("minimum number of genes set (%i) ", min_gene),
                                               "for PCA-based deconvolution\nSet the min_gene parameter to a lower integer."))
                     }
                     if (weighted){
                       # Get gene weights
                       ctw = markers_ct[, w_col]; names(ctw) = ctm; ctw = ctw[common_markers]
                       ng = GetNumGenesWeighted(expr, ctw, min_gene, max_gene) # Number of markers for PCA
                       expr = expr[,1:ng,drop=FALSE]
                       if (verbose){
                         base::cat(base::sprintf("Using %i genes for cell type %s; ", ng, ct))
                       }
                       ret = EstimatePCACellTypeProportions(expr, weighted=TRUE, w=ctw[1:ng])
                     }
                     else{
                       ng = GetNumGenes(expr, min_gene, max_gene)
                       expr = expr[,1:ng,drop=FALSE]
                       if (verbose){
                         base::cat(base::sprintf("Using %i genes for cell type %s; ", ng, ct))
                       }
                       ret = EstimatePCACellTypeProportions(expr)
                     }
                     # Flip the sign of the first PC if negatively correlated with most genes
                     cors = cor(expr, ret$pcs[,1])
                     n_pos = sum(cors[,1] > 0)
                     if (n_pos/base::length(cors[,1]) < (0.5)) ret$pcs[,1] = ret$pcs[,1] * -1
                     if (verbose){
                       cors = cor(expr, ret$pcs[,1]); n_pos = sum(cors[,1] > 0)
                       pct <- as.character(as.integer(100 * round(n_pos/base::length(cors[,1]), 2)))
                       clen <- as.character(base::length(cors[,1]))
                       base::cat(base::paste0(pct, "% of ", clen, " marker genes correlate positively with PC1 for cell type ", ct, "\n"))
                     }
                     return(ret)
                   })
  return(ctp)
}

