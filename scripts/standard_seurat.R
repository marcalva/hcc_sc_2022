
library(Seurat)

#==========================================================
# Normalization
#==========================================================

remove_multiplet <- function(x, lower=FALSE){
	nUMI <- x@meta.data[,"nCount_RNA"]
	upper_bound <- 10^(mean(log10(nUMI)) + 2*sd(log10(nUMI)))
	keep <- nUMI < upper_bound
	if (lower){
		lower_bound <- 10^(mean(log10(nUMI)) - 2*sd(log10(nUMI)))
		keep <- keep & (nUMI > lower_bound)
	}
	cells2keep <- rownames(x@meta.data)[keep]
	x <- subset(x, cells=cells2keep)
	return(x)
}

seurat_filt <- function(x,
                        subset_mt = TRUE,
                        subset_hb = TRUE,
                        mt_thresh = 5,
                        hb_thresh = 5,
                        hb_genes = c("HBA1", "HBA2", "HBB"), 
                        rm_multiplet=TRUE,
                        lower=FALSE){
	mt_genes <- grep(pattern="^mt-", rownames(x), value=TRUE, ignore.case=TRUE)
	# Note this is for human HB genes

	# Subset by MT %
	if (subset_mt){
	    x <- PercentageFeatureSet(x, features=mt_genes, col.name = "pct.mt")
		keep <- (x@meta.data[,"pct.mt"] < mt_thresh)
		cells2keep <- rownames(x@meta.data)[keep]
        if (length(cells2keep) == 0) return(NA)
		x <- subset(x, cells=cells2keep)
	}
    if (ncol(x) < 2) return(x)
	# Subset by hemoglobin
	if (subset_hb){
        x <- PercentageFeatureSet(x, features=hb_genes, col.name = "pct.hb")
		keep <- (x@meta.data[,"pct.hb"] < hb_thresh)
		cells2keep <- rownames(x@meta.data)[keep]
        if (length(cells2keep) == 0) return(NA)
		x <- subset(x, cells=cells2keep)
	}

	if (ncol(x) < 2) return(x)
	if (rm_multiplet) x <- remove_multiplet(x, lower=lower)
    return(x)
}

# Normalize data by either 
# SCTransform or 
# log normalizing, finding variable features, and scaling the data.
seurat_norm <- function(x, 
						vars.to.regress = NULL, 
                        batch_var = NULL, 
						method = "LogNormalize", 
						scale.factor = NULL, 
                        nvar = 3000, 
                        verbose = FALSE){
    if (method == "none") return(x)
	if (method == "sctransform"){
        if (is.null(batch_var)){
            x <- SCTransform(x, 
                             variable.features.n = nvar, 
                             vars.to.regress = vars.to.regress, 
                             verbose = verbose)
        } else {
            x <- SCTransform(x, 
                             variable.features.n = nvar, 
                             vars.to.regress = vars.to.regress, 
                             batch_var = batch_var, 
                             verbose = verbose)
        }
    } else {
        if (is.character(scale.factor)){
            f <- get(scale.factor)
            scale.factor <- f(Matrix::colSums(x@assays$RNA@counts))
        }
	    x <- NormalizeData(x, normalization.method=method, scale.factor=scale.factor)
	    x <- FindVariableFeatures(x, selection.method="vst", nfeatures=nvar)
	    x <- ScaleData(x, vars.to.regress = vars.to.regress)
    }
	return(x)
}

#==========================================================
# Integration
#==========================================================

# Seurat integration for use with Seurat v3.1.2
seurat_integrate <- function(x.list, 
                             n_intg_features=3000, 
                             run.umap=TRUE, 
                             python_lib="/u/local/apps/python/3.6.1-shared/bin/python3.6", 
                             run.cluster=TRUE, 
                             n_dim=30, 
                             resolution=0.8, 
                             verbose=TRUE){
    features <- SelectIntegrationFeatures(object.list = x.list, nfeatures = n_intg_features)
	anchors <- FindIntegrationAnchors(object.list = x.list, anchor.features = features, verbose=verbose)
	integrated <- IntegrateData(anchorset = anchors, verbose=verbose)
	integrated <- RunPCA(object = integrated, verbose = verbose)

    if (run.umap){
    	require(reticulate)
    	reticulate::use_python(python_lib)
    	integrated <- RunUMAP(object = integrated, dims = 1:n_dim)
    }

    if (run.cluster){
        integrated <- FindNeighbors(integrated, dims = 1:n_dim, verbose = verbose)
        integrated <- FindClusters(integrated, resolution=resolution, verbose = verbose)
    }

	return(integrated)
}

#==========================================================
# Cluster
#==========================================================

seurat_cluster <- function(x, 
                           n_dim=30, 
                           resolution = 0.8, 
                           python_lib="/u/local/apps/python/3.6.1-shared/bin/python3.6"){
	require(reticulate)
	reticulate::use_python(python_lib)
	x <- FindNeighbors(x, dims = 1:n_dim, verbose = TRUE)
	x <- FindClusters(x, resolution = resolution, verbose = TRUE)
	x <- RunUMAP(x, dims = 1:n_dim, reduction = "pca", verbose = TRUE)
	return(x)
}

#==========================================================
# Pipelines
#==========================================================

#' @param seur List of seurat objects to integrate
#' @param n_pcs Number of PCs to use
#' @param res Resolution(s) to use for Seurat clustering
log_cca <- function(seur.list, n_pcs = 30, res = c(1, 1.5, 2)){

    seur.list <- lapply(seur.list, function(x){
                        x@active.assay = "RNA"
                        x <- RunPCA(x, npcs = n_pcs, verbose = FALSE)
                        return(x) })


    anchors <- FindIntegrationAnchors(object.list = seur.list, dims = 1:n_pcs)
    seur <- IntegrateData(anchorset = anchors, dims = 1:n_pcs)

    seur <- ScaleData(seur, verbose = FALSE)
    seur <- RunPCA(seur, npcs = n_pcs, verbose = FALSE)

    return(seur)
}

#' @param seur List of seurat objects to integrate
#' @param n_pcs Number of PCs to use
#' @param res Resolution(s) to use for Seurat clustering
sct_cca <- function(seur.list, n_pcs = 30, nfeatures = 3e3){

    for (i in 1:length(seur.list)){
        seur.list[[i]]@active.assay <- "SCT"
        seur.list[[i]] <- RunPCA(seur.list[[i]], npcs = n_pcs, verbose = FALSE)
    }

    features <- SelectIntegrationFeatures(object.list = seur.list, 
                                          nfeatures = nfeatures)

    # This places residuals for features derived from 
    # SelectIntegrationFeatures
    seur.list <- PrepSCTIntegration(seur.list, anchor.features = features)

    anchors <- FindIntegrationAnchors(object.list = seur.list, reduction = "cca", 
                                      dims = 1:30, normalization.method = "SCT", 
                                      anchor.features = features)

    seur <- IntegrateData(anchorset = anchors, dims = 1:30, 
                          normalization.method = "SCT")
    seur <- RunPCA(seur, npcs = n_pcs, verbose = FALSE)

    return(seur)
}

