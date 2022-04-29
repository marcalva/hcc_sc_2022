
# export meta data and umap

setwd("../../")

library(Seurat)

fn <- "data/processed/sharma_aiz/liver.int_rand.rds"
seur <- readRDS(fn)

um <- seur@reductions$umap@cell.embeddings
md <- seur@meta.data

mdum <- cbind(md, um[rownames(md),])

out_fn <- "data/processed/sharma_aiz/liver.int_rand.md_umap.txt.gz"
gzf <- gzfile(out_fn, open = "w")
write.table(mdum, gzf, row.names=TRUE, col.names=NA, 
            quote=FALSE, sep='\t')
close(gzf)

