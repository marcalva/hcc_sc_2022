
setwd("../../")

library(ggplot2)
# source("scripts/go_enr/R/enrich.R")
# source("scripts/go_enr/R/read.R")
# source("scripts/go_enr/R/plot.R")
source("~/scripts/go_enr/R/enrich.R")
source("~/scripts/go_enr/R/plot.R")

#=============================================
#=============================================

infn <- "exp/d_sct_cca_all/merge/markers/seur.CellType.markers.txt"
m <- read.table(infn, header = TRUE, check.names = FALSE, 
                stringsAsFactors = FALSE)

infn <- "exp/d_sct_cca_all/merge/markers/exprd_gn.CellType.txt"
e <- read.table(infn, header = TRUE, check.names = FALSE, 
                stringsAsFactors = FALSE)

go_dat <- readRDS("~/scripts/go_enr/data/human/Reactome/pathway.rds")
td <- go_dat$term_description

cell_types <- unique(m$cluster)
go_l <- lapply(cell_types, function(x){
               message(x)
               genes <- m[m$cluster == x, "gene_name"]
               bg <- e[,"Name"]
               go_df <- fet_enr(genes = genes, 
                                genes2terms = go_dat$genes2terms, 
                                terms2genes = go_dat$terms2genes, 
                                bg_genes = bg)
               terms_i <- intersect(rownames(go_df), rownames(td))
               go_df <- cbind(td[terms_i,], go_df[terms_i,])
               return(go_df)
                })
names(go_l) <- cell_types
dirout <- "exp/d_sct_cca_all/merge/enr/CellType/Reactome/"
dir.create(dirout, showWarnings = FALSE, recursive = TRUE)

for (ct in cell_types){
    p <- plot_enr(go_l[[ct]])
    outfn <- paste0(dirout, "reactome.", ct, ".pdf")
    ggsave(outfn, width = 8, height = 7, dpi = 300)
    outfn <- paste0(dirout, "reactome.", ct, ".txt")
    write.table(go_l[[ct]], outfn, row.names = TRUE, 
                col.names = NA, quote = FALSE, sep = "\t")
}

# merge
go_l.s <- lapply(names(go_l), function(ct) {
                 x <- go_l[[ct]]
                 x <- subset(x, p_adj < 0.05)
                 x[,"CellType"] <- ct
                 return(x) })

merged <- do.call(rbind, go_l.s)

outfn <- paste0(dirout, "reactome.merged.sig.txt")
write.table(merged, outfn, 
            row.names = TRUE, col.names = NA, 
            quote = FALSE, sep = "\t")



