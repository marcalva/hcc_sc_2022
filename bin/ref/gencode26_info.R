
setwd("../../")

# Get map of ensembl to gene names

gtf_path <- "data/ref/gencode26/gencode.v26.annotation.gtf"

gtf <- read.table(gtf_path,
                  sep="\t",
                  header = FALSE,
                  stringsAsFactors=FALSE)
keep <- gtf[,3] == "gene"
gtf <- gtf[keep,,drop=FALSE]
gene_ids <- sapply(gtf[,9], function(g){
                   s <- strsplit(g, " ")[[1]]
                   for (i in 1:length(s)){
                       if (s[i] == "gene_id"){
                           sr = sub(";", "", s[i + 1])
                           return(sr)
                       }
                   }
                   return(NA)
                  })


gene_names <- sapply(gtf[,9], function(g){
                   s <- strsplit(g, " ")[[1]]
                   for (i in 1:length(s)){
                       if (s[i] == "gene_name"){
                           sr = sub(";", "", s[i + 1])
                           return(sr)
                       }
                   }
                   return(NA)
                  })


gene_types <- sapply(gtf[,9], function(g){
                   s <- strsplit(g, " ")[[1]]
                   for (i in 1:length(s)){
                       if (s[i] == "gene_type"){
                           sr = sub(";", "", s[i + 1])
                           return(sr)
                       }
                   }
                   return(NA)
                  })

datf <- data.frame(gtf[,c(1,4,5,7)])
colnames(datf) <- c("Chrm", "Start", "End", "Strand")
rownames(datf) <- gene_ids
datf[,"Name"] <- gene_names
datf[,"Type"] <- gene_types

dir_out <- "data/ref/gencode26/"
write.table(datf, 
            paste0(dir_out, "gencode.v26.annotation.txt"),
            row.names = TRUE, 
            col.names = NA, 
            sep = "\t", 
            quote = FALSE)

