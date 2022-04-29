
# add and format meta data, such as tumor, patient, viral 
# status for integrated data

setwd("../../")

library(Seurat)

#=========================================
# Functions
#=========================================

create_dir <- function(p){
    dir.create(p, showWarnings=FALSE, recursive=TRUE)
}

#=========================================
#=========================================

n_threads <- 8

# Gene data
gene_info <- read.table("data/ref/gencode26/gencode.v26.annotation.txt", 
                        header = TRUE, 
                        row.names = 1, 
                        stringsAsFactors = FALSE, 
                        sep = "\t")
gene_info[,"Name"] <- make.unique(gene_info[,"Name"])
symb2ens <- rownames(gene_info)
names(symb2ens) <- gene_info[,"Name"]

fn <- "data/processed/sharma_aiz/liver.int_rand.rds"
integrated <- readRDS(fn)

integrated@assays$SCT <- NULL
gc()

# annotate tumor or healthy
tum_stat <- integrated$Tumor
k <- integrated$source == "nash_hcc" & integrated$Tumor == "Tumor"
tum_stat[which(k)] <- "Tumor"
k <- integrated$source == "nash_hcc" & integrated$Tumor == "NonTumor"
tum_stat[which(k)] <- "NonTumor"
k <- integrated$source == "Sharma" & integrated$NormalvsTumor == 1
tum_stat[which(k)] <- "Tumor"
k <- integrated$source == "Sharma" & integrated$NormalvsTumor== 0
tum_stat[which(k)] <- "NonTumor"
k <- integrated$source == "Aizarani"
tum_stat[which(k)] <- "NonTumor"

integrated$TumorStat <- factor(tum_stat)

# annotate patient
pat <- as.character(integrated$Patient)
k <- integrated$source == "nash_hcc"
pat[which(k)] <- paste0("NASH_", as.character(integrated@meta.data[which(k), "Patient"]))
k <- integrated$source == "Sharma"
pat[which(k)] <- paste0("Sharma_", as.character(integrated@meta.data[which(k), "patientno"]))
k <- integrated$source == "Aizarani"
pat[which(k)] <- "Aizarani"

integrated$PatientStat <- factor(pat)

# annotate viral v non-viral
k <- is.na(integrated$ViralvsNonViral)
integrated@meta.data[which(k), "ViralvsNonViral"] <- 0

# annotate tumor type status
integrated$TumorType <- "Adjacent"
k <- which(integrated$PatientStat == "Aizarani" | integrated$PatientStat == "Sharma_0")
integrated@meta.data[k,"TumorType"] <- "Healthy"
k <- which(integrated$TumorStat == "Tumor")
integrated@meta.data[k,"TumorType"] <- "Core"
k <- which(integrated$source == "Sharma" & integrated$PNC == 2)
integrated@meta.data[k,"TumorType"] <- "Peripheral"

integrated@meta.data[,"TumorType"] <- factor(integrated@meta.data[,"TumorType"])

# remove some columns
to_rm <- c("SCT_snn_res.0.8", "SCT_snn_res.0.2", "SCT_snn_res.0.5", "SCT_snn_res.1", 
           "seurat_clusters", "Cluster", "ClusterProb", "score.debris")
k <- which(! colnames(integrated@meta.data) %in% to_rm )
integrated@meta.data <- integrated@meta.data[,k]

fn <- "data/processed/sharma_aiz/liver.int_rand.rds"
saveRDS(integrated, fn)


