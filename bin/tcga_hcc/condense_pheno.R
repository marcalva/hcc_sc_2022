
setwd("../../")

# read in data
out_dir <- "data/processed/tcga_hcc/expr/"
fn <- paste0(out_dir, "tcga.lihc.htseq_counts.rds")
counts <- readRDS(fn)
fn <- paste0(out_dir, "tcga.gencodev26.rds")
gencode <- readRDS(fn)


out_dir <- "data/processed/tcga_hcc/sample/"

fn <- paste0(out_dir, "tcga.clinical.rds")
clinical <- readRDS(fn)
fn <- paste0(out_dir, "tcga.sample_sheet.rds")
sample_sheet <- readRDS(fn)
fn <- paste0(out_dir, "tcga.analyte_rna.rds")
analyte.rna <- readRDS(fn)
fn <- paste0(out_dir, "tcga.slide.rds")
slide <- readRDS(fn)
fn <- paste0(out_dir, "tcga.bio_sample.rds")
bio_sample <- readRDS(fn)
fn <- paste0(out_dir, "tcga.cdr.rds")
cdr <- readRDS(fn)
fn <- paste0(out_dir, "tcga.ExtraEndpoints.rds")
ee <- readRDS(fn)


#===============================================================
# subset counts to only those with HCC (removes 10 cases)
# subset to samples that are primary tumor or normal tissue
# remove recurrent as these are from same person
#===============================================================

cases.k <- rownames(clinical)[clinical[,"histologic_diagnosis"] == "Hepatocellular Carcinoma"]
clinical <- clinical[cases.k,]

bio_sample <- bio_sample[colnames(counts),]

k <- bio_sample[,"case"] %in% cases.k
bio_sample <- bio_sample[k,]

k <- bio_sample[,"sample_type"] != "Recurrent Tumor"
bio_sample <- bio_sample[k,]

samples.k <- bio_sample[,"bcr_sample_barcode"]
cases.k <- unique(bio_sample[,"case"])

clinical <- clinical[cases.k,]
sample_sheet <- sample_sheet[samples.k,-c(1:5)]
analyte.rna <- analyte.rna[samples.k,]
slide <- slide[samples.k,-c(1:5)]
bio_sample <- bio_sample[samples.k,-c(1:3)]
colnames(cdr) <- paste0("cdr.", colnames(cdr))
cdr <- cdr[cases.k,-c(1)]
colnames(ee) <- paste0("ee.", colnames(ee))
ee <- ee[cases.k,-c(1)]

samples.pheno <- do.call(cbind, 
                         list(sample_sheet, analyte.rna, slide, bio_sample))
cases.pheno <- do.call(cbind, 
                       list(clinical, cdr, ee))

samples.pheno[,"sample_type2"] <- "Non-tumor"
ktum <- samples.pheno[,"sample_type"] == "Primary Tumor"
samples.pheno[ktum, "sample_type2"] <- "Tumor"
samples.pheno[,"sample_type2"] <- factor(samples.pheno[,"sample_type2"], levels = c("Tumor", "Non-tumor"))

# get tumor high and low by stage
tumor.low.high <- cases.pheno[,"ajcc_pathologic_tumor_stage"]
tumor.low.k <- which(tumor.low.high == "Stage I" | 
                     tumor.low.high == "Stage II")
tumor.high.k <- which(tumor.low.high == "Stage III" | 
                      tumor.low.high == "Stage IIIA" |
                      tumor.low.high == "Stage IIIB" |
                      tumor.low.high == "Stage IIIC" |
                      tumor.low.high == "Stage IV" |
                      tumor.low.high == "Stage IVA" |
                      tumor.low.high == "Stage IVB")
tumor.na <- which(tumor.low.high == "[Discrepancy]" | 
                  tumor.low.high == "[Not Available]")
tumor.low.high[tumor.low.k] <- "Low"
tumor.low.high[tumor.high.k] <- "High"
tumor.low.high[tumor.na] <- NA
cases.pheno[,"stage_low_high"] <- tumor.low.high

cases.pheno[cases.pheno[,"race"] == "[Not Available]" | 
         cases.pheno[,"race"] == "[Not Evaluated]" |
         cases.pheno[,"race"] == "[Unknown]", "race"] <- NA

cases.pheno[cases.pheno[,"age_at_diagnosis"] == "[Not Available]", "age_at_diagnosis"] <- NA

# HBV and HCV

# risk factors
rf <- cases.pheno[,"history_hepato_carcinoma_risk_factors"]
rf_spl <- lapply(rf, function(x) strsplit(x, split="\\|")[[1]])
rf_spl <- unique(do.call(c, rf_spl))

k.na1 <- grep("Not Available", cases.pheno[,"history_hepato_carcinoma_risk_factors"])
k.na2 <- grep("Unknown", cases.pheno[,"history_hepato_carcinoma_risk_factors"])
k.na <- union(k.na1, k.na2)

rfs2get <- c("HBV_risk" = "Hepatitis B", 
             "HCV_risk" = "Hepatitis C", 
             "NAFLD_risk" = "Non-Alcoholic Fatty Liver Disease", 
             "Alcohol_risk" = "Alcohol consumption", 
             "Hemochromatosis_risk" = "Hemochromatosis", 
             "SERPINA1_risk" = "Alpha-1 Antitrypsin Deficiency")

for (rfget in names(rfs2get)){
    k <- grep(rfs2get[rfget], cases.pheno[,"history_hepato_carcinoma_risk_factors"])
    cases.pheno[,rfget] <- 0
    cases.pheno[k,rfget] <- 1
    cases.pheno[k.na,rfget] <- NA
}

# serum test

cases.pheno[,"HBV"] <- 0
cases.pheno[,"HBV"] <- 0

hcnames <- cases.pheno[,"viral_hepatitis_serology"]
hcnames_spl <- lapply(hcnames, function(x) strsplit(x, split="\\|")[[1]])
hcnames_spl <- unique(do.call(c, hcnames_spl))
# sort(hcnames_spl)

k.na1 <- grep("Not Available", cases.pheno[,"viral_hepatitis_serology"])
# k.na2 <- grep("Unknown", cases.pheno[,"viral_hepatitis_serology"])
# k.na <- union(k.na1, k.na2)
k.na <- k.na1

cases.pheno[,"AntiHBc"] <- 0
k <- grep("HBV Core Antibody", cases.pheno[,"viral_hepatitis_serology"])
cases.pheno[k,"AntiHBc"] <- 1
cases.pheno[k.na,"AntiHBc"] <- NA

cases.pheno[,"AntiHBs"] <- 0
k <- grep("HBV Surface Antibody", cases.pheno[,"viral_hepatitis_serology"])
cases.pheno[k,"AntiHBs"] <- 1
cases.pheno[k.na,"AntiHBs"] <- NA

cases.pheno[,"HBsAg"] <- 0
k <- grep("Hepatitis B Surface Antigen", cases.pheno[,"viral_hepatitis_serology"])
cases.pheno[k,"HBsAg"] <- 1
cases.pheno[k.na,"HBsAg"] <- NA

cases.pheno[,"HBV_DNA"] <- 0
k <- grep("HBV DNA", cases.pheno[,"viral_hepatitis_serology"])
cases.pheno[k,"HBV_DNA"] <- 1
cases.pheno[k.na,"HBV_DNA"] <- NA

s2get <- c("HCV_Ab" = "Hepatitis  C Antibody", 
           "HCV_geno" = "HCV Genotype", 
           "HCV_RNA" = "Hepatitis C Virus RNA")

for (sget in names(s2get)){
    cases.pheno[,sget] <- 0
    k <- grep(s2get[sget], cases.pheno[,"viral_hepatitis_serology"])
    cases.pheno[k,sget] <- 1
    cases.pheno[k.na,sget] <- NA
}

cases.pheno[,"HBV_infect"] <- 0
k <- cases.pheno[,"AntiHBc"] == 1 & cases.pheno[,"HBsAg"] == 1
cases.pheno[which(k),"HBV_infect"] <- 1
cases.pheno[k.na,"HBV_infect"] <- NA

cases.pheno[,"HBV_immune"] <- 0
k <- cases.pheno[,"AntiHBc"] == 1 & cases.pheno[,"AntiHBs"] == 1
cases.pheno[which(k),"HBV_immune"] <- 1
cases.pheno[k.na,"HBV_immune"] <- NA

# HBV positive either immune or infectious
cases.pheno[,"ser_HBV_pos"] <- 0
k <- cases.pheno[,"HBV_infect"] == 1 | cases.pheno[,"HBV_immune"] == 1
cases.pheno[which(k), "ser_HBV_pos"] <- 1
cases.pheno[k.na, "ser_HBV_pos"] <- NA


# HCV positive either immune or infectious
cases.pheno[,"ser_HCV_pos"] <- 0
k <- cases.pheno[,"HCV_Ab"] == 1 | cases.pheno[,"HCV_RNA"] == 1 | cases.pheno[,"HCV_geno"] == 1
cases.pheno[which(k),"ser_HCV_pos"] <- 1
cases.pheno[k.na,"ser_HCV_pos"] <- NA


# save
out_dir <- "data/processed/tcga_hcc/sample/"

outfn <- paste0(out_dir, "samples.hcc.410.merged.rds")
saveRDS(samples.pheno, outfn)
outfn <- paste0(out_dir, "cases.hcc.361.merged.rds")
saveRDS(cases.pheno, outfn)

outfn <- paste0(out_dir, "samples.hcc.410.merged.txt")
write.table(samples.pheno, outfn, row.names=TRUE, col.names=NA, 
            quote=FALSE, sep="\t")
outfn <- paste0(out_dir, "cases.hcc.361.merged.txt")
write.table(cases.pheno, outfn, row.names=TRUE, col.names=NA, 
            quote=FALSE, sep="\t")

#==============================================================================
# Process CNV data
#==============================================================================

gist_dir <- "data/raw/tcga_hcc/gdac.LIHC-TP.CopyNumber_Gistic2/"

fn <- paste0(gist_dir, "all_lesions.conf_99.txt")
fn <- "data/raw/tcga_hcc/gdac.LIHC-TP.CopyNumber_Gistic2/all_lesions.conf_99.txt"
cnvs <- read.table(fn, header = TRUE, stringsAsFactors = FALSE, sep = "\t", check.names=FALSE)

fn <- "data/raw/tcga_hcc/gdac.LIHC-TP.CopyNumber_Gistic2/all_thresholded.by_genes.txt"
genes.band <- read.table(fn, header = TRUE, stringsAsFactors = FALSE, sep = "\t", check.names=FALSE)
genes.band <- genes.band[,c(1,3)]

k.val <- grep("values", cnvs[,1])
k.thr <- grep("values", cnvs[,1], invert = TRUE)

cnvs.val <- cnvs[k.val,]
cnvs.thr <- cnvs[k.thr,]

k.amp <- grep("Amp", cnvs.val[,1])
cnvs.val.amp <- cnvs.val[k.amp,]
k.del <- grep("Del", cnvs.val[,1])
cnvs.val.del <- cnvs.val[k.del,]

k.amp <- grep("Amp", cnvs.thr[,1])
cnvs.thr.amp <- cnvs.thr[k.amp,]
k.del <- grep("Del", cnvs.thr[,1])
cnvs.thr.del <- cnvs.thr[k.del,]

rownames(cnvs.val.amp) <- make.unique(trimws(paste0("Amp_", cnvs.val.amp[,2])), sep = "_")
rownames(cnvs.val.del) <- make.unique(trimws(paste0("Del_", cnvs.val.del[,2])), sep = "_")
rownames(cnvs.thr.amp) <- make.unique(trimws(paste0("Amp_", cnvs.thr.amp[,2])), sep = "_")
rownames(cnvs.thr.del) <- make.unique(trimws(paste0("Del_", cnvs.thr.del[,2])), sep = "_")

# add genes description
cnv.desc <- rbind(cnvs.val.amp[,2:7], cnvs.val.del[,2:7])
cnv.desc[1:4] <- apply(cnv.desc[,1:4], c(1,2), trimws)
cnv.desc[,"n_genes"] <- NA
cnv.desc[,"genes"] <- NA
for (i in rownames(cnv.desc)){
    k <- genes.band[,"Cytoband"] == cnv.desc[i, "Descriptor"]
    gs <- genes.band[k, "Gene Symbol"]
    gs <- sapply(gs, function(x) sub("\\|.*", "", x))
    cnv.desc[i,"n_genes"] <- length(gs)
    gs <- paste0(gs, collapse=",")
    cnv.desc[i,"genes"] <- gs
}


cnvs.val.amp <- cnvs.val.amp[,-c(1:9)]
cnvs.val.del <- cnvs.val.del[,-c(1:9)]
cnvs.thr.amp <- cnvs.thr.amp[,-c(1:9)]
cnvs.thr.del <- cnvs.thr.del[,-c(1:9)]

s.tum <- rownames(samples.pheno)[samples.pheno[,"sample_type"] == "Primary Tumor"]
cnvs.l <- list(cnvs.val.amp, cnvs.val.del, cnvs.thr.amp, cnvs.thr.del)
names(cnvs.l) <- c("val.amp", "val.del", "thr.amp", "thr.del")
for (i in 1:length(cnvs.l)){
    colnames(cnvs.l[[i]]) <- sapply(colnames(cnvs.l[[i]]), 
                                    function(x) 
                                        paste(strsplit(x, "-")[[1]][1:4], collapse="-"))
    sk <- intersect(s.tum, colnames(cnvs.l[[i]]))
    cnvs.l[[i]] <- cnvs.l[[i]][,sk]
}

out_dir <- "data/processed/tcga_hcc/cnv/"
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)
for (i in names(cnvs.l)){
    fn <- paste0(out_dir, "tcga.lihc.", i, ".354_tum.txt")
    write.table(cnvs.l[[i]], fn, row.names=TRUE, col.names=NA, quote=FALSE, sep = "\t")
}

fn <- paste0(out_dir, "tcga.lihc.cnv_descr.txt")
write.table(cnv.desc, fn, row.names=TRUE, col.names=NA, quote=FALSE, sep = "\t")


#==============================================================================
# Process mutation data
#==============================================================================

fn <- "data/raw/tcga_hcc/gdac.LIHC-TP.MutSigNozzleReport2CV/LIHC-TP.final_analysis_set.maf"
maf <- read.table(fn, sep="\t", stringsAsFactors=FALSE, quote="", fill=TRUE, header=TRUE)

bcs <- sapply(maf[,"Tumor_Sample_Barcode"], 
              function(x){
                  paste(strsplit(x, "-")[[1]][1:4], collapse="-") })
maf[,"Tumor_Sample_Barcode"] <- bcs

maf_samples <- unique(bcs)

# map variant types more generally
vt_names <- unique(maf[,"Variant_Classification"])
vt_map <- rep(NA, length(vt_names))
names(vt_map) <- vt_names

vt_map["Silent"] <- "Synonymous"
vt_map["Frame_Shift_Del"] <- "Frame_Shift"
vt_map["Frame_Shift_Ins"] <- "Frame_Shift"
vt_map["De_novo_Start_OutOfFrame"] <- "Frame_Shift"
vt_map["In_Frame_Del"] <- "In_Frame"
vt_map["In_Frame_Del"] <- "In_Frame"
vt_map["In_Frame_Ins"] <- "In_Frame"
vt_map["Missense_Mutation"] <- "Missense"
vt_map["Nonsense_Mutation"] <- "Nonsense"
vt_map["Splice_Site"] <- "Splice_Site"

maf[,"VariantTypeShort"] <- vt_map[maf[,"Variant_Classification"]]

k <- !is.na(maf[,"VariantTypeShort"])
maf <- maf[k,]

# keep significant genes only
fn <- "data/raw/tcga_hcc/gdac.LIHC-TP.MutSigNozzleReport2CV/sig_genes.txt"
sig_genes <- read.table(fn, sep="\t", quote="", stringsAsFactors=FALSE, header=TRUE)
k <- sig_genes[,"q"] < .1
sig_genes <- sig_genes[k,]
g <- sig_genes[,"gene"]
k <- maf[,"Hugo_Symbol"] %in% g
maf <- maf[k,]

# freq table
gene_type <- paste0(maf[,"Hugo_Symbol"], ":", maf[,"VariantTypeShort"])
gene_type.u <- unique(gene_type)

mut.type.f <- matrix(0, 
                nrow = length(gene_type.u), 
                ncol = length(maf_samples), 
                dimnames = list(gene_type.u, maf_samples))

for (i in 1:nrow(maf)){
    rn <- gene_type[i]
    cn <- maf[i,"Tumor_Sample_Barcode"]
    mut.type.f[rn,cn] = mut.type.f[rn,cn] + 1
}

mut.genes <- sapply(rownames(mut.type.f), function(x){
                    strsplit(x, ":")[[1]][1] })

mut.gene.f <- by(mut.type.f, mut.genes, FUN = function(datf){
                 colSums(datf) })
mut.gene.f <- do.call(rbind, mut.gene.f)


# gene-sample type
mut.tc.f <- matrix("None", nrow = nrow(mut.gene.f), ncol = ncol(mut.gene.f), 
                   dimnames = dimnames(mut.gene.f))

mut_types <- c("Synonymous", "In_Frame", "Missense", 
               "Splice_Site", "Frame_Shift", "Nonsense")

for (gene in sig_genes[,"gene"]){
    for (mut_type in mut_types){
        rn <- paste0(gene, ":", mut_type)
        if (rn %in% rownames(mut.type.f))
            mut.tc.f[gene, mut.type.f[rn,] == 1] <- mut_type
    }
}


s.tum <- rownames(samples.pheno)[samples.pheno[,"sample_type"] == "Primary Tumor"]
k <- intersect(s.tum, colnames(mut.type.f))
# mut.type.f <- mut.type.f[,k]
# mut.gene.f <- mut.gene.f[,k]
# mut.tc.f <- mut.tc.f[,k]

out_dir <- "data/processed/tcga_hcc/mut/"
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

fn <- paste0(out_dir, "tcga.lihc.mut.gene.type_t.txt")
write.table(mut.type.f, fn, row.names=TRUE, col.names=NA, quote=FALSE, sep="\t")

fn <- paste0(out_dir, "tcga.lihc.mut.gene.any.txt")
write.table(mut.gene.f, fn, row.names=TRUE, col.names=NA, quote=FALSE, sep="\t")

fn <- paste0(out_dir, "tcga.lihc.mut.gene.type.txt")
write.table(mut.tc.f, fn, row.names=TRUE, col.names=NA, quote=FALSE, sep="\t")

#==============================================================================
# Process broad CNV data
#==============================================================================





