
setwd("../../")

# annotations
fn <- "data/ref/gencode26/gencode.v26.annotation.txt"
gencode <- read.table(fn, header=TRUE, row.names=1, sep="\t", stringsAsFactors=FALSE)
gencode[,"gene_id"] <- sapply(rownames(gencode), function(s) strsplit(s, "\\.")[[1]][1])
k <- !duplicated(gencode[,"gene_id"])
gencode <- gencode[k,]
rownames(gencode) <- gencode[,"gene_id"]


raw_dir <- "data/raw/tcga_hcc/"

# sample sheet
fn <- paste0(raw_dir, "gdc_sample_sheet.2021-01-26.tsv")
sample_sheet <- read.table(fn, header=TRUE, stringsAsFactors=FALSE, sep="\t", 
                           check.name=FALSE)
rownames(sample_sheet) <- sample_sheet[,"Sample ID"]
sample2case <- sample_sheet[,"Case ID"]
names(sample2case) <- sample_sheet[,"Sample ID"]

# clinical patient data
fn <- "data/raw/tcga_hcc/Clinical/88b30bf9-b937-4381-8684-42bbcce98fa0/nationwidechildrens.org_clinical_patient_lihc.txt"
clin_pat <- read.table(fn, header = TRUE, row.names = 2, quote = '', 
                       stringsAsFactors = FALSE, sep = '\t')
clin_pat <- clin_pat[-c(1,2),]

# clinical follow-up (leave along for now)
fn <- "data/raw/tcga_hcc/Clinical/be89630c-6155-4e59-b472-868f367e26bf/nationwidechildrens.org_clinical_follow_up_v4.0_lihc.txt"
clin_pat.f <- read.table(fn, header = TRUE, quote = '', 
                         stringsAsFactors = FALSE, sep = '\t')
clin_pat.f <- clin_pat.f[-c(1,2),]

# endpoint
fn <- "data/raw/tcga_hcc/TCGA-CDR/TCGA-CDR.tsv"
cdr <- read.table(fn, header = TRUE, row.names = 1, sep = "\t", stringsAsFactors = FALSE)
rownames(cdr) <- cdr[,"bcr_patient_barcode"]
cdr <- cdr[cdr[,'type'] == "LIHC",,drop=FALSE]

fn <- "data/raw/tcga_hcc/TCGA-CDR/ExtraEndpoints.tsv"
ee <- read.table(fn, header = TRUE, row.names = 1, sep = "\t", stringsAsFactors = FALSE)
rownames(ee) <- ee[,"bcr_patient_barcode"]
ee <- ee[ee[,'type'] == "LIHC",,drop=FALSE]

# biospecimen analyte (for RIN)
fn <- "data/raw/tcga_hcc/Biospecimen/79f7fb36-5d9c-414f-beea-98556e26813d/nationwidechildrens.org_biospecimen_analyte_lihc.txt"
analyte <- read.table(fn, header = TRUE, quote = '', sep = "\t", 
                      na.string = "[Not Available]", stringsAsFactors = FALSE)
analyte.rna <- subset(analyte, analyte_type == "RNA")
k <- !duplicated(analyte.rna[,"bcr_sample_barcode"])
analyte.rna <- analyte.rna[k,]
rownames(analyte.rna) <- analyte.rna[,"bcr_sample_barcode"]

# biospecimen slide (for tumor percent, infiltration, etc.)
fn <- "data/raw/tcga_hcc/Biospecimen/50e9401a-3a6e-4083-babc-978b96acef12/nationwidechildrens.org_biospecimen_slide_lihc.txt"
slide <- read.table(fn, header = TRUE, quote = '', sep = "\t", 
                      na.string = "[Not Available]", stringsAsFactors = FALSE)
slide <- slide[-1,]
k <- !duplicated(slide[,"bcr_sample_barcode"])
slide <- slide[k,]
rownames(slide) <- slide[,"bcr_sample_barcode"]

# biospecimen sample (for primary, recurrent, or normal)
fn <- "data/raw/tcga_hcc/Biospecimen/c38f932e-a975-4551-838c-75851b2190a9/nationwidechildrens.org_biospecimen_sample_lihc.txt"
bio_sample <- read.table(fn, header = TRUE, quote = '', sep = "\t",
                         na.string = "[Not Available]", stringsAsFactors = FALSE)
bio_sample <- bio_sample[-1,]
k <- bio_sample[,'sample_type'] != "Blood Derived Normal"
bio_sample <- bio_sample[k,]
rownames(bio_sample) <- bio_sample[,"bcr_sample_barcode"]
case_id <- sapply(bio_sample$bcr_sample_barcode, function(x) {
                  x = strsplit(x, "-")[[1]]; 
                  paste(x[1:3], collapse='-') })
bio_sample[,"case"] <- case_id




# read counts
counts.l <- lapply(1:nrow(sample_sheet), function(x){
    message(x)
    line <- sample_sheet[x,]
    fn <- paste0(raw_dir, line[1], "/", line[2])
    cnt <- read.table(fn, header = FALSE, row.names = 1, colClasses = c("character", "integer"))
    colnames(cnt) <- line[,"Sample ID"]
    return(cnt)
})
counts.all <- do.call(cbind, counts.l)

# keep genes in gencode annotations used in snRNA-seq
gene_ids <- sapply(rownames(counts.all), function(s) strsplit(s, "\\.")[[1]][1])
rownames(counts.all) <- gene_ids
gk <- rownames(counts.all) %in% gencode[,"gene_id"]
counts.all <- counts.all[gk,]

gk <- rownames(gencode) %in% rownames(counts.all)
gencode <- gencode[gk,]

out_dir <- "data/processed/tcga_hcc/expr/"
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

outfn <- paste0(out_dir, "tcga.lihc.htseq_counts.rds")
saveRDS(counts.all, outfn)
outfn <- paste0(out_dir, "tcga.gencodev26.rds")
saveRDS(gencode, outfn)

out_dir <- "data/processed/tcga_hcc/sample/"
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)
outfn <- paste0(out_dir, "tcga.clinical.rds")
saveRDS(clin_pat, outfn)
outfn <- paste0(out_dir, "tcga.sample_sheet.rds")
saveRDS(sample_sheet, outfn)
outfn <- paste0(out_dir, "tcga.analyte_rna.rds")
saveRDS(analyte.rna, outfn)
outfn <- paste0(out_dir, "tcga.slide.rds")
saveRDS(slide, outfn)
outfn <- paste0(out_dir, "tcga.bio_sample.rds")
saveRDS(bio_sample, outfn)
outfn <- paste0(out_dir, "tcga.cdr.rds")
saveRDS(cdr, outfn)
outfn <- paste0(out_dir, "tcga.ExtraEndpoints.rds")
saveRDS(ee, outfn)

