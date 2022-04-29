
# map case IDs and sample IDs. Get paired samples
# Get tumor and non-tumor case IDs

setwd("../../")

fn <- "data/processed/tcga_hcc/sample/samples.hcc.410.merged.rds"
samples <- readRDS(fn)
samples[,"sample_type2"] <- factor(samples[,"sample_type2"], levels = c("Tumor", "Non-tumor"))

# separate samples by tumor and non-tumor
k_pt <- samples[,"sample_type2"] == "Tumor"
samples.pt <- samples[k_pt,]
k_nt <- samples[,"sample_type2"] == "Non-tumor"
samples.nt <- samples[k_nt,]

# case to sample map for tumor and non-tumor separately
cid2sid_tum <- rownames(samples.pt)
names(cid2sid_tum) <- samples.pt[,"case"]
cid2sid_nt <- rownames(samples.nt)
names(cid2sid_nt) <- samples.nt[,"case"]

# get case IDs with both tumor and non-tumor samples
cid_pair <- intersect(names(cid2sid_tum), names(cid2sid_nt))

# paired sample IDs
cid2sid_tum_pair <- cid2sid_tum[cid_pair]
cid2sid_nt_pair <- cid2sid_nt[cid_pair]

# get sample IDs with both tumor and non-tumor samples
smaps <- list("cid2sid_tum" = cid2sid_tum, 
              "cid2sid_nt" = cid2sid_nt,
              "cid2sid_tum_pair" = cid2sid_tum_pair, 
              "cid2sid_nt_pair" = cid2sid_nt_pair)
out_fn <-  "data/processed/tcga_hcc/sample/sam_case_map.rds"
saveRDS(smaps, out_fn)

