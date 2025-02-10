library(dplyr)
source("library/library.R")

asv_tab = read.table("~/FEDRR_all_2024/dada2_core/ASVs_counts.tsv", header = T)
asv_tab = read.table("data/FEDRR_all_2024/ASVs_counts.tsv", header = T)
ncol(asv_tab)
nrow(asv_tab)
colnames(asv_tab)
rownames(asv_tab)
#mod names to get rid of Illumina index ID
colnames(asv_tab) = sub("_S\\d+", "", colnames(asv_tab),perl = T)


metadata = read.csv("~/FEDRR_all_2024/collated_metadata_NY-OH-IA_02102025.csv")
metadata = read.csv("data/2024_metadata/collated/collated_metadata_NY-OH-IA_02102025.csv")
head(metadata)
nrow(metadata)

baseline_ids = metadata %>% filter(BaselineType == "BLO") %>% pull(SequenceID)
baseline_ids %>% length
pcrNeg_ids = metadata %>% filter(BaselineType == "PCRnegative") %>% pull(SequenceID)
pcrNeg_ids %>% length
samples_ids = metadata %>% filter(BaselineType != "BLO" & BaselineType != "PCRnegative" & !SequenceID %in% c("NY.42", "NY.35")) %>% pull(SequenceID)
samples_ids %>% length
samps_and_baseline_ids = metadata %>% filter(BaselineType != "PCRnegative") %>% pull(SequenceID)
samps_and_baseline_ids %>% length

######################
#trap wash baseline
asv_tab.baseline = asv_tab[,colnames(asv_tab) %in% baseline_ids]
ncol(asv_tab.baseline) == length(baseline_ids)
ncol(asv_tab.baseline)
nrow(asv_tab.baseline)
#how many samples with gt 1k seqs
sum(colSums(asv_tab.baseline) > 999)
#17
colSums(asv_tab.baseline) %>% sort
#
#how many nonzero asvs
sum(rowSums(asv_tab.baseline) > 0)
#2305
#
#remove 0 count asvs and samples with lt 1k seqs and transpose. Note we will also output the un filtered/unrarefied samples for purposes of *detection*
#note need to filter samples before ASVs otherwise may retain some more zero count ASVs
asv_tab.baseline = asv_tab.baseline[, colSums(asv_tab.baseline) > 999]
asv_tab.baseline = asv_tab.baseline[rowSums(asv_tab.baseline) > 0, ] %>% t() 

######################
#PCR negatives
asv_tab.pcr= asv_tab[,colnames(asv_tab) %in% pcrNeg_ids]
ncol(asv_tab.pcr) == length(pcrNeg_ids)
pcrNeg_ids[!pcrNeg_ids  %in% colnames(asv_tab.pcr)]
# the three samples that are missing did not have any sequences that passed QC

#how many samples with gt 1k seqs
asv_tab.pcr[colSums(asv_tab.pcr) > 999] %>% ncol()
#5
colSums(asv_tab.pcr) %>% sort
#there are two further samples with zero sequences
#We will filter at 400 instead of 1k to retain three more samples (and the drop is 417 to 46)
#how many samples with gt 1k seqs
asv_tab.pcr[colSums(asv_tab.pcr) > 400] %>% ncol()
#9
#how many nonzero asvs
asv_tab.pcr[rowSums(asv_tab.pcr) > 0,] %>% nrow()
#145
#
#remove 0 count asvs and samples with lt 1k seqs and transpose. Note we will also output the un filtered/unrarefied samples for purposes of *detection*
asv_tab.pcr = asv_tab.pcr[, colSums(asv_tab.pcr) > 400]
asv_tab.pcr = asv_tab.pcr[rowSums(asv_tab.pcr) > 0, ] %>% t()
ncol(asv_tab.pcr)
#132
sum(colSums(asv_tab.pcr) > 0)
#132
nrow(asv_tab.pcr)
#9

######################
#trap catch samples (i.e., the real deal)
asv_tab.samps = asv_tab[,colnames(asv_tab) %in% samples_ids]
ncol(asv_tab.samps) == length(samples_ids)
ncol(asv_tab.samps)
#317
samples_ids[!samples_ids  %in% colnames(asv_tab.samps)]
# four samples with no reads passing QC (reducing to 317 potential samples)

#how many samples with lt 1k seqs
asv_tab.samps[,colSums(asv_tab.samps) > 999] %>% colnames 
asv_tab.samps[,colSums(asv_tab.samps) > 999] %>% ncol
#313 (of 317 that had good seqs)
colSums(asv_tab.samps) %>% sort

#
#how many nonzero asvs
sum(rowSums(asv_tab.samps) > 0)
#24389
sum(rowSums(asv_tab.samps) > 1)
#24226
#
#remove 0 count asvs and samples with lt 1k seqs and transpose. Note we will also output the un filtered/unrarefied samples for purposes of *detection*
asv_tab.samps = asv_tab.samps[, colSums(asv_tab.samps) > 999]
asv_tab.samps = asv_tab.samps[rowSums(asv_tab.samps) > 0,] %>% t() 
ncol(asv_tab.samps)
nrow(asv_tab.samps)

######################
#trap catch AND baseline samples 
asv_tab.samps_and_baseline = asv_tab[,colnames(asv_tab) %in% samps_and_baseline_ids]
ncol(asv_tab.samps_and_baseline) == length(samps_and_baseline_ids)
samps_and_baseline_ids[!samps_and_baseline_ids  %in% colnames(asv_tab.samps_and_baseline)]
# four samples with no reads passing QC (reducing to 313 potential samples)

#how many samples with lt 1k seqs
asv_tab.samps_and_baseline[,colSums(asv_tab.samps_and_baseline) > 999] %>% colnames %>% length
#331 (313 + 28 = 331)
colSums(asv_tab.samps_and_baseline) %>% sort

#
#how many nonzero asvs
sum(rowSums(asv_tab.samps_and_baseline) > 0)
#24816
sum(rowSums(asv_tab.samps_and_baseline) > 1)
#24701
#
#remove 0 count asvs and samples with lt 1k seqs and transpose. Note we will also output the un filtered/unrarefied samples for purposes of *detection*
asv_tab.samps_and_baseline = asv_tab.samps_and_baseline[, colSums(asv_tab.samps_and_baseline) > 999]
asv_tab.samps_and_baseline = asv_tab.samps_and_baseline[rowSums(asv_tab.samps_and_baseline) > 0,] %>% t() 
ncol(asv_tab.samps_and_baseline)
nrow(asv_tab.samps_and_baseline)


max(rowSums(asv_tab.baseline))/min(rowSums(asv_tab.baseline)) 
#904
max(rowSums(asv_tab.pcr))/min(rowSums(asv_tab.pcr)) 
#98
max(rowSums(asv_tab.samps))/min(rowSums(asv_tab.samps)) 
#1179
max(rowSums(asv_tab.samps_and_baseline))/min(rowSums(asv_tab.samps_and_baseline))
#2786

#1000 rarefactions should be OK to get a representative sample, but the baseline and samps ratios are around 1k...
#will definitely want to work with the unfiltered tables for purposes of *detection*


###############
#perform rarefactions
#
write.csv(asv_tab.baseline, "data/FEDRR_all_2024/asv_tab.baseline.csv", quote = F)
write.csv(asv_tab.pcr, "data/FEDRR_all_2024/asv_tab.pcr.csv", quote = F)
write.csv(asv_tab.samps, "data/FEDRR_all_2024/asv_tab.samps.csv", quote = F)
write.csv(asv_tab.samps_and_baseline, "data/FEDRR_all_2024/asv_tab.samps_and_baseline.csv", quote = F)



min_seqs.baseline = min(rowSums(asv_tab.baseline)) 
rarefactions_list.baseline = multiple_subsamples(x = asv_tab.baseline, depth = min_seqs.baseline, iterations = 1000) 
saveRDS(rarefactions_list.baseline, "data/FEDRR_11062024/processed_tables/rarefactions.baseline.rds")

min_seqs.pcr = min(rowSums(asv_tab.pcr)) 
rarefactions_list.pcr = multiple_subsamples(x = asv_tab.pcr, depth = min_seqs.pcr, iterations = 1000) 
saveRDS(rarefactions_list.pcr, "data/FEDRR_11062024/processed_tables/rarefactions.pcr.rds")

min_seqs.samps = min(rowSums(asv_tab.samps)) 
rarefactions_list.samps = multiple_subsamples(x = asv_tab.samps, depth = min_seqs.samps, iterations = 1000) 
saveRDS(rarefactions_list.samps, "data/FEDRR_11062024/processed_tables/rarefactions.samps.rds")

min_seqs.samps_and_baseline = min(rowSums(asv_tab.samps_and_baseline)) 
rarefactions_list.samps_and_baseline = multiple_subsamples(x = asv_tab.samps_and_baseline, depth = min_seqs.samps_and_baseline, iterations = 1000) 
saveRDS(rarefactions_list.samps_and_baseline, "data/FEDRR_11062024/processed_tables/rarefactions.samps_and_baseline.rds")




########################
########################
# also running rarefaction of all samples at 400 seqs per sample
# for direct comparisons of samps to controls
#

asv_tab.min = asv_tab[, colSums(asv_tab) > 400]
asv_tab.min = asv_tab.min[rowSums(asv_tab.min) > 0, ] %>% t()
nrow(asv_tab.min)
#244
ncol(asv_tab.min)
#22883
min_seqs.min = min(rowSums(asv_tab.min)) 

# this is mem dumping, would need to run on premise, but can just combine the other tables
rarefactions_list.full_tab_min450 = multiple_subsamples(x = asv_tab.min, depth = min_seqs.min, iterations = 500) 

saveRDS(rarefactions_list.full_tab_min450, "data/FEDRR_11062024/processed_tables/rarefactions.full_tab_min450.rds")


