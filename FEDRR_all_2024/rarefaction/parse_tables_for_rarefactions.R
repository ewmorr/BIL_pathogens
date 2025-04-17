library(dplyr)
source("library/library.R")


asv_tab = read.table("data/FEDRR_all_2024/ASVs_counts.tsv", header = T)
ncol(asv_tab)
#538
nrow(asv_tab)
# 33678
colnames(asv_tab)
#rownames(asv_tab)
#right away get rid of Eli spore trap samples. She is doing that analysis
ST_samps = colnames(asv_tab)[grep("^X", colnames(asv_tab))]
length(ST_samps)
#34 ST samps
asv_tab = asv_tab[,!colnames(asv_tab) %in% ST_samps]
ncol(asv_tab)
#504

#mod names to get rid of Illumina index ID
colnames(asv_tab) = sub("_S\\d+", "", colnames(asv_tab),perl = T)

#metadata
metadata = read.csv("data/2024_metadata/collated/collated_metadata_NY-OH-IA-NH_02102025.csv")
head(metadata)
nrow(metadata)
#512

################################################
################################################
# we pull sample IDs based on the planned comparisons
# samples to PCR negatives
# samples to baselines (all BLO)
# NH samples to BLF and BLO baselines
# NH BLF and BLO baselines to each other

#sample type IDs
samples_ids = metadata %>% filter(dateOrBaselineType != "BLO" & dateOrBaselineType != "BLF" & dateOrBaselineType != "negative" & !SequenceID %in% c("NY.42", "NY.35")) %>% pull(SequenceID) # NY.35 and 42 are a duplicate (sample run twice) and an untraceable sample ID, respectively
samples_ids %>% length
#411

samps_and_baseline_ids = metadata %>% filter(dateOrBaselineType != "negative" & dateOrBaselineType != "BLF") %>% pull(SequenceID)
samps_and_baseline_ids %>% length
#478

baseline_ids = metadata %>% filter(dateOrBaselineType == "BLO" | dateOrBaselineType == "BLF") %>% pull(SequenceID)
baseline_ids %>% length
#66

samps_and_pcrNeg_ids = metadata %>% filter(dateOrBaselineType != "BLO" & dateOrBaselineType != "BLF") %>% pull(SequenceID)
samps_and_pcrNeg_ids %>% length
#446

samps_and_baseline_ids.NH = metadata %>% filter(dateOrBaselineType != "negative" & State == "NH") %>% pull(SequenceID)
samps_and_baseline_ids.NH %>% length
#131

baseline_ids.NH = metadata %>% filter((dateOrBaselineType == "BLO" | dateOrBaselineType == "BLF") & State == "NH") %>% pull(SequenceID)
baseline_ids.NH %>% length
#30
baseline_ids.NH

#sample type IDs
samples_ids.NH = metadata %>% filter(dateOrBaselineType != "BLO" & dateOrBaselineType != "BLF" & dateOrBaselineType != "negative" & State == "NH") %>% pull(SequenceID) 
samples_ids.NH %>% length
#101


######################
######################
#trap catch samples (i.e., the real deal)
asv_tab.samps = asv_tab[,colnames(asv_tab) %in% samples_ids]
ncol(asv_tab.samps) == length(samples_ids)
ncol(asv_tab.samps)
#406
samples_ids[!samples_ids  %in% colnames(asv_tab.samps)]
# five samples with no reads passing QC (reducing to 406 potential samples)

#how many samples with lt 1k seqs
asv_tab.samps[,colSums(asv_tab.samps) > 999] %>% colnames 
asv_tab.samps[,colSums(asv_tab.samps) > 999] %>% ncol
#403 (of 406 that had good seqs)
colSums(asv_tab.samps) %>% sort

#
#how many nonzero asvs
sum(rowSums(asv_tab.samps) > 0)
#31236
sum(rowSums(asv_tab.samps) > 1)
#30726
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
# five samples with no reads passing QC (reducing to 313 potential samples)

#how many samples with lt 1k seqs
asv_tab.samps_and_baseline[,colSums(asv_tab.samps_and_baseline) > 999] %>% colnames %>% length
#445 
colSums(asv_tab.samps_and_baseline) %>% sort

#
#how many nonzero asvs
sum(rowSums(asv_tab.samps_and_baseline) > 0)
#32127
sum(rowSums(asv_tab.samps_and_baseline) > 1)
#31741
#
#remove 0 count asvs and samples with lt 1k seqs and transpose. Note we will also output the un filtered/unrarefied samples for purposes of *detection*
asv_tab.samps_and_baseline = asv_tab.samps_and_baseline[, colSums(asv_tab.samps_and_baseline) > 999]
asv_tab.samps_and_baseline = asv_tab.samps_and_baseline[rowSums(asv_tab.samps_and_baseline) > 0,] %>% t() 
ncol(asv_tab.samps_and_baseline)
nrow(asv_tab.samps_and_baseline)

######################
#baseline samples  only
asv_tab.baseline = asv_tab[,colnames(asv_tab) %in% baseline_ids]
ncol(asv_tab.baseline) == length(baseline_ids)
baseline_ids[!baseline_ids  %in% colnames(asv_tab.baseline)]
# five samples with no reads passing QC (reducing to 313 potential samples)

#how many samples with lt 1k seqs
asv_tab.baseline[,colSums(asv_tab.baseline) > 999] %>% colnames %>% length
#41 
colSums(asv_tab.baseline) %>% sort

#
#how many nonzero asvs
sum(rowSums(asv_tab.baseline) > 0)
#4646
sum(rowSums(asv_tab.baseline) > 1)
#4305
#
#remove 0 count asvs and samples with lt 500 seqs and transpose. Note we will also output the un filtered/unrarefied samples for purposes of *detection*
asv_tab.baseline = asv_tab.baseline[, colSums(asv_tab.baseline) > 499]
asv_tab.baseline = asv_tab.baseline[rowSums(asv_tab.baseline) > 0,] %>% t() 
ncol(asv_tab.baseline)
nrow(asv_tab.baseline)

######################
#trap catch AND PCR negatives 
asv_tab.samps_and_pcrNeg = asv_tab[,colnames(asv_tab) %in% samps_and_pcrNeg_ids]
ncol(asv_tab.samps_and_pcrNeg) == length(samps_and_pcrNeg_ids)
samps_and_pcrNeg_ids[!samps_and_pcrNeg_ids  %in% colnames(asv_tab.samps_and_pcrNeg)]
# nine samples with no reads passing QC 

#how many samples with lt 1k seqs
asv_tab.samps_and_pcrNeg[,colSums(asv_tab.samps_and_pcrNeg) > 999] %>% colnames %>% length
#409 
colSums(asv_tab.samps_and_pcrNeg) %>% sort

#
#how many nonzero asvs
sum(rowSums(asv_tab.samps_and_pcrNeg) > 0)
#31239
sum(rowSums(asv_tab.samps_and_pcrNeg) > 1)
#30730
#
#remove 0 count asvs and samples with lt 1k seqs and transpose. Note we will also output the un filtered/unrarefied samples for purposes of *detection*
asv_tab.samps_and_pcrNeg = asv_tab.samps_and_pcrNeg[, colSums(asv_tab.samps_and_pcrNeg) > 999]
asv_tab.samps_and_pcrNeg = asv_tab.samps_and_pcrNeg[rowSums(asv_tab.samps_and_pcrNeg) > 0,] %>% t() 
ncol(asv_tab.samps_and_pcrNeg)
nrow(asv_tab.samps_and_pcrNeg)


######################
#trap catch AND baseline samples IN NH
asv_tab.samps_and_baseline.NH = asv_tab[,colnames(asv_tab) %in% samps_and_baseline_ids.NH]
ncol(asv_tab.samps_and_baseline.NH) == length(samps_and_baseline_ids.NH)
samps_and_baseline_ids.NH[!samps_and_baseline_ids.NH  %in% colnames(asv_tab.samps_and_baseline.NH)]
# one samples with no reads passing QC 

#how many samples with lt 1k seqs
asv_tab.samps_and_baseline.NH[,colSums(asv_tab.samps_and_baseline.NH) > 999] %>% colnames %>% length
#114 
colSums(asv_tab.samps_and_baseline.NH) %>% sort
asv_tab.samps_and_baseline.NH[,colSums(asv_tab.samps_and_baseline.NH) > 499] %>% colnames %>% length
#120

#
#how many nonzero asvs
sum(rowSums(asv_tab.samps_and_baseline.NH) > 0)
#18230
sum(rowSums(asv_tab.samps_and_baseline.NH) > 1)
#16758
#
#remove 0 count asvs and samples with lt 1k seqs and transpose. Note we will also output the un filtered/unrarefied samples for purposes of *detection*
asv_tab.samps_and_baseline.NH = asv_tab.samps_and_baseline.NH[, colSums(asv_tab.samps_and_baseline.NH) > 499]
asv_tab.samps_and_baseline.NH = asv_tab.samps_and_baseline.NH[rowSums(asv_tab.samps_and_baseline.NH) > 0,] %>% t() 
ncol(asv_tab.samps_and_baseline.NH)
nrow(asv_tab.samps_and_baseline.NH)


######################
#baseline samples IN NH
asv_tab.baseline.NH = asv_tab[,colnames(asv_tab) %in% baseline_ids.NH ]
ncol(asv_tab.baseline.NH) == length(baseline_ids.NH )
baseline_ids.NH [!baseline_ids.NH   %in% colnames(asv_tab.baseline.NH)]
# all samples passing QC 

#how many samples with lt 1k seqs
asv_tab.baseline.NH[,colSums(asv_tab.baseline.NH) > 999] %>% colnames %>% length
#14
colSums(asv_tab.baseline.NH) %>% sort
# the samples that don't pass the 1k filter account for more than half of total, 
# and we would like to compare these. Let's reduce the filter to min 500
# We will do the same above for the NH samps plus bl
asv_tab.baseline.NH[,colSums(asv_tab.baseline.NH) > 499] %>% colnames %>% length
#20

#
#how many nonzero asvs
sum(rowSums(asv_tab.baseline.NH) > 0)
#1580
sum(rowSums(asv_tab.baseline.NH) > 1)
#1419
#
#remove 0 count asvs and samples with lt 1k seqs and transpose. Note we will also output the un filtered/unrarefied samples for purposes of *detection*
asv_tab.baseline.NH = asv_tab.baseline.NH[, colSums(asv_tab.baseline.NH) > 499]
asv_tab.baseline.NH = asv_tab.baseline.NH[rowSums(asv_tab.baseline.NH) > 0,] %>% t() 
ncol(asv_tab.baseline.NH)
nrow(asv_tab.baseline.NH)


######################
#trap catch samples IN NH (for comparison to insects)
asv_tab.samps.NH = asv_tab[,colnames(asv_tab) %in% samples_ids.NH]
ncol(asv_tab.samps.NH) == length(samples_ids.NH)
samples_ids.NH[!samples_ids.NH  %in% colnames(asv_tab.samps.NH)]
# one samples with no reads passing QC 

#how many samples with lt 1k seqs
asv_tab.samps.NH[,colSums(asv_tab.samps.NH) > 999] %>% colnames %>% length
#100 
colSums(asv_tab.samps.NH) %>% sort
asv_tab.samps.NH[,colSums(asv_tab.samps.NH) > 499] %>% colnames %>% length
#100

#
#how many nonzero asvs
sum(rowSums(asv_tab.samps.NH) > 0)
#17966
sum(rowSums(asv_tab.samps.NH) > 1)
#16463
#
#remove 0 count asvs and samples with lt 1k seqs and transpose. Note we will also output the un filtered/unrarefied samples for purposes of *detection*
asv_tab.samps.NH = asv_tab.samps.NH[, colSums(asv_tab.samps.NH) > 499]
asv_tab.samps.NH = asv_tab.samps.NH[rowSums(asv_tab.samps.NH) > 0,] %>% t() 
ncol(asv_tab.samps.NH)
nrow(asv_tab.samps.NH)

#####################################
# The tables:

# asv_tab.samps
# asv_tab.samps_and_baseline
# asv_tab.samps_and_pcrNeg
# asv_tab.samps_and_baseline.NH
# asv_tab.baseline.NH
# asv_tab.samps.NH



max(rowSums(asv_tab.samps))/min(rowSums(asv_tab.samps)) 
#2795.356
max(rowSums(asv_tab.samps_and_baseline))/min(rowSums(asv_tab.samps_and_baseline)) 
#2844.312
max(rowSums(asv_tab.samps_and_pcrNeg))/min(rowSums(asv_tab.samps_and_pcrNeg)) 
#2795.356
max(rowSums(asv_tab.samps_and_baseline.NH))/min(rowSums(asv_tab.samps_and_baseline.NH))
#1762.574
max(rowSums(asv_tab.baseline.NH))/min(rowSums(asv_tab.baseline.NH))
#236.9521
max(rowSums(asv_tab.samps.NH))/min(rowSums(asv_tab.samps.NH))
#917

#1000 rarefactions should be OK to get a representative sample, but the baseline and samps ratios are around 2-3k...
#will definitely want to work with the unfiltered tables for purposes of *detection*


###############
#perform rarefactions
#
write.csv(asv_tab.samps, "data/FEDRR_all_2024/asv_tab.samps.csv", quote = F)
write.csv(asv_tab.samps_and_baseline, "data/FEDRR_all_2024/asv_tab.samps_and_baseline.csv", quote = F)
write.csv(asv_tab.samps_and_pcrNeg, "data/FEDRR_all_2024/asv_tab.samps_and_pcrNeg.csv", quote = F)
write.csv(asv_tab.baseline, "data/FEDRR_all_2024/asv_tab.baseline.csv", quote = F)
write.csv(asv_tab.samps_and_baseline.NH, "data/FEDRR_all_2024/asv_tab.samps_and_baseline.NH.csv", quote = F)
write.csv(asv_tab.baseline.NH, "data/FEDRR_all_2024/asv_tab.baseline.NH.csv", quote = F)
write.csv(asv_tab.samps.NH, "data/FEDRR_all_2024/asv_tab.samps.NH.csv", quote = F)



