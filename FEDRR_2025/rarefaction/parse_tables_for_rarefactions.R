library(dplyr)
source("library/library.R")


asv_tab = read.table("data/FEDRR_2025_11182025/ASVs_counts.tsv", header = T)
ncol(asv_tab)
#182
nrow(asv_tab)
# 10187
colnames(asv_tab)
#rownames(asv_tab)

#mod names to get rid of Illumina index ID
colnames(asv_tab) = sub("_S\\d+", "", colnames(asv_tab),perl = T)

#metadata
metadata = read.csv("data/2025_metadata/metadata_P1P2.csv")
head(metadata)
nrow(metadata)
#192

################################################
################################################
# we pull sample IDs based on the planned comparisons
# samples to PCR negatives
# samples to baselines (all BLO)
# NH samples to BLF and BLO baselines
# NH BLF and BLO baselines to each other

#sample type IDs
samples_ids = metadata %>% 
    filter(dateOrBaselineType != "BLF" & dateOrBaselineType != "negative" & !SequenceID %in% c("FEDRR25p2.140")) %>% pull(SequenceID) # FEDRR25p2.140 is a duplicate (sample run twice) 
samples_ids %>% length
#167

samps_and_baseline_ids = metadata %>% filter(dateOrBaselineType != "negative" & !SequenceID %in% c("FEDRR25p2.140")) %>% pull(SequenceID)
samps_and_baseline_ids %>% length
#182

baseline_ids = metadata %>% filter(dateOrBaselineType == "BLF") %>% pull(SequenceID)
baseline_ids %>% length
#15

samps_and_pcrNeg_ids = metadata %>% filter(dateOrBaselineType != "BLF" &  !SequenceID %in% c("FEDRR25p2.140")) %>% pull(SequenceID)
samps_and_pcrNeg_ids %>% length
#176

samps_and_baseline_ids.NH = metadata %>% filter(dateOrBaselineType != "negative" & State == "NH") %>% pull(SequenceID)
samps_and_baseline_ids.NH %>% length
#75

baseline_ids.NH = metadata %>% filter((dateOrBaselineType == "BLO" | dateOrBaselineType == "BLF") & State == "NH") %>% pull(SequenceID)
baseline_ids.NH %>% length
#15
baseline_ids.NH

#sample type IDs
samples_ids.NH = metadata %>% filter(dateOrBaselineType != "BLO" & dateOrBaselineType != "BLF" & dateOrBaselineType != "negative" & State == "NH") %>% pull(SequenceID) 
samples_ids.NH %>% length
#60


######################
######################
#trap catch samples (i.e., the real deal)
asv_tab.samps = asv_tab[,colnames(asv_tab) %in% samples_ids]
ncol(asv_tab.samps) == length(samples_ids)
length(samples_ids)
ncol(asv_tab.samps)
#166
samples_ids[!samples_ids  %in% colnames(asv_tab.samps)]
# FEDRR25p1.29 one samples with no reads passing QC (reducing to 166 potential samples)

#how many samples with lt 1k seqs
asv_tab.samps[,colSums(asv_tab.samps) > 999] %>% colnames 
asv_tab.samps[,colSums(asv_tab.samps) > 999] %>% ncol
#158 (of 166 that had good seqs)
colSums(asv_tab.samps) %>% sort

#
#how many nonzero asvs
sum(rowSums(asv_tab.samps) > 0)
#10149
sum(rowSums(asv_tab.samps) > 1)
#10129
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
# one samples with no reads passing QC (reducing to 181 potential samples)

#how many samples with lt 1k seqs
asv_tab.samps_and_baseline[,colSums(asv_tab.samps_and_baseline) > 999] %>% colnames %>% length
#172 (out of 181)
colSums(asv_tab.samps_and_baseline) %>% sort

#
#how many nonzero asvs
sum(rowSums(asv_tab.samps_and_baseline) > 0)
#10187
sum(rowSums(asv_tab.samps_and_baseline) > 1)
#10179
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
# 0 samples with no reads passing QC 

#how many samples with lt 1k seqs
asv_tab.baseline[,colSums(asv_tab.baseline) > 999] %>% colnames %>% length
#14 
colSums(asv_tab.baseline) %>% sort
# one nonpassing

#
#how many nonzero asvs
sum(rowSums(asv_tab.baseline) > 0)
#2081
sum(rowSums(asv_tab.baseline) > 1)
#1886
#
#remove 0 count asvs and samples with lt 1k seqs and transpose. Note we will also output the un filtered/unrarefied samples for purposes of *detection*
asv_tab.baseline = asv_tab.baseline[, colSums(asv_tab.baseline) > 999]
asv_tab.baseline = asv_tab.baseline[rowSums(asv_tab.baseline) > 0,] %>% t() 
ncol(asv_tab.baseline)
nrow(asv_tab.baseline)

######################
#trap catch AND PCR negatives 
asv_tab.samps_and_pcrNeg = asv_tab[,colnames(asv_tab) %in% samps_and_pcrNeg_ids]
ncol(asv_tab.samps_and_pcrNeg) == length(samps_and_pcrNeg_ids)
samps_and_pcrNeg_ids[!samps_and_pcrNeg_ids  %in% colnames(asv_tab.samps_and_pcrNeg)]
# 10 samples with no reads passing QC 

#how many samples with lt 1k seqs
asv_tab.samps_and_pcrNeg[,colSums(asv_tab.samps_and_pcrNeg) > 999] %>% colnames %>% length
#158 
colSums(asv_tab.samps_and_pcrNeg) %>% sort
#8 nonpassing

#
#how many nonzero asvs
sum(rowSums(asv_tab.samps_and_pcrNeg) > 0)
#10149
sum(rowSums(asv_tab.samps_and_pcrNeg) > 1)
#10129
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
#74 
colSums(asv_tab.samps_and_baseline.NH) %>% sort
asv_tab.samps_and_baseline.NH[,colSums(asv_tab.samps_and_baseline.NH) > 499] %>% colnames %>% length
#74

#
#how many nonzero asvs
sum(rowSums(asv_tab.samps_and_baseline.NH) > 0)
#7727
sum(rowSums(asv_tab.samps_and_baseline.NH) > 1)
#7252
#
#remove 0 count asvs and samples with lt 1k seqs and transpose. Note we will also output the un filtered/unrarefied samples for purposes of *detection*
asv_tab.samps_and_baseline.NH = asv_tab.samps_and_baseline.NH[, colSums(asv_tab.samps_and_baseline.NH) > 999]
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
# 1 nonpassing
asv_tab.baseline.NH[,colSums(asv_tab.baseline.NH) > 999] %>% colnames %>% length
#14

#
#how many nonzero asvs
sum(rowSums(asv_tab.baseline.NH) > 0)
#2081
sum(rowSums(asv_tab.baseline.NH) > 1)
#1886
#
#remove 0 count asvs and samples with lt 1k seqs and transpose. Note we will also output the un filtered/unrarefied samples for purposes of *detection*
asv_tab.baseline.NH = asv_tab.baseline.NH[, colSums(asv_tab.baseline.NH) > 999]
asv_tab.baseline.NH = asv_tab.baseline.NH[rowSums(asv_tab.baseline.NH) > 0,] %>% t() 
ncol(asv_tab.baseline.NH)
nrow(asv_tab.baseline.NH)



#####################################
# The tables:

# asv_tab.samps
# asv_tab.samps_and_baseline
# asv_tab.samps_and_pcrNeg
# asv_tab.samps_and_baseline.NH
# asv_tab.baseline.NH
# # asv_tab.baseline




max(rowSums(asv_tab.samps))/min(rowSums(asv_tab.samps)) 
#604.6748
max(rowSums(asv_tab.samps_and_baseline))/min(rowSums(asv_tab.samps_and_baseline)) 
#604.6748
max(rowSums(asv_tab.samps_and_pcrNeg))/min(rowSums(asv_tab.samps_and_pcrNeg)) 
#604.6748
max(rowSums(asv_tab.baseline))/min(rowSums(asv_tab.baseline)) 
#149.4465
max(rowSums(asv_tab.samps_and_baseline.NH))/min(rowSums(asv_tab.samps_and_baseline.NH))
#287.2201
max(rowSums(asv_tab.baseline.NH))/min(rowSums(asv_tab.baseline.NH))
#149.4465

#1000 rarefactions should be OK to get a representative sample
#will definitely want to work with the unfiltered tables for purposes of *detection*


###############
#perform rarefactions
#
write.csv(asv_tab.samps, "data/FEDRR_2025_11182025//asv_tab.samps.csv", quote = F)
write.csv(asv_tab.samps_and_baseline, "data/FEDRR_2025_11182025/asv_tab.samps_and_baseline.csv", quote = F)
write.csv(asv_tab.samps_and_pcrNeg, "data/FEDRR_2025_11182025/asv_tab.samps_and_pcrNeg.csv", quote = F)
write.csv(asv_tab.baseline, "data/FEDRR_2025_11182025/asv_tab.baseline.csv", quote = F)
write.csv(asv_tab.samps_and_baseline.NH, "data/FEDRR_2025_11182025/asv_tab.samps_and_baseline.NH.csv", quote = F)
write.csv(asv_tab.baseline.NH, "data/FEDRR_2025_11182025/asv_tab.baseline.NH.csv", quote = F)




