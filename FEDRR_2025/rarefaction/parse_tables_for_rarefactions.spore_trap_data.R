library(dplyr)
source("library/library.R")

#asv_tab = read.table("~/FEDRR_all_2024/dada2_core/ASVs_counts.tsv", header = T)
asv_tab = read.table("data/FEDRR_all_2024/ASVs_counts.tsv", header = T)
ncol(asv_tab)
nrow(asv_tab)
colnames(asv_tab)
rownames(asv_tab)
#mod names to get rid of Illumina index ID
colnames(asv_tab) = sub("_S\\d+", "", colnames(asv_tab),perl = T)


#funnel trap metadata
metadata = read.csv("data/2024_metadata/collated/collated_metadata_NY-OH-IA-NH_02102025.csv")
head(metadata)
nrow(metadata)

# we do not have the complete metadata for the spore traps at this point
# I.e., we don't even have the sequence-ssampleID mapping
# but we can just filter these based on the sample names (all began with a 
# number so now begin with X)
# Our goal is just to output a table of all the "real" samples to send to Eli
# and for rarefaction/averaging, so we won't bother with the ST metadata for now

metadata$dateOrBaselineType %>% unique()

baseline_ids = metadata %>% filter(dateOrBaselineType == "BLO" | dateOrBaselineType == "BLF") %>% pull(SequenceID)
baseline_ids %>% length
pcrNeg_ids = metadata %>% filter(dateOrBaselineType == "negative") %>% pull(SequenceID)
pcrNeg_ids %>% length
samples_ids = metadata %>% filter(dateOrBaselineType != "BLO" & dateOrBaselineType != "PCRnegative" & !SequenceID %in% c("NY.42", "NY.35")) %>% pull(SequenceID)
samples_ids %>% length
samples_ids
#note the samps does not include Eli's samples
#we need to negative filter based on BL and PCR neg
# also add Ny.42 and NY.35
toss_these = c(baseline_ids, pcrNeg_ids, "NY.42", "NY.35")
length(toss_these)
ncol(asv_tab) - length(toss_these)
nrow(metadata) - length(toss_these) + length( grep("^X", colnames(asv_tab) ))
#444
######################
#all samples
asv_tab.keeps = asv_tab[,!colnames(asv_tab) %in% toss_these]
ncol(asv_tab.keeps) == ncol(asv_tab) - length(toss_these)
ncol(asv_tab.keeps)
nrow(asv_tab.keeps)

metadata %>% filter(num_reads < 1000)

# the ncol is 4 greater than expected from subtracting the lengths bc there are 
# four PCR negatives that do not appear in the ASV tab bc they had no reads at all
# there are also four samples missing based on the 444 above. These are samps 
# that had zero reads
metadata %>% filter(num_reads == 0 & !dateOrBaselineType %in% c("BLO", "BLF", "negative"))

#how many samples with gt 1k seqs
sum(colSums(asv_tab.keeps) > 999)
#437
colSums(asv_tab.keeps) %>% sort
#
#how many nonzero asvs
sum(rowSums(asv_tab.keeps) > 0)
#32855

#Checking just the ST samples
asv_tab.ST = asv_tab[,grep("^X", colnames(asv_tab)) ]
ncol(asv_tab.ST)
#34
sum(colSums(asv_tab.ST) > 999)
#34
colSums(asv_tab.ST) %>% sort
#
#how many nonzero asvs
sum(rowSums(asv_tab.ST) > 0)
#13355
sum(rowSums(asv_tab.ST) > 0)/sum(rowSums(asv_tab.keeps) > 0)
#0.406483
asv_tab.nST = asv_tab[,-grep("^X", colnames(asv_tab)) ]
ncol(asv_tab.nST)
#504
sum(colSums(asv_tab.nST) > 999)
#451
colSums(asv_tab.nST) %>% sort
#
#how many nonzero asvs
sum(rowSums(asv_tab.nST) > 0)
#32135
sum(rowSums(asv_tab.nST) > 0)/sum(rowSums(asv_tab.keeps) > 0)
#0.9780855

#
#remove 0 count asvs and samples with lt 1k seqs and transpose. Note we will also output the un filtered/unrarefied samples for purposes of *detection*
#note need to filter samples before ASVs otherwise may retain some more zero count ASVs
asv_tab.keeps = asv_tab.keeps[, colSums(asv_tab.keeps) > 999]
asv_tab.keeps = asv_tab.keeps[rowSums(asv_tab.keeps) > 0, ] %>% t() 

nrow(asv_tab.keeps)
ncol(asv_tab.keeps)




max(rowSums(asv_tab.keeps))/min(rowSums(asv_tab.keeps)) 
#2795.356

#1000 rarefactions should be OK to get a representative sample, but the baseline and samps ratios are around 1k...
#will definitely want to work with the unfiltered tables for purposes of *detection*


###############
#perform rarefactions
#
write.csv(asv_tab.keeps, "data/FEDRR_all_2024/asv_tab.funnel_traps_and_spore_traps.csv", quote = F)

