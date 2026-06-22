library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)
library(lubridate)
source("library/library.R")

asv_tab = read.csv("data/FEDRR_all_2024/asv_tab.funnel_traps_and_spore_traps.csv", row.names = 1)
rownames(asv_tab)
asvs_tax = read.table("data/FEDRR_all_2024/ASVs_taxonomy.tsv", header = T)
tax_bs = read.table("data/FEDRR_all_2024/ASVs_taxonomy_bootstrapVals.tsv", header = T)
head(tax_bs < 90)
asvs_tax[tax_bs < 80] = NA

BRFA_asv = asvs_tax %>% filter(Genus == "g__Bretziella" & Species == "s__fagacearum") %>% row.names()

BRFA_tab = asv_tab[BRFA_asv]
rownames(BRFA_tab) = rownames(asv_tab)
BRFA_tab

rownames(BRFA_tab)[rownames(BRFA_tab) %>% grep("^X", .)]
write.csv(rownames(BRFA_tab)[rownames(BRFA_tab) %>% grep("^X", .)], "data/2024_metadata/sporeTrapSamples.csv")

write.csv(BRFA_tab, "data/2024_metadata/sporeTrapSamples.csv", row.names = T)

#The Bear brook insect trapping sites are 13-20
# positives at this site from metabarcoding include 15, 17,
# sample 1-6 are OH all of which are negative
# The remaining positives are all from Ryan's nitidulid sites