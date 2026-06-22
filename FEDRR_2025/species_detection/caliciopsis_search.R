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

caliciopsis_asvs = asvs_tax %>% filter(Genus == "g__Caliciopsis") %>% row.names()

asvs_tax[caliciopsis_asvs,]
tax_bs[caliciopsis_asvs,]

cali_tab = asv_tab[caliciopsis_asvs]
rownames(cali_tab) = rownames(asv_tab)
BRFcali_tabA_tab

rownames(cali_tab)[rownames(cali_tab) %>% grep("^X", .)]
write.csv(rownames(cali_tab)[rownames(cali_tab) %>% grep("^X", .)], "data/2024_metadata/caliciopsis_samples.csv")

write.csv(cali_tab, "data/2024_metadata/caliciopsis_samples.csv", row.names = T)

#