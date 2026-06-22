library(dplyr)
library(tidyr)
library(ggplot2)
library(vegan)

sp_tab = read.csv("data/2024_insect_data/species_tab.csv")
head(sp_tab)
str(sp_tab)
# there is a missing value with all zeros 
sp_tab = sp_tab[!is.na(sp_tab$Finest.ID),]

# set rownames to finest.ID for t() 
rownames(sp_tab) = sp_tab$Finest.ID

sp_tab

sp_tab.t = t(sp_tab %>% select(where(is.numeric)))
head(sp_tab.t)

rowSums(sp_tab.t)
sum(rowSums(sp_tab.t))

sp_nmds = metaMDS(sp_tab.t[rowSums(sp_tab.t) > 0,], distance = "bray", binary = F, try = 20, trymax = 100)

plot(sp_nmds)

metadata = read.csv("data/2024_metadata/collated/collated_metadata_NY-OH-IA-NH_02102025.csv")
head(metadata)

row.names(sp_tab.t) %>% sub("\\.\\.", " ", .) %>% gsub("\\.", "-", .) -> table_ids

metadata = metadata %>% filter(sampleID %in% table_ids)
head(metadata)
