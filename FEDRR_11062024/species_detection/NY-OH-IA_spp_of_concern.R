library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)
source("library/library.R")

asv_tab.bl = read.csv("data/FEDRR_all_2024/asv_tab.baseline.csv", row.names = 1)
#asv_tab.bl$SequenceID = rownames(asv_tab.bl)
asv_tab.samps = read.csv("data/FEDRR_all_2024/asv_tab.samps.csv", row.names = 1)
#asv_tab.samps$SequenceID = rownames(asv_tab.samps)
asvs_tax = read.table("data/FEDRR_all_2024/ASVs_taxonomy.tsv", header = T)
tax_bs = read.table("data/FEDRR_all_2024/ASVs_taxonomy_bootstrapVals.tsv", header = T)
head(tax_bs < 90)
asvs_tax[tax_bs < 80] = NA

asv_tab.samps = t(asv_tab.samps) %>% data.frame
asv_tab.samps$ASV = rownames(asv_tab.samps)
asv_tab.bl = t(asv_tab.bl) %>% data.frame
asv_tab.bl$ASV = rownames(asv_tab.bl)
asv_tab.samps_bl = full_join(asv_tab.samps, asv_tab.bl, by = "ASV") %>% mutate(across(everything(), ~replace_na(., 0))) %>% data.frame

metadata = read.csv("data/2024_metadata/collated/collated_metadata_NY-OH-IA_02102025.csv")
sp_search_terms = read.table("data/2024_metadata/NY_spp_concern_search_terms.txt", header = F)

asvs_gen_sp_search = data.frame(
    ASV = rownames(asvs_tax),
    Genus = sub("g__", "", asvs_tax$Genus),
    Species = sub("s__", "", asvs_tax$Species),
    stringsAsFactors = F
)


# formatting up search terms

test_split = function(x){
    if(grepl("_", x)){
        return( str_split_i(x, "_", 1) )
    } else {
        return(x)
    }
}
        
ny_gen_search = sapply(sp_search_terms$V1, test_split) %>% unname

ny_sp_search = mapply(\(x, y) sub(x, "", y), ny_gen_search, sp_search_terms$V1) %>% sub("^_", "", .) %>% unname 
ny_sp_search[ny_sp_search == ""] = NA

ny_search_terms = data.frame(
    Genus = ny_gen_search,
    Species = ny_sp_search,
    stringsAsFactors = F
)
#
#

#########################
# setting up search tables
# wel will have one for genus
# and one for species
# and mark >0 samples as yes
# 

asv_tab.samps.tax = left_join(asv_tab.samps_bl, asvs_gen_sp_search, by = "ASV")
head(asv_tab.samps.tax)


asv_tab.samps.gen = asv_tab.samps.tax %>% 
    filter(Genus %in% ny_search_terms$Genus) %>% 
    group_by(Genus) %>% 
    summarize(across(where(is.numeric), sum) )
asv_tab.samps.spp = asv_tab.samps.tax %>% 
    filter(Genus %in% ny_search_terms$Genus & Species %in% ny_search_terms$Species & !is.na(Species)) %>% 
    group_by(Genus, Species) %>% 
    summarize(across(where(is.numeric), sum) )


asv_tab.samps.gen.char = asv_tab.samps.gen %>%
    mutate(across(where(is.numeric), .fns = ~ifelse(. > 0, "genus", "absent"))) %>%
    data.frame()

asv_tab.samps.spp.char = asv_tab.samps.spp %>%
    mutate(across(where(is.numeric), .fns = ~ifelse(. > 0, "species", "absent"))) %>%
    data.frame()

#join with the genus character indicators first
ny_search_terms.samps = left_join(
    ny_search_terms,
    asv_tab.samps.gen.char,
    by = "Genus"
) %>% data.frame

#then test for where we found a species and change the value where applicable
for(i in 1:nrow(asv_tab.samps.spp.char) ){
    for(j in 3:ncol(asv_tab.samps.spp.char) ){
        if(asv_tab.samps.spp.char[i,j] == "species"){
            ny_search_terms.samps[ny_search_terms.samps$Genus == asv_tab.samps.spp.char[i,1] & ny_search_terms.samps$Species == asv_tab.samps.spp.char[i,2], j] = "species"
        }
    }
}
sum(ny_search_terms.samps == "species", na.rm = T)
#69
sum(ny_search_terms.samps == "genus", na.rm = T)
#241

#replace NAs with absent
ny_search_terms.samps.rp_na = ny_search_terms.samps %>%
    mutate(across(c(3:ncol(ny_search_terms.samps)), ~replace_na(., "absent")) ) %>%
    data.frame()
sum(ny_search_terms.samps.rp_na == "species", na.rm = T)
sum(ny_search_terms.samps.rp_na == "genus", na.rm = T)

ny_search_terms.long = left_join(
    ny_search_terms.samps.rp_na %>% 
        pivot_longer(cols = c(-Genus, - Species), names_to = "SequenceID", values_to = "gen_sp"),
    metadata %>% 
        select(SequenceID, State, Site, Risk, date, Lure, BaselineType, samplePeriod),
    by = "SequenceID"
) 
ny_search_terms.long$spec_name = paste(ny_search_terms.long$Genus, ny_search_terms.long$Species)
ny_search_terms.long$spec_name = sub(" NA", "", ny_search_terms.long$spec_name)
ny_search_terms.long$spec_name = sub("_", " ", ny_search_terms.long$spec_name)
ny_search_terms.long$spec_name = sub("-occidentalis", "", ny_search_terms.long$spec_name)
#ny_search_terms.long$spec_name = sub("Ceratocystis", "Breziella", ny_search_terms.long$spec_name)

ny_search_terms.long[ny_search_terms.long$BaselineType == 'BLO', "date"] = "BL"


#uniq_dates = ny_search_terms.long$date %>% unique() %>% str_sort(., numeric = T)


p1 = ggplot(ny_search_terms.long %>% filter(State == "NY" & spec_name != "Ophiostoma"), 
       aes(x = factor(samplePeriod, levels = c("BL", 1,2,3,4,5,6)), y = Lure, fill = factor(gen_sp, levels = c("species", "genus", "absent")) )
) +
    geom_tile() + 
    scale_fill_manual(
        values = c("red", "yellow", "grey"), 
        labels = c("species present", "genus present", "not detected")
    ) +
    labs(fill = "Presence/absence:", x = "Sample period") +
    facet_grid(Site~spec_name, scales = "free_x", labeller = label_wrap_gen(12) ) +
    my_gg_theme +
    theme(
        axis.text.x = element_text(angle = 55, hjust = 1),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 15),
        legend.position= "bottom"
    )
p1

pdf("figures/FEDRR_all_2024/NY_spp_concern.pdf", width = 34, height = 11)
p1
dev.off()


p2 = ggplot(ny_search_terms.long %>% filter(State == "OH" & spec_name != "Ophiostoma"), 
       aes(x = factor(samplePeriod, levels = c("BL", 1,2,3,4,5,6)), y = Lure, fill = factor(gen_sp, levels = c("species", "genus", "absent")) )
) +
    geom_tile() + 
    scale_fill_manual(
        values = c("red", "yellow", "grey"), 
        labels = c("species present", "genus present", "not detected")
    ) +
    labs(fill = "Presence/absence:", x = "Sample period") +
    facet_grid(Site~spec_name, scales = "free_x", labeller = label_wrap_gen(12) ) +
    my_gg_theme +
    theme(
        axis.text.x = element_text(angle = 55, hjust = 1),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 15),
        legend.position= "bottom"
    )
p2

pdf("figures/FEDRR_all_2024/OH_spp_concern.pdf", width = 34, height = 11)
p2
dev.off()

p3 = ggplot(ny_search_terms.long %>% filter(State == "IA" & spec_name != "Ophiostoma"), 
       aes(x = factor(samplePeriod, levels = c("BL", 1,2,3,4,5)), y = Lure, fill = factor(gen_sp, levels = c("species", "genus", "absent")) )
) +
    geom_tile() + 
    scale_fill_manual(
        values = c("yellow", "grey"), 
        labels = c("genus present", "not detected")
    ) +
    labs(fill = "Presence/absence:", x = "Sample period") +
    facet_grid(Site~spec_name, scales = "free_x", labeller = label_wrap_gen(12) ) +
    my_gg_theme +
    theme(
        axis.text.x = element_text(angle = 55, hjust = 1),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 15),
        legend.position= "bottom"
    )
p3

pdf("figures/FEDRR_all_2024/IA_spp_concern.pdf", width = 34, height = 11)
p3
dev.off()
