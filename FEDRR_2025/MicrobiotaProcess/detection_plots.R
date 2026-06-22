library(phyloseq)
library(MicrobiotaProcess)
library(ggplot2)
source("detection_plot_processing/detection_plots_lib.R")

# We will use MicroBiotaProcess package MPSE format 
# so it is easy to integrate with other code
# 
#Loading ASV calling results
#otus <- as.matrix(read.csv("data/FEDRR_2025_11182025/asv_tab.samps_and_baseline.csv", row.names = 1))
otus <- as.matrix(read.csv("data/FEDRR_2025_11182025/ASVs_counts.tsv", row.names = 1, sep = "\t"))
tax <- as.matrix(read.csv("data/FEDRR_2025_11182025/ASVs_taxonomy.tsv", row.names = 1, sep = "\t"))
sam <- as.data.frame(read.csv("data/2025_metadata/metadata_P1P2.csv", row.names = 1))
colnames(otus) = sub("_.*", "", colnames(otus))

sam$mdy = sam$dateOrBaselineType
sam$mdy[sam$mdy == "BLF"] = NA
sam$mdy = lubridate::mdy(sam$mdy)
sam$mdy %>% class()
sam$mdy %>% as.Date() %>% class()
###########################
###########################
# setting up MPSE obj
# 
## Basing from Eli's code

# filter extra taxa out to match ASV tab
#tax[row.names(tax) %in% colnames(otus),] -> tax

## Generating phyloseq objects from dada2 results
OTU = otu_table(otus, taxa_are_rows = T)
TAX = phyloseq::tax_table(tax)
SAM = sample_data(sam)
st_physeq = phyloseq(OTU, TAX,SAM)
st_physeq

## Converting phyloseq object to mpse
mpse2 <- st_physeq %>% as.MPSE()
mpse2

#filter out baseline samples
mpse2 %>% dplyr::pull(dateOrBaselineType) %>% unique()
mpse2 <- mpse2 %>% dplyr::filter(dateOrBaselineType != "BLF")
    
# Filter features without annotation in Phylum or Domain
mpse2 %>% mp_extract_feature(addtaxa=T) %>% dplyr::pull(Kingdom) %>% unique()
mpse2 <- mpse2 %>% dplyr::filter(!Kingdom %in% c("k__Unknown", "k__Viridiplantae","k__Rhizaria", "k__Metazoa", "k__Stramenopila"))
mpse2 %>% mp_extract_feature(addtaxa=T) %>% dplyr::pull(Kingdom) %>% unique()
mpse2 %>% mp_extract_feature(addtaxa=T) %>% dplyr::pull(Phylum) %>% unique()
mpse2 <- mpse2 %>% dplyr::filter(!Phylum %in% c("p__un_k__Fungi"))
mpse2 %>% mp_extract_feature(addtaxa=T) %>% dplyr::pull(Phylum) %>% unique()

mpse2


#########################################################
#########################################################
#########################################################
# run the detect alg
# 
# taxonomic ID database and species of interest table
dbPath <- "data/taxon_detect_test_data/sh_general_release_dynamic_s_all_19.02.2025.Bretziella.fasta"
spp_search_terms <- read.csv("data/taxon_detect_test_data/species_search_terms.csv")


taxon_detect_list.State_Site = taxon_detect(mpse=mpse2, groupVar = c("State", "Site"), abdVar = Abundance, dbPath = dbPath, speciesSearch = spp_search_terms)
taxon_detect_list.State_Site_Lure = taxon_detect(mpse=mpse2, groupVar = c("State", "Site", "Lure"), abdVar = Abundance, dbPath = dbPath, speciesSearch = spp_search_terms)
taxon_detect_list.Site_Lure_Date = taxon_detect(mpse=mpse2, groupVar = c("State", "Site", "Lure", "mdy"), abdVar = Abundance, dbPath = dbPath, speciesSearch = spp_search_terms)

head(taxon_detect_list.Site_Lure_Date)


library(ggplot2)

# color palette for different detection types
my_detect_pal = c(
    "not detected" = "light grey",
    "sp. not in db; no unkn. spp. match gen." = "dark grey",
    "unknown spp. in gen. detected; sp. not in db" = "yellow",
    "present <0.01%" = "#fee5d9",
    "present 0.01-0.1%" = "#fcae91",
    "present 0.1-1%" = "#fb6a4a",
    "present 1-10%" = "#de2d26",
    "present >10%" = "#a50f15"
    )


#species of different priority levels
taxon_detect_list.State_Site$db_matches %>% 
    filter(priority == "high") %>% 
    pull(species_name) -> high_priority_spp
taxon_detect_list.State_Site$db_matches %>% 
    filter(priority == "moderate") %>% 
    pull(species_name) -> med_priority_spp
taxon_detect_list.State_Site$db_matches %>% 
    filter(priority == "low") %>% 
    pull(species_name) -> low_priority_spp

# pull high priority species
state_site_detections.high_priority <- taxon_detect_list.State_Site$detections_tab %>%
    filter(Species %in% high_priority_spp)

state_site_lure_detections.high_priority <- taxon_detect_list.State_Site_Lure$detections_tab %>%
    filter(Species %in% high_priority_spp)
# pull med priority species
state_site_detections.med_priority <- taxon_detect_list.State_Site$detections_tab %>%
    filter(Species %in% med_priority_spp)

state_site_lure_detections.med_priority <- taxon_detect_list.State_Site_Lure$detections_tab %>%
    filter(Species %in% med_priority_spp)
# pull low priority species
state_site_detections.low_priority <- taxon_detect_list.State_Site$detections_tab %>%
    filter(Species %in% low_priority_spp)

state_site_lure_detections.low_priority <- taxon_detect_list.State_Site_Lure$detections_tab %>%
    filter(Species %in% low_priority_spp)


taxon_detect_list.Site_Lure_Date.high <- taxon_detect_list.Site_Lure_Date$detections_tab %>%
    filter(Species %in% high_priority_spp)
taxon_detect_list.Site_Lure_Date.med <- taxon_detect_list.Site_Lure_Date$detections_tab %>%
    filter(Species %in% med_priority_spp)
taxon_detect_list.Site_Lure_Date.low <- taxon_detect_list.Site_Lure_Date$detections_tab %>%
    filter(Species %in% low_priority_spp)


# two group vars
state_site_detections.high_priority %>%
ggplot(aes(x = Site, y = Species, fill = Detections)) +
    facet_grid(~State, scales = "free_x", space = "free_x") +
    geom_tile() +
    scale_fill_manual(values = my_detect_pal) +
    scale_y_discrete(limits = rev) +
    theme_bw() +
    theme(
        legend.position = "bottom",
        axis.text.x = element_text(angle = 35, hjust = 1),
        axis.text.y = element_text(face = "italic")
    ) -> p1
p1
state_site_detections.med_priority %>%
ggplot(aes(x = Site, y = Species, fill = Detections)) +
    facet_grid(~State, scales = "free_x", space = "free_x") +
    geom_tile() +
    scale_fill_manual(values = my_detect_pal) +
    scale_y_discrete(limits = rev) +
    theme_bw() +
    theme(
        legend.position = "bottom",
        axis.text.x = element_text(angle = 35, hjust = 1),
        axis.text.y = element_text(face = "italic")
    ) -> p2
p2
state_site_detections.low_priority %>%
ggplot(aes(x = Site, y = Species, fill = Detections)) +
    facet_grid(~State, scales = "free_x", space = "free_x") +
    geom_tile() +
    scale_fill_manual(values = my_detect_pal) +
    scale_y_discrete(limits = rev) +
    theme_bw() +
    theme(
        legend.position = "bottom",
        axis.text.x = element_text(angle = 35, hjust = 1),
        axis.text.y = element_text(face = "italic")
    ) -> p3
p3

pdf("figures/FEDRR_2025_11182025/MicrobiotaProcess/detections_state_site.pdf", width = 11, height = 4)
p1
p2
p3
dev.off()



# three group vars
state_site_lure_detections.high_priority %>%
ggplot(aes(x = Lure, y = Site, fill = Detections)) +
    facet_grid(State~Species, scales = "free_y", space = "free_y") +
    geom_tile() +
    scale_fill_manual(values = my_detect_pal) +
    theme_bw() +
    theme(
        legend.position = "bottom",
        axis.text.x = element_text(angle = 35, hjust = 1)
    ) -> p1
p1
state_site_lure_detections.med_priority %>%
ggplot(aes(x = Lure, y = Site, fill = Detections)) +
    facet_grid(State~Species, scales = "free_y", space = "free_y") +
    geom_tile() +
    scale_fill_manual(values = my_detect_pal) +
    theme_bw() +
    theme(
        legend.position = "bottom",
        axis.text.x = element_text(angle = 35, hjust = 1)
    ) -> p2
p2
state_site_lure_detections.low_priority %>%
ggplot(aes(x = Lure, y = Site, fill = Detections)) +
    facet_grid(State~Species, scales = "free_y", space = "free_y") +
    geom_tile() +
    scale_fill_manual(values = my_detect_pal) +
    theme_bw() +
    theme(
        legend.position = "bottom",
        axis.text.x = element_text(angle = 35, hjust = 1)
    ) -> p3
p3

pdf("figures/FEDRR_2025_11182025/MicrobiotaProcess/detections_state_site_lure.pdf", width = 16, height = 8)
p1
p2
p3
dev.off()


taxon_detect_list.Site_Lure_Date.high %>%
    filter(State == "NH") %>%
ggplot(aes(x = Lure, y = mdy, fill = Detections)) +
    facet_grid(Site~Species, scales = "free_y", space = "free_y") +
        geom_tile() +
    scale_fill_manual(values = my_detect_pal) +
    theme_bw() +
    theme(
        legend.position = "bottom",
        axis.text.x = element_text(angle = 35, hjust = 1)
    )

taxon_detect_list.Site_Lure_Date.med %>%
    filter(State == "NH") %>%
ggplot(aes(x = Lure, y = mdy, fill = Detections)) +
    facet_grid(Site~Species, scales = "free_y", space = "free_y") +
        geom_tile() +
    scale_fill_manual(values = my_detect_pal) +
    theme_bw() +
    theme(
        legend.position = "bottom",
        axis.text.x = element_text(angle = 35, hjust = 1)
    )

taxon_detect_list.Site_Lure_Date.low %>%
    filter(State == "NH") %>%
ggplot(aes(x = Lure, y = mdy, fill = Detections)) +
    facet_grid(Site~Species, scales = "free_y", space = "free_y") +
        geom_tile() +
    scale_fill_manual(values = my_detect_pal) +
    theme_bw() +
    theme(
        legend.position = "bottom",
        axis.text.x = element_text(angle = 35, hjust = 1)
    )
