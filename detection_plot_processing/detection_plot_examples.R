library(phyloseq)
library(MicrobiotaProcess)
source("detection_plot_processing/detection_plots_lib.R")

# We will use MicroBiotaProcess package MPSE format 
# so it is easy to integrate with other code
# 
#Loading ASV calling results
otus <- as.matrix(read.csv("data/taxon_detect_test_data/asv_tab.samps.NH.csv", row.names = 1))
tax <- as.matrix(read.csv("data/taxon_detect_test_data/ASVs_taxonomy.tsv", row.names = 1, sep = "\t"))
sam <- as.data.frame(read.csv("data/taxon_detect_test_data/metadata.samps.NH.csv", row.names = 1))
# taxonomic ID database and species of interest table
dbPath <- "data/taxon_detect_test_data/sh_general_release_dynamic_s_all_19.02.2025.Bretziella.fasta"
spp_search_terms <- read.csv("data/taxon_detect_test_data/species_search_terms.csv")

###########################
###########################
# setting up MPSE obj
# 
## Basing from Eli's code

# filter extra taxa out to match ASV tab
tax[row.names(tax) %in% colnames(otus),] -> tax

## Generating phyloseq objects from dada2 results
OTU = otu_table(otus, taxa_are_rows = F)
TAX = phyloseq::tax_table(tax)
SAM = sample_data(sam)
st_physeq = phyloseq(OTU, TAX,SAM)
st_physeq

## Converting phyloseq object to mpse
mpse2 <- st_physeq %>% as.MPSE()
mpse2

# Filter features without annotation in Phylum or Domain
mpse2 %>% mp_extract_feature(addtaxa=T) %>% dplyr::pull(Kingdom) %>% unique()
mpse2 <- mpse2 %>% dplyr::filter(!Kingdom %in% c("k__Unknown", "k__Viridiplantae","k__Rhizaria", "k__Metazoa", "k__Stramenopila"))
mpse2 %>% mp_extract_feature(addtaxa=T) %>% dplyr::pull(Kingdom) %>% unique()
mpse2 %>% mp_extract_feature(addtaxa=T) %>% dplyr::pull(Phylum) %>% unique()
mpse2 <- mpse2 %>% dplyr::filter(!Phylum %in% c("p__un_k__Fungi"))
mpse2 %>% mp_extract_feature(addtaxa=T) %>% dplyr::pull(Phylum) %>% unique()

mpse2

# the taxon detect function allows choice of which column to use for relative
# abundance calculation. We recommend using observed sequence abundance before 
# any subsampling (rarefaction) to avoid loss of sensitivity (i.e., artificially
# excluding detected species). However, if desired rarefied abundance can be 
# calculated (e.g., with mp_cal_alpha()) and input as abundance variable.

#Calculate rarefied abundance
#mpse2.ad <- mpse2 %>% mp_rrarefy() %>% mp_cal_alpha(.abundance = RareAbundance) %>% mp_cal_rarecurve(.abundance = RareAbundance)


################################################################################
################################################################################
# Run taxa detection. groupVar indicates which variables in the MPSE to use as
# grouping variables. We have tested with 1-3 categorical variables. 
# The abdVar variable indicates which column to use for relative abundance calc. 
# This is typically Abundance. Note that because of internal variable handing the
# grouping variables should be quoted and the abdVar variable should not. The
# speciesSearch is a data.frame and must contain a column named "species_name" 
# which is used as the species search terms. The search terms are case sensitive 
# and underscore separated (e.g., Bretziella_fagacearum). The species list may
# contain an arbitrary number of additional columns.
# The function returns a two element list. the first element is the original
# speciesSearch table with two logical columns appended indicating whether the
# species name and genus name were found in the reference database. The second 
# table indicates presence or absence of the search species in the provided 
# groups. "not detected" indicates the species was present in the ref db but
# was not detected in sample/group, "present (x)%" indicates detection and the
# % relative abundance, "genus present/absent (sp. not in db)" indicates that 
# the search term was not found in the ref db but the genus was present or 
# absent
# 
taxon_detect_list.Site = taxon_detect(mpse=mpse2, groupVar = c("Site"), abdVar = Abundance, dbPath = dbPath, speciesSearch = spp_search_terms)
taxon_detect_list.Site_Lure = taxon_detect(mpse=mpse2, groupVar = c("Site", "Lure"), abdVar = Abundance, dbPath = dbPath, speciesSearch = spp_search_terms)
taxon_detect_list.Site_Lure_Date = taxon_detect(mpse=mpse2, groupVar = c("Site", "Lure", "date_m.dd"), abdVar = Abundance, dbPath = dbPath, speciesSearch = spp_search_terms)

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
taxon_detect_list.Site$db_matches %>% 
    filter(priority == "high") %>% 
    pull(species_name) -> high_priority_spp
taxon_detect_list.Site$db_matches %>% 
    filter(priority == "moderate") %>% 
    pull(species_name) -> med_priority_spp
taxon_detect_list.Site$db_matches %>% 
    filter(priority == "low") %>% 
    pull(species_name) -> low_priority_spp

# pull high priority species
site_detections.high_priority <- taxon_detect_list.Site$detections_tab %>%
    filter(Species %in% high_priority_spp)

site_lure_detections.high_priority <- taxon_detect_list.Site_Lure$detections_tab %>%
    filter(Species %in% high_priority_spp)

site_lure_date_detections.high_priority <- taxon_detect_list.Site_Lure_Date$detections_tab %>%
    filter(Species %in% high_priority_spp)

########################################################################
########################################################################
# plot

# single group var
ggplot(site_detections.high_priority, aes(x = Site, y = Species, fill = Detections)) +
    geom_tile() +
    scale_fill_manual(values = my_detect_pal) +
    theme_bw()
# two group vars
ggplot(site_lure_detections.high_priority, aes(x = Site, y = Species, fill = Detections)) +
    facet_wrap(~Lure) +
    geom_tile() +
    scale_fill_manual(values = my_detect_pal) +
    theme_bw() +
    theme(
        legend.position = "bottom",
        axis.text.x = element_text(angle = 35, hjust = 1)
    )
# three group vars
ggplot(site_lure_date_detections.high_priority, aes(x = Site, y = date_m.dd, fill = Detections)) +
    facet_grid(Lure~Species) +
    geom_tile() +
    scale_fill_manual(values = my_detect_pal) +
    theme_bw() +
    theme(
        legend.position = "bottom",
        axis.text.x = element_text(angle = 35, hjust = 1)
    )

# all spp single group var
ggplot(taxon_detect_list.Site$detections_tab, aes(x = Site, y = Species, fill = Detections)) +
    geom_tile() +
    scale_fill_manual(values = my_detect_pal) +
    theme_bw()

