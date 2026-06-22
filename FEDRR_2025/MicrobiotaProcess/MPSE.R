library(phyloseq)
library(MicrobiotaProcess)
library(ggplot2)
library(dplyr)

# We will use MicroBiotaProcess package MPSE format 
# so it is easy to integrate with other code
# 
#Loading ASV calling results
#otus <- as.matrix(read.csv("data/FEDRR_2025_11182025/asv_tab.samps_and_baseline.csv", row.names = 1))
otus <- as.matrix(read.csv("data/FEDRR_2025_11182025/ASVs_counts.tsv", row.names = 1, sep = "\t"))
tax <- as.matrix(read.csv("data/FEDRR_2025_11182025/ASVs_taxonomy.tsv", row.names = 1, sep = "\t"))
sam <- as.data.frame(read.csv("data/2025_metadata/metadata_P1P2.csv", row.names = 1))
# taxonomic ID database and species of interest table
# spp_search_terms <- read.csv("data/taxon_detect_test_data/species_search_terms.csv")


sam %>% filter(dateOrBaselineType != "BLF" & dateOrBaselineType != "negative") %>%
    group_by(State, Site) %>%
    summarize(n = n())
sam %>% filter(dateOrBaselineType != "BLF" & dateOrBaselineType != "negative") %>%
    group_by(State) %>%
    summarize(n = n())
sam %>% filter(dateOrBaselineType != "BLF" & dateOrBaselineType != "negative") %>%
    select(State, Site) %>%
    distinct() %>%
    group_by(State) %>%
    summarize(n = n())
sam %>% filter(dateOrBaselineType != "BLF" & dateOrBaselineType != "negative") %>% nrow()


head(tax)
nrow(tax)
nrow(otus)
ncol(otus)
ncol(sam)
nrow(sam)
rownames(sam) %in% colnames(otus)
colnames(otus)
colnames(otus) = sub("_.*", "", colnames(otus))

sam$mdy = sam$dateOrBaselineType
sam$mdy[sam$mdy == "BLF"] = NA
sam$mdy = lubridate::mdy(sam$mdy)
sam$month = lubridate::month(sam$mdy)
unique(sam$month)
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


## Rarefaction curves
#Calculate Alpha diversity and rarefaction curves
mpse2.ad <- mpse2 %>% mp_rrarefy(raresize = 10000) %>% mp_cal_alpha(.abundance = RareAbundance) %>% mp_cal_rarecurve(.abundance = RareAbundance)

#Rarefaction Plot ungrouped
mpse2.rcp1 <- mpse2.ad %>% mp_plot_rarecurve(.rare=RareAbundanceRarecurve, .alpha = Observe)
mpse2.rcp1 <- mpse2.rcp1 + ylim(c(0, NA)) + xlim(c(0, NA))+ theme_bw()
mpse2.rcp1 + labs(x = "Number of reads", y = "Species") + theme(legend.position = "none") 

pdf("figures/FEDRR_2025_11182025/MicrobiotaProcess/rarefaction_by_sample.pdf", width = 9.5)
mpse2.rcp1 + labs(x = "Number of reads", y = "Species") + theme(legend.position = "none")
dev.off()


#Rarefaction Plot grouped by State
mpse2.rcp2 <- mpse2.ad %>% mp_plot_rarecurve(.rare=RareAbundanceRarecurve, .group = State, .alpha = Observe)
mpse2.rcp2 <- mpse2.rcp2 + ylim(c(0, NA)) + theme_bw() + scale_fill_manual(values=c("#000080","#FF2400", "#61d04f", "#2297e6","#DF536B","#28e2e5", "#cd0bbc", "#f5c710", "#9e939e", "#20A486FF", "#481F70F5","#FCA636FF"))  + scale_color_manual(values=c("#000080","#FF2400", "#61d04f", "#2297e6","#DF536B","#28e2e5", "#cd0bbc", "#f5c710", "#9e939e", "#20A486FF", "#481F70F5", "#FCA636FF"))
mpse2.rcp2 + labs(x = "Number of reads", y = "Species")

pdf("figures/FEDRR_2025_11182025/MicrobiotaProcess/rarefaction_state.pdf", width = 10)
mpse2.rcp2 + labs(x = "Number of reads", y = "Species")
dev.off()

## Alpha diversity
#Grouped by state

mpse2.ap.p.st <- mpse2.ad %>% mp_plot_alpha(.group = State, .alpha=c(Observe, Shannon, Simpson),step_increase = 0.09) + scale_fill_manual(values=c("#000080","#FF2400", "#61d04f", "#2297e6","#DF536B","#28e2e5", "#cd0bbc", "#f5c710", "#9e939e", "#20A486FF", "#481F70F5","#FCA636FF" ), guide="none") + scale_color_manual(values=c("#000080","#FF2400", "#61d04f", "#2297e6","#DF536B","#28e2e5", "#cd0bbc", "#f5c710", "#9e939e", "#20A486FF", "#481F70F5","#FCA636FF"))
mpse2.ap.p.st #Too many states to plot

pdf("figures/FEDRR_2025_11182025/MicrobiotaProcess/alpha_div_10k.pdf", width = 12)
mpse2.ap.p.st
dev.off()

#Export data by state
# Not run
#mpse2.ad.p.st <- mp_extract_sample(mpse2.ad)
#mpse2.ad.p.st <-as.data.frame(mpse2.ad.p.st) %>% select(-RareAbundanceRarecurve) %>% group_by(State) %>% dplyr::summarise(across(where(is.numeric), list(mean = mean, sd = sd)))
#write.csv(mpse2.ad.p.st, file = "/users/PAS0247/torres704/Spore_traps_project/2025/Data_analysis/Alpha_diversity_by_state.csv")


#Export data by state_site
# not run
#mpse2.ad.p.ss <- mp_extract_sample(mpse2.ad)
#mpse2.ad.p.ss <-as.data.frame(mpse2.ad.p.st) %>% select(-RareAbundanceRarecurve) %>% group_by(Group) %>% summarise(across(where(is.numeric), list(mean = mean, sd = sd)))
#write.csv(mpse2.ad.p.ss, file = "/users/PAS0247/torres704/Spore_traps_project/Eric_morrison_pipeline/Analysis/Alpha_diversity_by_group.csv")


# Taxonomy abundance (top 30 taxa) grouped by location (Group)
mpse2.ad$State_Site = paste(mpse2.ad$State, mpse2.ad$Site, sep = "_")
tax2.ab <- mpse2.ad %>% mp_cal_abundance( .abundance = RareAbundance) %>% mp_cal_abundance(.abundance=RareAbundance, .group=State_Site)

p <- tax2.ab%>% mp_plot_abundance(.abundance=RareAbundance, .group=State_Site, taxa.class = Phylum, topn = 30, plot.group = TRUE, geom = "bar") + theme(legend.position="right",legend.key.size = unit(0.5, 'cm'), legend.text = element_text(size = 10, face = "italic")) +     guides(color = guide_legend(ncol = 1), fill = guide_legend(ncol = 1))
p

g <- tax2.ab%>% mp_plot_abundance(.abundance=RareAbundance, .group=State_Site, taxa.class = Genus, topn = 30, plot.group = TRUE, geom = "bar") + theme(legend.position="right",legend.key.size = unit(0.2, 'cm'), legend.text = element_text(size = 8, face = "italic"), axis.text.x = element_text(angle = 60, hjust = 1)) +     guides(color = guide_legend(ncol = 1), fill = guide_legend(ncol = 1))
g

pdf("figures/FEDRR_2025_11182025/MicrobiotaProcess/top_abun_phyla.pdf", width = 12)
p
dev.off()

pdf("figures/FEDRR_2025_11182025/MicrobiotaProcess/top_abun_genera.pdf", width = 12)
g
dev.off()



# PCoAs for the dataset

#Standardize community data
bd.st <- tax2.ab %>% mp_decostand(.abundance=Abundance)
#Calculate the distance between the samples
dist.st <- bd.st %>% mp_cal_dist(.abundance=hellinger, distmethod="bray")
#Calculate principal coordinates analysis
pco.st <- dist.st %>% mp_cal_pcoa(.abundance=hellinger, distmethod="bray") %>% mp_cal_nmds(.abundance=hellinger, distmethod="bray")



pco.st$month[pco.st$month == 4] = "April"
pco.st$month[pco.st$month == 5] = "May"
pco.st$month[pco.st$month == 6] = "June"
pco.st$month[pco.st$month == 7] = "July"
pco.st$month = factor(pco.st$month, levels = c("April", "May", "June", "July"))
#pco.st$date_fac = as.factor(pco.st$month)
class(pco.st$month)

#Calculate adonis OH and NH locations
ado.st <- pco.st %>% mp_adonis(.abundance=hellinger, .formula=~State/Site+month+Lure, distmethod="bray", permutations=9999, action="add", by="terms")
ado.st %>% mp_extract_internal_attr(name='adonis') 


ado.plot_state <- pco.st %>% mp_adonis(.abundance=hellinger, .formula=~State, distmethod="bray", permutations=999, action="add", by="terms")
ado.plot_state %>% mp_extract_internal_attr(name='adonis') 

ado.plot_month <- pco.st %>% mp_adonis(.abundance=hellinger, .formula=~month, distmethod="bray", permutations=999, action="add", by="terms")
ado.plot_month %>% mp_extract_internal_attr(name='adonis') 

ado.plot_lure <- pco.st %>% mp_adonis(.abundance=hellinger, .formula=~Lure, distmethod="bray", permutations=999, action="add", by="terms")
ado.plot_lure %>% mp_extract_internal_attr(name='adonis') 

#Plot Pcoa
elis_palette = c("#000080","#FF2400", "#61d04f", "#2297e6","#DF536B","#28e2e5", "#cd0bbc", "#f5c710", "#9e939e", "#20A486FF", "#481F70F5","#FCA636FF")

pcoa.plot_state <- ado.plot_state %>% 
    mp_plot_ord(.ord = pcoa, .group = State, .color = State, .size =1.2,.alpha = 1, ellipse = TRUE, show.legend = FALSE, show.adonis=TRUE) +
    scale_fill_manual(values=elis_palette, guide="none") + 
    scale_color_manual(values=elis_palette, guide="none")
pcoa.plot_state + theme(plot.title = element_text(hjust = 0.5), plot.margin = margin(10, 20, 10, 10))

pcoa.plot_month <- ado.plot_month %>% 
    mp_plot_ord(.ord = pcoa, .group = month, .color = month, .size = 1.2,.alpha = 1, ellipse = TRUE, show.legend = FALSE, show.adonis=TRUE) +
    scale_fill_manual(values=elis_palette, guide="none") + 
    scale_color_manual(values=elis_palette, guide="none")
pcoa.plot_month + theme(plot.title = element_text(hjust = 0.5), plot.margin = margin(10, 20, 10, 10))

pcoa.plot_lure <- ado.plot_lure %>% 
    mp_plot_ord(.ord = pcoa, .group = Lure, .color = Lure, .size = 1.2,.alpha = 1, ellipse = TRUE, show.legend = FALSE, show.adonis=TRUE) +
    scale_fill_manual(values=elis_palette, guide="none") + 
    scale_color_manual(values=elis_palette, guide="none")
pcoa.plot_lure + theme(plot.title = element_text(hjust = 0.5), plot.margin = margin(10, 20, 10, 10))


pdf("figures/FEDRR_2025_11182025/MicrobiotaProcess/pcoa_state_month_lure.pdf",width = 9)
pcoa.plot_state + theme(plot.title = element_text(hjust = 0.5), plot.margin = margin(10, 20, 10, 10))
pcoa.plot_month + theme(plot.title = element_text(hjust = 0.5), plot.margin = margin(10, 20, 10, 10))
pcoa.plot_lure + theme(plot.title = element_text(hjust = 0.5), plot.margin = margin(10, 20, 10, 10))
dev.off()

#NMDS
nmds.plot_state <- ado.plot_state %>% 
    mp_plot_ord(.ord = nmds, .group = State, .color = State, .size =1.2,.alpha = 1, ellipse = TRUE, show.legend = FALSE, show.adonis=TRUE) +
    scale_fill_manual(values=elis_palette, guide="none") + 
    scale_color_manual(values=elis_palette, guide="none")
nmds.plot_state + theme(plot.title = element_text(hjust = 0.5), plot.margin = margin(10, 20, 10, 10))

nmds.plot_month <- ado.plot_month %>% 
    mp_plot_ord(.ord = nmds, .group = month, .color = month, .size = 1.2,.alpha = 1, ellipse = TRUE, show.legend = FALSE, show.adonis=TRUE) +
    #ggordpoint(., ) +
    scale_fill_manual(values=elis_palette, guide="none") + 
    scale_color_manual(values=elis_palette, guide="none")
nmds.plot_month + theme(plot.title = element_text(hjust = 0.5), plot.margin = margin(10, 20, 10, 10))

nmds.plot_lure <- ado.plot_lure %>% 
    mp_plot_ord(.ord = nmds, .group = Lure, .color = Lure, .size = 1.2,.alpha = 1, ellipse = TRUE, show.legend = FALSE, show.adonis=TRUE) +
    scale_fill_manual(values=elis_palette, guide="none") + 
    scale_color_manual(values=elis_palette, guide="none")
nmds.plot_lure + theme(plot.title = element_text(hjust = 0.5), plot.margin = margin(10, 20, 10, 10))


pdf("figures/FEDRR_2025_11182025/MicrobiotaProcess/nmds_state_month_lure.pdf",width = 9)
nmds.plot_state + theme(plot.title = element_text(hjust = 0.5), plot.margin = margin(10, 20, 10, 10))
nmds.plot_month + theme(plot.title = element_text(hjust = 0.5), plot.margin = margin(10, 20, 10, 10))
nmds.plot_lure + theme(plot.title = element_text(hjust = 0.5), plot.margin = margin(10, 20, 10, 10))
dev.off()

