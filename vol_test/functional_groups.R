library(dplyr)
library(tidyr)
library(ggplot2)
source("R_scripts/library.R")

top_hits = read.csv("data/vol_test/dada2_out/taxonomy_vsearch_global.top_hits.tsv", sep = "\t", header = F)
colnames(top_hits) = c("ASV", "annot", "ident", "len", "mism", "opens")
top_hits$ident.floor = floor(top_hits$ident)
#split taxon names
tax_string = strsplit(top_hits$annot, split = ";")
names(tax_string) = top_hits$ASV
tax_string.t = bind_rows(tax_string) %>% t() %>% data.frame
head(tax_string.t)
names(tax_string.t)
row.names(tax_string.t)
#bind to df
tax_string.t.asv = data.frame(
    ASV = row.names(tax_string.t),
    tax_string.t
)
colnames(tax_string.t.asv) = c("ASV", "hit", "phylum", "class", "order", "family", "genus", "species")
#now split the hit col based on the pipe
hit_string = strsplit(tax_string.t.asv$hit, split = "\\|")
head(hit_string)
names(hit_string) = tax_string.t.asv$ASV
hit_string.t = bind_rows(hit_string) %>% t() %>% data.frame
hit_string.t.asv = data.frame(
    ASV = row.names(hit_string.t),
    hit_string.t
)
colnames(hit_string.t.asv) = c("ASV", "lowest_rank", "ref_seq_ID", "species_hyp", "ref_seq_type", "kingdom")
head(hit_string.t.asv)

#join all tax info
all_tax_inf = left_join(tax_string.t.asv %>% select(-hit), hit_string.t.asv, by = "ASV") %>%
    left_join(., top_hits %>% select(ASV, ident.floor), by = "ASV") 
##################
#DONE taxonomic info ... MAKE THIS A FUNCTION!

##################
#polme traits
polme.csv = read.csv("~/repo/FungalTraits_Polme/Polme_FungalTraits_1.2_ver_16Dec_2020.csv")
trophic.p_s = data.frame(
    Genus.join = polme.csv$GENUS, #named as such bc this is the column to join by
    primary = polme.csv$primary_lifestyle,
    secondary = polme.csv$Secondary_lifestyle,
    morphology = polme.csv$Growth_form_template,
    stringsAsFactors = F
)

all_tax_inf$Genus.join = sub("g__", "", all_tax_inf$genus)

all_tax_inf.traits = left_join(all_tax_inf, trophic.p_s)

source("R_scripts/vol_test/read_data_and_split_objects.R")

#read rarefied counts
vol.asv = read.csv("data/vol_test/vol.asv_tab.rarefaction_avg.csv", row.names = 1)

#mean ASV abd of samples vol test samples
vol.post_ids = metadata.vol %>%
    filter(volumeOrWash != 0) %>%
    pull(sequenceID)
length(vol.post_ids)
nrow(metadata.vol)

vol.asv.post = vol.asv[rownames(vol.asv) %in% vol.post_ids]
vol.asv.post = vol.asv.post[colSums(vol.asv.post) != 0]
ncol(vol.asv) #2758
ncol(vol.asv.post) #2430
vol.asv.post$sample = rownames(vol.asv.post)

vol.asv.post.long = vol.asv.post %>% 
    pivot_longer(cols = -sample, names_to = "ASV", values_to = "count")

vol.asv.post.long.summarize = vol.asv.post.long %>%
    group_by(ASV) %>%
    summarize(mean = mean(count), SE = sd(count)/sqrt(n()))
vol.asv.post.long.summarize %>% filter(mean > 1) %>% nrow #94
vol.asv.post.long.summarize %>% filter(mean > 10) %>% nrow #23

#plot
vol.taxa.summary = vol.asv.post.long.summarize %>%
    #filter(mean > 10) %>%
    left_join(., all_tax_inf.traits, by = "ASV") 
vol.taxa.summary$label = sub("s__", "", vol.taxa.summary$species)

vol.taxa.summary.top_tax = vol.asv.post.long.summarize %>%
    filter(mean > 10) %>%
    left_join(., all_tax_inf.traits, by = "ASV") 
vol.taxa.summary$label = sub("s__", "", vol.taxa.summary$species)


vol.taxa.summary.primary = vol.taxa.summary %>% group_by(primary) %>% summarize(n_asvs = n())
vol.taxa.summary.primary %>% print(n = Inf)
vol.taxa.summary.morphology = vol.taxa.summary %>% group_by(morphology) %>% summarize(n_asvs = n())
vol.taxa.summary.morphology %>% print(n = Inf)

p1 = ggplot(vol.taxa.summary.primary, aes(x = reorder(primary, n_asvs), y = n_asvs)) +
    geom_col() + 
    my_gg_theme +
    labs(y = "ASVs count") +
    #scale_x_discrete(breaks = vol.taxa.summary$ASV, labels = vol.taxa.summary$lowest_rank) +
    theme(
        axis.text.x = element_text(angle = 65, hjust = 1, size = 16),
        axis.title.x = element_blank(),
        legend.title = element_text(size = 14)
    )
p1

p2 = ggplot(vol.taxa.summary.top_tax, aes(x = reorder(ASV, mean), y = mean, fill = primary)) +
    geom_col(color = "black") + 
    geom_errorbar(aes(ymin=mean-SE,ymax=mean+SE), width = 0.2) +
    my_gg_theme +
    labs(y = "Sequence count", fill = "Primary guild") +
    scale_fill_brewer(palette = "Paired") +
    scale_x_discrete(breaks = vol.taxa.summary$ASV, labels = vol.taxa.summary$lowest_rank) +
    theme(
        axis.text.x = element_text(angle = 65, hjust = 1, size = 12),
        axis.title.x = element_blank(),
        legend.title = element_text(size = 14)
    )
p2


pdf("figures/vol_test/vol_traits_asv_counts.pdf", width = 12, height = 8)
p1
p2
dev.off()

