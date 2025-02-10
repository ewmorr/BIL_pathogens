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

source("R_scripts/vol_test/read_data_and_split_objects.R")

#read rarefied counts
wash.asv = read.csv("data/vol_test/wash.asv_tab.rarefaction_avg.csv", row.names = 1)
head(wash.asv)
vol.asv = read.csv("data/vol_test/vol.asv_tab.rarefaction_avg.csv", row.names = 1)

#mean ASV abd of samples pot bleach wash
wash.post_ids = metadata.wash %>%
    filter(volumeOrWash == "post") %>%
    pull(sequenceID)

wash.asv.post = wash.asv[rownames(wash.asv) %in% wash.post_ids]
wash.asv.post = wash.asv.post[colSums(wash.asv.post) != 0]
ncol(wash.asv) #1555
ncol(wash.asv.post) #519
wash.asv.post$sample = rownames(wash.asv.post)

wash.asv.post.long = wash.asv.post %>% 
    pivot_longer(cols = -sample, names_to = "ASV", values_to = "count")

wash.asv.post.long.summarize = wash.asv.post.long %>%
    group_by(ASV) %>%
    summarize(mean = mean(count), SE = sd(count)/sqrt(n()))
wash.asv.post.long.summarize %>% filter(mean > 1) %>% nrow #23

#plot
wash.taxa.summary = wash.asv.post.long.summarize %>%
    filter(mean > 1) %>%
    left_join(., all_tax_inf, by = "ASV") 
wash.taxa.summary$label = sub("s__", "", wash.taxa.summary$species)

p1 = ggplot(wash.taxa.summary, aes(x = reorder(ASV, mean), y = mean, fill = as.factor(ident.floor))) +
        geom_col() + 
    geom_errorbar(aes(ymin=mean-SE,ymax=mean+SE), width = 0.2) +
    my_gg_theme +
    labs(y = "Sequence count", fill = "Match sequence\nsimilarity (%)") +
    scale_fill_brewer(palette = "Paired") +
    scale_x_discrete(breaks = wash.taxa.summary$ASV, labels = wash.taxa.summary$lowest_rank) +
    theme(
        axis.text.x = element_text(angle = 65, hjust = 1),
        axis.title.x = element_blank(),
        legend.title = element_text(size = 14)
    )
p1

pdf("figures/vol_test/wash_post_bleach_top_taxa.pdf", width = 12, height = 8)
p1
dev.off()



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
    filter(mean > 10) %>%
    left_join(., all_tax_inf, by = "ASV") 
vol.taxa.summary$label = sub("s__", "", vol.taxa.summary$species)


p2 = ggplot(vol.taxa.summary, aes(x = reorder(ASV, mean), y = mean, fill = as.factor(ident.floor))) +
    geom_col() + 
    geom_errorbar(aes(ymin=mean-SE,ymax=mean+SE), width = 0.2) +
    my_gg_theme +
    labs(y = "Sequence count", fill = "Match sequence\nsimilarity (%)") +
    scale_fill_brewer(palette = "Paired") +
    scale_x_discrete(breaks = vol.taxa.summary$ASV, labels = vol.taxa.summary$lowest_rank) +
    theme(
        axis.text.x = element_text(angle = 65, hjust = 1, size = 12),
        axis.title.x = element_blank(),
        legend.title = element_text(size = 14)
    )
p2

pdf("figures/vol_test/vol_post_bleach_top_taxa.pdf", width = 12, height = 8)
p2
dev.off()

