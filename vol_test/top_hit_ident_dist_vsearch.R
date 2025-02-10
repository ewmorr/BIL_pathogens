library(dplyr)
library(ggplot2)

top_hits = read.csv("data/vol_test/dada2_out/taxonomy_vsearch_global.top_hits.tsv", sep = "\t", header = F)

head(top_hits)
ncol(top_hits)
colnames(top_hits) = c("ASV", "annot", "ident", "len", "mism", "opens")
head(top_hits)

top_hits$ident.floor = floor(top_hits$ident)
head(top_hits)

top_hits %>% nrow()

top_hits %>%
    group_by(ident.floor) %>% 
    summarise(n = n()) %>%
    pull(n) %>% 
    sum

top_hits.n = top_hits %>%
    group_by(ident.floor) %>% 
    summarise(n = n())

top_hits.n$perc = top_hits.n$n/nrow(top_hits)

p1 = ggplot(top_hits.n, aes(x = ident.floor, y = perc*100)) +
    geom_col() +
    labs(x = "Sequence similarity (%)", y = "ASVs assigned (%)") +
    my_gg_theme +
    theme(
        axis.text = element_text(size = 14)
    )
p1

##################
#seq counts

source("R_scripts/vol_test/read_data_and_split_objects.R")
asv_counts = asv_tab %>% 
    rowSums
head(asv_counts)
length(asv_counts)
sum(asv_counts == 0)


asv_sim.counts = left_join(
    data.frame(ASV = names(asv_counts), counts = asv_counts),
    top_hits,
    by = "ASV"
)
head(asv_sim.counts)

top_hits.sum = asv_sim.counts %>%
    group_by(ident.floor) %>% 
    summarise(sum = sum(counts))
head(top_hits.sum)

total_seqs = sum(asv_counts)
total_seqs
top_hits.sum$perc = top_hits.sum$sum/total_seqs
head(top_hits.sum)

p2 = ggplot(top_hits.sum, aes(x = ident.floor, y = perc*100)) +
    geom_col() +
    labs(x = "Sequence similarity (%)", y = "Sequences assigned (%)") +
    my_gg_theme +
    theme(
        axis.text = element_text(size = 14)
    )
p2

pdf("figures/vol_test/all_sequences.seq_sim_assigned.pdf", width = 12, height = 6)
p1
p2
dev.off()


