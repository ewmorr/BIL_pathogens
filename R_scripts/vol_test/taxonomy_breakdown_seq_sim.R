library(dplyr)
library(ggplot2)
source("R_scripts/library.R")

top_hits = read.csv("data/vol_test/dada2_out/taxonomy_vsearch_global.top_hits.tsv", sep = "\t", header = F)

source("R_scripts/vol_test/read_data_and_split_objects.R")
asv_counts = data.frame(
    ASV = row.names(asv_tab),
    seqs = asv_tab %>% rowSums
)

head(top_hits)
ncol(top_hits)
colnames(top_hits) = c("ASV", "annot", "ident", "len", "mism", "opens")
head(top_hits)

top_hits$ident.floor = floor(top_hits$ident)
head(top_hits)

top_hits %>% nrow()

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
head(tax_string.t.asv)

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
    left_join(., top_hits %>% select(ASV, ident.floor), by = "ASV") %>%
    left_join(., asv_counts, by = "ASV")
head(all_tax_inf)
all_tax_inf$ASV.percent = 1/nrow(all_tax_inf)
all_tax_inf$seqs.percent = all_tax_inf$seqs/sum(all_tax_inf$seqs)
    


#counts for plots
all_tax_inf %>%
    group_by(ident.floor) %>% 
    summarise(n = n()) %>%
    pull(n) %>% 
    sum

top_hits.species = all_tax_inf %>%
    filter(ident.floor > 98) %>%
    group_by(
        kingdom,
        phylum,
        class,
        order,
        family,
        genus,
        species
    ) %>% 
    summarise(
        n_ASVs = n(), 
        n_seqs = sum(seqs), 
        percent_ASVs = sum(ASV.percent), 
        percent_seqs = sum(seqs.percent)
    )
top_hits.genus = all_tax_inf %>%
    filter(ident.floor > 96) %>%
    group_by(
        kingdom,
        phylum,
        class,
        order,
        family,
        genus
    ) %>% 
    summarise(
        n_ASVs = n(), 
        n_seqs = sum(seqs), 
        percent_ASVs = sum(ASV.percent), 
        percent_seqs = sum(seqs.percent)
    )

top_hits.family = all_tax_inf %>%
    filter(ident.floor > 93) %>%
    group_by(
        kingdom,
        phylum,
        class,
        order,
        family
    ) %>% 
    summarise(
        n_ASVs = n(), 
        n_seqs = sum(seqs), 
        percent_ASVs = sum(ASV.percent), 
        percent_seqs = sum(seqs.percent)
)

top_hits.order = all_tax_inf %>%
    filter(ident.floor > 90) %>%
    group_by(
        kingdom,
        phylum,
        class,
        order
    ) %>% 
    summarise(
        n_ASVs = n(), 
        n_seqs = sum(seqs), 
        percent_ASVs = sum(ASV.percent), 
        percent_seqs = sum(seqs.percent)
    )

sum(all_tax_inf$seqs.percent)
sum(top_hits.family$percent_ASVs)
sum(top_hits.family$percent_seqs)
sum(top_hits.order$percent_ASVs)
sum(top_hits.order$percent_seqs)

#write files
write.csv(all_tax_inf[order(all_tax_inf$seqs, decreasing = T),], "data/vol_test/ASV_top_hits.csv", row.names = F)

write.csv(top_hits.species[order(top_hits.species$n_ASVs, decreasing = T),], "data/vol_test/species.ASV_top_hits.csv", row.names = F)

write.csv(top_hits.genus[order(top_hits.genus$n_ASVs, decreasing = T),], "data/vol_test/genusr.ASV_top_hits.csv", row.names = F)

write.csv(top_hits.family[order(top_hits.family$n_ASVs, decreasing = T),], "data/vol_test/family.ASV_top_hits.csv", row.names = F)

write.csv(top_hits.order[order(top_hits.order$n_ASVs, decreasing = T),], "data/vol_test/order.ASV_top_hits.csv", row.names = F)


top_hits.n.ordered = top_hits.n[order(top_hits.n$perc, decreasing = T),]



p1 = ggplot(top_hits.n.ordered[1:10,], aes(x = reorder(order, perc), y = n)) +
    geom_col(width = 0.5) +
    labs(x = "Order", y = "ASVs assigned") +
    my_gg_theme +
    theme(
        axis.text = element_text(size = 14),
        axis.text.x = element_text(angle = 65, hjust = 1)
    )
p1

pdf("figures/vol_test/top_10_orders_gt_97_perc_sim.pdf", width = 10, height = 8)
p1
dev.off()
