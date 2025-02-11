library(dplyr)
library(tidyr)

#you can use this script to summarize sequence counts by taxonomic groups
#note that here we are reading in the asv_counts.tsv table that is the
# original output from dada2. Depending on the analysis you are performing
# you may instead wish to perform these operations on a table of average counts
# after rarefaction, e.g., the asv_tab.rarefaction_avg.csv table produced by
# avg_dist_and_div.R script.

taxonomy = read.table("data/ASVs_taxonomy.tsv", header = T)
head(taxonomy)
taxonomy$ASV = rownames(taxonomy)
asv_tab = read.table("data/ASVs_counts.tsv", header = T)
asv_tab$ASV = rownames(asv_tab)
asv_tab.tax = left_join(taxonomy, asv_tab, by = "ASV")

asv_tab.tax.long = asv_tab.tax %>% 
    pivot_longer(
        names_to = c("sample", "MID"), 
        names_sep = "_", 
        values_to = "seq_count",
        cols = where(is.numeric)
    )

#asv_tab.tax.long$sample = sub("X","", asv_tab.tax.long$sample) #this line deals with sample names that begin with "X" from importing numeric names
asv_tab.tax.long$MID = NULL

asv_tab.tax.long.taxon_summary = asv_tab.tax.long %>% 
    group_by(sample, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
    summarize(seq_count = sum(seq_count))

asv_tab.tax.long.taxon_summary
asv_tab.tax.long.taxon_sum = asv_tab.tax.long.taxon_summary %>%
    group_by(Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
    summarize(total_count = sum(seq_count))

asv_tab.taxon_summary = asv_tab.tax.long.taxon_summary %>%
    pivot_wider(names_from = "sample", values_from = seq_count)
head(asv_tab.taxon_summary)

asv_tab.taxon_summary.tots = left_join(asv_tab.taxon_summary, asv_tab.tax.long.taxon_sum)

#calculate rel abd
seqs_per_sample = colSums(asv_tab.taxon_summary.tots[,8:ncol(asv_tab.taxon_summary.tots)])
asv_tab.taxon_summary.RA = cbind(asv_tab.taxon_summary.tots[,1:7],
    apply(
        asv_tab.taxon_summary.tots[,8:ncol(asv_tab.taxon_summary.tots)], 
        1, 
        FUN = function(x) x/seqs_per_sample
    ) %>% t
)
colSums(asv_tab.taxon_summary.RA[,8:ncol(asv_tab.taxon_summary.RA)])


write.csv(asv_tab.taxon_summary.tots[order(asv_tab.taxon_summary.tots$total_count, decreasing = T),], "data/taxon_summary_table.csv", row.names = F, quote = F)
write.csv(asv_tab.taxon_summary.RA[order(asv_tab.taxon_summary.RA$total_count, decreasing = T),], "data/taxon_summary_table.rel_abd.csv", row.names = F, quote = F)

