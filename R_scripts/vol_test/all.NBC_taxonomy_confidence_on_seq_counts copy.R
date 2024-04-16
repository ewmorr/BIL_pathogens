library(dplyr)
library(tidyr)
library(ggplot2)
source("R_scripts/library.R")

asv_tax = read.table("data/vol_test/dada2_out/ASVs_taxonomy.tsv", sep = "\t", header = T)
head(asv_tax)
colnames(asv_tax)[1] = "ASV"
asv_boot = read.table("data/vol_test/dada2_out/ASVs_taxonomy_bootstrapVals.tsv", sep = "\t", header = T)
head(asv_boot)
colnames(asv_boot)[1] = "ASV"

#read metadata
source("R_scripts/vol_test/read_data_and_split_objects.R")
asv_counts = asv_tab %>% 
    rowSums
head(asv_counts)
length(asv_counts)
sum(asv_counts == 0)

length(asv_counts.wash)

asv_boot.counts = left_join(
    data.frame(ASV = names(asv_counts), counts = asv_counts),
    asv_boot,
    by = "ASV"
)
head(asv_boot.counts)

asv_boot.long = asv_boot.counts %>% 
    pivot_longer(
        cols = c(-ASV, -counts), 
        names_to = "tax_level", 
        values_to = "bootstrap"
    ) %>% group_by(tax_level, bootstrap) %>%
    summarize(counts = sum(counts))
head(asv_boot.long)



asv_boot.long.group_counts = within(asv_boot.long, {
       tax_level_seqs = ave(counts, tax_level,  FUN = sum)
    }
)
asv_boot.long.group_counts %>%
    filter(tax_level == "Species") %>%
    pull(counts) %>% sum()
asv_boot.long.group_counts %>%
    filter(tax_level == "Species" & bootstrap > 90) %>%
    pull(counts) %>% sum()
2333757/4936857
#0.473
asv_boot.long.group_counts$tax_level_seqs %>% unique

####################################
#Plots

#volume
p1 = ggplot(
    asv_boot.long.group_counts, 
    aes(x = bootstrap, y = counts/tax_level_seqs*100)
) + 
    #scale_x_binned(n.breaks = 20) +
    scale_x_binned(
        breaks = seq(0,100,by=5),
        labels = c("","5","","","","25","","","","","50","","","","","75","","","","95","")
    ) +
    geom_col() +
    facet_wrap(
        ~factor(tax_level, 
                levels = c(
                    "Kingdom", "Phylum", "Class", "Order", 
                    "Family", "Genus", "Species"
                )
        ),
        nrow = 1
    ) +
    scale_y_continuous(limits = c(0,100)) +
    labs(x = "Bootstrap score", y = "Sequences assigned (%)", title = "5% bins") +
    my_gg_theme +
    theme(
        axis.text = element_text(size = 14)
    )
p1

p2 = ggplot(
    asv_boot.long.group_counts, 
    aes(x = bootstrap, y = counts/tax_level_seqs*100)
) + 
    #scale_x_binned(n.breaks = 20) +
    scale_x_binned(
        breaks = seq(0,100,by=10),
        labels = c("","10","","30","","50","","70","","90","")
    ) +
    geom_col() +
    facet_wrap(
        ~factor(tax_level, 
                levels = c(
                    "Kingdom", "Phylum", "Class", "Order", 
                    "Family", "Genus", "Species"
                )
        ),
        nrow = 1
    ) +
    scale_y_continuous(limits = c(0,100)) +
    labs(x = "Bootstrap score", y = "Sequences assigned (%)", title = "10% bins") +
    my_gg_theme +
    theme(
        axis.text = element_text(size = 14)
    )
p2

p3 = ggplot(
    asv_boot.long.group_counts, 
    aes(x = bootstrap, y = counts/tax_level_seqs*100)
) + 
    #scale_x_binned(n.breaks = 20) +
    scale_x_binned(
        breaks = seq(0,100,by=20), 
        labels = c("",20,40,60,80,"")
    ) +
    geom_col() +
    facet_wrap(
        ~factor(tax_level, 
                levels = c(
                    "Kingdom", "Phylum", "Class", "Order", 
                    "Family", "Genus", "Species"
                )
        ),
        nrow = 1
    ) +
    scale_y_continuous(limits = c(0,100)) +
    labs(x = "Bootstrap score", y = "Sequences assigned (%)", title = "20% bins") +
    my_gg_theme +
    theme(
        axis.text = element_text(size = 14)
    )
p3

pdf("figures/vol_test/all_sequences.RDP-NBC_boostrap_by_tax_level.pdf", width = 12, height = 6)
p1
p2
p3
dev.off()

