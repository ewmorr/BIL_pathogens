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



p1 = ggplot(
    asv_boot %>% 
        pivot_longer(
            cols = c(-ASV), 
            names_to = "tax_level", 
            values_to = "bootstrap"
        ), 
    aes(x = bootstrap)
) + 
    geom_histogram(aes(y = after_stat(count/sum(count))), bins = 20 ) +
    facet_wrap(
        ~factor(tax_level, 
                levels = c(
                    "Kingdom", "Phylum", "Class", "Order", 
                    "Family", "Genus", "Species"
            )
        ),
        nrow = 1
    ) +
    labs(x = "Bootstrap score", y = "Count (%)", title = "5% bins") +
    scale_y_continuous(
        labels = scales::percent, 
        limits = c(0,0.11),
        breaks = c(0,0.05,0.1)
    ) +
    my_gg_theme +
    theme(
        axis.text = element_text(size = 14)
    )

p2 = ggplot(
    asv_boot %>% 
        pivot_longer(
            cols = c(-ASV), 
            names_to = "tax_level", 
            values_to = "bootstrap"
        ), 
    aes(x = bootstrap)
) + 
    geom_histogram(aes(y = after_stat(count/sum(count))), bins = 10 ) +
    #geom_histogram() +
    facet_wrap(
        ~factor(tax_level, 
                levels = c(
                    "Kingdom", "Phylum", "Class", "Order", 
                    "Family", "Genus", "Species"
                )
        ),
        nrow = 1
    ) +
    labs(x = "Bootstrap score", y = "Count (%)", title = "10% bins") +
    scale_y_continuous(
        labels = scales::percent, 
        limits = c(0,0.11),
        breaks = c(0,0.05,0.1)
    ) +
    my_gg_theme +
    theme(
        axis.text = element_text(size = 14)
    )

p3 = ggplot(
    asv_boot %>% 
        pivot_longer(
            cols = c(-ASV), 
            names_to = "tax_level", 
            values_to = "bootstrap"
        ), 
    aes(x = bootstrap)
) + 
    geom_histogram(aes(y = after_stat(count/sum(count))), bins = 5 ) +
    facet_wrap(
        ~factor(tax_level, 
                levels = c(
                    "Kingdom", "Phylum", "Class", "Order", 
                    "Family", "Genus", "Species"
                )
        ),
        nrow = 1
    ) +
    labs(x = "Bootstrap score", y = "Count (%)", title = "20% bins") +
    scale_y_continuous(
        labels = scales::percent, 
        limits = c(0,0.11),
        breaks = c(0,0.05,0.1)
    ) +
    scale_x_continuous(breaks = c(0,25,50,75,100)) +
    my_gg_theme +
    theme(
        axis.text = element_text(size = 14)
    )
p1
p2
p3

pdf("figures/vol_test/RDP-NBC_boostrap_by_tax_level.pdf", width = 12, height = 6)
p1
p2
p3
dev.off()
