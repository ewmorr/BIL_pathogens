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

#asv_boot$n_taxa = nrow(asv_boot)
head(asv_boot)
is.na(asv_boot) %>% sum



asv_boot.conf_by_level = asv_boot %>% 
    pivot_longer(
        cols = c(-ASV), 
        names_to = "tax_level", 
        values_to = "bootstrap"
    ) %>% 
    group_by(tax_level, bootstrap) %>%
    summarize(n_boot = n())

n_taxa = nrow(asv_boot)
asv_boot.conf_by_level$perc_boot = asv_boot.conf_by_level$n_boot/n_taxa

asv_boot.conf_by_level %>%
    filter(tax_level == "Species") %>%
    pull(perc_boot) %>%
    sum

p1 = ggplot(
    asv_boot.conf_by_level, 
    aes(x = bootstrap, y = perc_boot*100)
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
    labs(x = "Bootstrap score", y = "ASVs assigned (%)", title = "5% bins") +
    my_gg_theme +
    theme(
        axis.text = element_text(size = 14)
    )
p1

p2 = ggplot(
    asv_boot.conf_by_level, 
    aes(x = bootstrap, y = perc_boot*100)
) + 
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
    labs(x = "Bootstrap score", y = "ASVs assigned (%)", title = "10% bins") +
    my_gg_theme +
    theme(
        axis.text = element_text(size = 14)
    )

p3 = ggplot(
    asv_boot.conf_by_level, 
    aes(x = bootstrap, y = perc_boot*100)
) + 
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
    labs(x = "Bootstrap score", y = "ASVs assigned (%)", title = "20% bins") +
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
