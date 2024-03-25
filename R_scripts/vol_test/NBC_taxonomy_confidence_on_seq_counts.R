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
asv_counts.vol = asv_tab.vol %>% 
    t() %>%
    rowSums
head(asv_counts.vol)
length(asv_counts.vol)
sum(asv_counts.vol == 0)

asv_counts.wash = asv_tab.wash %>% 
    t() %>%
    rowSums
head(asv_counts.wash)
length(asv_counts.wash)

asv_boot.vol = left_join(
    data.frame(ASV = names(asv_counts.vol), counts = asv_counts.vol),
    asv_boot,
    by = "ASV"
)
head(asv_boot.vol)
asv_boot.wash = left_join(
    data.frame(ASV = names(asv_counts.wash), counts = asv_counts.wash),
    asv_boot,
    by = "ASV"
)
asv_boot.vol.long = asv_boot.vol %>% 
    pivot_longer(
        cols = c(-ASV, -counts), 
        names_to = "tax_level", 
        values_to = "bootstrap"
    ) %>% group_by(tax_level, bootstrap) %>%
    summarize(counts = sum(counts))
head(asv_boot.vol.long)

asv_boot.wash.long = asv_boot.wash %>% 
    pivot_longer(
        cols = c(-ASV, -counts), 
        names_to = "tax_level", 
        values_to = "bootstrap"
    ) %>% group_by(tax_level, bootstrap) %>%
    summarize(counts = sum(counts))

asv_boot.vol.long.group_counts = within(asv_boot.vol.long, {
       tax_level_seqs = ave(counts, tax_level,  FUN = sum)
    }
)
asv_boot.vol.long.group_counts %>%
    filter(tax_level == "Species") %>%
    pull(counts) %>% sum()
asv_boot.vol.long.group_counts %>%
    filter(tax_level == "Species" & bootstrap > 90) %>%
    pull(counts) %>% sum()
2193344/4604148
#0.476

ggplot(asv_boot.vol.long.group_counts, aes(x = bootstrap, y = counts/tax_level_seqs)) + 
    scale_x_binned(n.breaks = 20) +
    geom_col() +
    facet_wrap(
        ~factor(tax_level, 
            levels = c(
                "Kingdom", "Phylum", "Class", "Order", 
                    "Family", "Genus", "Species"
            )
        ),
        nrow = 1
    )
    

    labs(x = "Bootstrap score", y = "Count (%)", title = "5% bins") +
    scale_y_continuous(
        labels = scales::percent#, 
        #limits = c(0,0.11),
        #breaks = c(0,0.05,0.1)
    ) +
    my_gg_theme +
    theme(
        axis.text = element_text(size = 14)
    )
p1
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
