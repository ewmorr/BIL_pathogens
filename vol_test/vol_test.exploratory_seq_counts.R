library(dplyr)
library(ggplot2)
source("R_scripts/library.R")

theme_bw_noGrid = theme_bw() +
    theme(
        panel.grid = element_blank()
    )

#read the data
source("R_scripts/vol_test/read_data_and_split_objects.R")

#plot sequence counts by volume
p1 = ggplot(
    metadata.vol, 
    aes(x = volumeOrWash, y = seq_cts)
) +
    geom_point() +
    geom_smooth(method = "lm", se = F) +
    facet_wrap(~collectionID) +
    scale_y_log10() +
    my_gg_theme +
    labs(
        x = "Volume (ml) (0 equals blank)",
        y = "Total sequences"
) +
    theme(
        panel.grid = element_blank()
    )
p1

pdf("figures/vol_test/seq_counts_over_volume.pdf", width = 8, height = 6)
p1
dev.off()


#plot sequence counts by wash stage
p2 = ggplot(
    metadata.wash, 
    aes(x = factor(volumeOrWash, levels = c("pre", "mid", "post")), y = seq_cts)
) +
    geom_point() +
    geom_smooth(method = "lm", se = F) +
    facet_wrap(~collectionID) +
    scale_y_log10() +
    my_gg_theme +
    labs(
        x = "Wash stage",
        y = "Total sequences"
    )

p3 = ggplot(
    metadata.wash %>%
        group_by(volumeOrWash) %>%
        summarize(
            seq_cts.mean = mean(seq_cts),
            seq_cts.se = sd(seq_cts)/sqrt(n())
        ),
    aes(
        x = factor(volumeOrWash, levels = c("pre", "mid", "post")), 
        y = seq_cts.mean
    )
) +
    geom_point(size = 2) +
    geom_errorbar(
        aes(
            ymin = seq_cts.mean - seq_cts.se, 
            ymax = seq_cts.mean + seq_cts.se
        ),
        width = 0.1
    ) +
    my_gg_theme +
    labs(x = "Wash step", y = "Total sequences") 
p3

p4 = ggplot(metadata.wash, 
            aes(
                x = factor(volumeOrWash, levels = c("pre", "mid", "post")), 
                y = seq_cts
            )
) +
    geom_point(aes(color = collectionID), size = 3) +
    geom_smooth(method = "lm") +
    scale_color_manual(values = cbPalette) +
    scale_y_log10(labels = fancy_scientific) +
    my_gg_theme +
    labs(x = "Wash step", y = "Total sequences")
p4


pdf("figures/vol_test/seq_counts_over_wash_stage.pdf", width = 8, height = 4)
p2
p4
p3
dev.off()
