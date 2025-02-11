library(dplyr)
library(vegan)
library(ggplot2)
library(lubridate)
library(gridExtra)
library(stringr)
source("library/library.R")

metadata = read.csv("data/2024_metadata/collated/collated_metadata_NY-OH_11062024.csv")
head(metadata)

bl_v_samps = read.csv("data/FEDRR_11062024/processed_tables/avg_div/diversity.samps_and_baseline.csv")
colnames(bl_v_samps)[1] = "SequenceID"

bl_v_samps.meta = left_join(bl_v_samps, metadata)
bl_aov = aov(richness ~ BaselineType, data = bl_v_samps.meta)
summary(bl_aov)

bl_v_samps.meta[bl_v_samps.meta$BaselineType == "BLO", "date"] = "BL"
uniq_dates = bl_v_samps.meta$date %>% unique() %>% str_sort(., numeric =T)

p1 = ggplot(bl_v_samps.meta %>% filter(State == "NY"), 
       aes(
           x = factor(date, levels = c(uniq_dates[c(20,1:19)])),
           y = richness
       )
) +
    geom_boxplot() +
    geom_point(aes(color = Lure)) +
    facet_wrap(~Site, scales = "free_x", labeller = label_wrap_gen(17)) +
    scale_color_brewer(palette = "Set2") +
    labs(x = "Sample date/type", y = "ASV richness") +
    my_gg_theme +
    theme(
        legend.position = "bottom",
        axis.text.x = element_text(angle = 55, hjust = 1)
    )
p1

bl_v_samps.meta$Site = sub(" Coridor-12", "", bl_v_samps.meta$Site)
p2 = ggplot(bl_v_samps.meta %>% filter(State == "OH"), 
            aes(
                x = factor(date, levels = c(uniq_dates[c(20,1:19)])),
                y = richness
            )
) +
    geom_boxplot() +
    geom_point(aes(color = Lure)) +
    facet_wrap(~Site, scales = "free_x", labeller = label_wrap_gen(17)) +
    scale_color_brewer(palette = "Set2") +
    labs(x = "Sample date/type", y = "ASV richness") +
    my_gg_theme +
    theme(
        legend.position = "bottom",
        axis.text.x = element_text(angle = 55, hjust = 1)
    )
p2

pdf("figures/FEDRR_11062024/richness.NY-OH.pdf", width = 20, height = 8)
grid.arrange(p1,p2,ncol = 2)
dev.off()

p1 = ggplot(bl_v_samps.meta %>% filter(State == "NY"), 
            aes(
                x = factor(date, levels = c(uniq_dates[c(20,1:19)])),
                y = 1-simpson
            )
) +
    geom_boxplot() +
    geom_point(aes(color = Lure)) +
    facet_wrap(~Site, scales = "free_x", labeller = label_wrap_gen(17)) +
    scale_color_brewer(palette = "Set2") +
    labs(x = "Sample date/type", y = "Simpson D") +
    my_gg_theme +
    theme(
        legend.position = "bottom",
        axis.text.x = element_text(angle = 55, hjust = 1)
    )
p1

#bl_v_samps.meta$Site = sub(" Coridor-12", "", bl_v_samps.meta$Site)
p2 = ggplot(bl_v_samps.meta %>% filter(State == "OH"), 
            aes(
                x = factor(date, levels = c(uniq_dates[c(20,1:19)])),
                y = 1-simpson
            )
) +
    geom_boxplot() +
    geom_point(aes(color = Lure)) +
    facet_wrap(~Site, scales = "free_x", labeller = label_wrap_gen(17)) +
    scale_color_brewer(palette = "Set2") +
    labs(x = "Sample date/type", y = "Simpson D") +
    my_gg_theme +
    theme(
        legend.position = "bottom",
        axis.text.x = element_text(angle = 55, hjust = 1)
    )
p2

pdf("figures/FEDRR_11062024/simpson.NY-OH.pdf", width = 20, height = 8)
grid.arrange(p1,p2,ncol = 2)
dev.off()


p1 = ggplot(bl_v_samps.meta %>% filter(State == "NY"), 
            aes(
                x = factor(date, levels = c(uniq_dates[c(20,1:19)])),
                y = shannon
            )
) +
    geom_boxplot() +
    geom_point(aes(color = Lure)) +
    facet_wrap(~Site, scales = "free_x", labeller = label_wrap_gen(17)) +
    scale_color_brewer(palette = "Set2") +
    labs(x = "Sample date/type", y = "Shannon H'") +
    my_gg_theme +
    theme(
        legend.position = "bottom",
        axis.text.x = element_text(angle = 55, hjust = 1)
    )
p1

#bl_v_samps.meta$Site = sub(" Coridor-12", "", bl_v_samps.meta$Site)
p2 = ggplot(bl_v_samps.meta %>% filter(State == "OH"), 
            aes(
                x = factor(date, levels = c(uniq_dates[c(20,1:19)])),
                y = shannon
            )
) +
    geom_boxplot() +
    geom_point(aes(color = Lure)) +
    facet_wrap(~Site, scales = "free_x", labeller = label_wrap_gen(17)) +
    scale_color_brewer(palette = "Set2") +
    labs(x = "Sample date/type", y = "Shannon H'") +
    my_gg_theme +
    theme(
        legend.position = "bottom",
        axis.text.x = element_text(angle = 55, hjust = 1)
    )
p2

pdf("figures/FEDRR_11062024/shannon.NY-OH.pdf", width = 20, height = 8)
grid.arrange(p1,p2,ncol = 2)
dev.off()
