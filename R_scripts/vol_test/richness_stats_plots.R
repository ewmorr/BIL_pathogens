library(dplyr)
library(ggplot2)
library(RColorBrewer)
source("R_scripts/library.R")

#read data
vol.div = read.csv("data/vol_test/vol.diversity.csv")
wash.div = read.csv("data/vol_test/wash.diversity.csv")
head(vol.div)

#pull in metadata etc
source("R_scripts/vol_test/read_data_and_split_objects.R")
head(metadata.vol)
head(metadata.wash)

###
#Because we excluded some samples due to low sequence counts we will calculate
# richness/diversity based on the unrarefied samples (i.e., using all the sequences)

######### This is no longer necessary because we now calculated with and
# without excluding the minimums
# But we have converted to using the no_min tables (gives the same)
low_count_samples = all_vol[is.na(all_vol$richness) , "sequenceID"]
low_count_samples.div = data.frame(
    sequenceID = richness_calc(asv_tab.vol[low_count_samples, ]) %>% names,
    shannon = vegan::diversity(asv_tab.vol[low_count_samples, ], "shannon"),
    simpson = vegan::diversity(asv_tab.vol[low_count_samples, ], "simpson"),
    richness = richness_calc(asv_tab.vol[low_count_samples, ])
)
#join to the div table
vol.div.all = rbind(vol.div, low_count_samples.div)
###

#Join with metadata
all_vol = left_join(metadata.vol, vol.div.all, by = "sequenceID")
head(all_vol)
sum(all_vol$richness %>% is.na)

all_wash = left_join(metadata.wash, wash.div, by = "sequenceID")
head(all_wash)

###############
###############
#Wash exp 

#stats
wash_mod1 = aov(richness ~ volumeOrWash , data = all_wash)
plot(residuals(wash_mod1))
qqnorm(residuals(wash_mod1))

summary(wash_mod1)
TukeyHSD(wash_mod1)

wash_mod2 = aov(richness ~ volumeOrWash + collectionID, data = all_wash)
plot(residuals(wash_mod2))
qqnorm(residuals(wash_mod2))

summary(wash_mod2)
TukeyHSD(wash_mod2)

anova(wash_mod1, wash_mod2)

#plots

p1 = ggplot(all_wash, 
    aes(
        x = factor(volumeOrWash, levels = c("pre", "mid", "post")), 
        y = richness
    )
) +
    geom_point(aes(color = collectionID), size = 3) +
    geom_smooth(method = "lm") +
    scale_color_manual(values = cbPalette) +
    my_gg_theme +
    labs(x = "Wash step", y = "ASV richness")
p1

p2 = ggplot(
    all_wash %>% 
           group_by(volumeOrWash) %>%
           summarize(
               richness.mean = mean(richness), 
               richness.se = sd(richness)/sqrt(n()) 
           ),
    aes(
        x = factor(volumeOrWash, levels = c("pre", "mid", "post")), 
        y = richness.mean
    )
) +
    geom_point(size = 2) +
    geom_errorbar(
        aes(
            ymin = richness.mean - richness.se, 
            ymax = richness.mean + richness.se
        ),
        width = 0.1
    ) +
    my_gg_theme +
    labs(x = "Wash step", y = "ASV richness") +
    scale_y_continuous(limits = c(0, 150))
p2



pdf("figures/vol_test/wash_richness.pdf")
p1
p2
dev.off()

################


###############
###############
#Vol exp stats

all_vol %>% filter(volumeOrWash == 0)

vol_mod1 = lm(richness ~ volumeOrWash, data = all_vol)
vol_mod2 = lm(richness ~ volumeOrWash + I(volumeOrWash^2), data = all_vol)
vol_mod3 = lm(richness ~ volumeOrWash + I(volumeOrWash^2) + I(volumeOrWash^3), data = all_vol)

anova(vol_mod1, vol_mod2)
anova(vol_mod2, vol_mod3)

plot(residuals(vol_mod2))
qqnorm(residuals(vol_mod2))

summary(vol_mod2)

p1 = ggplot(
    all_vol, 
    aes(x = volumeOrWash, y = richness)
) +
    geom_line() +
    geom_point(aes(color = seq_cts), size = 2) +
    scale_color_gradient(trans = "log10") +
#    geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = F) +
#    geom_smooth(se = F, span = 1.35) +
    facet_wrap(~collectionID, ncol = 4) +
    my_gg_theme +
    labs(
        x = "Volume (ml) (0 equals blank)",
        y = "ASV richness",
        color = "Total seqs"
    ) +
    theme(
        panel.grid = element_blank(),
        legend.title = element_text(size = 17),
        legend.text = element_text(size = 16)
    )
p1

p2 = ggplot(
    all_vol, 
    aes(x = volumeOrWash, y = simpson)
) +
    geom_line() +
    geom_point(aes(color = seq_cts)) +
    scale_color_gradient(trans = "log10") +
    #    geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = F) +
    #    geom_smooth(se = F, span = 1.35) +
    facet_wrap(~collectionID, ncol = 4) +
    my_gg_theme +
    labs(
        x = "Volume (ml) (0 equals blank)",
        y = "ASV richness",
        color = "Total seqs"
    ) +
    theme(
        panel.grid = element_blank(),
        legend.title = element_text(size = 18)
    )
p2
p1

pdf("figures/vol_test/volume_richness.pdf", width = 8, height = 6)
p1
dev.off()

