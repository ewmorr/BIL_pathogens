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
#low_count_samples = all_vol[is.na(all_vol$richness) , "sequenceID"]
#low_count_samples.div = data.frame(
#    sequenceID = richness_calc(asv_tab.vol[low_count_samples, ]) %>% names,
#    shannon = vegan::diversity(asv_tab.vol[low_count_samples, ], "shannon"),
#    simpson = vegan::diversity(asv_tab.vol[low_count_samples, ], "simpson"),
#    richness = richness_calc(asv_tab.vol[low_count_samples, ])
#)
#join to the div table
#vol.div.all = rbind(vol.div, low_count_samples.div)
###

#Join with metadata
all_vol = left_join(metadata.vol, vol.div, by = "sequenceID")
head(all_vol)
sum(all_vol$richness %>% is.na)
all_vol %>% filter(volumeOrWash == 0)
all_vol.no_blanks = all_vol %>%
    filter(volumeOrWash != 0 & seq_cts >= 1900)
nrow(all_vol)
nrow(all_vol.no_blanks)

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
# p = 0.0007
TukeyHSD(wash_mod1)
#              diff        lwr       upr     p adj
#post-mid -92.97763 -142.29250 -43.66277 0.0028226
#pre-mid   24.55977  -24.75510  73.87463 0.3439894
#pre-post 117.53740   68.22253 166.85227 0.0008154


wash_mod2 = aov(richness ~ volumeOrWash + collectionID, data = all_wash)
plot(residuals(wash_mod2))
qqnorm(residuals(wash_mod2))

summary(wash_mod2)
#             Df Sum Sq Mean Sq F value  Pr(>F)   
#volumeOrWash  2  23063   11532  27.023 0.00475 **
#collectionID  2    618     309   0.724 0.53903   
#Residuals     4   1707     427                   

TukeyHSD(wash_mod2)
#same as above

anova(wash_mod1, wash_mod2)
# p = 0.539

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

p3 = all_wash %>% 
    group_by(volumeOrWash) %>%
    summarise(mean = mean(richness), SE = sd(richness)/n()) %>%
    ggplot(., 
            aes(
                x = factor(volumeOrWash, levels = c("pre", "mid", "post")), 
                y = mean
            )
) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = mean - SE, ymax = mean + SE), width = 0.1) +
    my_gg_theme +
    labs(x = "Wash step", y = "ASV richness")
p3


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
#p = 0.13
anova(vol_mod2, vol_mod3)
#p = 0.34

plot(residuals(vol_mod2))
qqnorm(residuals(vol_mod2))

summary(vol_mod2)

p1 = ggplot(
    all_vol.no_blanks, 
    aes(x = volumeOrWash, y = richness)
) +
    geom_line() +
    geom_point(aes(color = seq_cts), size = 2) +
    scale_color_gradient(trans = "log10") +
    scale_x_continuous(breaks = c(5,20,40,60,80), limits = c(0,90)) +
#    geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = F) +
#    geom_smooth(se = F, span = 1.35) +
    facet_wrap(~collectionID, ncol = 4) +
    my_gg_theme +
    labs(
        x = "Volume (ml)",
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
    all_vol.no_blanks, 
    aes(x = volumeOrWash, y = 1-simpson)
) +
    geom_line() +
    geom_point(aes(color = seq_cts)) +
    scale_color_gradient(trans = "log10") +
    scale_x_continuous(breaks = c(5,20,40,60,80), limits = c(0,90)) +
    #    geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = F) +
    #    geom_smooth(se = F, span = 1.35) +
    facet_wrap(~collectionID, ncol = 4) +
    my_gg_theme +
    labs(
        x = "Volume (ml)",
        y = "Simpson dominance (D)",
        color = "Total seqs"
    ) +
    theme(
        panel.grid = element_blank(),
        legend.title = element_text(size = 18)
    )
p2

p3 = ggplot(
    all_vol.no_blanks, 
    aes(x = volumeOrWash, y = shannon)
) +
    geom_line() +
    geom_point(aes(color = seq_cts)) +
    scale_color_gradient(trans = "log10") +
    scale_x_continuous(breaks = c(5,20,40,60,80), limits = c(0,90)) +
    #    geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = F) +
    #    geom_smooth(se = F, span = 1.35) +
    facet_wrap(~collectionID, ncol = 4) +
    my_gg_theme +
    labs(
        x = "Volume (ml)",
        y = "Shannon diversity (H')",
        color = "Total seqs"
    ) +
    theme(
        panel.grid = element_blank(),
        legend.title = element_text(size = 18)
    )
p3

pdf("figures/vol_test/volume_richness.pdf", width = 10, height = 6)
p1
p2
p3
dev.off()

