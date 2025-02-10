library(dplyr)
library(vegan)
library(ggplot2)
library(lubridate)
library(gridExtra)
source("library/library.R")

metadata = read.csv("data/2024_metadata/collated/collated_metadata_NY-OH_11062024.csv")
head(metadata)

#full rarefied table
#
bray_bin.dist = readRDS("data/FEDRR_11062024/processed_tables/avg_dist/bray-binary.full_tab_min450.rds")
bray_log.dist = readRDS("data/FEDRR_11062024/processed_tables/avg_dist/bray-logCounts.full_tab_min450.rds")

set.seed(12345)
##################################
#perform NMDS
bray_bin.nmds = metaMDS(bray_bin.dist, try = 20, trymax = 100, k = 2, weakties = T)
#bray_bin.nmds = metaMDS(bray_bin.dist, previous.best = bray_bin.nmds)
bray_bin.nmds

bray_log.nmds = metaMDS(bray_log.dist, try = 20, trymax = 100, k = 2)
#bray_log.nmds = metaMDS(bray_log.dist, previous.best = bray_log.nmds)
bray_log.nmds


nmds.metadata = left_join(
    data.frame(
        SequenceID = rownames(scores(bray_bin.nmds)),
        bray_bin.NMDS1 = scores(bray_bin.nmds)[,1],
        bray_bin.NMDS2 = scores(bray_bin.nmds)[,2]#,
#        bray_bin.NMDS3 = scores(bray_bin.nmds)[,3]
    ),
    data.frame(
        SequenceID = rownames(scores(bray_bin.nmds)),
        bray_log.NMDS1 = scores(bray_log.nmds)[,1],
        bray_log.NMDS2 = scores(bray_log.nmds)[,2]
    ),
    by = "SequenceID"
) %>% left_join(
    .,
    metadata,
    by = "SequenceID"
)
str(nmds.metadata)


nmds.metadata$Lure %>% unique
nmds.metadata$trapID %>% unique
nmds.metadata %>% filter(State == "NY") %>% pull(trapID) %>% unique
nmds.metadata$BaselineType %>% unique()

nmds.metadata$mdy = mdy(nmds.metadata$date)


#note if wish to set breaks manually can use date_breaks function
#e.g., scale_fill_gradientn(colours = rev(five_cols_gradient_palette), breaks = date_breaks)
#see lib and vol_nmds.r script
#
date_breaks(as.integer(nmds.metadata %>% filter(!is.na(mdy)) %>% pull(mdy)))
five_cols_gradient_palette = c('#ca0020','#f4a582','#f7f7f7','#92c5de','#0571b0')
five_cols_gradient_palette = c('#ca0020','#ef8a62','#f7f7f7','#67a9cf','#2166ac')
five_cols_gradient_palette = c('#b2182b','#d6604d','#f4a582','white','#2166ac')

nmds.metadata$state.type = paste(nmds.metadata$State, nmds.metadata$BaselineType, sep = " ")

nmds.metadata[nmds.metadata$BaselineType == "BLO", "mdy"] = NA

#fill by date
p1 = ggplot(nmds.metadata %>% filter(BaselineType != "PCRnegative"), 
       aes(x = bray_bin.NMDS1, y = bray_bin.NMDS2, shape = state.type, fill = as.integer(mdy))
) +
    geom_point(size = 3) +
    scale_shape_manual(values = c(25,24,23,22)) +
    scale_fill_gradientn(colours = rev(five_cols_gradient_palette), breaks = date_breaks, na.value = "#525252") +
    guides(shape = guide_legend(order = 1)) +
    labs(x="NMDS1", y = "NMDS2") +
    my_gg_theme +
    theme(
        legend.position = "inside",
        legend.position.inside = c(0.125,0.78) 
    )
p1

#baseline v samples
p2 = ggplot(nmds.metadata %>% filter(BaselineType != "PCRnegative"), 
       aes(x = bray_bin.NMDS1, y = bray_bin.NMDS2, fill = state.type, shape = BaselineType)
) +
    geom_point(size = 3) +
    scale_shape_manual(values = c(25,21)) +
    scale_fill_brewer(palette = "Paired") +
    labs(x="NMDS1", y = "NMDS2") +
    #ylim(-0.4, NA) +
    annotate(geom = "text", label = "stress = 0.217", x = min(nmds.metadata$bray_bin.NMDS1), y = -0.4, hjust = 0, size = 5.5) +
    guides(shape = "none", fill = guide_legend(override.aes = list(shape = c(25,21,25,21)))) +
    my_gg_theme +
    theme(
        legend.position = "inside",
        legend.position.inside = c(0.125,0.9) 
    )
p2

pdf("figures/FEDRR_11062024/baseline_comparisons.NY-OH.pdf", width = 16, height = 7)
grid.arrange(p2,p1,ncol = 2)
dev.off()

nmds.metadata$trapID %>% unique() %>% length()
nmds.metadata %>% filter(State == "OH") %>% pull(trapID) %>% unique %>% length

p3 = ggplot(nmds.metadata %>% filter(BaselineType != "PCRnegative" & State == "OH"), 
            aes(x = bray_bin.NMDS1, y = bray_bin.NMDS2, shape = state.type, fill = trapID)
) +
    geom_point(size = 3) +
    scale_shape_manual(values = c(25,21)) +
    scale_fill_manual(values = c25) +
    guides(shape = guide_legend(order = 1), fill = guide_legend(override.aes = list(shape = 21)) ) +
    labs(x="NMDS1", y = "NMDS2") +
    my_gg_theme
p3

p4 = ggplot(nmds.metadata %>% filter(BaselineType != "PCRnegative" & State == "OH"), 
            aes(x = bray_bin.NMDS1, y = bray_bin.NMDS2, shape = state.type, fill = Site)
) +
    geom_point(size = 3) +
    scale_shape_manual(values = c(25,21)) +
    scale_fill_brewer(palette = "Set1") +
    guides(shape = guide_legend(order = 1), fill = guide_legend(override.aes = list(shape = 21)) ) +
    labs(x="NMDS1", y = "NMDS2") +
    my_gg_theme
p4

pdf("figures/FEDRR_11062024/baseline_comparisons.OH.trapID-site.pdf", width = 20, height = 7)
grid.arrange(p4,p3,ncol = 2,widths = c(0.55, 0.45))
dev.off()

pdf("figures/FEDRR_11062024/baseline_comparisons.OH.site.pdf", width = 11, height = 7)
p4
dev.off()


