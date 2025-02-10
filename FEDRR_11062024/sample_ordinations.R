library(dplyr)
library(vegan)
library(ggplot2)
library(lubridate)
library(gridExtra)
library(gtable)
source("library/library.R")
set.seed(12345)

metadata = read.csv("data/2024_metadata/collated/collated_metadata_NY-OH_11062024.csv")
head(metadata)

#full rarefied table
#
bray_bin.dist = readRDS("data/FEDRR_11062024/processed_tables/avg_dist/bray-binary.samps.rds")
bray_log.dist = readRDS("data/FEDRR_11062024/processed_tables/avg_dist/bray-logCounts.samps.rds")


##################################
#perform NMDS
bray_bin.nmds = metaMDS(bray_bin.dist, try = 20, trymax = 100)
#bray_bin.nmds = metaMDS(bray_bin.dist, previous.best = bray_bin.nmds)
bray_bin.nmds
#stress = 0.163

bray_log.nmds = metaMDS(bray_log.dist, try = 20, trymax = 100)
#bray_log.nmds = metaMDS(bray_log.dist, previous.best = bray_log.nmds)
bray_log.nmds
#stress = 0.1998

nmds.metadata = left_join(
    data.frame(
        SequenceID = rownames(scores(bray_bin.nmds)),
        bray_bin.NMDS1 = scores(bray_bin.nmds)[,1],
        bray_bin.NMDS2 = scores(bray_bin.nmds)[,2]
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
five_cols_gradient_palette = c('#b2182b','#d6604d','#f4a582','white','#2166ac')

five_cols_gradient_palette = c('#ca0020','#ef8a62','#f7f7f7','#67a9cf','#2166ac')




#fill by date shape is state
p1 = ggplot(nmds.metadata, 
       aes(x = bray_bin.NMDS1, y = bray_bin.NMDS2, shape = State, fill = as.integer(mdy))
) +
    geom_point(size = 3) +
    scale_shape_manual(values = c(21,24)) +
    scale_fill_gradientn(colours = rev(five_cols_gradient_palette), breaks = date_breaks) +
    annotate(geom = "text", label = "stress = 0.163", x = min(nmds.metadata$bray_bin.NMDS1), y = min(nmds.metadata$bray_bin.NMDS2), hjust = 0, size = 5.5) +
    guides(shape = guide_legend(order = 1), fill = guide_colorbar(order = 2)) +
    labs(x="NMDS1", y = "NMDS2") +
    my_gg_theme +
    #annotate(geom = "segment", x = -0.275, xend = -0.275, y = 0.15, yend = Inf) +
    #annotate(geom = "segment", x = -Inf, xend = -0.275, y = 0.15, yend = 0.15) +
    theme(
        legend.position = "inside",
        legend.position.inside = c(0.075, 0.825)
    )
p1
pdf("figures/FEDRR_11062024/date.NY-OH.pdf", width = 8, height = 7)
p1
dev.off()

 
#separate the states
p7 = ggplot(nmds.metadata %>% filter(State == "NY"), 
            aes(x = bray_bin.NMDS1, y = bray_bin.NMDS2, color = as.integer(mdy), shape = Site)
) +
    geom_point(size = 3) +
    scale_color_gradientn(colours = rev(five_cols_gradient_palette), breaks = date_breaks, guide = guide_colorbar(title = element_blank(), order = 2)) +
    labs(x="NMDS1", y = "NMDS2", title = "New York") +
    my_gg_theme +
    theme(
        legend.title = element_text(size = 16)
    )
p7
p8 = ggplot(nmds.metadata %>% filter(State == "OH"), 
            aes(x = bray_bin.NMDS1, y = bray_bin.NMDS2, color = as.integer(mdy), shape = Site)
) +
    geom_point(size = 3) +
    scale_color_gradientn(colours = rev(five_cols_gradient_palette), breaks = date_breaks, guide = guide_colorbar(title = element_blank(), order = 2)) +
    labs(x="NMDS1", y = "NMDS2", color = "Risk level", title = "Ohio") +
    my_gg_theme +
    theme(
        legend.title = element_text(size = 16)
    )
p8

pdf("figures/FEDRR_11062024/date.NY-OH_facet.pdf", width = 19, height = 6)
grid.arrange(p7,p8,ncol = 2, widths = c(0.48,0.52))
dev.off()

####################################################
####################################################
#Site risk
#

p3 = ggplot(nmds.metadata, 
            aes(x = bray_bin.NMDS1, y = bray_bin.NMDS2, fill = Risk, shape = State)
) +
    geom_point(size = 3) +
    scale_shape_manual(values = c(21,24)) +
    annotate(geom = "text", label = "stress = 0.163", x = min(nmds.metadata$bray_bin.NMDS1), y = min(nmds.metadata$bray_bin.NMDS2), hjust = 0, size = 5.5) +
    scale_fill_brewer(palette = "Dark2") +
    guides(
        shape = guide_legend(order = 2, title = NULL), 
        fill = guide_legend(order = 1, override.aes = list(shape = 21))
    ) +
    labs(x="NMDS1", y = "NMDS2", fill = "Risk level") +
    my_gg_theme +
    theme(
        legend.position = "inside",
        legend.position.inside = c(0.0775, 0.865),
        legend.title = element_text(size = 15), 
        legend.margin = margin(t = -2)
    )
p3
pdf("figures/FEDRR_11062024/risk.NY-OH.pdf", width = 8, height = 7)
p3
dev.off()


#separate the states
p5 = ggplot(nmds.metadata %>% filter(State == "NY"), 
       aes(x = bray_bin.NMDS1, y = bray_bin.NMDS2, color = Risk, shape = Site)
) +
    geom_point(size = 3) +
#    scale_shape_manual(values = c(21,24)) +
    scale_color_brewer(palette = "Dark2") +
    labs(x="NMDS1", y = "NMDS2", color = "Risk level", title = "New York") +
    my_gg_theme +
    theme(
        legend.title = element_text(size = 16)
    )
p5
p6 = ggplot(nmds.metadata %>% filter(State == "OH"), 
            aes(x = bray_bin.NMDS1, y = bray_bin.NMDS2, color = Risk, shape = Site)
) +
    geom_point(size = 3) +
    #    scale_shape_manual(values = c(21,24)) +
    scale_color_brewer(palette = "Dark2") +
    labs(x="NMDS1", y = "NMDS2", color = "Risk level", title = "Ohio") +
    my_gg_theme +
    theme(
        legend.title = element_text(size = 16)
    )
p6

pdf("figures/FEDRR_11062024/risk.NY-OH_facet.pdf", width = 19, height = 6)
grid.arrange(p5,p6,ncol = 2, widths = c(0.48,0.52))
dev.off()



####################################################
####################################################
#Lure type
#
print("Alpha-pinene")
nmds.metadata$Lure[nmds.metadata$Lure == "Alpha-pinene_EtOH"] = "Alpha-pinene"

p9 = ggplot(nmds.metadata, 
            aes(x = bray_bin.NMDS1, y = bray_bin.NMDS2, fill = Lure, shape = State)
) +
    geom_point(size = 3) +
    annotate(geom = "text", label = "stress = 0.163", x = min(nmds.metadata$bray_bin.NMDS1), y = min(nmds.metadata$bray_bin.NMDS2), hjust = 0, size = 5.5) +
    scale_shape_manual(values = c(21,24)) +
    scale_fill_brewer(palette = "Set2") +
    guides(
        shape = guide_legend(order = 2, title = NULL), 
        fill = guide_legend(order = 1, override.aes = list(shape = 21))
    ) +
    labs(x="NMDS1", y = "NMDS2", fill = "Lure type") +
    my_gg_theme +
    theme(
        legend.position = "inside",
        legend.position.inside = c(0.13, 0.865),
        legend.title = element_text(size = 15), 
        legend.margin = margin(t = -2)
    )
p9
pdf("figures/FEDRR_11062024/lure.NY-OH.pdf", width = 8, height = 7)
p9
dev.off()


#separate the states
p10 = ggplot(nmds.metadata %>% filter(State == "NY"), 
            aes(x = bray_bin.NMDS1, y = bray_bin.NMDS2, color = Lure, shape = Site)
) +
    geom_point(size = 3) +
    #    scale_shape_manual(values = c(21,24)) +
    scale_color_brewer(palette = "Set2") +
    labs(x="NMDS1", y = "NMDS2", color = "Lure type", title = "New York") +
    my_gg_theme +
    theme(
        legend.title = element_text(size = 16)
    )
p10
p11 = ggplot(nmds.metadata %>% filter(State == "OH"), 
            aes(x = bray_bin.NMDS1, y = bray_bin.NMDS2, color = Lure, shape = Site)
) +
    geom_point(size = 3) +
    #    scale_shape_manual(values = c(21,24)) +
    scale_color_brewer(palette = "Set2") +
    labs(x="NMDS1", y = "NMDS2", color = "Lure type", title = "Ohio") +
    my_gg_theme +
    theme(
        legend.title = element_text(size = 16)
    )
p11

pdf("figures/FEDRR_11062024/lure.NY-OH_facet.pdf", width = 19, height = 6)
grid.arrange(p10,p11,ncol = 2, widths = c(0.48,0.52))
dev.off()

