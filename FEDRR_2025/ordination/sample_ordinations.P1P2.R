library(dplyr)
library(vegan)
library(ggplot2)
library(lubridate)
library(gridExtra)
library(gtable)
source("library/library.R")
set.seed(12345)

metadata = read.csv("data/2025_metadata/metadata_P1P2.csv")
head(metadata)

#full rarefied table
#
bray_bin.dist = readRDS("data/FEDRR_2025_11182025/rarefaction_avgs/bray-binary.samps.rds")
bray_log.dist = readRDS("data/FEDRR_2025_11182025/rarefaction_avgs/bray-logCounts.samps.rds")


##################################
#perform NMDS
bray_bin.nmds = metaMDS(bray_bin.dist, try = 20, trymax = 100)
bray_bin.nmds = metaMDS(bray_bin.dist, trymax = 200, previous.best = bray_bin.nmds)
bray_bin.nmds
#stress = 0.2029699

bray_log.nmds = metaMDS(bray_log.dist, try = 20, trymax = 100)
bray_log.nmds = metaMDS(bray_log.dist, previous.best = bray_log.nmds)
bray_log.nmds
#stress = 0.2083004

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
nmds.metadata$sampleID %>% unique
#nmds.metadata %>% filter(State == "NY") %>% pull(trapID) %>% unique
nmds.metadata$dateOrBaselineType %>% unique()
nmds.metadata %>% filter(dateOrBaselineType == "BLO")

nmds.metadata$date %>% class()
mdy(nmds.metadata$date)

nmds.metadata$mdy = mdy(nmds.metadata$date)
nmds.metadata$mdy %>% unique()
nmds.metadata$mdy %>% class()
nmds.metadata$mdy %>% as.character()


####################################################
# quick stats. We are not using how() 
# bc the sampling design is not compete yet
adonis2(bray_bin.dist ~ State, data = nmds.metadata)
#          Df SumOfSqs      R2      F Pr(>F)    
#Model      5    9.443 0.15406 5.5362  0.001 ***
#Residual 152   51.854 0.84594                  
#Total    157   61.297 1.00000                  

adonis2(bray_bin.dist ~ State + mdy, data = nmds.metadata, by = "terms")
#          Df SumOfSqs      R2      F Pr(>F)    
#State      5    9.443 0.15406 5.6914  0.001 ***
#mdy        1    1.746 0.02848 5.2607  0.001 ***
#Residual 151   50.108 0.81746                  
#Total    157   61.297 1.00000                  

adonis2(bray_bin.dist ~ State + Lure, data = nmds.metadata, by = "terms")
#          Df SumOfSqs      R2      F Pr(>F)    
#State      5    9.443 0.15406 5.5555  0.001 ***
#Lure       3    1.200 0.01957 1.1762  0.049 *  
#Residual 149   50.654 0.82637                  
#Total    157   61.297 1.00000                

adonis2(bray_bin.dist ~ State + mdy + Lure, data = nmds.metadata, by = "terms")
#          Df SumOfSqs      R2      F Pr(>F)    
#State      5    9.443 0.15406 5.7064  0.001 ***
#mdy        1    1.746 0.02848 5.2746  0.001 ***
#Lure       3    1.125 0.01835 1.1331  0.115    
#Residual 148   48.983 0.79911                  
#Total    157   61.297 1.00000               
adonis2(bray_bin.dist ~ State + mdy + Lure, data = nmds.metadata, by = "margin")
# this only changes the R^2 marginally from above

####################################################

p0 = ggplot(nmds.metadata, 
            aes(x = bray_bin.NMDS1, y = bray_bin.NMDS2, color = State, shape = State)
) +
    geom_point(size = 3) +
    #scale_shape_manual(values = c(21,24,23,22)) +
    annotate(
        geom = "text", 
        label = expr(paste("stress = ",!!round(bray_bin.nmds$stress,2))), 
        x = min(nmds.metadata$bray_bin.NMDS1), y = min(nmds.metadata$bray_bin.NMDS2), 
        hjust = 0, size = 5.5
    ) +
    scale_color_manual(values = rev(cbPalette)) +
    #scale_color_brewer(palette = "Set2") +
    #guides(
    #    shape = guide_legend(order = 2, title = NULL), 
    #    fill = guide_legend(order = 1, override.aes = list(shape = 21))
    #) +
    labs(x="NMDS1", y = "NMDS2") +
    my_gg_theme +
    theme(
        legend.position = "inside",
        legend.position.inside = c(0.95, 0.875),
        legend.title = element_text(size = 15), 
        legend.margin = margin(t = -2)
    )
p0

pdf("figures/FEDRR_2025_11182025/state.pdf", width = 9, height = 8)
p0
dev.off()

#note if wish to set breaks manually can use date_breaks function
#e.g., scale_fill_gradientn(colours = rev(five_cols_gradient_palette), breaks = date_breaks)
#see lib and vol_nmds.r script
#
date_breaks(as.integer(nmds.metadata %>% filter(!is.na(mdy)) %>% pull(mdy)))
five_cols_gradient_palette = c('#ca0020','#f4a582','#f7f7f7','#92c5de','#0571b0')
five_cols_gradient_palette = c('#b2182b','#d6604d','#f4a582','white','#2166ac')

five_cols_gradient_palette = c('#ca0020','#ef8a62','#f7f7f7','#67a9cf','#2166ac')

#fill by date shape is state

class(nmds.metadata$mdy)
year(nmds.metadata$mdy)

p1 = ggplot(nmds.metadata, 
       aes(x = bray_bin.NMDS1, y = bray_bin.NMDS2, shape = State, color = mdy)
) +
    geom_point(size = 3) +
    #scale_shape_manual(values = c(21,24,23,22)) +
    scale_color_gradientn(colours = rev(five_cols_gradient_palette), breaks = date_breaks) +
    annotate(geom = "text", label = expr(paste("stress = ",!!round(bray_bin.nmds$stress,2))), x = min(nmds.metadata$bray_bin.NMDS1), y = min(nmds.metadata$bray_bin.NMDS2), hjust = 0, size = 5.5) +
    guides(shape = guide_legend(order = 2), color = guide_colorbar(order = 1)) +
    labs(x="NMDS1", y = "NMDS2") +
    my_gg_theme +
    #annotate(geom = "segment", x = -0.275, xend = -0.275, y = 0.15, yend = Inf) +
    #annotate(geom = "segment", x = -Inf, xend = -0.275, y = 0.15, yend = 0.15) +
    theme(
        legend.position = "inside",
        #legend.position.inside = c(0.0665, 0.8182)
        legend.position.inside = c(0.942, 0.775)
    ) 
p1

pdf("figures/FEDRR_2025_11182025/date.pdf", width = 9, height = 8)
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
p9 = ggplot(nmds.metadata %>% filter(State == "IA"), 
            aes(x = bray_bin.NMDS1, y = bray_bin.NMDS2, color = as.integer(mdy), shape = Site)
) +
    scale_shape_manual(values = c(16,17,15,3,7,8,5)) +
    geom_point(size = 3) +
    scale_color_gradientn(colours = rev(five_cols_gradient_palette), breaks = date_breaks, guide = guide_colorbar(title = element_blank(), order = 2)) +
    labs(x="NMDS1", y = "NMDS2", color = "Risk level", title = "Iowa") +
    my_gg_theme +
    theme(
        legend.title = element_text(size = 16)
    )
p9
p10 = ggplot(nmds.metadata %>% filter(State == "NH"), 
            aes(x = bray_bin.NMDS1, y = bray_bin.NMDS2, color = as.integer(mdy), shape = Site)
) +
    scale_shape_manual(values = c(16,17,15,3,7,8,5)) +
    geom_point(size = 3) +
    scale_color_gradientn(colours = rev(five_cols_gradient_palette), breaks = date_breaks, guide = guide_colorbar(title = element_blank(), order = 2)) +
    labs(x="NMDS1", y = "NMDS2", color = "Risk level", title = "New Hampshire") +
    my_gg_theme +
    theme(
        legend.title = element_text(size = 16)
    )
p10

gp9 = ggplotGrob(p9)
gp10 = ggplotGrob(p10)
gp7 = ggplotGrob(p7)
gp8 = ggplotGrob(p8)

gpLeft = rbind(gp9, gp7)
gpRight = rbind(gp10, gp8)

pdf("figures/FEDRR_2025_11182025/date.NY-OH-IA-NH_facet.pdf", width = 17, height = 10)
grid.arrange(gpLeft, gpRight,ncol = 2, widths = c(0.49,0.51))
dev.off()

####################################################
####################################################
#Site risk
#

p3 = ggplot(nmds.metadata, 
            aes(x = bray_bin.NMDS1, y = bray_bin.NMDS2, fill = factor(Risk, levels = c("Low", "High")), shape = State)
) +
    geom_point(size = 3) +
    scale_shape_manual(values = c(21,24,23,22)) +
    annotate(geom = "text", label = expr(paste("stress = ",!!round(bray_bin.nmds$stress,2))), x = min(nmds.metadata$bray_bin.NMDS1), y = min(nmds.metadata$bray_bin.NMDS2), hjust = 0, size = 5.5) +
    scale_fill_brewer(palette = "Dark2") +
    guides(
        shape = guide_legend(order = 2, title = NULL), 
        fill = guide_legend(order = 1, override.aes = list(shape = 21))
    ) +
    labs(x="NMDS1", y = "NMDS2", fill = "Risk level") +
    my_gg_theme +
    guides(color = guide_legend(order = 1)) +
    theme(
        legend.position = "inside",
        legend.position.inside = c(0.9325, 0.15),
        legend.title = element_text(size = 15), 
        legend.margin = margin(t = -2)
    )
p3
pdf("figures/FEDRR_2025_11182025/risk.NY-OH-IA-NH.pdf", width = 9, height = 8)
p3
dev.off()


#separate the states
p5 = ggplot(nmds.metadata %>% filter(State == "IA"), 
       aes(x = bray_bin.NMDS1, y = bray_bin.NMDS2, color = factor(Risk, levels = c("Low", "High")), shape = Site)
) +
    geom_point(size = 3) +
    scale_shape_manual(values = c(16,17,15,3,7,8,5)) +
    scale_color_brewer(palette = "Dark2") +
    labs(x="NMDS1", y = "NMDS2", color = "Risk level", title = "Iowa") +
    my_gg_theme +
    guides(color = guide_legend(order = 1)) +
    theme(
        legend.title = element_text(size = 16)
    )
p5
p6 = ggplot(nmds.metadata %>% filter(State == "NH"), 
            aes(x = bray_bin.NMDS1, y = bray_bin.NMDS2, color = factor(Risk, levels = c("Low", "High")), shape = Site)
) +
    geom_point(size = 3) +
    #    scale_shape_manual(values = c(21,24)) +
    scale_color_brewer(palette = "Dark2") +
    labs(x="NMDS1", y = "NMDS2", color = "Risk level", title = "New Hampshire") +
    my_gg_theme +
    guides(color = guide_legend(order = 1)) +
    theme(
        legend.title = element_text(size = 16)
    )
p6
p7 = ggplot(nmds.metadata %>% filter(State == "NY"), 
            aes(x = bray_bin.NMDS1, y = bray_bin.NMDS2, color = factor(Risk, levels = c("Low", "High")), shape = Site)
) +
    geom_point(size = 3) +
    #    scale_shape_manual(values = c(21,24)) +
    scale_color_brewer(palette = "Dark2") +
    labs(x="NMDS1", y = "NMDS2", color = "Risk level", title = "New York") +
    my_gg_theme +
    guides(color = guide_legend(order = 1)) +
    theme(
        legend.title = element_text(size = 16)
    )
p7
p8 = ggplot(nmds.metadata %>% filter(State == "OH"), 
            aes(x = bray_bin.NMDS1, y = bray_bin.NMDS2, color = factor(Risk, levels = c("Low", "High")), shape = Site)
) +
    geom_point(size = 3) +
    #    scale_shape_manual(values = c(21,24)) +
    scale_color_brewer(palette = "Dark2") +
    labs(x="NMDS1", y = "NMDS2", color = "Risk level", title = "Ohio") +
    my_gg_theme +
    guides(color = guide_legend(order = 1)) +
    theme(
        legend.title = element_text(size = 16)
    )
p8

gp5 = ggplotGrob(p5)
gp6 = ggplotGrob(p6)
gp7 = ggplotGrob(p7)
gp8 = ggplotGrob(p8)

gpLeft = rbind(gp5, gp7)
gpRight = rbind(gp6, gp8)

pdf("figures/FEDRR_2025_11182025/risk.NY-OH-IA-NH_facet.pdf", width = 17, height = 10)
grid.arrange(gpLeft,gpRight,ncol = 2, widths = c(0.49,0.51))
dev.off()



####################################################
####################################################
#Lure type
#
print("Alpha-pinene")
nmds.metadata$Lure %>% unique()
#nmds.metadata$Lure[nmds.metadata$Lure == "Alpha-pinene_EtOH"] = "Alpha-pinene"
#nmds.metadata$Lure[nmds.metadata$Lure == "Alpha-pinene"] = "AP"

p9 = ggplot(nmds.metadata, 
            aes(x = bray_bin.NMDS1, y = bray_bin.NMDS2, color = Lure, shape = State)
) +
    geom_point(size = 3) +
    annotate(geom = "text", label = expr(paste("stress = ",!!round(bray_bin.nmds$stress,2))), x = min(nmds.metadata$bray_bin.NMDS1), y = min(nmds.metadata$bray_bin.NMDS2), hjust = 0, size = 5.5) +
    #scale_shape_manual(values = c(21,24,23,22)) +
    scale_color_brewer(palette = "Set2") +
    guides(
        shape = guide_legend(order = 2, title = NULL), 
        fill = guide_legend(order = 1, override.aes = list(shape = 21))
    ) +
    labs(x="NMDS1", y = "NMDS2", fill = "Lure type") +
    my_gg_theme +
    theme(
        legend.position = "inside",
        legend.position.inside = c(0.95, 0.8),
        legend.title = element_text(size = 15), 
        legend.margin = margin(t = -2)
    )
p9
pdf("figures/FEDRR_2025_11182025/lure.pdf", width = 9, height = 8)
p9
dev.off()


#separate the states
p5 = ggplot(nmds.metadata %>% filter(State == "IA"), 
            aes(x = bray_bin.NMDS1, y = bray_bin.NMDS2, color = Lure, shape = Site)
) +
    geom_point(size = 3) +
    scale_shape_manual(values = c(16,17,15,3,7,8,5)) +
    scale_color_brewer(palette = "Set2") +
    labs(x="NMDS1", y = "NMDS2", color = "Lure type", title = "Iowa") +
    guides(color = guide_legend(order = 1)) +
    my_gg_theme +
    theme(
        legend.title = element_text(size = 16)
    )
p5
p6 = ggplot(nmds.metadata %>% filter(State == "NH"), 
            aes(x = bray_bin.NMDS1, y = bray_bin.NMDS2, color = Lure, shape = Site)
) +
    geom_point(size = 3) +
    #    scale_shape_manual(values = c(21,24)) +
    scale_color_brewer(palette = "Set2") +
    labs(x="NMDS1", y = "NMDS2", color = "Lure type", title = "New Hampshire") +
    guides(color = guide_legend(order = 1)) +
    my_gg_theme +
    theme(
        legend.title = element_text(size = 16)
    )
p6
p7 = ggplot(nmds.metadata %>% filter(State == "NY"), 
            aes(x = bray_bin.NMDS1, y = bray_bin.NMDS2, color = Lure, shape = Site)
) +
    geom_point(size = 3) +
    #    scale_shape_manual(values = c(21,24)) +
    scale_color_brewer(palette = "Set2") +
    labs(x="NMDS1", y = "NMDS2", color = "Lure type", title = "New York") +
    guides(color = guide_legend(order = 1)) +
    my_gg_theme +
    theme(
        legend.title = element_text(size = 16)
    )
p7
p8 = ggplot(nmds.metadata %>% filter(State == "OH"), 
            aes(x = bray_bin.NMDS1, y = bray_bin.NMDS2, color = Lure, shape = Site)
) +
    geom_point(size = 3) +
    #    scale_shape_manual(values = c(21,24)) +
    scale_color_brewer(palette = "Set2") +
    labs(x="NMDS1", y = "NMDS2", color = "Lure type", title = "Ohio") +
    guides(color = guide_legend(order = 1)) +
    my_gg_theme +
    theme(
        legend.title = element_text(size = 16)
    )
p8

gp5 = ggplotGrob(p5)
gp6 = ggplotGrob(p6)
gp7 = ggplotGrob(p7)
gp8 = ggplotGrob(p8)

gpLeft = rbind(gp5, gp7)
gpRight = rbind(gp6, gp8)


pdf("figures/FEDRR_2025_11182025/lure.NY-OH-IA-NH_facet.pdf", width = 17, height = 10)
grid.arrange(gpLeft,gpRight,ncol = 2, widths = c(0.49,0.51))
dev.off()

#gpTop = cbind(gp5, gp6)
#gpBottom = cbind(gp7, gp8)
#gpAll = rbind(gpTop, gpBottom)

#pdf("figures/FEDRR_2025_11182025/lure.NY-OH-IA-NH_facet.pdf", width = 17, height = 10)
#plot(gpAll)
#dev.off()



