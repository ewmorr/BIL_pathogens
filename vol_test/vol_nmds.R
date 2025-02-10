library(vegan)
library(dplyr)
library(ggplot2)
library(lubridate)
source("R_scripts/library.R")

#read distance
bray_dist = readRDS("data/vol_test/vol.no_min.bray.rds")
brayLog_dist = readRDS("data/vol_test/vol.no_min.bray-logCounts.rds")
brayBin_dist = readRDS("data/vol_test/vol.no_min.bray-binary.rds")

#read metadata
source("R_scripts/vol_test/read_data_and_split_objects.R")
#add a var for blanks
metadata.vol$blank = ifelse(metadata.vol$volumeOrWash == 0, "blank", "sample")

#####################################
#Run NMDS and join points to metadata
set.seed(1234)

bray.nmds = metaMDS(comm = bray_dist, autotransform = F)
brayLog.nmds = metaMDS(comm = brayLog_dist, autotransform = F)
brayBin.nmds = metaMDS(comm = brayBin_dist, autotransform = F)

str(bray.nmds)

vol.ord = data.frame(
    sequenceID = row.names(bray.nmds$points),
    NMDS1.bray = bray.nmds$points[,1],
    NMDS2.bray = bray.nmds$points[,2],
    NMDS1.brayLog = brayLog.nmds$points[,1],
    NMDS2.brayLog = brayLog.nmds$points[,2],
    NMDS1.brayBin = brayBin.nmds$points[,1],
    NMDS2.brayBin = brayBin.nmds$points[,2]
)

vol.nmds.metadata = left_join(
    vol.ord,
    metadata.vol,
    by = "sequenceID"
)

#removing all of the blanks plus TR 5-11 80 ml
# The blanks all had low sequence counts except two that had PCR contamination 
# (evident when reexamining gels; each had one rep + but not the other)
# TR 5-11 80 ml also had low seq counts (see notes in run_and_save_rarefaction.R)

#filter blanks out

sample_names = vol.nmds.metadata %>%
    filter(blank != "blank" & sampleID != "TR 5-11 80 ml") %>%
    pull("sequenceID")


bray_dist.no_blank = as.matrix(bray_dist)[sample_names,sample_names] %>% as.dist
brayLog_dist.no_blank = as.matrix(brayLog_dist)[sample_names,sample_names] %>% as.dist
brayBin_dist.no_blank = as.matrix(brayBin_dist)[sample_names,sample_names] %>% as.dist

metadata.vol.no_blank = metadata.vol %>%
    filter(blank != "blank" & sampleID != "TR 5-11 80 ml")

bray.nmds = metaMDS(comm = bray_dist.no_blank, autotransform = F)
brayLog.nmds = metaMDS(comm = brayLog_dist.no_blank, autotransform = F)
brayBin.nmds = metaMDS(comm = brayBin_dist.no_blank, autotransform = F)

vol.ord = data.frame(
    sequenceID = row.names(bray.nmds$points),
    NMDS1.bray = bray.nmds$points[,1],
    NMDS2.bray = bray.nmds$points[,2],
    NMDS1.brayLog = brayLog.nmds$points[,1],
    NMDS2.brayLog = brayLog.nmds$points[,2],
    NMDS1.brayBin = brayBin.nmds$points[,1],
    NMDS2.brayBin = brayBin.nmds$points[,2]
)

vol.nmds.metadata.no_blank = left_join(
    vol.ord,
    metadata.vol.no_blank,
    by = "sequenceID"
)


#####################################

#############
#adonis tests

#blanks
adonis2(bray_dist ~ blank, data = vol.nmds.metadata)
#p = 0.047
adonis2(brayLog_dist ~ blank, data = vol.nmds.metadata)
#p = 0.009
adonis2(brayBin_dist ~ blank, data = vol.nmds.metadata)
#p = 0.002

#site*date

#these are pseudo repped because of the volume...
adonis2(bray_dist.no_blank ~ site+mdy(date)+volumeOrWash, data = vol.nmds.metadata.no_blank, by = "margin")
#site p = 0.001
#date p = 0.001
adonis2(bray_dist.no_blank ~ site*mdy(date) + volumeOrWash^2, data = vol.nmds.metadata.no_blank, by = "margin")
#site:date p = 0.001
adonis2(bray_dist.no_blank ~ site*mdy(date)*volumeOrWash, data = vol.nmds.metadata.no_blank, by = "terms")
#the tests of marginal effects have little effect on the p-val compared to the sequential effects

adonis2(brayLog_dist.no_blank ~ site*mdy(date)*volumeOrWash, data = vol.nmds.metadata.no_blank, by = "terms")
#site p = 0.001
#date p = 0.001
#site:date p = 0.001
adonis2(brayBin_dist.no_blank ~ (site*mdy(date)*volumeOrWash)^2, data = vol.nmds.metadata.no_blank, by = "terms")
#site p = 0.001
#date p = 0.001
#site:date p = 0.001

#############


#############
#Plots

#blanks comps
p1 = ggplot(
    vol.nmds.metadata,
    aes(x = NMDS1.bray, y = NMDS2.bray, shape = blank, size = blank, color = collectionID)
) +
    geom_point() +
    scale_color_brewer(palette = "Dark2") +
    scale_shape_manual(values = c(7, 16)) +
    scale_size_manual(values = c("blank" = 3.5, "sample" = 2.5)) +
    my_gg_theme +
    labs(
        x = "NMDS1",
        y = "NMDS2",
        title = "Bray-Curtis"
    )
p1

p2 = ggplot(
    vol.nmds.metadata,
    aes(x = NMDS1.brayLog, y = NMDS2.brayLog, shape = blank, color = collectionID, size = blank)
) +
    geom_point() +
    scale_color_brewer(palette = "Dark2") +
    scale_shape_manual(values = c(7, 16)) +
    scale_size_manual(values = c("blank" = 3.5, "sample" = 2.5)) +
    my_gg_theme +
    labs(
        x = "NMDS1",
        y = "NMDS2",
        title = "Bray-Curtis (log+1)"
    )
p2

p3 = ggplot(
    vol.nmds.metadata,
    aes(x = NMDS1.brayBin, y = NMDS2.brayBin, shape = blank, color = collectionID, size = blank)
) +
    geom_point() +
    scale_color_brewer(palette = "Dark2") +
    scale_shape_manual(values = c(7, 16)) +
    scale_size_manual(values = c("blank" = 3.5, "sample" = 2.5)) +
    my_gg_theme +
    labs(
        x = "NMDS1",
        y = "NMDS2",
        title = "Binary Bray-Curtis"
    )
p3

p4 = ggplot(
    vol.nmds.metadata,
    aes(x = NMDS1.bray, y = NMDS2.bray, shape = blank, size = as.factor(volumeOrWash), color = collectionID)
) +
    geom_point() +
    scale_color_brewer(palette = "Dark2") +
    scale_shape_manual(values = c(7, 16)) +
    scale_size_manual(values = c("0" = 3.5, "5" = 1, "10" = 1.5, "20" = 2, "40" = 2.5, "80" = 3)) +
    my_gg_theme +
    labs(
        x = "NMDS1",
        y = "NMDS2",
        title = "Bray-Curtis"
    )
p4

p5 = ggplot(
    vol.nmds.metadata,
    aes(x = NMDS1.brayLog, y = NMDS2.brayLog, shape = blank, color = collectionID, size = as.factor(volumeOrWash))
) +
    geom_point() +
    scale_color_brewer(palette = "Dark2") +
    scale_shape_manual(values = c(7, 16)) +
    scale_size_manual(values = c("0" = 3.5, "5" = 1, "10" = 1.5, "20" = 2, "40" = 2.5, "80" = 3)) +
    my_gg_theme +
    labs(
        x = "NMDS1",
        y = "NMDS2",
        title = "Bray-Curtis (log+1)"
    )
p5

p6 = ggplot(
    vol.nmds.metadata,
    aes(x = NMDS1.brayBin, y = NMDS2.brayBin, shape = blank, color = collectionID, size = as.factor(volumeOrWash))
) +
    geom_point() +
    scale_color_brewer(palette = "Dark2") +
    scale_shape_manual(values = c(7, 16)) +
    scale_size_manual(values = c("0" = 3.5, "5" = 1, "10" = 1.5, "20" = 2, "40" = 2.5, "80" = 3)) +
    my_gg_theme +
    labs(
        x = "NMDS1",
        y = "NMDS2",
        title = "Binary Bray-Curtis"
    )
p6

pdf("figures/vol_test/vol_nmds.blanks.pdf", width = 8, height = 6)
p1
p2
p3
p4
p5
p6
dev.off()


####
# No blanks
# TR 5-11 80 ml also removed do to low seq depth (see above)


####
#shape size is vol

p1 = ggplot(
    vol.nmds.metadata.no_blank,
    aes(x = NMDS1.bray, y = NMDS2.bray, size = as.factor(volumeOrWash), color = collectionID)
) +
    geom_point() +
    scale_color_brewer(palette = "Dark2") +
    scale_size_manual(values = c("5" = 1.5, "10" = 2, "20" = 2.5, "40" = 3, "80" = 3.5)) +
    #scale_size_manual(values = c("5" = 1.5, "10" = 2.5, "20" = 3.5, "40" = 4.5, "80" = 5.5)) +
    my_gg_theme +
    labs(
        x = "NMDS1",
        y = "NMDS2",
        title = "Bray-Curtis",
        size = "Volume (ml)",
        color = "Collection ID"
    ) +
    theme(
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16)
    )
p1

p2 = ggplot(
    vol.nmds.metadata.no_blank,
    aes(x = NMDS1.brayLog, y = NMDS2.brayLog, color = collectionID, size = as.factor(volumeOrWash))
) +
    geom_point() +
    scale_color_brewer(palette = "Dark2") +
    scale_size_manual(values = c("5" = 1.5, "10" = 2, "20" = 2.5, "40" = 3, "80" = 3.5)) +
    #scale_size_manual(values = c("5" = 1.5, "10" = 2.5, "20" = 3.5, "40" = 4.5, "80" = 5.5)) +
    my_gg_theme +
    labs(
        x = "NMDS1",
        y = "NMDS2",
        title = "Bray-Curtis (log+1)",
        size = "Volume (ml)",
        color = "Collection ID"
    ) +
    theme(
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16)
    )
p2

p3 = ggplot(
    vol.nmds.metadata.no_blank,
    aes(x = NMDS1.brayBin, y = NMDS2.brayBin, color = collectionID, size = as.factor(volumeOrWash))
) +
    geom_point() +
    scale_color_brewer(palette = "Dark2") +
    scale_shape_manual(values = c(7, 16)) +
    scale_size_manual(values = c("5" = 1.5, "10" = 2, "20" = 2.5, "40" = 3, "80" = 3.5)) +
    #scale_size_manual(values = c("5" = 1.5, "10" = 2.5, "20" = 3.5, "40" = 4.5, "80" = 5.5)) +
    my_gg_theme +
    labs(
        x = "NMDS1",
        y = "NMDS2",
        title = "Binary Bray-Curtis",
        size = "Volume (ml)",
        color = "Collection ID"
    ) +
    theme(
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16)
    )
p3


#site : date comps
vol.nmds.metadata.no_blank$mdy = mdy(vol.nmds.metadata.no_blank$date) 


#labeller function for color scale
#this can be used if feeding date into scale_x_gradient directly
Date_F <- function(x){
    month.name[as.numeric(format(as.Date(x, format = "%m-%d-%Y"),"%m"))]
}
Date_F(vol.nmds.metadata.no_blank$mdy)

#function to convert integer dates back to breaks and date labels
#use with scale_x_gradient2 if it is desired to define a midpoint 
date_breaks <- function(x){
    #breaks <- c(min(x),median(x),max(x))
    breaks <- quantile(x, probs = seq(0, 1, 0.25))
    attr(breaks,"labels") <- month.name[as.numeric(format(as.Date(breaks, origin="1970-01-01"),"%m"))]
    names(breaks) <- attr(breaks,"labels")
    return(breaks)
}

date_breaks(as.integer(vol.nmds.metadata.no_blank$mdy))

#April    May August 
#19460  19502  19586 
#mean of min-max is 19523

#April   April     May    June  August 
#19460.0 19477.5 19502.0 19516.0 19586.0 

five_cols_gradient_palette = c('#ca0020','#f4a582','#f7f7f7','#92c5de','#0571b0')
five_cols_gradient_palette = c('#b2182b','#d6604d','#f4a582','white','#2166ac')

p4 = ggplot(
    vol.nmds.metadata.no_blank,
    aes(x = NMDS1.bray, y = NMDS2.bray, shape = site, fill = as.integer(mdy))
) +
    geom_point(size = 3.5) +
    scale_shape_manual(values = c(21,22,24)) +
    #scale_fill_gradient2(high = "#b2182b", low = "#2166ac", mid = "white", midpoint = as.integer(as.Date("2023-06-01")), breaks = date_breaks) +
    scale_fill_gradientn(colours = rev(five_cols_gradient_palette), breaks = date_breaks) +
    #scale_fill_gradient(high = "#ca0020", low = "#0571b0", labels = Date_F) +
    my_gg_theme +
    labs(
        fill = "Date",
        shape = "Site",
        x = "NMDS1",
        y = "NMDS2",
        title = "Bray-Curtis"
    ) +
    theme(
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16)
    )
p4

p5 = ggplot(
    vol.nmds.metadata.no_blank,
    aes(x = NMDS1.brayLog, y = NMDS2.brayLog, shape = site, fill = mdy)
) +
    geom_point(size = 3.5) +
    scale_shape_manual(values = c(21,22,24)) +
    scale_fill_gradientn(colours = rev(five_cols_gradient_palette), breaks = date_breaks) +
    my_gg_theme +
    labs(
        fill = "Date",
        shape = "Site",
        x = "NMDS1",
        y = "NMDS2",
        title = "Bray-Curtis (log+1)"
    ) +
    theme(
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16)
    )
p5

#site and date
p6 = ggplot(
    vol.nmds.metadata.no_blank,
    aes(x = NMDS1.brayBin, y = NMDS2.brayBin, shape = site, fill = mdy)
) +
    geom_point(size = 3.5) +
    scale_shape_manual(values = c(21,22,24)) +
    scale_fill_gradientn(colours = rev(five_cols_gradient_palette), breaks = date_breaks) +
    my_gg_theme +
    labs(
        fill = "Date",
        shape = "Site",
        x = "NMDS1",
        y = "NMDS2",
        title = "Binary Bray-Curtis"
    ) +
    theme(
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16)
    )
p6

#date and site with size indcating vol

p7 = ggplot(
    vol.nmds.metadata.no_blank,
    aes(
        x = NMDS1.bray, 
        y = NMDS2.bray, 
        shape = site, 
        fill = as.integer(mdy), 
        size = as.factor(volumeOrWash)
    )
) +
    geom_point() +
    scale_shape_manual(values = c(21,22,24)) +
    #scale_size_manual(values = c("5" = 1.5, "10" = 2, "20" = 2.5, "40" = 3.5, "80" = 5)) +
    #scale_size_manual(values = c("5" = 1.5, "10" = 2.5, "20" = 3.5, "40" = 4.5, "80" = 5.5)) +
    scale_size_manual(values = c("5" = 1.5, "10" = 2, "20" = 2.5, "40" = 3, "80" = 3.5)) +
    scale_fill_gradientn(colours = rev(five_cols_gradient_palette), breaks = date_breaks) +
    my_gg_theme +
    labs(
        fill = "Date",
        shape = "Site",
        size = "Volume (ml)",
        x = "NMDS1",
        y = "NMDS2",
        title = "Bray-Curtis"
    ) +
    theme(
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16)
    ) +
    guides(
        shape = guide_legend(override.aes = list(size = 3))
    ) 
p7

p8 = ggplot(
    vol.nmds.metadata.no_blank,
    aes(x = NMDS1.brayLog, y = NMDS2.brayLog, shape = site, fill = mdy, size = as.factor(volumeOrWash))
) +
    geom_point() +
    scale_shape_manual(values = c(21,22,24)) +
    scale_size_manual(values = c("5" = 1.5, "10" = 2, "20" = 2.5, "40" = 3, "80" = 3.5)) +
    scale_fill_gradientn(colours = rev(five_cols_gradient_palette), breaks = date_breaks) +
    my_gg_theme +
    labs(
        fill = "Date",
        shape = "Site",
        size = "Volume (ml)",
        x = "NMDS1",
        y = "NMDS2",
        title = "Bray-Curtis (log+1)"
    ) +
    theme(
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16)
    ) +
    guides(
        shape = guide_legend(override.aes = list(size = 3))
    ) 
p8

#site and date
p9 = ggplot(
    vol.nmds.metadata.no_blank,
    aes(x = NMDS1.brayBin, y = NMDS2.brayBin, shape = site, fill = mdy, size = as.factor(volumeOrWash))
) +
    geom_point() +
    scale_shape_manual(values = c(21,22,24)) +
    scale_size_manual(values = c("5" = 1.5, "10" = 2, "20" = 2.5, "40" = 3, "80" = 3.5)) +
    scale_fill_gradientn(colours = rev(five_cols_gradient_palette), breaks = date_breaks) +
    my_gg_theme +
    labs(
        fill = "Date",
        shape = "Site",
        size = "Volume (ml)",
        x = "NMDS1",
        y = "NMDS2",
        title = "Binary Bray-Curtis"
    ) +
    theme(
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16)
    ) +
    guides(
        shape = guide_legend(override.aes = list(size = 3))
    ) 
p9

pdf("figures/vol_test/vol_nmds.no_blanks.pdf", width = 8, height = 6)
p1
p2
p3
p4
p5
p6
p7
p8
p9
dev.off()
