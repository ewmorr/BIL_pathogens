library(vegan)
library(dplyr)
library(ggplot2)
library(lubridate)
source("R_scripts/library.R")

#read distance
bray_dist = readRDS("data/vol_test/wash.bray.rds")
brayLog_dist = readRDS("data/vol_test/wash.bray-logCounts.rds")
brayBin_dist = readRDS("data/vol_test/wash.bray-binary.rds")

#read metadata
source("R_scripts/vol_test/read_data_and_split_objects.R")

#####################################
#Run NMDS and join points to metadata
set.seed(1234)

bray.nmds = metaMDS(comm = bray_dist, autotransform = F)
brayLog.nmds = metaMDS(comm = brayLog_dist, autotransform = F)
brayBin.nmds = metaMDS(comm = brayBin_dist, autotransform = F)

str(bray.nmds)

#warnings messages with log+1 and binary metrics: stress is nearly zero
# J Oksanen suggests using metric scaling (dbrda/capscale) instead of nmds in this case
?dbrda
bray.pcoa = capscale(bray_dist ~ 1)
brayLog.pcoa = capscale(brayLog_dist ~ 1)
brayBin.pcoa = capscale(brayBin_dist ~ 1)

str(bray.pcoa)

#eigenvalues
bray.pcoa.eig = bray.pcoa$CA$eig/sum(bray.pcoa$CA$eig)
brayLog.pcoa.eig = brayLog.pcoa$CA$eig/sum(brayLog.pcoa$CA$eig)
brayBin.pcoa.eig = brayBin.pcoa$CA$eig/sum(brayBin.pcoa$CA$eig)

#add NMDS in case want to look
wash.ord = data.frame(
    sequenceID = row.names(bray.nmds$points),
    NMDS1.bray = bray.nmds$points[,1],
    NMDS2.bray = bray.nmds$points[,2],
    NMDS1.brayLog = brayLog.nmds$points[,1],
    NMDS2.brayLog = brayLog.nmds$points[,2],
    NMDS1.brayBin = brayBin.nmds$points[,1],
    NMDS2.brayBin = brayBin.nmds$points[,2]
)

wash.nmds.metadata = left_join(
    wash.ord,
    metadata.wash,
    by = "sequenceID"
)

#then add pcoa axes
#add NMDS in case want to look
wash.ord = data.frame(
    sequenceID = row.names(bray.pcoa$CA$u),
    MDS1.bray = bray.pcoa$CA$u[,1],
    MDS2.bray = bray.pcoa$CA$u[,2],
    MDS1.brayLog = brayLog.pcoa$CA$u[,1],
    MDS2.brayLog = brayLog.pcoa$CA$u[,2],
    MDS1.brayBin = brayBin.pcoa$CA$u[,1],
    MDS2.brayBin = brayBin.pcoa$CA$u[,2]
)

wash.pcoa.metadata = left_join(
    metadata.wash,
    wash.ord,
    by = "sequenceID"
)

#####################################

#############
#adonis tests

#blanks
adonis2(bray_dist ~ volumeOrWash, data = wash.nmds.metadata)
#p = 0.029
adonis2(brayLog_dist ~ volumeOrWash, data = wash.nmds.metadata)
#p = 0.02
adonis2(brayBin_dist ~ volumeOrWash, data = wash.nmds.metadata)
#p = 0.026


#############


#############
#Plots

p1 = ggplot(
    wash.nmds.metadata,
    aes(x = NMDS1.bray, y = NMDS2.bray, color = factor(volumeOrWash, levels = c("pre", "mid", "post")))
) +
    geom_point(size = 3) +
    scale_color_manual(values = cbPalette) +
    my_gg_theme +
    labs(
        x = "NMDS1",
        y = "NMDS2",
        title = "Bray-Curtis"
    )
p1

p2 = ggplot(
    wash.nmds.metadata,
    aes(x = NMDS1.brayLog, y = NMDS2.brayLog, color = factor(volumeOrWash, levels = c("pre", "mid", "post")))
) +
    geom_point(size = 3) +
    scale_color_manual(values = cbPalette) +
    my_gg_theme +
    labs(
        x = "NMDS1",
        y = "NMDS2",
        title = "Bray-Curtis (log+1)"
    )
p2

p3 = ggplot(
    wash.nmds.metadata,
    aes(x = NMDS1.brayBin, y = NMDS2.brayBin, color = factor(volumeOrWash, levels = c("pre", "mid", "post")))
) +
    geom_point(size = 3) +
    scale_color_manual(values = cbPalette) +
    my_gg_theme +
    labs(
        x = "NMDS1",
        y = "NMDS2",
        title = "Binary Bray-Curtis"
    )
p3

p4 = ggplot(
    wash.pcoa.metadata,
    aes(x = MDS1.bray, y = MDS2.bray, color = factor(volumeOrWash, levels = c("pre", "mid", "post")))
) +
    geom_point(size = 3) +
    scale_color_manual(values = cbPalette) +
    my_gg_theme +
    labs(
        x = paste(paste("PCoA1 (", round(bray.pcoa.eig[1]*100, 0), sep = " "), "% variance)", sep = ""),
        y = paste(paste("PCoA2 (", round(bray.pcoa.eig[2]*100, 0), sep = " "), "% variance)", sep = ""),
        title = "Bray-Curtis"
    )
p4

p5 = ggplot(
    wash.pcoa.metadata,
    aes(x = MDS1.brayLog, y = MDS2.brayLog, color = factor(volumeOrWash, levels = c("pre", "mid", "post")))
) +
    geom_point(size = 3) +
    scale_color_manual(values = cbPalette) +
    my_gg_theme +
    labs(
        x = paste(paste("PCoA1 (", round(brayLog.pcoa.eig[1]*100, 0), sep = " "), "% variance)", sep = ""),
        y = paste(paste("PCoA2 (", round(brayLog.pcoa.eig[2]*100, 0), sep = " "), "% variance)", sep = ""),
        title = "Bray-Curtis (log+1)"
    )
p5

p6 = ggplot(
    wash.pcoa.metadata,
    aes(x = MDS1.brayBin, y = MDS2.brayBin, color = factor(volumeOrWash, levels = c("pre", "mid", "post")))
) +
    geom_point(size = 3) +
    scale_color_manual(values = cbPalette) +
    my_gg_theme +
    labs(
        x = paste(paste("PCoA1 (", round(brayBin.pcoa.eig[1]*100, 0), sep = " "), "% variance)", sep = ""),
        y = paste(paste("PCoA2 (", round(brayBin.pcoa.eig[2]*100, 0), sep = " "), "% variance)", sep = ""),
        title = "Binary Bray-Curtis"
    )
p6

pdf("figures/vol_test/wash_nmds.pdf", width = 8, height = 6)
p1
p2
p3
p4
p5
p6
dev.off()

