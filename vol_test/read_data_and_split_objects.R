library(dplyr)

#read the data
asv_tab = read.table("data/vol_test/dada2_out/ASVs_counts.tsv", header = T)

metadata = read.csv("data/vol_test/metadata.csv")

#format
asv_tab.t = t(asv_tab) %>% data.frame()

#seq counts
sample_counts = data.frame(
    seq_cts = rowSums(asv_tab.t),
    sequenceID = rownames(asv_tab.t) 
    
)
metadata.cts = left_join(metadata, sample_counts, by = "sequenceID")

# filter metadata by volume or wash stage
metadata.vol = metadata.cts %>%
    filter(experiment == "volume_test")
metadata.vol$volumeOrWash = as.numeric(metadata.vol$volumeOrWash)
metadata.wash = metadata.cts %>%
    filter(experiment == "trap_wash")

#split asv_tab
asv_tab.vol = asv_tab.t[rownames(asv_tab.t) %in% metadata.vol$sequenceID,]
#nrow(asv_tab.vol)
#ncol(asv_tab.vol)
asv_tab.vol = asv_tab.vol[colSums(asv_tab.vol) > 0]
#ncol(asv_tab.vol)

asv_tab.wash = asv_tab.t[rownames(asv_tab.t) %in% metadata.wash$sequenceID,]
#nrow(asv_tab.wash)
#ncol(asv_tab.wash)
asv_tab.wash = asv_tab.wash[colSums(asv_tab.wash) > 0]
#ncol(asv_tab.wash)
#colnames(asv_tab.wash)
