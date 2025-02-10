library(vegan)
library(dplyr)
source("R_scripts/library.R")

#read the data
source("R_scripts/vol_test/read_data_and_split_objects.R")

#look at read count distribution

#volume test
metadata.vol.cts = metadata.vol[c("collectionID", "volumeOrWash", "seq_cts")]
metadata.vol.cts[order(metadata.vol.cts$seq_cts),]
metadata.vol.cts[order(metadata.vol.cts$seq_cts),] %>% 
    filter(volumeOrWash == 0)
# subsampling at 1900 (or 1913) will remove four samples total
# three of the samples are blanks (vol == 0) and one is TR 5-11 80 ml
# five blanks remaining
# max depth is 437445
max(metadata.vol.cts$seq_cts)/1913
#228.7
median(metadata.vol.cts$seq_cts)
#75419.5

#also look at richness of low count samples 
left_join(metadata.vol,
    data.frame(
        sequenceID = rownames(asv_tab.vol),
        richness = richness_calc(asv_tab.vol)
    )
) %>% filter(seq_cts < 1900)
#the blanks have 2-6 asvs, the 80 ml sample has 60 (also has similar number total seqs as the blank with 9 asvs)

#wash test
metadata.wash.cts = metadata.wash[c("collectionID", "volumeOrWash", "seq_cts")]
metadata.wash.cts[order(metadata.wash.cts$seq_cts),]

# The post wash samples all have lowest seq count
# min is 501 and max is 1795
# sample at 501 to retain all 
# the max seq count is 247K so we will do lot's of iterations to get full sample
max(metadata.wash.cts$seq_cts)/min(metadata.wash.cts$seq_cts)

#also look at total richness
wash.rich = left_join(metadata.wash,
          data.frame(
              sequenceID = rownames(asv_tab.wash),
              richness = richness_calc(asv_tab.wash)
          )
) 
mean(wash.rich$seq_cts)
mean(wash.rich$richness)

#perform rarefaction on vol and wash tables
rarefactions_list.vol = multiple_subsamples(x = asv_tab.vol, depth = 1913, iterations = 10000) 
rarefactions_list.vol.no_min = multiple_subsamples_no_min(x = asv_tab.vol, depth = 1913, iterations = 10000) 
rarefactions_list.wash = multiple_subsamples(x = asv_tab.wash, depth = 501, iterations = 10000) 

saveRDS(rarefactions_list.vol, "data/vol_test/rarefactions.vol.rds")
saveRDS(rarefactions_list.vol.no_min, "data/vol_test/rarefactions.vol.no_min.rds")
saveRDS(rarefactions_list.wash, "data/vol_test/rarefactions.wash.rds")

