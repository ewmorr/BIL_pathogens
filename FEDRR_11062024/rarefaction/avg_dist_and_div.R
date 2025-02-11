library(vegan)
library(dplyr)
source("library/library.R")

######################
# full table 500 rarefactions 450 min
rarefactions_list = readRDS("data/FEDRR_11062024/processed_tables/rarefactions.full_tab_min450.rds")

#calc distance and alpha-div
## take avg
## save
## rm objs
## 
## binary bray
bray_binary_list = lapply(rarefactions_list, vegdist, method = "bray", binary = T, diag = T, upper = T)
bray_binary_avg = avg_matrix_list(bray_binary_list) %>% as.dist()
saveRDS(bray_binary_avg, "data/FEDRR_11062024/processed_tables/avg_dist/bray-binary.full_tab_min450.rds")
rm(bray_binary_list)
rm(bray_binary_avg)
gc()
        
#log bray
bray_logCts_list = lapply(rarefactions_list, log_dist, method = "bray")
bray_logCts_avg = avg_matrix_list(bray_logCts_list) %>% as.dist()
saveRDS(bray_logCts_avg, "data/FEDRR_11062024/processed_tables/avg_dist/bray-logCounts.full_tab_min450.rds")
rm(bray_logCts_list)
rm(bray_logCts_avg)
gc()



#alpha divs
shannon_list = lapply(rarefactions_list, diversity, index = "shannon")
simpson_list = lapply(rarefactions_list, diversity, index = "simpson")
richness_list = lapply(rarefactions_list, richness_calc)

div_avg = data.frame(
    sample = rownames(avg_matrix_list(shannon_list)),
    shannon = avg_matrix_list(shannon_list),
    simpson = avg_matrix_list(simpson_list),
    richness = avg_matrix_list(richness_list),
    stringsAsFactors = F
)
write.csv(div_avg, "data/FEDRR_11062024/processed_tables/avg_div/diversity.full_tab_min450.csv", row.names = F)
rm(shannon_list)
rm(simpson_list)
rm(richness_list)
rm(div_avg)
gc()

###########################################
#Also average the counts table
avg_counts = avg_matrix_list(rarefactions_list)

rm(rarefactions_list)
gc()

write.csv(avg_counts, "data/FEDRR_11062024/processed_tables/asv_tab.rarefaction_avg.full_tab_min450.csv")

rm(avg_counts)
gc()


############################################
############################################
############################################
############################################
############################################
# samps
rarefactions_list = readRDS("data/FEDRR_11062024/processed_tables/rarefactions.samps.rds")

#calc distance and alpha-div
## take avg
## save
## rm objs
## 
## binary bray
bray_binary_list = lapply(rarefactions_list, vegdist, method = "bray", binary = T, diag = T, upper = T)
bray_binary_avg = avg_matrix_list(bray_binary_list) %>% as.dist()
saveRDS(bray_binary_avg, "data/FEDRR_11062024/processed_tables/avg_dist/bray-binary.samps.rds")
rm(bray_binary_list)
rm(bray_binary_avg)
gc()

#log bray
bray_logCts_list = lapply(rarefactions_list, log_dist, method = "bray")
bray_logCts_avg = avg_matrix_list(bray_logCts_list) %>% as.dist()
saveRDS(bray_logCts_avg, "data/FEDRR_11062024/processed_tables/avg_dist/bray-logCounts.samps.rds")
rm(bray_logCts_list)
rm(bray_logCts_avg)
gc()



#alpha divs
shannon_list = lapply(rarefactions_list, diversity, index = "shannon")
simpson_list = lapply(rarefactions_list, diversity, index = "simpson")
richness_list = lapply(rarefactions_list, richness_calc)

div_avg = data.frame(
    sample = rownames(avg_matrix_list(shannon_list)),
    shannon = avg_matrix_list(shannon_list),
    simpson = avg_matrix_list(simpson_list),
    richness = avg_matrix_list(richness_list),
    stringsAsFactors = F
)
write.csv(div_avg, "data/FEDRR_11062024/processed_tables/avg_div/diversity.samps.csv", row.names = F)
rm(shannon_list)
rm(simpson_list)
rm(richness_list)
rm(div_avg)
gc()

###########################################
#Also average the counts table
avg_counts = avg_matrix_list(rarefactions_list)

rm(rarefactions_list)
gc()

write.csv(avg_counts, "data/FEDRR_11062024/processed_tables/asv_tab.rarefaction_avg.samps.csv")

rm(avg_counts)
gc()

############################################
############################################
############################################
############################################
############################################
# baseline
rarefactions_list = readRDS("data/FEDRR_11062024/processed_tables/rarefactions.baseline.rds")

#calc distance and alpha-div
## take avg
## save
## rm objs
## 
## binary bray
bray_binary_list = lapply(rarefactions_list, vegdist, method = "bray", binary = T, diag = T, upper = T)
bray_binary_avg = avg_matrix_list(bray_binary_list) %>% as.dist()
saveRDS(bray_binary_avg, "data/FEDRR_11062024/processed_tables/avg_dist/bray-binary.baseline.rds")
rm(bray_binary_list)
rm(bray_binary_avg)
gc()

#log bray
bray_logCts_list = lapply(rarefactions_list, log_dist, method = "bray")
bray_logCts_avg = avg_matrix_list(bray_logCts_list) %>% as.dist()
saveRDS(bray_logCts_avg, "data/FEDRR_11062024/processed_tables/avg_dist/bray-logCounts.baseline.rds")
rm(bray_logCts_list)
rm(bray_logCts_avg)
gc()



#alpha divs
shannon_list = lapply(rarefactions_list, diversity, index = "shannon")
simpson_list = lapply(rarefactions_list, diversity, index = "simpson")
richness_list = lapply(rarefactions_list, richness_calc)

div_avg = data.frame(
    sample = rownames(avg_matrix_list(shannon_list)),
    shannon = avg_matrix_list(shannon_list),
    simpson = avg_matrix_list(simpson_list),
    richness = avg_matrix_list(richness_list),
    stringsAsFactors = F
)
write.csv(div_avg, "data/FEDRR_11062024/processed_tables/avg_div/diversity.baseline.csv", row.names = F)
rm(shannon_list)
rm(simpson_list)
rm(richness_list)
rm(div_avg)
gc()

###########################################
#Also average the counts table
avg_counts = avg_matrix_list(rarefactions_list)

rm(rarefactions_list)
gc()

write.csv(avg_counts, "data/FEDRR_11062024/processed_tables/asv_tab.rarefaction_avg.baseline.csv")

rm(avg_counts)
gc()


############################################
############################################
############################################
############################################
############################################
# pcr
rarefactions_list = readRDS("data/FEDRR_11062024/processed_tables/rarefactions.pcr.rds")

#calc distance and alpha-div
## take avg
## save
## rm objs
## 
## binary bray
bray_binary_list = lapply(rarefactions_list, vegdist, method = "bray", binary = T, diag = T, upper = T)
bray_binary_avg = avg_matrix_list(bray_binary_list) %>% as.dist()
saveRDS(bray_binary_avg, "data/FEDRR_11062024/processed_tables/avg_dist/bray-binary.pcr.rds")
rm(bray_binary_list)
rm(bray_binary_avg)
gc()

#log bray
bray_logCts_list = lapply(rarefactions_list, log_dist, method = "bray")
bray_logCts_avg = avg_matrix_list(bray_logCts_list) %>% as.dist()
saveRDS(bray_logCts_avg, "data/FEDRR_11062024/processed_tables/avg_dist/bray-logCounts.pcr.rds")
rm(bray_logCts_list)
rm(bray_logCts_avg)
gc()



#alpha divs
shannon_list = lapply(rarefactions_list, diversity, index = "shannon")
simpson_list = lapply(rarefactions_list, diversity, index = "simpson")
richness_list = lapply(rarefactions_list, richness_calc)

div_avg = data.frame(
    sample = rownames(avg_matrix_list(shannon_list)),
    shannon = avg_matrix_list(shannon_list),
    simpson = avg_matrix_list(simpson_list),
    richness = avg_matrix_list(richness_list),
    stringsAsFactors = F
)
write.csv(div_avg, "data/FEDRR_11062024/processed_tables/avg_div/diversity.pcr.csv", row.names = F)
rm(shannon_list)
rm(simpson_list)
rm(richness_list)
rm(div_avg)
gc()

###########################################
#Also average the counts table
avg_counts = avg_matrix_list(rarefactions_list)

rm(rarefactions_list)
gc()

write.csv(avg_counts, "data/FEDRR_11062024/processed_tables/asv_tab.rarefaction_avg.pcr.csv")

rm(avg_counts)
gc()



############################################
############################################
############################################
############################################
############################################
# samps_and_baseline
rarefactions_list = readRDS("data/FEDRR_11062024/processed_tables/rarefactions.samps_and_baseline.rds")

#calc distance and alpha-div
## take avg
## save
## rm objs
## 
## binary bray
bray_binary_list = lapply(rarefactions_list, vegdist, method = "bray", binary = T, diag = T, upper = T)
bray_binary_avg = avg_matrix_list(bray_binary_list) %>% as.dist()
saveRDS(bray_binary_avg, "data/FEDRR_11062024/processed_tables/avg_dist/bray-binary.samps_and_baseline.rds")
rm(bray_binary_list)
rm(bray_binary_avg)
gc()

#log bray
bray_logCts_list = lapply(rarefactions_list, log_dist, method = "bray")
bray_logCts_avg = avg_matrix_list(bray_logCts_list) %>% as.dist()
saveRDS(bray_logCts_avg, "data/FEDRR_11062024/processed_tables/avg_dist/bray-logCounts.samps_and_baseline.rds")
rm(bray_logCts_list)
rm(bray_logCts_avg)
gc()



#alpha divs
shannon_list = lapply(rarefactions_list, diversity, index = "shannon")
simpson_list = lapply(rarefactions_list, diversity, index = "simpson")
richness_list = lapply(rarefactions_list, richness_calc)

div_avg = data.frame(
    sample = rownames(avg_matrix_list(shannon_list)),
    shannon = avg_matrix_list(shannon_list),
    simpson = avg_matrix_list(simpson_list),
    richness = avg_matrix_list(richness_list),
    stringsAsFactors = F
)
write.csv(div_avg, "data/FEDRR_11062024/processed_tables/avg_div/diversity.samps_and_baseline.csv", row.names = F)
rm(shannon_list)
rm(simpson_list)
rm(richness_list)
rm(div_avg)
gc()

###########################################
#Also average the counts table
avg_counts = avg_matrix_list(rarefactions_list)

rm(rarefactions_list)
gc()

write.csv(avg_counts, "data/FEDRR_11062024/processed_tables/asv_tab.rarefaction_avg.samps_and_baseline.csv")

rm(avg_counts)
gc()
