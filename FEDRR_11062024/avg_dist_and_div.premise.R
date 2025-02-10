library(vegan)
library(dplyr)
source("~/repo/BIL_pathogens/library/library.R")


############################################
############################################
# samps
rarefactions_list = readRDS("~/FEDRR_all_2024/rarefactions/rarefactions.samps.rds")

#calc distance and alpha-div
## take avg
## save
## rm objs
## 
## binary bray
bray_binary_list = lapply(rarefactions_list, vegdist, method = "bray", binary = T, diag = T, upper = T)
bray_binary_avg = avg_matrix_list(bray_binary_list) %>% as.dist()
saveRDS(bray_binary_avg, "~/FEDRR_all_2024/rarefactions/bray-binary.samps.rds")
rm(bray_binary_list)
rm(bray_binary_avg)
gc()

#log bray
bray_logCts_list = lapply(rarefactions_list, log_dist, method = "bray")
bray_logCts_avg = avg_matrix_list(bray_logCts_list) %>% as.dist()
saveRDS(bray_logCts_avg, "~/FEDRR_all_2024/rarefactions/bray-logCounts.samps.rds")
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
write.csv(div_avg, "~/FEDRR_all_2024/rarefactions/diversity.samps.csv", row.names = F)
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

write.csv(avg_counts, "~/FEDRR_all_2024/rarefactions/asv_tab.rarefaction_avg.samps.csv")

rm(avg_counts)
gc()

############################################
############################################
############################################
############################################
############################################
# baseline
rarefactions_list = readRDS("~/FEDRR_all_2024/rarefactions/rarefactions.baseline.rds")

#calc distance and alpha-div
## take avg
## save
## rm objs
## 
## binary bray
bray_binary_list = lapply(rarefactions_list, vegdist, method = "bray", binary = T, diag = T, upper = T)
bray_binary_avg = avg_matrix_list(bray_binary_list) %>% as.dist()
saveRDS(bray_binary_avg, "~/FEDRR_all_2024/rarefactions/bray-binary.baseline.rds")
rm(bray_binary_list)
rm(bray_binary_avg)
gc()

#log bray
bray_logCts_list = lapply(rarefactions_list, log_dist, method = "bray")
bray_logCts_avg = avg_matrix_list(bray_logCts_list) %>% as.dist()
saveRDS(bray_logCts_avg, "~/FEDRR_all_2024/rarefactions/bray-logCounts.baseline.rds")
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
write.csv(div_avg, "~/FEDRR_all_2024/rarefactions/diversity.baseline.csv", row.names = F)
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

write.csv(avg_counts, "~/FEDRR_all_2024/rarefactions/asv_tab.rarefaction_avg.baseline.csv")

rm(avg_counts)
gc()


############################################
############################################
############################################
############################################
############################################
# pcr
rarefactions_list = readRDS("~/FEDRR_all_2024/rarefactions/rarefactions.pcr.rds")

#calc distance and alpha-div
## take avg
## save
## rm objs
## 
## binary bray
bray_binary_list = lapply(rarefactions_list, vegdist, method = "bray", binary = T, diag = T, upper = T)
bray_binary_avg = avg_matrix_list(bray_binary_list) %>% as.dist()
saveRDS(bray_binary_avg, "~/FEDRR_all_2024/rarefactions/bray-binary.pcr.rds")
rm(bray_binary_list)
rm(bray_binary_avg)
gc()

#log bray
bray_logCts_list = lapply(rarefactions_list, log_dist, method = "bray")
bray_logCts_avg = avg_matrix_list(bray_logCts_list) %>% as.dist()
saveRDS(bray_logCts_avg, "~/FEDRR_all_2024/rarefactions/bray-logCounts.pcr.rds")
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
write.csv(div_avg, "~/FEDRR_all_2024/rarefactions/diversity.pcr.csv", row.names = F)
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

write.csv(avg_counts, "~/FEDRR_all_2024/rarefactions/asv_tab.rarefaction_avg.pcr.csv")

rm(avg_counts)
gc()



############################################
############################################
############################################
############################################
############################################
# samps_and_baseline
rarefactions_list = readRDS("~/FEDRR_all_2024/rarefactions/rarefactions.samps_and_baseline.rds")

#calc distance and alpha-div
## take avg
## save
## rm objs
## 
## binary bray
bray_binary_list = lapply(rarefactions_list, vegdist, method = "bray", binary = T, diag = T, upper = T)
bray_binary_avg = avg_matrix_list(bray_binary_list) %>% as.dist()
saveRDS(bray_binary_avg, "~/FEDRR_all_2024/rarefactions/bray-binary.samps_and_baseline.rds")
rm(bray_binary_list)
rm(bray_binary_avg)
gc()

#log bray
bray_logCts_list = lapply(rarefactions_list, log_dist, method = "bray")
bray_logCts_avg = avg_matrix_list(bray_logCts_list) %>% as.dist()
saveRDS(bray_logCts_avg, "~/FEDRR_all_2024/rarefactions/bray-logCounts.samps_and_baseline.rds")
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
write.csv(div_avg, "~/FEDRR_all_2024/rarefactions/diversity.samps_and_baseline.csv", row.names = F)
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

write.csv(avg_counts, "~/FEDRR_all_2024/rarefactions/asv_tab.rarefaction_avg.samps_and_baseline.csv")

rm(avg_counts)
gc()
