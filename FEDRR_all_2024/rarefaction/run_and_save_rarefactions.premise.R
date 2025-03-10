library(dplyr)
source("~/repo/BIL_pathogens/library/library.R")

#asv_tab.baseline = read.csv("~/FEDRR_all_2024/rarefactions/asv_tab.baseline.csv", header = T, row.names = 1)
#asv_tab.pcr = read.csv("~/FEDRR_all_2024/rarefactions/asv_tab.pcr.csv", header = T, row.names = 1)
asv_tab.samps = read.csv("~/FEDRR_all_2024/rarefactions/asv_tab.funnel_traps_and_spore_traps.csv", header = T, row.names = 1)
#asv_tab.samps_and_baseline = read.csv("~/FEDRR_all_2024/rarefactions/asv_tab.samps_and_baseline.csv", header = T, row.names = 1)



#min_seqs.baseline = min(rowSums(asv_tab.baseline))
#rarefactions_list.baseline = multiple_subsamples(x = asv_tab.baseline, depth = min_seqs.baseline, iterations = 1000)
#saveRDS(rarefactions_list.baseline, "~/FEDRR_all_2024/rarefactions/rarefactions.baseline.rds")

#min_seqs.pcr = min(rowSums(asv_tab.pcr))
#rarefactions_list.pcr = multiple_subsamples(x = asv_tab.pcr, depth = min_seqs.pcr, iterations = 1000)
#saveRDS(rarefactions_list.pcr, "~/FEDRR_all_2024/rarefactions/rarefactions.pcr.rds")

min_seqs.samps = min(rowSums(asv_tab.samps)) 
rarefactions_list.samps = multiple_subsamples(x = asv_tab.samps, depth = min_seqs.samps, iterations = 1000) 
saveRDS(rarefactions_list.samps, "~/FEDRR_all_2024/rarefactions/rarefactions.funnel_traps_and_spore_traps.rds")

#min_seqs.samps_and_baseline = min(rowSums(asv_tab.samps_and_baseline))
#rarefactions_list.samps_and_baseline = multiple_subsamples(x = asv_tab.samps_and_baseline, depth = min_seqs.samps_and_baseline, iterations = 1000)
#saveRDS(rarefactions_list.samps_and_baseline, "~/FEDRR_all_2024/rarefactions/rarefactions.samps_and_baseline.rds")
