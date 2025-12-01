library(dplyr)
source("~/repo/BIL_pathogens/library/library.R")

asv_tab.samps = read.csv("~/FEDRR_project/FEDRR_2025_11182025/rarefactions/asv_tab.samps.csv", header = T, row.names = 1)
asv_tab.samps_and_baseline = read.csv("~/FEDRR_project/FEDRR_2025_11182025/rarefactions/asv_tab.samps_and_baseline.csv", header = T, row.names = 1)
asv_tab.samps_and_pcrNeg = read.csv("~/FEDRR_project/FEDRR_2025_11182025/rarefactions/asv_tab.samps_and_pcrNeg.csv", header = T, row.names = 1)
asv_tab.samps_and_baseline.NH = read.csv("~/FEDRR_project/FEDRR_2025_11182025/rarefactions/asv_tab.samps_and_baseline.NH.csv", header = T, row.names = 1)
asv_tab.baseline.NH = read.csv("~/FEDRR_project/FEDRR_2025_11182025/rarefactions/asv_tab.baseline.NH.csv", header = T, row.names = 1)
asv_tab.baseline = read.csv("~/FEDRR_project/FEDRR_2025_11182025/rarefactions/asv_tab.baseline.csv", header = T, row.names = 1)




min_seqs.samps = min(rowSums(asv_tab.samps))
rarefactions_list.samps = multiple_subsamples(x = asv_tab.samps, depth = min_seqs.samps, iterations = 1000)
saveRDS(rarefactions_list.samps, "~/FEDRR_project/FEDRR_2025_11182025/rarefactions/rarefactions.samps.rds")
rm(rarefactions_list.samps)
gc()

min_seqs.samps_and_baseline = min(rowSums(asv_tab.samps_and_baseline))
rarefactions_list.samps_and_baseline = multiple_subsamples(x = asv_tab.samps_and_baseline, depth = min_seqs.samps_and_baseline, iterations = 1000)
saveRDS(rarefactions_list.samps_and_baseline, "~/FEDRR_project/FEDRR_2025_11182025/rarefactions/rarefactions.samps_and_baseline.rds")
rm(rarefactions_list.samps_and_baseline)
gc()


min_seqs.samps_and_pcrNeg = min(rowSums(asv_tab.samps_and_pcrNeg))
rarefactions_list.samps_and_pcrNeg = multiple_subsamples(x = asv_tab.samps_and_pcrNeg, depth = min_seqs.samps_and_pcrNeg, iterations = 1000)
saveRDS(rarefactions_list.samps_and_pcrNeg, "~/FEDRR_project/FEDRR_2025_11182025/rarefactions/rarefactions.samps_and_pcrNeg.rds")
rm(rarefactions_list.samps_and_pcrNeg)
gc()


min_seqs.samps_and_baseline.NH = min(rowSums(asv_tab.samps_and_baseline.NH))
rarefactions_list.samps_and_baseline.NH = multiple_subsamples(x = asv_tab.samps_and_baseline.NH, depth = min_seqs.samps_and_baseline.NH, iterations = 1000)
saveRDS(rarefactions_list.samps_and_baseline.NH, "~/FEDRR_project/FEDRR_2025_11182025/rarefactions/rarefactions.samps_and_baseline.NH.rds")
rm(rarefactions_list.samps_and_baseline.NH)
gc()

min_seqs.baseline.NH = min(rowSums(asv_tab.baseline.NH))
rarefactions_list.baseline.NH = multiple_subsamples(x = asv_tab.baseline.NH, depth = min_seqs.baseline.NH, iterations = 1000)
saveRDS(rarefactions_list.baseline.NH, "~/FEDRR_project/FEDRR_2025_11182025/rarefactions/rarefactions.baseline.NH.rds")
rm(rarefactions_list.baseline.NH)
gc()

min_seqs.baseline = min(rowSums(asv_tab.baseline))
rarefactions_list.baseline = multiple_subsamples(x = asv_tab.baseline, depth = min_seqs.baseline, iterations = 1000)
saveRDS(rarefactions_list.baseline, "~/FEDRR_project/FEDRR_2025_11182025/rarefactions/rarefactions.baseline.rds")
rm(rarefactions_list.baseline)
gc()

#min_seqs.samps = min(rowSums(asv_tab.samps))
#rarefactions_list.samps = multiple_subsamples(x = asv_tab.samps, depth = min_seqs.samps, iterations = 1000)
#saveRDS(rarefactions_list.samps, "~/FEDRR_project/FEDRR_2025_11182025/rarefactions/rarefactions.funnel_traps_and_spore_traps.rds")
