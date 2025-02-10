library(dplyr)

metadata.NY = read.csv("data/2024_metadata/collated/metadata.NY.csv")
metadata.OH = read.csv("data/2024_metadata/collated/metadata.OH.csv")
metadata.IA = read.csv("data/2024_metadata/collated/metadata.IA.csv")
metadata = rbind(metadata.NY, metadata.OH, metadata.IA)
ncol(metadata)
id_mapping = read.csv("data/2024_metadata/collated/sampleID_mapping_NY-OH-IA_02102025.csv")
reads_per_sample = read.table("data/2024_metadata/collated/reads_per_sample_NY-OH-IA_02102025.txt", header = F)
colnames(reads_per_sample) = c("SequenceID", "num_reads")
head(reads_per_sample)


full_metadata = full_join(id_mapping, metadata, by = "trapID") %>%
    full_join(., reads_per_sample, by = "SequenceID")
head(full_metadata)

sampleID.split = strsplit(full_metadata$sampleID, " ", fixed = T)
full_metadata$dateOrType = lapply(sampleID.split, function(x) x[2]) %>% unlist

write.csv(full_metadata, "data/2024_metadata/collated/collated_metadata_NY-OH-IA_02102025.csv", row.names = F)
