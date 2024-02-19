library(dplyr)
library(ggplot2)

lens = read.table("data/vol_test/itsx_test/R1.lens.txt")
lens_mod = read.table("data/vol_test/itsx_mod_test/R1.lens.txt")

comb_lens = left_join(lens[c("V1","V3")], lens_mod[c("V1","V3")], by = "V1") 
head(comb_lens)
cor(comb_lens$V3.x, comb_lens$V3.y)

ggplot(comb_lens, aes(V3.x, V3.y)) +
    geom_point()

lens = read.table("data/vol_test/itsx_test/R2.lens.txt")
lens_mod = read.table("data/vol_test/itsx_mod_test/R2.lens.txt")

comb_lens = left_join(lens[c("V1","V3")], lens_mod[c("V1","V3")], by = "V1") 
head(comb_lens)

ggplot(comb_lens, aes(V3.x, V3.y)) +
    geom_point()

cor(comb_lens$V3.x, comb_lens$V3.y)
