library(dplyr)
library(ggplot2)

lens = read.table("data/vol_test/all.lens.txt")
lens_mod = read.table("data/vol_test/all_mod.lens.txt")

comb_lens = left_join(lens, lens_mod, by = "V1") 
head(comb_lens)

ggplot(comb_lens, aes(V2.x, V2.y)) +
    geom_point()
