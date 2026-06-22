library(dplyr)

sites = read.csv("data/2025_metadata/database_sites_04202026.csv")
head(sites)
sites %>%
    group_by(State, Site.name) %>% 
    distinct() %>%
    group_by(State) %>%
    summarize(n = n())
