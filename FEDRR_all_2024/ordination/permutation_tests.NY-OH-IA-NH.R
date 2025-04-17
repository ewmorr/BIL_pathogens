library(dplyr)
library(vegan)
library(lubridate)
set.seed(12345)

metadata = read.csv("data/2024_metadata/collated/collated_metadata_NY-OH-IA-NH_02102025.csv")
head(metadata)

metadata$Lure %>% unique
metadata$trapID %>% unique
metadata %>% filter(State == "NY") %>% pull(trapID) %>% unique
metadata$BaselineType %>% unique()

metadata$BaselineType = metadata$dateOrBaselineType
metadata$BaselineType[grep("^[0-9].*", metadata$BaselineType)] = "sample"

metadata$date = paste0(metadata$dateOrBaselineType, "-2024")
metadata$date[grep("^BL.", metadata$date)] = NA
metadata$date[grep("^negative.*", metadata$date)] = NA
metadata$mdy = dmy(metadata$date)

########################################
########################################
########################################
#samps and baseline dist
#
bray_bin.dist = readRDS("data/FEDRR_all_2024/rarefaction_avgs/bray-binary.samps_and_baseline.rds")

###########################################
# Setting up permutations and perform tests

#checking that metadata matches dist obj for adonis2
length(rownames(as.matrix(bray_bin.dist))) == nrow(metadata)
rownames(metadata) = metadata$SequenceID

metadata.ordered = metadata[rownames(as.matrix(bray_bin.dist)),]
nrow(metadata.ordered)
sum(rownames(metadata.ordered) == rownames(as.matrix(bray_bin.dist)))

#for testing an effect that varies within subject (like baseline or not OR date)
h1 = with(metadata.ordered, how(nperm = 999, blocks = trapID))
# the above is equivalent to setting strata in adonis(2)
# 
# 
# supposedly this is the correct way to do it using the condition variable, see https://fromthebottomoftheheap.net/slides/advanced-vegan-webinar-2020/advanced-vegan#72
baseline.dbrda = dbrda(bray_bin.dist ~ BaselineType + Condition(trapID), permutations = h1, data = metadata.ordered)
anova(baseline.dbrda, permutations = h1, by ="terms", parallel = 4)
#              Df SumOfSqs      F Pr(>F)    
#BaselineType   2    4.415 6.3107  0.001 ***
#Residual     368  128.722                                        
#
anova(baseline.dbrda, permutations = h1, by ="terms", parallel = 4, model = "direct") #margins doesn't accomplish anmything with only one fixed effect
#
baseline.dbrda = dbrda(bray_bin.dist ~ BaselineType, permutations = h1, data = metadata.ordered)
anova(baseline.dbrda, permutations = h1, by ="terms")
#              Df SumOfSqs      F Pr(>F)    
#BaselineType   2    3.926 4.5619  0.001 ***
#Residual     442  190.204                  
adonis2(bray_bin.dist ~ BaselineType, permutations = h1, data = metadata.ordered)      
#          Df SumOfSqs      R2      F Pr(>F)    
#Model      2    3.926 0.02022 4.5619  0.001 ***
#Residual 442  190.204 0.97978                  
#Total    444  194.130 1.00000        


########################################
########################################
########################################
#main tests (date, lure, risk)
#
bray_bin.dist = readRDS("data/FEDRR_all_2024/rarefaction_avgs/bray-binary.samps.rds")

###########################################
# Setting up permutations and perform tests

#checking that metadata matches dist obj for adonis2
length(rownames(as.matrix(bray_bin.dist))) == nrow(metadata)
rownames(metadata) = metadata$SequenceID

metadata.ordered = metadata[rownames(as.matrix(bray_bin.dist)),]
nrow(metadata.ordered)
sum(rownames(metadata.ordered) == rownames(as.matrix(bray_bin.dist)))

metadata.ordered$mdy

#date and lure
#
metadata.ordered %>% colnames
h1 = with(metadata.ordered,
    how(
        nperm = 999,
        blocks = Site, 
        plots = Plots(strata = trapID, type = "free"),
        within = Within(type = "series")
    )
)
h1.dbrda = dbrda(bray_bin.dist ~ Lure + mdy + Condition(Site), data = metadata.ordered)
anova(h1.dbrda, permutations = h1)
# Error in check(sn, control = control, quietly = quietly) : 
# Design must be balanced if permuting 'strata'.

##########################################
# JUST DATE. COMPARE TO ADONIS2 BELOW
h1 = with(metadata.ordered,
          how(
              nperm = 999,
              blocks = trapID, 
              within = Within(type = "series")
          )
)
h1.dbrda = dbrda(bray_bin.dist ~ mdy + Condition(trapID), data = metadata.ordered)
anova(h1.dbrda, permutations = h1)
#Model: dbrda(formula = bray_bin.dist ~ mdy + Condition(trapID), data = metadata.ordered)
#          Df SumOfSqs      F Pr(>F)    
#Model      1     4.72 14.627  0.001 ***
#Residual 327   105.53                  
#---
4.72/sum(4.72,105.53)
# 0.04281179 
# this agrees with some of the adonis models below testing for mdy non nested

#########################################
# subsampling table for tests
metadata.ordered %>%
    group_by(State, Site, Lure, mdy) %>%
    summarize(n = n()) %>%
    filter(n > 1)

metadata.ordered %>%
    group_by(trapID, State, Site, Lure) %>%
    summarize(n = n()) %>%
    print(n = Inf)

metadata.ordered %>%
    group_by(trapID, State, Site, Lure) %>%
    summarize(n = n()) %>%
    filter(n <= 4)

# for test of Lure we are permuting traps within sites.
# five of the IA sites have n == 4
# downsample everything to 4 (throw out one IA site with 3)

min_n = 4

metadata.ordered %>%
    group_by(State, Site, Lure, trapID) %>%
    summarize(n = n()) %>%
    as.data.frame() -> trap_n
metadata.ordered %>%
    group_by(State, Site, Lure, trapID) %>%
    summarize(n = n()) %>%
    filter(n < min_n) %>%
    pull(trapID) -> low_n_trap
low_n_trap

lure_test.samples_exclude = metadata.ordered %>% filter(trapID %in% low_n_trap) %>% pull(SequenceID)
for(i in 1:nrow(trap_n)){
    if(trap_n$n[i] > min_n){
        tmpIDS = metadata.ordered %>% filter(trapID == trap_n$trapID[i]) %>% pull(SequenceID)
        n_samples_to_rm = trap_n$n[i] - min_n
        lure_test.samples_exclude = c(lure_test.samples_exclude, sample(tmpIDS, n_samples_to_rm, replace = F))
    }
}
length(lure_test.samples_exclude) 


bray_bin.mat = as.matrix(bray_bin.dist)
nrow(bray_bin.mat)
bray_bin.mat.lure_test = bray_bin.mat[!rownames(bray_bin.mat) %in% lure_test.samples_exclude, !colnames(bray_bin.mat) %in% lure_test.samples_exclude]
nrow(bray_bin.mat.lure_test)
bray_bin.dist.lure_test = as.dist(bray_bin.mat.lure_test)
metadata.lure_test = metadata.ordered[rownames(bray_bin.mat.lure_test),]
metadata.lure_test %>%
    group_by(trapID, State, Site, Lure) %>%
    summarize(n = n()) %>%
    print(n = Inf)
h2 = with(metadata.lure_test,
    how(
        nperm = 999,
        blocks = Site,
        plots = Plots(strata = trapID,type = "free"),
        within = Within(type = "none")
    )
)
h2.dbrda = dbrda(bray_bin.dist.lure_test ~ Lure + Condition(Site), data = metadata.lure_test)
anova(h2.dbrda, permutations = h2)
#Model      2    0.992 1.4488  0.006 **
#Residual 269   92.124                 
anova(h2.dbrda, permutations = h2, by ="terms")
#          Df SumOfSqs      F Pr(>F)   
#Lure       2    0.992 1.4488  0.002 **
#Residual 269   92.124                 
0.992/sum(0.992,92.124)
#0.01065338
# this agrees with some of the adonis models below testing for Lure non nested
h2.dbrda = dbrda(bray_bin.dist.lure_test ~ trapID + Condition(Site), data = metadata.lure_test)

#for Risk we are permuting sites within states (block = state, strata = plots, within = none
metadata.ordered %>%
    group_by(State, Site) %>%
    summarize(n = n()) %>%
    print(n = Inf)
# we don't want to throw out a whole site so we will downsample everything to the lowest n (12)
metadata.ordered %>%
    group_by(State, Site) %>%
    summarize(n = n()) -> site_n

n_samples = 12
#just testing the indexing
for(i in 1:10){
    print(paste(i*n_samples-n_samples+1, i*n_samples))
}
#
risk_test.samples_include = vector(mode = "character", length = n_samples*nrow(site_n))
for(i in 1:nrow(site_n)){
    tmpIDS = metadata.ordered %>% filter(Site == site_n$Site[i]) %>% pull(SequenceID)
    risk_test.samples_include[(i*n_samples-n_samples+1):(i*n_samples)] = sample(tmpIDS, n_samples, replace = F)
}

bray_bin.mat.risk_test = bray_bin.mat[rownames(bray_bin.mat) %in% risk_test.samples_include, colnames(bray_bin.mat) %in% risk_test.samples_include]
bray_bin.dist.risk_test = as.dist(bray_bin.mat.risk_test)
metadata.risk_test = metadata.ordered[rownames(bray_bin.mat.risk_test),]
nrow(bray_bin.mat.risk_test) == nrow(metadata.risk_test)

h3 = with(metadata.risk_test,
          how(
              nperm = 999,
              blocks = State,
              plots = Plots(strata = Site,type = "free"),
              within = Within(type = "none")
          )
)
h3.dbrda = dbrda(bray_bin.dist.risk_test ~ Risk + Condition(State), data = metadata.risk_test)
anova(h3.dbrda, permutations = h3)
factorial(6)^2
#          Df SumOfSqs      F Pr(>F)
#Model      1     0.83 2.1481  0.707
#Residual 295   114.04              
0.83/sum(0.83,114.04)
#0.007225559
# this r^2 agrees with some of the adonis models below testing for risk non nested


metadata.risk_test %>% group_by(State, Site, Risk) %>% summarize(n = n())

# to get a sense of variance explained by location we test trap within site
# and site within state
h4.dbrda = dbrda(bray_bin.dist.lure_test ~ trapID, data = metadata.lure_test)
anova(h4.dbrda, permutations = h2)
#Model     73   52.994 2.1644      1
#Residual 222   74.459              
52.994/sum(52.994,74.459)
# 0.4157925
h5.dbrda = dbrda(bray_bin.dist.risk_test ~ Site, data = metadata.risk_test)
anova(h5.dbrda, permutations = h3)
#          Df SumOfSqs      F Pr(>F)
#Model     24   35.166 4.2855      1
#Residual 275   94.024              
35.166/sum(35.166,94.024)
#0.2722037


adonis2(bray_bin.dist ~ State/Site/trapID, data = metadata.ordered, by = "terms", permutations = 1999)
#                   Df SumOfSqs      R2       F Pr(>F)    
#State               3   19.101 0.11027 18.9431  5e-04 ***
#State:Site         21   24.721 0.14272  3.5023  5e-04 ***
#State:Site:trapID  50   19.151 0.11056  1.1395  5e-04 ***
#Residual          328  110.246 0.63645                   
#Total             402  173.219 1.00000                   


###############################################################################

adonis2(bray_bin.dist ~ Lure + mdy, strata= metadata.ordered$Site, data = metadata.ordered)
#          Df SumOfSqs     R2      F Pr(>F)    
#Model      3    6.479 0.0374 5.1682  0.001 ***
#Residual 399  166.739 0.9626                  
#Total    402  173.219 1.0000                  

adonis2(bray_bin.dist ~ Lure + mdy, strata= metadata.ordered$Site, data = metadata.ordered, by = "terms")
#          Df SumOfSqs     R2      F Pr(>F)    
#Lure       2    1.129 0.00652  1.351  0.001 ***
#mdy        1    5.350 0.03089 12.803  0.001 ***
#Residual 399  166.739 0.96260                  
#Total    402  173.219 1.00000                  

adonis2(bray_bin.dist ~ Lure + mdy, strata= metadata.ordered$Site, data = metadata.ordered, by = "margin")
#          Df SumOfSqs      R2       F Pr(>F)    
#Lure       2     1.13 0.00652  1.3517  0.001 ***
#mdy        1     5.35 0.03089 12.8026  0.001 ***
#Residual 399   166.74 0.96260                   
#Total    402   173.22 1.00000                   

adonis2(bray_bin.dist ~ Site + Lure + mdy, strata = metadata.ordered$State, data = metadata.ordered)
#          Df SumOfSqs      R2      F Pr(>F)    
#Model     27   49.775 0.28736 5.6004  0.001 ***
#Residual 375  123.443 0.71264                  
#Total    402  173.219 1.00000                  

adonis2(bray_bin.dist ~ Site/Lure/mdy, strata = metadata.ordered$State, data = metadata.ordered) ####################
#Model    149  101.667 0.58693 2.4127  0.001 ***
#Residual 253   71.551 0.41307                  
#Total    402  173.219 1.00000                  

adonis2(bray_bin.dist ~ Site + Lure + mdy + Site/Lure/mdy, strata = metadata.ordered$State, data = metadata.ordered)
#          Df SumOfSqs      R2      F Pr(>F)    
#Model    149  101.667 0.58693 2.4127  0.001 ***
#Residual 253   71.551 0.41307                  
#Total    402  173.219 1.00000                  

adonis2(bray_bin.dist ~ Site/Lure/mdy, strata = metadata.ordered$State, data = metadata.ordered, by = "terms")
#               Df SumOfSqs      R2      F Pr(>F)    
#Site           24   43.822 0.25299 6.4563  0.001 ***
#Site:Lure      50   19.151 0.11056 1.3543  0.001 ***
#Site:Lure:mdy  75   38.694 0.22338 1.8243  0.001 ***
#Residual      253   71.551 0.41307                  
#Total         402  173.219 1.00000                  

adonis2(bray_bin.dist ~ Site + Lure + mdy + Site/Lure/mdy, strata = metadata.ordered$State, data = metadata.ordered, by = "terms")
#               Df SumOfSqs      R2       F Pr(>F)    
#Site           24   43.822 0.25299  6.4563  0.001 ***
#Lure            2    1.120 0.00647  1.9809  0.001 ***
#mdy             1    4.833 0.02790 17.0888  0.001 ***
#Site:Lure      48   17.918 0.10344  1.3199  0.001 ***
#Site:Lure:mdy  74   33.974 0.19613  1.6234  0.001 ***
#Residual      253   71.551 0.41307                   
#Total         402  173.219 1.00000                   
#
#Running this second way pulls a bit of variance out of Site:Lure and a bit more out of Site:Lure:mdy. 
# These can be attributed to the main effect presumably, i.e., consistent variation in these across the dataset?
# 
# 

adonis2(bray_bin.dist ~ Site/trapID/mdy + Risk + Lure + mdy, strata = metadata.ordered$State, data = metadata.ordered, by = "terms")

adonis2(bray_bin.dist ~ Risk + Lure + Site/trapID/mdy + mdy, strata = metadata.ordered$State, data = metadata.ordered, by = "terms")
#                 Df SumOfSqs      R2       F Pr(>F)    
#Risk              1    1.014 0.00585  3.5860  0.001 ***
#Lure              2    1.129 0.00652  1.9957  0.001 ***
#Site             23   42.800 0.24708  6.5798  0.001 ***
#mdy               1    4.833 0.02790 17.0888  0.001 ***
#Site:trapID      48   17.918 0.10344  1.3199  0.001 ***
#Site:trapID:mdy  74   33.974 0.19613  1.6234  0.001 ***
#Residual        253   71.551 0.41307                   
#Total           402  173.219 1.00000                   

#interstingly we get a sig effect of risk level here but the correct way to perform the sig test is using the restricted perms above

# need to do some more digging to fully understand the nested effects structure and if we should use the r2 values from one of the nested adonis models
