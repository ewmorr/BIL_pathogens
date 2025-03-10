library(dplyr)
library(vegan)
library(lubridate)
set.seed(12345)

metadata = read.csv("data/2024_metadata/collated/collated_metadata_NY-OH-IA_02102025.csv")
head(metadata)

########################################
########################################
########################################
#samps and baseline dist
#
bray_bin.dist = readRDS("data/FEDRR_all_2024/avg_dist/bray-binary.samps_and_baseline.rds")

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
#BaselineType   1    2.619 7.2497  0.001 ***
#Residual     273   98.620                                   
#
anova(baseline.dbrda, permutations = h1, by ="terms", parallel = 4, model = "direct") #margins doesn't accomplish anmything with only one fixed effect
#
baseline.dbrda = dbrda(bray_bin.dist ~ BaselineType, permutations = h1, data = metadata.ordered)
anova(baseline.dbrda, permutations = h1, by ="terms")
#              Df SumOfSqs      F Pr(>F)    
#BaselineType   1    2.584 5.9668  0.001 ***
#Residual     329  142.478                            
adonis2(bray_bin.dist ~ BaselineType, permutations = h1, data = metadata.ordered)
#          Df SumOfSqs      R2      F Pr(>F)    
#Model      1    2.584 0.01781 5.9668  0.001 ***
#Residual 329  142.478 0.98219                  
#Total    330  145.062 1.00000                    


########################################
########################################
########################################
#main tests (date, lure, risk)
#
bray_bin.dist = readRDS("data/FEDRR_all_2024/avg_dist/bray-binary.samps.rds")

###########################################
# Setting up permutations and perform tests

#checking that metadata matches dist obj for adonis2
length(rownames(as.matrix(bray_bin.dist))) == nrow(metadata)
rownames(metadata) = metadata$SequenceID

metadata.ordered = metadata[rownames(as.matrix(bray_bin.dist)),]
nrow(metadata.ordered)
sum(rownames(metadata.ordered) == rownames(as.matrix(bray_bin.dist)))

metadata.ordered$mdy = mdy(metadata.ordered$date)

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
#Df SumOfSqs      F Pr(>F)    
#Model      1    4.671 14.076  0.001 ***
# Residual 245   81.303   
#---
4.671/sum(4.671,81.303)
# 0.05433038 
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
# for test of Lure we are permuting traps within sites.
# should throw out the n of <=4 and downsample the rest to five (at random)
metadata.ordered %>%
    group_by(State, Site, Lure, trapID) %>%
    summarize(n = n()) %>%
    as.data.frame() -> trap_n
metadata.ordered %>%
    group_by(State, Site, Lure, trapID) %>%
    summarize(n = n()) %>%
    filter(n <= 4) %>%
    pull(trapID) -> low_n_trap

lure_test.samples_exclude = metadata.ordered %>% filter(trapID %in% low_n_trap) %>% pull(SequenceID)
for(i in 1:nrow(trap_n)){
    if(trap_n$n[i] > 5){
        tmpIDS = metadata.ordered %>% filter(trapID == trap_n$trapID[i]) %>% pull(SequenceID)
        lure_test.samples_exclude = c(lure_test.samples_exclude, sample(tmpIDS, 1, replace = F))
    }
}
length(lure_test.samples_exclude) == (trap_n %>% filter(n == 6) %>% nrow)+3


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
#Model      2    1.037 1.4639  0.002 **
#Residual 224   79.331  
anova(h2.dbrda, permutations = h2, by ="terms")
#          Df SumOfSqs      F Pr(>F)   
#Lure       2    1.037 1.4639  0.001 ***
#Residual 224   79.331                               
1.037/sum(1.037,79.331)
#0.01290315
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

#just testing the indexing
for(i in 1:10){
    print(paste(i*12-12+1, i*12))
}
#
risk_test.samples_include = vector(mode = "character", length = 12*nrow(site_n))
for(i in 1:nrow(site_n)){
    tmpIDS = metadata.ordered %>% filter(Site == site_n$Site[i]) %>% pull(SequenceID)
    risk_test.samples_include[(i*12-12+1):(i*12)] = sample(tmpIDS, 12, replace = F)
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
#Model      1    0.989 2.5048  0.379
#Residual 224   88.486      
0.989/sum(0.989,88.486)
#0.01105337
# this r^2 agrees with some of the adonis models below testing for risk non nested


metadata.risk_test %>% group_by(State, Site, Risk) %>% summarize(n = n())

# to get a sense of variance explained by location we test trap within site
# and site within state
h4.dbrda = dbrda(bray_bin.dist.lure_test ~ trapID, data = metadata.lure_test)
anova(h4.dbrda, permutations = h2)
#Model     48   37.255 2.2125      1
#Residual 196   68.758      
37.255/sum(37.255,68.758)
# 0.3514192
h5.dbrda = dbrda(bray_bin.dist.risk_test ~ Site, data = metadata.risk_test)
anova(h5.dbrda, permutations = h3)
#          Df SumOfSqs      F Pr(>F)
#Model     18   24.138 3.7714      1
#Residual 209   74.314              
24.138/sum(24.138,74.314)
#0.2451753


adonis2(bray_bin.dist ~ State/Site/trapID, data = metadata.ordered, by = "terms", permutations = 1999)
Df SumOfSqs      R2       F Pr(>F)    
#State               2   11.465 0.08755 16.4024  5e-04 ***
#State:Site         16   18.121 0.13838  3.2407  5e-04 ***
#State:Site:trapID  38   15.394 0.11755  1.1591  5e-04 ***


###############################################################################

adonis2(bray_bin.dist ~ Lure + mdy, strata= metadata.ordered$Site, data = metadata.ordered)
#          Df SumOfSqs     R2      F Pr(>F)    
#Model      3     6.67 0.05093 5.3489  0.001 ***
#Residual 299   124.28 0.94907                  
#Total    302   130.95 1.00000                  
adonis2(bray_bin.dist ~ Lure + mdy, strata= metadata.ordered$Site, data = metadata.ordered, by = "terms")
#Df SumOfSqs      R2      F Pr(>F)    
#Lure       2    1.239 0.00946  1.4905  0.001 ***
#mdy        1    5.431 0.04147 13.0659  0.001 ***
#Residual 299  124.284 0.94907                   
#Total    302  130.954 1.00000                   
#--
adonis2(bray_bin.dist ~ Lure + mdy, strata= metadata.ordered$Site, data = metadata.ordered, by = "margin")
#          Df SumOfSqs      R2       F Pr(>F)    
#Lure       2    1.244 0.00950  1.4969  0.001 ***
#mdy        1    5.431 0.04147 13.0659  0.001 ***
#Residual 299  124.284 0.94907                   
#Total    302  130.954 1.00000                   

adonis2(bray_bin.dist ~ Site + Lure + mdy, strata = metadata.ordered$State, data = metadata.ordered)
#          Df SumOfSqs      R2      F Pr(>F)    
#Model     14   22.383 0.25947 4.8052  0.001 ***
#Residual 192   63.882 0.74053                  
#Total    206   86.265 1.00000                  
adonis2(bray_bin.dist ~ Site/Lure/mdy, strata = metadata.ordered$State, data = metadata.ordered)
#Model     21   35.585 0.27174 4.9929  0.001 ***
#Residual 281   95.369 0.72826                  
#Total    302  130.954 1.00000                  
adonis2(bray_bin.dist ~ Site + Lure + mdy + Site/Lure/mdy, strata = metadata.ordered$State, data = metadata.ordered)
#          Df SumOfSqs      R2      F Pr(>F)    
#Model    113   75.402 0.57579 2.2702  0.001 ***
#Residual 189   55.552 0.42421                  
#Total    302  130.954 1.00000                  
adonis2(bray_bin.dist ~ Site/Lure/mdy, strata = metadata.ordered$State, data = metadata.ordered, by = "terms")
#               Df SumOfSqs      R2      F Pr(>F)    
#Site           18   29.586 0.22593 5.5922  0.001 ***
#Site:Lure      38   15.394 0.11755 1.3782  0.001 ***
#Site:Lure:mdy  57   30.422 0.23231 1.8158  0.001 ***
#Residual      189   55.552 0.42421                  
#Total         302  130.954 1.00000                  
#
adonis2(bray_bin.dist ~ Site + Lure + mdy + Site/Lure/mdy, strata = metadata.ordered$State, data = metadata.ordered, by = "terms")
#Df SumOfSqs      R2       F Pr(>F)    
#Site           18   29.586 0.22593  5.5922  0.001 ***
#Lure            2    1.218 0.00930  2.0715  0.001 ***
#mdy             1    4.781 0.03651 16.2664  0.001 ***
#Site:Lure      36   14.066 0.10741  1.3293  0.001 ***
#Site:Lure:mdy  56   25.751 0.19664  1.5645  0.001 ***
#Residual      189   55.552 0.42421                   
#Total         302  130.954 1.00000                   
#
#Running this second way pulls a bit of variance out of Site:Lure and a bit more out of Site:Lure:mdy. 
# It seems the first way is more appropriate though
# 
# 

adonis2(bray_bin.dist ~ Site/trapID/mdy + Risk + Lure + mdy, strata = metadata.ordered$State, data = metadata.ordered, by = "terms")

adonis2(bray_bin.dist ~ Risk + Lure + Site/trapID/mdy + mdy, strata = metadata.ordered$State, data = metadata.ordered, by = "terms")
#Df SumOfSqs      R2       F Pr(>F)    
#Risk              1    1.109 0.00847  3.7727  0.001 ***
#Lure              2    1.239 0.00946  2.1078  0.001 ***
#Site             17   28.456 0.21730  5.6949  0.001 ***
#mdy               1    4.781 0.03651 16.2664  0.001 ***
#Site:trapID      36   14.066 0.10741  1.3293  0.001 ***
#Site:trapID:mdy  56   25.751 0.19664  1.5645  0.001 ***
#Residual        189   55.552 0.42421                   
#Total           302  130.954 1.00000                   

#interstingly we get a sig effect of risk level here but the correct way to perform the sig test is using the restricted perms above

# need to do some more digging to fully understand the nested effects structure and if we should use the r2 values from one of the nested adonis models
