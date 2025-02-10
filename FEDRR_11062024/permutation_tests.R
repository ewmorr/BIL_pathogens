library(dplyr)
library(vegan)
library(lubridate)
set.seed(12345)

metadata = read.csv("data/2024_metadata/collated/collated_metadata_NY-OH_11062024.csv")
head(metadata)

########################################
########################################
########################################
#samps and baseline dist
#
bray_bin.dist = readRDS("data/FEDRR_11062024/processed_tables/avg_dist/bray-binary.samps_and_baseline.rds")

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
#BaselineType   1    2.627 7.0445  0.001 ***
#Residual     197   73.469                  
#
anova(baseline.dbrda, permutations = h1, by ="terms", parallel = 4, model = "direct") #margins doesn't accomplish anmything with only one fixed effect
#
baseline.dbrda = dbrda(bray_bin.dist ~ BaselineType, permutations = h1, data = metadata.ordered)
anova(baseline.dbrda, permutations = h1, by ="terms")
#              Df SumOfSqs      F Pr(>F)    
#BaselineType   1    2.666 6.2471  0.001 ***
#    Residual     232   98.999                  
adonis2(bray_bin.dist ~ BaselineType, permutations = h1, data = metadata.ordered)
#          Df SumOfSqs      R2      F Pr(>F)    
#Model      1    2.666 0.02622 6.2471  0.001 ***
#Residual 232   98.999 0.97378                  


########################################
########################################
########################################
#main tests (date, lure, risk)
#
bray_bin.dist = readRDS("data/FEDRR_11062024/processed_tables/avg_dist/bray-binary.samps.rds")

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
#Model      1    3.638 11.065  0.001 ***
#    Residual 168   55.230                  
#---
3.638/sum(3.638,55.230)
# 0.06179928 
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
# should throw out the n of 3 and downsample the rest to five (at random)
metadata.ordered %>%
    group_by(State, Site, Lure, trapID) %>%
    summarize(n = n()) %>%
    as.data.frame() -> trap_n

lure_test.samples_exclude = metadata.ordered %>% filter(trapID == 1825) %>% pull(SequenceID)
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
anova(h2.dbrda, permutations = h2, by ="terms")
#          Df SumOfSqs      F Pr(>F)   
#Lure       2    0.869 1.2479  0.002 **
#Residual 161   56.088                 
0.869/sum(0.869,56.088)
#0.01525712
# this agrees with some of the adonis models below testing for Lure non nested
h2.dbrda = dbrda(bray_bin.dist.lure_test ~ trapID + Condition(Site), data = metadata.lure_test)

#for Risk we are permuting sites within states (block = state, strata = plots, within = none
metadata.ordered %>%
    group_by(State, Site) %>%
    summarize(n = n()) %>%
    print(n = Inf)
# we don't want to throw out a whole site so we will downsample everything to the lowest n (15)
metadata.ordered %>%
    group_by(State, Site) %>%
    summarize(n = n()) -> site_n

#just testing the indexing
for(i in 1:10){
    print(paste(i*15-15+1, i*15))
}
#
risk_test.samples_include = vector(mode = "character", length = 15*nrow(site_n))
for(i in 1:nrow(site_n)){
    tmpIDS = metadata.ordered %>% filter(Site == site_n$Site[i]) %>% pull(SequenceID)
    risk_test.samples_include[(i*15-15+1):(i*15)] = sample(tmpIDS, 15, replace = F)
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
#Model      1    1.063 2.7003  0.487
#Residual 177   69.700              
1.063/sum(1.063,69.700)
#0.015
# this r^2 agrees with some of the adonis models below testing for risk non nested


metadata.risk_test %>% group_by(State, Site, Risk) %>% summarize(n = n())

# to get a sense of variance explained by location we test trap within site
# and site within state
h4.dbrda = dbrda(bray_bin.dist.lure_test ~ trapID, data = metadata.lure_test)
anova(h4.dbrda, permutations = h2)
#Model     34   24.651 2.0931      1
#Residual 140   48.494              
24.651/sum(24.651,48.494)
# 0.3370155
h5.dbrda = dbrda(bray_bin.dist.risk_test ~ Site, data = metadata.risk_test)
anova(h5.dbrda, permutations = h3)
#          Df SumOfSqs      F Pr(>F)
#Model     11   15.792 4.0767      1
#Residual 168   59.161              
15.792/sum(15.792,59.161)
#0.210692


adonis2(bray_bin.dist ~ State/Site/trapID, data = metadata.ordered, by = "terms", permutations = 1999)
Df SumOfSqs      R2       F Pr(>F)    
#State               1    4.711 0.05509 13.5247 0.0005 ***
#    State:Site         10   13.081 0.15296  3.7553 0.0005 ***
#    State:Site:trapID  24    8.857 0.10357  1.0594 0.0790 .  


###############################################################################

adonis2(bray_bin.dist ~ Lure + mdy, strata= metadata.ordered$Site, data = metadata.ordered)
#          Df SumOfSqs     R2      F Pr(>F)    
#Model      3    6.082 0.0705 5.1327  0.001 ***
#Residual 203   80.183 0.9295                  
#Total    206   86.265 1.0000                  
adonis2(bray_bin.dist ~ Lure + mdy, strata= metadata.ordered$Site, data = metadata.ordered, by = "terms")
#Df SumOfSqs      R2      F Pr(>F)    
#Lure       2    0.867 0.01005  1.097  0.021 *  
#    mdy        1    5.215 0.06046 13.204  0.001 ***
#    Residual 203   80.183 0.92950                  
#Total    206   86.265 1.00000                  
#--
adonis2(bray_bin.dist ~ Lure + mdy, strata= metadata.ordered$Site, data = metadata.ordered, by = "margin")
#          Df SumOfSqs      R2       F Pr(>F)    
#    Lure       2    0.869 0.01008  1.1005  0.034 *  
#    mdy        1    5.215 0.06046 13.2040  0.001 ***
#    Residual 203   80.183 0.92950                   
#    Total    206   86.265 1.00000  

adonis2(bray_bin.dist ~ Site + Lure + mdy, strata = metadata.ordered$State, data = metadata.ordered)
#          Df SumOfSqs      R2      F Pr(>F)    
#Model     14   22.383 0.25947 4.8052  0.001 ***
#Residual 192   63.882 0.74053                  
#Total    206   86.265 1.00000                  
adonis2(bray_bin.dist ~ Site/Lure/mdy, strata = metadata.ordered$State, data = metadata.ordered)
#          Df SumOfSqs      R2      F Pr(>F)    
#Model     71   46.215 0.53573 2.1941  0.001 ***
#Residual 135   40.050 0.46427                  
#Total    206   86.265 1.00000                  
adonis2(bray_bin.dist ~ Site + Lure + mdy + Site/Lure/mdy, strata = metadata.ordered$State, data = metadata.ordered)
#          Df SumOfSqs      R2      F Pr(>F)    
#Model     71   46.215 0.53573 2.1941  0.001 ***
#Residual 135   40.050 0.46427                  
#Total    206   86.265 1.00000                  
adonis2(bray_bin.dist ~ Site/Lure/mdy, strata = metadata.ordered$State, data = metadata.ordered, by = "terms")
#               Df SumOfSqs      R2      F Pr(>F)    
#Site           11   17.813 0.20649 5.4585  0.001 ***
#Site:Lure      24    8.913 0.10332 1.2518  0.001 ***
#Site:Lure:mdy  36   19.489 0.22592 1.8248  0.001 ***
#Residual      135   40.050 0.46427                  
#Total         206   86.265 1.00000                  
#
adonis2(bray_bin.dist ~ Site + Lure + mdy + Site/Lure/mdy, strata = metadata.ordered$State, data = metadata.ordered, by = "terms")
#Df SumOfSqs      R2       F Pr(>F)    
#Site           11   17.813 0.20649  5.4585  0.001 ***
#    Lure            2    0.852 0.00988  1.4364  0.005 ** 
#    mdy             1    3.718 0.04310 12.5316  0.001 ***
#    Site:Lure      22    7.989 0.09261  1.2240  0.001 ***
#    Site:Lure:mdy  35   15.844 0.18366  1.5259  0.001 ***
#    Residual      135   40.050 0.46427                   
#Total         206   86.265 1.00000                   
#
#Running this second way pulls a bit of variance out of Site:Lure and a bit more out of Site:Lure:mdy. 
# It seems the first way is more appropriate though
# 
# 

adonis2(bray_bin.dist ~ Site/trapID/mdy + Risk + Lure + mdy, strata = metadata.ordered$State, data = metadata.ordered, by = "terms")

adonis2(bray_bin.dist ~ Risk + Lure + Site/trapID/mdy + mdy, strata = metadata.ordered$State, data = metadata.ordered, by = "terms")
#Df SumOfSqs      R2       F Pr(>F)    
#Risk              1    1.152 0.01336  3.8841  0.001 ***
#Lure              2    0.864 0.01002  1.4569  0.007 ** 
#Site             10   16.648 0.19299  5.6118  0.001 ***
#mdy               1    3.718 0.04310 12.5316  0.001 ***
#Site:trapID      22    7.989 0.09261  1.2240  0.001 ***
#Site:trapID:mdy  35   15.844 0.18366  1.5259  0.001 ***
#Residual        135   40.050 0.46427                   
#Total           206   86.265 1.00000                   


