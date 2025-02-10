library(dplyr)
library(vegan)
library(lubridate)
library(stringr)


metadata = read.csv("data/2024_metadata/collated/collated_metadata_NY-OH_11062024.csv")
head(metadata)

#full rarefied table
#
bray_bin.dist = readRDS("data/FEDRR_11062024/processed_tables/avg_dist/bray-binary.full_tab_min450.rds")



########################################
# Perform adonis

#checking that metadata matches dist obj for adonis2
length(rownames(as.matrix(bray_bin.dist))) == nrow(metadata)
rownames(metadata) = metadata$SequenceID

metadata.ordered = metadata[rownames(as.matrix(bray_bin.dist)),]
nrow(metadata.ordered)
sum(rownames(metadata.ordered) == rownames(as.matrix(bray_bin.dist)))

metadata.ordered.noPCR = metadata.ordered
metadata.ordered.noPCR[metadata.ordered.noPCR$BaselineType == "PCRnegative",]

metadata.ordered.noPCR$BaselineType = ifelse(metadata.ordered.noPCR$BaselineType == "PCRnegative", NA, metadata.ordered.noPCR$BaselineType)
metadata.ordered.noPCR
metadata.ordered.noPCR[metadata.ordered.noPCR$trapID == "PCRnegative",]

#no PCR negative
adonis2(bray_bin.dist ~ BaselineType, data = metadata.ordered.noPCR, na.action = na.omit) 
#with PCR negative
adonis2(bray_bin.dist ~ BaselineType, data = metadata.ordered, na.action = na.omit) 

#with PCR negative add state
adonis2(bray_bin.dist ~ BaselineType + State, data = metadata.ordered.noPCR, na.action = na.omit) 
adonis2(bray_bin.dist ~ BaselineType + Site, data = metadata.ordered.noPCR, na.action = na.omit) 

metadata.ordered.noPCR$mdy = mdy(metadata.ordered.noPCR$date)
str(metadata.ordered.noPCR)
adonis2(bray_bin.dist ~ BaselineType + Site + mdy, 
        data = metadata.ordered.noPCR, 
        na.action = na.omit, 
        strata = metadata.ordered.noPCR %>% filter(!is.na(BaselineType)) %>% pull(trapID)
) 

adonis2(bray_bin.dist ~ BaselineType + Site + mdy, 
        data = metadata.ordered.noPCR, 
        na.action = na.omit#, 
        #        strata = metadata.ordered.noPCR %>% filter(!is.na(BaselineType)) %>% pull(trapID)
) 
adonis2(bray_bin.dist ~ BaselineType + trapID + mdy, 
        data = metadata.ordered.noPCR, 
        na.action = na.omit
) 


adonis2(bray_bin.dist ~ BaselineType + Site + mdy, 
        data = metadata.ordered.noPCR, 
        na.action = na.omit,
        by = "margin"
) 
adonis2(bray_bin.dist ~ BaselineType + trapID + mdy, 
        data = metadata.ordered.noPCR, 
        na.action = na.omit,
        by = "margin"
) 
adonis2(bray_bin.dist ~ BaselineType + Site + trapID + mdy, 
        data = metadata.ordered.noPCR, 
        na.action = na.omit,
        by = "margin"
) 

adonis2(bray_bin.dist ~ BaselineType + State + Site + trapID + mdy, 
        data = metadata.ordered.noPCR, 
        na.action = na.omit,
        by = "margin"
) 
#setting by = terms with the correct ordering essentially runs a nested model
# by first assessing significance of 
adonis2(bray_bin.dist ~ BaselineType + State + Site + trapID + mdy, 
        data = metadata.ordered.noPCR, 
        na.action = na.omit,
        by = "terms"
) 
#nested design 
adonis2(bray_bin.dist ~ BaselineType + State/Site/trapID + mdy, 
        data = metadata.ordered.noPCR, 
        na.action = na.omit,
        by = "terms"
) 
#nested design 
adonis2(bray_bin.dist ~ BaselineType + State/Site/trapID + mdy, 
        data = metadata.ordered.noPCR, 
        na.action = na.omit,
        by = "margin"
) 

adonis2(bray_bin.dist ~ BaselineType + State/Site/trapID + mdy, 
        data = metadata.ordered.noPCR, 
        na.action = na.omit
) 

# see https://github.com/vegandevs/vegan/discussions/600 for complex designs
# especially gavin simpson answer about howto permute the second test versus the
# first test when testing an effect need to think about where the "treatment"
# varies, if it does not vary within subects (e.g., pre and post treatment)
# but instead varies between subjects, need to shuffle subjects and NOT shuffle
# samples within subjects (i.e., repeated measures)
# if the treatment instead varies within subjects use subject as block (i.e., permutations
# only allowed within blocks)
#
#its all about permute::how .... (and use dbRDA with Condition() for partialling)
#

#for testing an effect that varies within subject (like baseline or not OR date)
# where subject == trap ID
h1 = with(your_data, how(nperm = 999, blocks = subject))
# the above is equivalent to setting strata in vegan

#for testing an effect that varies between subjects
h2 = with(your_data, 
          how(nperm = 999,
              plots = Plots(strata = subject, type = "free"),
              within = Within(type = "none")
          )
)
# However, note that the number of observations within each plot must be equal!

dbrda(bray_bin.dist ~ within_trt + Condition(subject), permutations = h1, data = your_data)

dbrda(bray_bin.dist ~ between_trt + Condition(subject), permutations = h2, data = your_data)

# Gavin Simpson:
#In both cases I don't think it (the conditioning) is needed for the purposes of the permutation test, but if you ignore it, you might inflate the inertia ("variance", sums of squares, etc) attributable to the treatment effect when you report the inertia explained value.

# the key is to consider the level of (in)dependence. E.g., because multiple
# measurements occur on units we must restrict permutation to those units
# (either within or between)
# This also guides choice of permutation unit with multiple levels of dependence, 
# i.e., nested or hierarchical designs in this case the level to be restricted
# should be determined by the next level of dependence. For example in our design
# with traps within sites within states, to test the effect of state we can
# permute sites (i.e., the whole block of samples constituting a site) but 
# not permute samples *within* sites. To test the effect of a trap level 
# treatment where multiple samples are performed on each trap, we can permute 
# the samples from each trap as a block (as is done in h2). To test the effect
# of a site level treatment we can permute whole blocks of samples at the site 
# level. Here we can (should) condition on state to account for the variation due
# to state. We could consider restricting permutations to within state using
# the block level of how (which sits above plot
# 

h3 = with(your_data, 
          how(nperm = 999,
              blocks = State,
              plots = Plots(strata = Site, type = "free"),
              within = Within(type = "none")
          )
)

dbrda(bray_bin.dist ~ between_site_trt + Condition(State), permutations = h3, data = your_data)


###############################################################################
###############################################################################
###############################################################################
# Permutations examples !!!
# 
# 

# Setting up a nested design. This could be thought of as a split plot or a 
# geogrpahically/spatially nested design

#2 blocks (e.g., two fields) A and B
# each block with 2 plots C to F
# each plot has 3 subjects (1 to 18), these could be thought of as individual
# plants, animals, traps (e.g., an insect trap), a secondary spatial subdivision
# (e.g., subplot or transect)
# Each subplot with four repeated measures (e.g., measurements repeated through 
# a growing season OR individuals with pre treatment and post treatment measures
# 
# n = 48
# 
# Note this is similar to our USDA-FS BIL design where each state (block) has
# six sites (plots) each with three traps (subjects) for which multiple measures 
# are made through the growing season (repeated measures)
block_v = rep(LETTERS[1:2], each = 24)
#plot_v = rep(LETTERS[3:6], each = 12)
plot_v = rep(LETTERS[1:4], each = 12)
#subject_v = rep(letters[7:18], each = 4)
subject_v = rep(rep(letters[1:3], each = 4), 4)
#rep_mesr_v = paste0("_", 1:48)
rep_mesr_v = paste0("_", 1:4)

block_design = data.frame(
    block = block_v,
    plot = paste0(block_v, plot_v), 
    #subject = stringr::str_sort(rep(paste0("s", 1:12), 48/12), numeric = T), #stringr for natural sort
    subject = paste0(block_v, plot_v, subject_v),
    rep = paste0(block_v, plot_v, subject_v, rep_mesr_v)
)

block_design

plot_lvl_trt = rep( rep(letters[1:2], each = 12), 2)
subj_lvl_trt = rep( rep(letters[1:3], each = 4), 4)
repeat_lvl_trt = rep( rep(letters[1:2], each = 2), 12) 

# treatment applied to repeated measures of subjects
# OR could be something like examining the effect of time (e.g., the measure effect itself)
h1 = with(block_design, how(nperm = 99, blocks = subject))
shuffleSet(n = 48, nset = 99, control = h1)
h1.p = shuffle(48, control = h1)

block_design.Tr = block_design

block_design.Tr$trt = repeat_lvl_trt
block_design.Tr
block_design.Tr[h1.p,]
block_design.Tr[h1.p, ] == block_design.Tr

#for testing an effect that varies between subjects
# e.g., this could be trap lure
# here we will want to use plot as blocvk
h2 = with(block_design, 
          how(nperm = 99,
              blocks = plot,
              plots = Plots(strata = subject, type = "free"),
              within = Within(type = "none")
          )
)
h2
s = shuffleSet(n = 48, nset = 99, control = h2)
str(s)

h2.p = shuffle(48, control = h2)

block_design.Tr = block_design
block_design.Tr
block_design.Tr$trt = subj_lvl_trt
block_design.Tr[h2.p,]
block_design.Tr[h2.p, ] == block_design.Tr


#for testing an effect that varies between plots
# e.g., this could be risk type
# here we will use block as block
h3 = with(block_design, 
          how(nperm = 999,
              complete = T,
              blocks = block,
              plots = Plots(strata = plot, type = "free"),
              within = Within(type = "none")
          )
)
h3
h3.p = shuffle(48, control = h3)

block_design.Tr = block_design
block_design.Tr$trt = plot_lvl_trt
block_design.Tr
block_design.Tr[h3.p,]
block_design.Tr[h3.p, ] == block_design.Tr


#perform a test
sp_dat = data.frame(replicate(50,sample(0:10,48,rep=TRUE)))
adonis2(sp_dat ~ trt, data = block_design.Tr, permutations = h3)
#NOTE! If the number of *possible* permutations is less than the requested
# number of permutations there is no message thrown by how(), instead the 
# message is thrown at the stage of testing i.e., by adonis2 or anova.cca
# 
foo = dbrda(sp_dat ~ trt, data = block_design.Tr)
anova(foo, permutations = h3)
#
#
#OR the number of perms can be checked by running through shuffleSet
shuffleSet(n = 48, nset = 99, control = h3)

#NOTE the number of possible permutations is given by factorial(s)^p where 
# s = the number of subjects to permute and p = the number of blocks of s
# i.e., restricted units in which s can be permuted within but not between
# for example in h2 above
h2 = with(block_design, 
          how(nperm = 99,
              blocks = plot,
              plots = Plots(strata = subject, type = "free"),
              within = Within(type = "none")
          )
)
# plots are 4 and subjects per plot are 3 so
factorial(3)^4
# = 1296
# note that             
s = shuffleSet(n = 48, nset = 99, control = h2)
s
str(s)
# $maxperm gives the maximum possible permutations (when less than minperm)


###################################################
###################################################
# Not sure what the implementation of this interms of design,
# but it is also possible to permute both plots and within plots
# Let's see what it looks like
# We run with series for time series/transect type permutes

h4 = with(block_design,
    how(nperm = 99,
        blocks = plot,
        plots = Plots(strata = subject, type = "free"),
        within = Within(type = "series")
    )
)
h4
shuffleSet(n = 48, nset = 99, control = h4)
              
h4.p = shuffle(48, control = h4)

block_design.Tr = block_design
block_design.Tr$trt = subj_lvl_trt
block_design.Tr
block_design.Tr[h4.p,]
block_design.Tr[h4.p, ] == block_design.Tr
# as expected this shuffles reps within subjects and moves groups of subjects 
# within plots (i.e., the subjects are not moved between plots, and subject reps stay together)
# this could be used to test treatments that are appplied at the subject and within subject level while
# respecting the dependence of groupsof subjects within plots and plots within 
# blocks. A similar analysis could be apllied one leve up (i.e, block/plot level
# 
h5 = with(block_design,
          how(nperm = 99,
              blocks = block,
              plots = Plots(strata = plot, type = "free"),
              within = Within(type = "free")
          )
)
h5
shuffleSet(n = 48, nset = 99, control = h5)

h5.p = shuffle(48, control = h5)

block_design.Tr = block_design
block_design.Tr$trt = plot_lvl_trt
block_design.Tr
block_design.Tr[h5.p,]
block_design.Tr[h5.p, ] == block_design.Tr

#note that with this design it is the *REPS*that are permuted at the lowest
# level, This is because the algorithm permutes each row of subject individually
# not recognizing that there are groups of subjects. This means that within a 
# plot ALL of the reps can be permuted freely.

h6 = with(block_design,
          how(nperm = 99,
              blocks = block,
              plots = Plots(strata = plot, type = "free"),
              within = Within(type = "series")
          )
)
h6
shuffleSet(n = 48, nset = 99, control = h6)

h6.p = shuffle(48, control = h6)

block_design.Tr = block_design
block_design.Tr$trt = plot_lvl_trt
block_design.Tr
block_design.Tr[h6.p,]
block_design.Tr[h6.p, ] == block_design.Tr
#applying series retains the dependence between reps in that they are 
# ordered according to the original ordering. This is the way to go probably
# even though the dependence between subject is not respected
# 
#It may be possible to perform a two level permutation by first permuting within
# lower levels and then permuting higher levels (on the permuted df) while
# indicating wihtin = "none"







