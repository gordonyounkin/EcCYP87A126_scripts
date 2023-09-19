####### Script for paired t-tests for Myzus persicae choice assay on CYP87A126 mutants########
####### Gordon Younkin, June 23, 2023 #######
####### Script written for CYP87A126 paper, figure 5C #######

# load M. persicae choice data
mp <- read.csv("5C1_Mp_choice.csv")

##### Rearrange and prepare for analysis #####
# make unique ID for each leaf and for each plant
mp$leafID <- factor(paste(mp$genotype, mp$ID))
# remove samples where one or both leaves were limp/dried out
mp.c <- mp[mp$notes=="", ]
# reorder genotype levels
mp.c$genotype <- factor(mp.c$genotype, levels=c("x2b16", "x2e1"))

#### Paired t-test (for each genotype separately)
t2b16 <- t.test(mp.c[mp.c$genotype=="x2b16", "n_WT"], mp.c[mp.c$genotype=="x2b16", "n_mut"], paired = TRUE, alternative = "two.sided")
t2b16
# 2b16 p-value: 0.12
t2e1 <- t.test(mp.c[mp.c$genotype=="x2e1", "n_WT"], mp.c[mp.c$genotype=="x2e1", "n_mut"], paired = TRUE, alternative = "two.sided")
t2e1
# 2e1 p-value = 0.08