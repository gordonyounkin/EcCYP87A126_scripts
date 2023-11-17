####### Script for paired t-tests for Brevicoryne brassicae choice assay on CYP87A126 mutants########
####### Gordon Younkin, September 27, 2023 #######
####### Script written for CYP87A3 paper, figure 5E #######

##### load B. brassicae choice data
bb <- read.csv("./20230927_cabbage_aphid_choice.csv")

##### Rearrange and prepare for analysis #####
# make unique ID for each leaf and for each plant
bb$leafID <- factor(paste(bb$genotype, bb$ID))
# remove samples where one or both leaves were limp/dried out
bb.c <- bb[bb$notes=="", ]
# reorder genotype levels
bb.c$genotype <- factor(bb.c$genotype, levels=c("cyp87a126-1", "cyp87a126-2"))

#### Paired t-test (for each genotype separately)
t2b16 <- t.test(bb.c[bb.c$genotype=="cyp87a126-1", "n_WT"], bb.c[bb.c$genotype=="cyp87a126-1", "n_mut"], paired = TRUE, alternative = "two.sided")
t2b16
# 2b16 p-value: 0.0.004
t2e1 <- t.test(bb.c[bb.c$genotype=="cyp87a126-2", "n_WT"], bb.c[bb.c$genotype=="cyp87a126-2", "n_mut"], paired = TRUE, alternative = "two.sided")
t2e1
# 2e1 p-value = 0.009