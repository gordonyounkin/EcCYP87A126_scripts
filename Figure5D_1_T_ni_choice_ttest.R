####### Script to do paired t-tests for T. ni choice assay on CYP87A126 mutants########
####### Gordon Younkin, June 7, 2023 #######
####### Script written for CYP87A126 paper, figure 5D #######

##### load T. ni choice data #####
tn <- read.csv("5D1_T_ni_choice")

##### Reformat and prepare for analysis #####
# make unique ID for each leaf and for each plant
tn$leafID <- factor(paste(tn$genotype, tn$ID))
# remove caterpillars with no data or marked as to be removed
tn.c <- tn[!is.na(tn$mut_damage_area), ]
# reorder genotype levels
tn.c$genotype <- factor(tn.c$genotype, levels=c("cyp87a126-1", "cyp87a126-2"))

##### Paired t-test (for each genotype separately) #####
t2b16 <- t.test(tn.c[tn.c$genotype=="cyp87a126-1", "WT_damage_area"], tn.c[tn.c$genotype=="cyp87a126-1", "mut_damage_area"], paired = TRUE, alternative = "two.sided")
t2b16
# 2b16 (cyp87a126-1) p-value: 0.00013
t2e1 <- t.test(tn.c[tn.c$genotype=="cyp87a126-2", "WT_damage_area"], tn.c[tn.c$genotype=="cyp87a126-2", "mut_damage_area"], paired = TRUE, alternative = "two.sided")
t2e1
# 2e1 (cyp87a126-2) p-value = 0.01299