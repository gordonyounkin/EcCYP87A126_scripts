####### Script to do ANOVA for M. persicae colony growth on CYP87A126 mutants########
####### Gordon Younkin, May 23, 2023 #######
####### Script written for CYP87A126 paper, figure 5C #######

library(multcompView)

##### Load data #####
# M. persicae colony counts
mp <- read.csv("5C2_Mp_growth.csv")

##### Rearrange and prepare for analysis #####
# reformat where needed
mp$genotype <- factor(mp$genotype, levels=c("WT", "cyp87a126-1", "cyp87a126-2"))
mp.bygenotype <- aggregate(total_aphids ~ genotype, data=mp, function(x) c(mean=mean(x), sd=sd(x)))
mp.bygenotype <- data.frame(genotype=mp.bygenotype[,1], mp.bygenotype$total_aphids)


### ANOVA
aov <- aov(total_aphids ~ genotype, data=mp)
summary(aov) 
# p = 0.767, no differences between groups