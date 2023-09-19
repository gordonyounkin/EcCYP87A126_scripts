####### Script for one-way ANOVA for T. ni weights on CYP87A126 mutants########
####### Gordon Younkin, May 23, 2023 #######
####### Script written for CYP87A126 paper, figure 5D #######

library(multcompView)

##### load data on T. ni growth #####
tn <- read.csv("5D2_T_ni_growth.csv")

##### Reformat and prepare for statistical analysis #####
# make unique ID for each leaf and for each plant
tn$leafID <- factor(paste(tn$genotype, tn$flat, tn$plant, tn$leaf))
tn$plantID <- factor(paste(tn$genotype, tn$flat, tn$plant))
# for weights after 12 days: remove caterpillars with no data or marked as to be removed
tn.w <- tn[!is.na(tn$weight_5_23), ]
tn.w <- tn.w[!grepl("remove", tn.w$notes_5_23),]
# remove caterpillars on leaf position 6--this was an error in bagging caterpillars
tn.w <- tn.w[tn.w$leaf <= 5, ]

# aggregate (take means)
# reorder genotype levels
tn.w$genotype <- factor(tn.w$genotype, levels=c("WT", "2b16", "x2e1"))
tn.w$leaf <- factor(tn.w$leaf)
tn.w$flat <- factor(tn.w$flat)
tn.w$plantID <- factor(tn.w$plantID)
tn.byleaf <- aggregate(weight_5_23 ~ genotype + flat + plant + leaf, data=tn.w, mean)
tn.byplant <- aggregate(weight_5_23 ~ genotype + flat + plant, data=tn.byleaf, mean)
tn.bygenotype <- aggregate(weight_5_23 ~ genotype, data=tn.byplant, function(x) c(mean=mean(x), sd=sd(x)))
tn.bygenotype <- data.frame(genotype=tn.bygenotype[,1], tn.bygenotype$weight_5_23)

### group by genotype and leaf position
tn.byleaf.geno <- aggregate(weight_5_23 ~ genotype + leaf, data=tn.w, mean)
tn.byleaf$geno_leaf <- paste(tn.byleaf$genotype, tn.byleaf$leaf)
tn.byleaf$flat <- factor(tn.byleaf$flat)

### Normalize by average weight of caterpillar on WT leaf of same age
# get WT means
wtmeans <- tn.byleaf.geno[tn.byleaf.geno$genotype=="WT", "weight_5_23"]
names(wtmeans) <- tn.byleaf.geno[tn.byleaf.geno$genotype=="WT", "leaf"]
# normalize weights to average on that leaf age
tn.byleaf$weight_norm <- tn.byleaf$weight_5_23 / wtmeans[tn.byleaf$leaf]
# aggregate by genotype
tn.bygeno.norm <- aggregate(weight_norm ~ genotype, data=tn.byleaf, function(x) c(mean=mean(x), sd=sd(x)))
tn.bygeno.norm <- data.frame(genotype=tn.bygeno.norm[,1], tn.bygeno.norm$weight_norm)


###################################################
##### ANOVA for normalized caterpillar growth #####
aov.norm <- aov(weight_norm ~ genotype, data=tn.byleaf)
summary(aov.norm) 
# p=0.0187, continue with post-hoc Tukey's HSD
tukey <- TukeyHSD(aov.norm)
tukey
cld <- multcompLetters4(aov.norm, tukey)