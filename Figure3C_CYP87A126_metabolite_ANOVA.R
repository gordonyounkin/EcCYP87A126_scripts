####### Script to do ANOVA for CYP87A3 vs WT cardiac glycoside and glucosinolate abundances########
####### Gordon Younkin, May 19, 2023 #######
####### Script written for CYP87A126 paper, figure 3C #######

library(multcompView)

###############################################################################
##### CARDENOLIDES ############################################################
##### Load data ######
cg <- read.csv("3C_1_CYP87A126_cardenolides.csv")

##### Rearrange and aggregate data to prepare for statistical analysis ######
# reorder genotype levels
cg$genotype <- factor(cg$genotype, levels=c("WT", "4g2680.2b16", "4g2680.2e1"))
# Normalize to average abundance in WT
cg$total_norm <- cg$total / mean(cg[cg$genotype=="WT", "total"], na.rm=TRUE)
# aggregate by genotype
cg.bygenotype <- aggregate(total_norm ~ genotype, data=cg, function(x) c(mean=mean(x), sd=sd(x)))
cg.bygenotype <- data.frame(genotype=cg.bygenotype[,1], cg.bygenotype$total_norm)

### one-way ANOVA ###
# ANOVA and post-hoc test on log10(ion count)
cg$genotype <- factor(cg$genotype, levels=c("WT", "4g2680.2b16", "4g2680.2e1"))
cg.aov <- aov(total ~ genotype, data=cg)
# get ANOVA summary with p-value
summary(cg.aov) # anova p-value <<< 0.001

# significant ANOVA, now post-hoc (Tukey HSD)
tukey <- TukeyHSD(cg.aov)
tukey
cld <- multcompLetters4(cg.aov, tukey)

###############################################################################
###############################################################################

###############################################################################
##### Glucosinolates ##########################################################
##### Load data #####
gsl <- read.csv("3C_2_CYP87A126_glucosinolates.csv")

##### Rearrange and aggregate data to prepare for statistical analysis ######
# reorder genotype levels
gsl$genotype <- factor(gsl$genotype, levels=c("WT", "4g2680.2b16", "4g2680.2e1"))
# normalize such that WT average = 1, for aliphatic and indole glucosinolates separately
gsl$aliphatic.norm <- gsl$aliphatic / mean(gsl[gsl$genotype=="WT", "aliphatic"])
gsl$indole.norm <- gsl$indole / mean(gsl[gsl$genotype=="WT", "indole"])
# aggregate by genotype
gsl.bygenotype <- aggregate(cbind(aliphatic.norm, indole.norm, total) ~ genotype, data=gsl, function(x) c(mean=mean(x), sd=sd(x)))
gsl.bygenotype <- data.frame(genotype=gsl.bygenotype[,1], gsl.bygenotype$aliphatic.norm, gsl.bygenotype$indole.norm, gsl.bygenotype$total)
names(gsl.bygenotype) <- c("genotype", "aliphatic.mean", "aliphatic.sd", "indole.mean", "indole.sd", "total.mean", "total.sd")

### ANOVA ###
# ANOVA and post-hoc test on log10(ion count)
gsl$genotype <- factor(gsl$genotype, levels=c("WT", "4g2680.2b16", "4g2680.2e1"))
# ALIPHATIC
ali.aov <- aov(aliphatic ~ genotype, data=gsl)
summary(ali.aov) # p = 0.149
# aliphatic glucosinolate ANOVA not significant, no post-hoc test

# INDOLE
ind.aov <- aov(indole ~ genotype, data=gsl)
summary(ind.aov) # p = 0.23
# indole glucosinolate ANOVA not significant, no post-hoc test


