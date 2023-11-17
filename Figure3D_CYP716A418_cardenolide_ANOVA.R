####### Script for ANOVA of CYP716A418 vs WT cardiac glycoside abundances ########
####### Gordon Younkin, May 19, 2023 #######
####### Script written for CYP87A126 paper, figure 3D #######

library(multcompView)

# load data
cg <- read.csv("3D_CYP716A418_cardenolides.csv")

##### Rearrange and aggregate data to prepare for statistical analysis ######
# reorder genotype levels
cg$genotype <- factor(cg$genotype, levels=c("WT", "cyp716a418-1", "cyp716a418-2"))
# normalize such that WT average = 1 for each genin
cg$digitoxigenin.norm <- cg$digitoxigenin / mean(cg[cg$genotype=="WT", "digitoxigenin"])
cg$cannogenol.norm <- cg$cannogenol / mean(cg[cg$genotype=="WT", "cannogenol"])
cg$cannogenin.norm <- cg$cannogenin / mean(cg[cg$genotype=="WT", "cannogenin"])
cg$strophanthidin.norm <- cg$strophanthidin / mean(cg[cg$genotype=="WT", "strophanthidin"])
# aggregate by genotype
cg.bygenotype <- aggregate(cbind(digitoxigenin.norm, cannogenol.norm, cannogenin.norm, strophanthidin.norm) ~ genotype, data=cg, function(x) c(mean=mean(x), sd=sd(x)))
cg.bygenotype <- data.frame(genotype=cg.bygenotype[,1], cg.bygenotype$digitoxigenin.norm, cg.bygenotype$cannogenol.norm, cg.bygenotype$cannogenin.norm, cg.bygenotype$strophanthidin.norm)
names(cg.bygenotype) <- c("genotype", "digitoxigenin.mean", "digitoxigenin.sd", "cannogenol.mean", "cannogenol.sd", "cannogenin.mean", "cannogenin.sd",
                          "strophanthidin.mean", "strophanthidin.sd")

##### ANOVA #####
# ANOVA and post-hoc test on normalized ion count, separate for each genin
cg$genotype <- factor(cg$genotype, levels=c("WT", "cyp716a418-1", "cyp716a418-2"))
# DIGITOXIGENIN
dig.aov <- aov(digitoxigenin ~ genotype, data=cg)
summary(dig.aov) # p < 0.001
dig.tukey <- TukeyHSD(dig.aov)
dig.cld <- multcompLetters4(dig.aov, dig.tukey)
dig.letters <- dig.cld$genotype$Letters
# CANNOGENOL
canol.aov <- aov(cannogenol ~ genotype, data=cg)
summary(canol.aov) # p < 0.001
canol.tukey <- TukeyHSD(canol.aov)
canol.cld <- multcompLetters4(canol.aov, canol.tukey)
canol.letters <- canol.cld$genotype$Letters
# CANNOGENIN
canin.aov <- aov(cannogenin ~ genotype, data=cg)
summary(canin.aov)
canin.tukey <- TukeyHSD(canin.aov) # p < 0.001
canin.cld <- multcompLetters4(canin.aov, canin.tukey)
canin.letters <- canin.cld$genotype$Letters
summary(canin.aov)
# STROPHANTHIDIN
stro.aov <- aov(strophanthidin ~ genotype, data=cg)
summary(stro.aov) # p < 0.001
stro.tukey <- TukeyHSD(stro.aov)
stro.cld <- multcompLetters4(stro.aov, stro.tukey)
stro.letters <- stro.cld$genotype$Letters

