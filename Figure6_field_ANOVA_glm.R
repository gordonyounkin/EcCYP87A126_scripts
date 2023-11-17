####### Script for conducting statistical analyses on field data associated with cyp87a126 knockout lines  ########
####### Gordon Younkin, October 26, 2023 #######
####### Script written for CYP87A126 paper, figure 6 #######

library(boot)
library(multcomp)
library(multcompView)
library(emmeans)
library(effects)

##### Load Data #####
# load field observations
obs <- read.csv("20231003_Erysimum_field_observations_wide.csv")
obs$block <- as.factor(obs$block)

# visitors observed at >25 plants over 4 censuses:
#[1] "flea_beetle_yellow" "other_flea_beetle"  "turnip_aphid"       "pieris_rapae_egg"   "leafhopper"         "snail"              "spider"            
#[8] "leaf_miner"


################################################################################
###############   linear model for Phyllotreta striolata   #####################
################################################################################

##### get P. striolata (flea_beetle_yellow) #####
herbivore <- "flea_beetle_yellow"
# get data for that herbivore
herb.df <- cbind(obs[,c("ID", "block", "plant", "genotype", "flowering_visit_2", "height_3", "height_4")], obs[,grep(paste0(herbivore, "_"), names(obs))])
# get presence/absence
herb.df$herb.sum <-  rowSums(herb.df[,grep(paste0(herbivore, "_"), names(herb.df))])
# remove NAs
herb.df <- herb.df[!is.na(herb.df$herb.sum), ]
# remove plants shorter than 30cm on visit 3--these had been topped but had axillary shoots growing and were much smaller
# and generally had fewer visitors
herb.df <- herb.df[herb.df$height_3 > 30,]
herb.df$genotype <- factor(herb.df$genotype, levels=c("WT","cyp87a126-1","cyp87a126-2"))


### one-way ANOVA count of P. striolata by genotype ###
aov <- aov(flea_beetle_yellow_4 ~ genotype, data=herb.df)
summary(aov)
# p=0.113, no post-hoc test


################################################################################
##############   logistic model for all other visitors   #######################
################################################################################
# make vector with species found on >25 plants during the four censuses
glm_specs <- c("other_flea_beetle", "turnip_aphid", "pieris_rapae_egg", "leafhopper", "snail", "spider", "leaf_miner")
# make empty list to store models and estimates
models <- as.list(glm_specs)
names(models) <- glm_specs
estimates <- as.list(glm_specs)
names(estimates) <- glm_specs
cis <- as.list(glm_specs)
names(cis) <- glm_specs

# loop through species for which to generate models
for(herbivore in glm_specs) {
  # get data for that herbivore
  herb.df <- cbind(obs[,c("ID", "block", "plant", "genotype", "flowering_visit_2", "height_3", "height_4")], obs[,grep(paste0(herbivore, "_"), names(obs))])
  # get presence/absence
  herb.df$present <-  rowSums(herb.df[,grep(paste0(herbivore, "_"), names(herb.df))]) >= 1
  # remove NAs
  herb.df <- herb.df[!is.na(herb.df$present), ]
  herb.df$genotype <- factor(herb.df$genotype, levels=c("WT","cyp87a126-1","cyp87a126-2"))
  
  ##### fit logistic model (glm + family=binomial) #####
  herb.glm <- glm(present ~ genotype, data=herb.df, family="binomial")
  models[[herbivore]] <- summary(herb.glm)
  
  ##### comparisons between genotype estimated marginal means #####
  marginal <- emmeans(herb.glm, ~ genotype)
  estimates[[herbivore]] <- pairs(marginal)
  
  ##### confidence interval for estimated marginal means, back-transformed #####
  cis[[herbivore]] <- as.data.frame(Effect("genotype", herb.glm))
  
  }
