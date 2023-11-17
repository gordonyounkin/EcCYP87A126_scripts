####### Script for chi-squared for P. rapae oviposition and feeding on CYP87A126 mutants########
####### Gordon Younkin, June 28, 2023 #######
####### Script written for CYP87A126 paper, figure 5F #######

################################################################################
#################### OVIPOSITION ###############################################
##### load P. rapae oviposition data
pr <- read.csv("5F1_Pr_oviposition.csv")

##### Reformat and prepare for analysis
pr <- pr[pr$notes != "no eggs", ]

##### Chi-squared test for aggregate number of eggs laid on WT vs mutant, separate for each genotype #####
cyp87a126.2 <- chisq.test(colSums(pr[pr$genotype=="cyp87a126-2", c("n_WT", "n_mut")]), correct=FALSE)
cyp87a126.2 # cyp87a126-2: p<0.001
cyp87a126.1 <- chisq.test(colSums(pr[pr$genotype=="cyp87a126-1", c("n_WT", "n_mut")]), correct=FALSE)
cyp87a126.1 # cyp87a126-1: p<0.001

################################################################################
################################################################################

################################################################################
#################### CATERPILLAR FEEDING #######################################
#### load P. rapae feeding data
pr_damage <- read.csv("5F2_Pr_feeding.csv")

##### tabulate number of damaged leaves per genotype
table(pr_damage[pr_damage$leaf_damaged=="N", "genotype"])
undamaged <- c(14,1,0) # number of undamaged leaves for WT, cyp87a126-1, cyp87a126-2
table(pr_damage[pr_damage$leaf_damaged=="Y", "genotype"], useNA="ifany")
damaged <- c(0,16,13) # number of damaged leaves for WT, cyp87a126-1, cyp87a126-2
pr_dam2 <- data.frame("damaged"=damaged, "undamaged"=undamaged)
row.names(pr_dam2) <- c("WT","cyp87a126-1","cyp87a126-2")

### Chi squared test for caterpillar commencing feeding and causing damage to leaf, by genotype
chisq.test(pr_dam2, correct=FALSE) # p<0.001