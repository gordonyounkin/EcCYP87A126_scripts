####### Script to do Kruskal-Wallis test for B. brassicae colony growth on CYP87A126 mutants########
####### Gordon Younkin, October 23, 2023 #######
####### Script written for CYP87A126 paper, figure 5E #######

##### Load data #####
# load B. brassicae colony counts
bb <- read.csv("20231020_B_brassicae_colony_growth.csv")

##### Rearrange and prepare for analysis #####
bb$genotype <- factor(bb$genotype, levels=c("WT", "cyp87a126-2", "cyp87a126-2"))


### ANOVA (Use Kruskal-Wallis because WT is almost entirely 0==not normally distributed)
kwt <- kruskal.test(total_aphids ~ genotype, data=bb)
kwt
# p < 0.001, do post-hoc test Pairwise wilcox test
bb.wilcox <- pairwise.wilcox.test(bb$total_aphids, bb$genotype, p.adjust.method="bonferroni")