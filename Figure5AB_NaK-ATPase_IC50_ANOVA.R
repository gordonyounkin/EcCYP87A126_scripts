####### Script for estimating Na/K ATPase IC50 and ANOVA for CYP87A126 mutant extracts########
####### Gordon Younkin, May 23, 2023 #######
####### Script written for CYP87A126 paper, figure 5A and B #######

library(nlme)
library(multcompView)

##### Load Na/K inhibition data #####
d.nak <- read.csv("5AB_NaK_inhibition")

##### Reformat data to prepare for analysis
d.nak$Plate<- as.factor(d.nak$Plate)
# log transform dilution
d.nak$log.conc<-log(d.nak$Dilution)
# subset data by plant
d.nak.ery<-subset(d.nak, Type=="Plant")
# make unique identifier for each dilution series
d.nak.ery$ID<-factor(paste(d.nak.ery$Species, d.nak.ery$Plant, d.nak.ery$Plate))
# make assay plate and genotype ("Species") a factor
d.nak.ery$Plate<-as.factor(d.nak.ery$Plate)
d.nak.ery$Species<-as.factor(d.nak.ery$Species)
# make groupedData object
model.data.ery<-groupedData(R.BG ~ log.conc | ID, outer=~Plate, inner=~Species, data=d.nak.ery)
# set up parameters
myControl<-nlmeControl(maxIter=1000, msMaxIter=500, msWarnNoConv=F )

#########################################################################
############ model enzyme inhibition for each sample ####################
#########################################################################

###### First round of model fitting #######
### fit nlsList to generate starting value
test.fit2 <- nlsList(R.BG ~ SSfpl(log.conc, A, B, xmid, scal), model.data.ery)
### fit nlme model
ery.m1 <- nlme(test.fit2, control=myControl)

### plot residuals and AIC values to evaluate model fits
plot(resid(ery.m1)~fitted(ery.m1));abline(h=0, lty=2)
text(0.25, 0.02, round(AIC(ery.m1), 0))

###### Second round of model fitting #######
# define pseudo-4parameter-logistic model with optimal fixed scal parameter
logis4<-function(x, A, B, xmid){
  A+(B-A)      / ((1+ exp((xmid-x)/1)))}
# fit nlsList to generate starting value
test.fit2 <- nlsList(R.BG ~ SSfpl(log.conc, A, B, xmid, scal), model.data.ery)
# fit nlme model
ery.m2 <- nlme(test.fit2, control=myControl)

# add plate-specific variances
ery.m3<-update(ery.m2, weights=varIdent(form=~1|Plate))

### compare model fits
anova(ery.m2, ery.m3) # p=0.307 --> no weights parameter required


#########################################################################
############ calculate IC50 for each sample #############################
#########################################################################

### make object to receive concentration values
nak.concentrations<-data.frame(ID=levels(d.nak.ery$ID))

### extract species names from variable ID
nak.concentrations$Species<-  as.factor(unlist(strsplit(as.character(nak.concentrations$ID), " "))[seq(1, 3*length(nak.concentrations$ID), by=3)])

### get sample dilution at 50% inhibition
for(i in 1:nrow(nak.concentrations)){
  nak.concentrations$dilution50[i]<-  exp(fixef(ery.m2)[3] + ranef(ery.m2)[which(rownames(ranef(ery.m2))==nak.concentrations$ID[i]),3])}


### ANOVA ###
# one-way ANOVA and post-hoc test on log10(IC50)
nak.aov <- aov(log10(dilution50) ~ Species, data=nak.concentrations)
# shapiro-wilk test for normality of residuals
shapiro.test(nak.aov$residuals)
# get ANOVA summary with p-value
summary(nak.aov) # anova p-value < 0.001

# significant anova, now post-hoc (Tukey HSD)
tukey <- TukeyHSD(nak.aov)
tukey
cld <- multcompLetters4(nak.aov, tukey)


##############################################################################
######### For plot in figure 5A, make one model per Species/Line #############
##############################################################################

### make groupedData object
model.data.all<-groupedData(R.BG ~ log.conc | Species, data=d.nak)

### fit nlsList to generate starting value
test.fit4 <- nlsList(R.BG ~ SSfpl(log.conc, A, B, xmid, scal), model.data.all)

### fit nlme model
std.m4<-nlme(test.fit4, control=myControl)

### evaluate residuals
plot(std.m4)
