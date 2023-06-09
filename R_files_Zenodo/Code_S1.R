################################################################################
## Non random plant selection for medicine using Benin medicinal Flora        ##
## Orou Gaoue & Yessoufou Kowiyou,  March 8, 2021 | ogaoue@utk.edu            ##
## University of Tennessee Koxville, University of Johannesburg               ##
################################################################################

################################################################################
## PART A: Moerman analysis corrected                                         ##
################################################################################
## I. Data: Complied from the following documents: Benin medicinal plants compendium (Adjanonhou et al 1986) and Benin Flora (Akeogninou et al 2005)

## II. Research questions:
## 1. Are local people using medicinal plants in a random way? 
## 1.1. Are most families under- or over-represented in the local pharmacopoeia? Are the same kind of families over- under-represented as in other studies? 
## 1.2. From a methodological perspective, would using poisson or negbin error structure (given that the response variable here is a count data), rather than the classical general linear model (withoiut transformation: Moerman 1996 and others, or with log+1 transformation: Ford and Gaoue 2017), yield different results, thereby demonstrating a classifaction bias associated with incorrect use of the model?

## 2. What are the drivers of over-representation? Biochemical?
## 2.1. Is there any phylogenetic signal for over-utilized families? in other words, are over-utilized families clustered, thereby signaling an evolutionary underlying mechanisms for their overuse in indigenous medicine?
## 2.2. Is there any phylogenetic signal for organs used and is that consistent with prediction from the life history theory and optimal defense theory? Are species groupped by life form and organs used? Are some families more used for a given organ than other? Does it demonstrate any link between the organs that are most likely to be vital for species fitness (life history theory) and therefore more protected by secondary compounds (optimal defense theory)?

## 3. Can we use phylogeny to identify the most likely source of medicine?
## Are there nodes or clades that are more medicinal than other? What are teh characteristics of these nodes and how does it compare with global prediction from State of Wolrd Plants 2017?

rm(list=ls(all=TRUE))  ## Clear all

install_packages("devtools")
# install_packages("R2admb")
# devtools::install_github("bbolker/glmmadmb")

###################################################################################
## A. Dataset
###################################################################################
setwd("/Users/ogaoue/Dropbox/PROJECTS/Benin Med Flora/Data/CSV")
benin_med_flora<-read.csv("Benin_Med_Plants_APGIV_July20_2017Final.csv", header=T)
## Use Data: Dataset_S1.csv *******

## Sort data in ascending order to facilitate ploting credible intervals later in the code.
benin_med_flora<-benin_med_flora[order(benin_med_flora$tot_species),]

## Checking the nature of each variable
str(benin_med_flora) ## Data has 180 families (observations)

# 'data.frame':	180 obs. of  4 variables:
#   $ old_family  : Factor w/ 179 levels "Acanthaceae",..: 15 23 24 25 27 29 33 35 39 42 ...
# $ apgIV_family: Factor w/ 178 levels "Acanthaceae",..: 15 23 24 25 27 29 33 35 38 41 ...
# $ tot_species : int  1 1 1 1 1 1 1 1 1 1 ...
# $ med_species : int  0 0 0 0 0 0 0 0 0 1 ...

attach(benin_med_flora)
head(benin_med_flora)

###################################################################################
## B1 - Moerman approach: Linear regression
###################################################################################

moer<-lm(med_species~tot_species-1)
summary(moer)

# Residuals:
#   Min       1Q   Median       3Q      Max 
# -25.8284  -0.3172  -0.1586   0.8414  13.6479 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# tot_species 0.158624   0.005695   27.85   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 3.508 on 179 degrees of freedom
# Multiple R-squared:  0.8125,	Adjusted R-squared:  0.8115 
# F-statistic: 775.8 on 1 and 179 DF,  p-value: < 2.2e-16

## residuals
resid_lm<-residuals(moer)
predict_lm<-predict(moer)

## S3 method for class 'glm'
# predict(object, newdata = NULL,
#         type = c("link", "response", "terms"),
#         se.fit = FALSE, dispersion = NULL, terms = NULL,
#         na.action = na.pass, ...)
# type: The type of prediction required. The default is on the scale of the linear predictors; the alternative "response" is on the scale of the response variable. Thus for a default binomial model the default predictions are of log-odds (probabilities on logit scale) and type = "response" gives the predicted probabilities. The "terms" option returns a matrix giving the fitted values of each term in the model formula on the linear predictor scale.

## A1. plot: Without transformation
op<-par(mfrow=c(2,2))
plot(tot_species, med_species, ylab="Medicinal species per family", xlab="Number of species per family", pch=16, col="orange", xlim=c(0,500))
abline(moer, col="red")  ## adding the fit line
## install.packages("calibrate")
library(calibrate) ## add label to data point
textxy(tot_species, med_species, apgIV_family) ## add labels to  plot's points

## close up plot [0-70]
plot(tot_species, med_species, ylab="Medicinal species per family", xlab="Number of species per family", pch=16, col="grey", ylim=c(0,20), xlim=c(0,70))
abline(moer, col="red")  ## adding the fit line
library(calibrate) ## add label to data point
textxy(tot_species, med_species, apgIV_family) ## add labels to  plot's points

## close up plot [0-20]
plot(tot_species, med_species, ylab="Medicinal species per family", xlab="Number of species per family", pch=16, col="orange", ylim=c(0,7), xlim=c(0,20))
abline(moer, col="red")  ## adding the fit line
library(calibrate) ## add label to data point
textxy(tot_species, med_species, apgIV_family, offset=0.6) ## add labels to  plot's points
par(op)

## Diagnostic plots
op<-par(mfrow=c(2,2))
plot(moer)
par(op)

## The diagnostic suggest a slight departure from normality and slightly increasing variance as function of fitted values, even though we have small sample size for families with larger number of species. We will log (x+1) transform the data to deal with non-normallity as done in Moerman papers. We added +1 given that some families may have no medicinal plant species.

###################################################################################
## B2 - Moerman approach: Linear regression after log-log transformation
###################################################################################

lmoer<-lm(log(med_species+1)~log(tot_species+1)-1)

## Note here that to aid interpretation we fit the model through the origin (intercetp is 0), that is why I added -1 to the model

## Diagnostic plots
op<-par(mfrow=c(2,2))
plot(lmoer)
par(op)

## This yields better diagnostic plot and suggests that the normality and homogeneity of variance assumptions have been met. Let's then look at the results and extact the residuals

summary(lmoer)

# Residuals:
#   Min      1Q  Median      3Q     Max 
# -1.2860 -0.5096 -0.2987  0.3521  1.5999 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# log(tot_species + 1)  0.46383    0.01946   23.84   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.5725 on 179 degrees of freedom
# Multiple R-squared:  0.7605,	Adjusted R-squared:  0.7592 
# F-statistic: 568.4 on 1 and 179 DF,  p-value: < 2.2e-16

## The summary suggests significant positive relationship between the number of species per family and number of medicinal species. The slope of the relationship beta=0.46383 +/- 0.01946 (t= 23.84, p<0.0001). In average, one would expect nearly half (0.46383 on a log scale) of the total number of species in a family used for medicinal purpose.

predict_loglm<-predict(lmoer)
resid_loglm<-residuals(lmoer)
resid_family<-data.frame(apgIV_family, resid_loglm)
write.csv(resid_family, "Residual_by_Family.csv")

## ## A2. plot: With log transformation  *********
## using calibrate package: Within the calibrate package, the textxy() function can be used to label a plot's data points. The textxy() function accepts the following arugments ("Label points in a plot," n.d.).

#plot(total_species, med_species, log="xy", ylab="Medicinal species per family", xlab="Number of species per family")

op<-par(mfrow=c(2,2))
plot(log(tot_species+1), log(med_species+1), ylab="Medicinal species per family (log)", xlab="Number of species per family (log)", ylim=c(0,5), xlim=c(0.5,7), pch=20, col=hcl(h=130, c=50, l=60, alpha=1), cex = 1.5)
abline(lmoer, col=hcl(h=10, c=1, l=55, alpha=1), lwd=0.75)  ## adding the fit line
library(calibrate) ## add label to data point
textxy(log(tot_species+1), log(med_species+1), apgIV_family, offset=0.6) ## add labels to the preexisting plot's points

## closeup [0-3]
plot(log(tot_species+1), log(med_species+1), ylab="Medicinal species per family (log)", xlab="Number of species per family (log)", ylim=c(0,1.5), xlim=c(0.5,3), pch=20, col=hcl(h=130, c=50, l=60, alpha=1), cex = 1.5)
abline(lmoer, col=hcl(h=10, c=1, l=55, alpha=1), lwd=0.75)  ## adding the fit line
library(calibrate) ## add label to data point
textxy(log(tot_species+1), log(med_species+1), apgIV_family, offset=0.6) ## add labels to the preexisting plot's points

hist(resid_loglm, prob=TRUE, ylab="", xlab="Residuals",col=hcl(h=230, c=50, l=70, alpha=1), main="", ylim=c(0,1)); box()
par(op)


###################################################################################
## C - Alternative approch: poisson and negative binomial
###################################################################################
library(MASS)

## We fit a poisson model first
poim<-glm(med_species~tot_species, poisson)
summary(poim)

## there is a strong overdispersion, deviance/df= 693.57/178 = 3.89 >>>>1. We will use negative binomial error structure instead to account for overdispersion

poi_nb<-glm.nb(med_species~tot_species-1)
summary(poi_nb)

# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -3.3279  -1.1525  -0.2183   0.2667   1.9986  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# tot_species 0.029158   0.001883   15.48   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for Negative Binomial(0.7614) family taken to be 1)
# 
# Null deviance: 430.67  on 180  degrees of freedom
# Residual deviance: 171.61  on 179  degrees of freedom
# AIC: 637.53

## This suggest a significant positive link between total family size and number of medicinal plant species. beta=0.029158 +/- 0.001883 (Z=15.48, p<0.00001, AIC: 637.53). The mean number of medicinal species is 1.029587=exp(0.029158) per family.

## The residuals are:
predict_nb<-predict(poi_nb,  type="response")
resid_nb<-residuals(poi_nb)
resid_nb_family<-data.frame(apgIV_family, resid_nb)
write.csv(resid_nb_family, "Residual_NegBin.csv")

## Figure 1
ttsp<-range(tot_species)
tt_species<-seq(min(ttsp), max(ttsp), 1)

dev.off()
#op<-par(mfrow=c(2,1), pty="s", las=1)
#par(pty="s", las=1)
op<-par(mfrow=c(1,2), mgp=c(2,0.5,0),mar=c(3,3,1,1),oma=c(1,2,0,0),
        cex=0.9 ,cex.lab=1.1,cex.axis=1.2) #Four graphs
## Fig1-A
plot(tot_species, med_species, pch=16, col=hcl(h=230, c=150, l=70, alpha=1), xlim=c(0,600), ylim=c(0,90), xlab="", ylab="")
#points(total_species,poi_nb$fitted.values, col="red", type = "l")
points(tt_species, exp(poi_nb$coefficients*tt_species), col="red", type="l", lwd=2)
textxy(tot_species, med_species, apgIV_family, offset=0.6) ## add labels to the preexisting plot's points
abline(lm(med_species~tot_species-1), lwd=2, col="blue", lty=1) ## linear model
points(tt_species, exp(lmoer$coefficients*log(tt_species)), col="purple", type="l", lwd=2)  ## adding the fit line
mtext("(a)", side=3, adj=0)

## Fig1-B
#par(mai=c(.35, .9, .2, 0.25)) ## c(bottom, left, top, right)
plot(tot_species, med_species, pch=16, col=hcl(h=230, c=150, l=70, alpha=1), xlim=c(0,100), ylim=c(0,20), xlab="", ylab="")
#points(total_species,poi_nb$fitted.values, col="red", type = "l")
points(tt_species, exp(poi_nb$coefficients*tt_species), col="red", lwd=2, type="l")
textxy(tot_species, med_species, apgIV_family, offset=0.6) ## add labels to the preexisting plot's points
abline(lm(med_species~tot_species-1), lwd=2, col="blue", lty=1) ## linear model
points(tt_species, exp(lmoer$coefficients*log(tt_species)), lwd=2, col="purple", type="l")  ## adding the fit line
mtext("(b)", side=3, adj=0)
legend(0,20, c("Negative bionomial", "Log-transformed", "Linear"), lty=c(1,1,1), col=c("red","purple", "blue"), lwd=c(2,2,2), cex=0.9, bty="n")

par(oma=c(4, 5, 1, 0)) ## outer margin (bottom, left, top, right)
title(xlab="Number of plant species per family", ylab="Number of medicinal species per family", outer=T, cex.lab=1.25)
par(op)

### Fig 2
op<-par(mfrow=c(1,2))
hist(resid_loglm, prob=TRUE, ylab="Probability density", xlab="Residuals",col=hcl(h=230, c=50, l=70, alpha=1), main="Log normal", ylim=c(0,1)); box()

hist(resid_nb, prob=TRUE, ylab="", xlab="Residuals",col=hcl(h=230, c=50, l=70, alpha=1), main="Negative binomial", ylim=c(0,1)); box()
par(op)

## Comparing the residuals of linear models (moerman) with the negative binomial model
resids_family<-data.frame(apgIV_family, resid_nb, predict_nb, resid_loglm, predict_loglm, resid_lm, predict_lm)
write.csv(resids_family, "Residuals_all.csv")

#########################################################################
### RCEAP: residuals by hand
#########################################################################
mean_negbin<-exp(poi_nb$coefficients*tot_species)
mean_loglm<-exp(lmoer$coefficients*log(tot_species)) ## backtransformed
mean_lm<-moer<-moer$coefficients*tot_species  ## linear without transformation

resid_negbin<-med_species-mean_negbin
resid_loglm<-med_species-mean_loglm
resid_lm<-med_species-mean_lm
residuals_all<-data.frame(apgIV_family, med_species, tot_species, resid_negbin,resid_loglm,resid_lm)
write.csv(residuals_all, "Residuals_for_3_models.csv")

########################################################################
### Figure 2
########################################################################
overunder<-read.csv("/Users/ogaoue/Dropbox/PROJECTS/Benin Med Flora/R Script/over_under_family_Long.csv", header = T)
## Use Data: TableS6.csv *******
str(overunder)

library(ggplot2)

#install.packages("ggpubr")
# Install
#if(!require(devtools)) install.packages("devtools")
#devtools::install_github("kassambara/ggpubr")

install.packages("ggsci")

library(ggpubr) ## allows you to use ggarrange() to organize plots
library("ggsci") ## allows you to use the jco pallette of colors
over<-subset(overunder,overunder$Classification=="Overutilized")
under<-subset(overunder,overunder$Classification=="Underutilized")

# A bar graph
ggplot(data=overunder, aes(x=reorder(Family, Residuals), y=Residuals, fill=Model)) + 
  facet_wrap(~Classification,  scales ="free_y") +
  geom_bar(colour="black", stat="identity",
           position=position_dodge(),
           size=.3) +                  # Thinner lines
  scale_fill_hue(name="Model") +       # Set legend title
  scale_fill_jco()+
  scale_color_jco()+
  xlab("") + ylab("Residuals") + # Set axis labels
  ggtitle("") +                        # Set title
  theme_bw() + 
  coord_flip()


nb<-subset(overunder,overunder$Model=="NegBin")
ggplot(data=nb, aes(x=reorder(Family, Residuals), y=Residuals)) + 
  facet_wrap(~Classification,  scales ="free_y") +
  geom_bar(colour="black", stat="identity",
           position=position_dodge(),
           size=.3) +                  # Thinner lines
  scale_fill_hue(name="Model") +       # Set legend title
  scale_fill_jco()+
  xlab("") + ylab("Residuals") + # Set axis labels
  ggtitle("") +                        # Set title
  theme_bw() + 
  coord_flip()

################################################################################
## PART B: Phylogenetics signal analysis                                      ##
################################################################################
## Phylogenetic signal tests

# First read in the arguments listed at the command line
args=(commandArgs(TRUE))
library(picante)
library(ape)
library(adephylo)
library(ade4)
library(phylobase)
library(geiger)
library(phytools)
library(dplyr)

tre_med <- read.nexus("Table_S7.txt") ## phylogeny
#plot(tre_med,show.tip.label=F) # read phylogenetic tree

# data that contain medicinal variables
data_med <- read.table("Medicinal uses_2.txt",header=TRUE) 
## Data; Table_S5

attach(data_med)
names(data_med)

t1 <- tre_med
m <- data_med
s1 <- data.frame(species=t1$tip.label)

colnames(m)
m1 <- unique(m[, c(1,7)]) # to select the corresponding column to the plant organ to represent on the phylogeny, here column 7

colnames(m1)[[1]] <- "species"

d1 <- merge(m1, s1, by="species")
head(d1)

z <- as.data.frame(d1)
rownames(z) <- z[,1]
z[,1] <- NULL

#write.csv(z, "SA/Species_with_specimens_SA.csv")
subphy <- match.phylo.data(t1, z)$phy
subdat <- match.phylo.data(t1, z)$data

names(subdat)[1] <- "x"
y <- as.numeric(subdat$x)
names(y) <- rownames(subdat)
b <- di2multi(subphy) # ensuring the tree is fully resolved
z1 <- contMap(b, y, type="fan", outline=T, fsize=0.0000000001,main="root", legend=FALSE,lwd=2)

## To calculate D --> phylogenetic signal test in organ selection
library(caper)

data_med <- read.table("Medicinal uses_2.txt",header=TRUE)
## Data; Table_S5
attach(data_med)
names(data_med)
rownames(data_med) <- data_med$species

med_comp <- comparative.data(tre_med, data_med, species)
d_root <- phylo.d(med_comp, binvar=root) # to get the D value and p value for root selection
d_bark <- phylo.d(med_comp, binvar=bark)
d_leave <- phylo.d(med_comp, binvar=leave)
d_fruit <- phylo.d(med_comp, binvar=fruit)
d_stem <- phylo.d(med_comp, binvar=stem)
d_use <- phylo.d(med_comp, binvar=utilization_status)

print(d_root) # to show D and p values
print(d_bark)
print(d_leave)
print(d_fruit)
print(d_stem)
print(d_use)
