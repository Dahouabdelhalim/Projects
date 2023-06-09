library("ctv")
install.views("Phylogenetics")
library(ape)
library(caper)
library(geiger)
library(nlme)
library(phangorn)
library(phytools)
library(phylolm)
library(corHMM)
library(gee)
library(car)
library(emmeans)
library(tidyverse)
library(broom)

#Reading in the rooted tree#
RallidTree<-read.tree("Data/Rallid_Tree.tre")
RallidTree
plotTree(RallidTree,ftype="i",fsize=0.6,lwd=1)
Ntip(RallidTree)

#Making a slightly different tree to remove species without duet information
CertainTree <- drop.tip(RallidTree, c("Rallina_canningi", "Laterallus_levraudi",
                                      "Laterallus_ruber", "Habroptila_wallacii",
                                      "Hypotaenidia_owstoni", 
                                      "Pardirallus_maculatus", 
                                      "Hapalocrex_flaviventer", 
                                      "Laterallus_spilopterus", "Porzana_fluminea",
                                      "Zapornia_akool", "Zapornia_bicolor",
                                      "Paragallinula_angulata",
                                      "Tribonyx_ventralis", "Fulica_cornuta"))
CertainTree
plotTree(CertainTree,ftype="i",fsize=0.6,lwd=1)
Ntip(CertainTree)

#Checking for polytomies#
is.binary.tree(CertainTree)

#Checking for ultramectricity
is.ultrametric(CertainTree)
UltraTree<-force.ultrametric(CertainTree)
is.ultrametric(UltraTree)

#Phylogenetic generalized least-squares regressions of discrete data in APE and caper
Rallid_Reduced<-read.csv("Data/RallidReducedData.csv",header=TRUE,row.names=1)
Rallid_Reduced

#Comparing species in tree to data
NameCheck<-name.check(CertainTree,Rallid_Reduced)
NameCheck


#Phylogenetic generalized linear models in phylolm for binary traits

#Get variables from continuous to discrete factors
Rallid_Reduced$Habitat<-as.factor(Rallid_Reduced$Habitat)
Rallid_Reduced$SocialBond<-as.factor(Rallid_Reduced$SocialBond)
Rallid_Reduced$Territory<-as.factor(Rallid_Reduced$Territory)
Rallid_Reduced$CoopBreed<-as.factor(Rallid_Reduced$CoopBreed)
Rallid_Reduced$Migrant<-as.factor(Rallid_Reduced$Migrant)
Rallid_Reduced$Bioregion<-as.factor(Rallid_Reduced$Bioregion)

#Checking for multicollinearity
vf.lm = lm(Duets~Habitat+SocialBond+Territory+CoopBreed+Migrant+
                                   Bioregion+Latitude+SSD, data=Rallid_Reduced)
vif(mod=vf.lm)
vf.lm2 = lm(Duets~Habitat+SocialBond+Territory+CoopBreed+Migrant+
             Latitude+SSD, data=Rallid_Reduced)
vif(mod=vf.lm2)

#Running the models and comparing AIC and log likelihood values
pglmm.duets0 = phyloglm(Duets~Habitat + SocialBond + Territory + CoopBreed +
                        Migrant + Latitude + SSD, phy=CertainTree,
                        method = "logistic_MPLE", data=Rallid_Reduced, btol=30,
                        boot = 1000)
pglmm.duets0
summary(pglmm.duets0)
AIC(pglmm.duets0)
plot(pglmm.duets0)

pglmm.duets1 = phyloglm(Duets~Habitat + SocialBond + Territory + CoopBreed +
                        Migrant + Latitude, phy=CertainTree,
                        method = "logistic_MPLE", data=Rallid_Reduced, btol=30, 
                        boot = 1000)
pglmm.duets1
summary(pglmm.duets1)
AIC(pglmm.duets1)
plot(pglmm.duets1)

LRTteststat1 <- -2 * (-23.07--22.99)
p.val1 <- pchisq(LRTteststat1, df = 1, lower.tail = FALSE)

pglmm.duets2 = phyloglm(Duets~Habitat + CoopBreed + Territory + 
                          Migrant + Latitude, phy=CertainTree,
                        method = "logistic_MPLE", data=Rallid_Reduced, btol=30,
                        boot = 1000)
pglmm.duets2
summary(pglmm.duets2)
AIC(pglmm.duets2)
plot(pglmm.duets2)

LRTteststat2 <- -2 * (-23.43--23.07)
p.val2 <- pchisq(LRTteststat2, df = 1, lower.tail = FALSE)

pglmm.duets3 = phyloglm(Duets~Habitat + CoopBreed + Territory + Migrant, 
                        phy=CertainTree, 
                        method = "logistic_MPLE", data=Rallid_Reduced, btol=30, 
                        boot = 1000)
pglmm.duets3
summary(pglmm.duets3)
AIC(pglmm.duets3)
plot(pglmm.duets3)

LRTteststat3 <- -2 * (-23.56--23.43)
p.val3 <- pchisq(LRTteststat3, df = 1, lower.tail = FALSE)

pglmm.duets4 = phyloglm(Duets~Habitat + Migrant + Territory, 
                          phy=CertainTree, 
                        method = "logistic_MPLE", data=Rallid_Reduced, btol=30, 
                        boot = 1000)
pglmm.duets4
summary(pglmm.duets4)
AIC(pglmm.duets4)
plot(pglmm.duets4)

LRTteststat4 <- -2 * (-23.42--23.56)
p.val4 <- pchisq(LRTteststat4, df = 1, lower.tail = FALSE)

pglmm.duets5 = phyloglm(Duets~Habitat + Territory, phy=CertainTree, 
                        method = "logistic_MPLE", data=Rallid_Reduced, btol=30, 
                        boot = 1000)
pglmm.duets5
summary(pglmm.duets5)
AIC(pglmm.duets5)
plot(pglmm.duets5)

LRTteststat5 <- -2 * (-27.26--23.42)
p.val5 <- pchisq(LRTteststat5, df = 1, lower.tail = FALSE)

pglmm.duets6 = phyloglm(Duets~Habitat + Migrant, phy=CertainTree, 
                        method = "logistic_MPLE", data=Rallid_Reduced, btol=30, 
                        boot = 1000)
pglmm.duets6
summary(pglmm.duets6)
AIC(pglmm.duets6)
plot(pglmm.duets6)

LRTteststat6 <- -2 * (-30.46--23.42)
p.val6 <- pchisq(LRTteststat6, df = 1, lower.tail = FALSE)

pglmm.duets7 = phyloglm(Duets~Migrant + Territory, phy=CertainTree, 
                        method = "logistic_MPLE", data=Rallid_Reduced, btol=30, 
                        boot = 1000)
pglmm.duets7
summary(pglmm.duets7)
AIC(pglmm.duets7)
plot(pglmm.duets7)

LRTteststat7 <- -2 * (-26.13--23.42)
p.val7 <- pchisq(LRTteststat7, df = 1, lower.tail = FALSE)