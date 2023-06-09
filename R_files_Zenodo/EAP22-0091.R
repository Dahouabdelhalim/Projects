#Set the working directory
setwd("D:/Rdata/wan")

bio.wf <- read.csv("D:/Rdata/wan/native_biomass.csv")
bio.wf<-bio.wf[(!is.na(bio.wf$native_biomass)),] #?Þ³?na????
View(bio.wf)

library(lattice)
library(nlme)
library(stats)
library(base)
library(dplyr)
library(knitr)
library(MuMIn)
library(Rmisc)
library(multcomp)
library(ggpubr)
library(funrar)
library(lme4)


str(bio.wf)
bio.wf$community_diversity=as.factor(bio.wf$community_diversity)
bio.wf$community_density=as.factor(bio.wf$community_density)
bio.wf$nutrient=as.factor(bio.wf$nutrient)
bio.wf$severed=as.factor(bio.wf$severed)
str(bio.wf)
#Try including variance structure
vf1<-varIdent(form=~1|nutrient)

bio.wf <- bio.wf %>%
          filter(!is.na(bio.wf$alien_root_ratio)) %>%
          droplevels()
        plot.resid<-function(model,breaks=30){
          plot(resid(model)~fitted(model), xlab = "Fitted values", ylab="Residuals")
          hist(resid(model),breaks=breaks)
          qqnorm(resid(model))
          qqline(resid(model))
        }


#----------------------------------------------end-----------------------------------------------------------------------#


#data translate
#origin
u<-gl(1,1,length(bio.wf$alien_root_ratio))
LCVlme01 <- lme(alien_root_ratio~community_diversity+community_density+nutrient+severed+community_diversity:community_density+community_diversity:nutrient+community_diversity:severed+community_density:nutrient+community_density:severed+nutrient:severed+community_diversity:community_density:nutrient+community_diversity:community_density:severed+community_density:nutrient:severed+community_diversity:nutrient:severed+community_diversity:community_density:nutrient:severed,random=list(u=pdBlocked(list(
                      pdIdent(form=~alien_identity-1),
                      pdIdent(form=~native_community-1)
                    ))), weights=vf1, data=bio.wf,method="REML")
summary(LCVlme01)
r.squaredGLMM(LCVlme01)
opm1a <- par(mfrow = c(2,3), mar = c(4,4,2,2))
plot.resid(LCVlme01)
hist(bio.wf$alien_root_ratio)
#sqrt
u<-gl(1,1,length(bio.wf$alien_root_ratio))
LCVlme01_1 <- lme(sqrt(alien_root_ratio)~community_diversity+community_density+nutrient+severed+community_diversity:community_density+community_diversity:nutrient+community_diversity:severed+community_density:nutrient+community_density:severed+nutrient:severed+community_diversity:community_density:nutrient+community_diversity:community_density:severed+community_density:nutrient:severed+community_diversity:nutrient:severed+community_diversity:community_density:nutrient:severed,random=list(u=pdBlocked(list(
                      pdIdent(form=~alien_identity-1),
                      pdIdent(form=~native_community-1)
                    ))), weights=vf1, data=bio.wf,method="REML")
summary(LCVlme01_1)
r.squaredGLMM(LCVlme01_1)
opm1a <- par(mfrow = c(2,3), mar = c(4,4,2,2))
plot.resid(LCVlme01_1)
hist(sqrt(bio.wf$alien_root_ratio))
#cubic
u<-gl(1,1,length(bio.wf$alien_root_ratio))
LCVlme01_2 <- lme((alien_root_ratio)^(1/3)~community_diversity+community_density+nutrient+severed+community_diversity:community_density+community_diversity:nutrient+community_diversity:severed+community_density:nutrient+community_density:severed+nutrient:severed+community_diversity:community_density:nutrient+community_diversity:community_density:severed+community_density:nutrient:severed+community_diversity:nutrient:severed+community_diversity:community_density:nutrient:severed,random=list(u=pdBlocked(list(
                      pdIdent(form=~alien_identity-1),
                      pdIdent(form=~native_community-1)
                    ))), weights=vf1, data=bio.wf,method="REML")
summary(LCVlme01_2)
r.squaredGLMM(LCVlme01_2)
opm1a <- par(mfrow = c(2,3), mar = c(4,4,2,2))
plot.resid(LCVlme01_2)
hist((bio.wf$alien_root_ratio)^(1/3))
#log
u<-gl(1,1,length(bio.wf$alien_root_ratio))
LCVlme01_3 <- lme(log(alien_root_ratio)~community_diversity+community_density+nutrient+severed+community_diversity:community_density+community_diversity:nutrient+community_diversity:severed+community_density:nutrient+community_density:severed+nutrient:severed+community_diversity:community_density:nutrient+community_diversity:community_density:severed+community_density:nutrient:severed+community_diversity:nutrient:severed+community_diversity:community_density:nutrient:severed,random=list(u=pdBlocked(list(
                      pdIdent(form=~alien_identity-1),
                      pdIdent(form=~native_community-1)
                    ))), weights=vf1, data=bio.wf,method="REML")
summary(LCVlme01_3)
r.squaredGLMM(LCVlme01_3)
opm1a <- par(mfrow = c(2,3), mar = c(4,4,2,2))
plot.resid(LCVlme01_3)
hist(log(bio.wf$alien_root_ratio))
#logit
u<-gl(1,1,length(bio.wf$alien_root_ratio))
LCVlme01_4 <- lme(logit(alien_root_ratio)~community_diversity+community_density+nutrient+severed+community_diversity:community_density+community_diversity:nutrient+community_diversity:severed+community_density:nutrient+community_density:severed+nutrient:severed+community_diversity:community_density:nutrient+community_diversity:community_density:severed+community_density:nutrient:severed+community_diversity:nutrient:severed+community_diversity:community_density:nutrient:severed,random=list(u=pdBlocked(list(
                      pdIdent(form=~alien_identity-1),
                      pdIdent(form=~native_community-1)
                    ))), weights=vf1, data=bio.wf,method="REML")
summary(LCVlme01_4)
r.squaredGLMM(LCVlme01_4)
opm1a <- par(mfrow = c(2,3), mar = c(4,4,2,2))
plot.resid(LCVlme01_4)
hist(logit(bio.wf$alien_root_ratio))

# Now I choose cubic translate, but if necessary, we can delete some data,then choose logit translate 
LCVlme02 <- lme((alien_root_ratio)^(1/3)~community_diversity+community_density+nutrient+severed+community_diversity:community_density+community_diversity:nutrient+community_diversity:severed+community_density:nutrient+community_density:severed+nutrient:severed+community_diversity:community_density:nutrient+community_diversity:community_density:severed+community_density:nutrient:severed+community_diversity:nutrient:severed+community_diversity:community_density:nutrient:severed,random=list(u=pdBlocked(list(
                      pdIdent(form=~alien_identity-1),
                      pdIdent(form=~native_community-1)
                    ))), data=bio.wf,method="REML",control=list(maxIter=1000,msMaxIter=1000,niterEM=1000,opt="optim")) 

summary(LCVlme02)



#Use AIC to select best model
anova(LCVlme01_2,LCVlme02)
#LCVlme01 seems to be best, but actually not to be different with LCVlme02.

#read the data file into R:

################################################
#test fixed effects on native_biomass, use method="ML"#
################################################
LCV2lme01a =lme((alien_root_ratio)^(1/3)~community_diversity+community_density+nutrient+severed+community_diversity:community_density+community_diversity:nutrient+community_diversity:severed+community_density:nutrient+community_density:severed+nutrient:severed+community_diversity:community_density:nutrient+community_diversity:community_density:severed+community_diversity:nutrient:severed+community_density:nutrient:severed+community_diversity:community_density:nutrient:severed,random=list(u=pdBlocked(list(
                      pdIdent(form=~alien_identity-1),
                      pdIdent(form=~native_community-1)
                    ))), weights=vf1, data=bio.wf,method="ML",control=list(maxIter=1000,msMaxIter=1000,niterEM=1000,opt="optim"))
summary(LCV2lme01a)

#community_diversity:community_density:nutrient:severed
LCV2lme01b =lme((alien_root_ratio)^(1/3)~community_diversity+community_density+nutrient+severed+community_diversity:community_density+community_diversity:nutrient+community_diversity:severed+community_density:nutrient+community_density:severed+nutrient:severed+community_diversity:community_density:nutrient+community_diversity:community_density:severed+community_diversity:nutrient:severed+community_density:nutrient:severed,random=list(u=pdBlocked(list(pdIdent(form=~alien_identity-1),
                     pdIdent(form=~native_community-1)
                    ))), weights=vf1, data=bio.wf,method="ML",control=list(maxIter=1000,msMaxIter=1000,niterEM=1000,opt="optim"))
anova(LCV2lme01a,LCV2lme01b)   

#community_density:nutrient:severed
LCV2lme01c =lme((alien_root_ratio)^(1/3)~community_diversity+community_density+nutrient+severed+community_diversity:community_density+community_diversity:nutrient+community_diversity:severed+community_density:nutrient+community_density:severed+nutrient:severed+community_diversity:community_density:nutrient+community_diversity:community_density:severed+community_diversity:nutrient:severed, random=list(u=pdBlocked(list(pdIdent(form=~alien_identity-1),
                     pdIdent(form=~native_community-1)
                    ))), weights=vf1, data=bio.wf,method="ML",control=list(maxIter=1000,msMaxIter=1000,niterEM=1000,opt="optim"))
anova(LCV2lme01b,LCV2lme01c)   

#community_diversity:nutrient:severed
LCV2lme01d =lme((alien_root_ratio)^(1/3)~community_diversity+community_density+nutrient+severed+community_diversity:community_density+community_diversity:nutrient+community_diversity:severed+community_density:nutrient+community_density:severed+nutrient:severed+community_diversity:community_density:nutrient+community_diversity:community_density:severed,random=list(u=pdBlocked(list(pdIdent(form=~alien_identity-1),
                     pdIdent(form=~native_community-1)
                    ))), weights=vf1, data=bio.wf,method="ML",control=list(maxIter=1000,msMaxIter=1000,niterEM=1000,opt="optim"))
anova(LCV2lme01c,LCV2lme01d) 

#community_diversity:community_density:severed
LCV2lme01e =lme((alien_root_ratio)^(1/3)~community_diversity+community_density+nutrient+severed+community_diversity:community_density+community_diversity:nutrient+community_diversity:severed+community_density:nutrient+community_density:severed+nutrient:severed+community_diversity:community_density:nutrient,random=list(u=pdBlocked(list(pdIdent(form=~alien_identity-1),
                     pdIdent(form=~native_community-1)
                    ))), weights=vf1, data=bio.wf,method="ML",control=list(maxIter=1000,msMaxIter=1000,niterEM=1000,opt="optim"))
anova(LCV2lme01d,LCV2lme01e) 

#community_diversity:community_density:nutrient
LCV2lme01f =lme((alien_root_ratio)^(1/3)~community_diversity+community_density+nutrient+severed+community_diversity:community_density+community_diversity:nutrient+community_diversity:severed+community_density:nutrient+community_density:severed+nutrient:severed,random=list(u=pdBlocked(list(pdIdent(form=~alien_identity-1),
                     pdIdent(form=~native_community-1)
                    ))), weights=vf1, data=bio.wf,method="ML",control=list(maxIter=1000,msMaxIter=1000,niterEM=1000,opt="optim")) 
anova(LCV2lme01e,LCV2lme01f) 

#nutrient:severed
LCV2lme01g =lme((alien_root_ratio)^(1/3)~community_diversity+community_density+nutrient+severed+community_diversity:community_density+community_diversity:nutrient+community_diversity:severed+community_density:nutrient+community_density:severed,random=list(u=pdBlocked(list(pdIdent(form=~alien_identity-1),
                     pdIdent(form=~native_community-1)
                    ))), weights=vf1, data=bio.wf,method="ML",control=list(maxIter=1000,msMaxIter=1000,niterEM=1000,opt="optim"))
anova(LCV2lme01f,LCV2lme01g) 

#community_density:severed
LCV2lme01h =lme((alien_root_ratio)^(1/3)~community_diversity+community_density+nutrient+severed+community_diversity:community_density+community_diversity:nutrient+community_diversity:severed+community_density:nutrient,random=list(u=pdBlocked(list(pdIdent(form=~alien_identity-1),
                     pdIdent(form=~native_community-1)
                    ))), weights=vf1, data=bio.wf,method="ML",control=list(maxIter=1000,msMaxIter=1000,niterEM=1000,opt="optim")) 
anova(LCV2lme01g,LCV2lme01h) 

#community_density:nutrient
LCV2lme01i =lme((alien_root_ratio)^(1/3)~community_diversity+community_density+nutrient+severed+community_diversity:community_density+community_diversity:nutrient+community_diversity:severed,random=list(u=pdBlocked(list(pdIdent(form=~alien_identity-1),
                     pdIdent(form=~native_community-1)
                    ))), weights=vf1, data=bio.wf,method="ML",control=list(maxIter=1000,msMaxIter=1000,niterEM=1000,opt="optim"))
anova(LCV2lme01h,LCV2lme01i)

#community_diversity:severed
LCV2lme01j =lme((alien_root_ratio)^(1/3)~community_diversity+community_density+nutrient+severed+community_diversity:community_density+community_diversity:nutrient,random=list(u=pdBlocked(list(pdIdent(form=~alien_identity-1),
                     pdIdent(form=~native_community-1)
                    ))), weights=vf1, data=bio.wf,method="ML",control=list(maxIter=1000,msMaxIter=1000,niterEM=1000,opt="optim")) 
anova(LCV2lme01i,LCV2lme01j) 

#community_diversity:nutrient
LCV2lme01k =lme((alien_root_ratio)^(1/3)~community_diversity+community_density+nutrient+severed+community_diversity:community_density,random=list(u=pdBlocked(list(pdIdent(form=~alien_identity-1),
                     pdIdent(form=~native_community-1)
                    ))), weights=vf1, data=bio.wf,method="ML",control=list(maxIter=1000,msMaxIter=1000,niterEM=1000,opt="optim")) 
anova(LCV2lme01j,LCV2lme01k) 

#community_diversity:community_density
LCV2lme01l =lme((alien_root_ratio)^(1/3)~community_diversity+community_density+nutrient+severed,random=list(u=pdBlocked(list(pdIdent(form=~alien_identity-1),
                     pdIdent(form=~native_community-1)
                    ))), weights=vf1, data=bio.wf,method="ML",control=list(maxIter=1000,msMaxIter=1000,niterEM=1000,opt="optim"))
anova(LCV2lme01k,LCV2lme01l) 

#severed
LCV2lme01m =lme((alien_root_ratio)^(1/3)~community_diversity+community_density+nutrient,random=list(u=pdBlocked(list(pdIdent(form=~alien_identity-1),
                     pdIdent(form=~native_community-1)
                    ))), weights=vf1, data=bio.wf,method="ML",control=list(maxIter=1000,msMaxIter=1000,niterEM=1000,opt="optim"))
anova(LCV2lme01l,LCV2lme01m) 

#nutrient
LCV2lme01n=lme((alien_root_ratio)^(1/3)~community_diversity+community_density+severed,random=list(u=pdBlocked(list(pdIdent(form=~alien_identity-1),
                     pdIdent(form=~native_community-1)
                    ))), weights=vf1, data=bio.wf,method="ML",control=list(maxIter=1000,msMaxIter=1000,niterEM=1000,opt="optim"))
anova(LCV2lme01l,LCV2lme01n) 

#community_density
LCV2lme01o=lme((alien_root_ratio)^(1/3)~community_diversity+nutrient+severed,random=list(u=pdBlocked(list(pdIdent(form=~alien_identity-1),
                     pdIdent(form=~native_community-1)
                    ))),weights=vf1, data=bio.wf,method="ML",control=list(maxIter=1000,msMaxIter=1000,niterEM=1000,opt="optim"))
anova(LCV2lme01l,LCV2lme01o) 

#community_diversity
LCV2lme01p=lme((alien_root_ratio)^(1/3)~community_density+nutrient+severed,random=list(u=pdBlocked(list(pdIdent(form=~alien_identity-1),pdIdent(form=~native_community-1)))), weights=vf1, data=bio.wf,method="ML",control=list(maxIter=1000,msMaxIter=1000,niterEM=1000,opt="optim")) 
anova(LCV2lme01l,LCV2lme01p)