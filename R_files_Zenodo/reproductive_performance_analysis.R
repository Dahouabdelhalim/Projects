# Script for:
# Fitness and fur colouration - testing the camouflage and thermoregulation hypotheses in an Arctic mammal
# Cecilia Di Bernardi, Anne-Mathilde Thierry, Nina E. Eide, Diana E. Bowler, Lars Rød-Eriksen, Stefan Blumentrath, Lukas Tietgen, Brett K. Sandercock, Øystein Flagstad, Arild Landa

#This R script reads in the associated data files on dryad and includes the code for the statistical analysis of reproductive performance (breeding propensity, litter size)

library(dplyr)


########## 1) BREEDING PROPENSITY

#import dataset
br = read.csv2("breeding.csv") 
br$rodent <- as.factor(br$rodent)

#GAMs
library(mgcv)
library(bbmle)
library(wiqid)
library(MuMIn)

#models GAM Binomial

gam_null<-gam(detected_breeding~ 1,select=T,data=br,family = "binomial")

#rodent+colour
gamA<-gam(detected_breeding~rodent+colour+s(age,k=6),select=T,data=br,family = "binomial")

#rodent*colour
gamAA<-gam(detected_breeding~ rodent*colour+s(age,k=6),select=T,data=br,family = "binomial")

#climate+rodent+color*sex
gamB<-gam(detected_breeding~ snow1+rodent+colour*sex+s(age,k=6),select=T,data=br,family = "binomial")

gamBB<-gam(detected_breeding~ snow2+rodent+colour*sex+s(age,k=6),select=T,data=br,family = "binomial")

gamBBB<-gam(detected_breeding~ avgwin +rodent+colour*sex+s(age,k=6),select=T,data=br,family = "binomial")

#climate+colour*sex
gamC<-gam(detected_breeding~ snow1+colour*sex+s(age,k=6),select=T,data=br,family = "binomial")

gamCC<-gam(detected_breeding~ snow2+colour*sex+s(age,k=6),select=T,data=br,family = "binomial")

gamCCCC<-gam(detected_breeding~ avgwin+colour*sex+s(age,k=6),select=T,data=br,family = "binomial")

#colour
gam_D <- gam(detected_breeding~colour+ s(age,k=6),select=T,data=br,family = "binomial")

#climate+rodent+origin+colour
gamY<-gam(detected_breeding~ snow1+rodent+origin+colour+s(age,k=6),select=T,data=br,family = "binomial")

gamYY<-gam(detected_breeding~ snow2+rodent+origin+colour+s(age,k=6),select=T,data=br,family = "binomial")

gamYYY<-gam(detected_breeding~ avgwin+rodent+origin+colour+s(age,k=6),select=T,br,family = "binomial")

#climate+rodent
gamJ <-gam(detected_breeding~ snow1+rodent+s(age,k=6),select=T,data=br,family = "binomial")

gamJJ <-gam(detected_breeding~ snow2+rodent+s(age,k=6),select=T,data=br,family = "binomial")

gamJJJ <-gam(detected_breeding~ avgwin+rodent+s(age,k=6),select=T,data=br,family = "binomial")

#rodent+colour*climate+sex+origin
gamT<-gam(detected_breeding~ rodent+colour*snow1+sex+origin+s(age,k=6),select=T,data=br,family = "binomial")

gamTT<-gam(detected_breeding~ rodent+colour*snow2+sex+origin+s(age,k=6),select=T,data=br,family = "binomial")

gamTTTT<-gam(detected_breeding~ rodent+colour*avgwin+sex+origin+s(age,k=6),select=T,data=br,family = "binomial")

#climate+rodent+colour
gamK<-gam(detected_breeding~ snow1+rodent+colour+s(age,k=6),select=T,data=br,family = "binomial")

gamKK<-gam(detected_breeding~ snow2+rodent+colour+s(age,k=6),select=T,data=br,family = "binomial")

gamKKK<-gam(detected_breeding~ avgwin+rodent+colour+s(age,k=6),select=T,data=br,family = "binomial")


#model selection: AICc
AICc_values_revision1 <- MuMIn::AICc(gam_null,gamA,gamAA,gam_D,gamY,gamYY,gamYYY,gamJ,gamJJ,gamJJJ,gamT,gamTT,gamTTTT,
                                     gamK,gamKK,gamKKK,gamB,gamBB,gamBBB,gamC,gamCC,gamCCCC) 
AICc_tavola <- AICtable(AICc_values_revision1, digits =3)




########## 2) LITTER SIZE

#subset the dataset to narrow down to the breeding events
litter <- br %>% filter(detected_breeding==1)

#models GAMM Poisson

#null model
gamm_null <- gamm(sizelitter ~1,random=list(breeding_id=~1),data=litter,family = "poisson")

#colour+origin+sex+rodent
gamm1 <-gamm(sizelitter~colour+origin+sex+rodent+s(age,k=4),random=list(breeding_id=~1),data=litter,family = "poisson")

#colour*climate+origin+sex+rodent
gamm2 <-gamm(sizelitter~colour*snow1+origin+sex+rodent+s(age,k=4),random=list(breeding_id=~1),data=litter,family = "poisson")

gamm22 <-gamm(sizelitter~colour*snow2+origin+sex+rodent+s(age,k=4),random=list(breeding_id=~1),data=litter,family = "poisson")

gamm222 <-gamm(sizelitter~colour*avgwin+origin+sex+rodent+s(age,k=4),random=list(breeding_id=~1),data=litter,family = "poisson")

#colour*climate+origin+rodent
gamm3 <-gamm(sizelitter~colour*snow1+origin+rodent+s(age,k=4),random=list(breeding_id=~1),data=litter,family = "poisson")

gamm33 <-gamm(sizelitter~colour*snow2+origin+rodent+s(age,k=4),random=list(breeding_id=~1),data=litter,family = "poisson")

gamm333 <-gamm(sizelitter~colour*avgwin+origin+rodent+s(age,k=4),random=list(breeding_id=~1),data=litter,family = "poisson")

#colour*climate+rodent
gamm4 <-gamm(sizelitter~colour*snow1+rodent+s(age,k=4),random=list(breeding_id=~1),data=litter,family = "poisson")

gamm44 <-gamm(sizelitter~colour*snow2+rodent+s(age,k=4),random=list(breeding_id=~1),data=litter,family = "poisson")

gamm444 <-gamm(sizelitter~colour*avgwin+rodent+s(age,k=4),random=list(breeding_id=~1),data=litter,family = "poisson")

#climate+rodent
gamm5 <-gamm(sizelitter~snow1+rodent+s(age,k=4),random=list(breeding_id=~1),data=litter,family = "poisson")

gamm55 <-gamm(sizelitter~snow2+rodent+s(age,k=4),random=list(breeding_id=~1),data=litter,family = "poisson")

gamm555 <-gamm(sizelitter~avgwin+rodent+s(age,k=4),random=list(breeding_id=~1),data=litter,family = "poisson")

#climate+rodent+colour*sex
gamm7 <-gamm(sizelitter~snow1+rodent+colour*sex+s(age,k=4),random=list(breeding_id=~1),data=litter,family = "poisson")

gamm77 <-gamm(sizelitter~snow2+rodent+colour*sex+s(age,k=4),random=list(breeding_id=~1),data=litter,family = "poisson")

gamm777 <-gamm(sizelitter~avgwin+rodent+colour*sex+s(age,k=4),random=list(breeding_id=~1),data=litter,family = "poisson")

#climate+rodent+colour+origin
gamm9 <-gamm(sizelitter~snow1+rodent+colour+origin+s(age,k=4),random=list(breeding_id=~1),data=litter,family = "poisson")

gamm99 <-gamm(sizelitter~snow2+rodent+colour+origin+s(age,k=4),random=list(breeding_id=~1),data=litter,family = "poisson")

gamm999 <-gamm(sizelitter~avgwin+rodent+colour+origin+s(age,k=4),random=list(breeding_id=~1),data=litter,family = "poisson")

#colour
gamm100 <-gamm(sizelitter~colour+s(age,k=4),random=list(breeding_id=~1),data=litter,family = "poisson")

#colour+rodent
gamm101 <-gamm(sizelitter~colour+rodent+s(age,k=4),random=list(breeding_id=~1),data=litter,family = "poisson")


#model selection: AICc
AICc_values_litter <- MuMIn::AICc(gamm1,gamm2,gamm22,gamm222,gamm3,gamm33,gamm333,gamm4,gamm44,gamm444,
                                  gamm5,gamm55,gamm555,gamm7,gamm77,gamm777,gamm9,gamm99,gamm999,gamm100,gamm101, gamm_null)

AICc_table_litter <- AICtable(AICc_values_litter) 




