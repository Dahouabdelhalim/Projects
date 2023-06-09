############################    SETTTING UP THE ANALYSIS     ####################################################

#load libraries
library(lme4)
library(ggplot2)
library(reshape2)
library(MASS)
library(car)
library(Hmisc)
library(corrplot) 
library(gridExtra)
library(extrafont)
library(glmmADMB) 
library(influence.ME)
library(DHARMa)

#Load dataset
pan<- read.csv("~/PhD/Journal of Applied Ecology/Dryaddata/beebowldata.R.na.csv")
str(pan) #view variables

#y variables, rename/characterize
bee.ab<-pan$bee.ab
bee.div<-pan$bee.div
bee.rich<-pan$bee.rich
cwm<-pan$CWM
lasioglossum<-pan$lasioglossum
miner<-pan$miner
renter<-pan$renter

#collection variables, rename/characterize
year<-as.factor(pan$year)
season<-as.factor(pan$season)
site<-pan$site
trt<-as.factor(pan$treat)
hood<-as.factor(pan$hood)
n.pans<-pan$n.pans

#explanatory, rename/characterize
lpi.g<-pan$lpi.green
ENN<-pan$enn.mn
gs.size<-pan$gs.area
gs.ca<-pan$gs.ca

bloom.ab<-as.integer(pan$bloom.ab)
bloom.ar<-pan$bloom.ar
vdiv<-pan$vdiv
vheight<-pan$vheight
biomass<-pan$biomass

#################MULTICOLINEARITY CHECK #######################################
require(corrplot)
my_data<-cbind(bloom.ab,bloom.ar,vdiv,vheight,biomass,ENN,gs.size,lpi.g,gs.ca) #create data subset

print(my_data)
na.exclude(my_data)
res1<-rcorr(as.matrix(my_data))
res1
res1$r
res1$P
#view variable correlations with a plot
corrplot(res1$r, type="upper", order="hclust", p.mat = res1$P, sig.level = 0.01, insig = "blank")

#check Variance Inflation Factors
vif(lm(bee.ab ~ ENN+lpi.g+gs.size+bloom.ab+bloom.ar+vdiv+vheight+biomass+gs.ca))

#take out gs.ca which has a too high VIF and is correlated with ENN (see plot)
vif(lm(bee.ab ~ ENN+lpi.g+gs.size+bloom.ab+bloom.ar+vdiv+vheight+biomass))
#all VIF < 2, therefore okay to include in analysis
#correlation between biomass and vegheight is 0.54 but ok

#new plot
my_data<-cbind(bloom.ab,bloom.ar,vdiv,vheight,biomass,ENN,gs.size,lpi.g,hood) 

print(my_data)
na.exclude(my_data)
res1<-rcorr(as.matrix(my_data))
res1
res1$r
res1$P
corrplot(res1$r, type="upper", order="hclust", p.mat = res1$P, sig.level = 0.01, insig = "blank")

############# centering,scaling,creating data ###################################################################
newdata<-cbind(lpi.g,ENN,gs.size,bloom.ab,bloom.ar,vdiv,vheight,biomass,bee.div)
mydataedited<-as.data.frame(newdata)

#center and scale variables
mydataed<-scale(mydataedited, center=TRUE, scale=TRUE)
mydata.edited<-as.data.frame(mydataed)
str(mydata.edited)

#rename variables
lpi.g<-mydata.edited$lpi.g 
ENN<-mydata.edited$ENN
gs.size<-mydata.edited$gs.size
bloom.ab<-mydata.edited$bloom.ab
bloom.ar<-mydata.edited$bloom.ar
vdiv<-mydata.edited$vdiv
vheight<-mydata.edited$vheight
biomass<-mydata.edited$biomass
bee.div<-mydata.edited$bee.div

#create data set for analysis with centered/scaled explanatory variables
newdata2<-cbind(bee.ab,bee.div,bee.rich,cwm,miner,renter,lasioglossum,lpi.g,ENN,gs.size,bloom.ab,bloom.ar,vdiv,vheight,biomass,year,season,trt,n.pans,site,hood)
mydata<-as.data.frame(newdata2)
mydata<-transform(mydata,hood=factor(hood))
mydata<-transform(mydata,year=factor(year))
mydata<-transform(mydata,season=factor(season))
str(mydata) #check variables

mydata <- na.exclude(mydata)
head(mydata)
na.action(mydata)

#########################OVERDISPERSION CHECK#######################################################################

#created a dispersion check function that will determine dispersion parameters for most models
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval) }

## Testing distributions for full models across all response variables

### Response variable: Bee ABUNDANCE
poissontest<-glmer(bee.ab~ENN+lpi.g+gs.size+bloom.ar+bloom.ab+vdiv+vheight+biomass+season+year+(1|hood)+offset(log(n.pans)),family=poisson,data=mydata)
summary(poissontest)
overdisp_fun(poissontest) #dispersion parameter is 429, overdispersed

nbtest<-glmer.nb(bee.ab~ENN+lpi.g+gs.size+bloom.ar+bloom.ab+vdiv+vheight+biomass+season+year+(1|hood)+offset(log(n.pans)),data=mydata)
summary(nbtest)
overdisp_fun(nbtest) #dispersion parameter is 1.03, use negative binomial distribution

### Response variable: Bee RICHNESS
poissontest<-glmer(bee.rich~ENN+lpi.g+gs.size+bloom.ar+bloom.ab+vdiv+vheight+biomass+season+year+(1|hood)+offset(log(n.pans)),family=poisson,data=mydata)
summary(poissontest)
overdisp_fun(poissontest) #dispersion parameter is 187, overdispersed

nbtest<-glmer.nb(bee.rich~ENN+lpi.g+gs.size+bloom.ar+bloom.ab+vdiv+vheight+biomass+season+year+(1|hood)+offset(log(n.pans)),data=mydata)
summary(nbtest)
overdisp_fun(nbtest) #dispersion parameter is 1.06, use negative binomial distribution

### Response variable: RENTER Bees (cavity nesting)
poissontest<-glmer(renter~ENN+lpi.g+gs.size+bloom.ar+bloom.ab+vdiv+vheight+biomass+season+year+(1|hood)+offset(log(n.pans)),family=poisson,data=mydata)
summary(poissontest)
overdisp_fun(poissontest) #dispersion parameter is 281, overdispersed

nbtest<-glmer.nb(renter~ENN+lpi.g+gs.size+bloom.ar+bloom.ab+vdiv+vheight+biomass+season+year+(1|hood)+offset(log(n.pans)),data=mydata)
summary(nbtest)
overdisp_fun(nbtest)#dispersion parameter is 0.96, use negative binomial distribution

### Response variable: MINER Bees (ground nesting)
poissontest<-glmer(miner~ENN+lpi.g+gs.size+bloom.ar+bloom.ab+vdiv+vheight+biomass+season+year+(1|hood)+offset(log(n.pans)),family=poisson,data=mydata)
summary(poissontest)
overdisp_fun(poissontest) #dispersion parameter is 256, overdispersed

nbtest<-glmer.nb(miner~ENN+lpi.g+gs.size+bloom.ar+bloom.ab+vdiv+vheight+biomass+year+season+(1|hood)+offset(log(n.pans)),data=mydata)
summary(nbtest)
overdisp_fun(nbtest) #dispersion parameter is 1.08, use negative binomial distribution

### Response variable: Bee DIVERSITY
#continuous response variable so used lmer instead of glmer
lmtest<-lmer(bee.div~ENN+lpi.g+gs.size+bloom.ar+bloom.ab+vdiv+vheight+biomass+season+year+(1|hood)+offset(log(n.pans)),data=mydata)
summary(lmtest)
plot(lmtest) #normality looks ok

### Response variable: CWM or community weighted mean of bee body size
####### subsetting data for examining bee body size (CWM) to exclude all 0s or absences of bees
cwmdata<-subset(mydata, cwm > 0)

#continuous response variable so used lmer instead of glmer
lmtest<-lmer(cwm~ENN+lpi.g+gs.size+bloom.ar+bloom.ab+vdiv+vheight+year+biomass+season+(1|hood)+offset(log(n.pans)),data=cwmdata)
summary(lmtest) 
plot(lmtest) #normality looks good

## Response Variable: Lasioglossum abundances
## lasioglossum
poissontest<-glmer(lasioglossum~ENN+lpi.g+gs.size+bloom.ar+bloom.ab+vdiv+vheight+biomass+season+year+(1|hood)+offset(log(n.pans)),family=poisson,data=mydata)
summary(poissontest)
overdisp_fun(poissontest) #dispersion parameter is 237, overdispersed

nbtest<-glmer.nb(lasioglossum~ENN+lpi.g+gs.size+bloom.ar+bloom.ab+vdiv+vheight+biomass+season+year+(1|hood)+offset(log(n.pans)),data=mydata)
summary(nbtest)
overdisp_fun(nbtest)#dispersion parameter is 1.03, use negative binomial distribution

################ TREATMENT T-tests ###########################################
#Model 1: Bee Abundance
treatment.ab<-glmer.nb(bee.ab~trt+(1|hood)+offset(log(n.pans)),data=mydata)
control.ab<-glmer.nb(bee.ab~1+(1|hood)+offset(log(n.pans)),data=mydata)
anova(treatment.ab,control.ab) # no significant difference

#model 2: Bee Richness
treatment.rich<-glmer.nb(bee.rich~trt+(1|hood)+offset(log(n.pans)),data=mydata)
control.rich<-glmer.nb(bee.rich~1+(1|hood)+offset(log(n.pans)),data=mydata)
anova(treatment.rich,control.rich) #no significant difference

#model3: Renter bees
treatment.rent<-glmer.nb(renter~trt+(1|hood)+offset(log(n.pans)),data=mydata)
control.rent<-glmer.nb(renter~1+(1|hood)+offset(log(n.pans)),data=mydata)
anova(treatment.rent,control.rent) #no significant difference

#model4: Miner bees
treatment.mine<-glmer.nb(miner~trt+(1|hood)+offset(log(n.pans)),data=mydata)
control.mine<-glmer.nb(miner~1+(1|hood)+offset(log(n.pans)),data=mydata)
anova(treatment.mine,control.mine) #no significant difference

#model5: Bee diversity
treatment.div<-lmer(bee.div~trt+(1|hood)+offset(log(n.pans)),data=mydata)
control.div<-lmer(bee.div~1+(1|hood)+offset(log(n.pans)),data=mydata)
anova(treatment.div,control.div) #no significant difference

#model6: CWM
treatment.cwm<-lmer(cwm~trt+(1|hood)+offset(log(n.pans)),data=cwmdata)
control.cwm<-lmer(cwm~1+(1|hood)+offset(log(n.pans)),data=cwmdata)
anova(treatment.cwm,control.cwm) #no significant difference

#model7: Lasioglossum
treatment.las<-lmer(lasioglossum~trt+(1|hood)+offset(log(n.pans)),data=mydata)
control.las<-lmer(lasioglossum~1+(1|hood)+offset(log(n.pans)),data=mydata)
anova(treatment.las,control.las) #no significant difference

############# BACKWARDS MODEL SELECTION########################################

############# FIRST RESPONSE VARIABLE : BEE ABUNDANCE
#first step
full.model<-glmer.nb(bee.ab~ENN+lpi.g+gs.size+bloom.ar+bloom.ab+vdiv+vheight+year+season+biomass+(1|hood)+offset(log(n.pans)),data = mydata) 
summary(full.model) # remove bloom.ar variable

modelb<-glmer.nb(bee.ab~ENN+lpi.g+gs.size+bloom.ab+vdiv+vheight+year+season+biomass+(1|hood)+offset(log(n.pans)),data = mydata) 
anova(modelb,full.model) #no significant difference

#secondstep
summary(modelb) #remove vdiv variable

modelc<-glmer.nb(bee.ab~ENN+lpi.g+gs.size+bloom.ab+vheight+year+season+biomass+(1|hood)+offset(log(n.pans)),data = mydata) 
anova(modelc,modelb) #no significant difference

#thirdstep
summary(modelc) #remove ENN variable

modeld<-glmer.nb(bee.ab~lpi.g+gs.size+bloom.ab+vheight+year+season+biomass+(1|hood)+offset(log(n.pans)),data = mydata)
anova(modelc,modeld) #no significant difference

#fourthstep
summary(modeld) #remove bloom.ab variable

modele<-glmer.nb(bee.ab~lpi.g+gs.size+vheight+year+season+biomass+(1|hood)+offset(log(n.pans)),data = mydata)
anova(modeld,modele) #no significant difference

#fifthstep
summary(modele) #remove gs.size variable

modelf<-glmer.nb(bee.ab~lpi.g+vheight+year+season+biomass+(1|hood)+offset(log(n.pans)),data = mydata)
anova(modele,modelf) #no significant difference

#sixthstep
summary(modelf) #remove lpi.g variable

modelg<-glmer.nb(bee.ab~vheight+year+season+biomass+(1|hood)+offset(log(n.pans)),data = mydata)
anova(modelf,modelg) #no significant difference

#seventhstep
summary(modelg) #all variables significant

modelh<-glmer.nb(bee.ab~1+(1|hood)+offset(log(n.pans)),data = mydata)#intercept only model
anova(modelg,modelh) #significant difference, keep model g

### BEST Bee abundance model
abundancemodel<-glmer.nb(bee.ab~vheight+year+season+biomass+(1|hood)+offset(log(n.pans)),data = mydata)
summary(abundancemodel)
plot(abundancemodel)#looks ok
inf=influence(abundancemodel,obs=T)
plot(inf,which="cook")##no points violate cooks distance of >0.5

#residual diagnostics, confirming model performance okay
fittedModel<-glmer.nb(bee.ab~vheight+year+season+biomass+(1|hood)+offset(log(n.pans)),data = mydata)
simulationOutput <- simulateResiduals(fittedModel = fittedModel)
residuals(simulationOutput)
plot(simulationOutput)
testDispersion(simulationOutput)

############# SECOND RESPONSE VARIABLE : BEE RICHNESS
#first step
full.model<-glmer.nb(bee.rich~ENN+lpi.g+gs.size+bloom.ar+bloom.ab+vdiv+vheight+biomass+year+season+(1|hood)+offset(log(n.pans)),data = mydata) 
summary(full.model) # remove vdiv variable

modelb<-glmer.nb(bee.rich~ENN+lpi.g+gs.size+bloom.ar+bloom.ab+vheight+year+season+biomass+(1|hood)+offset(log(n.pans)),data = mydata) 
anova(modelb,full.model) #no significant difference

#secondstep
summary(modelb) #remove ENN variable

modelc<-glmer.nb(bee.rich~lpi.g+gs.size+bloom.ar+bloom.ab+vheight+year+season+biomass+(1|hood)+offset(log(n.pans)),data = mydata) 
anova(modelc,modelb) #no significant difference

#thirdstep
summary(modelc) #remove bloom.ab variable

modeld<-glmer.nb(bee.rich~lpi.g+gs.size+bloom.ar+vheight+year+season+biomass+(1|hood)+offset(log(n.pans)),data = mydata) 
anova(modelc,modeld) #no significant difference

#fourthstep
summary(modeld) #remove gs.size variable

modele<-glmer.nb(bee.rich~lpi.g+bloom.ar+vheight+year+season+biomass+(1|hood)+offset(log(n.pans)),data = mydata) 
anova(modeld,modele) #no significant difference

#fifthstep
summary(modele) #remove bloom.ar variable

modelf<-glmer.nb(bee.rich~lpi.g+vheight+year+season+biomass+(1|hood)+offset(log(n.pans)),data = mydata) 
anova(modele,modelf) #no significant difference

#sixthstep
summary(modelf) #all variables significant

modelg<-glmer.nb(bee.rich~1+(1|hood)+offset(log(n.pans)),data = mydata) #intercept only model
anova(modelf,modelg) #significant difference, model better than intercept only model

##BEST BEE RICHNESS MODEL:
richnessmodel<-glmer.nb(bee.rich~lpi.g+vheight+year+season+biomass+(1|hood)+offset(log(n.pans)),data = mydata) 
summary(richnessmodel)
plot(richnessmodel)
inf=influence(richnessmodel,obs=T)
plot(inf,which="cook")#no points violate cooks distance of >0.5

#residual diagnostics, confirming model performance okay
fittedModel<-glmer.nb(bee.rich~lpi.g+vheight+year+season+biomass+(1|hood)+offset(log(n.pans)),data = mydata) 
simulationOutput <- simulateResiduals(fittedModel = fittedModel)
residuals(simulationOutput)
plot(simulationOutput)
testDispersion(simulationOutput)

############# THIRD RESPONSE VARIABLE : RENTER BEES

#first step
full.model<-glmer.nb(renter~ENN+lpi.g+gs.size+bloom.ar+bloom.ab+vdiv+vheight+year+season+biomass+(1|hood)+offset(log(n.pans)),data = mydata) 
summary(full.model) # remove gs.size variable

modelb<-glmer.nb(renter~ENN+lpi.g+bloom.ar+bloom.ab+vdiv+vheight+year+season+biomass+(1|hood)+offset(log(n.pans)),data = mydata) 
anova(modelb,full.model) #no significant difference

#secondstep
summary(modelb) #remove biomass variable

modelc<-glmer.nb(renter~ENN+lpi.g+bloom.ar+bloom.ab+vdiv+vheight+year+season+(1|hood)+offset(log(n.pans)),data = mydata) 
anova(modelc,modelb) #no significant difference

#thirdstep
summary(modelc) #remove year variable 

modeld<-glmer.nb(renter~ENN+lpi.g+bloom.ar+bloom.ab+vdiv+vheight+season+(1|hood)+offset(log(n.pans)),data = mydata) 
anova(modelc,modeld) #no significant difference

#fourthstep
summary(modeld) #remove bloom.area variable

modele<-glmer.nb(renter~ENN+lpi.g+bloom.ab+vdiv+vheight+season+(1|hood)+offset(log(n.pans)),data = mydata) 
anova(modeld,modele) #no significant difference

#fifthstep
summary(modele) #remove lpi.g variable

modelf<-glmer.nb(renter~ENN+bloom.ab+vdiv+vheight+season+(1|hood)+offset(log(n.pans)),data = mydata) 
anova(modele,modelf) #no significant difference

#sixthstep
summary(modelf) #remove veg height variable 

modelg<-glmer.nb(renter~ENN+bloom.ab+vdiv+season+(1|hood)+offset(log(n.pans)),data = mydata) 
anova(modelf,modelg) #no significant difference

#seventhstep
summary(modelg) #remove vdiv variable

modelh<-glmer.nb(renter~ENN+bloom.ab+season+(1|hood)+offset(log(n.pans)),data = mydata) 
anova(modelg,modelh) #no significant difference

#eighth step
summary(modelh) #remove ENN

modeli<-glmer.nb(renter~bloom.ab+season+(1|hood)+offset(log(n.pans)),data = mydata) 
anova(modelh,modeli) #no significant difference

#ninth step
summary(modeli) #remove bloom.ab

modelj<-glmer.nb(renter~season+(1|hood)+offset(log(n.pans)),data = mydata) 
anova(modeli,modelj) #no significant difference

#tentrh step
summary(modelj) #everything significant

intercept<-glmer.nb(renter~1+(1|hood)+offset(log(n.pans)),data = mydata) 
anova(modelj,intercept) #significant difference, keep model j

### BEST Bee renter model
rentermodel<-glmer.nb(renter~season+(1|hood)+offset(log(n.pans)),data = mydata)
summary(rentermodel)
plot(rentermodel) #looks ok
inf=influence(rentermodel,obs=T)
plot(inf,which="cook")##no points violate cooks distance of >0.5

#residual diagnostics, confirming model performance okay
fittedModel<-glmer.nb(renter~season+(1|hood)+offset(log(n.pans)),data = mydata)
simulationOutput <- simulateResiduals(fittedModel = fittedModel)
residuals(simulationOutput)
plot(simulationOutput)
testDispersion(simulationOutput)

############# FOURTH RESPONSE VARIABLE : miner BEES

#first step
full.model<-glmer.nb(miner~ENN+lpi.g+gs.size+bloom.ar+bloom.ab+vdiv+vheight+year+season+biomass+(1|hood)+offset(log(n.pans)),data = mydata) 
summary(full.model) # remove bloom area variable

modelb<-glmer.nb(miner~ENN+lpi.g+gs.size+bloom.ab+vdiv+vheight+year+season+biomass+(1|hood)+offset(log(n.pans)),data = mydata) 
anova(modelb,full.model) #no significant difference

#secondstep
summary(modelb) #remove ENN variable

modelc<-glmer.nb(miner~lpi.g+gs.size+bloom.ab+vdiv+vheight+year+season+biomass+(1|hood)+offset(log(n.pans)),data = mydata) 
anova(modelc,modelb) #no significant difference

#thirdstep
summary(modelc) #remove bloom.ab variable

modeld<-glmer.nb(miner~lpi.g+gs.size+vdiv+vheight+year+season+biomass+(1|hood)+offset(log(n.pans)),data = mydata) 
anova(modelc,modeld) #no significant difference

#fourthstep
summary(modeld) #remove vdiv variable

modele<-glmer.nb(miner~lpi.g+gs.size+vheight+year+season+biomass+(1|hood)+offset(log(n.pans)),data = mydata) 
anova(modeld,modele) #no significant difference

#fifthstep
summary(modele) #remove biomass variable

modelf<-glmer.nb(miner~lpi.g+gs.size+vheight+year+season+(1|hood)+offset(log(n.pans)),data = mydata) 
anova(modele,modelf) #no significant difference

#sixthstep
summary(modelf) #remove vheight variable

modelg<-glmer.nb(miner~lpi.g+gs.size+year+season+(1|hood)+offset(log(n.pans)),data = mydata) 
anova(modelf,modelg) #no significant difference

#seventhstep
summary(modelg) #remove gs.size variable

modelh<-glmer.nb(miner~lpi.g+year+season+(1|hood)+offset(log(n.pans)),data = mydata) 
anova(modelg,modelh) #no significant difference

#eighth step
summary(modelh) #remove lpi.g

modeli<-glmer.nb(miner~year+season+(1|hood)+offset(log(n.pans)),data = mydata) 
anova(modelh,modeli) #no significant difference

#ninthstep
summary(modeli)# all variables significant

intercept<-glmer.nb(miner~1+(1|hood)+offset(log(n.pans)),data = mydata) #intercept only model
anova(modeli,intercept) #significant difference, keep model i

### BEST Bee miner model
minermodel<-glmer.nb(miner~year+season+(1|hood)+offset(log(n.pans)),data = mydata)
summary(minermodel)
plot(minermodel) #looks ok
inf=influence(minermodel,obs=T)
plot(inf,which="cook") ###no points violate cooks distance of >0.5

#residual diagnostics, confirming model performance okay
fittedModel<-glmer.nb(miner~year+season+(1|hood)+offset(log(n.pans)),data = mydata)
simulationOutput <- simulateResiduals(fittedModel = fittedModel)
residuals(simulationOutput)
plot(simulationOutput)
testDispersion(simulationOutput)

############# FIFTH RESPONSE VARIABLE : cwm BEES
cwmdata<-subset(mydata, cwm > 0)
#first step
full.model<-lmer(cwm~ENN+lpi.g+gs.size+bloom.ar+bloom.ab+vdiv+vheight+year+season+biomass+(1|hood)+offset(log(n.pans)),data = cwmdata) 
summary(full.model) # remove ENN variable

modelb<-lmer(cwm~lpi.g+gs.size+bloom.ar+bloom.ab+vdiv+vheight+year+season+biomass+(1|hood)+offset(log(n.pans)),data = cwmdata) 
anova(modelb,full.model) #no significant difference

#secondstep
summary(modelb) #remove vdiv variable

modelc<-lmer(cwm~lpi.g+gs.size+bloom.ar+bloom.ab+vheight+year+season+biomass+(1|hood)+offset(log(n.pans)),data = cwmdata) 
anova(modelc,modelb) #no significant difference

#thirdstep
summary(modelc) #remove season variable

modeld<-lmer(cwm~lpi.g+gs.size+bloom.ar+bloom.ab+vheight+year+biomass+(1|hood)+offset(log(n.pans)),data = cwmdata) 
anova(modelc,modeld) #no significant difference

#fourthstep
summary(modeld) #remove gs.size variable

modele<-lmer(cwm~lpi.g+bloom.ar+bloom.ab+vheight+year+biomass+(1|hood)+offset(log(n.pans)),data = cwmdata) 
anova(modeld,modele) #no significant difference

#fifthstep
summary(modele) #remove bloom.ab variable

modelf<-lmer(cwm~lpi.g+bloom.ar+vheight+year+biomass+(1|hood)+offset(log(n.pans)),data = cwmdata) 
anova(modele,modelf) #no significant difference

#sixthstep
summary(modelf) #remove biomass variable

modelg<-lmer(cwm~lpi.g+bloom.ar+vheight+year+(1|hood)+offset(log(n.pans)),data = cwmdata) 
anova(modelf,modelg) #significant difference, keep model f

#sixthstep
summary(modelg) #all significant

intercept<-lmer(cwm~1+(1|hood)+offset(log(n.pans)),data = cwmdata) 
anova(modelg,intercept) #significant difference, keep model g

### BEST Bee cwm model
cwmmodel<-lmer(cwm~lpi.g+bloom.ar+vheight+year+(1|hood)+offset(log(n.pans)),data = cwmdata)
summary(cwmmodel)

#normality check
plot(cwmmodel) # looks great

#checking for influential points
inf=influence(cwmmodel,obs=T)
plot(inf,which="cook",sort=TRUE,cexaxis=1.5) #no points violate cooks distance of >0.5

############# SIXTH RESPONSE VARIABLE : Bee Diversity

#first step
full.model<-lmer(bee.div~ENN+lpi.g+gs.size+bloom.ar+bloom.ab+vdiv+vheight+year+season+biomass+(1|hood)+offset(log(n.pans)),data = mydata) 
summary(full.model) # remove vdiv variable

modelb<-lmer(bee.div~ENN+lpi.g+gs.size+bloom.ar+bloom.ab+vheight+year+season+biomass+(1|hood)+offset(log(n.pans)),data = mydata) 
anova(modelb,full.model) #no significant difference

#secondstep
summary(modelb) #remove bloom.ab variable

modelc<-lmer(bee.div~ENN+lpi.g+gs.size+bloom.ar+vheight+year+season+biomass+(1|hood)+offset(log(n.pans)),data = mydata) 
anova(modelc,modelb) #no significant difference

#thirdstep
summary(modelc) #remove ENN variable

modeld<-lmer(bee.div~lpi.g+gs.size+bloom.ar+vheight+year+season+biomass+(1|hood)+offset(log(n.pans)),data = mydata) 
anova(modelc,modeld) #no significant difference

#fourthstep
summary(modeld) #remove gs.size variable

modele<-lmer(bee.div~lpi.g+bloom.ar+vheight+year+season+biomass+(1|hood)+offset(log(n.pans)),data = mydata) 
anova(modeld,modele) #no significant difference

#fifthstep
summary(modele) #remove bloom area variable

modelf<-lmer(bee.div~lpi.g+vheight+year+season+biomass+(1|hood)+offset(log(n.pans)),data = mydata) 
anova(modele,modelf) #no significant difference

#sixthstep
summary(modelf) #remove biomass variable

modelg<-lmer(bee.div~lpi.g+vheight+year+season+(1|hood)+offset(log(n.pans)),data = mydata) 
anova(modelf,modelg) #no significant difference

#seventhstep
summary(modelg) #remove vegheight

modelh<-lmer(bee.div~lpi.g+year+season+(1|hood)+offset(log(n.pans)),data = mydata) 
anova(modelg,modelh) #no significant difference

#eighth step
summary(modelh) #remove season

modeli<-lmer(bee.div~lpi.g+year+(1|hood)+offset(log(n.pans)),data = mydata) 
anova(modelh,modeli) #significant difference, keep model h

#intercept check
intercept<-lmer(bee.div~1+(1|hood)+offset(log(n.pans)),data = mydata) 
anova(modelh,intercept) #keep model h

### BEST Bee bee.div model
bee.divmodel<-lmer(bee.div~lpi.g+year+season+(1|hood)+offset(log(n.pans)),data = mydata) 
summary(bee.divmodel)

#normality check
plot(bee.divmodel) #looks good

#checking for influential points
inf=influence(bee.divmodel,obs=T)
plot(inf,which="cook") #does not violate cook's distance

############# SEVENTH RESPONSE VARIABLE : Lasioglossum BEES

#first step
full.model<-glmer.nb(lasioglossum~ENN+lpi.g+gs.size+bloom.ar+bloom.ab+vdiv+vheight+year+biomass+season+(1|hood)+offset(log(n.pans)),data = mydata) 
summary(full.model) # remove vdiv variable

modelb<-glmer.nb(lasioglossum~ENN+lpi.g+gs.size+bloom.ar+bloom.ab+vheight+year+season+biomass+(1|hood)+offset(log(n.pans)),data = mydata) 
anova(modelb,full.model) #no significant difference

#secondstep
summary(modelb) #remove ENN variable

modelc<-glmer.nb(lasioglossum~lpi.g+gs.size+bloom.ar+bloom.ab+vheight+year+season+biomass+(1|hood)+offset(log(n.pans)),data = mydata) 
anova(modelc,modelb) #no significant difference

#thirdstep
summary(modelc) #remove bloom.ab variable

modeld<-glmer.nb(lasioglossum~lpi.g+gs.size+bloom.ar+vheight+year+season+biomass+(1|hood)+offset(log(n.pans)),data = mydata) 
anova(modelc,modeld) #no significant difference

#fourthstep
summary(modeld) #remove bloom.ar variable

modele<-glmer.nb(lasioglossum~lpi.g+gs.size+vheight+year+season+biomass+(1|hood)+offset(log(n.pans)),data = mydata) 
anova(modeld,modele) #no significant difference

#fifthstep
summary(modele) #remove gs.size variable

modelf<-glmer.nb(lasioglossum~lpi.g+vheight+year+season+biomass+(1|hood)+offset(log(n.pans)),data = mydata) 
anova(modele,modelf) #no significant difference

#sixthstep
summary(modelf) #remove lpi.g variable

modelg<-glmer.nb(lasioglossum~vheight+year+season+biomass+(1|hood)+offset(log(n.pans)),data = mydata) 
anova(modelf,modelg) #no significant difference

#seventhstep
summary(modelg) #all variables significant

intercept<-glmer.nb(lasioglossum~1+(1|hood)+offset(log(n.pans)),data = mydata) 
anova(modelg,intercept) #signfigant, keep model g

### BEST Bee lasioglossum model
lasioglossummodel<-glmer.nb(lasioglossum~vheight+biomass+year+season+(1|hood)+offset(log(n.pans)),data = mydata)
summary(lasioglossummodel)
plot(lasioglossummodel) #looks ok
inf=influence(lasioglossummodel,obs=T)
plot(inf,which="cook") #no points violate cooks distance of >0.5

#residual diagnostics, confirming model performance okay
fittedModel<-glmer.nb(lasioglossum~vheight+biomass+year+season+(1|hood)+offset(log(n.pans)),data = mydata)
simulationOutput <- simulateResiduals(fittedModel = fittedModel)
residuals(simulationOutput)
plot(simulationOutput)
testDispersion(simulationOutput)

########################### Obtaining partial residuals for plotting #############################

## RESPONSE VARIABLE ONE: BEE ABUNDANCE
abfit<-glmmadmb(bee.ab~vheight+year+season+biomass+(1|hood)+offset(log(n.pans)),family="nbinom",data=mydata)
X <- model.matrix(~vheight+year+season+biomass,data=mydata)
beta <- fixef(abfit)

## multiply each column of X by the corresponding beta 
beta_X <- sweep(X,MARGIN=2,STATS=beta,FUN="*")

## add residuals (columnwise addition works automatically, ## although we could use sweep with MARGIN=1) 
p_resid <- sweep(beta_X,MARGIN=1,STATS=residuals(abfit),FUN="+")
partialresiduals<-as.data.frame(p_resid)
partialresiduals

abvheightpresids<-partialresiduals$vheight
abbiomasspresids<-partialresiduals$biomass

## RESPONSE VARIABLE TWO: BEE RICHNESS
richfit <- glmmadmb(bee.rich~lpi.g+vheight+year+season+biomass+(1|hood)+offset(log(n.pans)),family="nbinom",data=mydata)
summary(richfit)
X <- model.matrix(~lpi.g+vheight+year+season+biomass,data=mydata)
beta <- fixef(richfit)

## multiply each column of X by the corresponding beta 
beta_X <- sweep(X,MARGIN=2,STATS=beta,FUN="*")

## add residuals (columnwise addition works automatically, ## although we could use sweep with MARGIN=1) 
p_resid <- sweep(beta_X,MARGIN=1,STATS=residuals(richfit),FUN="+")
partialresiduals<-as.data.frame(p_resid)
partialresiduals

rlpi.gpresids<-partialresiduals$lpi.g
rvheightpresids<-partialresiduals$vheight
rbiomasspresids<-partialresiduals$biomass

## FIFTH RESPONSE VARIABLE SIX: CWM bee size
cwmdata<-subset(mydata, cwm > 0)
cwmfit<-lmer(cwm~lpi.g+bloom.ar+vheight+year+(1|hood)+offset(log(n.pans)),data = cwmdata)
     
summary(cwmfit)
X <- model.matrix(~lpi.g+bloom.ar+vheight+year,data=cwmdata)
beta <- fixef(cwmfit)

## multiply each column of X by the corresponding beta 
beta_X <- sweep(X,MARGIN=2,STATS=beta,FUN="*")

## add residuals (columnwise addition works automatically, ## although we could use sweep with MARGIN=1) 
p_resid <- sweep(beta_X,MARGIN=1,STATS=residuals(cwmfit),FUN="+")
partialresiduals<-as.data.frame(p_resid)
partialresiduals

clpi.gpresids<-partialresiduals$lpi.g
cvheightpresids<-partialresiduals$vheight
cbloompresids<-partialresiduals$bloom.ar

### SIXTH RESPONSE VARIABLE, LASIOGLOSSUM
lasfit <- glmmadmb(lasioglossum~vheight+biomass+year+season+(1|hood)+offset(log(n.pans)),family="nbinom",data=mydata)
summary(lasfit)
X <- model.matrix(~vheight+biomass+year+season,data=mydata)
beta <- fixef(lasfit)

## multiply each column of X by the corresponding beta 
beta_X <- sweep(X,MARGIN=2,STATS=beta,FUN="*")

## add residuals (columnwise addition works automatically, ## although we could use sweep with MARGIN=1) 
p_resid <- sweep(beta_X,MARGIN=1,STATS=residuals(lasfit),FUN="+")
partialresiduals<-as.data.frame(p_resid)
partialresiduals

lvheightpresids<-partialresiduals$vheight
lbiomasspresids<-partialresiduals$biomass

#Reload dataset
pan<- read.csv("~/PhD/Journal of Applied Ecology/Dryaddata/beebowldata.R.na.csv")
str(pan) #view variables

#y variables, rename/characterize
bee.ab<-pan$bee.ab
bee.rich<-pan$bee.rich
cwm<-pan$CWM
lasioglossum<-pan$lasioglossum

#collection variables, rename/characterize
year<-as.factor(pan$year)
season<-as.factor(pan$season)
hood<-as.factor(pan$hood)
n.pans<-pan$n.pans

#explanatory, rename/characterize
lpi.g<-pan$lpi.green
bloom.ar<-pan$bloom.ar
vheight<-pan$vheight
biomass<-pan$biomass

#create subsetted dataset for plotting with original, not scaled features
plotdata<-cbind(bee.ab,bee.rich,cwm,lasioglossum,lpi.g,bloom.ar,vheight,biomass,year,season,n.pans,hood)
plotdata<-as.data.frame(plotdata)
plotdata<-transform(plotdata,hood=factor(hood))
plotdata<-transform(plotdata,year=factor(year))
plotdata<-transform(plotdata,season=factor(season))
str(plotdata) #check variables

plotdata <- na.exclude(plotdata)
na.action(plotdata)
cwmplotdata<-subset(plotdata, cwm > 0)

###   BEE ABUNDANCE PLOTS
abundvheight<-ggplot(plotdata, aes(x=vheight, y=abvheightpresids))+ geom_smooth(method=glm , col=colors()[153],fill=colors()[142],se=TRUE) + 
  geom_point(alpha=0.35,size=1.5)+theme_classic() +labs(x="\\nVegetation Height\\n", y="\\nPartial Residuals\\n") + 
  theme(axis.text=element_text(size=18,family="Calibri"), axis.title=element_text(size=20,family="Calibri",lineheight=.5), plot.title=element_text(size=22,face="bold",family="Calibri"))

abundbiomass<-ggplot(plotdata, aes(x=vheight, y=abbiomasspresids))+ geom_smooth(method=glm , col=colors()[153],fill=colors()[142],se=TRUE) + 
  geom_point(alpha=0.35,size=1.5)+theme_classic() +labs(x="\\nBiomass\\n", y="\\nPartial Residuals\\n") + 
  theme(axis.text=element_text(size=18,family="Calibri"), axis.title=element_text(size=20,family="Calibri",lineheight=.5), plot.title=element_text(size=22,face="bold",family="Calibri"))

###   BEE RICHNESS PLOTS
richlpi<-ggplot(plotdata, aes(x=lpi.g, y=rlpi.gpresids))+ geom_smooth(method=glm , col=colors()[153],fill=colors()[142],se=TRUE) + 
  geom_point(alpha=0.35,size=1.5)+theme_classic() +labs(x="\\nLargest Patch Size (LPI)\\n", y="\\nPartial Residuals\\n") + 
  theme(axis.text=element_text(size=18,family="Calibri"), axis.title=element_text(size=20,family="Calibri",lineheight=.5), plot.title=element_text(size=22,face="bold",family="Calibri"))

richvheight<-ggplot(plotdata, aes(x=vheight, y=rvheightpresids))+ geom_smooth(method=glm , col=colors()[153],fill=colors()[142],se=TRUE) + 
  geom_point(alpha=0.35,size=1.5)+theme_classic() +labs(x="\\nVegetation Height\\n", y="\\nPartial Residuals\\n") + 
  theme(axis.text=element_text(size=18,family="Calibri"), axis.title=element_text(size=20,family="Calibri",lineheight=.5), plot.title=element_text(size=22,face="bold",family="Calibri"))

richbiomass<-ggplot(plotdata, aes(x=vheight, y=rbiomasspresids))+ geom_smooth(method=glm , col=colors()[153],fill=colors()[142],se=TRUE) + 
  geom_point(alpha=0.35,size=1.5)+theme_classic() +labs(x="\\nBiomass\\n", y="\\nPartial Residuals\\n") + 
  theme(axis.text=element_text(size=18,family="Calibri"), axis.title=element_text(size=20,family="Calibri",lineheight=.5), plot.title=element_text(size=22,face="bold",family="Calibri"))

## CWM PLOTS
cwmlpi<-ggplot(cwmplotdata, aes(x=lpi.g, y=clpi.gpresids))+ geom_smooth(method=glm , col=colors()[153],fill=colors()[142],se=TRUE) + 
  geom_point(alpha=0.35,size=1.5)+theme_classic() +labs(x="\\nLargest Patch Size (LPI)\\n", y="\\nPartial Residuals\\n") + 
  theme(axis.text=element_text(size=18,family="Calibri"), axis.title=element_text(size=20,family="Calibri",lineheight=.5), plot.title=element_text(size=22,face="bold",family="Calibri"))

cwmvheight<-ggplot(cwmplotdata, aes(x=vheight, y=cvheightpresids))+ geom_smooth(method=glm , col=colors()[153],fill=colors()[142],se=TRUE) + 
  geom_point(alpha=0.35,size=1.5)+theme_classic() +labs(x="\\nVegetation Height\\n", y="\\nPartial Residuals\\n") + 
  theme(axis.text=element_text(size=18,family="Calibri"), axis.title=element_text(size=20,family="Calibri",lineheight=.5), plot.title=element_text(size=22,face="bold",family="Calibri"))

cwmbloom<-ggplot(cwmplotdata, aes(x=bloom.ar, y=cbloompresids))+ geom_smooth(method=glm , col=colors()[153],fill=colors()[142],se=TRUE) + 
  geom_point(alpha=0.35,size=1.5)+theme_classic() +labs(x="\\nBloom Area\\n", y="\\nPartial Residuals\\n") + 
  theme(axis.text=element_text(size=18,family="Calibri"), axis.title=element_text(size=20,family="Calibri",lineheight=.5), plot.title=element_text(size=22,face="bold",family="Calibri"))

##Lasioglossum plots
lasvheight<-ggplot(plotdata, aes(x=vheight, y=lvheightpresids))+ geom_smooth(method=glm , col=colors()[153],fill=colors()[142],se=TRUE) + 
  geom_point(alpha=0.35,size=1.5)+theme_classic() +labs(x="\\nVegetation Height\\n", y="\\nPartial Residuals\\n") + 
  theme(axis.text=element_text(size=18,family="Calibri"), axis.title=element_text(size=20,family="Calibri",lineheight=.5), plot.title=element_text(size=22,face="bold",family="Calibri"))

lasbiomass<-ggplot(plotdata, aes(x=vheight, y=lbiomasspresids))+ geom_smooth(method=glm , col=colors()[153],fill=colors()[142],se=TRUE) + 
  geom_point(alpha=0.35,size=1.5)+theme_classic() +labs(x="\\nBiomass\\n", y="\\nPartial Residuals\\n") + 
  theme(axis.text=element_text(size=18,family="Calibri"), axis.title=element_text(size=20,family="Calibri",lineheight=.5), plot.title=element_text(size=22,face="bold",family="Calibri"))

###### GRID TOGETHER

grid.arrange(richlpi,richvheight,richbiomass,
             cwmlpi, cwmvheight,cwmbloom,
             abundvheight,abundbiomass, nichevheight,
             lasvheight,lasbiomass,
             ncol=3, nrow=4)
