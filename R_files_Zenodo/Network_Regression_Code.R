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
library(effects)

#Load dataset
pan<- read.csv("~/PhD/Journal of Applied Ecology/Dryaddata/4. Networks Regression Analysis/Network Regression Data.csv")
names(pan)
str(pan)

#yvariables
general<-pan$generality
niche<-pan$niche

#collection variables
year<-as.factor(pan$year)
season<-as.factor(pan$season)
trt<-as.factor(pan$treat)
hood<-as.factor(pan$hood)

#explanatory
lpi.g<-pan$lpi.green #
ENN<-pan$enn.mn
gs.size<-pan$gs.area
gs.ca<-pan$gs.ca

bloom.ab<-as.integer(pan$bloom.ab)
bloom.ar<-pan$bloom.ar
vdiv<-pan$vdiv
vheight<-pan$vheight
biomass<-pan$biomass

#################MULTICOLINEARITY CHECK #######################################

my_data<-cbind(bloom.ab,bloom.ar,vdiv,vheight,biomass,ENN,gs.size,lpi.g,gs.ca,hood) 

print(my_data)
na.exclude(my_data)
res1<-rcorr(as.matrix(my_data))
res1
res1$r
res1$P
corrplot(res1$r, type="upper", order="hclust", p.mat = res1$P, sig.level = 0.01, insig = "blank")

#ENN correlates with gs,ca (total greenspace)

vif(lm(niche ~ ENN+lpi.g+gs.size+bloom.ab+bloom.ar+vdiv+vheight+biomass+gs.ca))
#gs.ca is 6.7, violates VIF < 3 thus take out

vif(lm(niche ~ ENN+lpi.g+gs.size+bloom.ab+bloom.ar+vdiv+vheight+biomass))
#bloom area violates VIF < 3 thus take out

vif(lm(niche ~ ENN+lpi.g+gs.size+bloom.ab+vdiv+vheight+biomass))
#all variables  VIF < 2 ok to proceed

my_data<-cbind(bloom.ab,vdiv,vheight,biomass,ENN,gs.size,lpi.g,hood) 

print(my_data)
na.exclude(my_data)
res1<-rcorr(as.matrix(my_data))
res1
res1$r
res1$P
corrplot(res1$r, type="upper", order="hclust", p.mat = res1$P, sig.level = 0.01, insig = "blank")

############# centering,scaling,creating data ###################################################################
newdata<-cbind(lpi.g,ENN,gs.size,bloom.ab,vdiv,vheight,biomass)
mydataedited<-as.data.frame(newdata)

mydataed<-scale(mydataedited, center=TRUE, scale=TRUE)
mydata.edited<-as.data.frame(mydataed)
str(mydata.edited)

lpi.g<-mydata.edited$lpi.g #
ENN<-mydata.edited$ENN
gs.size<-mydata.edited$gs.size
bloom.ab<-mydata.edited$bloom.ab
vdiv<-mydata.edited$vdiv
vheight<-mydata.edited$vheight
biomass<-mydata.edited$biomass

newdata2<-cbind(niche,general,lpi.g,ENN,gs.size,bloom.ab,vdiv,vheight,biomass,year,season,trt,hood)
mydata<-as.data.frame(newdata2)
str(mydata)

mydata<-transform(mydata,trt=as.factor(trt))
mydata<-transform(mydata,season=as.factor(season))
mydata<-transform(mydata,year=as.factor(year))
mydata<-transform(mydata,hood=as.factor(hood))
str(mydata)
print(mydata)
mydata<-na.exclude(mydata)
print(mydata)

######################### NORMALITY CHECK#######################################################################

#niche overlap
lmtest<-lmer(niche~ENN+lpi.g+gs.size+bloom.ab+vdiv+vheight+biomass+(1|hood),data=mydata)
summary(lmtest)
plot(lmtest) #normality looks good

#generality
lmtest<-lmer(general~ENN+lpi.g+gs.size+bloom.ab+vdiv+vheight+biomass+(1|hood),data=mydata)
summary(lmtest)
plot(lmtest) #normality looks good

####################  simple t -tests of treatments #######################################33
#niche overlap
treatment.n<-lmer(niche~trt*season+(1|hood),data=mydata)
summary(treatment.n) #trt tends towards significant
control.n<-lmer(niche~1+(1|hood),data=mydata)
summary(control.n)
anova(treatment.n,control.n)

#generality
treatment.g<-lmer(general~(1|hood)+trt*season,data=mydata)
summary(treatment.g) #trt IS significant
control.g<-lmer(general~1+(1|hood),data=mydata)
summary(control.g)
anova(treatment.g,control.g) #keep treatment*season in the model

##PLOTTING SEASON*TRT EFFECT
#one way to display without confidence interval bars
plot(predictorEffect("season", treatment.g),lines=list(multiline=TRUE))

#second way to display with confidence interval bars
e1.mm1 <- predictorEffect("season", treatment.g, ylevels = list(trt=c(1,5)))
e1.mm1
plot(e1.mm1, lines=list(multiline=TRUE), confint=list(style="bars"))


### LOCAL AND LANDSCAPE LMER

#generality
full.model<-lmer(general~ENN+lpi.g+gs.size+vdiv+vheight+bloom.ab+biomass+(1|hood), data = mydata) 
summary(full.model) #remove variable: biomass

modelb<-lmer(general~ENN+lpi.g+gs.size+vdiv+vheight+bloom.ab+(1|hood), data = mydata) 
anova(full.model,modelb) #no difference

#### STEP ONE
summary(modelb) #remove variable: ENN

modelc<-lmer(general~lpi.g+gs.size+vdiv+vheight+bloom.ab+(1|hood), data = mydata)

anova(modelb,modelc) #no difference

##### STEP 2
summary(modelc) #remove variable: lpi.g

modeld<-lmer(general~gs.size+vdiv+vheight+bloom.ab+(1|hood), data = mydata) 

anova(modelc,modeld) #no difference

### STEP 3
summary(modeld) #remove variable: vheight

modele<-lmer(general~gs.size+vdiv+bloom.ab+(1|hood), data = mydata) 

anova(modeld, modele) #no difference

## STEP 4
summary(modele) #remove variable: gs.size

modelf<-lmer(general~vdiv+bloom.ab+(1|hood), data = mydata) 

anova(modele,modelf)

## STEP 5
summary(modelf) #remove variable: vdiv

modelg<-lmer(general~bloom.ab+(1|hood), data = mydata)

anova(modelf,modelg) #no difference

### STEP 6
summary(modelg) #remove variable: bloom.ab

intercept<-lmer(general~1+(1|hood),data = mydata)

anova(intercept,modelg) #no difference

##best model for generality below
genmodel<-lmer(general~1+(1|hood),data = mydata)

############niche overlap model selection
full.model<-lmer(niche~ENN+lpi.g+gs.size+vdiv+vheight+bloom.ab+biomass+(1|hood), data = mydata) 
summary(full.model) #remove variable: vdiv

modelb<-lmer(niche~ENN+lpi.g+gs.size+vheight+bloom.ab+biomass+(1|hood), data = mydata) 
anova(full.model,modelb) #no difference

#### STEP ONE
summary(modelb) #remove variable: ENN

modelc<-lmer(niche~lpi.g+gs.size+vheight+bloom.ab+biomass+(1|hood), data = mydata) 
anova(modelb,modelc) #no difference

##### STEP 2
summary(modelc) #remove variable: biomass

modeld<-lmer(niche~lpi.g+gs.size+vheight+bloom.ab+(1|hood), data = mydata) #boundary is singular warning

anova(modelc,modeld) #no difference

### STEP 3
summary(modeld) #remove variable: LPI.g

modele<-lmer(niche~gs.size+bloom.ab+vheight+(1|hood), data = mydata) #boundary is singular warning

anova(modeld, modele) #no difference

## STEP 4
summary(modele) #remove variable: gs.size

modelf<-lmer(niche~bloom.ab+vheight+(1|hood), data = mydata) #boundary is singular warning

anova(modele,modelf) 

## STEP 5
summary(modelf) #remove variable: bloom.ab

modelg<-lmer(niche~vheight+(1|hood), data = mydata) ##boundary is singular warning
anova(modelf,modelg) #no difference

## STEP 7 
summary(modelg) #remove variable: vheight

intercept<-lmer(niche~1+(1|hood),data = mydata) ##boundary is singular warning
anova(modelg,intercept) #significant difference, keep model g

### BEST MODEL: niche overlap
nichefit<-lmer(niche~vheight+(1|hood), data = mydata)

#######Getting partial residuals for plotting
nichefit<-lmer(niche~vheight+(1|hood), data = mydata)
X <- model.matrix(~vheight,data=mydata)
beta <- fixef(nichefit)
## multiply each column of X by the corresponding beta 
beta_X <- sweep(X,MARGIN=2,STATS=beta,FUN="*")
## add residuals (columnwise addition works automatically, ## although we could use sweep with MARGIN=1) 
p_resid <- sweep(beta_X,MARGIN=1,STATS=residuals(nichefit),FUN="+")
partialresiduals<-as.data.frame(p_resid)
partialresiduals

nvheightpresids<-partialresiduals$vheight

########## PLOTTING NICHE, VEG HEIGHT EFFECT

#Reload dataset
pan<- read.csv("~/PhD/Journal of Applied Ecology/Dryaddata/4. Networks Regression Analysis/MalisaData.csv")

#y variable, rename/characterize
niche<-pan$niche

#explanatory, rename/characterize
vheight<-pan$vheight
hood<-as.factor(pan$hood)

#create subsetted dataset for plotting with original, not scaled features
plotdata<-cbind(niche,vheight,hood)
plotdata<-as.data.frame(plotdata)
plotdata<-transform(plotdata,hood=factor(hood))
str(plotdata) #check variables
# exclude NAS
plotdata <- na.exclude(plotdata)
na.action(plotdata)

## CREATE PLOT
nichevheight<-ggplot(plotdata, aes(x=vheight, y=nvheightpresids))+ geom_smooth(method=glm , col=colors()[153],fill=colors()[142],se=TRUE) + 
  +   geom_point(alpha=0.35,size=1.5)+theme_classic() +labs(x="\\nVegetation Height (cm)\\n", y="\\nPartial Residuals\\n") + 
  +   theme(axis.text=element_text(size=18,family="Calibri"), axis.title=element_text(size=20,family="Calibri",lineheight=.5), plot.title=element_text(size=22,face="bold",family="Calibri"))
nichevheight

