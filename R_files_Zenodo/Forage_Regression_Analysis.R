#####load libraries
library(MASS)
library(blmeco)
library(glmmADMB)

####Load dataset
forage.0<- read.csv("~/PhD/Journal of Applied Ecology/Dryaddata/forage.csv")
str(forage.0)
forage = forage.0[-c(92,132,152),] #take out NAs which represent missing data from site vandalism of flowers
str(forage)

forage<- transform(forage,hood=factor(hood))
forage<- transform(forage,year=factor(year))
forage<- transform(forage,season=factor(season))
forage<- transform(forage,trt=factor(trt))

str(forage)

#########################OVERDISPERSION function#######################################################################

#created a dispersion check function that will work with most models
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval) }

########################## RED CLOVER ###########################
#Checking for overdisperson and comparing model distributions for full models
fullp.red<-glmmadmb(red.bee~as.numeric(red.area)+season+year+(1|hood), family="poisson", data=forage)
summary(fullp.red)
overdisp_fun(fullp.red) #very overdispersed, dispersion parameter = 600

fullnb.red<-glmmadmb(red.bee~as.numeric(red.area)+season+year+(1|hood), family="nbinom", data=forage)
overdisp_fun(fullnb.red) # dispersion parameter 0.84--> go with a negative binomial model
summary(fullnb.red)

#Compare a full negative binomial model vs. reduced

full.red<-glmmadmb(red.bee~as.numeric(red.area)+season+year+(1|hood), family="nbinom",data=forage)
summary(full.red)#red clover bloom area sig. driver of number of bees foraging on red clover

reduced.red<-glmmadmb(red.bee~as.numeric(red.area)+(1|hood),family="nbinom", data=forage)
summary(reduced.red)
overdisp_fun(reduced.red) #dispersion parameter = 1

anova(reduced.red,full.red) #significant difference between full and reduced models so include the full model

############################ QUEEN ANNES LACE #########################################
##Checking for overdisperson and comparing model distributions for full models
fullp.QAL<-glmmadmb(QAL.bee~as.numeric(QAL.area)+season+year+(1|hood), family="poisson", data=forage)
summary(fullp.QAL)
overdisp_fun(fullp.QAL) #18 dispersion parameter, overdispersed

fullnb.QAL<-glmmadmb(QAL.bee~as.numeric(QAL.area)+season+year+(1|hood), family="nbinom", data=forage)
overdisp_fun(fullnb.QAL) #0.8 dispersion parameter, neg. bin. better than poisson
summary(fullnb.QAL)

#Compare a full negative binomial model vs. reduced

full.QAL<-glmmadmb(QAL.bee~as.numeric(QAL.area)+season+year+(1|hood), family="nbinom", data=forage)
summary(full.QAL) #queen anne's lace sig. driver of number of bees foraging on QAL

reduced.QAL<-glmmadmb(QAL.bee~as.numeric(QAL.area)+(1|hood),family="nbinom", data=forage)
summary(reduced.QAL)
overdisp_fun(reduced.QAL) #dispersion parameter = 0.6

anova(reduced.QAL,full.QAL) #sig diff, go with FULL Model

################### CHICORY ######################################################
##Checking for overdisperson and comparing model distributions for full models
fullp.chicory<-glmmadmb(chicory.bee~as.numeric(chicory.area)+season+year+(1|hood), family="poisson", data=forage)
summary(fullp.chicory)
overdisp_fun(fullp.chicory) #dispersion parameter 640

fullnb.chicory<-glmmadmb(chicory.bee~as.numeric(chicory.area)+season+year+(1|hood), family="nbinom", data=forage)
overdisp_fun(fullnb.chicory) #dispersion parameter= 0.57, better than poisson model
summary(fullnb.chicory)

#Compare a full negative binomial model vs. reduced

full.chicory<-glmmadmb(chicory.bee~as.numeric(chicory.area)+season+year+(1|hood), family="nbinom", data=forage)
summary(full.chicory)

reduced.chicory<-glmmadmb(chicory.bee~as.numeric(chicory.area)+(1|hood),family="nbinom", data=forage)
summary(reduced.chicory) #chicory area not significant driver

anova(reduced.chicory,full.chicory) #no sig difference, use reduced model