####Condition script

##Using this script to look at mass/tarsus condition to see if it varies across years differently for different age classes


####Load Data####

library(here)

##Load condition data from database for 2015-2019 seasons
cond = read.csv(here::here("Input files","RBFW Condition scores 15_19.csv"))

##Load in age data from database
ages = read.csv(here::here("Input files","All RBFW ind histories for ages.csv"))


library(dplyr)
library(reshape2)
library(ggplot2)
library(cowplot)
library(lme4)
library(DHARMa)
library(lubridate)
library(viridis)



####Set up condition dataframe####

##Filter for males
ages.m = ages %>% filter(Sex=="M")
cond.m = cond %>% filter(Fwnumber %in% ages.m$FWno)

#Melt ages twice
ages.years = ages %>% select(contains("Year")) #get year data
ages.ages = ages %>% select(1:4,contains("Age")) #get age data
ages.years.m = melt(ages.years,id.vars = NULL) #melt both
ages.ages.m = melt(ages.ages,id.vars = 1:4) #melt both
ages.ages.m = ages.ages.m %>% select(-variable) 
ages.py = data.frame(ages.ages.m,ages.years.m$value) #combine
colnames(ages.py)[5:6] = c("Age","Year") #rename columns
ages.py = ages.py[!is.na(ages.py$Year),] #remove rows with no data
ages.py = ages.py %>% select(FWno,Age,Age.exact,Year) %>% distinct(FWno,Year,.keep_all=T)

##Add ages to condition data
cond.ma = merge(cond.m,ages.py,by.x=c("Fwnumber","Year"),by.y=c("FWno","Year"))
cond.ma = cond.ma %>% filter(!Age.exact=="") %>% filter(Age>0)

##Remove males at 1min or 2min, but keep in males with 3min and above - know they're old
cond.ma = cond.ma %>% filter(!(Age<3 & Age.exact=="min"))

##Add age categories column
cond.ma$Age.cat = NA
for (i in 1:nrow(cond.ma)) {if(cond.ma$Age[i]>=3) {cond.ma$Age.cat[i]="3+"} else {cond.ma$Age.cat[i]=cond.ma$Age[i]}} 

##Get capture dates in julian date form
cond.ma$Date = as.Date(cond.ma$Date,"%m/%d/%y")
cond.ma$jdate = yday(cond.ma$Date)


####Calculate mass-tarsus residuals####

##First see when captures are from
ggplot(data=cond.ma,aes(x=jdate)) + geom_histogram(binwidth = 5,color="black",fill="gray") + theme_cowplot()
ggplot(data=cond.ma,aes(x=jdate)) + geom_histogram(binwidth = 5,color="black",fill="gray") + theme_cowplot() + facet_wrap(facets="Year") + 
  geom_vline(aes(xintercept=200))

##Select for early captures - 200 is in mid-July - July 19th - and week 29 which is before most of the 2-year-olds began molting
#in the dry years. 
cond.mae = cond.ma %>% filter(jdate<=200)

##Remove individuals missing mass or tarsus
cond.mae = cond.mae %>% filter(Mass!="",Tarsus!="")

##Plot mass by tarsus
ggplot(data=cond.mae,aes(x=Tarsus,y=Mass)) + geom_point() + geom_smooth(method="lm") 
ggplot(data=cond.mae,aes(x=Tarsus,y=Mass)) + geom_point(aes(color=as.factor(Year)),size=2) + geom_smooth(method="lm") + 
  scale_color_viridis(discrete = T,name="Year") + theme_cowplot()
ggplot(data=cond.mae,aes(x=Tarsus,y=Mass)) + geom_point() + geom_smooth(method="lm") + facet_wrap(facets="Year")

##Calculate residuals
cond.mae$mtresid = resid(lm(Mass~Tarsus,data=cond.mae))

##Plot residuals by age class by year
ggplot(data=cond.mae,aes(x=Age.cat,y=mtresid)) + geom_hline(aes(yintercept=0),size=0.5,lty=2) + geom_boxplot() + 
  facet_wrap(facets="Year") + theme_cowplot()

##Group by wet and dry years and plot
cond.mae$WetorDry = NA
for(i in 1:nrow(cond.mae)) {if(cond.mae$Year[i] %in% c("2016","2018","2019")) {cond.mae$WetorDry[i]="Wet"} else 
  {cond.mae$WetorDry[i]="Dry"}}
#Plot
ggplot(data=cond.mae,aes(x=Age.cat,y=mtresid)) + geom_hline(aes(yintercept=0),size=0.5,lty=2) + geom_boxplot(width=0.5) + 
  facet_wrap(facets="WetorDry") + theme_cowplot() + geom_point(size=1.2,position = position_jitter(height=0,width=0.1))

##Bring in raincloud plots
setwd("/Users/Joe/Documents/Research/R Packages/RainCloudPlots-master/tutorial_R")
source("R_rainclouds.R")
ggplot(data=cond.mae,aes(x=Age.cat,y=mtresid)) + geom_hline(aes(yintercept=0),size=0.6,lty=1) + 
  geom_boxplot(width=0.4,position = position_nudge(x=0),aes()) + 
  facet_wrap(facets="WetorDry") + theme_cowplot() + 
  geom_point(size=1.2,position = position_jitter(width=0.08),alpha=0.6) +
  geom_flat_violin(position = position_nudge(x = .25, y = 0),adjust=0.8, trim = FALSE,fill="gray",linetype=0) + 
  xlab("Age") + ylab("Residual condition mass/tarsus")

##Rename Wet/Dry to Wetter/Drier
cond.mae$WetorDry2 = NA
for(i in 1:nrow(cond.mae)) {if(cond.mae$WetorDry[i]=="Wet") {cond.mae$WetorDry2[i]="Wetter"} else
  cond.mae$WetorDry2[i]="Drier"} 

##Plot to compare within age classes better
ggplot(data=cond.mae,aes(x=WetorDry2,y=mtresid)) + geom_hline(aes(yintercept=0),size=0.6,lty=1) + 
  geom_boxplot(width=0.4,position = position_nudge(x=0),aes()) + 
  facet_wrap(facets="Age.cat") + theme_cowplot() + 
  geom_point(size=1.2,position = position_jitter(width=0.08),alpha=0.6) +
  geom_flat_violin(position = position_nudge(x = .25, y = 0),adjust=0.8, trim = FALSE,fill="gray",linetype=0) + 
  xlab("Non-breeding season rainfall conditions") + ylab("Residual condition mass/tarsus")




####Test to see if differences among years are significant####

#Use fwnumber as random effect since individuals can be captured twice

library(lme4)
library(lmerTest)

##1-year-olds
cond.mae1 = cond.mae %>% filter(Age.cat=="1")
lm1 = lmer(mtresid~WetorDry + (1|Fwnumber),data=cond.mae1)
summary(lm1)
#LRT
lm1b = lmer(mtresid~1 + (1|Fwnumber),data=cond.mae1)
anova(lm1,lm1b)
#Paper: (X2=12.63, p<0.001)

##2-year-olds
cond.mae2 = cond.mae %>% filter(Age.cat=="2")
lm2 = lmer(mtresid~WetorDry + (1|Fwnumber),data=cond.mae2)
summary(lm2)
#LRT
lm2b = lmer(mtresid~1 + (1|Fwnumber),data=cond.mae2)
anova(lm2,lm2b)
#Paper: (X2=5.72, p=0.017)

##3-year-olds
cond.mae3 = cond.mae %>% filter(Age.cat=="3+")
lm3 = lmer(mtresid~WetorDry + (1|Fwnumber),data=cond.mae3)
summary(lm3)
#LRT
lm3b = lmer(mtresid~1 + (1|Fwnumber),data=cond.mae3)
anova(lm3,lm3b)
#Paper: (X2=5.78, p=0.016)

##Paper: All age classes were in better condition in wet years than in dry years





####Look at Fat and Keel as well####

##Fat
cond.maeF = cond.mae %>% filter(Fat<=3) #There's one error score of 5

ggplot(data=cond.maeF,aes(x=WetorDry2,y=Fat)) + 
  geom_boxplot(width=0.4,position = position_nudge(x=0),aes()) + 
  facet_wrap(facets="Age.cat") + theme_cowplot() + 
  geom_point(size=1.2,position = position_jitter(width=0.08),alpha=1) +
  geom_flat_violin(position = position_nudge(x = .25, y = 0),adjust=0.8, trim = FALSE,fill="gray",linetype=0) + 
  xlab("Non-breeding season rainfall conditions") + ylab("Fat score")


##1-year-olds
cond.maeF1 = cond.maeF %>% filter(Age.cat=="1")
lm4 = lmer(Fat~WetorDry + (1|Fwnumber),data=cond.maeF1)
summary(lm4)
#LRT
lm4b = lmer(Fat~1 + (1|Fwnumber),data=cond.maeF1)
anova(lm4,lm4b)
#Paper: (X2=2.91, p=0.088)

##2-year-olds
cond.maeF2 = cond.maeF %>% filter(Age.cat=="2")
lm5 = lmer(Fat~WetorDry + (1|Fwnumber),data=cond.maeF2)
summary(lm5)
#LRT
lm5b = lmer(Fat~1 + (1|Fwnumber),data=cond.maeF2)
anova(lm5,lm5b)
#Paper: (X2=0.38, p=0.537)

##3-year-olds
cond.maeF3 = cond.maeF %>% filter(Age.cat=="3+")
lm6 = lmer(Fat~WetorDry + (1|Fwnumber),data=cond.maeF3)
summary(lm6)
#LRT
lm6b = lmer(Fat~1 + (1|Fwnumber),data=cond.maeF3)
anova(lm6,lm6b)
#Paper: (X2=1.04, p=0.308)



##Keel
ggplot(data=cond.mae,aes(x=WetorDry2,y=Keel)) + 
  geom_boxplot(width=0.4,position = position_nudge(x=0),aes()) + 
  facet_wrap(facets="Age.cat") + theme_cowplot() + 
  geom_point(size=1.2,position = position_jitter(width=0.08),alpha=1) +
  geom_flat_violin(position = position_nudge(x = .25, y = 0),adjust=0.8, trim = FALSE,fill="gray",linetype=0) + 
  xlab("Non-breeding season rainfall conditions") + ylab("Muscle score")


##1-year-olds
lm7 = lmer(Keel~WetorDry + (1|Fwnumber),data=cond.mae1)
summary(lm7)
#LRT
lm7b = lmer(Keel~1 + (1|Fwnumber),data=cond.mae1)
anova(lm7,lm7b)
#Paper: (X2=9.11, p=0.003)

##2-year-olds
lm8 = lmer(Keel~WetorDry + (1|Fwnumber),data=cond.mae2)
summary(lm8)
#LRT
lm8b = lmer(Keel~1 + (1|Fwnumber),data=cond.mae2)
anova(lm8,lm8b)
#Paper: (X2=1.61, p=0.204)

##3-year-olds
lm9 = lmer(Keel~WetorDry + (1|Fwnumber),data=cond.mae3)
summary(lm9)
#LRT
lm9b = lmer(Keel~1 + (1|Fwnumber),data=cond.mae3)
anova(lm9,lm9b)
#Paper: (X2=1.99, p=0.158)

##Paper: 1-year-olds did differ in keel score across wet and dry years, but no other social classes did. 



##Write conditions scores to .csv 
#write.csv(cond.mae,here::here("Output files","Condition scores 15_19.csv"),row.names=F)







####Condition and molt date####


##Load molt date dataframe
pag2 = read.csv(here::here("Output files","twoyo_moltdates_dataframe.csv"))

##Get only males that were seen during molt
pag2.sdm = pag2 %>% filter(molt.time1>0)
pag2.sdm$Year = as.factor(pag2.sdm$Year)

##Merge with condition scores 
cond.mae2.short = cond.mae2 %>% select(Fwnumber,Year,mtresid,jdate)
pag2.sdmc = merge(pag2.sdm,cond.mae2.short,by.x=c("Fwnumber","Year"),by.y=c("Fwnumber","Year"))

##Get time to event column - days between condition measurement and molt date
pag2.sdmc$tte = pag2.sdmc$molt.date - pag2.sdmc$jdate

##Model
moltdatem = lmer(tte~mtresid + (1|Year) + (1|Fwnumber),data=pag2.sdmc)
summary(moltdatem)

#Check residuals
library(car)
hist(resid(moltdatem))
qqPlot(resid(moltdatem))
#Look good

##Likelihood ratio test for condition
moltdatemnoc = lmer(tte~1 + (1|Year) + (1|Fwnumber),data=pag2.sdmc)
anova(moltdatem,moltdatemnoc)

##Paper: Condition was an important predictor of when two-year-old males molted into ornamented plumage, with males in higher 
#condition molting sooner after the date their condition score was measured (X2=5.81,p=0.016). 



####Plot model prediction####


#Function to get confidence intervals from Ben Bolker
easyPredCI <- function(model,newdata=NULL,alpha=0.05) {
  ## baseline prediction, on the linear predictor (logit) scale:
  pred0 <- predict(model,re.form=NA,newdata=newdata)
  ## fixed-effects model matrix for new data
  X <- model.matrix(formula(model,fixed.only=TRUE)[-2],newdata)
  beta <- fixef(model) ## fixed-effects coefficients
  V <- vcov(model)     ## variance-covariance matrix of beta
  pred.se <- sqrt(diag(X %*% V %*% t(X))) ## std errors of predictions
  ## inverse-link function
  linkinv <- family(model)$linkinv
  ## construct 95% Normal CIs on the link scale and
  ##  transform back to the response (probability) scale:
  crit <- -qnorm(alpha/2)
  linkinv(cbind(conf.low=pred0-crit*pred.se,
                conf.high=pred0+crit*pred.se))
}


#Get prediction dataframe - need to have columns for each of the fixed effects here
cond.predict.df = data.frame(mtresid=seq(from=-0.80,to=0.90,by=0.01))

#Predict values 
cond.predict.p = predict(moltdatem,newdata=cond.predict.df,type="response",re.form=NA)

#Get confidence intervals
cond.predict.ci = easyPredCI(moltdatem,newdata=cond.predict.df)

#Combine prediction dataframes
cond.predict.df2 = cbind(cond.predict.df,tte = cond.predict.p,cond.predict.ci)

#Plot
ggplot(data=cond.predict.df2,aes(x=mtresid,y=tte)) + 
  geom_ribbon(aes(x=mtresid,y=tte,ymin=conf.low,ymax=conf.high),alpha=0.1) +
  geom_point(data=pag2.sdmc,aes(x=mtresid,y=tte,color=Year),size=2.5) +
  geom_line(color="black",size=1.2) + theme_cowplot() + xlab("Residual condition (mass/tarsus)") +
  ylab("Days remaining before molt into\\nornamented plumage") + scale_color_viridis(discrete = T)

ggplot(data=cond.predict.df2,aes(x=mtresid,y=tte)) + 
  geom_ribbon(aes(x=mtresid,y=tte,ymin=conf.low,ymax=conf.high),alpha=0.1) +
  geom_point(data=pag2.sdmc,aes(x=mtresid,y=tte),size=2.5) +
  geom_line(color="black",size=1.2) + theme_cowplot() + xlab("Residual condition (mass/tarsus)") +
  ylab("Days remaining before molt into\\nornamented plumage") 

#Plot with facet so size is same as other graph in Figure 3 - take out header in Affinity Designer
cond.predict.df2$age = 2

ggplot(data=cond.predict.df2,aes(x=mtresid,y=tte)) + 
  geom_ribbon(aes(x=mtresid,y=tte,ymin=conf.low,ymax=conf.high),alpha=0.1) +
  geom_point(data=pag2.sdmc,aes(x=mtresid,y=tte),size=2.5) +
  geom_line(color="black",size=1.2) + theme_cowplot() + xlab("Residual condition (mass/tarsus)") +
  ylab("Days remaining before molt into\\nornamented plumage") +
  facet_wrap(facets="age")



