library(mgcv)
#read in data
both<-rbind(read.csv("histCnem.csv",row.names=1),read.csv("modCnem.csv",row.names=1))
#set reference level of habitat factor to "wood"
both$habitat<-relevel(both$habitat,ref="wood")
#calculate counts for analyses
#First set total. 2 ways to do this. One must be commented out. If you need the other one, change the commenting out (#).
#This one uses the original totals as entered. A few are wrong. This should give analyses as in PLoSONE paper (and corrigendum).
both$tot<-both$origtot
#This one adds up the totals from the individual morph counts and ignores the entered totals. 
#Gives numerically slightly different results from paper
#both$tot<-both$P0+both$P3+both$P5+both$Y0+both$Y3+both$Y5+both$B0+both$B3+both$B5
#Now other counts
both$y<-both$Y0+both$Y3+both$Y5
both$ny<-both$tot-both$y
both$unb<-both$Y0+both$P0+both$B0
both$nunb<-both$tot-both$unb
both$mid<-both$Y3+both$P3+both$B3
both$nmid<-both$tot-both$mid
both$nmidbutbd<-both$nunb-both$mid
#set up subset with data for banded only
both.justbanded<-subset(both,nunb>0)
#main GAM fitting to produce more of the numbers in Table 2
GAM.yellow<-gam(cbind(y, ny) ~ s(DEM.m) + habitat + modhist + s(jan_precip) + 
    s(july_precip) + s(jan_min, july_max) + s(long, lat, by = modhist) + 
    year,family="quasibinomial",data=both)
summary(GAM.yellow)
GAM.midbband<-gam(cbind(mid, nmidbutbd) ~ s(DEM.m) + habitat + modhist + s(jan_precip) + 
    s(july_precip) + s(jan_min, july_max) + s(long, lat, by = modhist) + 
    year,family="quasibinomial",data=both.justbanded)
summary(GAM.midbband)
GAM.unb<-gam(cbind(unb,nunb) ~ s(DEM.m) + habitat + modhist + s(jan_precip) + 
    s(july_precip) + s(jan_min, july_max) + s(long, lat, by = modhist) + 
    year,family="quasibinomial",data=both)
summary(GAM.unb)
#do F tests for entire habitat factor
anova(GAM.yellow, update(GAM.yellow,~.-habitat),test="F")
anova(GAM.midbband, update(GAM.midbband,~.-habitat),test="F")
anova(GAM.unb, update(GAM.unb,~.-habitat),test="F")
#check year effects with temperature out of model
summary(update(GAM.yellow,~.-s(jan_min,july_max)))
summary(update(GAM.midbband,~.-s(jan_min,july_max)))
summary(update(GAM.unb,~.-s(jan_min,july_max)))
#check year effect on midbanded on historical only
summary(update(GAM.midbband,~.-modhist,data=subset(both.justbanded,modhist=="hist")))
#check year effect on midbanded with dune samples omitted
summary(update(GAM.midbband,data=subset(both.justbanded,habitat!="sand")))
#obtaining corrected (by model) morph frequencies and SEs for Fig 3
#Predictions at “typical” values which are for year 1970 and have mean values across whole 
#data set for climate variables, long, lat, altitude, and have modhist=hist. 
#They differ only in terms of habitat.
#First set up the data:
temp<-data.frame(DEM.m=mean(both$DEM.m),habitat=factor("wood"),
  modhist=factor("hist"),jan_precip=mean(both$jan_precip),
  july_precip=mean(both$july_precip),
  jan_min=mean(both$jan_min),
  july_max=mean(both$july_max),
  long=mean(both$long),lat=mean(both$lat),
  year=1970)
habitat.artif.data<-rbind(temp,temp,temp,temp)
habitat.artif.data$year<-c(1970,1970,1970,1970)
habitat.artif.data$habitat<-factor(c("wood","hedge","grass","sand"))
#Now calculate the numbers:
predict(GAM.yellow, habitat.artif.data,type="response",se.fit=T)
predict(GAM.midbband, habitat.artif.data,type="response",se.fit=T)
predict(GAM.unb, habitat.artif.data,type="response",se.fit=T)
#Correlation between predicted changes in midband and temp change 1950-2000
#set up data subset for modern data only
mod<-subset(both,modhist=="mod")
#take out sites where no banded
mod.justbanded<-subset(mod,nunb>0)
#take out sites where temperature diff is NA 
#(because that means no 1950 temperature data available)
mod.justbanded.fortemp<-na.exclude(mod.justbanded)
#get model predictions for modern midbanded
realpreds<-predict(GAM.midbband,type="response",newdata=mod.justbanded.fortemp)
#set up artificial old data for modern sites
artificial.hist.justbanded<-mod.justbanded.fortemp
artificial.hist.justbanded$modhist<-factor("hist")
artificial.hist.justbanded$year<-1950
#get model predictions for artificial old data
artpreds<-predict(GAM.midbband,type="response",newdata=artificial.hist.justbanded)
#differences
preddiffs<-realpreds-artpreds
tempdiff.edited<-unlist(mod.justbanded.fortemp$tempdiff)
cor.test(preddiffs,tempdiff.edited)
