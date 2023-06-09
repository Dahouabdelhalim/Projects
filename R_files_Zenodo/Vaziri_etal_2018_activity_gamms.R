library(dplyr)
library(lubridate)
library(ggplot2)
library(tidyverse)
library(MuMIn)

setwd()

hsact<- read.csv("Vaziri_etal_2018_clean_data_ready_for_analysis.csv")
names(hsact)
head(hsact)


#make a column called date_time for the time an observation was made and make 
#that column and the column "start_time" into datetime objects
#round both into 30 minute intervals
hsact$date_time<-paste(hsact$date, hsact$time, sep = " ")

hsact$date_time<-mdy_hms(hsact$date_time, tz= "America/Chicago")
hsact$date_time<-round_date(hsact$date_time,unit = "minute")

hsact$start_time<-mdy_hm(hsact$start_time, tz= "America/Chicago")
hsact$start_time<-round_date(hsact$start_time, unit = "minute")


#cut this down so that I only have what I need
b <- hsact %>%
  select(bird_id_round, social_treat, lps_treat, power, date_time, start_time)
head(b)

##find the interval of time in hours between when a bird was treated and when a measurement was recorded

b$time_since_trt<-interval(b$start_time,b$date_time)
b$time_since_trt<-as.duration(b$time_since_trt)
dur_hr<-duration (minutes = 60)
dur_min<-duration (minutes = 1)


#this makes a column that has minutes and fraction of minutes since a bird was injected
b$min_since_inj<-b$time_since_trt/ dur_min
str(b)
#round that column to the nearest minute
b$min_since_inj<-round(b$min_since_inj)
head(b)

class(b$min_since_inj)

#find the number of hours in the interval
require(plyr)
b$half_hr_since_inj<-b$time_since_trt/ dur_hr
head(b)
#round that column to the nearest half hour
b$half_hr_since_inj<-round_any(b$half_hr_since_inj, 0.5, ceiling)

names(b)
unique(b$bird_id_round)
unique(b$min_since_inj)

#make a dataframe with average power per minute per bird
require(dplyr)
c<-dplyr::group_by(b, bird_id_round, min_since_inj, half_hr_since_inj)%>%
  dplyr::summarize(powmin=mean(power, na.rm=T))

head(c)


#show the number of minutes that have elapsed since the previous reading, by bird
c$elapsed_time<- ave(c$min_since_inj, c$bird_id_round,FUN = function(x) c(0, diff(x)))
head(c)
#make a column detnoting whether the number of minutes elapsed beween one minute and the next was one or not one
c$elapsed_min_1<-ifelse(c$elapsed_time==1, 1,NA)
head(c)
c[100:200,]
unique(c$elapsed_min_1)


#make a column that shows the difference in power from one minute to the next
c$power_diff<- ave(c$powmin, c$bird_id_round, FUN = function(x) c(0, diff(x)))
head(c)
dim(c)

#make a column that evaluates whether the previous observation occured one minute ago, and if so,
#whether the bird was active between a minute and the previous minute
c$active<-ifelse(is.na(c$elapsed_min_1 == T), NA,
                   ifelse(c$power_diff >= 4 |c$power_diff <= -4, 1,0))


#look at the count of minutes recorded for each hour each bird
test<-c(1,2,5,NA,NA,4)
length(test)
length(test[which(is.na(test)!=T)])
length.noNA<-function(x){length(x[which(is.na(x)!=T)])}
length.noNA(test)

obs_count<-tapply(c$active, list(c$bird_id_round,c$half_hr_since_inj), length.noNA )
dim(obs_count)
head(obs_count)
obs_count<-as.vector(obs_count)

#look at the sum of minutes active for each hour each bird
obs_active_sum<-tapply(c$active, list(c$bird_id_round,c$half_hr_since_inj), sum, na.rm=T)
head(obs_active_sum)
length(obs_active_sum)
dim(obs_active_sum)
# taking row names and repeating 48 times
bird_id_round<-rep(dimnames(obs_active_sum)[[1]],98)  
length(bird_id_round)
# taking col names repeating each one the same # of time as there were rows in powmin
half_hr_since_inj<-rep((dimnames(obs_active_sum)[[2]]),each=nrow(obs_active_sum))        
length(half_hr_since_inj)
#make obs_active_sum a vector
obs_active_sum<-as.vector(obs_active_sum)
head(obs_active_sum)
length(obs_active_sum)
#make a new dataframe with each bird and each time interval and the proportion of that time interval that the bird was active
obs<-na.omit(data.frame(bird_id_round, half_hr_since_inj, obs_active_sum, obs_count))
head(obs)
obs$propact<-obs$obs_active_sum/obs$obs_count
range(obs$propact)
head(obs)

#get rid of any rows that have fewer than 10 observations for each half hour
obs<-obs %>%
  filter(!obs_count < 11)

#get rid of columns I don't care about
obs<-obs %>%
  select(bird_id_round, half_hr_since_inj, propact)
obs

#match the treatment info and start time back to obs from c
obs$inj_time<-b$start_time[match(obs$bird_id_round,b$bird_id_round)]
#add in social treatment information for each bird
obs$soc<-hsact$social_treat[match(obs$bird_id_round,hsact$bird_id_round)]

#add in lps treatment information for each bird
obs$lps<-hsact$lps_treat[match(obs$bird_id_round,hsact$bird_id_round)]

names(obs)

#make column for combined treatment
obs$comb_trt<-as.factor(paste(obs$soc, obs$lps, sep = "_"))
str(obs)

ggplot(data=obs, aes(x=half_hr_since_inj, y = propact)) +
  geom_point()+
  geom_line()+
  facet_wrap(~comb_trt)

# getting transmitter on time from original data base

head(hsact)

hsact$on_time_decimal<-hour(hsact$start_time)+(minute(hsact$start_time)/60)
bybird<-group_by(hsact,bird_id_round)
birdstart <- dplyr::summarize(bybird, count = n(), on_time = mean(on_time_decimal, na.rm = T))
obs$on_time_decimal<-birdstart$on_time[match(obs$bird_id_round,birdstart$bird_id_round)]
str(obs)

#rename the treatments so Control is first alphabetically
levels(obs$comb_trt)
levels(obs$comb_trt)<-c("Treat-All", "Control-Mix", "Treat-Mix")

obs$comb_trt<-relevel(obs$comb_trt, ref = "Treat-All")

class(obs$half_hr_since_inj)
obs$half_hr_since_inj<-as.character(obs$half_hr_since_inj)
obs$half_hr_since_inj<-as.numeric(obs$half_hr_since_inj)

head(obs$half_hr_since_inj)
range(obs$half_hr_since_inj)

obs<-obs %>%
  filter(half_hr_since_inj < 8.5)
range(obs$half_hr_since_inj)

#write.csv(obs, file = "Vaziri_etal_2018_data_for_activity_figure.csv", row.names = F)

#do the gamms
hosp.act.gamm1<-gamm(propact~comb_trt+on_time_decimal+s(half_hr_since_inj, by = as.factor(comb_trt)),
                      random=list(bird_id_round=~1),data=obs)

hosp.act.gamm1.Exp<-gamm(propact~comb_trt+on_time_decimal+s(half_hr_since_inj, by = as.factor(comb_trt)),
                          random=list(bird_id_round=~1),data=obs,corr=corExp(form=~half_hr_since_inj|bird_id_round))
hosp.act.gamm1.Exp.nug<-gamm(propact~comb_trt+on_time_decimal+s(half_hr_since_inj, by = as.factor(comb_trt)),
                              random=list(bird_id_round=~1),data=obs,corr=corExp(form=~half_hr_since_inj|bird_id_round), nugget = T)
hosp.act.gamm1.Gaus<-gamm(propact~comb_trt+on_time_decimal+s(half_hr_since_inj, by = as.factor(comb_trt)),
                           random=list(bird_id_round=~1),data=obs,corr=corGaus(form=~half_hr_since_inj|bird_id_round))
hosp.act.gamm1.Gaus.nug<-gamm(propact~comb_trt+on_time_decimal+s(half_hr_since_inj, by = as.factor(comb_trt)),
                               random=list(bird_id_round=~1),data=obs,corr=corGaus(form=~half_hr_since_inj|bird_id_round), nugget = T)
hosp.act.gamm1.Lin<-gamm(propact~comb_trt+on_time_decimal+s(half_hr_since_inj, by = as.factor(comb_trt)),
                          random=list(bird_id_round=~1),data=obs,corr=corLin(form=~half_hr_since_inj|bird_id_round))
hosp.act.gamm1.Lin.nug<-gamm(propact~comb_trt+on_time_decimal+s(half_hr_since_inj, by = as.factor(comb_trt)),
                              random=list(bird_id_round=~1),data=obs,corr=corLin(form=~half_hr_since_inj|bird_id_round), nugget = T)
hosp.act.gamm1.Ratio<-gamm(propact~comb_trt+on_time_decimal+s(half_hr_since_inj, by = as.factor(comb_trt)),
                            random=list(bird_id_round=~1),data=obs,corr=corRatio(form=~half_hr_since_inj|bird_id_round))
hosp.act.gamm1.Ratio.nug<-gamm(propact~comb_trt+on_time_decimal+s(half_hr_since_inj, by = as.factor(comb_trt)),
                                random=list(bird_id_round=~1),data=obs,corr=corRatio(form=~half_hr_since_inj|bird_id_round), nugget = T)
hosp.act.gamm1.Spher<-gamm(propact~comb_trt+on_time_decimal+s(half_hr_since_inj, by = as.factor(comb_trt)),
                            random=list(bird_id_round=~1),data=obs,corr=corSpher(form=~half_hr_since_inj|bird_id_round))
hosp.act.gamm1.Spher.nug<-gamm(propact~comb_trt+on_time_decimal+s(half_hr_since_inj, by = as.factor(comb_trt)),
                                random=list(bird_id_round=~1),data=obs,corr=corSpher(form=~half_hr_since_inj|bird_id_round), nugget = T)

summary(hosp.act.gamm1$lme)

anova(hosp.act.gamm1$lme,hosp.act.gamm1.Exp$lme,hosp.act.gamm1.Exp.nug$lme,
      hosp.act.gamm1.Gaus$lme,hosp.act.gamm1.Gaus.nug$lme,
      hosp.act.gamm1.Lin$lme,hosp.act.gamm1.Lin.nug$lme,
      hosp.act.gamm1.Ratio$lme,hosp.act.gamm1.Ratio.nug$lme,
      hosp.act.gamm1.Spher$lme,hosp.act.gamm1.Spher.nug$lme)


AICc(hosp.act.gamm1,hosp.act.gamm1.Exp.nug,hosp.act.gamm1.Exp,
     hosp.act.gamm1.Gaus,hosp.act.gamm1.Gaus.nug,
     hosp.act.gamm1.Lin.nug,
     hosp.act.gamm1.Ratio,hosp.act.gamm1.Ratio.nug,
     hosp.act.gamm1.Spher,hosp.act.gamm1.Spher.nug)
#Exp best

#Exp seems to be best correlation structure
plot(Variogram (hosp.act.gamm1$lme,form=~half_hr_since_inj,maxDist=10,resType="n",robust=F))
plot(Variogram (hosp.act.gamm1.Exp.nug$lme,form=~half_hr_since_inj,maxDist=10,resType="n",robust=F))

summary(hosp.act.gamm1.Exp$lme)

# Fixed effects: y ~ X - 1 
# Value Std.Error  DF   t-value p-value
# X(Intercept)                                             0.7018935 0.5384411 283  1.303566  0.1934
# Xcomb_trtControl-Mix                                     0.2478185 0.0496436  19  4.991955  0.0001
# Xcomb_trtTreat-Mix                                      -0.0392375 0.0459797  19 -0.853367  0.4041
# Xon_time_decimal                                        -0.0288264 0.0411898  19 -0.699844  0.4925
# Xs(half_hr_since_inj):as.factor(comb_trt)Treat-AllFx1   -0.6922385 0.1488365 283 -4.650999  0.0000
# Xs(half_hr_since_inj):as.factor(comb_trt)Control-MixFx1 -0.0285301 0.0290709 283 -0.981397  0.3272
# Xs(half_hr_since_inj):as.factor(comb_trt)Treat-MixFx1   -0.3930189 0.1541834 283 -2.549036  0.0113


par(mfrow=c(2,2))
plot (hosp.act.gamm1.Exp$gam, select =c(1),main="A_T")
plot (hosp.act.gamm1.Exp$gam, select =c(2),main="M_C")
plot (hosp.act.gamm1.Exp$gam, select =c(3),main="M_T")

# diff smoothers for each social group
hosp.act.gamm2.Exp<-gamm(propact~comb_trt+on_time_decimal+s(half_hr_since_inj, by = as.factor(soc)),
                         random=list(bird_id_round=~1),data=obs,corr=corExp(form=~half_hr_since_inj|bird_id_round))

par(mfrow=c(2,2))
plot (hosp.act.gamm2.Exp$gam, select =c(1),main="All")
plot (hosp.act.gamm2.Exp$gam, select =c(2),main="Mix")


summary(hosp.act.gamm2.Exp$lme)
# Fixed effects: y ~ X - 1 
# Value Std.Error  DF   t-value p-value
# X(Intercept)                              0.7523166 0.5342318 284  1.408221  0.1602
# Xcomb_trtControl-Mix                      0.2427304 0.0494121  19  4.912373  0.0001
# Xcomb_trtTreat-Mix                       -0.0331837 0.0456070  19 -0.727601  0.4757
# Xon_time_decimal                         -0.0326904 0.0408686  19 -0.799890  0.4337
# Xs(half_hr_since_inj):as.factor(soc)AFx1 -0.6897812 0.1500446 284 -4.597175  0.0000
# Xs(half_hr_since_inj):as.factor(soc)MFx1 -0.3169239 0.1135262 284 -2.791637  0.0056

# diff smoothers for each LPS group:
hosp.act.gamm3.Exp<-gamm(propact~comb_trt+on_time_decimal+s(half_hr_since_inj, by = as.factor(lps)),
                         random=list(bird_id_round=~1),data=obs,corr=corExp(form=~half_hr_since_inj|bird_id_round))
par(mfrow=c(2,2))
plot (hosp.act.gamm3.Exp$gam, select =c(1),main="Control")
plot (hosp.act.gamm3.Exp$gam, select =c(2),main="LPS")

summary(hosp.act.gamm3.Exp$lme)
# Fixed effects: y ~ X - 1 
                                            # Value Std.Error  DF   t-value p-value
# X(Intercept)                              0.6945037 0.5363748 284  1.294810  0.1964
# Xcomb_trtControl-Mix                      0.2488300 0.0495030  19  5.026563  0.0001
# Xcomb_trtTreat-Mix                       -0.0366954 0.0456482  19 -0.803873  0.4314
# Xon_time_decimal                         -0.0283265 0.0410267  19 -0.690442  0.4983
# Xs(half_hr_since_inj):as.factor(lps)CFx1 -0.0287914 0.0291803 284 -0.986673  0.3246
# Xs(half_hr_since_inj):as.factor(lps)TFx1 -0.6603821 0.1294217 284 -5.102563  0.0000


#only 1 smoother:
hosp.act.gamm4.Exp<-gamm(propact~comb_trt+on_time_decimal+s(half_hr_since_inj),
                         random=list(bird_id_round=~1),data=obs,corr=corExp(form=~half_hr_since_inj|bird_id_round))

par(mfrow=c(1,2))
plot (hosp.act.gamm4.Exp$gam, select =c(1),main="All birds")
summary(hosp.act.gamm4.Exp$lme)
# Fixed effects: y ~ X - 1 
# Value Std.Error  DF   t-value p-value
# X(Intercept)              0.7112453 0.5346001 285  1.330425  0.1844
# Xcomb_trtControl-Mix      0.2417928 0.0493597  19  4.898587  0.0001
# Xcomb_trtTreat-Mix       -0.0373622 0.0455047  19 -0.821063  0.4218
# Xon_time_decimal         -0.0293313 0.0408911  19 -0.717302  0.4819
# Xs(half_hr_since_inj)Fx1 -0.5648036 0.1125644 285 -5.017604  0.0000

# no smoother for time since inj
hosp.act.gamm5.Exp<-gamm(propact~comb_trt+on_time_decimal,
                           random=list(bird_id_round=~1),data=obs,corr=corExp(form=~half_hr_since_inj|bird_id_round))

summary(hosp.act.gamm5.Exp$lme)
# Fixed effects: y ~ X - 1 
# Value Std.Error  DF   t-value p-value
# X(Intercept)          0.7566174 0.6376107 286  1.186645  0.2364
# Xcomb_trtControl-Mix  0.2262408 0.0590549  19  3.831025  0.0011
# Xcomb_trtTreat-Mix   -0.0405019 0.0542567  19 -0.746487  0.4645
# Xon_time_decimal     -0.0317040 0.0487693  19 -0.650082  0.5234

require(MuMIn)
AICc(hosp.act.gamm1.Exp,hosp.act.gamm2.Exp,hosp.act.gamm3.Exp,hosp.act.gamm4.Exp,hosp.act.gamm5.Exp)

# 
#                     df       AICc
# hosp.act.gamm1.Exp   13 -137.56472
# hosp.act.gamm2.Exp   11 -135.38538
# hosp.act.gamm3.Exp   11 -149.89822
# hosp.act.gamm4.Exp    9 -147.32093
# hosp.act.gamm5.Exp    7  -98.16312

summary(hosp.act.gamm3.Exp$lme)

# Fixed effects: y ~ X - 1 
#                                             Value   Std.Error  DF   t-value p-value
# X(Intercept)                              0.6945037 0.5363748 284  1.294810  0.1964
# Xcomb_trtControl-Mix                      0.2488300 0.0495030  19  5.026563  0.0001
# Xcomb_trtTreat-Mix                       -0.0366954 0.0456482  19 -0.803873  0.4314
# Xon_time_decimal                         -0.0283265 0.0410267  19 -0.690442  0.4983
# Xs(half_hr_since_inj):as.factor(lps)CFx1 -0.0287914 0.0291803 284 -0.986673  0.3246
# Xs(half_hr_since_inj):as.factor(lps)TFx1 -0.6603821 0.1294216 284 -5.102563  0.0000

anova(hosp.act.gamm3.Exp$lme)
      
      