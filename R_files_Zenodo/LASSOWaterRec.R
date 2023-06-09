setwd("H:/BigLake/Paper/LASSOLimitedDataH2OAct")
alldata <- read.csv("LagosDec2019.csv",na.strings=c("NA"),header=TRUE)

library(tidyverse)
library(glmnet)
library(mpath)

QBasic <- alldata %>% select(summerpud,
                             lake_area_ha,
                             hu12_prism_ppt_30yr_normal_800mm2_annual_mean,
                             hu12_prism_tmax_30yr_normal_800mm2_annual_mean,
                             hu12_prism_tmean_30yr_normal_800mm2_annual_mean,
                             hu12_prism_tmin_30yr_normal_800mm2_annual_mean,
                             hu12_streamdensity_streams_density_mperha,
                             hu12_wl_allwetlandsundissolved_avgsize_ha,
                             hu12_wl_allwetlandsundissolved_overlapping_area_pct,
                             hu12_wl_allwetlandsundissolved_count,
                             hu12_wl_forestedwetlandsundissolved_avgsize_ha,
                             hu12_wl_forestedwetlandsundissolved_overlapping_area_pct,
                             hu12_wl_forestedwetlandsundissolved_count,
                             hu12_wl_scrubshrubwetlandsundissolved_avgsize_ha,
                             hu12_wl_scrubshrubwetlandsundissolved_overlapping_area_pct,
                             hu12_wl_scrubshrubwetlandsundissolved_count,
                             hu12_wl_openwaterwetlandsundissolved_avgsize_ha,
                             hu12_wl_openwaterwetlandsundissolved_overlapping_area_pct,
                             hu12_wl_openwaterwetlandsundissolved_count,
                             boatlaunch,
                             beach,
                             hotels,
                             shelter,
                             toilets,
                             picnic,
                             bbq,
                             marina,
                             chlameasure9514,
                             no2no3measure9514,
                             tpmeasure9514)

QWq <- alldata %>% select(maxdepth,avgsecchi9514,nsecchigreater19514)

QLand <- alldata %>% select(lakes4ha_buffer500m_nlcd2011_pct_11,
                            lakes4ha_buffer500m_nlcd2011_pct_21,
                            lakes4ha_buffer500m_nlcd2011_pct_22,
                            lakes4ha_buffer500m_nlcd2011_pct_23,
                            lakes4ha_buffer500m_nlcd2011_pct_24,
                            lakes4ha_buffer500m_nlcd2011_pct_31,
                            lakes4ha_buffer500m_nlcd2011_pct_41,
                            lakes4ha_buffer500m_nlcd2011_pct_42,
                            lakes4ha_buffer500m_nlcd2011_pct_43,
                            lakes4ha_buffer500m_nlcd2011_pct_52,
                            lakes4ha_buffer500m_nlcd2011_pct_71,
                            lakes4ha_buffer500m_nlcd2011_pct_81,
                            lakes4ha_buffer500m_nlcd2011_pct_82,
                            lakes4ha_buffer500m_nlcd2011_pct_90,
                            lakes4ha_buffer500m_nlcd2011_pct_95)

QLandSum <- rowSums(QLand, na.rm = FALSE, dims = 1)
QLandSumCheck <- matrix(0,51107,1)
QLandSumCheck[QLandSum>95,]<-1
QLand[QLandSumCheck==0,] <- NA


Forest <- QLand$lakes4ha_buffer500m_nlcd2011_pct_41+
  QLand$lakes4ha_buffer500m_nlcd2011_pct_42+
  QLand$lakes4ha_buffer500m_nlcd2011_pct_43
Developed <- QLand$lakes4ha_buffer500m_nlcd2011_pct_21+
  QLand$lakes4ha_buffer500m_nlcd2011_pct_22+
  QLand$lakes4ha_buffer500m_nlcd2011_pct_23+
  QLand$lakes4ha_buffer500m_nlcd2011_pct_24
Agriculture <- QLand$lakes4ha_buffer500m_nlcd2011_pct_81+
  QLand$lakes4ha_buffer500m_nlcd2011_pct_82
Wetlands <- QLand$lakes4ha_buffer500m_nlcd2011_pct_90+
  QLand$lakes4ha_buffer500m_nlcd2011_pct_95
OtherLand <- QLand$lakes4ha_buffer500m_nlcd2011_pct_31+
  QLand$lakes4ha_buffer500m_nlcd2011_pct_52+
  QLand$lakes4ha_buffer500m_nlcd2011_pct_71

QFull <- cbind(QBasic,QWq,Forest,Agriculture,Developed,Wetlands,OtherLand)

colnames(QFull) <- c("pud","area","ppt","tmax","tmean","tmin","streamden",
                     "wetlanddis","wetlandarea","wetlandcount",
                     "fwetlanddis","fwetlandarea","fwetlandcount",
                     "swetlanddis","swetlandarea","swetlandcount",
                     "owetlanddis","owetlandarea","owetlandcount",
                     "boat","beach","hotel","shelter","toilets","picnic","bbq","marina",
                     "chlameasure","no2no3measure","tpmeasure",
                     "maxdepth","avgsecchi","secching1",
                     "Forest","Agriculture","Developed","Wetlands","OtherLand")

P <- alldata %>% select(HU12PopSqKM,DistCBSAmeters)
colnames(P) <- c("popdenhu12","DistCBSA")

QOtherw <- alldata %>% select(hu12_lakes_lakes4ha_count,
                              hu12_lakes_lakes4ha_overlapping_area_pct,
                              avgdistnear5lakes,avgsizenear5lakes,avgdistnear5lakeswsecchimeasure,
                              lagsecchi,lagboatlaunch,lagbeach,laghotels,
                              lagshelter,lagtoilets,lagpicnic,lagbbq,
                              lagmarina)
colnames(QOtherw) <- c("hu12lakecount","hu12lakeareapct","distnear5lakes","sizenear5lakes","distlkwsecchi","lagsecchi","lagboat",
                       "lagbeach","laghotels","lagshelter","lagtoilets","lagpicnic",
                       "lagbbq","lagmarina")

POther <- alldata %>% select(avgcbsadistnear5lakes)

Control <- alldata %>% select(lowag,highag,nhd_lat,nhd_long) # otherag is the omitted region.

Demo <- alldata %>% select(HU12PercBach,HU12PercPoor,HU12PercNonHispWhite,
                           HU12medianAge,HU12avgMedHHinc)

colnames(Demo) <- c("PercBachelors","PercPoor","PerWhite","MedianAge","MedHHInc")

lagy <- alldata$lagy

Z <- cbind(QFull,QOtherw,P,POther,Control,Demo,lagy)
Z_temp <- Z[which(Z$boat>0 | Z$beach>0 | Z$marina>0),]
Z_complete <- Z_temp[complete.cases(Z_temp), ]


set.seed(1)
train=sample(1:nrow(Z_complete),nrow(Z_complete)/2)
test=(-train)
Z_complete_train <- Z_complete[train, ]
Z_complete_test <- Z_complete[test, ]
write.csv(Z_complete_train,"Z_complete_train.csv")
write.csv(Z_complete_test,"Z_complete_test.csv")

#### LASSO #####

grid=10^seq(10,-2,length=300)

x <- model.matrix(pud~.,Z_complete_train)[,-1]
y <- Z_complete_train$pud

lasso <- glmnet(x,y,family="poisson",alpha=1,maxit = 1000000,thresh=1e-06,lambda=grid)
print(lasso)

train=sample(1:nrow(x),nrow(x)/2)

lasso.train.cv <- cv.glmnet(x[train,],y[train],alpha=1,family="poisson",maxit = 1000000,
                            thresh=1e-05,type.measure="mse")
plot(lasso.train.cv)
print(lasso.train.cv)

bestlam=lasso.train.cv$lambda.min

test=(-train)

predbestlamcoef <- as.data.frame(predict(lasso,family="poisson",type="coefficients",
                                         s=bestlam,newx=x[test,])[1:65,])
predbestlamyvalues <- exp(predict(lasso,family="poisson",s=bestlam,newx=x[test,]))
cvm <- mean((predbestlamyvalues-y[test])^2)

write.csv(predbestlamcoef,file="predbestlamcoefit.csv")


#Step 1

Z_complete_train_nosecchi <-  Z_complete_train[-c(32)]

x.s1 <- model.matrix(pud~.,Z_complete_train_nosecchi)[,-1]
y.s1 <- Z_complete_train_nosecchi$pud

lasso.s1 <- glmnet(x.s1,y.s1,family="poisson",alpha=1,maxit = 1000000,thresh=1e-06,lambda=grid)
print(lasso.s1)

lasso.s1.train.cv <- cv.glmnet(x.s1[train,],y.s1[train],alpha=1,family="poisson",maxit = 1000000,
                               thresh=1e-05,type.measure="mse")
plot(lasso.s1.train.cv)
print(lasso.s1.train.cv)

bestlam.s1=lasso.s1.train.cv$lambda.min

predbestlamcoef <- as.data.frame(predict(lasso.s1,family="poisson",type="coefficients",s=bestlam.s1,newx=x[test,])[1:64,])
predbestlamyvalues <- exp(predict(lasso.s1,family="poisson",s=bestlam,newx=x.s1[test,]))
cvm.s1 <- mean((predbestlamyvalues-y.s1[test])^2)

write.csv(predbestlamcoef,file="predbestlamcoefs1it.csv")


# Step 2

Z_complete_train_nopud <-  Z_complete_train[-c(1)]
x.s2 <- model.matrix(avgsecchi~.,Z_complete_train_nopud)[,-1]
y.s2 <- Z_complete_train_nopud$avgsecchi

lasso.s2 <- glmnet(x.s2,y.s2,family="poisson",alpha=1,maxit = 1000000,thresh=1e-06,lambda=grid)
print(lasso.s2)

lasso.s2.train.cv <- cv.glmnet(x.s2[train,],y.s2[train],alpha=1,family="poisson",maxit = 1000000,
                               thresh=1e-05,type.measure="mse")
plot(lasso.s2.train.cv)
print(lasso.s2.train.cv)

bestlam.s2=lasso.s2.train.cv$lambda.min


predbestlamcoef <- as.data.frame(predict(lasso.s2,family="poisson",type="coefficients",s=bestlam.s2,newx=x[test,])[1:64,])
predbestlamyvalues <- exp(predict(lasso.s2,family="poisson",s=bestlam,newx=x.s2[test,]))
cvm.s2 <- mean((predbestlamyvalues-y.s2[test])^2)

write.csv(predbestlamcoef,file="predbestlamcoefs2it.csv")

teststat <- c(cvm,cvm.s1,cvm.s2,bestlam,bestlam.s1,bestlam.s2)
write.csv(teststat,file="teststatit.csv")
