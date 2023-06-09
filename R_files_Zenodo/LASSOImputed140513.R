setwd("H:/BigLake/Paper/LASSOImputedData0514")
alldata <- read.csv("LagosDec2019withImputed.csv",na.strings=c("NA"),header=TRUE)

library(tidyverse)
library(glmnet)
library(mpath)

QBasic <- alldata %>% select(summerpud,
                             lake_area_ha,
                             hu12_prism_ppt_30yr_normal_800mm,
                             hu12_prism_tmax_30yr_normal_800m,
                             hu12_prism_tmean_30yr_normal_800,
                             hu12_prism_tmin_30yr_normal_800m,
                             hu12_streamdensity_streams_densi,
                             hu12_wl_allwetlandsundissolved_a,
                             hu12_wl_allwetlandsundissolved_o,
                             hu12_wl_allwetlandsundissolved_c,
                             hu12_wl_forestedwetlandsundissol,v40,v41,
                             hu12_wl_scrubshrubwetlandsundiss,v43,v44,
                             hu12_wl_openwaterwetlandsundisso,v46,v47,
                             boatlaunch,beach,hotels,shelter,toilets,picnic,bbq,marina,
                             chlameasure0514,no2no3measure0514,tpmeasure0514)

colnames(QBasic) <- c("pud","area","ppt","tmax","tmean","tmin","streamden",
                      "wetlanddis","wetlandarea","wetlandcount",
                      "fwetlanddis","fwetlandarea","fwetlandcount",
                      "swetlanddis","swetlandarea","swetlandcount",
                      "owetlanddis","owetlandarea","owetlandcount",
                      "boatlaunch","beach","hotel",
                      "shelter","toilets","picnic","bbq","marina",
                      "chlameasure","no2no3measure","tpmeasure")

QWq <- alldata %>% select(X_14_avgsecchi0514,X_14_maxdepth)
colnames(QWq) <- c("secchi","maxdepth")                             

QLand <- alldata %>% select(lakes4ha_buffer500m_nlcd2011_pct,
                            v125,v126,v127,v128,
                            v129,
                            v130,v131,v132,
                            v133,v134,
                            v135,v136,
                            v137,v138)
QLandSum <- rowSums(QLand, na.rm = FALSE, dims = 1)
QLandSumCheck <- matrix(0,51107,1)
QLandSumCheck[QLandSum>95,]<-1
QLand[QLandSumCheck==0,] <- NA

Forest <- QLand$v130+QLand$v131+QLand$v132
Developed <- QLand$v125+QLand$v126+QLand$v127+QLand$v128
Agriculture <- QLand$v135+QLand$v136
Wetlands <- QLand$v137+QLand$v138
OtherLand <- QLand$v129+QLand$v133+QLand$v134

QLandFinal <- cbind(Forest,Developed,Agriculture,Wetlands,OtherLand)

P <- alldata %>% select(hu12popsqkm,distcbsameters)
colnames(P) <- c("popdenhu12","DistCBSA")

QOtherw <- alldata %>% select(hu12_lakes_lakes4ha_count,
                              hu12_lakes_lakes4ha_overlapping_,
                              avgdistnear5lakes,avgsizenear5lakes,avgdistnear5lakeswsecchimeasure,
                              lagsecchi,lagboatlaunch,lagbeach,laghotels,
                              lagshelter,lagtoilets,lagpicnic,lagbbq,
                              lagmarina)
colnames(QOtherw) <- c("hu12lakecount","hu12lakeareapct","distnear5lakes","sizenear5lakes","distlkwsecchi","lagsecchi","lagboat",
                       "lagbeach","laghotels","lagshelter","lagtoilets","lagpicnic",
                       "lagbbq","lagmarina")

POther <- alldata %>% select(avgcbsadistnear5lakes)

Control <- alldata %>% select(lowag,highag,nhd_lat,nhd_long) # otherag is the omitted region.

Demo <- alldata %>% select(hu12percbach,hu12percpoor,hu12percnonhispwhite,
                           hu12medianage,hu12avgmedhhinc)

colnames(Demo) <- c("PercBachelors","PercPoor","PerWhite","MedianAge","MedHHInc")

lagy <- alldata$lagy            

Z <- cbind(QBasic,QWq,QOtherw,QLandFinal,P,POther,Control,Demo,lagy)
Z_complete <- Z[complete.cases(Z), ]

set.seed(1)
train=sample(1:nrow(Z_complete),nrow(Z_complete)/2)
test=(-train)
Z_complete_train <- Z_complete[train, ]
Z_complete_test <- Z_complete[test, ]
#write.csv(Z_complete_train,"Z_complete_train_14.csv")
#write.csv(Z_complete_test,"Z_complete_test_14.csv")

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
                                         s=bestlam,newx=x[test,])[1:64,])
predbestlamyvalues <- exp(predict(lasso,family="poisson",s=bestlam,newx=x[test,]))
cvm <- mean((predbestlamyvalues-y[test])^2)

write.csv(predbestlamcoef,file="predbestlamcoefit_14.csv")


#Step 1

Z_complete_train_nosecchi <-  Z_complete_train[-c(31)]

x.s1 <- model.matrix(pud~.,Z_complete_train_nosecchi)[,-1]
y.s1 <- Z_complete_train_nosecchi$pud

lasso.s1 <- glmnet(x.s1,y.s1,family="poisson",alpha=1,maxit = 1000000,thresh=1e-06,lambda=grid)
print(lasso.s1)

lasso.s1.train.cv <- cv.glmnet(x.s1[train,],y.s1[train],alpha=1,family="poisson",maxit = 1000000,
                               thresh=1e-05,type.measure="mse")
plot(lasso.s1.train.cv)
print(lasso.s1.train.cv)

bestlam.s1=lasso.s1.train.cv$lambda.min

predbestlamcoef <- as.data.frame(predict(lasso.s1,family="poisson",type="coefficients",s=bestlam.s1,newx=x[test,])[1:63,])
predbestlamyvalues <- exp(predict(lasso.s1,family="poisson",s=bestlam,newx=x.s1[test,]))
cvm.s1 <- mean((predbestlamyvalues-y.s1[test])^2)

write.csv(predbestlamcoef,file="predbestlamcoefs1it_14.csv")


# Step 2

Z_complete_train_nopud <-  Z_complete_train[-c(1)]
x.s2 <- model.matrix(secchi~.,Z_complete_train_nopud)[,-1]
y.s2 <- Z_complete_train_nopud$secchi

lasso.s2 <- glmnet(x.s2,y.s2,family="poisson",alpha=1,maxit = 1000000,thresh=1e-06,lambda=grid)
print(lasso.s2)

lasso.s2.train.cv <- cv.glmnet(x.s2[train,],y.s2[train],alpha=1,family="poisson",maxit = 1000000,
                               thresh=1e-05,type.measure="mse")
plot(lasso.s2.train.cv)
print(lasso.s2.train.cv)

bestlam.s2=lasso.s2.train.cv$lambda.min


predbestlamcoef <- as.data.frame(predict(lasso.s2,family="poisson",type="coefficients",s=bestlam.s2,newx=x[test,])[1:63,])
predbestlamyvalues <- exp(predict(lasso.s2,family="poisson",s=bestlam,newx=x.s2[test,]))
cvm.s2 <- mean((predbestlamyvalues-y.s2[test])^2)

write.csv(predbestlamcoef,file="predbestlamcoefs2it_14.csv")

teststat <- c(cvm,cvm.s1,cvm.s2,bestlam,bestlam.s1,bestlam.s2)
write.csv(teststat,file="teststatit_14.csv")


