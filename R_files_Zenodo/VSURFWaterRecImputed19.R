setwd("~/BigLake/Paper/VSURFImputedDataH2OAct")
#setwd("H:/BigLake/Paper/VSURFImputedDataH2OAct")

alldata <- read.csv("LagosDec2019withImputed.csv",na.strings=c("NA"),header=TRUE)
library(tidyverse)
library(VSURF)

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
                             chlameasure9514,no2no3measure9514,tpmeasure9514)

colnames(QBasic) <- c("pud","area","ppt","tmax","tmean","tmin","streamden",
                      "wetlanddis","wetlandarea","wetlandcount",
                      "fwetlanddis","fwetlandarea","fwetlandcount",
                      "swetlanddis","swetlandarea","swetlandcount",
                      "owetlanddis","owetlandarea","owetlandcount",
                      "boatlaunch","beach","hotel",
                      "shelter","toilets","picnic","bbq","marina",
                      "chlameasure","no2no3measure","tpmeasure")

QWq <- alldata %>% select(X_19_avgsecchi9514,X_19_maxdepth)
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
Z_temp <- Z[which(Z$boatlaunch>0 | Z$beach>0 | Z$marina>0),]
Z_complete <- Z_temp[complete.cases(Z_temp), ]


set.seed(1)
train=sample(1:nrow(Z_complete),nrow(Z_complete)/2)
test=(-train)
Z_complete_train <- Z_complete[train, ]
Z_complete_test <- Z_complete[test, ]
#write.csv(Z_complete_train,"Z_complete_train_19.csv")
#write.csv(Z_complete_test,"Z_complete_test_19.csv")

v <- VSURF(pud ~ ., data = Z_complete_train, ntree=100, parallel = TRUE, ncores = 16)

threshvar <- as.data.frame(v[["varselect.thres"]])
interpvar <- as.data.frame(v[["varselect.interp"]])
predvar <- as.data.frame(v[["varselect.pred"]])

write.csv(threshvar,"threshvar_19.csv")
write.csv(interpvar,"interpvar_19.csv")
write.csv(predvar,"predvar_19.csv")

