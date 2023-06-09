#setwd("~/BigLake/Paper/VSURFLimitedData")
setwd("H:/BigLake/Paper/VSURFLimitedDataH2OAct")

alldata <- read.csv("LagosDec2019.csv",na.strings=c("NA"),header=TRUE)
library(tidyverse)
library(VSURF)

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

v <- VSURF(pud ~ ., data = Z_complete_train, ntree=100, parallel = TRUE, ncores = 16)

threshvar <- as.data.frame(v[["varselect.thres"]])
interpvar <- as.data.frame(v[["varselect.interp"]])
predvar <- as.data.frame(v[["varselect.pred"]])

write.csv(threshvar,"threshvarlimiteddatawithoutSDSecchi.csv")
write.csv(interpvar,"interpvarlimiteddatawithoutSDSecchi.csv")
write.csv(predvar,"predvarlimiteddatawithoutSDSecchi.csv")

