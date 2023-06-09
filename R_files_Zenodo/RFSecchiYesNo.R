setwd("H:/BigLake/Paper/Table2")

alldata <- read.csv("LagosDec2019.csv",na.strings=c("NA"),header=TRUE)
library(tidyverse)

QBasic <- alldata %>% select(summerpud,lake_area_ha,
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
                             lakes4ha_buffer500m_nlcd2011_pct_11,
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
                             lakes4ha_buffer500m_nlcd2011_pct_95,
                             boatlaunch,
                             beach,
                             hotels,
                             shelter,
                             toilets,
                             picnic,
                             bbq,
                             marina,
                             secchimeasure9514,
                             chlameasure9514,
                             no2no3measure9514,
                             tpmeasure9514,
                             maxdepthmeasure)

colnames(QBasic) <- c("pud","area","ppt","tmax","tmean","tmin","streamden",
                     "wetlanddis","wetlandarea","wetlandcount",
                     "fwetlanddis","fwetlandarea","fwetlandcount",
                     "swetlanddis","swetlandarea","swetlandcount",
                     "owetlanddis","owetlandarea","owetlandcount",
                     "lc11","lc21","lc22","lc23","lc24","lc31",
                     "lc41","lc42","lc43","lc52","lc71","lc81",
                     "lc82","lc90","lc95","boat","beach","hotel",
                     "shelter","toilets","picnic","bbq","marina",
                     "secchimeasure","chlameasure","no2no3measure","tpmeasure",
                     "maxdepthmeasure")

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

Z <- cbind(QBasic,QOtherw,P,POther,Control,Demo,lagy)
Z_complete <- Z[complete.cases(Z), ]

summersecchicount <- sum(Z_complete$secchimeasure)/length(Z_complete$secchimeasure)

library(randomForest)
library(caret)

train_index <- sample(1:nrow(Z_complete), 0.8 * nrow(Z_complete))
test_index <- setdiff(1:nrow(Z_complete), train_index)

train <- Z_complete[train_index,]
test <- Z_complete[test_index,]

secchi<-randomForest(formula=as.factor(secchimeasure) ~ .,
                data = train, ntree=500)

print(secchi)

test$secchipred <- predict(secchi,test)
test$secchipred <- as.factor(test$secchipred)
test$secchimeasure <- as.factor(test$secchimeasure)
confusionMatrix(test$secchipred,test$secchimeasure)

pred <- predict(object=secchi,newdata=test)
library(Metrics)
auc(actual=test$secchimeasure,predicted=pred)


theta <- summersecchicount
randomanalysis <- matrix(0,500,4)
N <- length(test_index)

for (j in 1:500)
{
flips <- rbinom(n = N, 
                size = 1, 
                prob = theta)

flips <-as.factor(flips)

cm <- confusionMatrix(flips,test$secchimeasure)

randomanalysis[j,] <- c(cm[["table"]][1],cm[["table"]][2],cm[["table"]][3],cm[["table"]][4]) 
}

sensitivity <- randomanalysis[,1]/(randomanalysis[,1] + randomanalysis[,2])
specificity <- randomanalysis[,4]/(randomanalysis[,3] + randomanalysis[,4])
accuracy <- (randomanalysis[,1]+randomanalysis[,4])/(randomanalysis[,1] + randomanalysis[,2] + randomanalysis[,3] + randomanalysis[,4])
mean(accuracy)
mean(sensitivity)
mean(specificity)


