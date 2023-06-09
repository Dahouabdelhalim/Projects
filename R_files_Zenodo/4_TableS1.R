# load packages and data
rm(list = ls())
packages <- c("dplyr","lubridate","xtable","tidyr","qdapTools","gghighlight",
              "hexbin","ggpubr","censReg","VGAM","stargazer","gridExtra","tree", "glmnet", "tidyverse", 
              "ISLR", "randomForest", "matrixStats")
lapply(packages, require, character.only = TRUE)
setwd("C:/Users/szabomorvai.agnes/Dropbox/Research/MASZK/Lasso")
data <- read.csv2("MASZK_dataCore.csv", header=T, sep=",", dec = ".")

IdVacc = data %>% select(c(Id, VaccineType)) 
IdVacc$NoVacc = as.factor(ifelse(IdVacc$VaccineType == "NoVaccine", 1, 0))
IdVacc$TypePf = as.factor(ifelse(IdVacc$VaccineType == "Pfizer", 1, 0))
IdVacc$TypeMo = as.factor(ifelse(IdVacc$VaccineType == "Moderna", 1, 0))
IdVacc$TypeAs = as.factor(ifelse(IdVacc$VaccineType == "Astrazeneca", 1, 0))
IdVacc$TypeSp = as.factor(ifelse(IdVacc$VaccineType == "Sputnik", 1, 0))
IdVacc$TypeSi = as.factor(ifelse(IdVacc$VaccineType == "Sinopharm", 1, 0))
#IdVacc$TypeJa = as.factor(ifelse(IdVacc$VaccineType == 6, 1, 0))

CharShort = data %>% select(c(Id,Age,IsFemale,CitySize,Schooling,Smoking,ChronicIll,AcuteIll,WorkLastWeek,WealthPreCovid,WealthNow))
Id <- data %>% select(c(Id))
PrefHet <- data %>% select(c(PfizerRating, ModernaRating, AstraRating, SputRating, SinoRating)) 

PrefHet <- data.matrix(PrefHet) 
PrefVar <- rowVars(PrefHet, na.rm = TRUE)
PrefMean <- rowMeans(PrefHet, na.rm = TRUE)
PrefMin <- rowMins(PrefHet, na.rm = TRUE)
PrefMax <- rowMaxs(PrefHet, na.rm = TRUE)
PrefHetMat <- cbind(Id, PrefVar, PrefMean, PrefMin, PrefMax)
PrefHetNew <- as.data.frame(PrefHetMat)
PrefHetNew$UnaccAny <- PrefMin
PrefHetNew[PrefHetNew$UnaccAny > 1,]$UnaccAny = 0

dataListen <- data %>% select(c(Id, VaccineDoctor, VaccineScientist, VaccineSceptical, VaccinePolitician, VaccineFamily, VaccineFriends, VaccineJournalist, VaccineCelebrity, SourceWebNews, SourceSocialNetwork, SourcePress, SourceRadio, SourceTV, SourceFamily, SourceFriends, SourceOther))

LassoVars = CharShort %>% 
  left_join(dataListen, by="Id") %>%    
  left_join(PrefHetNew, by="Id")
# eliminate missing variables
dim(LassoVars)
LassoVars = na.omit(LassoVars)
dim(LassoVars)
# Lasso
# outcome: PrefVar
x = model.matrix(PrefVar~ . -Id -PrefMean -PrefMax -PrefMin ,LassoVars )[,-1]
y = LassoVars$PrefVar

#outcome: PrefMean
#x = model.matrix(PrefMean  ~ . -Id -PrefVar -PrefMax -PrefMin ,LassoVars )[,-1]
#y = LassoVars$PrefMean

#outcome: PrefMin
#x = model.matrix(PrefMin  ~ . -Id -PrefVar -PrefMax -PrefMean ,LassoVars )[,-1]
#y = LassoVars$PrefMin

#outcome: PrefMax
#x = model.matrix(PrefMax  ~ . -Id -PrefVar -PrefMin -PrefMean ,LassoVars )[,-1]
#y = LassoVars$PrefMax

#outcome: UnaccAny
#x = model.matrix(UnaccAny  ~ . -Id -PrefVar -PrefMin -PrefMean -PrefMax ,LassoVars )[,-1]
#y = LassoVars$UnaccAny

grid =10^seq(10,-2, length =100)
set.seed (1)
train=sample (1: nrow(x), nrow(x)/2)

test=(- train )
y.test=y[test]

# Lasso
lasso.mod =glmnet (x[train ,],y[train],alpha =1, lambda =grid)
plot(lasso.mod)


set.seed (1)
cv.out =cv.glmnet (x[train ,],y[train],alpha =1)
plot(cv.out)
bestlam =cv.out$lambda.min
lasso.pred=predict (lasso.mod ,s=bestlam ,newx=x[test ,])
out=glmnet (x,y,alpha =1, lambda =grid)
lasso.coef=predict (out ,type = "coefficients",s=bestlam )[1:20 ,]
lasso.coef
lasso.coef[lasso.coef !=0]

