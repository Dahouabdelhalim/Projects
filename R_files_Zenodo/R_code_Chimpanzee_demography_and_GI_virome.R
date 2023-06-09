#R code for statistical analyses reported in the manuscript entitled ...

## "Demography, life history trade-offs, and the gastrointestinal virome of wild chimpanzees"

### Authors: Jacob D. Negrey, Melissa Emery Thompson, Kevin E. Langergraber, Zarin P. Machanda, John C. Mitani, Martin N. Muller, Emily Otali, Leah A. Owens, Richard W. Wrangham, and Tony L. Goldberg




#Open packages
library("dplyr")
library("tidyr")
library("car")
library("broom")
library("purrr")


#Import Datasets
virus <- read.csv("Data_Aging_Kibale_Chimpanzees_Viral_Load.csv")

virusPres <- read.csv("Data_Aging_Kibale_Chimpanzees_Viral_Presence.csv")






#Reformat to dataframes
virus<-tbl_df(virus)

virusPres<-tbl_df(virusPres)





#LINEAR MODELS FOR VIRAL RICHNESS

richness1<-lm(Richness ~ AgeJul16*sex+Community, data=virusPres)
summary(richness1)


shapiro.test(residuals(richness1))
qqPlot(residuals(richness1))





#Removing viruses most contributing to age*sex effect of richness

richnessPost1<-lm(Richness_No_Sali ~ AgeJul16*sex+Community, data=virusPres)
summary(richnessPost1)


richnessPost2<-lm(Richness_No_Sali_Chisa ~ AgeJul16*sex+Community, data=virusPres)
summary(richnessPost2)


richnessPost3<-lm(Richness_No_Sali_Chisa_Porpris1 ~ AgeJul16*sex+Community, data=virusPres)
summary(richnessPost3)





#Simulations to test richness results - Random selection of viruses and linear model x 1,000

Randomized = function(i) { 
  virusPres2 = virusPres[7:18] 
  Random = virusPres2[, sample(ncol(virusPres2), 3, replace=FALSE)] 
  virusPres$ToRemove = apply(Random,MARGIN=1,FUN=sum,na.rm=TRUE) 
  virusPres$TestRich = virusPres$Richness - virusPres$ToRemove 
  richnessRandom = lm(TestRich ~ AgeJul16*sex+Community, data=virusPres)
  richnessRandom
}

result<-lapply(1:1000,Randomized)

sapply(result, coef)

summaries <- lapply(result, summary)

SimulationCoefs <- lapply(summaries, function(x) x$coefficients[, c(1)])
SimulationPvalues <- lapply(summaries, function(x) x$coefficients[, c(4)])


library("xlsx")
write.xlsx(SimulationCoefs, "SimulationCoefs.xlsx") 
write.xlsx(SimulationPvalues, "SimulationPvalues.xlsx") 


#Descriptive stats were calculated in Excel. Please see supplementary data file for P values from simulations.







#LINEAR MODELS FOR TOTAL VIRAL LOAD


load1a<-lm(Total_Load ~ AgeJul16*sex+Community, data=virus)
summary(load1a)


shapiro.test(residuals(load1a))
qqPlot(residuals(load1a))


#Residuals are not normally distributed; transforming and running again


virus1<-filter(virus,Total_Load>0)


powerTransform(virus1$Total_Load,family="bcPower")
# 0.3584495



load1b<-lm((Total_Load^0.3584495)/0.3584495 ~ AgeJul16*sex+Community, data=virus1)
summary(load1b)

vif(load1b)

shapiro.test(residuals(load1b))
qqPlot(residuals(load1b))





#GENERALIZED LINEAR MODELS FOR VIRAL OCCURRENCE


#Astrovirus
glm1a<-glm(Astrovirus ~ AgeJul16*sex+Community, data=virusPres, family=binomial)

#Separation occurred; reducing model complexity
glm1b<-glm(Astrovirus ~ AgeJul16*sex, data=virusPres, family=binomial)
summary(glm1b)





#Bufavirus
glm2<-glm(Bufavirus ~ AgeJul16*sex+Community, data=virusPres, family=binomial)
summary(glm2)




#Circular ssDNA virus
glm3a<-glm(Circular_ssDNA_virus ~ AgeJul16*sex+Community, data=virusPres, family=binomial)

#Separation occurred; reducing model complexity

glm3b<-glm(Circular_ssDNA_virus ~ AgeJul16*sex, data=virusPres, family=binomial)
summary(glm3b)




#Chisavirus
glm4<-glm(Chisavirus ~ AgeJul16*sex+Community, data=virusPres, family=binomial)
summary(glm4)




#Picobirna-Like Virus
glm5a<-glm(Picobirna_like_virus ~ AgeJul16*sex+Community, data=virusPres, family=binomial)

#Separation occurred; reducing model complexity
glm5b<-glm(Picobirna_like_virus ~ AgeJul16*sex, data=virusPres, family=binomial)
summary(glm5b)




#Porprismacovirus 1
glm6<-glm(Porprismacovirus_1 ~ AgeJul16*sex+Community, data=virusPres, family=binomial)
summary(glm6)




#Porprismacovirus 2
glm7<-glm(Porprismacovirus_2 ~ AgeJul16*sex+Community, data=virusPres, family=binomial)
summary(glm7)




#Porprismacovirus 3
glm8<-glm(Porprismacovirus_3 ~ AgeJul16*sex+Community, data=virusPres, family=binomial)
summary(glm8)




#Porprismacovirus 4
glm9<-glm(Porprismacovirus_4 ~ AgeJul16*sex+Community, data=virusPres, family=binomial)
summary(glm9)




#Porprismacovirus 5
glm10<-glm(Porprismacovirus_5 ~ AgeJul16*sex+Community, data=virusPres, family=binomial)
summary(glm10)




#Porprismacovirus 6
glm11<-glm(Porprismacovirus_6 ~ AgeJul16*sex+Community, data=virusPres, family=binomial)
summary(glm11)




#Salivirus
glm12a<-glm(Salivirus ~ AgeJul16*sex+Community, data=virusPres, family=binomial)
summary(glm12a)

#Separation occurred; reducing model complexity

glm12b<-glm(Salivirus ~ AgeJul16*sex, data=virusPres, family=binomial)

#Separation occurred; reducing model complexity

glm12c<-glm(Salivirus ~ AgeJul16+sex, data=virusPres, family=binomial)
summary(glm12c)






#LINEAR MODELS FOR LOAD PER VIRUS


#Bufavirus
Bufav<-dplyr::select(virus,AgeJul16,sex,Community,Bufavirus)
Bufav<-filter(Bufav,Bufavirus>0)

lm1<-lm(Bufavirus ~ AgeJul16*sex+Community, data=Bufav)
summary(lm1)

qqPlot(residuals(lm1))
shapiro.test(residuals(lm1))




#Porprismacovirus 1
Porpris1<-dplyr::select(virus,AgeJul16,sex,Community,Porprismacovirus_1)
Porpris1<-filter(Porpris1,Porprismacovirus_1>0)

lm2<-lm(Porprismacovirus_1 ~ AgeJul16*sex+Community, data=Porpris1)
summary(lm2)

qqPlot(residuals(lm2))
shapiro.test(residuals(lm2))




#Porprismacovirus 2
Porpris2<-dplyr::select(virus,AgeJul16,sex,Community,Porprismacovirus_2)
Porpris2<-filter(Porpris2,Porprismacovirus_2>0)

lm3<-lm(Porprismacovirus_2 ~ AgeJul16*sex+Community, data=Porpris2)
summary(lm3)

qqPlot(residuals(lm3))
shapiro.test(residuals(lm3))




#Porprismacovirus 4
Porpris4<-dplyr::select(virus,AgeJul16,sex,Community,Porprismacovirus_4)
Porpris4<-filter(Porpris4,Porprismacovirus_4>0)

lm4<-lm(Porprismacovirus_4 ~ AgeJul16*sex+Community, data=Porpris4)
summary(lm4)

qqPlot(residuals(lm4))
shapiro.test(residuals(lm4))




#Porprismacovirus 5
Porpris5<-dplyr::select(virus,AgeJul16,sex,Community,Porprismacovirus_5)
Porpris5<-filter(Porpris5,Porprismacovirus_5>0)

lm5<-lm(Porprismacovirus_5 ~ AgeJul16*sex+Community, data=Porpris5)
summary(lm5)


qqPlot(residuals(lm5))
shapiro.test(residuals(lm5))




#Porprismacovirus 6
Porpris6<-dplyr::select(virus,AgeJul16,sex,Community,Porprismacovirus_6)
Porpris6<-filter(Porpris6,Porprismacovirus_6>0)

lm6<-lm(Porprismacovirus_6 ~ AgeJul16*sex+Community, data=Porpris6)
summary(lm6)


qqPlot(residuals(lm6))
shapiro.test(residuals(lm6))







#NOTE: P values were grouped by predictor and subsequently adjusted using the Benjamini-Hochberg method in the function p.adjust()

## The following form was used:

### p.adjust(PredictorX,method="BH")