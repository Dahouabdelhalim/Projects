library(lavaan)
library(MVN)

d_full <- as.data.frame(gone_wild)
d <- subset(d_full, Cup_Number != 11 & Cup_Number != 20)
subd <- d[,c(1,6,10,13,16,23,28)] #keeping wanted columns

#Standardize the data (putting everything into units of standard deviation)
#different variables have different units and you don't want that to affect the model 
d1 <- data.frame(lapply(2:7, function(x,...){(subd[,x]-mean(subd[,x],na.rm = T))/sd(subd[,x],na.rm = T)}))

#addID column back in 
d2<- cbind(subd[,1],d1)
#putting column names back in 
colnames(d2) <- colnames(subd[1:7])

#check for multivariate normality. 
mvn(d2, mvnTest = c("mardia"), covariance = T)


#### Original_MetaModel ####
m<-
  'Percent_Total_Consumed_Lesion ~ Lesion_Penetrometer + Lesion_Phenolics + Lesion_C_N
Lesion_Phenolics ~ Plant_Lesion_Cover
Lesion_Penetrometer ~ Plant_Lesion_Cover'

fit_1 <- sem(m, data=d2[,2:7], missing = "ml")
summary(fit_1, standardized=T, rsquare=TRUE)

#### Original_minus_C_N (attempt 1 at correcting model fit) ####
m2<-
  'Percent_Total_Consumed_Lesion ~ Lesion_Penetrometer + Lesion_Phenolics
Lesion_Phenolics ~ Plant_Lesion_Cover
Lesion_Penetrometer ~ Plant_Lesion_Cover
Lesion_C_N ~ Lesion_C_N'

fit_2 <- sem(m2, data=d2[,2:7], missing = "ml")
summary(fit_2, standardized=T, rsquare=TRUE)

#### Original_minus_C_N (attempt 2 at correcting model fit) ####
m2<-
  'Percent_Total_Consumed_Lesion ~ Lesion_Penetrometer + Lesion_Phenolics
Lesion_Phenolics ~ Plant_Lesion_Cover
Lesion_Penetrometer ~ Plant_Lesion_Cover'

fit_2 <- sem(m2, data=d2[,2:7], missing = "ml")
summary(fit_2, standardized=T, rsquare=TRUE)
AIC(fit_2)

#### Model_Fits_minus_Penetrometer (removing first insignificant path) ####
m2<-
  'Percent_Total_Consumed_Lesion ~ Lesion_Phenolics
Lesion_Phenolics ~ Plant_Lesion_Cover
Lesion_Penetrometer ~ Plant_Lesion_Cover'

fit_2 <- sem(m2, data=d2[,2:7], missing = "ml")
summary(fit_2, standardized=T, rsquare=TRUE)
AIC(fit_2)

#### Model_Fits_minus_Penetrometer (removing second insignificant path) ####
m2<-
  'Percent_Total_Consumed_Lesion ~ Lesion_Phenolics
Lesion_Phenolics ~ Plant_Lesion_Cover
Lesion_Penetrometer ~ Lesion_Penetrometer'

fit_2 <- sem(m2, data=d2[,2:7], missing = "ml")
summary(fit_2, standardized=T, rsquare=TRUE)
AIC(fit_2)

#### Just Phenolics #### 
#lesion cover --> phenolics --> consumption
mphenolics<-
  'Percent_Total_Consumed_Lesion ~ Lesion_Phenolics 
Lesion_Phenolics ~ Plant_Lesion_Cover'

fit_phenolics <- sem(mphenolics, data=d2[,2:7], missing = "ml")
summary(fit_phenolics, standardized=T, rsquare=TRUE)

#### Adding Clip Width ####
m2<-
  'Percent_Total_Consumed_Lesion ~ Lesion_Phenolics + Lesion_Clip_Width
Lesion_Phenolics ~ Plant_Lesion_Cover
Lesion_Penetrometer ~ Lesion_Penetrometer'

fit_2 <- sem(m2, data=d2[,2:7], missing = "ml")
summary(fit_2, standardized=T, rsquare=TRUE)
AIC(fit_2)


#### Checking AIC with width in the model matrix ####
m2<-
  'Percent_Total_Consumed_Lesion ~ Lesion_Phenolics
Lesion_Phenolics ~ Plant_Lesion_Cover
Lesion_Penetrometer ~ Lesion_Penetrometer
Lesion_Clip_Width ~ Lesion_Clip_Width'

fit_2 <- sem(m2, data=d2[,2:7], missing = "ml")
summary(fit_2, standardized=T, rsquare=TRUE)
AIC(fit_2)
