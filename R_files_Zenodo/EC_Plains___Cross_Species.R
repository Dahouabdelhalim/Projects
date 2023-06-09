##########################################################################################
##########################################################################################
##################### PLAINS ZEBRA EPIGENETIC CLOCKS AND VALIDATION ######################

#### Set environment ####

rm(list=ls())
#setwd('')
par(mfrow=c(1,1))

library(plyr)
library(glmnet)

#### Load sample sheet - select sample sets

sample_sheet <- read.csv("SampleSheetAgeN06final.csv")
sample_sheet <- sample_sheet[which(sample_sheet$CanBeUsedForAgingStudies=="yes"),] 
sample_sheet <- sample_sheet[which(sample_sheet$ConfidenceInAgeEstimate >= 90),] 
nrow(sample_sheet) #96 # all samples (biopsy + blood)

#transform age #optional 
sample_sheet$sqrtAge <- sqrt(sample_sheet$ageAtSamYrsR)

###################################### BLOOD CLOCK #######################################

#subset sample info
blood_samples <- sample_sheet[which(sample_sheet$Tissue=="Blood"),] #blood only
nrow(blood_samples) #76

#load normalized methylation data
load(file = "Zebra_probes_sesame_96BLOODsamples.Rdata")
betas <- t(normalized_betas_sesame) 

#subset blood methylation data to the 76 high confidence sampleshead(blood_samples)
blood_betas <- betas[which(row.names(betas)%in%blood_samples$Basename),]
count(row.names(blood_betas)==as.character(blood_samples$Basename)) #check if samples are in the same order

#initialize new dataframe to receive data from linear loo regression - use sqrt of age
#ages <- as.data.frame(blood_samples$ageAtSamYrsR) #use this for untransformed ages
ages <- as.data.frame(blood_samples$sqrtAge) #use this for sqrt transformed ages

colnames(ages)<- c("Age")
rownames(ages) <- rownames(blood_betas)
ages$sampleID <- as.character(blood_samples$SampleID)

ages$Prediction <- NA

#train the clock on plains zebra blood samples and test on same plains zebra blood samples. It leaves one sample out each iteration # a warning such as Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold is ok
for (i in 1:nrow(ages)){
  cvfit = cv.glmnet(blood_betas[-i,], ages$Age[-i], nfolds = 75, alpha= .5) 
  ages$Prediction[i] = predict(cvfit, newx = blood_betas[i,,drop=FALSE], type = "response", s = "lambda.min")
}

#results as # sqrt transformed ages  
cor.test(ages$Age,ages$Prediction)# 0.956, ci 0.931-0.971 

#plot predicted age against chronological age (retransform back to true age by squaring)
ages$Age2 <- ages$Age^2
plot(ages$Prediction^2 ~ ages$Age2, pch=16) 
abline(lm(Prediction^2 ~ Age2, data = ages), col='blue')
abline(0,1, lty='dashed')

#output the age predictions of the LOO model
ages$Prediction
write.csv(ages[c(1:3)], file="AgePrediction_LOO_BLOOD_sqrtTrans.csv")

#run the full cross-validated model (trained with all samples) and get model coefficients reported in Supplementary Data 1
#we only show this for the blood clock as it was used to predict age in other Equus species # you can adapt this code to get the final models for biopsy and combined

cvfit = cv.glmnet(blood_betas, blood_samples$sqrtAge, nfolds = 76, alpha= .5)   # a warning message about fold < 3 is normal

Final_coeficients <-as.matrix(coef(cvfit, s= "lambda.min")) # extract coefficients of the best model
Final_coeficients <-as.data.frame(Final_coeficients[which(Final_coeficients!=0),])
Final_coeficients$CG <- rownames(Final_coeficients)
Final_coeficients$lambda <- cvfit$lambda.min
blood_zebra_clock <- cvfit #this will be the clock we use to predict age in the other species

#output the coefficients of the full model - we call it the Final clock and use for evaluation of cross-species age prediction
write.csv(Final_coeficients, file="Coefficients_Final_Clock_Zebra_BLOOD_sqrtTransformed.csv") #EC_sqrtTrans_blood tab of Supplementary data 1
save(blood_zebra_clock, file="Final_Clock_Zebra_BLOOD_sqrtTransformed.Rdata")

##################################### BIOPSY CLOCK #######################################

biopsy_samples <- sample_sheet[which(sample_sheet$Tissue=="Biopsy"),] #biopsy only
nrow(biopsy_samples) #20

#load normalized methylation data
load(file = "Zebra_probes_sesame_24BIOPSYsamples.Rdata")
betas <- t(normalized_betas_sesame) 

#subset biopsy methylation data to the 20 high confidence samples
biopsy_betas <- betas[which(row.names(betas)%in%biopsy_samples$Basename),]
count(row.names(biopsy_betas)==as.character(biopsy_samples$Basename)) #check if samples are in the same order

#initialize new dataframe to receive data from linear loo regression - use sqrt of age
ages <- as.data.frame(biopsy_samples$sqrtAge) #use this for sqrt transformed ages

colnames(ages)<- c("Age")
rownames(ages) <- rownames(biopsy_betas)
ages$sampleID <- as.character(biopsy_samples$SampleID)

ages$Prediction <- NA

#train the clock on plains zebra biopsy samples and test on same plains zebra biopsy samples. Run the glm leaving one sample out each time # a warning such as Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold is ok
for (i in 1:nrow(ages)){
  cvfit = cv.glmnet(biopsy_betas[-i,], ages$Age[-i], nfolds = 19, alpha= .5) 
  ages$Prediction[i] = predict(cvfit, newx = biopsy_betas[i,,drop=FALSE], type = "response", s = "lambda.min")
}

#results as # sqrt transformed ages  
cor.test(ages$Age,ages$Prediction)## 0.615, ci 0.237-0.831 # 

#plot predicted age against chronological age (retransform back to true age by squaring)
ages$Age2 <- ages$Age^2
plot(ages$Prediction^2 ~ ages$Age2, pch=16, ylim=c(0,25))
abline(lm(Prediction^2 ~ Age2, data = ages), col='blue')
abline(0,1, lty='dashed')

#output the age predictions of the LOO model
write.csv(ages[c(1:3)], file="AgePrediction_LOO_20BIOPSYsamples_sqrtTransformation.csv")

#################################### COMBINED CLOCK ######################################

#load normalized methylation data
load(file = "Zebra_probes_IN-COMMON_sesame_120samplesMASK.Rdata")
betas <- t(betas_both)

#subset methylation data to high confidence samples
betas <- betas[which(row.names(betas)%in%sample_sheet$Basename),]
count(row.names(betas)==as.character(sample_sheet$Basename)) # check if samples are in the same order

#initialize new dataframe to receive data from linear loo regression - use sqrt of age
ages <-as.data.frame(sample_sheet$sqrtAge)

colnames(ages)<- c("Age")
rownames(ages) <- rownames(betas)
ages$sampleID <- as.character(sample_sheet$SampleID)
ages$Tissue <- as.factor(sample_sheet$Tissue)

#train the clock on plains zebra both tissue types samples and test on both tissue types. It leaves one sample out each iteration # a warning such as Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold is ok
ages$Prediction <- NA
for (i in 1:nrow(ages)){
  cvfit = cv.glmnet(betas[-i,], ages$Age[-i], nfolds = 95, alpha= .5) 
  ages$Prediction[i] = predict(cvfit, newx = betas[i,,drop=FALSE], type = "response", s = "lambda.min")
}

#results #  sqrt transformed ages 
cor.test(ages$Age,ages$Prediction) # 0.890, ci 0.840-0.926

#plot predicted age against chronological age (retransform back to true age by squaring)
ages$Age2 <- ages$Age^2
plot(ages$Prediction^2 ~ ages$Age2, pch=16, col=ages$Tissue, xlim=c(0,25), ylim=c(0,25))
abline(lm(Prediction^2 ~ Age2, data = ages), col='blue')
abline(0,1, lty='dashed')

#write out file with ages and predicted ages
write.csv(ages[c(1:3)], file="AgePrediction_LOO_COMBINED_sqrtTrans.csv")


##########################################################################################
############################## CROSS SPECIES AGE PREDICTION ##############################
############# Use final zebra blood clock to estimate age from other Equids ##############


######################################## Horses ##########################################

## reload blood_zebra_clock if need be
#load("Final_Clock_Zebra_BLOOD_sqrtTransformed.Rdata")

# load Horse samples and methylation data
horse_samples <- read.csv('/Users/ren/Dropbox/N23.2019-9032HorsesCarrieFinno/SampleSheetAgeN23final.csv')

horse_samples <- horse_samples[which(horse_samples$ConfidenceInAgeEstimate>90),]
horse_samples <- horse_samples[which(horse_samples$TissueDetailed=="Blood -Buffy Coat"),]
# remove outliers
horse_samples <- horse_samples[-c(which(horse_samples$Basename %in% c("203203210134_R05C01", "203203210107_R04C02"))),]
nrow(horse_samples) # 188

# sqrt transform horse ages
horse_samples$Age_sqrt <- sqrt(horse_samples$Age)

load(file = "/Users/ren/Dropbox/Zebras_Epi_collaboration/Horse_probes_sesame_190BLOODsamples.Rdata")
betas_Horse <- t(normalized_betas_sesame)
betas_Horse <- betas_Horse[which(row.names(betas_Horse)%in%horse_samples$Basename),]
count(row.names(betas_Horse)==as.character(horse_samples$Basename))
rm(normalized_betas_sesame)

#subset blood and horse betas to a set of common markers
ncol(betas_Horse)# 31836
ncol(blood_betas)# 31836
count(colnames(betas_Horse)==colnames(blood_betas))

horse_samples$Prediction <- NA
horse_samples$Prediction = predict(blood_zebra_clock, newx = betas_Horse, type = "response", s = "lambda.min")

cor.test(horse_samples$Age_sqrt,horse_samples$Prediction) #0.928
summary(lm(Prediction~Age_sqrt, data = horse_samples))

write.csv(horse_samples[,c("Basename","Age_sqrt", "Prediction")], file="AgePrediction_PlainsZebra_BLOOD_on_Horse.csv")

#### Plot Horse Predictions ####
plot(horse_samples$Prediction^2 ~ horse_samples$Age, pch=16)
abline(lm(Prediction^2~Age, data = horse_samples), col='blue')
abline(0,1, lty='dashed')


########################### Grevy's zebras and Somali Asses ##############################
### a little different because the small amount of data is provided in a different format

data <- read.csv("~/Dropbox/N06.2018-9323ZebrasBloodRenLarison/DataEquusFromN27/dat0SeSaMe.csv")

cgs <- data[,1]
data <- t(data[,-1])
colnames(data)<- cgs
count(colnames(data) %in% colnames(blood_betas))# 5718 cpgs are not in the zebra markers
data <- data[,which(colnames(data) %in% colnames(blood_betas))]
count(colnames(data) == colnames(blood_betas))# same order

data_info <- read.csv("~/Dropbox/N06.2018-9323ZebrasBloodRenLarison/DataEquusFromN27/datSample.csv")
nrow(data_info)# 11

Equus_grevyi_info <- data_info[which(data_info$SpeciesLatinName=="Equus grevyi"),]
Equus_grevyi <- data[which(rownames(data)%in% c(paste("X", Equus_grevyi_info$Basename, sep = ""))),]
count(rownames(Equus_grevyi)== paste("X",Equus_grevyi_info$Basename, sep = ""))

Equus_africanus_info <- data_info[which(data_info$SpeciesLatinName=="Equus africanus somaliensis"),]
Equus_africanus <- data[which(rownames(data)%in% c(paste("X", Equus_africanus_info$Basename, sep = ""))),]
count(rownames(Equus_africanus)== paste("X",Equus_africanus_info$Basename, sep = ""))


# Equus grevyi
Equus_grevyi_info$Age_sqrt <- sqrt(Equus_grevyi_info$Age)
Equus_grevyi_info$Prediction <- NA
Equus_grevyi_info$Prediction = predict(blood_zebra_clock, newx = Equus_grevyi, type = "response", s = "lambda.min")

cor.test(Equus_grevyi_info$Age_sqrt,Equus_grevyi_info$Prediction) #0.92
summary(lm(Prediction~Age_sqrt, data = Equus_grevyi_info))

write.csv(Equus_grevyi_info[,c("Basename","Age_sqrt", "Prediction")], file="AgePrediction_PlainsZebra_BLOOD_on_Equus_grevyi.csv")

#Plot
plot(Equus_grevyi_info$Prediction^2 ~ Equus_grevyi_info$Age, pch=16, xlim=c(0,19), ylim=c(0,19))
abline(lm(Prediction^2~Age, data = Equus_grevyi_info), col='blue')
abline(0,1, lty='dashed')

# Equus africanus
Equus_africanus_info$Age_sqrt <- sqrt(Equus_africanus_info$Age)
Equus_africanus_info$Prediction <- NA
Equus_africanus_info$Prediction = predict(blood_zebra_clock, newx = Equus_africanus, type = "response", s = "lambda.min")

cor.test(Equus_africanus_info$Age_sqrt,Equus_africanus_info$Prediction) #0.92
summary(lm(Prediction~Age_sqrt, data = Equus_africanus_info))

write.csv(Equus_africanus_info[,c("Basename","Age_sqrt", "Prediction")], file="AgePrediction_PlainsZebra_BLOOD_on_Equus_africanus.csv")

#Plot
plot(Equus_africanus_info$Prediction^2 ~ Equus_africanus_info$Age, pch=16, xlim=c(0,10), ylim=c(0,10))
abline(lm(Prediction^2~Age, data = Equus_africanus_info), col='blue')
abline(0,1, lty='dashed')




