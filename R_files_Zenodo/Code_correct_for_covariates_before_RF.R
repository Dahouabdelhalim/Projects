## Script to correct phenotypes and genotypes for potential confounding factors
## before conducting Random Forest analyses on CESRF Chinook

## For each phenotypic trait, I will correct both the phenotype and the genotypes using Linear Models using 
## all potential confounding covariates: Origin (i.e. hatchery line), sex, age, coordinate on Principal Component 1, and generation (i.e year)
## I'll essentially perform 6 corrections of the genotypes and phenotypes, one for each trait.

# Read in genotype data with missing values imputed, since RF can't take missing values.
# I'll use the "No Founders" data set since there is only one group of founders, and excluding it
# provides a most balanced experimental design in terms of equal representation of INT and SEG groups.

# Note: Files names and formats can be used in the R package GAPIT as well

myGD <- read.table("GD_no_founders.txt", header = TRUE)

# Read in phenotypes (including missing values). Potential covariates to include (Origin, Generation, sex, age, PC1) are in the first few columns
myY <- read.table("RF_phenotypes_no_founders.txt", header = TRUE)

# Combine genotypes and phenotypes
all_data_with_missing <- merge(myY,myGD,by="Individual")

###################################################################################################################################################################################################################################
###################################################################################################################################################################################################################################

# Now correct the phenotypes and genotypes for covariates, and use the residual from the model as the corrected phenotype and genotypes


# Age at maturity. Since age is a discrete trait, and correction will transform it into a more continuous trait,
# we will only correct the genotypes. Keeping the response variable discrete is supported by my simulation work.

all_data_age <- all_data_with_missing[which(all_data_with_missing$age!="NA"),c(1:6,12:9119)] # exclude individuals with missing phenotypes and omit other phenotypes
write.csv(all_data_age, file="RF_age_uncorrected.csv") # write the raw, uncorrected data to a file

age_corrected_all_covariates <- all_data_age # create another data frame to which the uncorrected phenotype and corrected genotypes (residuals) will be written
age_corrected_all_covariates[,7:9114] <- NA  # Keep the first few columns the same, but then replace with NA's over which you can write the residuals. I'm keeping the first few columns out of simplicity for the loop[i].

for (i in 7:ncol(all_data_age)){
  LM_SNP_i <- lm(all_data_age[,i] ~ factor(all_data_age$Origin) + factor(all_data_age$Generation) + factor(all_data_age$sex) + all_data_age$PC1, na.action=na.omit) 
  age_corrected_all_covariates[,i] <- LM_SNP_i$residuals  #This step gives the same result as Marine's method of using the predict function and subtracting it from the observed value
  colnames(age_corrected_all_covariates)[i]<-colnames(all_data_age)[i] 
  if(i%%50==0) print(i)
  # checked the output of the last SNP and YES! The new data set contains the residuals of the model, so it worked! 
}

# Verify that the residuals have been written to the data frame properly, using the last column as an example
age_corrected_all_covariates[,9114]-LM_SNP_i$residuals  #Should all be zero if correct

write.csv(age_corrected_all_covariates[,c(1,6:9114)],file="Age_corrected_all_covariates.csv",row.names=FALSE) #exclude covariate columns


###################################################################################################################################################################################################################################
###################################################################################################################################################################################################################################

# Daily Growth Coefficient. 

all_data_DGC <- all_data_with_missing[which(all_data_with_missing$Roza_CESRF_DGC!="NA"),c(1:6,11:9119)] # exclude individuals with missing phenotypes and omit other phenotypes
write.csv(all_data_DGC, file="RF_DGC_uncorrected.csv")

all_data_DGC2 <- all_data_DGC[which(all_data_DGC$age!="NA"),] #The model won't work if we have missing covariate values either
DGC_corrected_all_covariates <- all_data_DGC2 # create another data frame for the corrected phenotypes and genotypes (residuals)
DGC_corrected_all_covariates[,7:9115] <- NA  # Keep the first few columns the same, but then replace with NA's over which you can write the residuals. I'm keeping the first few columns out of simplicity for the loop[i].

for (i in 7:ncol(all_data_DGC2)){
  LM_SNP_i <- lm(all_data_DGC2[,i] ~ factor(all_data_DGC2$Origin) + factor(all_data_DGC2$Generation) + factor(all_data_DGC2$sex) + all_data_DGC2$age + all_data_DGC2$PC1, na.action=na.omit) 
  DGC_corrected_all_covariates[,i] <- LM_SNP_i$residuals  #This step gives the same result as Marine's method of using the predict function and subtracting it from the observed value
  colnames(DGC_corrected_all_covariates)[i]<-colnames(all_data_DGC2)[i] 
  if(i%%50==0) print(i)
  # checked the output of the last SNP and YES! The new data set contains the residuals of the model, so it worked! 
}

# Verify that the residuals have been written to the data frame properly, using the last column as an example
DGC_corrected_all_covariates[,9115]-LM_SNP_i$residuals  #Should all be zero if correct

write.csv(DGC_corrected_all_covariates[,c(1,7:9115)],file="DGC_corrected_all_covariates.csv",row.names=FALSE) #exclude covariate columns


###################################################################################################################################################################################################################################
###################################################################################################################################################################################################################################

# Forklength at Roza. 

all_data_forklength <- all_data_with_missing[which(all_data_with_missing$forklgth_Roza!="NA"),c(1:6,9,12:9119)] # exclude individuals with missing phenotypes and omit other phenotypes
write.csv(all_data_forklength, file="RF_forklength_uncorrected.csv")

all_data_forklength2 <- all_data_forklength[which(all_data_forklength$age!="NA"),] #Unfortunately the model won't work if we have missing covariate values either
forklength_corrected <- all_data_forklength2 # create another data frame for the corrected phenotypes and genotypes (residuals)
forklength_corrected[,7:9115] <- NA  # Keep the first few columns the same, but then replace with NA's over which you can write the residuals. I'm keeping the first few columns out of simplicity for the loop[i].

for (i in 7:ncol(all_data_forklength2)){
  LM_SNP_i <- lm(all_data_forklength2[,i] ~ factor(all_data_forklength2$Origin) + factor(all_data_forklength2$Generation) + factor(all_data_forklength2$sex) + all_data_forklength2$age + all_data_forklength2$PC1, na.action=na.omit) 
  forklength_corrected[,i] <- LM_SNP_i$residuals  #This step gives the same result as Marine's method of using the predict function and subtracting it from the observed value
  colnames(forklength_corrected)[i]<-colnames(all_data_forklength2)[i] 
  if(i%%50==0) print(i)
  # checked the output of the last SNP and YES! The new data set contains the residuals of the model, so it worked! 
}


# Verify that the residuals have been written to the data frame properly, using the last column as an example
forklength_corrected[,9115]-LM_SNP_i$residuals  #Should all be zero if correct

write.csv(forklength_corrected[,c(1,7:9115)],file="Forklength_corrected_all_covariates.csv",row.names=FALSE) #exclude covariate columns


###################################################################################################################################################################################################################################
###################################################################################################################################################################################################################################

# Weight at Roza. 

all_data_weight <- all_data_with_missing[which(all_data_with_missing$weight_Roza!="NA"),c(1:6,10,12:9119)] # exclude individuals with missing phenotypes and omit other phenotypes
all_data_weight2 <- all_data_weight[which(all_data_weight$age!="NA"),] #The model won't work if we have missing covariate values either

write.csv(all_data_weight2,file="RF_weight_uncorrected.csv",row.names=FALSE)

weight_corrected_all_covariates <- all_data_weight2 # create another data frame for the corrected phenotypes and genotypes (residuals)
weight_corrected_all_covariates[,7:9115] <- NA  # Keep the first few columns the same, but then replace with NA's over which you can write the residuals. I'm keeping the first few columns out of simplicity for the loop[i].

for (i in 7:ncol(all_data_weight2)){
  LM_SNP_i <- lm(all_data_weight2[,i] ~ factor(all_data_weight2$Origin) + factor(all_data_weight2$Generation) + factor(all_data_weight2$sex) + all_data_weight2$age + all_data_weight2$PC1, na.action=na.omit) 
  weight_corrected_all_covariates[,i] <- LM_SNP_i$residuals  #This step gives the same result as Marine's method of using the predict function and subtracting it from the observed value
  colnames(weight_corrected_all_covariates)[i]<-colnames(all_data_weight2)[i] 
  if(i%%50==0) print(i)
  # checked the output of the last SNP and YES! The new data set contains the residuals of the model, so it worked! 
}


# Verify that the residuals have been written to the data frame properly, using the last column as an example
weight_corrected_all_covariates[,9115]-LM_SNP_i$residuals #Should all be zero if correct

write.csv(weight_corrected_all_covariates[,c(1,7:9115)],file="weight_corrected_all_covariates.csv",row.names=FALSE) #exclude covariate columns


###################################################################################################################################################################################################################################
###################################################################################################################################################################################################################################

# Return time to Roza. 

all_data_return <- all_data_with_missing[which(all_data_with_missing$roza_day_of_year!="NA"),c(1:7,12:9119)] # exclude individuals with missing phenotypes and omit other phenotypes
all_data_return2 <- all_data_return[which(all_data_return$age!="NA"),] #The model won't work if we have missing covariate values either
write.csv(all_data_return2,file="RF_return_uncorrected.csv",row.names=FALSE)

return_corrected_all_covariates <- all_data_return2 # create another data frame for the corrected phenotypes and genotypes (residuals)
return_corrected_all_covariates[,7:9115] <- NA  # Keep the first few columns the same, but then replace with NA's over which you can write the residuals. I'm keeping the first few columns out of simplicity for the loop[i].

for (i in 7:ncol(all_data_return2)){
  LM_SNP_i <- lm(all_data_return2[,i] ~ factor(all_data_return2$Origin) + factor(all_data_return2$Generation) + factor(all_data_return2$sex) + all_data_return2$age + all_data_return2$PC1, na.action=na.omit) 
  return_corrected_all_covariates[,i] <- LM_SNP_i$residuals  #This step gives the same result as Marine's method of using the predict function and subtracting it from the observed value
  colnames(return_corrected_all_covariates)[i]<-colnames(all_data_return2)[i] 
  if(i%%50==0) print(i)
  # checked the output of the last SNP and YES! The new data set contains the residuals of the model, so it worked! 
}


# Verify that the residuals have been written to the data frame properly, using the last column as an example
return_corrected_all_covariates[,9115]-LM_SNP_i$residuals  #Should all be zero if correct

write.csv(return_corrected_all_covariates[,c(1,7:9115)],file="return_corrected_all_covariates.csv",row.names=FALSE) #exclude covariate columns


###################################################################################################################################################################################################################################
###################################################################################################################################################################################################################################

# Spawn timing at CESRF. 

all_data_spawn <- all_data_with_missing[which(all_data_with_missing$spawn_day_of_year!="NA"),c(1:6,8,12:9119)] # exclude individuals with missing phenotypes and omit other phenotypes
all_data_spawn2 <- all_data_spawn[which(all_data_spawn$age!="NA"),] #Unfortunately the model won't work if we have missing covariate values either

write.csv(all_data_spawn, file="RF_spawn_uncorrected.csv")

spawn_corrected_all_covariates <- all_data_spawn2 # create another data frame for the corrected phenotypes and genotypes (residuals)
spawn_corrected_all_covariates[,7:9115] <- NA  # Keep the first few columns the same, but then replace with NA's over which you can write the residuals. I'm keeping the first few columns out of simplicity for the loop[i].

for (i in 7:ncol(all_data_spawn2)){
  LM_SNP_i <- lm(all_data_spawn2[,i] ~ factor(all_data_spawn2$Origin) + factor(all_data_spawn2$Generation) + factor(all_data_spawn2$sex) + all_data_spawn2$age + all_data_spawn2$PC1, na.action=na.omit) 
  spawn_corrected_all_covariates[,i] <- LM_SNP_i$residuals  #This step gives the same result as Marine's method of using the predict function and subtracting it from the observed value
  colnames(spawn_corrected_all_covariates)[i]<-colnames(all_data_spawn2)[i] 
  if(i%%50==0) print(i)
  # checked the output of the last SNP and YES! The new data set contains the residuals of the model, so it worked! 
}


# Verify that the residuals have been written to the data frame properly, using the last column as an example
spawn_corrected_all_covariates[,9115]-LM_SNP_i$residuals  #Should all be zero if correct

write.csv(spawn_corrected_all_covariates[,c(1,7:9115)],file="spawn_corrected_all_covariates.csv",row.names=FALSE) #exclude covariate columns

