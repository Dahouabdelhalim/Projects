## Modelling red fox home range size ##
# April 2020 #

# Load the appropriate packages
library(lme4)
library(MuMIn)
library(arm)

# Import the data and check the structure
foxdf<-read.csv("fox_data.csv")
str(foxdf) 

# Create a correlation matrix to identify any predictor variables that may be auto-correlated
cor(foxdf[,c("NDVI_mean","NDVI_seas","HFI","Temp_mean","Temp_seas","Precip_mean","Precip_seas")])
# NDVI_seas has a moderately high correlation with Temp_mean and Temp_seas, so we should exclude it from the model

# Fit a global model containing the predictor variables and interactions of interest
m0<-lmer(HR_mean ~ 
						NDVI_mean * HFI +
						Temp_mean + 
						Temp_seas + 
						Precip_mean + 
						Precip_seas + (1| Study_ID), 
						data = foxdf, REML=FALSE)					
plot(m0) # Check the residual plot to see if the model is a good fit for the data (residual plot doesn't look good)

# Refit the model with the response variable log-transformed
m1<-lmer(log_HR_mean ~ 
						NDVI_mean * HFI +
						Temp_mean + 
						Temp_seas + 
						Precip_mean + 
						Precip_seas + (1| Study_ID), 
						data = foxdf, REML=FALSE)
plot(m1) # Check the residual plot (it looks a lot better)
m1_std<-standardize(m1) # Standardise the predictor variables using the arm package

# Use the dredge function from MuMIn to fit all possible combinations of predictors from the global model
options(na.action = "na.fail") 
modsel_1<-dredge(m1_std); modsel_1   
# Note that the 'boundary (singular) fit' warning occurs because Study_ID accounts for zero variance in some of the models. This is not a problem for model fit or intepretation, though.

# Get the parameter estimates and R2 values for the top model. Change the '1' inside get.models() to 2, 3, 4, to get values for the second, third and fourth models (etc).
summary(get.models(modsel_1, 1)[[1]]);confint(get.models(modsel_1, 1)[[1]]);r.squaredGLMM(get.models(modsel_1, 1)[[1]])

###################################
# Analyse the HFI components separately 

# Create a correlation matrix to exclude correlated variables (r > 0.5) from featuring in the same model
smat <- abs(cor(foxdf[, c("Built_areas", "Cropland", "Lights", "Nav_waterways", "Pasture", "Human_popn", "Railways","Roads")])) <= 0.5
smat[!lower.tri(smat)] <- NA
z.smat<-smat
colnames(z.smat) <- paste("z.", colnames(z.smat),sep="")
rownames(z.smat) <- paste("z.", rownames(z.smat),sep="")

# Fit the model
m2<-lmer(log_HR_mean ~ 	
						Built_areas +
						Cropland +
						Lights + 
						Nav_waterways + 
						Pasture + 
						Human_popn + 
						Railways +
						Roads + (1| Study_ID), 
						data = foxdf, REML=FALSE)
# Note that the 'boundary (singular) fit' warning occurs because Study_ID accounts for zero variance in the global model (see summary(m2)). The same applies when the dredge function is run below. This is not a problem for model fit or intepretation, though.

m2_std<-standardize(m2) # Standardize the predictor variables

# Use the dredge function to fit all possible combinations of predictors from the global model, excluding models containing highly correlated variables (subset term)
modsel2<-dredge(m2_std, subset=z.smat) 
modsel2_subset<-modsel2[-18,] # Exclude the model that did not converge
modsel2_subset<- model.sel(modsel2_subset); modsel2_subset

# Get the parameter estimates for the top model
summary(get.models(modsel2_subset, 1)[[1]]);confint(get.models(modsel2_subset, 1)[[1]]);r.squaredGLMM(get.models(modsel2_subset, 1)[[1]])

# Compare models that only contain HFI or human population density
m_hfp<-standardize(lmer(log_HR_mean ~ 	HFI  + (1| Study_ID), 	data = foxdf, REML=FALSE))
m_hpd<-standardize(lmer(log_HR_mean ~ 	Human_popn  + (1| Study_ID), 	data = foxdf, REML=FALSE))
summary(m_hfp);AICc(m_hfp);r.squaredGLMM(m_hfp)
summary(m_hpd);AICc(m_hpd);r.squaredGLMM(m_hpd)
