### Read Data #############################################################
data <- read.csv("CT_data.csv")
rptdata <- read.csv("Repeatability.csv")

### Packages ##############################################################
library(pbkrtest) # Enables Kenward-Roger F-tests
library(lmerTest) # Enables Kenward-Roger F-tests
library(FactoMineR) #PCA package
library(lme4) # Linear models
library(car) # Linear models
library(Hmisc) #Correlations
library(OutlierDetection) # Outlier testing
library(rptR) # Repeatability analysis
library(scales)# Re-scaling of response variable
library(predictmeans)
library(EMAtools)#Cohens d effect size
library(simr)
library(effectsize)
library(nlme)
library(MuMIn)
library(RColorBrewer)

### Repeatability #########################################################

UMe <- rpt(U.mean ~ (1|ID), grname = "ID", data= rptdata,
           datatype= "Gaussian", nboot=10000, npermut=10000)
UMa <- rpt(U.max ~ (1|ID), grname = "ID", data= rptdata,
           datatype= "Gaussian", nboot=10000, npermut=10000)
UMi <- rpt(U.min ~ (1|ID), grname = "ID", data= rptdata,
           datatype= "Gaussian", nboot=10000, npermut=10000)

LMe <- rpt(L.mean ~ (1|ID), grname = "ID", data= rptdata,
           datatype= "Gaussian", nboot=10000, npermut=10000)
LMa <- rpt(L.max ~ (1|ID), grname = "ID", data= rptdata,
           datatype= "Gaussian", nboot=10000, npermut=10000)
LMi <- rpt(L.min ~ (1|ID), grname = "ID", data= rptdata,
           datatype= "Gaussian", nboot=10000, npermut=10000)

print(UMe) #R=0.635, CI(0.147,0.856), p=0.00105
print(UMa) #R=0.604, CI(0.11,0.838), p=0.00199
print(UMi) #R=0, CI(0,0.4), p=0.5
print(LMe) #R=0.769, CI(0.368,0.914), p=2.54e-05
print(LMa) #R=0.741, CI(0.33,0.903), p=6.66e-05
print(LMi) #R=0, CI(0,0.406), p=1

# As minimum values for the upper and lower tract are non-repeatable, we will not use them in subsequent analyses.

#### Correlation Tests ####################################################
names(data)
vars=c("U.mean","U.max","U.min","L.mean","L.max","L.min")
Corr=data[vars]
mat=as.matrix(Corr)
mode(mat) = "numeric"
rcorr(mat,type=c("pearson"))

#### PCA and Data Compilation #############################################
names(data)
pc.subset <- data[,c(10:11,13:14)]
pc.comp <- PCA(pc.subset, scale.unit=TRUE, graph = TRUE)
data.pca <- cbind(data[,c(1:6)],pc.comp$ind$coord[,1:3]) # Adds loading values for the first three PC's

pc.comp$eig # PC eigenvalues
pc.comp$var # PC trait loadings

data.pca$Bias <- as.factor(data.pca$Bias)
data.pca$Conflict <- as.factor(data.pca$Conflict)
data.pca$Voltage <- as.factor(data.pca$Voltage)
data.pca$Population <- as.factor(data.pca$Population)

### Re-scale Dim.1 ###
#Data fails to converge with negative eigenvalue
data.pca$Dim.1 <- rescale(data.pca$Dim.1, to = c(0, 1)) 
data.pca$Dim.2 <- rescale(data.pca$Dim.2, to = c(0, 1)) 


#### Checking for Outliers################################################

names(data.pca)
outlier.data <- data.pca[, c(5, 11)]
outlier.data
dens(outlier.data, k=3, boottimes=50000, rnames=TRUE)


#### PC1 Model Reduction #################################################
pc1.model1 <- lmer(Dim.1 ~ Conflict+Bias+Weight+(1 + Weight * Conflict|Population) + (1|Voltage)
                   + Conflict * Bias
                   + Conflict * Weight
                   + Bias * Weight, data=data.pca, REML=TRUE)

pc1.model2 <- lmer(Dim.1 ~ Conflict+Bias+Weight+(1 + Weight * Conflict|Population) + (1|Voltage)
                  + Conflict * Bias
                  + Conflict * Weight, data=data.pca, REML=TRUE)

KR.A1 <- KRmodcomp(pc1.model1, pc1.model2)
summary(KR.A1) #Non-sig (remove)

pc1.model3 <- lmer(Dim.1 ~ Conflict+Bias+Weight+(1 + Weight * Conflict|Population) + (1|Voltage)
                   + Conflict * Bias, data=data.pca, REML=TRUE)

KR.A2 <- KRmodcomp(pc1.model2, pc1.model3)
summary(KR.A2) #Non-sig (remove)

pc1.model4 <- lmer(Dim.1 ~ Conflict+Bias+Weight+(1 + Weight * Conflict|Population) + (1|Voltage), data=data.pca, REML=TRUE)

KR.A3 <- KRmodcomp(pc1.model3, pc1.model4)# Non-sig (remove)
summary(KR.A3)

#Final model results#
anova(pc1.model4, type = c("III","II","I"), ddf="Kenward-Roger")
summary(pc1.model4, ddf="Kenward-Roger")


pc1.nlme <- lme(Dim.1 ~ Conflict + Bias + Weight, random = ~1 + Conflict|Population, data=data.pca, method = "REML")
anova(pc1.nlme)



# Means, CI's and Effect Size ##############################################
emmeans(pc1.model4, ~ Conflict) # Estimated means
emmeans(pc1.model4, ~ Bias)

Standardized.model <- standardize(pc1.model4)
summary(Standardized.model, type = c("III","II","I"), ddf="Kenward-Roger")

standardize_parameters(pc1.model4, method = "refit", ci=0.95) # Coefficients match those of above, but 95% CI are provided. 


  
# Standardization method: refit

#            | Coefficient (std.) |        95% CI
#------------------------------------------------
#(Intercept) |               0.65 | [-0.17, 1.46]
#ConflictL   |              -0.43 | [-0.92, 0.05]
#BiasMale    |              -0.28 | [-0.76, 0.21]
#Weight      |               0.49 | [ 0.20, 0.78]

