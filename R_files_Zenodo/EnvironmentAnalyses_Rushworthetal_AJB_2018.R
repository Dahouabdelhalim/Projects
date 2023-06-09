### R code for analyses from Rushworth et al. American Journal of Botany 
### "Ecological differentiation facilitates fine-scale coexistence of sexual and asexual Boechera"
### Please contact C. Rushworth with questions (catherine.rushworth@gmail.com)


library(DescTools)
library(adegenet)
library(plyr)
library(lsmeans)
library(lme4)
library(fields)
library(ape)

setwd()



#######################################################################################################
##### G-test for goodness of fit, pure asexual vs. pure sexual vs. mixed sexual/asexual (num pops) ####
#######################################################################################################


observed = c(55,58,36) # number of mixed pops, pure asex pops, pure sex pops (in that order)
expected = c(0.3333333333,0.3333333333,0.3333333333) # expected proportions
GTest(x=observed, p=expected, correct="none")
 


#######################################################################################################
##################################             DAPC           #########################################
#######################################################################################################


### Load in data for DAPC and MANCOVA
fulldata <- read.csv("FullEnvironmentalData_Rushworthetal_AJB_2018.csv", header=T) # FC1, FC2, FC4
FC3data <- read.csv("FC3FullEnvironmentalData_Rushworthetal_AJB_2018.csv", header=T) 
# FC3data contains only the four focal taxa

### Work with "fulldata" to get it ready for FC4
mixdata <- na.omit(fulldata) 
### removing NAs limits the dataset to pure asexual and mixed sexual/asexual populations
mixdata$AMS <- as.factor(mixdata$AMS)


### FC1: Asexual/sexual comparison

pop <- fulldata[,2] # population code numbers
AMS <- fulldata[,3] # mixed asexual/sexual vs. pure asexual pops: NA=pure sexual, 1=mixed, 2=pure asexual
RSP <- fulldata[,4] # reproductive system + ploidy; 1=diploid sexual, 2=diploid asexual, 3=triploid asexual
RS <- fulldata[,5] # reproductive system; asexual or sexual, binary
vars1 <- fulldata[,6:29] # raw environmental data

dapc.full <- dapc(vars1, grp=RS, n.pca=NULL, center=TRUE, scale=TRUE) # note that this is now standardized
scatter(dapc.full, scree.da=FALSE, bg = "white", pch = 17:18, col=mycol, solid=.5) 

dapc.full$pca.loadings # PCA loadings
dapc.full$tab # PC values; these were appended to the datafile for use in MANCOVA
dapc.full$pca.eig/sum(dapc.full$pca.eig)*100 # % variation explained by each pc
contribA <- loadingplot(dapc.full$var.contr, axis=1, thres=.07, lab.jitter=1) # variables with highest loadings

### FC2: Repeat with "RSP" as group (comparing diploid sexual, diploid asexual, triploid asexual)
### FC4: Repeat with "mix" as group (comparing pure asexual populations with mixed sexual/asexual)


### 2. FC3: Compare sexual retrofracta, sexual stricta, asexual retrofracta, asexual retrofracta x stricta

### SppNum contains codes for each species/hybrid combination; repeat above analyses with SppNum as group
### retrofracta=1, stricta=2, asexual retro=11, asexual hybrid=12

### Remove asexual stricta (code 22) prior to analyses
FC3data <- subset(FC3data, FC3data$SppNum != 22) 

spp <- FC3data[,2] # SppNum: species/hybrid codes
pop <- FC3data[,3] # population code number
vars <- FC3data[,8:31] # non-standardized environmental data

dapc.FC3 <- dapc(vars, grp=spp, n.pca=NULL, center=TRUE, scale=TRUE)




#######################################################################################################
##############################                  MANCOVA                  ##############################
#######################################################################################################


### PCs retained from DAPC (dapc.full$tab, dapc.FC3$tab), appended to data file
### Note that PCs are not included in data files provided because they vary with analysis!
myYs <- cbind(fulldata$PC1, fulldata$PC2, fulldata$PC3, fulldata$PC4, fulldata$PC5)
FC3Ys <- cbind(FC3data$PC1, FC3data$PC2, FC3data$PC3, FC3data$PC4, FC3data$PC5)
FC4Ys <- cbind(mixdata$PC1, mixdata$PC2, mixdata$PC3, mixdata$PC4, mixdata$PC5)

### MANCOVA for FC1 with population as a covariate
FC1model <- manova(myYs ~ RS + Pop,data = fulldata) # RS is reproductive system, sexual/asexual
summary(FC1model, test = "Wilks")
summary.aov(FC1model)

### FC2: For comparison of diploid sexual, diploid asexual, triploid asexual, use myYs ~ RSP
### RSP is reproductive system + ploidy (1=diploid sexual, 2=diploid asexual, 3=triploid asexual)
### pairwise comparisons as below
FC2model.pw1 <- manova(myYs ~ as.factor(RSP) + Pop, data=fulldata, subset=as.factor(RSP) %in% c("1","2"))
summary(FC2model.pw1,test="Wilks")

### FC3: sexual stricta, sexual retro, asexual retro, asexual retro x stricta
FC3model <- manova(FC3Ys ~ SppNum + Pop, data = FC3data)
summary(FC3model, test = "Wilks")

### pairwise comparisons as below
FC3model.pw1 <- manova (FC3Ys ~ SppNum + Pop, data=FC3data, subset=SppNum %in% c("1","2"))
summary(FC3model.pw1,test="Wilks")

### FC4: mixed sexual/asexual populations vs. pure asexual populations
FC4model <- manova(FC4Ys ~ AMS + Pop,data = mixdata)
summary(FC4model, test = "Wilks")

### correct for multiple comparisons using function p.adjust in package 'lsmeans' (method=holm)



#######################################################################################################
#############################                    GLMs                   ###############################
#######################################################################################################


### NOTE that variables selected for models were selected based on visual assessment of PCA, run in JMP
### environmental variables were re-scaled using function "scale" following data cleaning for each subsequent analysis

popdata <- read.csv("PopEnvironmentalData_Rushworthetal_AJB_2018.csv", header=T) # full data, population level
FC3popdata <- read.csv("FC3PopEnvironmentalData_Rushworthetal_AJB_2018.csv", header=T) # FC3 data, population level

### Work with "popdata" to get it ready for FC4
mixpopdata <- popdata[!is.na(popdata$AMS),] 
### removing NAs limits the dataset to pure asexual and mixed sexual/asexual populations
mixpopdata <- subset(mixpopdata, PopSize>3) # removing small populations
mixpopdata$AMS <- as.factor(mixpopdata$AMS)


### FC1: sexuals vs. asexuals 
### NumSex is the number of sexual individuals in a population
### (PopSize - NumSex) gives the number of asexuals in a population
### populations containing fewer than 4 individuals have already been removed
fullvars <- popdata[,10:33] # select environmental variables
stand.full <- scale(fullvars) 
### I found that appending the standardized variables to the dataset and renaming them was easiest

FC1glm <- glm(cbind(NumSex,PopSize-NumSex) ~ BIO1 + BIO7 +  BIO12 + BIO15 + Elev + Slope + Distroad, 
              data=popdata, 
              family=binomial())

### Moran's I to look for evidence of geographic structure in residuals; repeated for models below
ll.dists <- rdist.earth(cbind(popdata$Long,popdata$Lat),cbind(popdata$Long,popdata$Lat),miles=FALSE)
ll.dists.inv <- 1/ll.dists # calculate inverse distances
diag(ll.dists.inv) <- 0 # make diagonals zeros

obs <- (pops$NumSex/pops$PopSize) 
FC1resid <- (cbind(obs, fitted.values(FC1glm))) # fitted values are the predicted values
### residuals are obs - fitted values; did this by hand and saved as column Resids
Moran.I(FC1resid$Resids,ll.dists.inv)


### FC2: Remove triploids (limiting comparison to diploid sexuals and diploid asexuals)
### NumDips is the number of diploid individuals in a population
### (NumDips - NumSex) gives the number of diploid asexuals in a population

FC2popdata <- subset(popdata, NumDips>3) # remove small populations
### rescale environmental data here using code above (lines 141-142) prior to running model

FC2glm <- glm(cbind(NumSex,NumDips-NumSex) ~ BIO1 + BIO7 + BIO12 + BIO15 + Elev + Slope + Distroad, 
               data=FC2popdata, 
               family=binomial())


### FC3: sexual retrofracta, sexual stricta, asexual retrofracta, asexual retrofracta x stricta
### NumInd is the population size only including the four focal taxa
### (NumInd - NumSex) gives the number of asexuals in the population
### remove small populations (N<4)
FC3vars <- FC3popdata[,13:36]
stand.FC3 <- scale(FC3vars) # scale environmental variables

FC3glm <- glm(cbind(NumSex,NumInd-NumSex) ~ BIO1 + BIO3 + BIO4 + BIO15 + Elev + Slope + Distroad, 
               data=FC3popdata,
               family=binomial())


### FC4: Mixed sexual/asexual vs. pure asexual pops
### note that this dataset 'mixpops' is the dataset 'pops' with pure sexual populations removed (see above)
mixvars <- mixpopdata[,10:33]
stand.mix <- scale(mixvars) # rescale following removal of those populations

FC4glm <- glm(AMS ~ BIO3 + BIO7 + BIO12 + BIO15 + Slope + Distroad, 
              data=mixpopdata, 
              family=binomial())

### correction for multiple comparisons using function p.adjust in package 'lsmeans' (method=holm)


