#Helgeland house sparrow inbreeding depression analyses using INLA
#Estimating the effect of inbreeding on annual reproductive success (AR) and proportion of variance explained by F within each island
#######################################################################
#Alina Niskanen & Stefanie Muff
#alina.niskanen@gmail.com
#April 2020


library(nadiv)
library(INLA)
library(lme4)
library(MASS)
library(MCMCpack)
library(MasterBayes)
library(MCMCglmm)
library(lmerTest)
library(pedigreemm)
library(tidyr)
library(pscl)
require(arm)
library(SMisc)
library(parallel)
library(ggplot2)


#Import the phenotypic data and inbreeding estimates
data_LRS <- read.csv("Data.txt", header = T, stringsAsFactors = F, sep = "\\t")

#Import Best population size estimates combination of observational population size and genotyped adult sample size per island for all study years.
Pop_size_w_SNPs <- read.csv("Pop_size_1998_2013.csv", header = T, stringsAsFactors = F, sep = "\\t") 

###
###Subsetting the data
###

#Include only the individuals that have been observed as adults on one of the 8 study islands
data_adult_temp <- data_LRS[which(data_LRS$laflok %in% c("20","22","23","24","26","27","28","38")),]
table(data_adult_temp$laflok)
#Only birds that hatched on one of the 8 study islands
data_hatch_temp <- data_LRS[which(data_LRS$fiflok %in% c("20","22","23","24","26","27","28","38")),]
table(data_hatch_temp$fiflok)

###
###Per year recruit production
###

#Start choosing the correct data from "data_adult_temp" that includes all individuals that were observed as adults (laflok) on one of the 8 study islands
data_no_YR <- data_adult_temp[(is.na(data_adult_temp$X1998) & is.na(data_adult_temp$X1999) & is.na(data_adult_temp$X2000) & is.na(data_adult_temp$X2001) & is.na(data_adult_temp$X2002) & is.na(data_adult_temp$X2003) & is.na(data_adult_temp$X2004) & is.na(data_adult_temp$X2005) & is.na(data_adult_temp$X2006) & is.na(data_adult_temp$X2007) & is.na(data_adult_temp$X2008) & is.na(data_adult_temp$X2009) & is.na(data_adult_temp$X2010) & is.na(data_adult_temp$X2011) & is.na(data_adult_temp$X2012)),]

#Make the final yearly estimate dataset, includes only individuals with laflok on one of the 8 study islands and YR estimate available
data_YR_temp <- data_adult_temp[!data_adult_temp$id %in% data_no_YR$id,]

#Scale (center) FGRM
data_YR_temp$z_FGRM <- scale(data_YR_temp$FGRM, center = T, scale = T)
data_YR_temp$z_FROH <- scale(data_YR_temp$FROH, center = T, scale = T)

#Transpose the data to include each individual yearly estimate on its own row
data_YR <- gather(data_YR_temp, "Year", "recruits", 22:36)
data_YR$Year <- gsub("X","", data_YR$Year)

#Remove all rows with missing YR to reduce the file size
data_YR <- data_YR[!is.na(data_YR$recruits),]
table(data_YR$recruits)
hist(data_YR$recruits)
#Add age of the individual when having the offspring
data_YR$Year <- as.integer(data_YR$Year)
data_YR$age <- data_YR$Year-data_YR$all_hatchyears

#Remove temp files
rm(list=ls(pattern="temp"))


###
###Prepare the variables for the dataset that has correct YR data for adult islands
###

# make island factor
data_YR$f_laflok <- as.factor(data_YR$laflok)
#Combine age classes 6 and above
data_YR$age_comb <- data_YR$age
data_YR$age_comb[which(data_YR$age>=6)] <- 6
#Center combined age
data_YR$c_age_comb <- scale(data_YR$age_comb, center = T, scale = F)
data_YR$c_age_comb2 <- data_YR$c_age_comb^2

#Add combined variable for hatch year and adult island
data_YR$IAY <- as.factor(paste(data_YR$all_hatchyears,data_YR$laflok,sep="_"))
#Add combined variable for hatch year and adult island to be used in
data_YR$IAY2 <- data_YR$IAY


################
#Add population size information to the individual data
################

#Combine population size information with F-data
data_YR <- merge(data_YR, Pop_size_w_SNPs, by.x=c("Year","laflok"), by.y=c("Year","Flok"))

#Mean center the annual population size
data_YR$c_Correct_size <- scale(data_YR$Correct_size, center = T, scale = F)

#Plot YR vs population size
plot(data_YR$Correct_size,data_YR$recruits)
boxplot(data_YR$recruits~data_YR$Correct_size)


############################################################
### Animal model to account for relatedness among individuals
##############################################################

#Import the pedigree
#Column names: individual="id", dam="dam", sire="sire"
ped <- read.table("pedigree.txt", header = T, stringsAsFactors = F)

# need an ordered Pedigree for the inverseA() function:
ped <- orderPed(ped)

# introduce the "ID", a new identity for individuals, enumerated from 1 to number of ind.:
ped$ID <- 1:(nrow(ped))

# fathers and mothers are replaced by these new IDs
d.map <- ped[,c("id","ID")]
ped$dam.id <- d.map[match(ped$dam, d.map$id),"ID"]
ped$sire.id <- d.map[match(ped$sire, d.map$id),"ID"]

# compute A inverse, using the new IDs:
Cmatrix <- inverseA(ped[,c("ID","dam.id","sire.id")])$Ainv

# also need to add ID column to data file:
data_YR$ID <- d.map[match(data_YR$id, d.map$id), "ID"]
data_YR$ID2 <- data_YR$ID

data_YR$f_laflok2 <- data_YR$f_laflok
table(data_YR$laflok)
data_YR$f_all_hatchyears2 <- data_YR$f_all_hatchyears

##########
#Defining and running INLA models
##########

#"Final model" including centered age as quadratic term, FGRM:laflok and FGRM:hatchyear_laflok random interaction, FGRM scaled
formula_13_YR = recruits ~ z_FGRM + gen_sex + c_age_comb + c_age_comb2 +
  f(f_laflok,model="iid",param = c(0.1,0.01)) +
  f(f_laflok2,z_FGRM,model="iid",param = c(0.1,0.01)) +
  f(IAY, model = "iid",param = c(0.1,0.01)) +
  f(IAY2, z_FGRM, model = "iid",param = c(0.1,0.01)) +
    f(ID,model="generic0",
    Cmatrix=Cmatrix,
    hyper=list(
      prec=list( param = c(0.1,0.01))
    )) +
  f(ID2,model="iid",param = c(0.1,0.01),
    constr=TRUE)

#Run INLA
model.poi_YR_13 = inla(formula=formula_13_YR, family="poisson",
                       data=data_YR,
                       control.compute=list(dic=T, config=T))
summary(model.poi_YR_13)

inla_emarginal(model.poi_YR_13)
inla_mmarginal(model.poi_YR_13)


#"Full model" including centered age as quadratic term, sex:z_FGRM, island_type:z_FGRM and centered annual population size:z_FGRM interactions
#Random effects FGRM:laflok and FGRM:hatchyear_laflok interaction
formula_pop_size_YR_full2 = recruits ~ z_FGRM + c_age_comb + c_age_comb2 + gen_sex*z_FGRM + laflok_status*z_FGRM + c_Correct_size*z_FGRM +
  f(f_laflok,model="iid",param = c(0.1,0.01)) +
  f(f_laflok2,z_FGRM,model="iid",param = c(0.1,0.01)) +
  f(IAY, model = "iid",param = c(0.1,0.01)) +
  f(IAY2, z_FGRM, model = "iid",param = c(0.1,0.01)) +
  f(ID,model="generic0",
    Cmatrix=Cmatrix,
    hyper=list(
      prec=list( param = c(0.1,0.01))
    )) +
  f(ID2,model="iid",param = c(0.1,0.01),
    constr=TRUE)

#Run INLA
model.poi_pop_size_YR_full2 = inla(formula=formula_pop_size_YR_full2, family="poisson",
                                  data=data_YR,
                                  control.compute=list(dic=T, config=T))
summary(model.poi_pop_size_YR_full2)


#Table the variances of the random effects

### First make a function for transforming precisions to variances
inlaPosteriors_Poi_3 <- function(model){
  sigma.island = inla.tmarginal(function(x) 1/x,  model$marginals.hyperpar$"Precision for f_laflok")
  m.island=inla.mmarginal(sigma.island)
  e.island=inla.emarginal(function(x) x, sigma.island)

  sigma.island_2 = inla.tmarginal(function(x) 1/x,  model$marginals.hyperpar$"Precision for f_laflok2")
  m.island_2=inla.mmarginal(sigma.island_2)
  e.island_2=inla.emarginal(function(x) x, sigma.island_2)  

  sigma.ID = inla.tmarginal(function(x) 1/x,  model$marginals.hyperpar$"Precision for ID")
  m.ID=inla.mmarginal(sigma.ID)
  e.ID=inla.emarginal(function(x) x, sigma.ID)
  
  sigma.IAY = inla.tmarginal(function(x) 1/x,  model$marginals.hyperpar$"Precision for IAY")
  m.IAY=inla.mmarginal(sigma.IAY)
  e.IAY=inla.emarginal(function(x) x, sigma.IAY)

  sigma.IAY2 = inla.tmarginal(function(x) 1/x,  model$marginals.hyperpar$"Precision for IAY2")
  m.IAY2=inla.mmarginal(sigma.IAY2)
  e.IAY2=inla.emarginal(function(x) x, sigma.IAY2) 

  sigma.animal = inla.tmarginal(function(x) 1/x,  model$marginals.hyperpar$"Precision for ID2")
  m.animal=inla.mmarginal(sigma.animal)
  e.animal=inla.emarginal(function(x) x, sigma.animal)
  
  results_tab <- cbind(rbind(island = e.island, island_2 = e.island_2, Island_year=e.IAY, Island_year_2=e.IAY2, ID= e.ID, varA=e.animal),
                       rbind(island = m.island, island_2 = m.island_2, Island_year=m.IAY, Island_year_2=m.IAY2, ID= m.ID, varA=m.animal),
                       1/model$summary.hyperpar[,c(5,3)])
  
  names(results_tab) <- c("Mean","Mode", "2.5%", "97.5 %")
  return(results_tab)
}

###Posterior variances of the random effects of model 13 "Final model"
inlaPosteriors_Poi_3(model.poi_YR_13)

###Posterior variances of the random effects of "Full model" including population size and habitat type
inlaPosteriors_Poi_3(model.poi_pop_size_YR_full2)



###################################################################
### Calculate an approximate measure of variable importance
### Start by doing this for the entire population, irrespective of island, and only for the model labelled "model.poi_YR_13" (final model)
### To this end we use the function F_importance that we code here first
F_importance <- function(r.out,data){
  # Variances by the random effects given as the sum of all posterior models of the estimated variance components:
  var_random <- sum(inla_emarginal(r.out)) 
  
  # Variances by the two fixed effects
  var_F <- r.out$summary.fixed["z_FGRM",1]^2 * var(data$z_FGRM)
  var_sex <- r.out$summary.fixed["gen_sexm",1]^2 * var(as.numeric(as.factor(data$gen_sex)))
  var_age <- r.out$summary.fixed["c_age_comb",1]^2 * var(as.numeric(as.factor(data$c_age_comb)))
  var_age2 <- r.out$summary.fixed["c_age_comb2",1]^2 * var(as.numeric(as.factor(data$c_age_comb2)))
  
  # Sum of variances explained by the fixed effects
  var_fixed <-  var_F + var_sex + var_age + var_age2
  
  # And finally the proportion explained by FGRM
  importance <- var_F / (var_fixed + var_random)
  return(importance)
}


# Now use the function to calculate the proportion of variance explained by F in the full model
F_importance(model.poi_YR_13,data_YR)
var(data_YR$z_FGRM)

### Now do it island-wise
# Check number of data points per island:
data.frame(table(data_YR$laflok))
# Generate vector of islands
islands <- sort(unique(data_YR$laflok))

# Fit INLA-models island-wise and calculate proportion of variance explained by F:
# INLA formula (without laflok, because each model will be run only for one island)
formula = recruits ~ z_FGRM + gen_sex + c_age_comb + c_age_comb2 +
  f(IAY, model = "iid",param = c(0.1,0.01)) +
  f(ID,model="generic0",
    Cmatrix=Cmatrix,
    hyper=list(
      prec=list( param = c(0.1,0.01))
    )) +
  f(ID2,model="iid",param = c(0.1,0.01),
    constr=TRUE)

# Run island-wise analyses in parallel
F_impo_islands_AR <- mclapply(islands, function(i) {
  # Temporary dataset with i'th islands data
  data.tmp <- subset(data_YR,laflok==i)
  #Run INLA
  model.out = inla(formula=formula, family="poisson",
       data=data.tmp,
       control.compute=list(dic=T, config=T))
  
  # Return the fixed effects estimates and F-importances:
  return(list(model.out$summary.fixed,F_importance(model.out,data.tmp),inla_emarginal(model.out)))
}, mc.cores =2)

###
# Generate a temporary file that contains the island names and mean population sizes
tmp <- unique(data_YR[,c("laflok","mean.n")])
 

# The following is the list of proportions of variance explained by F for each island, plus the size of the islands
(d.variance_AR <- data.frame(
  island=islands,
  variance_explained_F = do.call(rbind, lapply(F_impo_islands_AR, `[[`, 2)),
  size = tmp[tmp$flok==islands,2]  #data.frame(table(data_YR$laflok))[,2]
))


ggplot(data=d.variance_AR,aes(x=log(size),y=log(variance_explained_F))) + geom_point(size=2) + 
  xlab("log(Size)")+ 
  ylab("log(Variance explained by F)") +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold")) +
  theme_bw()

summary(lm(log(variance_explained_F) ~ log(size), d.variance_AR))

# Check out the two smallest pops and what the difference is there
# Pops 20 and 22
F_impo_islands_AR[[1]]
F_impo_islands_AR[[2]]




