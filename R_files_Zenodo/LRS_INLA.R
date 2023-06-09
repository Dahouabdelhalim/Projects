#Helgeland house sparrow inbreeding depression analyses using INLA
#Estimating the effect of inbreeding on lifetime reproductive success (LRS) and proportion of variance explained by F within each island
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
require(arm)
library(parallel)
library(SMisc)
library(ggplot2)


#Import the phenotypic data and inbreeding estimates
data_LRS <- read.csv("Data.txt", header = T, stringsAsFactors = F, sep = "\\t")

#Import population size estimate that include only hatch year population sizes between 1997-2012 for inner islands and 2003-2012 for all outer islands
Pop_size_1997_2012 <- read.csv("Pop_size_1997_2012.csv", header = T, stringsAsFactors = F, sep = "\\t")

###
###Subsetting the data
###

#Include only the individuals that have been observed as adults on one of the 8 study islands
data_adult_temp <- data_LRS[which(data_LRS$laflok %in% c("20","22","23","24","26","27","28","38")),]
table(data_adult_temp$laflok)


###
###Lifetime recruit production
###

#Only those individuals that have LRS_OK estimate and have been adults on one of the study islands
data_ok <- data_adult_temp[!is.na(data_adult_temp$LRS_OK),] 


#Remove temp files
rm(list=ls(pattern="temp"))

###
###Prepare the variables for the dataset that has correct LRS data for adult islands
###

#Scale (center) FGRM and FROH
data_ok$c_FGRM <- scale(data_ok$FGRM, center = T, scale = F)
data_ok$c_FGRM2 <- data_ok$c_FGRM
data_ok$c_FGRM_q <- data_ok$c_FGRM^2
data_ok$c_FROH <- scale(data_ok$FROH, center = T, scale = F)
data_ok$z_FGRM <- scale(data_ok$FGRM, center = T, scale = T)
data_ok$z_FROH <- scale(data_ok$FROH, center = T, scale = T)

# make hatchyear and island factor covariates
data_ok$f_all_hatchyears <- as.factor(data_ok$all_hatchyears)
data_ok$f_laflok <- as.factor(data_ok$laflok)

#Add combined variable for hatch year and adult island
data_ok$IAY <- as.factor(paste(data_ok$all_hatchyears,data_ok$laflok,sep="_"))


#Import the pedigree
ped <- read.table("pedigree.txt", header = T, stringsAsFactors = F)



############################################################
### Animal model to account for relatedness among individuals
##############################################################

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
data_ok$ID <- d.map[match(data_ok$id, d.map$id), "ID"]

#Add new columns for random variables
data_ok$f_laflok2 <- data_ok$f_laflok
table(data_ok$laflok)
data_ok$f_all_hatchyears2 <- data_ok$f_all_hatchyears
data_ok$IAY2 <- data_ok$IAY


###
###Population size/density
###

#Combine the population size information for each individual for their hatch year and hatch island
data_ok_hatchyear <- merge(data_ok, Pop_size_1997_2012[,c(1,3:8)], by.x=c("all_hatchyears","fiflok"), by.y=c("Year","Flok"))

#Rename the columns for first adult year
colnames(data_ok_hatchyear)[59:63] <- c("BestPopEst_hatchy","g_pop_size_hatchy","Correct_size_hatchy","rel.n_hatchy","mean.n_hatchy")

#Mean center the mean population size of hatch year
data_ok_hatchyear$c_mean.n_hatchy <- scale(data_ok_hatchyear$mean.n_hatchy, center = T, scale = F)

#Mean center the population size of hatch year
data_ok_hatchyear$c_cor_size_hatchy <- scale(data_ok_hatchyear$Correct_size_hatchy, center = T, scale = F)

##########
#Defining and running INLA models
##########

#Final model, including FGRM:laflok and FGRM:year_island interactions in the random term, FGRM scaled
formula_7 = LRS_OK ~ z_FGRM + gen_sex +  
  f(f_laflok,model="iid",param = c(0.1,0.01)) +
  f(f_laflok2,z_FGRM,model="iid",param = c(0.1,0.01)) +
  f(IAY, model = "iid",param = c(0.1,0.01)) +
  f(IAY2, z_FGRM, model = "iid",param = c(0.1,0.01)) +
  f(ID,model="generic0",
    Cmatrix=Cmatrix,
    hyper=list(
      prec=list( param = c(0.1,0.01))
    )
  ) 

#Run INLA
model.zip_7 = inla(formula=formula_7, family="zeroinflatedpoisson1",
                   data=data_ok,
                   control.family=list(link='log'),
                   control.compute=list(dic=T, config=T), 
                   control.predictor=list(link=1, compute=TRUE)
)
summary(model.zip_7)

inla_emarginal(model.zip_7)
inla_mmarginal(model.zip_7)

#"Full model" including sex, habitat type and annual population size interactions with FGRM
#Centered annual population size of the parental hatch year is used
#And island*z_FGRM + island-year*z_FGRM interactions in the random term
formula_pop_size_8_full = LRS_OK ~ z_FGRM + gen_sex*z_FGRM + laflok_status*z_FGRM + c_cor_size_hatchy*z_FGRM +
  f(f_laflok,model="iid",param = c(0.1,0.01)) +
  f(f_laflok2,z_FGRM,model="iid",param = c(0.1,0.01)) +
  f(IAY, model = "iid",param = c(0.1,0.01)) +
  f(IAY2,z_FGRM, model = "iid",param = c(0.1,0.01)) +
  f(ID,model="generic0",
    Cmatrix=Cmatrix,
    hyper=list(
      prec=list( param = c(0.1,0.01))
    )
  )

#Run INLA
model.zip_pop8_full = inla(formula=formula_pop_size_8_full, family="zeroinflatedpoisson1",
                                   data=data_ok_hatchyear,
                                   control.compute=list(dic=T, config=T)
)
summary(model.zip_pop8_full)



#Table the variances of the random effects

### First make a function for transforming precisions to variances
inlaPosteriors_Zip_4 <- function(model){
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

  sigma.e = inla.tmarginal(function(x) 1/x,  model$marginals.hyperpar$"zero-probability parameter for zero-inflated poisson_1")
  m.e =inla.mmarginal(sigma.e)
  e.e =inla.emarginal(function(x) x, sigma.e)  

  results_tab <- cbind(rbind(zero_probability = e.e, island = e.island, island_2 = e.island_2, Island_year=e.IAY, Island_year_2=e.IAY2, ID= e.ID),
                       rbind(zero_probability = m.e, island = m.island, island_2 = m.island_2, Island_year=m.IAY, Island_year_2=m.IAY2, ID= m.ID),
                       1/model$summary.hyperpar[,c(5,3)])
  
  names(results_tab) <- c("Mean","Mode", "2.5%", "97.5 %")
  return(results_tab)
}


###Posterior variances of the random effects model 7 'Final model'
inlaPosteriors_Zip_4(model.zip_7)

###Posterior variances of the random effects of "Full model" including population size and habitat type
inlaPosteriors_Zip_4(model.zip_pop8_full)


###################################################################
### Calculate an approximate measure of variable importance
### Start by doing this for the entire population, irrespective of island, and only for the model labelled with "7" (final model)
### To this end we use the function F_importance that we code here first
F_importance <- function(r.out, data){
  # Variances by the random effects given as the sum of all posterior models of the estimated variance components:
  var_random <- sum(inla_emarginal(r.out)[-1]) 
  
  # Variances by the two fixed effects
  var_F <- r.out$summary.fixed["z_FGRM",1]^2 * var(data$z_FGRM)
  var_sex <- r.out$summary.fixed["gen_sexm",1]^2 * var(as.numeric(as.factor(data$gen_sex)))
  
  # Sum of variances explained by the fixed effects
  var_fixed <-  var_F + var_sex 
  
  # And finally the proportion explained by FGRM
  importance <- var_F / (var_fixed + var_random)
  return(importance)
}
# Now use the function to calculate the proportion of variance explained by F in the full model
F_importance(model.zip_7, data_ok)


### Now do it island-wise
# Check number of data points per island:
table(data_ok$laflok)
# Generate vector of islands
islands <- sort(unique(data_ok$laflok))

# Fit INLA-models island-wise and calculate proportion of variance explained by F:
# INLA formula (without laflok, because each model will be run only for one island)
formula = LRS_OK ~ z_FGRM + gen_sex +  
  f(IAY, model = "iid",param = c(0.1,0.01)) +
  f(ID,model="generic0",
    Cmatrix=Cmatrix,
    hyper=list(
      prec=list( param = c(0.1,0.01))
    )
  )

# Then run (in parallel) the island-specific models
F_impo_islands_LRS <- mclapply(islands, function(i) {
  data.tmp <- subset(data_ok,laflok==i)
 #Run INLA
  model.out = inla(formula=formula, family="zeroinflatedpoisson1",
                        data=data.tmp,
                        control.family=list(link='log'),
                        control.compute=list(dic=T, config=T), 
                        control.predictor=list(link=1, compute=TRUE),
                   num.threads=2
  )
  
  # Return the fixed effects estimates and F-importances:
  return(list(model.out$summary.fixed,F_importance(model.out,data.tmp),inla_mmarginal(model.out)))
}, mc.cores =2)


# Generate a temporary file that contains the island names and mean population sizes
tmp <- unique(Pop_size_1997_2012[,c("Flok","mean.n")])

# The following is the list of proportions of variance explained by F for each island:
(d.variance_LRS <- data.frame(
  island=islands,
  variance_explained_F = do.call(rbind, lapply(F_impo_islands_LRS, `[[`, 2)),
  size = tmp[tmp$Flok==islands,2]  # data.frame(table(data_ok$laflok))[,2]
))

ggplot(data=d.variance_LRS,aes(x=log(size),y=log(variance_explained_F))) + geom_point(size=2) + 
  xlab("log(Size)")+ 
  ylab("log(Variance explained by F)") +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold")) +
  theme_bw()

summary(lm(log(variance_explained_F) ~ log(size), d.variance_LRS))


# Check out the two smallest pops and what the difference is there
# Pops 20 and 22
F_impo_islands_LRS[[1]]
F_impo_islands_LRS[[2]]






