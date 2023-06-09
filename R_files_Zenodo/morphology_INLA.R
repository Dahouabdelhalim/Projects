#Helgeland house sparrow inbreeding depression analyses using INLA
#Estimating the effect of inbreeding on morphological traits
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


#Import the phenotypic data, inbreeding estimates and the morphology data, adjusted for May of age 1 year
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
###Include individuals with morphology information
###

#Only those individuals that have age1mass estimate and have been adults on one of the study islands
data_morph <- data_adult_temp[!is.na(data_adult_temp$age1mass),]
table(data_morph$laflok) 
table(data_morph$fiflok) 

#Remove temp files
rm(list=ls(pattern="temp"))

###
###Prepare the variables for the dataset that has correct morpohology data for adult islands
###

#Scale FGRM
data_morph$z_FGRM <- scale(data_morph$FGRM, center = T, scale = T)

# make island factor covariate
data_morph$f_laflok <- as.factor(data_morph$laflok)
data_morph$f_laflok2 <- data_morph$f_laflok

#Add combined variable for hatch year and island
data_morph$IAY <- as.factor(paste(data_morph$all_hatchyears,data_morph$laflok,sep="_"))
data_morph$IAY2 <- data_morph$IAY

############################################################
### Animal model to account for relatedness among individuals
##############################################################

#Import the pedigree
ped <- read.table("SNP_pedigree_Helgeland_05122017.txt", header = T, stringsAsFactors = F)

# # need an ordered Pedigree for the inverseA() function:
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
data_morph$ID <- d.map[match(data_morph$id, d.map$id), "ID"]


###
###Population size
###

#Combine here the population size information for each individual for their hatch year and hatch island
data_morph_hatchyear <- merge(data_morph, Pop_size_1997_2012, by.x=c("all_hatchyears","fiflok"), by.y=c("Year","Flok"))
#Rename the columns for first adult year
names(data_morph_hatchyear)[names(data_morph_hatchyear) == "Correct_size"] <- c("Correct_size_hatchy")

#Mean center the annual population size of hatch year
data_morph_hatchyear$c_Correct_size_hatchy <- scale(data_morph_hatchyear$Correct_size_hatchy, center = T, scale = F)


##########
#Defining and running INLA models
##########


#########
#Beak length
#########


#Model with FGRM:island and FGRM:island_year interactions, FGRM scaled, "Final model".
formula_6_bL = age1billL ~ z_FGRM + gen_sex +
  f(f_laflok,model="iid", param = c(0.1,0.01)) +
  f(f_laflok2, z_FGRM, model="iid",param = c(0.1,0.01)) +
  f(IAY, model = "iid",param = c(0.1,0.01)) +
  f(IAY2, z_FGRM, model = "iid", param = c(0.1,0.01)) +
  f(ID,model="generic0",
    Cmatrix=Cmatrix,
    hyper=list(
      prec=list( param = c(0.1,0.01))
    )
  )

#Run INLA
model.beak_L_6 = inla(formula=formula_6_bL, family="gaussian",
                      data=data_morph,
                      control.compute=list(dic=T, config=T)
)
summary(model.beak_L_6)

#Full model with centered annual population size for the parental hatch year and habitat type.
formula_bL_size_full2 = age1billL ~ z_FGRM + gen_sex*z_FGRM + laflok_status*z_FGRM + c_Correct_size_hatchy*z_FGRM +
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
model.beak_L_full2 = inla(formula=formula_bL_size_full2, family="gaussian",
                                 data=data_morph_hatchyear,
                                 control.compute=list(dic=T, config=T)
)
summary(model.beak_L_full2)


#Table the variances of the random effects

### First make a function for transforming precisions to variances
inlaPosteriors_morph_4 <- function(model){
  sigma.island = inla.tmarginal(function(x) 1/x,  model$marginals.hyperpar$"Precision for f_laflok")
  m.island=inla.mmarginal(sigma.island)
  e.island=inla.emarginal(function(x) x, sigma.island)

  sigma.island2 = inla.tmarginal(function(x) 1/x,  model$marginals.hyperpar$"Precision for f_laflok2")
  m.island2=inla.mmarginal(sigma.island2)
  e.island2=inla.emarginal(function(x) x, sigma.island2)  

  sigma.ID = inla.tmarginal(function(x) 1/x,  model$marginals.hyperpar$"Precision for ID")
  m.ID=inla.mmarginal(sigma.ID)
  e.ID=inla.emarginal(function(x) x, sigma.ID)
  
  sigma.IAY = inla.tmarginal(function(x) 1/x,  model$marginals.hyperpar$"Precision for IAY")
  m.IAY=inla.mmarginal(sigma.IAY)
  e.IAY=inla.emarginal(function(x) x, sigma.IAY)
  
  sigma.IAY2 = inla.tmarginal(function(x) 1/x,  model$marginals.hyperpar$"Precision for IAY2")
  m.IAY2=inla.mmarginal(sigma.IAY2)
  e.IAY2=inla.emarginal(function(x) x, sigma.IAY2) 

  sigma.e = inla.tmarginal(function(x) 1/x,  model$marginals.hyperpar$"Precision for the Gaussian observations")
  m.e =inla.mmarginal(sigma.e)
  e.e =inla.emarginal(function(x) x, sigma.e)
  
  results_tab <- cbind(rbind(VarE = e.e, island = e.island, island2 = e.island2, Island_year=e.IAY, Island_year_2=e.IAY2, ID= e.ID),
                       rbind(VarE = m.e, island = m.island, island2 = m.island2, Island_year=m.IAY, Island_year_2=m.IAY2, ID= m.ID),
                       1/model$summary.hyperpar[,c(5,3)])
  
  names(results_tab) <- c("Mean","Mode", "2.5%", "97.5 %")
  return(results_tab)
}


###Posterior variances  of the random effects model 6, "final model".
inlaPosteriors_morph_4(model.beak_L_6)

###Posterior variances of the random effects of "Full model" including population size and habitat type.
inlaPosteriors_morph_4(model.beak_L_full2)



#################
#Beak Depth
#################

#Model with FGRM:island and FGRM:island_year interactions, FGRM scaled, "Final model".
formula_6_bD = age1billD ~ z_FGRM + gen_sex +
  f(f_laflok,model="iid", param = c(0.1,0.01)) +
  f(f_laflok2, z_FGRM, model="iid",param = c(0.1,0.01)) +
  f(IAY, model = "iid",param = c(0.1,0.01)) +
  f(IAY2, z_FGRM, model = "iid", param = c(0.1,0.01)) +
  f(ID,model="generic0",
    Cmatrix=Cmatrix,
    hyper=list(
      prec=list( param = c(0.1,0.01))
    )
  )

  )

#Run INLA
model.beak_D_6 = inla(formula=formula_6_bD, family="gaussian",
                       data=data_morph,
                       control.compute=list(dic=T, config=T)
)
summary(model.beak_D_6)

#Full model with centered annual population size for the parental hatch year and habitat type.
formula_bD_size_full2 = age1billD ~ z_FGRM + gen_sex*z_FGRM + laflok_status*z_FGRM + c_Correct_size_hatchy*z_FGRM +
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
model.beak_D_full2 = inla(formula=formula_bD_size_full2, family="gaussian",
data=data_morph_hatchyear,
control.compute=list(dic=T, config=T)
)
summary(model.beak_D_full2)


#Table the variances of the random effects

###Model 6 beak depth "final model".
inlaPosteriors_morph_4(model.beak_D_6)

###Full model with centered annual population size for the parental hatch year and habitat type.
inlaPosteriors_morph_4(model.beak_D_full2)



##############################
#Tarsus length
##############################

# "Final model".
formula_6_T = age1tarsus ~ z_FGRM + gen_sex +
  f(f_laflok,model="iid", param = c(0.1,0.01)) +
  f(f_laflok2, z_FGRM, model="iid",param = c(0.1,0.01)) +
  f(IAY, model = "iid",param = c(0.1,0.01)) +
  f(IAY2, z_FGRM, model = "iid", param = c(0.1,0.01)) +
  f(ID,model="generic0",
    Cmatrix=Cmatrix,
    hyper=list(
      prec=list( param = c(0.1,0.01))
    )
  )

#Run INLA
model.tarsus_6 = inla(formula=formula_6_T, family="gaussian",
                       data=data_morph,
                       control.compute=list(dic=T, config=T)
)
summary(model.tarsus_6)

#Full model.
formula_T_size_full2 = age1tarsus ~ z_FGRM + gen_sex*z_FGRM + laflok_status*z_FGRM + c_Correct_size_hatchy*z_FGRM +
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
model.tarsus_full2 = inla(formula=formula_T_size_full2, family="gaussian",
                                 data=data_morph_hatchyear,
                                 control.compute=list(dic=T, config=T)
)
summary(model.tarsus_full2)


#Table the variances of the random effects
###Model 6 tarsus "final model"
inlaPosteriors_morph_4(model.tarsus_6)

###Full model with centered annual population size for the parental hatch year and habitat type.
inlaPosteriors_morph_4(model.tarsus_full2)



##############################
#Wing length
##############################


#Model with FGRM:island and FGRM:island_year interactions, FGRM scaled, "Final model".
formula_6_W = age1wing ~ z_FGRM + gen_sex +
  f(f_laflok,model="iid", param = c(0.1,0.01)) +
  f(f_laflok2, z_FGRM, model="iid",param = c(0.1,0.01)) +
  f(IAY, model = "iid",param = c(0.1,0.01)) +
  f(IAY2, z_FGRM, model = "iid", param = c(0.1,0.01)) +
  f(ID,model="generic0",
    Cmatrix=Cmatrix,
    hyper=list(
      prec=list( param = c(0.1,0.01))
    )
  )

model.wing_6 = inla(formula=formula_6_W, family="gaussian",
                    data=data_morph,
                    control.compute=list(dic=T, config=T)
)
summary(model.wing_6)

#Full model with centered annual population size for the parental hatch year and habitat type.
formula_W_size_full2 = age1wing ~ z_FGRM + gen_sex*z_FGRM + laflok_status*z_FGRM + c_Correct_size_hatchy*z_FGRM +
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

model.wing_full2 = inla(formula=formula_W_size_full2, family="gaussian",
                               data=data_morph_hatchyear,
                               control.compute=list(dic=T, config=T)
)
summary(model.wing_full2)

#Table the variances of the random effects

###Model 6 wing, FGRM:island and FGRM:island_year interaction, FGRM scaled, "Final model".
inlaPosteriors_morph_4(model.wing_6)

###Full model2 with centered annual population size for the parental hatch year and habitat type.
inlaPosteriors_morph_4(model.wing_full2)


##############################
#Mass
##############################


#Model with FGRM:island and FGRM:island_year interactions, FGRM scaled, "Final model".
formula_6_M = age1mass ~ z_FGRM + gen_sex +
  f(f_laflok,model="iid", param = c(0.1,0.01)) +
  f(f_laflok2, z_FGRM, model="iid",param = c(0.1,0.01)) +
  f(IAY, model = "iid",param = c(0.1,0.01)) +
  f(IAY2, z_FGRM, model = "iid", param = c(0.1,0.01)) +
  f(ID,model="generic0",
    Cmatrix=Cmatrix,
    hyper=list(
      prec=list( param = c(0.1,0.01))
    )
  )

model.mass_6 = inla(formula=formula_6_M, family="gaussian",
                     data=data_morph,
                     control.compute=list(dic=T, config=T)
)
summary(model.mass_6)

#Full model with centered annual population size for the parental hatch year and habitat type.
formula_M_size_full2 = age1mass ~ z_FGRM + gen_sex*z_FGRM + laflok_status*z_FGRM + c_Correct_size_hatchy*z_FGRM +
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

model.mass_full2 = inla(formula=formula_M_size_full2, family="gaussian",
                               data=data_morph_hatchyear,
                               control.compute=list(dic=T, config=T)
)
summary(model.mass_full2)



#Table the variances of the random effects

###Model 6 mass, FGRM:island and FGRM:island_year interaction, FGRM scaled, "Final model".
inlaPosteriors_morph_4(model.mass_6)

###Full model with centered annual population size for the parental hatch year and habitat type.
inlaPosteriors_morph_4(model.mass_full2)


##############################
#Total badge size, only males
##############################

#Make subset only for males
d <- data_morph[which(data_morph$gen_sex=="m"),]
#Make subset only for males including hatch year population sizes
d_hatchyear <- data_morph_hatchyear[which(data_morph_hatchyear$gen_sex=="m"),]
#Make subset only for males including first adult year population sizes
d_minad <- data_morph_minad[which(data_morph_minad$gen_sex=="m"),]

##Specify the models for INLA


#Model with FGRM:island and FGRM:island_year interactions, FGRM scaled, "Final model".
formula_6_totba = age1totba ~ z_FGRM +
  f(f_laflok,model="iid", param = c(0.1,0.01)) +
  f(f_laflok2, z_FGRM, model="iid",param = c(0.1,0.01)) +
  f(IAY, model = "iid",param = c(0.1,0.01)) +
  f(IAY2, z_FGRM, model = "iid", param = c(0.1,0.01)) +
  f(ID,model="generic0",
    Cmatrix=Cmatrix,
    hyper=list(
      prec=list( param = c(0.1,0.01))
    )
  )

model.totba_6 = inla(formula=formula_6_totba, family="gaussian",
                     data=d,
                     control.compute=list(dic=T, config=T)
)
summary(model.totba_6)

#Full model with centered annual population size for the parental hatch year and habitat type.
formula_totba_size_full2 = age1totba ~ z_FGRM + laflok_status*z_FGRM + c_Correct_size_hatchy*z_FGRM +
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

model.totba_full2 = inla(formula=formula_totba_size_full2, family="gaussian",
                                data=d_hatchyear,
                                control.compute=list(dic=T, config=T)
)
summary(model.totba_full2)


#Table the variances of the random effects

###Model 6 total badge size, "Final model".
inlaPosteriors_morph_4(model.totba_6)

###Full model with centered annual population size for the parental hatch year and habitat type.
inlaPosteriors_morph_4(model.totba_full2)


##############################
#Visible badge size, only males
##############################


#Model with FGRM:island and FGRM:island_year interactions, FGRM scaled, "Final model".
formula_6_visba = age1visba ~ z_FGRM +
  f(f_laflok,model="iid", param = c(0.1,0.01)) +
  f(f_laflok2, z_FGRM, model="iid",param = c(0.1,0.01)) +
  f(IAY, model = "iid",param = c(0.1,0.01)) +
  f(IAY2, z_FGRM, model = "iid", param = c(0.1,0.01)) +
  f(ID,model="generic0",
    Cmatrix=Cmatrix,
    hyper=list(
      prec=list( param = c(0.1,0.01))
    )
  )

model.visba_6 = inla(formula=formula_6_visba, family="gaussian",
                     data=d,
                     control.compute=list(dic=T, config=T)
)
summary(mode.visba_6)


#Full model with centered annual population size for the parental hatch year and habitat type.
formula_visba_size_full2 = age1visba ~ z_FGRM + laflok_status*z_FGRM + c_Correct_size_hatchy*z_FGRM +
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


model.visba_full2 = inla(formula=formula_visba_size_full2, family="gaussian",
                                data=d_hatchyear,
                                control.compute=list(dic=T, config=T)
)
summary(model.visba_full2)


#Table the variances of the random effects

###Model 6 visible badge size, "Final model".
inlaPosteriors_morph_4(model.visba_6)

###Full model with centered annual population size for the parental hatch year and habitat type.
inlaPosteriors_morph_4(model.visba_full2)


