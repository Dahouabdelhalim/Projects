#Permutation analyses for estimating the probability of the observed variance in inbreeding depression in different house sparrow traits explained by the interaction between FGRM and island or year nested within island
#Script originally by Yimen Araya-Ajoy 2018
#Modifications by Alina Niskanen in April 2020
#alina.niskanen@gmail.com

library(INLA)
library(MasterBayes)
library(MCMCglmm)
library(tidyr)

#Import the phenotypic data and inbreeding estimates
data_LRS <- read.csv("Data.txt", header = T, stringsAsFactors = F, sep = "\\t")

#Import the pedigree
ped <- read.table("pedigree.txt", header = T, stringsAsFactors = F)

## need an ordered Pedigree for the inverseA() function:
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
data_LRS$ID <- d.map[match(data_LRS$id, d.map$id), "ID"]
data_LRS$ID2 <- data_LRS$ID

###
###Subsetting the data
###

#Include only the individuals that have been observed as adults on one of the 8 study islands
data_adult_temp <- data_LRS[which(data_LRS$laflok %in% c("20","22","23","24","26","27","28","38")),]
table(data_adult_temp$laflok)


##########################
###Lifetime recruit production
##########################

#Only those individuals that have LRS_OK estimate and have been adults on one of the study islands
data_ok <- data_adult_temp[!is.na(data_adult_temp$LRS_OK),]

#Scale FGRM
data_ok$z_FGRM <- scale(data_ok$FGRM, center = T, scale = T)

# make island factor
data_ok$f_laflok <- as.factor(data_ok$laflok)
data_ok$f_laflok2 <- data_ok$f_laflok

#Add combined variable for hatch year and adult island
data_ok$IAY <- as.factor(paste(data_ok$all_hatchyears,data_ok$laflok,sep="_"))
data_ok$IAY2 <- data_ok$IAY


##########################
###Per year recruit production
##########################

#Start choosing the correct data from "data_adult_temp" that includes all individuals that were observed as adults (laflok) on one of the 8 study islands
data_no_YR <- data_adult_temp[(is.na(data_adult_temp$X1998) & is.na(data_adult_temp$X1999) & is.na(data_adult_temp$X2000) & is.na(data_adult_temp$X2001) & is.na(data_adult_temp$X2002) & is.na(data_adult_temp$X2003) & is.na(data_adult_temp$X2004) & is.na(data_adult_temp$X2005) & is.na(data_adult_temp$X2006) & is.na(data_adult_temp$X2007) & is.na(data_adult_temp$X2008) & is.na(data_adult_temp$X2009) & is.na(data_adult_temp$X2010) & is.na(data_adult_temp$X2011) & is.na(data_adult_temp$X2012)),]

#Make the final yearly estimate dataset, includes only individuals with laflok on one of the 8 study islands and YR estimate available
data_YR_temp <- data_adult_temp[!data_adult_temp$id %in% data_no_YR$id,]

#Scale FGRM
data_YR_temp$z_FGRM <- scale(data_YR_temp$FGRM, center = T, scale = T)

#Transpose the data to include each individual yearly estimate on its own row
data_YR <- gather(data_YR_temp, "Year", "recruits", 21:35)
data_YR$Year <- gsub("X","", data_YR$Year)

#Remove all rows with missing YR to reduce the file size
data_YR <- data_YR[!is.na(data_YR$recruits),]
#Add age of the individual when having the offspring
data_YR$Year <- as.integer(data_YR$Year)
data_YR$age <- data_YR$Year-data_YR$all_hatchyears

# make island factor
data_YR$f_laflok <- as.factor(data_YR$laflok)
data_YR$f_laflok2 <- data_YR$laflok
data_YR$f_age <- as.factor(data_YR$age)
#add squared age for quadratic effect
data_YR$age2 <- data_YR$age^2
#Center age
data_YR$c_age <- scale(data_YR$age, center = T, scale = F)
data_YR$c_age2 <- data_YR$c_age^2
#Combine age classes 6 and above
data_YR$age_comb <- data_YR$age
data_YR$age_comb[which(data_YR$age>=6)] <- 6
#Center combined age
data_YR$c_age_comb <- scale(data_YR$age_comb, center = T, scale = F)
data_YR$c_age_comb2 <- data_YR$c_age_comb^2
#Add combined variable for hatch year and adult island
data_YR$IAY <- as.factor(paste(data_YR$all_hatchyears,data_YR$laflok,sep="_"))
data_YR$IAY2 <- data_YR$IAY


###Function to transform precision to variance in lifetime-reproductive success model
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


###Function to transform precision to variance in annual reproductive success model
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


##########################
###Simulations for lifetime recruit production
##########################

#Some of the permutation rounds introduce errors, so the number of permutations is 1500 and only 1000 of the successful rounds will be used to estimate the probability of the observed variance

n.sim=100 #1500
n.ran=6

data_ok$y=data_ok$LRS_OK

formula = y ~ z_FGRM + gen_sex +  
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

model = inla(formula=formula, family="zeroinflatedpoisson1",
                   data=data_ok,
                   control.compute=list(dic=T, config=T)
)

res.obs<-inlaPosteriors_Zip_4(model)

null.V.mean.LRS<-matrix(NA, n.sim, n.ran)
null.V2.mean.LRS<-matrix(NA, n.sim, n.ran)

null.V.mode.LRS<-matrix(NA, n.sim, n.ran)
null.V2.mode.LRS<-matrix(NA, n.sim, n.ran)
error<-rep(NA, n.sim)


#More warnings about different Hessian Eigenvalues being negative, so changed tolerance from 1e-3 to 1e-6, but either of those values helped the model to converge, so went back to the default value
for(i in 1:n.sim){
   data_ok2<-data_ok
   R<-sample(1:nrow(data_ok), nrow(data_ok))
   data_ok2$f_laflok<-data_ok$f_laflok[R]
   data_ok2$f_laflok2<-data_ok2$f_laflok
   data_ok2$IAY<-data_ok$IAY[R]
   data_ok2$IAY2<-data_ok2$IAY
   model.null= inla(formula=formula, family="zeroinflatedpoisson1",
                         data=data_ok2,
                         control.compute=list(dic=T, config=T))

   res<-inlaPosteriors_Zip_4(model.null)
     null.V2.mean.LRS[i,]<-as.numeric(0<(res[,1]-res.obs[,1]))
     null.V.mean.LRS[i,]<-res[,1]
     error[i]<-model.null$mode$mode.status
     null.V2.mode.LRS[i,]<-as.numeric(0<(res[,2]-res.obs[,2]))
     null.V.mode.LRS[i,]<-res[,2]
     }

colnames(null.V.mode.LRS)<-rownames(res.obs)
colnames(null.V.mean.LRS)<-rownames(res.obs)
colnames(null.V2.mode.LRS)<-rownames(res.obs)
colnames(null.V2.mean.LRS)<-rownames(res.obs)

#Include only successful simulation rounds
error2<-error==0
null.V.mode.LRS2<-null.V.mode.LRS[error2,]
null.V.mean.LRS2<-null.V.mean.LRS[error2,]
null.V2.mode.LRS2<-null.V2.mode.LRS[error2,]
null.V2.mean.LRS2<-null.V2.mean.LRS[error2,]

#Print outputs for successful simulations
write.table(null.V.mode.LRS2, file="null.V.mode.LRS2.txt", row.names = F, col.names = T, quote = F)
write.table(null.V.mean.LRS2, file="null.V.mean.LRS2.txt", row.names = F, col.names = T, quote = F)
write.table(null.V2.mode.LRS2, file="null.V2.mode.LRS2.txt", row.names = F, col.names = T, quote = F)
write.table(null.V2.mean.LRS2, file="null.V2.mean.LRS2.txt", row.names = F, col.names = T, quote = F)

#Print outputs for mean and mode differences between observed and simulated results
sink("Variance_probabilities_LRS.txt")
cat("Probabilities of mean and mode differences between observed data and 1000 successful simulations.")
cat("\\n")
apply(null.V2.mean.LRS2[1:1000,],2,sum)/1000
apply(null.V2.mode.LRS2[1:1000,],2,sum)/1000
sink()


##########################
###Simulations for per year recruit production
##########################

#Some of the permutation rounds introduce errors, so the number of permutations is 1500 and only 1000 of the successful rounds will be used to estimate the probability of the observed variance

n.sim=1500
n.ran=6

data_YR$y=data_YR$recruits

formula = y ~ z_FGRM + gen_sex + c_age_comb + c_age_comb2 +
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

model = inla(formula=formula, family="poisson",
                          data=data_YR,
                          control.compute=list(dic=T, config=T))
   
#Transform precision to variance
res.obs<-inlaPosteriors_Poi_3(model)   


null.V.mean.YR<-matrix(NA, n.sim, n.ran)
null.V2.mean.YR<-matrix(NA, n.sim, n.ran)

null.V.mode.YR<-matrix(NA, n.sim, n.ran)
null.V2.mode.YR<-matrix(NA, n.sim, n.ran)
error<-rep(NA, n.sim)


for(i in 1:n.sim){
   data_YR2<-data_YR
   R<-sample(1:nrow(data_YR), nrow(data_YR))
   data_YR2$f_laflok<-data_YR$f_laflok[R]
   data_YR2$f_laflok2<-data_YR2$f_laflok
   data_YR2$IAY<-data_YR$IAY[R]
   data_YR2$IAY2<-data_YR2$IAY
   model.null= inla(formula=formula, family="poisson",
                         data=data_YR2,
                         control.compute=list(dic=T, config=T))
   res<-inlaPosteriors_Poi_3(model.null)
     null.V2.mean.YR[i,]<-as.numeric(0<(res[,1]-res.obs[,1]))
     null.V.mean.YR[i,]<-res[,1]
     error[i]<-model.null$mode$mode.status
     null.V2.mode.YR[i,]<-as.numeric(0<(res[,2]-res.obs[,2]))
     null.V.mode.YR[i,]<-res[,2]
     }



colnames(null.V.mode.YR)<-rownames(res.obs)
colnames(null.V.mean.YR)<-rownames(res.obs)
colnames(null.V2.mode.YR)<-rownames(res.obs)
colnames(null.V2.mean.YR)<-rownames(res.obs)

error2<-error==0
null.V.mode.YR2<-null.V.mode.YR[error2,]
null.V.mean.YR2<-null.V.mean.YR[error2,]
null.V2.mode.YR2<-null.V2.mode.YR[error2,]
null.V2.mean.YR2<-null.V2.mean.YR[error2,]


#Print outputs for successful simulations

write.table(null.V.mode.YR2, file="null.V.mode.YR2.txt", row.names = F, col.names = T, quote = F)
write.table(null.V.mean.YR2, file="null.V.mean.YR2.txt", row.names = F, col.names = T, quote = F)
write.table(null.V2.mode.YR2, file="null.V2.mode.YR2.txt", row.names = F, col.names = T, quote = F)
write.table(null.V2.mean.YR2, file="null.V2.mean.YR2.txt", row.names = F, col.names = T, quote = F)

#Print outputs for mean and mode differences
sink("Variance_probabilities_AR.txt")
cat("Probabilities of mean and mode differences between observed data and 1000 successful simulations.")
cat("\\n")
cat("\\n")
apply(null.V2.mean.YR2[1:1000,],2,sum)/1000
apply(null.V2.mode.YR2[1:1000,],2,sum)/1000
cat("\\n")
cat("\\n")
sink()

