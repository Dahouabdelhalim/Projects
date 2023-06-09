#Permutation analyses for estimating the probability of the observed variance in inbreeding (FGRM) explained by island and year nested within island
#Script originally by Yimen Araya-Ajoy 2018
#Modifications by Alina Niskanen in April 2020
#alina.niskanen@gmail.com

library(INLA)
library(MasterBayes)
library(MCMCglmm)

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

#Add presence of an island in 2013 for which no offspring data is available based on SNPs
data_adult_temp$X2013 <- NA
data_adult_temp$X2013[which(data_adult_temp$minadyear<=2013 & data_adult_temp$maxadyear >=2013)] <- 1

#Start choosing the correct data from "data_adult_temp" that includes all individuals that were observed as adults (laflok) on one of the 8 study islands
data_absence <- data_adult_temp[(is.na(data_adult_temp$X1998) & is.na(data_adult_temp$X1999) & is.na(data_adult_temp$X2000) & is.na(data_adult_temp$X2001) & is.na(data_adult_temp$X2002) & is.na(data_adult_temp$X2003) & is.na(data_adult_temp$X2004) & is.na(data_adult_temp$X2005) & is.na(data_adult_temp$X2006) & is.na(data_adult_temp$X2007) & is.na(data_adult_temp$X2008) & is.na(data_adult_temp$X2009) & is.na(data_adult_temp$X2010) & is.na(data_adult_temp$X2011) & is.na(data_adult_temp$X2012) & is.na(data_adult_temp$X2013)),]

#Make the final yearly estimate dataset, includes only individuals with laflok on one of the 8 study islands and adults on one of the 8 islands on 1998-2013
data_inbreeding <- data_adult_temp[!data_adult_temp$id %in% data_absence$id,]

# make island factor
data_inbreeding$f_laflok <- as.factor(data_inbreeding$laflok)

#Add combined variable for hatch year and adult island
data_inbreeding$IAY <- as.factor(paste(data_inbreeding$all_hatchyears,data_inbreeding$laflok,sep="_"))
#Add combined variable for hatch year and adult island
data_inbreeding$IAY2 <- data_inbreeding$IAY

#Scale FGRM
data_inbreeding$z_FGRM <- scale(data_inbreeding$FGRM, center = T, scale = T)


#Function for transforming precisions to variances
inlaPosteriors_inbr_0 <- function(model){
  sigma.island = inla.tmarginal(function(x) 1/x,  model$marginals.hyperpar$"Precision for f_laflok")
  m.island=inla.mmarginal(sigma.island)
  e.island=inla.emarginal(function(x) x, sigma.island)
  
  sigma.IAY = inla.tmarginal(function(x) 1/x,  model$marginals.hyperpar$"Precision for IAY")
  m.IAY=inla.mmarginal(sigma.IAY)
  e.IAY=inla.emarginal(function(x) x, sigma.IAY)
  
  sigma.e = inla.tmarginal(function(x) 1/x,  model$marginals.hyperpar$"Precision for the Gaussian observations")
  m.e =inla.mmarginal(sigma.e)
  e.e =inla.emarginal(function(x) x, sigma.e)
  
  results_tab <- cbind(rbind(VarE = e.e, island = e.island, Island_year=e.IAY),
                       rbind(VarE = m.e, island = m.island, Island_year=m.IAY),
                       1/model$summary.hyperpar[,c(5,3)])
  
  names(results_tab) <- c("Mean","Mode", "2.5%", "97.5 %")
  return(results_tab)
}


###Simulations for inbreeding differences between islands and island-years, including sex as fixed factor

n.sim=1000
n.ran=3

data_inbreeding$y=data_inbreeding$z_FGRM

formula_0.1_inbr = z_FGRM ~ gen_sex +
  f(f_laflok,model="iid",param = c(0.1,0.01)) +
  f(IAY, model = "iid",param = c(0.1,0.01))


model = inla(formula=formula_0.1_inbr, family="gaussian",
             data=data_inbreeding,
             control.compute=list(dic=T, config=T)
)

res.obs <- inlaPosteriors_inbr_0(model)


null.V.mean.inbreeding_1<-matrix(NA, n.sim, n.ran)
null.V2.mean.inbreeding_1<-matrix(NA, n.sim, n.ran)

null.V.mode.inbreeding_1<-matrix(NA, n.sim, n.ran)
null.V2.mode.inbreeding_1<-matrix(NA, n.sim, n.ran)



for(i in 1:n.sim){
  data_inbreeding2<-data_inbreeding
  R<-sample(1:nrow(data_inbreeding), nrow(data_inbreeding))
  data_inbreeding2$f_laflok<-data_inbreeding$f_laflok[R]
  data_inbreeding2$f_laflok2<-data_inbreeding2$f_laflok
  data_inbreeding2$IAY<-data_inbreeding$IAY[R]
  data_inbreeding2$IAY2<-data_inbreeding2$IAY
  model.null= inla(formula=formula_0.1_inbr, family="gaussian",
                   data=data_inbreeding2,
                   control.compute=list(dic=T, config=T))
  res<-inlaPosteriors_inbr_0(model.null)
  null.V2.mean.inbreeding_1[i,]<-as.numeric(0<(res[,1]-res.obs[,1]))
  null.V.mean.inbreeding_1[i,]<-res[,1]
  
  null.V2.mode.inbreeding_1[i,]<-as.numeric(0<(res[,2]-res.obs[,2]))
  null.V.mode.inbreeding_1[i,]<-res[,2]
}


colnames(null.V.mode.inbreeding_1)<-rownames(res.obs)
colnames(null.V.mean.inbreeding_1)<-rownames(res.obs)
colnames(null.V2.mode.inbreeding_1)<-rownames(res.obs)
colnames(null.V2.mean.inbreeding_1)<-rownames(res.obs)

apply(null.V2.mean.inbreeding_1,2,sum)/n.sim
apply(null.V2.mode.inbreeding_1,2,sum)/n.sim


#Print outputs for simulations
write.table(null.V.mean.inbreeding_1, file="null.V.mean.inbreeding_1.txt", row.names = F, col.names = T, quote = F)
write.table(null.V2.mean.inbreeding_1, file="null.V2.mean.inbreeding_1.txt", row.names = F, col.names = T, quote = F)
write.table(null.V.mode.inbreeding_1, file="null.V.mode.inbreeding_1.txt", row.names = F, col.names = T, quote = F)
write.table(null.V2.mode.inbreeding_1, file="null.V2.mode.inbreeding_1.txt", row.names = F, col.names = T, quote = F)

#Print outputs for mean and mode differences
sink("variance_probabilities_inbreeding_1.txt")
cat("Probabilities of mean and mode differences between observed data and 1000 simulations for FGRM, model including sex as fixed factor.")
cat("\\n")
cat("\\n")
apply(null.V2.mean.inbreeding_1,2,sum)/n.sim
apply(null.V2.mode.inbreeding_1,2,sum)/n.sim
cat("\\n")
cat("\\n")
sink()


    
    ###Simulations for inbreeding differences between islands and island-years, including sex and habitat type as fixed factors

n.sim=1000
n.ran=3

data_inbreeding$y=data_inbreeding$z_FGRM


formula_3.1_inbr = z_FGRM ~ gen_sex + laflok_status + 
  f(f_laflok,model="iid",param = c(0.1,0.01)) +
  f(IAY, model = "iid",param = c(0.1,0.01))


model = inla(formula=formula_3.1_inbr, family="gaussian",
                       data=data_inbreeding,
                       control.compute=list(dic=T, config=T)
)

res.obs <- inlaPosteriors_inbr_0(model)


null.V.mean.inbreeding<-matrix(NA, n.sim, n.ran)
null.V2.mean.inbreeding<-matrix(NA, n.sim, n.ran)

null.V.mode.inbreeding<-matrix(NA, n.sim, n.ran)
null.V2.mode.inbreeding<-matrix(NA, n.sim, n.ran)



for(i in 1:n.sim){
   data_inbreeding2<-data_inbreeding
   R<-sample(1:nrow(data_inbreeding), nrow(data_inbreeding))
   data_inbreeding2$f_laflok<-data_inbreeding$f_laflok[R]
   data_inbreeding2$f_laflok2<-data_inbreeding2$f_laflok
   data_inbreeding2$IAY<-data_inbreeding$IAY[R]
   data_inbreeding2$IAY2<-data_inbreeding2$IAY
   model.null= inla(formula=formula_3.1_inbr, family="gaussian",
                         data=data_inbreeding2,
                         control.compute=list(dic=T, config=T))
   res<-inlaPosteriors_inbr_0(model.null)
     null.V2.mean.inbreeding[i,]<-as.numeric(0<(res[,1]-res.obs[,1]))
     null.V.mean.inbreeding[i,]<-res[,1]

     null.V2.mode.inbreeding[i,]<-as.numeric(0<(res[,2]-res.obs[,2]))
     null.V.mode.inbreeding[i,]<-res[,2]
     }


colnames(null.V.mode.inbreeding)<-rownames(res.obs)
colnames(null.V.mean.inbreeding)<-rownames(res.obs)
colnames(null.V2.mode.inbreeding)<-rownames(res.obs)
colnames(null.V2.mean.inbreeding)<-rownames(res.obs)

apply(null.V2.mean.inbreeding,2,sum)/n.sim
apply(null.V2.mode.inbreeding,2,sum)/n.sim


#Print outputs for simulations
write.table(null.V.mean.inbreeding, file="null.V.mean.inbreeding_04052020.txt", row.names = F, col.names = T, quote = F)
write.table(null.V2.mean.inbreeding, file="null.V2.mean.inbreeding_04052020.txt", row.names = F, col.names = T, quote = F)
write.table(null.V.mode.inbreeding, file="null.V.mode.inbreeding_04052020.txt", row.names = F, col.names = T, quote = F)
write.table(null.V2.mode.inbreeding, file="null.V2.mode.inbreeding_04052020.txt", row.names = F, col.names = T, quote = F)

#Print outputs for mean and mode differences
sink("variance_probabilities_inbreeding_04052020.txt")
cat("Probabilities of mean and mode differences between observed data and 1000 simulations for FGRM, model including sex and habitat type as fixed factors.")
cat("\\n")
cat("\\n")
apply(null.V2.mean.inbreeding,2,sum)/n.sim
apply(null.V2.mode.inbreeding,2,sum)/n.sim
cat("\\n")
cat("\\n")
sink()

