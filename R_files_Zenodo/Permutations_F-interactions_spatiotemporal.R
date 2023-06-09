#Permutation analyses for estimating the probability of the observed variance in inbreeding depression in different house sparrow traits explained by the interaction between FGRM and island or year nested within island
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

#Only those individuals that have age1mass estimate and have been adults on one of the study islands
data_morph <- data_adult_temp[!is.na(data_adult_temp$age1mass),]

#Scale FGRM
data_morph$z_FGRM <- scale(data_morph$FGRM, center = T, scale = T)

# make island factor covariate
data_morph$f_laflok <- as.factor(data_morph$laflok)
data_morph$f_laflok2 <- data_morph$f_laflok

#Add combined variable for hatch year and island
data_morph$IAY <- as.factor(paste(data_morph$all_hatchyears,data_morph$laflok,sep="_"))
data_morph$IAY2 <- data_morph$IAY

#Make a separate dataset for only males to use in badge size tests
d <- data_morph[which(data_morph$gen_sex=="m"),]


#Function for transforming precisions to variances
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


###Simulations for mass

n.sim=1000
n.ran=6

data_morph$y=data_morph$age1mass

formula = y ~ z_FGRM + gen_sex +
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

model = inla(formula=formula, family="gaussian",
                       data=data_morph,
                       control.compute=list(dic=T, config=T)
)


res.obs<-inlaPosteriors_morph_4(model)

null.V.mean.mass<-matrix(NA, n.sim, n.ran)
null.V2.mean.mass<-matrix(NA, n.sim, n.ran)

null.V.mode.mass<-matrix(NA, n.sim, n.ran)
null.V2.mode.mass<-matrix(NA, n.sim, n.ran)



for(i in 1:n.sim){
   data_morph2<-data_morph
   R<-sample(1:nrow(data_morph), nrow(data_morph))
   data_morph2$f_laflok<-data_morph$f_laflok[R]
   data_morph2$f_laflok2<-data_morph2$f_laflok
   data_morph2$IAY<-data_morph$IAY[R]
   data_morph2$IAY2<-data_morph2$IAY
   model.null= inla(formula=formula, family="gaussian",
                         data=data_morph2,
                         control.compute=list(dic=T, config=T))

   res<-inlaPosteriors_morph_4(model.null)
     null.V2.mean.mass[i,]<-as.numeric(0<(res[,1]-res.obs[,1]))
     null.V.mean.mass[i,]<-res[,1]

     null.V2.mode.mass[i,]<-as.numeric(0<(res[,2]-res.obs[,2]))
     null.V.mode.mass[i,]<-res[,2]
     }


colnames(null.V.mode.mass)<-rownames(res.obs)
colnames(null.V.mean.mass)<-rownames(res.obs)
colnames(null.V2.mode.mass)<-rownames(res.obs)
colnames(null.V2.mean.mass)<-rownames(res.obs)

apply(null.V2.mean.mass,2,sum)/n.sim
apply(null.V2.mode.mass,2,sum)/n.sim


###Simulations for beak length

n.sim=1000
n.ran=6

data_morph$y=data_morph$age1billL

formula = y ~ z_FGRM + gen_sex +
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

model = inla(formula=formula, family="gaussian",
                       data=data_morph,
                       control.compute=list(dic=T, config=T)
)


res.obs<-inlaPosteriors_morph_4(model)

null.V.mean.billL<-matrix(NA, n.sim, n.ran)
null.V2.mean.billL<-matrix(NA, n.sim, n.ran)

null.V.mode.billL<-matrix(NA, n.sim, n.ran)
null.V2.mode.billL<-matrix(NA, n.sim, n.ran)



for(i in 1:n.sim){
   data_morph2<-data_morph
   R<-sample(1:nrow(data_morph), nrow(data_morph))
   data_morph2$f_laflok<-data_morph$f_laflok[R]
   data_morph2$f_laflok2<-data_morph2$f_laflok
   data_morph2$IAY<-data_morph$IAY[R]
   data_morph2$IAY2<-data_morph2$IAY
   model.null= inla(formula=formula, family="gaussian",
                         data=data_morph2,
                         control.compute=list(dic=T, config=T))

   res<-inlaPosteriors_morph_4(model.null)
     null.V2.mean.billL[i,]<-as.numeric(0<(res[,1]-res.obs[,1]))
     null.V.mean.billL[i,]<-res[,1]

     null.V2.mode.billL[i,]<-as.numeric(0<(res[,2]-res.obs[,2]))
     null.V.mode.billL[i,]<-res[,2]
     }


colnames(null.V.mode.billL)<-rownames(res.obs)
colnames(null.V.mean.billL)<-rownames(res.obs)
colnames(null.V2.mode.billL)<-rownames(res.obs)
colnames(null.V2.mean.billL)<-rownames(res.obs)

apply(null.V2.mean.billL,2,sum)/n.sim
apply(null.V2.mode.billL,2,sum)/n.sim


###Simulations for beak depth

n.sim=1000
n.ran=6

data_morph$y=data_morph$age1billD


formula = y ~ z_FGRM + gen_sex +
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

model = inla(formula=formula, family="gaussian",
                       data=data_morph,
                       control.compute=list(dic=T, config=T)
)


res.obs<-inlaPosteriors_morph_4(model)

null.V.mean.billD<-matrix(NA, n.sim, n.ran)
null.V2.mean.billD<-matrix(NA, n.sim, n.ran)

null.V.mode.billD<-matrix(NA, n.sim, n.ran)
null.V2.mode.billD<-matrix(NA, n.sim, n.ran)



for(i in 1:n.sim){
   data_morph2<-data_morph
   R<-sample(1:nrow(data_morph), nrow(data_morph))
   data_morph2$f_laflok<-data_morph$f_laflok[R]
   data_morph2$f_laflok2<-data_morph2$f_laflok
   data_morph2$IAY<-data_morph$IAY[R]
   data_morph2$IAY2<-data_morph2$IAY
   model.null= inla(formula=formula, family="gaussian",
                         data=data_morph2,
                         control.compute=list(dic=T, config=T))

   res<-inlaPosteriors_morph_4(model.null)
     null.V2.mean.billD[i,]<-as.numeric(0<(res[,1]-res.obs[,1]))
     null.V.mean.billD[i,]<-res[,1]

     null.V2.mode.billD[i,]<-as.numeric(0<(res[,2]-res.obs[,2]))
     null.V.mode.billD[i,]<-res[,2]
     }


colnames(null.V.mode.billD)<-rownames(res.obs)
colnames(null.V.mean.billD)<-rownames(res.obs)
colnames(null.V2.mode.billD)<-rownames(res.obs)
colnames(null.V2.mean.billD)<-rownames(res.obs)

apply(null.V2.mean.billD,2,sum)/n.sim
apply(null.V2.mode.billD,2,sum)/n.sim


###Simulations for wing length 

n.sim=1000
n.ran=6

data_morph$y=data_morph$age1wing


formula = y ~ z_FGRM + gen_sex +
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

model_1 = inla(formula=formula, family="gaussian",
             data=data_morph,
             control.compute=list(dic=T, config=T)
)

res.obs<-inlaPosteriors_morph_4(model_1)

null.V.mean.wing_2<-matrix(NA, n.sim, n.ran)
null.V2.mean.wing_2<-matrix(NA, n.sim, n.ran)

null.V.mode.wing_2<-matrix(NA, n.sim, n.ran)
null.V2.mode.wing_2<-matrix(NA, n.sim, n.ran)


for(i in 1:n.sim){
   data_morph2<-data_morph
   R<-sample(1:nrow(data_morph), nrow(data_morph))
   data_morph2$f_laflok<-data_morph$f_laflok[R]
   data_morph2$f_laflok2<-data_morph2$f_laflok
   data_morph2$IAY<-data_morph$IAY[R]
   data_morph2$IAY2<-data_morph2$IAY
   model.null= inla(formula=formula, family="gaussian",
                         data=data_morph2,
                         control.compute=list(dic=T, config=T))

   res<-inlaPosteriors_morph_4(model.null)
     null.V2.mean.wing_2[i,]<-as.numeric(0<(res[,1]-res.obs[,1]))
     null.V.mean.wing_2[i,]<-res[,1]
     
     null.V2.mode.wing_2[i,]<-as.numeric(0<(res[,2]-res.obs[,2]))
     null.V.mode.wing_2[i,]<-res[,2]
     }


colnames(null.V.mode.wing_2)<-rownames(res.obs)
colnames(null.V.mean.wing_2)<-rownames(res.obs)
colnames(null.V2.mode.wing_2)<-rownames(res.obs)
colnames(null.V2.mean.wing_2)<-rownames(res.obs)

apply(null.V2.mean.wing_2,2,sum)/n.sim
apply(null.V2.mode.wing_2,2,sum)/n.sim


###Simulations for tarsus length

n.sim=1000
n.ran=6

data_morph$y=data_morph$age1tarsus

data_morph$ob<-1:nrow(data_morph)

formula = y ~ z_FGRM + gen_sex +
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

model = inla(formula=formula, family="gaussian",
                       data=data_morph,
                       control.compute=list(dic=T, config=T)
)


res.obs<-inlaPosteriors_morph_4(model)

null.V.mean.tarsus<-matrix(NA, n.sim, n.ran)
null.V2.mean.tarsus<-matrix(NA, n.sim, n.ran)

null.V.mode.tarsus<-matrix(NA, n.sim, n.ran)
null.V2.mode.tarsus<-matrix(NA, n.sim, n.ran)



for(i in 1:n.sim){
   data_morph2<-data_morph
   R<-sample(1:nrow(data_morph), nrow(data_morph))
   data_morph2$f_laflok<-data_morph$f_laflok[R]
   data_morph2$f_laflok2<-data_morph2$f_laflok
   data_morph2$IAY<-data_morph$IAY[R]
   data_morph2$IAY2<-data_morph2$IAY
   model.null= inla(formula=formula, family="gaussian",
                         data=data_morph2,
                         control.compute=list(dic=T, config=T))

   res<-inlaPosteriors_morph_4(model.null)
     null.V2.mean.tarsus[i,]<-as.numeric(0<(res[,1]-res.obs[,1]))
     null.V.mean.tarsus[i,]<-res[,1]

     null.V2.mode.tarsus[i,]<-as.numeric(0<(res[,2]-res.obs[,2]))
     null.V.mode.tarsus[i,]<-res[,2]
     }


colnames(null.V.mode.tarsus)<-rownames(res.obs)
colnames(null.V.mean.tarsus)<-rownames(res.obs)
colnames(null.V2.mode.tarsus)<-rownames(res.obs)
colnames(null.V2.mean.tarsus)<-rownames(res.obs)

apply(null.V2.mean.tarsus,2,sum)/n.sim
apply(null.V2.mode.tarsus,2,sum)/n.sim


###Simulations for total badge size
#Only males

n.sim=1000
n.ran=6

d$y=d$age1totba

d$ob<-1:nrow(d)

formula = y ~ z_FGRM +
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

model = inla(formula=formula, family="gaussian",
             data=d,
             control.compute=list(dic=T, config=T)
)


res.obs<-inlaPosteriors_morph_4(model)

null.V.mean.totba<-matrix(NA, n.sim, n.ran)
null.V2.mean.totba<-matrix(NA, n.sim, n.ran)

null.V.mode.totba<-matrix(NA, n.sim, n.ran)
null.V2.mode.totba<-matrix(NA, n.sim, n.ran)



for(i in 1:n.sim){
  d2<-d
  R<-sample(1:nrow(d), nrow(d))
  d2$f_laflok<-d$f_laflok[R]
  d2$f_laflok2<-d2$f_laflok
  d2$IAY<-d$IAY[R]
  d2$IAY2<-d2$IAY
  model.null= inla(formula=formula, family="gaussian",
                   data=d2,
                   control.compute=list(dic=T, config=T))
  
  res<-inlaPosteriors_morph_4(model.null)
  null.V2.mean.totba[i,]<-as.numeric(0<(res[,1]-res.obs[,1]))
  null.V.mean.totba[i,]<-res[,1]
  
  null.V2.mode.totba[i,]<-as.numeric(0<(res[,2]-res.obs[,2]))
  null.V.mode.totba[i,]<-res[,2]
}


colnames(null.V.mode.totba)<-rownames(res.obs)
colnames(null.V.mean.totba)<-rownames(res.obs)
colnames(null.V2.mode.totba)<-rownames(res.obs)
colnames(null.V2.mean.totba)<-rownames(res.obs)

apply(null.V2.mean.totba,2,sum)/n.sim
apply(null.V2.mode.totba,2,sum)/n.sim


###Simulations for visible badge size
#Only males

n.sim=1000
n.ran=6

d$y=d$age1visba

d$ob<-1:nrow(d)

formula = y ~ z_FGRM +
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

model = inla(formula=formula, family="gaussian",
             data=d,
             control.compute=list(dic=T, config=T)
)


res.obs<-inlaPosteriors_morph_4(model)

null.V.mean.visba<-matrix(NA, n.sim, n.ran)
null.V2.mean.visba<-matrix(NA, n.sim, n.ran)

null.V.mode.visba<-matrix(NA, n.sim, n.ran)
null.V2.mode.visba<-matrix(NA, n.sim, n.ran)



for(i in 1:n.sim){
  d2<-d
  R<-sample(1:nrow(d), nrow(d))
  d2$f_laflok<-d$f_laflok[R]
  d2$f_laflok2<-d2$f_laflok
  d2$IAY<-d$IAY[R]
  d2$IAY2<-d2$IAY
  model.null= inla(formula=formula, family="gaussian",
                   data=d2,
                   control.compute=list(dic=T, config=T))
  
  res<-inlaPosteriors_morph_4(model.null)
  null.V2.mean.visba[i,]<-as.numeric(0<(res[,1]-res.obs[,1]))
  null.V.mean.visba[i,]<-res[,1]
  
  null.V2.mode.visba[i,]<-as.numeric(0<(res[,2]-res.obs[,2]))
  null.V.mode.visba[i,]<-res[,2]
}


colnames(null.V.mode.visba)<-rownames(res.obs)
colnames(null.V.mean.visba)<-rownames(res.obs)
colnames(null.V2.mode.visba)<-rownames(res.obs)
colnames(null.V2.mean.visba)<-rownames(res.obs)

apply(null.V2.mean.visba,2,sum)/n.sim
apply(null.V2.mode.visba,2,sum)/n.sim


#Print outputs for simulations
write.table(null.V.mean.mass, file="null.V.mean.mass.txt", row.names = F, col.names = T, quote = F)
write.table(null.V2.mean.mass, file="null.V2.mean.mass.txt", row.names = F, col.names = T, quote = F)
write.table(null.V.mode.mass, file="null.V.mode.mass.txt", row.names = F, col.names = T, quote = F)
write.table(null.V2.mode.mass, file="null.V2.mode.mass.txt", row.names = F, col.names = T, quote = F)

write.table(null.V.mean.billL, file="null.V.mean.billL.txt", row.names = F, col.names = T, quote = F)
write.table(null.V2.mean.billL, file="null.V2.mean.billL.txt", row.names = F, col.names = T, quote = F)
write.table(null.V.mode.billL, file="null.V.mode.billL.txt", row.names = F, col.names = T, quote = F)
write.table(null.V2.mode.billL, file="null.V2.mode.billL.txt", row.names = F, col.names = T, quote = F)

write.table(null.V.mean.billD, file="null.V.mean.billD.txt", row.names = F, col.names = T, quote = F)
write.table(null.V2.mean.billD, file="null.V2.mean.billD.txt", row.names = F, col.names = T, quote = F)
write.table(null.V.mode.billD, file="null.V.mode.billD.txt", row.names = F, col.names = T, quote = F)
write.table(null.V2.mode.billD, file="null.V2.mode.billD.txt", row.names = F, col.names = T, quote = F)

write.table(null.V.mean.wing_2, file="null.V.mean.wing.txt", row.names = F, col.names = T, quote = F)
write.table(null.V2.mean.wing_2, file="null.V2.mean.wing.txt", row.names = F, col.names = T, quote = F)
write.table(null.V.mode.wing_2, file="null.V.mode.wing.txt", row.names = F, col.names = T, quote = F)
write.table(null.V2.mode.wing_2, file="null.V2.mode.wing.txt", row.names = F, col.names = T, quote = F)

write.table(null.V.mean.tarsus, file="null.V.mean.tarsus.txt", row.names = F, col.names = T, quote = F)
write.table(null.V2.mean.tarsus, file="null.V2.mean.tarsus.txt", row.names = F, col.names = T, quote = F)
write.table(null.V.mode.tarsus, file="null.V.mode.tarsus.txt", row.names = F, col.names = T, quote = F)
write.table(null.V2.mode.tarsus, file="null.V2.mode.tarsus.txt", row.names = F, col.names = T, quote = F)

write.table(null.V.mean.totba, file="null.V.mean.totba.txt", row.names = F, col.names = T, quote = F)
write.table(null.V2.mean.totba, file="null.V2.mean.totba.txt", row.names = F, col.names = T, quote = F)
write.table(null.V.mode.totba, file="null.V.mode.totba.txt", row.names = F, col.names = T, quote = F)
write.table(null.V2.mode.totba, file="null.V2.mode.totba.txt", row.names = F, col.names = T, quote = F)

write.table(null.V.mean.visba, file="null.V.mean.visba.txt", row.names = F, col.names = T, quote = F)
write.table(null.V2.mean.visba, file="null.V2.mean.visba.txt", row.names = F, col.names = T, quote = F)
write.table(null.V.mode.visba, file="null.V.mode.visba.txt", row.names = F, col.names = T, quote = F)
write.table(null.V2.mode.visba, file="null.V2.mode.visba.txt", row.names = F, col.names = T, quote = F)


#Print outputs for mean and mode differences
sink("variance_probabilities_morphology.txt")
cat("Probabilities of mean and mode differences between observed data and 1000 simulations.")
cat("\\n")
cat("\\n")
cat("Mass", sep="\\n")
cat("\\n")
apply(null.V2.mean.mass,2,sum)/n.sim
apply(null.V2.mode.mass,2,sum)/n.sim
cat("\\n")
cat("\\n")
cat("Bill length", sep="\\n")
cat(sep="/")
apply(null.V2.mean.billL,2,sum)/n.sim
apply(null.V2.mode.billL,2,sum)/n.sim
cat("\\n")
cat("\\n")
cat("Bill depth", sep="\\n")
cat("\\n")
apply(null.V2.mean.billD,2,sum)/n.sim
apply(null.V2.mode.billD,2,sum)/n.sim
cat("\\n")
cat("\\n")
cat("Wing length", sep="\\n")
cat("\\n")
apply(null.V2.mean.wing_2,2,sum)/n.sim
apply(null.V2.mode.wing_2,2,sum)/n.sim
cat("\\n")
cat("\\n")
cat("Tarsus length", sep="\\n")
cat("\\n")
apply(null.V2.mean.tarsus,2,sum)/n.sim
apply(null.V2.mode.tarsus,2,sum)/n.sim
cat("\\n")
cat("\\n")
cat("Total badge size", sep="\\n")
apply(null.V2.mean.totba,2,sum)/n.sim
apply(null.V2.mode.totba,2,sum)/n.sim
cat("\\n")
cat("\\n")
cat("Visual badge size", sep="\\n")
cat("\\n")
apply(null.V2.mean.visba,2,sum)/n.sim
apply(null.V2.mode.visba,2,sum)/n.sim
cat("\\n")
sink()

