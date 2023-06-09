#Install packages
install.packacges("ecospat")
install.packacges("dplyr")
install.packacges("adehabitatMA")
install.packacges("adehabitatHR")

library(ecospat)
library(dplyr)
library(adehabitatMA)
library(adehabitatHR)

###########################
#### data preparation #####
###########################

#Import csv with occurrences and background climatic information
bioclimdata_all<-read.csv("PathTo/bioclimdata_all.csv")
names(bioclimdata_all)
bioclimdata_all<-bioclimdata_all[complete.cases(bioclimdata_all),]

#occurence sites for the species (column names should be x,y)
occ<-bioclimdata_all[bioclimdata_all$obs_type == "occ",1:3]
bg<-bioclimdata_all[bioclimdata_all$obs_type == "background",1:3]

# species list
sp.list<-levels(occ[,1])
sp.nbocc<-c()
sp.nbg<-c()

for (i in 1:length(sp.list)){
  sp.nbocc<-c(sp.nbocc,length(which(occ[,1] == sp.list[i])))
  sp.nbg<-c(sp.nbg,length(which(bg[,1] == sp.list[i])))
} #calculate the nb of occurences per species
sp.list<-sp.list[sp.nbocc>4] # remove species with less than 5 occurences
nb.sp<-length(sp.list) #nb of species

################################
#### niche quantifications #####
################################

# selection of species to analyze
sp.choice<- c("amf","atf") #CHOOSE THE SET OF SPECIES FOR PAIRWISE ANALYSES

sp.combn<-combn(sp.choice,2) 

sp.combn<-combn(sp.list,2) 

x<-dim(sp.combn)[2]


# storage matrices
results <- matrix(nrow=x,ncol=4)
colnames(results) <- c("comb","D","N_Sim","Equiv")


# loop of niche quantification for each combination of species
t1=Sys.time()

i=1

for(i in 1:ncol(sp.combn)) { 
  
  spa<-sp.combn[1,i] #name of species a
  spb<-sp.combn[2,i] #name of species b
  
  results[i,"comb"]<-paste(spa,"vs",spb)
  
  occ.sp1<- bioclimdata_all[bioclimdata_all$sp == spa &
                              bioclimdata_all$obs_type == "occ",]
  
  occ.sp1<-occ.sp1[complete.cases(occ.sp1),]
  occ.sp1<-occ.sp1[,-c(1,2,3,11)]
  names(occ.sp1)
  occ.sp2<-bioclimdata_all[bioclimdata_all$sp == spb & 
                             bioclimdata_all$obs_type == "occ",]
  occ.sp2<-occ.sp2[complete.cases(occ.sp2),]
  occ.sp2<-occ.sp2[,-c(1,2,3,11)]
  
  clim1<- bioclimdata_all[bioclimdata_all$sp == spa & 
                            bioclimdata_all$obs_type == "background",]
  clim1<-clim1[complete.cases(clim1),]
  clim1<-clim1[,-c(1,2,3,11)]
  
  clim2<- bioclimdata_all[bioclimdata_all$sp == spb & 
                            bioclimdata_all$obs_type == "background",]
  clim2<-clim2[complete.cases(clim2),]
  clim2<-clim2[,-c(1,2,3,11)]
  
  clim12<-rbind(clim1,clim2)

  ## ANALYSIS - selection of parameters
  
  # selection of the type of analysis.
  # If PROJ =F, the models are calibrated on both ranges.
  # If PROJ =T, the models are calibrated on species 1 range only and projected to range 2.
  PROJ = F
  #number of interation for the tests of equivalency and similarity
  iterations<-100
  #resolution of the gridding of the climate space
  R=100
  
  # selection of variables to include in the analyses
  length(names(clim12))
  Xvar<-c(1:11)
  nvar<-length(Xvar)
  
  ################################
  #### niche quantifications #####
  ################################
  
  #calibration of PCA-env 
  pca.env <-dudi.pca(clim12, center = T, scale = T, scannf = F, nf = 2)
  
  # predict the scores on the PCA axes
  scores.clim12<- pca.env$li
  scores.clim1<- suprow(pca.env,clim1)$lisup
  scores.clim2<- suprow(pca.env,clim2)$lisup
  scores.occ.sp1<- suprow(pca.env,occ.sp1)$lisup
  scores.occ.sp2<- suprow(pca.env,occ.sp2)$lisup

  # calculation of occurence density
  za<- ecospat.grid.clim.dyn(scores.clim12,scores.clim1,scores.occ.sp1,100) # Back especifico para cada sp
  zb<- ecospat.grid.clim.dyn(scores.clim12,scores.clim2,scores.occ.sp2,100) # Back especifico para cada sp
  
  # overlap corrected by availabilty of background conditions
  x<-ecospat.niche.overlap(za,zb,cor=T) 
  results[i,"D"]<- round(x$D,3)
  
  
  # test of niche equivalency and similarity
  equ<-ecospat.niche.equivalency.test(za,zb,rep=100, alternative = "lower") #put at least 100 for final analyses
  results[i,"Equiv"]<-round(equ$p.D,3)
  
  nsim<-ecospat.niche.similarity.test(za,zb,rep=1000,alternative = "lower",
                                    rand.type = 1)
  
  ecospat.plot.overlap.test(nsim,"D","Niche Conservatims")
  results[i,"N_Sim"]<-round(nsim$p.D,3)
  
  #Figures
  pdf(file=paste('08ENMeval/broennimann/Figuras/',
                 paste(spa,"_",spb),".pdf", sep=""),width=14,height=8)
  
  layout(matrix(c(1,2,3,4,5,6), nrow=3, byrow=TRUE))
  ecospat.plot.contrib(pca.env$co,pca.env$eig)
  ecospat.plot.niche (za, title=spa, name.axis1="PC1", name.axis2="PC2", cor=F)
  ecospat.plot.niche (zb, title=spb, name.axis1="PC1", name.axis2="PC2", cor=F)
  ecospat.plot.overlap.test(equ,"D","Equivalency")
  ecospat.plot.overlap.test(nsim,"D","Niche Similarity Test")
  dev.off()
    
  #counter
  cat(paste(i))
}

t2=Sys.time()

t2-t1 #Time difference 
t1=10.23

#write results
write.csv(results,"08ENMeval/broennimann/bronn_res_eq_lower.csv",row.names = F)

