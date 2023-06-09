require(OutFLANK)
require(vcfR)
require(adegenet)
require(gdm)
require(gradientForest)
require(foreach)
require(doParallel)
require(pbapply)
require(gdata)
require(data.table)
require(PresenceAbsence)
require(ROCR)
require(modEvA)
require(dplyr)
require(grid)
require(gridExtra)
require(gtools)
require(stringr)
require(reshape2)
require(hierfstat)
require(clusterGeneration)
#require(ggplot2)

#setwd("/Users/akijarl/Desktop/PostDoc/TTT_Offset_Vulnerability_GF_Sims/")
#setwd("/Users/akijarl/Desktop/TTT_Offset_Vulnerability_GF_Sims/")
#setwd("E:/Research_AJL/TTT_Offset_Vulnerability_GF_Sims/")
#load("results/R_results/1602454886424.RData")

#f1<-list.files("results/SLiM_output/Sim_sum/",pattern=".txt")[11]
# (seed<-substr(f1, start=0, stop=13))

#details=file.info(list.files("results/SLiM_output/Sim_sum/",full.names = T))
#(f1<-row.names(details[with(details,order(as.POSIXct(mtime),decreasing = T)),])[1])
#(seed<-substr(f1, start=30, stop=42))

seed=commandArgs(trailingOnly=T)
fit<-read.table(paste("results/SLiM_output/Sim_sum/",seed,"_Freq_SL_WF_Antag.txt",sep=""))
#fit<-read.table(paste("results/SLiM_output/Sim_sum/",f1,sep=""))
#fit<-read.table(f1)

fit_nam <- NULL
for(i in 1:100){
  fit_nam <- c(fit_nam,paste("P",i,"_fit",sep=""))
}

freq_nam <- NULL
for(i in 1:100){
  freq_nam <- c(freq_nam,paste("P",i,"_freq",sep=""))
}

env_nam <- NULL
for(i in 1:100){
  env_nam <- c(env_nam,paste("P",i,"_env",sep=""))
}

colnames(fit)<-c("m","n","u","r","Env_rate","Burnin","Env_shift", "Generation", fit_nam, freq_nam, env_nam)

N<-data.frame(fit$n[1]*100)
colnames(N)<-"N"
(specs<-data.frame(seed,fit[1,1:2],N,fit[1,3:7]))

#write.table(specs,"R_results/output_metadata.txt",append=F,quote=F,sep=",",row.names=F)

plotTitle <- paste(colnames(specs)[1],":",specs[[1]],", ", colnames(specs)[2],":",specs[2],", ",colnames(specs)[3],":",specs[3],", ",colnames(specs)[4],":",specs[4],",",colnames(specs)[5],":",1e-7,", ",colnames(specs)[6],":",specs[6],", ",colnames(specs)[7],":",specs[7],", ",colnames(specs)[9],":",specs[9],sep="")

gen_nam <- paste("Gen",fit$Generation,sep="")

fitt<-data.frame(t(fit[,-1:-8]))

colnames(fitt)<-gen_nam

#fitt$Location <- factor(rep(paste("P",seq(1,10,1),sep=""),30),levels = paste("P",seq(1,10,1),sep=""))
fitt$Location <- factor(rep(paste("P",1:100,sep=""),3), levels = unique(rep(paste("P",1:100,sep=""),3)))
fitt$Environment <- factor(rep(paste("E",1:10,sep=""),30), levels = unique(rep(paste("E",1:10,sep=""),30)))

fitt$Type <- as.factor(c(rep("Fit",100),rep("Freq",100),rep("Env",100)))

#VCF files are filtered with vcftools, as it is much faster than R. The filtering for MAF > 0.05 is accomplished with the following code:
#vcftools --vcf T1.vcf --maf 0.05 --recode --recode-INFO-all --out T1

# number of cores to use for parallel processing
# cores <- 3
# cl <- makeCluster(cores)
# registerDoParallel(cl)

#VCF files are filtered with vcftools, as it is much faster than R. The filtering for MAF > 0.05 is accomplished with the following code:
#vcf1 <- read.vcfR(paste("results/SLiM_output/VCF_files/T1_",seed,".vcf",sep=""))
# geno1 <- vcf1@gt[,-1] # Remove 1st column, which is 'Format'
# position1 <- getPOS(vcf1) # Positions in bp
# chromosome1 <- getCHROM(vcf1) # Chromosome information
# 
# rm(vcf1)
# gc()

# vcf2 <- read.vcfR(paste("results/SLiM_output/VCF_files/T2_",seed,".vcf",sep=""))
# geno2 <- vcf2@gt[,-1] # Remove 1st column, which is 'Format'
# position2 <- getPOS(vcf2) # Positions in bp
# chromosome2 <- getCHROM(vcf2) # Chromosome information
# 
# rm(vcf2)
# gc()

#Population values
Pop <- NULL
for(i in 1:100){
  Pop <- c(Pop,paste("P",i,sep=""))
}

vcf1_filt <- read.vcfR(paste("results/SLiM_output/VCF_files/T2_",seed,"_filtered_subset.recode.vcf",sep=""))
#vcf1_filt <- read.vcfR("results/SLiM_output/VCF_files/T2_1621040440900_filtered_subset.recode.vcf")
geno1_filt <- vcf1_filt@gt[,-1] # Remove 1st column, which is 'Format'
position1_filt <- as.numeric(getPOS(vcf1_filt)) # Positions in bp
chromosome1_filt <- as.numeric(getCHROM(vcf1_filt)) # Chromosome information
positionM2<-paste("M",position1_filt[which(position1_filt==224999)],sep="")

No_A<-unname(dim(vcf1_filt)[1]) #Get the number of filtered alleles

#vcf1_m2<-vcf1_filt[which(vcf1_filt@fix[,2]==224999),]
#FST_tot1<-pairwise_genetic_diff(vcf1_filt,pops = Pop)

rm(vcf1_filt)

# vcf2_filt <- read.vcfR(paste("results/SLiM_output/VCF_files/T2_",seed,".recode.vcf",sep=""))
# geno2_filt <- vcf2_filt@gt[,-1] # Remove 1st column, which is 'Format'
# position2_filt <- as.numeric(getPOS(vcf2_filt)) # Positions in bp
# chromosome2_filt <- as.numeric(getCHROM(vcf2_filt)) # Chromosome information
# 
# rm(vcf2_filt)
# gc()

### Convert VCF to 012 format ####
# Character matrix containing the genotypes
# individuals in columns

#Create Genotype matrix
G1f <- matrix(NA, nrow = nrow(geno1_filt), ncol = ncol(geno1_filt))
G1f[geno1_filt %in% c("0/0", "0|0")] <- 0
G1f[geno1_filt %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
G1f[geno1_filt %in% c("1/1", "1|1")] <- 2

#Create Genotype matrix
# G2f <- matrix(NA, nrow = nrow(geno2_filt), ncol = ncol(geno2_filt))
# G2f[geno2_filt %in% c("0/0", "0|0")] <- 0
# G2f[geno2_filt %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
# G2f[geno2_filt %in% c("1/1", "1|1")] <- 2

#stopCluster(cl)

#Check number of duplicate positions
sum(duplicated(position1_filt))
position1_filt[duplicated(position1_filt)]
#sum(duplicated(position2_filt))

for(i in 1:length(position1_filt)){
  if(duplicated(position1_filt)[i]){
    position1_filt[i]<-position1_filt[i]+0.5
    #print(position1_filt[i])
  }
}

#Create Genotype matrix
# G1 <- matrix(NA, nrow = nrow(geno1), ncol = ncol(geno1))
# G1[geno1 %in% c("0/0", "0|0")] <- 0
# G1[geno1  %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
# G1[geno1 %in% c("1/1", "1|1")] <- 2

Start<-seq(1,ncol(geno1_filt),10)
Stop<-seq(10,ncol(geno1_filt),10)

#No subsample (neutral / size variation)
# St=1
# demes=as.numeric(unname(specs[,9:18]))
# for(i in 1:9){
#   St=c(St,St[i]+demes[i])
# }
# 
# Start<-NULL
# for(i in seq(0,9000,1000)){
#   Start<-c(Start,St+i)
# }
# 
# Sp=demes[1]
# for(i in 1:9){
#   Sp=c(Sp,Sp[i]+demes[i+1])
# }
# 
# Stop<-NULL
# for(i in seq(0,9000,1000)){
#   Stop<-c(Stop,Sp+i)
# }

Pop_afreq1<-NULL
for(i in 1:100){
  Pop_afreq1<-rbind(Pop_afreq1,rowSums(G1f[,Start[i]:Stop[i]])/(2*ncol(G1f[,Start[i]:Stop[i]])))
}
Pop_afreq1<-data.frame(Pop_afreq1)
colnames(Pop_afreq1)<-paste("M",position1_filt,sep="")

# Pop_afreq2<-NULL
# for(i in 1:100){
#   Pop_afreq2<-rbind(Pop_afreq2,rowSums(G2f[,Start[i]:Stop[i]])/(2*ncol(G2f[,Start[i]:Stop[i]])))
# }
# Pop_afreq2<-data.frame(Pop_afreq2)
# colnames(Pop_afreq2)<-paste("M",position2_filt,sep="")

#Pop_afreq2_inP1<-Pop_afreq2[colnames(Pop_afreq2)%in%colnames(Pop_afreq1)]

#Combine allele frequency, genomic position, & meta data into one data frame
#PreN<-data.frame(colMeans(Pop_afreq1),position1_filt,PO_pre_filt,GO_pre_filt,MT_pre_filt,row.names=MID_pre_filt,stringsAsFactors = F)
#colnames(PreN)<-c("AF_l","PP","PO","GO","LT")
#PreN[grep("MT=2",PreN$ID1),]

#PostN<-data.frame(a_freq2,position2,ID2)
#PostN[grep("MT=2",ID2),]

#Subset the environmental variables for the generation you're considering (make sure the M2 AF and environmental data are not being compare across generations)
# envPop<-data.frame(fitt[fitt$Type=="Env",gen_nam[length(gen_nam)-30]]) #300 years prior to the end of the simulation is taken as the "before environmental shift" time
# names(envPop) <- "envSelect"


#######################################################
# Set up environmental vartiables for Gradient Forest
#######################################################
cg<-read.table(paste("results/SLiM_output/CG_files/",seed,"_fitnessmat_pop.txt",sep=""),header=F)

colnames(cg)<-rep(paste("H",seq(1:100),sep=""))
cg$Transplant<-paste("T",seq(1:100),sep="")

cg_df <- melt(cg,  id.vars ="Transplant",  variable.name = "Home",
              value.name = "Fitness")

cg_df$Transplant<-factor(cg_df$Transplant,levels=unique(cg_df$Transplant))

#Subset the environmental variables for the generation you're considering (make sure the M2 AF and environmental data are not being compare across generations)
envPop<-data.frame(fitt[fitt$Type=="Env",gen_nam[length(gen_nam)]]) #300 years prior to the end of the simulation is taken as the "before environmental shift" time
names(envPop) <- "envPop"

fakeEnv <- envPop$envPop + rnorm(nrow(envPop),0,1.3)

cor(envPop$envPop, fakeEnv)

Pop <- levels(fitt$Location)
Popsenv <- data.frame(Pop, envPop,fakeEnv)

nfake <- 10

Popsenv[,4:(3+nfake)] <- NA

cov1 <- genPositiveDefMat(nfake,covMethod="unifcorrmat" )
head(cov1)

a<- mvrnorm(nrow(Popsenv),mu=rep(0, nfake), Sigma=cov1$Sigma)

Popsenv[,4:(3+nfake)] <- a

# head(Popsenv)
# tail(Popsenv)

sel_env_cols <- 2
all_env_cols <- 2:ncol(Popsenv)
sel_env_cols_plus2 <- 2:3

cov_allEnv <- cov(Popsenv[,all_env_cols])
cov_selEnv <- cov(Popsenv[,c(sel_env_cols,sel_env_cols)])
cov_selEnv_plus2 <- cov(Popsenv[,sel_env_cols_plus2])

start_time <- Sys.time()
for (i in 1:nrow(cg_df)){
  # get the row in PopsenvStnd for the common garden genotype
  row1 = which(Popsenv==gsub("T","P",as.character(cg_df$Transplant[i])))
  # get the row in PopsenvStnd for the genotype source
  row2 = which(Popsenv==gsub("H","P",as.character(cg_df$Home[i])))
  
  # Look up the envi (all var)
  (envpop1_all <- Popsenv[row1,all_env_cols])
  (envpop2_all <- Popsenv[row2,all_env_cols])
  
  # Look up the envi selected only
  (envpop1_sel <- Popsenv[row1,sel_env_cols])
  (envpop2_sel <- Popsenv[row2,sel_env_cols])
  
  # Look up the envi selected plus 2 correlated env
  (envpop1_sel_plus2 <- Popsenv[row1,sel_env_cols_plus2])
  (envpop2_sel_plus2 <- Popsenv[row2,sel_env_cols_plus2])  
  
  
  ### Calculate the environmental distance between the two rows
  
  # Euclidean distance for selective environments
  cg_df$EdSelEnv[i] <- dist(rbind(envpop1_sel,
                                  envpop2_sel))
  
  # Mahalanobis distance for selective environments
  #cg_df$MdSelEnv[i] <- mahalanobis(sapply(envpop1_sel,as.numeric),
  #                                 sapply(envpop2_sel,as.numeric),
  #                                 cov_selEnv)
  
  # Euclidean distance for ALL environments
  cg_df$EdAllEnv[i] <- dist(rbind(envpop1_all,
                                  envpop2_all))
  
  # Mahalanobis distance for ALL environments
  cg_df$MdAllEnv[i] <- mahalanobis(as.numeric(unlist(envpop1_all)),
                                   as.numeric(unlist(envpop2_all)),
                                   cov_allEnv)
  
  # Euclidean distance for selective environments + 2 env
  cg_df$EdSelEnvPlus2[i] <- dist(rbind(envpop1_sel_plus2,
                                       envpop2_sel_plus2))
  
  # Mahalanobis distance for selective environments + 2 env
  cg_df$MdSelEnvPlus2[i] <- mahalanobis(sapply(envpop1_sel_plus2,as.numeric), 
                                         sapply(envpop2_sel_plus2,as.numeric), 
                                         cov_selEnv_plus2) 
}
end_time <- Sys.time()
writeLines("Environmental distance calculation:")
print(end_time - start_time)
# this take a few minutes

# envPop.shift<-data.frame(fitt[fitt$Type=="Env",gen_nam[length(gen_nam)]])
# names(envPop.shift) <- "envSelect"

envTab <- Popsenv[,-1]

#Merge the population specific allele frequencies of all neutral (M1) alleles with the population specific frequency of the selected (M2) allele
alFreq<-cbind(data.frame(Pop_afreq1)) 

##############################################
# Chunk to fit GF models to minor allele frequencies at the level of
# populations
# GF is fit to each SNP individually to 
# ease computational / memory burden

#cl <- makeCluster(cores)
#registerDoParallel(cl)

##### added by MCF, running all loci in on model ##########
vars <- colnames(envTab)

gfMod <- gradientForest(data=data.frame(envTab[, vars], alFreq),
                        predictor.vars=vars,
                        response.vars=colnames(alFreq),
                        corr.threshold=0.5, 
                        ntree=500, 
                        trace=T)

#cImp <- cumimp(gfMod, "envSelect", type="Species")
#cImp <- data.frame(rbindlist(cImp, idcol="allele"))

#stopCluster(cl)

# Calculate genomic offset
# note that I am doing this for the avearge across all alleles since 
# GF was fit to all alleles simultaneously
# The more correct way is to calculate offset for adaptive alleles only,
# either individually or for a model fit to just those alleles.

plot(gfMod, plot.type="Overall.Importance")

GF_dist_gen <- NULL

for(i in 1:ncol(envTab)){
  tmp <- data.frame(envTab[,i])
  colnames(tmp) <- colnames(envTab)[i]
  
  writeLines(paste("Predicting ",colnames(envTab)[i]))
  gftmp <- predict(gfMod, tmp)
  
  GF_dist_gen <- data.frame(c(GF_dist_gen,gftmp))
}

# gfTrans1 <- predict(gfMod, envPop)
# colnames(gfTrans1)<-"C.Imp_genome_before"
# 
# gfTrans2 <- predict(gfMod, envPop.shift)
# colnames(gfTrans2)<-"C.Imp_genome_after"

# offset needs to be considered using absolute values ()
# offset <- gfTrans2-gfTrans1
# offset <- gfTrans2-gfTrans1
# colnames(offset)<-"D_C.Imp_genome"

M2<-data.frame(Pop_afreq1$M224999)
colnames(M2)<-"M2"

# Selected allele only
gfM2 <- gradientForest(data=data.frame(envTab[,vars], M2),
                       predictor.vars=vars,
                       response.vars=colnames(M2),
                       corr.threshold=0.5, 
                       ntree=500, 
                       trace=T)

GF_dist_caus <- NULL

for(i in 1:ncol(envTab)){
  tmp <- data.frame(envTab[,i])
  colnames(tmp) <- colnames(envTab)[i]
  
  writeLines(paste("Predicting ",colnames(envTab)[i]))
  gftmp <- predict(gfM2, tmp)
  
  GF_dist_caus <- data.frame(c(GF_dist_caus,gftmp))
}

# gfM2Trans1 <- predict(gfM2, envPop)
# colnames(gfM2Trans1)<-"C.Imp_M2_before"
# 
# gfM2Trans2 <- predict(gfM2, envPop.shift)
# colnames(gfM2Trans2)<-"C.Imp_M2_after"
# 
# M2offset <- gfM2Trans2-gfM2Trans1
# colnames(M2offset)<-"D_C.Imp_M2"

# gfM1 <- gradientForest(data=data.frame(envPop, 1-M2),
#                        predictor.vars=colnames(envPop),
#                        response.vars=colnames(M2),
#                        corr.threshold=0.5, 
#                        ntree=500, 
#                        trace=T)
# 
# gfM1Trans1 <- predict(gfM1, envPop)
# colnames(gfM1Trans1)<-"C.Imp_M1_before"

##############################################################################
#Get  Weir & Cockerham F_ST values from the VCF files and use population data

#Create an object listing every population in the whole dataset
PopsALL <- NULL
for(j in rep(1:100)){
  for(i in rep(j,10)){
    PopsALL <- c(PopsALL,i)
  }
}

#Create an object splitting a single population into a Pre ("T1") and Post ("T2") "population"
PopsP <- c(rep("T1",100),rep("T2",100))

# cores<-3
# cl <- makeCluster(cores)
# registerDoParallel(cl)

# G2 <- matrix(NA, nrow = nrow(geno2), ncol = ncol(geno2))
# G2[geno2 %in% c("0/0", "0|0")] <- 0
# G2[geno2  %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
# G2[geno2 %in% c("1/1", "1|1")] <- 2
# 
# #Calculate allele frequencies across the whole meta population
# a_freq1 <- rowSums(G1)/(2*ncol(G1))
# a_freq2 <- rowSums(G2)/(2*ncol(G2))

#Prepare the Pre (G1) and Post (G2) G matrices for FST calculation
# Gt1<-t(G1)
# rownames(Gt1)<-PopsALL
# #colnames(Gt1)<-MID_pre
# colnames(Gt1)<-paste("M",position1,sep="")
# 
# # Gt2<-t(G2)
# # rownames(Gt2)<-PopsALL
# # #colnames(Gt2)<-MID_post
# # colnames(Gt2)<-paste("M",position2,sep="")
# 
# #Add loop to iterate across each TRUE population, not just x-location
# #Filter by population (same pop. before and after env. shift)
# listGt<-list()
# for(i in 1:100){
#   Gt1_i<-data.frame(Gt1[rownames(Gt1)==i,])
#   Gt2_i<-data.frame(Gt2[rownames(Gt2)==i,])
#   Gt1_i[setdiff(names(Gt2_i), names(Gt1_i))] <- 0
#   Gt2_i[setdiff(names(Gt1_i), names(Gt2_i))] <- 0
#   listGt[[i]]<-rbind(Gt1_i,Gt2_i)
# }
# 
# #Filter the files based on Major and Minor AF filtration, then seed the files so that all variants are present at both time points to be compared
# a_freq<-list()
# a_freq_filt<-list()
# listGt_filt<-list()
# for(i in 1:100){
#   a_freq[[i]] <- colSums(listGt[[i]])/(2*nrow(listGt[[i]]))
#   a_freq_filt[[i]] <- a_freq[[i]][a_freq[[i]]>0.05 & a_freq[[i]]<0.95]
#   listGt_filt[[i]]<-listGt[[i]][colnames(listGt[[i]])%in%names(a_freq_filt[[i]])]
# }
# 
# #Calculate per locus FST values
# listPfst<-list()
# for(i in 1:100){
#   listPfst[[i]]<-MakeDiploidFSTMat(SNPmat = listGt_filt[[i]], locusNames = colnames(listGt_filt[[i]]), popNames = PopsP)
# }
# 
# #Filter out NA values
# listPfst_noNa<-list()
# for(i in 1:100){
#   listPfst_noNa<-lapply(listPfst,function(x) x[!is.na(x$FST),])
# }
# 
# #Calculate FST values for each populations
# FST_genome_pop<-NULL
# for(i in 1:100){
#   FST_genome_pop<-c(FST_genome_pop,mean(listPfst_noNa[[i]]$T1)/mean(listPfst_noNa[[i]]$T2))
# }

getFSTs_diploids = function(popNameList, SNPDataColumn){  
  #eliminating the missing data for this locus
  popnames=unlist(as.character(popNameList))
  popNameTemp=popnames[which(SNPDataColumn!=9)]
  snpDataTemp=SNPDataColumn[SNPDataColumn!=9]
  
  HetCounts <- tapply(snpDataTemp, list(popNameTemp,snpDataTemp), length)
  HetCounts[is.na(HetCounts)] = 0
  
  #Case: all individuals are genetically identical at this locus
  if(dim(HetCounts)[2]==1){
    return (list(He=NA,FST=NA, T1=NA, T2=NA,FSTNoCorr=NA, T1NoCorr=NA, T2NoCorr=NA,meanAlleleFreq = NA))
  }
  
  if(dim(HetCounts)[2]==2){
    if(paste(colnames(HetCounts),collapse="")=="01"){HetCounts=cbind(HetCounts,"2"=0)}
    if(paste(colnames(HetCounts),collapse="")=="12"){HetCounts=cbind("0"=0,HetCounts)} 
    if(paste(colnames(HetCounts),collapse="")=="02"){HetCounts=cbind(HetCounts[,1],"1"=0, HetCounts[,2])}
  }
  
  out = WC_FST_Diploids_2Alleles(HetCounts)	
  return(out)
}

MakeDiploidFSTMat_2<-function(SNPmat,locusNames,popNames){
  locusname <- unlist(locusNames)
  popname <- unlist(popNames)
  snplevs <- levels(as.factor(unlist(SNPmat)))
  if(any(!(snplevs%in%c(0,1,2,9)))==TRUE) {
    print("Error: Your snp matrix has a character other than 0,1,2 or 9")
    break
  }
  if (dim(SNPmat)[1] != length(popname)) {
    print("Error: your population names do not match your SNP matrix")
    break
  }
  if (dim(SNPmat)[2] != length(locusname)) {
    print("Error:  your locus names do not match your SNP matrix")
    break
  }
  writeLines("Calculating FSTs, may take a few minutes...")
  nloci <- length(locusname)
  FSTmat <- matrix(NA, nrow = nloci, ncol = 8)
  for (i in 1:nloci) {
    FSTmat[i, ] = unlist(getFSTs_diploids(popname, SNPmat[,i]))
    if (i%%10000 == 0) {
      print(paste(i, "done of", nloci))
    }
  }
  outTemp = as.data.frame(FSTmat)
  outTemp = cbind(locusname, outTemp)
  colnames(outTemp) = c("LocusName", "He", "FST", "T1", "T2", 
                        "FSTNoCorr", "T1NoCorr", "T2NoCorr", "meanAlleleFreq")
  return(outTemp)
}

########################################################
#Calculate per locus FST values pre environmental shift

Gt1f<-t(G1f)
colnames(Gt1f)<-paste("M",position1_filt,sep="")

# Gt2f<-t(G2f)
# colnames(Gt2f)<-paste("M",position2_filt,sep="")

Pre_geno<-data.frame(PopsALL,Gt1f)
colnames(Pre_geno)[1]<-"Locality"

# Post_geno<-data.frame(PopsALL,Gt2f)
# colnames(Post_geno)[1]<-"Locality"

#Pull out only the M2 allele

# Gt1_m2<-matrix(Gt1f[,positionM2])
# colnames(Gt1_m2)<-positionM2
# rownames(Gt1_m2)<-PopsALL

# Gt2_m2<-matrix(Gt2[,positionM2])
# colnames(Gt2_m2)<-positionM2
# rownames(Gt2_m2)<-PopsALL
# 
# listGt_m2<-list()
# for(i in 1:100){
#   Gt1_m2i<-matrix(Gt1_m2[rownames(Gt1_m2)==i,])
#   colnames(Gt1_m2i)<-positionM2
#   rownames(Gt1_m2i)<-PopsALL[which(PopsALL==i)]
#   Gt2_m2i<-matrix(Gt2_m2[rownames(Gt2_m2)==i,])
#   colnames(Gt2_m2i)<-positionM2
#   rownames(Gt2_m2i)<-PopsALL[which(PopsALL==i)]
#   listGt_m2[[i]]<-rbind(Gt1_m2i,Gt2_m2i)
# }

#Calculate M2 locus FST values
# list_M2_Pfst<-list()
# for(i in 1:100){
#   list_M2_Pfst[[i]]<-MakeDiploidFSTMat_2(SNPmat = listGt_m2[[i]], locusNames = colnames(listGt_m2[[i]]), popNames = PopsP)
# }

#Filter out NA values
# list_M2_Pfst_noNa<-list()
# for(i in 1:100){
#   list_M2_Pfst_noNa<-lapply(list_M2_Pfst,function(x) x[!is.na(x$FST),])
# }

#Calculate M2 FST values for each populations
# FST_M2_pop<-NULL
# for(i in 1:100){
#   FST_M2_pop<-c(FST_M2_pop,mean(list_M2_Pfst_noNa[[i]]$T1)/mean(list_M2_Pfst_noNa[[i]]$T2))
# }
# 
# for(i in 1:length(FST_M2_pop)){
#   if(is.na(FST_M2_pop[i])){
#     FST_M2_pop[i]<-0
#   }
# }

# G_m2<-data.frame(row.names(Gt2_m2))
# G_m2$M224999<-data.frame(Gt2_m2)
# colnames(G_m2)<-c("Locality","M224999")
# 
# FST_m2<-pairwise.WCfst(G_m2,diploid=TRUE)
# 
# write.table(FST_m2,file=paste("~/Desktop/",seed,"_Pre_FST_allele.txt",sep=""),sep=",",col.names = T,row.names = T,quote = F)
# write.table(Summary_Pop[,c("Env_before","Env_after","F_ST_M2_bef.aft.")],file=paste("~/Desktop/",seed,"_sum_Pop_subset.txt",sep=""),sep=",",col.names = T,row.names = F,quote = F)


CI<-NULL
for(i in 1:nrow(GF_dist_gen)){
  for(j in 1:nrow(GF_dist_gen)){
    CI<-c(CI,dist(rbind(GF_dist_gen[i,],GF_dist_gen[j,])))
    #CI<-c(CI,CI_bf[j]-CI_bf[i])
  }
}

cg_df$D_CI_allEnv<-CI

GFcaus_dist_gen<-matrix(GF_dist_gen[,1])
CI_causal<-NULL
for(i in 1:nrow(GFcaus_dist_gen)){
  for(j in 1:nrow(GFcaus_dist_gen)){
    CI_causal<-c(CI_causal,dist(rbind(GFcaus_dist_gen[i,],GFcaus_dist_gen[j,])))
    #CI<-c(CI,CI_bf[j]-CI_bf[i])
  }
}

cg_df$D_CI_causEnv<-CI_causal


CI_sel<-NULL
for(i in 1:nrow(GF_dist_caus)){
  for(j in 1:nrow(GF_dist_caus)){
    CI_sel<-c(CI_sel,dist(rbind(GF_dist_caus[i,],GF_dist_caus[j,])))
    #CI<-c(CI,CI_bf[j]-CI_bf[i])
  }
}

cg_df$D_CI_caus_allEnv<-CI_sel

GFcaus_dist_caus<-matrix(GF_dist_caus[,1])
CIcaus_causal<-NULL
for(i in 1:nrow(GFcaus_dist_caus)){
  for(j in 1:nrow(GFcaus_dist_caus)){
    CIcaus_causal<-c(CIcaus_causal,dist(rbind(GFcaus_dist_caus[i,],GFcaus_dist_caus[j,])))
    #CI<-c(CI,CI_bf[j]-CI_bf[i])
  }
}

cg_df$D_CI_caus_causEnv<-CIcaus_causal


options(scipen = 0)

########################################################
#Calculate per locus FST values pre environmental shift for Edge and cCore populations

#Whole genome
Gt1f<-t(G1f)
colnames(Gt1f)<-paste("M",position1_filt,sep="")

Pre_geno <- data.frame(PopsALL,Gt1f)
colnames(Pre_geno)[1]<-"Locality"

Edges <- c(1,5,6,10,12,19,41,50,51,60,82,89,91,95,96,100)
Cores <- c(34, 35, 36, 37, 44, 45, 46, 47, 54, 55, 56, 57, 64, 65, 66, 67)

Pre_geno_comb <- Pre_geno[Pre_geno$Locality%in%c(Edges,Cores),]

start_time <- Sys.time()
Pre_FST_comb<-pairwise.WCfst(Pre_geno_comb,diploid = T)
end_time <- Sys.time()
print(end_time - start_time)

mean(Pre_FST_comb, na.rm=T)

Gt1f_sel<-matrix(Gt1f[,colnames(Gt1f)%in%positionM2])
colnames(Gt1f_sel)<-positionM2
Pre_geno_sel<-data.frame(PopsALL,Gt1f_sel)
colnames(Pre_geno_sel)[1]<-"Locality"

Pre_geno_sel_comb <- Pre_geno_sel[Pre_geno_sel$Locality%in%c(Edges,Cores),]

start_time <- Sys.time()
Pre_FST_sel_comb<-genet.dist(Pre_geno_sel_comb,diploid = T,method="Fst")
end_time <- Sys.time()
print(end_time - start_time)

Pre_FST_Edge <- Pre_FST_comb[row.names(Pre_FST_comb)%in%Edges,]
Pre_FST_Cores <- Pre_FST_comb[row.names(Pre_FST_comb)%in%Cores,]

Pre_FST_sel_comb<-as.matrix(Pre_FST_sel_comb)
row.names(Pre_FST_sel_comb)<-row.names(Pre_FST_comb)

Pre_FST_sel_Edge <- Pre_FST_sel_comb[row.names(Pre_FST_sel_comb)%in%Edges,]
Pre_FST_sel_Cores <- Pre_FST_sel_comb[row.names(Pre_FST_sel_comb)%in%Cores,]


Pre_FST_Cores[is.na(Pre_FST_Cores)]<-0
Pre_FST_Edge[is.na(Pre_FST_Edge)]<-0
Pre_FST_sel_Cores[is.na(Pre_FST_sel_Cores)]<-0
Pre_FST_sel_Edge[is.na(Pre_FST_sel_Edge)]<-0

cg_df_Edges <- cg_df[as.numeric(substr(as.character(cg_df$Home),2,nchar(as.character(cg_df$Home))))%in%c(Edges,Cores),]
cg_df_Edges <- cg_df_Edges[as.numeric(substr(as.character(cg_df_Edges$Transplant),2,nchar(as.character(cg_df_Edges$Transplant))))%in%Edges,]

cg_df_Edges$Pre_FST <- as.vector(Pre_FST_Edge)
cg_df_Edges$Pre_FSt_sel <- as.vector(Pre_FST_sel_Edge)

cg_df_Cores <- cg_df[as.numeric(substr(as.character(cg_df$Home),2,nchar(as.character(cg_df$Home))))%in%c(Edges,Cores),]
cg_df_Cores <- cg_df_Cores[as.numeric(substr(as.character(cg_df_Cores$Transplant),2,nchar(as.character(cg_df_Cores$Transplant))))%in%Cores,]

cg_df_Cores$Pre_FST <- as.vector(Pre_FST_Cores)
cg_df_Cores$Pre_FSt_sel <- as.vector(Pre_FST_sel_Cores)


#compare edge pops
#plot(x=cg_df_Edges$D_CI_allEnv, y=cg_df_Edges$Fitness)
dat_sum<-cor.test(x=cg_df_Edges$D_CI_allEnv, y=cg_df_Edges$Fitness, method = "spearman")
paste("r = ",round(dat_sum$estimate[[1]],3),"\\nslope = ", round(dat_sum$statistic[[1]],3),"\\np-value = ", signif(dat_sum$p.value,3),sep="")
#summary(dat_sum)
GF_off_genome_allEnv_Edges <- round(dat_sum$estimate[[1]],3)

dat_sum<-cor.test(x=cg_df_Edges$D_CI_causEnv, y=cg_df_Edges$Fitness, method = "spearman")
paste("r = ",round(dat_sum$estimate[[1]],3),"\\nslope = ", round(dat_sum$statistic[[1]],3),"\\np-value = ", signif(dat_sum$p.value,3),sep="")
GF_off_genome_causEnv_Edges <- round(dat_sum$estimate[[1]],3)

dat_sum<-cor.test(x=cg_df_Edges$D_CI_caus_allEnv, y=cg_df_Edges$Fitness, method = "spearman")
paste("r = ",round(dat_sum$estimate[[1]],3),"\\nslope = ", round(dat_sum$statistic[[1]],3),"\\np-value = ", signif(dat_sum$p.value,3),sep="")
GF_off_causal_allEnv_Edges <- round(dat_sum$estimate[[1]],3)

dat_sum<-cor.test(x=cg_df_Edges$D_CI_caus_causEnv, y=cg_df_Edges$Fitness, method = "spearman")
paste("r = ",round(dat_sum$estimate[[1]],3),"\\nslope = ", round(dat_sum$statistic[[1]],3),"\\np-value = ", signif(dat_sum$p.value,3),sep="")
GF_off_causal_causEnv_Edges <- round(dat_sum$estimate[[1]],3)

dat_sum<-cor.test(x=cg_df_Edges$EdSelEnv, y=cg_df_Edges$Fitness, method = "spearman")
paste("r = ",round(dat_sum$estimate[[1]],3),"\\nslope = ", round(dat_sum$statistic[[1]],3),"\\np-value = ", signif(dat_sum$p.value,3),sep="")
Env_ED_caus_env_Edges <- round(dat_sum$estimate[[1]],3)

dat_sum<-cor.test(x=cg_df_Edges$MdSelEnv, y=cg_df_Edges$Fitness, method = "spearman")
paste("r = ",round(dat_sum$estimate[[1]],3),"\\nslope = ", round(dat_sum$statistic[[1]],3),"\\np-value = ", signif(dat_sum$p.value,3),sep="")
Env_MD_caus_env_Edges <- round(dat_sum$estimate[[1]],3)

dat_sum<-cor.test(x=cg_df_Edges$EdAllEnv, y=cg_df_Edges$Fitness, method = "spearman")
paste("r = ",round(dat_sum$estimate[[1]],3),"\\nslope = ", round(dat_sum$statistic[[1]],3),"\\np-value = ", signif(dat_sum$p.value,3),sep="")
Env_ED_all_env_Edges <- round(dat_sum$estimate[[1]],3)

dat_sum<-cor.test(x=cg_df_Edges$MdAllEnv, y=cg_df_Edges$Fitness, method = "spearman")
paste("r = ",round(dat_sum$estimate[[1]],3),"\\nslope = ", round(dat_sum$statistic[[1]],3),"\\np-value = ", signif(dat_sum$p.value,3),sep="")
Env_MD_all_env_Edges <- round(dat_sum$estimate[[1]],3)

dat_sum<-cor.test(x=cg_df_Edges$EdSelEnvPlus2, y=cg_df_Edges$Fitness, method = "spearman")
paste("r = ",round(dat_sum$estimate[[1]],3),"\\nslope = ", round(dat_sum$statistic[[1]],3),"\\np-value = ", signif(dat_sum$p.value,3),sep="")
Env_ED_causP2_env_Edges <- round(dat_sum$estimate[[1]],3)

dat_sum<-cor.test(x=cg_df_Edges$MdSelEnvPlus2, y=cg_df_Edges$Fitness, method = "spearman")
paste("r = ",round(dat_sum$estimate[[1]],3),"\\nslope = ", round(dat_sum$statistic[[1]],3),"\\np-value = ", signif(dat_sum$p.value,3),sep="")
Env_MD_causP2_env_Edges <- round(dat_sum$estimate[[1]],3)

dat_sum<-cor.test(x=cg_df_Edges$Pre_FST, y=cg_df_Edges$Fitness, method = "spearman")
paste("r = ",round(dat_sum$estimate[[1]],3),"\\nslope = ", round(dat_sum$statistic[[1]],3),"\\np-value = ", signif(dat_sum$p.value,3),sep="")
FST_genome_Edges <- round(dat_sum$estimate[[1]],3)

dat_sum<-cor.test(x=cg_df_Edges$Pre_FSt_sel, y=cg_df_Edges$Fitness, method = "spearman")
paste("r = ",round(dat_sum$estimate[[1]],3),"\\nslope = ", round(dat_sum$statistic[[1]],3),"\\np-value = ", signif(dat_sum$p.value,3),sep="")
FST_causal_Edges <- round(dat_sum$estimate[[1]],3)

#compare core pops
dat_sum<-cor.test(x=cg_df_Cores$D_CI_allEnv, y=cg_df_Cores$Fitness, method = "spearman")
paste("r = ",round(dat_sum$estimate[[1]],3),"\\nslope = ", round(dat_sum$statistic[[1]],3),"\\np-value = ", signif(dat_sum$p.value,3),sep="")
GF_off_genome_allEnv_Cores <- round(dat_sum$estimate[[1]],3)

dat_sum<-cor.test(x=cg_df_Cores$D_CI_causEnv, y=cg_df_Cores$Fitness, method = "spearman")
paste("r = ",round(dat_sum$estimate[[1]],3),"\\nslope = ", round(dat_sum$statistic[[1]],3),"\\np-value = ", signif(dat_sum$p.value,3),sep="")
GF_off_genome_causEnv_Cores <- round(dat_sum$estimate[[1]],3)

dat_sum<-cor.test(x=cg_df_Cores$D_CI_caus_allEnv, y=cg_df_Cores$Fitness, method = "spearman")
paste("r = ",round(dat_sum$estimate[[1]],3),"\\nslope = ", round(dat_sum$statistic[[1]],3),"\\np-value = ", signif(dat_sum$p.value,3),sep="")
GF_off_causal_allEnv_Cores <- round(dat_sum$estimate[[1]],3)

dat_sum<-cor.test(x=cg_df_Cores$D_CI_caus_causEnv, y=cg_df_Cores$Fitness, method = "spearman")
paste("r = ",round(dat_sum$estimate[[1]],3),"\\nslope = ", round(dat_sum$statistic[[1]],3),"\\np-value = ", signif(dat_sum$p.value,3),sep="")
GF_off_causal_causEnv_Cores <- round(dat_sum$estimate[[1]],3)

dat_sum<-cor.test(x=cg_df_Cores$EdSelEnv, y=cg_df_Cores$Fitness, method = "spearman")
paste("r = ",round(dat_sum$estimate[[1]],3),"\\nslope = ", round(dat_sum$statistic[[1]],3),"\\np-value = ", signif(dat_sum$p.value,3),sep="")
Env_ED_caus_env_Cores <- round(dat_sum$estimate[[1]],3)

dat_sum<-cor.test(x=cg_df_Cores$MdSelEnv, y=cg_df_Cores$Fitness, method = "spearman")
paste("r = ",round(dat_sum$estimate[[1]],3),"\\nslope = ", round(dat_sum$statistic[[1]],3),"\\np-value = ", signif(dat_sum$p.value,3),sep="")
Env_MD_caus_env_Cores <- round(dat_sum$estimate[[1]],3)

dat_sum<-cor.test(x=cg_df_Cores$EdAllEnv, y=cg_df_Cores$Fitness, method = "spearman")
paste("r = ",round(dat_sum$estimate[[1]],3),"\\nslope = ", round(dat_sum$statistic[[1]],3),"\\np-value = ", signif(dat_sum$p.value,3),sep="")
Env_ED_all_env_Cores <- round(dat_sum$estimate[[1]],3)

dat_sum<-cor.test(x=cg_df_Cores$MdAllEnv, y=cg_df_Cores$Fitness, method = "spearman")
paste("r = ",round(dat_sum$estimate[[1]],3),"\\nslope = ", round(dat_sum$statistic[[1]],3),"\\np-value = ", signif(dat_sum$p.value,3),sep="")
Env_MD_all_env_Cores <- round(dat_sum$estimate[[1]],3)

dat_sum<-cor.test(x=cg_df_Cores$EdSelEnvPlus2, y=cg_df_Cores$Fitness, method = "spearman")
paste("r = ",round(dat_sum$estimate[[1]],3),"\\nslope = ", round(dat_sum$statistic[[1]],3),"\\np-value = ", signif(dat_sum$p.value,3),sep="")
Env_ED_causP2_env_Cores <- round(dat_sum$estimate[[1]],3)

dat_sum<-cor.test(x=cg_df_Cores$MdSelEnvPlus2, y=cg_df_Cores$Fitness, method = "spearman")
paste("r = ",round(dat_sum$estimate[[1]],3),"\\nslope = ", round(dat_sum$statistic[[1]],3),"\\np-value = ", signif(dat_sum$p.value,3),sep="")
Env_MD_causP2_env_Cores <- round(dat_sum$estimate[[1]],3)

dat_sum<-cor.test(x=cg_df_Cores$Pre_FST, y=cg_df_Cores$Fitness, method = "spearman")
paste("r = ",round(dat_sum$estimate[[1]],3),"\\nslope = ", round(dat_sum$statistic[[1]],3),"\\np-value = ", signif(dat_sum$p.value,3),sep="")
FST_genome_Cores <- round(dat_sum$estimate[[1]],3)

#plot(x=cg_df_Cores$Pre_FSt_sel, y=cg_df_Cores$Fitness)
dat_sum<-cor.test(x=cg_df_Cores$Pre_FSt_sel, y=cg_df_Cores$Fitness, method = "spearman")
paste("r = ",round(dat_sum$estimate[[1]],3),"\\nslope = ", round(dat_sum$statistic[[1]],3),"\\np-value = ", signif(dat_sum$p.value,3),sep="")
FST_causal_Cores <- round(dat_sum$estimate[[1]],3)

case <- data.frame(GF_off_genome_allEnv_Edges, 
                   GF_off_genome_causEnv_Edges, 
                   GF_off_causal_allEnv_Edges, 
                   GF_off_causal_causEnv_Edges, 
                   Env_ED_caus_env_Edges,
                   Env_MD_caus_env_Edges,
                   Env_ED_all_env_Edges,
                   Env_MD_all_env_Edges,
                   Env_ED_causP2_env_Edges,
                   Env_MD_causP2_env_Edges, 
                   FST_genome_Edges, 
                   FST_causal_Edges,
                   GF_off_genome_allEnv_Cores, 
                   GF_off_genome_causEnv_Cores, 
                   GF_off_causal_allEnv_Cores, 
                   GF_off_causal_causEnv_Cores, 
                   Env_ED_caus_env_Cores,
                   Env_MD_caus_env_Cores,
                   Env_ED_all_env_Cores,
                   Env_MD_all_env_Cores,
                   Env_ED_causP2_env_Cores,
                   Env_MD_causP2_env_Cores, 
                   FST_genome_Cores, 
                   FST_causal_Cores)

write.csv(case, paste("results/R_results/",seed,"_SL_corr_allEnv.csv",sep=""),row.names = F)

# write.table(Post_FST,file="~/Desktop/Post_FST_2.txt",sep=",",col.names = F,row.names = F,quote = F)
########################################################

#Get per population FST pre and post environmental shift
Pfst_pre_filt<-MakeDiploidFSTMat(SNPmat = Gt1f, locusNames = colnames(Gt1f), popNames = PopsALL)

#Filter out NA values
Pfst_pre_noNa<-Pfst_pre_filt[!is.na(Pfst_pre_filt$FST),]

#ink_bef<-data.frame(cor(Gt1_m2,Gt1f[,colnames(Gt1f)!=M2_MID],method="pearson"))

#Heterozgosity per allele before env. shift
Het_bef<-Pfst_pre_filt$He

#Calculate FST value for each allele pre environmental shift
F_ST_ll1<-Pfst_pre_noNa$T1/Pfst_pre_noNa$T2

#Calculate FST values averaged across each allele pre environmental shift
F_ST_l1<-mean(Pfst_pre_noNa$T1)/mean(Pfst_pre_noNa$T2)

# Gt2f<-t(G2f)
# rownames(Gt2f)<-PopsALL
# #colnames(Gt2f)<-MID_post_filt
# colnames(Gt2f)<-paste("M",position2_filt,sep="")
# #Calculate per locus FST values post environmental shift
# Pfst_post<-MakeDiploidFSTMat(SNPmat = Gt2f, locusNames = colnames(Gt2f), popNames = PopsALL)
# 
# #Filter out NA values
# Pfst_post_noNa<-Pfst_post[!is.na(Pfst_post$FST),]
# 
# #Heterozgosity per allele after env. shift
# Het_aft<-Pfst_post_noNa$He
# 
# #Calculate FST value for each allele pre environmental shift
# F_ST_ll2<-Pfst_post_noNa$T1/Pfst_post_noNa$T2
# 
# #Calculate FST values for each  allele pre environmental shift
# F_ST_l2<-mean(Pfst_post_noNa$T1)/mean(Pfst_post_noNa$T2)
# 
# #Pop_afreq2 is not filtered for MAF in order to properly compare all AF shifts from Pop_afreq1
# Pop_afreq2<-NULL
# for(i in 1:100){
#   Pop_afreq2<-rbind(Pop_afreq2,rowSums(G2[,Start[i]:Stop[i]])/(2*ncol(G2[,Start[i]:Stop[i]])))
# }
# 
# Pop_afreq2<-data.frame(Pop_afreq2)
# colnames(Pop_afreq2)<-paste("M",position2,sep="")


#stopCluster(cl)

##################################
#Population specific summary stats
##################################

#Location values
Loc <- NULL
for(j in 1:10){
  for(i in 1:10){
    Loc <- c(Loc,paste("A",i,sep=""))
  }
}
Loc<-factor(Loc,levels=Loc[1:10])

# #Population values
# Pop <- NULL
# for(i in 1:100){
#   Pop <- c(Pop,paste("P",i,sep=""))
# }

X <- NULL
for(j in 1:10){
  for(i in 1:10){
    X <- c(X,i)
  }
}

Y <- NULL
for(j in 1:10){
  Y<-c(Y,rep(j,10))
}

Env_before<-envPop$envPop

#Env_after<-envPop.shift$envSelect

#Diff_env<-Env_after-Env_before

#Env_range<- envPop$envSelect%in%round(envPop.shift$envSelect,1)&round(envPop.shift$envSelect,1)%in%envPop$envSelect

M2_AF_before<-data.frame(Pop_afreq1[,c(which(colnames(Pop_afreq1)==positionM2))])
colnames(M2_AF_before)<-"M2_AF_before"

# M2_AF_after<-data.frame(Pop_afreq2[,c(which(colnames(Pop_afreq2)==positionM2))])
# colnames(M2_AF_after)<-"M2_AF_after"

# M2_AF_diff<-M2_AF_after-M2_AF_before
# colnames(M2_AF_diff)<-"M2_AF_diff"

M1_AF_before_all<-Pop_afreq1[,-which(colnames(Pop_afreq1)==positionM2)]

#M1_AF_after_all<-Pop_afreq2[,-which(colnames(Pop_afreq2)==positionM2)]

#If neutral
#M1_AF_before_all<-Pop_afreq1
#M1_AF_after_all<-Pop_afreq2

# M1_AF_before_shared<-M1_AF_before_all[colnames(M1_AF_before_all)%in%colnames(M1_AF_after_all)]
# M1_AF_after_shared<-M1_AF_after_all[colnames(M1_AF_after_all)%in%colnames(M1_AF_before_all)]

# M1_AF_before<-data.frame(rowMeans(M1_AF_before_shared))
# colnames(M1_AF_before)<-"M1_AF_before"
# 
# M1_AF_after<-data.frame(rowMeans(M1_AF_after_shared))
# colnames(M1_AF_after)<-"M1_AF_after"
# 
# M1_AF_diff<-M1_AF_after-M1_AF_before
# colnames(M1_AF_diff)<-"M1_AF_diff"
# 
# F_ST_genome_bef.aft.<-data.frame(FST_genome_pop)
# colnames(F_ST_genome_bef.aft.)<-"F_ST_genome_bef.aft."

# F_ST_M2_bef.aft.<-data.frame(FST_M2_pop)
# colnames(F_ST_M2_bef.aft.)<-"F_ST_M2_bef.aft."

Rel_Fit_before <- data.frame(fitt[fitt$Type=="Fit",gen_nam[length(gen_nam)]])
colnames(Rel_Fit_before)<-"Rel_Fit_before"

# Rel_Fit_after <- data.frame(fitt[fitt$Type=="Fit",gen_nam[length(gen_nam)]])
# colnames(Rel_Fit_after)<-"Rel_Fit_after"

# Rel_Fit_diff<-Rel_Fit_after-Rel_Fit_before
# colnames(Rel_Fit_diff)<-"Rel_Fit_diff"

 Pop_size<-rep(specs$n,100)
 Summary_Pop<-cbind(Pop,X,Y,Env_before,GF_dist_gen,GF_dist_caus,CI,M2_AF_before,Rel_Fit_before)
 #Summary_Pop<-cbind(Pop,X,Y,Env_before,Env_after,Diff_env,Env_range,gfTrans1,gfTrans2,offset,gfM2Trans1,gfM2Trans2,M2offset,M2_AF_before,M2_AF_after,M2_AF_diff,M1_AF_before,M1_AF_after,M1_AF_diff,F_ST_genome_bef.aft.,F_ST_M2_bef.aft.,Rel_Fit_before,Rel_Fit_after,Rel_Fit_diff)
 
#If neutral:
# Pop_size<-unlist(rep(unname(specs[9:18]),10))
# Summary_Pop<-cbind(Pop,X,Y,Pop_size,Env_before,Env_after,Diff_env,Env_range,gfTrans1,gfTrans2,offset,gfM2Trans1=0,gfM2Trans2=0,M2offset=0,M2_AF_before=0,M2_AF_after=0,M2_AF_diff=0,M1_AF_before,M1_AF_after,M1_AF_diff,F_ST_genome_bef.aft.,F_ST_M2_bef.aft.=0,Rel_Fit_before,Rel_Fit_after,Rel_Fit_diff)

#Summary_Pop<-cbind(Pop,X,Y,Env_before,Env_after,Diff_env,Env_range,gfTrans1,gfTrans2,offset,gfM2Trans1=0,gfM2Trans2=0,M2offset=0,M2_AF_before=0,M2_AF_after=0,M2_AF_diff=0,M1_AF_before,M1_AF_after,M1_AF_diff,F_ST_genome_bef.aft.,F_ST_M2_bef.aft.=0,Rel_Fit_before,Rel_Fit_after,Rel_Fit_diff)
write.csv(Summary_Pop,file=paste("results/R_results/",seed,"_summary_Pop.csv",sep=""),row.names=F)
#Summary_Pop<-read.csv("Summary_Pop_1576675870126.csv")

##################################
# Allele specific summary stats
##################################

R2<-data.frame(gfMod$result) #Those allelese which had an R2 value > 0
colnames(R2)<-"R2"

R0<-colnames(alFreq[,!colnames(alFreq)%in%names(gfMod$result)]) #Get all alleles, regardless of R2 value
R0<-data.frame(rep(0,length(R0)),row.names = colnames(alFreq[,!colnames(alFreq)%in%names(gfMod$result)])) #Filter out those who we already have saved in R2
colnames(R0)<-"R2"
R2_all<-rbind(R2,R0) #merge them so we have all alleles with accompanying R2 values 
R2_all$Pos <- as.numeric(substring(row.names(R2_all),2))

Rho_Env<-cor(as.matrix(alFreq),envPop,method = "spearman")

temp<-merge(Rho_Env,R2_all,by=0)
colnames(temp)<-c("Locus","Rho_Env","R2","Position")

Linkage<-NULL
position1_filt_scaled<-NULL
for(i in 1:length(position1_filt)){
  if(position1_filt[i]>0 & position1_filt[i]<50001){
    Linkage<-c(Linkage,1)
    position1_filt_scaled<-c(position1_filt_scaled,position1_filt[i])
  }
  if(position1_filt[i]>50000 & position1_filt[i]<100001){
    Linkage<-c(Linkage,2)
    position1_filt_scaled<-c(position1_filt_scaled,position1_filt[i]-50000)
  }
  if(position1_filt[i]>100000 & position1_filt[i]<150001){
    Linkage<-c(Linkage,3)
    position1_filt_scaled<-c(position1_filt_scaled,position1_filt[i]-100000)
  }
  if(position1_filt[i]>150000 & position1_filt[i]<200001){
    Linkage<-c(Linkage,4)
    position1_filt_scaled<-c(position1_filt_scaled,position1_filt[i]-150000)
  }
  if(position1_filt[i]>200000 & position1_filt[i]<250001){
    Linkage<-c(Linkage,5)
    position1_filt_scaled<-c(position1_filt_scaled,position1_filt[i]-200000)
  }
  if(position1_filt[i]>250000 & position1_filt[i]<300001){
    Linkage<-c(Linkage,6)
    position1_filt_scaled<-c(position1_filt_scaled,position1_filt[i]-250000)
  }
  if(position1_filt[i]>300000 & position1_filt[i]<350001){
    Linkage<-c(Linkage,7)
    position1_filt_scaled<-c(position1_filt_scaled,position1_filt[i]-300000)
  }
  if(position1_filt[i]>350000 & position1_filt[i]<400001){
    Linkage<-c(Linkage,8)
    position1_filt_scaled<-c(position1_filt_scaled,position1_filt[i]-350000)
  }
  if(position1_filt[i]>400000 & position1_filt[i]<450001){
    Linkage<-c(Linkage,9)
    position1_filt_scaled<-c(position1_filt_scaled,position1_filt[i]-400000)
  }
  if(position1_filt[i]>450000 & position1_filt[i]<500001){
    Linkage<-c(Linkage,10)
    position1_filt_scaled<-c(position1_filt_scaled,position1_filt[i]-450000)
  }
}

#PreN$LG<-Linkage
#PreN[rownames(PreN)!=Pfst_pre_filt$LocusName,]
#PreN$FST<-Pfst_pre_filt$FST
#PreN$DistM2<-abs(PreN$PP-PreN$PP[PreN$LT=="M2"])

Summary_Locus<-cbind(temp[order(temp$Position),],Linkage)

write.csv(Summary_Locus,file=paste("results/R_results/",seed,"_summary_Loc.csv",sep=""),row.names=F)

##################################
# Simulation specific summary stats
##################################

#No_A<-unname(dim(vcf1_filt)[1])
PR2<-gfMod$species.pos.rsq/No_A
#PR2_UL
F_ST_l1
#F_ST_l2
Rel.Fit_l1<-mean(fitt[fitt$Type=="Fit",gen_nam[length(gen_nam)]])
#Rel.Fit_l2<-mean(fitt[fitt$Type=="Fit",gen_nam[length(gen_nam)]])
Rho_EnvR2<-cor(temp$Rho_Env,temp$R2,method = "spearman")

Summary_Sim<-cbind(seed,No_A,PR2,F_ST_l1,Rel.Fit_l1,Rho_EnvR2)
#Summary_Sim<-cbind(seed,No_A,PR2,F_ST_l1,F_ST_l2,Rel.Fit_l1,Rel.Fit_l2,Rho_EnvR2)

write.csv(Summary_Sim,file=paste("results/R_results/",seed,"_summary_Sim.csv",sep=""),row.names=F)

###############################################
#Visualize R^2>0 compared to Spearman correltation of alFreq to each env.
###############################################

#EnvCor<-cor(as.matrix(R2),envPop,method = "spearman")
#EnvCor<-data.frame(rownames(EnvCor),EnvCor)
#colnames(EnvCor)<-c("MID","rho")

#R2MID<-unique(data.frame(impDat$allele,impDat$r2))
#colnames(R2MID)<-c("MID","R2")

#Comp<-merge(EnvCor,R2MID, by="MID")

#Link<-NULL
#for(i in 1:length(Comp$MID)){
#  if(Comp$MID[i]%in%linked_MID){
#    Link<-c(Link,"Linked")
#    }
#  else{
#    Link<-c(Link,"Unlinked")
#  }
#}

#Comp$Link<-Link

#r2<-data.frame(gfMod$result)
#colnames(r2)<-"r2"

#save.image("~/Desktop/PostDoc/SLiMstuff/SLiM_output/250K/VCF_output/300K/1727520158823.RData")
#write.table(gfTrans1$C.Imp_genome_before,paste("results/R_results/",seed,"_CI",sep=""),row.names = F,col.names = F)
save.image(paste("results/R_results/",seed,".RData",sep=""))

# dat_sum<-cor.test(x=cg_df$D_CI_caus_causEnv, y=cg_df$Fitness, method = "spearman")
# plot(x=cg_df$D_CI_caus_causEnv, y=cg_df$Fitness)
# paste("rho = ",round(dat_sum$estimate[[1]],3),"\\np-value = ", signif(dat_sum$p.value,3),sep="")
# 
# 
# cg_subset<-cg_df[cg_df$Transplant=="T1" | cg_df$Transplant=="T5" | cg_df$Transplant=="T10",]
# cg_subset$Environment<-as.factor(rep(1:10,each=3,times=10))
# 
# dat_sum<-cor.test(x=cg_subset$D_CI_caus_allEnv, y=cg_subset$Fitness, method = "spearman")
# plot(x=cg_subset$D_CI_caus_causEnv, y=cg_subset$Fitness)
# ggplot(data=cg_subset)+
#   geom_point(aes(x=D_CI_caus_causEnv, y=Fitness, col=Transplant, shape=Environment))+
#   theme_classic()
# 
# paste("rho = ",round(dat_sum$estimate[[1]],3),"\\np-value = ", signif(dat_sum$p.value,3),sep="")
# 
# demes <-unique(cg_df$Transplant)
# CI_plots <- vector('list', length(demes))
# max(cg_df[!is.na(cg_df$D_CI_caus_causEnv),]$D_CI_caus_causEnv)
# min(cg_df[!is.na(cg_df$D_CI_caus_causEnv),]$D_CI_caus_causEnv)
# for (i in seq_along(demes)) {
#   dat_filt<-cg_df[cg_df$Transplant%in%demes[i],]
#   dat_sum<-cor.test(x=dat_filt$D_CI_caus_causEnv, y=dat_filt$Fitness, method = "pearson")
#   CI_plots[[i]] <- local({
#     i <- i
#     p1 <- ggplot(data=dat_filt)+
#       geom_smooth(aes(x=D_CI_caus_causEnv, y=Fitness), method="lm",se = F)+
#       annotate("text",x=0.025,y=0.9,label=paste("r = ",round(dat_sum$estimate[[1]],3),"\\np-value = ", signif(dat_sum$p.value,3),sep=""))+
#       geom_point(aes(x=D_CI_caus_causEnv, y=Fitness, col=Home), size=3)+
#       labs(shape="Home", colour="Common\\nGarden\\nLocation")+
#       scale_color_manual(values = colfunc(100), guide=F)+
#        ylab("Average relative fitness in common garden")+
#        xlab("GF Offset (genome)")+
#       #ylab("")+
#       #xlab("")+
#       #ylim(0,1)+
#       #xlim(0,0.05)+
#       ggtitle(paste("Transplant to",str_sub(demes[i],2,4)))+
#       theme_classic()+
#       theme(axis.text=element_text(size=12),axis.title=element_text(size=14)
#       )
#     #print(p1)
#   })
# }
# 
# plot_grid(CI_plots[[1]],CI_plots[[2]],CI_plots[[3]],CI_plots[[4]],CI_plots[[5]],CI_plots[[6]],CI_plots[[7]],CI_plots[[8]],CI_plots[[9]],CI_plots[[10]],
#           CI_plots[[11]],CI_plots[[12]],CI_plots[[13]],CI_plots[[14]],CI_plots[[15]],CI_plots[[16]],CI_plots[[17]],CI_plots[[18]],CI_plots[[19]],CI_plots[[20]],
#           CI_plots[[21]],CI_plots[[22]],CI_plots[[23]],CI_plots[[24]],CI_plots[[25]],CI_plots[[26]],CI_plots[[27]],CI_plots[[28]],CI_plots[[29]],CI_plots[[30]],
#           CI_plots[[31]],CI_plots[[32]],CI_plots[[33]],CI_plots[[34]],CI_plots[[35]],CI_plots[[36]],CI_plots[[37]],CI_plots[[38]],CI_plots[[39]],CI_plots[[40]],
#           CI_plots[[41]],CI_plots[[42]],CI_plots[[43]],CI_plots[[44]],CI_plots[[45]],CI_plots[[46]],CI_plots[[47]],CI_plots[[48]],CI_plots[[49]],CI_plots[[50]],
#           nrow=5)
# 
# plot_grid(CI_plots[[51]],CI_plots[[52]],CI_plots[[53]],CI_plots[[54]],CI_plots[[55]],CI_plots[[56]],CI_plots[[57]],CI_plots[[58]],CI_plots[[59]],CI_plots[[60]],
#           CI_plots[[66]],CI_plots[[62]],CI_plots[[63]],CI_plots[[64]],CI_plots[[65]],CI_plots[[66]],CI_plots[[67]],CI_plots[[68]],CI_plots[[69]],CI_plots[[70]],
#           CI_plots[[71]],CI_plots[[77]],CI_plots[[73]],CI_plots[[74]],CI_plots[[75]],CI_plots[[76]],CI_plots[[77]],CI_plots[[78]],CI_plots[[79]],CI_plots[[80]],
#           CI_plots[[81]],CI_plots[[82]],CI_plots[[88]],CI_plots[[84]],CI_plots[[85]],CI_plots[[86]],CI_plots[[87]],CI_plots[[88]],CI_plots[[89]],CI_plots[[90]],
#           CI_plots[[91]],CI_plots[[92]],CI_plots[[93]],CI_plots[[94]],CI_plots[[95]],CI_plots[[96]],CI_plots[[97]],CI_plots[[98]],CI_plots[[99]],CI_plots[[100]],
#           nrow=5)
# 
# 
# dat_sum<-cor.test(x=cg_df_Edges$D_CI_allEnv, y=cg_df_Edges$Fitness, method = "spearman")
# plot(x=cg_df_Edges$D_CI_allEnv, y=cg_df_Edges$Fitness)
# paste("rho = ",round(dat_sum$estimate[[1]],3),"\\np-value = ", signif(dat_sum$p.value,3),sep="")
# 
# dat_sum<-cor.test(x=cg_df_Edges$D_CI_causEnv, y=cg_df_Edges$Fitness, method = "spearman")
# plot(x=cg_df_Edges$D_CI_causEnv, y=cg_df_Edges$Fitness)
# paste("rho = ",round(dat_sum$estimate[[1]],3),"\\np-value = ", signif(dat_sum$p.value,3),sep="")
# 
# dat_sum<-cor.test(x=cg_df_Edges$D_CI_caus_allEnv, y=cg_df_Edges$Fitness, method = "spearman")
# plot(x=cg_df_Edges$D_CI_caus_allEnv, y=cg_df_Edges$Fitness)
# paste("rho = ",round(dat_sum$estimate[[1]],3),"\\np-value = ", signif(dat_sum$p.value,3),sep="")
# 
# dat_sum<-cor.test(x=cg_df_Edges$D_CI_caus_causEnv, y=cg_df_Edges$Fitness, method = "spearman")
# plot(x=cg_df_Edges$D_CI_caus_causEnv, y=cg_df_Edges$Fitness)
# paste("rho = ",round(dat_sum$estimate[[1]],3),"\\np-value = ", signif(dat_sum$p.value,3),sep="")
# 
# dat_sum<-cor.test(x=cg_df_Cores$D_CI_allEnv, y=cg_df_Cores$Fitness, method = "spearman")
# plot(x=cg_df_Cores$D_CI_allEnv, y=cg_df_Cores$Fitness)
# paste("rho = ",round(dat_sum$estimate[[1]],3),"\\np-value = ", signif(dat_sum$p.value,3),sep="")
# 
# dat_sum<-cor.test(x=cg_df_Cores$D_CI_causEnv, y=cg_df_Cores$Fitness, method = "spearman")
# plot(x=cg_df_Cores$D_CI_causEnv, y=cg_df_Cores$Fitness)
# paste("rho = ",round(dat_sum$estimate[[1]],3),"\\np-value = ", signif(dat_sum$p.value,3),sep="")
# 
# dat_sum<-cor.test(x=cg_df_Cores$D_CI_caus_allEnv, y=cg_df_Cores$Fitness, method = "spearman")
# plot(x=cg_df_Cores$D_CI_caus_allEnv, y=cg_df_Cores$Fitness)
# paste("rho = ",round(dat_sum$estimate[[1]],3),"\\np-value = ", signif(dat_sum$p.value,3),sep="")
# 
# dat_sum<-cor.test(x=cg_df_Cores$D_CI_caus_causEnv, y=cg_df_Cores$Fitness, method = "spearman")
# plot(x=cg_df_Cores$D_CI_caus_causEnv, y=cg_df_Cores$Fitness)
# paste("rho = ",round(dat_sum$estimate[[1]],3),"\\np-value = ", signif(dat_sum$p.value,3),sep="")
# 
# #Visualize environments
# ggplot()+
#   geom_point(aes(cg_df_Edges$D_CI_allEnv,cg_df_Edges$Pre_FST),col="red")+
#   geom_point(aes(cg_df_Edges$D_CI_causEnv,cg_df_Edges$Pre_FST),col="darkred")+
#   theme_classic()
# 
# ggplot()+
#   geom_point(aes(cg_df_Edges$D_CI_allEnv,cg_df_Edges$Pre_FST),col="red")+
#   geom_point(aes(cg_df_Edges$D_CI_causEnv,cg_df_Edges$Pre_FST),col="darkred")+
#   theme_classic()
# 
# ggplot()+  
#   geom_point(aes(cg_df_Edges$D_CI_allEnv,cg_df_Edges$Pre_FSt_sel),col="red")+
#   geom_point(aes(cg_df_Edges$D_CI_causEnv,cg_df_Edges$Pre_FSt_sel),col="darkred")+
#   theme_classic()
# 
# ggplot()+  
#   geom_point(aes(cg_df_Edges$D_CI_caus_allEnv,cg_df_Edges$Pre_FST),size=2,col="orange")+
#   geom_point(aes(cg_df_Edges$D_CI_caus_causEnv,cg_df_Edges$Pre_FST),col="brown")+
#   theme_classic()
# 
# ggplot()+ 
#   geom_point(aes(cg_df_Edges$D_CI_caus_allEnv,cg_df_Edges$Pre_FSt_sel),size=2,col="orange")+
#   geom_point(aes(cg_df_Edges$D_CI_caus_causEnv,cg_df_Edges$Pre_FSt_sel),col="brown")+
#   theme_classic()
# 
# plot(cg_df_Cores$D_CI_allEnv,cg_df_Cores$Pre_FST)
# plot(cg_df_Cores$D_CI_allEnv,cg_df_Cores$Pre_FSt_sel)
# plot(cg_df_Cores$D_CI_causEnv,cg_df_Cores$Pre_FST)
# plot(cg_df_Cores$D_CI_causEnv,cg_df_Cores$Pre_FSt_sel)
# plot(cg_df_Cores$D_CI_caus_allEnv,cg_df_Cores$Pre_FST)
# plot(cg_df_Cores$D_CI_caus_allEnv,cg_df_Cores$Pre_FSt_sel)
# plot(cg_df_Cores$D_CI_caus_causEnv,cg_df_Cores$Pre_FST)
# plot(cg_df_Cores$D_CI_caus_causEnv,cg_df_Cores$Pre_FSt_sel)
# 
# ggplot()+
#   geom_point(aes(cg_df_Edges$D_CI_allEnv,cg_df_Edges$Fitness),col="red")+
#   theme_classic()
# 
# ggplot()+
#   geom_point(aes(cg_df_Cores$D_CI_causEnv,cg_df_Cores$Fitness),col="darkred")+
#   theme_classic()
# 
# ggplot()+
#   geom_point(aes(cg_df_Edges$D_CI_caus_allEnv,cg_df_Edges$Fitness),col="red")+
#   theme_classic()
# 
# ggplot()+
#   geom_point(aes(cg_df_Cores$D_CI_caus_causEnv,cg_df_Cores$Fitness),col="darkred")+
#   theme_classic()
