#library.path <- c("C:/Users/18084/Documents/R/win-library/3.6","E:/R-3.6.2/library")
library.path <- .libPaths()

require(OutFLANK, lib.loc=library.path)
require(vcfR, lib.loc=library.path)
require(adegenet, lib.loc=library.path)
require(gdm, lib.loc=library.path)
require(gradientForest, lib.loc=library.path)
require(foreach, lib.loc=library.path)
require(doParallel, lib.loc=library.path)
require(pbapply, lib.loc=library.path)
require(gdata, lib.loc=library.path)
require(data.table, lib.loc=library.path)
require(PresenceAbsence, lib.loc=library.path)
require(ROCR, lib.loc=library.path)
require(modEvA, lib.loc=library.path)
require(dplyr, lib.loc=library.path)
require(grid, lib.loc=library.path)
require(gridExtra, lib.loc=library.path)
require(gtools, lib.loc=library.path)
require(stringr, lib.loc=library.path)
require(reshape2, lib.loc=library.path)
require(hierfstat)

#setwd("/Users/akijarl/Desktop/PostDoc/TTT_Offset_Vulnerability_GF_Sims/")
#setwd("/Users/akijarl/Desktop/TTT_Offset_Vulnerability_GF_Sims/")
setwd("E:/Research_AJL/TTT_Offset_Vulnerability_GF_Sims/")
#seed=commandArgs(trailingOnly=T)

# seed = 3393663285300
# seed = 2207643041387
# seed = 6359250998124
# seed = 1940819310024
# seed = 6040998529642
# seed = 6417941976338
# seed = 3323115928098
# seed = 6232745500281
# seed = 5643197868019
# seed = 1278335993644

# seed = 2519699755125 
# seed = 5518873473655 
# seed = 5536135188855 
# seed = 5555300357555 
# seed = 5597604322855 
# seed = 5526691746355 
# seed = 5561122192855 
# seed = 5551393883655 
# seed = 5521073151055 
# seed = 5583941596055 

# seed = 1589041996747 
# seed = 5385408036410 
# seed = 6842444356530 
# seed = 4621752902865 
# seed = 3809137545526 
# seed = 8081899648532 
# seed = 2876069253310 
# seed = 4280735488981 
# seed = 1254488294944 
# seed = 4698739288374 

#Read in all files with the seed values listed in 'seed_r.txt' and run them through GradientForest
r.seeds <- read.table("seeds_Neutral.txt")
options(scipen = 999)
for(k in 1:length(r.seeds$V1)){
  paste("seed = ",r.seeds$V3[k])

#f1<-list.files("results/SLiM_output/Sim_sum")[grep(".txt",list.files("results/SLiM_output/Sim_sum"))][5]
#(seed<-substr(f1, start=1, stop=13))
#fit<-read.table(paste("Fit_SP_100_",seed,".txt",sep=""),fill=T) 
fit<-read.table(paste("results/SLiM_output/Sim_sum/",r.seeds$V3[k],"_Freq_NL.txt",sep=""))

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

 colnames(fit)<-c("m","n1","n2","n3","n4","n5","n6","n7","n8","n9","n10","u","r","Env_rate","Burnin","Env_shift", "Generation", fit_nam, freq_nam, env_nam)


N<-data.frame(sum(fit$n1[1],fit$n2[1],fit$n3[1],fit$n4[1],fit$n5[1],fit$n6[1],fit$n7[1],fit$n8[1],fit$n9[1],fit$n10[1])*10)
colnames(N)<-"N"
(specs<-data.frame(r.seeds$V1[k],N,fit[1,c(1,12:16)],fit[1,c(2:11)]))

#write.table(specs,"R_results/output_metadata.txt",append=F,quote=F,sep=",",row.names=F)

plotTitle <- paste(colnames(specs)[1],":",specs[[1]],", ", colnames(specs)[2],":",specs[2],", ",colnames(specs)[3],":",specs[3],", ",colnames(specs)[4],":",1e-7,",",colnames(specs)[5],":",specs[5],", ",colnames(specs)[6],":",specs[6],",\\n ",colnames(specs)[7],":",specs[7],", ",colnames(specs)[8],":",specs[8],sep="")

gen_nam <- paste("Gen",fit$Generation,sep="")

fitt<-data.frame(t(fit[,-1:-17])) # For neutral simulation

colnames(fitt)<-gen_nam

fitt$Location <- as.factor(rep(paste("A",seq(1,10,1),sep=""),30))
fitt$Location <- factor(fitt$Location, levels = unique(fitt$Location))

fitt$Type <- as.factor(c(rep("Fit",100),rep("Freq",100),rep("Env",100)))

#VCF files are filtered with vcftools, as it is much faster than R. The filtering for MAF > 0.05 is accomplished with the following code:
# vcf1 <- read.vcfR(paste("results/SLiM_output/VCF_files/T1_",seed,"_unfiltered_subset.recode.vcf",sep=""))
# geno1 <- vcf1@gt[,-1] # Remove 1st column, which is 'Format'
# position1 <- getPOS(vcf1) # Positions in bp
# chromosome1 <- getCHROM(vcf1) # Chromosome information
# 
# rm(vcf1)
# gc()
# 
# vcf2 <- read.vcfR(paste("results/SLiM_output/VCF_files/T2_",seed,"_unfiltered_subset.recode.vcf",sep=""))
# geno2 <- vcf2@gt[,-1] # Remove 1st column, which is 'Format'
# position2 <- getPOS(vcf2) # Positions in bp
# chromosome2 <- getCHROM(vcf2) # Chromosome information
# 
# rm(vcf2)
# gc()

vcf1_filt <- read.vcfR(paste("results/SLiM_output/VCF_files/T1_",r.seeds$V3[k],"_filtered_subset.recode.vcf",sep=""))
geno1_filt <- vcf1_filt@gt[,-1] # Remove 1st column, which is 'Format'
position1_filt <- as.numeric(getPOS(vcf1_filt)) # Positions in bp
chromosome1_filt <- as.numeric(getCHROM(vcf1_filt)) # Chromosome information

No_A<-unname(dim(vcf1_filt)[1]) #Get the number of filtered alleles

rm(vcf1_filt)

# vcf2_filt <- read.vcfR(paste("results/SLiM_output/VCF_files/T2_",seed,"_filtered_subset.recode.vcf",sep=""))
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

#Create Genotype matrix
# G1 <- matrix(NA, nrow = nrow(geno1), ncol = ncol(geno1))
# G1[geno1 %in% c("0/0", "0|0")] <- 0
# G1[geno1  %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
# G1[geno1 %in% c("1/1", "1|1")] <- 2

#Start<-seq(1,ncol(geno1_filt),100)
#Stop<-seq(100,ncol(geno1_filt),100)

#Start and Stop for subsampled vcf file
Start<-seq(1,ncol(geno1_filt),10)
Stop<-seq(10,ncol(geno1_filt),10)

#No subsample (neutral / size variation)
# St=1
# demes=as.numeric(unname(specs[,9:18]/10))
# for(i in 1:9){
#   St=c(St,St[i]+demes[i])
# }
# 
# Start<-NULL
# for(i in seq(0,900,100)){
#   Start<-c(Start,St+i)
# }
# 
# Sp=demes[1]
# for(i in 1:9){
#   Sp=c(Sp,Sp[i]+demes[i+1])
# }
# 
# Stop<-NULL
# for(i in seq(0,900,100)){
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
envPop<-data.frame(fitt[fitt$Type=="Env",gen_nam[length(gen_nam)-30]]) #300 years prior to the end of the simulation is taken as the "before environmental shift" time
names(envPop) <- "envSelect"

# envPop.shift<-data.frame(fitt[fitt$Type=="Env",gen_nam[length(gen_nam)]])
# names(envPop.shift) <- "envSelect"

#Merge the population specific allele frequencies of all neutral (M1) alleles with the population specific frequency of the selected (M2) allele
alFreq<-data.frame(cbind(Pop_afreq1))
#alFreq<-Pop_afreq1

##############################################
# Chunk to fit GF models to minor allele frequencies at the level of
# populations
# GF is fit to each SNP individually to 
# ease computational / memory burden

#cl <- makeCluster(cores)
#registerDoParallel(cl)

##### added by MCF, running all loci in on model ##########

gfMod <- gradientForest(data=data.frame(envPop, alFreq),
                        predictor.vars=colnames(envPop),
                        response.vars=colnames(alFreq),
                        corr.threshold=0.5, 
                        ntree=500, 
                        trace=T)

#stopCluster(cl)

# Calculate genomic offset
# note that I am doing this for the avearge across all alleles since 
# GF was fit to all alleles simultaneously
# The more correct way is to calculate offset for adaptive alleles only,
# either individually or for a model fit to just those alleles.
gfTrans1 <- predict(gfMod, envPop)
colnames(gfTrans1)<-"C.Imp_genome_before"


CI<-NULL
for(i in 1:nrow(gfTrans1)){
  for(j in 1:nrow(gfTrans1)){
    CI<-c(CI,dist(rbind(gfTrans1[i,],gfTrans1[j,])))
    #CI<-c(CI,CI_bf[j]-CI_bf[i])
  }
}

#plot(gfMod, plot.type="Overall.Importance")
#plot(gfMod, plot.type="C", show.species=T)

#gfTrans2 <- predict(gfMod, envPop.shift)
#colnames(gfTrans2)<-"C.Imp_genome_after"

# offset needs to be considered using absolute values ()
#offset <- gfTrans2-gfTrans1
#colnames(offset)<-"D_C.Imp_genome"

##############################################################################
#Get  Weir & Cockerham F_ST values from the VCF files and use population data

#Create an object listing every population in the whole dataset
# PopsALL <- NULL
# for(j in rep(1:100)){
#   for(i in rep(j,100)){
#     PopsALL <- c(PopsALL,i)
#   }
# }
PopsALL <- NULL
for(j in rep(1:100)){
  for(i in rep(j,10)){
    PopsALL <- c(PopsALL,i)
  }
}

#Create an object splitting a single population into a Pre ("T1") and Post ("T2") "population"
# PopsP <- c(rep("T1",100),rep("T2",100))
PopsP <- c(rep("T1",10),rep("T2",10))

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
# Gt2<-t(G2)
# rownames(Gt2)<-PopsALL
# #colnames(Gt2)<-MID_post
# colnames(Gt2)<-paste("M",position2,sep="")
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

########################################################
#Calculate per locus FST values pre environmental shift for Edge and cCore populations

#Whole genome
Gt1f<-t(G1f)
colnames(Gt1f)<-paste("M",position1_filt,sep="")

Pre_geno <- data.frame(PopsALL,Gt1f)
colnames(Pre_geno)[1]<-"Locality"

Edges <- c(1,5,6,10,12,19,41,50,51,60,82,89,91,95,96,100)

Pre_geno_Edge <- Pre_geno[Pre_geno$Locality%in%Edges,]

start_time <- Sys.time()
Pre_FST_Edge<-pairwise.WCfst(Pre_geno_Edge,diploid = T)
end_time <- Sys.time()
print(end_time - start_time)

mean(Pre_FST_Edge, na.rm=T)

Cores <- c(34, 35, 36, 37, 44, 45, 46, 47, 54, 55, 56, 57, 64, 65, 66, 67)

Pre_geno_Cores <- Pre_geno[Pre_geno$Locality%in%Cores,]

start_time <- Sys.time()
Pre_FST_Cores<-pairwise.WCfst(Pre_geno_Cores,diploid = T)
end_time <- Sys.time()
print(end_time - start_time)

mean(Pre_FST_Cores, na.rm=T)

Pre_FST_Cores[is.na(Pre_FST_Cores)]<-0
Pre_FST_Edge[is.na(Pre_FST_Edge)]<-0

cg_df <- data.frame(rep(1:nrow(gfTrans1),each=100),rep(1:nrow(gfTrans1),100),CI)
colnames(cg_df)<-c("Home","Transplant","D_CI")


cg_df_Edges <- cg_df[as.numeric(cg_df$Home)%in%Edges,]
cg_df_Edges <- cg_df_Edges[as.numeric(cg_df_Edges$Transplant)%in%Edges,]

cg_df_Edges$FST <- as.vector(Pre_FST_Edge)

cg_df_Cores <- cg_df[as.numeric(cg_df$Home)%in%Cores,]
cg_df_Cores <- cg_df_Cores[as.numeric(cg_df_Cores$Transplant)%in%Cores,]

cg_df_Cores$FST <- as.vector(Pre_FST_Edge)


dat_sum<-cor.test(x=cg_df_Edges$D_CI, y=cg_df_Edges$FST, method = "spearman")
paste("r = ",round(dat_sum$estimate[[1]],3),"\\nslope = ", round(dat_sum$statistic[[1]],3),"\\np-value = ", signif(dat_sum$p.value,3),sep="")
GF_off_genome_Edges <- round(dat_sum$estimate[[1]],3)

dat_sum<-cor.test(x=cg_df_Cores$D_CI, y=cg_df_Cores$FST, method = "spearman")
paste("r = ",round(dat_sum$estimate[[1]],3),"\\nslope = ", round(dat_sum$statistic[[1]],3),"\\np-value = ", signif(dat_sum$p.value,3),sep="")
GF_off_genome_Cores <- round(dat_sum$estimate[[1]],3)

x <- melt(data.frame(GF_off_genome_Edges,GF_off_genome_Cores))

ggplot(x, aes(variable, value, fill=variable)) + 
  geom_bar(position="dodge",stat="identity")+
  theme_classic()+
  ylab("Spearman's rho \\nbetween FST and GF Offset")+ 
  theme(legend.position = "none")

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
# cores <- 7
# cl <- makeCluster(cores)
# registerDoParallel(cl)
# 
# start_time <- Sys.time()
# Pre_FST<-pairwise.WCfst(Pre_geno,diploid = T)
# end_time <- Sys.time()
# print(paste("Run time:",end_time - start_time))
# stopCluster(cl)
# 
# start_time <- Sys.time()
# Post_FST<-pairwise.WCfst(Post_geno,diploid = T)
# end_time <- Sys.time()
# print(paste("Run time:",end_time - start_time))
# 
# write.table(Pre_FST,file="~/Desktop/Pre_FST_2.txt",sep=",",col.names = F,row.names = F,quote = F)
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
#colnames(Gt2f)<-MID_post_filt
# colnames(Gt2f)<-paste("M",position2_filt,sep="")
#Calculate per locus FST values post environmental shift
# Pfst_post<-MakeDiploidFSTMat(SNPmat = Gt2f, locusNames = colnames(Gt2f), popNames = PopsALL)

#Filter out NA values
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

#Pop_afreq2 is not filtered for MAF in order to properly compare all AF shifts from Pop_afreq1
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

#Population values
Pop <- NULL
for(i in 1:100){
  Pop <- c(Pop,paste("P",i,sep=""))
}

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

Env_before<-envPop$envSelect

#Env_after<-envPop.shift$envSelect

#Diff_env<-Env_after-Env_before

#Env_range<- envPop$envSelect%in%round(envPop.shift$envSelect,1)&round(envPop.shift$envSelect,1)%in%envPop$envSelect

# M2_AF_before<-data.frame(Pop_afreq1[,c(which(colnames(Pop_afreq1)==positionM2))])
# colnames(M2_AF_before)<-"M2_AF_before"
# 
# M2_AF_after<-data.frame(Pop_afreq2[,c(which(colnames(Pop_afreq2)==positionM2))])
# colnames(M2_AF_after)<-"M2_AF_after"
# 
# M2_AF_diff<-M2_AF_after-M2_AF_before
# colnames(M2_AF_diff)<-"M2_AF_diff"

# M1_AF_before_all<-Pop_afreq1[,-which(colnames(Pop_afreq1)==positionM2)]
# M1_AF_after_all<-Pop_afreq2[,-which(colnames(Pop_afreq2)==positionM2)]
# 
#If neutral
M1_AF_before_all<-Pop_afreq1
#M1_AF_after_all<-Pop_afreq2

# M1_AF_before_shared<-M1_AF_before_all[colnames(M1_AF_before_all)%in%colnames(M1_AF_after_all)]
# M1_AF_after_shared<-M1_AF_after_all[colnames(M1_AF_after_all)%in%colnames(M1_AF_before_all)]

# M1_AF_before<-data.frame(rowMeans(M1_AF_before_shared))
# colnames(M1_AF_before)<-"M1_AF_before"

# M1_AF_after<-data.frame(rowMeans(M1_AF_after_shared))
# colnames(M1_AF_after)<-"M1_AF_after"
# 
# M1_AF_diff<-M1_AF_after-M1_AF_before
# colnames(M1_AF_diff)<-"M1_AF_diff"

# F_ST_genome_bef.aft.<-data.frame(FST_genome_pop)
# colnames(F_ST_genome_bef.aft.)<-"F_ST_genome_bef.aft."

# F_ST_M2_bef.aft.<-data.frame(FST_M2_pop)
# colnames(F_ST_M2_bef.aft.)<-"F_ST_M2_bef.aft."

Rel_Fit_before <- data.frame(fitt[fitt$Type=="Fit",gen_nam[length(gen_nam)-30]])
colnames(Rel_Fit_before)<-"Rel_Fit_before"

# Rel_Fit_after <- data.frame(fitt[fitt$Type=="Fit",gen_nam[length(gen_nam)]])
# colnames(Rel_Fit_after)<-"Rel_Fit_after"
# 
# Rel_Fit_diff<-Rel_Fit_after-Rel_Fit_before
# colnames(Rel_Fit_diff)<-"Rel_Fit_diff"

# Pop_size<-rep(specs$n,100)
Summary_Pop<-cbind(Pop,X,Y,Env_before,Env_after,Diff_env,Env_range,gfTrans1,gfTrans2,offset,gfM2Trans1,gfM2Trans2,M2offset,M2_AF_before,M2_AF_after,M2_AF_diff,M1_AF_before,M1_AF_after,M1_AF_diff,F_ST_genome_bef.aft.,F_ST_M2_bef.aft.,Rel_Fit_before,Rel_Fit_after,Rel_Fit_diff)

#If neutral:
Pop_size<-unlist(rep(unname(specs[9:18]),10))
#Summary_Pop<-cbind(Pop,X,Y,Pop_size,Env_before,Env_after,Diff_env,Env_range,gfTrans1,gfTrans2,offset,gfM2Trans1=0,gfM2Trans2=0,M2offset=0,M2_AF_before=0,M2_AF_after=0,M2_AF_diff=0,M1_AF_before,M1_AF_after,M1_AF_diff,F_ST_genome_bef.aft.,F_ST_M2_bef.aft.=0,Rel_Fit_before,Rel_Fit_after,Rel_Fit_diff)
#Summary_Pop<-cbind(Pop,X,Y,Pop_size,Env_before,Env_after,Diff_env,Env_range,gfTrans1,M1_AF_before,Rel_Fit_before,Rel_Fit_after,Rel_Fit_diff)

#Summary_Pop<-cbind(Pop,X,Y,Env_before,gfTrans1,M1_AF_before_all,Rel_Fit_before,)
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

write.csv(Summary_Locus,file=paste("results/R_results/",r.seeds$V1[k],"_summary_Loc.csv",sep=""),row.names=F)

##################################
# Simulation specific summary stats
##################################

#No_A<-unname(dim(vcf1_filt)[1])
PR2<-gfMod$species.pos.rsq/No_A
#PR2_UL
F_ST_l1
#F_ST_l2
Rel.Fit_l1<-mean(fitt[fitt$Type=="Fit",gen_nam[length(gen_nam)-30]])
#Rel.Fit_l2<-mean(fitt[fitt$Type=="Fit",gen_nam[length(gen_nam)]])
Rho_EnvR2<-cor(temp$Rho_Env,temp$R2,method = "spearman")

#Summary_Sim<-cbind(seed,No_A,PR2,F_ST_l1,F_ST_l2,Rel.Fit_l1,Rel.Fit_l2,Rho_EnvR2)
Summary_Sim<-cbind(r.seeds$V1[k],No_A,PR2,F_ST_l1,Rel.Fit_l1,Rho_EnvR2)

write.csv(Summary_Sim,file=paste("results/R_results/",r.seeds$V1[k],"_summary_Sim.csv",sep=""),row.names=F)

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
write.table(gfTrans1$C.Imp_genome_before,paste("results/R_results/",r.seeds$V1[k],"_CI",sep=""),row.names = F,col.names = F)
save.image(paste("results/R_results/",r.seeds$V1[k],".RData",sep=""))

#rm(list=ls())
}
