StMichel



#RSCRIPT TO CALCULATE BETAS WITH BOOTSTRAPPING

####packages
require(tidyr)
require(dplyr)
library(vcfR)
library(adegenet)
require(hierfstat)


#set wd where I have the files
setwd("/scratch/wally/FAC/FBM/DEE/mchapuis/default/afontcub/POPGEN")

##vcf
selysi.asocial4<-read.vcfR("selysi_compromise.Asocial.vcf.gz")


## Metadata for asocial compromise 152 indvs
geod_filt4<-read.table('metadata_good_filt_Compromise.csv', header=T, sep=",", dec= ".")


#create genind object
asocial4.gd<- vcfR2genind(selysi.asocial4)
## Define Pop
pop(asocial4.gd)<-geod_filt4$Pop_social

#SAmpling pops
dplyr::count(geod_filt4,Pop_social)
# I will hve to remvoe (n =1) Aubenas_P, Dalaas_M, Leuk_P, Riddes_M. 
#and be careful with (n=3) Bussets_P, Hauderes_P, Luette_P, Riddes_P, 

##filter small subpops
asocial4.gd2<-asocial4.gd[!(asocial4.gd@pop=="Aubenas_P"|asocial4.gd@pop=="Dalaas_M"|asocial4.gd@pop=="Leuk_P"|asocial4.gd@pop=="Riddes_M")]# rm pops < 3
#asocial4.gd2<-asocial4.gd[!(asocial4.gd@pop=="Aubenas_P"|asocial4.gd@pop=="Dalaas_M"|asocial4.gd@pop=="Leuk_P"|asocial4.gd@pop=="Riddes_M"|asocial4.gd@pop=="Bussets_P"|asocial4.gd@pop=="Hauderes_P"|asocial4.gd@pop=="Luette_P"|asocial4.gd@pop=="Riddes_P")]#rm pops < 5
geno<-genind2hierfstat(asocial4.gd2) #convert to hierfstat
dim(geno)

#filter metadata
meta<-geod_filt4[!(geod_filt4$Pop_social=="Aubenas_P"|geod_filt4$Pop_social=="Dalaas_M"|geod_filt4$Pop_social=="Leuk_P"|geod_filt4$Pop_social=="Riddes_M"),]%>%select( seqIDvcf, Population, Pop_social)
#meta<-geod_filt4[!(geod_filt4$Pop_social=="Aubenas_P"|geod_filt4$Pop_social=="Dalaas_M"|geod_filt4$Pop_social=="Leuk_P"|geod_filt4$Pop_social=="Riddes_M"|geod_filt4$Pop_social=="Bussets_P"|geod_filt4$Pop_social=="Hauderes_P"|geod_filt4$Pop_social=="Luette_P"|geod_filt4$Pop_social=="Riddes_P"),]%>%select( seqIDvcf, Population, Pop_social)
names(meta)<-c("ID", "pop","subpop"); dim(meta)

geno$pop<-meta$subpop#make sure pop factor in hierfstat object is pop social



#FUNCTION to get matrix of rarefaction sample sizes
pops <- unique(meta$pop)
subpops <- unique(meta$subpop)
pops_combinations <- t(combn(pops,2))
pops_resample_sizes <- matrix(data=NA,nrow=length(subpops),ncol=length(subpops))
rownames(pops_resample_sizes) <- subpops
colnames(pops_resample_sizes) <- subpops

for(row in 1:nrow(pops_combinations)){
  pop1 <- pops_combinations[row,1]
  pop2 <- pops_combinations[row,2]

  pops_combination_subset <- as.data.frame(subset(meta, pop == pop1 | pop == pop2))

  subpop_sizes <- c()
  subpop_names <- unique(pops_combination_subset$subpop)
  for(subpop_i in subpop_names){
    subpop_subset <- subset(pops_combination_subset, subpop == subpop_i)
    n <- nrow(subpop_subset)
    subpop_sizes <- c(subpop_sizes, n)
  }

  n_min <- min(subpop_sizes)

  subpops_combinations <- t(combn(subpop_names,2))
  for(row in 1:nrow(subpops_combinations)){
    subpop1 <- as.character(subpops_combinations[row,1])
    subpop2 <-as.character(subpops_combinations[row,2])
    pops_resample_sizes[subpop1,subpop2] <- n_min
    pops_resample_sizes[subpop2,subpop1] <- n_min
  }
}

print(pops_resample_sizes)

####FUNCTION to Calculate betas  with BOOTSRTAP iterations
#(i changed slightly the resampling function using base R instead of dplyr but does the same)
#test_geno<-geno[,1:21]
#geno<-test_geno #removed once tested
nb_ite=100# nb of iterations
subpop_names<-unique(geno$pop)
subpop_combis=t(combn(subpop_names,2))

pairfst <- matrix(data=NA,nrow=length(subpop_names),ncol=length(subpop_names))#matrix for mean beta
rownames(pairfst) <- subpop_names; colnames(pairfst) <- subpop_names
pairSD <- matrix(data=NA,nrow=length(subpop_names),ncol=length(subpop_names)) #matrix for SD of beta
rownames(pairSD) <- subpop_names; colnames(pairSD) <- subpop_names
pairMedianFst <- matrix(data=NA,nrow=length(subpop_names),ncol=length(subpop_names)) #matrix for SD of beta
rownames(pairMedianFst) <- subpop_names; colnames(pairMedianFst) <- subpop_names

#nrow=25

for (nrow in 1:nrow(subpop_combis)) {
    subpop1<-as.character(subpop_combis[nrow,1])
    subpop2<-as.character(subpop_combis[nrow,2])
      subpop_subset<-geno[(geno$pop==subpop1|geno$pop==subpop2),]
      n<-pops_resample_sizes[subpop1,subpop2]
      beta_vec<-NULL
      #with manual iteration (bootrstraping)
      for (iteration in 1:nb_ite){
        subpop_subset$pop<-factor(subpop_subset$pop,levels=as.character(unique(subpop_subset$pop)))#set correct nb of levels
        resample<-lapply(split(subpop_subset,subpop_subset$pop),function (x) x[sample(1:nrow(x),n),])
        resample_subset<-do.call("rbind",resample) #unlist and make 1 df
        beta<-betas(resample_subset)$betaW #calculate beta
        beta_vec<-c(beta_vec,beta)#store in vector
        }
      mean_beta<-mean(beta_vec,na.rm=T)#take the mean across nb_ite iterations
      sd_beta<-sd(beta_vec,na.rm=T)
      median_beta<-median(beta_vec,na.rm=T)#take median

      pairfst[subpop1,subpop2]<- pairfst[subpop2,subpop1]<-mean_beta
      pairSD[subpop1,subpop2] <-  pairSD[subpop2,subpop1]<-sd_beta
      pairMedianFst[subpop1,subpop2]<-pairMedianFst[subpop2,subpop1]<-median_beta

      print(paste(nrow,as.character(subpop_combis[nrow,])))
      print(mean_beta)
      print(sd_beta)
      print(median_beta)

    }

print(pairfst)
print(pairSD)
print(pairMedianFst)


sink("Boot_MEANbeta_asocial4_Wg17_Popsocial.txt")
print(pairfst)
sink()

sink("Boot_SDbeta_asocial4_Wg17_Popsocial.txt")
print(pairSD)
sink()

sink("Boot_MEDIANbeta_asocial4_Wg17_Popsocial.txt")
print(pairMedianFst)
sink()
