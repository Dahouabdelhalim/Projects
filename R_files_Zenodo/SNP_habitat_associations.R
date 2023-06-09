#SNP analysis 
#We want to know if microhabitat where fish were captured (like depth and flow rate) are predicted by SNPs
#this tells us if fish have habitat preferences influenced by genotype

#load packages 
{
  BiocManager::install("qvalue", version = "3.8")
  library(AER)
  library(MASS)
  library(qvalue) #for pval corrections
  library(VGAM)
  library(tidyr)
  library(ggplot2)
  library(car)
  library(tidyverse)
  }


#load in data
metadata<-read.csv("a_Metadata.pooledgenotyping_80pcompDP8.csv") # shows sequence, habitat, watershed
genotypes<-read.csv("RAD_genotypes_reordered copy.csv") # this is the full file- 840:78225
trap<-read.csv("07.ecoloTrap.csv") # fish IDs as rows 
RADseq<-read.csv("06.ddRADmeta.csv") # has fishID and fishID universal


#eliminates duplicates in trap data
trap <- trap[order(trap$fishID.univ),]
trap <- trap[-c(seq(1,39, by = 2)),]

#add in vegetation and substrate (to trap)
substrate.matrix <- data.frame(A =  grepl("A", trap$substrate),  B =  grepl("B", trap$substrate), C =  grepl("C", trap$substrate),
                               D =  grepl("D", trap$substrate),  E =  grepl("E", trap$substrate),  F =  grepl("F", trap$substrate),
                               G =  grepl("G", trap$substrate),  H =  grepl("H", trap$substrate),  I =  grepl("I", trap$substrate), 
                               J =  grepl("J", trap$substrate))
trap$substratePCA1 <- prcomp(substrate.matrix, scale = T)$x[,1] # First 2 PCA axes
trap$substratePCA2 <- prcomp(substrate.matrix, scale = T)$x[,2] # First 2 PCA axes

vegetation.matrix <- data.frame(A =  grepl("A", trap$veg_structure),  D =  grepl("D", trap$veg_structure), E =  grepl("E", trap$veg_structure),
                                F =  grepl("F", trap$veg_structure),  Bh =  grepl("Bh", trap$veg_structure),  Cg =  grepl("Cg", trap$veg_structure),
                                Bg =  grepl("Bg", trap$veg_structure),  Blp =  grepl("Blp", trap$veg_structure),  Clp =  grepl("Clp", trap$veg_structure), 
                                Cros =  grepl("Cros", trap$veg_structure), Bspp =  grepl("Bspp", trap$veg_structure), Bsc =  grepl("Bsc", trap$veg_structure), 
                                Calg =  grepl("Calg", trap$veg_structure), Ch =  grepl("Ch", trap$veg_structure))
trap$vegetationPCA1 <- prcomp(vegetation.matrix, scale = T)$x[,1] # First PCA axis
trap$vegetationPCA2 <- prcomp(vegetation.matrix, scale = T)$x[,2] # 2nd PC



#############################################################################################################################
#### MERGE ALL DATA TOGETHER

#need to eventually merge SNP data (genotypes) and habitat data (trap) but they don't have any columns in common
#metadata has fishID, which is in the same order as genotype data, so first merge those
#then using fishID in metadata, which we rename fish.id to match RADseq, we merge those data frames
#RADseq has fish.id.univ which allows us to match habitat data in trap to SNPs

metadata$fish.id<-metadata$fishID #add column with new name to match RADseq
dim(metadata) #both have same number of fish (rows) 840:6 fish were genotyped
dim(genotypes) #840:78225 (so 78224 are SNPS, first was ID #)
merged<-as.data.frame(cbind(metadata, genotypes)) #this merges the two data frames together, because they are same number of rows
dim(merged) #this has both data frames together now, with fish ID included, 840: 78231
with.fish.id<-merge(merged, RADseq, by = "fish.id", all.x = T) #this now puts SNP data with RADseq so you have fish.ID.univ
trapreduced<-trap[,c(1, 12, 13, 25, 27, 86, 88)] #only selecting columns I need from trap data
trapSNP<-merge(with.fish.id, trapreduced, by = "fishID.univ", all.x = T) #merge together SNP data and habitat data
dim(trapSNP) #840:78251

columns.to.analyze <- c(9:78232) #all the columns with SNP data in trapSNP


####################################################################################################################################
#if you want to run lakes and streams together, you need to include a habitat term

full <- trapSNP[(c(-97:-144, -673:-696)),] #eliminate marine rows from habitat.x, kept all columns
SNP_full_ID <-data.frame(full[,c(1,5,6, 9:78232, 78248:78251)]) 
dim(SNP_full_ID) 
SNP_full_ID[4:78231] <- lapply(SNP_full_ID[4:78231], as.numeric) 
#ABOVE:everything is an integer, convert to numeric 
#but leave out fish ID, watershed, habitat
focalSNP<-colnames(SNP_full_ID[,(4:78227)]) #these are all the columns with SNPs

####################################################################################################################################
#DEPTH in lakes and streams

results.storage.stream.AND.lake.depth <- as.data.frame(matrix(nrow = length(focalSNP), ncol = 11)) 

for (i in 1:length(focalSNP)){
  print(focalSNP[i])
  SNP <- SNP_full_ID[,i+3] #SNP is the column name, skip first 3 columns
  SNPfreq <- tapply(SNP, SNP_full_ID$watershed.x, mean, na.rm = T)/2
  SNPfreq[is.nan(SNPfreq)] <- NA
  Npops_polymorph <- sum(SNPfreq > 0, na.rm  = T)
  Mean_SNP_freq <- mean(SNPfreq, na.rm = T)
  if(Npops_polymorph >= 5 & Mean_SNP_freq > 0.01 & Mean_SNP_freq < 0.99){
    model<-lm(log(depth.conv.cm) ~  watershed.x*habitat.x + SNP, data=SNP_full_ID, na.action = na.omit) #SNP is second now
    sum_model<-Anova(model, type = "II")
    results.storage.stream.AND.lake.depth[i, 2:5]<-sum_model$`F value`[1:4]
    results.storage.stream.AND.lake.depth[i, 6:9]<-sum_model$`Pr(>F)`[1:4]
    results.storage.stream.AND.lake.depth[i,1]<-focalSNP[i]
  results.storage.stream.AND.lake.depth[i,10] <-Mean_SNP_freq
  results.storage.stream.AND.lake.depth[i,11] <-Npops_polymorph
  }
}


colnames(results.storage.stream.AND.lake.depth)<-c("SNP_pos", "watershed_F", "habitat_F", "SNP_F",
                                          "watershedhabitat_F", "watershed_P",
                                          "habitat_P", "SNP_P","watershedhabitat_P", "Mean_SNP_freq", "Npops_polymorph") #reordered because now doing Anova not anova

dim(results.storage.stream.AND.lake.depth) #78224:11
#lots of NAs from the SNPs that did not fit the frequency criteria- elimate those
results.storage.stream.AND.lake.depth.clean<-na.omit(results.storage.stream.AND.lake.depth, cols="SNP_F") #remove rows where SNP_F is NA
dim(results.storage.stream.AND.lake.depth.clean) #29240:11, these are the number of SNPs that fit criteria above
#there are 1164 SNPs with pvalue < 0.05, this is 4% of the SNPs that fit the criteria

#########################################
sum(depthresults$SNP_P<0.05) #1164 SNPs are considered significant
SNP_pvalue<-depthresults$SNP_P #the p values for SNP and depth
SNPq<-qvalue(p=SNP_pvalue, fdr.level = 0.05) #this sets the false discovery rate to 5%, so should get ~1100 SNPs remaining significant after correction (this is still a lot!)
sum(SNPq$significant=="TRUE") #asks how many are still significant after qvalue, says zero

#correct the pvalues 
SNPqvals<-qvalue(results.storage.stream.AND.lake.depth.clean$SNP_P, fdr.level = 0.05) #5% false discovery rate
results.storage.stream.AND.lake.depth.clean$SNPsignificance<-SNPqvals$significant #true is significant


####################################################################################################################################
#VEGETATION in lakes and streams
#cube root transformation improves it a bit (makes more normal), so I did that before doing lm

results.storage.stream.AND.lake.vegetation <- as.data.frame(matrix(nrow = length(focalSNP), ncol = 11)) #33 by 4, so will need 132 columns. 

for (i in 1:length(focalSNP)){
  print(focalSNP[i])
  SNP <- SNP_full_ID[,i+3] #SNP is the column name, skip first 3 columns
  SNPfreq <- tapply(SNP, SNP_full_ID$watershed.x, mean, na.rm = T)/2
  SNPfreq[is.nan(SNPfreq)] <- NA
  Npops_polymorph <- sum(SNPfreq > 0, na.rm  = T)
  Mean_SNP_freq <- mean(SNPfreq, na.rm = T)
  if(Npops_polymorph >= 5 & Mean_SNP_freq > 0.01 & Mean_SNP_freq < 0.99){
    model<-lm((vegetationPCA1)^(1/3) ~  watershed.x*habitat.x + SNP, data=SNP_full_ID, na.action = na.omit) #SNP is second now
    sum_model<-Anova(model, type = "II")
    results.storage.stream.AND.lake.vegetation[i, 2:5]<-sum_model$`F value`[1:4]
    results.storage.stream.AND.lake.vegetation[i, 6:9]<-sum_model$`Pr(>F)`[1:4]
    results.storage.stream.AND.lake.vegetation[i,1] <-focalSNP[i]
    results.storage.stream.AND.lake.vegetation[i,10] <-Mean_SNP_freq
    results.storage.stream.AND.lake.vegetation[i,11] <-Npops_polymorph
  }
}


colnames(results.storage.stream.AND.lake.vegetation)<-c("SNP_pos", "watershed_F", "habitat_F", "SNP_F",
                                                        "watershedhabitat_F", "watershed_P",
                                                        "habitat_P", "SNP_P","watershedhabitat_P", "Mean_SNP_freq", "Npops_polymorph") #reordered because now doing Anova not anova

dim(results.storage.stream.AND.lake.vegetation) #78224:11
results.storage.stream.AND.lake.vegetation.clean<-na.omit(results.storage.stream.AND.lake.vegetation, cols="SNP_F")
dim(results.storage.stream.AND.lake.vegetation.clean) #28249: 11

#correct the pvalues 
SNPqvals<-qvalue(results.storage.stream.AND.lake.vegetation.clean$SNP_P, fdr.level = 0.05) #5% false discovery rate
results.storage.stream.AND.lake.vegetation.clean$SNPsignificance<-SNPqvals$significant #true is significant


############################################################################################################################
#substrate is not normal

#SUBSTRATE in lakes and streams

results.storage.stream.AND.lake.substrate <- as.data.frame(matrix(nrow = length(focalSNP), ncol = 11)) 

for (i in 1:length(focalSNP)){
  print(focalSNP[i])
  SNP <- SNP_full_ID[,i+3] #SNP is the column name, skip first 3 columns
  SNPfreq <- tapply(SNP, SNP_full_ID$watershed.x, mean, na.rm = T)/2
  SNPfreq[is.nan(SNPfreq)] <- NA
  Npops_polymorph <- sum(SNPfreq > 0, na.rm  = T)
  Mean_SNP_freq <- mean(SNPfreq, na.rm = T)
  if(Npops_polymorph >= 5 & Mean_SNP_freq > 0.01 & Mean_SNP_freq < 0.99){
    model<-lm((substratePCA1^(1/3)) ~  watershed.x*habitat.x + SNP, data=SNP_full_ID, na.action = na.omit) #SNP is second now
    sum_model<-Anova(model, type = "II")
    results.storage.stream.AND.lake.substrate[i, 2:5]<-sum_model$`F value`[1:4]
    results.storage.stream.AND.lake.substrate[i, 6:9]<-sum_model$`Pr(>F)`[1:4]
    results.storage.stream.AND.lake.substrate[i,1] <-focalSNP[i]
    results.storage.stream.AND.lake.substrate[i,10] <-Mean_SNP_freq
    results.storage.stream.AND.lake.substrate[i,11] <-Npops_polymorph
  }
}


colnames(results.storage.stream.AND.lake.substrate)<-c("SNP_pos", "watershed_F", "habitat_F", "SNP_F",
                                                       "watershedhabitat_F", "watershed_P",
                                                       "habitat_P", "SNP_P","watershedhabitat_P", "Mean_SNP_freq", "Npops_polymorph") #reordered because now doing Anova not anova

dim(results.storage.stream.AND.lake.substrate) #78224:11
results.storage.stream.AND.lake.substrate.clean<-na.omit(results.storage.stream.AND.lake.substrate, cols="SNP_F")
dim(results.storage.stream.AND.lake.substrate.clean) #29161:11

#correct the pvalues
SNPqvals<-qvalue(results.storage.stream.AND.lake.substrate.clean$SNP_P, fdr.level = 0.05) #5% false discovery rate
results.storage.stream.AND.lake.substrate.clean$SNPsignificance<-SNPqvals$significant #true is significant



#now we want to separate the STREAM data out to run a model only on SNP data from the stream fish because flowrate was only in streams
############################################################################################################################
#STREAM flow

stream <- trapSNP[trapSNP$habitat.x == "stream",]
SNP_stream_ID <-data.frame(stream[,c(1,5,6, 9:78232, 78248:78251)]) #all the SNPs taken from stream fish and Fish ID (1) and habitat characteristics
dim(SNP_stream_ID) 
SNP_stream_ID[4:78231] <- lapply(SNP_stream_ID[4:78231], as.numeric) 
focalSNP<-colnames(SNP_stream_ID[,(4:78227)])
SNP_stream_ID$newflowrate<-SNP_stream_ID$flowrate
SNP_stream_ID$newflowrate[SNP_stream_ID$newflowrate>0]<-1 #replace all nonzero flowrates with 1 so it becomes binary (flow, no flow)

results.storage.stream.flowrate <- as.data.frame(matrix(nrow = length(focalSNP), ncol = 4))

#May have to run in sections due to NAs

for (i in 1:length(focalSNP)){
  print(focalSNP[i])
  SNP <- SNP_stream_ID[,i+3] #SNP is the column name, skip first 3 columns
  SNPfreq <- tapply(SNP, SNP_stream_ID$watershed.x, mean, na.rm = T)/2
  SNPfreq[is.nan(SNPfreq)] <- NA
  Npops_polymorph <- sum(SNPfreq > 0, na.rm  = T)
  Mean_SNP_freq <- mean(SNPfreq, na.rm = T)
  if(Npops_polymorph >= 5 & Mean_SNP_freq > 0.01 & Mean_SNP_freq < 0.99) {
    model<-glm(newflowrate ~ watershed.x + SNP, data=SNP_stream_ID, family=binomial, na.action = na.omit) #newflowrate is 0 or 1, with NAs
    sum_model<-summary(model)
    results.storage.stream.flowrate[i, 2]<-coef(sum_model)['SNP', 'Pr(>|z|)'] 
    results.storage.stream.flowrate[i,1] <-focalSNP[i]
    results.storage.stream.flowrate[i,3] <-Mean_SNP_freq
    results.storage.stream.flowrate[i,4] <-Npops_polymorph
}
}

colnames(results.storage.stream.flowrate)<-c("SNP_pos", "SNP_P", "Mean_SNP_freq", "Npops_polymorph") #not saving watershed because they are all separated
                                                        
dim(results.storage.stream.flowrate) #78224:4
results.storage.stream.flowrate.clean<-na.omit(results.storage.stream.flowrate, cols="SNP_F")
dim(results.storage.stream.flowrate.clean) #22714:11

#correct the pvalues and save significance
SNPqvals<-qvalue(results.storage.stream.flowrate.clean$SNP_P, fdr.level = 0.05) #5% false discovery rate
results.storage.stream.flowrate.clean$SNPsignificance<-SNPqvals$significant #true is significant

####################################################################################################################################

