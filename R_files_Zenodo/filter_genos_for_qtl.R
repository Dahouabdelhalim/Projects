setwd("/Users/cpeichel/Desktop")

library(data.table) #load the big data file reading library
SNP.sorted.merged <-  fread("output_Roberts_genotypes.tab") #read in the file 
#output_Pye_genotypes.tab
#output_Misty_fam1.tab
#output_Misty_fam2.tab
#output_Boot_genotypes.tab; used for Boot fam 2
#output_Boot_combined.tab; used for Boot fam 1 (added resequencing of two F1s and 3 grandparents)

SNP.sorted.merged[1:10,1:10]
temp <- apply(SNP.sorted.merged, 1, function(x) length(unique(x))) #count number of unique character in a row, this will be used to filter
unique(temp) #number of unique characters among loci 
# 3 4 5 6
SNP.sorted.merged2 <- cbind(SNP.sorted.merged,temp) #bind the counts used for filtering with the data 
dim(SNP.sorted.merged2)
#Roberts: 929819 237
#Pye: 1049411 383
#Misty: 1151554     567
#Misty_fam1: 860000	287
#Misty_fam2: 781225    297
#Boot: 1404482     593
#Boot_combined: 1726951     595

SNP.sorted.merged3 <- subset(SNP.sorted.merged2, SNP.sorted.merged2$temp > 4) #makes is so that all rows have at least 2 of the three possible genotypes. i.e. we want variation
SNP.sorted.merged3$temp <- NULL #delete the count column
dim(SNP.sorted.merged3) #check the dimentions of the dataframe to determine number of SNPs still in action
#Roberts: 85861 236
#Pye: 93961   382
#Misty: 92497   566
#Misty: 64307   272
#Misty_fam 1: 76352  286
#Misty_fam2: 54463   296
#Boot: 99588   592
#Boot_combined: 101741    594

SNP.sorted.merged3[SNP.sorted.merged3 == "NN"] <- NA #replace all NN values with NAs for filtering 
SNP.sorted.merged3$combined <- paste(SNP.sorted.merged3$CHROM, SNP.sorted.merged3$POS,sep="_") #make combined locus column 

head(SNP.sorted.merged3) #to figure out which columns the F1s are in
dim(SNP.sorted.merged3) 
#Roberts: 85861 237
#Pye: 93961   383
#Misty: 92497   567
#Misty_fam1: 76352  287
#Misty_fam2: 54463   297
#Boot: 99588   593
#Boot_combined: 101741    595
  
#look at sites that are hets in all F1 parents 
#all combinations of heterozygotes  
#  AC   AG   AT   CA     CG   CT   GA   GC   GT   TA   TC   TG  
#F1 parents
#Rob_LxS_f3xm3_F1_f1 (225)
#Rob_LxS_f3xm3_F1_f4 (226)
#Rob_LxS_f3xm3_F1_m1 (227)
#Rob_LxS_f3xm3_F1_m4 (228)
#Rob_LxS_f5xm5_F1_f1 (229)
#Rob_LxS_f5xm5_F1_m1 (230)
#Rob_LxS_f5xm5_F1_m3 (231)
#Rob_LxS_f5xm5_F1_m4 (232)
# 8 (ok)

#Pye_LxS_f1xm1_F1_f2
#Pye_LxS_f1xm1_F1_m2
#Pye_LxS_f1xm1_F1_m3
#Pye_LxS_f2xm2_F1_f1
#Pye_LxS_f2xm2_F1_m1
#Pye_LxS_f4xm4_F1_f1
#Pye_LxS_f4xm4_F1_f2
#Pye_LxS_f4xm4_F1_m1
#Pye_LxS_f4xm4_F1_m3
#Pye_LxS_f4xm4_F1_m4
# 10 (ok)

#Misty_fam1 file
#Misty_LxS_f1xm1_F1_f2 (273); EXCLUDE!!!
#Misty_LxS_f1xm1_F1_f3 (274)
#Misty_LxS_f1xm1_F1_f4 (275)
#Misty_LxS_f1xm1_F1_m2 (276)
#Misty_LxS_f1xm1_F1_m3 (277)
#Misty_LxS_f1xm1_F1_m4 (278)
# 6 (Misty_LxS_f1xm1_F1_f1 was not sequenced)

#Misty_fam2 file
#Misty_LxS_f2xm2_F1_f1 (289)
#Misty_LxS_f2xm2_F1_f4 (290)
#Misty_LxS_f2xm2_F1_f7 (291)
#Misty_LxS_f2xm2_F1_m11 (292)
#Misty_LxS_f2xm2_F1_m12 (293)
#Misty_LxS_f2xm2_F1_m4 (294)
# 6 (ok)

#For Boot fam 1, used Boot combined genotypes file
#Boot_LXS_f1xm1_F1_f2 (3); this is the new data but was not combined because of typo in name
#Boot_LxS_f1xm1_F1_f2 (580); this is the original data for this F1
#Boot_LxS_f1xm1_F1_f4 (581)
#Boot_LxS_f1xm1_F1_f5 (582)
#Boot_LxS_f1xm1_F1_m2 (583); this male was also resequenced but data was combined
#Boot_LxS_f1xm1_F1_m4 (584)
#Boot_LxS_f1xm1_F1_m5 (585)
#Boot_LxS_f2xm2_F1_f2 (586)
#Boot_LxS_f2xm2_F1_f3 (587)
#Boot_LxS_f2xm2_F1_f5 (588)
#Boot_LxS_f2xm2_F1_m2 (589)
#Boot_LxS_f2xm2_F1_m3 (590)
#Boot_LxS_f2xm2_F1_m5 (591)
#12 (ok)

#For Boot fam 2, used original Boot genotypes file
#Boot_LxS_f1xm1_F1_f2 (579)
#Boot_LxS_f1xm1_F1_f4 (580)
#Boot_LxS_f1xm1_F1_f5 (581)
#Boot_LxS_f1xm1_F1_m2 (582)
#Boot_LxS_f1xm1_F1_m4 (583)
#Boot_LxS_f1xm1_F1_m5 (584)
#Boot_LxS_f2xm2_F1_f2 (585)
#Boot_LxS_f2xm2_F1_f3 (586)
#Boot_LxS_f2xm2_F1_f5 (587)
#Boot_LxS_f2xm2_F1_m2 (588)
#Boot_LxS_f2xm2_F1_m3 (589)
#Boot_LxS_f2xm2_F1_m5 (590)


f1 <- SNP.sorted.merged3[,c(1:2,156:159)] #subset out F1 parents and meta data the column numbers for the F1s will vary depending on the watershed 
#Roberts_all: f1<-SNP.sorted.merged3[,c(1:2,225:232)] 
#Roberts_fam3: f1<-SNP.sorted.merged3[,c(1:2,225:228)]
#Roberts_fam5: f1<-SNP.sorted.merged3[,c(1:2,229:232)]
#Pye_all: f1<-SNP.sorted.merged3[,c(1:2,367:376)]
#Pye_fam1: f1<-SNP.sorted.merged3[,c(1:2,367:369)] 
#Pye_fam2: f1<-SNP.sorted.merged3[,c(1:2,370:371)]
#Pye_fam4: f1<-SNP.sorted.merged3[,c(1:2,372:376)]
#Misty_all: f1 <- SNP.sorted.merged3[,c(1:2,553:564)]
#Misty_fam1: f1 <- SNP.sorted.merged3[,c(1:2,274:278)]
#Misty_fam2: f1 <- SNP.sorted.merged3[,c(1:2,289:294)]
#Boot_all: f1<-SNP.sorted.merged3[,c(1:2,579:590)]
#Boot_fam1: f1<-SNP.sorted.merged3[,c(1:2,579:584)]
#Boot_fam2: f1<-SNP.sorted.merged3[,c(1:2,585:590)]
#Boot_fam1_combined: f1<-SNP.sorted.merged3[,c(1:3,581:585)]
#Boot_fam1_combined: f1<-SNP.sorted.merged3[,c(1:2,581,584)]; only used these two F1s for filtering because too many NAs in the other 4 F1s (see below)

head(f1) #to check that got the right columns
dim(f1)
#Roberts_all: 85861 10
#Roberts_fam3: 85861 6
#Roberts_fam5: 85861 6
#Pye_all: 93961    12
#Pye_fam1: 93961 5
#Pye_fam2: 93961     4
#Pye_fam4: 93961     7
#Misty_all: 92497    14 
#Misty_fam1: 76352     7
#Misty_fam2: 54463     8
#Boot_all: 99588    14
#Boot_fam1: 99588    8
#Boot_fam2: 99588    8 
#Boot_fam1_combined: 101741      8 (all F1s)
#Boot_fam1_combined: 101741      4 (2 good F1s)

f1_2 <- f1[-which(rowMeans(is.na(f1)) > 0.05, )] #this drops all the sites with more than 5% missing data 
#Roberts: this should be all sites with at least 1 NA in an F1
#Pye: this should be all sites with at least 1 NA in an F1
#Misty: if I do this for Misty_all, I lose almost all sites; a lot of missing data in Misty F1s?
#Misty_fam1: if I use >0.05 I get 378 markers, if I use >0.2 (1 fish ok with NA), I get 3592 markers
#Misty_fam2: this should be all sites with at least 1 NA in an F1
#colMeans(is.na(f1_2))
#             CHROM                   POS Misty_LxS_f1xm1_F1_f2 Misty_LxS_f1xm1_F1_f3 Misty_LxS_f1xm1_F1_f4 
#          0.000000000           0.000000000           0.795935412           0.003340757           0.052616927 
#Misty_LxS_f1xm1_F1_m2 Misty_LxS_f1xm1_F1_m3 Misty_LxS_f1xm1_F1_m4 
#          0.034521158           0.004732739           0.003619154 
#individual Misty_LxS_f1xm1_F1_f2 is mostly NAs; remake the file for filtering without this individual

#Boot_fam1: if I do this, I lose almost all sites, a lot of missing data in Boot F1s?; if I use >0.25 (2 fish with NA), I get 370 markers, add in resequencing data
#Boot_fam1_combined
#look at all f1s
#f1<-SNP.sorted.merged3[,c(1:3,580:591)]
#colMeans(is.na(f1))
# 			CHROM                  POS 				Boot_LXS_f1xm1_F1_f2 Boot_LxS_f1xm1_F1_f2 Boot_LxS_f1xm1_F1_f4 Boot_LxS_f1xm1_F1_f5 
#           0.0000000            0.0000000            0.9941027            0.9999509            0.7347382            0.9833400 
#Boot_LxS_f1xm1_F1_m2 Boot_LxS_f1xm1_F1_m4 Boot_LxS_f1xm1_F1_m5 Boot_LxS_f2xm2_F1_f2 Boot_LxS_f2xm2_F1_f3 Boot_LxS_f2xm2_F1_f5 
#           0.9753688            0.8561838            0.9651861            0.7987832            0.9019766            0.8147355 
#Boot_LxS_f2xm2_F1_m2 Boot_LxS_f2xm2_F1_m3 Boot_LxS_f2xm2_F1_m5 
#           0.8358675            0.6472907            0.8170256 
#data for Boot_LXS_f1xm1_F1_f2 and Boot_LxS_f1xm1_F1_m2 still bad
#f1_2 <- f1[-which(rowMeans(is.na(f1)) > 0.05, )]
#Boot_fam1_combined 53  8 (>0.05)
#Boot_fam1_combined 665  8 (>0.25)
#this is not enough, just use the two good F1s and sort later

dim(f1_2) #look at number of sites
#Roberts_all: 3634 10
#Roberts_fam3: 6820 6 
#Roberts_fam5: 12068     6
#Pye_all: 465 12
#Pye_fam1: 14029     5
#Pye_fam2: 4018 4
#Pye_fam4: 12663 7
#Misty_all: 1 14
#Misty_fam1: 378   8 (>0.05)
#Misty_fam1: 3592  8 (>0.20)
#Misty_fam1: 3237	7 (when I exclude problematic F1)
#Misty_fam2: 1954    8
#Boot_fam1: 0  8
#Boot_fam2: 4616    8 
#Boot_fam1_combined: 10977     4 


head(f1_2)
#this is the numbering for filtering 
f1_2[f1_2 == "AC"] <- 1 
f1_2[f1_2 == "AG"] <- 1
f1_2[f1_2 == "AT"] <- 1 
f1_2[f1_2 == "CA"] <- 1 
f1_2[f1_2 == "CG"] <- 1  
f1_2[f1_2 == "CT"] <- 1 
f1_2[f1_2 == "GA"] <- 1 
f1_2[f1_2 == "GC"] <- 1
f1_2[f1_2 == "GT"] <- 1 
f1_2[f1_2 == "TA"] <- 1 
f1_2[f1_2 == "TC"] <- 1  
f1_2[f1_2 == "TG"] <- 1 
f1_2[f1_2 == "TT"] <- 0
f1_2[f1_2 == "AA"] <- 0
f1_2[f1_2 == "GG"] <- 0
f1_2[f1_2 == "CC"] <- 0


str(f1_2[1:10,])#depending upon number of F1s
#need to makesure all the relevant columns are formatted to be numeric - will have to change the column names to match all the relevant F1s:
f1_2$Rob_LxS_f3xm3_F1_f1<- as.numeric(f1_2$Rob_LxS_f3xm3_F1_f1)
f1_2$Rob_LxS_f3xm3_F1_f4<- as.numeric(f1_2$Rob_LxS_f3xm3_F1_f4)
f1_2$Rob_LxS_f3xm3_F1_m1<- as.numeric(f1_2$Rob_LxS_f3xm3_F1_m1)
f1_2$Rob_LxS_f3xm3_F1_m4<- as.numeric(f1_2$Rob_LxS_f3xm3_F1_m4)

f1_2$Rob_LxS_f5xm5_F1_f1<- as.numeric(f1_2$Rob_LxS_f5xm5_F1_f1)
f1_2$Rob_LxS_f5xm5_F1_m1<- as.numeric(f1_2$Rob_LxS_f5xm5_F1_m1)
f1_2$Rob_LxS_f5xm5_F1_m3<- as.numeric(f1_2$Rob_LxS_f5xm5_F1_m3)
f1_2$Rob_LxS_f5xm5_F1_m4<- as.numeric(f1_2$Rob_LxS_f5xm5_F1_m4)

f1_2$Pye_LxS_f1xm1_F1_f2<- as.numeric(f1_2$Pye_LxS_f1xm1_F1_f2)
f1_2$Pye_LxS_f1xm1_F1_m2<- as.numeric(f1_2$Pye_LxS_f1xm1_F1_m2)
f1_2$Pye_LxS_f1xm1_F1_m3<- as.numeric(f1_2$Pye_LxS_f1xm1_F1_m3)

f1_2$Pye_LxS_f2xm2_F1_f1<- as.numeric(f1_2$Pye_LxS_f2xm2_F1_f1)
f1_2$Pye_LxS_f2xm2_F1_m1<- as.numeric(f1_2$Pye_LxS_f2xm2_F1_m1)

f1_2$Pye_LxS_f4xm4_F1_f1<- as.numeric(f1_2$Pye_LxS_f4xm4_F1_f1)
f1_2$Pye_LxS_f4xm4_F1_f2<- as.numeric(f1_2$Pye_LxS_f4xm4_F1_f2)
f1_2$Pye_LxS_f4xm4_F1_m1<- as.numeric(f1_2$Pye_LxS_f4xm4_F1_m1)
f1_2$Pye_LxS_f4xm4_F1_m3<- as.numeric(f1_2$Pye_LxS_f4xm4_F1_m3)
f1_2$Pye_LxS_f4xm4_F1_m4<- as.numeric(f1_2$Pye_LxS_f4xm4_F1_m4)

#f1_2$Misty_LxS_f1xm1_F1_f2 <- as.numeric(f1_2$Misty_LxS_f1xm1_F1_f2)
f1_2$Misty_LxS_f1xm1_F1_f3 <- as.numeric(f1_2$Misty_LxS_f1xm1_F1_f3)
f1_2$Misty_LxS_f1xm1_F1_f4 <- as.numeric(f1_2$Misty_LxS_f1xm1_F1_f4)
f1_2$Misty_LxS_f1xm1_F1_m2 <- as.numeric(f1_2$Misty_LxS_f1xm1_F1_m2)
f1_2$Misty_LxS_f1xm1_F1_m3 <- as.numeric(f1_2$Misty_LxS_f1xm1_F1_m3)
f1_2$Misty_LxS_f1xm1_F1_m4 <- as.numeric(f1_2$Misty_LxS_f1xm1_F1_m4)

f1_2$Misty_LxS_f2xm2_F1_f1 <- as.numeric(f1_2$Misty_LxS_f2xm2_F1_f1)
f1_2$Misty_LxS_f2xm2_F1_f4 <- as.numeric(f1_2$Misty_LxS_f2xm2_F1_f4)
f1_2$Misty_LxS_f2xm2_F1_f7 <- as.numeric(f1_2$Misty_LxS_f2xm2_F1_f7)
f1_2$Misty_LxS_f2xm2_F1_m11 <- as.numeric(f1_2$Misty_LxS_f2xm2_F1_m11)
f1_2$Misty_LxS_f2xm2_F1_m12 <- as.numeric(f1_2$Misty_LxS_f2xm2_F1_m12)
f1_2$Misty_LxS_f2xm2_F1_m4 <- as.numeric(f1_2$Misty_LxS_f2xm2_F1_m4)

#f1_2$Boot_LxS_f1xm1_F1_f2 <- as.numeric(f1_2$Boot_LxS_f1xm1_F1_f2)
f1_2$Boot_LxS_f1xm1_F1_f4 <- as.numeric(f1_2$Boot_LxS_f1xm1_F1_f4)
#f1_2$Boot_LxS_f1xm1_F1_f5 <- as.numeric(f1_2$Boot_LxS_f1xm1_F1_f5)
#f1_2$Boot_LxS_f1xm1_F1_m2 <- as.numeric(f1_2$Boot_LxS_f1xm1_F1_m2)
f1_2$Boot_LxS_f1xm1_F1_m4 <- as.numeric(f1_2$Boot_LxS_f1xm1_F1_m4)
#f1_2$Boot_LxS_f1xm1_F1_m5 <- as.numeric(f1_2$Boot_LxS_f1xm1_F1_m5)

f1_2$Boot_LxS_f2xm2_F1_f2 <- as.numeric(f1_2$Boot_LxS_f2xm2_F1_f2)
f1_2$Boot_LxS_f2xm2_F1_f3 <- as.numeric(f1_2$Boot_LxS_f2xm2_F1_f3)
f1_2$Boot_LxS_f2xm2_F1_f5 <- as.numeric(f1_2$Boot_LxS_f2xm2_F1_f5)
f1_2$Boot_LxS_f2xm2_F1_m2 <- as.numeric(f1_2$Boot_LxS_f2xm2_F1_m2)
f1_2$Boot_LxS_f2xm2_F1_m3 <- as.numeric(f1_2$Boot_LxS_f2xm2_F1_m3)
f1_2$Boot_LxS_f2xm2_F1_m5 <- as.numeric(f1_2$Boot_LxS_f2xm2_F1_m5)

temp <- rowSums(f1_2[,c(3:10)]) #count numbers of hets -change to match dimentions of dataframe 
f1_3 <- cbind(f1_2,temp)
#only interested in loci with a temp which matches the number of F1s you have - in this example there were 6 reference F1s 
#Rob_all should be 8 F1s  
#Pye_all should be 10 F1s
#Misty_all should be 12 F1s
#Misty_fam1 shoudl be 5 F1s
#Misty_fam2 should be 6 F1s
#Boot_fam2 should be 6 F1s
#Boot_fam1 should be 2 F1s

f1_loci <- subset(f1_3,f1_3$temp == 8) #change to match number of F1s
dim(f1_loci)
#Roberts_all: 8 F1s are het in 220 loci
#Roberts_fam3: 4 F1s are het in 1395 loci
#Roberts_fam5: 4 F1s are het in 2427 loci
#Pye_all: 10 F1s are het in 65 loci
#Pye_fam1: 3 F1s are het in 3688 loci
#Pye_fam2: 2 F1s are het in 1731 loci
#Pye_fam4: 5 F1s are het in 2469 loci
#Misty_all: 12 F1s are het in 0 loci
#Misty_fam1: 5 F1s are het in 752 loci 
#Misty_fam2: 6 F1s are het in 683 loci
#Boot_fam2: 6 F1s are het in 994 loci
#Boot_fam1: 2 F2s are het in 3998 loci

head(f1_loci)
f1_loci$combined <- paste(f1_loci$CHROM,f1_loci$POS,sep="_") #make combined locus column
keeps <- f1_loci$combined #these are the loci we want to look at in the F2s
#now we want to filter these out of the bigger dataset 
SNP_for_qtl <- SNP.sorted.merged3[(SNP.sorted.merged3$combined%in% keeps),]
dim(SNP_for_qtl)
#Misty_fam1: 752 287
#Misty_fam2: 683 297

SNP_for_qtl_fam <- SNP_for_qtl[,c(1:75,225:237)]#this will only include F2s from one family, need to change for each file depending on number of F2s in each family
#Roberts_fam3: SNP_for_qtl_fam <- SNP_for_qtl[,c(1:75,225:237)]
#Roberts_fam5: SNP_for_qtl_fam <- SNP_for_qtl[,c(1:2,76:237)]
#Pye_fam1: SNP_for_qtl_fam <- SNP_for_qtl[,c(1:75,367:383)]
#Pye_fam2: SNP_for_qtl_fam <- SNP_for_qtl[,c(1:2,76:168,367:383)]
#Pye_fam4: SNP_for_qtl_fam <- SNP_for_qtl[,c(1:2,169:383)]
#Boot_fam2: SNP_for_qtl_fam <- SNP_for_qtl[,c(1:2,283:593)]
#Boot_fam1: SNP_for_qtl_fam <- SNP_for_qtl[,c(1:2,4:283,578:595)]

dim(SNP_for_qtl_fam)
#Roberts_fam3: 1395 88
#Roberts_fam5: 2427 164
#Pye_fam1: 3688   92
#Pye_fam2: 1731  112
#Pye_fam4: 2469  217
#Boot_fam2: 994  313
#Boot_fam1: 3998  300

head(SNP_for_qtl_fam) 
#to check that I removed the correct F2s

#count individuals with more than 50% missing data
#z<-colMeans(is.na(SNP_for_qtl))
#z<-colMeans(is.na(SNP_for_qtl_fam))
#write.table(z,"Xwatershed_famX_NAs.txt", sep="\\t", col.names=TRUE, quote=FALSE,) #to record percent missing data for each individual at these SNPs
#length(which(colMeans(is.na(SNP_for_qtl)) > 0.50 ))
#length(which(colMeans(is.na(SNP_for_qtl_fam)) > 0.50 ))
#which(colMeans(is.na(SNP_for_qtl)) > 0.5)
#which(colMeans(is.na(SNP_for_qtl_fam)) > 0.5)
#Roberts_fam3: 7 (this includes 4 GPs) = 3 bad F2s/73 F2s = 70 good F2s
#Roberts_fam5: 87 (this includes 1 F1, 4 GPs) = 82 bad F2s/149 F2s = 67 good F2s
#Pye_fam1: 42 (this includes 6 F1s from other fams, 5 GPs = 11) = 31 bad F2s/73 F2s = 42 good F2s
#Pye_fam2: 10 (this includes 2 F1s from other fams, 3 GPs = 5) = 5 bad F2s/93 F2s = 88 good F2s
#Pye_fam4: 56 (this includes 4 F1s from other fams, 4 GPs = 8) = 48 bad F2s/198 F2s = 150 good F2s
#Misty_fam1: 143 (this includes 6 F1s from other fam, 1 bad F1 this fam, 4 GPs = 11) = 132 bad F2s/268 F2s = 136 good F2s
#Misty_fam2: 85 (this includes 6 F1s from other fam, 2 GPs = 8) = 77 bad F2s/280 F2s = 203 good F2s
#Boot_fam2: 63 (this includes 4 F1s from other fam, 4 GPs = 8) = 55 bad F2s/294 F2s = 239 good F2s
#Boot_fam1: 99 (this includes 4 F1s this fam, 1 F1 other fam, 5 GPs = 10) = 89 bad F2s/280 F2s = 191 good F2s

#length(which(colMeans(is.na(SNP_for_qtl)) < 0.50 ))
#length(which(colMeans(is.na(SNP_for_qtl_fam)) < 0.50 ))
#Roberts_fam3: 81 (this includes 8 F1s, plus 3 locus columns) = 70 good F2s
#Roberts_fam5: 77 (this includes 7 F1s, plus 3 locus columens) = 67 good F2s
#Pye_fam1: 50 (this includes 3 F1s from this fam, 1 F1 other fam, 1 GP this fam, plus 3 locus columns) = 42 good F2s
#Pye_fam2: 102 (this includes 9 F1s, 2 GPs, plus 3 locus columns) = 88 good F2s
#Pye_fam4: 161 (this incldes 6 F1s, 2 GPs, plus 3 locus columns) = 150 good F2s
#Misty_fam1: 144, this includes 5 F1s this fam, plus 3 locus columns = 136 good F2s 
#Misty_fam2: 212, this includes 6 F1s this fam, plus 3 locus columns = 203 good F2s
#Boot_fam2: 250, this includes 6 F1s this fam, 2 F1s other fam, plus 3 locus columns = 239 good F2s
#Boot_fam1: 200, this includes 2 F1s this fam, 5 F1s other fam, plus 3 locus columns = 190 good F2s

#SNP_for_qtl_filt <- SNP_for_qtl[,-which(,colMeans(is.na(SNP_for_qtl)) > 0.50 )]
#this doesn't work for some reason - returns NULL!!!!
#Diana said to remove metadata columns first
#decided to keep these individuals for filtering later in JoinMap or R/qtl

SNP_for_qtl_filt <- SNP_for_qtl[-which(rowMeans(is.na(SNP_for_qtl)) > 0.10, )] #this gives us all the sites with less than 10% missing data
#SNP_for_qtl_filt <- SNP_for_qtl_fam[-which(rowMeans(is.na(SNP_for_qtl_fam)) > 0.10, )]
#SNP_for_qtl_filt <- SNP_for_qtl_filt_ind[-which(rowMeans(is.na(SNP_for_qtl)) > 0.10, )]

SNP_for_qtl_filt <- SNP_for_qtl[-which(rowMeans(is.na(SNP_for_qtl)) > 0.10, )]
dim(SNP_for_qtl_filt)
#Roberts_all: 20 loci have less than 10% missing data
#Roberts_all: 85 loci have less than 20% missing data
#Roberts_all: 136 loci have less than 30% missing data
#Roberts_all: 169 loci have less than 40% missing data
#Roberts_all: 203 loci have less than 50% missing data

SNP_for_qtl_filt <- SNP_for_qtl_fam[-which(rowMeans(is.na(SNP_for_qtl_fam)) > 0.10, )]
dim(SNP_for_qtl_filt)
#Roberts_fam3: 97 loci have less than 10% missing data 
#Roberts_fam3: 470 loci have less than 15% missing data 
#Roberts_fam3: 713 loci have less than 20% missing data
#Roberts_fam3: 905 loci have less than 25% missing data; use this 
#Roberts_fam3: 1035 loci have less than 30% missing data 
#Roberts_fam3: 1216 loci have less than 40% missing data 
#Roberts_fam3: 1323 loci have less than 50% missing data

SNP_for_qtl_filt <- SNP_for_qtl_fam[-which(rowMeans(is.na(SNP_for_qtl_fam)) > 0.10, )]
dim(SNP_for_qtl_filt)
#Roberts_fam5: 83 loci have less than 10% missing data
#Roberts_fam5: 233 loci have less than 15% missing data 
#Roberts_fam5: 418 loci have less than 20% missing data
#Roberts_fam5: 588 loci have less than 25% missing data; use this
#Roberts_fam5: 725 loci have less than 30% missing data
#Roberts_fam5: 984 loci have less than 40% missing data
#Roberts_fam5: 1267 loci have less than 50% missing data

SNP_for_qtl_filt <- SNP_for_qtl[-which(rowMeans(is.na(SNP_for_qtl)) > 0.10, )]
dim(SNP_for_qtl_filt)
#Pye_all: 2 loci have less than 10% missing data
#Pye_all: 26 loci have less than 20% missing data
#Pye_all: 54 loci have less than 30% missing data
#Pye_all: 62 loci have less than 40% missing data
#Pye_all: 0 loci have less than 50% missing data

SNP_for_qtl_filt <- SNP_for_qtl_fam[-which(rowMeans(is.na(SNP_for_qtl_fam)) > 0.10, )]
dim(SNP_for_qtl_filt)
#Pye_fam1: 33 loci have less than 10% missing data
#Pye_fam1: 394 loci have less than 20% missing data
#Pye_fam1: 661 loci have less than 25% missing data; use this
#Pye_fam1: 840 loci have less than 30% missing data
#Pye_fam1: 1273 loci have less than 40% missing data
#Pye_fam1: 2853 loci have less than 50% missing data

SNP_for_qtl_filt <- SNP_for_qtl_fam[-which(rowMeans(is.na(SNP_for_qtl_fam)) > 0.10, )]
dim(SNP_for_qtl_filt)
#Pye_fam2: 216 loci have less than 10% missing data
#Pye_fam2: 848 loci have less than 20% missing data
#Pye_fam2: 1077 loci have less than 25% missing data; use this
#Pye_fam2: 1197 loci have less than 30% missing data
#Pye_fam2: 1437 loci have less than 40% missing data
#Pye_fam2: 1602 loci have less than 50% missing data

SNP_for_qtl_filt <- SNP_for_qtl_fam[-which(rowMeans(is.na(SNP_for_qtl_fam)) > 0.10, )]
dim(SNP_for_qtl_filt)
#Pye_fam4: 26 loci have less than 10% missing data
#Pye_fam4: 228 loci have less than 20% missing data
#Pye_fam4: 731 loci have less than 25% missing data; use this
#Pye_fam4: 1252 loci have less than 30% missing data
#Pye_fam4: 1988 loci have less than 40% missing data
#Pye_fam4: 2363 loci have less than 50% missing data

SNP_for_qtl_filt <- SNP_for_qtl[-which(rowMeans(is.na(SNP_for_qtl)) > 0.10, )]
dim(SNP_for_qtl_filt)
#Misty_fam1: 4 loci have less than 10% missing data
#Misty_fam1: 94 loci have less than 20% missing data
#Misty_fam1: 130 loci have less than 25% missing data
#Misty_fam1: 190 loci have less than 30% missing data
#Misty_fam1: 271 loci have less than 40% missing data
#Misty_fam1: 382 loci have less than 50% missing data; use this

SNP_for_qtl_filt <- SNP_for_qtl[-which(rowMeans(is.na(SNP_for_qtl)) > 0.10, )]
dim(SNP_for_qtl_filt)
#Misty_fam2: 18 loci have less than 10% missing data
#Misty_fam2: 87 loci have less than 20% missing data
#Misty_fam2: 195 loci have less than 25% missing data
#Misty_fam2: 277 loci have less than 30% missing data
#Misty_fam2: 431 loci have less than 40% missing data
#Misty_fam2: 530 loci have less than 50% missing data; use this

SNP_for_qtl_filt <- SNP_for_qtl_fam[-which(rowMeans(is.na(SNP_for_qtl_fam)) > 0.10, )]
dim(SNP_for_qtl_filt)
#Boot_fam2: 29 loci have less than 10% missing data
#Boot_fam2: 343 loci have less than 20% missing data
#Boot_fam2: 468 loci have less than 25% missing data; use this
#Boot_fam2: 575 loci have less than 30% missing data
#Boot_fam2: 792 loci have less than 40% missing data
#Boot_fam2: 907 loci have less than 50% missing data

SNP_for_qtl_filt <- SNP_for_qtl_fam[-which(rowMeans(is.na(SNP_for_qtl_fam)) > 0.10, )]
dim(SNP_for_qtl_filt)
#Boot_fam1: 46 loci have less than 10% missing data
#Boot_fam1: 706 loci have less than 20% missing data
#Boot_fam1: 1102 loci have less than 25% missing data; use this but filter further in Excel to find markers segregating 1:2:1 in all F2 subfamilies
#Boot_fam1: 1435 loci have less than 30% missing data
#Boot_fam1: 2063 loci have less than 40% missing data
#Boot_fam1: 2578 loci have less than 50% missing data
#in excel, calculate freq of hets in each subfamily; delete markers with <0.3 hets or >0.7 hets in that family
#785 markers remain; overall freq of hets from 0.438-0.644

write.table(SNP_for_qtl_filt,"Xwatershed_famX_qtl_SNPs.txt", sep="\\t", row.names=FALSE, quote=FALSE,)
#Further filtering performed in JoinMap 4.1
