#load file of genotypes that is generated in the "PCA.R" script
#each row is an individual
setwd("~/Documents")
df <- read.csv("majorcounts_April15-2020.csv",header=T, stringsAsFactors = F)

#name sample column, then keep plants from the 24 populations
colnames(df)[1]<- "Sample"
df <- df[1:167,]
ids <- df[,1]

#Genotyper converts into BiA format
#for(y in 1:nrow(df)){
#  for(z in 1:ncol(df)){ 
#      if((is.na(df[y,z]))>0) { df[y,z] <- -9}
#      if((!is.na(df[y,z])>0)&((df[y,z]==0)>0)) { df[y,z] <- "AAAA"}
#      if((!is.na(df[y,z])>0)&((df[y,z]==1)>0)) { df[y,z] <- "AAAB"}
#      if((!is.na(df[y,z])>0)&((df[y,z]==2)>0)) { df[y,z] <- "AABB"}
#      if((!is.na(df[y,z])>0)&((df[y,z]==3)>0)) { df[y,z] <- "ABBB"}
#      if((!is.na(df[y,z])>0)&((df[y,z]==4)>0)) { df[y,z] <- "BBBB"}
#    }}

#force read proportions into genotype calls (0, 0.25, 0.50, 0.75, and 1)
df[,2:ncol(df)] <- df[,2:ncol(df)]/4
df[is.na(df)] <- -9

#Sample column is a factor
df$Sample <- as.factor(ids)

#set of 24
pops <- c("AL2012", "AL79", "ALBG", "AR125", "AR56", "IA10", "IN46", "IN77", "KS60", "KY51", "MI126", "MI127", "MN117", "MN118", "MO115", "MO116", "MO49", "MO57", "OH119", "OH64", "OK61", "PA27", "TN19", "WI128")
coder <- ids
for(i in 1:length(ids)){
  name <- pops[pmatch(pops,ids[i], nomatch=NA)>0]
  coder[i] <-name[!is.na(name)]
}

coder <- coder[1:161]
df<- df[1:161,]

#insert Pop column
df$Pop <- coder
df$Pop <- as.factor(df$Pop)

#ploidy numeric vector
df$Ploidy <- 4

#Genotype format is biallelic (BiA) or allele frequency (freq)
df$Format <- "freq"
df$Format <- as.factor(df$Format)

#move columns so ordered as in stampp: Sample,Pop,Ploidy,Format,nLoci
df2 <- df[,c(1,5531:5533,2:5530)]
df2 <- df2[,c(1:4,5:5533)]

write.csv(df2, "stampp_input.csv")
#upload this file now as df2, with stringsAsFactors=T, then peel off the first column
df2 <- read.csv("stampp_input.csv", stringsAsFactors = T)
df2 <- df2[,-1]
#subset df2 based on pops with missing data

#eliminate sites with >95% missing data
df2[df2==-9.0] <- NA
df2 <- df2[,!sapply(df2, function(x) mean(is.na(x)))>0.95]

#use stAMMP package to calculate pairwise Fst using Weir and Cockerham (1984) method
library(StAMPP)
#convert genotypes to dataframe of allele frequencies
freq <- stamppConvert(df2,"r")
#compute pairwise Fst based on population allele frequencies 
Fst <- stamppFst(freq, nboots=100, percent=95, nclusters=1)
write.csv(Fst$Fsts, "stampp_Fsts.csv")

###Compare Weir and Cockerham (1984) from StAMPP to modified Weir and Hill (2002) from Bedassle
setwd("C:/Users/Jeremiah/Documents")
df <- read.csv("pairwise_FST_May13.csv", header=T, stringsAsFactors=F)
df <- df[,-1]
rownames(df)<- colnames(df)

have <- read.csv("pairwise_haversine_dist.csv", header=T, stringsAsFactors = F)
have <- have[,-1]
rownames(have)<-colnames(have)

stamp <- read.csv("stampp_Fsts.csv", header=T, stringsAsFactors = F)
stamp <- stamp[,-1]
rownames(stamp) <- colnames(stamp)
for(i in 1:ncol(stamp)) {for(j in 1:i) {stamp[j,i]=stamp[i,j] }}
for(i in 1:nrow(stamp)) {stamp[i,i] <- 0}

pops <- c("AL2012",  "AL79", "ALBG", "AR125", "AR56", "IA10", "IN46", "IN77", "KS60", "KY51", "MI126", "MI127", "MN117", "MN118", "MO115", "MO116", "MO49", "MO57", "OH119", "OH64", "OK61", "PA27", "TN19", "WI128")
#df <- df[,unlist(lapply(pops, function(x) grep(x,colnames(df))))]
#df <- df[rownames(df) %in% pops,]
have <- have[,unlist(lapply(pops, function(x) grep(x,colnames(have))))] 
have <- have[rownames(have) %in% pops,]
stamp <- stamp[,unlist(lapply(pops, function(x) grep(x,colnames(stamp))))]
stamp <- stamp[rownames(stamp) %in% pops,]

#converting matrices to dataframe
new.names <- matrix(apply(combn(ncol(stamp), 2), 2, function(x) c(colnames(stamp)[x[1]], colnames(stamp)[x[2]])),ncol=2,byrow=T)
#pair_df <- data.frame( t(combn(names(df),2)) , dist=t(df)[lower.tri(df)])
stamp_df <- data.frame( t(combn(names(stamp),2)) , dist=t(stamp)[lower.tri(stamp)])
have_df <- data.frame( t(combn(names(have),2)) , dist=t(have)[lower.tri(have)])
new.df <- data.frame(cbind(new.names,as.numeric(stamp_df[,3])))
#new.df <- data.frame(cbind(new.df,as.numeric(stamp_df[,3])))
new.df <- data.frame(cbind(new.df,as.numeric(have_df[,3])))
colnames(new.df) <- c("pop1","pop2","stamp_fst","distance")
write.csv(new.df, "FST_2020.csv")

#plotting
df <- read.csv("FST_2020.csv")
df <- df[,-1]
#plot(df$fst ~ df$stamp_fst, xlab="W+C(1984) Fst", ylab="W+H(2002) Fst")
#plot(df$fst ~ df$distance)

#creat matrix (table s1)
pops <- c("AL2012", "AL79", "ALBG", "AR125", "AR56", "IA10", "IN46", "IN77", "KS60", "KY51", "MI126", "MI127", "MN117", "MN118", "MO115", "MO116", "MO49", "MO57", "OH119", "OH64", "OK61", "PA27", "TN19", "WI128")
stamp[upper.tri(stamp)] <- have[upper.tri(have)]

plot(df$stamp_fst ~ df$distance, xlab="geographic distance (km)", ylab="W+C(1984) Fst")
df <- df[order(df$stamp_fst),]

##IBD plotting by watershed
ibd <- df
#insert df$flag column, with colors written in the character column 
east <- c("AL2012", "AL79", "ALBG",  "IN46", "KY51", "OH119", "OH64", "PA27", "TN19")
west <- c("AR125", "AR56", "IA10", "IN77", "KS60", "MI126", "MI127", "MN117", "MN118", "MO115", "MO116", "MO49", "MO57", "OK61", "WI128")
df$reg1 <- 0
df$reg2 <- 0
df[df$pop1 %in% east,]$reg1 <- 1
df[df$pop2 %in% east,]$reg2 <- 1
df$reg1 <- as.numeric(df$reg1, na.rm=T)
df$reg2 <- as.numeric(df$reg2, na.rm=T)
df$compare <- df$reg1+df$reg2
df$compare <- as.factor(df$compare)

#average population pairwise fst 
#who's the most isolated?
results <- data.frame(pop=character(), avgf=numeric(), stringsAsFactors = F)
for(i in unique(pops)){
  working <- subset(df, (df$pop1==i)|(df$pop2==i))
  results[nrow(results)+1,]<-c(i, mean(working$fst))
}

results$avgf <- as.numeric(results$avgf)
results <- results[order(-results$avgf),]

#color points in IBD plot by regional (compare) or wisc factors
plot(df$distance, df$stamp_fst, pch=19, col=c("black", "green", "red")[df$compare], xlab="pairwise geographic distance (km)", ylab="W+C(1984) Fst")
df <- df[order(-df$stamp_fst),]

#Mantel tests
library(ecodist)
mantel(df$stamp_fst~df$distance, mrank=F, nboot=1000, pboot=0.9, cboot=0.95)

eastonly <- subset(df, df$compare==2)
westonly <- subset(df, df$compare==0)
across <- subset(df, df$compare==1)
within <- subset(df, df$compare!=1)
mantel(eastonly$stamp_fst~eastonly$distance, mrank=F, nboot=1000, pboot=0.9, cboot=0.95)
mantel(westonly$stamp_fst~westonly$distance, mrank=F, nboot=1000, pboot=0.9, cboot=0.95)

plot(eastonly$distance, eastonly$stamp_fst, pch=19, xlab="pairwise distance", ylab="pairwise fst")
plot(westonly$distance, westonly$stamp_fst, pch=19, xlab="pairwise distance", ylab="pairwise fst")
plot(across$distance, across$fst, pch=19, xlab="pairwise distance", ylab="pairwise fst")

