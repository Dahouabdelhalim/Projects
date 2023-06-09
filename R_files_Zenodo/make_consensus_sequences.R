#this R script generates consensus sequences for each individual in our dataset
library(dplyr)

# generate .phy files from df4
df <- read.csv("~/df4_0.50_revised.csv")
#df <- read.csv("df4_0.50_revised.csv")

#paste together pop and ind
df$plant <- paste(df$Pop, df$Ind, sep=":")

#sort by plant
df <- df[order(df$plant),]

# unique individuals
plant <- unique(df$plant)
loci <- unique(df$Locus)

#counter for locus, because later we will be writing and cannot use colons
genome <- 0

#summary on variability across all loci
summary <- data.frame(locus=character(), S=numeric(), ntaxa=numeric(), stringsAsFactors=F)

# generate .phy file at each locus
for(i in loci)
{
  genome <- genome+1
  #consider one locus at a time
  results <- data.frame(ID=character(),seq=character(), stringsAsFactors = F)
  #character vector to remove SNPs that are not biallelic or monomorphic
  scrub <- vector()
  working <- subset(df, df$Locus==i)
  phy <- data.frame(working$ID, working$Seq)
  phy[,1] <- as.character(phy[,1])
  phy[,2] <- as.character(phy[,2])
  
  #count the length of sequences at a locus
  nseq <- length(strsplit(phy[1,2], split="")[[1]])
  
  #generate consensus for each individual, including IUPAC ambiguity
  #first, all haplotypes for an individual have the same name
  phy[,1] <- gsub("\\\\..*","", phy[,1] )

  for(cz in unique(plant))
  {
    phy[grep(cz, phy[,1]),1] <- cz 
  }

   #consider one plant at a time, determing consensus for each individual
  for(j in unique(phy[,1]))
  {
    
    #single individual slice of dataset
    slice <- subset(phy, phy[,1]==j)
    
    #new df for doing math
    snps <- data.frame(matrix(NA, nrow = nrow(slice), ncol = nseq))
    
    #snp array for a single individual
    for(k in 1:nrow(snps))
    {
      for(a in 1:ncol(snps))
      {
            snps[k,a] <- strsplit(slice[k,2], split="")[[1]][a]
      }
    }
    
    #eliminate taxa with 5%+ NNNs
    lastcol <- ncol(snps)+1
    snps$consensus <- 0
    for(l in 1:nrow(snps))
    {
      snps[l,lastcol] <- sum(snps[l,]=="N")/nseq
    }
    snps <- subset(snps, snps$consensus<0.05)
    snps <- snps[,-lastcol]
    
    #logic to find the most common SNP in each position, including ambiguity for heterozygotes
    #ambiguous IUPAC code: 0.25 < frequency(minor allele) < 0.50
    
    if(nrow(snps)>0)
    {
    major <- nrow(snps)+1
    for(b in 1:ncol(snps))
    {
          
        if(dim(table(snps[,b]))==1)
        {
          snps[major,b] <- names(which.max(table(snps[,b])))
        }
      
          maxnum <- sum(snps[,b]==names(which.max(table(snps[,b]))),na.rm = T)
          minnum <- sum(snps[,b]==names(which.min(table(snps[,b]))),na.rm = T)

        if(dim(table(snps[,b]))>1&minnum/(minnum+maxnum)<=0.25)
        {
          snps[major,b] <- names(which.max(table(snps[,b])))
        }
          
        if(dim(table(snps[,b]))>1&minnum/(minnum+maxnum)>0.25)
        {
          if ( sort(unique(snps[,b]))[1]=="A" & sort(unique(snps[,b]))[2]=="C" )
          {snps[major,b] <- "M"}
          
          if ( sort(unique(snps[,b]))[1]=="A" & sort(unique(snps[,b]))[2]=="G" )
          {snps[major,b] <- "R"}
          
          if ( sort(unique(snps[,b]))[1]=="A" & sort(unique(snps[,b]))[2]=="T" )
          {snps[major,b] <- "W"}
          
          if ( sort(unique(snps[,b]))[1]=="C" & sort(unique(snps[,b]))[2]=="G" )
          {snps[major,b] <- "S"}
          
          if ( sort(unique(snps[,b]))[1]=="C" & sort(unique(snps[,b]))[2]=="T")
          {snps[major,b] <- "Y"}
          
          if ( sort(unique(snps[,b]))[1]=="G" & sort(unique(snps[,b]))[2]=="T" )
          {snps[major,b] <- "K"}
        }
    }
    
    ########turn all "NA" into "N"
    snps[is.na(snps)] <- "N"
    
    #now I just want to take a single line for each individual
    newrow <- nrow(results)+1
    results[newrow,1]<- j
    results[newrow,2]<- paste(snps[nrow(snps),], sep = "", collapse="")    
  }}
  
  #determine number of positions that are variable across all individual consensus sequences
  if(nrow(results)>0)
  {
  filter <- data.frame(matrix(NA, nrow = nrow(results), ncol = nseq))
  
  #snp array across populations
  for(ggg in 1:nrow(filter))
  {
    for(hhh in 1:ncol(filter))
    {
      filter[ggg,hhh] <- strsplit(results[ggg,2], split="")[[1]][hhh]
    }
  }
  
  #counter for segregating sites and ensuring sites are at most biallelic
  counter <- 0
  for(x in 1:ncol(filter))
  {
    #identify segregating sites
    if(dim(table(filter[,x]))>1)
    {
      counter <- counter+1
      
      #identify sites with >3 IUPAC codes for later removal
      if(length(unique(filter[,x]))>3)
      {
        scrub[length(scrub)+1]<- c(x)
        #ignore the position by filling with N
        filter[,x]<-"N"
      }
      
      #ensure that only 2 SNPs are segregating -- check IUPAC set
      if(length(unique(filter[,x]))>1)
      {
        test <- unique(filter[,x])
        #identify IUPAC ambiguity code with variable poly
        poly <- test[which(test %in% c("M","R","W","S","Y","K"))]
        #now remove polymorphic code(s) from set
        test <- test[-which(test %in% c("M","R","W","S","Y","K"))]
       
        #having >1 polymorphic code means the position is not biallelic
        if(length(unique(poly))>1)
        {
          scrub[length(scrub)+1]<- c(x)
          #ignore the position by filling with N
          filter[,x]<-"N"
        }
        #look at positions with a poly code, check to make sure it comports with actual nucleotides
        if(length(unique(poly))==1)
        {
          #M is A/C
          if(poly=="M")
          {
            if(all(test %in% c("A","C"))==FALSE){
            scrub[length(scrub)+1]<- c(x)
            #ignore the position by filling with N
            filter[,x]<-"N"}
          }
          #R is A/G
          if(poly=="R")
          {
            if(all(test %in% c("A","G"))==FALSE){
            scrub[length(scrub)+1]<- c(x)
            #ignore the position by filling with N
            filter[,x]<-"N"}
          }
          #W is A/T
          if(poly=="W")
          {
            if(all(test %in% c("A","T"))==FALSE){
            scrub[length(scrub)+1]<- c(x)
            #ignore the position by filling with N
            filter[,x]<-"N"}
          }
          #S is C/G
          if(poly=="S")
          {
            if(all(test %in% c("C","G"))==FALSE){
            scrub[length(scrub)+1]<- c(x)
            #ignore the position by filling with N
            filter[,x]<-"N"}
          }
          #Y is C/T
          if(poly=="Y")
          {
            if(all(test %in% c("C","T"))==FALSE){
            scrub[length(scrub)+1]<- c(x)
            #ignore the position by filling with N
            filter[,x]<-"N"}
          }
          #K is G/T
          if(poly=="K")
          {
            if(all(test %in% c("G","T"))==FALSE){
            scrub[length(scrub)+1]<- c(x)
            #ignore the position by filling with N
            filter[,x]<-"N"}
          }
        }
      }
    }
    
    ###must collapse down the filter array so the CATG string is a single field
    for(ii in 1:nrow(results))
    {
    results[ii,2]<- paste(filter[ii,], sep = "", collapse="")    
    }
    
  }
    
  #update the summary file for segregating sites, remembering to remove positions that are not biallelic
  counter <- counter - length(scrub)
  ntaxa <- nrow(results)
  summary[nrow(summary)+1,] <- c(as.character(i),counter,ntaxa)

  ###-----------------###
  #prepare file for writing
  #test 1: counter >0 & counter<5
  #test 2: counter >0 & counter>=5
  
  if(nrow(results)>0)
  {
  ntaxa <- nrow(results)
  #write the file 
  colnames(results) <- NULL
  rownames(results) <- c()
  opts <- options(useFancyQuotes = FALSE)
  on.exit(options(opts))
  h1 <- paste(c(ntaxa))
  h2 <- paste(c(nseq))
  header <- c(h1, h2)
  results <- rbind(header, results)
  write.table(results, paste(c(genome),".phy",sep=""), append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)  
  }
  }}
#summary of segregating sites for each RAD locus
write.csv(summary, "S.csv")

library(phylotools)

#set working directory
setwd("C:/Users/Jeremiah/Desktop/phy")

#take summary file and then do some QC
df <- read.csv("S.csv")

#keep only those RADloci with >2/3 representation of individuals
df <- subset(df, df$ntaxa>115)

noms <- unique(df$X)
for(i in 1:length(unique(noms)))
{
  noms[i] <- paste(noms[i], ".phy", sep="")
}

#concatenate .phy files into a single .phy file
#files <- list.files(pattern = ".phy")

#make fasta files
for(i in 1:length(noms))
{
  tempy <- read.phylip(noms[i])
  nouveau <- gsub("\\\\..*","", noms[i])
  dat2fasta(tempy, outfile=paste(nouveau, ".fasta", sep="", collapse=""))
}

#list fasta files and concatenate loci
lister <- list.files(pattern=".fasta")
supermat(infiles= lister, outfile="supermat.phy", partition.file="partitions.txt")

