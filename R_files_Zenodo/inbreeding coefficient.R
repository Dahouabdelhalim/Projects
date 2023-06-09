# Upload vcf
df <- read.csv("bellmerge.vcf", header=T, sep="\\t", stringsAsFactors=F)

# Drop extraneous columns
# ID,REF,ALT,QUAL,FILTER,INFO,FORMAT
# Format = CATG. Strip off "GT:DP:"
df <- df[,-(3:9)]

# Simplify formatting for allele calls
df <- as.data.frame(apply(df,2,function(x) sub(".*:","",x)),stringsAsFactors=F)
df$POS <- as.numeric(df$POS)
rownames(df) = make.names(paste(df$CHROM,df$POS,sep="-"), unique=TRUE)

# Subset df to include 24 focal populations
pops <- c("AL2012", "AL79", "ALBG", "AR125", "AR56", "IA10", "IN46", "IN77", "KS60", "KY51", "MI126", "MI127", "MN117", "MN118", "MO115", "MO116", "MO49", "MO57", "OH119", "OH64", "OK61", "PA27", "TN19", "WI128", "MD5", "WV72", "VA73", "NC36", "NC90", "VA71")
df <- df[,unlist(lapply(pops, function(x) grep(x,colnames(df))))] 

# Decompose combined CATG strings into separate calls 
string.cutter <- function(string.vector,position) {
  result <- lapply(string.vector,function(y) strsplit(y, split=",")[[1]][position])
  return(as.numeric(result))
}

# Make a 3 dimensional array of reads.
# Arrays are referenced by [row,column,depth]
df.array.names <- c("C","A","T","G")
df.array.vector <- unlist(lapply(1:4, function (pos) apply(df, 2, function(x) string.cutter(x,pos))))
df.array.dims <- c(nrow(df),ncol(df),4)
df.array <- array(df.array.vector,dim=df.array.dims,dimnames=list(rownames(df),colnames(df),df.array.names))

# Clean up. Remove large objects to free up memory
rm(df.array.vector)

# Eliminate errors (read count < 5%) -- applied at the individual level
for(i in 1:ncol(df.array))
{
  df.array[,i,] <- df.array[,i,] * (df.array[,i,]/rowSums(df.array[,i,]) > 0.05)
}

# Reduce the array down to two dimensions to get read depth at each site. 
df.read.depth <- df.array[,,1] + df.array[,,2] + df.array[,,3] + df.array[,,4]

#remove monomorphic sites and those with >2 nucleotides
counter <- apply(df.array, 1, function(x) sum(colSums(x, na.rm=TRUE)>0))
df.array <- df.array[counter==2,,]
df.read.depth <- df.read.depth[counter==2,]

# If total read depth is below a threshold for an individual, do not use that site
df.read.depth[df.read.depth <= 20] <- NA
df.array[df.read.depth <= 20] <- NA

#convert from 3D to 2D, by using MAJOR/MINOR rule for SNPs
#Major nucleotide at a polymorphic site is kept for all individuals
ref <- data.frame(C=apply(df.array[,,1],1,sum, na.rm=TRUE),A=apply(df.array[,,2],1,sum,na.rm=TRUE),T=apply(df.array[,,3],1,sum,na.rm=TRUE),G=apply(df.array[,,4],1,sum, na.rm=TRUE)) 
ref$major.pos <- max.col(ref, ties.method = "random")
ref$major.call <- colnames(ref)[ref$major.pos]
df.major <- data.frame(apply(df.array,2,function(x) x[cbind(seq_along(ref$major.pos), ref$major.pos)]),row.names=rownames(df.array[,,1]))

#proportions expressed as proportion of the major read
props <- df.major/df.read.depth

#generate matrix of tetraploid genotypes for individuals (2 = major allele)
#for loops take a while 
genotype <- data.frame(matrix(27, nrow=nrow(props), ncol=ncol(props)))
dimnames(genotype) <- dimnames(props)
for(i in 1:nrow(genotype)) {
  for(k in 1:ncol(genotype)) {
    if((is.na(props[i,k]))>0) { genotype[i,k] <- NA}
    if((!is.na(props[i,k])>0)&((props[i,k]<=0.05)>0)) { genotype[i,k] <- 0}
    if((!is.na(props[i,k])>0)&((props[i,k]>0.05)>0)&((props[i,k]<=0.33)>0)) { genotype[i,k] <- 1}
    if((!is.na(props[i,k])>0)&((props[i,k]>0.33)>0)&((props[i,k]<=0.67)>0)) { genotype[i,k] <- 2}
    if((!is.na(props[i,k])>0)&((props[i,k]>0.67)>0)&((props[i,k]<=0.95)>0)) { genotype[i,k] <- 3}
    if((!is.na(props[i,k])>0)&((props[i,k]>0.95)>0)) { genotype[i,k] <- 4}
  }
}

#individuals in rows, with alleles in columns
major.counts <- as.matrix(t(genotype))
#save genotypic data
write.csv(major.counts, "genotype_array.csv")

#each row is an individual
df <- read.csv("genotype_array.csv",header=T, stringsAsFactors = F)

#name sample column, then keep plants from the 24 populations
df <- df[1:167,]

#force read proportions into genotype calls (0, 0.25, 0.50, 0.75, and 1)
df <- as.data.frame(df)
ids <- df[,1]
df <- df[,-1]
df[] <- lapply(df, as.numeric)
df[,1:ncol(df)] <- df[,1:ncol(df)]/4
df <- t(df)
colnames(df) <- ids

#subset into populations
pops <- c("AL2012", "AL79", "ALBG", "AR125", "AR56", "IA10", "IN46", "IN77", "KS60", "KY51", "MI126", "MI127", "MN117", "MN118", "MO115", "MO116", "MO49", "MO57", "OH119", "OH64", "OK61", "PA27", "TN19", "WI128")

#generate dataframe of allele frequencies for each population
pq <- data.frame(matrix(NA, ncol = length(pops), nrow = nrow(df))) 
colnames(pq)<-pops
z<-0
for(i in pops) {
z<-z+1
slice <- df[ , grepl(i , colnames( df ) ) ]
Ho <- rowMeans(slice, na.rm=T)
Ho[is.nan(Ho)]<- NA
pq[,z]<-Ho
}

#generate information on sample size for each population
samples <- data.frame(matrix(NA, ncol = length(pops), nrow = nrow(df))) 
colnames(samples)<-pops
z<-0
for(i in pops) {
  z<-z+1
  slice <- df[ , grepl(i , colnames( df ) ) ]
  numbers<- rowSums(!is.na(slice))
  samples[,z]<-numbers
}

#convert p,q information for expected heterozygosity
#see Meirmans 2018 for defense of "gene diversity" approach in a polyploid
he <- pq
for(i in 1:nrow(pq)){
  for(j in 1:ncol(pq)){
   he[i,j]<-1-(pq[i,j]^2)-(1-pq[i,j])^2
  }}

#generate dataframe of observed heterozygosities (Hs)
#Moody et al. 1993 for gametic heterozygosity of tetraploid genotypes
#alleles sampled without replacement during movement into gametes
hs <- df
for(i in 1:nrow(df)){
  for(j in 1:ncol(df)){
    if(!is.na(df[i,j])&(df[i,j]==0.5)){hs[i,j]<-0.67}    
    if(!is.na(df[i,j])&(df[i,j]==0.75)){hs[i,j]<-0.5}
    if(!is.na(df[i,j])&(df[i,j]==0.25)){hs[i,j]<-0.5}
    if(!is.na(df[i,j])&(df[i,j]==1)){hs[i,j]<-0}
    if(!is.na(df[i,j])&(df[i,j]==0)){hs[i,j]<-0}
  }}

#generate dataframe of observed heterozygosities
hobs <- data.frame(matrix(NA, ncol = length(pops), nrow = nrow(df))) 
colnames(hobs)<-pops
z<-0
for(i in pops) {
  z<-z+1
  slice <- hs[ , grepl(i , colnames( hs ) ) ]
  Ho <- rowMeans(slice, na.rm=T)
  Ho[is.nan(Ho)]<- NA
  hobs[,z]<-Ho
}

#now go over the dataframe to get He, Ho
#conditioned on sample size of 5+ plants
results <- data.frame(pop=character(), He=numeric(), Ho=numeric(), stringsAsFactors=F)
z<-0
for(i in pops){
  z<-z+1
  workhe <- he[(samples[,z]>5),]
  workhobs <- hobs[(samples[,z]>5),]
  expected <- mean(workhe[,z], na.rm=T)
  observed <- mean(workhobs[,z], na.rm=T)
  results[nrow(results)+1,]<- c(i,expected,observed)
}
results$He <- as.numeric(results$He)
results$Ho <- as.numeric(results$Ho)
results$F <- 1-(results$Ho/results$He)
  
#examine per locus inbreeding coefficients, then average
fis <- 1-(workhobs/workhe)
fis[sapply(fis, is.nan)] <- NA

z<-0
for(i in pops){
  z<-z+1
  meanFis <- mean(fis[,z], na.rm=T)
  results[z,5]<-meanFis
}
colnames(results)[5]<-"Fis"

#save the output
write.csv(results, "inbreeding_coefficients.csv")

