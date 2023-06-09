#Script to quantify population-genetic diversity within populations
#load vcf
df <- read.csv("bellmerge.vcf", header=T, sep="\\t", stringsAsFactors=F)

# Simplify formatting for allele calls
df <- as.data.frame(apply(df,2,function(x) sub(".*:","",x)),stringsAsFactors=F)
df$POS <- as.numeric(df$POS)
colnames(df)[1]<- "CHROM"

# Subset df to include 24 focal populations
# Fed one population at a time
pops <- c("AL2012", "AL79", "ALBG", "AR125", "AR56", "IA10", "IN46", "IN77", "KS60", "KY51", "MI126", "MI127", "MN117", "MN118", "MO115", "MO116", "MO49", "MO57", "OH119", "OH64", "OK61", "PA27", "TN19", "WI128")
df <- cbind(df[,1:2],df[,unlist(lapply(pops, function(x) grep(x,colnames(df))))])

# Decompose combined CATG strings into separate calls 
string.cutter <- function(string.vector) {
  split <- strsplit(as.character(string.vector), split=",")
  return(split)
}

# Make a 3 dimensional array of reads.
# Arrays are referenced by [row,column,depth]
df.array.cols <- colnames(df[-c(1:2)])
df.array.rows <- rownames(df)
df.array.names <- c("C","A","T","G")
df.array.vector <- as.numeric(unlist(apply(df[,-c(1:2)], 1, string.cutter)))
df.array.dims <- c(4,ncol(df)-2,nrow(df))
df.array <- array(df.array.vector,dim=df.array.dims,dimnames=list(df.array.names,df.array.cols,df.array.rows))
df.array <- aperm(df.array)

# Identifies each row in df.array by contig and position
key <- data.frame(contig=df$CHROM,position=df$POS,stringsAsFactors=F)

# Eliminate errors (read count < 5 percent) -- applied at the individual level
for(i in 1:ncol(df.array))
{
  df.array[,i,] <- df.array[,i,] * (df.array[,i,]/rowSums(df.array[,i,]) > 0.05)
}

# Reduce the array down to two dimensions to get read depth at each site. 
df.read.depth <- df.array[,,1] + df.array[,,2] + df.array[,,3] + df.array[,,4]

# If total read depth is below a threshold for individuals, do not use that site
df.read.depth[df.read.depth <= 20] <- NA

# This function calculates the length of sequence, L, given a population-contig pair
getL <- function(contig) {
  working <- key[key$contig==contig,]
  length <- max(working$position)-min(working$position)
  return(length)
}

# This function calculates Watterson's correction factor a, assuming tetraploidy
getA <- function(popContig) {
  
  if(sum(!is.na(popContig))<1) { return(NA) }
  
  # Add number of individuals with a non-zero C,T,A, or G
  popContig <- apply(popContig, c(1,2), sum, na.rm=T)
  
  k <- rowSums(popContig>0, na.rm=T)
  index <- 4*k
  index[index<1]<- NA
  index <- index[!is.na(index)]
  working<-data.frame(k=numeric(), a=numeric(), stringsAsFactors=F)
  a <- 0
  for(i in 1:length(index)) { 
    working[nrow(working)+1,] <- c(index[i],0)    
  }
  for(i in 1:nrow(working)){
    a<-0
    top <- working[i,1]-1
    for(j in 1:top){
      a=1/j+a
      working[i,2]<-a
    }
  }
  a<-mean(working[,2], na.rm=T)
  return(a)
}

getS <- function(popContig) {
  
  # Make sure popContig is not an empty vector
  if(is.null(nrow(popContig))) { return(NA) }
  
  # Add read depth for individuals together across all sites
  # results in a 2D matrix where rows are sites within the contig 
  # and columns are aggregate CATG read depth counts
  S <- apply(popContig,c(1,3),sum,na.rm=T)
  
  # Count the number of nucleotide positions with data for each site
  S <- rowSums(S>0, na.rm=T)
  
  # How many of those are polymorphic?
  S <- sum(S>1, na.rm=T)
  
  # Returns an integer representing the number of polymorphic sites
  # at the given contig
  return(S)
}

# This function calculates Tajima's theta
# need x1 and x2 across each SNP
tajima <- function(popContig, readDepth) {
  
  popContig[,,1] <- popContig[,,1]/as.matrix(readDepth)
  popContig[,,2] <- popContig[,,2]/as.matrix(readDepth)
  popContig[,,3] <- popContig[,,3]/as.matrix(readDepth)
  popContig[,,4] <- popContig[,,4]/as.matrix(readDepth)
  
  #### Genotyping #################
  
  # Within a single individual we need to identify how many copies of each allele
  # exist at an specific site. To do this we create allele frequency bins:
  bin.1 <- 0.05
  bin.2 <- 0.33
  bin.3 <- 0.67
  bin.4 <- 0.95
  
  # Now we use these bins to infer allele copy number. Because these are
  # tetraploids, alleles can either be present in 0,1,2,3, or 4 copies.
  popContig[popContig > bin.4] <- 4
  popContig[popContig > bin.3 & popContig <= bin.4] <- 3
  popContig[popContig > bin.2 & popContig <= bin.3] <- 2
  popContig[popContig > bin.1 & popContig <= bin.2] <- 1
  popContig[popContig <= bin.1] <- 0
  
  ################################
  
  ### Set up a data frame to calculate pi ###
  df_pi <- data.frame(rowSums(popContig,na.rm=T))
  popContig[is.na(popContig)] <- 0
  df_pi$k <- apply(popContig,1,function(x) max(colSums(x)))
  colnames(df_pi)[1]<-"n"
  df_pi$pi_c <- 2*df_pi$k*(df_pi$n-df_pi$k)/(df_pi$n*(df_pi$n-1))
  pi <- sum(df_pi$pi, na.rm=T)
  return(pi)
}

# This function gets a row of data
getData <- function(currentPop, currentContig) {
 
  # Working on a population-contig pair
  slice <- df.array[key$contig==currentContig,grep(currentPop,colnames(df.array)),]
  
  # If there are no sites that contain data fora particular population / contig pair 
  # the array will be reduced by a dimension. If so bail early.
  if(length(dim(slice)) < 3 | nrow(slice) < 3) return(c(NA,NA,NA,NA,NA,NA))
  
  rd.depth <- rowSums(slice, dims=2, na.rm=T)
  slice[,,1] <- replace(slice[,,1],rd.depth < 20, NA)
  slice[,,2] <- replace(slice[,,2],rd.depth < 20, NA)
  slice[,,3] <- replace(slice[,,3],rd.depth < 20, NA)
  slice[,,4] <- replace(slice[,,4],rd.depth < 20, NA)
  
  # Remove sites with more than 2 nucleotides segregating in a population
  test <- apply(slice, 1, function(x) sum(colSums(x)>0, na.rm=T))
  
  slice <- slice[test<3,,]
  rd.depth <- rowSums(slice, dims=2, na.rm=T)
  
  L <- getL(currentContig)
  S <- getS(slice)
  a <- getA(slice)
  
  if(is.na(L)) return(c(NA,NA,NA,NA,NA,NA))
  meanRD <- mean(rd.depth,na.rm=T)
  thetaW <- S/(L*a)
  thetaPi <- tajima(slice, rd.depth)/L
  return(c(meanRD, L, S, a, thetaW, thetaPi))
} 

results <- data.frame( pop=rep(pops, each=length(unique(df$CHROM))), contig=rep(unique(df$CHROM),length(pops)) )
d <- cbind(results,t(apply(results,1, function(x) getData(x[1],x[2]))))
colnames(d) <- c("pop", "contig", "meanRD","L", "S", "a", "thetaW", "thetaPi")

#save the results
write.csv(d, file="results.csv",row.names=F)

###---process the results to generate population-level estimates---###

#to go straight to results:
df <- read.csv("results.csv", header=T, stringsAsFactors=F)
#df <- subset(df,df$pop!="WI128")

#load file and name columns
df <- df[order(df$pop),]

#eliminate contigs with missing data
df <- na.omit(df)

#generate diff for numerator of Tajima's D
df$D <- df$thetaPi - df$thetaW

#find denominator of Tajima's D for each population
df$tajima_D <- 0
for(i in unique(df$pop))
{
  working <- subset(df, df$pop==i)
  denom <- sqrt(var(working$D, na.rm=T))
  df[df$pop==i,]$tajima_D <- working$D/denom
}

##FILTERING
#make a table on how many times a specific value occurs in dataset (here it is ID)
b <- as.data.frame(table(df$contig), stringsAsFactors=F)
hist(b$Freq, breaks=20, xlab="Number of Populations Present", ylab="Contig Frequency", main="Contig Distribution")

#retain contigs present in some threshold of populations, say 18
b <- subset(b,b$Freq>18)
hist(b$Freq, breaks=20, xlab="Number of Populations Present", ylab="Contig Frequency", main="Contig Distribution")

#in output, keep only those IDs that pass above criterion
df <- subset(df, df$contig %in% b$Var1)

#require 50+ reads 
df <- subset(df, df$meanRD>50)

#this 'L' is the distance when considered at the VCF level
df <- subset(df, df$L<1000)

#generating population means
results <- data.frame(pop=character(), contigs=numeric(), meanRD=numeric(), length=numeric(), a=numeric(), thetaw=numeric(), thetapi=numeric(), tajimaD=numeric(), stringsAsFactors=F)

for(reg in unique(df$pop))
{
  working <- subset(df, df$pop==reg)
  contigs <- nrow(working)
  results[nrow(results)+1,] <- c(reg, contigs, mean(working$meanRD), mean(working$L), mean(working$a), mean(working$thetaW, na.rm=T), mean(working$thetaPi, na.rm=T), mean(working$tajima_D, na.rm=T))		
}

#save the population-level estimates
write.csv(results, "Genetic summaries.csv")