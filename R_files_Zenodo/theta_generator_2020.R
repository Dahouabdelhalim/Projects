# load vcf
#setwd("C:/Users/Nate/Dropbox/Academia/WSU/Campanula/")
setwd("C:/Users/Jeremiah/Documents")

df <- read.csv("bellflower.vcf", header=T, sep="\\t", stringsAsFactors=F)

df$POS <- as.numeric(df$POS)
colnames(df)[1] = "CHROM"
df$CHROM <- gsub(":", "_", df$CHROM)

# Simplify formatting for allele calls
df <- as.data.frame(apply(df,2,function(x) sub(".*:","",x)),stringsAsFactors=F) 

# Subset df to include 24 focal populations
# Fed one population at a time
pops <- c("AL2012", "AL79", "ALBG", "AR125", "AR56", "IA10", "IN46", "IN77", "KS60", "KY51", "MI126", "MI127", "MN117", "MN118", "MO115", "MO116", "MO49", "MO57", "OH119", "OH64", "OK61", "PA27", "TN19", "WI128")
#, "GA22", "IA45", "IA8", "ILCB", "MINorth", "MS55", "NE59")

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
#length <- max(working$position)-min(working$position)
#correction now that L is 85 for all contigs

getL <- function(contig) {
  working <- key[key$contig==contig,]
  length <- 85
  return(length)
}

# This function calculates Watterson's correction factor a, assuming tetraploidy
getA <- function(popContig) {
  
  if(sum(!is.na(popContig))<1) { return(NA) }
  
  # Add number of individuals with a non-zero C,T,A, or G
  n <- NA
  memContig <- popContig
  popContig <- apply(popContig, c(1,2), sum, na.rm=T)
  k <- rowSums(popContig>0, na.rm=T)
  index <- 4*k
  index[index<1]<- NA
  index <- index[!is.na(index)]
  n <- mean(index, na.rm=T)
  
  if(!is.na(n)){
    #these depend on n
    a <- sum(1/1:(n-1)) }
  return(a)
}

getS <- function(popContig) {
  
  # Make sure popContig isn't an empty vector
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
  pi <- sum(df_pi$pi_c, na.rm=T)
  return(pi)
}

###### tajima's D ########
get_tajimasD <- function(popContig){
  
  # Make sure popContig isn't an empty vector
  if(!is.null(nrow(popContig))) { 
    
    # Add number of individuals with a non-zero C,T,A, or G
    n <- NA
    memContig <- popContig
    popContig <- apply(popContig, c(1,2), sum, na.rm=T)
    k <- rowSums(popContig>0, na.rm=T)
    index <- 4*k
    index[index<1]<- NA
    index <- index[!is.na(index)]
    n <- mean(index, na.rm=T)
    
    if(!is.na(n)){
      #these depend on n
      a1 <- sum(1/1:(n-1))
      a2 <- sum(1/(1:(n-1))^2)
      b1 <- (n+1)/(3*(n-1))
      b2 <- 2*(n^2+n+3)/(9*n*(n-1))
      #these depend on summaries of n
      c1 <- b1 - 1/a1
      c2 <- b2 - ((n+2)/(a1*n))+(a2/a1^2)
      e1 <- c1/a1
      e2 <- c2/(a1^2+a2)
    }
  } 
  
  #getS
  if(!is.null(nrow(memContig))) { 
    
    # Add read depth for individuals together across all sites
    # results in a 2D matrix where rows are sites within the contig 
    # and columns are aggregate CATG read depth counts
    S <- apply(memContig,c(1,3),sum,na.rm=T)
    
    # Count the number of nucleotide positions with data for each site
    S <- rowSums(S>0, na.rm=T)
    
    # How many of those are polymorphic?
    S <- sum(S>1, na.rm=T)
  } 
  tajimasd <- NA
  if(!is.null(nrow(memContig))){
    if(S>0){
    tajimasd <- sqrt(e1*S + e2*S*(S-1))}} 
  return(tajimasd)
}

############################################################################
# This function gets a row of data
getData <- function(currentPop, currentContig) {
  # Working on a population-contig pair
  #print(c(currentPop,currentContig))
  slice <- df.array[key$contig==currentContig,grep(currentPop,colnames(df.array)),]
  #slice <- aperm(slice, c(3,2,1))
  
  # print(head(slice))
  # Error checking. If there are no sites that contain data for
  # a particular population / contig pair the array will be reduced by
  # by a dimension. If so bail early.
  if(length(dim(slice)) < 3 | nrow(slice) < 3) return(c(NA,NA,NA,NA,NA,NA,NA))
  
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
  
  if(is.na(L)) return(c(NA,NA,NA,NA,NA,NA,NA))
  meanRD <- mean(rd.depth,na.rm=T)
  thetaW <- S/(L*a)
  thetaPi <- tajima(slice, rd.depth)/L
  #new code to calculate tajima's D
  tajsd <- get_tajimasD(slice)
  #Tajima's D covers the whole locus, not each site
  tD <- 85*(thetaPi-thetaW)/tajsd
  return(c(meanRD, L, S, a, thetaW, thetaPi, tD))
} 

results <- data.frame( pop=rep(pops, each=length(unique(df$CHROM))), contig=rep(unique(df$CHROM),length(pops)) )
d <- cbind(results,t(apply(results,1, function(x) getData(x[1],x[2]))))
colnames(d) <- c("pop", "contig", "meanRD","L", "S", "a", "thetaW", "thetaPi", "tajimasD")
write.csv(d, "Apr27_20_pi.csv",row.names=F)

###Now generate population means


#to go straight to results:
df <- read.csv("Apr27_20_pi.csv", header=T, stringsAsFactors=F)
#df <- subset(df,df$pop!="WI128")

#load file and name columns
df <- df[order(df$pop),]

#eliminate contigs with missing data
#df <- na.omit(df)

##FILTERING
#make a table on how many times a specific value occurs in dataset (here it is ID)
#clear plot history
b <- as.data.frame(table(df$contig), stringsAsFactors=F)
hist(b$Freq, breaks=20, xlab="Number of Populations Present", ylab="Contig Frequency", main="Contig Distribution")

#retain contigs present in some threshold of populations
b <- subset(b,b$Freq>23)
hist(b$Freq, breaks=20, xlab="Number of Populations Present", ylab="Contig Frequency", main="Contig Distribution")

#in output, keep only those IDs that pass above criterion
df <- subset(df, df$contig %in% b$Var1)

#generating population means
results <- data.frame(pop=character(), contigs=numeric(), meanRD=numeric(), length=numeric(), a=numeric(), thetaw=numeric(), thetapi=numeric(), tajimaD=numeric(), stringsAsFactors=F)

for(reg in unique(df$pop))
{
  working <- subset(df, df$pop==reg)
  contigs <- nrow(working)
  results[nrow(results)+1,] <- c(reg, contigs, mean(working$meanRD), mean(working$L), mean(working$a), mean(working$thetaW, na.rm=T), mean(working$thetaPi, na.rm=T), mean(working$tajimasD, na.rm=T))		
}

#add other variables such as distance from Appalachian origin, distance from Miss. origin, latitude, longitude, autonomy, and outcrossing rate (from other papers)
results$origin <- c(574.4121,625.2450,378.0876,830.5776,709.0799,1036.8446,346.6879,331.4669,1075.0405,121.0806,610.5049,606.2039,1194.9849,1366.4832,775.4209,691.0662,687.4743,774.9982,321.7093,487.1012,1058.3564,533.5685,432.0201,879.9210)
results$origin2 <- c(775.6419, 630.9584, 520.4695, 358.6388, 272.8686, 514.0102, 321.9157, 262.5629, 484.9015, 505.0224, 586.8635, 486.4209, 769.9644, 884.1728, 185.2947, 215.2572, 97.6708, 204.9848, 542.3751, 784.2143, 645.0784, 895.3698, 340.8064, 527.5753)
results$long <- c(-85.95084, -88.20820, -86.51697, -92.71624, -91.38649, -93.6725, -86.39908, -87.00707, -95.52129, -84.24557, -85.34209, -86.58336, -93.19248, -95.88795, -92.01166, -91.26709, -91.1036, -92.21757, -83.9967, -81.5181, -94.5669, -80.08326, -88.0687, -90.04498)
results$lat <- c(32.26718, 32.9293, 34.65581, 36.03251, 36.22686, 42.0728, 39.14864, 38.10728, 38.97137, 37.93824, 42.32093, 41.92333, 44.90078, 45.02566, 38.92989, 36.76468, 38.4711, 37.85376, 39.88547, 41.1147, 33.9464, 41.00795, 35.7583, 43.14964)
results$auto <- c(0.54176231, 0.41465841, 0.33812854, 0.70378968, 0.4801948, 0.75522123, 0.50690206, 0.46653918, 0.59512501, 0.44181021, 0.62231644, 0.7642207, 0.60510913, 0.74235381, 0.56691715, 0.58925164, 0.51668144, 0.57978036, 0.51983586, 0.47999632, 0.57208208, 0.35789006, 0.56829612, 0.64201384)
results$tm <- c(0.797, 0.706, 0.818, 0.889, NA, 0.564, 0.67, 0.806, 0.543, 0.894, 0.851, 0.729, 0.787, 0.947, 0.695, 0.888, 0.823, 0.883, 0.882, 0.908, 0.911, 0.947, 0.785, 0.976) 


