##R script to convert allelic information into genotype
##Follows MAC-PR method from Esselink et al. 2004
##allele size from PeakScanner is followed by area under the peak
##for data input format, see "microsatellite data.csv" file

#### File import and cleaning section
# Read in .csv
df <- read.csv("path/filename.csv")

# Check for factors
str(df)

# Set any NA to 0
df[is.na(df)] <- 0

#### Scoring allelic richness within a population
df.size<-df[,grep("size",colnames(df))]
x <- list(1:4,5:8,9:12,13:16)
rich <- do.call(cbind,lapply(x,function(i) rowSums(df.size[, i]!=0)))
df <- cbind(df,rich=rich)

#### 1 allele - copy it over to all four slots (homozygote)
rich.1 <- (rich[,rep(1:ncol(rich),each=4)]==1)*df.size
empty <- c(2:4,6:8,10:12,14:16)
first <- rep(c(1,5,9,13),each=3)
df[,colnames(rich.1[,empty])] <- df[,colnames(rich.1[,empty])] + rich.1[,first]


#### 2 alleles - 1) find ratio of locusX.1.size / (locusX.1.size + locusX.2.size)
# specify delta centered around 0.5 to determine bins
bin <- 1/3

# create some working arrays
rich.2 <- (rich[,rep(1:ncol(rich),each=4)]==2)*df.size
area.2 <- (rich[,rep(1:ncol(rich),each=4)]==2)*df.area
ratios <- do.call(cbind,lapply(x,function(i) area.2[,i]/rowSums(area.2[, i])))
ratios[is.na(ratios)] <- 0

# make an array that identifies which allele to choose. Decision logic is all in the second line below
z <- ratios[,c(1,5,9,13)]
z <- cbind(ifelse(z>bin,1,2), ifelse(z>(2*bin),1,2))
z <- z[,order(colnames(z))]

# choose arrays from the rich.2 table based on the alleles identified in x
result <- data.frame(1:nrow(z))
t<- rich.2[,c(1,2,5,6,9,10,13,14)]

for(i in seq(1,8,by=2)) {
  a <- z[,i] + (i-1)
  b <- z[,i+1] + (i-1)
  result <- cbind(result, t[cbind(seq_along(a), a)] , t[cbind(seq_along(b), b)] )
}

result<-result[,-1]
colnames(result) <- colnames(rich.2[,c(3,4,7,8,11,12,15,16)])
df[,colnames(result)] <- result + df[,colnames(result)]

#### 3 alleles - copy the allele with the highest area to the empty slot
# create some working arrays
rich.3 <- (rich[,rep(1:ncol(rich),each=4)]==3)*df.size
df.area<-df[,grep("area",colnames(df))]
area.3 <- (rich[,rep(1:ncol(rich),each=4)]==3)*df.area

# identify position of maximum in each subset of the vector
x <- list(1:4,5:8,9:12,13:16)
maximums <- do.call(cbind,lapply(x,function(i) max.col(area.3[,i])))

# choose the size that matches the max indices found in area
rich.3$locus1.4.size <- rich.3[,1:4][cbind(seq_along(maximums[,1]), maximums[,1])]
rich.3$locus2.4.size <- rich.3[,5:8][cbind(seq_along(maximums[,2]), maximums[,2])]
rich.3$locus3.4.size <- rich.3[,9:12][cbind(seq_along(maximums[,3]), maximums[,3])]
rich.3$locus4.4.size <- rich.3[,13:16][cbind(seq_along(maximums[,4]), maximums[,4])]

# insert these values into the appropriate place in the larger dataset
rich.3 <- rich.3[,grep("\\\\.4",colnames(rich.3))]
df[,colnames(rich.3)] <- rich.3 + df[,colnames(rich.3)]

#remove area information, since not relevant now to final genotype
df.final <- df[,-grep("area",colnames(df))]

head(df.final)
write.csv(df.final, file="genotypes.csv")