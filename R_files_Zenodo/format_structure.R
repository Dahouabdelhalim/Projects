## Takes fineRADstructure input file and transforms it into STRUCTURE input file format

# Set working directory to location of fine_revised_MADed_0.005.in
setwd("~/Desktop/")

# Read in fine.in and transpose it so that individuals are rows loci are columns
data <- t(read.csv("fine_revised_MAFed_0.005.in", sep = '\\t', header = T, stringsAsFactors = F))

# Make a dataframe with four rows per individual to hold genotypes
require(gdata)
data2 <- interleave(data,data)
data3 <- interleave(data2,data2)

range <- (1:nrow(data) * 4) - 3

# Logic to split strings and only keep desired genotype position
for( i in range) {
  for( j in 1:ncol(data3)) {
    data3[i,j] <- strsplit(data3[i,j], '/')[[1]][1]
    data3[i+1,j] <- strsplit(data3[i+1,j], '/')[[1]][2]
    data3[i+2,j] <- strsplit(data3[i+2,j], '/')[[1]][3]
    data3[i+3,j] <- strsplit(data3[i+3,j], '/')[[1]][4]
  }
}

# Breaks apart phased locus into independant alleles
data4 <- as.data.frame(t(as.data.frame(strsplit(data3[,1],""),stringsAsFactors=F)), stringsAsFactors = F)
for(i in 2:ncol(data3)) {
  data4 <- cbind(data4, as.data.frame(t(as.data.frame(strsplit(data3[,i],""),stringsAsFactors=F)), stringsAsFactors=F))
}

# Replace CATG with integer values
data4[data4=="C"] <- 1
data4[data4=="A"] <- 2
data4[data4=="T"] <- 3
data4[data4=="G"] <- 4
data4[is.na(data4)] <- -9
data4[data4==""] <- -9

# Make column of individual and population names
data5 <- cbind(rep(rownames(data),each=4),rep(rownames(data),each=4),data4)
data5[,2] <- gsub("\\\\..*$","",data5[,2])
data5[,2] <- gsub("\\\\-.*$","",data5[,2])

# Get a list of unique popnames without dots or dashes
popnames <- unique(data5[,2])

# Change population names to a number
for(i in 1:length(popnames)) {
  data5[data5[,2]==popnames[i],2] <- i
}

# Save file as space separated text file with no row or col names. No quotes.
write.table(data5, "structure_0.05cutoff_infile.csv", row.names=F, col.names=F, sep=" ", quote=F)
