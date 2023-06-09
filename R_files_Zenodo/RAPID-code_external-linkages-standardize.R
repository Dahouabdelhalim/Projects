##############
###Title: RAPID Code to Standardize External Linkages
###Author: Katelin D. Pearson
###Date: February 4, 2021
#############

library(splitstackshape)

#load the dataset of specimen records
#the default associatedSequences field in our dataset is named associatedSequences_gbifP
bats <- read.csv("PATH")

#make a subset of the dataset for which the associatedSequences field is not empty
#(to get a simple count)
assocSeq <- subset(bats, bats$associatedSequences_gbifP!="")

#make a new column for the re-formatted sequences
#DO NOT RUN THIS if you have already run the GenBank-Specimen Linkage code on your dataset
#bats$associatedSequences_rapid <- NA

bats$associatedSequences_rapid <- as.character(bats$associatedSequences_rapid)

#populate the new (or previously existing) column
#we are using the NCBI URL as recommended by the Darwin Core
for(i in 1:dim(bats)[1]){
  if(bats$associatedSequences_gbifP[i]==""){
    next
  } else {
    if(is.na(bats$associatedSequences_gbifP[i])){
      next
    } else {
      #transfer over the sequences that are already correctly formatted with the nuccure URL
      if(grepl("nuccore",as.character(bats$associatedSequences_gbifP[i]))){
        bats$associatedSequences_rapid[i] <- as.character(bats$associatedSequences_gbifP[i])
      } else {
        #reformat the records that are only an NCBI number
        if(nchar(as.character(bats$associatedSequences_gbifP[i]))<10){
          bats$associatedSequences_rapid[i] <- as.character(paste("http://www.ncbi.nlm.nih.gov/nuccore/",as.character(bats$associatedSequences_gbifP[i]),sep=""))
        }
      }
    }
  }
}

#make sure it works
t <- table(bats$associatedSequences_rapid)
t <- as.data.frame(t)

#now to deal with columns that have two sequences
#preserve the original column
bats$associatedSequences_gbifP_orig <- bats$associatedSequences_gbifP

#split the column into two columns, assuming the sequences are separated by a |
bats <- cSplit(bats,"associatedSequences_gbifP",sep="|",type.convert = FALSE)

#isolate each incorrect URL and apply the correct URL, then join them together
#into the associatedSequences_rapid column
for(i in 1:dim(bats)[1]){
  if(!is.na(bats$associatedSequences_gbifP_2[i])){
    first <- sub(".*term=", "", as.character(bats$associatedSequences_gbifP_1[i]))
    second <- sub(".*term=", "", as.character(bats$associatedSequences_gbifP_2[i]))
    bats$associatedSequences_rapid[i] <- paste(paste("http://www.ncbi.nlm.nih.gov/nuccore/",first, sep = ""), paste("http://www.ncbi.nlm.nih.gov/nuccore/",second,sep = ""),sep = "|")
    print(i)
  }
}

#remove erroneous associatedSequences
for(i in 1:dim(bats)[1]){
  if(is.na(bats$associatedSequences_rapid[i])){
    next
  } else {
    if(nchar(as.character(bats$associatedSequences_rapid[i]))<10){
      bats$associatedSequences_rapid[i] <- NA
    }
  }
}