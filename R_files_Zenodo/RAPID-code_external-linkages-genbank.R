##############
###Title: RAPID Code to Find GenBank Sequences
###Author: Katelin D. Pearson
###Date: February 4, 2021
#############

library(splitstackshape)
library(tidyverse)
library(DataCombine)

##Part 1: Extract Catalog Numbers from GenBank Results

#Read in dataset (a gff3 file converted to csv from NCBI Nuccore)
sample <- read_csv("PATH",
  col_types = cols(
  text = col_character()
))

#Make a dataframe to hold each line that contains a catalog number
positives <- matrix(NA,nrow=100000)
positives <- as.data.frame(positives)
colnames(positives) <- c("text")

#Go through the dataset and extract only lines with catalog numbers
#In GenBank, these are values in the "specimen-voucher" field
k = 1
#ptm <- proc.time() #for if you want to time the loop
for(i in 1:dim(sample)[1]){
  if(sample$text[i]==""){
    next
  } else {
    if(is.na(sample$text[i])){
      next
    } else {
      if(grepl("specimen-voucher=",sample$text[i])){
        positives$text[k] <- as.character(sample$text[i])
        k <- k + 1
      }
    }
  }
}
#ptm - proc.time() #for if you want to time the loop

#Remove any blank rows
positives <- positives %>% filter_all(any_vars(!is.na(.)))

#Split the results into two columns
positives <- cSplit(positives, "text", sep = "+")

#Because the first column always has the GenBankID at the beginning followed by a tab
#We can extract this value using cSplit and put it into its own column
positives <- cSplit(positives, "text_1", sep = "/t")
colnames(positives)[2] <- "GenBankID"

#Make a new catalogNumber field
positives$catalogNumber <- NA

#Extract the catalogNumbers from each specimen record
for(i in 1:dim(positives)[1]){
  cn <- sub(".*specimen-voucher=","", positives$text_2[i])
  positives$catalogNumber[i] <- cn
}

#Remove garbage on the other side of the catalog number
positives <- cSplit(positives,"catalogNumber", sep = ";")

#Make a cleaned dataframe of only catalog number and GenBankID
cleaned <- cbind(catalogNumber = as.character(positives$catalogNumber_1), GenBankID = as.character(positives$GenBankID))
cleaned <- as.data.frame(cleaned)

#standardize ROM catalog numbers that say ROM MAM or ROM: into just ROM
tofrom <- matrix(ncol = 2, nrow = 1)
tofrom <- as.data.frame(tofrom)
colnames(tofrom) <- c("from","to")
tofrom$from[1] <- "ROM MAM"
tofrom$to[1] <- "ROM"
cleaned <- FindReplace(cleaned, "catalogNumber", tofrom, from = "from", to = "to", exact = FALSE)
tofrom$from[1] <- "ROM:"
tofrom$to[1] <- "ROM "
cleaned <- FindReplace(cleaned, "catalogNumber", tofrom, from = "from", to = "to", exact = FALSE)

#standardize Western Australian Museum records
for(i in 1:dim(cleaned)[1]){
  if(grepl("Western Australian", cleaned$catalogNumber[i])){
    this1 <- sub(" Western Australian Museum.*", "", as.character(cleaned$catalogNumber[i]))
    this2 <- paste("WAM ",this1)
    cleaned$catalogNumber[i] <- this2
  } else {
    next
  }
}

write.csv(cleaned, "PATH")

##Part 2: Unify Catalog Number in Specimen Dataset
#The GenBank specimen voucher numbers are generally the institution code followed by the number
#Often the catalog number in the specimen dataset is only the catalog number
#We need to make a catalog number field that matches the GenBank format for matching purposes

#load GenBank dataset if need be 
#cleaned <- read.csv("PATH")

specs <- read.csv("PATH")

specs$unifiedCatNum <- paste(as.character(specs$institutionCode_gbifP),as.character(specs$catalogNumber_gbifP),sep = " ")

#look at what else needs to be cleaned or standardized:

#look at values in the catalogNumber field for the GenBank data
# table(substr(as.character(cleaned$catalogNumber), 1, 4))
#look at values in the unifiedCatNum field for the specimen data
# table(substr(specs$institutionCode_gbifP, 1, 4))

##Part 3: Compare Catalog Number in Specimen Dataset to Catalog Numbers for GenBank sequences

#only uncomment the first time you are making this column
# specs$associatedSequences_rapid <- NA

for(i in 1:dim(specs)[1]){
  matches <- which(cleaned$catalogNumber %in% specs$unifiedCatNum[i])
  if(length(matches)==0){
    next
  } else {
    specs$associatedSequences_rapid[i] <- paste(paste("http://www.ncbi.nlm.nih.gov/nuccore/",as.character(cleaned$GenBankID[matches]),sep=""), collapse = "|")
  }
}

View(specs$associatedSequences_rapid)
t <- table(specs$associatedSequences_rapid)
View(t)

write.csv(specs, "PATH")