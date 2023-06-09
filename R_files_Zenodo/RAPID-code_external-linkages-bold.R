##############
###Title: RAPID Code to Find and Link Barcode of Life External Linkages
###Author: Katelin D. Pearson
###Date: February 12, 2021
#############

library(DataCombine)

#load in all Hipposideridae and Rhinolophidae sequences from BOLD
boldbats <- read.csv("C:/Users/Katie/Documents/FSU_RAPID/bold_hippo_rhino.csv")

#exclude records that were already mined from NCBI
boldbats <- subset(boldbats,boldbats$institutionCode!="Mined from GenBank, NCBI")

#standardize the ROM and MZB catalog numbers
#no other catalog numbers were identified as needing standardization
tofrom <- matrix(ncol = 2, nrow = 2)
tofrom <- as.data.frame(tofrom)
colnames(tofrom) <- c("from","to")
tofrom$from[1] <- "ROM:MAM:"
tofrom$to[1] <- "ROM "
tofrom$from[2] <- "MZB_"
tofrom$to[2] <- "MZB "
boldbats <- FindReplace(boldbats, "catalogNumber", tofrom, from = "from", to = "to", exact = FALSE)
boldbats$catalogNumber <- as.character(boldbats$catalogNumber)

#need to add ROM when belongs to Royal Ontario Museum (some were missing acronym in catalog number)
for(i in 1:dim(boldbats)[1]){
  if(grepl("Royal Ontario",boldbats$institutionCode[i])){
    if(is.na(boldbats$catalogNumber[i])){
      next
    } else {
      if(grepl("ROM",boldbats$catalogNumber[i])==FALSE){
        boldbats$catalogNumber[i] <- paste("ROM",boldbats$catalogNumber[i])
      }
    }
  } else {
    next
  }
}

#load the previous specimen data
specs <- read.csv("C:/Users/Katie/Documents/FSU_RAPID/FinalParsedFiles/final_standardized_assocSeq.csv")

#add a flag for specimens that already have sequences from NCBI
#these are expected to be dupliates, not new sequences
specs$inBOLD <- NA
specs$associatedSequences_rapid <- as.character(specs$associatedSequences_rapid)

#find BOLD sequences corresponding to the specimens in our dataset
for(i in 1:dim(specs)[1]){
  matches <- which(boldbats$catalogNumber %in% specs$unifiedCatNum[i])
  if(length(matches)==0){
    next
  } else {
    if(!is.na(specs$associatedSequences_rapid[i])){
     specs$inBOLD[i] <- TRUE
    } else {
      specs$associatedSequences_rapid[i] <- paste(paste("http://www.barcodinglife.org/index.php/Public_RecordView?processid=",as.character(boldbats$occurrenceID[matches]),sep=""), collapse = "|")
    }
  }
}

#how many specimens were from BOLD?
boldIDs <- subset(specs,grepl("barcoding",specs$associatedSequences_rapid))

#how many specimens were both in NCBI and BOLD?
boldAndNCBI <- subset(specs, specs$inBOLD==TRUE)

#output data file
write.csv(specs,"C:/Users/Katie/Documents/FSU_RAPID/associatedSequences_after_bold.csv")
