#####
#title: ExtractDuplicateGeoreferences, version 2.0
#author: Katelin (Katie) Pearson
#contact: katelin.d.pearson24@gmail.com
#date: March 25, 2021
#details: This script was developed for the California Phenology Network (capturingcaliforniasflowers.org) and made possible by National Science Foundation Awards 1802312 and 2001500. Any opinions, findings, and conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the National Science Foundation.
####

#load packages
library(dplyr)
library(stringr)
library(gtools)
library(DataCombine)

#load the table of specimens (occids) that have duplicates
#in most Symbiota instances, this is the omoccurduplicatelink table
duplinks <- read.csv("PATH")

#load the table of georeference data from the database
#we used a modified version of omoccurrences that only includes the columns occid,
#decimalLatitude, decimalLongitude, geodeticDatum, coordinateUncertaintyInMeters, footprintWKT,
#coordinatePrecision, georeferencedBy, georeferenceSources, and georeferenceRemarks
#for specimens for which decimalLatitude IS NOT NULL
others <- read.csv("PATH")

#load the table of data from the collection for which you are importing data
#make sure that the first column contains the occurrence ID and is called "id"
mycoll <- read.csv("PATH")

#make a new data frame that will hold all the new georeferences
#note that the maximum is 50,000 records here; increase size as needed
newMycoll <- matrix(ncol=length(colnames(others)),nrow=50000)
colnames(newMycoll) <- colnames(others)
newMycoll <- as.data.frame(newMycoll)
f <- 1

#look at each record in the mycoll dataframe
for(i in 1:dim(mycoll)[1]){
  print(i)
  #if desired, uncomment the line below to track progress (it will slow down the code, though)
  #print(paste((i/dim(mycoll)[1])*100,"% complete"))
  
  #uncomment lines 45-47 and line 103 if your input data (from your collection) contains georeferenced specimens
  #if the record in mycoll already has georeference data, go to the next record
  #if(!is.na(mycoll$decimalLatitude[i])){
  #  next
  #} else {
  
  #does the mycoll record have a duplicate?
  m <- match(mycoll$id[i],duplinks$occid)
  
  #if no duplicate found, go to the next record
  if(is.na(m)){
    next
  }  
  else {
    this <- mycoll$id[i]
    #find the id number of the duplicate cluster to which this record belongs
    dupid <- duplinks$duplicateid[duplinks$occid==this]
    #make a temporary dataframe of all the records that belong to this cluster
    dupes <- subset(duplinks,duplicateid==dupid)
    #from the temporary dataframe, remove the record that we are searching against (the original record)
    dupes <- subset(dupes, occid!=mycoll$id[i])
    #if there are no other duplicate records (because of a mistake), go to the next record in mycoll
    if(dim(dupes)[1]<1){
      next
    } else {
      #for every duplicate in the temporary dataframe
      for(j in 1:dim(dupes)[1]){
        #find the occid for that duplicate
        thisdup <- dupes$occid[j]
        #if the duplicate does not exist within the full dataset (i.e., no georeferenced duplicate exists), go to next
        if(thisdup %in% others$occid){
          
          #uncomment lines 77-79 and line 98 if your "others" dataset contains specimens that are not georeferenced
          #if the duplicate does not have georeference data, go to the next duplicate in the temporary dataframe
          #if(is.na(others$decimalLatitude[others$occid==thisdup])){
          #  next
          #} else {
          
          #if the duplicate DOES have georeference data, add all the data to the "newMycoll" dataframe
          newMycoll$occid[f]=mycoll$id[i]
          newMycoll$decimalLatitude[f] <- as.character(others$decimalLatitude[others$occid==thisdup])
          newMycoll$decimalLongitude[f] <- as.character(others$decimalLongitude[others$occid==thisdup])
          newMycoll$geodeticDatum[f] <- as.character(others$geodeticDatum[others$occid==thisdup])
          newMycoll$coordinateUncertaintyInMeters[f] <- as.character(others$coordinateUncertaintyInMeters[others$occid==thisdup])
          newMycoll$footprintWKT[f] <- as.character(others$footprintWKT[others$occid==thisdup])
          newMycoll$coordinatePrecision[f]=as.character(others$coordinatePrecision[others$occid==thisdup])
          newMycoll$georeferencedBy[f] <- as.character(others$georeferencedBy[others$occid==thisdup])
          newMycoll$georeferenceSources[f] <- as.character(others$georeferenceSources[others$occid==thisdup])
          newMycoll$georeferenceRemarks[f] <- paste("copied from duplicate collid=",as.character(others$collid[others$occid==thisdup])," ",as.character(others$catalogNumber[others$occid==thisdup]),"; ", as.character(others$georeferenceRemarks[others$occid==thisdup]),sep="")
          newMycoll$catalogNumber[f] <- as.character(mycoll$catalogNumber[i])
          newMycoll$otherCatalogNumbers[f] <- as.character(mycoll$otherCatalogNumbers[i])
          newMycoll$recordedBy[f] <- as.character(mycoll$recordedBy[i])
          newMycoll$recordNumber[f] <- as.character(mycoll$recordNumber[i])
          f <- f+1
        }
        #}
      }
    }
  }
}
#}

#Export result to new file
write.csv(newMycoll,"PATH")


#in the code above, the georeferenceRemarks fields is populated with the statement "copied from duplicate"
#followed by the collid of the duplicate specimen, the catalog number of the duplicate specimen, and the georeference remarks from that specimen
#the code below replaces the collid values with the acronyms of their respective institutions

#to use this line of code, you will need a file that contains a column of collection ID numbers in the format "collid=123" (for example)
#and a column of the collection acronyms corresponding to each collid.
#You can create this file from the omcollections table in your Symbiota portal.

#NOTE THAT THIS FILE NEEDS TO BE SORTED IN DESCENDING ORDER (in order of collid) TO AVOID ERRONEOUS REPLACEMENTS
replaces <- read.csv("PATH")
tempMycoll <- FindReplace(newMycoll, "georeferenceRemarks", replaces, from = "CollidEquals", to = "Acronym", exact = FALSE)
write.csv(tempMycoll,"PATH")