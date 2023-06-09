

### Computational script-------------------------------------------------------------------------------------------------------------------------------------

# This script aims to harvest links to related outputs from the Scholexplorer API, on the basis of the dataset DOIs available in the file "Data_DataCite_harvest_Flemish_universities.txt". The harvest took place the 26th of December 2020. The output of this harvesting process can be found in the file "Data_Scholexplorer_harvest_Flemish_universities", available in txt- and csv-formats. 


## Load libraries--------------------------------------------------------------------------------------------------------------------------------------------

.libPaths("C:/R/library") # Modify the path to where the R libraries are available. If the packages are not yet installed, please install them first.

library(httr)
library(stringr)
library(jsonlite)
library(dplyr)


## Load data------------------------------------------------------------------------------------------------------------------------------------------------

data<-read.table(file=file.choose(), header=TRUE, sep="\\t", comment.char="", fill=TRUE, na.strings="character", quote="", row.names=NULL, stringsAsFactors = FALSE, strip.white = TRUE, encoding = "utf-8", blank.lines.skip = TRUE) # input the file "Data_DataCite_harvest_Flemish_universities.txt"

data_2<-filter(data[,1:15], Inclusion == "yes" & Versioning_manual == "primary")
data<-data_2

DOI_data<-as.data.frame(data[,3])


## Request data from Scholexplorer API-----------------------------------------------------------------------------------------------------------------------

scholix_api <- function(query) {
  query_2<-str_replace_all(query, "/", "%2F")
  url <- paste("http://api.scholexplorer.openaire.eu/v2/Links?sourcePid=", query_2, sep = "")
  GET(url)
  }

scholix_resp<-list()

for (i in 1:length(DOI_data[,1])) {
  
  temporary<-DOI_data[i,]
  scholix_resp[[i]]<- scholix_api(temporary)
  
}

scholix_parsed<-list()

for (i in 1:length(scholix_resp)) {
  
  if (scholix_resp[[i]]$status_code == 200) {
  
  temporary_2 <- scholix_resp[[i]]
  scholix_parsed[[i]] <- content(temporary_2, "parsed")
    
  } else {
    
    scholix_parsed <- "Wrong status code"
    }}

scholix_results<-list()

for (i in 1:length(scholix_parsed)) {

  if (length(scholix_parsed[[i]]$result) == 0) {
    
    temporary_3 <- scholix_parsed[[i]]
    scholix_results[[i]] <- "No scholix links"
    
  } else {
    
    scholix_results[[i]]<-list()
    number_of_results <- length(scholix_parsed[[i]]$result)
    
    for (j in 1:(length(number_of_results))) {
    
      DOI_of_source_dataset<-DOI_data[i,1]
      
      RelationshipType<-scholix_parsed[[i]]$result[[j]]$RelationshipType$Name
      
      RelationshipSubtype<-scholix_parsed[[i]]$result[[j]]$RelationshipType$SubType
      
      Type_of_research_output<-scholix_parsed[[i]]$result[[j]]$target$Type
      
      Persistent_identifier<-scholix_parsed[[i]]$result[[j]]$target$Identifier[[1]]$IDURL
      
      Results_total<-data.frame(DOI_of_source_dataset, RelationshipType, RelationshipSubtype, Type_of_research_output, Persistent_identifier)
    
      scholix_results[[i]][[j]]<-Results_total
      
    }
    
    next
    
    }}
  


scholix_number_links<-vector()

for (i in 1:length(scholix_results)) {
  
  scholix_number_links<-append(scholix_number_links, ifelse(scholix_results[[i]] == "No scholix links", 0, length(scholix_results[[i]])))
}

scholix_number_links[scholix_number_links > 1]

scholix_link_to<-vector()

for (i in 1:length(scholix_results)) {
  
  scholix_link_to<-append(scholix_link_to, ifelse(scholix_results[[i]] == "No scholix links", "No_scholix_links", as.character(scholix_results[[i]][[1]]$Type_of_research_output)))
  }


scholix_DOI<-vector()
for (i in 1:length(scholix_results)) {
  
  scholix_DOI<-append(scholix_DOI, ifelse(scholix_results[[i]] == "No scholix links", "No_DOI", as.character(scholix_results[[i]][[1]]$Persistent_identifier)))
  
}


scholix_table_end<-data.frame(scholix_link_to, scholix_DOI)

write.table(scholix_table_end, file="Scholix_info.txt", sep = '\\t', row.names = FALSE, col.names = TRUE, quote=FALSE, eol="\\n") # create output file

##-----------------------------------------------------------------------------------------------------------------------------------------------------------

rm(list=ls())





