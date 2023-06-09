
### Computational script-------------------------------------------------------------------------------------------------------------------------------------

# This script aims to harvest dataset metadata (... - 2020) related to one or more of the Flemish universities (Katholieke Universiteit Leuven, Universiteit Gent, Universiteit Hasselt, Universiteit Antwerpen), according to predetermined queries. The harvest took place the 26th of December 2020. The output of this harvesting process can be found in the file "Data_DataCite_harvest_Flemish_universities", available in txt- and csv-formats. 


## Load libraries--------------------------------------------------------------------------------------------------------------------------------------------

.libPaths("C:/R/library") # Modify the path to where the R libraries are available. If the packages are not yet installed, please install them first.

library("rdatacite")
library("dplyr")
library("stringr")
library("reshape")


## DataCite - Flemish universities---------------------------------------------------------------------------------------------------------------------------


# Request data from DataCite API
queries <- c("Vrije AND Universiteit AND Brussel", "https://ror.org/006e5kg04", "Universiteit AND Gent", "UGent", "https://ror.org/00cv9y106", "KU AND Leuven", "Katholieke AND Universiteit AND Leuven", "KULeuven", "Catholic AND University AND Leuven", "University AND Leuven", "https://ror.org/05f950310", "Universiteit AND Antwerpen", "UAntwerpen", "UAntwerp", "University AND Antwerp", "https://ror.org/008x57b05", "Universiteit AND Hasselt", "UHasselt", "https://ror.org/04nbhqj75") # list of search items: different entities referring to the five Flemish universities

results_list <- list()
for (i in 1:length(queries)) { # loop over the search items, to obtain the data
  
  response<-dc_dois(query = queries[i], resource_type_id = "dataset", limit = 500)
  data<-response[["data"]]

  results_list[[i]]<-data
}

response_Ghent_University<-dc_dois(query = "Ghent AND University AND NOT (Global AND Biodiversity AND Information AND Facility)", resource_type_id = "dataset", limit = 500, random = TRUE) # separate search to exclude hits from the Global Biodiversity Information Facility 
results_list[[(length(results_list)+1)]]<-response_Ghent_University[["data"]]

response_Hasselt_University<-dc_dois(query = "University AND Hasselt AND NOT (Global AND Biodiversity AND Information AND Facility)", resource_type_id = "dataset", limit = 500, random = TRUE) # separate search to exclude hits from the Global Biodiversity Information Facility
results_list[[(length(results_list)+1)]]<-response_Hasselt_University[["data"]]

queries_2 <- c("Vrije AND Universiteit AND Brussel", "https://ror.org/006e5kg04", "Universiteit AND Gent", "UGent", "https://ror.org/00cv9y106", "KU AND Leuven", "Katholieke AND Universiteit AND Leuven", "KULeuven", "Catholic AND University AND Leuven", "University AND Leuven", "https://ror.org/05f950310", "Universiteit AND Antwerpen", "UAntwerpen", "UAntwerp", "University AND Antwerp", "https://ror.org/008x57b05", "Universiteit AND Hasselt", "UHasselt", "https://ror.org/04nbhqj75", "Ghent AND University AND NOT (Global AND Biodiversity AND Information AND Facility)", "University AND Hasselt AND NOT (Global AND Biodiversity AND Information AND Facility)") # overview of all the queries

                                                        
# Transform selection of requested data into dataframe
final_results_list<-list()

for (i in 1:length(results_list)) {

  if (length(results_list[[i]]) != 0) {
  
  response_data <- results_list[[i]]
  
  response_data_URL<-select(response_data$attributes, url) # url
  response_data_DOI<-select(response_data$attributes, doi) # doi
  response_query<-rep(queries_2[i], length(response_data_DOI[,1])) # query
  
  response_descriptions<-select(response_data$attributes, descriptions)$descriptions # descriptions
  response_descriptions_2<-unlist(lapply(response_descriptions, function(x) ifelse(length(x) == 0, "no_description", x$description)))
  
  response_data_publisher<-select(response_data$attributes, publisher) # publisher
  response_data_client<-select(response_data$relationships, client)$client$data$id # client id
  response_data_publicationYear<-select(response_data$attributes, publicationYear) # publicationYear
  
  response_data_relatedIdentifiers<-select(response_data$attributes, relatedIdentifiers) # identifiers of related research outputs 
  response_data_relatedIdentifiers_2<-unlist(lapply(response_data_relatedIdentifiers[[1]], function(x) ifelse(length(x) != 0, paste(x, collapse=" ; "), "no_identifier")))
  
  response_data_creators<-select(response_data$attributes, creators) # data creators
  response_data_creators_names<-unlist(lapply(response_data_creators[[1]], function(x) paste(x$name, collapse=" / "))) # name data creators
  
  response_data_creators_nameIdentifiers<-lapply(response_data_creators[[1]], function(x) ifelse(length(x$nameIdentifiers) != 0, paste(unlist(x), collapse=" ; "), "no_identifier")) # identifiers of data creators
  response_data_creators_nameIdentifiers_2<-unlist(lapply(response_data_creators_nameIdentifiers, function(x) ifelse(x == "list()" | x == "no_identifier", "no_identifier", paste(unlist(str_extract_all(unlist(x[[1]]), "https://orcid.org/...................")), collapse=" ; "))))
  response_data_creators_nameIdentifiers_3<-ifelse(response_data_creators_nameIdentifiers_2 == "", "no_identifier", response_data_creators_nameIdentifiers_2)
  
  response_data_creators_affiliation<-unlist(lapply(response_data_creators[[1]], function(x) paste(x$affiliation, collapse=" / "))) # affiliation of data creators
  response_data_creators_affiliation_2<-unlist(lapply(response_data_creators_affiliation, function(x) ifelse(grepl("list()", x, fixed = FALSE), "no_affiliation", x)))
  
  response_results<-data.frame(response_query, response_data_DOI, response_data_URL, response_descriptions_2, response_data_publisher, response_data_client, response_data_publicationYear, response_data_relatedIdentifiers_2, response_data_creators_names, response_data_creators_nameIdentifiers_3, response_data_creators_affiliation_2) # resulting data frame with all information
  
  final_results_list[[i]]<-response_results
}}

final_results_list = final_results_list[-which(sapply(final_results_list, is.null))] # filter out the null results

headers_for_columns<-c("Query", "DOI", "URL", "Description", "Publisher", "Client_ID", "PublicationYear", "Identifiers_related_outputs", "Names_data_creators", "Identifiers_data_creators", "Affiliations_data_creators")

total_results_2<-lapply(final_results_list, setNames, nm = headers_for_columns )
total_results_3<-merge_all(total_results_2)

total_results_3$Query <- str_replace_all(total_results_3$Query, '\\n', '') # eliminate newlines in the data
total_results_3$DOI <- str_replace_all(total_results_3$DOI, '\\n', '')
total_results_3$URL <- str_replace_all(total_results_3$URL, '\\n', '')
total_results_3$Description <- str_replace_all(total_results_3$Description, '\\\\$', '')
total_results_3$Description <- str_replace_all(total_results_3$Description, '\\n\\n', '')
total_results_3$Description <- str_replace_all(total_results_3$Description, '\\r\\n', '')
total_results_3$Description <- str_replace_all(total_results_3$Description, '\\n', '')
total_results_3$Publisher <- str_replace_all(total_results_3$Publisher, '\\n', '')
total_results_3$Client_ID <- str_replace_all(total_results_3$Client_ID, '\\n', '')
total_results_3$PublicationYear <- str_replace_all(total_results_3$PublicationYear, '\\n', '')
total_results_3$Identifiers_related_outputs <- str_replace_all(total_results_3$Identifiers_related_outputs, '\\n', '')
total_results_3$Names_data_creators <- str_replace_all(total_results_3$Names_data_creators, '\\n', '')
total_results_3$Identifiers_data_creators <- str_replace_all(total_results_3$Identifiers_data_creators, '\\n', '')
total_results_3$Affiliations_data_creators <- str_replace_all(total_results_3$Affiliations_data_creators, '\\n', '')


# Detect which metadata fields contain affiliation information. Of course, since the cities (Ghent etc.) are used in the queries, there will obviously be a certain number of false positives. Therefore, this procedure also had to be manually checked/validated afterwards.
search_affiliation<-c("Brussel", "brussel", "Gent", "gent", "Ghent", "ghent", "Antwerpen", "antwerpen", "Antwerp", "antwerp", "Hasselt", "hasselt", "Leuven", "leuven") # search terms to spot affiliation information

affiliation_hits_total<-vector()
for (i in 1:length(search_affiliation)) {
affiliation_hits<-sapply(colnames(total_results_3)[2:ncol(total_results_3)], function(x) grep(search_affiliation[i], total_results_3[,x]))
affiliation_hits_vector<-unlist(affiliation_hits)
affiliation_hits_total<-append(affiliation_hits_total, affiliation_hits_vector)
}
names(affiliation_hits_total)<-str_replace_all(names(affiliation_hits_total), "[:digit:]", "")

affiliation_information_fields<-vector()
for (i in 1:nrow(total_results_3)) {
  all_matches<-affiliation_hits_total[which(affiliation_hits_total==i)]
  if (length(all_matches) != 0) {
  unified<-paste(sort(unique(names(all_matches))), collapse="+")}
  else {unified <- "no_field"}
  affiliation_information_fields<-append(affiliation_information_fields, unified)
}
total_results_3$affiliation_information_fields <- affiliation_information_fields


# Remove duplicate DOIs
DOI_duplicates<-as.character(filter(data.frame(table(total_results_3$DOI)), Freq>1)$Var1)

rows_to_be_deleted<-vector()
for (i in 1:length(DOI_duplicates)) {
  position_duplicates<-which(total_results_3$DOI == DOI_duplicates[i])
  number_to_be_deleted<-length(position_duplicates)-1
  rows_to_be_deleted<-append(rows_to_be_deleted, position_duplicates[1:number_to_be_deleted]) 
}
total_results_4<-total_results_3[-rows_to_be_deleted,]


# Indicate DOIs that are versions or parts of main dataset DOI
potential_doubles<-filter(as.data.frame(table(total_results_4$Client_ID, total_results_4$PublicationYear, total_results_4$Names_data_creators)), Freq > 1)

positions_doubles<-list()
for (i in 1:length(potential_doubles[,1])) {
  potential_doubles_2<-potential_doubles[i,]
  positions_doubles[[i]]<-which(total_results_4$Client_ID == potential_doubles_2$Var1 & total_results_4$PublicationYear == potential_doubles_2$Var2 & total_results_4$Names_data_creators == potential_doubles_2$Var3)
  }

non_primary<-vector() # indices of rows with data points that do not refer to the main dataset DOIs
for (i in 1:length(positions_doubles)) {
  temporary_data_frame<-total_results_4[positions_doubles[[i]],]
  total_rows_temporary<-nrow(temporary_data_frame)
  
  harvard_file <- str_detect(temporary_data_frame$URL, "https://dataverse.harvard.edu/file") # Harvard dataverse URLs contain "file" to indicate separate files
  version_in_DOI <- str_ends(temporary_data_frame$DOI, "\\\\.v[:digit:]") # detect for example "10.6084/m9.figshare.10003385.v1", with ending .v1 (= first version)
  version_of <- str_detect(temporary_data_frame$Identifiers_related_outputs, "VersionOf") # detect related identifiers with relation  type "VersionOf"
  part_of <- str_detect(temporary_data_frame$Identifiers_related_outputs, "PartOf") # detect related identifiers with relation  type "PartOf"
  
  if (TRUE %in% harvard_file & length(which(harvard_file == TRUE)) < total_rows_temporary){ # the second condition makes sure that not all DOIs referring to the same data package are marked as "non primary", if all of them are separate files 
      non_primary<-append(non_primary, positions_doubles[[i]][which(harvard_file == TRUE)])
      
      } else if (length(which(harvard_file == TRUE)) == total_rows_temporary){    # if all DOIs refer to separate files: select one as "primary"
      non_primary<-append(non_primary, positions_doubles[[i]][2:total_rows_temporary])
      
      } else if (TRUE %in% version_in_DOI & length(which(version_in_DOI == TRUE)) < total_rows_temporary){
      non_primary<-append(non_primary, positions_doubles[[i]][which(version_in_DOI == TRUE)])
      
      } else if (TRUE %in% version_of & length(which(version_of == TRUE)) < total_rows_temporary){
      non_primary<-append(non_primary, positions_doubles[[i]][which(version_of == TRUE)])
      
      } else if (TRUE %in% part_of & length(which(part_of == TRUE)) < total_rows_temporary){
        non_primary<-append(non_primary, positions_doubles[[i]][which(part_of == TRUE)])
    
        # detect overlapping descriptions of substantial size
      } else if (min(nchar(temporary_data_frame$Description)) > 40 & length(unique(temporary_data_frame$Description)) == 1){ 
        non_primary<-append(non_primary, positions_doubles[[i]][2:total_rows_temporary])
      }}

Versioning<-rep("primary", nrow(total_results_4))
for (i in 1:length(Versioning)) {
  if (i %in% unique(non_primary)){
    Versioning[i]<-"non_primary"
  }}
total_results_4$Versioning <- Versioning

total_results_5<-filter(total_results_4, Versioning == "primary")
potential_doubles_2<-filter(as.data.frame(table(total_results_5$Client_ID, total_results_5$PublicationYear, total_results_5$Names_data_creators)), Freq > 1) # After the automatic analysis, there are still some remaining cases, listed in potential_doubles_2, that have to be checked manually (cf. the variable Versioning_manual in file "Data_DataCite_harvest_Flemish_universities").

write.table(total_results_4, file="Data_DataCite_harvest_Flemish_universities.txt", sep = '/t', row.names = FALSE, col.names = TRUE, quote=FALSE, eol="\\n") # create output file

##-----------------------------------------------------------------------------------------------------------------------------------------------------------

rm(list=ls())

