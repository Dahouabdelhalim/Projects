
### Computational script-------------------------------------------------------------------------------------------------------------------------------------

# This script aims to harvest dataset DOIs from DataCite in a random fashion. Results are restricted to repositories currently used by researchers affiliated to Flemish universities. For each repository, a sample of 1000 metadata records is collected (if possible), limited to the publication years 2019 and 2020. Next, the raw data are trimmed down to samples of 450 records per repository (when possible). 


## Load libraries--------------------------------------------------------------------------------------------------------------------------------------------

.libPaths("C:/R/library") # Modify the path to where the R libraries are available. If the packages are not yet installed, please install them first.

library("rdatacite")
library("dplyr")
library("tidyr")
library("stringr")
library("reshape")


## DataCite - repositories used by researchers affiliated to Flemish universities----------------------------------------------------------------------------

data<-read.table(file=file.choose(), header=TRUE, sep="\\t", comment.char="", fill=TRUE, na.strings="character", quote="", row.names=NULL, stringsAsFactors = FALSE, strip.white = TRUE, encoding = "utf-8", blank.lines.skip = TRUE) # cf. the file "Data_DataCite_harvest_Flemish_universities.txt"

data_2<-filter(data[,1:15], Inclusion == "yes" & Versioning_manual == "primary")
data<-data_2

client_id_repo<-unique(data$Client_ID) # repositories used by Flemish researchers
results_repo_list <- list()

for (j in 1:length(client_id_repo)) {

response_repo<-dc_dois(query = "publicationYear:2019 OR publicationYear:2020", client_id = client_id_repo[j], resource_type_id = "dataset", sample_size = "1000", random = TRUE, limit = 1000) # sample of 1000 hits, for publication years 2019 and 2020
data_repo<-response_repo[["data"]]

results_repo_list[[j]]<-data_repo
}


# Transform selection of requested data into dataframe       
final_results_repo_list<-list()

for (i in 1:length(results_repo_list)) {
  
  if (length(results_repo_list[[i]]) != 0) {
    
    response_data_repo <- results_repo_list[[i]]
    
    response_data_repo_DOI<-select(response_data_repo$attributes, doi) #doi
    
    response_data_repo_publisher<-select(response_data_repo$attributes, publisher) #publisher
    response_data_repo_client<-select(response_data_repo$relationships, client)$client$data$id # client id
    response_data_repo_publicationYear<-select(response_data_repo$attributes, publicationYear) #publicationYear
    
    response_data_repo_relatedIdentifiers<-select(response_data_repo$attributes, relatedIdentifiers) #identifiers of related research outputs 
    response_data_repo_relatedIdentifiers_2<-unlist(lapply(response_data_repo_relatedIdentifiers[[1]], function(x) ifelse(length(x) != 0, paste(x, collapse=" ; "), "no_identifier")))
    
    response_data_repo_creators<-select(response_data_repo$attributes, creators) #data creators
    response_data_repo_creators_names<-unlist(lapply(response_data_repo_creators[[1]], function(x) paste(x$name, collapse=" / "))) # name data creators
    
    response_data_repo_creators_nameIdentifiers<-lapply(response_data_repo_creators[[1]], function(x) ifelse(length(x$nameIdentifiers) != 0, paste(unlist(x), collapse=" ; "), "no_identifier"))# identifiers of data creators
    response_data_repo_creators_nameIdentifiers_2<-unlist(lapply(response_data_repo_creators_nameIdentifiers, function(x) ifelse(x == "list()" | x == "no_identifier", "no_identifier", paste(unlist(str_extract_all(unlist(x[[1]]), "https://orcid.org/...................")), collapse=" ; "))))
    response_data_repo_creators_nameIdentifiers_3<-ifelse(response_data_repo_creators_nameIdentifiers_2 == "", "no_identifier", response_data_repo_creators_nameIdentifiers_2)
    
    response_data_repo_creators_affiliation<-unlist(lapply(response_data_repo_creators[[1]], function(x) paste(x$affiliation, collapse=" / "))) # affiliation of data creators
    response_data_repo_creators_affiliation_2<-unlist(lapply(response_data_repo_creators_affiliation, function(x) ifelse(grepl("list()", x, fixed = FALSE), "no_affiliation", x)))
    
    response_results_repo<-data.frame(response_data_repo_DOI, response_data_repo_publisher, response_data_repo_client, response_data_repo_publicationYear, response_data_repo_relatedIdentifiers_2, response_data_repo_creators_names, response_data_repo_creators_nameIdentifiers_3, response_data_repo_creators_affiliation_2) # resulting data frame with all information
    
    final_results_repo_list[[i]]<-response_results_repo
  }}


headers_for_columns<-c("DOI", "Publisher", "Client_ID", "PublicationYear", "Identifiers_related_outputs", "Names_data_creators", "Identifiers_data_creators", "Affiliations_data_creators")

total_results_repo_2<-lapply(final_results_repo_list, setNames, nm = headers_for_columns )
total_results_repo_3<-merge_all(total_results_repo_2)

total_results_repo_3$DOI <- str_replace_all(total_results_repo_3$DOI, '\\n', '')
total_results_repo_3$Publisher <- str_replace_all(total_results_repo_3$Publisher, '\\n', '')
total_results_repo_3$Client_ID <- str_replace_all(total_results_repo_3$Client_ID, '\\n', '')
total_results_repo_3$PublicationYear <- str_replace_all(total_results_repo_3$PublicationYear, '\\n', '')
total_results_repo_3$Identifiers_related_outputs <- str_replace_all(total_results_repo_3$Identifiers_related_outputs, '\\n', '')
total_results_repo_3$Names_data_creators <- str_replace_all(total_results_repo_3$Names_data_creators, '\\n', '')
total_results_repo_3$Identifiers_data_creators <- str_replace_all(total_results_repo_3$Identifiers_data_creators, '\\n', '')
total_results_repo_3$Affiliations_data_creators <- str_replace_all(total_results_repo_3$Affiliations_data_creators, '\\n', '')


potential_doubles_repo<-filter(as.data.frame(table(total_results_repo_3$Client_ID, total_results_repo_3$PublicationYear, total_results_repo_3$Names_data_creators)), Freq > 1)

positions_doubles_repo<-list()
for (i in 1:length(potential_doubles_repo[,1])) {
  potential_doubles_repo_2<-potential_doubles_repo[i,]
  positions_doubles_repo[[i]]<-which(total_results_repo_3$Client_ID == potential_doubles_repo_2$Var1 & total_results_repo_3$PublicationYear == potential_doubles_repo_2$Var2 & total_results_repo_3$Names_data_creators == potential_doubles_repo_2$Var3)
}

non_primary_repo<-vector()
for (i in 1:length(positions_doubles_repo)) {
  temporary_data_frame_repo<-total_results_repo_3[positions_doubles_repo[[i]],]
  total_rows_temporary_repo<-nrow(temporary_data_frame_repo)
  non_primary_repo<-append(non_primary_repo, positions_doubles_repo[[i]][2:total_rows_temporary_repo])
  }

Versioning_repo<-rep("primary", nrow(total_results_repo_3))
for (i in 1:length(Versioning_repo)) {
  if (i %in% unique(non_primary_repo)){
    Versioning_repo[i]<-"non_primary"
  }}
total_results_repo_3$Versioning <- Versioning_repo



write.table(total_results_repo_3, file="Data_DataCite_harvest_random.txt", sep = '\\t', row.names = FALSE, col.names = TRUE, quote=FALSE, eol="\\n") #create output file "Data_DataCite_harvest_random.txt"

total_results_repo_3<-read.table(file=file.choose(), header=TRUE, sep="\\t", comment.char="", fill=TRUE, na.strings="character", quote="", row.names=NULL, stringsAsFactors = FALSE, strip.white = TRUE, encoding = "utf-8", blank.lines.skip = TRUE) # cf. the file "Data_DataCite_harvest_random.txt"

total_results_repo_4<-filter(total_results_repo_3, Versioning == "primary")
(Distribution_publisher<-as.data.frame(table(total_results_repo_4$Client_ID))[order(as.data.frame(table(total_results_repo_4$Client_ID))[,2], decreasing = TRUE),])


# Recoding of categories
total_results_repo_4$Identifiers_data_creators[total_results_repo_4$Identifiers_data_creators!="no_identifier"]<-"identifier" # everything that is not "no_identifier", is recoded as "identifier"
total_results_repo_4$Affiliations_data_creators[total_results_repo_4$Affiliations_data_creators!="no_affiliation"]<-"affiliation" # everything that is not "no_affiliation", is recoded as "affiliation"

for (i in 1:nrow(total_results_repo_4)) {  # hierarchical attribution
  Supplement <- str_detect(total_results_repo_4$Identifiers_related_outputs[i], "Supplement")
  Reference <- str_detect(total_results_repo_4$Identifiers_related_outputs[i], "Reference")
  Cite <- str_detect(total_results_repo_4$Identifiers_related_outputs[i], "Cite")
  Document <- str_detect(total_results_repo_4$Identifiers_related_outputs[i], "Document")
  
  if (Supplement == TRUE){
    total_results_repo_4$Identifiers_related_outputs_3[i]<-"Supplement" # related identifiers in supplement-relation
  
  } else if (Cite == TRUE | Reference == TRUE){
    total_results_repo_4$Identifiers_related_outputs_3[i]<-"Cite/Reference" # related identifiers in cite/reference-relation
    
  } else if (Document == TRUE){
    total_results_repo_4$Identifiers_related_outputs_3[i]<-"Document" # related identifiers in document-relation
    
  } else if (total_results_repo_4$Identifiers_related_outputs[i] != "no_identifier") { # other types of related identifiers
    total_results_repo_4$Identifiers_related_outputs_3[i]<-"other_identifier"
  
  } else {
    total_results_repo_4$Identifiers_related_outputs_3[i]<-"no_identifier" # no identifier
    
  }}


## Further processing----------------------------------------------------------------------------------------------------------------------------------------

# In this section, the raw harvested data are further trimmed down: samples of 450 hits per repository when possible

data_repo_final<-data.frame(total_results_repo_4$Client_ID, total_results_repo_4$Identifiers_data_creators, total_results_repo_4$Affiliations_data_creators, total_results_repo_4$Identifiers_related_outputs_3)
colnames(data_repo_final)<-c("Client_ID", "Identifiers_data_creators", "Affiliations_data_creators", "Identifiers_related_outputs")

selected_repo<-c("dryad.dryad", "figshare.ars", "bl.mendeley", "cern.zenodo", "gesis.icpsr", "delft.data4tu", "ieee.dataport", "pangaea.repository", "gdcc.harvard-dv") # selected repositories to be examined in this analysis
data_repo_final<-filter(data_repo_final, Client_ID %in% selected_repo)

list_equal_sizes<-list()
for (i in 1:length(selected_repo)){
  temporary <- filter(data_repo_final, Client_ID == selected_repo[i])
  temporary_2 <- sample_n(temporary, 450, replace = FALSE) # sample of 450 hits per repository
  list_equal_sizes[[i]]<-temporary_2
}

headers_for_columns<-c("Client_ID", "Identifiers_data_creators", "Affiliations_data_creators", "Identifiers_related_outputs")

data_repo_final_bis<-lapply(list_equal_sizes, setNames, nm = headers_for_columns )
data_repo_final<-merge_all(data_repo_final_bis)

data_repo_final<- data_repo_final %>% mutate_if(is.factor,as.character)

data_repo_final$Identifiers_data_creators[data_repo_final$Identifiers_data_creators =="no_identifier"]<-"no_data_crea_id"
data_repo_final$Identifiers_data_creators[data_repo_final$Identifiers_data_creators =="identifier"]<-"data_crea_id"
data_repo_final$Identifiers_related_outputs[data_repo_final$Identifiers_related_outputs =="no_identifier"]<- "no/other_related_output_id"
data_repo_final$Identifiers_related_outputs[data_repo_final$Identifiers_related_outputs =="other_identifier"]<- "no/other_related_output_id"
data_repo_final$Identifiers_related_outputs[data_repo_final$Identifiers_related_outputs =="Cite/Reference"]<- "supplement_and_other_id"
data_repo_final$Identifiers_related_outputs[data_repo_final$Identifiers_related_outputs =="Document"]<- "supplement_and_other_id"
data_repo_final$Identifiers_related_outputs[data_repo_final$Identifiers_related_outputs =="Supplement"]<- "supplement_and_other_id"

data_repo_final<- data_repo_final %>% mutate_if(is.character,as.factor)
summary(data_repo_final)

write.table(data_repo_final, file="Data_DataCite_processed_random.txt", sep = '\\t', row.names = FALSE, col.names = TRUE, quote=FALSE, eol="\\n") # create output file


##-----------------------------------------------------------------------------------------------------------------------------------------------------------

rm(list=ls())
