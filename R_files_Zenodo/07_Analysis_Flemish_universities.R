
### Computational script-------------------------------------------------------------------------------------------------------------------------------------

# This script aims to analyze the data in the file "Data_DataCite_harvest_Flemish_universities.txt". 


## Load libraries--------------------------------------------------------------------------------------------------------------------------------------------

.libPaths("C:/R/library") # Modify the path to where the R libraries are available. If the packages are not yet installed, please install them first.

library("ggplot2")
library("scales")
library("dplyr")


## Load data------------------------------------------------------------------------------------------------------------------------------------------------

data<-read.table(file=file.choose(), header=TRUE, sep="\\t", comment.char="", fill=TRUE, na.strings="character", quote="", row.names=NULL, stringsAsFactors = FALSE, strip.white = TRUE, encoding = "utf-8", blank.lines.skip = TRUE) # cf. the file "Data_DataCite_harvest_Flemish_universities.txt"

data_2<-filter(data[,1:15], Inclusion == "yes" & Versioning_manual == "primary")
data<-data_2


## Analysis--------------------------------------------------------------------------------------------------------------------------------------------------


# Univariate testing

(Distribution_query<-as.data.frame(table(data$Query))[order(as.data.frame(table(data$Query))[,2], decreasing = TRUE),])
(Distribution_publisher<-as.data.frame(table(data$Publisher))[order(as.data.frame(table(data$Publisher))[,2], decreasing = TRUE),])
(Distribution_publicationyear<-as.data.frame(table(data$PublicationYear))[order(as.data.frame(table(data$PublicationYear))[,2], decreasing = TRUE),])
(Distribution_Affiliation_field<-as.data.frame(table(data$Affiliation_information_fields))[order(as.data.frame(table(data$Affiliation_information_fields))[,2], decreasing = TRUE),])

data_2<-data # collapse certain values into broader categories
data_2$Identifiers_related_outputs[data_2$Identifiers_related_outputs!="no_identifier"]<-"identifier"
data_2$Identifiers_data_creators[data_2$Identifiers_data_creators!="no_identifier"]<-"identifier"
data_2$Affiliations_data_creators[data_2$Affiliations_data_creators!="no_affiliation"]<-"affiliation"

(Distribution_related_outputs<-as.data.frame(table(data_2$Identifiers_related_outputs))[order(as.data.frame(table(data_2$Identifiers_related_outputs))[,2], decreasing = TRUE),])
(Distribution_identifiers_creators<-as.data.frame(table(data_2$Identifiers_data_creators))[order(as.data.frame(table(data_2$Identifiers_data_creators))[,2], decreasing = TRUE),])
(Distribution_affiliations_creators<-as.data.frame(table(data_2$Affiliations_data_creators))[order(as.data.frame(table(data_2$Affiliations_data_creators))[,2], decreasing = TRUE),])


# Multivariate testing

table(data_2$Identifiers_data_creators, data_2$PublicationYear)
table(data_2$Affiliation_information_fields, data_2$PublicationYear)
publisher_orcid<-as.data.frame(table(data_2$Publisher, data_2$Identifiers_data_creators))


## Overview of plots----------------------------------------------------------------------------------------------------------------------------------------


# Plot 1
Number_of_DOIs_per_year_and_affiliationfield<-aggregate(DOI ~ PublicationYear + Affiliation_information_fields, data_2, function(x) length(unique(x))) 

Number_of_DOIs_per_year_and_affiliationfield$PublicationYear<-as.numeric(Number_of_DOIs_per_year_and_affiliationfield$PublicationYear)
Number_of_DOIs_per_year_and_affiliationfield$DOI<-as.numeric(Number_of_DOIs_per_year_and_affiliationfield$DOI)
Number_of_DOIs_per_year_and_affiliationfield$Affiliation_information_fields<-as.factor(Number_of_DOIs_per_year_and_affiliationfield$Affiliation_information_fields)

data_points<-as.data.frame(table(Number_of_DOIs_per_year_and_affiliationfield$Affiliation_information_fields))
affiliation_info_major<-filter(data_points, Freq > 2)

Number_of_DOIs_per_year_and_affiliationfield<-filter(Number_of_DOIs_per_year_and_affiliationfield, Affiliation_information_fields %in% affiliation_info_major[,1])

plot_1<-ggplot(Number_of_DOIs_per_year_and_affiliationfield, aes(x=PublicationYear, y= DOI, color = Affiliation_information_fields)) + theme_minimal() + geom_jitter(size=3, width = 0.1) + geom_line(size = 1, alpha = 0.7) + scale_x_continuous(breaks = scales::pretty_breaks(n = 10), limits=c(2005.75,2020.25)) + scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) + theme(axis.text.x=element_text(angle=60, hjust=1)) + xlab("publication year of the dataset") + ylab("number of DOIs") + labs(color = "Affiliation info") + scale_color_manual(values=c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#A65628", "#999999")) # zoom in!

jpeg(file = "Number_of_DOIs_per_year_and_affiliationfield.jpeg", width = 748, height = 634, quality = 100, res = 100)
plot_1
dev.off() 


# Plot 2
Number_of_DOIs_per_year_and_ORCID<-aggregate(DOI ~ PublicationYear + Identifiers_data_creators, data_2, function(x) length(unique(x)))

Number_of_DOIs_per_year_and_ORCID$PublicationYear<-as.numeric(Number_of_DOIs_per_year_and_ORCID$PublicationYear)
Number_of_DOIs_per_year_and_ORCID$DOI<-as.numeric(Number_of_DOIs_per_year_and_ORCID$DOI)
Number_of_DOIs_per_year_and_ORCID$Identifiers_data_creators<-as.factor(Number_of_DOIs_per_year_and_ORCID$Identifiers_data_creators)

plot_2<-ggplot(Number_of_DOIs_per_year_and_ORCID, aes(x=PublicationYear, y= DOI, color = Identifiers_data_creators)) + theme_minimal() + geom_point(size=2) + geom_line(size = 1) + scale_x_continuous(breaks = scales::pretty_breaks(n = 10), limits=c(2005.75,2020.25)) + scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) + theme(axis.text.x=element_text(angle=60, hjust=1)) + xlab("publication year of the dataset") + ylab("number of DOIs") + labs(color = "ORCID") + scale_color_manual(values=c("#e41a1c", "#377eb8")) # zoom in!

jpeg(file = "Number_of_DOIs_per_year_and_ORCID.jpeg", width = 748, height = 634, quality = 100, res = 100)
plot_2
dev.off()


# Plot 3
Number_of_DOIs_per_publisher_and_identifier<-aggregate(DOI ~ Client_ID + Identifiers_data_creators, data_2, function(x) length(unique(x)))

Number_of_DOIs_per_publisher_and_identifier$DOI<-as.numeric(Number_of_DOIs_per_publisher_and_identifier$DOI)
Number_of_DOIs_per_publisher_and_identifier$Client_ID<-as.factor(Number_of_DOIs_per_publisher_and_identifier$Client_ID)
Number_of_DOIs_per_publisher_and_identifier$Identifiers_data_creators<-as.factor(Number_of_DOIs_per_publisher_and_identifier$Identifiers_data_creators)

plot_3<-ggplot(Number_of_DOIs_per_publisher_and_identifier, aes(x= Client_ID, y= DOI))  + geom_bar(stat="identity") + facet_grid(Identifiers_data_creators ~ ., labeller = labeller(Identifiers_data_creators = c(`identifier` = "ORCID", `no_identifier` = "no ORCID"))) + scale_x_discrete(label = function(x) stringr::str_trunc(x, 15)) + coord_flip() + theme_minimal() + theme(strip.background = element_rect(colour="black", fill="white", size=1, linetype="solid")) + xlab("archiving organization") + ylab("number of DOIs")
  
jpeg(file = "Number_of_DOIs_per_publisher_and_identifier.jpeg", width = 748, height = 634, quality = 100, res = 100)
plot_3
dev.off() 


# Plot 4
Number_of_DOIs_per_publisher_and_related_output<-aggregate(DOI ~ Client_ID + Identifiers_related_outputs, data_2, function(x) length(unique(x)))

Number_of_DOIs_per_publisher_and_related_output$DOI<-as.numeric(Number_of_DOIs_per_publisher_and_related_output$DOI)
Number_of_DOIs_per_publisher_and_related_output$Client_ID<-as.factor(Number_of_DOIs_per_publisher_and_related_output$Client_ID)
Number_of_DOIs_per_publisher_and_related_output$Identifiers_related_outputs<-as.factor(Number_of_DOIs_per_publisher_and_related_output$Identifiers_related_outputs)

plot_4<-ggplot(Number_of_DOIs_per_publisher_and_related_output, aes(x= Client_ID, y= DOI))  + geom_bar(stat="identity") + facet_grid(Identifiers_related_outputs ~ ., labeller = labeller(Identifiers_related_outputs = c(`identifier` = "related outputs detected", `no_identifier` = "no related outputs detected"))) + scale_x_discrete(label = function(x) stringr::str_trunc(x, 15)) + coord_flip() + theme_minimal() + theme(strip.background = element_rect(colour="black", fill="white", size=1, linetype="solid")) + xlab("archiving organization") + ylab("number of DOIs")

jpeg(file = "Number_of_DOIs_per_publisher_and_related_output.jpeg", width = 748, height = 634, quality = 100, res = 100)
plot_4
dev.off() 


## Overview of tables----------------------------------------------------------------------------------------------------------------------------------------


# Table 3
Distribution_publicationyear<-as.data.frame(table(data$PublicationYear))[order(as.data.frame(table(data$PublicationYear))[,1], decreasing = FALSE),]
View(Distribution_publicationyear)
sum(Distribution_publicationyear[,2])

# Table 4
Distribution_publisher<-as.data.frame(table(data$Client_ID))[order(as.data.frame(table(data$Client_ID))[,2], decreasing = TRUE),]
View(Distribution_publisher)
sum(Distribution_publisher[,2])


##-----------------------------------------------------------------------------------------------------------------------------------------------------------

rm(list=ls())







