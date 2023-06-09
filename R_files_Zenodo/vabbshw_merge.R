# This R script imports the files as data frames and offers some suggestions how to merge the different components of the VABB dataset

# Installation of packages

list.of.packages <- c("data.table", "dplyr", "plyr","stringr","rio")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library(dplyr)
library(data.table)
library(plyr)
library(stringr)
library(rio)

# SET WORK DIRECTORY

myworkdirectory <- ""

# IMPORT FILES (Copy the files in your work directory or change path)

vabb_publications <- read.csv(file=paste0(myworkdirectory,"/","vabb_publications.csv"),sep=",",header=TRUE,encoding = "UTF-8")
vabb_containers <- read.csv(file=paste0(myworkdirectory,"/","vabb_containers.csv"),sep=",",header=TRUE,encoding = "UTF-8")
vabb_additionalid <- read.csv(file=paste0(myworkdirectory,"/","vabb_additionalid.csv"),sep=",",header=TRUE,encoding = "UTF-8")
vabb_authors <- read.csv(file=paste0(myworkdirectory,"/","vabb_authors.csv"),sep=",",header=TRUE,encoding = "UTF-8")
vabb_disciplines <- read.csv(file=paste0(myworkdirectory,"/","vabb_disciplines.csv"),sep=",",header=TRUE,encoding = "UTF-8")
vabb_journals <- read.csv(file=paste0(myworkdirectory,"/","vabb_journals.csv"),sep=",",header=TRUE,encoding = "UTF-8")
vabb_publications_disciplines <- read.csv(file=paste0(myworkdirectory,"/","vabb_publications_disciplines.csv"),sep=",",header=TRUE,encoding = "UTF-8")
vabb_publishers <- read.csv(file=paste0(myworkdirectory,"/","vabb_publishers.csv"),sep=",",header=TRUE,encoding = "UTF-8")
vabb_publications_standardnumbers <- read.csv(file=paste0(myworkdirectory,"/","vabb_publications_standardnumbers.csv"),sep=",",header=TRUE,encoding = "UTF-8")

## MERGE FILES

# Add authors to publications

mytable_publications_authors <- merge(vabb_publications,vabb_authors,by="vabb_publications_id",all.y=TRUE)

# Add disciplines to publications

mytable_publications_disciplines <- merge(vabb_publications,vabb_publications_disciplines,by="vabb_publications_id")
mytable_publications_disciplines <- merge(mytable_publications_disciplines,vabb_disciplines,by="vabb_disciplines_id")

# Add standard numbers to publications

mytable_publications_standardnumbers <- merge(vabb_publications,vabb_publications_standardnumbers,by="vabb_publications_id",all.x=TRUE)

# Split up into a file with journal articles (ISSN) and a file with book publications (ISBN)

mytable_publications_journalarticles <- mytable_publications_standardnumbers %>% filter(vabb_standardnumbers_type == "ISSN")
mytable_publications_bookpublications <- mytable_publications_standardnumbers %>% filter(vabb_standardnumbers_type == "ISBN")

## JOURNALS

# Prepare journals file for merge

journals_a <- vabb_journals %>% select(vabb_standardnumbers_number,vabb_journals_description,vabb_journals_statuspeerreview)
journals_b <- vabb_journals %>% select(vabb_journals_alternativeissn,vabb_journals_description,vabb_journals_statuspeerreview)
setnames(journals_b, old=c("vabb_journals_alternativeissn"),new=c("vabb_standardnumbers_number"))
journals_b <- journals_b %>% filter(!is.na(vabb_standardnumbers_number) & vabb_standardnumbers_number != "")
vabb_journals <- rbind(journals_a,journals_b)
vabb_journals <- unique(vabb_journals)
rm(journals_a)
rm(journals_b)

mytable_publications_journalarticles$vabb_standardnumbers_number <- gsub("([[:digit:]]{4})[[:blank:]]([[:digit:]]{3,})","\\\\1-\\\\2",mytable_publications_journalarticles$vabb_standardnumbers_number)
mytable_publications_journalarticles$vabb_standardnumbers_number <- gsub("([[:digit:]]{4})[[:blank:]]-[[:blank:]]([[:digit:]]{3,})","\\\\1-\\\\2",mytable_publications_journalarticles$vabb_standardnumbers_number)

# Merge journal articles with journal names

mytable_publications_journalarticles <- merge(mytable_publications_journalarticles,vabb_journals,by="vabb_standardnumbers_number",all.x=TRUE)

## BOOKS

mytable_publications_bookpublications <- merge(mytable_publications_bookpublications,vabb_publishers,by="vabb_standardnumbers_prefix",all.x=TRUE)
