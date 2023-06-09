# This script produces all the files required to deploy an SNSF data story. 
# 
# Data story template: https://github.com/snsf-data/datastory_template
# 
# By running this file, the following components of a data story are generated
# and stored in the output directory: 
# 
# 1) a HTML file (self-contained), which contains all visualizations and 
#   images in encoded form, one for every specified language.
# 2) one file "metadata.json", which contains the metadata essential for 
#   the story (including all language versions in one file).
# 
# The files are stored in output/xxx, where xxx stands for the title of the
# data story in English, how it can also be used for the vanity URL to the 
# story, that means: no special characters, only lowercase.

# Unique name of this data story in English (all lowercase, underscore as 
# white-space, no special characters etc.)
datastory_name <- "70_years_snsf"

# Language-specific names, do adapt! (used for vanity URL! Format: all 
# lowercase, minus as white-space (!) and no special characters, no special 
# characters etc.)
datastory_name_de <- "70-jahre-snf-von-wenigen-millionen-zur-milliarde"
datastory_name_en <- "70-years-snsf-from-a-few-million-to-a-billion"
datastory_name_fr <- "70-ans-du-fns-un-budget-en-forte-hausse"

# English title and lead of the story (Mandatory, even if no EN version)
title_en <- "70 years of the SNSF: from a few million to a billion" 
lead_en <- "The SNSF has been awarding grants to research projects across all scientific disciplines since 1952. How has the amount of funding evolved over time?" 
# German title and lead of the story (Mandatory, even if no DE version)
title_de <- "70 Jahre SNF: Von wenigen Millionen zur Milliarde" 
lead_de <- "Der SNF vergibt seit 1952 Beiträge an Forschungsprojekte in allen wissenschaftlichen Disziplinen. Wie hat sich das Fördervolumen während dieser Zeit entwickelt?" 
# French title and lead of the story (Mandatory, even if no FR version)
title_fr <- "70 ans du FNS : un budget en forte hausse"
lead_fr <- "Depuis 1952, le FNS finance des projets de recherche dans toutes les disciplines scientifiques. Comment le volume des subsides a-t-il évolué au fil des années ?" 
# Contact persons, always (first name + last name)
contact_person <- c("This Rutishauser", "Daniel Schnyder", "Julius Mattern")
# Mail address to be displayed as contact persons, put "datastories@snf.ch" for 
# every name of a contact person listed above.
contact_person_mail <- c("datastories@snf.ch", "datastories@snf.ch", "datastories@snf.ch")
# One of the following categories:  "standard", "briefing", "techreport", 
# "policybrief", "flagship", "figure". Category descriptions are 
datastory_category <- "standard"
# Date, after which the story should be published. Stories not displayed if the 
# date lies in the future. 
publication_date <- "2022-08-10 00:01:00"
# Available language versions in lowercase, possible: "en", "de", "fr".
languages <- c("en", "de", "fr") 
# Whether this story should be a "Feature Story" story
feature_story <- FALSE
# DOI of the story (optional)
doi_url <- "https://doi.org/10.46446/datastory.70-years-snsf-from-a-few-million-to-a-billion" 
# URL to Github page (optional)
github_url <- "https://github.com/snsf-data/70_years_snsf" 
# Put Tag IDs here. Only choose already existing tags. 
tags_ids <- c(260, # research and development
              160) # value of research

# IMPORTANT: Put a title image (as .jpg) into the output directory.
# example: "output/datastory-template.jpg"

# Install snf.datastory package if not available, otherwise load it
if (!require("snf.datastory")) {
  if (!require("devtools")) {
    install.packages("devtools")
    library(devtools)
  }
  install_github("snsf-data/snf.datastory")
  library(snf.datastory)
}

# Load packages
library(tidyverse)
library(scales)
library(conflicted)
library(jsonlite)
library(here)

# Conflict preferences
conflict_prefer("filter", "dplyr")

# Function to validate a mandatory parameter value
is_valid <- function(param_value) {
  if (is.null(param_value)) {
    return(FALSE)
  }
  if (is.na(param_value)) {
    return(FALSE)
  }
  if (str_trim(param_value) == "") {
    return(FALSE)
  }
  return(TRUE)
}

# Validate parameters and throw error message when not correctly filled
if (any(
  !is_valid(datastory_name),
  !is_valid(title_en),
  !is_valid(title_de),
  !is_valid(title_fr),
  !is_valid(datastory_category),
  !is_valid(publication_date),
  length(languages) < 1,
  !is_valid(lead_en),
  !is_valid(lead_de),
  !is_valid(lead_fr)
)) {
  stop("Incorrect value for at least one of the mandatory metadata values.")
}

# Create output directory in main directory
if (!dir.exists(here("output"))) {
  dir.create(here("output"))
}

# Create story directory in output directory
if (!dir.exists(here("output", datastory_name))) {
  dir.create(here("output", datastory_name))
}

# Create a JSON file with the metadata and save it in the output directory
tibble(
  title_en = title_en,
  title_de = title_de,
  title_fr = title_fr,
  author = paste(contact_person, collapse = ";"),
  datastory_category = datastory_category,
  publication_date = publication_date,
  languages = paste(languages, collapse = ";"),
  short_desc_en = lead_en,
  short_desc_de = lead_de,
  short_desc_fr = lead_fr,
  tags = paste(paste0("T", tags_ids, "T"), collapse = ","),
  author_url = paste(contact_person_mail, collapse = ";"),
  top_story = feature_story,
  github_url = github_url,
  doi = doi_url
) %>%
  toJSON() %>%
  write_lines(here("output", datastory_name, "metadata.json"))

# Knit HTML output for each language version
for (idx in seq_len(length(languages))) {
  current_lang <- languages[idx]
  output_file <- here(
    "output", datastory_name,
    paste0(
      str_replace_all(
        get(paste0("datastory_name_", current_lang)), "_", "-"
      ),
      "-", current_lang, ".html"
    )
  )
  print(paste0("Generating output for ", current_lang, " version..."))
  rmarkdown::render(
    input = here(paste0(current_lang, ".Rmd")),
    output_file = output_file,
    params = list(
      title = get(paste0("title_", current_lang)),
      publication_date = publication_date,
      doi = doi_url
    ),
    envir = new.env()
  )
  
  # # Generate raw text version of each story for the search feature
  # readLines(output_file)
  
  # # inputFile <- "
  # con  <- file(output_file, open = "r")
  # while (length(oneLine <- readLines(con, n = 1)) > 0) {
  #   myLine <- unlist((strsplit(oneLine, ",")))
  #   print(myLine)
  # }
  # close(con)
}