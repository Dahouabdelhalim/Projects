
#### Script 1 - Data Preparation 
#### Incorporate external guidelines scores into journal scores to account for journal referral to external guidelines using hierarchical approach 

# Load packages
library(tidyr)
library(dplyr)

# Import journal data file (Relative Path)
journals<-read.csv("JournalData.csv", stringsAsFactors=FALSE, na.strings = "")%>%
  select(complete:Reference...Notes.10)%>% # omit NA columns that get imported 
  filter(!is.na(journ_tit))

# Omit excluded journals
journals<-journals[!(journals$complete=="exclude"),]
journals<-journals[!(journals$complete=="invert-exclude"),]

# Separate so that there is 1 column called 'external level' (columns ex_1, ex_2, ex_3) and 1 columnn for 'external guideline' (the guideline that matches that column)
journals_expanded<-gather(journals, "external_level", "external_guideline", ex_1:ex_3)%>%
  mutate(external_guideline=sub("\\\\s+$", "", external_guideline), # Remove trailing blanks from external_guideline
         external_level=substr(external_level, 4, 4))%>%  # Turn "ex_1, ex_2, ex_3" into "1, 2, 3"
  separate_rows(external_guideline, sep=";")  # Separate out cases where >1 external guideline listed in the same cell 


# Import External guidelines
external<-read.csv("C:\\\\Users\\\\field\\\\Dropbox\\\\Grad School\\\\Journal Scoring Project\\\\Data\\\\Newest\\\\ExternalGuidelines.csv", stringsAsFactors=FALSE,na.strings = "")%>%
  filter(!is.na(ex_tit))
external<-external[!(external$complete=="exclude"),]

# Join the tables
journals_and_external<-left_join(journals_expanded, external, by=c("external_guideline"="ex_tit")) #left_join keeps all value from "left" data.frame ("journals expanded"), and any matching values on "right" data.frame ("external"); by=c tells it which values to join by (first value is for left table, second for right)

journals_and_external[journals_and_external=="n/a"]<-0 # Replace everything imputed as "n/a" with 0

journals_and_external[is.na(journals_and_external)]<-0 # Replace blanks with 0 

# For each external guideline, find 'path of least resistance'
journal_and_external_combined<-mutate(journals_and_external, 
                                      ex_nat_least_resistance=pmin(external_level, ex_nat),
                                      ex_ac_approve_least_resistance=pmin(external_level, ex_ac_approve),
                                      ex_ac_inst_least_resistance=pmin(external_level, ex_ac_inst),
                                      ex_ac_num_least_resistance=pmin(external_level, ex_ac_num),
                                      ex_state_leg_least_resistance=pmin(external_level, ex_state_leg),
                                      ex_state_nat_least_resistance=pmin(external_level, ex_state_nat),
                                      ex_any_ac_state_YN_least_resistance=pmin(external_level, ex_any_ac_state_YN),
                                      ex_pain_dis_suff_YN_least_resistance=pmin(external_level, ex_pain_dis_suff_YN),
                                      ex_field_YN_least_resistance=pmin(external_level, ex_field_YN),
                                      ex_three_r_least_resistance=pmin(external_level, ex_three_r),
                                      ex_leg_least_resistance=pmin(external_level, ex_leg))

# In journals with >1 external guideline, find the *strictest* of the 'paths of least resistance' for each given criterion
journal_and_external_strictest_external<-journal_and_external_combined%>%
  group_by(journ_tit)%>%  # within each journal,
  summarize(ex_state_nat_strictest_least_path=max(ex_state_nat_least_resistance),
            ex_state_leg_strictest_least_path=max(ex_state_leg_least_resistance),
            ex_ac_num_strictest_least_path=max(ex_ac_num_least_resistance),
            ex_ac_inst_strictest_least_path=max(ex_ac_inst_least_resistance),
            ex_ac_approve_strictest_least_path=max(ex_ac_approve_least_resistance),
            ex_any_ac_state_YN_strictest_least_path=max(ex_any_ac_state_YN_least_resistance),
            ex_nat_strictest_least_path=max(ex_nat_least_resistance),
            ex_pain_dis_suff_YN_strictest_least_path=max(ex_pain_dis_suff_YN_least_resistance),
            ex_field_YN_strictest_least_path=max(ex_field_YN_least_resistance),
            ex_three_r_strictest_least_path=max(ex_three_r_least_resistance),
            ex_leg_strictest_least_path=max(ex_leg_least_resistance))

# Merge incorporated external scorings back to original journal spreadsheet
journal_with_sorted_out_external<-left_join(journals, journal_and_external_strictest_external, by="journ_tit")

# Use strictest out of internal journal criterion vs. external criterion ("final" scoring)
journal_and_external_strictest_overall<-mutate(journal_with_sorted_out_external,
                                               final_state_nat=pmax(ex_state_nat_strictest_least_path, state_nat),
                                               final_state_leg=pmax(ex_state_leg_strictest_least_path, state_leg),
                                               final_ac_num=pmax(ex_ac_num_strictest_least_path, ac_num),
                                               final_ac_inst=pmax(ex_ac_inst_strictest_least_path, ac_inst),
                                               final_ac_approve=pmax(ex_ac_approve_strictest_least_path, ac_approve),
                                               final_nat=pmax(ex_nat_strictest_least_path, nat),
                                               final_leg=pmax(ex_leg_strictest_least_path, leg),
                                               final_any_ac_state_YN=pmax(ex_any_ac_state_YN_strictest_least_path, any_ac_state_YN),
                                               final_pain_dis_suff_YN=pmax(ex_pain_dis_suff_YN_strictest_least_path, pain_dis_suff_YN),
                                               final_field_YN=pmax(ex_field_YN_strictest_least_path, field_YN),
                                               final_three_r=pmax(ex_three_r_strictest_least_path, three_r))%>%
  plyr::rename(c("oa"="oa_YN", "interdis"="interdis_YN", "conserv"="conserv_YN"))%>%
  select(journ_tit, if., country, oa_YN, interdis_YN, conserv_YN, final_ac_num, final_ac_inst, final_ac_approve, final_any_ac_state_YN, final_field_YN, final_three_r, condition)
  
regression_data <- journal_and_external_strictest_overall

# Assign animal welfare law score to each country according to Animal Protection Index 

regression_data$welfare_law <- NA
regression_data$welfare_law <- regression_data$country
regression_data$welfare_law[regression_data$welfare_law=="Brazil"] <- 0
regression_data$welfare_law[regression_data$welfare_law=="Canada"] <- 1
regression_data$welfare_law[regression_data$welfare_law=="China"] <- 0
regression_data$welfare_law[regression_data$welfare_law=="France"] <- 0
regression_data$welfare_law[regression_data$welfare_law=="Mexico"] <- 0
regression_data$welfare_law[regression_data$welfare_law=="Pakistan"] <- 0
regression_data$welfare_law[regression_data$welfare_law=="Russia"] <- 0
regression_data$welfare_law[regression_data$welfare_law=="South Africa"] <- 0
regression_data$welfare_law[regression_data$welfare_law=="Spain"] <- 0
regression_data$welfare_law[regression_data$welfare_law=="USA"] <- 0

regression_data$welfare_law[regression_data$welfare_law=="Australia"] <- 1
regression_data$welfare_law[regression_data$welfare_law=="Austria"] <- 1
regression_data$welfare_law[regression_data$welfare_law=="Belgium"] <- 1
regression_data$welfare_law[regression_data$welfare_law=="Bulgaria"] <- 1
regression_data$welfare_law[regression_data$welfare_law=="Chile"] <- 1
regression_data$welfare_law[regression_data$welfare_law=="Colombia"] <- 1
regression_data$welfare_law[regression_data$welfare_law=="Czech Republic"] <- 1
regression_data$welfare_law[regression_data$welfare_law=="Denmark"] <- 1
regression_data$welfare_law[regression_data$welfare_law=="England"] <- 1
regression_data$welfare_law[regression_data$welfare_law=="Finland"] <- 1
regression_data$welfare_law[regression_data$welfare_law=="Germany"] <- 1
regression_data$welfare_law[regression_data$welfare_law=="Hungary"] <- 1
regression_data$welfare_law[regression_data$welfare_law=="Italy"] <- 1
regression_data$welfare_law[regression_data$welfare_law=="Japan"] <- 1
regression_data$welfare_law[regression_data$welfare_law=="Kenya"] <- 1
regression_data$welfare_law[regression_data$welfare_law=="Netherlands"] <- 1
regression_data$welfare_law[regression_data$welfare_law=="New Zealand"] <- 1
regression_data$welfare_law[regression_data$welfare_law=="Poland"] <- 1
regression_data$welfare_law[regression_data$welfare_law=="Romania"] <- 1
regression_data$welfare_law[regression_data$welfare_law=="Singapore"] <- 1
regression_data$welfare_law[regression_data$welfare_law=="Slovakia"] <- 1
regression_data$welfare_law[regression_data$welfare_law=="Switzerland"] <- 1
regression_data$welfare_law[regression_data$welfare_law=="Taiwan"] <- 1
regression_data$welfare_law[regression_data$welfare_law=="Turkey"] <- 1


