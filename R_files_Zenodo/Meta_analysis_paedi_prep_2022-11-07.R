#### 1. Check and set directory #### 
getwd()
#### 2. Name file ####
##Name the file: 

####3. Packages and functions ####
library(readxl)
library(tidyverse)
# library(meta)
# library(metafor)
library(janitor)

## short function for checking if all variables are numeric
is_all_numeric <- function(x){
  !any(is.na(suppressWarnings(as.numeric(na.omit(x))))) & is.character(x)
}

## abbreviations
## AE = adverse event


#### 4. Data set ####

GTT_paed <- read_excel("./Extraction_table_pediatric_AEs_systematic_review_2022-11-18.xlsx")

QAT <- read_excel("./Summary_QAT_2022-11-18.xlsx")


# 5. Data cleaning --------------------------------------------------------
Paed <- clean_names(GTT_paed)
QAT <- clean_names(QAT)
rm(GTT_paed)

## Exclude all rows with irrelevant data
Paed1 <- Paed %>%
  filter(!combinewithinstudy_1_no_2_sum_3_mean_99_exclude == 99)

## Merge QAT PICO information with Paed1
Paed1_qat <- left_join(Paed1, 
                  QAT %>% dplyr::select(pico1_2, study_id), 
                  by = "study_id")

## Matlow is now double because of 2 PICOs
# Selection of the correct PICO for the corresponding sample
Paed2 <- Paed1_qat %>%
  filter(!(study_id == 25 & subgroups %in% "ICU" & pico1_2 %in% "PICO 1")) %>%
  filter(!(study_id == 25 & subgroups %in% "Normal level of inpatient care" & pico1_2 %in% "PICO 2")) %>%
  select(study_id, first_author_s_last_name, publication_year, subgroups, pico1_2, everything())

rm(Paed1, Paed1_qat)

## Exclude Tartaglia, not enough reported data (Number of admissions is not correct) [04.05.2022]
Paed2 <- Paed2 %>%
  filter(!(study_id == 40))

## Change publication year for Stockwell (because we have excluded the other Stockwell with the same year)
Paed2 <- Paed2 %>%
  mutate(publication_year = replace(publication_year, study_id ==36, 2018))

## Correct data from Study ID 19 Kirkendal
# 62 Admissions with AEs (25.8%) instead of 80.
Paed2 <- Paed2 %>%
  mutate(number_of_admissions_with_1_ae = replace(number_of_admissions_with_1_ae, study_id == 19, 62),
         percent_of_admissions_with_1_ae = replace(percent_of_admissions_with_1_ae, study_id ==19, 25.8))

## change Solevag to Solevåg
Paed2 <- Paed2 %>%
  mutate(first_author_s_last_name = replace(first_author_s_last_name, study_id ==32, "Solevåg"))

## add Number of admissions with AEs from Solevag, due to additional author information
Paed2 <- Paed2 %>%
  mutate(number_of_admissions_with_1_ae = replace(number_of_admissions_with_1_ae, study_id ==32, 41))

         
# 5.1 Replacing NA --------------------------------------------------------
Paed2 <- Paed2 %>%
  na_if("NA")

## Variables are saved as characters, change into numeric for numerical variables
Paed3 <- Paed2 %>%
  mutate_if(is_all_numeric, as.numeric)

# reordering the data frame
Paed3 <- Paed3 %>%
  select(study_id, first_author_s_last_name, publication_year, subgroups, pico1_2,
         sample_size, inclusion_of_more_than_one_ae_per_patient, number_of_hospital_days, number_of_a_es,
         number_of_admissions_with_1_ae, no_of_preventable_a_es, number_of_admissions_with_1_preventable_ae, everything())

## adding variable for reporting information
Paed3$Reporting_tag <- NA
Paed3$Reporting_tag2 <- NA
Paed3$Reporting_tag3 <- NA


# 5.1.1 means / sums ------------------------------------------------------
test <- Paed3 %>%
  filter(combinewithinstudy_1_no_2_sum_3_mean_99_exclude == 2)
### there are no studies to calculate the mean
### there are three studies to calculate the sum

## create the sums of all numeric variables (there are only relevant variables extracted, otherwise it would not work this way)
test2 <- test %>%
  group_by(study_id) %>%
  summarise_if(is.numeric, sum)

## adding reporting tag info
test2 <- test2 %>%
  mutate(Reporting_tag = "sub_sum")

## add lost variables due to non.numeric status
test3 <- bind_rows(test, test2)

## filling up relevant information for the sum rows.
test4 <- test3 %>%
  group_by(study_id) %>%
  mutate(first_author_s_last_name= replace(first_author_s_last_name, Reporting_tag %in% "sub_sum", first_author_s_last_name),
         publication_year = replace(publication_year, Reporting_tag %in% "sub_sum", publication_year),
         pico1_2 = replace(pico1_2, Reporting_tag %in% "sub_sum", pico1_2),
         subgroups = replace(subgroups, Reporting_tag %in% "sub_sum", "paediatric sum"),
         inclusion_of_more_than_one_ae_per_patient = replace(inclusion_of_more_than_one_ae_per_patient, Reporting_tag %in% "sub_sum", inclusion_of_more_than_one_ae_per_patient),
         combinewithinstudy_1_no_2_sum_3_mean_99_exclude = replace(combinewithinstudy_1_no_2_sum_3_mean_99_exclude, Reporting_tag %in% "sub_sum", 1))
## there are several warnings but it worked properly.

## selection of final sum score information
test5 <- test4 %>%
  filter(Reporting_tag %in% "sub_sum")

## add sum score studies to data set (keep original study strata)
Paed4 <- bind_rows(Paed3, test5)
rm(test, test2, test3, test4, test5)

## reduce data set to summed studies (get rid of original study strata)
Paed4 <- Paed4 %>%
  filter(combinewithinstudy_1_no_2_sum_3_mean_99_exclude == 1)

# 5.2 Calculations of primary & secondary outcomes --------------------------------------------------------

# 5.2.1 AEs per 100 admissions -----------------------------------------------
## Overview on variables

## Total number of AEs
summary(Paed4$number_of_a_es)
table(Paed4$number_of_a_es)
## 3 NA's [07.11.2022]

## Total number of admissions
summary(Paed4$sample_size)
table(Paed4$sample_size)
# 0 NA's [07.11.2022]

## AEs per 100 admissions
summary(Paed4$a_es_per_100_admissions)
## 23 out of 33 NA's [07.11.2022]

# 5.2.1.1 Calculation Number of AEs-------------------------------------
# AE/100adm * number of patient admissions / 100 = total AE detected

# select all studies that didn't report total AE detected
test <- Paed4 %>%
  filter(is.na(number_of_a_es) | is.na(sample_size)) ## three studies

# select all studies that a calculation is possible
test2 <- test %>%
  filter(!(is.na(sample_size) | is.na(a_es_per_100_admissions))) ## one study (Study ID = 34)

# calculate number of AEs
test3 <- test2 %>%
  mutate(number_of_a_es = replace(number_of_a_es, is.na(number_of_a_es), values = ceiling(a_es_per_100_admissions * sample_size / 100))) %>%
  mutate(Reporting_tag = replace(Reporting_tag, is.na(Reporting_tag), value = "Cal_AE"))

## create full data set without study that was manipulated
test4 <- Paed4 %>%
  filter(!study_id == 34)

## add manipulated study information to full data set
Paed5 <- union(test4, test3)
rm(test, test2, test3, test4)

summary(Paed5$number_of_a_es)
## two remaining studies that will not be included in the calculation for AE/100 admissions [07.11.2022]
## 5 Brennan (will be added in chapter 5.2.1.2)
## 48 Zegers


# 5.2.1.2 Calculation of AEs for Brennan ----------------------------------
## Additional data for Brennan from the Tables:
# Leape et al 1991 (Study ID 22) Table 4
# Rates of AEs for Age group 0-15 total: 12.91 AEs per 1000 discharges --> 1.291 per 100 discharges

## N(AEs) = AE rate * Sample size / 100
(1.291*6661)/100
# = 85.9935

# Brennan et al 1991 (Study ID 5) Table 3
# Rates of AEs according to age
# Newborn, N admissions = 3,595; Crude AE rate 0.6% +/- 0.1 (SE)
# < 15 yr, N admissions = 3,066, Crude AE rate 2.1% +/- 0.4 (SE)

0.006*3595
# = 21.57

0.021*3066
# 64.386

0.006*3595 + 0.021*3066
# 85.956

## Both calculations yield 86 AEs for this sample.

## Add the calculated information in data set
Paed5 <- Paed5 %>%
  mutate(number_of_a_es= replace(number_of_a_es, study_id==5, value = 86))%>%
  mutate(number_of_admissions_with_1_ae = replace(number_of_admissions_with_1_ae, study_id==5, value = 86)) %>%
  mutate(Reporting_tag2 = replace(Reporting_tag2, Reporting_tag %in% "sub_sum" & study_id ==5, "Cal_AE")) %>%
  mutate(Reporting_tag3 = replace(Reporting_tag3, Reporting_tag2 %in% "Cal_AE", "Cal_Adm_AE"))


# 5.2.2 AEs per 1000 patient days -----------------------------------------------
## Overview on variables

## Total number of patient days
summary(Paed5$number_of_hospital_days)
table(Paed5$number_of_hospital_days)
# 10 NA's [07.11.2022]

## AEs per 1000 patient days
summary(Paed5$ae_per_1000_patient_days)
## 20 NA's [07.11.2022]

# 5.2.2.1 Calculation of Number of hospital days --------------------------
# Total AE detected / AE/1000 * 1000 = Total patient days

# select all studies that didn't report total patient days
test <- Paed5 %>%
  filter(is.na(number_of_hospital_days)) ## 10 studies

# select all studies that a calculation is possible
test2 <- test %>%
  filter(!(is.na(ae_per_1000_patient_days)  | is.na(number_of_a_es))) ## 2 studies

# Study ID 35
# 54.9 harms per 1000 patient days (95% CI:48.3–62.3).

# Study ID 36
# 19.0 AEs (95% confidence interval [CI] 17.3–21.0) per 1000 patient days

## calculating total patient days
test3 <- test2 %>%
  mutate(number_of_hospital_days = replace(number_of_hospital_days, is.na(number_of_hospital_days), value= round((number_of_a_es / ae_per_1000_patient_days * 1000),0))) %>%
  mutate(Reporting_tag = replace(Reporting_tag, is.na(Reporting_tag), value="Cal_HospDays"))

## reducing full data set by original studies
test4 <- Paed5 %>%
  filter(!study_id == c(35,36))

## add manipulated studies to full data set
Paed6 <- bind_rows(test4, test3)
rm(test, test2, test3, test4)


summary(Paed6$number_of_hospital_days)
## 8 remaining studies that will not be included in the analysis of 1000 patient days. [07.11.2022]
# 5 Brennan
# 6 Champan
# 23 Letaief
# 27 Maziero
# 33 Sommella
# 34 Soop
# 47 Woods
# 48 Zegers


# 5.2.3 % of admissions with AEs -----------------------------------------------
## Overview on variables

## Total number of admissions with AEs
summary(Paed6$number_of_admissions_with_1_ae)
table(Paed6$number_of_admissions_with_1_ae)
# 3 NA's [07.11.2022]

## Percentage of admissions with AEs
summary(Paed6$percent_of_admissions_with_1_ae)
# 16 NA's [07.11.2022]

# 5.2.3.1 Calculation Number of Admissions with AEs -----------------------
# % of Admissions with AE * Sample size/ 100 = total number of admissions with AEs


# select all studies that didn't report total number of admissions with at least 1 AE
test <- Paed6 %>%
  filter(is.na(number_of_admissions_with_1_ae)) ## 3 studies

# select all studies that a calculation is possible
test2 <- test %>%
  filter(!(is.na(percent_of_admissions_with_1_ae) | is.na(sample_size))) ## 1 study

# calculate total number of admissions with at least 1 AE
test3 <- test2 %>%
  mutate(number_of_admissions_with_1_ae = replace(number_of_admissions_with_1_ae, is.na(number_of_admissions_with_1_ae), value = (percent_of_admissions_with_1_ae * sample_size / 100))) %>%
  mutate(Reporting_tag = replace(Reporting_tag, is.na(Reporting_tag), values = "Cal_Adm_AE"))

## round number of admissions with at least one AE
test3$number_of_admissions_with_1_ae <- round(test3$number_of_admissions_with_1_ae, digits = 0)

## reduce full data set by original study
test4 <- Paed6 %>%
  filter(!study_id == 11)

## add manipulated study to full data set
Paed7 <- bind_rows(test4, test3)

rm(test, test2, test3, test4)

summary(Paed7$number_of_admissions_with_1_ae)
## 1 study will be excluded from the calculation of % of admissions with AEs.
# 31 Sharek
# 34 Soop --> AE = admissions (see 5.2.3.2)


# 5.2.3.2 Adding number of admissions for Soop (ID 34) ----------------------------
## Each admissions has max 1 AE, therefore number of admissions with 1 AE = number of AEs.
Paed7 <- Paed7 %>%
  mutate(number_of_admissions_with_1_ae= replace(number_of_admissions_with_1_ae, study_id==34, value = number_of_a_es[study_id==34]))%>%
  mutate(Reporting_tag2 = replace(Reporting_tag2, Reporting_tag %in% "Cal_AE" & study_id ==34, "Cal_Adm_AE"))

# 5.2.4 % of preventable AEs from total AEs detected-----------------------------------------------
## Total number of preventable AEs
summary(Paed7$no_of_preventable_a_es)
table(Paed7$no_of_preventable_a_es)
# 19 NA's [07.11.2022]

## Percentage of preventable AEs
summary(Paed7$percent_of_preventable_a_es_out_of_all_a_es)
# 21 NA's [07.11.2022]


# 5.2.4.1 Calculation of Number of preventable AEs ------------------------
# % of prev. AEs * Total AE detected / 100 = total number of preventable AEs

# select all studies that didn't report total preventable AE
test <- Paed7 %>%
  filter(is.na(no_of_preventable_a_es)) ## 19 studies

# select all studies that a calculation is possible
test2 <- test %>%
  filter(!(is.na(number_of_a_es) | is.na(percent_of_preventable_a_es_out_of_all_a_es))) ## 1 study

# calculate total preventable AE detected
test3 <- test2 %>%
  mutate(no_of_preventable_a_es = replace(no_of_preventable_a_es, is.na(no_of_preventable_a_es), value= (percent_of_preventable_a_es_out_of_all_a_es *number_of_a_es/100 ))) %>%
  mutate(Reporting_tag = replace(Reporting_tag, is.na(Reporting_tag), value = "Cal_preAE"))

## rounding
test3$no_of_preventable_a_es <- round(test3$no_of_preventable_a_es,0)

## reduce full data set by original study
test4 <- Paed7 %>%
  filter(!study_id == 43)

## add manipulated study to full data
Paed8 <- bind_rows(test4, test3)

summary(Paed8$no_of_preventable_a_es)
## 18 (Matlow 2x) studies will be excluded from the calculation of % of preventable AEs.

rm(test, test2, test3, test4)


# 5.2.4.2 Calculation preventable AEs for Study Soop ----------------------
## Soop included only 1 AE per patient.
## They reported preventable AEs / 100 admissions.
## Therefore we can calculate the total number of preventable AEs.

test <- Paed8 %>%
  filter(study_id ==34)

# prev AEs/100adm * number of patient admissions / 100 = total preventable AE detected
test3 <- test %>%
  mutate(no_of_preventable_a_es = replace(no_of_preventable_a_es, is.na(no_of_preventable_a_es), values = ceiling(preventable_a_es_per_100_admissions * sample_size / 100))) %>%
  mutate(Reporting_tag3 = replace(Reporting_tag3, is.na(Reporting_tag3), value = "Cal_preAE"))

## reduce full data set by original study
test4 <- Paed8 %>%
  filter(!study_id == 34)

## add manipulated study to full data
Paed8a <- bind_rows(test4, test3)
Paed8 <- Paed8a ## after check that everything was calculated correctly

rm(test, test3, test4, Paed8a)


# 5.2.5 % admissions with preventable AEs-----------------------------------------------
## Overview on variables

## Total number of admissions with preventable AEs
summary(Paed8$number_of_admissions_with_1_preventable_ae)
table(Paed8$number_of_admissions_with_1_preventable_ae)
# 22 NA's [07.11.2022]

##  Percentage of admissions with preventable AEs
summary(Paed8$percent_of_admissions_with_1_preventable_ae)
# 29 NA's [07.11.2022]


# 5.2.5.1 Calculation Number of Admissions with prev AEs -----------------------
# % of Admissions with prev AE * Sample size/ 100 = total number of admissions with prev AEs

# select all studies that didn't report total number of admissions with at least 1 prev AE
test <- Paed8 %>%
  filter(is.na(number_of_admissions_with_1_preventable_ae)) ## 22 studies

# select all studies that a calculation is possible
test2 <- test %>%
  filter(!(is.na(percent_of_admissions_with_1_preventable_ae) | is.na(sample_size)))
# 0 studies where number of admissions with preventable AEs can be calculated

## Only 10 studies reporting on admissions with preventable AEs:
# 2 Agarwal
# 8 Davis
# 17 Jorro-Baron
# 21 Larsen
# 25 Matlow
# 36 Stockwell
# 42 Unbeck
# 44 Verlaat
# 46 Wilson
# 47 Woods

rm(test, test2)

# 5.2.5.2 Adding number of admissions with prev AEs for Soop ----------------------------
Paed8$Reporting_tag4 <- NA

## Number of admissions with preventable AEs = Number of preventable AEs.
Paed8b <- Paed8 %>%
  mutate(number_of_admissions_with_1_preventable_ae= replace(number_of_admissions_with_1_preventable_ae, 
                                                             study_id==34, value = no_of_preventable_a_es[study_id==34]))%>%
  mutate(Reporting_tag4 = replace(Reporting_tag4, study_id ==34, "Cal_Adm_prevAE"))

Paed8 <- Paed8b # first check if calculated correctly

rm(Paed3, Paed4, Paed5, Paed6, Paed7, Paed8b)


# 5.3 Adding information from QATs ----------------------------------------
## Exclude Tartaglia (ID 40)
QAT <- QAT %>%
  filter(!(study_id == 40)) 

# 5.3.1 Adding corrections to RoB -----------------------------------------
QAT2 <- QAT %>%
  mutate(outcome_risk = replace(outcome_risk, study_id == 19, "Low"),
         outcome_risk = replace(outcome_risk, study_id == 28, "High"),
         outcome_risk = replace(outcome_risk, study_id == 32, "Low"))
## Corrections after final version of QAT assessments [24.05.2022]


# 5.3.2 Merging QAT and Paed ----------------------------------------------

## Merge QAT results with Paed8
Paed9 <- left_join(Paed8, 
                       QAT2 %>% dplyr::select(patient_selection_risk, patient_selection_concern, 
                                             reviewer_risk, reviewer_concern,
                                             rr_process_risk, rr_process_concern,
                                             outcome_risk, outcome_concern,
                                             flow_timing_risk,
                                             overall_risk_of_bias, overall_concern, 
                                             study_id), 
                       by = "study_id")

## Matlow are doubled again
Paed9 <- Paed9 %>%
  filter(!(study_id == 25 & subgroups %in% "ICU" & is.na(overall_concern))) %>%
  filter(!(study_id == 25 & subgroups %in% "Normal level of inpatient care" & is.na(overall_concern))) %>%
  select(study_id, first_author_s_last_name, publication_year, subgroups, pico1_2, everything())

rm(Paed8, QAT)

# 5.4 Differentiation into different analysis groups -------------------------------------------------

# 5.4.1 Inclusion of >= 1 AE ------------------------------------------------
table(Paed9$inclusion_of_more_than_one_ae_per_patient)

## Overview of all studies applying GTT/TT methodology
Paed_GTT <- Paed9 %>%
  filter(inclusion_of_more_than_one_ae_per_patient == "Yes")

## Overview of all studies applying HMPS methodology
Paed_HMPS <- Paed9 %>%
  filter(inclusion_of_more_than_one_ae_per_patient == "No")


# 5.4.2. Population (general care vs. intensive care) ---------------------------------------------------------------
table(Paed9$pico1_2)
table(Paed_GTT$pico1_2)
table(Paed_HMPS$pico1_2)

## Relevel the character variables for display
Paed9$pico1_2 <- factor(Paed9$pico1_2, levels = c("PICO 1", "PICO 2"))

## Renaming variable pico1_2 to Population
Paed9 <- Paed9 %>%
  rename(Population = pico1_2)

## Rename variables of the population
Paed9 <- Paed9 %>%
  mutate(Population=recode(Population, 
                  "PICO 1"="General care",
                  "PICO 2"="Intensive care"))


# 5.5 Adding asterisk -----------------------------------------------------
## For displaying the asterisk in the forest plots, those are added as variables

## Asterisk for % of admissions with AEs
Paed9$Tag1_AEadm <- NA
Paed9$Tag2_AEadm <- NA

## Adding Asterisk for a sum score manipulation and for the calculation of admissions with AEs
Paed9a <- Paed9 %>%
  mutate(Tag1_AEadm = replace(Tag1_AEadm, Reporting_tag %in% "sub_sum", "\\u0024" ),
         Tag1_AEadm = replace(Tag1_AEadm, Reporting_tag %in% "Cal_Adm_AE", "\\u00A5"),
         Tag1_AEadm = replace(Tag1_AEadm, Reporting_tag2 %in% "Cal_Adm_AE", "\\u00A5"),
         Tag2_AEadm = replace(Tag2_AEadm, Reporting_tag3 %in% "Cal_Adm_AE", "\\u00A5"),
         Tag1_AEadm = replace(Tag1_AEadm, study_id == 8 | study_id == 46, "\\u00A2"))

## Asterisk for AEs per 100 admissions
Paed9a$Tag1_AE100 <- NA
Paed9a$Tag2_AE100 <- NA

## Adding Asterisk for sum score and calculation of total number of AEs.
Paed9b <- Paed9a %>%
  mutate(Tag1_AE100 = replace(Tag1_AE100, Reporting_tag %in% "Cal_AE", "\\u0023"),
         Tag1_AE100 = replace(Tag1_AE100, Reporting_tag %in% "sub_sum", "\\u0024"),
         Tag2_AE100 = replace(Tag2_AE100, Reporting_tag2 %in% "Cal_AE", "\\u0023"),
         Tag1_AE100 = replace(Tag1_AE100, study_id == 8 | study_id == 46, "\\u00A2"))

## Asterisk for AEs per 1000 patient days
Paed9b$Tag1_AE1000 <- NA
Paed9b$Tag2_AE1000 <- NA

## Adding Asterisk for sum score, calculation of total number of AEs, calculation of total number of hospital days.
Paed9c <- Paed9b %>%
  mutate(Tag1_AE1000 = replace(Tag1_AE1000, Reporting_tag %in% "Cal_AE", "\\u0023"),
         Tag1_AE1000 = replace(Tag1_AE1000, Reporting_tag %in% "sub_sum", "\\u0024"),
         Tag2_AE1000 = replace(Tag2_AE1000, Reporting_tag2 %in% "Cal_AE", "\\u0023"),
         Tag1_AE1000 = replace(Tag1_AE1000, Reporting_tag %in% "Cal_HospDays", "\\u00A3"),
         Tag1_AE1000 = replace(Tag1_AE1000, study_id == 8 | study_id == 46, "\\u00A2"))

## Asterisk for % of preventable AEs
Paed9c$Tag1_prevAE <- NA
Paed9c$Tag2_prevAE <- NA

## Asterisk for sum score, calculation of AEs, calculation of preventable AEs
Paed9d <- Paed9c %>%
  mutate(Tag1_prevAE = replace(Tag1_prevAE, Reporting_tag %in% "Cal_AE", "\\u0023"),
         Tag1_prevAE = replace(Tag1_prevAE, Reporting_tag %in% "sub_sum", "\\u0024"),
         Tag2_prevAE = replace(Tag2_prevAE, Reporting_tag2 %in% "Cal_AE", "\\u0023"),
         Tag1_prevAE = replace(Tag1_prevAE, Reporting_tag %in% "Cal_preAE", "\\u0026"),
         Tag1_prevAE = replace(Tag1_prevAE, study_id == 8 | study_id == 46, "\\u00A2"),
         Tag2_prevAE = replace(Tag2_prevAE, Reporting_tag3 %in% "Cal_preAE", "\\u0026"))

## Asterisk for admissions with preventable AEs
Paed9d$Tag1_prevAdmAE <- NA

## Asterisk for sum score, calculation of admissions with preventable AEs
Paed9e <- Paed9d %>%
  mutate(Tag1_prevAdmAE = replace(Tag1_prevAdmAE, Reporting_tag4 %in% "Cal_Adm_prevAE", "\\u00B6"),
         Tag1_prevAdmAE = replace(Tag1_prevAdmAE, Reporting_tag %in% "sub_sum", "\\u0024"),
         Tag1_prevAdmAE = replace(Tag1_prevAdmAE, study_id == 8 | study_id == 46, "\\u00A2"))

## Uniting the asterisks and removing all the remaining variables
Paed10 <- Paed9e %>%
  unite(Study_asterisk_AE100, Tag1_AE100, Tag2_AE100, sep = "", remove = T, na.rm = T) %>%
  unite(Study_asterisk_AE1000, Tag1_AE1000, Tag2_AE1000, sep = "", remove = T, na.rm = T) %>%
  unite(Study_asterisk_AEadm, Tag1_AEadm, Tag2_AEadm, sep = "", remove = T, na.rm = T) %>%
  unite(Study_asterisk_prevAE, Tag1_prevAE, Tag2_prevAE, sep = "", remove = T, na.rm = T)

rm(Paed9a, Paed9b, Paed9c, Paed9d, Paed9e)  


# 5.6 Calculation of outcome measures -------------------------------------
Paed10 <- Paed10 %>%
  mutate(AEadm_cal = number_of_admissions_with_1_ae/sample_size*100,
         AE_100_cal = number_of_a_es/sample_size*100,
         AE_1000_cal = number_of_a_es/number_of_hospital_days*1000,
         AEprev_cal = no_of_preventable_a_es/number_of_a_es*100,
         AEadmprev_cal = number_of_admissions_with_1_preventable_ae/sample_size*100)


# 5.7 Adding and revising Information on RRR Method -----------------------
## Information on the RRR method was only recorded on one row of the studies.
## In step 5. Data cleaning currently (25.08.22) Line 53, this information is excluded.
## Some of the information was updated after the analysis, which will here be corrected.
table(Paed10$rrr_method)

Paed11 <- Paed10 %>%
  mutate(rrr_method = replace(rrr_method, study_id == 5, "HMPS"), ## was lost during selection
         rrr_method = replace(rrr_method, study_id == 7, "Global Trigger tool (GTT)"),  ## changed from TT to GTT
         rrr_method = replace(rrr_method, study_id == 8, "HMPS"),  ## was lost during selection
         rrr_method = replace(rrr_method, study_id == 23, "HMPS"),  ## was lost during selection
         rrr_method = replace(rrr_method, study_id == 25, "Trigger tool (TT)"), ## was lost during selection
         rrr_method = replace(rrr_method, study_id == 28, "Global Trigger tool (GTT)"),  ## changed from TT to GTT
         rrr_method = replace(rrr_method, study_id == 33, "HMPS"),  ## was lost during selection
         rrr_method = replace(rrr_method, study_id == 34, "HMPS"),  ## was lost during selection
         rrr_method = replace(rrr_method, study_id == 46, "HMPS"),  ## was lost during selection
         rrr_method = replace(rrr_method, study_id == 48, "HMPS")  ## was lost during selection
         )

# rrr <- Paed11 %>%
#   select(study_id, first_author_s_last_name, rrr_method)
table(Paed11$rrr_method)

table(Paed11$rrr_method, Paed11$Population)


# 6. Exporting final adapted data ----------------------------------------
Paed11a <- Paed11 %>%
  select(study_id, first_author_s_last_name, publication_year, Population, sample_size,
         inclusion_of_more_than_one_ae_per_patient, number_of_hospital_days, number_of_a_es,
         number_of_admissions_with_1_ae, number_of_admissions_with_1_preventable_ae, no_of_preventable_a_es,
         AEadm_cal, AE_100_cal, AE_1000_cal, AEprev_cal, AEadmprev_cal, percent_of_admissions_with_1_ae,
         a_es_per_100_admissions, ae_per_1000_patient_days, percent_of_admissions_with_1_preventable_ae,
         percent_of_preventable_a_es_out_of_all_a_es, rrr_method)

# write.csv(Paed11a, file = "./11_GTT_Paedi_Data/GTT_paedi_overview_2022-11-18.csv")

rm(Paed11a)
