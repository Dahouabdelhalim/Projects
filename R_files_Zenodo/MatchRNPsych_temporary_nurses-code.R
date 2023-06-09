## Title: The deployment of temporary nurses and its association with permanent nurses' outcomes in Swiss psychiatric hospitals: A secondary analysis.
## Author: Leonel Oliveira
## Co-Authors:  B. Gehri, M. Simon
## Data from MatchRN Psychiatry (https://matchrnpsychiatrie.nursing.unibas.ch/)
## Protocol: https://www.researchprotocols.org/2021/8/e26700
## University of Basel, Faculty of Medicine, Department of Public Health, Institute of Nursing Science
## Date and Version: 20.11.2022 Version 1.1

### load libraries and data set ====
# load libraries
library(tidyverse)
library(gtsummary)
library(psych)
library(lme4)
library(sjPlot)

# Setwd and place files there

setwd("")

# load data set
MRNPsych_Unit<-read.csv("MRNP_UNIT.csv", sep = ",", header = T)
MRNPsych_Nurses<-read.csv("MRNP_NURSES.csv", sep = ",", header = T)

### data preparation ====
## create original variables
MRNPsych_Nurses<-MRNPsych_Nurses %>%
  mutate(across(everything(), list(orig = ~.))) %>%
  select(-unit_code_orig)

MRNPsych_Unit<-MRNPsych_Unit %>%
  mutate(U.E.3_orig = U.E.3)

## correct variables 
# age, somatic, turnover, total registered nurses and frequency of temporary nurses use
MRNPsych_Nurses<-MRNPsych_Nurses %>%
  mutate(H.2 = case_when(H.2 == 0 ~ NA_integer_,                                # age 0 not possible
                         T ~ H.2),
         D.2.9.c = case_when(D.2.8 >= 1 & is.na(D.2.9.c) ~ as.integer(0),       # if registered nurse was responsible for more than one patient, we assume that a missing somatic answer is considered as none = 0
                             D.2.9.c >= 40 ~ NA_integer_,                       # if more than 40 patients are somatic we considered it as a wrong entry = NA
                             T ~ as.integer(D.2.9.c)),
         D.2.9.e = case_when(D.2.8 >= 1 & is.na(D.2.9.e) ~ as.integer(0),       # if registered nurse was responsible for more than one patient, we assume that a missing admission/discharge answer is considered as none = 0
                             T ~ as.integer(D.2.9.e)),
         D.2.11 = case_when(D.2.11 == 0 & is.na(D.2.11) ~ NA_integer_,          # if total registered nurses on most recent shift was 0, we coded it as NA.
                            T ~ as.integer(D.2.11)),
         across(A.8.1:A.8.5, ~ case_when(. == 1 ~ 5, . == 2 ~ 4, . == 3 ~ 3, . == 4 ~ 2, . == 5 ~ 1)))    # reverse burnout for better understanding and interpretation of the results  

MRNPsych_Unit<-MRNPsych_Unit %>%
  mutate(U.E.3 = case_when(U.E.3 == 6 ~ as.integer(0),
                           T ~ as.integer(U.E.3)))

## create numeric variables 
# employment percentage, patient to nurse ratio, somatic ratio, turnover ratio, mean score burnout and mean score leadership
MRNPsych_Nurses<-MRNPsych_Nurses %>%
mutate(Perc = case_when(is.na(H.8.1) & H.8 == 1 ~ as.integer(100),              # combine the two answers to employment percentage, since working as full-time is considered to be 100%. If both questions were answered, we considred the second answer as the employment percentage.
                        is.na(H.8.1) & H.8 == 2 ~ NA_integer_,
                        is.na(H.8.1) & is.na(H.8) ~ NA_integer_,
                        !is.na(H.8.1) & H.8 == 1 ~ as.integer(H.8.1),
                        !is.na(H.8.1) & H.8 == 2 ~ as.integer(H.8.1),
                        !is.na(H.8.1) & is.na(H.8) ~ as.integer(H.8.1)),
       PNR = case_when((D.2.10/D.2.11) == Inf ~ NA_real_,
                       T ~ D.2.10/D.2.11),
       Skill = (D.2.11/(D.2.11+D.2.12)*100),
       Somatic = case_when((D.2.9.c/D.2.8) == Inf | (D.2.9.c/D.2.8) > 1 ~ NA_real_,    # infinite and ratio values over 1 (since not possible) coded as NA
                           T ~ (D.2.9.c/D.2.8)*100),
       Turnover = case_when((D.2.9.e/D.2.8) == Inf | (D.2.9.e/D.2.8) > 1 ~ NA_real_,   # infinite and ratio values over 1 (since not possible) coded as NA
                            T ~ (D.2.9.e/D.2.8)*100),
       Burnout = rowMeans(across(c(A.8.1:A.8.5)),na.rm = T),
       Leadership = rowMeans(across(c(A.1.2,A.1.6,A.1.8,A.1.12)),na.rm = T))

## create factor variables 
# change shift and employment percentage into factors
MRNPsych_Nurses<-MRNPsych_Nurses %>%
  mutate(Shift = factor(D.2.1, levels = 1:3, labels = c("Early shift",
                                                        "Late shift",
                                                        "Night shift")),
         Sex = factor(H.1, levels = 1:2, labels = c("Female",
                                                    "Male")),
         Perc = case_when(Perc <= 60 ~ "<60%",
                          Perc >= 61 & Perc <= 95 ~ "61%-95%",
                          Perc >= 96 ~ "96%-100%",
                          T ~ NA_character_),
         Perc = fct_relevel(Perc, "<60%"))

## create dichotomous variables 
# for temporary nurses' frequency
MRNPsych_Unit<-MRNPsych_Unit %>%
  mutate(TempFreq = case_when(U.E.3 >= 3 ~ "Frequent",
                              U.E.3 <= 2 ~ "Occasional"),
         TempFreq = fct_relevel(TempFreq, "Occasional"))

## create remaining variables
# job satisfaction, intention to leave organisation, intention to leave profession and age
MRNPsych_Nurses<-MRNPsych_Nurses %>%
  mutate(JobSat = A.2,
         ITLJob = A.3,
         ITLProf = A.4,
         Age = H.2)

## select just registred nurses
MRNPsych_Nurses<-MRNPsych_Nurses %>%
  filter(H.3 <= 3)

## select just complete cases
MRNPsych_Unit<-MRNPsych_Unit %>%
  drop_na(TempFreq)

### adjusted staffing ====
## calculating adjusted staffing
# create data set for the adjusted staffing model
MRNPsych_Model<-MRNPsych_Nurses %>%
  select(unit_code,PNR,Skill,Shift,Somatic,Turnover) %>%
  drop_na(PNR,Skill,Shift,Somatic,Turnover) %>%
  left_join(., MRNPsych_Unit %>% select(unit_code,TempFreq), by = c("unit_code"))

# adjusted staffing model
Staffing_MLM<-lmer(PNR ~ Skill + Shift + Somatic + Turnover + (1 | unit_code), data = MRNPsych_Model)

tab_model(Staffing_MLM,
          pred.labels = c("Intercept","Skill & Grade Mix","Late shift","Night shift","Somatic diagnoses ratio","Turnover ratio"),
          dv.labels = c("Patient-to-nurse ratio"),
          string.pred = "Coefficient",
          string.ci = "CI (95%)",
          p.style = "stars",
          digits = 3,
          file = "adjusted_staffing_model.docx")

# testing adjusted staffing by adding temporary nurses' frequency
Staffing_MLM_Temp<-lmer(PNR ~ Skill + Shift + Somatic + Turnover + TempFreq + (1 | unit_code), data = MRNPsych_Model)

tab_model(Staffing_MLM_Temp,
          pred.labels = c("Intercept","Skill & Grade Mix","Late shift","Night shift","Somatic diagnoses ratio","Turnover ratio","Frequent deployment of temporary"),
          dv.labels = c("Patient-to-nurse ratio"),
          string.pred = "Coefficient",
          string.ci = "CI (95%)",
          p.style = "stars",
          digits = 3,
          file = "adjusted_staffing_model_tempFreq.doc")

## calculate average PNR for Early shift with average predictors.
round(Staffing_MLM@beta[1] + (Staffing_MLM@beta[2]  * mean(MRNPsych_Model$Skill)) + (Staffing_MLM@beta[5] * mean(MRNPsych_Model$Somatic)) + (Staffing_MLM@beta[6] * mean(MRNPsych_Model$Turnover)),2)

Staffing_Example<-paste0(round(Staffing_MLM@beta[1]+(Staffing_MLM@beta[2]*mean(Staffing_MLM@frame$Skill))+(Staffing_MLM@beta[5]*mean(Staffing_MLM@frame$Somatic))+(Staffing_MLM@beta[6]*mean(Staffing_MLM@frame$Turnover)),2),
                         " (95%-CI = ",
                         round(confint(Staffing_MLM)[3,1]+(confint(Staffing_MLM)[4,1]*mean(Staffing_MLM@frame$Skill))+(confint(Staffing_MLM)[7,1]*mean(Staffing_MLM@frame$Somatic))+(confint(Staffing_MLM)[8,1]*mean(Staffing_MLM@frame$Turnover)),2),
                         " - ",
                         round(confint(Staffing_MLM)[3,2]+(confint(Staffing_MLM)[4,2]*mean(Staffing_MLM@frame$Skill))+(confint(Staffing_MLM)[7,2]*mean(Staffing_MLM@frame$Somatic))+(confint(Staffing_MLM)[8,2]*mean(Staffing_MLM@frame$Turnover)),2),
                         ")")

## extract random effect for adjusted staffing variable
MRNPsych_Unit<-MRNPsych_Unit %>%
  left_join(.,as.data.frame(ranef(Staffing_MLM)), by = c("unit_code" = "grp")) %>%
  select(-c(grpvar,term,condsd)) %>%
  rename(RandEff = "condval") %>%
  mutate(Intercept = fixef(Staffing_MLM)[1],
         adjStaff = rowSums(across(c(Intercept,RandEff)))) %>%
  select(-c(Intercept,RandEff))

### create final data set ====
## merge nurse and unit data set together
MRNPsych_Merge<-MRNPsych_Nurses %>%
  left_join(., MRNPsych_Unit, by = c("unit_code")) %>%
  drop_na(TempFreq)

## create final data set
MRNPsych_Final<-MRNPsych_Merge %>%
  select(unit_code,
         JobSat,ITLJob,ITLProf,Burnout,
         Leadership,
         Shift,PNR,Skill,Somatic,Turnover,
         Sex,Age,Perc,
         TempFreq,adjStaff) %>%
  drop_na(JobSat,Burnout,ITLJob,ITLProf,
          Leadership,
          adjStaff,Age,Sex,Perc)

### descriptive analysis ====
## descriptive analyses nurse level
MRNPsych_Final %>%
  select(JobSat,Burnout,ITLJob,ITLProf,Leadership,Age,Sex,Perc) %>%
  tbl_summary(label = list(JobSat ~ "Job satisfaction",
                           ITLJob ~ "Intention to leave organization",
                           ITLProf ~ "Intention to leave profession",
                           Sex ~ "Gender",
                           Perc ~ "Employment percentage"),
              type = list(JobSat ~ 'continuous2',
                          ITLJob ~ 'continuous2',
                          ITLProf ~ 'continuous2',
                          Burnout ~ 'continuous2',
                          Leadership ~ 'continuous2',
                          Age ~ 'continuous2'),
              statistic = all_continuous2() ~ c("{mean} ({sd})",
                                                "{median} [{p25}, {p75}]",
                                                "{min}, {max}"),
              digits = all_continuous() ~ 1) %>%
  as_flex_table() %>% 
  flextable::save_as_docx(path = "overall_nurses.docx")

## descriptive analyses unit level
MRNPsych_Final %>%
  filter(!duplicated(unit_code)) %>%
  select(TempFreq,adjStaff) %>%
  tbl_summary(label = list(TempFreq ~ "Temporary nurses use",
                           adjStaff ~ "Adjusted staffing"),
              type = list(adjStaff ~ 'continuous2'),
              statistic = all_continuous2() ~ c("{mean} ({sd})",
                                                "{median} [{p25}, {p75}]",
                                                "{min}, {max}"),
              digits = all_continuous() ~ 1) %>%
  as_flex_table() %>% 
  flextable::save_as_docx(path = "overall_units.docx")

## descriptive analysis nurses of the adjusted staffing model
MRNPsych_Model %>%
  select(PNR,Skill,Shift,Somatic,Turnover) %>%
  tbl_summary(label = list(PNR ~ "Patient-to-nurse ratio",
                           Skill ~ "Skill % grade mix in %",
                           Shift ~ "Shift",
                           Somatic ~ "Somtic diagnosis ratio in %",
                           Turnover ~ "Turnover ratio in %"),
              type = list(PNR ~ 'continuous2',
                          Skill ~ 'continuous2',
                          Somatic ~ 'continuous2',
                          Turnover ~ 'continuous2'),
              statistic = all_continuous2() ~ c("{mean} ({sd})",
                                                "{median} [{p25}, {p75}]",
                                                "{min}, {max}"),
              digits = all_continuous() ~ 1) %>%
  as_flex_table() %>% 
  flextable::save_as_docx(path = "overall_staffing.docx")

### inferential analysis ====
## final multilevel analysis
# Job satisfaction
Final_MLM_JobSat_wo<-lmer(JobSat ~ Age + Sex + Perc + adjStaff + TempFreq + (1 | unit_code), data = MRNPsych_Final)

Final_MLM_JobSat_w<-lmer(JobSat ~ Age + Sex + Perc + Leadership + adjStaff + TempFreq + (1 | unit_code), data = MRNPsych_Final)

# Bournout
Final_MLM_Burnout_wo<-lmer(Burnout ~  Age + Sex + Perc + adjStaff + TempFreq + (1 | unit_code), data = MRNPsych_Final)

Final_MLM_Burnout_w<-lmer(Burnout ~  Age + Sex + Perc + Leadership + adjStaff + TempFreq + (1 | unit_code), data = MRNPsych_Final)

# Intention to leave
Final_MLM_ITLProf_wo<-lmer(ITLProf ~ Age + Sex + Perc + adjStaff + TempFreq + (1 | unit_code), data = MRNPsych_Final)
Final_MLM_ITLJob_wo<-lmer(ITLJob ~ Age + Sex + Perc + adjStaff + TempFreq + (1 | unit_code), data = MRNPsych_Final)

Final_MLM_ITLProf_w<-lmer(ITLProf ~ Age + Sex + Perc + Leadership + adjStaff + TempFreq + (1 | unit_code), data = MRNPsych_Final)
Final_MLM_ITLJob_w<-lmer(ITLJob ~ Age + Sex + Perc + Leadership + adjStaff + TempFreq + (1 | unit_code), data = MRNPsych_Final)

# overall final models
tab_model(Final_MLM_JobSat_wo,Final_MLM_Burnout_wo,Final_MLM_ITLJob_wo,Final_MLM_ITLProf_wo,
          #pred.labels = c("Intercept","Age","Sex [Male]","Employment percentage [61%-95%]","Employment percentage [96%-100%]","Leadership","Adjusted staffing","Frequent use of temporary nurses [Frequent]"),
          dv.labels = c("Job satisfaction","Burnout","Intention to leave organization","Intention to leave profession"),
          string.pred = "Coefficient",
          string.ci = "CI (95%)",
          p.style = "stars",
          digits = 3,
          file = "final_model_wo.doc")

tab_model(Final_MLM_JobSat_w,Final_MLM_Burnout_w,Final_MLM_ITLJob_w,Final_MLM_ITLProf_w,
          #pred.labels = c("Intercept","Age","Sex [Male]","Employment percentage [61%-95%]","Employment percentage [96%-100%]","Leadership","Adjusted staffing","Frequent use of temporary nurses [Frequent]"),
          dv.labels = c("Job satisfaction","Burnout","Intention to leave organization","Intention to leave profession"),
          string.pred = "Coefficient",
          string.ci = "CI (95%)",
          p.style = "stars",
          digits = 3,
          file = "final_model_w.doc")