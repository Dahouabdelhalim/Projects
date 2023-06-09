# library
library(psych)
library(lme4)
library(lmerTest)
library(data.table)

####################
######### Sample 1 (TDL)
####################
# set the working directory
setwd("~/switchdrive/Boris/Recherche/Ambizione/Expe/2019/Eye tracker/Results_DPT_tobii/Sample 1/")

# load the eprime data
Data_dpt_sample1 <- read.csv("dpt_sample1_27_05_2019.csv",sep= ";")

# delete the rows associated with the training procedure
Data_dpt_sample1 <- Data_dpt_sample1[which((Data_dpt_sample1$Procedure!="PracticeProc")), ]

# delete the trials with wrong answer
Data_dpt_sample1 <- Data_dpt_sample1[which((Data_dpt_sample1$Stimulus.ACC!="0")), ]

## select reaction times from 200 to 1500 ms
limit_RT <- c(200,1500)
## indicate RT as numeric 
Data_dpt_sample1$RT <- as.numeric(as.character(Data_dpt_sample1$Stimulus.RT))
## selected RT as numeric 
Data_dpt_sample1$selected <- Data_dpt_sample1$RT > limit_RT [1] & Data_dpt_sample1$RT < limit_RT [2]
# Keep RT between 200-1500
Data_dpt_sample1 <- Data_dpt_sample1[Data_dpt_sample1$selected,]

## keep only useful variables
#keepVars <- c("subject", "Stimuli","Action", "RT", "images")
#data_main <- data_main [, keepVars]

# create stimuli category (PA versus SB) on the left versus rigth side of the screen
Data_dpt_sample1$stimuli_left <- NA
Data_dpt_sample1$stimuli_left[is.element(Data_dpt_sample1$Image5,c("s1.jpg", "s2.jpg", "s3.jpg",
                                                     "s4.jpg", "s5.jpg"))] <- "PA"

Data_dpt_sample1$stimuli_left[is.element(Data_dpt_sample1$Image5,c("l1.jpg", "l2.jpg", "l3.jpg",
                                                  "l4.jpg", "l5.jpg"))] <- "SB"

Data_dpt_sample1$stimuli_right <- NA
Data_dpt_sample1$stimuli_right[is.element(Data_dpt_sample1$Image6,c("s1.jpg", "s2.jpg", "s3.jpg",
                                                   "s4.jpg", "s5.jpg"))] <- "PA"

Data_dpt_sample1$stimuli_right[is.element(Data_dpt_sample1$Image6,c("l1.jpg", "l2.jpg", "l3.jpg",
                                                   "l4.jpg", "l5.jpg"))] <- "SB"

# create a variable indicating where thre dot appeared on the screen
Data_dpt_sample1$dot_side <- Data_dpt_sample1$Image8
Data_dpt_sample1$dot_side <- factor(Data_dpt_sample1$dot_side,
                            levels = c("white.jpg",
                                       "black_dot.jpg"),
                            labels = c("left", "right"))

# Create a variable indicating if the dot appears on the side of PA or SB
Data_dpt_sample1$dot_stim <- ifelse(Data_dpt_sample1$dot_side == "right",
                            Data_dpt_sample1$stimuli_right, Data_dpt_sample1$stimuli_left)
Data_dpt_sample1$dot_stim <- factor(Data_dpt_sample1$dot_stim,
                            levels = c("SB", "PA"),
                            labels = c("Dot on SB side",
                                       "Dot on PA side"))

Data_dpt_sample1$sample <- "Sample1"

# to create the right name for the subject (the "S" become "s")
Data_dpt_sample1$Subject
Data_dpt_sample1$subject <- as.numeric(as.character(Data_dpt_sample1$Subject))
Data_dpt_sample1$subject <- Data_dpt_sample1$subject + 38000

# recode correct value for the participants
Data_dpt_sample1$subject <- ifelse(Data_dpt_sample1$subject ==47533, 38533, Data_dpt_sample1$subject)
Data_dpt_sample1$subject
Data_dpt_sample1$subject <- as.factor(as.character(Data_dpt_sample1$subject))

# load  the sample 1 (TDL) self-report data
setwd("~/switchdrive/Boris/Recherche/Ambizione/Expe/2019/TDL_motivation/Results/")
load("data_all_SR.RData")

### create a dataset with only meaninful variables
data_SR_sample1 <- data.frame (subject =data_all_SR$subject, age = data_all_SR$age, BMI = data_all_SR$BMI,
                                         gender= data_all_SR$gender,
                                         Intention_PA_all =data_all_SR$Intention_PA_all,
                                         Intention_limit_SB_all=data_all_SR$Intention_limit_SB_all,
                                         IPAQ_cat = data_all_SR$IPAQ_cat,IPAQ_cat_2 = data_all_SR$IPAQ_cat_2,
                                         total_time_sed = data_all_SR$total_time_sed, 
                                         total_time_PA=data_all_SR$total_time_PA, total_time_MVPA_free_time=data_all_SR$total_time_MVPA_free_time, 
                                         Sit_free_time = data_all_SR$Sit_free_time,
                                         add_sev_all = data_all_SR$add_sev_all, 
                                         add_conti_all = data_all_SR$add_conti_all, 
                                         add_tol_all = data_all_SR$add_tol_all, 
                                         add_lack_cont_all = data_all_SR$add_lack_cont_all, 
                                         add_red_act_all = data_all_SR$add_red_act_all, 
                                         add_int_all = data_all_SR$add_int_all, 
                                         add_time_all = data_all_SR$add_time_all,
                                         add_total_all = data_all_SR$add_total_all )

# Merge behavioral and self-reported data
dpt_all_sample1 <- merge(Data_dpt_sample1, data_SR_sample1, by="subject", all.x =TRUE)

#test mixed models
m1_test <- lmer(RT   ~  1 + dot_stim + (1|subject) + (1|Image5) + (1|Image6), data=dpt_all_sample1,  REML=FALSE)
summary(m1_test)


##############
### Sample 2 (inhibition)
###############
# set the working directory
setwd("~/switchdrive/Boris/Recherche/Ambizione/Expe/2019/Eye tracker/Results_DPT_tobii/Sample 2/")

# load the eprime data
Data_dpt_sample2 <- read.csv("dpt_sample2_27_05_2019.csv",sep= ";")

# delete the rows associated with the training procedure
Data_dpt_sample2 <- Data_dpt_sample2[which((Data_dpt_sample2$Procedure!="PracticeProc")), ]

# delete the trials with wrong answer
Data_dpt_sample2 <- Data_dpt_sample2[which((Data_dpt_sample2$Stimulus.ACC!="0")), ]

## select reaction times from 200 to 1500 ms
limit_RT <- c(200,1500)
## indicate RT as numeric 
Data_dpt_sample2$RT <- as.numeric(as.character(Data_dpt_sample2$Stimulus.RT))
## selected RT as numeric 
Data_dpt_sample2$selected <- Data_dpt_sample2$RT > limit_RT [1] & Data_dpt_sample2$RT < limit_RT [2]
# Keep RT between 200-1500
Data_dpt_sample2 <- Data_dpt_sample2[Data_dpt_sample2$selected,]

## keep only useful variables
#keepVars <- c("subject", "Stimuli","Action", "RT", "images")
#data_main <- data_main [, keepVars]

# create stimuli category (PA versus SB) on the left versus rigth side of the screen
Data_dpt_sample2$stimuli_left <- NA
Data_dpt_sample2$stimuli_left[is.element(Data_dpt_sample2$Image5,c("s1.jpg", "s2.jpg", "s3.jpg",
                                                                   "s4.jpg", "s5.jpg"))] <- "PA"

Data_dpt_sample2$stimuli_left[is.element(Data_dpt_sample2$Image5,c("l1.jpg", "l2.jpg", "l3.jpg",
                                                                   "l4.jpg", "l5.jpg"))] <- "SB"

Data_dpt_sample2$stimuli_right <- NA
Data_dpt_sample2$stimuli_right[is.element(Data_dpt_sample2$Image6,c("s1.jpg", "s2.jpg", "s3.jpg",
                                                                    "s4.jpg", "s5.jpg"))] <- "PA"

Data_dpt_sample2$stimuli_right[is.element(Data_dpt_sample2$Image6,c("l1.jpg", "l2.jpg", "l3.jpg",
                                                                    "l4.jpg", "l5.jpg"))] <- "SB"

# create a variable indicating where thre dot appeared on the screen
Data_dpt_sample2$dot_side <- Data_dpt_sample2$Image8
Data_dpt_sample2$dot_side <- factor(Data_dpt_sample2$dot_side,
                                    levels = c("white.jpg",
                                               "black_dot.jpg"),
                                    labels = c("left", "right"))

# Create a variable indicating if the dot appears on the side of PA or SB
Data_dpt_sample2$dot_stim <- ifelse(Data_dpt_sample2$dot_side == "right",
                                    Data_dpt_sample2$stimuli_right, Data_dpt_sample2$stimuli_left)
Data_dpt_sample2$dot_stim <- factor(Data_dpt_sample2$dot_stim,
                                    levels = c("SB", "PA"),
                                    labels = c("Dot on SB side",
                                               "Dot on PA side"))
Data_dpt_sample2$sample <- "Sample2"

# to create the right name for the subject (the "S" become "s")
Data_dpt_sample2$subject <- NA
Data_dpt_sample2$subject <- as.numeric(as.character(Data_dpt_sample2$Subject))
Data_dpt_sample2$subject <- as.factor(as.character(Data_dpt_sample2$subject))

# for the sample 2 (Inhibition)
setwd("~/switchdrive/Boris/Recherche/Ambizione/Expe/2019/Inhibition/Results/")
load("data_SR1_inhib.RData")

### create a dataset with only meaninful variables
data_SR_sample2 <- data.frame (subject =data_SR1_inhib$subject, age = data_SR1_inhib$age, BMI = data_SR1_inhib$BMI,
                                         gender= data_SR1_inhib$gender,
                                         Intention_PA_all =data_SR1_inhib$Intention_PA_all,
                                         Intention_limit_SB_all=data_SR1_inhib$Intention_limit_SB_all,
                                         IPAQ_cat = data_SR1_inhib$IPAQ_cat,IPAQ_cat_2 = data_SR1_inhib$IPAQ_cat_2,
                                         total_time_sed = data_SR1_inhib$total_time_sed, 
                                         total_time_PA=data_SR1_inhib$total_time_PA, total_time_MVPA_free_time=data_SR1_inhib$total_time_MVPA_free_time, 
                                         Sit_free_time = data_SR1_inhib$Sit_free_time,
                                         add_sev_all = data_SR1_inhib$add_sev_all, 
                                         add_conti_all = data_SR1_inhib$add_conti_all, 
                                         add_tol_all = data_SR1_inhib$add_tol_all, 
                                         add_lack_cont_all = data_SR1_inhib$add_lack_cont_all, 
                                         add_red_act_all = data_SR1_inhib$add_red_act_all, 
                                         add_int_all = data_SR1_inhib$add_int_all, 
                                         add_time_all = data_SR1_inhib$add_time_all,
                                         add_total_all = data_SR1_inhib$add_total_all)

# Merge behavioral and self-reported data
dpt_all_sample2 <- merge(Data_dpt_sample2, data_SR_sample2, by="subject", all.x =TRUE)

#test mixed models
m2_test <- lmer(RT   ~  1 + dot_stim + (1|subject) + (1|Image5) + (1|Image6), data=dpt_all_sample2,  REML=FALSE)
summary(m2_test)

########
### Merge the two dataset together
#############
names(dpt_all_sample1)
names(dpt_all_sample2)
# check same columns
all.equal(names(dpt_all_sample1), names(dpt_all_sample2))
dpt_both_sample <- rbind(dpt_all_sample1, dpt_all_sample2)

### save the new data set
setwd("~/switchdrive/Boris/Recherche/Ambizione/Expe/2019/Eye tracker/")
save(dpt_both_sample, file="dpt_behavioral_both_sample.RData")

m3_test <- lmer(RT   ~  1 + dot_stim + (1|Subject) + (1|Image5) + (1|Image6), data=dpt_both_sample,  REML=FALSE)
summary(m3_test)

