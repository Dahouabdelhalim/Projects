#libraries
library(data.table)

############################
########### Sample 1 (TDL)
############################
# setwd
setwd("~/switchdrive/Boris/Recherche/Ambizione/Expe/2019/Eye tracker/Results_DPT_tobii/Sample 1/")

### Gaze data
### Read loop
gazeDat <- NULL
toRead <- list.files(pattern = "gazedata")
for (file in toRead) { 
  gazeDat <- rbind(gazeDat, fread(file = file))
}

gazeDat[, cue := {
  rows <- which(.SD$CurrentObject == "cue")
  cutpoints <- c(min(rows),
                 rows[which(diff(rows) > 1)],
                 max(rows))
  cueGroups <- cut(rows, cutpoints,
                   labels = 
                     paste0("cue",
                            c(1:(length(cutpoints) - 1))),
                   include.lowest = TRUE)
  ret <- factor(rep(NA_integer_, .N),
                levels = levels(cueGroups))
  ret[rows] <- cueGroups
  
  ret
}, .(Subject)]

gazeDat[, durationStim := c(diff(.SD$TETTime), NA),
        .(Subject)]

gazeDat[!is.na(cue), LookingAt := {
  ifelse(AOI == "centre", "cross",
         ifelse(AOI == "droite", MaskRight,
                ifelse(AOI == "gauche", MaskLeft, "nothing")))
}, .(Subject, cue)]

gazeDat[, LookingAt := {
  tmp <- gsub("s.\\\\.jpg", "Activity", LookingAt)
  tmp <- gsub("l.\\\\.jpg", "Sedentary", tmp)
  tmp
}]

# Time looking at stimuli
gazeDat[!is.na(cue), FirstLook := {
  LookingAt[AOI %in% c("gauche", "droite")][1]
}, .(Subject, cue)]

# Split area of interests (AOI) up within cues
gazeDat[!is.na(cue), AOISpell := {
  cutpoints <- which(diff(as.integer(factor(.SD$AOI))) != 0)
  cutpoints <- c(0, cutpoints, length(.SD$AOI))
  ret <- cut(c(1:length(.SD$AOI)),
             cutpoints,
             labels = paste0("Spell",
                             c(1:(length(cutpoints) - 1))),
             include.lowest = TRUE)
  levels(ret)[ret]
}, .(Subject, cue)]

gazeDat[!is.na(cue),
        timeInSpell := sum(.SD$durationStim),
        by = .(Subject, cue, LookingAt, AOISpell)]

gazeDat[!is.na(cue),
        totTimeStim := sum(.SD$durationStim),
        by = .(Subject, cue, LookingAt)]

gazeDat[!is.na(cue), 
        c("totTimeSed",
          "totTimeAct",
          "totTimeNeut") :=
          list(.SD[LookingAt == "Sedentary"]$totTimeStim[1],
               .SD[LookingAt == "Activity"]$totTimeStim[1],
               .SD[LookingAt == "cross"]$totTimeStim[1]),
        by = .(Subject, cue)]

### Pupil stuff

## Take the mean of pupil dilation for both eyes
## DiameterPupilLeft/RightEye
gazeDat[, DiameterPupilMean := rowMeans(
  gazeDat[, .(DiameterPupilLeftEye,
              DiameterPupilRightEye)])]

## Maximum public dilation per cue
gazeDat[!is.na(cue),
        MaxDiameterPupilCue := max(DiameterPupilMean),
        by = .(Subject, cue)]

gazeDat[!is.na(cue),
        c("MaxDiameterPupilAct",
          "MaxDiameterPupilSed",
          "MaxDiameterPupilNeut") :=
          list(if (.SD[LookingAt == "Sedentary", .N] > 0) {
            max(.SD[LookingAt == "Sedentary"]$DiameterPupilMean)
          } else {
            as.double(NA)
          },
          if (.SD[LookingAt == "Activity", .N] > 0) {
            max(.SD[LookingAt == "Activity"]$DiameterPupilMean)
          } else {
            as.double(NA)
          },
          if (.SD[LookingAt == "cross", .N] > 0) {
            max(.SD[LookingAt == "cross"]$DiameterPupilMean)
          } else {
            as.double(NA)
          }),
        by = .(Subject, cue)]

gazeDat[!is.na(cue),
        c("MedianDiameterPupilAct",
          "MedianDiameterPupilSed",
          "MedianDiameterPupilNeut") :=
          list(if (.SD[LookingAt == "Sedentary", .N] > 0) {
            median(.SD[LookingAt == "Sedentary"]$DiameterPupilMean)
          } else {
            as.double(NA)
          },
          if (.SD[LookingAt == "Activity", .N] > 0) {
            median(.SD[LookingAt == "Activity"]$DiameterPupilMean)
          } else {
            as.double(NA)
          },
          if (.SD[LookingAt == "cross", .N] > 0) {
            median(.SD[LookingAt == "cross"]$DiameterPupilMean)
          } else {
            as.double(NA)
          }),
        by = .(Subject, cue)]

gazeDat$sample <- "Sample1"

gazeDat$subject <- gazeDat$Subject

gazeDatCue <- gazeDat[, .SD[1, ], by = .(Subject, cue),
                         .SDcols = c("FirstLook", "MaskLeft", "MaskRight","CueLeft","CueRight","RT",
                                     "totTimeSed", "totTimeAct", "totTimeNeut",
                                     "DiameterPupilMean", "MaxDiameterPupilCue",
                                     "MaxDiameterPupilAct", "MaxDiameterPupilSed",
                                     "MaxDiameterPupilNeut", "MedianDiameterPupilAct",
                                     "MedianDiameterPupilSed", "MedianDiameterPupilNeut",
                                     "sample", "subject", "ACC")]

gazeDatCue <- gazeDatCue[!is.na(cue)]

# create the code for the behavioral performance 
# create stimuli category (PA versus SB) on the left versus rigth side of the screen
gazeDatCue$stimuli_left <- NA
gazeDatCue$stimuli_left[is.element(gazeDatCue$MaskLeft,c("s1.jpg", "s2.jpg", "s3.jpg",
                                                                   "s4.jpg", "s5.jpg"))] <- "PA"

gazeDatCue$stimuli_left[is.element(gazeDatCue$MaskLeft,c("l1.jpg", "l2.jpg", "l3.jpg",
                                                                   "l4.jpg", "l5.jpg"))] <- "SB"

gazeDatCue$stimuli_right <- NA
gazeDatCue$stimuli_right[is.element(gazeDatCue$MaskRight,c("s1.jpg", "s2.jpg", "s3.jpg",
                                                                    "s4.jpg", "s5.jpg"))] <- "PA"

gazeDatCue$stimuli_right[is.element(gazeDatCue$MaskRight,c("l1.jpg", "l2.jpg", "l3.jpg",
                                                                    "l4.jpg", "l5.jpg"))] <- "SB"

# create a variable indicating where thre dot appeared on the screen
gazeDatCue$dot_side <- gazeDatCue$CueRight
gazeDatCue$dot_side <- factor(gazeDatCue$dot_side,
                                    levels = c("white.jpg",
                                               "black_dot.jpg"),
                                    labels = c("left", "right"))

# Create a variable indicating if the dot appears on the side of PA or SB
gazeDatCue$dot_stim <- ifelse(gazeDatCue$dot_side == "right",
                              gazeDatCue$stimuli_right, gazeDatCue$stimuli_left)
gazeDatCue$dot_stim <- factor(gazeDatCue$dot_stim,
                                    levels = c("SB", "PA"),
                                    labels = c("Dot on SB side",
                                               "Dot on PA side"))

## Set total time looking at SED and ACT to 0 if NA and other isn't NA
gazeDatCue$totTimeAct <- ifelse(!is.na(gazeDatCue$totTimeSed) & is.na(gazeDatCue$totTimeAct),
                                0, gazeDatCue$totTimeAct)
gazeDatCue$totTimeSed <- ifelse(is.na(gazeDatCue$totTimeSed) & !is.na(gazeDatCue$totTimeAct),
                                0, gazeDatCue$totTimeSed)
gazeDatCue$DiameterPupilMean[gazeDatCue$DiameterPupilMean == -1] <- NA 

gazeDatCue$totTimeSum <- apply(gazeDatCue[, .(totTimeSed, totTimeAct, totTimeNeut)], MARGIN = 1,
                               sum, na.rm = TRUE)

# difference on time PA versus SB
gazeDatCue$actTimeDif_PASB <- gazeDatCue$totTimeAct - gazeDatCue$totTimeSed

# difference on time PA versus SB in comparison to time total
gazeDatCue$actTimeDif_PASB_Ratio <- (gazeDatCue$totTimeAct - gazeDatCue$totTimeSed)/gazeDatCue$totTimeSum

############
### Pupil
###########
describe(gazeDatCue$MaxDiameterPupilAct)
hist(gazeDatCue$MaxDiameterPupilAct)

describe(gazeDatCue$MaxDiameterPupilSed)
hist(gazeDatCue$MaxDiameterPupilSed)

describe(gazeDatCue$MaxDiameterPupilNeut)
hist(gazeDatCue$MaxDiameterPupilNeut)

## select pupil between 2 mm to 4.5 mm 
gazeDatCue$MaxDiameterPupilAct_clean <- ifelse(gazeDatCue$MaxDiameterPupilAct>2,gazeDatCue$MaxDiameterPupilAct, NA)
gazeDatCue$MaxDiameterPupilAct_clean <- ifelse(gazeDatCue$MaxDiameterPupilAct_clean<4,gazeDatCue$MaxDiameterPupilAct_clean, NA)
describe(gazeDatCue$MaxDiameterPupilAct_clean)

gazeDatCue$MaxDiameterPupilNeut_clean <- ifelse(gazeDatCue$MaxDiameterPupilNeut>2,gazeDatCue$MaxDiameterPupilNeut, NA)
gazeDatCue$MaxDiameterPupilNeut_clean <- ifelse(gazeDatCue$MaxDiameterPupilNeut_clean<4,gazeDatCue$MaxDiameterPupilNeut_clean, NA)
describe(gazeDatCue$MaxDiameterPupilNeut_clean)

gazeDatCue$MaxDiameterPupilSed_clean <- ifelse(gazeDatCue$MaxDiameterPupilSed>2,gazeDatCue$MaxDiameterPupilSed, NA)
gazeDatCue$MaxDiameterPupilSed_clean <- ifelse(gazeDatCue$MaxDiameterPupilSed_clean<4,gazeDatCue$MaxDiameterPupilSed_clean, NA)
describe(gazeDatCue$MaxDiameterPupilSed_clean)

# difference on max pupil PA versus SB
gazeDatCue$MaxPupilDif_PASB <- gazeDatCue$MaxDiameterPupilAct - gazeDatCue$MaxDiameterPupilSed

describe(gazeDatCue$MaxPupilDif_PASB)
hist(gazeDatCue$MaxPupilDif_PASB)

gazeDatCue$MaxPupilDif_PASB_clean <- ifelse(gazeDatCue$MaxPupilDif_PASB>-1,gazeDatCue$MaxPupilDif_PASB, NA)
gazeDatCue$MaxPupilDif_PASB_clean <- ifelse(gazeDatCue$MaxPupilDif_PASB_clean<1,gazeDatCue$MaxPupilDif_PASB_clean, NA)
describe(gazeDatCue$MaxPupilDif_PASB_clean)
hist(gazeDatCue$MaxPupilDif_PASB_clean)

# to create the right name for the subject
gazeDatCue$subject
gazeDatCue$subject <- as.numeric(as.character(gazeDatCue$subject))
gazeDatCue$subject <- gazeDatCue$subject + 38000

# recode correct value for the participants
gazeDatCue$subject <- ifelse(gazeDatCue$subject ==47533, 38533, gazeDatCue$subject)
gazeDatCue$subject
gazeDatCue$subject <- as.factor(as.character(gazeDatCue$subject))

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
data_SR_sample1$subject <- as.factor(as.character(data_SR_sample1$subject ))
dpt_all_sample1_EyeTrack <- merge(gazeDatCue, data_SR_sample1, by="subject", all.x =TRUE)
table(dpt_all_sample1_EyeTrack$subject)


############################
########### Sample 2 (Inhibition)
############################
# setwd
setwd("~/switchdrive/Boris/Recherche/Ambizione/Expe/2019/Eye tracker/Results_DPT_tobii/Sample 2/")

### Gaze data
### Read loop
gazeDat <- NULL
toRead <- list.files(pattern = "gazedata")
for (file in toRead) { 
  gazeDat <- rbind(gazeDat, fread(file = file))
}

gazeDat[, cue := {
  rows <- which(.SD$CurrentObject == "cue")
  cutpoints <- c(min(rows),
                 rows[which(diff(rows) > 1)],
                 max(rows))
  cueGroups <- cut(rows, cutpoints,
                   labels = 
                     paste0("cue",
                            c(1:(length(cutpoints) - 1))),
                   include.lowest = TRUE)
  ret <- factor(rep(NA_integer_, .N),
                levels = levels(cueGroups))
  ret[rows] <- cueGroups
  
  ret
}, .(Subject)]

gazeDat[, durationStim := c(diff(.SD$TETTime), NA),
        .(Subject)]

gazeDat[!is.na(cue), LookingAt := {
  ifelse(AOI == "centre", "cross",
         ifelse(AOI == "droite", MaskRight,
                ifelse(AOI == "gauche", MaskLeft, "nothing")))
}, .(Subject, cue)]

gazeDat[, LookingAt := {
  tmp <- gsub("s.\\\\.jpg", "Activity", LookingAt)
  tmp <- gsub("l.\\\\.jpg", "Sedentary", tmp)
  tmp
}]

# Time looking at stimuli
gazeDat[!is.na(cue), FirstLook := {
  LookingAt[AOI %in% c("gauche", "droite")][1]
}, .(Subject, cue)]

# Split area of interests (AOI) up within cues
gazeDat[!is.na(cue), AOISpell := {
  cutpoints <- which(diff(as.integer(factor(.SD$AOI))) != 0)
  cutpoints <- c(0, cutpoints, length(.SD$AOI))
  ret <- cut(c(1:length(.SD$AOI)),
             cutpoints,
             labels = paste0("Spell",
                             c(1:(length(cutpoints) - 1))),
             include.lowest = TRUE)
  levels(ret)[ret]
}, .(Subject, cue)]

gazeDat[!is.na(cue),
        timeInSpell := sum(.SD$durationStim),
        by = .(Subject, cue, LookingAt, AOISpell)]

gazeDat[!is.na(cue),
        totTimeStim := sum(.SD$durationStim),
        by = .(Subject, cue, LookingAt)]

gazeDat[!is.na(cue), 
        c("totTimeSed",
          "totTimeAct",
          "totTimeNeut") :=
          list(.SD[LookingAt == "Sedentary"]$totTimeStim[1],
               .SD[LookingAt == "Activity"]$totTimeStim[1],
               .SD[LookingAt == "cross"]$totTimeStim[1]),
        by = .(Subject, cue)]

### Pupil stuff
## Take the mean of pupil dilation for both eyes
## DiameterPupilLeft/RightEye
gazeDat[, DiameterPupilMean := rowMeans(
  gazeDat[, .(DiameterPupilLeftEye,
              DiameterPupilRightEye)])]

## Maximum public dilation per cue
gazeDat[!is.na(cue),
        MaxDiameterPupilCue := max(DiameterPupilMean),
        by = .(Subject, cue)]

gazeDat[!is.na(cue),
        c("MaxDiameterPupilAct",
          "MaxDiameterPupilSed",
          "MaxDiameterPupilNeut") :=
          list(if (.SD[LookingAt == "Sedentary", .N] > 0) {
            max(.SD[LookingAt == "Sedentary"]$DiameterPupilMean)
          } else {
            as.double(NA)
          },
          if (.SD[LookingAt == "Activity", .N] > 0) {
            max(.SD[LookingAt == "Activity"]$DiameterPupilMean)
          } else {
            as.double(NA)
          },
          if (.SD[LookingAt == "cross", .N] > 0) {
            max(.SD[LookingAt == "cross"]$DiameterPupilMean)
          } else {
            as.double(NA)
          }),
        by = .(Subject, cue)]

gazeDat[!is.na(cue),
        c("MedianDiameterPupilAct",
          "MedianDiameterPupilSed",
          "MedianDiameterPupilNeut") :=
          list(if (.SD[LookingAt == "Sedentary", .N] > 0) {
            median(.SD[LookingAt == "Sedentary"]$DiameterPupilMean)
          } else {
            as.double(NA)
          },
          if (.SD[LookingAt == "Activity", .N] > 0) {
            median(.SD[LookingAt == "Activity"]$DiameterPupilMean)
          } else {
            as.double(NA)
          },
          if (.SD[LookingAt == "cross", .N] > 0) {
            median(.SD[LookingAt == "cross"]$DiameterPupilMean)
          } else {
            as.double(NA)
          }),
        by = .(Subject, cue)]


gazeDat$sample <- "Sample2"

gazeDat$subject <- gazeDat$Subject

gazeDatCue <- gazeDat[, .SD[1, ], by = .(Subject, cue),
                      .SDcols = c("FirstLook", "MaskLeft", "MaskRight","CueLeft","CueRight", "RT",
                                  "totTimeSed", "totTimeAct", "totTimeNeut",
                                  "DiameterPupilMean", "MaxDiameterPupilCue",
                                  "MaxDiameterPupilAct", "MaxDiameterPupilSed",
                                  "MaxDiameterPupilNeut", "MedianDiameterPupilAct",
                                  "MedianDiameterPupilSed", "MedianDiameterPupilNeut",
                                  "sample", "subject", "ACC")]

gazeDatCue <- gazeDatCue[!is.na(cue)]

# create the code for the behavioral performance 
# create stimuli category (PA versus SB) on the left versus rigth side of the screen
gazeDatCue$stimuli_left <- NA
gazeDatCue$stimuli_left[is.element(gazeDatCue$MaskLeft,c("s1.jpg", "s2.jpg", "s3.jpg",
                                                         "s4.jpg", "s5.jpg"))] <- "PA"

gazeDatCue$stimuli_left[is.element(gazeDatCue$MaskLeft,c("l1.jpg", "l2.jpg", "l3.jpg",
                                                         "l4.jpg", "l5.jpg"))] <- "SB"

gazeDatCue$stimuli_right <- NA
gazeDatCue$stimuli_right[is.element(gazeDatCue$MaskRight,c("s1.jpg", "s2.jpg", "s3.jpg",
                                                           "s4.jpg", "s5.jpg"))] <- "PA"

gazeDatCue$stimuli_right[is.element(gazeDatCue$MaskRight,c("l1.jpg", "l2.jpg", "l3.jpg",
                                                           "l4.jpg", "l5.jpg"))] <- "SB"

# create a variable indicating where thre dot appeared on the screen
gazeDatCue$dot_side <- gazeDatCue$CueRight
gazeDatCue$dot_side <- factor(gazeDatCue$dot_side,
                              levels = c("white.jpg",
                                         "black_dot.jpg"),
                              labels = c("left", "right"))

# Create a variable indicating if the dot appears on the side of PA or SB
gazeDatCue$dot_stim <- ifelse(gazeDatCue$dot_side == "right",
                              gazeDatCue$stimuli_right, gazeDatCue$stimuli_left)
gazeDatCue$dot_stim <- factor(gazeDatCue$dot_stim,
                              levels = c("SB", "PA"),
                              labels = c("Dot on SB side",
                                         "Dot on PA side"))

## Set total time looking at SED and ACT to 0 if NA and other isn't NA
gazeDatCue$totTimeAct <- ifelse(!is.na(gazeDatCue$totTimeSed) & is.na(gazeDatCue$totTimeAct),
                                0, gazeDatCue$totTimeAct)
gazeDatCue$totTimeSed <- ifelse(is.na(gazeDatCue$totTimeSed) & !is.na(gazeDatCue$totTimeAct),
                                0, gazeDatCue$totTimeSed)
gazeDatCue$DiameterPupilMean[gazeDatCue$DiameterPupilMean == -1] <- NA 

gazeDatCue$totTimeSum <- apply(gazeDatCue[, .(totTimeSed, totTimeAct, totTimeNeut)], MARGIN = 1,
                               sum, na.rm = TRUE)

# difference on time PA versus SB
gazeDatCue$actTimeDif_PASB <- gazeDatCue$totTimeAct - gazeDatCue$totTimeSed

# difference on time PA versus SB in comparison to time total
gazeDatCue$actTimeDif_PASB_Ratio <- (gazeDatCue$totTimeAct - gazeDatCue$totTimeSed)/gazeDatCue$totTimeSum

############
### Pupil
###########
describe(gazeDatCue$MaxDiameterPupilAct)
hist(gazeDatCue$MaxDiameterPupilAct)

describe(gazeDatCue$MaxDiameterPupilSed)
hist(gazeDatCue$MaxDiameterPupilSed)

describe(gazeDatCue$MaxDiameterPupilNeut)
hist(gazeDatCue$MaxDiameterPupilNeut)

## select pupil between 2 mm to 4.5 mm 
gazeDatCue$MaxDiameterPupilAct_clean <- ifelse(gazeDatCue$MaxDiameterPupilAct>2,gazeDatCue$MaxDiameterPupilAct, NA)
gazeDatCue$MaxDiameterPupilAct_clean <- ifelse(gazeDatCue$MaxDiameterPupilAct_clean<4,gazeDatCue$MaxDiameterPupilAct_clean, NA)
describe(gazeDatCue$MaxDiameterPupilAct_clean)

gazeDatCue$MaxDiameterPupilNeut_clean <- ifelse(gazeDatCue$MaxDiameterPupilNeut>2,gazeDatCue$MaxDiameterPupilNeut, NA)
gazeDatCue$MaxDiameterPupilNeut_clean <- ifelse(gazeDatCue$MaxDiameterPupilNeut_clean<4,gazeDatCue$MaxDiameterPupilNeut_clean, NA)
describe(gazeDatCue$MaxDiameterPupilNeut_clean)

gazeDatCue$MaxDiameterPupilSed_clean <- ifelse(gazeDatCue$MaxDiameterPupilSed>2,gazeDatCue$MaxDiameterPupilSed, NA)
gazeDatCue$MaxDiameterPupilSed_clean <- ifelse(gazeDatCue$MaxDiameterPupilSed_clean<4,gazeDatCue$MaxDiameterPupilSed_clean, NA)
describe(gazeDatCue$MaxDiameterPupilSed_clean)

# difference on max pupil PA versus SB
gazeDatCue$MaxPupilDif_PASB <- gazeDatCue$MaxDiameterPupilAct - gazeDatCue$MaxDiameterPupilSed

describe(gazeDatCue$MaxPupilDif_PASB)
hist(gazeDatCue$MaxPupilDif_PASB)

gazeDatCue$MaxPupilDif_PASB_clean <- ifelse(gazeDatCue$MaxPupilDif_PASB>-1,gazeDatCue$MaxPupilDif_PASB, NA)
gazeDatCue$MaxPupilDif_PASB_clean <- ifelse(gazeDatCue$MaxPupilDif_PASB_clean<1,gazeDatCue$MaxPupilDif_PASB_clean, NA)
describe(gazeDatCue$MaxPupilDif_PASB_clean)
hist(gazeDatCue$MaxPupilDif_PASB_clean)


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
data_SR_sample2$subject <- as.factor(as.character(data_SR_sample2$subject))
gazeDatCue$subject <- as.factor(as.character(gazeDatCue$subject))
dpt_all_sample2_EyeTrack <- merge(gazeDatCue, data_SR_sample2, by="subject", all.x =TRUE)

# Merge the tow sample together 
dpt_both_sample_EyeTrack <- rbind(dpt_all_sample1_EyeTrack, dpt_all_sample2_EyeTrack)

### save the new data set
setwd("~/switchdrive/Boris/Recherche/Ambizione/Expe/2019/Eye tracker/")
save(dpt_both_sample_EyeTrack, file="dpt_both_sample_EyeTrack.RData")
