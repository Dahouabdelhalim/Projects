#########################################
###   ANALYSING TRENDS IN SEX RATIO   ###
#########################################

library(MuMIn)

# Whole-island ------------------------------------------------------------

whole_island_sex_ratio <- read.csv("SI_ASR_total_island.csv")

#First assessing whether a linear or non-linear temporal trend is most appropriate

year1 <-glm(cbind(Num_males,Num_females)~poly(Year, 1), data=whole_island_sex_ratio, family = "binomial")

year2 <-glm(cbind(Num_males,Num_females)~poly(Year, 2), data=whole_island_sex_ratio, family = "binomial")

year3 <-glm(cbind(Num_males,Num_females)~poly(Year, 3), data=whole_island_sex_ratio, family = "binomial")

model.sel(year1, year2, year3)

#Now comparing trends in the different age categories

whole_island_m1 <-glm(cbind(Num_males,Num_females)~1, data=whole_island_sex_ratio, family = "binomial")

whole_island_m2 <-glm(cbind(Num_males,Num_females)~poly(Year, 1), data=whole_island_sex_ratio, family = "binomial")

whole_island_m3 <-glm(cbind(Num_males,Num_females)~poly(Year, 1) + Age_category, data=whole_island_sex_ratio, family = "binomial")

whole_island_m4 <-glm(cbind(Num_males,Num_females)~poly(Year, 1)*Age_category, data=whole_island_sex_ratio, family = "binomial")

model.sel(whole_island_m1, whole_island_m2, whole_island_m3, whole_island_m4)

summary(whole_island_m4)

# Subdivision level -------------------------------------------------------

subdivision_sex_ratio <- read.csv("SI_ASR_subdivisions.csv")

subdivision_m1 <-glm(cbind(Num_males,Num_females)~1, data=subdivision_sex_ratio, family = "binomial")

subdivision_m2 <-glm(cbind(Num_males,Num_females)~ Year, data=subdivision_sex_ratio, family = "binomial")

subdivision_m3 <-glm(cbind(Num_males,Num_females)~Year + Age_category, data=subdivision_sex_ratio, family = "binomial")

subdivision_m4 <-glm(cbind(Num_males,Num_females)~Year + Age_category + Subdivision, data=subdivision_sex_ratio, family = "binomial")

subdivision_m5 <-glm(cbind(Num_males,Num_females)~ Year*Age_category + Subdivision, data=subdivision_sex_ratio, family = "binomial")

subdivision_m6 <-glm(cbind(Num_males,Num_females)~Year + Age_category*Subdivision, data=subdivision_sex_ratio, family = "binomial")

subdivision_m7 <-glm(cbind(Num_males,Num_females)~Year*Age_category+Year*Subdivision, data=subdivision_sex_ratio, family = "binomial")

subdivision_m8 <-glm(cbind(Num_males,Num_females)~Year*Age_category+Age_category*Subdivision, data=subdivision_sex_ratio, family = "binomial")

subdivision_m9 <-glm(cbind(Num_males,Num_females)~Year*Age_category*Subdivision, data=subdivision_sex_ratio, family = "binomial")

model.sel(subdivision_m1, subdivision_m2, subdivision_m3, subdivision_m4, subdivision_m5, subdivision_m6, subdivision_m7, subdivision_m8, subdivision_m9)

summary(subdivision_m5)
