# C mineralization Note
# Data Analysis
# 1-30-20 CV Updated 5-22-20

# set working directory
setwd("C:/Users/carme/Desktop")

# load in data for comparing C mineralization through time across treatments
dry<-read.csv("cmin_drying.csv")
dry<-na.omit(dry)

# two-way ANCOVA
modeld<-aov(dry$CO2~dry$Treatment*dry$Day+dry$Moisture)
summary(modeld)
plot(dry$CO2~dry$Moisture)

# look at residuals and check assumptions
hist(resid(modeld))
shapiro.test(resid(modeld))
bartlett.test(dry$CO2~dry$Treatment)
boxplot(dry$CO2~dry$Treatment)

# load in data for Franzluebbers assay
franz<-read.csv("Franzluebbers.csv")

# run ANOVA
modelf<-aov(franz$CO2~franz$Treatment)
summary(modelf)

# look at residuals and check assumptions
hist(resid(modelf))
shapiro.test(resid(modelf))
bartlett.test(franz$CO2~franz$Treatment)
boxplot(franz$CO2~franz$Treatment)

# load in corrected data
corrected<-read.csv("corrected.csv")

# see if treatment effect disappears after we correct for it
modelc<-aov(corrected$Corrected_Cmin~corrected$Treatment)
summary(modelc)
hist(resid(modelc))
shapiro.test(resid(modelc))
bartlett.test(corrected$Corrected_Cmin~corrected$Treatment)
boxplot(corrected$Corrected_Cmin~corrected$Treatment)


