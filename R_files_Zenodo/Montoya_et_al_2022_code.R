#### MONTOYA ET AL 2022 ANALYSES & FIGURES ####


#libraries
library(asreml) #need a license for this 
library(ggplot2)
library(ggpubr) 
library(varhandle)
library(tidyr)
library(stringr)
library(lme4)
library(emmeans)
library(tibble)
library(plyr)
library(dplyr)
library(gridExtra)
library(grid)



#### IMPORT SYMBIOSIS PHENOTYPING DATA ####
# "full" Data has two observations per pot (for 2-strain inocula: each strain as focal (strain.o) and as competitor (strain.comp.o))
D1= read.csv("symbiosis_phenotype_data.csv")

# seperate metadata from data
meta.D1 = select(D1, "Meta.Data")
D1=select(D1,-"Meta.Data")
str(D1) #D1 has 1040 obs. of 25 variables

# make relevant columns into factors
D1$Unique.ID = as.factor(D1$Unique.ID)
D1$strain.label = as.factor(D1$strain.label)
D1$strain.o = as.factor(D1$strain.o)
D1$strain.comp.o = as.factor(D1$strain.comp.o)
D1$Block = as.factor(D1$Block)
D1$inoc.type = as.factor(D1$inoc.type)
D1$marker = as.factor(D1$marker)

#### Data processing for 1) Partner Choice and Sanctions

# add 1-strain mean shoot mass for focal strain
ss.means.strain.p = aggregate(shoot.mass ~ strain.o, FUN=mean, data=subset(D1, inoc.type=='single'))
names(ss.means.strain.p)[2] = 'shoot.mass.ss.o'
D1 = merge(D1, ss.means.strain.p, all.x=T)

# add 1-strain mean shoot mass for competitor strain
ss.means.comp.strain.p = ss.means.strain.p
names(ss.means.comp.strain.p) = c('strain.comp.o', 'shoot.mass.ss.comp')
D1 = merge(D1, ss.means.comp.strain.p, all.x=T)

# Nodule number for 1-strain inocula, mean and se
ss.means.strain.nod = aggregate((nod.sum+ambiguous.nods) ~ strain.o, FUN=mean, data=subset(D1, inoc.type=='single'))
ss.means.strain.nod$nod.se = aggregate((nod.sum+ambiguous.nods) ~ strain.o, FUN=function(x) sd(x)/sqrt(length(x)), data=subset(D1, inoc.type=='single'))[,2]
colnames(ss.means.strain.nod)=c("strain.o", "nod.sum", "nod.sum.se")

# add mean nodule count to 1-strain inocula mean shoot mass df
ss.means = merge(ss.means.strain.p, ss.means.strain.nod, by='strain.o')

# add nodule proportion to ss.means data frame
ss.means$nod.prop = aggregate(nod.sum/(nod.sum+comp.nod.sum) ~ strain.o, FUN=mean, data=subset(D1, inoc.type=='mix'))[,2]
ss.means$nod.prop.se = aggregate(nod.sum ~ strain.o, FUN=function(x) sd(x)/sqrt(length(x)), data=subset(D1, inoc.type=='single'))[,2]

# separate 2-strain inocula, 1-strain inocula, and Neg  (D1 now has 1040 obs. of 27 variables)
D1m = D1[which(D1$inoc.type=="mix"),] # 896 obs. of  27 variables
str(D1m)
D1s = D1[which(D1$inoc.type=="single"),] # 128 obs. of  27 variables
str(D1s)
D1n = D1[which(D1$inoc.type=="NEG"),] # 16 obs. of  27 variables #need this for neg controls on plots
str(D1n)

# create new nodule proportion column for mixed strains that includes 0s 
D1m$nod.prop = with(D1m, nod.sum/nods.tot)
D1m$comp.nod.prop = with(D1m, 1-nod.prop)
# check NAs for nodule proportion (ie any with plants with no nodules in mixed inoc)
D1m[which(is.na(D1m$nod.prop)),]
D1m[which(is.na(D1m$comp.nod.prop)),]

# make 1-strain inocula and NEG nod.prop & nod.comp.prop = NA
D1s$nod.prop = NA
D1s$comp.nod.prop = NA
D1n$nod.prop = NA
D1n$comp.nod.prop = NA


# add mixed back to single and neg (all observations)
D1 = rbind(D1m, D1s, D1n) # 1040 obs. of  29 variables 
str(D1) 


#### Data processing for 2) Absolute and Conditional Components of Partner Choice

# use D1 dataset (symbiosis phenotype data (ie "full" dataset))
# create "half" data set from "full" data set
# Odd blocks
do=rbind((D1m[which(D1m$Block=="B01"),]),
         (D1m[which(D1m$Block=="B03"),]),
         (D1m[which(D1m$Block=="B05"),]),
         (D1m[which(D1m$Block=="B07"),]))

# odd blocks focal strain is red
do2=rbind(do[which(do$marker == "mScarlet"),])

# Even blocks
de=rbind((D1m[which(D1m$Block=="B02"),]),
         (D1m[which(D1m$Block=="B04"),]),
         (D1m[which(D1m$Block=="B06"),]),
         (D1m[which(D1m$Block=="B08"),]))

# even blocks focal strain is green 
de2=rbind(de[which(de$marker == "sfGFP"),])

# combine to make "half" dataset
D1.half=rbind(do2,de2)
str(D1.half) #D1.half has 448 obs. of 29 variables


#### Data processing for 3) Benefits of Partner choice

#make a copy of D1m here to use in the supplemental analysis 
SD1m = D1m

# use D1m dataset (symbiosis phenotype data for 2-strain inocula)
# add a column to ID more beneficial strain in a pot (used for checks)
D1m$strain.best = apply(D1m, 1, function(x) ifelse(x['shoot.mass.ss.o'] > x['shoot.mass.ss.comp'], x['strain.o'], x['strain.comp.o'])) 
D1m$strain.worst = apply(D1m, 1, function(x) ifelse(x['shoot.mass.ss.o'] < x['shoot.mass.ss.comp'], x['strain.o'], x['strain.comp.o']))

# Calculate the Deviation from a neutral expectation (dne) metric: (W_s1s2 - mean(W_s1 + W_s2))/mean(W_s1 + W_s2) 
# (ie the neutral expectation for host fitness in the absence of partner choice)

# calculate the neutral expectation of shoot mass under no partner choice (ie the mean shoot mass of the 2 strains in 1-strain inocula)
D1m$MeanShoot2Strains = (D1m$shoot.mass.ss.comp + D1m$shoot.mass.ss.o)/2 #used for check
# calculate the difference between actual shoot mass in 2-strain inocula and expected shoot mass
D1m$Dif = D1m$shoot.mass - ((D1m$shoot.mass.ss.comp + D1m$shoot.mass.ss.o)/2) #used for check
# calculate the shoot mass deviation from the neutral expectation (use in model)
D1m$dne = (D1m$shoot.mass - ((D1m$shoot.mass.ss.comp + D1m$shoot.mass.ss.o)/2))/((D1m$shoot.mass.ss.comp + D1m$shoot.mass.ss.o)/2)

# add new column "NodPropBestStrain" to assign nod prop of more beneficial strain
D1m$NodPropBestStrain = apply(D1m, 1, function(x) ifelse(x['shoot.mass.ss.o'] > x['shoot.mass.ss.comp'], x['nod.prop'], x['comp.nod.prop'])) 
head(D1m)

# Check data
str(D1m) #D1m has 896 obs. of  35 variables
# change to numeric
D1m$NodPropBestStrain = as.numeric(D1m$NodPropBestStrain)
# droplevels (no neg)
D1m=droplevels(D1m)
levels(D1m$strain.o) #good
levels(D1m$strain.comp.o) #good

# Use only one observation (row) per plant pot from D1m
D1m2 = D1m[!duplicated(D1m[, "Unique.ID"]),] 
str(D1m2) # D1m2 has 448 obs. of  35 variables

#Drop rows where nodpropbeststrain=na, so models are fitted to same data set despite na (ex. na for nodpropbeststrain but not for marker-results in diff data sets))
D1m2=D1m2[complete.cases(D1m2[,35]),] #35 is column for nodpropbeststrain 


#### IMPORT CFU DATA ####
DC = read.csv("cfu_data.csv", header = T)

# separate metadata from data
meta.cfu = select(DC, "Meta.Data")
DC = select(DC, -"Meta.Data")

str(DC)
#change to factor 
DC$C.conc = as.factor(DC$C.conc)
DC$E.conc = as.factor(DC$E.conc)
DC$F.conc = as.factor(DC$F.conc)
DC$G.conc = as.factor(DC$G.conc)
DC$H.conc = as.factor(DC$H.conc)
DC$I.conc = as.factor(DC$I.conc)

#### Convert saturated assays ("s") to 9999

# replace saturated assays ("s") with 9999-- easily recognizable number
# this lets us "unfactor" the count data so that they can all be treated as numbers

# C dilution
# calculate num, the number of levels in this column
# the last level in this column is probably "s", but we will check this
num = DC$C.conc %>% levels() %>% length()
levels(DC$C.conc)
levels(DC$C.conc)[num]
# "s"
levels(DC$C.conc)[num] = 9999
# converts to numeric
DC$C.conc = unfactor(DC$C.conc)

# E dilution
# calculate num, the number of levels in this column
# the last level in this column is probably "s", but we will check this
num = DC$E.conc %>% levels() %>% length()
levels(DC$E.conc)
levels(DC$E.conc)[num]
# "s"
levels(DC$E.conc)[num] = 9999
# converts to numeric
DC$E.conc = unfactor(DC$E.conc)

# F dilution
# calculate num, the number of levels in this column
# the last level in this column is probably "s", but we will check this
num = DC$F.conc %>% levels() %>% length()
levels(DC$F.conc)
levels(DC$F.conc)[num]
# "s"
levels(DC$F.conc)[num] = 9999
# converts to numeric
DC$F.conc = unfactor(DC$F.conc)

# G dilution
# calculate num, the number of levels in this column
# the last level in this column is probably "s", but we will check this
num = DC$G.conc %>% levels() %>% length()
levels(DC$G.conc)
levels(DC$G.conc)[num]
# "s"
levels(DC$G.conc)[num] = 9999
# converts to numeric
DC$G.conc = unfactor(DC$G.conc)

# H dilution
# calculate num, the number of levels in this column
# the last level in this column is probably "s", but we will check this
num = DC$H.conc %>% levels() %>% length()
levels(DC$H.conc)
levels(DC$H.conc)[num]
# "s"
levels(DC$H.conc)[num] = 9999
# converts to numeric
DC$H.conc = unfactor(DC$H.conc)

# I dilution
# calculate num, the number of levels in this column
# the last level in this column is probably "s", but we will check this
num = DC$I.conc %>% levels() %>% length()
levels(DC$I.conc)
levels(DC$I.conc)[num]
# "s"
levels(DC$I.conc)[num] = 9999
# converts to numeric
DC$I.conc = unfactor(DC$I.conc)


#### Calculate CFU from non-saturated assays 

# make vector of dilution factors corresponding to dilutions C through I
df = c(125, 3125, 15625, 78125, 390625, 1953125)


# calculate cfu for each dilution level in each replicate
# do not use counts of 9999 (i.e., saturated)
# do not use counts outside 1-50 range
DC$C.cfu = ifelse(DC$C.conc>50,"",
                   ifelse(DC$C.conc<1,"",df[1]*DC$C.conc))
DC$C.cfu = as.numeric(DC$C.cfu)

DC$E.cfu = ifelse(DC$E.conc>50,"",
                   ifelse(DC$E.conc<1,"",df[2]*DC$E.conc))
DC$E.cfu = as.numeric(DC$E.cfu)

DC$F.cfu = ifelse(DC$F.conc>50,"",
                   ifelse(DC$F.conc<1,"",df[3]*DC$F.conc))
DC$F.cfu = as.numeric(DC$F.cfu)

DC$G.cfu = ifelse(DC$G.conc>50,"",
                   ifelse(DC$G.conc<1,"",df[4]*DC$G.conc))
DC$G.cfu = as.numeric(DC$G.cfu)

DC$H.cfu = ifelse(DC$H.conc>50,"",
                   ifelse(DC$H.conc<1,"",df[5]*DC$H.conc))
DC$H.cfu = as.numeric(DC$H.cfu)

DC$I.cfu = ifelse(DC$I.conc>50,"",
                   ifelse(DC$I.conc<1,"",df[5]*DC$I.conc))
DC$I.cfu = as.numeric(DC$I.cfu)


# add up cfu across all six dilutions within each replicate
DC$cfu.without.s = rowSums(DC[,c("C.cfu", "E.cfu", "F.cfu", "G.cfu", "H.cfu", "I.cfu")], na.rm = TRUE)

# count up how many dilutions contributed useable data to cfu.without.s within each replicate
DC$C.num = ifelse(is.na(DC$C.cfu),0,1)
DC$E.num = ifelse(is.na(DC$E.cfu),0,1)
DC$F.num = ifelse(is.na(DC$F.cfu),0,1)
DC$G.num = ifelse(is.na(DC$G.cfu),0,1)
DC$H.num = ifelse(is.na(DC$H.cfu),0,1)
DC$I.num = ifelse(is.na(DC$I.cfu),0,1)
DC$number.without.s = DC$C.num + DC$E.num + DC$F.num + DC$G.num + DC$H.num + DC$I.num


# remove some of the processing columns that we don't need anymore
DC = select(DC, -c("C.cfu", "E.cfu", "F.cfu", "G.cfu", "H.cfu", "I.cfu",
                    "C.num", "E.num", "F.num", "G.num", "H.num", "I.num"))


#### Calculate CFU from saturated assays 

# we will only use saturated assays for replicates where the only data was zero or saturated
# then, we will assign a count of 1 to the highest dilution in which growth was saturated

# flag replicates in which count data has only values of 0 and/or 9999 and/or NA (i.e., no growth or saturated)
# column will have value TRUE for these rows of interest
DC$Flag = (DC$C.conc==0 | DC$C.conc==9999 | is.na(DC$C.conc)) &
  (DC$E.conc==0 | DC$E.conc==9999 | is.na(DC$E.conc)) &
  (DC$F.conc==0 | DC$F.conc==9999 | is.na(DC$F.conc)) &
  (DC$G.conc==0 | DC$G.conc==9999 | is.na(DC$G.conc)) &
  (DC$H.conc==0 | DC$H.conc==9999 | is.na(DC$H.conc)) &
  (DC$I.conc==0 | DC$I.conc==9999 | is.na(DC$I.conc))

# for rows where Flag = TRUE, if highest dilution (I) is saturated, assign a cfu value of 1*df for I
# if dilution I is not saturated, check dilution H, and so on
DC$cfu.using.s = ifelse(DC$Flag==FALSE,"",
                         ifelse(DC$I.conc==9999,df[6],
                                ifelse(DC$H.conc==9999,df[5],
                                       ifelse(DC$G.conc==9999,df[4],
                                              ifelse(DC$F.conc==9999,df[3],
                                                     ifelse(DC$E.conc==9999,df[2],
                                                            ifelse(DC$C.conc==9999,df[1],"")))))))

DC$cfu.using.s = as.numeric(DC$cfu.using.s)      


# count up how many dilutions contributed useable data to cfu.using.s within each replicate
# basically, just record 1 for each row where cfu.using.s has a value (is not NA)
DC$number.using.s = ifelse(is.na(DC$cfu.using.s),"",1)
DC$number.using.s = as.numeric(DC$number.using.s)   



#### Add up CFU counts from saturated and non-saturated assays, within each nodule and colony color 
D2 = DC %>%
  group_by(Unique.ID, Unique.nod.ID, Nodule.Type, Nodule.Desc, Drop.plate.colony.color, marker) %>%
  summarize(cfu.without.s = sum(cfu.without.s, na.rm = TRUE),
            cfu.using.s = sum(cfu.using.s, na.rm = TRUE),
            number.without.s = sum(number.without.s, na.rm = TRUE),
            number.using.s = sum(number.using.s, na.rm = TRUE))

#### Average the CFU counts within each nodule and colony color
D2$Mean.CFU = (D2$cfu.without.s + D2$cfu.using.s) / 
  (D2$number.without.s + D2$number.using.s)
D2$Mean.CFU = as.numeric(D2$Mean.CFU)
# some rows have a value of NaN (not a number) because of dividing by zero; set these to NA
D2$Mean.CFU[is.nan(D2$Mean.CFU)] = NA

# remove columns we don't need anymore
D2 = select(D2, -c("cfu.without.s", "cfu.using.s", "number.without.s", "number.using.s"))
head(D2)
# Spread out data in wide format so that there is one row per nodule, with columns for red and green CFU
D2 = pivot_wider(D2, names_from = "Drop.plate.colony.color", values_from = "Mean.CFU")
head(D2)
names(D2)
names(D2)[6:7] = c("Green.CFU", "Red.CFU")

D2 = as.data.frame(D2) 

# Merge processed cfu data with D1 (symbiosis phenotype dataset)
cfu=merge(D2, D1, by= c("Unique.ID", "marker"), all.x=T, all.y=F)

# Make new cfu column per row (ie green colonies assigned to sfGFP marked strain and vice versa)
cfu$CFU=ifelse(cfu$marker == "sfGFP", cfu$Green.CFU, cfu$Red.CFU)
head(cfu)

# Make cfu.total column for each unique nod ID (will be used to calculate cfu proportion in mixed nodules)
cfu.total=aggregate(CFU ~ Unique.nod.ID, data = cfu, FUN = sum)
colnames(cfu.total)[2]= "CFU.total"

cfu=merge(cfu, cfu.total, by="Unique.nod.ID", all.x=T, all.y=F)



#### 1. PARTNER CHOICE & SANCTIONS ####

#### **Partner choice ####
# Do more beneficial strains: 1) form more nodules in 1-strain inocula 
M1=lm(nod.sum ~ shoot.mass.ss.o, data = ss.means)

summary(M1)
#Coefficients:
#                Estimate Std. Error t value Pr(>|t|)   
#(Intercept)       -7.182      9.605  -0.748   0.4829   
#shoot.mass.ss.o 1096.466    239.742   4.574   0.0038 **
#Residual standard error: 7.897 on 6 degrees of freedom
#Multiple R-squared:  0.7771,	Adjusted R-squared:  0.7399 
#F-statistic: 20.92 on 1 and 6 DF,  p-value: 0.003796

# plot residuals 
par(mfrow=c(2,2))
plot(M1)

# test effect of shoot mass in 1-strain inocula on nodule number in 1-strain inocula
drop1(M1, test = "F")
#                Df Sum of Sq     RSS    AIC F value   Pr(>F)   
#<none>                        374.19 34.763                    
#shoot.mass.ss.o  1    1304.5 1678.68 44.771  20.917 0.003796 **

# More beneficial strains form more nodules in 1-strain inocula



# Do more beneficial strains: 2) form a higher proportion of nodules in 2-strain inocula
M2=lm(nod.prop ~ shoot.mass.ss.o, data = ss.means)

summary(M2)
#Coefficients:
#                Estimate Std. Error t value Pr(>|t|)    
#(Intercept)      -0.3677     0.1279  -2.876 0.028194 *  
#shoot.mass.ss.o  22.4296     3.1911   7.029 0.000414 ***
#Residual standard error: 0.1051 on 6 degrees of freedom
#Multiple R-squared:  0.8917,	Adjusted R-squared:  0.8737 
#F-statistic:  49.4 on 1 and 6 DF,  p-value: 0.0004142

#plot residuals
par(mfrow=c(2,2))
plot(M2)

# test effect of shoot mass in 1-strain inocula on nodule proportion in 2-strain inocula
drop1(M2, test = "F")
#                Df Sum of Sq     RSS     AIC F value    Pr(>F)    
#<none>                       0.06630 -34.345                      
#shoot.mass.ss.o  1   0.54587 0.61217 -18.562  49.404 0.0004142 ***

# More beneficial strains form a greater proportion of nodules in 2-strain inocula



#### **Sanctions ####

# Subset data for mixed inoc and nodule.type "full" nodules (nodules occupied fully by one strain in 2-strain inocula)
cfu.full = subset(cfu, inoc.type=='mix' & Nodule.Desc=="full")

# Genotypic means for "full" nodules
cfu.means.full = aggregate(CFU ~ strain.o, FUN=mean, data=cfu.full)
cfu.means.full$se = aggregate(CFU ~ strain.o, FUN=function(x) sd(x)/sqrt(length(x)), data=subset(cfu.full))[,2]
names(cfu.means.full)[2] = 'mean.cfu'
names(cfu.means.full)[3] = 'cfu.se'

# combine CFU "full" nodules with single strain inoc shoot mass means 
cfu.full.shoot.means=cbind(ss.means, cfu.means.full)
cfu.full.shoot.means=cfu.full.shoot.means[,-7]

# Do more beneficial strains: 1) form a greater number of CFUs in nodules containing a single strain (i.e., among nodule sanctions)
M3=lm(mean.cfu ~ shoot.mass.ss.o, data = cfu.full.shoot.means)

summary(M3)
#Coefficients:
#                Estimate Std. Error t value Pr(>|t|)  
#(Intercept)      -526080     641421   -0.82   0.4435  
#shoot.mass.ss.o 46752641   16009152    2.92   0.0266 *
#Residual standard error: 527300 on 6 degrees of freedom
#Multiple R-squared:  0.587,	Adjusted R-squared:  0.5182 
#F-statistic: 8.529 on 1 and 6 DF,  p-value: 0.02662


# plot residuals
par(mfrow=c(2,2))
plot(M3)

# test effect of shoot mass in 1-strain inocula on cfu per nodule in 2-strain inocula
drop1(M3, test = "F")
#                Df  Sum of Sq        RSS    AIC F value  Pr(>F)  
#<none>                        1.6685e+12 212.51                  
#shoot.mass.ss.o  1 2.3717e+12 4.0403e+12 217.58  8.5286 0.02662 *

# more beneficial strains form greater cfu's per nodule in 2-strain inocula


# Subset cfu data for mixed inoc and "mixed" nodules
cfu.mixed=subset(cfu, inoc.type=='mix' & Nodule.Desc=='mixed')

# add cfu proportion within mixed nodules to data frame
cfu.mixed$cfu.prop=cfu.mixed$CFU/cfu.mixed$CFU.total

# now drop any cfu.prop=1 because a strain didn't grow 
# only used data for mixed nodules where both strains grew
cfu.mixed$cfu.prop=ifelse(cfu.mixed$cfu.prop==1, NA, cfu.mixed$cfu.prop)

# Genotypic means of cfu proportion 
cfu.means.mixed.prop = aggregate(cfu.prop ~ strain.o, FUN=mean, data=cfu.mixed)
cfu.means.mixed.prop$se = aggregate(cfu.prop ~ strain.o, FUN=function(x) sd(x)/sqrt(length(x)), data=cfu.mixed)[,2]
names(cfu.means.mixed.prop)[2] = 'mean.cfu.prop'
names(cfu.means.mixed.prop)[3] = 'cfu.prop.se'

# combine cfu.prop means and shoot.mass means 
cfu.prop.mixed.shoot.means=cbind(ss.means,cfu.means.mixed.prop)
cfu.prop.mixed.shoot.means=cfu.prop.mixed.shoot.means[,-7]

# Do more beneficial strains: 2) form a greater proportion of CFUs in nodules containing two strains (i.e., intra-nodule sanctions).
M4=lm(mean.cfu.prop ~ shoot.mass.ss.o, data = cfu.prop.mixed.shoot.means)

summary(M4)
#Coefficients:
#                Estimate Std. Error t value Pr(>|t|)  
#(Intercept)      0.02679    0.12694   0.211   0.8398  
#shoot.mass.ss.o 10.56877    3.16834   3.336   0.0157 *
#Residual standard error: 0.1044 on 6 degrees of freedom
#Multiple R-squared:  0.6497,	Adjusted R-squared:  0.5913 
#F-statistic: 11.13 on 1 and 6 DF,  p-value: 0.01569

# plot residuals
par(mfrow=c(2,2))
plot(M4)

# test effect of shoot mass in 1-strain inocula on cfu proportion in mixed nodules
drop1(M4, test = "F")
#                Df Sum of Sq      RSS     AIC F value  Pr(>F)  
#<none>                       0.065353 -34.459                  
#shoot.mass.ss.o  1    0.1212 0.186552 -28.068  11.127 0.01569 *

# more beneficial strains form a greater proportion of cfu's in mixed nodules in 2-strain inocula


#### 2. ABSOLUTE & CONDITIONAL COMPONENTS OF PARTNER CHOICE ####
# use D1 "full" dataset for nodule number
# used D1.half "half" dataset for nodule proportion and shoot mass

#### **Nodule number ####

# model nodule number use D1 with plant pot as random effect
M.nod = asreml(fixed= nod.sum ~ marker,
                random=~ Unique.ID + strain.comp.o:strain.o + strain.o + strain.comp.o + Block,
                residual = ~idv(units),
                data = subset(D1, inoc.type=='mix'),
                na.action = na.method(x="omit", y="omit"))

M.nod$converge  #TRUE

# check model assumptions 
plot(M.nod) # slight heterogeneity of variance

# summary of fixed effect
summary(M.nod,coef = TRUE)$coef.fixed
#                solution std error   z.ratio
#marker_mScarlet  0.00000        NA        NA
#marker_sfGFP    -1.06250 0.7902569 -1.344499
#(Intercept)     18.63728 4.8145416  3.871039

# summary of random effects
summary(M.nod)$varcomp
#                       component std.error   z.ratio bound %ch
#strain.o               118.37067 66.148706  1.789463     P   0
#strain.comp.o           50.04776 29.644278  1.688277     P   0
#Block                   10.29737  6.425759  1.602514     P   0
#strain.comp.o:strain.o  26.25581  8.036613  3.267024     P   0
#Unique.ID               26.48986  8.248706  3.211396     P   0
#units!units            139.88966  9.668900 14.468002     P   0
#units!R                  1.00000        NA        NA     F   0

# test fixed effects
wald.asreml(M.nod) 
#              Df Sum of Sq Wald statistic Pr(Chisq)    
#(Intercept)    1   14.2387        14.2387  0.000161 ***
#marker         1    1.8077         1.8077  0.178787 # marker not significant
#residual (MS)       1.0000 

# Models for LRTs

# model GxG removed
M.nodb = asreml(fixed= nod.sum ~ marker,
                 random=~ Unique.ID + strain.o + strain.comp.o  + Block,
                 residual = ~idv(units),
                 data = subset(D1, inoc.type=='mix'),
                 na.action = na.method(x="omit", y="omit"))

# model DGE removed
M.nodc = asreml(fixed= nod.sum ~ marker,
                 random=~Unique.ID + strain.comp.o:strain.o + strain.comp.o  + Block,
                 residual = ~idv(units),
                 data = subset(D1, inoc.type=='mix'),
                 na.action = na.method(x="omit", y="omit"))

# model main SGE removed
M.nodd = asreml(fixed= nod.sum ~ marker,
                 random=~ Unique.ID + strain.comp.o:strain.o + strain.o  + Block,
                 residual = ~idv(units),
                 data = subset(D1, inoc.type=='mix'),
                 na.action = na.method(x="omit", y="omit"))

# model Block removed
M.node = asreml(fixed= nod.sum ~ marker,
                 random=~ Unique.ID + strain.comp.o:strain.o + strain.o + strain.comp.o,
                 residual = ~idv(units),
                 data = subset(D1, inoc.type=='mix'),
                 na.action = na.method(x="omit", y="omit"))

# model plant pot removed
M.nodf = asreml(fixed= nod.sum ~marker,
                 random=~ strain.comp.o:strain.o + strain.o + strain.comp.o + Block,
                 residual = ~idv(units),
                 data = subset(D1, inoc.type=='mix'),
                 na.action = na.method(x="omit", y="omit"))

# Use LRT to calculate chi-squared stat and p-value

# GxG SGE (strain.o:strain.comp.o)
lrt(M.nod, M.nodb)
#               df LR-statistic Pr(Chisq)    
#M.nod/M.nodb  1       49.536 9.739e-13 ***

# DGE (strain.o)
lrt(M.nod, M.nodc)
#               df LR-statistic Pr(Chisq)    
#M.nod/M.nodc  1       49.119 1.204e-12 ***

# main SGE (strain.comp.o)
lrt(M.nod, M.nodd)
#               df LR-statistic Pr(Chisq)    
#M.nod/M.nodd  1       25.605 2.095e-07 ***

# Block
lrt(M.nod, M.node)
#               df LR-statistic Pr(Chisq)    
#M.nod/M.node  1       26.316  1.45e-07 ***
  
# Plant pot (Unique.ID)
lrt(M.nod, M.nodf)
#               df LR-statistic Pr(Chisq)    
#M.nod/M.nodf  1       10.726 0.0005282 ***

# change variance ratios into variance components

# DGE (strain.o)
dge = vpredict(M.nod, hA ~ V1/(V1 + V2 + V3 + V4 + V5 + V6))
rownames(dge) = "DGE"
dge

# Main SGE (strain.comp.o)
main = vpredict(M.nod, hA ~ V2/(V1 + V2 + V3 + V4 + V5 + V6))
rownames(main) = "Main SGE"
main

# GxG SGE (strain.o:strain.comp.o) 
gxg = vpredict(M.nod, hA ~ V4/(V1 + V2 + V3 + V4 + V5 + V6))
rownames(gxg) = "GxG SGE"
gxg

# Block 
block = vpredict(M.nod, hA ~ V3/(V1 + V2 + V3 + V4 + V5 + V6))
rownames(block) = "Block"
block

# Plant pot (Unique.ID)
pot = vpredict(M.nod, hA ~ V5/(V1 + V2 + V3 + V4 + V5 +V6))
rownames(pot) = "Plant pot"
pot

# residual variance 
resid = vpredict(M.nod, hA ~ V6/(V1 + V2 + V3 + V4 + V5+ V6))
rownames(resid) = "residual"
resid

# combine variance dataframes
tnod = rbind(dge, main, gxg, block, pot, resid)
# change estimates and SE to percent and round to 1 decimal place
tnod$percentVar = round(tnod$Estimate*100, 1)
tnod$percentSE = round(tnod$SE*100, 1)
tnod

#            Estimate         SE percentVar percentSE
#DGE       0.31875674 0.12476714       31.9      12.5
#Main SGE  0.13477206 0.07358025       13.5       7.4
#GxG SGE   0.07070345 0.02483226        7.1       2.5
#Block     0.02772946 0.01769027        2.8       1.8
#Plant pot 0.07133373 0.02586277        7.1       2.6
#residual  0.37670456 0.07712636       37.7       7.7

#### **Nodule proportion ####
# model nodule proportion use D1.half dataset and variance constraints

# only strain.o and strain.comp.o constrained to be equal (Vfocal=Vcompetitor) with correlation=-1
# initial 'run' to generate starting parameter file so we can set constraints
M.nod.prop = asreml(fixed= nod.prop ~ marker,
                  random=~ strain.comp.o:strain.o + Block + str(~strain.o + strain.comp.o,
                                                               ~corgh(2, init = c(0.1,0.1,0.1)):id(strain.o)), #init specifies the correlation then the variance, id() needs levels specified but struggling here)
                  residual = ~idv(units),
                  data = D1.half,
                  maxiter = 20,
                  equate.levels=c('strain.o','strain.comp.o'),
                  na.action = na.method(x="omit", y="omit"),
                  start.values = TRUE)

# now edit table to impose constraints needed  
M.nod.prop$vparameters.table

constraints = M.nod.prop$vparameters.table
constraints[3,2] = -0.9999     #replace foc-comp correlation with (almost) -1. 
constraints[3,3] = "F"         #fix to this starting value, so when you run the model the correlation is imposed not estimated 

# note third component in table is the correlation which we have already fixed to -1.
# set the constraint the fourth (Focal strain) and fifth (Competitor strain) variance components are equal 
M = as.matrix(data.frame(V1=c(1,2,3,4,4,5,6), V2=c(1,1,1,1,1,1,1)))  
dimnames(M)[[1]] = constraints$Component 

# now actually run model, fitting to data to estimate variance with constraints 
# use G.param=constraints and vcc=M 
M.nod.prop = asreml(fixed= nod.prop ~ marker,
                  random=~strain.comp.o:strain.o + Block + str(~strain.o + strain.comp.o,
                                                               ~corgh(2, init = c(0.1,0.1,0.1)):id(strain.o)), #init specifies the correlation then the variance, id() specifies level
                  residual = ~idv(units),
                  data = D1.half,
                  maxiter = 20,
                  equate.levels=c('strain.o','strain.comp.o'),
                  na.action = na.method(x="omit", y="omit"),
                  G.param=constraints,
                  vcc=M)

M.nod.prop$converge  #True

# check model assumptions 
plot(M.nod.prop)

# summary of fixed effect
summary(M.nod.prop,coef = TRUE)$coef.fixed
#                    solution  std error    z.ratio
#marker_mScarlet  0.000000000         NA         NA
#marker_sfGFP    -0.005560554 0.01786858 -0.3111917
#(Intercept)      0.493257297 0.01642989 30.0219544

# summary of random effects
summary(M.nod.prop)$varcomp
#                                                      component    std.error   z.ratio bound %ch
#strain.comp.o:strain.o                             0.0061097931 0.0016068344  3.802379     P   0
#Block                                              0.0003989788 0.0003684938  1.082729     P   0
#strain.o+strain.comp.o!corgh(2)!2:!corgh(2)!1.cor -0.9999000000           NA        NA     F   0
#strain.o+strain.comp.o!corgh(2)_1                  0.0599270710 0.0322287915  1.859427     P   0
#strain.o+strain.comp.o!corgh(2)_2                  0.0599270710 0.0322287915  1.859427     C   0
#units!units                                        0.0131205873 0.0009574000 13.704394     P   0
#units!R                                            1.0000000000           NA        NA     F   0

# test fixed effects
wald.asreml(M.nod.prop) 
#              Df Sum of Sq Wald statistic Pr(Chisq)    
#(Intercept)    1    1261.4         1261.4    <2e-16 ***
#marker         1       0.1            0.1    0.7557 # marker not significant
#residual (MS)          1.0  


# Models for LRTs

# initial model GxG removed
M.nod.propb = asreml(fixed= nod.prop ~ marker,
                   random=~ Block + str(~strain.o + strain.comp.o,
                                        ~corgh(2, init = c(0.1,0.1,0.1)):id(strain.o)), #init specifies the correlation then the variance, id() needs levels specified but struggling here)
                   residual = ~idv(units),
                   data = D1.half,
                   maxiter = 20,
                   equate.levels=c('strain.o','strain.comp.o'),
                   na.action = na.method(x="omit", y="omit"),
                   start.values = TRUE)


# now edit table to impose constraints needed  
M.nod.propb$vparameters.table

constraints = M.nod.propb$vparameters.table
constraints[2,2] = -0.9999     #replace foc-comp correlation with (almost) -1. 
constraints[2,3] = "F"         #fix to this starting value, so when you run the model the correlation is imposed not estimated 

M = as.matrix(data.frame(V1=c(1,2,3,3,4,5), V2=c(1,1,1,1,1,1)))  
dimnames(M)[[1]] = constraints$Component 

# now actually run model, fitting to data to estimate variance with constraints 
# use G.param=constraints and vcc=M 
M.nod.propb = asreml(fixed= nod.prop ~ marker,
                   random=~ Block + str(~strain.o + strain.comp.o,
                                        ~corgh(2, init = c(0.1,0.1,0.1)):id(strain.o)), #init specifies the correlation then the variance, id() needs levels specified but struggling here)
                   residual = ~idv(units),
                   data = D1.half,
                   maxiter = 20,
                   equate.levels=c('strain.o','strain.comp.o'),
                   na.action = na.method(x="omit", y="omit"),
                   G.param=constraints,
                   vcc=M)

 
# initial model Block removed
M.nod.propc = asreml(fixed= nod.prop ~ marker,
                   random=~ strain.o:strain.comp.o + str(~strain.o + strain.comp.o,
                                                         ~corgh(2, init = c(0.1,0.1,0.1)):id(strain.o)), #init specifies the correlation then the variance, id() needs levels specified but struggling here)
                   residual = ~idv(units),
                   data = D1.half,
                   maxiter = 20,
                   equate.levels=c('strain.o','strain.comp.o'),
                   na.action = na.method(x="omit", y="omit"),
                   start.values = TRUE)


# now edit table to impose constraints needed  
M.nod.propc$vparameters.table

constraints = M.nod.propc$vparameters.table
constraints[2,2] = -0.9999     #replace foc-comp correlation with (almost) -1. 
constraints[2,3] = "F"         #fix to this starting value, so when you run the model the correlation is imposed not estimated 

M = as.matrix(data.frame(V1=c(1,2,3,3,4,5), V2=c(1,1,1,1,1,1)))  
dimnames(M)[[1]] = constraints$Component 

# now actually run model, fitting to data to estimate variance with constraints 
# use G.param=constraints and vcc=M 
M.nod.propc = asreml(fixed= nod.prop ~ marker,
                   random=~ strain.o:strain.comp.o + str(~strain.o + strain.comp.o,
                                                         ~corgh(2, init = c(0.1,0.1,0.1)):id(strain.o)), #init specifies the correlation then the variance, id() needs levels specified but struggling here)
                   residual = ~idv(units),
                   data = D1.half,
                   maxiter = 20,
                   equate.levels=c('strain.o','strain.comp.o'),
                   na.action = na.method(x="omit", y="omit"),
                   G.param=constraints,
                   vcc=M)


# to test significance of DGE (strain.o) and main SGE (strain.comp.o) when they are contrained to be equal, remove both terms and their correlation structure
# model strain.o + strain.comp.o and correlation structure removed
M.nod.propd = asreml(fixed= nod.prop ~ marker,
                   random=~ strain.o:strain.comp.o + Block,  
                   residual = ~idv(units),
                   data = D1.half,
                   maxiter = 20,
                   na.action = na.method(x="omit", y="omit"))


#Use LRT to calculate chi-squared stat and p-value

# GxG SGE (strain.o:strain.comp.o)
lrt(M.nod.prop, M.nod.propb)
#                   df LR-statistic Pr(Chisq)    
#M.nod.prop/M.nod.propb  1       69.587 < 2.2e-16 ***

# Block
lrt(M.nod.prop, M.nod.propc)
#                   df LR-statistic Pr(Chisq)  
#M.nod.prop/M.nod.propc  1       4.0383   0.02224 *

# DGE (strain.o) and main SGE(strain.comp.o)
lrt(M.nod.prop, M.nod.propd)
#                   df LR-statistic Pr(Chisq)    
#M.nod.prop/M.nod.propd  1       120.58 < 2.2e-16 ***


# change variance ratios into variance components
# DGE (strain.o)
dge = vpredict(M.nod.prop, hA ~ V4/(V1 + V2 + V4 + V5 + V6))
rownames(dge) = "DGE"
dge

# Main SGE (strain.comp.o)
main = vpredict(M.nod.prop, hA ~ V5/(V1 + V2 + V4 + V5 + V6))
rownames(main) = "Main SGE"
main

# GxG SGE (strain.o:strain.comp.o) 
gxg = vpredict(M.nod.prop, hA ~ V1/(V1 + V2 + V4 + V5 + V6))
rownames(gxg) = "GxG SGE"
gxg

# Block 
block = vpredict(M.nod.prop, hA ~ V2/(V1 + V2 + V4 + V5 + V6))
rownames(block) = "Block"
block

# residual variance 
resid = vpredict(M.nod.prop, hA ~ V6/(V1 + V2 + V4 + V5 + V6))
rownames(resid) = "residual"
resid

# combine variance dataframes
tnodprop = rbind(dge, main, gxg, block, resid)
# change estimates and SE to percent and round to 1 decimal place
tnodprop$percentVar = round(tnodprop$Estimate*100, 1)
tnodprop$percentSE = round(tnodprop$SE*100, 1)
tnodprop

#            Estimate          SE percentVar percentSE
#DGE      0.429635552 0.033018462       43.0       3.3
#Main SGE 0.429635552 0.033018462       43.0       3.3
#GxG SGE  0.043802981 0.023090143        4.4       2.3
#Block    0.002860402 0.002948279        0.3       0.3
#residual 0.094065514 0.043936335        9.4       4.4



#### **Shoot mass ####
# model shoot mass use D1.half data set and variance constraints

# only strain.o and strain.comp.o constrained to be equal (Vfocal=Vcompetitor) with correlation=+1
# initial 'run' to generate starting parameter file so we can set constraints 
M.shoot = asreml(fixed= log(shoot.mass) ~ marker,
                  random=~ strain.comp.o:strain.o + Block + str(~strain.o + strain.comp.o,
                                                               ~corgh(2, init = c(0.1,0.1,0.1)):id(strain.o)),  
                  residual = ~idv(units),
                  data = D1.half,
                  maxiter = 20,
                  equate.levels=c('strain.o','strain.comp.o'),
                  na.action = na.method(x="omit", y="omit"),
                  start.values = TRUE)

# now edit table to impose constraints needed  
M.shoot$vparameters.table

constraints = M.shoot$vparameters.table
constraints[3,2] = 0.9999     #replace foc-comp correlation with (almost) 1. 
constraints[3,3] = "F"         #fix to this starting value, so when you run the model the correlation is imposed not estimated 

# note third component in table is the correlation which we have already fixed to +1.
# set the constraint the fourth (Focal strain) and fifth (Competitor strain) variance components are equal 
M = as.matrix(data.frame(V1=c(1,2,3,4,4,5,6), V2=c(1,1,1,1,1,1,1)))  
dimnames(M)[[1]] = constraints$Component 

# now actually run model, fitting to data to estimate variance with constraints 
# use G.param=constraints and vcc=M 
M.shoot = asreml(fixed= log(shoot.mass) ~ marker,
                  random=~ strain.comp.o:strain.o + Block + str(~strain.o + strain.comp.o,
                                                               ~corgh(2, init = c(0.1,0.1,0.1)):id(strain.o)), #init specifies the correlation then the variance, id() specifies levels
                  residual = ~idv(units),
                  data = D1.half,
                  maxiter = 20,
                  equate.levels=c('strain.o','strain.comp.o'),
                  na.action = na.method(x="omit", y="omit"),
                  G.param=constraints,
                  vcc=M)
M.shoot$converge  #True

# check model assumptions 
plot(M.shoot)

#summary of fixed effect
summary(M.shoot,coef = TRUE)$coef.fixed
#                   solution std error     z.ratio
#marker_mScarlet  0.00000000        NA          NA
#marker_sfGFP     0.08809501 0.1561780   0.5640682
#(Intercept)     -3.33786006 0.1359665 -24.5491328

# summary of random effects
summary(M.shoot)$varcomp
#                                                   component   std.error   z.ratio bound %ch
#strain.comp.o:strain.o                            0.02484806 0.013145260  1.890268     P   0
#Block                                             0.04333663 0.028167676  1.538523     P   0
#strain.o+strain.comp.o!corgh(2)!2:!corgh(2)!1.cor 0.99990000          NA        NA     F   0
#strain.o+strain.comp.o!corgh(2)_1                 0.01169539 0.009120089  1.282376     P   0
#strain.o+strain.comp.o!corgh(2)_2                 0.01169539 0.009120089  1.282376     C   0
#units!units                                       0.30500301 0.021983075 13.874447     P   0
#units!R                                           1.00000000          NA        NA     F   0

# test fixed effects
wald.asreml(M.shoot) 
#              Df Sum of Sq Wald statistic Pr(Chisq)    
#(Intercept)    1    875.71         875.71    <2e-16 ***
#marker         1      0.32           0.32    0.5727 # marker not significant  
#residual (MS)         1.00 



# Models for LRTs

# initial model GxG removed
M.shootb = asreml(fixed= log(shoot.mass) ~ marker,
                   random=~ Block + str(~strain.o + strain.comp.o,
                                        ~corgh(2, init = c(0.1,0.1,0.1)):id(strain.o)), 
                   residual = ~idv(units),
                   data = D1.half,
                   maxiter = 20,
                   equate.levels=c('strain.o','strain.comp.o'),
                   na.action = na.method(x="omit", y="omit"),
                   start.values = TRUE)


# now edit table to impose constraints needed  
M.shootb$vparameters.table

constraints = M.shootb$vparameters.table
constraints[2,2] = 0.9999     #replace foc-comp correlation with (almost) -1. 
constraints[2,3] = "F"         #fix to this starting value, so when you run the model the correlation is imposed not estimated 

M = as.matrix(data.frame(V1=c(1,2,3,3,4,5), V2=c(1,1,1,1,1,1)))  
dimnames(M)[[1]] = constraints$Component 

# now actually run model, fitting to data to estimate variance with constraints 
# use G.param=constraints and vcc=M 
M.shootb = asreml(fixed= log(shoot.mass) ~ marker,
                   random=~ Block + str(~strain.o + strain.comp.o,
                                        ~corgh(2, init = c(0.1,0.1,0.1)):id(strain.o)), #init specifies the correlation then the variance, id() needs levels specified but struggling here)
                   residual = ~idv(units),
                   data = D1.half,
                   maxiter = 20,
                   equate.levels=c('strain.o','strain.comp.o'),
                   na.action = na.method(x="omit", y="omit"),
                   G.param=constraints,
                   vcc=M)


# initial model Block removed
M.shootc = asreml(fixed= log(shoot.mass) ~ marker,
                   random=~ strain.o:strain.comp.o + str(~strain.o + strain.comp.o,
                                                         ~corgh(2, init = c(0.1,0.1,0.1)):id(strain.o)), #init specifies the correlation then the variance, id() needs levels specified but struggling here)
                   residual = ~idv(units),
                   data = D1.half,
                   maxiter = 20,
                   equate.levels=c('strain.o','strain.comp.o'),
                   na.action = na.method(x="omit", y="omit"),
                   start.values = TRUE)


# now edit table to impose constraints needed  
M.shootc$vparameters.table

constraints = M.shootc$vparameters.table
constraints[2,2] = 0.9999     #replace foc-comp correlation with (almost) -1. 
constraints[2,3] = "F"         #fix to this starting value, so when you run the model the correlation is imposed not estimated 

M = as.matrix(data.frame(V1=c(1,2,3,3,4,5), V2=c(1,1,1,1,1,1)))  
dimnames(M)[[1]] = constraints$Component 

# now actually run model, fitting to data to estimate variance with constraints 
# use G.param=constraints and vcc=M 
M.shootc = asreml(fixed= log(shoot.mass) ~ marker,
                   random=~ strain.o:strain.comp.o + str(~strain.o + strain.comp.o,
                                                         ~corgh(2, init = c(0.1,0.1,0.1)):id(strain.o)), #init specifies the correlation then the variance, id() needs levels specified but struggling here)
                   residual = ~idv(units),
                   data = D1.half,
                   maxiter = 20,
                   equate.levels=c('strain.o','strain.comp.o'),
                   na.action = na.method(x="omit", y="omit"),
                   G.param=constraints,
                   vcc=M)

# to test significance of DGE (strain.o) and main SGE (strain.comp.o) when they are contrained to be equal, remove both terms and their correlation structure
# model strain.o + strain.comp.o with correlation structure removed
M.shootd = asreml(fixed= log(shoot.mass) ~ marker,
                   random=~ strain.o:strain.comp.o + Block,  
                   residual = ~idv(units),
                   data = D1.half,
                   maxiter = 20,
                   na.action = na.method(x="omit", y="omit"))


# Use LRT to calculate chi-squared stat and p-value

# GxG SGE (strain.o:strain.comp.o)
lrt(M.shoot, M.shootb)
#                   df LR-statistic Pr(Chisq)   
#M.shoot/M.shootb  1       6.1168  0.006695 **

# Block
lrt(M.shoot, M.shootc)
#                   df LR-statistic Pr(Chisq)    
#M.shoot/M.shootc  1       31.889 8.161e-09 ***
  
# DGE (strain.o) and main SGE (strain.comp.o)
lrt(M.shoot, M.shootd)
#                   df LR-statistic Pr(Chisq)   
#M.shoot/M.shootd  1       5.5283  0.009356 **


# change variance ratios into variance components
# DGE (strain.o) 
dge = vpredict(M.shoot, hA ~ V4/(V1 + V2 + V4 + V5 + V6))
rownames(dge)="DGE"
dge

# Main SGE (strain.comp.o)
main = vpredict(M.shoot, hA ~ V5/(V1 + V2 + V4 + V5 + V6))
rownames(main) = "Main SGE"
main

# GxG SGE (strain.o:strain.comp.o) 
gxg = vpredict(M.shoot, hA ~ V1/(V1 + V2 + V4 + V5 + V6))
rownames(gxg) = "GxG SGE"
gxg

# Block 
block = vpredict(M.shoot, hA ~ V2/(V1 + V2 + V4 + V5 + V6))
rownames(block) = "Block"
block

# residual variance 
resid = vpredict(M.shoot, hA ~ V6/(V1 + V2 + V4 + V5 + V6))
rownames(resid) = "residual"
resid

# combine variance dataframes
tshoot = rbind(dge, main, gxg, block, resid)
# change estimates and SE to percent and round to 1 decimal place
tshoot$percentVar = round(tshoot$Estimate*100, 1)
tshoot$percentSE = round(tshoot$SE*100, 1)
tshoot

#           Estimate         SE percentVar percentSE
#DGE      0.02949073 0.02191896        2.9       2.2
#Main SGE 0.02949073 0.02191896        2.9       2.2
#GxG SGE  0.06265609 0.03271044        6.3       3.3
#Block    0.10927630 0.06382434       10.9       6.4
#residual 0.76908616 0.07075047       76.9       7.1



#### 3. BENEFIT OF PARTNER CHOICE ####
# use D1m2 dataset (processed symbiosis phenotype data for 2-strain inocula with deviation from neutral expecation (dne) calculation)

# Benefit of partner choice to the host model
M5 = lmer(dne ~ NodPropBestStrain + marker + (1|strain.o) + (1|strain.comp.o) + (1|Block), data=D1m2, na.action=na.exclude)

summary(M5)
#Random effects:
# Groups        Name        Variance Std.Dev.
# Block         (Intercept) 0.059026 0.2430  
# strain.comp.o (Intercept) 0.005898 0.0768  
# strain.o      (Intercept) 0.000000 0.0000  # singular fit (addressed below)
# Residual                  0.426170 0.6528  
#Fixed effects:
#                  Estimate Std. Error t value
#(Intercept)       -0.16852    0.14504  -1.162
#NodPropBestStrain  0.43702    0.12952   3.374
#markersfGFP        0.02371    0.06235   0.380

# check model assumptions
plot(M5) #heterogeneity of variance (addressed below)
qqnorm(resid(M5))
qqline(resid(M5))

# test fixed effects
drop1(M5, test = "Chisq") 
#                  npar    AIC     LRT  Pr(Chi)   
#<none>                 903.14                    
#NodPropBestStrain    1 911.48 10.3439 0.001299 **
#marker               1 901.28  0.1453 0.703029 # marker not significant

# partner choice benefits host 


# Reran model M5 with log(dne) to improve residuals and singular fit and overall result unchanged
M5.1=lmer(log(dne +1) ~  NodPropBestStrain + marker + (1|strain.o) + (1|strain.comp.o) + (1|Block), data=D1m2, na.action=na.omit)

summary(M5.1) #fit is no longer singular
#Random effects:
# Groups        Name        Variance  Std.Dev.
# Block         (Intercept) 0.0380225 0.19499 
# strain.comp.o (Intercept) 0.0032447 0.05696 
# strain.o      (Intercept) 0.0003029 0.01741 
# Residual                  0.3228418 0.56819 
#Fixed effects:
#                  Estimate Std. Error t value
#(Intercept)       -0.30858    0.12261  -2.517
#NodPropBestStrain  0.40128    0.11343   3.538
#markersfGFP       -0.02255    0.05427  -0.416

# check model assumptions
plot(M5.1) # improves residuals
qqnorm(resid(M5.1))
qqline(resid(M5.1))

# test fixed effects
drop1(M5.1, test = "Chisq")
#                  npar    AIC     LRT   Pr(Chi)    
#<none>                 779.65                      
#NodPropBestStrain    1 788.69 11.0471 0.0008883 ***
#marker               1 777.82  0.1718 0.6785172 # marker not significant

#Increasing homogeneity of variance in this model using a log transformation on the response variable did not impact the results 


#### MAIN FIGURES #### 

#** FIG 1: Partner choice ####
#a,b,c are nodule photos
#***** d. Shoot mass in 1-strain inocula (ie indicates how beneficial each strain is to the host) ####
#Shoot mass 1-strain inoc standard error (mean calculated above for models) **probably move back above but won't use til here
ss.means.strain.p$se = aggregate(shoot.mass ~ strain.o, FUN=function(x) sd(x)/sqrt(length(x)), data=subset(D1, inoc.type=='single'))[,2]
names(ss.means.strain.p)[3] = 'shoot.se'

#Reorder strains: least to most beneficial (based on ss.inoc shoot mass)
ss.means.strain.p$strain.o = factor(ss.means.strain.p$strain.o, levels = ss.means.strain.p[order(ss.means.strain.p$shoot.mass.ss.o), 'strain.o'])

#plot
ggplot(ss.means.strain.p, aes(x=strain.o, y=shoot.mass.ss.o)) +
  geom_bar(stat='identity') + 
  geom_errorbar(aes(ymin=shoot.mass.ss.o-shoot.se, ymax = shoot.mass.ss.o+shoot.se), width=0, size=1) + 
  geom_hline(aes(yintercept = mean(D1n$shoot.mass) + sd(D1n$shoot.mass)/sqrt(length(D1n$shoot.mass))), linetype='dashed', size=1) +
  geom_hline(aes(yintercept = mean(D1n$shoot.mass) - sd(D1n$shoot.mass)/sqrt(length(D1n$shoot.mass))), linetype='dashed', size=1) + 
  xlab('Rhizobium strain') + ylab('Shoot mass (g)') +
  theme(panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.border = element_blank(),axis.line = element_line(colour = "black"),axis.ticks = element_line(colour = "black", size = 1.5),
        axis.ticks.length = unit(0.2, "cm"),text = element_text(size=32), axis.text.x = element_text(size=28, vjust=0.5, angle=90), 
        axis.text.y = element_text(size=28,hjust=0.95), plot.title = element_text(hjust = 0.5))


#***** e. Nodule number by shoot mass in 1-strain inocula ####

#plot of single strain nodule number predicted by shoot mass 
ggplot(ss.means, aes(x=shoot.mass.ss.o, y=nod.sum)) + geom_point(size=6) + stat_smooth(method='lm', size=3, color='black') + 
  theme_bw() + xlab('Shoot mass (g) (1-strain)') + ylab('Nodule number (1-strain)')  +
  theme(panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),axis.line = element_line(colour = "black"),axis.ticks = element_line(colour = "black", size = 1.5), axis.ticks.length = unit(0.25, "cm"),text = element_text(size=32), axis.text.x = element_text(size=28, vjust=0.95), axis.text.y = element_text(size=28,hjust=0.95)) +
  stat_cor(label.y = 60, label.x=0.02,size=10, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))


#***** f. Nodule proportion in 2-strain inocula by shoot mass in 1-strain inocula ####

# plot of nodule proportion predicted by shoot mass 
ggplot(ss.means, aes(x=shoot.mass.ss.o, y=nod.prop)) + geom_point(size=6) + stat_smooth(method='lm', size=3, color='black') + 
  xlab('Shoot mass (g) (1-strain)') + ylab('Nodule proportion (2-strains)') + ylim(-.1,1) +
  theme(panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),axis.line = element_line(colour = "black"),axis.ticks = element_line(colour = "black", size = 1.5), axis.ticks.length = unit(0.25, "cm"),text = element_text(size=32), axis.text.x = element_text(size=28, vjust=0.95), axis.text.y = element_text(size=28,hjust=0.95)) +
  stat_cor(label.y = 1, label.x=0.0205,size=10, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))



#** FIG 2: Absolute & Conditional components of partner choice ####

#*** Percent of variance explained by DGE, Main SGE, GxG SGE ####
#***** a. Nodule number ####
# use tnod table

# Keep rows DGE, Main SGE, GxG SGE 
fnod = tnod[c("DGE", "Main SGE", "GxG SGE"),]

# change row names to a column
fnod = tibble::rownames_to_column(fnod, "component")

# reorder components
fnod$component = factor(fnod$component, levels=c("DGE", "Main SGE", "GxG SGE"))

#plot
fig2a = ggplot(fnod, aes(x=component, y=percentVar)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=percentVar-percentSE, ymax = percentVar+percentSE), width=0, size=1)+
  ylab('Percent of variance') + ylim(0, 48)+
  theme(panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black", size = 0.5), axis.ticks.length = unit(0.2, "cm"),
        text = element_text(size=20), axis.text.x = element_text(size=20,vjust=0.95),
        axis.text.y = element_text(size=20,hjust=0.95),
        axis.title.y = element_text(size=20, hjust=0.5),
        axis.title.x = element_blank()) +
  ggtitle("Nodule number") +
  theme(plot.title = element_text(hjust = 0.5))

fig2a

#***** b. Nodule proportion ####
#use tnodprop

# Keep rows DGE, Main SGE, GxG SGE 
fnodprop = tnodprop[c("DGE", "Main SGE", "GxG SGE"),]

# change row names to a column
fnodprop = tibble::rownames_to_column(fnodprop, "component")

# reorder components
fnodprop$component = factor(fnodprop$component, levels=c("DGE", "Main SGE", "GxG SGE"))

#plot
fig2b= ggplot(fnodprop, aes(x=component, y=percentVar)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=percentVar-percentSE, ymax = percentVar+percentSE), width=0, size=1)+
  ylab('Percent of variance') + ylim(0, 48)+
  theme(panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black", size = 0.5), axis.ticks.length = unit(0.2, "cm"), 
        text = element_text(size=20), axis.text.x = element_text(size=20, vjust=0.95),
        axis.text.y = element_text(size=20,hjust=0.95),
        axis.title.y = element_text(size=20, hjust=0.5),
        axis.title.x = element_blank())+
  ggtitle("Nodule proportion") +
  theme(plot.title = element_text(hjust = 0.5))

fig2b

#***** c. Shoot mass ####
# use tshoot

# Keep rows DGE, Main SGE, GxG SGE 
fshoot = tshoot[c("DGE", "Main SGE", "GxG SGE"),]

# change row names to a column
fshoot = tibble::rownames_to_column(fshoot, "component")

# reorder components
fshoot$component = factor(fshoot$component, levels=c("DGE", "Main SGE", "GxG SGE"))

#plot
fig2c = ggplot(fshoot, aes(x=component, y=percentVar)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=percentVar-percentSE, ymax = percentVar+percentSE), width=0, size=1)+
  xlab('Component') + ylab('Percent of variance') + ylim(0, 48) +
  theme(panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black", size = 0.5), axis.ticks.length = unit(0.2, "cm"),
        text = element_text(size=20), axis.text.x = element_text(size=20, vjust=0.95), 
        axis.text.y = element_text(size=20,hjust=0.95),
        axis.title.y = element_text(size=20, hjust=0.5),
        axis.title.x = element_blank())+
  ggtitle("Shoot mass") +
  theme(plot.title = element_text(hjust = 0.5))

fig2c

# print all together (nodule number, nodule proportion, shoot mass)
grid.arrange(fig2a, fig2b, fig2c, nrow=1)



#*** Heatmaps #### 

# Create subsets for heatmaps (nodule number, nodule proportion, and shootmass)

# Subset competitor strain C394B
C394B.comp=rbind(
  (D1[which(D1$strain.label=="C394B_C264A"),]), 
  (D1[which(D1$strain.label=="C277A_C394B"),]), 
  (D1[which(D1$strain.label=="C394B_C395A"),]), 
  (D1[which(D1$strain.label=="C394B_C412B"),]), 
  (D1[which(D1$strain.label=="C067A_C394B"),]), 
  (D1[which(D1$strain.label=="C066B_C394B"),]),
  (D1[which(D1$strain.label=="C399B_C394B"),]))

# remove C394B from subset
C394B.comp=C394B.comp[-which(C394B.comp$strain.o=="C394B"),]

# Drop levels after subsetting
C394B.comp = droplevels(C394B.comp)

# Subset competitor strain C066B
C066B.comp=rbind(
  (D1[which(D1$strain.label=="C066B_C067A"),]), 
  (D1[which(D1$strain.label=="C066B_C264A"),]), 
  (D1[which(D1$strain.label=="C066B_C277A"),]), 
  (D1[which(D1$strain.label=="C066B_C412B"),]), 
  (D1[which(D1$strain.label=="C066B_C395A"),]), 
  (D1[which(D1$strain.label=="C066B_C394B"),]),
  (D1[which(D1$strain.label=="C066B_C399B"),]))

# remove C066B from subset
C066B.comp=C066B.comp[-which(C066B.comp$strain.o=="C066B"),]

# Drop levels after subsetting
C066B.comp = droplevels(C066B.comp)

# Subset competitor strain C395A 
C395A.comp=rbind(
  (D1[which(D1$strain.label=="C395A_C264A"),]), 
  (D1[which(D1$strain.label=="C277A_C395A"),]), 
  (D1[which(D1$strain.label=="C394B_C395A"),]), 
  (D1[which(D1$strain.label=="C395A_C412B"),]), 
  (D1[which(D1$strain.label=="C067A_C395A"),]), 
  (D1[which(D1$strain.label=="C066B_C395A"),]),
  (D1[which(D1$strain.label=="C399B_C395A"),]))

# remove C395A from subset
C395A.comp=C395A.comp[-which(C395A.comp$strain.o=="C395A"),]

# Drop levels after subsetting
C395A.comp = droplevels(C395A.comp)

# Subset competitor strain C277A 
C277A.comp=rbind(
  (D1[which(D1$strain.label=="C277A_C264A"),]), 
  (D1[which(D1$strain.label=="C277A_C395A"),]), 
  (D1[which(D1$strain.label=="C277A_C394B"),]), 
  (D1[which(D1$strain.label=="C277A_C412B"),]), 
  (D1[which(D1$strain.label=="C067A_C277A"),]), 
  (D1[which(D1$strain.label=="C066B_C277A"),]),
  (D1[which(D1$strain.label=="C399B_C277A"),]))

# remove C277A from subset
C277A.comp=C277A.comp[-which(C277A.comp$strain.o=="C277A"),]

# Drop levels after subsetting
C277A.comp = droplevels(C277A.comp)

# Subset competitor strain C067A 
C067A.comp=rbind(
  (D1[which(D1$strain.label=="C067A_C264A"),]), 
  (D1[which(D1$strain.label=="C067A_C395A"),]), 
  (D1[which(D1$strain.label=="C067A_C394B"),]), 
  (D1[which(D1$strain.label=="C067A_C412B"),]), 
  (D1[which(D1$strain.label=="C067A_C399B"),]), 
  (D1[which(D1$strain.label=="C066B_C067A"),]),
  (D1[which(D1$strain.label=="C067A_C277A"),]))

# remove C067A from subset
C067A.comp=C067A.comp[-which(C067A.comp$strain.o=="C067A"),]

# Drop levels after subsetting
C067A.comp = droplevels(C067A.comp)

# Subset competitor strain C399B (C399B and S399B)
C399B.comp=rbind(
  (D1[which(D1$strain.label=="C399B_C264A"),]), 
  (D1[which(D1$strain.label=="C399B_C395A"),]), 
  (D1[which(D1$strain.label=="C399B_C394B"),]), 
  (D1[which(D1$strain.label=="C399B_C412B"),]), 
  (D1[which(D1$strain.label=="C067A_C399B"),]), 
  (D1[which(D1$strain.label=="C066B_C399B"),]),
  (D1[which(D1$strain.label=="C399B_C277A"),]))

# remove C399B from subset
C399B.comp=C399B.comp[-which(C399B.comp$strain.o=="C399B"),]

# Drop levels after subsetting
C399B.comp = droplevels(C399B.comp)

# Subset competitor strain C264A (C264A and S264A)
C264A.comp=rbind(
  (D1[which(D1$strain.label=="C264A_C412B"),]), 
  (D1[which(D1$strain.label=="C277A_C264A"),]), 
  (D1[which(D1$strain.label=="C394B_C264A"),]), 
  (D1[which(D1$strain.label=="C395A_C264A"),]), 
  (D1[which(D1$strain.label=="C067A_C264A"),]), 
  (D1[which(D1$strain.label=="C066B_C264A"),]),
  (D1[which(D1$strain.label=="C399B_C264A"),]))

# remove C264A from subset
C264A.comp=C264A.comp[-which(C264A.comp$strain.o=="C264A"),]

# Drop levels after subsetting
C264A.comp = droplevels(C264A.comp)

# Subset competitor strain C412B 
C412B.comp=rbind(
  (D1[which(D1$strain.label=="C264A_C412B"),]), 
  (D1[which(D1$strain.label=="C277A_C412B"),]), 
  (D1[which(D1$strain.label=="C394B_C412B"),]), 
  (D1[which(D1$strain.label=="C395A_C412B"),]), 
  (D1[which(D1$strain.label=="C067A_C412B"),]), 
  (D1[which(D1$strain.label=="C066B_C412B"),]),
  (D1[which(D1$strain.label=="C399B_C412B"),]))

# remove C412B from subset
C412B.comp=C412B.comp[-which(C412B.comp$strain.o=="C412B"),]

# Drop levels after subsetting
C412B.comp = droplevels(C412B.comp)


#***** d. Nodule number ####

# Mean nodule number for each strain against the 8 competitors 
# change column names for heatmap code. Each dataframe will show each focal strains mean nodule number against a competitor

# Competitor C394B
ns1=aggregate(nod.sum ~ strain.o, C394B.comp, mean)
ns1
colnames(ns1)=c("focal.strain","C394B")

# Competitor C066B
ns2=aggregate(nod.sum ~ strain.o, C066B.comp, mean)
ns2
colnames(ns2)=c("focal.strain","C066B")

# Competitor C395A 
ns3=aggregate(nod.sum ~ strain.o, C395A.comp, mean)
ns3
colnames(ns3)=c("focal.strain","C395A")

# Competitor C277A 
ns4=aggregate(nod.sum ~ strain.o, C277A.comp, mean)
ns4
colnames(ns4)=c("focal.strain","C277A")

# Competitor C067A 
ns5=aggregate(nod.sum ~ strain.o, C067A.comp, mean)
ns5
colnames(ns5)=c("focal.strain","C067A") 

# Competitor C399B
ns6=aggregate(nod.sum ~ strain.o, C399B.comp, mean)
ns6
colnames(ns6)=c("focal.strain","C399B")

# Competitor C264A
ns7=aggregate(nod.sum ~ strain.o, C264A.comp, mean)
ns7
colnames(ns7)=c("focal.strain","C264A")

# Competitor C412B
ns8=aggregate(nod.sum ~ strain.o, C412B.comp, mean)
ns8
colnames(ns8)=c("focal.strain","C412B")

# Merge all competitor data frames
aa=merge(ns1, ns2, by = "focal.strain", all.x = T, all.y = T)
aa

bb= merge(aa, ns3, by = "focal.strain", all.x = T, all.y = T)
bb

cc= merge(bb, ns4, by = "focal.strain", all.x = T, all.y = T)
cc

dd= merge(cc, ns5, by = "focal.strain", all.x = T, all.y = T)
dd

ee= merge(dd, ns6, by = "focal.strain", all.x = T, all.y = T)
ee

ff= merge(ee, ns7, by = "focal.strain", all.x = T, all.y = T)
ff

gg= merge(ff, ns8, by = "focal.strain", all.x = T, all.y = T)

# add mean nodule number in 1-strain inocula halved along diagonal of heat map (ie to compare to how many nodules a strain formed in competiton)
# mean nodule number already calculated (ss.means). subset and calculate half of mean nodule number in 1-strain inocula
nod.ss = ss.means[,c("strain.o", "nod.sum")]
nod.ss$nod.sum.halved = nod.ss$nod.sum/2
nod.ss
#  strain.o nod.sum nod.sum.halved
#1    C066B 38.7500       19.37500
#2    C067A 48.9375       24.46875
#3    C264A  9.1250        4.56250
#4    C277A 41.4375       20.71875
#5    C394B 49.6250       24.81250
#6    C395A 39.7500       19.87500
#7    C399B 38.7500       19.37500
#8    C412B 12.4375        6.21875

#  add values that will be on heatmap diagonal
gg1= gg %>% 
    mutate(C394B=ifelse(focal.strain=='C394B', 24.81250, C394B),
         C066B=ifelse(focal.strain=='C066B', 19.37500, C066B),
         C395A=ifelse(focal.strain=='C395A', 19.87500, C395A),
         C277A=ifelse(focal.strain=='C277A', 20.71875, C277A),
         C067A=ifelse(focal.strain=='C067A', 24.46875, C067A),
         C399B=ifelse(focal.strain=='C399B', 19.37500, C399B),
         C264A=ifelse(focal.strain=='C264A', 4.56250, C264A),
         C412B=ifelse(focal.strain=='C412B', 6.21875, C412B))

# Now organize data
nod.sum.long = pivot_longer(data=gg1,
                            cols = -("focal.strain"),
                            names_to = "Competitor", 
                            values_to = "Nodule.number")

# Change order of competitor strains: least beneficial to most beneficial
nod.sum.long$Competitor=factor(nod.sum.long$Competitor, levels=c("C412B", "C264A", "C399B", "C067A", "C277A", "C395A", "C066B", "C394B"))

# Change order of focal strains: least beneficial to most beneficial
nod.sum.long$focal.strain=factor(nod.sum.long$focal.strain, levels=c("C412B", "C264A", "C399B", "C067A", "C277A", "C395A", "C066B", "C394B"))


# Nodule number heatmap   
fig2d = ggplot(data = nod.sum.long, mapping = aes(x = focal.strain, y = Competitor, fill = Nodule.number)) +
  geom_tile(color = "gray") +
  scale_fill_continuous(high = "#132B43", low = "#56B1F7", na.value="white") +
  scale_x_discrete(labels=c('412B','264A','399B','067A','277A','395A','066B','394B')) +
  scale_y_discrete(labels=c('412B','264A','399B','067A','277A','395A','066B','394B')) +
  xlab(label = "Focal strain") +
  ylab(label = "Competitor strain")+
theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
        axis.ticks = element_blank(), legend.position="right", legend.justification="top", legend.title=element_text(), 
        text = element_text(size=30), axis.text.x = element_text(size=28, vjust=0.5, angle=90), axis.text.y = element_text(size=28,hjust=0.95))+
  labs(fill='Focal
strain 
nodule 
number')

fig2d




#***** e. Nodule proportion ####

# Mean nodule proportion for each strain against the 8 competitors 
# change column names for heatmap code. Each dataframe will show each focal strains mean nodule proportion against a competitor

# Competitor C394B
np1=aggregate(nod.prop ~ strain.o, C394B.comp, mean)
np1
colnames(np1)=c("focal.strain","C394B")

# Competitor C066B
np2=aggregate(nod.prop ~ strain.o, C066B.comp, mean)
np2
colnames(np2)=c("focal.strain","C066B")

# Competitor C395A 
np3=aggregate(nod.prop ~ strain.o, C395A.comp, mean)
np3
colnames(np3)=c("focal.strain","C395A")

# Competitor C277A 
np4=aggregate(nod.prop ~ strain.o, C277A.comp, mean)
np4
colnames(np4)=c("focal.strain","C277A")

# Competitor C067A 
np5=aggregate(nod.prop ~ strain.o, C067A.comp, mean)
np5
colnames(np5)=c("focal.strain","C067A") 

# Competitor C399B
np6=aggregate(nod.prop ~ strain.o, C399B.comp, mean)
np6
colnames(np6)=c("focal.strain","C399B")

# Competitor C264A
np7=aggregate(nod.prop ~ strain.o, C264A.comp, mean)
np7
colnames(np7)=c("focal.strain","C264A")

# Competitor C412B
np8=aggregate(nod.prop ~ strain.o, C412B.comp, mean)
np8
colnames(np8)=c("focal.strain","C412B")

# Merge all competitor data frames
aa2=merge(np1, np2, by = "focal.strain", all.x = T, all.y = T)
aa2

bb2= merge(aa2, np3, by = "focal.strain", all.x = T, all.y = T)
bb2

cc2= merge(bb2, np4, by = "focal.strain", all.x = T, all.y = T)
cc2

dd2= merge(cc2, np5, by = "focal.strain", all.x = T, all.y = T)
dd2

ee2= merge(dd2, np6, by = "focal.strain", all.x = T, all.y = T)
ee2

ff2= merge(ee2, np7, by = "focal.strain", all.x = T, all.y = T)
ff2

gg2= merge(ff2, np8, by = "focal.strain", all.x = T, all.y = T)

# Now organize data
nod.prop.long = pivot_longer(data=gg2,
                            cols = -("focal.strain"),
                            names_to = "Competitor", 
                            values_to = "Nodule.prop")

# Change order of competitor strains: least beneficial to most beneficial
nod.prop.long$Competitor=factor(nod.prop.long$Competitor, levels=c("C412B", "C264A", "C399B", "C067A", "C277A", "C395A", "C066B", "C394B"))

# Change order of focal strains: least beneficial to most beneficial
nod.prop.long$focal.strain=factor(nod.prop.long$focal.strain, levels=c("C412B", "C264A", "C399B", "C067A", "C277A", "C395A", "C066B", "C394B"))

# Nodule proportion heatmap  
fig2e = ggplot(data = nod.prop.long, mapping = aes(x = focal.strain, y = Competitor, fill = Nodule.prop)) +
  geom_tile(color="gray") +
  scale_fill_continuous(high = "#132B43", low = "#56B1F7", na.value="white") +
  scale_x_discrete(labels=c('412B','264A','399B','067A','277A','395A','066B','394B')) +
  scale_y_discrete(labels=c('412B','264A','399B','067A','277A','395A','066B','394B')) +
  xlab(label = "Focal strain") +
  ylab(label = "Competitor strain")+
  theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
        axis.ticks = element_blank(), legend.position="right", legend.justification="top", legend.title=element_text(), 
        text = element_text(size=30), axis.text.x = element_text(size=28, vjust=0.5, angle=90), axis.text.y = element_text(size=28,hjust=0.95))+
  labs(fill='Focal 
strain
nodule 
proportion')

fig2e

#***** f. Shoot mass ####

# Mean shoot mass of each 2-strain combination 
# change column names for heatmap code. Each dataframe will show the mean shoot mass for every 2-strain combination

# Competitor C394B
sm1=aggregate(shoot.mass ~ strain.o, C394B.comp, mean)
sm1
colnames(sm1)=c("focal.strain","C394B")

# Competitor C066B
sm2=aggregate(shoot.mass ~ strain.o, C066B.comp, mean)
sm2
colnames(sm2)=c("focal.strain","C066B")

# Competitor C395A 
sm3=aggregate(shoot.mass ~ strain.o, C395A.comp, mean)
sm3
colnames(sm3)=c("focal.strain","C395A")

# Competitor C277A 
sm4=aggregate(shoot.mass ~ strain.o, C277A.comp, mean)
sm4
colnames(sm4)=c("focal.strain","C277A")

# Competitor C067A 
sm5=aggregate(shoot.mass ~ strain.o, C067A.comp, mean)
sm5
colnames(sm5)=c("focal.strain","C067A") 

# Competitor C399B
sm6=aggregate(shoot.mass ~ strain.o, C399B.comp, mean)
sm6
colnames(sm6)=c("focal.strain","C399B")

# Competitor C264A
sm7=aggregate(shoot.mass ~ strain.o, C264A.comp, mean)
sm7
colnames(sm7)=c("focal.strain","C264A")

# Competitor C412B
sm8=aggregate(shoot.mass ~ strain.o, C412B.comp, mean)
sm8
colnames(sm8)=c("focal.strain","C412B")

# Merge all competitor data frames
aa3=merge(sm1, sm2, by = "focal.strain", all.x = T, all.y = T)
aa3

bb3= merge(aa3, sm3, by = "focal.strain", all.x = T, all.y = T)
bb3

cc3= merge(bb3, sm4, by = "focal.strain", all.x = T, all.y = T)
cc3

dd3= merge(cc3, sm5, by = "focal.strain", all.x = T, all.y = T)
dd3

ee3= merge(dd3, sm6, by = "focal.strain", all.x = T, all.y = T)
ee3

ff3= merge(ee3, sm7, by = "focal.strain", all.x = T, all.y = T)
ff3

gg3= merge(ff3, sm8, by = "focal.strain", all.x = T, all.y = T)

# add mean shoot mass in 1-strain inocula along diagonal of heat map (ie to compare to shoot mass in strain competiton)
# mean shoot mass already calculated (ss.means)
shoot.ss = ss.means[,c("strain.o", "shoot.mass.ss.o")]
shoot.ss
#  strain.o shoot.mass.ss.o
#1    C066B      0.05160000
#2    C067A      0.03918750
#3    C264A      0.02037500
#4    C277A      0.04092500
#5    C394B      0.05241250
#6    C395A      0.04495000
#7    C399B      0.03723125
#8    C412B      0.02000625

# add values that will be on heatmap diagonal
gg3= gg3 %>% 
    mutate(C394B=ifelse(focal.strain=='C394B', 0.05241250, C394B),
         C066B=ifelse(focal.strain=='C066B', 0.05160000, C066B),
         C395A=ifelse(focal.strain=='C395A', 0.04495000, C395A),
         C277A=ifelse(focal.strain=='C277A', 0.04092500, C277A),
         C067A=ifelse(focal.strain=='C067A', 0.03918750, C067A),
         C399B=ifelse(focal.strain=='C399B', 0.03723125, C399B),
         C264A=ifelse(focal.strain=='C264A', 0.02037500, C264A),
         C412B=ifelse(focal.strain=='C412B', 0.02000625, C412B))

# Now organize data
shoot.mass.long = pivot_longer(data=gg3,
                               cols = -("focal.strain"),
                               names_to = "Competitor", 
                               values_to = "Shoot.mass")

# Change order of competitor strains: least beneficial to most beneficial
shoot.mass.long$Competitor=factor(shoot.mass.long$Competitor, levels=c("C412B", "C264A", "C399B", "C067A", "C277A", "C395A", "C066B", "C394B"))

# Change order of focal strains: least beneficial to most beneficial
shoot.mass.long$focal.strain=factor(shoot.mass.long$focal.strain, levels=c("C412B", "C264A", "C399B", "C067A", "C277A", "C395A", "C066B", "C394B"))

# remove top left half of heatmap (redundant data shown) #returns warning but OK
shoot.mass.longf= shoot.mass.long %>% 
    mutate(Shoot.mass=ifelse(focal.strain=='C066B' & Competitor=='C394B', NA, Shoot.mass),
         Shoot.mass=ifelse(focal.strain=='C395A' & Competitor==c('C394B','C066B'), NA, Shoot.mass),
         Shoot.mass=ifelse(focal.strain=='C277A' & Competitor==c('C394B','C066B','C395A'), NA, Shoot.mass),
         Shoot.mass=ifelse(focal.strain=='C067A' & Competitor==c('C394B','C066B','C395A','C277A'), NA, Shoot.mass),
         Shoot.mass=ifelse(focal.strain=='C399B' & Competitor==c('C394B','C066B','C395A','C277A','C067A'), NA, Shoot.mass),
         Shoot.mass=ifelse(focal.strain=='C264A' & Competitor!=c('C264A','C412B'), NA, Shoot.mass),
         Shoot.mass=ifelse(focal.strain=='C412B' & Competitor!='C412B', NA, Shoot.mass))

# shoot mass heatmap  
fig2f = ggplot(data = shoot.mass.longf, mapping = aes(x = focal.strain, y = Competitor, fill = Shoot.mass)) +
  geom_tile(color = "white") +
  scale_fill_continuous(high = "#132B43", low = "#56B1F7", na.value="white") +
  scale_x_discrete(labels=c('412B','264A','399B','067A','277A','395A','066B','394B')) +
  scale_y_discrete(position = "right",labels=c('412B','264A','399B','067A','277A','395A','066B','394B')) +
  xlab(label = "Focal strain") +
  ylab(label = "Competitor strain")+
  theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
        axis.ticks = element_blank(), legend.position="right", legend.justification="top", legend.title=element_text(), 
        text = element_text(size=30), axis.text.x = element_text(size=28, vjust=0.5, angle=90), axis.text.y = element_text(size=28,hjust=0.95), axis.title.y.right = element_text(vjust=1.5))+
  labs(fill='Shoot 
mass (g) 
per pot')

fig2f

#** FIG 3: Benefit of partner choice ####

# Benefit to host figure prep
D1m2$combo = paste(D1m2$strain.o, D1m2$strain.comp.o,sep=':')

# aggregate the means nodpropbeststrain and "dne" for the plot
D1m3 = aggregate(cbind(NodPropBestStrain, dne) ~ combo + strain.o + strain.comp.o, FUN=mean, data=D1m2)
# add standard error
D1m3 = cbind(D1m3, aggregate(cbind(NodPropBestStrain, dne) ~ combo + strain.o + strain.comp.o, FUN=function(x) sd(x)/length(x), data=D1m2)[4:5])
head(D1m3)
# Rename se columns
names(D1m3)[6:7] = paste(names(D1m3)[6:7],'se',sep='.')

## Fig 3: partner choice benefits the host
ggplot(D1m3, aes(x=NodPropBestStrain, y=dne, alpha=.6)) + 
  geom_point(size=4) + 
  geom_errorbar(aes(ymin=dne-dne.se, ymax=dne+dne.se)) +
  geom_errorbarh(aes(xmin=NodPropBestStrain-NodPropBestStrain.se, xmax=NodPropBestStrain+NodPropBestStrain.se)) +
  stat_smooth(method='lm', alpha=.2, size=1, color='black') +
  geom_hline(yintercept = 0, linetype='longdash', size=1) + #plant fitness under neutral nod occupancy
  geom_vline(xintercept = 0.5, linetype='dashed', color='black', size=1) + #partner choice:better strain occupies >50% of nodules
  theme(legend.position='none') + xlab('Nodule proportion of more beneficial strain')+ ylab('Shoot mass (g) deviation from neutral expectation') +
  theme(panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),axis.line = element_line(colour = "black"),axis.ticks = element_line(colour = "black", size = 0.5), axis.ticks.length = unit(0.2, "cm"),text = element_text(size=22), axis.text.x = element_text(size=22, vjust=0.95), axis.text.y = element_text(size=22,hjust=0.95)) 




#** FIG 4: Host sanctions (2-strain inocula) ####

#***** a. single strain nodules: cfu per nodule (2-strain) by shoot mass (1-strain) ####
ggplot(cfu.full.shoot.means, aes(x=shoot.mass.ss.o, y=mean.cfu)) + geom_point(size=6) + stat_smooth(method='lm', size=3, color='black') + 
  theme_set(theme_bw(base_size=6)) + xlab('Shoot mass (g) (1-strain)') + ylab('CFU per nodule') +
  theme(panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),
  axis.line = element_line(colour = "black"),axis.ticks = element_line(colour = "black", size = 1.5), axis.ticks.length = unit(0.25, "cm"),
  text = element_text(size=32), axis.text.x = element_text(size=32, vjust=0.95), axis.text.y = element_text(size=32,hjust=0.95)) +
  stat_cor(label.y = 2700000, label.x=0.02,size=10, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
    ggtitle("Single strain nodules") +
  theme(plot.title = element_text(size=32, hjust = 0.5))


#***** b. mixed nodules: cfu proportion (2-strain) by shoot mass (1-strain) ####
ggplot(cfu.prop.mixed.shoot.means, aes(x=shoot.mass.ss.o, y=mean.cfu.prop)) + geom_point(size=6) + stat_smooth(method='lm', size=3, color="black") + 
  theme_set(theme_bw(base_size=6))+xlab('Shoot mass (g) (1-strain)') +
  ylab('CFU proportion per nodule') +
  theme(panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),
  axis.line = element_line(colour = "black"),axis.ticks = element_line(colour = "black", size = 1.5), axis.ticks.length = unit(0.25, "cm"),
  text = element_text(size=32), axis.text.x = element_text(size=32, vjust=0.95), axis.text.y = element_text(size=32,hjust=0.95)) +
  stat_cor(label.y = 0.8, label.x=0.02, size=10,aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
      ggtitle("Mixed nodules") +
  theme(plot.title = element_text(size=32, hjust = 0.5))



#### SUPPLEMENTARY FIGURES & ANALYSES #### 

#** FIG S2: Verification of unbiased half dataset ####
# nodule number used the full data and is not presented here
# nodule proportion and shoot mass used the half dataset

#***** Nodule proportion ####
# NOTE: can skip to "Read in file from the output" since code takes a while to run


# run iterations of model using randomized datasets 

# make empty df for storing output of each model iteration
# each row will store output for one iteration
# Columns with "v" prefixes store variance percentages
# Columns with "p" prefixes store p-values from LRTs

iter.df1 = data.frame(v.GxG = numeric(),
                      v.strain.o = numeric(),
                      v.strain.comp.o = numeric(),
                      v.Block = numeric(),
                      v.Residual = numeric(),
                      p.GxG = numeric(),
                      p.strain.o = numeric(),
                      p.strain.comp.o = numeric(),
                      p.Block = numeric())



# NOTE: this takes a long time to run 
for (i in 1:10000) {
  
  # create mixed data set with only one observation per pot (choose between strain and strain.comp at random)
  D1m.rand = ddply(D1m, .(Unique.ID), function(x) x[sample(nrow(x),1),])
  
  # run model on randomized dataset
  M7 = lmer(nod.prop ~ (1|strain.o:strain.comp.o) + (1|strain.o) + (1|strain.comp.o) + (1|Block), data=D1m.rand, na.action=na.omit)
  
  # partition the variance from the model
  var.compM7 = as.data.frame(summary(M7)$varcor)
  var.compM7$variance = var.compM7$sdcor^2
  var.compM7$percent = 100*(var.compM7$variance/sum(var.compM7$variance))
  
  # transpose the small variance dataframe so that variance estimates are in a single row
  temp = select(var.compM7, percent) %>% t()
  
  # add the row of variance estimates to the ith row of the iteration dataframe
  iter.df1[i,1:5] = temp[1:5]
  
  # test significance of GxG and add p-value to dataframe
  M7x = update(M7, .~. -(1|strain.o:strain.comp.o))
  results = anova(M7, M7x)
  pvalue = results$`Pr(>Chisq)`[2]
  iter.df1[i,]$p.GxG = pvalue
  
  # test significance of strain.o and add p-value to dataframe
  M7x = update(M7, .~. -(1|strain.o))
  results = anova(M7, M7x)
  pvalue = results$`Pr(>Chisq)`[2]
  iter.df1[i,]$p.strain.o = pvalue
  
  # test significance of strain.comp.o and add p-value to dataframe
  M7x = update(M7, .~. -(1|strain.comp.o))
  results = anova(M7, M7x)
  pvalue = results$`Pr(>Chisq)`[2]
  iter.df1[i,]$p.strain.comp.o = pvalue
  
  # test significance of Block and add p-value to dataframe
  M7x = update(M7, .~. -(1|Block))
  results = anova(M7, M7x)
  pvalue = results$`Pr(>Chisq)`[2]
  iter.df1[i,]$p.Block = pvalue
  
}

# write the iteration dataframe to a file for safekeeping, since it takes a while to make
# write.csv(x = iter.df1, file = paste("noduleproportion_models_randomized.csv"), row.names = F)

#****** Read in file from the output ####
iter.df1=read.csv("noduleproportion_models_randomized.csv")


# summarize data from model iterations

# look at histograms of iteration df

# variance percentages
par(mfrow = c(2,2))
hist(iter.df1$v.GxG)
hist(iter.df1$v.strain.o)
hist(iter.df1$v.strain.comp.o)
hist(iter.df1$v.Block)


# make an empty dataframe to store summary statistics from the iteration dataframe
summary.df = data.frame(Source = factor(levels = c("v.GxG", 
                                                    "v.strain.o",
                                                    "v.strain.comp.o")),
                         Median = numeric(),
                         Mean = numeric(),
                         LCL = numeric(),
                         UCL = numeric())

# fill in Source column in summary.df
summary.df[1:3,]$Source = c("v.GxG", 
                             "v.strain.o",
                             "v.strain.comp.o")

# calculate summary statistics for each column of the iteration df
# summary statistics include:
# lower confidence interval (LCL), which is the 2.5% quantile
# upper confidence interval (UCL), which is the 97.5% quantile
# Median (i.e., 50% quantile)
# Mean

for (i in 1:3) {
  temp = quantile(iter.df1[,i], probs = seq(0, 1, 0.025))
  LCL = temp["2.5%"]
  Median = temp["50%"]
  Mean = mean(iter.df1[,i])
  UCL = temp["97.5%"]
  
  # Add the summary statistics to the ith row of the summary dataframe
  summary.df[i,]$Median = Median
  summary.df[i,]$Mean = Mean
  summary.df[i,]$LCL = LCL
  summary.df[i,]$UCL = UCL
}

# view the output
summary.df

# plot mean variance percentages

# just pull out variance data
usedata = summary.df[1:3,]

# make labels and order them 
usedata$Source.label = ifelse(usedata$Source=="v.GxG", "GxG SGE",
                               ifelse(usedata$Source=="v.strain.o", "DGE",
                                      ifelse(usedata$Source=="v.strain.comp.o", "Main SGE", "")))

usedata$Source.label =ordered(usedata$Source.label, levels = c("DGE", "Main SGE", "GxG SGE"))

# plot
p = ggplot(usedata, aes(x=Source.label, y=Mean)) + 
  geom_col() +
  geom_errorbar(aes(ymin=LCL, ymax=UCL, width = 0.1)) +
  xlab("") +
  ylab("Percent of variance") + ylim(0, 50)+
  theme(panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black", size = 0.5), axis.ticks.length = unit(0.2, "cm"),
        text = element_text(size=20), axis.text.x = element_text(size=20, vjust=0.95), 
        axis.text.y = element_text(size=20,hjust=0.95),
        axis.title.y = element_text(size=20, hjust=0.5),
        axis.title.x = element_blank())
p1 = arrangeGrob(p, top = textGrob("Nodule Proportion", x = unit(0, "npc"),
                                     y = unit(1, "npc"), just=c("left","top"),
                                     gp=gpar(col="black", fontsize=24)))
grid.arrange(p1)



#***** Shoot mass ####
# NOTE: can skip to "Read in file from the output" since code takes a while to run

# *** run iterations of model using randomized datasets

# make empty df for storing output of each model iteration
# each row will store output for one iteration
# Columns with "v" prefixes store variance percentages
# Columns with "p" prefixes store p-values from LRTs

iter.df2 = data.frame(v.GxG = numeric(),
                     v.strain.o = numeric(),
                     v.strain.comp.o = numeric(),
                     v.Block = numeric(),
                     v.Residual = numeric(),
                     p.GxG = numeric(),
                     p.strain.o = numeric(),
                     p.strain.comp.o = numeric(),
                     p.Block = numeric())


# NOTE: this takes a long time to run
for (i in 1:10000) {

# create mixed data set with only one observation per pot (choose between strain and strain.comp at random)
D1m.rand = ddply(D1m, .(Unique.ID), function(x) x[sample(nrow(x),1),])

# run model on randomized dataset
M8 = lmer(log(shoot.mass) ~ (1|strain.o:strain.comp.o) + (1|strain.o) + (1|strain.comp.o) + (1|Block), data=D1m.rand, na.action=na.omit)

# partition the variance from the model
var.compM8 = as.data.frame(summary(M8)$varcor)
var.compM8$variance = var.compM8$sdcor^2
var.compM8$percent = 100*(var.compM8$variance/sum(var.compM8$variance))

# transpose the small variance dataframe so that variance estimates are in a single row
temp = select(var.compM8, percent) %>% t()

# add the row of variance estimates to the ith row of the iteration dataframe
iter.df2[i,1:5] = temp[1:5]

# test significance of GxG and add p-value to dataframe
M8x = update(M8, .~. -(1|strain.o:strain.comp.o))
results = anova(M8, M8x)
pvalue = results$`Pr(>Chisq)`[2]
iter.df2[i,]$p.GxG = pvalue

# test significance of strain.o and add p-value to dataframe
M8x = update(M8, .~. -(1|strain.o))
results = anova(M8, M8x)
pvalue = results$`Pr(>Chisq)`[2]
iter.df2[i,]$p.strain.o = pvalue

# test significance of strain.comp.o and add p-value to dataframe
M8x = update(M8, .~. -(1|strain.comp.o))
results = anova(M8, M8x)
pvalue = results$`Pr(>Chisq)`[2]
iter.df2[i,]$p.strain.comp.o = pvalue

# test significance of Block and add p-value to dataframe
M8x = update(M8, .~. -(1|Block))
results = anova(M8, M8x)
pvalue = results$`Pr(>Chisq)`[2]
iter.df2[i,]$p.Block = pvalue

}

# write the iteration dataframe to a file for safekeeping, since it takes a while to make
# write.csv(x = iter.df2, file = paste("shootmass_models_randomized.csv"), row.names = F)

#****** Read in file from the output ####
iter.df2 = read.csv("shootmass_models_randomized.csv")

# *** summarize data from model iterations ####

# look at histograms of iteration df

# variance percentages
par(mfrow = c(2,2))
hist(iter.df2$v.GxG)
hist(iter.df2$v.strain.o)
hist(iter.df2$v.strain.comp.o)
hist(iter.df2$v.Block)

# make an empty dataframe to store summary statistics from the iteration dataframe
summary.df = data.frame(Source = factor(levels = c("v.GxG", 
                                                    "v.strain.o",
                                                    "v.strain.comp.o")),
                         Median = numeric(),
                         Mean = numeric(),
                         LCL = numeric(),
                         UCL = numeric())

# fill in Source column in summary.df
summary.df[1:3,]$Source = c("v.GxG", 
                             "v.strain.o",
                             "v.strain.comp.o")

# calculate summary statistics for each column of the iteration df
# summary statistics include:
# lower confidence interval (LCL), which is the 2.5% quantile
# upper confidence interval (UCL), which is the 97.5% quantile
# Median (i.e., 50% quantile)
# Mean

for (i in 1:3) {
  temp = quantile(iter.df2[,i], probs = seq(0, 1, 0.025))
  LCL = temp["2.5%"]
  Median = temp["50%"]
  Mean = mean(iter.df2[,i])
  UCL = temp["97.5%"]
  
  # Add the summary statistics to the ith row of the summary dataframe
  summary.df[i,]$Median = Median
  summary.df[i,]$Mean = Mean
  summary.df[i,]$LCL = LCL
  summary.df[i,]$UCL = UCL
}

# view the output
summary.df

# plot mean variance percentages

# just pull out variance data
usedata = summary.df[1:3,]

# make labels and order them 
usedata$Source.label = ifelse(usedata$Source=="v.GxG", "GxG SGE",
                               ifelse(usedata$Source=="v.strain.o", "DGE",
                                      ifelse(usedata$Source=="v.strain.comp.o", "Main SGE", "")))

usedata$Source.label =ordered(usedata$Source.label, levels = c("DGE", "Main SGE", "GxG SGE"))

# plot
p = ggplot(usedata, aes(x=Source.label, y=Mean)) + 
  geom_col() +
  geom_errorbar(aes(ymin=LCL, ymax=UCL, width = 0.1)) +
  xlab("") +
  ylab("Percent of variance") + ylim(0, 50)+
  theme(panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black", size = 0.5), axis.ticks.length = unit(0.2, "cm"),
        text = element_text(size=20), axis.text.x = element_text(size=20, vjust=0.95), 
        axis.text.y = element_text(size=20,hjust=0.95),
        axis.title.y = element_text(size=20, hjust=0.5),
        axis.title.x = element_blank())
p2 = arrangeGrob(p, top = textGrob("Shoot mass", x = unit(0, "npc"),
                                     y = unit(1, "npc"), just=c("left","top"),
                                     gp=gpar(col="black", fontsize=24)))
grid.arrange(p2)

#
# print all figures together (nodule proportion, shoot mass)
grid.arrange(p1, p2, nrow = 1)



#** FIG S3: Nodule number per shoot mass in 1-strain inocula ####

# Is the trend we observe that more beneficial strains form more nodules simply driven by larger plants forming more nodules?  

# calculate nodule number per shoot in 1-strain inocula
D1s.sub = D1s[,c(2,8,17,21)]
head(D1s.sub)
D1s.sub$nod.per.shoot = (D1s.sub$nod.sum+D1s.sub$ambiguous.nods)/D1s.sub$shoot.mass

mean.nodpershoot = aggregate(nod.per.shoot ~ strain.o, FUN=mean, data=D1s.sub)
mean.nodpershoot$nod.per.shoot.se = aggregate(nod.per.shoot ~ strain.o, FUN=function(x) sd(x)/sqrt(length(x)), data=D1s.sub)[,2]

# add number of nodules per shoot mass in 1-strain inocula to ss.means
ss.means = merge(ss.means, mean.nodpershoot, by="strain.o")

head(ss.means)

# Do more beneficial strains form more nodules per shoot mass in 1-strain inocula 
M1.1 = lm(nod.per.shoot ~ shoot.mass.ss.o, data = ss.means)

summary(M1.1)
#Coefficients:
#                Estimate Std. Error t value Pr(>|t|)  
#(Intercept)        269.7      299.4   0.901   0.4025  
#shoot.mass.ss.o  19982.9     7473.3   2.674   0.0368 *
#Residual standard error: 246.2 on 6 degrees of freedom
#Multiple R-squared:  0.5437,	Adjusted R-squared:  0.4677 
#F-statistic:  7.15 on 1 and 6 DF,  p-value: 0.03683

# plot residuals 
par(mfrow=c(2,2))
plot(M1.1)

# test effect of shoot mass in 1-strain inocula on nodule number per shoot mass in 1-strain inocula
drop1(M1.1, test = "F")
#                Df Sum of Sq    RSS    AIC F value  Pr(>F)  
#<none>                       363605 89.795                  
#shoot.mass.ss.o  1    433280 796885 94.072  7.1498 0.03683 *

# This indicates the trend we found that more beneficial strains form more nodules in 1-strain inocula is not just driven by larger plants forming more nodules 


#plot of single strain nodule number per shoot mass predicted by shoot mass 
ggplot(ss.means, aes(x=shoot.mass.ss.o, y=nod.per.shoot)) + geom_point(size=6) + stat_smooth(method='lm', size=3, color='black') + 
  theme_bw() + xlab('Shoot mass (g) (1-strain)') +
  ylab(expression(atop('Nodules per shoot mass',paste((g^{-1}),"(1-strain)")))) +
  theme(panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),axis.line = element_line(colour = "black"),axis.ticks = element_line(colour = "black", size = 1.5), axis.ticks.length = unit(0.25, "cm"),text = element_text(size=32), axis.text.x = element_text(size=28, vjust=0.95), axis.text.y = element_text(size=28,hjust=0.95)) +
  stat_cor(label.y = 1600, label.x=0.02,size=10, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))


#** FIG S4: Shoot mass per nodule in 1-strain inocula ####

#calculate the benefit each strain provides to the host
#subtract the mean of the neg control from mean shoot mass
# mean of neg = 0.020125 g
ss.means$shoot.temp = ss.means$shoot.mass.ss.o - 0.020125
ss.means
#shoot mass per nodule means
ss.means$shoot.per.nod = ss.means$shoot.temp/ss.means$nod.sum

# Do more beneficial strains provide more benefit to the host on a per nodule basis in 1-strain inocula? 
M1.2 = lm(shoot.per.nod ~ shoot.mass.ss.o, data = ss.means)

summary(M1.2)
#Coefficients:
#                  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)     -4.496e-04  7.683e-05  -5.852   0.0011 ** 
#shoot.mass.ss.o  2.294e-02  1.918e-03  11.962 2.07e-05 ***
#Residual standard error: 6.317e-05 on 6 degrees of freedom
#Multiple R-squared:  0.9598,	Adjusted R-squared:  0.953 
#F-statistic: 143.1 on 1 and 6 DF,  p-value: 2.068e-05

# plot residuals 
par(mfrow=c(2,2))
plot(M1.2)

# test effect of shoot mass in 1-strain inocula on shoot mass per nodule number in 1-strain inocula
drop1(M1.2, test = "F")
#                Df  Sum of Sq        RSS     AIC F value    Pr(>F)    
#<none>                        2.3940e-08 -153.02                      
#shoot.mass.ss.o  1 5.7094e-07 5.9488e-07 -129.31  143.09 2.068e-05 ***

# More beneficial strains provide more benefit to the host per nodule

#plot of single strain shoot per nodule predicted by shoot mass 
ggplot(ss.means, aes(x=shoot.mass.ss.o, y=shoot.per.nod)) + geom_point(size=6) + stat_smooth(method='lm', size=3, color='black') + 
  theme_bw() + xlab('Shoot mass (g) (1-strain)') + ylab('Shoot mass (g) increase 
per nodule (1-strain)')  +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  theme(panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),axis.line = element_line(colour = "black"),axis.ticks = element_line(colour = "black", size = 1.5), axis.ticks.length = unit(0.25, "cm"),text = element_text(size=32), axis.text.x = element_text(size=28, vjust=0.95), axis.text.y = element_text(size=28,hjust=0.95)) +
  stat_cor(label.y = 0.001, label.x=0.02,size=10, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))



#** Growth rate vs. nodule number (1-strain) or proportion (2-strains) ####

# growth data
growth=read.csv("growth_data.csv")
head(growth)
#remove meta.data
growth=growth[,2:4]

#change strain column name to match other files
colnames(growth)[1]= "strain.o"

# calc mean growth per strain  
mean.growth = aggregate(OD72hr ~ strain.o, FUN=mean, data=growth)
head(mean.growth)

# merge growth means with ss.means data
ss.means = merge(ss.means, mean.growth, by="strain.o")

# model nod number in 1 strain inocula predicted by growth rate
Mgrow=lm(nod.sum~OD72hr, data=ss.means)

summary(Mgrow)
#Coefficients:
#            Estimate Std. Error t value Pr(>|t|)  
#(Intercept)    24.70      12.19   2.025   0.0892 .
#OD72hr         22.53      24.12   0.934   0.3864  
#Residual standard error: 15.63 on 6 degrees of freedom
#Multiple R-squared:  0.1269,	Adjusted R-squared:  -0.01861 
#F-statistic: 0.8721 on 1 and 6 DF,  p-value: 0.3864


# plot residuals 
par(mfrow=c(2,2))
plot(Mgrow)

# test fixed effect
drop1(Mgrow, test = "F")
#       Df Sum of Sq    RSS    AIC F value Pr(>F)
#<none>              1465.7 45.685               
#OD72hr  1    213.03 1678.7 44.771  0.8721 0.3864

# growth has no relationship with nodule number in 1-strain inocula 

# model nod prop in 2 strain inocula predicted by growth rate
Mgrow.prop=lm(nod.prop~OD72hr, data=ss.means)

summary(Mgrow.prop)
#Coefficients:
#            Estimate Std. Error t value Pr(>|t|)
#(Intercept)   0.2024     0.2110   0.959    0.374
#OD72hr        0.6429     0.4174   1.540    0.174
#Residual standard error: 0.2704 on 6 degrees of freedom
#Multiple R-squared:  0.2834,	Adjusted R-squared:  0.1639 
#F-statistic: 2.372 on 1 and 6 DF,  p-value: 0.1744


# plot residuals 
par(mfrow=c(2,2))
plot(Mgrow.prop)

# test fixed effect
drop1(Mgrow.prop, test = "F")
#       Df Sum of Sq     RSS     AIC F value Pr(>F)
#<none>              0.43871 -19.227               
#OD72hr  1   0.17346 0.61217 -18.561  2.3723 0.1744

# growth has no relationship with nodule proportion in 2-strain inocula 


#** Main and GxG SGE examples ####

#subset dataframes
# subset focal strain C394B
Rev.C394B=rbind(
  (D1[which(D1$strain.label=="C394B_C264A"),]), 
  (D1[which(D1$strain.label=="C277A_C394B"),]), 
  (D1[which(D1$strain.label=="C394B_C395A"),]), 
  (D1[which(D1$strain.label=="C394B_C412B"),]), 
  (D1[which(D1$strain.label=="C067A_C394B"),]), 
  (D1[which(D1$strain.label=="C066B_C394B"),]),
  (D1[which(D1$strain.label=="C399B_C394B"),]))

# remove other focal strains besides C394B from subset
Rev.C394B=Rev.C394B[-which(Rev.C394B$strain.o!="C394B"),]

# Drop levels after subsetting
Rev.C394B = droplevels(Rev.C394B)


# subset focal strain C066B 
Rev.C066B=rbind(
  (D1[which(D1$strain.label=="C066B_C067A"),]), 
  (D1[which(D1$strain.label=="C066B_C264A"),]), 
  (D1[which(D1$strain.label=="C066B_C277A"),]), 
  (D1[which(D1$strain.label=="C066B_C412B"),]), 
  (D1[which(D1$strain.label=="C066B_C395A"),]), 
  (D1[which(D1$strain.label=="C066B_C394B"),]),
  (D1[which(D1$strain.label=="C066B_C399B"),]))

# remove other focal strains besides C066B from subset
Rev.C066B=Rev.C066B[-which(Rev.C066B$strain.o!="C066B"),]

# Drop levels after subsetting
Rev.C066B = droplevels(Rev.C066B)


# subset focal strain C395A 
Rev.C395A=rbind(
  (D1[which(D1$strain.label=="C395A_C264A"),]), 
  (D1[which(D1$strain.label=="C277A_C395A"),]), 
  (D1[which(D1$strain.label=="C394B_C395A"),]), 
  (D1[which(D1$strain.label=="C395A_C412B"),]), 
  (D1[which(D1$strain.label=="C067A_C395A"),]), 
  (D1[which(D1$strain.label=="C066B_C395A"),]),
  (D1[which(D1$strain.label=="C399B_C395A"),]))

# remove other focal strains besides C395A from subset
Rev.C395A=Rev.C395A[-which(Rev.C395A$strain.o!="C395A"),]

# Drop levels after subsetting
Rev.C395A = droplevels(Rev.C395A)


# subset focal strain C277A 
Rev.C277A=rbind(
  (D1[which(D1$strain.label=="C277A_C264A"),]), 
  (D1[which(D1$strain.label=="C277A_C395A"),]), 
  (D1[which(D1$strain.label=="C277A_C394B"),]), 
  (D1[which(D1$strain.label=="C277A_C412B"),]), 
  (D1[which(D1$strain.label=="C067A_C277A"),]), 
  (D1[which(D1$strain.label=="C066B_C277A"),]),
  (D1[which(D1$strain.label=="C399B_C277A"),]))

# remove other focal strains besides C277A from subset
Rev.C277A=Rev.C277A[-which(Rev.C277A$strain.o!="C277A"),]

# Drop levels after subsetting
Rev.C277A = droplevels(Rev.C277A)


# subset focal strain C067A 
Rev.C067A=rbind(
  (D1[which(D1$strain.label=="C067A_C264A"),]), 
  (D1[which(D1$strain.label=="C067A_C395A"),]), 
  (D1[which(D1$strain.label=="C067A_C394B"),]), 
  (D1[which(D1$strain.label=="C067A_C412B"),]), 
  (D1[which(D1$strain.label=="C067A_C399B"),]), 
  (D1[which(D1$strain.label=="C066B_C067A"),]),
  (D1[which(D1$strain.label=="C067A_C277A"),]))

# remove other focal strains besides C067A from subset
Rev.C067A=Rev.C067A[-which(Rev.C067A$strain.o!="C067A"),]

# Drop levels after subsetting
Rev.C067A = droplevels(Rev.C067A)


# subset focal strain C399B 
Rev.C399B=rbind(
  (D1[which(D1$strain.label=="C399B_C264A"),]), 
  (D1[which(D1$strain.label=="C399B_C395A"),]), 
  (D1[which(D1$strain.label=="C399B_C394B"),]), 
  (D1[which(D1$strain.label=="C399B_C412B"),]), 
  (D1[which(D1$strain.label=="C067A_C399B"),]), 
  (D1[which(D1$strain.label=="C066B_C399B"),]),
  (D1[which(D1$strain.label=="C399B_C277A"),]))

# remove other focal strains besides C399B from subset
Rev.C399B=Rev.C399B[-which(Rev.C399B$strain.o!="C399B"),]

# Drop levels after subsetting
Rev.C399B = droplevels(Rev.C399B)

# subset focal strain C264A 
Rev.C264A=rbind(
  (D1[which(D1$strain.label=="C264A_C412B"),]), 
  (D1[which(D1$strain.label=="C277A_C264A"),]), 
  (D1[which(D1$strain.label=="C394B_C264A"),]), 
  (D1[which(D1$strain.label=="C395A_C264A"),]), 
  (D1[which(D1$strain.label=="C067A_C264A"),]), 
  (D1[which(D1$strain.label=="C066B_C264A"),]),
  (D1[which(D1$strain.label=="C399B_C264A"),]))

# remove other focal strains besides C264A from subset
Rev.C264A=Rev.C264A[-which(Rev.C264A$strain.o!="C264A"),]

# Drop levels after subsetting
Rev.C264A = droplevels(Rev.C264A)

# subset focal strain C412B 
Rev.C412B=rbind(
  (D1[which(D1$strain.label=="C264A_C412B"),]), 
  (D1[which(D1$strain.label=="C277A_C412B"),]), 
  (D1[which(D1$strain.label=="C394B_C412B"),]), 
  (D1[which(D1$strain.label=="C395A_C412B"),]), 
  (D1[which(D1$strain.label=="C067A_C412B"),]), 
  (D1[which(D1$strain.label=="C066B_C412B"),]),
  (D1[which(D1$strain.label=="C399B_C412B"),]))


# remove other focal strains besides C412B from subset
Rev.C412B=Rev.C412B[-which(Rev.C412B$strain.o!="C412B"),]

# Drop levels after subsetting
Rev.C412B = droplevels(Rev.C412B)

# nodule number and nodule proportion models for each focal strain  

# focal strain 394B nod number and nod prop models

# calc mean nod number of focal strain explained by competitor strain
Rn394B=aggregate(nod.sum ~ strain.comp.o, Rev.C394B, mean)
Rn394B$se = aggregate(nod.sum ~ strain.comp.o, FUN=function(x) sd(x)/sqrt(length(x)), data=Rev.C394B)[,2]
Rn394B

Rn394B$Strain = "394B"
colnames(Rn394B)[1]="Competitor"
Rn394B


# calc mean nod proportion of focal strain explained by competitor strain
Rp394B=aggregate(nod.prop ~ strain.comp.o, Rev.C394B, mean)
Rp394B$se = aggregate(nod.prop ~ strain.comp.o, FUN=function(x) sd(x)/sqrt(length(x)), data=Rev.C394B)[,2]
Rp394B

Rp394B$Strain = "394B"
colnames(Rp394B)[1]="Competitor"
Rp394B


# nod number model
Mn394B= lm(nod.sum ~ strain.comp.o, data = Rev.C394B)
summary(Mn394B)

# plot residuals 
par(mfrow=c(2,2))
plot(Mn394B)

# test fixed effect 
drop1(Mn394B, test = "F")

#              Df Sum of Sq   RSS    AIC F value  Pr(>F)  
#<none>                     22578 608.30                  
#strain.comp.o  6    3492.3 26071 612.41  2.7068 0.01747 *

# emmeans
Mn394B %>% emmeans(~ strain.comp.o) %>%
  contrast(., method = "pairwise") %>%
  summary(by = NULL, adjust = "holm") ->  temp #holm adjustment for multiple comparisons

# list just comparisons that are significant
temp[which(temp$p.value<=0.05),]

#[1] contrast estimate SE       df       t.ratio  p.value 
#<0 rows> (or 0-length row.names)


# nod prop model 
Mp394B= lm(nod.prop ~ strain.comp.o, data = Rev.C394B)
summary(Mp394B)

# plot residuals 
par(mfrow=c(2,2))
plot(Mp394B)

# test fixed effect 
drop1(Mp394B, test = "F")

#              Df Sum of Sq     RSS     AIC F value    Pr(>F)    
#<none>                     0.82099 -536.56                      
#strain.comp.o  6    1.6971 2.51808 -423.04  36.175 < 2.2e-16 **

# emmeans
Mp394B %>% emmeans(~ strain.comp.o) %>%
  contrast(., method = "pairwise") %>%
  summary(by = NULL, adjust = "holm") ->  temp #holm adjustment for multiple comparisons

# list just comparisons that are significant
temp[which(temp$p.value<=0.05),]

#        contrast   estimate        SE  df   t.ratio      p.value
#1  C066B - C067A -0.1094125 0.0312629 105 -3.499756 6.161757e-03
#2  C066B - C264A -0.2674807 0.0312629 105 -8.555852 1.665072e-12
#6  C066B - C412B -0.2688696 0.0312629 105 -8.600278 1.409890e-12
#7  C067A - C264A -0.1580682 0.0312629 105 -5.056096 2.187474e-05
#8  C067A - C277A  0.1529382 0.0312629 105  4.892004 3.977883e-05
#9  C067A - C395A  0.1181645 0.0312629 105  3.779704 2.610025e-03
#10 C067A - C399B  0.1008265 0.0312629 105  3.225116 1.343141e-02
#11 C067A - C412B -0.1594571 0.0312629 105 -5.100522 1.964447e-05
#12 C264A - C277A  0.3110064 0.0312629 105  9.948100 1.603142e-15
#13 C264A - C395A  0.2762327 0.0312629 105  8.835800 4.469037e-13
#14 C264A - C399B  0.2588947 0.0312629 105  8.281212 5.901845e-12
#18 C277A - C412B -0.3123953 0.0312629 105 -9.992527 1.337791e-15
#20 C395A - C412B -0.2776216 0.0312629 105 -8.880226 3.755667e-13
#21 C399B - C412B -0.2602836 0.0312629 105 -8.325638 5.045506e-12



# focal strain 066B nod number and nod prop models

# calc mean nod number of focal strain explained by competitor strain
Rn066B=aggregate(nod.sum ~ strain.comp.o, Rev.C066B, mean)
Rn066B$se = aggregate(nod.sum ~ strain.comp.o, FUN=function(x) sd(x)/sqrt(length(x)), data=Rev.C066B)[,2]
Rn066B

Rn066B$Strain = "066B"
colnames(Rn066B)[1]="Competitor"
Rn066B

# calc mean nod proportion of focal strain explained by competitor strain
Rp066B=aggregate(nod.prop ~ strain.comp.o, Rev.C066B, mean)
Rp066B$se = aggregate(nod.prop ~ strain.comp.o, FUN=function(x) sd(x)/sqrt(length(x)), data=Rev.C066B)[,2]
Rp066B

Rp066B$Strain = "066B"
colnames(Rp066B)[1]="Competitor"
Rp066B

# nod number model 
Mn066B= lm(nod.sum ~ strain.comp.o, data = Rev.C066B)
summary(Mn066B)

# plot residuals 
par(mfrow=c(2,2))
plot(Mn066B)

# test fixed effect 
drop1(Mn066B, test = "F")

#              Df Sum of Sq   RSS    AIC F value    Pr(>F)    
#<none>                     28523 634.47                      
#strain.comp.o  6     15924 44446 672.16  9.7699 1.539e-08 ***

# emmeans
Mn066B %>% emmeans(~ strain.comp.o) %>%
  contrast(., method = "pairwise") %>%
  summary(by = NULL, adjust = "holm") ->  temp #holm adjustment for multiple comparisons

# list just comparisons that are significant
temp[which(temp$p.value<=0.05),]

#        contrast  estimate       SE  df   t.ratio      p.value
#2  C067A - C277A  21.31250 5.827137 105  3.657457 6.001407e-03
#3  C067A - C394B  31.40625 5.827137 105  5.389654 8.729723e-06
#4  C067A - C395A  24.21875 5.827137 105  4.156201 1.125130e-03
#8  C264A - C394B  25.65625 5.827137 105  4.402891 4.639893e-04
#9  C264A - C395A  18.46875 5.827137 105  3.169438 2.802041e-02
#15 C277A - C412B -23.68750 5.827137 105 -4.065032 1.487560e-03
#18 C394B - C412B -33.78125 5.827137 105 -5.797230 1.502786e-06
#20 C395A - C412B -26.59375 5.827137 105 -4.563776 2.601597e-04
#21 C399B - C412B -18.03125 5.827137 105 -3.094358 3.286591e-02

# nod prop model
Mp066B= lm(nod.prop ~ strain.comp.o, data = Rev.C066B)
summary(Mp066B)

# plot residuals 
par(mfrow=c(2,2))
plot(Mp066B)

# test fixed effect 
drop1(Mp066B, test = "F")

#              Df Sum of Sq    RSS     AIC F value    Pr(>F)    
#<none>                     1.4438 -467.99                      
#strain.comp.o  6     7.429 8.8728 -278.45  89.185 < 2.2e-16 ***

# emmeans
Mp066B %>% emmeans(~ strain.comp.o) %>%
  contrast(., method = "pairwise") %>%
  summary(by = NULL, adjust = "holm") ->  temp #holm adjustment for multiple comparisons

# list just comparisons that are significant
temp[which(temp$p.value<=0.05),]

#        contrast   estimate         SE  df    t.ratio      p.value
#1  C067A - C264A -0.2577849 0.04165791 104  -6.188139 1.222045e-07
#2  C067A - C277A  0.2968997 0.04165791 104   7.127092 1.674217e-09
#3  C067A - C394B  0.4733454 0.04165791 104  11.362679 9.248056e-19
#4  C067A - C395A  0.1785266 0.04234652 104   4.215851 2.962089e-04
#5  C067A - C399B  0.2213189 0.04165791 104   5.312771 4.326770e-06
#6  C067A - C412B -0.2577849 0.04165791 104  -6.188139 1.222045e-07
#7  C264A - C277A  0.5546847 0.04165791 104  13.315231 6.284518e-23
#8  C264A - C394B  0.7311304 0.04165791 104  17.550818 1.533771e-31
#9  C264A - C395A  0.4363115 0.04234652 104  10.303364 1.977341e-16
#10 C264A - C399B  0.4791039 0.04165791 104  11.500911 5.176470e-19
#12 C277A - C394B  0.1764457 0.04165791 104   4.235587 2.962089e-04
#13 C277A - C395A -0.1183731 0.04234652 104  -2.795346 2.470134e-02
#15 C277A - C412B -0.5546847 0.04165791 104 -13.315231 6.284518e-23
#16 C394B - C395A -0.2948188 0.04234652 104  -6.962057 3.422509e-09
#17 C394B - C399B -0.2520265 0.04165791 104  -6.049908 1.849391e-07
#18 C394B - C412B -0.7311304 0.04165791 104 -17.550818 1.533771e-31
#20 C395A - C412B -0.4363115 0.04234652 104 -10.303364 1.977341e-16
#21 C399B - C412B -0.4791039 0.04165791 104 -11.500911 5.176470e-19


# focal strain 395A nod number and nod prop models

# calc mean nod number of focal strain explained by competitor
Rn395A=aggregate(nod.sum ~ strain.comp.o, Rev.C395A, mean)
Rn395A$se = aggregate(nod.sum ~ strain.comp.o, FUN=function(x) sd(x)/sqrt(length(x)), data=Rev.C395A)[,2]
Rn395A

Rn395A$Strain = "395A"
colnames(Rn395A)[1]="Competitor"
Rn395A

# calc mean nod proportion of focal strain explained by competitor strain 
Rp395A=aggregate(nod.prop ~ strain.comp.o, Rev.C395A, mean)
Rp395A$se = aggregate(nod.prop ~ strain.comp.o, FUN=function(x) sd(x)/sqrt(length(x)), data=Rev.C395A)[,2]
Rp395A

Rp395A$Strain = "395A"
colnames(Rp395A)[1]="Competitor"
Rp395A

# nod number model
Mn395A= lm(nod.sum ~ strain.comp.o, data = Rev.C395A)
summary(Mn395A)

# plot residuals 
par(mfrow=c(2,2))
plot(Mn395A)

# test fixed effect 
drop1(Mn395A, test = "F")

#              Df Sum of Sq   RSS    AIC F value   Pr(>F)   
#<none>                     14695 560.20                    
#strain.comp.o  6    3421.9 18117 571.65  4.0749 0.001029 **

# emmeans 
Mn395A %>% emmeans(~ strain.comp.o) %>%
  contrast(., method = "pairwise") %>%
  summary(by = NULL, adjust = "holm") ->  temp #holm adjustment for multiple comparisons

# list just comparisons that are significant
temp[which(temp$p.value<=0.05),]

#        contrast  estimate       SE  df   t.ratio     p.value
#2  C066B - C264A -15.03125 4.182633 105 -3.593729 0.009462785
#12 C264A - C277A  15.62500 4.182633 105  3.735685 0.006094505
#13 C264A - C394B  16.12500 4.182633 105  3.855227 0.004190473

# nod prop model 
Mp395A= lm(nod.prop ~ strain.comp.o, data = Rev.C395A)
summary(Mp395A)

# plot residuals 
par(mfrow=c(2,2))
plot(Mp395A)

# test fixed effect 
drop1(Mp395A, test = "F")

#              Df Sum of Sq    RSS     AIC F value    Pr(>F)    
#<none>                     1.9709 -433.45                      
#strain.comp.o  6    7.5964 9.5673 -270.08  66.807 < 2.2e-16 ***

# emmeans
Mp395A %>% emmeans(~ strain.comp.o) %>%
  contrast(., method = "pairwise") %>%
  summary(by = NULL, adjust = "holm") ->  temp #holm adjustment for multiple comparisons

# list just comparisons that are significant
temp[which(temp$p.value<=0.05),]

#contrast   estimate         SE  df    t.ratio      p.value
#1  C066B - C067A -0.2927386 0.04947552 104  -5.916837 4.245111e-07
#2  C066B - C264A -0.5532718 0.04947552 104 -11.182738 2.475003e-18
#4  C066B - C394B  0.1586900 0.04947552 104   3.207444 1.151590e-02
#5  C066B - C399B -0.1354173 0.04947552 104  -2.737057 2.187760e-02
#6  C066B - C412B -0.5576121 0.04947552 104 -11.270464 1.678992e-18
#7  C067A - C264A -0.2605332 0.04867099 104  -5.352947 4.156367e-06
#8  C067A - C277A  0.3005550 0.04867099 104   6.175239 1.556745e-07
#9  C067A - C394B  0.4514285 0.04867099 104   9.275106 4.197652e-14
#10 C067A - C399B  0.1573213 0.04867099 104   3.232342 1.151590e-02
#11 C067A - C412B -0.2648735 0.04867099 104  -5.442123 3.172387e-06
#12 C264A - C277A  0.5610882 0.04867099 104  11.528186 4.769381e-19
#13 C264A - C394B  0.7119617 0.04867099 104  14.628052 1.092453e-25
#14 C264A - C399B  0.4178545 0.04867099 104   8.585289 1.238487e-12
#16 C277A - C394B  0.1508736 0.04867099 104   3.099867 1.245580e-02
#17 C277A - C399B -0.1432337 0.04867099 104  -2.942897 1.604028e-02
#18 C277A - C412B -0.5654284 0.04867099 104 -11.617362 3.196089e-19
#19 C394B - C399B -0.2941073 0.04867099 104  -6.042764 2.627645e-07
#20 C394B - C412B -0.7163020 0.04867099 104 -14.717228 7.475871e-26
#21 C399B - C412B -0.4221947 0.04867099 104  -8.674464 8.468695e-13



# focal strain 277A nod number and nod prop models

# calc mean nod number of focal strain explained by competitor strain 
Rn277A=aggregate(nod.sum ~ strain.comp.o, Rev.C277A, mean)
Rn277A$se = aggregate(nod.sum ~ strain.comp.o, FUN=function(x) sd(x)/sqrt(length(x)), data=Rev.C277A)[,2]
Rn277A

Rn277A$Strain = "277A"
colnames(Rn277A)[1]="Competitor"
Rn277A

# calc mean nod proportion of focal strain explained by competitor strain 
Rp277A=aggregate(nod.prop ~ strain.comp.o, Rev.C277A, mean)
Rp277A$se = aggregate(nod.prop ~ strain.comp.o, FUN=function(x) sd(x)/sqrt(length(x)), data=Rev.C277A)[,2]
Rp277A

Rp277A$Strain = "277A"
colnames(Rp277A)[1]="Competitor"
Rp277A

# nod number model
Mn277A= lm(nod.sum ~ strain.comp.o, data = Rev.C277A)
summary(Mn277A)

# plot residuals 
par(mfrow=c(2,2))
plot(Mn277A)

# test fixed effect 
drop1(Mn277A, test = "F")

# emmeans
Mn277A %>% emmeans(~ strain.comp.o) %>%
  contrast(., method = "pairwise") %>%
  summary(by = NULL, adjust = "holm") ->  temp #holm adjustment for multiple comparisons

# list just comparisons that are significant
temp[which(temp$p.value<=0.05),]

#        contrast estimate       SE  df   t.ratio      p.value
#2  C066B - C264A -16.4375 5.404623 105 -3.041378 4.758860e-02
#6  C066B - C412B -18.9375 5.404623 105 -3.503945 1.147635e-02
#12 C264A - C394B  24.6875 5.404623 105  4.567849 2.694550e-04
#13 C264A - C395A  22.1250 5.404623 105  4.093717 1.504558e-03
#18 C394B - C412B -27.1875 5.404623 105 -5.030415 4.264766e-05
#20 C395A - C412B -24.6250 5.404623 105 -4.556284 2.694550e-04

# nod prop model 
Mp277A= lm(nod.prop ~ strain.comp.o, data = Rev.C277A)
summary(Mp277A)

# plot residuals 
par(mfrow=c(2,2))
plot(Mp277A)

# test fixed effect 
drop1(Mp277A, test = "F")

#              Df Sum of Sq    RSS     AIC F value    Pr(>F)    
#<none>                     1.4283 -474.55                      
#strain.comp.o  6    6.1602 7.5885 -299.49  75.476 < 2.2e-16 ***

# emmeans
Mp277A %>% emmeans(~ strain.comp.o) %>%
  contrast(., method = "pairwise") %>%
  summary(by = NULL, adjust = "holm") ->  temp #holm adjustment for multiple comparisons

# list just comparisons that are significant
temp[which(temp$p.value<=0.05),]

#        contrast   estimate         SE  df    t.ratio      p.value
#1  C066B - C067A -0.2294654 0.04123547 105  -5.564759 1.820752e-06
#2  C066B - C264A -0.4453153 0.04123547 105 -10.799327 1.773457e-17
#3  C066B - C394B  0.2422894 0.04123547 105   5.875752 5.015146e-07
#6  C066B - C412B -0.4453153 0.04123547 105 -10.799327 1.773457e-17
#7  C067A - C264A -0.2158499 0.04123547 105  -5.234569 6.827831e-06
#8  C067A - C394B  0.4717548 0.04123547 105  11.440511 6.909003e-19
#9  C067A - C395A  0.2126453 0.04123547 105   5.156854 7.138426e-06
#10 C067A - C399B  0.1891807 0.04123547 105   4.587814 6.221630e-05
#11 C067A - C412B -0.2158499 0.04123547 105  -5.234569 6.827831e-06
#12 C264A - C394B  0.6876047 0.04123547 105  16.675079 5.910296e-30
#13 C264A - C395A  0.4284952 0.04123547 105  10.391423 1.295796e-16
#14 C264A - C399B  0.4050305 0.04123547 105   9.822383 2.149694e-15
#16 C394B - C395A -0.2591095 0.04123547 105  -6.283657 8.414819e-08
#17 C394B - C399B -0.2825741 0.04123547 105  -6.852697 6.131630e-09
#18 C394B - C412B -0.6876047 0.04123547 105 -16.675079 5.910296e-30
#20 C395A - C412B -0.4284952 0.04123547 105 -10.391423 1.295796e-16
#21 C399B - C412B -0.4050305 0.04123547 105  -9.822383 2.149694e-15


# focal strain 067A nod number and nod prop models

# calc mean nod number of focal strain explained by competitor strain
Rn067A=aggregate(nod.sum ~ strain.comp.o, Rev.C067A, mean)
Rn067A$se = aggregate(nod.sum ~ strain.comp.o, FUN=function(x) sd(x)/sqrt(length(x)), data=Rev.C067A)[,2]
Rn067A

Rn067A$Strain = "067A"
colnames(Rn067A)[1]="Competitor"
Rn067A

# calc mean nod proportion of focal strain explained by competitor strain
Rp067A=aggregate(nod.prop ~ strain.comp.o, Rev.C067A, mean)
Rp067A$se = aggregate(nod.prop ~ strain.comp.o, FUN=function(x) sd(x)/sqrt(length(x)), data=Rev.C067A)[,2]
Rp067A

Rp067A$Strain = "067A"
colnames(Rp067A)[1]="Competitor"
Rp067A


# nod number model
Mn067A= lm(nod.sum ~ strain.comp.o, data = Rev.C067A)
summary(Mn067A)

# plot residuals 
par(mfrow=c(2,2))
plot(Mn067A)

# test fixed effect 
drop1(Mn067A, test = "F")

#              Df Sum of Sq   RSS    AIC F value    Pr(>F)    
#<none>                     36431 661.88                      
#strain.comp.o  6     21384 57815 701.61  10.272 6.273e-09 ***
  
# emmeans
Mn067A %>% emmeans(~ strain.comp.o) %>%
  contrast(., method = "pairwise") %>%
  summary(by = NULL, adjust = "holm") ->  temp #holm adjustment for multiple comparisons

# list just comparisons that are significant
temp[which(temp$p.value<=0.05),]

#        contrast  estimate       SE  df   t.ratio      p.value
#1  C066B - C264A -25.75000 6.585594 105 -3.910050 1.966286e-03
#6  C066B - C412B -27.62500 6.585594 105 -4.194762 8.013721e-04
#7  C264A - C277A  30.06250 6.585594 105  4.564888 2.181161e-04
#8  C264A - C394B  31.87500 6.585594 105  4.840110 8.184669e-05
#9  C264A - C395A  31.59375 6.585594 105  4.797403 9.072368e-05
#10 C264A - C399B  27.15625 6.585594 105  4.123584 9.721660e-04
#15 C277A - C412B -31.93750 6.585594 105 -4.849600 8.184669e-05
#18 C394B - C412B -33.75000 6.585594 105 -5.124822 2.862826e-05
#20 C395A - C412B -33.46875 6.585594 105 -5.082115 3.266831e-05
#21 C399B - C412B -29.03125 6.585594 105 -4.408296 3.786107e-04

# nod prop model
Mp067A= lm(nod.prop ~ strain.comp.o, data = Rev.C067A)
summary(Mp067A)

# plot residuals 
par(mfrow=c(2,2))
plot(Mp067A)

# test fixed effect
drop1(Mp067A, test = "F")

#              Df Sum of Sq     RSS     AIC F value    Pr(>F)    
#<none>                      0.9884 -515.78                      
#strain.comp.o  6    13.271 14.2589 -228.84  234.97 < 2.2e-16 **

# emmeans
Mp067A %>% emmeans(~ strain.comp.o) %>%
  contrast(., method = "pairwise") %>%
  summary(by = NULL, adjust = "holm") ->  temp #holm adjustment for multiple comparisons

# list just comparisons that are significant
temp[which(temp$p.value<=0.05),]

#        contrast    estimate         SE  df    t.ratio      p.value
#1  C066B - C264A -0.74221505 0.03430179 105 -21.637796 2.920018e-39
#3  C066B - C394B  0.09832785 0.03430179 105   2.866552 4.513323e-02
#6  C066B - C412B -0.73874542 0.03430179 105 -21.536646 4.109621e-39
#7  C264A - C277A  0.78415011 0.03430179 105  22.860329 2.768782e-41
#8  C264A - C394B  0.84054290 0.03430179 105  24.504348 6.521619e-44
#9  C264A - C395A  0.72905012 0.03430179 105  21.253999 1.193550e-38
#10 C264A - C399B  0.69297995 0.03430179 105  20.202445 7.551595e-37
#15 C277A - C412B -0.78068048 0.03430179 105 -22.759179 3.867096e-41
#16 C394B - C395A -0.11149278 0.03430179 105  -3.250349 1.549339e-02
#17 C394B - C399B -0.14756296 0.03430179 105  -4.301903 4.188224e-04
#18 C394B - C412B -0.83707327 0.03430179 105 -24.403198 8.992459e-44
#20 C395A - C412B -0.72558049 0.03430179 105 -21.152849 1.673447e-38
#21 C399B - C412B -0.68951031 0.03430179 105 -20.101295 1.060670e-36



# focal strain 399B nod number and nod prop models

# calc mean nod number of focal strain explained by competitor strain
Rn399B=aggregate(nod.sum ~ strain.comp.o, Rev.C399B, mean)
Rn399B$se = aggregate(nod.sum ~ strain.comp.o, FUN=function(x) sd(x)/sqrt(length(x)), data=Rev.C399B)[,2]
Rn399B

Rn399B$Strain = "399B"
colnames(Rn399B)[1]="Competitor"
Rn399B

# calc mean nod proportion of focal strain explained by competitor strain
Rp399B=aggregate(nod.prop ~ strain.comp.o, Rev.C399B, mean)
Rp399B$se = aggregate(nod.prop ~ strain.comp.o, FUN=function(x) sd(x)/sqrt(length(x)), data=Rev.C399B)[,2]
Rp399B

Rp399B$Strain = "399B"
colnames(Rp399B)[1]="Competitor"
Rp399B

# nod number model
Mn399B= lm(nod.sum ~ strain.comp.o, data = Rev.C399B)
summary(Mn399B)

# plot residuals 
par(mfrow=c(2,2))
plot(Mn399B)

# test fixed effect 
drop1(Mn399B, test = "F")

#              Df Sum of Sq   RSS    AIC F value    Pr(>F)    
#<none>                     21104 600.73                      
#strain.comp.o  6     10752 31856 634.85  8.9165 7.299e-08 ***
  
# emmeans
Mn399B %>% emmeans(~ strain.comp.o) %>%
  contrast(., method = "pairwise") %>%
  summary(by = NULL, adjust = "holm") ->  temp #holm adjustment for multiple comparisons

# list just comparisons that are significant
temp[which(temp$p.value<=0.05),]

#        contrast  estimate      SE  df   t.ratio      p.value
#6  C066B - C412B -15.90625 5.01231 105 -3.173437 3.360042e-02
#13 C264A - C394B  25.43750 5.01231 105  5.075005 3.198073e-05
#14 C264A - C395A  24.31250 5.01231 105  4.850558 7.723392e-05
#18 C277A - C412B -15.84375 5.01231 105 -3.160968 3.360042e-02
#20 C394B - C412B -27.21875 5.01231 105 -5.430380 7.673487e-06
#21 C395A - C412B -26.09375 5.01231 105 -5.205933 1.929792e-05

# nodule prop model 
Mp399B= lm(nod.prop ~ strain.comp.o, data = Rev.C399B)
summary(Mp399B)

# plot residuals 
par(mfrow=c(2,2))
plot(Mp399B)

# test fixed effect 
drop1(Mp399B, test = "F")

#              Df Sum of Sq    RSS     AIC F value    Pr(>F)    
#<none>                     1.3251 -482.95                      
#strain.comp.o  6    8.3495 9.6745 -272.29  110.27 < 2.2e-16 ***

# emmeans
Mp399B %>% emmeans(~ strain.comp.o) %>%
  contrast(., method = "pairwise") %>%
  summary(by = NULL, adjust = "holm") ->  temp #holm adjustment for multiple comparisons

# list just comparisons that are significant
temp[which(temp$p.value<=0.05),]

#        contrast   estimate         SE  df    t.ratio      p.value
#1  C066B - C067A -0.2138761 0.03971709 105  -5.384988 3.118121e-06
#2  C066B - C264A -0.5186836 0.03971709 105 -13.059455 1.479385e-22
#4  C066B - C394B  0.2188203 0.03971709 105   5.509475 2.065897e-06
#6  C066B - C412B -0.5183451 0.03971709 105 -13.050933 1.479385e-22
#7  C067A - C264A -0.3048075 0.03971709 105  -7.674467 1.076880e-10
#8  C067A - C277A  0.2879494 0.03971709 105   7.250013 7.363231e-10
#9  C067A - C394B  0.4326964 0.03971709 105  10.894463 7.841731e-18
#10 C067A - C395A  0.2647088 0.03971709 105   6.664859 1.134412e-08
#11 C067A - C412B -0.3044690 0.03971709 105  -7.665945 1.076880e-10
#12 C264A - C277A  0.5927569 0.03971709 105  14.924480 1.986887e-26
#13 C264A - C394B  0.7375039 0.03971709 105  18.568930 1.252228e-33
#14 C264A - C395A  0.5695163 0.03971709 105  14.339326 3.012742e-25
#16 C277A - C394B  0.1447470 0.03971709 105   3.644451 2.092335e-03
#18 C277A - C412B -0.5924184 0.03971709 105 -14.915958 1.986887e-26
#19 C394B - C395A -0.1679876 0.03971709 105  -4.229604 3.010084e-04
#20 C394B - C412B -0.7371654 0.03971709 105 -18.560409 1.252228e-33
#21 C395A - C412B -0.5691778 0.03971709 105 -14.330804 3.012742e-25


# focal strain 264A nod number and nod prop models

# calc mean nod number of focal strain explained by competitor strain
Rn264A=aggregate(nod.sum ~ strain.comp.o, Rev.C264A, mean)
Rn264A$se = aggregate(nod.sum ~ strain.comp.o, FUN=function(x) sd(x)/sqrt(length(x)), data=Rev.C264A)[,2]
Rn264A

Rn264A$Strain = "264A"
colnames(Rn264A)[1]="Competitor"
Rn264A

# calc mean nod proportion of focal strain explained by competitor strain
Rp264A=aggregate(nod.prop ~ strain.comp.o, Rev.C264A, mean)
Rp264A$se = aggregate(nod.prop ~ strain.comp.o, FUN=function(x) sd(x)/sqrt(length(x)), data=Rev.C264A)[,2]
Rp264A

Rp264A$Strain = "264A"
colnames(Rp264A)[1]="Competitor"
Rp264A

# nod number model
Mn264A= lm(nod.sum ~ strain.comp.o, data = Rev.C264A)
summary(Mn264A)

# plot residuals 
par(mfrow=c(2,2))
plot(Mn264A)

# test fixed effect  
drop1(Mn264A, test = "F")

#              Df Sum of Sq    RSS     AIC F value  Pr(>F)  
#<none>                     19.094 -184.14                  
#strain.comp.o  6    2.1205 21.214 -184.35  1.9435 0.08051 

# emmeans
Mn264A %>% emmeans(~ strain.comp.o) %>%
  contrast(., method = "pairwise") %>%
  summary(by = NULL, adjust = "holm") ->  temp #holm adjustment for multiple comparisons

# list just comparisons that are significant
temp[which(temp$p.value<=0.05),]

#[1] contrast estimate SE       df       t.ratio  p.value 
#<0 rows> (or 0-length row.names)

# nod prop model
Mp264A= lm(nod.prop ~ strain.comp.o, data = Rev.C264A)
summary(Mp264A)

# plot residuals 
par(mfrow=c(2,2))
plot(Mp264A)

# test fixed effect  
drop1(Mp264A, test = "F")

#              Df Sum of Sq    RSS     AIC F value    Pr(>F)    
#<none>                     1.3898 -434.78                      
#strain.comp.o  6   0.70222 2.0920 -404.25  8.1685 3.707e-07 ***
  
# emmeans
Mp264A %>% emmeans(~ strain.comp.o) %>%
  contrast(., method = "pairwise") %>%
  summary(by = NULL, adjust = "holm") ->  temp #holm adjustment for multiple comparisons

# list just comparisons that are significant
temp[which(temp$p.value<=0.05),]

#        contrast   estimate         SE df   t.ratio      p.value
#6  C066B - C412B -0.3104167 0.05183109 97 -5.989006 7.505168e-07
#11 C067A - C412B -0.3104167 0.05183109 97 -5.989006 7.505168e-07
#15 C277A - C412B -0.3104167 0.05183109 97 -5.989006 7.505168e-07
#18 C394B - C412B -0.3090278 0.05183109 97 -5.962209 7.505168e-07
#20 C395A - C412B -0.3000000 0.05183109 97 -5.788032 1.399708e-06
#21 C399B - C412B -0.3082041 0.05183109 97 -5.946318 7.505168e-07



# focal strain 412B nod number and nod prop models

# calc mean nod number of focal strain explained by competitor strain
Rn412B=aggregate(nod.sum ~ strain.comp.o, Rev.C412B, mean)
Rn412B$se = aggregate(nod.sum ~ strain.comp.o, FUN=function(x) sd(x)/sqrt(length(x)), data=Rev.C412B)[,2]
Rn412B

Rn412B$Strain = "412B"
colnames(Rn412B)[1]="Competitor"
Rn412B

# calc mean nod proportion of focal strain explained by competitor strain
Rp412B=aggregate(nod.prop ~ strain.comp.o, Rev.C412B, mean)
Rp412B$se = aggregate(nod.prop ~ strain.comp.o, FUN=function(x) sd(x)/sqrt(length(x)), data=Rev.C412B)[,2]
Rp412B

Rp412B$Strain = "412B"
colnames(Rp412B)[1]="Competitor"
Rp412B

# nod number model
Mn412B= lm(nod.sum ~ strain.comp.o, data = Rev.C412B)
summary(Mn412B)

# plot residuals 
par(mfrow=c(2,2))
plot(Mn412B)

# test fixed effect 
drop1(Mn412B, test = "F")

#              Df Sum of Sq    RSS    AIC F value    Pr(>F)    
#<none>                     172.91 62.636                      
#strain.comp.o  6    49.308 222.21 78.736  4.9905 0.0001557 ***

# emmeans
Mn412B %>% emmeans(~ strain.comp.o) %>%
  contrast(., method = "pairwise") %>%
  summary(by = NULL, adjust = "holm") ->  temp #holm adjustment for multiple comparisons

# list just comparisons that are significant
temp[which(temp$p.value<=0.05),]

#        contrast estimate        SE  df   t.ratio      p.value
#2  C066B - C264A -1.96875 0.4536968 105 -4.339352 0.0006923385
#7  C067A - C264A -1.78125 0.4536968 105 -3.926080 0.0026287113
#12 C264A - C277A  1.96875 0.4536968 105  4.339352 0.0006923385
#13 C264A - C394B  1.96875 0.4536968 105  4.339352 0.0006923385
#14 C264A - C395A  1.75000 0.4536968 105  3.857202 0.0031702601
#15 C264A - C399B  1.84375 0.4536968 105  4.063837 0.0016809227

# nod prop model
Mp412B= lm(nod.prop ~ strain.comp.o, data = Rev.C412B)
summary(Mp412B)

# plot residuals 
par(mfrow=c(2,2))
plot(Mp412B)

# test fixed effect 
drop1(Mp412B, test = "F")

#              Df Sum of Sq    RSS     AIC F value    Pr(>F)    
#<none>                     1.3746 -435.93                      
#strain.comp.o  6    3.4916 4.8661 -316.46  41.065 < 2.2e-16 ***

# emmeans
Mp412B %>% emmeans(~ strain.comp.o) %>%
  contrast(., method = "pairwise") %>%
  summary(by = NULL, adjust = "holm") ->  temp #holm adjustment for multiple comparisons

# list just comparisons that are significant
temp[which(temp$p.value<=0.05),]

#        contrast   estimate         SE df   t.ratio      p.value
#2  C066B - C264A -0.6895833 0.05154629 97 -13.37794 1.986221e-22
#7  C067A - C264A -0.6861137 0.05154629 97 -13.31063 2.210888e-22
#12 C264A - C277A  0.6895833 0.05154629 97  13.37794 1.986221e-22
#13 C264A - C394B  0.6895833 0.05154629 97  13.37794 1.986221e-22
#14 C264A - C395A  0.6835069 0.05154629 97  13.26006 2.644176e-22
#15 C264A - C399B  0.6870323 0.05154629 97  13.32845 2.151548e-22


#** FIG S5 a-h: nodule number  ####

# 394B nod number plot
Rn394B$Competitor=factor(Rn394B$Competitor, levels=c("C412B", "C264A", "C399B", "C067A", "C277A", "C395A", "C066B"))

ggplot(Rn394B, aes(x= Competitor, y= nod.sum)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=nod.sum-se, ymax = nod.sum+se), width=0, size=1, position=position_dodge(.9)) +
  geom_hline(aes(yintercept = mean(nod.sum)), linetype='dashed', size=1, color="red") +
  xlab('Rhizobium Competitor') + ylab('Focal strain
nodule number') +
  ylim(0,52)+
  scale_x_discrete(labels=c('412B','264A','399B','067A','277A','395A','066B')) +
  theme(panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.border = element_blank(),axis.line = element_line(colour = "black"),axis.ticks = element_line(colour = "black", size = 0.5),
        axis.ticks.length = unit(0.2, "cm"),text = element_text(size=32), axis.text.x = element_text(size=28, vjust=0.5, angle=90), 
        axis.text.y = element_text(size=28,hjust=0.95), plot.title = element_text(hjust = 0.5,size=28)) +
  ggtitle("Focal strain 394B")

# 066B nod number plot
Rn066B$Competitor=factor(Rn066B$Competitor, levels=c("C412B", "C264A", "C399B","C067A", "C277A", "C395A", "C394B"))

ggplot(Rn066B, aes(x= Competitor, y= nod.sum)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=nod.sum-se, ymax = nod.sum+se), width=0, size=1, position=position_dodge(.9)) +
  geom_hline(aes(yintercept = mean(nod.sum)), linetype='dashed', size=1, color="red") +
  xlab('Rhizobium Competitor') + ylab('Focal strain
nodule number') +
  ylim(0,52)+
  scale_x_discrete(labels=c('412B','264A','399B','067A','277A','395A','394B')) +
  theme(panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.border = element_blank(),axis.line = element_line(colour = "black"),axis.ticks = element_line(colour = "black", size = 0.5),
        axis.ticks.length = unit(0.2, "cm"),text = element_text(size=32), axis.text.x = element_text(size=28, vjust=0.5, angle=90), 
        axis.text.y = element_text(size=28,hjust=0.95), plot.title = element_text(hjust = 0.5, size=28)) +
  annotate("text", x = 1, y = 52, label = "a", size = 7) +
  annotate("text", x = 2, y = 45, label = "abc", size = 7) +
  annotate("text", x = 3, y = 32, label = "bcd", size = 7) +
  annotate("text", x = 4, y = 52, label = "ab", size = 7) +
  annotate("text", x = 5, y = 25, label = "cd", size = 7) +
  annotate("text", x = 6, y = 23, label = "d", size = 7) +
  annotate("text", x = 7, y = 15, label = "d", size = 7) +
  ggtitle("Focal strain 066B")


# 395A nod number plot 
Rn395A$Competitor=factor(Rn395A$Competitor, levels=c("C412B", "C264A", "C399B", "C067A", "C277A", "C066B", "C394B"))

ggplot(Rn395A, aes(x= Competitor, y= nod.sum)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=nod.sum-se, ymax = nod.sum+se), width=0, size=1, position=position_dodge(.9)) +
  geom_hline(aes(yintercept = mean(nod.sum)), linetype='dashed', size=1, color="red") +
  xlab('Rhizobium Competitor') + ylab('Focal strain
nodule number') +
  ylim(0,52) +
  scale_x_discrete(labels=c('412B','264A', '399B','067A','277A','066B','394B')) +
  theme(panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.border = element_blank(),axis.line = element_line(colour = "black"),axis.ticks = element_line(colour = "black", size = 0.5),
        axis.ticks.length = unit(0.2, "cm"),text = element_text(size=32), axis.text.x = element_text(size=28, vjust=0.5, angle=90), 
        axis.text.y = element_text(size=28,hjust=0.95), plot.title = element_text(hjust = 0.5, size=28)) +
  annotate("text", x = 1, y = 25, label = "ab", size = 7) +
  annotate("text", x = 2, y = 34, label = "a", size = 7) +
  annotate("text", x = 3, y = 23, label = "ab", size = 7) +
  annotate("text", x = 4, y = 28, label = "ab", size = 7) +
  annotate("text", x = 5, y = 20, label = "b", size = 7) +
  annotate("text", x = 6, y = 20, label = "b", size = 7) +
  annotate("text", x = 7, y = 20, label = "b", size = 7) +
  ggtitle("Focal strain 395A")


# 277A nod number plot
Rn277A$Competitor=factor(Rn277A$Competitor, levels=c("C412B", "C264A", "C399B", "C067A", "C395A", "C066B", "C394B"))

ggplot(Rn277A, aes(x= Competitor, y= nod.sum)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=nod.sum-se, ymax = nod.sum+se), width=0, size=1, position=position_dodge(.9)) +
  geom_hline(aes(yintercept = mean(nod.sum)), linetype='dashed', size=1, color="red") +
  xlab('Rhizobium Competitor') + ylab('Focal strain
nodule number') +
  ylim(0,52)+
  scale_x_discrete(labels=c('412B','264A', '399B','067A','395A','066B','394B')) +
  theme(panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.border = element_blank(),axis.line = element_line(colour = "black"),axis.ticks = element_line(colour = "black", size = 0.5),
        axis.ticks.length = unit(0.2, "cm"),text = element_text(size=32), axis.text.x = element_text(size=28, vjust=0.5, angle=90), 
        axis.text.y = element_text(size=28,hjust=0.95), plot.title = element_text(hjust = 0.5, size=28)) +
  annotate("text", x = 1, y = 50, label = "a", size = 7) +
  annotate("text", x = 2, y = 47, label = "a", size = 7) +
  annotate("text", x = 3, y = 37, label = "ab", size = 7) +
  annotate("text", x = 4, y = 37, label = "ab", size = 7) +
  annotate("text", x = 5, y = 22, label = "b", size = 7) +
  annotate("text", x = 6, y = 30, label = "b", size = 7) +
  annotate("text", x = 7, y = 20, label = "b", size = 7) +
  ggtitle("Focal strain 277A")



# 067A nod number plot
Rn067A$Competitor=factor(Rn067A$Competitor, levels=c("C412B", "C264A", "C399B", "C277A", "C395A", "C066B", "C394B"))

ggplot(Rn067A, aes(x= Competitor, y= nod.sum)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=nod.sum-se, ymax = nod.sum+se), width=0, size=1, position=position_dodge(.9)) +
  geom_hline(aes(yintercept = mean(nod.sum)), linetype='dashed', size=1, color="red") +
  xlab('Rhizobium Competitor') + ylab('Focal strain
nodule number') +
  ylim(0,52)+
  scale_x_discrete(labels=c('412B','264A','399B','277A','395A','066B','394B')) +
  theme(panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.border = element_blank(),axis.line = element_line(colour = "black"),axis.ticks = element_line(colour = "black", size = 0.5),
        axis.ticks.length = unit(0.2, "cm"),text = element_text(size=32), axis.text.x = element_text(size=28, vjust=0.5, angle=90), 
        axis.text.y = element_text(size=28,hjust=0.95), plot.title = element_text(hjust = 0.5, size=28)) +
  annotate("text", x = 1, y = 50, label = "a", size = 7) +
  annotate("text", x = 2, y = 52, label = "a", size = 7) +
  annotate("text", x = 3, y = 22, label = "b", size = 7) +
  annotate("text", x = 4, y = 22, label = "b", size = 7) +
  annotate("text", x = 5, y = 22, label = "b", size = 7) +
  annotate("text", x = 6, y = 22, label = "b", size = 7) +
  annotate("text", x = 7, y = 22, label = "b", size = 7) +
  ggtitle("Focal strain 067A")


# 399B nod number plot
Rn399B$Competitor=factor(Rn399B$Competitor, levels=c("C412B", "C264A", "C067A", "C277A", "C395A", "C066B", "C394B"))

ggplot(Rn399B, aes(x= Competitor, y= nod.sum)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=nod.sum-se, ymax = nod.sum+se), width=0, size=1, position=position_dodge(.9)) +
  geom_hline(aes(yintercept = mean(nod.sum)), linetype='dashed', size=1, color="red") +
  xlab('Rhizobium Competitor') + ylab('Focal strain
nodule number') +
  ylim(0,52) +
  scale_x_discrete(labels=c('412B','264A','067A','277A','395A','066B','394B')) +
  theme(panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.border = element_blank(),axis.line = element_line(colour = "black"),axis.ticks = element_line(colour = "black", size = 0.5),
        axis.ticks.length = unit(0.2, "cm"),text = element_text(size=32), axis.text.x = element_text(size=28, vjust=0.5, angle=90), 
        axis.text.y = element_text(size=28,hjust=0.95), plot.title = element_text(hjust = 0.5, size=28)) +
  annotate("text", x = 1, y = 45, label = "a", size = 7) +
  annotate("text", x = 2, y = 45, label = "ab", size = 7) +
  annotate("text", x = 3, y = 30, label = "abc", size = 7) +
  annotate("text", x = 4, y = 30, label = "bc", size = 7) +
  annotate("text", x = 5, y = 18, label = "c", size = 7) +
  annotate("text", x = 6, y = 30, label = "bc", size = 7) +
  annotate("text", x = 7, y = 18, label = "c", size = 7) +
  ggtitle("Focal strain 399B")


# 264A nod number plot
Rn264A$Competitor=factor(Rn264A$Competitor, levels=c("C412B", "C399B","C067A", "C277A", "C395A","C066B", "C394B"))

ggplot(Rn264A, aes(x= Competitor, y= nod.sum)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=nod.sum-se, ymax = nod.sum+se), width=0, size=1, position=position_dodge(.9)) +
  geom_hline(aes(yintercept = mean(nod.sum)), linetype='dashed', size=1, color="red") +
  xlab('Rhizobium Competitor') + ylab('Focal strain
nodule number') +
  ylim(0,52)+
  scale_x_discrete(labels=c('412B','399B','067A','277A','395A','066B','394B')) +
  theme(panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.border = element_blank(),axis.line = element_line(colour = "black"),axis.ticks = element_line(colour = "black", size = 0.5),
        axis.ticks.length = unit(0.2, "cm"),text = element_text(size=32), axis.text.x = element_text(size=28, vjust=0.5, angle=90), 
        axis.text.y = element_text(size=28,hjust=0.95), plot.title = element_text(hjust = 0.5, size=28)) +
  ggtitle("Focal strain 264A")


# 412B nod number plot
Rn412B$Competitor=factor(Rn412B$Competitor, levels=c("C264A", "C399B","C067A", "C277A", "C395A","C066B", "C394B"))

ggplot(Rn412B, aes(x= Competitor, y= nod.sum)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=nod.sum-se, ymax = nod.sum+se), width=0, size=1, position=position_dodge(.9)) +
  geom_hline(aes(yintercept = mean(nod.sum)), linetype='dashed', size=1, color="red") +
  xlab('Rhizobium Competitor') + ylab('Focal strain
nodule number') +
  ylim(0,52)+
  scale_x_discrete(labels=c('264A','399B','067A','277A','395A','066B','394B')) +
  theme(panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.border = element_blank(),axis.line = element_line(colour = "black"),axis.ticks = element_line(colour = "black", size = 0.5),
        axis.ticks.length = unit(0.2, "cm"),text = element_text(size=32), axis.text.x = element_text(size=28, vjust=0.5, angle=90), 
        axis.text.y = element_text(size=28,hjust=0.95), plot.title = element_text(hjust = 0.5, size=28)) +
  annotate("text", x = 1, y = 8, label = "a", size = 7) +
  annotate("text", x = 2, y = 8, label = "b", size = 7) +
  annotate("text", x = 3, y = 8, label = "b", size = 7) +
  annotate("text", x = 4, y = 8, label = "b", size = 7) +
  annotate("text", x = 5, y = 8, label = "b", size = 7) +
  annotate("text", x = 6, y = 8, label = "b", size = 7) +
  annotate("text", x = 7, y = 8, label = "b", size = 7) +
  ggtitle("Focal strain 412B")



#** Figure S6: Heatmap: focal strain nodule number with printed values  ####

# add values that will be on heatmap diagonal
ss.means
#add 1-strain inoc nod number on diag
Sgg= gg %>% 
    mutate(C394B=ifelse(focal.strain=='C394B', 49.6250, C394B),
         C066B=ifelse(focal.strain=='C066B', 38.7500, C066B),
         C395A=ifelse(focal.strain=='C395A', 39.7500, C395A),
         C277A=ifelse(focal.strain=='C277A', 41.4375, C277A),
         C067A=ifelse(focal.strain=='C067A', 48.9375, C067A),
         C399B=ifelse(focal.strain=='C399B', 38.7500, C399B),
         C264A=ifelse(focal.strain=='C264A', 9.1250, C264A),
         C412B=ifelse(focal.strain=='C412B', 12.4375, C412B))
# Now organize data
Snod.sum.long = pivot_longer(data=Sgg,
                            cols = -("focal.strain"),
                            names_to = "Competitor", 
                            values_to = "Nodule.number")

# Change order of competitor strains: least beneficial to most beneficial
Snod.sum.long$Competitor=factor(Snod.sum.long$Competitor, levels=c("C412B", "C264A", "C399B", "C067A", "C277A", "C395A", "C066B", "C394B"))

# Change order of focal strains: least beneficial to most beneficial
Snod.sum.long$focal.strain=factor(Snod.sum.long$focal.strain, levels=c("C412B", "C264A", "C399B", "C067A", "C277A", "C395A", "C066B", "C394B"))


ggplot(data = Snod.sum.long, mapping = aes(x = focal.strain, y = Competitor, fill = Nodule.number)) +
  geom_tile(color = "gray") +
  geom_text(aes(label=round(Nodule.number,0)),color= "white")+
  scale_fill_continuous(high = "#132B43", low = "#56B1F7", na.value="white") +
  scale_x_discrete(labels=c('412B','264A','399B','067A','277A','395A','066B','394B')) +
  scale_y_discrete(labels=c('412B','264A','399B','067A','277A','395A','066B','394B')) +
  xlab(label = "Focal strain") +
  ylab(label = "Competitor strain")+
theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
        axis.ticks = element_blank(), legend.position="right", legend.justification="top", legend.title=element_text(), 
        text = element_text(size=30), axis.text.x = element_text(size=28, vjust=0.5, angle=90), axis.text.y = element_text(size=28,hjust=0.95))+
  labs(fill='Focal 
strain
nodule 
number')





#***FIG S7 a-h: nodule proportion ####

# 394B nod prop plot
Rp394B$Competitor=factor(Rp394B$Competitor, levels=c("C412B", "C264A", "C399B", "C067A", "C277A", "C395A", "C066B"))

ggplot(Rp394B, aes(x= Competitor, y= nod.prop)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=nod.prop-se, ymax = nod.prop+se), width=0, size=1, position=position_dodge(.9)) +
  geom_hline(aes(yintercept = mean(nod.prop)), linetype='dashed', size=1, color="red") +
  xlab('Rhizobium Competitor') + ylab('Focal strain 
nodule proportion') +
  ylim(0,1.05)+
  scale_x_discrete(labels=c('412B','264A','399B','067A','277A','395A','066B')) +
  theme(panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.border = element_blank(),axis.line = element_line(colour = "black"),axis.ticks = element_line(colour = "black", size = 0.5),
        axis.ticks.length = unit(0.2, "cm"),text = element_text(size=32), axis.text.x = element_text(size=28, vjust=0.5, angle=90), 
        axis.text.y = element_text(size=28,hjust=0.95), plot.title = element_text(hjust = 0.5,size=28)) +
  annotate("text", x = 1, y = 1.05, label = "a", size = 7) +
  annotate("text", x = 2, y = 1.05, label = "a", size = 7) +
  annotate("text", x = 3, y = 0.87, label = "c", size = 7) +
  annotate("text", x = 4, y = 0.95, label = "b", size = 7) +
  annotate("text", x = 5, y = 0.87, label = "c", size = 7) +
  annotate("text", x = 6, y = 0.87, label = "c", size = 7) +
  annotate("text", x = 7, y = 0.87, label = "c", size = 7) +
  ggtitle("Focal strain 394B")

# 066B nod prop plot
Rp066B$Competitor=factor(Rp066B$Competitor, levels=c("C412B", "C264A", "C399B","C067A", "C277A", "C395A", "C394B"))

ggplot(Rp066B, aes(x= Competitor, y= nod.prop)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=nod.prop-se, ymax = nod.prop+se), width=0, size=1, position=position_dodge(.9)) +
  geom_hline(aes(yintercept = mean(nod.prop)), linetype='dashed', size=1, color="red") +
  xlab('Rhizobium Competitor') + ylab('Focal strain 
nodule proportion') +
  ylim(0,1.05) +
  scale_x_discrete(labels=c('412B','264A','399B','067A','277A','395A','394B')) +
  theme(panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.border = element_blank(),axis.line = element_line(colour = "black"),axis.ticks = element_line(colour = "black", size = 0.5),
        axis.ticks.length = unit(0.2, "cm"),text = element_text(size=32), axis.text.x = element_text(size=28, vjust=0.5, angle=90), 
        axis.text.y = element_text(size=28,hjust=0.95), plot.title = element_text(hjust = 0.5, size=28)) +
  annotate("text", x = 1, y = 1.05, label = "a", size = 7) +
  annotate("text", x = 2, y = 1.05, label = "a", size = 7) +
  annotate("text", x = 3, y = 0.7, label = "cd", size = 7) +
  annotate("text", x = 4, y = 0.85, label = "b", size = 7) +
  annotate("text", x = 5, y = 0.58, label = "d", size = 7) +
  annotate("text", x = 6, y = 0.7, label = "c", size = 7) +
  annotate("text", x = 7, y = 0.39, label = "e", size = 7) +
  ggtitle("Focal strain 066B")


# 395A nod prop plot
Rp395A$Competitor=factor(Rp395A$Competitor, levels=c("C412B", "C264A", "C399B", "C067A", "C277A", "C066B", "C394B"))

ggplot(Rp395A, aes(x= Competitor, y= nod.prop)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=nod.prop-se, ymax = nod.prop+se), width=0, size=1, position=position_dodge(.9)) +
  geom_hline(aes(yintercept = mean(nod.prop)), linetype='dashed', size=1, color="red") +
  xlab('Rhizobium Competitor') + ylab('Focal strain 
nodule proportion') +
  ylim(0, 1.05)+
  scale_x_discrete(labels=c('412B','264A', '399B','067A','277A','066B','394B')) +
  theme(panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.border = element_blank(),axis.line = element_line(colour = "black"),axis.ticks = element_line(colour = "black", size = 0.5),
        axis.ticks.length = unit(0.2, "cm"),text = element_text(size=32), axis.text.x = element_text(size=28, vjust=0.5, angle=90), 
        axis.text.y = element_text(size=28,hjust=0.95), plot.title = element_text(hjust = 0.5, size=28)) +
  annotate("text", x = 1, y = 1.05, label = "a", size = 7) +
  annotate("text", x = 2, y = 1.05, label = "a", size = 7) +
  annotate("text", x = 3, y = 0.72, label = "c", size = 7) +
  annotate("text", x = 4, y = 0.85, label = "b", size = 7) +
  annotate("text", x = 5, y = 0.58, label = "d", size = 7) +
  annotate("text", x = 6, y = 0.58, label = "d", size = 7) +
  annotate("text", x = 7, y = 0.43, label = "e", size = 7) +
  ggtitle("Focal strain 395A")


# 277A nod prop plot
Rp277A$Competitor=factor(Rp277A$Competitor, levels=c("C412B", "C264A", "C399B", "C067A", "C395A", "C066B", "C394B"))

ggplot(Rp277A, aes(x= Competitor, y= nod.prop)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=nod.prop-se, ymax = nod.prop+se), width=0, size=1, position=position_dodge(.9)) +
  geom_hline(aes(yintercept = mean(nod.prop)), linetype='dashed', size=1, color="red") +
  xlab('Rhizobium Competitor') + ylab('Focal strain 
nodule proportion') +
  ylim(0,1.05)+
  scale_x_discrete(labels=c('412B','264A', '399B','067A','395A','066B','394B')) +
  theme(panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.border = element_blank(),axis.line = element_line(colour = "black"),axis.ticks = element_line(colour = "black", size = 0.5),
        axis.ticks.length = unit(0.2, "cm"),text = element_text(size=32), axis.text.x = element_text(size=28, vjust=0.5, angle=90), 
        axis.text.y = element_text(size=28,hjust=0.95), plot.title = element_text(hjust = 0.5, size=28)) +
  annotate("text", x = 1, y = 1.05, label = "a", size = 7) +
  annotate("text", x = 2, y = 1.05, label = "a", size = 7) +
  annotate("text", x = 3, y = 0.74, label = "c", size = 7) +
  annotate("text", x = 4, y = 0.9, label = "b", size = 7) +
  annotate("text", x = 5, y = 0.74, label = "c", size = 7) +
  annotate("text", x = 6, y = 0.74, label = "c", size = 7) +
  annotate("text", x = 7, y = 0.43, label = "d", size = 7) +
  ggtitle("Focal strain 277A")

# 067A nod prop plot
Rp067A$Competitor=factor(Rp067A$Competitor, levels=c("C412B", "C264A", "C399B", "C277A", "C395A", "C066B", "C394B"))

ggplot(Rp067A, aes(x= Competitor, y= nod.prop)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=nod.prop-se, ymax = nod.prop+se), width=0, size=1, position=position_dodge(.9)) +
  geom_hline(aes(yintercept = mean(nod.prop)), linetype='dashed', size=1, color="red") +
  xlab('Rhizobium Competitor') + ylab('Focal strain 
nodule proportion') +
  scale_x_discrete(labels=c('412B','264A','399B','277A','395A','066B','394B')) +
  ylim(0,1.05)+
  theme(panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.border = element_blank(),axis.line = element_line(colour = "black"),axis.ticks = element_line(colour = "black", size = 0.5),
        axis.ticks.length = unit(0.2, "cm"),text = element_text(size=32), axis.text.x = element_text(size=28, vjust=0.5, angle=90), 
        axis.text.y = element_text(size=28,hjust=0.95), plot.title = element_text(hjust = 0.5, size=28)) +
  annotate("text", x = 1, y = 1.05, label = "a", size = 7) +
  annotate("text", x = 2, y = 1.05, label = "a", size = 7) +
  annotate("text", x = 3, y = 0.41, label = "b", size = 7) +
  annotate("text", x = 4, y = 0.35, label = "bc", size = 7) +
  annotate("text", x = 5, y = 0.41, label = "b", size = 7) +
  annotate("text", x = 6, y = 0.41, label = "b", size = 7) +
  annotate("text", x = 7, y = 0.28, label = "c", size = 7) +
  ggtitle("Focal strain 067A")

# 399B nod prop plot
Rp399B$Competitor=factor(Rp399B$Competitor, levels=c("C412B", "C264A", "C067A", "C277A", "C395A", "C066B", "C394B"))

ggplot(Rp399B, aes(x= Competitor, y= nod.prop)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=nod.prop-se, ymax = nod.prop+se), width=0, size=1, position=position_dodge(.9)) +
  geom_hline(aes(yintercept = mean(nod.prop)), linetype='dashed', size=1, color="red") +
  xlab('Rhizobium Competitor') + ylab('Focal strain 
nodule proportion') +
  ylim(0,1.05)+
  scale_x_discrete(labels=c('412B','264A','067A','277A','395A','066B','394B')) +
  theme(panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.border = element_blank(),axis.line = element_line(colour = "black"),axis.ticks = element_line(colour = "black", size = 0.5),
        axis.ticks.length = unit(0.2, "cm"),text = element_text(size=32), axis.text.x = element_text(size=28, vjust=0.5, angle=90), 
        axis.text.y = element_text(size=28,hjust=0.95), plot.title = element_text(hjust = 0.5, size=28)) +
  annotate("text", x = 1, y = 1.05, label = "a", size = 7) +
  annotate("text", x = 2, y = 1.05, label = "a", size = 7) +
  annotate("text", x = 3, y = 0.83, label = "b", size = 7) +
  annotate("text", x = 4, y = 0.58, label = "c", size = 7) +
  annotate("text", x = 5, y = 0.58, label = "c", size = 7) +
  annotate("text", x = 6, y = 0.58, label = "c", size = 7) +
  annotate("text", x = 7, y = 0.38, label = "d", size = 7) +
  ggtitle("Focal strain 399B")


# 264A nod prop plot
Rp264A$Competitor=factor(Rp264A$Competitor, levels=c("C412B", "C399B","C067A", "C277A", "C395A","C066B", "C394B"))

#plot with stats
ggplot(Rp264A, aes(x= Competitor, y= nod.prop)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=nod.prop-se, ymax = nod.prop+se), width=0, size=1, position=position_dodge(.9)) +
  geom_hline(aes(yintercept = mean(nod.prop)), linetype='dashed', size=1, color="red") +
  xlab('Rhizobium Competitor') + ylab('Focal strain 
nodule proportion') +
  ylim(0, 1.05)+
  scale_x_discrete(labels=c('412B','399B','067A','277A','395A','066B','394B')) +
  theme(panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.border = element_blank(),axis.line = element_line(colour = "black"),axis.ticks = element_line(colour = "black", size = 0.5),
        axis.ticks.length = unit(0.2, "cm"),text = element_text(size=32), axis.text.x = element_text(size=28, vjust=0.5, angle=90), 
        axis.text.y = element_text(size=28,hjust=0.95), plot.title = element_text(hjust = 0.5, size=28)) +
  annotate("text", x = 1, y = 0.55, label = "a", size = 7) +
  annotate("text", x = 2, y = 0.17, label = "b", size = 7) +
  annotate("text", x = 3, y = 0.17, label = "b", size = 7) +
  annotate("text", x = 4, y = 0.17, label = "b", size = 7) +
  annotate("text", x = 5, y = 0.17, label = "b", size = 7) +
  annotate("text", x = 6, y = 0.17, label = "b", size = 7) +
  annotate("text", x = 7, y = 0.17, label = "b", size = 7) +
  ggtitle("Focal strain 264A")

# 412 nod prop plot
Rp412B$Competitor=factor(Rp412B$Competitor, levels=c("C264A", "C399B","C067A", "C277A", "C395A","C066B", "C394B"))

#plot with stats
ggplot(Rp412B, aes(x= Competitor, y= nod.prop)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=nod.prop-se, ymax = nod.prop+se), width=0, size=1, position=position_dodge(.9)) +
  geom_hline(aes(yintercept = mean(nod.prop)), linetype='dashed', size=1, color="red") +
  xlab('Rhizobium Competitor') + ylab('Focal strain 
nodule proportion') +
  ylim(0,1.05)+
  scale_x_discrete(labels=c('264A','399B','067A','277A','395A','066B','394B')) +
  theme(panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.border = element_blank(),axis.line = element_line(colour = "black"),axis.ticks = element_line(colour = "black", size = 0.5),
        axis.ticks.length = unit(0.2, "cm"),text = element_text(size=32), axis.text.x = element_text(size=28, vjust=0.5, angle=90), 
        axis.text.y = element_text(size=28,hjust=0.95), plot.title = element_text(hjust = 0.5, size=28)) +
  annotate("text", x = 1, y = 0.95, label = "a", size = 7) +
  annotate("text", x = 2, y = 0.18, label = "b", size = 7) +
  annotate("text", x = 3, y = 0.18, label = "b", size = 7) +
  annotate("text", x = 4, y = 0.18, label = "b", size = 7) +
  annotate("text", x = 5, y = 0.18, label = "b", size = 7) +
  annotate("text", x = 6, y = 0.18, label = "b", size = 7) +
  annotate("text", x = 7, y = 0.18, label = "b", size = 7) +
  ggtitle("Focal strain 412B")






#** FIG S8: Shoot mass per nodule benefit of partner choice model  ####

# use SD1m dataset 
# calculate shoot mass per nodule in 2-strain inocula
SD1m$plant.benefit = SD1m$shoot.mass - 0.020125
SD1m$shootpernodtot= SD1m$plant.benefit/SD1m$nods.tot

# add shoot per nod in 1-strain innocula to SD1m
spn = ss.means[,c(1,10)]
names(spn)[2] = 'shootpernod.ss.o'
SD1m = merge(SD1m, spn, all.x=T)

spn.comp = spn
names(spn.comp) = c('strain.comp.o', 'shootpernod.ss.comp')
SD1m = merge(SD1m, spn.comp, all.x=T)

# add a column to ID more beneficial strain in a pot (used for checks)
SD1m$strain.best = apply(SD1m, 1, function(x) ifelse(x['shootpernod.ss.o'] > x['shootpernod.ss.comp'], x['strain.o'], x['strain.comp.o'])) 
SD1m$strain.worst = apply(SD1m, 1, function(x) ifelse(x['shootpernod.ss.o'] < x['shootpernod.ss.comp'], x['strain.o'], x['strain.comp.o']))

# Calculate the Deviation from a neutral expectation (dne) metric: (W_s1s2 - mean(W_s1 + W_s2))/mean(W_s1 + W_s2) 
# (ie the neutral expectation for host fitness in the absence of partner choice)

# calculate the neutral expectation of shoot mass under no partner choice (ie the mean shoot mass of the 2 strains in 1-strain inocula)
SD1m$MeanShoot2Strains = (SD1m$shootpernod.ss.comp + SD1m$shootpernod.ss.o)/2 #used for check
# calculate the difference between actual shoot mass in 2-strain inocula and expected shoot mass
SD1m$Dif = SD1m$shootpernodtot - ((SD1m$shootpernod.ss.comp + SD1m$shootpernod.ss.o)/2) #used for check
# calculate the shoot mass deviation from the neutral expectation (use in model)
SD1m$dne = (SD1m$shootpernodtot - ((SD1m$shootpernod.ss.comp + SD1m$shootpernod.ss.o)/2))/((SD1m$shootpernod.ss.comp + SD1m$shootpernod.ss.o)/2)

# add new column "NodPropBestStrain" to assign nod prop of more beneficial strain
SD1m$NodPropBestStrain = apply(SD1m, 1, function(x) ifelse(x['shootpernod.ss.o'] > x['shootpernod.ss.comp'], x['nod.prop'], x['comp.nod.prop'])) 

# Check data
str(SD1m) #SD1m has 896 obs. of  39 variables
# change to numeric
SD1m$NodPropBestStrain = as.numeric(SD1m$NodPropBestStrain)
# droplevels (no neg)
SD1m=droplevels(SD1m)
levels(SD1m$strain.o) #good
levels(SD1m$strain.comp.o) #good

# Use only one observation (row) per plant pot from SD1m
SD1m2 = SD1m[!duplicated(SD1m[, "Unique.ID"]),] 
str(SD1m2) # 439 obs. of 39 variables 

#Drop rows where nodpropbeststrain=na, so models are fitted to same data set despite na (ex. na for nodpropbeststrain but not for marker-results in diff data sets))
SD1m2=SD1m2[complete.cases(SD1m2[,39]),] #39 is column for nodpropbeststrain 


# Benefit of partner choice to the host model
M5.2 = lmer(dne ~ NodPropBestStrain + marker + (1|strain.o) + (1|strain.comp.o) + (1|Block), data=SD1m2, na.action=na.exclude)

summary(M5.2)
#Random effects:
# Groups        Name        Variance Std.Dev.
# Block         (Intercept)    0.0    0.00   
# strain.comp.o (Intercept)  101.3   10.07   
# strain.o      (Intercept)  113.0   10.63   
# Residual                  1557.1   39.46   
#Fixed effects:
#                  Estimate Std. Error t value
#(Intercept)        -46.743      9.970  -4.688
#NodPropBestStrain   62.586      9.383   6.670
#markersfGFP         -6.921      3.772  -1.835

# check model assumptions
plot(M5.2) #heterogeneity of variance 
qqnorm(resid(M5.2))
qqline(resid(M5.2))

# test fixed effects
drop1(M5.2, test = "Chisq") 
#                  npar    AIC     LRT   Pr(Chi)    
#<none>                 4504.3                      
#NodPropBestStrain    1 4529.8 27.5913 1.498e-07 ***
#marker               1 4505.7  3.4046   0.06501 .



# Benefit to host figure prep
SD1m2$combo = paste(SD1m2$strain.o, SD1m2$strain.comp.o,sep=':')

# aggregate the means nodpropbeststrain and "dne" for the plot
SD1m3 = aggregate(cbind(NodPropBestStrain, dne) ~ combo + strain.o + strain.comp.o, FUN=mean, data=SD1m2)
# add standard error
SD1m3 = cbind(SD1m3, aggregate(cbind(NodPropBestStrain, dne) ~ combo + strain.o + strain.comp.o, FUN=function(x) sd(x)/length(x), data=SD1m2)[4:5])
head(SD1m3)
# Rename se columns
names(SD1m3)[6:7] = paste(names(SD1m3)[6:7],'se',sep='.')

## *****a. plot benefit of partner choice with influential point ####
ggplot(SD1m3, aes(x=NodPropBestStrain, y=dne, alpha=.6)) + 
  geom_point(size=4) + 
  geom_errorbar(aes(ymin=dne-dne.se, ymax=dne+dne.se)) +
  geom_errorbarh(aes(xmin=NodPropBestStrain-NodPropBestStrain.se, xmax=NodPropBestStrain+NodPropBestStrain.se)) +
  stat_smooth(method='lm', alpha=.2, size=1, color='black') +
  geom_hline(yintercept = 0, linetype='longdash', size=1) + #plant fitness under neutral nod occupancy
  geom_vline(xintercept = 0.5, linetype='dashed', color='black', size=1) + #partner choice:better strain occupies >50% of nodules
  theme(legend.position='none') + xlab('Nodule proportion of more beneficial strain')+ ylab('Shoot mass (g) per nodule 
deviation from neutral expectation') +
  theme(panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),axis.line = element_line(colour = "black"),axis.ticks = element_line(colour = "black", size = 0.5), axis.ticks.length = unit(0.2, "cm"),text = element_text(size=22), axis.text.x = element_text(size=22, vjust=0.95), axis.text.y = element_text(size=22,hjust=0.95)) 



### drop bad strain combo and rerun model
SD1m2temp=SD1m2[SD1m2$strain.label !="C264A_C412B",]

M5.3=lmer(dne ~  NodPropBestStrain + marker + (1|strain.o) + (1|strain.comp.o) + (1|Block), data=SD1m2temp, na.action=na.omit)


summary(M5.3) #fit is still singular
#Random effects:
# Groups        Name        Variance Std.Dev.
# Block         (Intercept) 0.3589   0.5991  
# strain.comp.o (Intercept) 0.1632   0.4039  
# strain.o      (Intercept) 0.0000   0.0000  
# Residual                  6.4304   2.5358  
#Fixed effects:
#                   Estimate Std. Error t value
#(Intercept)       -0.630412   0.557296  -1.131
#NodPropBestStrain  1.692057   0.564099   3.000
#markersfGFP       -0.004277   0.244386  -0.018

# check model assumptions
plot(M5.3) # residuals are better than previous model
qqnorm(resid(M5.3))
qqline(resid(M5.3))

# test fixed effects
drop1(M5.3, test = "Chisq")
#                   npar    AIC    LRT  Pr(Chi)   
#<none>                 2052.9                   
#NodPropBestStrain    1 2057.8 6.8955 0.008641 **
#marker               1 2050.9 0.0003 0.985630


# Benefit to host figure prep
SD1m2temp$combo = paste(SD1m2temp$strain.o, SD1m2temp$strain.comp.o,sep=':')

# aggregate the means nodpropbeststrain and "dne" for the plot
SD1m3temp = aggregate(cbind(NodPropBestStrain, dne) ~ combo + strain.o + strain.comp.o, FUN=mean, data=SD1m2temp)
# add standard error
SD1m3temp = cbind(SD1m3temp, aggregate(cbind(NodPropBestStrain, dne) ~ combo + strain.o + strain.comp.o, FUN=function(x) sd(x)/length(x), data=SD1m2temp)[4:5])
head(SD1m3temp)
# Rename se columns
names(SD1m3temp)[6:7] = paste(names(SD1m3temp)[6:7],'se',sep='.')


## *****b. plot benefit of partner choice without influential point ####
ggplot(SD1m3temp, aes(x=NodPropBestStrain, y=dne, alpha=.6)) + 
  geom_point(size=4) + 
  geom_errorbar(aes(ymin=dne-dne.se, ymax=dne+dne.se)) +
  geom_errorbarh(aes(xmin=NodPropBestStrain-NodPropBestStrain.se, xmax=NodPropBestStrain+NodPropBestStrain.se)) +
  stat_smooth(method='lm', alpha=.2, size=1, color='black') +
  geom_hline(yintercept = 0, linetype='longdash', size=1) + #plant fitness under neutral nod occupancy
  geom_vline(xintercept = 0.5, linetype='dashed', color='black', size=1) + #partner choice:better strain occupies >50% of nodules
  theme(legend.position='none') + xlab('Nodule proportion of more beneficial strain')+ ylab('Shoot mass (g) per nodule 
deviation from neutral expectation') +
  theme(panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),axis.line = element_line(colour = "black"),axis.ticks = element_line(colour = "black", size = 0.5), axis.ticks.length = unit(0.2, "cm"),text = element_text(size=22), axis.text.x = element_text(size=22, vjust=0.95), axis.text.y = element_text(size=22,hjust=0.95)) 




#** Fig S9: cfu heatmaps ####

#***** a. cfu: among nodule sanctions ####


#### Heatmap of cfu per nodule in "full" nodules

#  Create subset for heatmap 

# Subset competitor strain C394B
C394B.comp.full=rbind(
  (cfu.full[which(cfu.full$strain.label=="C394B_C264A"),]), 
  (cfu.full[which(cfu.full$strain.label=="C277A_C394B"),]), 
  (cfu.full[which(cfu.full$strain.label=="C394B_C395A"),]), 
  (cfu.full[which(cfu.full$strain.label=="C394B_C412B"),]), 
  (cfu.full[which(cfu.full$strain.label=="C067A_C394B"),]), 
  (cfu.full[which(cfu.full$strain.label=="C066B_C394B"),]),
  (cfu.full[which(cfu.full$strain.label=="C399B_C394B"),]))

# remove C394B from subset
C394B.comp.full=C394B.comp.full[-which(C394B.comp.full$strain.o=="C394B"),]

# Drop levels after subsetting
C394B.comp.full = droplevels(C394B.comp.full)

# Subset competitor strain C066B
C066B.comp.full=rbind(
  (cfu.full[which(cfu.full$strain.label=="C066B_C067A"),]), 
  (cfu.full[which(cfu.full$strain.label=="C066B_C264A"),]), 
  (cfu.full[which(cfu.full$strain.label=="C066B_C277A"),]), 
  (cfu.full[which(cfu.full$strain.label=="C066B_C412B"),]), 
  (cfu.full[which(cfu.full$strain.label=="C066B_C395A"),]), 
  (cfu.full[which(cfu.full$strain.label=="C066B_C394B"),]),
  (cfu.full[which(cfu.full$strain.label=="C066B_C399B"),]))

# remove C066B from subset
C066B.comp.full=C066B.comp.full[-which(C066B.comp.full$strain.o=="C066B"),]

# Drop levels after subsetting
C066B.comp.full = droplevels(C066B.comp.full)

# Subset competitor strain C395A 
C395A.comp.full=rbind(
  (cfu.full[which(cfu.full$strain.label=="C395A_C264A"),]), 
  (cfu.full[which(cfu.full$strain.label=="C277A_C395A"),]), 
  (cfu.full[which(cfu.full$strain.label=="C394B_C395A"),]), 
  (cfu.full[which(cfu.full$strain.label=="C395A_C412B"),]), 
  (cfu.full[which(cfu.full$strain.label=="C067A_C395A"),]), 
  (cfu.full[which(cfu.full$strain.label=="C066B_C395A"),]),
  (cfu.full[which(cfu.full$strain.label=="C399B_C395A"),]))

# remove C395A from subset
C395A.comp.full=C395A.comp.full[-which(C395A.comp.full$strain.o=="C395A"),]

# Drop levels after subsetting
C395A.comp.full = droplevels(C395A.comp.full)


# Subset competitor strain C277A 
C277A.comp.full=rbind(
  (cfu.full[which(cfu.full$strain.label=="C277A_C264A"),]), 
  (cfu.full[which(cfu.full$strain.label=="C277A_C395A"),]), 
  (cfu.full[which(cfu.full$strain.label=="C277A_C394B"),]), 
  (cfu.full[which(cfu.full$strain.label=="C277A_C412B"),]), 
  (cfu.full[which(cfu.full$strain.label=="C067A_C277A"),]), 
  (cfu.full[which(cfu.full$strain.label=="C066B_C277A"),]),
  (cfu.full[which(cfu.full$strain.label=="C399B_C277A"),]))

# remove C277A from subset
C277A.comp.full=C277A.comp.full[-which(C277A.comp.full$strain.o=="C277A"),]

# Drop levels after subsetting
C277A.comp.full = droplevels(C277A.comp.full)

# Subset competitor strain C067A 
C067A.comp.full=rbind(
  (cfu.full[which(cfu.full$strain.label=="C067A_C264A"),]), 
  (cfu.full[which(cfu.full$strain.label=="C067A_C395A"),]), 
  (cfu.full[which(cfu.full$strain.label=="C067A_C394B"),]), 
  (cfu.full[which(cfu.full$strain.label=="C067A_C412B"),]), 
  (cfu.full[which(cfu.full$strain.label=="C067A_C399B"),]), 
  (cfu.full[which(cfu.full$strain.label=="C066B_C067A"),]),
  (cfu.full[which(cfu.full$strain.label=="C067A_C277A"),]))

# remove C067A from subset
C067A.comp.full=C067A.comp.full[-which(C067A.comp.full$strain.o=="C067A"),]

# Drop levels after subsetting
C067A.comp.full = droplevels(C067A.comp.full)

# Subset competitor strain C399B (C399B and S399B)
C399B.comp.full=rbind(
  (cfu.full[which(cfu.full$strain.label=="C399B_C264A"),]), 
  (cfu.full[which(cfu.full$strain.label=="C399B_C395A"),]), 
  (cfu.full[which(cfu.full$strain.label=="C399B_C394B"),]), 
  (cfu.full[which(cfu.full$strain.label=="C399B_C412B"),]), 
  (cfu.full[which(cfu.full$strain.label=="C067A_C399B"),]), 
  (cfu.full[which(cfu.full$strain.label=="C066B_C399B"),]),
  (cfu.full[which(cfu.full$strain.label=="C399B_C277A"),]))

# remove C399B from subset
C399B.comp.full=C399B.comp.full[-which(C399B.comp.full$strain.o=="C399B"),]

# Drop levels after subsetting
C399B.comp.full = droplevels(C399B.comp.full)

# Subset competitor strain C264A (C264A and S264A)
C264A.comp.full=rbind(
  (cfu.full[which(cfu.full$strain.label=="C264A_C412B"),]), 
  (cfu.full[which(cfu.full$strain.label=="C277A_C264A"),]), 
  (cfu.full[which(cfu.full$strain.label=="C394B_C264A"),]), 
  (cfu.full[which(cfu.full$strain.label=="C395A_C264A"),]), 
  (cfu.full[which(cfu.full$strain.label=="C067A_C264A"),]), 
  (cfu.full[which(cfu.full$strain.label=="C066B_C264A"),]),
  (cfu.full[which(cfu.full$strain.label=="C399B_C264A"),]))

# remove C264A from subset
C264A.comp.full=C264A.comp.full[-which(C264A.comp.full$strain.o=="C264A"),]

# Drop levels after subsetting
C264A.comp.full = droplevels(C264A.comp.full)

# Subset competitor strain C412B 
C412B.comp.full=rbind(
  (cfu.full[which(cfu.full$strain.label=="C264A_C412B"),]), 
  (cfu.full[which(cfu.full$strain.label=="C277A_C412B"),]), 
  (cfu.full[which(cfu.full$strain.label=="C394B_C412B"),]), 
  (cfu.full[which(cfu.full$strain.label=="C395A_C412B"),]), 
  (cfu.full[which(cfu.full$strain.label=="C067A_C412B"),]), 
  (cfu.full[which(cfu.full$strain.label=="C066B_C412B"),]),
  (cfu.full[which(cfu.full$strain.label=="C399B_C412B"),]))

# remove C412B from subset
C412B.comp.full=C412B.comp.full[-which(C412B.comp.full$strain.o=="C412B"),]

# Drop levels after subsetting
C412B.comp.full = droplevels(C412B.comp.full)


# Mean cfu per nodule for each strain against the 8 competitors 
# change column names for heatmap code. Each dataframe will show each focal strains mean cfu per nodule against a competitor

# Competitor C394B
fc1=aggregate(CFU ~ strain.o, C394B.comp.full, mean)
fc1
colnames(fc1)=c("focal.strain","C394B")

# Competitor C066B
fc2=aggregate(CFU ~ strain.o, C066B.comp.full, mean)
fc2
colnames(fc2)=c("focal.strain","C066B")

# Competitor C395A 
fc3=aggregate(CFU ~ strain.o, C395A.comp.full, mean)
fc3
colnames(fc3)=c("focal.strain","C395A")

# Competitor C277A 
fc4=aggregate(CFU ~ strain.o, C277A.comp.full, mean)
fc4
colnames(fc4)=c("focal.strain","C277A")

# Competitor C067A 
fc5=aggregate(CFU ~ strain.o, C067A.comp.full, mean)
fc5
colnames(fc5)=c("focal.strain","C067A") 

# Competitor C399B
fc6=aggregate(CFU ~ strain.o, C399B.comp.full, mean)
fc6
colnames(fc6)=c("focal.strain","C399B")

# Competitor C264A
fc7=aggregate( CFU ~ strain.o, C264A.comp.full, mean)
fc7
colnames(fc7)=c("focal.strain","C264A")

# Competitor C412B
fc8=aggregate(CFU ~ strain.o, C412B.comp.full, mean)
fc8
colnames(fc8)=c("focal.strain","C412B")

# Merge all competitor data frames
aa4=merge(fc1, fc2, by = "focal.strain", all.x = T, all.y = T)
aa4

bb4= merge(aa4, fc3, by = "focal.strain", all.x = T, all.y = T)
bb4

cc4= merge(bb4, fc4, by = "focal.strain", all.x = T, all.y = T)
cc4

dd4= merge(cc4, fc5, by = "focal.strain", all.x = T, all.y = T)
dd4

ee4= merge(dd4, fc6, by = "focal.strain", all.x = T, all.y = T)
ee4

ff4= merge(ee4, fc7, by = "focal.strain", all.x = T, all.y = T)
ff4

gg4= merge(ff4, fc8, by = "focal.strain", all.x = T, all.y = T)


# add 1-strain inocula mean cfu values along the heatmap diagonal
# mean cfu in 1-strain inocula and standard error
cfu.means = aggregate(CFU ~ strain.o, FUN=mean, data=subset(cfu,inoc.type=='single'))
names(cfu.means)[2] = 'ss.mean.cfu'
cfu.means

#  strain.o ss.mean.cfu
#1    C066B   1542559.5
#2    C067A    890145.8
#3    C264A      4275.0
#4    C277A    755825.0
#5    C394B    752546.9
#6    C395A   1494791.7
#7    C399B    127267.9
#8    C412B   1219867.6

# add values that will be on heatmap diagonal
gg4= gg4 %>% 
    mutate(C394B=ifelse(focal.strain=='C394B', 752546.9, C394B),
         C066B=ifelse(focal.strain=='C066B', 1542559.5, C066B),
         C395A=ifelse(focal.strain=='C395A', 1494791.7, C395A),
         C277A=ifelse(focal.strain=='C277A', 755825.0, C277A),
         C067A=ifelse(focal.strain=='C067A', 890145.8, C067A),
         C399B=ifelse(focal.strain=='C399B', 127267.9, C399B),
         C264A=ifelse(focal.strain=='C264A', 4275.0, C264A),
         C412B=ifelse(focal.strain=='C412B', 1219867.6, C412B))

#Now organize data
CFU.long = pivot_longer(data=gg4,
                        cols = -("focal.strain"),
                        names_to = "Competitor", 
                        values_to = "CFU")

# Change order of competitor strains: least beneficial to most beneficial
CFU.long$Competitor=factor(CFU.long$Competitor, levels=c("C412B", "C264A", "C399B", "C067A", "C277A", "C395A", "C066B", "C394B"))

# Change order of focal strains: least beneficial to most beneficial
CFU.long$focal.strain=factor(CFU.long$focal.strain, levels=c("C412B", "C264A", "C399B", "C067A", "C277A", "C395A", "C066B", "C394B"))

#heatmap
cfu.full.heatmap = ggplot(data = CFU.long, mapping = aes(x = focal.strain, y = Competitor, fill = CFU)) +
  geom_tile(color = "gray") +
  scale_fill_continuous(high = "#132B43", low = "#56B1F7", na.value="white") +
  scale_x_discrete(labels=c('412B','264A','399B','067A','277A','395A','066B','394B')) +
  scale_y_discrete(labels=c('412B','264A','399B','067A','277A','395A','066B','394B')) +
  xlab(label = "Focal strain") +
  ylab(label = "Competitor strain")+
  theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
        axis.ticks = element_blank(), legend.position="right", legend.justification="top", legend.title=element_text(), 
        text = element_text(size=30), axis.text.x = element_text(size=28, vjust=0.5, angle=90), axis.text.y = element_text(size=28,hjust=0.95), axis.title.y.right = element_text(vjust=1.5))+
  labs(fill='CFU per 
nodule')

cfu.full.heatmap

#note: gray squares were added in powerpoint

#***** b. cfu: Intra-nodule sanctions ####
#### Heatmap of cfu proportion in mixed nodules 

#  Create subset for heatmap 

# Subset competitor strain C394B
C394B.comp.mixed=rbind(
  (cfu.mixed[which(cfu.mixed$strain.label=="C394B_C264A"),]), 
  (cfu.mixed[which(cfu.mixed$strain.label=="C277A_C394B"),]), 
  (cfu.mixed[which(cfu.mixed$strain.label=="C394B_C395A"),]), 
  (cfu.mixed[which(cfu.mixed$strain.label=="C394B_C412B"),]), 
  (cfu.mixed[which(cfu.mixed$strain.label=="C067A_C394B"),]), 
  (cfu.mixed[which(cfu.mixed$strain.label=="C066B_C394B"),]),
  (cfu.mixed[which(cfu.mixed$strain.label=="C399B_C394B"),]))

# remove C394B from subset
C394B.comp.mixed=C394B.comp.mixed[-which(C394B.comp.mixed$strain.o=="C394B"),]

# Drop levels after subsetting
C394B.comp.mixed = droplevels(C394B.comp.mixed)

# Subset competitor strain C066B
C066B.comp.mixed=rbind(
  (cfu.mixed[which(cfu.mixed$strain.label=="C066B_C067A"),]), 
  (cfu.mixed[which(cfu.mixed$strain.label=="C066B_C264A"),]), 
  (cfu.mixed[which(cfu.mixed$strain.label=="C066B_C277A"),]), 
  (cfu.mixed[which(cfu.mixed$strain.label=="C066B_C412B"),]), 
  (cfu.mixed[which(cfu.mixed$strain.label=="C066B_C395A"),]), 
  (cfu.mixed[which(cfu.mixed$strain.label=="C066B_C394B"),]),
  (cfu.mixed[which(cfu.mixed$strain.label=="C066B_C399B"),]))

# remove C066B from subset
C066B.comp.mixed=C066B.comp.mixed[-which(C066B.comp.mixed$strain.o=="C066B"),]

# Drop levels after subsetting
C066B.comp.mixed = droplevels(C066B.comp.mixed)

# Subset competitor strain C395A 
C395A.comp.mixed=rbind(
  (cfu.mixed[which(cfu.mixed$strain.label=="C395A_C264A"),]), 
  (cfu.mixed[which(cfu.mixed$strain.label=="C277A_C395A"),]), 
  (cfu.mixed[which(cfu.mixed$strain.label=="C394B_C395A"),]), 
  (cfu.mixed[which(cfu.mixed$strain.label=="C395A_C412B"),]), 
  (cfu.mixed[which(cfu.mixed$strain.label=="C067A_C395A"),]), 
  (cfu.mixed[which(cfu.mixed$strain.label=="C066B_C395A"),]),
  (cfu.mixed[which(cfu.mixed$strain.label=="C399B_C395A"),]))

# remove C395A from subset
C395A.comp.mixed=C395A.comp.mixed[-which(C395A.comp.mixed$strain.o=="C395A"),]

# Drop levels after subsetting
C395A.comp.mixed = droplevels(C395A.comp.mixed)


# Subset competitor strain C277A 
C277A.comp.mixed=rbind(
  (cfu.mixed[which(cfu.mixed$strain.label=="C277A_C264A"),]), 
  (cfu.mixed[which(cfu.mixed$strain.label=="C277A_C395A"),]), 
  (cfu.mixed[which(cfu.mixed$strain.label=="C277A_C394B"),]), 
  (cfu.mixed[which(cfu.mixed$strain.label=="C277A_C412B"),]), 
  (cfu.mixed[which(cfu.mixed$strain.label=="C067A_C277A"),]), 
  (cfu.mixed[which(cfu.mixed$strain.label=="C066B_C277A"),]),
  (cfu.mixed[which(cfu.mixed$strain.label=="C399B_C277A"),]))

# remove C277A from subset
C277A.comp.mixed=C277A.comp.mixed[-which(C277A.comp.mixed$strain.o=="C277A"),]

# Drop levels after subsetting
C277A.comp.mixed = droplevels(C277A.comp.mixed)

# Subset competitor strain C067A 
C067A.comp.mixed=rbind(
  (cfu.mixed[which(cfu.mixed$strain.label=="C067A_C264A"),]), 
  (cfu.mixed[which(cfu.mixed$strain.label=="C067A_C395A"),]), 
  (cfu.mixed[which(cfu.mixed$strain.label=="C067A_C394B"),]), 
  (cfu.mixed[which(cfu.mixed$strain.label=="C067A_C412B"),]), 
  (cfu.mixed[which(cfu.mixed$strain.label=="C067A_C399B"),]), 
  (cfu.mixed[which(cfu.mixed$strain.label=="C066B_C067A"),]),
  (cfu.mixed[which(cfu.mixed$strain.label=="C067A_C277A"),]))

# remove C067A from subset
C067A.comp.mixed=C067A.comp.mixed[-which(C067A.comp.mixed$strain.o=="C067A"),]

# Drop levels after subsetting
C067A.comp.mixed = droplevels(C067A.comp.mixed)

# Subset competitor strain C399B (C399B and S399B)
C399B.comp.mixed=rbind(
  (cfu.mixed[which(cfu.mixed$strain.label=="C399B_C264A"),]), 
  (cfu.mixed[which(cfu.mixed$strain.label=="C399B_C395A"),]), 
  (cfu.mixed[which(cfu.mixed$strain.label=="C399B_C394B"),]), 
  (cfu.mixed[which(cfu.mixed$strain.label=="C399B_C412B"),]), 
  (cfu.mixed[which(cfu.mixed$strain.label=="C067A_C399B"),]), 
  (cfu.mixed[which(cfu.mixed$strain.label=="C066B_C399B"),]),
  (cfu.mixed[which(cfu.mixed$strain.label=="C399B_C277A"),]))

# remove C399B from subset
C399B.comp.mixed=C399B.comp.mixed[-which(C399B.comp.mixed$strain.o=="C399B"),]

# Drop levels after subsetting
C399B.comp.mixed = droplevels(C399B.comp.mixed)

# Subset competitor strain C264A (C264A and S264A)
C264A.comp.mixed=rbind(
  (cfu.mixed[which(cfu.mixed$strain.label=="C264A_C412B"),]), 
  (cfu.mixed[which(cfu.mixed$strain.label=="C277A_C264A"),]), 
  (cfu.mixed[which(cfu.mixed$strain.label=="C394B_C264A"),]), 
  (cfu.mixed[which(cfu.mixed$strain.label=="C395A_C264A"),]), 
  (cfu.mixed[which(cfu.mixed$strain.label=="C067A_C264A"),]), 
  (cfu.mixed[which(cfu.mixed$strain.label=="C066B_C264A"),]),
  (cfu.mixed[which(cfu.mixed$strain.label=="C399B_C264A"),]))

# remove C264A from subset
C264A.comp.mixed=C264A.comp.mixed[-which(C264A.comp.mixed$strain.o=="C264A"),]

# Drop levels after subsetting
C264A.comp.mixed = droplevels(C264A.comp.mixed)

# Subset competitor strain C412B 
C412B.comp.mixed=rbind(
  (cfu.mixed[which(cfu.mixed$strain.label=="C264A_C412B"),]), 
  (cfu.mixed[which(cfu.mixed$strain.label=="C277A_C412B"),]), 
  (cfu.mixed[which(cfu.mixed$strain.label=="C394B_C412B"),]), 
  (cfu.mixed[which(cfu.mixed$strain.label=="C395A_C412B"),]), 
  (cfu.mixed[which(cfu.mixed$strain.label=="C067A_C412B"),]), 
  (cfu.mixed[which(cfu.mixed$strain.label=="C066B_C412B"),]),
  (cfu.mixed[which(cfu.mixed$strain.label=="C399B_C412B"),]))

# remove C412B from subset
C412B.comp.mixed=C412B.comp.mixed[-which(C412B.comp.mixed$strain.o=="C412B"),]

# Drop levels after subsetting
C412B.comp.mixed = droplevels(C412B.comp.mixed)


# Mean cfu proportion per mixed nodule for each strain against the 8 competitors 
# change column names for heatmap code. Each dataframe will show each focal strains mean cfu proportion per mixed nodule against a competitor

# Competitor C394B
mc1=aggregate(cfu.prop ~ strain.o, C394B.comp.mixed, mean)
mc1
colnames(mc1)=c("focal.strain","C394B")

# Competitor C066B
mc2=aggregate(cfu.prop ~ strain.o, C066B.comp.mixed, mean)
mc2
colnames(mc2)=c("focal.strain","C066B")

# Competitor C395A 
mc3=aggregate(cfu.prop ~ strain.o, C395A.comp.mixed, mean)
mc3
colnames(mc3)=c("focal.strain","C395A")

# Competitor C277A 
mc4=aggregate(cfu.prop ~ strain.o, C277A.comp.mixed, mean)
mc4
colnames(mc4)=c("focal.strain","C277A")

# Competitor C067A 
mc5=aggregate(cfu.prop ~ strain.o, C067A.comp.mixed, mean)
mc5
colnames(mc5)=c("focal.strain","C067A") 

# Competitor C399B
mc6=aggregate(cfu.prop ~ strain.o, C399B.comp.mixed, mean)
mc6
colnames(mc6)=c("focal.strain","C399B")

# Competitor C264A
mc7=aggregate( cfu.prop ~ strain.o, C264A.comp.mixed, mean)
mc7
colnames(mc7)=c("focal.strain","C264A")

# Competitor C412B
mc8=aggregate(cfu.prop ~ strain.o, C412B.comp.mixed, mean)
mc8
colnames(mc8)=c("focal.strain","C412B")

# Merge all competitor data frames
aa5=merge(mc1, mc2, by = "focal.strain", all.x = T, all.y = T)
aa5

bb5= merge(aa5, mc3, by = "focal.strain", all.x = T, all.y = T)
bb5

cc5= merge(bb5, mc4, by = "focal.strain", all.x = T, all.y = T)
cc5

dd5= merge(cc5, mc5, by = "focal.strain", all.x = T, all.y = T)
dd5

ee5= merge(dd5, mc6, by = "focal.strain", all.x = T, all.y = T)
ee5

ff5= merge(ee5, mc7, by = "focal.strain", all.x = T, all.y = T)
ff5

gg5= merge(ff5, mc8, by = "focal.strain", all.x = T, all.y = T)

# Now organize data
CFU.mixed.long = pivot_longer(data=gg5,
                        cols = -("focal.strain"),
                        names_to = "Competitor", 
                        values_to = "CFU.prop")


# Change order of competitor strains: least beneficial to most beneficial
CFU.mixed.long$Competitor=factor(CFU.mixed.long$Competitor, levels=c("C412B", "C264A", "C399B", "C067A", "C277A", "C395A", "C066B", "C394B"))

# Change order of focal strains: least beneficial to most beneficial
CFU.mixed.long$focal.strain=factor(CFU.mixed.long$focal.strain, levels=c("C412B", "C264A", "C399B", "C067A", "C277A", "C395A", "C066B", "C394B"))

# heatmap
cfu.mixed.heatmap = ggplot(data = CFU.mixed.long, mapping = aes(x = focal.strain, y = Competitor, fill = CFU.prop)) +
  geom_tile(color = "gray") +
  scale_fill_continuous(high = "#132B43", low = "#56B1F7", na.value="white") +
  scale_x_discrete(labels=c('412B','264A','399B','067A','277A','395A','066B','394B')) +
  scale_y_discrete(labels=c('412B','264A','399B','067A','277A','395A','066B','394B')) +
  xlab(label = "Focal strain") +
  ylab(label = "Competitor strain")+
  theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
        axis.ticks = element_blank(), legend.position="right", legend.justification="top", legend.title=element_text(), 
        text = element_text(size=30), axis.text.x = element_text(size=28, vjust=0.5, angle=90), axis.text.y = element_text(size=28,hjust=0.95), axis.title.y.right = element_text(vjust=1.5))+
  labs(fill='CFU 
proportion 
per nodule')

cfu.mixed.heatmap 

#note: gray squares were added in powerpoint

# print both cfu heatmaps
grid.arrange(cfu.full.heatmap, cfu.mixed.heatmap, nrow=1)


