####
# Christina Baer
# Begun: November 2018
# Last edited: December 20, 2019
# Trait variation and quantitative genetics for families of C. dilaticollis,
# dorsalis, and placida at different temperatures
# Analysis of final flat files for publication


# Set up workspace --------------------------------------------------------
# Working directory stem will need to be changed
setwd("C:/Users/cbaer/Dropbox/3_spp_quant_genetics_manuscript")
#library(readxl)
library(purrr)
library(dplyr)
library(stringr)
library(ggplot2)
library(scales)
library(lubridate)
library(survival)
library(coxme)
library(lme4)
library(nlme)
library(reshape2)
library(multcomp)


# Useful functions --------------------------------------------------------
# se: function for standard errors
se <- function(x) var(na.omit(x))/sqrt(length(na.omit(x)))

# theme_classic2 function -------------------------------------------------
# a modified theme_classic plot function to give grid-free but square plots
theme_classic2<-function (base_size = 11, base_family = "") 
{
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black",
                               size = 0.5),
      legend.key = element_blank(),
      strip.background = element_rect(
        fill = "white",
        colour = "black",
        size = 1
      ),
      aspect.ratio = 1,
      axis.text = element_text(size=10,family="sans"),
      complete = TRUE
    )
}

# MeanSESD ----------------------------------------------------------------
# MeanSESD: a function that will calculate the mean, mean+-standard error (SE),
# and mean+-standard deviation (SD) for use as a boxplot
MeanSESD <- function(x) {
  v <- c(mean(x) - sd(x), mean(x) - sd(x)/sqrt(length(x)), mean(x), mean(x) + sd(x)/sqrt(length(x)), mean(x) + sd(x))
  names(v) <- c("ymin", "lower", "middle", "upper", "ymax")
  v
}



# Import data -------------------------------------------------------------
# Cephaloleia dilaticollis data
dila<-read.csv("C_dilaticollis_flat_traits_2018-Nov-27.csv")

# Cephaloleia dorsalis data
dors<-read.csv("C_dorsalis_flat_traits_2018-Nov-27.csv")

# Cephaloleia placida data
plac<-read.csv("C_placida_flat_traits_2018-Nov-27.csv")

# check that all the column names are identical
colnames(dila)==colnames(dors)
colnames(dila)==colnames(plac)
# both true; by the transitive property, the column names for dors and plac are
# identical

# set temperatures and families as factors
dila$Pareja<-as.factor(dila$Pareja)
dila$Tratamiento_temperatura_C<-as.factor(dila$Tratamiento_temperatura_C)
dors$Pareja<-as.factor(dors$Pareja)
dors$Tratamiento_temperatura_C<-as.factor(dors$Tratamiento_temperatura_C)
plac$Pareja<-as.factor(plac$Pareja)
plac$Tratamiento_temperatura_C<-as.factor(plac$Tratamiento_temperatura_C)

# combine the three species into a single dataframe
larval_set1<-rbind(dila, dors, plac)

# Data fields -----------------------------------------------------------
# X: Row number from original dataset containing additional species and
# eggs that did not produce larvae (Numeric)
# sheet: Name of original Excel spreadsheet containing data for that population
# (Character)
# Consecutivo: Consecutive row number in the original Excel spreadsheet (Numeric)
# Especie: Beetle species (Categorical)
# Pareja: Parental pair that produced an individual (Categorical)
# Localidad: Locality where the parental beetles were collected (Categorical)
# Temperatura_pareja: Temperature at which the parental pair was kept
# (Constant, 22°C)
# Numero_huevo: Consecutive egg number for a parental pair (Numeric)
# Tratamiento_temperatura_C: Temperature at which an individual was reared
# (Categorical)
# Fecha_oviposition: Date the egg was laid (Date:MM/DD/YYYY)
# Fecha_eclosion: Date the egg hatched (Date:MM/DD/YYYY)
# Huevo_causa_muerte: Cause of egg death; since these datasets only include
# eggs that hatched, this is blank for all individuals (Character)
# Fecha_muerte_larva: Date an individual died as a larva, if applicable (Date:MM/DD/YYYY)
# Fecha_pupado: Date an indvidual pupated, if applicable (Date:MM/DD/YYYY)
# Fecha_nace_adulto: Date an individual emerged as an adult, if applicable (Date:MM/DD/YYYY)
# Tamano_dia_1_mm: Length of one-day-old larva in millimeters (Numeric)
# Tamano_dia_15_mm: Length of sixteen-day-old larva in millimeters, if
# applicable. This was originally intended to be the length of fifteen-day-old
# larvae, but the function for calculating the measurement date was set up as
# Fecha_eclosion+15, and this mistake was only detected midway through the
# experiment. (Numeric)
# Peso_pupa_g: Mass of pupa in grams, if applicable (Numeric)
# Tamano_pupa_mm: Length of pupa in millimeters, if applicable (Numeric)
# Fecha_muerte_pupa: Date an individual died as a pupa, if applicable (Date:MM/DD/YYYY)
# Peso_adulto_g: Mass of adult in grams, if applicable (Numeric)
# Tamano_adulto_mm: Length of adult in millimeters (Numeric)
# Sexo_adulto: Sex of adult, if applicable (Categorical)
# Fecha_muerte_adulto: Date an individual died as an adult, if applicable (Date:MM/DD/YYYY)
# Notas: Notes about an individual's condition, disappearance, etc., if applicable
# (Character)
# Fecha_mediacion: Date for sixteen-day-old larval measurement, if applicable (Date:MM/DD/YYYY)
# Comentarios.revisar.datos.Abril.2017: Any notes from a review of the data
# entered in April 2017 (Character)
# hatch: Whether an individual did (1) or did not (0) hatch (Logical)
# pupate: Whether an individual did (1) or did not (0) pupate (Logical)
# emerge: Whether an individual did (1) or did not (0) emerge as an adult (Logical)
# dev_time: Days from hatching to adult emergence, if applicable (Numeric)
# Asurv_time: Days an individual survived as an adult, if applicable (Numeric)
# hatch_yrwk: Year and week number that an individual hatched in (Date)
# death_yrwk: Year and week number that an individual died in (Date)
# emerge_yrwk: Year and week number that an individual emerged as an adult in,
# if applicable (Date)
# ldays: Days as a larva (Numeric)
# pdays: Days as a pupa, if applicable (Numeric)
# adays: Days as an adult, if applicable (Numeric)
# Checked_CB: Notes from when the record was checked by Christina S. Baer in November 2018
# prior to beginning data analysis (Character)



# Tables 1-3: Sample sizes, ANOVAs, and t-tests for developmental traits ------------------------------------
# Goal: Report sample sizes for developmental traits and test whether those 
# traits are affected by temperature.

# Following Garcia-Robledo and Horvitz (2012), lengths and masses will be 
# log-transformed for ANOVAs and t-tests.


# Sample sizes for developmental traits -----------------------------------
# sample sizes for length traits
aggregate(Tamano_dia_1_mm~Tratamiento_temperatura_C+Especie, data=larval_set1, na.action = na.omit, FUN=length)
aggregate(Tamano_dia_15_mm~Tratamiento_temperatura_C+Especie, data=larval_set1, na.action = na.omit, FUN=length)
aggregate(Tamano_pupa_mm~Tratamiento_temperatura_C+Especie, data=larval_set1, na.action = na.omit, FUN=length)
aggregate(Tamano_adulto_mm~Tratamiento_temperatura_C+Especie, data=larval_set1, na.action = na.omit, FUN=length)

# sample sizes for mass traits
aggregate(Peso_pupa_g~Tratamiento_temperatura_C+Especie, data=larval_set1, na.action = na.omit, FUN=length)
aggregate(Peso_adulto_g~Tratamiento_temperatura_C+Especie, data=larval_set1, na.action = na.omit, FUN=length)

# sample sizes for development time traits
aggregate(ldays~Tratamiento_temperatura_C+Especie, data=larval_set1[larval_set1$pupate==1,], na.action = na.omit, FUN=length)
aggregate(pdays~Tratamiento_temperatura_C+Especie, data=larval_set1[larval_set1$emerge==1,], na.action = na.omit, FUN=length)


# dilaticollis ANOVAs and t-tests (Table 1) -----------------
# Day 1 length
dila_1<-aov(log(Tamano_dia_1_mm)~Tratamiento_temperatura_C, 
            data=larval_set1[larval_set1$Especie=="C_dilaticollis",])
summary(dila_1)
TukeyHSD(dila_1)
# 30 and 15 are the same, but all other comparisons are different

# Day 16 length
dila_15<-aov(log(Tamano_dia_15_mm)~Tratamiento_temperatura_C, 
             data=larval_set1[larval_set1$Especie=="C_dilaticollis",])
summary(dila_15)
TukeyHSD(dila_15)
# all different

# # growth (controlling for length at Day 1), not included in table
# dila_g<-aov(log(Tamano_dia_15_mm-Tamano_dia_1_mm)~Tratamiento_temperatura_C, 
#             data=larval_set1[larval_set1$Especie=="C_dilaticollis",])
# summary(dila_g)
# TukeyHSD(dila_g)
# # all different

# pupal length
dila_p<-aov(log(Tamano_pupa_mm)~Tratamiento_temperatura_C, 
            data=larval_set1[larval_set1$Especie=="C_dilaticollis",])
summary(dila_p)
TukeyHSD(dila_p)
# 25 and 20 are different, 20 and 15 are borderline similar (p = 0.063), 
# 25 and 15 are the same

# adult length
# dila_a<-aov(log(Tamano_adulto_mm)~Tratamiento_temperatura_C, data=larval_set1[larval_set1$Especie=="C_dilaticollis",])
# summary(dila_a)
# TukeyHSD(dila_a)
# # only 25 and 20 are different

# April 22, 2019: Because there are so few 15C adults (only 5), this analysis
# should really be a t-test between 20C and 25C
t.test(log(dila$Tamano_adulto_mm[dila$Tratamiento_temperatura_C==20]),
       log(dila$Tamano_adulto_mm[dila$Tratamiento_temperatura_C==25]))
# significantly different

# pupal mass
dila_pp<-aov(log(Peso_pupa_g)~Tratamiento_temperatura_C, 
             data=larval_set1[larval_set1$Especie=="C_dilaticollis",])
summary(dila_pp)
TukeyHSD(dila_pp)
# the weights for 20 and 25 are different, but both overlap with 15 (2
# overlapping groups)

# adult mass
# dila_ap<-aov(log(Peso_adulto_g)~Tratamiento_temperatura_C, 
#              data=larval_set1[larval_set1$Especie=="C_dilaticollis",])
# summary(dila_ap)
# TukeyHSD(dila_ap)
# # only the weights for 20 and 25 are different, but both overlap with 15

# April 22, 2019: Because there are so few 15C adults (only 5), this analysis
# should really be a t-test between 20C and 25C
t.test(log(dila$Peso_adulto_g[dila$Tratamiento_temperatura_C==20]),
       log(dila$Peso_adulto_g[dila$Tratamiento_temperatura_C==25]))
# significantly different


# development times
# larval development times
dila_dl<-aov(ldays~Tratamiento_temperatura_C, 
             data=larval_set1[larval_set1$Especie=="C_dilaticollis" & larval_set1$pupate==1,])
summary(dila_dl)
TukeyHSD(dila_dl)
# all different

# pupal development times
# dila_dp<-aov(pdays~Tratamiento_temperatura_C, 
#              data=larval_set1[larval_set1$Especie=="C_dilaticollis" & larval_set1$emerge==1,])
# summary(dila_dp)
# TukeyHSD(dila_dp)
# # all different

# April 22, 2019: Because there are so few 15C adults (only 5), this analysis
# should really be a t-test between 20C and 25C
t.test(dila$pdays[dila$Tratamiento_temperatura_C==20],
       dila$pdays[dila$Tratamiento_temperatura_C==25])
# significantly different



# dorsalis ANOVAs and t-tests (Table 2) -----------------------------------
# Day 1 length
dors_1<-aov(log(Tamano_dia_1_mm)~Tratamiento_temperatura_C, 
            data=larval_set1[larval_set1$Especie=="C_dorsalis",])
summary(dors_1)
TukeyHSD(dors_1)
# 30 and 15 are the same, and 25 and 20 are the same (2 groups)

# Day 16 length
dors_15<-aov(log(Tamano_dia_15_mm)~Tratamiento_temperatura_C, 
             data=larval_set1[larval_set1$Especie=="C_dorsalis",])
summary(dors_15)
TukeyHSD(dors_15)
# 30 and 20 are the same, all other comparisons are different (3 groups)

# # growth (controlling for length at Day 1), not included in Table 2
# dors_g<-aov(log(Tamano_dia_15_mm-Tamano_dia_1_mm)~Tratamiento_temperatura_C,
#             data=larval_set1[larval_set1$Especie=="C_dorsalis",])
# summary(dors_g)
# TukeyHSD(dors_g)
# # 30 and 20 are the same, all other comparisons are different (3 groups)

# pupal length
dors_p<-aov(log(Tamano_pupa_mm)~Tratamiento_temperatura_C, 
            data=larval_set1[larval_set1$Especie=="C_dorsalis",])
summary(dors_p)
TukeyHSD(dors_p)
# all different

# adult length
# dors_a<-aov(log(Tamano_adulto_mm)~Tratamiento_temperatura_C, 
#             data=larval_set1[larval_set1$Especie=="C_dorsalis",])
# summary(dors_a)
# TukeyHSD(dors_a)
# # 20 and 25 are the same, but different from 15 (2 groups)

# April 22, 2019: Because there are so few 15C adults (only 4), this analysis
# should really be a t-test between 20C and 25C
t.test(log(dors$Tamano_adulto_mm[dors$Tratamiento_temperatura_C==20]),
       log(dors$Tamano_adulto_mm[dors$Tratamiento_temperatura_C==25]))
# not significantly different

# pupal mass
dors_pp<-aov(log(Peso_pupa_g)~Tratamiento_temperatura_C, 
             data=larval_set1[larval_set1$Especie=="C_dorsalis",])
summary(dors_pp)
TukeyHSD(dors_pp)
# 20 and 25 weights are the same, but different from 15 (2 groups)

# adult mass
# dors_ap<-aov(log(Peso_adulto_g)~Tratamiento_temperatura_C, 
#              data=larval_set1[larval_set1$Especie=="C_dorsalis",])
# summary(dors_ap)
# # no effect of temperature on adult weight

# April 22, 2019: Because there are so few 15C adults (only 4), this analysis
# should really be a t-test between 20C and 25C
t.test(log(dors$Peso_adulto_g[dors$Tratamiento_temperatura_C==20]),
       log(dors$Peso_adulto_g[dors$Tratamiento_temperatura_C==25]))
# not significantly different


# development times
# larval development time
dors_dl<-aov(ldays~Tratamiento_temperatura_C, 
             data=larval_set1[larval_set1$Especie=="C_dorsalis" & larval_set1$pupate==1,])
summary(dors_dl)
TukeyHSD(dors_dl)
# all different

# pupal development time
# dors_dp<-aov(pdays~Tratamiento_temperatura_C, data=larval_set1[larval_set1$Especie=="C_dorsalis" & larval_set1$emerge==1,])
# summary(dors_dp)
# TukeyHSD(dors_dp)
# # all different

# April 22, 2019: Because there are so few 15C adults (only 4), this analysis
# should really be a t-test between 20C and 25C
t.test(dors$pdays[dors$Tratamiento_temperatura_C==20],
       dors$pdays[dors$Tratamiento_temperatura_C==25])



# placida ANOVAs (Table 3) ------------------------------------
# Day 1 length
plac_1<-aov(log(Tamano_dia_1_mm)~Tratamiento_temperatura_C, 
            data=larval_set1[larval_set1$Especie=="C_placida",])
summary(plac_1)
TukeyHSD(plac_1)
# 30 and 15 are the same, and 25 and 20 are the same (2 groups)

# Day 16 length
plac_15<-aov(log(Tamano_dia_15_mm)~Tratamiento_temperatura_C, 
             data=larval_set1[larval_set1$Especie=="C_placida",])
summary(plac_15)
TukeyHSD(plac_15)
# all different

# growth (controlling for length at day 1), not included in Table 3
# this doesn't work because there is one individual (Pareja 6, Huevo 42) that 
# didn't grow at all (you can't take the log of 0). 
# plac_g<-aov(log(Tamano_dia_15_mm-Tamano_dia_1_mm)~Tratamiento_temperatura_C, 
#             data=larval_set1[larval_set1$Especie=="C_placida",])
# summary(plac_g)
# TukeyHSD(plac_g)

# pupal length
plac_p<-aov(log(Tamano_pupa_mm)~Tratamiento_temperatura_C, 
            data=larval_set1[larval_set1$Especie=="C_placida",])
summary(plac_p)
TukeyHSD(plac_p)
# 25 and 15 the same, 25 and 20 are different, and 20-15 p = 0.051 (2
# overlapping groups)

# adult length
plac_a<-aov(log(Tamano_adulto_mm)~Tratamiento_temperatura_C, 
            data=larval_set1[larval_set1$Especie=="C_placida",])
summary(plac_a)
TukeyHSD(plac_a)
# 25 and 20 different, but 15 overlaps with both (2 overlapping groups)

# pupal mass
plac_pp<-aov(log(Peso_pupa_g)~Tratamiento_temperatura_C, 
             data=larval_set1[larval_set1$Especie=="C_placida",])
summary(plac_pp)
TukeyHSD(plac_pp)
# 20 and 15 weights are the same, but different from 25 (2 groups)

# adult mass
plac_ap<-aov(log(Peso_adulto_g)~Tratamiento_temperatura_C, 
             data=larval_set1[larval_set1$Especie=="C_placida",])
summary(plac_ap)
TukeyHSD(plac_ap)
# 25 and 15 weights are the same, but different from 20 (2 groups)

# development times
# larval development time
plac_dl<-aov(ldays~Tratamiento_temperatura_C, 
             data=larval_set1[larval_set1$Especie=="C_placida" & larval_set1$pupate==1,])
summary(plac_dl)
TukeyHSD(plac_dl)
# all comparisons different

# adult development time
plac_dp<-aov(pdays~Tratamiento_temperatura_C, 
             data=larval_set1[larval_set1$Especie=="C_placida" & larval_set1$emerge==1,])
summary(plac_dp)
TukeyHSD(plac_dp)
# all different

# Fig. 1: Length traits ---------------------------------------------------
# Goal: Display mean +- SE +- SD larval (Day 1 and Day 16), pupal, and adult
# lengths by temperature for each species

# Individual dilaticollis length plots (Fig. 1A-D)----------------------------
# Day 1 length plot (Fig. 1A)
p_L1_dila<-ggplot(larval_set1[larval_set1$Especie=="C_dilaticollis",],
                  aes(x = Tratamiento_temperatura_C, y = Tamano_dia_1_mm, 
                      color = Tratamiento_temperatura_C))+ 
  stat_summary(fun.data=MeanSESD, geom="boxplot")+ theme_classic2() + 
  expand_limits(y=c(1.8,3))+
  scale_y_continuous(breaks=pretty_breaks(3)) +
  scale_x_discrete(limits=c("15","20","25","30"))+
  scale_color_manual(breaks=c(15,20,25,30),labels=c("15","20","25","30"),
                    values=c("turquoise", "yellow","orange","red"))+
  labs(x="Temperature (°C)", y="Day 1 length (mm)", title="dilaticollis", 
       family="sans")+
  guides(color=FALSE)

p_L1_dila

# Day 16 length plot (Fig. 1B)
p_L15_dila<-ggplot(larval_set1[larval_set1$Especie=="C_dilaticollis",],
                  aes(x = Tratamiento_temperatura_C, y = Tamano_dia_15_mm, 
                      color = Tratamiento_temperatura_C)) + 
  stat_summary(fun.data=MeanSESD, geom="boxplot")+ theme_classic2() + 
  expand_limits(y=c(2.5,6.5))+
  scale_y_continuous(breaks=pretty_breaks(4)) +
  scale_x_discrete(limits=c("15","20","25","30"))+
  scale_color_manual(breaks=c(15,20,25,30),labels=c("15","20","25","30"),
                     values=c("turquoise", "yellow","orange"))+
  labs(x="Temperature (°C)", y="Day 16 length (mm)", title="dilaticollis", 
       family="sans")+
  guides(color=FALSE)

p_L15_dila

# Pupal length plot (Fig. 1C)
p_LP_dila<-ggplot(larval_set1[larval_set1$Especie=="C_dilaticollis",],
                   aes(x = Tratamiento_temperatura_C, y = Tamano_pupa_mm, 
                       color = Tratamiento_temperatura_C))+ 
  stat_summary(fun.data=MeanSESD, geom="boxplot")+ theme_classic2() + 
  expand_limits(y=c(5.5,8))+
  scale_y_continuous(breaks=pretty_breaks(3)) +
  scale_x_discrete(limits=c("15","20","25","30"))+
  scale_color_manual(breaks=c(15,20,25,30),labels=c("15","20","25","30"),
                     values=c("turquoise", "yellow","orange","red"))+
  labs(x="Temperature (°C)", y="Pupal length (mm)", title="dilaticollis", 
       family="sans")+
  guides(color=FALSE)

p_LP_dila 

# Adult length plot (Fig. 1D)
p_LA_dila<-ggplot(larval_set1[larval_set1$Especie=="C_dilaticollis",],
                  aes(x = Tratamiento_temperatura_C, y = Tamano_adulto_mm, 
                      color = Tratamiento_temperatura_C))+ 
  stat_summary(fun.data=MeanSESD, geom="boxplot")+ theme_classic2() + 
  expand_limits(y=c(3.8,6.4))+
  scale_y_continuous(breaks=pretty_breaks(3)) +
  scale_x_discrete(limits=c("15","20","25","30"))+
  scale_color_manual(breaks=c(15,20,25,30),labels=c("15","20","25","30"),
                     values=c("turquoise", "yellow","orange","red"))+
  labs(x="Temperature (°C)", y="Adult length (mm)", title="dilaticollis", 
       family="sans")+
  guides(color=FALSE)

p_LA_dila 

# Individual dorsalis length plots (Fig. 1E-H)--------------------------------
# Day 1 length plot (Fig. 1E)
p_L1_dors<-ggplot(larval_set1[larval_set1$Especie=="C_dorsalis",],
                  aes(x = Tratamiento_temperatura_C, y = Tamano_dia_1_mm, 
                      color = Tratamiento_temperatura_C))+ 
  stat_summary(fun.data=MeanSESD, geom="boxplot")+ theme_classic2() + 
  expand_limits(y=c(1.8,3))+
  scale_y_continuous(breaks=pretty_breaks(3)) +
  scale_x_discrete(limits=c("15","20","25","30"))+
  scale_color_manual(breaks=c(15,20,25,30),labels=c("15","20","25","30"),
                     values=c("turquoise", "yellow","orange","red"))+
  labs(x="Temperature (°C)", y="Day 1 length (mm)", title="dorsalis", 
       family="sans")+
  guides(color=FALSE)

p_L1_dors 

# Day 16 length plot (Fig. 1F)
p_L15_dors<-ggplot(larval_set1[larval_set1$Especie=="C_dorsalis",],
                  aes(x = Tratamiento_temperatura_C, y = Tamano_dia_15_mm, 
                      color = Tratamiento_temperatura_C))+ 
  stat_summary(fun.data=MeanSESD, geom="boxplot")+ theme_classic2() + 
  expand_limits(y=c(2.5,6.5))+
  scale_y_continuous(breaks=pretty_breaks(4)) +
  scale_x_discrete(limits=c("15","20","25","30"))+
  scale_color_manual(breaks=c(15,20,25,30),labels=c("15","20","25","30"),
                     values=c("turquoise", "yellow","orange","red"))+
  labs(x="Temperature (°C)", y="Day 16 length (mm)", title="dorsalis", 
       family="sans")+
  guides(color=FALSE)

p_L15_dors 

# Pupal length plot (Fig. 1G)
p_LP_dors<-ggplot(larval_set1[larval_set1$Especie=="C_dorsalis",],
                  aes(x = Tratamiento_temperatura_C, y = Tamano_pupa_mm, 
                      color = Tratamiento_temperatura_C))+ 
  stat_summary(fun.data=MeanSESD, geom="boxplot")+ theme_classic2() + 
  expand_limits(y=c(5.5,8))+
  scale_y_continuous(breaks=pretty_breaks(3)) +
  scale_x_discrete(limits=c("15","20","25","30"))+
  scale_color_manual(breaks=c(15,20,25,30),labels=c("15","20","25","30"),
                     values=c("turquoise", "yellow","orange","red"))+
  labs(x="Temperature (°C)", y="Pupal length (mm)", title="dorsalis", 
       family="sans")+
  guides(color=FALSE)

p_LP_dors 

# Adult length plot (Fig. 1H)
p_LA_dors<-ggplot(larval_set1[larval_set1$Especie=="C_dorsalis",],
                  aes(x = Tratamiento_temperatura_C, y = Tamano_adulto_mm, 
                      color = Tratamiento_temperatura_C))+ 
  stat_summary(fun.data=MeanSESD, geom="boxplot")+ theme_classic2() + 
  expand_limits(y=c(3.8,6.4))+
  scale_y_continuous(breaks=pretty_breaks(3)) +
  scale_x_discrete(limits=c("15","20","25","30"))+
  scale_color_manual(breaks=c(15,20,25,30),labels=c("15","20","25","30"),
                     values=c("turquoise", "yellow","orange","red"))+
  labs(x="Temperature (°C)", y="Adult length (mm)", title="dorsalis", 
       family="sans")+
  guides(color=FALSE)

p_LA_dors 

# Individual placida length plots (Fig. 1I-L)---------------------------------
# Day 1 length plot (Fig. 1I)
p_L1_plac<-ggplot(larval_set1[larval_set1$Especie=="C_placida",],
                  aes(x = Tratamiento_temperatura_C, y = Tamano_dia_1_mm, 
                      color = Tratamiento_temperatura_C))+ 
  stat_summary(fun.data=MeanSESD, geom="boxplot")+ theme_classic2() + 
  expand_limits(y=c(1.8,3))+
  scale_y_continuous(breaks=pretty_breaks(3)) + 
  scale_x_discrete(limits=c("15","20","25","30"))+
  scale_color_manual(breaks=c(15,20,25,30),labels=c("15","20","25","30"),
                     values=c("turquoise", "yellow","orange","red"))+
  labs(x="Temperature (°C)", y="Day 1 length (mm)", title="placida", 
       family="sans")+
  guides(color=FALSE)

p_L1_plac 

# Day 16 length plot (Fig. 1J)
p_L15_plac<-ggplot(larval_set1[larval_set1$Especie=="C_placida",],
                   aes(x = Tratamiento_temperatura_C, y = Tamano_dia_15_mm, 
                       color = Tratamiento_temperatura_C))+ 
  stat_summary(fun.data=MeanSESD, geom="boxplot")+ theme_classic2() + 
  expand_limits(y=c(2.5,6.5))+
  scale_y_continuous(breaks=pretty_breaks(4)) +
  scale_x_discrete(limits=c("15","20","25","30"))+
  scale_color_manual(breaks=c(15,20,25,30),labels=c("15","20","25","30"),
                     values=c("turquoise", "yellow","orange","red"))+
  labs(x="Temperature (°C)", y="Day 16 length (mm)", title="placida", 
       family="sans")+
  guides(color=FALSE)

p_L15_plac 

# Pupal length plot (Fig. 1K)
p_LP_plac<-ggplot(larval_set1[larval_set1$Especie=="C_placida",],
                  aes(x = Tratamiento_temperatura_C, y = Tamano_pupa_mm, 
                      color = Tratamiento_temperatura_C))+ 
  stat_summary(fun.data=MeanSESD, geom="boxplot")+ theme_classic2() + 
  expand_limits(y=c(5.5,8))+
  scale_y_continuous(breaks=pretty_breaks(3)) +
  scale_x_discrete(limits=c("15","20","25","30"))+
  scale_color_manual(breaks=c(15,20,25,30),labels=c("15","20","25","30"),
                     values=c("turquoise", "yellow","orange","red"))+
  labs(x="Temperature (°C)", y="Pupal length (mm)", title="placida", 
       family="sans")+
  guides(color=FALSE)

p_LP_plac 

# Adult length plot (Fig. 1L)
p_LA_plac<-ggplot(larval_set1[larval_set1$Especie=="C_placida",],
                  aes(x = Tratamiento_temperatura_C, y = Tamano_adulto_mm, 
                      color = Tratamiento_temperatura_C)) +
  stat_summary(fun.data=MeanSESD, geom="boxplot")+ theme_classic2() + 
  expand_limits(y=c(3.8,6.4))+
  scale_y_continuous(breaks=pretty_breaks(3)) +
  scale_x_discrete(limits=c("15","20","25","30"))+
  scale_color_manual(breaks=c(15,20,25,30),labels=c("15","20","25","30"),
                     values=c("turquoise", "yellow","orange","red"))+
  labs(x="Temperature (°C)", y="Adult length (mm)", title="placida", 
       family="sans")+
  guides(color=FALSE)

p_LA_plac

# save Figure 1 length plots ----------------------------------------------
# ggsave("Fig1_new_dila_Day1.pdf",plot=p_L1_dila,device="pdf")
# ggsave("Fig1_new_dila_Day16.pdf",plot=p_L16_dila,device="pdf")
# ggsave("Fig1_new_dila_Pupae.pdf",plot=p_LP_dila,device="pdf")
# ggsave("Fig1_new_dila_Adults.pdf",plot=p_LA_dila,device="pdf")

# ggsave("Fig1_new_dors_Day1.pdf",plot=p_L1_dors,device="pdf")
# ggsave("Fig1_new_dors_Day16.pdf",plot=p_L16_dors,device="pdf")
# ggsave("Fig1_new_dors_Pupae.pdf",plot=p_LP_dors,device="pdf")
# ggsave("Fig1_new_dors_Adults.pdf",plot=p_LA_dors,device="pdf")

# ggsave("Fig1_new_plac_Day1.pdf",plot=p_L1_plac,device="pdf")
# ggsave("Fig1_new_plac_Day16.pdf",plot=p_L16_plac,device="pdf")
# ggsave("Fig1_new_plac_Pupae.pdf",plot=p_LP_plac,device="pdf")
# ggsave("Fig1_new_plac_Adults.pdf",plot=p_LA_plac,device="pdf")

# Fig. 2: Mass traits ---------------------------------------------------
# Goal: Display mean +- SE +- SD pupal and adult lengths by temperature for each
# species. 
# Masses are graphed in milligrams rather than grams

# Individual dilaticollis masses (Fig. 2A-B)----------------------------------
# Pupal mass  plot (Fig. 2A)
p_MP_dila<-ggplot(larval_set1[larval_set1$Especie=="C_dilaticollis",],
                  aes(x = Tratamiento_temperatura_C, y = Peso_pupa_g*1000, 
                      color = Tratamiento_temperatura_C))+ 
  stat_summary(fun.data=MeanSESD, geom="boxplot")+ theme_classic2() + 
  expand_limits(y=c(5,21))+
  scale_y_continuous(breaks=pretty_breaks(3)) +
  scale_x_discrete(limits=c("15","20","25"))+
  scale_color_manual(breaks=c(15,20,25),labels=c("15","20","25"),
                     values=c("turquoise", "yellow","orange"))+
  labs(x="Temperature (°C)", y="Pupal mass (mg)", title="dilaticollis", 
       family="sans")+
  guides(color=FALSE)

p_MP_dila 

# Adult mass  plot (Fig. 2B)
p_MA_dila<-ggplot(larval_set1[larval_set1$Especie=="C_dilaticollis",],
                  aes(x = Tratamiento_temperatura_C, y = Peso_adulto_g*1000, 
                      color = Tratamiento_temperatura_C))+ 
  stat_summary(fun.data=MeanSESD, geom="boxplot")+ theme_classic2() + 
  expand_limits(y=c(3.5,12.5))+
  scale_y_continuous(breaks=pretty_breaks(4)) +
  scale_x_discrete(limits=c("15","20","25"))+
  scale_color_manual(breaks=c(15,20,25),labels=c("15","20","25"),
                     values=c("turquoise", "yellow","orange"))+
  labs(x="Temperature (°C)", y="Adult mass (mg)", title="dilaticollis", 
       family="sans")+
  guides(color=FALSE)

p_MA_dila 


# Individual dorsalis masses (Fig. 2C-D)------------------------------------
# Pupal mass  plot (Fig. 2C)
p_MP_dors<-ggplot(larval_set1[larval_set1$Especie=="C_dorsalis",],
                  aes(x = Tratamiento_temperatura_C, y = Peso_pupa_g*1000, 
                      color = Tratamiento_temperatura_C))+ 
  stat_summary(fun.data=MeanSESD, geom="boxplot")+ theme_classic2() + 
  expand_limits(y=c(5,21))+
  scale_y_continuous(breaks=pretty_breaks(3)) +
  scale_x_discrete(limits=c("15","20","25"))+
  scale_color_manual(breaks=c(15,20,25),labels=c("15","20","25"),
                     values=c("turquoise", "yellow","orange"))+
  labs(x="Temperature (°C)", y="Pupal mass (mg)", title="dorsalis", 
       family="sans")+
  guides(color=FALSE)

p_MP_dors 

# Adult mass  plot (Fig. 2D)
p_MA_dors<-ggplot(larval_set1[larval_set1$Especie=="C_dorsalis",],
                  aes(x = Tratamiento_temperatura_C, y = Peso_adulto_g*1000, 
                      color = Tratamiento_temperatura_C))+ 
  stat_summary(fun.data=MeanSESD, geom="boxplot")+ theme_classic2() + 
  expand_limits(y=c(3.5,12.5))+
  scale_y_continuous(breaks=pretty_breaks(4)) +
  scale_x_discrete(limits=c("15","20","25"))+
  scale_color_manual(breaks=c(15,20,25),labels=c("15","20","25"),
                     values=c("turquoise", "yellow","orange"))+
  labs(x="Temperature (°C)", y="Adult mass (mg)", title="dorsalis", 
       family="sans")+
  guides(color=FALSE)

p_MA_dors 

# Individual placida masses (Fig. 2E-F)----------------------------------
# Pupal mass  plot (Fig. 2E)
p_MP_plac<-ggplot(larval_set1[larval_set1$Especie=="C_placida",],
                  aes(x = Tratamiento_temperatura_C, y = Peso_pupa_g*1000, 
                      color = Tratamiento_temperatura_C))+ 
  stat_summary(fun.data=MeanSESD, geom="boxplot")+ theme_classic2() + 
  expand_limits(y=c(5,21))+
  scale_y_continuous(breaks=pretty_breaks(3)) +
  scale_x_discrete(limits=c("15","20","25"))+
  scale_color_manual(breaks=c(15,20,25),labels=c("15","20","25"),
                     values=c("turquoise", "yellow","orange"))+
  labs(x="Temperature (°C)", y="Pupal mass (mg)", title="placida", 
       family="sans")+
  guides(color=FALSE)

p_MP_plac 

# Adult mass  plot (Fig. 2F)
p_MA_plac<-ggplot(larval_set1[larval_set1$Especie=="C_placida",],
                  aes(x = Tratamiento_temperatura_C, y = Peso_adulto_g*1000, 
                      color = Tratamiento_temperatura_C))+ 
  stat_summary(fun.data=MeanSESD, geom="boxplot")+ theme_classic2() + 
  expand_limits(y=c(3.5,12.5))+
  scale_y_continuous(breaks=pretty_breaks(4)) +
  scale_x_discrete(limits=c("15","20","25"))+
  scale_color_manual(breaks=c(15,20,25),labels=c("15","20","25"),
                     values=c("turquoise", "yellow","orange"))+
  labs(x="Temperature (°C)", y="Adult mass (mg)", title="placida", 
       family="sans")+
  guides(color=FALSE)

p_MA_plac 

# save Figure 2 mass plots ----------------------------------------------
# ggsave("Fig2_new_dila_Pupae.pdf",plot=p_MP_dila,device="pdf")
# ggsave("Fig2_new_dila_Adults.pdf",plot=p_MA_dila,device="pdf")

# ggsave("Fig2_new_dors_Pupae.pdf",plot=p_MP_dors,device="pdf")
# ggsave("Fig2_new_dors_Adults.pdf",plot=p_MA_dors,device="pdf")

# ggsave("Fig2_new_plac_Pupae.pdf",plot=p_MP_plac,device="pdf")
# ggsave("Fig2_new_plac_Adults.pdf",plot=p_MA_plac,device="pdf")
# Fig. 3: Development times -----------------------------------------------
# Goal: Display mean +- SE +- SD development times for larvae and pupae of each
# species.

# dilaticollis development times (Fig. 3A-B)---------------------------------
# larval development plot (Fig. 3A)
p_ld_dila<-ggplot(larval_set1[larval_set1$Especie=="C_dilaticollis" & larval_set1$pupate==1 &
                                larval_set1$Tratamiento_temperatura_C!="NA",],
                  aes(x = Tratamiento_temperatura_C, y = ldays))+ 
  stat_summary(fun.data=MeanSESD, geom="boxplot")+ theme_classic2() + 
  expand_limits(y=c(0,170)) + 
  scale_x_discrete(limits=c("15","20","25"))+
  labs(x="Temperature (°C)", y="Development time (days)", 
       title="dilaticollis larvae", family="sans")

p_ld_dila 

# pupal development plot (Fig. 3B)
p_pd_dila<-ggplot(larval_set1[larval_set1$Especie=="C_dilaticollis"& larval_set1$emerge==1,],
                  aes(x = Tratamiento_temperatura_C, y = pdays))+ 
  stat_summary(fun.data=MeanSESD, geom="boxplot")+ theme_classic2() + 
  expand_limits(y=c(0,90)) + 
  scale_x_discrete(limits=c("15","20","25"))+
  labs(x="Temperature (°C)", y="Development time (days)", 
       title="dilaticollis pupae", family="sans")

p_pd_dila 

# dorsalis development times (Fig. 3C-D)-------------------------------------
# larval development plot (Fig. 3C)
p_ld_dors<-ggplot(larval_set1[larval_set1$Especie=="C_dorsalis"& larval_set1$pupate==1,],
                  aes(x = Tratamiento_temperatura_C, y = ldays))+ 
  stat_summary(fun.data=MeanSESD, geom="boxplot")+ theme_classic2() + 
  expand_limits(y=c(0,170)) + 
  scale_x_discrete(limits=c("15","20","25"))+
  labs(x="Temperature (°C)", y="Development time (days)", 
       title="dorsalis larvae", family="sans")

p_ld_dors 

# pupal development plot (Fig. 3D)
p_pd_dors<-ggplot(larval_set1[larval_set1$Especie=="C_dorsalis"& larval_set1$emerge==1,],
                  aes(x = Tratamiento_temperatura_C, y = pdays))+ 
  stat_summary(fun.data=MeanSESD, geom="boxplot")+ theme_classic2() + 
  expand_limits(y=c(0,90)) + 
  scale_x_discrete(limits=c("15","20","25"))+
  labs(x="Temperature (°C)", y="Development time (days)", 
       title="dorsalis pupae", family="sans")

p_pd_dors 

# placida development times (Fig. 3E-F)--------------------------------------
# larval development plot (Fig. 3E)
p_ld_plac<-ggplot(larval_set1[larval_set1$Especie=="C_placida"& larval_set1$pupate==1,],
                  aes(x = Tratamiento_temperatura_C, y = ldays))+ 
  stat_summary(fun.data=MeanSESD, geom="boxplot")+ theme_classic2() + 
  expand_limits(y=c(0,170)) + 
  scale_x_discrete(limits=c("15","20","25"))+
  labs(x="Temperature (°C)", y="Development time (days)", 
       title="placida larvae", family="sans")

p_ld_plac 

# pupal development plot (Fig. 3F)
p_pd_plac<-ggplot(larval_set1[larval_set1$Especie=="C_placida"& larval_set1$emerge==1,],
                  aes(x = Tratamiento_temperatura_C, y = pdays))+ 
  stat_summary(fun.data=MeanSESD, geom="boxplot")+ theme_classic2() + 
  expand_limits(y=c(0,90)) + 
  scale_x_discrete(limits=c("15","20","25"))+
  labs(x="Temperature (°C)", y="Development time (days)", 
       title="placida pupae", family="sans")

p_pd_plac 

# save Figure 3 development time plots ----------------------------------------------
# ggsave("Fig3_new_dila_LDT.pdf",plot=p_ld_dila,device="pdf")
# ggsave("Fig3_new_dila_PDT.pdf",plot=p_pd_dila,device="pdf")

# ggsave("Fig3_new_dors_LDT.pdf",plot=p_ld_dors,device="pdf")
# ggsave("Fig3_new_dors_PDT.pdf",plot=p_pd_dors,device="pdf")

# ggsave("Fig3_new_plac_LDT.pdf",plot=p_ld_plac,device="pdf")
# ggsave("Fig3_new_plac_PDT.pdf",plot=p_pd_plac,device="pdf")

# Subset balanced families for quantitative genetics ---------------------
# We need to know which families per species have at least 5 individuals each
# from 20 and 25C (and potentially other temperatures) for each trait, and how
# many individuals there are. The number of temperatures and families included
# in an analysis will depend on the trait in question.

#  dilaticollis families for quantitative genetics ------------------------
# Day 1 lengths
# count number of individuals per family (Pareja) and temperature with Day 1
# lengths
fam_dila_1<-aggregate(Tamano_dia_1_mm~Tratamiento_temperatura_C+Pareja,
  data=larval_set1[larval_set1$Especie=="C_dilaticollis",], length)
# only consider temperature*family combinations with >= 5 individuals
fam_dila_1[fam_dila_1$Tamano_dia_1_mm>=5,]
# for each family, count how many temperatures have >= 5 individuals
(fam_dila_1[fam_dila_1$Tamano_dia_1_mm>=5,]%>%
    count(Pareja))->dila_T1_fam_list
# only consider families that have individuals in at least 2 treatments 
dila_T1_fam_list<-dila_T1_fam_list[dila_T1_fam_list$n>=2,]

# Day 16 lengths
# count number of individuals per family (Pareja) and temperature with Day 16
# lengths
fam_dila_15<-aggregate(Tamano_dia_15_mm~Tratamiento_temperatura_C+Pareja,
                      data=larval_set1[larval_set1$Especie=="C_dilaticollis",], length)
# only consider temperature*family combinations with >= 5 individuals
fam_dila_15[fam_dila_15$Tamano_dia_15_mm>=5,]
# for each family, count how many temperatures have >= 5 individuals
(fam_dila_15[fam_dila_15$Tamano_dia_15_mm>=5,]%>%
  count(Pareja))->dila_T15_fam_list
# only consider families that have individuals in at least 2 treatments 
dila_T15_fam_list<-dila_T15_fam_list[dila_T15_fam_list$n>=2,]

# pupal lengths
# count number of individuals per family (Pareja) and temperature with pupal
# lengths
fam_dila_p<-aggregate(Tamano_pupa_mm~Tratamiento_temperatura_C+Pareja,
                       data=larval_set1[larval_set1$Especie=="C_dilaticollis",], length)
# only consider temperature*family combinations with >= 5 individuals
fam_dila_p[fam_dila_p$Tamano_pupa_mm>=5,]
# for each family, count how many temperatures have >= 5 individuals
(fam_dila_p[fam_dila_p$Tamano_pupa_mm>=5,]%>%
    count(Pareja))->dila_Tp_fam_list
# only consider families that have individuals in at least 2 treatments 
dila_Tp_fam_list<-dila_Tp_fam_list[dila_Tp_fam_list$n>=2,]

# pupal masses
# count number of individuals per family (Pareja) and temperature with pupal
# masses
fam_dila_p2<-aggregate(Peso_pupa_g~Tratamiento_temperatura_C+Pareja,
                      data=larval_set1[larval_set1$Especie=="C_dilaticollis",], length)
# only consider temperature*family combinations with >= 5 individuals
fam_dila_p2[fam_dila_p2$Peso_pupa_g>=5,]
# for each family, count how many temperatures have >= 5 individuals
(fam_dila_p2[fam_dila_p2$Peso_pupa_g>=5,]%>%
    count(Pareja))->dila_Pp_fam_list
# only consider families that have individuals in at least 2 treatments 
dila_Pp_fam_list<-dila_Pp_fam_list[dila_Pp_fam_list$n>=2,]

# compare the family lists for pupal lengths and pupal masses
dila_Tp_fam_list==dila_Pp_fam_list 
# The lists are the same, which makes sense since these measurements are paired.
# Going forward, I will assume that I only need one list each for pupae and
# adults.

# adult lengths (and masses)
# count number of individuals per family (Pareja) and temperature with adults
# lengths
fam_dila_a<-aggregate(Tamano_adulto_mm~Tratamiento_temperatura_C+Pareja,
                      data=larval_set1[larval_set1$Especie=="C_dilaticollis",], length)
# only consider temperature*family combinations with >= 5 individuals
fam_dila_a[fam_dila_a$Tamano_adulto_mm>=5,]
# for each family, count how many temperatures have >= 5 individuals
(fam_dila_a[fam_dila_a$Tamano_adulto_mm>=5,]%>%
    count(Pareja))->dila_a_fam_list
# only consider families that have individuals in at least 2 treatments 
dila_a_fam_list<-dila_a_fam_list[dila_a_fam_list$n>=2,]
dila_a_fam_list
# There are only four families that are balanced for 20-25



# dorsalis families for quantitative genetics -----------------------------
# Day 1 lengths
# count number of individuals per family (Pareja) and temperature with Day 1
# lengths
fam_dors_1<-aggregate(Tamano_dia_1_mm~Tratamiento_temperatura_C+Pareja,
                      data=larval_set1[larval_set1$Especie=="C_dorsalis",], length)
# only consider temperature*family combinations with >= 5 individuals
fam_dors_1[fam_dors_1$Tamano_dia_1_mm>=5,]
# for each family, count how many temperatures have >= 5 individuals
(fam_dors_1[fam_dors_1$Tamano_dia_1_mm>=5,]%>%
    count(Pareja))->dors_T1_fam_list
# only consider families that have individuals in at least 2 treatments 
dors_T1_fam_list<-dors_T1_fam_list[dors_T1_fam_list$n>=2,]

# Day 16 lengths
# count number of individuals per family (Pareja) and temperature with Day 16
# lengths
fam_dors_15<-aggregate(Tamano_dia_15_mm~Tratamiento_temperatura_C+Pareja,
                       data=larval_set1[larval_set1$Especie=="C_dorsalis",], length)
# only consider temperature*family combinations with >= 5 individuals
fam_dors_15[fam_dors_15$Tamano_dia_15_mm>=5,]
# for each family, count how many temperatures have >= 5 individuals
(fam_dors_15[fam_dors_15$Tamano_dia_15_mm>=5,]%>%
    count(Pareja))->dors_T15_fam_list
# only consider families that have individuals in at least 2 treatments 
dors_T15_fam_list<-dors_T15_fam_list[dors_T15_fam_list$n>=2,]

# pupal lengths (and masses)
# count number of individuals per family (Pareja) and temperature with pupal
# lengths
fam_dors_p<-aggregate(Tamano_pupa_mm~Tratamiento_temperatura_C+Pareja,
                      data=larval_set1[larval_set1$Especie=="C_dorsalis",], length)
# only consider temperature*family combinations with >= 5 individuals
fam_dors_p[fam_dors_p$Tamano_pupa_mm>=5,]
# for each family, count how many temperatures have >= 5 individuals
(fam_dors_p[fam_dors_p$Tamano_pupa_mm>=5,]%>%
    count(Pareja))->dors_p_fam_list
# only consider families that have individuals in at least 2 treatments 
dors_p_fam_list<-dors_p_fam_list[dors_p_fam_list$n>=2,]

# adult lengths (and masses)
# count number of individuals per family (Pareja) and temperature with adult
# lengths
fam_dors_a<-aggregate(Tamano_adulto_mm~Tratamiento_temperatura_C+Pareja,
                      data=larval_set1[larval_set1$Especie=="C_dorsalis",], length)
# only consider temperature*family combinations with >= 5 individuals
fam_dors_a[fam_dors_a$Tamano_adulto_mm>=5,]
# for each family, count how many temperatures have >= 5 individuals
(fam_dors_a[fam_dors_a$Tamano_adulto_mm>=5,]%>%
    count(Pareja))->dors_a_fam_list
# only consider families that have individuals in at least 2 treatments 
dors_a_fam_list<-dors_a_fam_list[dors_a_fam_list$n>=2,]
dors_a_fam_list
# six families


# placida families for quantitative genetics ------------------------------
# Day 1 lengths
# count number of individuals per family (Pareja) and temperature with Day 1
# lengths
fam_plac_1<-aggregate(Tamano_dia_1_mm~Tratamiento_temperatura_C+Pareja,
                      data=larval_set1[larval_set1$Especie=="C_placida",], length)
# only consider temperature*family combinations with >= 5 individuals
fam_plac_1[fam_plac_1$Tamano_dia_1_mm>=5,]
# for each family, count how many temperatures have >= 5 individuals
(fam_plac_1[fam_plac_1$Tamano_dia_1_mm>=5,]%>%
    count(Pareja))->plac_T1_fam_list
# only consider families that have individuals in at least 2 treatments 
plac_T1_fam_list<-plac_T1_fam_list[plac_T1_fam_list$n>=2,]

#Day 16 lengths
# count number of individuals per family (Pareja) and temperature with Day 16
# lengths
fam_plac_15<-aggregate(Tamano_dia_15_mm~Tratamiento_temperatura_C+Pareja,
                       data=larval_set1[larval_set1$Especie=="C_placida",], length)
# only consider temperature*family combinations with >= 5 individuals
fam_plac_15[fam_plac_15$Tamano_dia_15_mm>=5,]
# for each family, count how many temperatures have >= 5 individuals
(fam_plac_15[fam_plac_15$Tamano_dia_15_mm>=5,]%>%
    count(Pareja))->plac_T15_fam_list
# only consider families that have individuals in at least 2 treatments 
plac_T15_fam_list<-plac_T15_fam_list[plac_T15_fam_list$n>=2,]

# pupal lengths (and masses)
# count number of individuals per family (Pareja) and temperature with pupal
# lengths
fam_plac_p<-aggregate(Tamano_pupa_mm~Tratamiento_temperatura_C+Pareja,
                      data=larval_set1[larval_set1$Especie=="C_placida",], length)
# only consider temperature*family combinations with >= 5 individuals
fam_plac_p[fam_plac_p$Tamano_pupa_mm>=5,]
# for each family, count how many temperatures have >= 5 individuals
(fam_plac_p[fam_plac_p$Tamano_pupa_mm>=5,]%>%
    count(Pareja))->plac_p_fam_list
# only consider families that have individuals in at least 2 treatments 
plac_p_fam_list<-plac_p_fam_list[plac_p_fam_list$n>=2,]

#adult lengths (and masses)
# count number of individuals per family (Pareja) and temperature with adult
# lengths
fam_plac_a<-aggregate(Tamano_adulto_mm~Tratamiento_temperatura_C+Pareja,
                      data=larval_set1[larval_set1$Especie=="C_placida",], length)
# only consider temperature*family combinations with >= 5 individuals
fam_plac_a[fam_plac_a$Tamano_adulto_mm>=5,]
# for each family, count how many temperatures have >= 5 individuals
(fam_plac_a[fam_plac_a$Tamano_adulto_mm>=5,]%>%
    count(Pareja))->plac_a_fam_list
# only consider families that have individuals in at least 2 treatments 
plac_a_fam_list<-plac_a_fam_list[plac_a_fam_list$n>=2,]
plac_a_fam_list


# Table and Fig. 4: Larval survivorship ------------------------------------
# Goal: Testing the effects of temperature and family on survivorship.

# Because the coding of pupate is 0=dead larva, 1=pupated, we can use
# 1-pupate to represent survivorship. 

# We need to use Cox proportional hazards here so that we can look at
# temperature x family interactions using a mixed effects model. This can be
# done using either the frailty function in coxph or using a new package, coxme.
# coxme is newer, but it doesn't have the ability to test for an interaction
# effect. We'll use coxph with frailty.

# The recommended method to compare multiple effects in Cox proportional hazards
# models built with coxph is to perform log-likelihood ratio (LLR) tests rather
# than just use the p-values from the summary (Therneau 2017).So we'll start by
# comparing the full model to the additive model, to test whether the
# temperature x family interaction is significant, then continue with the 
# main effects.

# A note on naming models and LLR tests: I am naming models based on which
# effects they include (add=additive, f=family-only, etc.). I am naming LLR
# tests based on the effect they test the significance of. So an LLR test
# comparing an additive model and a family-only model will be named LLR_t,
# because it tests the significance of the temperature effect.


# dilaticollis survivorship -----------------------------------------------
# The other question is whether we need to analyze a subset of the data or not.
# Because survivorship only requires that a family have produced >= 5 individuals
# per temperature, the list of suitable families is dila_T1_fam_list
dila_T1_fam_list

# Note: There is also 1 larva that hatched at 35C and died the same day. It
# needs to be omitted  because it's a treatment sampling problem.

# limit analysis to families with at least two N>=5 treatments at day 1
dila_lim<-dila[as.numeric(as.character(dila$Pareja)) %in%
                 as.numeric(as.character(dila_T1_fam_list$Pareja)) 
               # & as.numeric(as.character(dila_lim$Tratamiento_temperatura_C))==20 |
               #as.numeric(as.character(dila_lim$Tratamiento_temperatura_C))==25
               ,]

# # we need to get rid of empty levels in the factor Pareja
# dila_lim$Pareja<- factor(dila_lim$Pareja, exclude=c(1, 21, 26, 28, 34, 35))

# we also need to get rid of an empty level in the factor Tratamiento_temperatura_C
dila_lim$Tratamiento_temperatura_C<- factor(dila_lim$Tratamiento_temperatura_C, 
                                            exclude=35)

# full model for these data
cph_dilaLFULL<-coxph(Surv(ldays, event=(1-pupate), type='right') ~ 
                      Tratamiento_temperatura_C*frailty(Pareja),
                    data=dila_lim)
summary(cph_dilaLFULL)
logLik(cph_dilaLFULL)

# is the temeperature x family interaction significant?
# additive mixed model
cph_dilaLadd<-coxph(Surv(ldays, event=(1-pupate), type='right') ~ 
                    Tratamiento_temperatura_C+frailty(Pareja),
                  data=dila_lim)
summary(cph_dilaLadd)
logLik(cph_dilaLadd)

# compare the additive and full models
logLik(cph_dilaLFULL)
LLRLi<-(logLik(cph_dilaLadd)-logLik(cph_dilaLFULL))*-2
LLRLi
pchisq(LLRLi[1],9.930744,lower.tail = FALSE)
# No. 

# is family significant?
# Temperature-only model
cph_dilaLt<-coxph(Surv(ldays, event=(1-pupate), type='right') ~ 
                    Tratamiento_temperatura_C,
                  data=dila_lim)
summary(cph_dilaLt)
logLik(cph_dilaLt)
# compare additive and temperature-only models
LLRLf<-(logLik(cph_dilaLt)-logLik(cph_dilaLadd))*-2
LLRLf
pchisq(LLRLf[1],3,lower.tail = FALSE) 
# YES.

# is temperature significant?
# family-only model
cph_dilaLf<-coxph(Surv(ldays, event=(1-pupate), type='right') ~ 
                    frailty(Pareja),
                  data=dila_lim)
summary(cph_dilaLf)
logLik(cph_dilaLf)
# compare additive and family-only models
LLRLt<-(logLik(cph_dilaLf)-logLik(cph_dilaLadd))*-2
LLRLt
# df from summary read-out (df=9.81)
pchisq(LLRLt[1],9.81,lower.tail = FALSE) 
# Yes.

# The the best model of these data is the additive mixed model.

survdiff(Surv(ldays, event=(1-pupate), type='right') ~ Tratamiento_temperatura_C,
         data=dila_lim)

# pairwise comparisons
survdiff(Surv(ldays, event=(1-pupate), type='right') ~ Tratamiento_temperatura_C,
         data=dila_lim[dila_lim$Tratamiento_temperatura_C==15 |
                     dila_lim$Tratamiento_temperatura_C==30,])
# 15 and 30 are different

survdiff(Surv(ldays, event=(1-pupate), type='right') ~ Tratamiento_temperatura_C,
         data=dila_lim[dila_lim$Tratamiento_temperatura_C==15 |
                     dila_lim$Tratamiento_temperatura_C==25,])
# 15 and 25 are the same

survdiff(Surv(ldays, event=(1-pupate), type='right') ~ Tratamiento_temperatura_C,
         data=dila_lim[dila_lim$Tratamiento_temperatura_C==15 |
                     dila_lim$Tratamiento_temperatura_C==20,])
# 15 and 20 are different

survdiff(Surv(ldays, event=(1-pupate), type='right') ~ Tratamiento_temperatura_C,
         data=dila_lim[dila_lim$Tratamiento_temperatura_C==20 |
                     dila_lim$Tratamiento_temperatura_C==25,])
# 20 and 25 are different

survdiff(Surv(ldays, event=(1-pupate), type='right') ~ Tratamiento_temperatura_C,
         data=dila_lim[dila_lim$Tratamiento_temperatura_C==20 |
                         dila_lim$Tratamiento_temperatura_C==30,])
# 20 and 30 are different

survdiff(Surv(ldays, event=(1-pupate), type='right') ~ Tratamiento_temperatura_C,
         data=dila_lim[dila_lim$Tratamiento_temperatura_C==30 |
                         dila_lim$Tratamiento_temperatura_C==25,])
# 25 and 30 are different

# dilaticollis survival plot (Fig. 4A)----------------------------------------
# Instead of graphing this using the survival package, we'll calculate the
# proportion alive as larvae each day and plot it in ggplot. 

# limited more-balanced data set
# make a matrix for dilaticollis that includes temperature (15, 20 ,25, 30), day, and the number and proportion of
# live individuals for each temperature
# make matrix and fill in temperatures
l_dila_day_by_day_lim<-matrix(nrow=4*(max(dila_lim$ldays, na.rm=TRUE)+1),ncol=4)
colnames(l_dila_day_by_day_lim)<-c("Temp", "Day", "Live_inds", "Prop_alive")
l_dila_day_by_day_lim<-as.data.frame(l_dila_day_by_day_lim)
l_dila_day_by_day_lim[1:dim(l_dila_day_by_day_lim)[1]/4,1]<-15
l_dila_day_by_day_lim[((dim(l_dila_day_by_day_lim)[1]/4)+1):(2*(dim(l_dila_day_by_day_lim)[1]/4)),1]<-20
l_dila_day_by_day_lim[(2*(dim(l_dila_day_by_day_lim)[1]/4)+1):(3*(dim(l_dila_day_by_day_lim)[1]/4)),1]<-25
l_dila_day_by_day_lim[(3*(dim(l_dila_day_by_day_lim)[1]/4)+1):(4*(dim(l_dila_day_by_day_lim)[1]/4)),1]<-30

# fill in the day numbers
l_dila_day_by_day_lim[,2]<-c(0:max(dila_lim$ldays, na.rm=TRUE))


# for each day, calculate the number of live individuals and the proportion 
# alive
for (i in 1:dim(l_dila_day_by_day_lim)[1]){
  #i
  #print(l_dila_day_by_day_lim$Temp[i])# the number of live larvae is the sum of those that survive to pupate (live) 
  # and the ones that die after day i (notdeadyet)
  live<-sum(dila_lim$pupate[dila_lim$Tratamiento_temperatura_C==l_dila_day_by_day_lim$Temp[i]], na.rm=TRUE)
  #print(live)
  notdeadyet<-dim(dila_lim[dila_lim$Tratamiento_temperatura_C==l_dila_day_by_day_lim$Temp[i] &
                             dila_lim$pupate==0 &
                             dila_lim$ldays>l_dila_day_by_day_lim$Day[i],])[1]
  # number of live larvae
  l_dila_day_by_day_lim[i,3]<-live+notdeadyet
  # number of live larvae/total number of larvae in treatment
  l_dila_day_by_day_lim[i,4]<-l_dila_day_by_day_lim[i,3]/dim(dila_lim[as.numeric(as.character(dila_lim$Tratamiento_temperatura_C))==
                                                                as.numeric(as.character(l_dila_day_by_day_lim$Temp[i])),])[1]
}

# we want the lines to truncate after they reach their minimum
for (i in 1:4){
  t<-15+(i-1)*5
  cut<-min(which(l_dila_day_by_day_lim$Prop_alive[l_dila_day_by_day_lim$Temp==t]==
                   min(l_dila_day_by_day_lim$Prop_alive[l_dila_day_by_day_lim$Temp==t], na.rm=TRUE)))
  if(cut<max(l_dila_day_by_day_lim$Day, na.rm=TRUE)){
    l_dila_day_by_day_lim[l_dila_day_by_day_lim$Temp==t & l_dila_day_by_day_lim$Day>cut,3:4]<-NA
  }
}

l_dila_day_by_day_lim$Temp<-as.factor(l_dila_day_by_day_lim$Temp)

# plot balanced data set
p_l_dila_surv_lim<- ggplot(data=l_dila_day_by_day_lim, aes(x=Day, y=Prop_alive, group=Temp)) + geom_line(aes(color=Temp),size=1.2,alpha=0.6) + 
  scale_color_manual(breaks=c(15,20,25,30),labels=c("15","20","25","30"),
                     values=c("turquoise", "yellow","orange","red"))+
  expand_limits(x=c(0,240))+
  scale_x_continuous(breaks=pretty_breaks(2)) +
  scale_y_continuous(breaks=pretty_breaks(2)) +
  labs(x="Days", y="Proportion of larvae alive", title="balanced dilaticollis larval survival", family="sans", color="Temperature") +
  theme_classic2() 

p_l_dila_surv_lim 

# dorsalis survivorship ---------------------------------------------------
# limited to families with at least two N>=5 treatments at day 1
dors_lim<-dors[as.numeric(as.character(dors$Pareja)) %in%
                 as.numeric(as.character(dors_T1_fam_list$Pareja)) 
               # & as.numeric(as.character(dors_lim$Tratamiento_temperatura_C))==20 |
               #as.numeric(as.character(dors_lim$Tratamiento_temperatura_C))==25
               ,]

# # we need to get rid of empty levels in the factor Pareja
# dors_lim$Pareja<- factor(dors_lim$Pareja, exclude=c(1, 3,8:11,15,19, 26, 28, 33:34))

# full model for these data
cph_dorsFULL<-coxph(Surv(ldays, event=(1-pupate), type='right') ~ 
                      Tratamiento_temperatura_C*frailty(Pareja),
                    data=dors_lim)
summary(cph_dorsFULL)

# additive model
cph_dorsLadd<-coxph(Surv(ldays, event=(1-pupate), type='right') ~ 
                    Tratamiento_temperatura_C+frailty(Pareja),
                  data=dors_lim)
summary(cph_dorsLadd)
logLik(cph_dorsLadd)

# is temperature x family significant?
logLik(cph_dorsFULL)
LLRLi<-(logLik(cph_dorsLadd)-logLik(cph_dorsFULL))*-2
LLRLi
pchisq(LLRLi[1],5.881594,lower.tail = FALSE)
#No.

# is family significant?
# temperature-only model
cph_dorsLt<-coxph(Surv(ldays, event=(1-pupate), type='right') ~ 
                    Tratamiento_temperatura_C,
                  data=dors_lim)
summary(cph_dorsLt)
logLik(cph_dorsLt)
# compare to additive model
LLRLf<-(logLik(cph_dorsLt)-logLik(cph_dorsLadd))*-2
LLRLf
pchisq(LLRLf[1],3,lower.tail = FALSE)
# YES.

# is temperature significant?
# family-only model
cph_dorsLf<-coxph(Surv(ldays, event=(1-pupate), type='right') ~ 
                    frailty(Pareja),
                  data=dors_lim)
summary(cph_dorsLf)
logLik(cph_dorsLf)
#compare to additive model
LLRLt<-(logLik(cph_dorsLf)-logLik(cph_dorsLadd))*-2
LLRLt
# df from summary read-out (13.33)
pchisq(LLRLt[1],13.33,lower.tail = FALSE)
# YES.

# The additive mixed model is the best for these data.

# which temperature groups are different?
survdiff(Surv(ldays, event=(1-pupate), type='right') ~ Tratamiento_temperatura_C,
         data=dors_lim)

# pairwise comparisons
survdiff(Surv(ldays, event=(1-pupate), type='right') ~ Tratamiento_temperatura_C,
         data=dors_lim[dors_lim$Tratamiento_temperatura_C==30 |
                     dors_lim$Tratamiento_temperatura_C==15,])
# 30 and 15 are different

survdiff(Surv(ldays, event=(1-pupate), type='right') ~ Tratamiento_temperatura_C,
         data=dors_lim[dors_lim$Tratamiento_temperatura_C==15 |
                     dors_lim$Tratamiento_temperatura_C==25,])
# 15 and 25 are different

survdiff(Surv(ldays, event=(1-pupate), type='right') ~ Tratamiento_temperatura_C,
         data=dors_lim[dors_lim$Tratamiento_temperatura_C==20 |
                     dors_lim$Tratamiento_temperatura_C==25,])
# 20 and 25 are different

survdiff(Surv(ldays, event=(1-pupate), type='right') ~ Tratamiento_temperatura_C,
         data=dors_lim[dors_lim$Tratamiento_temperatura_C==20 |
                         dors_lim$Tratamiento_temperatura_C==30,])
# 20 and 30 are different

survdiff(Surv(ldays, event=(1-pupate), type='right') ~ Tratamiento_temperatura_C,
         data=dors_lim[dors_lim$Tratamiento_temperatura_C==25 |
                         dors_lim$Tratamiento_temperatura_C==30,])
# 25 and 30 are different


# dorsalis survivorship plot (Fig. 4B)---------------------------------------
# Instead of graphing this using the survival package, we'll calculate the
# proportion alive as larvae each day and plot it in ggplot. 
# Make a matrix for dorsalis that includes temperature (15, 20 ,25, 30), day,
# and the number and proportion of live individuals for each temperature.

# limited dataset
l_dors_day_by_day_lim<-matrix(nrow=4*(max(dors_lim$ldays, na.rm=TRUE)+1),ncol=4)
colnames(l_dors_day_by_day_lim)<-c("Temp", "Day", "Live_inds", "Prop_alive")
l_dors_day_by_day_lim<-as.data.frame(l_dors_day_by_day_lim)
l_dors_day_by_day_lim[1:dim(l_dors_day_by_day_lim)[1]/4,1]<-15
l_dors_day_by_day_lim[((dim(l_dors_day_by_day_lim)[1]/4)+1):(2*(dim(l_dors_day_by_day_lim)[1]/4)),1]<-20
l_dors_day_by_day_lim[(2*(dim(l_dors_day_by_day_lim)[1]/4)+1):(3*(dim(l_dors_day_by_day_lim)[1]/4)),1]<-25
l_dors_day_by_day_lim[(3*(dim(l_dors_day_by_day_lim)[1]/4)+1):(4*(dim(l_dors_day_by_day_lim)[1]/4)),1]<-30


l_dors_day_by_day_lim[,2]<-c(0:max(dors_lim$ldays, na.rm=TRUE))


for (i in 1:dim(l_dors_day_by_day_lim)[1]){
  #i
  #print(l_dors_day_by_day_lim$Temp[i])
  live<-sum(dors_lim$pupate[dors_lim$Tratamiento_temperatura_C==l_dors_day_by_day_lim$Temp[i]], na.rm=TRUE)
  #print(live)
  notdeadyet<-dim(dors_lim[dors_lim$Tratamiento_temperatura_C==l_dors_day_by_day_lim$Temp[i] &
                             dors_lim$pupate==0 &
                             dors_lim$ldays>l_dors_day_by_day_lim$Day[i],])[1]
  l_dors_day_by_day_lim[i,3]<-live+notdeadyet
  l_dors_day_by_day_lim[i,4]<-l_dors_day_by_day_lim[i,3]/dim(dors_lim[as.numeric(as.character(dors_lim$Tratamiento_temperatura_C))==
                                                                as.numeric(as.character(l_dors_day_by_day_lim$Temp[i])),])[1]
}

# we want the lines to truncate after they reach their minimum
for (i in 1:4){
  t<-15+(i-1)*5
  cut<-min(which(l_dors_day_by_day_lim$Prop_alive[l_dors_day_by_day_lim$Temp==t]==
                   min(l_dors_day_by_day_lim$Prop_alive[l_dors_day_by_day_lim$Temp==t], na.rm=TRUE)))
  if(cut<max(l_dors_day_by_day_lim$Day, na.rm=TRUE)){
    l_dors_day_by_day_lim[l_dors_day_by_day_lim$Temp==t & l_dors_day_by_day_lim$Day>cut,3:4]<-NA
  }
}

l_dors_day_by_day_lim$Temp<-as.factor(l_dors_day_by_day_lim$Temp)

p_l_dors_surv_lim<- ggplot(data=l_dors_day_by_day_lim, aes(x=Day, y=Prop_alive, group=Temp))+ 
  geom_line(aes(color=Temp),size=1.2,alpha=0.6) + 
  scale_color_manual(breaks=c(15,20,25,30),labels=c("15","20","25","30"),
                     values=c("turquoise", "yellow","orange","red"))+
  expand_limits(x=c(0,240))+
  scale_x_continuous(breaks=pretty_breaks(2)) +
  scale_y_continuous(breaks=pretty_breaks(2)) +
  labs(x="Days", y="Proportion of larvae alive", title="balanced dorsalis larval survival", family="sans", color="Temperature") +
  theme_classic2() 

p_l_dors_surv_lim 


# placida survivorship ----------------------------------------------------
# limited to families with at least two N>=5 treatments at day 1
# there are two larvae that hatched at 35C and died the same day. We should
# exclude these and get rid of 35C as a level in the factor Tratamiento_temperatura_C
plac_lim<-plac_lim[as.numeric(as.character(plac_lim$Tratamiento_temperatura_C))!=35,]
plac_lim$Tratamiento_temperatura_C<-factor(plac_lim$Tratamiento_temperatura_C, exclude=35)

# full model for these data
cph_placFULL<-coxph(Surv(ldays, event=(1-pupate), type='right') ~ 
                      Tratamiento_temperatura_C*frailty(Pareja),
                    data=plac_lim)
summary(cph_placFULL)

# additive model
cph_placLadd<-coxph(Surv(ldays, event=(1-pupate), type='right') ~ 
                    Tratamiento_temperatura_C+frailty(Pareja),
                  data=plac_lim)
summary(cph_placLadd)
logLik(cph_placLadd)

# is temperature x family significant?
logLik(cph_placFULL)
LLRLi<-(logLik(cph_placLadd)-logLik(cph_placFULL))*-2
LLRLi
pchisq(LLRLi[1],0.986526,lower.tail = FALSE)
# NO

# is family significant?
# temperature-only model
cph_placLt<-coxph(Surv(ldays, event=(1-pupate), type='right') ~ 
                    Tratamiento_temperatura_C,
                  data=plac_lim)
summary(cph_placLt)
logLik(cph_placLt)
# compare temperature-only and additive models
LLRLf<-(logLik(cph_placLt)-logLik(cph_placLadd))*-2
LLRLf
pchisq(LLRLf[1],3,lower.tail = FALSE)
# NO


# is temperature significant?
# family-only model
cph_placLf<-coxph(Surv(ldays, event=(1-pupate), type='right') ~ 
                    frailty(Pareja),
                  data=plac_lim)
summary(cph_placLf)
logLik(cph_placLf)
# compare family-only and additive models
LLRLt<-(logLik(cph_placLf)-logLik(cph_placLadd))*-2
LLRLt
# df from summary (9.1)
pchisq(LLRLt[1],9.1,lower.tail = FALSE)
# YES


# The temperature-only model is the best for these data.
# which temperatures are different?
survdiff(Surv(ldays, event=(1-pupate), type='right') ~ Tratamiento_temperatura_C,
         data=plac_lim[as.numeric(as.character(plac_lim$Tratamiento_temperatura_C))!=35,])

# pairwise comparisons
survdiff(Surv(ldays, event=(1-pupate), type='right') ~ Tratamiento_temperatura_C,
         data=plac_lim[as.numeric(as.character(plac_lim$Tratamiento_temperatura_C))==30 |
                     as.numeric(as.character(plac_lim$Tratamiento_temperatura_C))==15,])
# 30 and 15 are different

survdiff(Surv(ldays, event=(1-pupate), type='right') ~ Tratamiento_temperatura_C,
         data=plac_lim[as.numeric(as.character(plac_lim$Tratamiento_temperatura_C))==25 |
                     as.numeric(as.character(plac_lim$Tratamiento_temperatura_C))==15,])
# 15 and 25 are different

survdiff(Surv(ldays, event=(1-pupate), type='right') ~ Tratamiento_temperatura_C,
         data=plac_lim[as.numeric(as.character(plac_lim$Tratamiento_temperatura_C))==15 |
                         as.numeric(as.character(plac_lim$Tratamiento_temperatura_C))==20,])
#  15 and 20 are different

survdiff(Surv(ldays, event=(1-pupate), type='right') ~ Tratamiento_temperatura_C,
         data=plac_lim[as.numeric(as.character(plac_lim$Tratamiento_temperatura_C))==25 |
                     as.numeric(as.character(plac_lim$Tratamiento_temperatura_C))==20,])
# 20 and 25 are the same

survdiff(Surv(ldays, event=(1-pupate), type='right') ~ Tratamiento_temperatura_C,
         data=plac_lim[as.numeric(as.character(plac_lim$Tratamiento_temperatura_C))==30 |
                         as.numeric(as.character(plac_lim$Tratamiento_temperatura_C))==20,])
# 20 and 30 are different

survdiff(Surv(ldays, event=(1-pupate), type='right') ~ Tratamiento_temperatura_C,
         data=plac_lim[as.numeric(as.character(plac_lim$Tratamiento_temperatura_C))==30 |
                         as.numeric(as.character(plac_lim$Tratamiento_temperatura_C))==25,])
# 25 and 30 are different

# placida survivorship plot (Fig. 4C)---------------------------------------
# Instead of graphing this using the survival package, we'll calculate the
# proportion alive as larvae each day and plot it in ggplot. 
# Make a matrix for placalis that includes temperature (15, 20 ,25, 30), day,
# and the number and proportion of live individuals for each temperature.

# limited dataset
# make matrix and fill in temperatures
l_plac_day_by_day_lim<-matrix(nrow=4*(max(plac_lim$ldays, na.rm=TRUE)+1),ncol=4)
colnames(l_plac_day_by_day_lim)<-c("Temp", "Day", "Live_inds", "Prop_alive")
l_plac_day_by_day_lim<-as.data.frame(l_plac_day_by_day_lim)
l_plac_day_by_day_lim[1:dim(l_plac_day_by_day_lim)[1]/4,1]<-15
l_plac_day_by_day_lim[((dim(l_plac_day_by_day_lim)[1]/4)+1):(2*(dim(l_plac_day_by_day_lim)[1]/4)),1]<-20
l_plac_day_by_day_lim[(2*(dim(l_plac_day_by_day_lim)[1]/4)+1):(3*(dim(l_plac_day_by_day_lim)[1]/4)),1]<-25
l_plac_day_by_day_lim[(3*(dim(l_plac_day_by_day_lim)[1]/4)+1):(4*(dim(l_plac_day_by_day_lim)[1]/4)),1]<-30

# fill in the day numbers
l_plac_day_by_day_lim[,2]<-c(0:max(plac_lim$ldays, na.rm=TRUE))

# for each day, calculate the number of live individuals and the proportion 
# alive
for (i in 1:dim(l_plac_day_by_day_lim)[1]){
  #i
  #print(l_plac_day_by_day_lim$Temp[i])
  # the number of live larvae is the sum of those that survive to pupate (live) 
  # and the ones that die after day i (notdeadyet)
  live<-sum(plac_lim$pupate[plac_lim$Tratamiento_temperatura_C==l_plac_day_by_day_lim$Temp[i]], na.rm=TRUE)
  #print(live)
  notdeadyet<-dim(plac_lim[plac_lim$Tratamiento_temperatura_C==l_plac_day_by_day_lim$Temp[i] &
                             plac_lim$pupate==0 &
                             plac_lim$ldays>l_plac_day_by_day_lim$Day[i],])[1]
  # number of live larvae
  l_plac_day_by_day_lim[i,3]<-live+notdeadyet
  # number of live larvae/total number of larvae in treatment
  l_plac_day_by_day_lim[i,4]<-l_plac_day_by_day_lim[i,3]/dim(plac_lim[as.numeric(as.character(plac_lim$Tratamiento_temperatura_C))==
                                                                        as.numeric(as.character(l_plac_day_by_day_lim$Temp[i])),])[1]
}

# we want the lines to truncate after they reach their minimum
for (i in 1:4){
  t<-15+(i-1)*5
  cut<-min(which(l_plac_day_by_day_lim$Prop_alive[l_plac_day_by_day_lim$Temp==t]==
                   min(l_plac_day_by_day_lim$Prop_alive[l_plac_day_by_day_lim$Temp==t], na.rm=TRUE)))
  if(cut<max(l_plac_day_by_day_lim$Day, na.rm=TRUE)){
    l_plac_day_by_day_lim[l_plac_day_by_day_lim$Temp==t & l_plac_day_by_day_lim$Day>cut,3:4]<-NA
  }
}

l_plac_day_by_day_lim$Temp<-as.factor(l_plac_day_by_day_lim$Temp)

p_l_plac_surv_lim<- ggplot(data=l_plac_day_by_day_lim, aes(x=Day, y=Prop_alive, group=Temp))+ 
  geom_line(aes(color=Temp),size=1.2,alpha=0.6) + 
  scale_color_manual(breaks=c(15,20,25,30),labels=c("15","20","25","30"),
                     values=c("turquoise", "yellow","orange","red"))+
  expand_limits(x=c(0,240))+
  scale_x_continuous(breaks=pretty_breaks(2)) +
  scale_y_continuous(breaks=pretty_breaks(2)) +
  labs(x="Days", y="Proportion of larvae alive", title="balanced placida larval survival", family="sans", color="Temperature") +
  theme_classic2() 

p_l_plac_surv_lim 

# save Figure 4 larval survivorship plots ----------------------------------------------
# ggsave("Fig4_dila_lim.pdf",plot=p_l_dila_surv_lim,device="pdf")
# ggsave("Fig4_dors_lim.pdf",plot=p_l_dors_surv_lim,device="pdf")
# ggsave("Fig4_plac_lim.pdf",plot=p_l_plac_surv_lim,device="pdf")

# General genotype x environment trait data preparation function ---------
# This is a function that can prepare data from any dataset, species, and trait
# for genotype x environment analyses (mixed effect models, reaction norms, or
# correlations). It will only include families that have >= 5 individuals in
# both 20C and 25C. The output will be a dataframe with the mean trait value and
# sample size for each family x temperature combination.
FT_prep<- function (data1, pop, trait) {
  print(pop)
  print(trait)
  print(head(data1[data1$Especie==pop,colnames(data1)==trait]))
  
  # calculate the mean trait value for each family x temperature combination
  trait_data<-aggregate(data1[data1$Especie==pop,colnames(data1)==trait] ~ 
                          data1$Pareja[data1$Especie==pop] + 
                          data1$Tratamiento_temperatura_C[data1$Especie==pop],
                        FUN=mean,
                        na.remove=TRUE)
  colnames(trait_data)<-c("Pareja", "Temp", trait)
  trait_data$Pareja<-as.numeric(trait_data$Pareja)
  trait_data$Temp<-as.numeric(as.character(trait_data$Temp))
  print(trait_data)
  
  print(data1$Pareja[data1$Especie==pop & is.na(data1[,colnames(data1)==trait])==FALSE])
  # calculate the sample size for each family x temperature combination
  N_trait<-table(data1$Pareja[data1$Especie==pop & is.na(data1[,colnames(data1)==trait])==FALSE],
                 data1$Tratamiento_temperatura_C[data1$Especie== pop & 
                                                   is.na(data1[,colnames(data1)==trait])==FALSE])
  
  N_trait<-as.data.frame(N_trait)
  colnames(N_trait)<-c("Pareja","Temp","N")
  
  N_trait$Pareja<-as.numeric(N_trait$Pareja)
  N_trait$Temp<-as.numeric(as.character(N_trait$Temp))
  N_trait$N<-as.numeric(N_trait$N)
  print(N_trait)
  
  #print(dim(trait_data))
  #print(dim(N_trait))
  
  # put the means and trait values together in one matrix, and only include
  # families that have >= 5 individuals in both 20C and 25C
  trait_data2<-matrix(ncol=dim(trait_data)[2])
  trait_data2<-as.data.frame(trait_data2)
  colnames(trait_data2)<-colnames(trait_data)
  #trait_data2[,3]<-as.numeric(trait_data2)
  #print(trait_data2)
  for (i in 1:max(trait_data$Pareja)) {
    print(c("i =",i, N_trait))
    print(N_trait$N[N_trait$Pareja==i & N_trait$Temp==20])
    print(N_trait$N[N_trait$Pareja==i & N_trait$Temp==25])
    print(N_trait$N[N_trait$Pareja==i & N_trait$Temp==20]>=5 &
            N_trait$N[N_trait$Pareja==i & N_trait$Temp==25]>=5)
    #print(N_trait$N[N_trait$Pareja==i & N_trait$Temp==15]>=5)
    if (i > max(N_trait$Pareja[N_trait$Temp==20]) | i > 
        max(N_trait$Pareja[N_trait$Temp==25])) {}
    else if (any(trait_data$Pareja==i)==FALSE) {}
    else if (#N_trait$N[N_trait$Pareja==i & N_trait$Temp==15]>=5 &
      N_trait$N[N_trait$Pareja==i & N_trait$Temp==20]>=5 &
      N_trait$N[N_trait$Pareja==i & N_trait$Temp==25]>=5) {
      
      if (is.na(trait_data2[1,1])) {
        print(c("is numeric?",is.numeric(trait_data$Temp)))
        #print(trait_data[trait_data$Pareja==i & trait_data$Temp<30,])
        #print(dim(trait_data2))
        stub<-trait_data[trait_data$Pareja==i & 
                           trait_data$Temp<30,]
        print(stub)
        trait_data2[1:dim(stub)[1],1:3]<-stub
        print("start")
        #print(head(trait_data2))
      }
      else {
        #print(colnames(trait_data2))
        #print(colnames(trait_data))
        trait_data2<-rbind(trait_data2,trait_data[trait_data$Pareja==i & trait_data$Temp<30,])
        print("continue")}
      
    }
  }
  
  print(head(trait_data2))
  print(as.numeric(as.character(trait_data2$Temp[1:5])))
  
  # if too few families have enough individuals at 15 degrees, omit the 15C data
  if (length(N_trait$N[N_trait$Temp==15 & N_trait$N>=5])<5) {
    trait_data2<-trait_data2[as.numeric(as.character(trait_data2$Temp))>15,]
  } 
  # add Ns into the trait_data2 data frame
  trait_data2[,4]<-NA
  for (ii in 1:dim(trait_data2)[1]){
    print(c("ii=",ii))
    #print(N_trait[N_trait$Pareja==trait_data2$Pareja[ii] &
    #                N_trait$Temp==trait_data2$Temp[ii],3])
    trait_data2[ii,4]<-N_trait[N_trait$Pareja==trait_data2$Pareja[ii] &
                                 N_trait$Temp==trait_data2$Temp[ii],3]
    #print(trait_data2[ii,])
  }
  colnames(trait_data2)<-c(colnames(trait_data),"N")
  trait_data2$Pareja<-as.factor(trait_data2$Pareja)
  trait_data2$Temp<-as.factor(trait_data2$Temp)
  
  return(trait_data2)
  
}



# Table 5: Models for family x temperature effects on traits ---------------
# Goal: Test whether there are significant family x temperature (genotype x
# environment) interactions, as well as significant family or temperature
# effects on developmental traits.

# Dec. 12, 2018: The R Book (2nd ed.) indicates that you can code a grouping
# variable as both a fixed and a random effect (p. 686, 707) if you are
# interested in the estimates for the different groups.

# But because you'll want to compare models with different fixed effects, you
# need to use maximum likelihood instead of REML. From the best model, you would
# then report the varianace distribution across the different random effects,
# the model statistics, and the fixed parameters of interest.

# There are too many missing interactions to do this using lme, so we'll use
# lmer.

# As with the ANOVAs and t-tests, lengths and masses will be log-transformed.

# A note on naming models and LLR tests: I am naming models based on which
# effects they include (add=additive, f=family-only, etc.). I am naming LLR
# tests based on the effect they test the significance of. So an LLR test
# comparing an additive model and a family-only model will be named LLR_t,
# because it tests the significance of the temperature effect.

# Note: these models do not include the larval survival mixed models, which 
# are in the larval survivorship section

# dilaticollis mixed effect models-------------------------------------------
# Day 1 lengths
# Dec. 6, 2018: I am filtering the data to include only 20C and 25C individuals 
# from well-represented families. Explore whether this is the best way to
# do this

# Dec. 17, 2018: Go back and check the dila_T1_fam_list!
dila_T1_fam_list

# the fewest balanced families are for all four temperatures (n=6)
# we'll use the complete dila_lim data set for now

dila_fam_T1_full<-lmer(log(Tamano_dia_1_mm)~Tratamiento_temperatura_C*Pareja+(1|Pareja),
                       data=dila_lim,
                       na.action="na.omit",
                       REML=FALSE)
summary(dila_fam_T1_full)
anova(dila_fam_T1_full)

dila_fam_T1_add<-lmer(log(Tamano_dia_1_mm)~Tratamiento_temperatura_C+Pareja+(1|Pareja),
                       data=dila_lim,
                       na.action="na.omit",
                       REML=FALSE)
summary(dila_fam_T1_add)
anova(dila_fam_T1_add)

anova(dila_fam_T1_full,dila_fam_T1_add)
# so the interaction doesn't improve the additive model

# what happens if we remove Pareja as a fixed effect?
dila_fam_T1_Tonly<-lmer(log(Tamano_dia_1_mm)~Tratamiento_temperatura_C+(1|Pareja),
                      data=dila_lim,
                      na.action="na.omit",
                      REML=FALSE)
summary(dila_fam_T1_Tonly)
anova(dila_fam_T1_Tonly)
anova(dila_fam_T1_Tonly, dila_fam_T1_add)
# the additive model is better

# what about a random-only model?
dila_fam_T1_rand<-lmer(log(Tamano_dia_1_mm)~(1|Pareja),
                       data=dila_lim,
                       na.action="na.omit",
                       REML=FALSE)
summary(dila_fam_T1_rand)
anova(dila_fam_T1_rand)
anova(dila_fam_T1_rand, dila_fam_T1_Tonly)

# Dec. 19, 2018: But if Pareja is significant, then the significance of
# temperature should be tested by comparing the additive model to a Pareja-only
# model
dila_fam_T1_Ponly<-lmer(log(Tamano_dia_1_mm)~Pareja+(1|Pareja),
                       data=dila_lim,
                       na.action="na.omit",
                       REML=FALSE)
summary(dila_fam_T1_Ponly)
anova(dila_fam_T1_Ponly)
anova(dila_fam_T1_Ponly, dila_fam_T1_add)
# Temperature's still highly significant

# Day 16 length
dila_T15_fam_list
# there are enough three-way balanced families to include 15C

dila_fam_T15_full<-lmer(log(Tamano_dia_15_mm)~Tratamiento_temperatura_C*Pareja+(1|Pareja),
                       data=dila[as.numeric(as.character(dila$Pareja)) %in% 
                                   as.numeric(as.character(dila_T15_fam_list$Pareja)),],
                       na.action="na.omit",
                       REML=FALSE)
summary(dila_fam_T15_full)
anova(dila_fam_T15_full)

dila_fam_T15_add<-lmer(log(Tamano_dia_15_mm)~Tratamiento_temperatura_C+Pareja+(1|Pareja),
                      data=dila[as.numeric(as.character(dila$Pareja)) %in% 
                                  as.numeric(as.character(dila_T15_fam_list$Pareja)),],
                      na.action="na.omit",
                      REML=FALSE)
summary(dila_fam_T15_add)
anova(dila_fam_T15_add)

anova(dila_fam_T15_full,dila_fam_T15_add)
# so the interaction doesn't improve the additive model

# what happens if we remove Pareja as a fixed effect?
dila_fam_T15_Tonly<-lmer(log(Tamano_dia_15_mm)~Tratamiento_temperatura_C+(1|Pareja),
                        data=dila[as.numeric(as.character(dila$Pareja)) %in% 
                                    as.numeric(as.character(dila_T15_fam_list$Pareja)),],
                        na.action="na.omit",
                        REML=FALSE)
summary(dila_fam_T15_Tonly)
anova(dila_fam_T15_Tonly)
anova(dila_fam_T15_Tonly, dila_fam_T15_add)
# Including Pareja as a fixed effect doesn't improve the model

# what about a random-only model?
dila_fam_T15_rand<-lmer(log(Tamano_dia_15_mm)~(1|Pareja),
                        data=dila[as.numeric(as.character(dila$Pareja)) %in% 
                                    as.numeric(as.character(dila_T15_fam_list$Pareja)),],
                        na.action="na.omit",
                       REML=FALSE)
summary(dila_fam_T15_rand)
anova(dila_fam_T15_rand)
anova(dila_fam_T15_rand, dila_fam_T15_Tonly)
# Temperature's significant

# Pupal lengths
dila_Tp_fam_list
# now we're down to only 20-25
# we still need to limit these families to only the 20 and 25 data

dila_fam_Tp_full<-lmer(log(Tamano_pupa_mm)~Tratamiento_temperatura_C*Pareja+(1|Pareja),
                        data=dila[as.numeric(as.character(dila$Pareja)) %in% 
                                    as.numeric(as.character(dila_Tp_fam_list$Pareja)) &
                                   (as.numeric(as.character(dila$Tratamiento_temperatura_C))==20 |
                                      as.numeric(as.character(dila$Tratamiento_temperatura_C))==25),],
                        na.action="na.omit",
                        REML=FALSE)
summary(dila_fam_Tp_full)
anova(dila_fam_Tp_full)

dila_fam_Tp_add<-lmer(log(Tamano_pupa_mm)~Tratamiento_temperatura_C+Pareja+(1|Pareja),
                      data=dila[as.numeric(as.character(dila$Pareja)) %in% 
                                  as.numeric(as.character(dila_Tp_fam_list$Pareja)) &
                                  (as.numeric(as.character(dila$Tratamiento_temperatura_C))==20 |
                                     as.numeric(as.character(dila$Tratamiento_temperatura_C))==25),],
                      na.action="na.omit",
                       REML=FALSE)
summary(dila_fam_Tp_add)
anova(dila_fam_Tp_add)

anova(dila_fam_Tp_full,dila_fam_Tp_add)
# the interaction doesn't improve the additive model

# what happens if we remove Pareja as a fixed effect?
dila_fam_Tp_Tonly<-lmer(log(Tamano_pupa_mm)~Tratamiento_temperatura_C+(1|Pareja),
                        data=dila[as.numeric(as.character(dila$Pareja)) %in% 
                                    as.numeric(as.character(dila_Tp_fam_list$Pareja)) &
                                    (as.numeric(as.character(dila$Tratamiento_temperatura_C))==20 |
                                       as.numeric(as.character(dila$Tratamiento_temperatura_C))==25),],
                        na.action="na.omit",
                         REML=FALSE)
summary(dila_fam_Tp_Tonly)
anova(dila_fam_Tp_Tonly)
anova(dila_fam_Tp_Tonly, dila_fam_Tp_add)
# family is significant

# what about a random-only model?
dila_fam_Tp_rand<-lmer(log(Tamano_pupa_mm)~(1|Pareja),
                       data=dila[as.numeric(as.character(dila$Pareja)) %in% 
                                   as.numeric(as.character(dila_Tp_fam_list$Pareja)) &
                                   (as.numeric(as.character(dila$Tratamiento_temperatura_C))==20 |
                                      as.numeric(as.character(dila$Tratamiento_temperatura_C))==25),],
                       na.action="na.omit",
                        REML=FALSE)
summary(dila_fam_Tp_rand)
anova(dila_fam_Tp_rand)
anova(dila_fam_Tp_rand, dila_fam_Tp_Tonly)
# Temperature is significant

# Dec. 19, 2018: But if Pareja is significant, then the significance of
# temperature should be tested by comparing the additive model to a Pareja-only
# model
dila_fam_Tp_Ponly<-lmer(log(Tamano_pupa_mm)~Pareja+(1|Pareja),
                       data=dila[as.numeric(as.character(dila$Pareja)) %in% 
                                   as.numeric(as.character(dila_Tp_fam_list$Pareja)) &
                                   (as.numeric(as.character(dila$Tratamiento_temperatura_C))==20 |
                                      as.numeric(as.character(dila$Tratamiento_temperatura_C))==25),],
                       na.action="na.omit",
                       REML=FALSE)
summary(dila_fam_Tp_Ponly)
anova(dila_fam_Tp_Ponly)
anova(dila_fam_Tp_Ponly, dila_fam_Tp_add)
# Temperature is still highly significant

# Adult length
dila_a_fam_list
# there are only four balanced families

# We'll give it a shot:
dila_fam_Ta_full<-lmer(log(Tamano_adulto_mm)~Tratamiento_temperatura_C*Pareja+(1|Pareja),
                       data=dila[as.numeric(as.character(dila$Pareja)) %in% 
                                   as.numeric(as.character(dila_a_fam_list$Pareja)) &
                                   (as.numeric(as.character(dila$Tratamiento_temperatura_C))==20 |
                                      as.numeric(as.character(dila$Tratamiento_temperatura_C))==25),],
                       na.action="na.omit",
                       REML=FALSE)
summary(dila_fam_Ta_full)
anova(dila_fam_Ta_full)

dila_fam_Ta_add<-lmer(log(Tamano_adulto_mm)~Tratamiento_temperatura_C+Pareja+(1|Pareja),
                      data=dila[as.numeric(as.character(dila$Pareja)) %in% 
                                  as.numeric(as.character(dila_a_fam_list$Pareja)) &
                                  (as.numeric(as.character(dila$Tratamiento_temperatura_C))==20 |
                                     as.numeric(as.character(dila$Tratamiento_temperatura_C))==25),],
                      na.action="na.omit",
                      REML=FALSE)
summary(dila_fam_Ta_add)
anova(dila_fam_Ta_add)

anova(dila_fam_Ta_full,dila_fam_Ta_add)
# the interaction doesn't improve the additive model

# is family significant?
dila_fam_Ta_Tonly<-lmer(log(Tamano_adulto_mm)~Tratamiento_temperatura_C+(1|Pareja),
                      data=dila[as.numeric(as.character(dila$Pareja)) %in% 
                                  as.numeric(as.character(dila_a_fam_list$Pareja)) &
                                  (as.numeric(as.character(dila$Tratamiento_temperatura_C))==20 |
                                     as.numeric(as.character(dila$Tratamiento_temperatura_C))==25),],
                      na.action="na.omit",
                      REML=FALSE)
summary(dila_fam_Ta_Tonly)
anova(dila_fam_Ta_Tonly)

anova(dila_fam_Ta_add,dila_fam_Ta_Tonly)
# not quite

# random-only model
dila_fam_Ta_rand<-lmer(log(Tamano_adulto_mm)~(1|Pareja),
                        data=dila[as.numeric(as.character(dila$Pareja)) %in% 
                                    as.numeric(as.character(dila_a_fam_list$Pareja)) &
                                    (as.numeric(as.character(dila$Tratamiento_temperatura_C))==20 |
                                       as.numeric(as.character(dila$Tratamiento_temperatura_C))==25),],
                        na.action="na.omit",
                        REML=FALSE) 
summary(dila_fam_Ta_rand)
anova(dila_fam_Ta_rand)

anova(dila_fam_Ta_rand,dila_fam_Ta_Tonly)
# temperature is significant


# Pupal mass
dila_Pp_fam_list
# we need to limit these families to the 20 and 25C data

dila_fam_Pp_full<-lmer(Peso_pupa_g~Tratamiento_temperatura_C*Pareja+(1|Pareja),
                       data=dila[as.numeric(as.character(dila$Pareja)) %in% 
                                   as.numeric(as.character(dila_Pp_fam_list$Pareja))&
                                   (as.numeric(as.character(dila$Tratamiento_temperatura_C))==20 |
                                      as.numeric(as.character(dila$Tratamiento_temperatura_C))==25),],
                       na.action="na.omit",
                       REML=FALSE)
summary(dila_fam_Pp_full)
anova(dila_fam_Pp_full)

dila_fam_Pp_add<-lmer(Peso_pupa_g~Tratamiento_temperatura_C+Pareja+(1|Pareja),
                      data=dila[as.numeric(as.character(dila$Pareja)) %in% 
                                  as.numeric(as.character(dila_Pp_fam_list$Pareja))&
                                  (as.numeric(as.character(dila$Tratamiento_temperatura_C))==20 |
                                     as.numeric(as.character(dila$Tratamiento_temperatura_C))==25),],
                      na.action="na.omit",
                      REML=FALSE)
summary(dila_fam_Pp_add)
anova(dila_fam_Pp_add)

anova(dila_fam_Pp_full,dila_fam_Pp_add)
# the interaction doesn't improve the additive model

# what happens if we remove Pareja as a fixed effect?
dila_fam_Pp_Tonly<-lmer(Peso_pupa_g~Tratamiento_temperatura_C+(1|Pareja),
                        data=dila[as.numeric(as.character(dila$Pareja)) %in% 
                                    as.numeric(as.character(dila_Pp_fam_list$Pareja))&
                                    (as.numeric(as.character(dila$Tratamiento_temperatura_C))==20 |
                                       as.numeric(as.character(dila$Tratamiento_temperatura_C))==25),],
                        na.action="na.omit",
                        REML=FALSE)
summary(dila_fam_Pp_Tonly)
anova(dila_fam_Pp_Tonly)
anova(dila_fam_Pp_Tonly, dila_fam_Pp_add)
# family is significant

# what about a random-only model?
dila_fam_Pp_rand<-lmer(Peso_pupa_g~(1|Pareja),
                       data=dila[as.numeric(as.character(dila$Pareja)) %in% 
                                   as.numeric(as.character(dila_Pp_fam_list$Pareja))&
                                   (as.numeric(as.character(dila$Tratamiento_temperatura_C))==20 |
                                      as.numeric(as.character(dila$Tratamiento_temperatura_C))==25),],
                       na.action="na.omit",
                       REML=FALSE)
summary(dila_fam_Pp_rand)
anova(dila_fam_Pp_rand)
anova(dila_fam_Pp_rand, dila_fam_Pp_Tonly)
# temperature is significant

# Dec. 19, 2018: But if Pareja is significant, then the significance of
# temperature should be tested by comparing the additive model to a Pareja-only
# model
dila_fam_Pp_Ponly<-lmer(Peso_pupa_g~Pareja+(1|Pareja),
                       data=dila[as.numeric(as.character(dila$Pareja)) %in% 
                                   as.numeric(as.character(dila_Pp_fam_list$Pareja))&
                                   (as.numeric(as.character(dila$Tratamiento_temperatura_C))==20 |
                                      as.numeric(as.character(dila$Tratamiento_temperatura_C))==25),],
                       na.action="na.omit",
                       REML=FALSE)
summary(dila_fam_Pp_Ponly)
anova(dila_fam_Pp_Ponly)
anova(dila_fam_Pp_Ponly, dila_fam_Pp_add)
# Temperature still highly significant

# Adult mass
# again there only four balanced families

# We'll give it a try:
dila_fam_Pa_full<-lmer(Peso_adulto_g~Tratamiento_temperatura_C*Pareja+(1|Pareja),
                       data=dila[as.numeric(as.character(dila$Pareja)) %in% 
                                   as.numeric(as.character(dila_a_fam_list$Pareja))&
                                   (as.numeric(as.character(dila$Tratamiento_temperatura_C))==20 |
                                      as.numeric(as.character(dila$Tratamiento_temperatura_C))==25),],
                       na.action="na.omit",
                       REML=FALSE) 
summary(dila_fam_Pa_full)
anova(dila_fam_Pa_full)

dila_fam_Pa_add<-lmer(Peso_adulto_g~Tratamiento_temperatura_C+Pareja+(1|Pareja),
                      data=dila[as.numeric(as.character(dila$Pareja)) %in% 
                                  as.numeric(as.character(dila_a_fam_list$Pareja))&
                                  (as.numeric(as.character(dila$Tratamiento_temperatura_C))==20 |
                                     as.numeric(as.character(dila$Tratamiento_temperatura_C))==25),],
                      na.action="na.omit",
                      REML=FALSE)
summary(dila_fam_Pa_add)
anova(dila_fam_Pa_add)

anova(dila_fam_Pa_full,dila_fam_Pa_add)
# the interaction doesn't improve the additive model

# is family significant?
dila_fam_Pa_Tonly<-lmer(Peso_adulto_g~Tratamiento_temperatura_C+(1|Pareja),
                      data=dila[as.numeric(as.character(dila$Pareja)) %in% 
                                  as.numeric(as.character(dila_a_fam_list$Pareja))&
                                  (as.numeric(as.character(dila$Tratamiento_temperatura_C))==20 |
                                     as.numeric(as.character(dila$Tratamiento_temperatura_C))==25),],
                      na.action="na.omit",
                      REML=FALSE)
summary(dila_fam_Pa_Tonly)
anova(dila_fam_Pa_Tonly)

anova(dila_fam_Pa_add,dila_fam_Pa_Tonly)
# family is significant

# random-only model?
dila_fam_Pa_rand<-lmer(Peso_adulto_g~(1|Pareja),
                        data=dila[as.numeric(as.character(dila$Pareja)) %in% 
                                    as.numeric(as.character(dila_a_fam_list$Pareja))&
                                    (as.numeric(as.character(dila$Tratamiento_temperatura_C))==20 |
                                       as.numeric(as.character(dila$Tratamiento_temperatura_C))==25),],
                        na.action="na.omit",
                        REML=FALSE)
summary(dila_fam_Pa_rand)
anova(dila_fam_Pa_rand)

anova(dila_fam_Pa_rand,dila_fam_Pa_Tonly)
# Temperature's significant

# since family was significant, also try a Pareja-only model
dila_fam_Pa_Ponly<-lmer(Peso_adulto_g~Pareja+(1|Pareja),
                       data=dila[as.numeric(as.character(dila$Pareja)) %in% 
                                   as.numeric(as.character(dila_a_fam_list$Pareja))&
                                   (as.numeric(as.character(dila$Tratamiento_temperatura_C))==20 |
                                      as.numeric(as.character(dila$Tratamiento_temperatura_C))==25),],
                       na.action="na.omit",
                       REML=FALSE)
summary(dila_fam_Pa_Ponly)
anova(dila_fam_Pa_Ponly)

anova(dila_fam_Pa_Ponly,dila_fam_Pa_add)
# temperature is still significant

# larval development time
# the balanced families for this will be the same as those that had balanced
# pupal measurements
dila_Pp_fam_list
# only consider 20 and 25C

dila_fam_LDT_full<-lmer(ldays~Tratamiento_temperatura_C*Pareja+(1|Pareja),
                       data=dila[dila$pupate==1 &
                         as.numeric(as.character(dila$Pareja)) %in% 
                                   as.numeric(as.character(dila_Pp_fam_list$Pareja))&
                           (as.numeric(as.character(dila$Tratamiento_temperatura_C))==20 |
                              as.numeric(as.character(dila$Tratamiento_temperatura_C))==25),],
                       na.action="na.omit",
                       REML=FALSE)
summary(dila_fam_LDT_full)
anova(dila_fam_LDT_full)

dila_fam_LDT_add<-lmer(ldays~Tratamiento_temperatura_C+Pareja+(1|Pareja),
                       data=dila[dila$pupate==1 &
                                   as.numeric(as.character(dila$Pareja)) %in% 
                                   as.numeric(as.character(dila_Pp_fam_list$Pareja))&
                                   (as.numeric(as.character(dila$Tratamiento_temperatura_C))==20 |
                                      as.numeric(as.character(dila$Tratamiento_temperatura_C))==25),],
                       na.action="na.omit",
                      REML=FALSE)
summary(dila_fam_LDT_add)
anova(dila_fam_LDT_add)

anova(dila_fam_LDT_full,dila_fam_LDT_add)
# the interaction doesn't improve the additive model

# what happens if we remove Pareja as a fixed effect?
dila_fam_LDT_Tonly<-lmer(ldays~Tratamiento_temperatura_C+(1|Pareja),
                         data=dila[dila$pupate==1 &
                                     as.numeric(as.character(dila$Pareja)) %in% 
                                     as.numeric(as.character(dila_Pp_fam_list$Pareja))&
                                     (as.numeric(as.character(dila$Tratamiento_temperatura_C))==20 |
                                        as.numeric(as.character(dila$Tratamiento_temperatura_C))==25),],
                         na.action="na.omit",
                        REML=FALSE)
summary(dila_fam_LDT_Tonly)
anova(dila_fam_LDT_Tonly)
anova(dila_fam_LDT_Tonly, dila_fam_LDT_add)
# family isn't significant

# what about a random-only model?
dila_fam_LDT_rand<-lmer(ldays~(1|Pareja),
                        data=dila[dila$pupate==1 &
                                    as.numeric(as.character(dila$Pareja)) %in% 
                                    as.numeric(as.character(dila_Pp_fam_list$Pareja))&
                                    (as.numeric(as.character(dila$Tratamiento_temperatura_C))==20 |
                                       as.numeric(as.character(dila$Tratamiento_temperatura_C))==25),],
                        na.action="na.omit",
                       REML=FALSE)
summary(dila_fam_LDT_rand)
anova(dila_fam_LDT_rand)
anova(dila_fam_LDT_rand, dila_fam_LDT_Tonly)
# Temperature is significant

# Pupal development time
# again, there are only four balenced families available for this.
# But we'll give it a shot:
# only consider 20 and 25C

dila_fam_PDT_full<-lmer(pdays~Tratamiento_temperatura_C*Pareja+(1|Pareja),
                        data=dila[dila$emerge==1 &
                                    as.numeric(as.character(dila$Pareja)) %in% 
                                    as.numeric(as.character(dila_a_fam_list$Pareja))&
                                    (as.numeric(as.character(dila$Tratamiento_temperatura_C))==20 |
                                       as.numeric(as.character(dila$Tratamiento_temperatura_C))==25),],
                        na.action="na.omit",
                        REML=FALSE)
summary(dila_fam_PDT_full)
anova(dila_fam_PDT_full)

dila_fam_PDT_add<-lmer(pdays~Tratamiento_temperatura_C+Pareja+(1|Pareja),
                       data=dila[dila$emerge==1 &
                                   as.numeric(as.character(dila$Pareja)) %in% 
                                   as.numeric(as.character(dila_a_fam_list$Pareja))&
                                   (as.numeric(as.character(dila$Tratamiento_temperatura_C))==20 |
                                      as.numeric(as.character(dila$Tratamiento_temperatura_C))==25),],
                       na.action="na.omit",
                       REML=FALSE)
summary(dila_fam_PDT_add)
anova(dila_fam_PDT_add)

anova(dila_fam_PDT_full,dila_fam_PDT_add)
# interaction doesn't improve additive model

# is family significant?
dila_fam_PDT_Tonly<-lmer(pdays~Tratamiento_temperatura_C+(1|Pareja),
                       data=dila[dila$emerge==1 &
                                   as.numeric(as.character(dila$Pareja)) %in% 
                                   as.numeric(as.character(dila_a_fam_list$Pareja))&
                                   (as.numeric(as.character(dila$Tratamiento_temperatura_C))==20 |
                                      as.numeric(as.character(dila$Tratamiento_temperatura_C))==25),],
                       na.action="na.omit",
                       REML=FALSE)
summary(dila_fam_PDT_Tonly)
anova(dila_fam_PDT_Tonly)

anova(dila_fam_PDT_Tonly,dila_fam_PDT_add)
# no

# random-only model
dila_fam_PDT_rand<-lmer(pdays~(1|Pareja),
                         data=dila[dila$emerge==1 &
                                     as.numeric(as.character(dila$Pareja)) %in% 
                                     as.numeric(as.character(dila_a_fam_list$Pareja))&
                                     (as.numeric(as.character(dila$Tratamiento_temperatura_C))==20 |
                                        as.numeric(as.character(dila$Tratamiento_temperatura_C))==25),],
                         na.action="na.omit",
                         REML=FALSE)
summary(dila_fam_PDT_rand)
anova(dila_fam_PDT_rand)

anova(dila_fam_PDT_Tonly,dila_fam_PDT_rand)
# temperature's significant

# dorsalis mixed effects models ----------------------------------------------
# Day 1 length
View(dors_T1_fam_list)
# we should omit the 30C data

dors_fam_T1_full<-lmer(log(Tamano_dia_1_mm)~Tratamiento_temperatura_C*Pareja+(1|Pareja),
                       data=dors_lim[as.numeric(as.character(dors_lim$Tratamiento_temperatura_C))<30,],
                       na.action="na.omit",
                       REML=FALSE)
summary(dors_fam_T1_full)
anova(dors_fam_T1_full)

dors_fam_T1_add<-lmer(log(Tamano_dia_1_mm)~Tratamiento_temperatura_C+Pareja+(1|Pareja),
                      data=dors_lim[as.numeric(as.character(dors_lim$Tratamiento_temperatura_C))<30,],
                      na.action="na.omit",
                      REML=FALSE)
summary(dors_fam_T1_add)
anova(dors_fam_T1_add)

anova(dors_fam_T1_full,dors_fam_T1_add)
# so the interaction doesn't improve the additive model

# what happens if we remove Pareja as a fixed effect?
dors_fam_T1_Tonly<-lmer(log(Tamano_dia_1_mm)~Tratamiento_temperatura_C+(1|Pareja),
                        data=dors_lim[as.numeric(as.character(dors_lim$Tratamiento_temperatura_C))<30,],
                        na.action="na.omit",
                        REML=FALSE)
summary(dors_fam_T1_Tonly)
anova(dors_fam_T1_Tonly)
anova(dors_fam_T1_Tonly, dors_fam_T1_add)
# the additive model is better

# what about a random-only model?
dors_fam_T1_rand<-lmer(log(Tamano_dia_1_mm)~(1|Pareja),
                       data=dors_lim[as.numeric(as.character(dors_lim$Tratamiento_temperatura_C))<30,],
                       na.action="na.omit",
                       REML=FALSE)
summary(dors_fam_T1_rand)
anova(dors_fam_T1_rand)
anova(dors_fam_T1_rand, dors_fam_T1_Tonly)
# Temperature's significant

# Dec. 19, 2018: But if Pareja is significant, then the significance of
# temperature should be tested by comparing the additive model to a Pareja-only
# model
dors_fam_T1_Ponly<-lmer(log(Tamano_dia_1_mm)~Pareja+(1|Pareja),
                       data=dors_lim[as.numeric(as.character(dors_lim$Tratamiento_temperatura_C))<30,],
                       na.action="na.omit",
                       REML=FALSE) 
summary(dors_fam_T1_Ponly)
anova(dors_fam_T1_Ponly)
anova(dors_fam_T1_Ponly, dors_fam_T1_add)
# Temperature is still highly significant

# Day 16 lengths
dors_T15_fam_list
# only 4 families with 15C; just do 20 vs 25.

dors_fam_T15_full<-lmer(log(Tamano_dia_15_mm)~Tratamiento_temperatura_C*Pareja+(1|Pareja),
                        data=dors[as.numeric(as.character(dors$Pareja)) %in% 
                                    as.numeric(as.character(dors_T15_fam_list$Pareja)) &
                                    (as.numeric(as.character(dors$Tratamiento_temperatura_C))==20 |
                                       as.numeric(as.character(dors$Tratamiento_temperatura_C))==25),],
                        na.action="na.omit",
                        REML=FALSE)
summary(dors_fam_T15_full)
anova(dors_fam_T15_full)

dors_fam_T15_add<-lmer(log(Tamano_dia_15_mm)~Tratamiento_temperatura_C+Pareja+(1|Pareja),
                       data=dors[as.numeric(as.character(dors$Pareja)) %in% 
                                   as.numeric(as.character(dors_T15_fam_list$Pareja)) &
                                   (as.numeric(as.character(dors$Tratamiento_temperatura_C))==20 |
                                      as.numeric(as.character(dors$Tratamiento_temperatura_C))==25),],
                       na.action="na.omit",
                       REML=FALSE)
summary(dors_fam_T15_add)
anova(dors_fam_T15_add)

anova(dors_fam_T15_full,dors_fam_T15_add)
# the interaction is significant!

# what happens if we remove Pareja as a fixed effect?
dors_fam_T15_Tonly<-lmer(log(Tamano_dia_15_mm)~Tratamiento_temperatura_C+(1|Pareja),
                         data=dors[as.numeric(as.character(dors$Pareja)) %in% 
                                     as.numeric(as.character(dors_T15_fam_list$Pareja)) &
                                     (as.numeric(as.character(dors$Tratamiento_temperatura_C))==20 |
                                        as.numeric(as.character(dors$Tratamiento_temperatura_C))==25),],
                         na.action="na.omit",
                         REML=FALSE)
summary(dors_fam_T15_Tonly)
anova(dors_fam_T15_Tonly)
anova(dors_fam_T15_Tonly, dors_fam_T15_add)
# Including Pareja as a fixed effect isn't (quite) significant

# what about a random-only model?
dors_fam_T15_rand<-lmer(log(Tamano_dia_15_mm)~(1|Pareja),
                        data=dors[as.numeric(as.character(dors$Pareja)) %in% 
                                    as.numeric(as.character(dors_T15_fam_list$Pareja)) &
                                    (as.numeric(as.character(dors$Tratamiento_temperatura_C))==20 |
                                       as.numeric(as.character(dors$Tratamiento_temperatura_C))==25),],
                        na.action="na.omit",
                         REML=FALSE)
summary(dors_fam_T15_rand)
anova(dors_fam_T15_rand)
anova(dors_fam_T15_rand, dors_fam_T15_Tonly)
# Temperature's significant

# Dec. 19, 2018: But if Pareja is significant, then the significance of
# temperature should be tested by comparing the additive model to a Pareja-only
# model
dors_fam_T15_Ponly<-lmer(log(Tamano_dia_15_mm)~Pareja+(1|Pareja),
                        data=dors[as.numeric(as.character(dors$Pareja)) %in% 
                                    as.numeric(as.character(dors_T15_fam_list$Pareja)) &
                                    (as.numeric(as.character(dors$Tratamiento_temperatura_C))==20 |
                                       as.numeric(as.character(dors$Tratamiento_temperatura_C))==25),],
                        na.action="na.omit",
                        REML=FALSE)
summary(dors_fam_T15_Ponly)
anova(dors_fam_T15_Ponly)
anova(dors_fam_T15_Ponly, dors_fam_T15_add)
# Temperature's still significant

# Pupal length
dors_p_fam_list
# need to limit these families to 20 and 25C only

dors_fam_Tp_full<-lmer(log(Tamano_pupa_mm)~Tratamiento_temperatura_C*Pareja+(1|Pareja),
                       data=dors[as.numeric(as.character(dors$Pareja)) %in% 
                                   as.numeric(as.character(dors_p_fam_list$Pareja))&
                                   (as.numeric(as.character(dors$Tratamiento_temperatura_C))==20 |
                                      as.numeric(as.character(dors$Tratamiento_temperatura_C))==25),],
                       na.action="na.omit",
                       REML=FALSE)
summary(dors_fam_Tp_full)
anova(dors_fam_Tp_full)

dors_fam_Tp_add<-lmer(log(Tamano_pupa_mm)~Tratamiento_temperatura_C+Pareja+(1|Pareja),
                      data=dors[as.numeric(as.character(dors$Pareja)) %in% 
                                  as.numeric(as.character(dors_p_fam_list$Pareja))&
                                  (as.numeric(as.character(dors$Tratamiento_temperatura_C))==20 |
                                     as.numeric(as.character(dors$Tratamiento_temperatura_C))==25),],
                      na.action="na.omit",
                      REML=FALSE)
summary(dors_fam_Tp_add)
anova(dors_fam_Tp_add)

anova(dors_fam_Tp_full,dors_fam_Tp_add)
# the interaction doesn't improve the additive model

# what happens if we remove Pareja as a fixed effect?
dors_fam_Tp_Tonly<-lmer(log(Tamano_pupa_mm)~Tratamiento_temperatura_C+(1|Pareja),
                        data=dors[as.numeric(as.character(dors$Pareja)) %in% 
                                    as.numeric(as.character(dors_p_fam_list$Pareja))&
                                    (as.numeric(as.character(dors$Tratamiento_temperatura_C))==20 |
                                       as.numeric(as.character(dors$Tratamiento_temperatura_C))==25),],
                        na.action="na.omit",
                        REML=FALSE)
summary(dors_fam_Tp_Tonly)
anova(dors_fam_Tp_Tonly)
anova(dors_fam_Tp_Tonly, dors_fam_Tp_add)
# family isn't significant

# what about a random-only model?
dors_fam_Tp_rand<-lmer(log(Tamano_pupa_mm)~(1|Pareja),
                       data=dors[as.numeric(as.character(dors$Pareja)) %in% 
                                   as.numeric(as.character(dors_p_fam_list$Pareja))&
                                   (as.numeric(as.character(dors$Tratamiento_temperatura_C))==20 |
                                      as.numeric(as.character(dors$Tratamiento_temperatura_C))==25),],
                       na.action="na.omit",
                        REML=FALSE)
summary(dors_fam_Tp_rand)
anova(dors_fam_Tp_rand)
anova(dors_fam_Tp_rand, dors_fam_Tp_Tonly)
# temperature doesn't quite add anything to the model, either!

# Adult length
dors_a_fam_list
# need to limit these families to 20 and 25C data

dors_fam_Ta_full<-lmer(log(Tamano_adulto_mm)~Tratamiento_temperatura_C*Pareja+(1|Pareja),
                       data=dors[as.numeric(as.character(dors$Pareja)) %in% 
                                   as.numeric(as.character(dors_a_fam_list$Pareja))&
                                   (as.numeric(as.character(dors$Tratamiento_temperatura_C))==20 |
                                      as.numeric(as.character(dors$Tratamiento_temperatura_C))==25),],
                       na.action="na.omit",
                       REML=FALSE)
summary(dors_fam_Ta_full)
anova(dors_fam_Ta_full)

dors_fam_Ta_add<-lmer(log(Tamano_adulto_mm)~Tratamiento_temperatura_C+Pareja+(1|Pareja),
                      data=dors[as.numeric(as.character(dors$Pareja)) %in% 
                                  as.numeric(as.character(dors_a_fam_list$Pareja))&
                                  (as.numeric(as.character(dors$Tratamiento_temperatura_C))==20 |
                                     as.numeric(as.character(dors$Tratamiento_temperatura_C))==25),],
                      na.action="na.omit",
                      REML=FALSE)
summary(dors_fam_Ta_add)
anova(dors_fam_Ta_add)

anova(dors_fam_Ta_full,dors_fam_Ta_add)
# the interaction doesn't improve the additive model

# what happens if we remove Pareja as a fixed effect?
dors_fam_Ta_Tonly<-lmer(log(Tamano_adulto_mm)~Tratamiento_temperatura_C+(1|Pareja),
                        data=dors[as.numeric(as.character(dors$Pareja)) %in% 
                                    as.numeric(as.character(dors_a_fam_list$Pareja))&
                                    (as.numeric(as.character(dors$Tratamiento_temperatura_C))==20 |
                                       as.numeric(as.character(dors$Tratamiento_temperatura_C))==25),],
                        na.action="na.omit",
                        REML=FALSE)
summary(dors_fam_Ta_Tonly)
anova(dors_fam_Ta_Tonly)
anova(dors_fam_Ta_Tonly, dors_fam_Ta_add)
# family isn't (quite) significant

# what about the random-only model?
dors_fam_Ta_rand<-lmer(log(Tamano_adulto_mm)~(1|Pareja),
                       data=dors[as.numeric(as.character(dors$Pareja)) %in% 
                                   as.numeric(as.character(dors_a_fam_list$Pareja))&
                                   (as.numeric(as.character(dors$Tratamiento_temperatura_C))==20 |
                                      as.numeric(as.character(dors$Tratamiento_temperatura_C))==25),],
                       na.action="na.omit",
                        REML=FALSE)
summary(dors_fam_Ta_rand)
anova(dors_fam_Ta_rand)
anova(dors_fam_Ta_rand, dors_fam_Ta_Tonly)
# huh. temperature doesn't add anything to the model either(!)


# Pupal mass
dors_p_fam_list
# again, needs to be limited to 20 and 25C

dors_fam_Pp_full<-lmer(Peso_pupa_g~Tratamiento_temperatura_C*Pareja+(1|Pareja),
                       data=dors[as.numeric(as.character(dors$Pareja)) %in% 
                                   as.numeric(as.character(dors_p_fam_list$Pareja))&
                                   (as.numeric(as.character(dors$Tratamiento_temperatura_C))==20 |
                                      as.numeric(as.character(dors$Tratamiento_temperatura_C))==25),],
                       na.action="na.omit",
                       REML=FALSE)
summary(dors_fam_Pp_full)
anova(dors_fam_Pp_full)

dors_fam_Pp_add<-lmer(Peso_pupa_g~Tratamiento_temperatura_C+Pareja+(1|Pareja),
                      data=dors[as.numeric(as.character(dors$Pareja)) %in% 
                                  as.numeric(as.character(dors_p_fam_list$Pareja))&
                                  (as.numeric(as.character(dors$Tratamiento_temperatura_C))==20 |
                                     as.numeric(as.character(dors$Tratamiento_temperatura_C))==25),],
                      na.action="na.omit",
                      REML=FALSE)
summary(dors_fam_Pp_add)
anova(dors_fam_Pp_add)

anova(dors_fam_Pp_full,dors_fam_Pp_add)
# the interaction doesn't improve the additive model

# what happens if we remove Pareja as a fixed effect?
dors_fam_Pp_Tonly<-lmer(Peso_pupa_g~Tratamiento_temperatura_C+(1|Pareja),
                        data=dors[as.numeric(as.character(dors$Pareja)) %in% 
                                    as.numeric(as.character(dors_p_fam_list$Pareja))&
                                    (as.numeric(as.character(dors$Tratamiento_temperatura_C))==20 |
                                       as.numeric(as.character(dors$Tratamiento_temperatura_C))==25),],
                        na.action="na.omit",
                        REML=FALSE)
summary(dors_fam_Pp_Tonly)
anova(dors_fam_Pp_Tonly)
anova(dors_fam_Pp_Tonly, dors_fam_Pp_add)
# family is significant

# what about a random-only model?
dors_fam_Pp_rand<-lmer(Peso_pupa_g~(1|Pareja),
                       data=dors[as.numeric(as.character(dors$Pareja)) %in% 
                                   as.numeric(as.character(dors_p_fam_list$Pareja))&
                                   (as.numeric(as.character(dors$Tratamiento_temperatura_C))==20 |
                                      as.numeric(as.character(dors$Tratamiento_temperatura_C))==25),],
                       na.action="na.omit",
                        REML=FALSE)
summary(dors_fam_Pp_rand)
anova(dors_fam_Pp_rand)
anova(dors_fam_Pp_rand, dors_fam_Pp_Tonly)
# temperature isn't significant!

# Dec. 19, 2018: But if Pareja is significant, then the significance of
# temperature should be tested by comparing the additive model to a Pareja-only
# model
dors_fam_Pp_Ponly<-lmer(Peso_pupa_g~Pareja+(1|Pareja),
                       data=dors[as.numeric(as.character(dors$Pareja)) %in% 
                                   as.numeric(as.character(dors_p_fam_list$Pareja))&
                                   (as.numeric(as.character(dors$Tratamiento_temperatura_C))==20 |
                                      as.numeric(as.character(dors$Tratamiento_temperatura_C))==25),],
                       na.action="na.omit",
                       REML=FALSE)
summary(dors_fam_Pp_Ponly)
anova(dors_fam_Pp_Ponly)
anova(dors_fam_Pp_Ponly, dors_fam_Pp_add)
# temperature still isn't significant


# Adult mass
dors_a_fam_list
# need to limit data to 20 and 25C

dors_fam_Pa_full<-lmer(Peso_adulto_g~Tratamiento_temperatura_C*Pareja+(1|Pareja),
                       data=dors[as.numeric(as.character(dors$Pareja)) %in% 
                                   as.numeric(as.character(dors_a_fam_list$Pareja))&
                                   (as.numeric(as.character(dors$Tratamiento_temperatura_C))==20 |
                                      as.numeric(as.character(dors$Tratamiento_temperatura_C))==25),],
                       na.action="na.omit",
                       REML=FALSE)
summary(dors_fam_Pa_full)
anova(dors_fam_Pa_full)

dors_fam_Pa_add<-lmer(Peso_adulto_g~Tratamiento_temperatura_C+Pareja+(1|Pareja),
                      data=dors[as.numeric(as.character(dors$Pareja)) %in% 
                                  as.numeric(as.character(dors_a_fam_list$Pareja))&
                                  (as.numeric(as.character(dors$Tratamiento_temperatura_C))==20 |
                                     as.numeric(as.character(dors$Tratamiento_temperatura_C))==25),],
                      na.action="na.omit",
                      REML=FALSE)
summary(dors_fam_Pa_add)
anova(dors_fam_Pa_add)

anova(dors_fam_Pa_full,dors_fam_Pa_add)
# the interaction doesn't improve the additive model

# what happens if we remove Pareja as a fixed effect?
dors_fam_Pa_Tonly<-lmer(Peso_adulto_g~Tratamiento_temperatura_C+(1|Pareja),
                        data=dors[as.numeric(as.character(dors$Pareja)) %in% 
                                    as.numeric(as.character(dors_a_fam_list$Pareja))&
                                    (as.numeric(as.character(dors$Tratamiento_temperatura_C))==20 |
                                       as.numeric(as.character(dors$Tratamiento_temperatura_C))==25),],
                        na.action="na.omit",
                        REML=FALSE)
summary(dors_fam_Pa_Tonly)
anova(dors_fam_Pa_Tonly)
anova(dors_fam_Pa_Tonly, dors_fam_Pa_add)
# family isn't significant

# random-only model?
dors_fam_Pa_rand<-lmer(Peso_adulto_g~(1|Pareja),
                       data=dors[as.numeric(as.character(dors$Pareja)) %in% 
                                   as.numeric(as.character(dors_a_fam_list$Pareja))&
                                   (as.numeric(as.character(dors$Tratamiento_temperatura_C))==20 |
                                      as.numeric(as.character(dors$Tratamiento_temperatura_C))==25),],
                       na.action="na.omit",
                        REML=FALSE)
summary(dors_fam_Pa_rand)
anova(dors_fam_Pa_rand)
anova(dors_fam_Pa_rand, dors_fam_Pa_Tonly)
# temperature isn't significant either!

# larval development time
# the balanced families for this will be the same as those that had balanced
# pupal measurements
dors_p_fam_list

dors_fam_LDT_full<-lmer(ldays~Tratamiento_temperatura_C*Pareja+(1|Pareja),
                        data=dors[dors$pupate==1 &
                                    as.numeric(as.character(dors$Pareja)) %in% 
                                    as.numeric(as.character(dors_p_fam_list$Pareja))&
                                    (as.numeric(as.character(dors$Tratamiento_temperatura_C))==20 |
                                       as.numeric(as.character(dors$Tratamiento_temperatura_C))==25),],
                        na.action="na.omit",
                        REML=FALSE)
summary(dors_fam_LDT_full)
anova(dors_fam_LDT_full)

dors_fam_LDT_add<-lmer(ldays~Tratamiento_temperatura_C+Pareja+(1|Pareja),
                       data=dors[dors$pupate==1 &
                                   as.numeric(as.character(dors$Pareja)) %in% 
                                   as.numeric(as.character(dors_p_fam_list$Pareja))&
                                   (as.numeric(as.character(dors$Tratamiento_temperatura_C))==20 |
                                      as.numeric(as.character(dors$Tratamiento_temperatura_C))==25),],
                       na.action="na.omit",
                       REML=FALSE)
summary(dors_fam_LDT_add)
anova(dors_fam_LDT_add)

anova(dors_fam_LDT_full,dors_fam_LDT_add)
# the interaction isn't significant

# what happens if we remove Pareja as a fixed effect?
dors_fam_LDT_Tonly<-lmer(ldays~Tratamiento_temperatura_C+(1|Pareja),
                         data=dors[dors$pupate==1 &
                                     as.numeric(as.character(dors$Pareja)) %in% 
                                     as.numeric(as.character(dors_p_fam_list$Pareja))&
                                     (as.numeric(as.character(dors$Tratamiento_temperatura_C))==20 |
                                        as.numeric(as.character(dors$Tratamiento_temperatura_C))==25),],
                         na.action="na.omit",
                         REML=FALSE)
summary(dors_fam_LDT_Tonly)
anova(dors_fam_LDT_Tonly)
anova(dors_fam_LDT_Tonly, dors_fam_LDT_add)
# family is significant

# random-only model?
dors_fam_LDT_rand<-lmer(ldays~(1|Pareja),
                        data=dors[dors$pupate==1 &
                                    as.numeric(as.character(dors$Pareja)) %in% 
                                    as.numeric(as.character(dors_p_fam_list$Pareja))&
                                    (as.numeric(as.character(dors$Tratamiento_temperatura_C))==20 |
                                       as.numeric(as.character(dors$Tratamiento_temperatura_C))==25),],
                        na.action="na.omit",
                         REML=FALSE)
summary(dors_fam_LDT_rand)
anova(dors_fam_LDT_rand)
anova(dors_fam_LDT_rand, dors_fam_LDT_Tonly)
# temperature is significant

# Dec. 19, 2018: But if Pareja is significant, then the significance of
# temperature should be tested by comparing the additive model to a Pareja-only
# model
dors_fam_LDT_Ponly<-lmer(ldays~Pareja+(1|Pareja),
                        data=dors[dors$pupate==1 &
                                    as.numeric(as.character(dors$Pareja)) %in% 
                                    as.numeric(as.character(dors_p_fam_list$Pareja))&
                                    (as.numeric(as.character(dors$Tratamiento_temperatura_C))==20 |
                                       as.numeric(as.character(dors$Tratamiento_temperatura_C))==25),],
                        na.action="na.omit",
                        REML=FALSE)
summary(dors_fam_LDT_Ponly)
anova(dors_fam_LDT_Ponly)
anova(dors_fam_LDT_Ponly, dors_fam_LDT_add)
# temperature is still highly significant

# Pupal development time
# the balanced families for this will be the same as those that had balanced
# adult measurements
dors_a_fam_list

dors_fam_PDT_full<-lmer(pdays~Tratamiento_temperatura_C*Pareja+(1|Pareja),
                        data=dors[dors$emerge==1 &
                                    as.numeric(as.character(dors$Pareja)) %in% 
                                    as.numeric(as.character(dors_a_fam_list$Pareja))&
                                    (as.numeric(as.character(dors$Tratamiento_temperatura_C))==20 |
                                       as.numeric(as.character(dors$Tratamiento_temperatura_C))==25),],
                        na.action="na.omit",
                        REML=FALSE)
summary(dors_fam_PDT_full)
anova(dors_fam_PDT_full)

dors_fam_PDT_add<-lmer(pdays~Tratamiento_temperatura_C+Pareja+(1|Pareja),
                       data=dors[dors$emerge==1 &
                                   as.numeric(as.character(dors$Pareja)) %in% 
                                   as.numeric(as.character(dors_a_fam_list$Pareja))&
                                   (as.numeric(as.character(dors$Tratamiento_temperatura_C))==20 |
                                      as.numeric(as.character(dors$Tratamiento_temperatura_C))==25),],
                       na.action="na.omit",
                       REML=FALSE)
summary(dors_fam_PDT_add)
anova(dors_fam_PDT_add)

anova(dors_fam_PDT_full,dors_fam_PDT_add)
# the interaction isn't quite significant

# what happens if we remove Pareja as a fixed effect?
dors_fam_PDT_Tonly<-lmer(pdays~Tratamiento_temperatura_C+(1|Pareja),
                         data=dors[dors$emerge==1 &
                                     as.numeric(as.character(dors$Pareja)) %in% 
                                     as.numeric(as.character(dors_a_fam_list$Pareja))&
                                     (as.numeric(as.character(dors$Tratamiento_temperatura_C))==20 |
                                        as.numeric(as.character(dors$Tratamiento_temperatura_C))==25),],
                         na.action="na.omit",
                         REML=FALSE)
summary(dors_fam_PDT_Tonly)
anova(dors_fam_PDT_Tonly)
anova(dors_fam_PDT_Tonly, dors_fam_PDT_add)
# family isn't significant

# random-only model?
dors_fam_PDT_rand<-lmer(pdays~(1|Pareja),
                        data=dors[dors$emerge==1 &
                                    as.numeric(as.character(dors$Pareja)) %in% 
                                    as.numeric(as.character(dors_a_fam_list$Pareja))&
                                    (as.numeric(as.character(dors$Tratamiento_temperatura_C))==20 |
                                       as.numeric(as.character(dors$Tratamiento_temperatura_C))==25),],
                        na.action="na.omit",
                         REML=FALSE)
summary(dors_fam_PDT_rand)
anova(dors_fam_PDT_rand)
anova(dors_fam_PDT_rand, dors_fam_PDT_Tonly)
# temperature is significant

# placida mixed effect models-----------------------------------------------
# Day 1 length
plac_T1_fam_list
# only three families representing 30C, so we should only consider 15-25C

plac_fam_T1_full<-lmer(log(Tamano_dia_1_mm)~Tratamiento_temperatura_C*Pareja+(1|Pareja),
                       data=plac_lim[as.numeric(as.character(plac_lim$Tratamiento_temperatura_C))<30,],
                       na.action="na.omit",
                       REML=FALSE)
summary(plac_fam_T1_full)
anova(plac_fam_T1_full)

plac_fam_T1_add<-lmer(log(Tamano_dia_1_mm)~Tratamiento_temperatura_C+Pareja+(1|Pareja),
                      data=plac_lim[as.numeric(as.character(plac_lim$Tratamiento_temperatura_C))<30,],
                      na.action="na.omit",
                      REML=FALSE)
summary(plac_fam_T1_add)
anova(plac_fam_T1_add)

anova(plac_fam_T1_full,plac_fam_T1_add)
# so the interaction doesn't improve the additive model

# what happens if we remove Pareja as a fixed effect?
plac_fam_T1_Tonly<-lmer(log(Tamano_dia_1_mm)~Tratamiento_temperatura_C+(1|Pareja),
                        data=plac_lim[as.numeric(as.character(plac_lim$Tratamiento_temperatura_C))<30,],
                        na.action="na.omit",
                        REML=FALSE)
summary(plac_fam_T1_Tonly)
anova(plac_fam_T1_Tonly)
anova(plac_fam_T1_Tonly, plac_fam_T1_add)
# the additive model is better

# random-only model?
plac_fam_T1_rand<-lmer(log(Tamano_dia_1_mm)~(1|Pareja),
                        data=plac_lim[as.numeric(as.character(plac_lim$Tratamiento_temperatura_C))<30,],
                        na.action="na.omit",
                        REML=FALSE)
summary(plac_fam_T1_rand)
anova(plac_fam_T1_rand)
anova(plac_fam_T1_rand, plac_fam_T1_Tonly)
# Temperature is significant

# Dec. 19, 2018: But if Pareja is significant, then the significance of
# temperature should be tested by comparing the additive model to a Pareja-only
# model
plac_fam_T1_Ponly<-lmer(log(Tamano_dia_1_mm)~Pareja+(1|Pareja),
                       data=plac_lim[as.numeric(as.character(plac_lim$Tratamiento_temperatura_C))<30,],
                       na.action="na.omit",
                       REML=FALSE)
summary(plac_fam_T1_Ponly)
anova(plac_fam_T1_Ponly)
anova(plac_fam_T1_Ponly, plac_fam_T1_add)
# Temperature's still significant

# Day 16 length
plac_T15_fam_list
# again, omit 30C (only one family)

plac_fam_T15_full<-lmer(log(Tamano_dia_15_mm)~Tratamiento_temperatura_C*Pareja+(1|Pareja),
                        data=plac[as.numeric(as.character(plac$Pareja)) %in% 
                                    as.numeric(as.character(plac_T15_fam_list$Pareja)) &
                                    as.numeric(as.character(plac$Tratamiento_temperatura_C))!=30,],
                        na.action="na.omit",
                        REML=FALSE)
summary(plac_fam_T15_full)
anova(plac_fam_T15_full)

plac_fam_T15_add<-lmer(log(Tamano_dia_15_mm)~Tratamiento_temperatura_C+Pareja+(1|Pareja),
                       data=plac[as.numeric(as.character(plac$Pareja)) %in% 
                                   as.numeric(as.character(plac_T15_fam_list$Pareja)) &
                                   as.numeric(as.character(plac$Tratamiento_temperatura_C))!=30,],
                       na.action="na.omit",
                       REML=FALSE)
summary(plac_fam_T15_add)
anova(plac_fam_T15_add)

anova(plac_fam_T15_full,plac_fam_T15_add)
# the interaction is significant!

# what happens if we remove Pareja as a fixed effect?
plac_fam_T15_Tonly<-lmer(log(Tamano_dia_15_mm)~Tratamiento_temperatura_C+(1|Pareja),
                         data=plac[as.numeric(as.character(plac$Pareja)) %in% 
                                     as.numeric(as.character(plac_T15_fam_list$Pareja)) &
                                     as.numeric(as.character(plac$Tratamiento_temperatura_C))!=30,],
                         na.action="na.omit",
                         REML=FALSE)
summary(plac_fam_T15_Tonly)
anova(plac_fam_T15_Tonly)
anova(plac_fam_T15_Tonly, plac_fam_T15_add)
# Including Pareja as a fixed effect alone isn't significant

# random-only model?
plac_fam_T15_rand<-lmer(log(Tamano_dia_15_mm)~(1|Pareja),
                        data=plac[as.numeric(as.character(plac$Pareja)) %in% 
                                    as.numeric(as.character(plac_T15_fam_list$Pareja)) &
                                    as.numeric(as.character(plac$Tratamiento_temperatura_C))!=30,],
                        na.action="na.omit",
                         REML=FALSE)
summary(plac_fam_T15_rand)
anova(plac_fam_T15_rand)
anova(plac_fam_T15_rand, plac_fam_T15_Tonly)
# Temperature's significant

# Pupal length
plac_p_fam_list
# we'll leave in 15C

plac_fam_Tp_full<-lmer(log(Tamano_pupa_mm)~Tratamiento_temperatura_C*Pareja+(1|Pareja),
                       data=plac[as.numeric(as.character(plac$Pareja)) %in% 
                                   as.numeric(as.character(plac_p_fam_list$Pareja)),],
                       na.action="na.omit",
                       REML=FALSE)
summary(plac_fam_Tp_full)
anova(plac_fam_Tp_full)

plac_fam_Tp_add<-lmer(log(Tamano_pupa_mm)~Tratamiento_temperatura_C+Pareja+(1|Pareja),
                      data=plac[as.numeric(as.character(plac$Pareja)) %in% 
                                  as.numeric(as.character(plac_p_fam_list$Pareja)),],
                      na.action="na.omit",
                      REML=FALSE)
summary(plac_fam_Tp_add)
anova(plac_fam_Tp_add)

anova(plac_fam_Tp_full,plac_fam_Tp_add)
# the interaction doesn't (quite) improve the additive model

# what happens if we remove Pareja as a fixed effect?
plac_fam_Tp_Tonly<-lmer(log(Tamano_pupa_mm)~Tratamiento_temperatura_C+(1|Pareja),
                        data=plac[as.numeric(as.character(plac$Pareja)) %in% 
                                    as.numeric(as.character(plac_p_fam_list$Pareja)),],
                        na.action="na.omit",
                        REML=FALSE)
summary(plac_fam_Tp_Tonly)
anova(plac_fam_Tp_Tonly)
anova(plac_fam_Tp_Tonly, plac_fam_Tp_add)
# family isn't significant

# random-only model?
plac_fam_Tp_rand<-lmer(log(Tamano_pupa_mm)~(1|Pareja),
                        data=plac[as.numeric(as.character(plac$Pareja)) %in% 
                                    as.numeric(as.character(plac_p_fam_list$Pareja)),],
                        na.action="na.omit",
                        REML=FALSE)
summary(plac_fam_Tp_rand)
anova(plac_fam_Tp_rand)
anova(plac_fam_Tp_rand, plac_fam_Tp_Tonly)
# Temperature's significant

# Adult length
plac_a_fam_list
# we should only include the 20 and 25 treatments for this, I think (only 3 families
# have enough for 15C)

plac_fam_Ta_full<-lmer(log(Tamano_adulto_mm)~Tratamiento_temperatura_C*Pareja+(1|Pareja),
                       data=plac[as.numeric(as.character(plac$Tratamiento_temperatura_C))>15 &
                         as.numeric(as.character(plac$Pareja)) %in% 
                                   as.numeric(as.character(plac_a_fam_list$Pareja)),],
                       na.action="na.omit",
                       REML=FALSE)
summary(plac_fam_Ta_full)
anova(plac_fam_Ta_full)

plac_fam_Ta_add<-lmer(log(Tamano_adulto_mm)~Tratamiento_temperatura_C+Pareja+(1|Pareja),
                      data=plac[as.numeric(as.character(plac$Tratamiento_temperatura_C))>15 &
                                  as.numeric(as.character(plac$Pareja)) %in% 
                                  as.numeric(as.character(plac_a_fam_list$Pareja)),],
                      na.action="na.omit",
                      REML=FALSE)
summary(plac_fam_Ta_add)
anova(plac_fam_Ta_add)

anova(plac_fam_Ta_full,plac_fam_Ta_add)
# the interaction doesn't improve the additive model

# what happens if we remove Pareja as a fixed effect?
plac_fam_Ta_Tonly<-lmer(log(Tamano_adulto_mm)~Tratamiento_temperatura_C+(1|Pareja),
                        data=plac[as.numeric(as.character(plac$Tratamiento_temperatura_C))>15 &
                                    as.numeric(as.character(plac$Pareja)) %in% 
                                    as.numeric(as.character(plac_a_fam_list$Pareja)),],
                        na.action="na.omit",
                        REML=FALSE)
summary(plac_fam_Ta_Tonly)
anova(plac_fam_Ta_Tonly)
anova(plac_fam_Ta_Tonly, plac_fam_Ta_add)
# family isn't significant

# random-only model?
plac_fam_Ta_rand<-lmer(log(Tamano_adulto_mm)~(1|Pareja),
                        data=plac[as.numeric(as.character(plac$Tratamiento_temperatura_C))>15 &
                                    as.numeric(as.character(plac$Pareja)) %in% 
                                    as.numeric(as.character(plac_a_fam_list$Pareja)),],
                        na.action="na.omit",
                        REML=FALSE)
summary(plac_fam_Ta_rand)
anova(plac_fam_Ta_rand)
anova(plac_fam_Ta_rand, plac_fam_Ta_Tonly)
# Temperature's significant


# Pupal mass
plac_p_fam_list

plac_fam_Pp_full<-lmer(Peso_pupa_g~Tratamiento_temperatura_C*Pareja+(1|Pareja),
                       data=plac[as.numeric(as.character(plac$Pareja)) %in% 
                                   as.numeric(as.character(plac_p_fam_list$Pareja)),],
                       na.action="na.omit",
                       REML=FALSE)
summary(plac_fam_Pp_full)
anova(plac_fam_Pp_full)

plac_fam_Pp_add<-lmer(Peso_pupa_g~Tratamiento_temperatura_C+Pareja+(1|Pareja),
                      data=plac[as.numeric(as.character(plac$Pareja)) %in% 
                                  as.numeric(as.character(plac_p_fam_list$Pareja)),],
                      na.action="na.omit",
                      REML=FALSE)
summary(plac_fam_Pp_add)
anova(plac_fam_Pp_add)

anova(plac_fam_Pp_full,plac_fam_Pp_add)
# the interaction doesn't improve the additive model

# what happens if we remove Pareja as a fixed effect?
plac_fam_Pp_Tonly<-lmer(Peso_pupa_g~Tratamiento_temperatura_C+(1|Pareja),
                        data=plac[as.numeric(as.character(plac$Pareja)) %in% 
                                    as.numeric(as.character(plac_p_fam_list$Pareja)),],
                        na.action="na.omit",
                        REML=FALSE)
summary(plac_fam_Pp_Tonly)
anova(plac_fam_Pp_Tonly)
anova(plac_fam_Pp_Tonly, plac_fam_Pp_add)
# family isn't significant

# random-only model?
plac_fam_Pp_rand<-lmer(Peso_pupa_g~(1|Pareja),
                        data=plac[as.numeric(as.character(plac$Pareja)) %in% 
                                    as.numeric(as.character(plac_p_fam_list$Pareja)),],
                        na.action="na.omit",
                        REML=FALSE)
summary(plac_fam_Pp_rand)
anova(plac_fam_Pp_rand)
anova(plac_fam_Pp_rand, plac_fam_Pp_Tonly)
# Temperature's significant


# Adult mass
plac_a_fam_list
# again, we should probably limit this to 20 and 25C

plac_fam_Pa_full<-lmer(Peso_adulto_g~Tratamiento_temperatura_C*Pareja+(1|Pareja),
                       data=plac[as.numeric(as.character(plac$Tratamiento_temperatura_C))>15 &
                                   as.numeric(as.character(plac$Pareja)) %in% 
                                   as.numeric(as.character(plac_a_fam_list$Pareja)),],
                       na.action="na.omit",
                       REML=FALSE)
summary(plac_fam_Pa_full)
anova(plac_fam_Pa_full)

plac_fam_Pa_add<-lmer(Peso_adulto_g~Tratamiento_temperatura_C+Pareja+(1|Pareja),
                      data=plac[as.numeric(as.character(plac$Tratamiento_temperatura_C))>15 &
                                  as.numeric(as.character(plac$Pareja)) %in% 
                                  as.numeric(as.character(plac_a_fam_list$Pareja)),],
                      na.action="na.omit",
                      REML=FALSE)
summary(plac_fam_Pa_add)
anova(plac_fam_Pa_add)

anova(plac_fam_Pa_full,plac_fam_Pa_add)
# the interaction doesn't improve the additive model

# what happens if we remove Pareja as a fixed effect?
plac_fam_Pa_Tonly<-lmer(Peso_adulto_g~Tratamiento_temperatura_C+(1|Pareja),
                        data=plac[as.numeric(as.character(plac$Tratamiento_temperatura_C))>15 &
                                    as.numeric(as.character(plac$Pareja)) %in% 
                                    as.numeric(as.character(plac_a_fam_list$Pareja)),],
                        na.action="na.omit",
                        REML=FALSE)
summary(plac_fam_Pa_Tonly)
anova(plac_fam_Pa_Tonly)
anova(plac_fam_Pa_Tonly, plac_fam_Pa_add)
# family isn't significant

# random-only model?
plac_fam_Pa_rand<-lmer(Peso_adulto_g~(1|Pareja),
                        data=plac[as.numeric(as.character(plac$Tratamiento_temperatura_C))>15 &
                                    as.numeric(as.character(plac$Pareja)) %in% 
                                    as.numeric(as.character(plac_a_fam_list$Pareja)),],
                        na.action="na.omit",
                        REML=FALSE)
summary(plac_fam_Pa_rand)
anova(plac_fam_Pa_rand)
anova(plac_fam_Pa_rand, plac_fam_Pa_Tonly)
# Temperature's significant

# larval development time
# the balanced families for this will be the same as those that had balanced
# pupal measurements
plac_p_fam_list

plac_fam_LDT_full<-lmer(ldays~Tratamiento_temperatura_C*Pareja+(1|Pareja),
                        data=plac[plac$pupate==1 &
                                    as.numeric(as.character(plac$Pareja)) %in% 
                                    as.numeric(as.character(plac_p_fam_list$Pareja)),],
                        na.action="na.omit",
                        REML=FALSE)
summary(plac_fam_LDT_full)
anova(plac_fam_LDT_full)

plac_fam_LDT_add<-lmer(ldays~Tratamiento_temperatura_C+Pareja+(1|Pareja),
                       data=plac[plac$pupate==1 &
                                   as.numeric(as.character(plac$Pareja)) %in% 
                                   as.numeric(as.character(plac_p_fam_list$Pareja)),],
                       na.action="na.omit",
                       REML=FALSE)
summary(plac_fam_LDT_add)
anova(plac_fam_LDT_add)

anova(plac_fam_LDT_full,plac_fam_LDT_add)
# the interaction is highly significant!

# what happens if we remove Pareja as a fixed effect?
plac_fam_LDT_Tonly<-lmer(ldays~Tratamiento_temperatura_C+(1|Pareja),
                         data=plac[plac$pupate==1 &
                                     as.numeric(as.character(plac$Pareja)) %in% 
                                     as.numeric(as.character(plac_p_fam_list$Pareja)),],
                         na.action="na.omit",
                         REML=FALSE)
summary(plac_fam_LDT_Tonly)
anova(plac_fam_LDT_Tonly)
anova(plac_fam_LDT_Tonly, plac_fam_LDT_add)
# family is significant

# random-only model?
plac_fam_LDT_rand<-lmer(ldays~(1|Pareja),
                         data=plac[plac$pupate==1 &
                                     as.numeric(as.character(plac$Pareja)) %in% 
                                     as.numeric(as.character(plac_p_fam_list$Pareja)),],
                         na.action="na.omit",
                         REML=FALSE)
summary(plac_fam_LDT_rand)
anova(plac_fam_LDT_rand)
anova(plac_fam_LDT_rand, plac_fam_LDT_Tonly)
# Temperature's significant

# Dec. 19, 2018: But if Pareja is significant, then the significance of
# temperature should be tested by comparing the additive model to a Pareja-only
# model
plac_fam_LDT_Ponly<-lmer(ldays~Pareja+(1|Pareja),
                        data=plac[plac$pupate==1 &
                                    as.numeric(as.character(plac$Pareja)) %in% 
                                    as.numeric(as.character(plac_p_fam_list$Pareja)),],
                        na.action="na.omit",
                        REML=FALSE)
summary(plac_fam_LDT_Ponly)
anova(plac_fam_LDT_Ponly)
anova(plac_fam_LDT_Ponly, plac_fam_LDT_add)
# Temperature's significant

# Pupal development time
# the balanced families for this will be the same as those that had balanced
# adult measurements (only include 20 and 25C treatments)
plac_a_fam_list

plac_fam_PDT_full<-lmer(pdays~Tratamiento_temperatura_C*Pareja+(1|Pareja),
                        data=plac[plac$emerge==1 &
                                    as.numeric(as.character(plac$Tratamiento_temperatura_C))>15 &
                                    as.numeric(as.character(plac$Pareja)) %in% 
                                    as.numeric(as.character(plac_a_fam_list$Pareja)),],
                        na.action="na.omit",
                        REML=FALSE)
summary(plac_fam_PDT_full)
anova(plac_fam_PDT_full)

plac_fam_PDT_add<-lmer(pdays~Tratamiento_temperatura_C+Pareja+(1|Pareja),
                       data=plac[plac$emerge==1 &
                                   as.numeric(as.character(plac$Tratamiento_temperatura_C))>15 &
                                   as.numeric(as.character(plac$Pareja)) %in% 
                                   as.numeric(as.character(plac_a_fam_list$Pareja)),],
                       na.action="na.omit",
                       REML=FALSE)
summary(plac_fam_PDT_add)
anova(plac_fam_PDT_add)

anova(plac_fam_PDT_full,plac_fam_PDT_add)
# the interaction isn't significant

# what happens if we remove Pareja as a fixed effect?
plac_fam_PDT_Tonly<-lmer(pdays~Tratamiento_temperatura_C+(1|Pareja),
                         data=plac[plac$emerge==1 &
                                     as.numeric(as.character(plac$Tratamiento_temperatura_C))>15 &
                                     as.numeric(as.character(plac$Pareja)) %in% 
                                     as.numeric(as.character(plac_a_fam_list$Pareja)),],
                         na.action="na.omit",
                         REML=FALSE)
summary(plac_fam_PDT_Tonly)
anova(plac_fam_PDT_Tonly)
anova(plac_fam_PDT_Tonly, plac_fam_PDT_add)
# family isn't particularly significant

# random-only model?
plac_fam_PDT_rand<-lmer(pdays~(1|Pareja),
                         data=plac[plac$emerge==1 &
                                     as.numeric(as.character(plac$Tratamiento_temperatura_C))>15 &
                                     as.numeric(as.character(plac$Pareja)) %in% 
                                     as.numeric(as.character(plac_a_fam_list$Pareja)),],
                         na.action="na.omit",
                         REML=FALSE)
summary(plac_fam_PDT_rand)
anova(plac_fam_PDT_rand)
anova(plac_fam_PDT_rand, plac_fam_PDT_Tonly)
# Temperature's significant

# Figures 5-6. Reaction norms for traits with significant GxE interactions --------
# From the mixed effect models above, we get three cases with significant
# quantifable genotype x environment interactions: dorsalis Day 16 length,
# placida Day 16 length, and placida larval development time.

# Goal: To show what these interactions look like by plotting lines connecting
# family means across each pair of temperatures (reaction norms).
# dorsalis Day 16 length reaction norm (Fig. 5)------------------------------
FTP_FT15dors

p_dors15_rn<-ggplot(data=FTP_FT15dors, aes(x=as.numeric(as.character(Temp)), 
                                           y=Tamano_dia_15_mm, 
                                           group=Pareja))+
  geom_line()+theme_classic2()+ 
  labs(x= "Temperature", y= "Mean Day 16 Length (mm)", title= "dorsalis", 
       family="sans")+
  scale_x_continuous(breaks=c(20,25))+
  expand_limits(x=c(18.5,26.5))

p_dors15_rn

# placida Day 16 length reaction norms (Fig. 6A-C)----------------------------
FTP_FT15plac

# we need to do three pairwise comparisons here
# 15 vs. 20 (Fig. 6A)
p_plac15_rn15v20<-ggplot(data=FTP_FT15plac[as.numeric(as.character(FTP_FT15plac$Temp))==15 |
                                             as.numeric(as.character(FTP_FT15plac$Temp))==20,], 
                         aes(x=as.numeric(as.character(Temp)), 
                                           y=Tamano_dia_15_mm, 
                                           group=Pareja))+
  geom_line()+theme_classic2()+ 
  labs(x= "Temperature", y= "Mean Day 16 Length (mm)", title= "placida 15 v 20", 
       family="sans")+
  scale_x_continuous(breaks=c(15,20))+
  expand_limits(x=c(13.5,21.5),y=c(2.75,4.75))+
  scale_y_continuous(breaks=pretty_breaks(4))

p_plac15_rn15v20

# 15 vs. 25 (Fig. 6B)
p_plac15_rn15v25<-ggplot(data=FTP_FT15plac[as.numeric(as.character(FTP_FT15plac$Temp))==15 |
                                             as.numeric(as.character(FTP_FT15plac$Temp))==25,], 
                         aes(x=as.numeric(as.character(Temp)), 
                             y=Tamano_dia_15_mm, 
                             group=Pareja))+
  geom_line()+theme_classic2()+ 
  labs(x= "Temperature", y= "Mean Day 16 Length (mm)", title= "placida 15 v 25", 
       family="sans")+
  scale_x_continuous(breaks=c(15,25))+
  expand_limits(x=c(12,28),y=c(2.75,4.75))+
  scale_y_continuous(breaks=pretty_breaks(4))

p_plac15_rn15v25

# 20 vs. 25 (Fig. 6C)
p_plac15_rn20v25<-ggplot(data=FTP_FT15plac[as.numeric(as.character(FTP_FT15plac$Temp))==20 |
                                             as.numeric(as.character(FTP_FT15plac$Temp))==25,], 
                         aes(x=as.numeric(as.character(Temp)), 
                             y=Tamano_dia_15_mm, 
                             group=Pareja))+
  geom_line()+theme_classic2()+ 
  labs(x= "Temperature", y= "Mean Day 16 Length (mm)", title= "placida 20 v 25",
       family="sans")+
  scale_x_continuous(breaks=c(20,25))+
  expand_limits(x=c(18.5,26.5),y=c(2.75,4.75))+
  scale_y_continuous(breaks=pretty_breaks(4))

p_plac15_rn20v25

# placida larval development time reaction norms (Fig. 6D-F)------------------
FTP_FLDTplac

# 15 vs. 20 (Fig. 6D)
p_placLDT_rn15v20<-ggplot(data=FTP_FLDTplac[as.numeric(as.character(FTP_FLDTplac$Temp))==15 |
                                             as.numeric(as.character(FTP_FLDTplac$Temp))==20,], 
                         aes(x=as.numeric(as.character(Temp)), 
                             y=ldays, 
                             group=Pareja))+
  geom_line()+theme_classic2()+ 
  labs(x= "Temperature", y= "Mean Larval Development Time (days)", 
       title= "placida 15 v 20", family="sans")+
  scale_x_continuous(breaks=c(15,20))+
  expand_limits(x=c(13.5,21.5),y=c(40,200))+
  scale_y_continuous(breaks=pretty_breaks(4))

p_placLDT_rn15v20

# 15 vs. 25 (Fig. 6E)
p_placLDT_rn15v25<-ggplot(data=FTP_FLDTplac[as.numeric(as.character(FTP_FLDTplac$Temp))==15 |
                                             as.numeric(as.character(FTP_FLDTplac$Temp))==25,], 
                         aes(x=as.numeric(as.character(Temp)), 
                             y=ldays, 
                             group=Pareja))+
  geom_line()+theme_classic2()+ 
  labs(x= "Temperature", y= "Mean Larval Development Time (days)", 
       title= "placida 15 v 25", family="sans")+
  scale_x_continuous(breaks=c(15,25))+
  expand_limits(x=c(12,28),y=c(40,200))+
  scale_y_continuous(breaks=pretty_breaks(4))

p_placLDT_rn15v25

# 20 vs. 25 (Fig. 6F)
p_placLDT_rn20v25<-ggplot(data=FTP_FLDTplac[as.numeric(as.character(FTP_FLDTplac$Temp))==20 |
                                             as.numeric(as.character(FTP_FLDTplac$Temp))==25,], 
                         aes(x=as.numeric(as.character(Temp)), 
                             y=ldays, 
                             group=Pareja))+
  geom_line()+theme_classic2()+ 
  labs(x= "Temperature", y= "Mean Larval Development Time (days)", 
       title= "placida 20 v 25", family="sans")+
  scale_x_continuous(breaks=c(20,25))+
  expand_limits(x=c(18.5,26.5),y=c(40,200))+
  scale_y_continuous(breaks=pretty_breaks(4))

p_placLDT_rn20v25

# Save reaction norm plots ------------------------------------------------
# ggsave("Fig5_dors16length.pdf",plot=p_dors15_rn,device="pdf")

# ggsave("Fig6_plac16length_15-20.pdf",plot=p_plac15_rn15v20,device="pdf")
# ggsave("Fig6_plac16length_15-25.pdf",plot=p_plac15_rn15v25,device="pdf")
# ggsave("Fig6_plac16length_20-25.pdf",plot=p_plac15_rn20v25,device="pdf")

# ggsave("Fig6_placLDT_15-20.pdf",plot=p_placLDT_rn15v20,device="pdf")
# ggsave("Fig6_placLDT_15-25.pdf",plot=p_placLDT_rn15v25,device="pdf")
# ggsave("Fig6_placLDT_20-25.pdf",plot=p_placLDT_rn20v25,device="pdf")

# Fig. 7: Within-family trait correlations across temperatures --------------
# Goal: To estimate and display all the applicable pairwise within-family
# trait correlations across temperatures for each species.

# dilaticollis trait correlations--------------------------------------------
# make a data frame for collecting the correlations and their P-values
# the dummy rows will make nice blank spaces between each trait on the plot 
dila_corrs<-as.data.frame(matrix(nrow = 9*4,ncol = 5))
colnames(dila_corrs)<-c("Species","Trait","Comparison","r","P")
dila_corrs$Species<-"C_dilaticollis"
dila_corrs$Trait<-c("Day 1 L","Day 1 L","Day 1 L","Day 1 L",
                    "Day 15 L","Day 15 L","Day 15 L","Day 15 L",
                    "Pupae L","Pupae L","Pupae L","Pupae L",
                    "Adults L","Adults L","Adults L","Adults L",
                    "Pupae M","Pupae M","Pupae M","Pupae M",
                    "Adults M","Adults M","Adults M","Adults M",
                    "Larvae D","Larvae D","Larvae D","Larvae D",
                    "Pupae D","Pupae D","Pupae D","Pupae D",
                    "Larval S","Larval S","Larval S","Larval S")
dila_corrs$Comparison<-c("15 vs 20","15 vs 25","20 vs 25","dummy")
# the final dummy row will create a weird, extra-large space at that end of the 
# graph, so we'll drop it now.
dila_corrs<-dila_corrs[1:35,]

# dilaticollis Day 1 length
FTP_FT1dila<-FT_prep(larval_set1,"C_dilaticollis","Tamano_dia_1_mm")

any(FTP_FT1dila$N<5)
# No.

# what do the family x temperature correlations look like?
# 15 vs. 20
fill<-cor.test(log(FTP_FT1dila$Tamano_dia_1_mm[FTP_FT1dila$Temp==15]),
               log(FTP_FT1dila$Tamano_dia_1_mm[FTP_FT1dila$Temp==20 & as.numeric(as.character(FTP_FT1dila$Pareja))<16]))
# sample size
length(log(FTP_FT1dila$Tamano_dia_1_mm[FTP_FT1dila$Temp==15]))

dila_corrs[1,4]<-fill$estimate
dila_corrs[1,5]<-fill$p.value

# 15 vs. 25
fill<-cor.test(log(FTP_FT1dila$Tamano_dia_1_mm[FTP_FT1dila$Temp==15]),
               log(FTP_FT1dila$Tamano_dia_1_mm[FTP_FT1dila$Temp==25 & as.numeric(as.character(FTP_FT1dila$Pareja))<16]))
# sample size
length(log(FTP_FT1dila$Tamano_dia_1_mm[FTP_FT1dila$Temp==15]))

dila_corrs[2,4]<-fill$estimate
dila_corrs[2,5]<-fill$p.value

# 20 vs. 25
fill<-cor.test(log(FTP_FT1dila$Tamano_dia_1_mm[FTP_FT1dila$Temp==25]),
               log(FTP_FT1dila$Tamano_dia_1_mm[FTP_FT1dila$Temp==20]))
# sample size
length(log(FTP_FT1dila$Tamano_dia_1_mm[FTP_FT1dila$Temp==25]))

dila_corrs[3,4]<-fill$estimate
dila_corrs[3,5]<-fill$p.value

dila_corrs[1:4,]
# all significantly different from 0

# Day 16 length
FTP_FT15dila<-FT_prep(larval_set1,"C_dilaticollis","Tamano_dia_15_mm")

any(FTP_FT15dila$N<5)
# Yes.

# what do the family x temperature correlations look like?
# 15 vs. 20
fill<-cor.test(log(FTP_FT15dila$Tamano_dia_15_mm[FTP_FT15dila$Temp==15]),
               log(FTP_FT15dila$Tamano_dia_15_mm[FTP_FT15dila$Temp==20 & as.numeric(as.character(FTP_FT15dila$Pareja))<16]))
# sample size
length(log(FTP_FT15dila$Tamano_dia_15_mm[FTP_FT15dila$Temp==15]))

dila_corrs[5,4]<-fill$estimate
dila_corrs[5,5]<-fill$p.value

# 15 vs. 25
fill<-cor.test(log(FTP_FT15dila$Tamano_dia_15_mm[FTP_FT15dila$Temp==15]),
               log(FTP_FT15dila$Tamano_dia_15_mm[FTP_FT15dila$Temp==25 & as.numeric(as.character(FTP_FT15dila$Pareja))<16]))
# sample size
length(log(FTP_FT15dila$Tamano_dia_15_mm[FTP_FT15dila$Temp==15]))

dila_corrs[6,4]<-fill$estimate
dila_corrs[6,5]<-fill$p.value

# 20 vs. 25
fill<-cor.test(log(FTP_FT15dila$Tamano_dia_15_mm[FTP_FT15dila$Temp==25]),
               log(FTP_FT15dila$Tamano_dia_15_mm[FTP_FT15dila$Temp==20]))
# sample size
length(log(FTP_FT15dila$Tamano_dia_15_mm[FTP_FT15dila$Temp==25]))

dila_corrs[7,4]<-fill$estimate
dila_corrs[7,5]<-fill$p.value

dila_corrs[5:8,]
# none different from 0

# Pupal length
FTP_FTPdila<-FT_prep(larval_set1,"C_dilaticollis","Tamano_pupa_mm")

any(FTP_FTPdila$N<5)
# No.

# what do the family x temperature correlations look like?
# there are no suitable 15 degree families for pupae, so there are no 
# 15v20 and 15v25 comparisons
dila_corrs[9,4]<-NA
dila_corrs[9,5]<-NA

dila_corrs[10,4]<-NA
dila_corrs[10,5]<-NA

# 20 vs. 25
fill<-cor.test(log(FTP_FTPdila$Tamano_pupa_mm[FTP_FTPdila$Temp==25]),
               log(FTP_FTPdila$Tamano_pupa_mm[FTP_FTPdila$Temp==20]))
# sample size
length(log(FTP_FTPdila$Tamano_pupa_mm[FTP_FTPdila$Temp==25]))

dila_corrs[11,4]<-fill$estimate
dila_corrs[11,5]<-fill$p.value

dila_corrs[9:12,]
# not significantly different from 0

# Adult length
# only four families(!)
FTP_FTAdila<-FT_prep(larval_set1,"C_dilaticollis","Tamano_adulto_mm")

any(FTP_FTAdila$N<5)
# No.

# what do the family x temperature correlations look like?
# there are no suitable 15 degree families for adults, so there are no 
# 15v20 and 15v25 comparisons
dila_corrs[13,4]<-NA
dila_corrs[13,5]<-NA

dila_corrs[14,4]<-NA
dila_corrs[14,5]<-NA

# 20 vs. 25
fill<-cor.test(log(FTP_FTAdila$Tamano_adulto_mm[FTP_FTAdila$Temp==25]),
               log(FTP_FTAdila$Tamano_adulto_mm[FTP_FTAdila$Temp==20]))
# sample size
length(log(FTP_FTAdila$Tamano_adulto_mm[FTP_FTAdila$Temp==25]))

dila_corrs[15,4]<-fill$estimate
dila_corrs[15,5]<-fill$p.value

dila_corrs[13:16,]
# correlation not significantly different from 0 

# Pupal mass
FTP_FPPdila<-FT_prep(larval_set1,"C_dilaticollis","Peso_pupa_g")

any(FTP_FPPdila$N<5)
# No.

# what do the family x temperature correlations look like?
# there are no suitable 15 degree families for pupae, so there are no 
# 15v20 and 15v25 comparisons
dila_corrs[17,4]<-NA
dila_corrs[17,5]<-NA

dila_corrs[18,4]<-NA
dila_corrs[18,5]<-NA

# 20 vs. 25
fill<-cor.test(log(FTP_FPPdila$Peso_pupa_g[FTP_FPPdila$Temp==25]),
               log(FTP_FPPdila$Peso_pupa_g[FTP_FPPdila$Temp==20]))
# sample size
length(log(FTP_FPPdila$Peso_pupa_g[FTP_FPPdila$Temp==25]))

dila_corrs[19,4]<-fill$estimate
dila_corrs[19,5]<-fill$p.value

dila_corrs[17:20,]
# the 20-25 correlation is significantly different from 0

# Adult mass
FTP_FPAdila<-FT_prep(larval_set1,"C_dilaticollis","Peso_adulto_g")

any(FTP_FPAdila$N<5)
# No.

# what do the family x temperature correlations look like?
# there are no suitable 15 degree families for adults, so there are no 
# 15v20 and 15v25 comparisons
dila_corrs[21,4]<-NA
dila_corrs[21,5]<-NA

dila_corrs[22,4]<-NA
dila_corrs[22,5]<-NA

# 20 vs. 25
fill<-cor.test(log(FTP_FPAdila$Peso_adulto_g[FTP_FPAdila$Temp==25]),
               log(FTP_FPAdila$Peso_adulto_g[FTP_FPAdila$Temp==20]))
# sample size
length(log(FTP_FPAdila$Peso_adulto_g[FTP_FPAdila$Temp==25]))

dila_corrs[23,4]<-fill$estimate
dila_corrs[23,5]<-fill$p.value

dila_corrs[21:24,]
# the 20-25 correlation is not significantly different from 0

# larval development time
FTP_FLDTdila<-FT_prep(larval_set1[larval_set1$pupate==1,],
                      "C_dilaticollis","ldays")

any(FTP_FLDTdila$N<5)
# No.

# what do the family x temperature correlations look like?
# there are no suitable 15 degree families for larval development, so there are
# no 15v20 and 15v25 comparisons
dila_corrs[25,4]<-NA
dila_corrs[25,5]<-NA

dila_corrs[26,4]<-NA
dila_corrs[26,5]<-NA

# 20 vs. 25
fill<-cor.test(FTP_FLDTdila$ldays[FTP_FLDTdila$Temp==25],
               FTP_FLDTdila$ldays[FTP_FLDTdila$Temp==20])
# sample size
length(FTP_FLDTdila$ldays[FTP_FLDTdila$Temp==25])

dila_corrs[27,4]<-fill$estimate
dila_corrs[27,5]<-fill$p.value

dila_corrs[25:28,]
# not significantly different from 0

# pupal development time
FTP_FPDTdila<-FT_prep(larval_set1[larval_set1$emerge==1,],
                      "C_dilaticollis","pdays")

any(FTP_FPDTdila$N<5)
# No.

# what do the family x temperature correlations look like?
# there are no suitable 15 degree families for pupal development, so there are
# no 15v20 and 15v25 comparisons
dila_corrs[29,4]<-NA
dila_corrs[29,5]<-NA

dila_corrs[30,4]<-NA
dila_corrs[30,5]<-NA

# 20 vs. 25
fill<-cor.test(FTP_FPDTdila$pdays[FTP_FPDTdila$Temp==25],
               FTP_FPDTdila$pdays[FTP_FPDTdila$Temp==20])
# sample size
length(FTP_FPDTdila$pdays[FTP_FPDTdila$Temp==25])

dila_corrs[31,4]<-fill$estimate
dila_corrs[31,5]<-fill$p.value

dila_corrs[29:32,]
# not significantly different from 0


# larval survival
# to calculate this, we'll just take the average of "pupate" for all the larvae,
# since pupation (survival) = 1 and death equals 0
FTP_FLSdila<-FT_prep(larval_set1,
                     "C_dilaticollis","pupate")

any(FTP_FLSdila$N<5)
# No.

# what do the family x temperature correlations look like?
# 15 vs. 20
fill<-cor.test(asin(FTP_FLSdila$pupate[FTP_FLSdila$Temp==15]),
               asin(FTP_FLSdila$pupate[FTP_FLSdila$Temp==20 & 
                                         as.numeric(as.character(FTP_FLSdila$Pareja))<16]))
# sample size
length(asin(FTP_FLSdila$pupate[FTP_FLSdila$Temp==15]))

dila_corrs[33,4]<-fill$estimate
dila_corrs[33,5]<-fill$p.value

# 15 vs. 25
fill<-cor.test(asin(FTP_FLSdila$pupate[FTP_FLSdila$Temp==15]),
               asin(FTP_FLSdila$pupate[FTP_FLSdila$Temp==25 & as.numeric(as.character(FTP_FLSdila$Pareja))<16]))
# sample size
length(asin(FTP_FLSdila$pupate[FTP_FLSdila$Temp==15]))

dila_corrs[34,4]<-fill$estimate
dila_corrs[34,5]<-fill$p.value

# 20 vs. 25
fill<-cor.test(asin(FTP_FLSdila$pupate[FTP_FLSdila$Temp==25]),
               asin(FTP_FLSdila$pupate[FTP_FLSdila$Temp==20]))
# sample size
length(asin(FTP_FLSdila$pupate[FTP_FLSdila$Temp==25]))

dila_corrs[35,4]<-fill$estimate
dila_corrs[35,5]<-fill$p.value

dila_corrs[33:35,]
# the 20-25 comparison is significantly different from 0


# dilaticollis correlations plot (Fig. 7A)----------------------------------
# to display them properly, we to label them individually with a single label
# (rather than in subset combinations of trait x comparison)
dila_corrs[,6]<-paste(dila_corrs[,2],dila_corrs[,3])
colnames(dila_corrs)[6]<-"Trait_comparison"
dila_corrs$Trait_comparison<-as.factor(dila_corrs$Trait_comparison)
dila_corrs$Trait_comparison<-factor(dila_corrs$Trait_comparison,levels=rev(dila_corrs$Trait_comparison))

p_dila_corrs<-ggplot(data=dila_corrs, aes(x=Trait_comparison, y=r))+
  geom_col(position="dodge") + theme_classic2() + coord_flip() +
  theme(axis.text.y= element_text(size=7), 
        axis.title.y=element_text(angle=0, vjust=0.5), 
        panel.border = element_rect(colour = "black", fill=NA)) +
  geom_hline(yintercept=0) +
  expand_limits(y=c(-1,1)) +
  scale_x_discrete(expand=c(0.02,0.1) 
                   # note: this next line makes the plot really hard to
                   # interpret, but in the final figure, the traits are going to
                   # be be indicated using brackets and sub-brackets. It's
                   # easier to take out the traits with one line of code here
                   # than by changing each by hand in Illustrator. But be sure
                   # to do the labels while looking at a reference plot that has
                   # the full labels!
                   , labels=rev(dila_corrs$Comparison)
  ) +
  scale_y_continuous(breaks=pretty_breaks(5)) + 
  labs(x= "Trait", y= "r", fill="Comparison", title= "dilaticollis pairwise correlations", family="sans")

p_dila_corrs

dila_corrs

# dorsalis trait correlations -----------------------------------------------------
# make a data frame for collecting the correlations and their P-values
# the dummy rows will make nice blank spaces between each trait on the plot 
dors_corrs<-as.data.frame(matrix(nrow = 9*4,ncol = 5))
colnames(dors_corrs)<-c("Species","Trait","Comparison","r","P")
dors_corrs$Species<-"C_dorsalis"
dors_corrs$Trait<-c("Day 1 L","Day 1 L","Day 1 L","Day 1 L",
                    "Day 15 L","Day 15 L","Day 15 L","Day 15 L",
                    "Pupae L","Pupae L","Pupae L","Pupae L",
                    "Adults L","Adults L","Adults L","Adults L",
                    "Pupae M","Pupae M","Pupae M","Pupae M",
                    "Adults M","Adults M","Adults M","Adults M",
                    "Larvae D","Larvae D","Larvae D","Larvae D",
                    "Pupae D","Pupae D","Pupae D","Pupae D",
                    "Larval S","Larval S","Larval S","Larval S")
dors_corrs$Comparison<-c("15 vs 20","15 vs 25","20 vs 25","dummy")
# the final dummy row will create a weird, extra-large space at that end of the
# graph, so we'll drop it now.
dors_corrs<-dors_corrs[1:35,]

# Day 1 length
FTP_FT1dors<-FT_prep(larval_set1,"C_dorsalis","Tamano_dia_1_mm")

any(FTP_FT1dors$N<5)
# Yes.Some of the families have fewer than 5 in 15C
# the list to include for the 15C correlations is:
as.numeric(as.character(FTP_FT1dors$Pareja[FTP_FT1dors$Temp==15 & FTP_FT1dors$N>=5]))

# what do the family x temperature correlations look like?
# 15 vs. 20
fill<-cor.test(log(FTP_FT1dors$Tamano_dia_1_mm[FTP_FT1dors$Temp==15 &
                                                 as.character(FTP_FT1dors$Pareja) %in%
                                                 as.character(FTP_FT1dors$Pareja[FTP_FT1dors$Temp==15 & FTP_FT1dors$N>=5])]),
               log(FTP_FT1dors$Tamano_dia_1_mm[FTP_FT1dors$Temp==20 & as.character(FTP_FT1dors$Pareja) %in%
                                                 as.character(FTP_FT1dors$Pareja[FTP_FT1dors$Temp==15 & FTP_FT1dors$N>=5])]))
# sample size
length(log(FTP_FT1dors$Tamano_dia_1_mm[FTP_FT1dors$Temp==15 &
                                         as.character(FTP_FT1dors$Pareja) %in%
                                         as.character(FTP_FT1dors$Pareja[FTP_FT1dors$Temp==15 & FTP_FT1dors$N>=5])]))


dors_corrs[1,4]<-fill$estimate
dors_corrs[1,5]<-fill$p.value

# 15 vs. 25
fill<-cor.test(log(FTP_FT1dors$Tamano_dia_1_mm[FTP_FT1dors$Temp==15 &
                                                 as.character(FTP_FT1dors$Pareja) %in%
                                                 as.character(FTP_FT1dors$Pareja[FTP_FT1dors$Temp==15 & FTP_FT1dors$N>=5])]),
               log(FTP_FT1dors$Tamano_dia_1_mm[FTP_FT1dors$Temp==25 &
                                                 as.character(FTP_FT1dors$Pareja) %in%
                                                 as.character(FTP_FT1dors$Pareja[FTP_FT1dors$Temp==15 & FTP_FT1dors$N>=5])]))
# sample size
length(log(FTP_FT1dors$Tamano_dia_1_mm[FTP_FT1dors$Temp==15 &
                                         as.character(FTP_FT1dors$Pareja) %in%
                                         as.character(FTP_FT1dors$Pareja[FTP_FT1dors$Temp==15 & FTP_FT1dors$N>=5])]))

dors_corrs[2,4]<-fill$estimate
dors_corrs[2,5]<-fill$p.value

# 20 vs. 25
fill<-cor.test(log(FTP_FT1dors$Tamano_dia_1_mm[FTP_FT1dors$Temp==25]),
               log(FTP_FT1dors$Tamano_dia_1_mm[FTP_FT1dors$Temp==20]))
# sample size
length(log(FTP_FT1dors$Tamano_dia_1_mm[FTP_FT1dors$Temp==25]))

dors_corrs[3,4]<-fill$estimate
dors_corrs[3,5]<-fill$p.value

dors_corrs[1:4,]
# only the 20-25 comparison is significantly different from 0

# Day 16 length
FTP_FT15dors<-FT_prep(larval_set1,"C_dorsalis","Tamano_dia_15_mm")

any(FTP_FT15dors$N<5)
# No.

# what do the family x temperature correlations look like?
# not enough individuals for the 15C comparisons
dors_corrs[5,4]<-NA
dors_corrs[5,5]<-NA

dors_corrs[6,4]<-NA
dors_corrs[6,5]<-NA

# 20 vs. 25
fill<-cor.test(log(FTP_FT15dors$Tamano_dia_15_mm[FTP_FT15dors$Temp==25]),
               log(FTP_FT15dors$Tamano_dia_15_mm[FTP_FT15dors$Temp==20]))
# sample size
length(log(FTP_FT15dors$Tamano_dia_15_mm[FTP_FT15dors$Temp==25]))

dors_corrs[7,4]<-fill$estimate
dors_corrs[7,5]<-fill$p.value

dors_corrs[5:8,]
# not significantly different from 0

# Pupal length
FTP_FTPdors<-FT_prep(larval_set1,"C_dorsalis","Tamano_pupa_mm")

any(FTP_FTPdors$N<5)
# No.

# what do the family x temperature correlations look like?
# there are no suitable 15 degree families for pupae, so there are no 
# 15v20 and 15v25 comparisons
dors_corrs[9,4]<-NA
dors_corrs[9,5]<-NA

dors_corrs[10,4]<-NA
dors_corrs[10,5]<-NA

# 20 vs. 25
fill<-cor.test(log(FTP_FTPdors$Tamano_pupa_mm[FTP_FTPdors$Temp==25]),
               log(FTP_FTPdors$Tamano_pupa_mm[FTP_FTPdors$Temp==20]))
# sample size
length(log(FTP_FTPdors$Tamano_pupa_mm[FTP_FTPdors$Temp==25]))

dors_corrs[11,4]<-fill$estimate
dors_corrs[11,5]<-fill$p.value

dors_corrs[9:12,]
# not significantly different from 01

# Adult length
FTP_FTAdors<-FT_prep(larval_set1,"C_dorsalis","Tamano_adulto_mm")

any(FTP_FTAdors$N<5)
# No.

# what do the family x temperature correlations look like?
# there are no suitable 15 degree families for adults, so there are no 
# 15v20 and 15v25 comparisons
dors_corrs[13,4]<-NA
dors_corrs[13,5]<-NA

dors_corrs[14,4]<-NA
dors_corrs[14,5]<-NA

# 20 vs. 25
fill<-cor.test(log(FTP_FTAdors$Tamano_adulto_mm[FTP_FTAdors$Temp==25]),
               log(FTP_FTAdors$Tamano_adulto_mm[FTP_FTAdors$Temp==20]))
# sample size
length(log(FTP_FTAdors$Tamano_adulto_mm[FTP_FTAdors$Temp==25]))

dors_corrs[15,4]<-fill$estimate
dors_corrs[15,5]<-fill$p.value

dors_corrs[13:16,]
# correlation borderline similar to 0

# Pupal mass
FTP_FPPdors<-FT_prep(larval_set1,"C_dorsalis","Peso_pupa_g")

any(FTP_FPPdors$N<5)
# No.

# what do the family x temperature correlations look like?
# there are no suitable 15 degree families for pupae, so there are no 
# 15v20 and 15v25 comparisons
dors_corrs[17,4]<-NA
dors_corrs[17,5]<-NA

dors_corrs[18,4]<-NA
dors_corrs[18,5]<-NA

# 20 vs. 25
fill<-cor.test(log(FTP_FPPdors$Peso_pupa_g[FTP_FPPdors$Temp==25]),
               log(FTP_FPPdors$Peso_pupa_g[FTP_FPPdors$Temp==20]))
# sample size
length(log(FTP_FPPdors$Peso_pupa_g[FTP_FPPdors$Temp==25]))

dors_corrs[19,4]<-fill$estimate
dors_corrs[19,5]<-fill$p.value

dors_corrs[17:20,]
# not significantly different from 0

# Adult mass
FTP_FPAdors<-FT_prep(larval_set1,"C_dorsalis","Peso_adulto_g")

any(FTP_FPAdors$N<5)
# No.

# what do the family x temperature correlations look like?
# there are no suitable 15 degree families for adults, so there are no 
# 15v20 and 15v25 comparisons
dors_corrs[21,4]<-NA
dors_corrs[21,5]<-NA

dors_corrs[22,4]<-NA
dors_corrs[22,5]<-NA

# 20 vs. 25
fill<-cor.test(log(FTP_FPAdors$Peso_adulto_g[FTP_FPAdors$Temp==25]),
               log(FTP_FPAdors$Peso_adulto_g[FTP_FPAdors$Temp==20]))
# sample size
length(log(FTP_FPAdors$Peso_adulto_g[FTP_FPAdors$Temp==25]))

dors_corrs[23,4]<-fill$estimate
dors_corrs[23,5]<-fill$p.value

dors_corrs[21:24,]
# the 20-25 correlation is not significantly different from 0

# larval development time
FTP_FLDTdors<-FT_prep(larval_set1[larval_set1$pupate==1,],
                      "C_dorsalis","ldays")

any(FTP_FLDTdors$N<5)
# No.

# what do the family x temperature correlations look like?
# there are no suitable 15 degree families for larval development, so there are
# no 15v20 and 15v25 comparisons
dors_corrs[25,4]<-NA
dors_corrs[25,5]<-NA

dors_corrs[26,4]<-NA
dors_corrs[26,5]<-NA

# 20 vs. 25
fill<-cor.test(FTP_FLDTdors$ldays[FTP_FLDTdors$Temp==25],
               FTP_FLDTdors$ldays[FTP_FLDTdors$Temp==20])
# sample size
length(FTP_FLDTdors$ldays[FTP_FLDTdors$Temp==25])

dors_corrs[27,4]<-fill$estimate
dors_corrs[27,5]<-fill$p.value

dors_corrs[25:28,]
# significantly different from

# pupal development time
FTP_FPDTdors<-FT_prep(larval_set1[larval_set1$emerge==1,],
                      "C_dorsalis","pdays")

any(FTP_FPDTdors$N<5)
# No.

# what do the family x temperature correlations look like?
# there are no suitable 15 degree families for pupal development, so there are
# no 15v20 and 15v25 comparisons
dors_corrs[29,4]<-NA
dors_corrs[29,5]<-NA

dors_corrs[30,4]<-NA
dors_corrs[30,5]<-NA

# 20 vs. 25
fill<-cor.test(FTP_FPDTdors$pdays[FTP_FPDTdors$Temp==25],
               FTP_FPDTdors$pdays[FTP_FPDTdors$Temp==20])
# sample size
length(FTP_FPDTdors$pdays[FTP_FPDTdors$Temp==25])

dors_corrs[31,4]<-fill$estimate
dors_corrs[31,5]<-fill$p.value

dors_corrs[29:32,]
# not significantly different from 0


# larval survival
# to calculate this, we'll just take the average of "pupate" for all the larvae,
# since pupation (survival) = 1 and death equals 0
FTP_FLSdors<-FT_prep(larval_set1,
                     "C_dorsalis","pupate")

any(FTP_FLSdors$N<5)
# Yes.

# what do the family x temperature correlations look like?
# because some of the families have N<5 at 15C, we need a list of acceptable
# families for the 20C data
match_list<-as.numeric(as.character(FTP_FLSdors$Pareja[FTP_FLSdors$Temp==15 & FTP_FLSdors$N>=5]))

# 15 vs. 20
fill<-cor.test(asin(FTP_FLSdors$pupate[FTP_FLSdors$Temp==15 & FTP_FLSdors$N>=5]),
               asin(FTP_FLSdors$pupate[FTP_FLSdors$Temp==20 & as.numeric(as.character(FTP_FLSdors$Pareja)) %in% match_list]))
# sample size
length(asin(FTP_FLSdors$pupate[FTP_FLSdors$Temp==15 & FTP_FLSdors$N>=5]))

dors_corrs[33,4]<-fill$estimate
dors_corrs[33,5]<-fill$p.value

# 15 vs. 25
fill<-cor.test(asin(FTP_FLSdors$pupate[FTP_FLSdors$Temp==15 & FTP_FLSdors$N>=5]),
               asin(FTP_FLSdors$pupate[FTP_FLSdors$Temp==25 & as.numeric(as.character(FTP_FLSdors$Pareja))%in% match_list]))
# sample size
length(asin(FTP_FLSdors$pupate[FTP_FLSdors$Temp==15 & FTP_FLSdors$N>=5]))

dors_corrs[34,4]<-fill$estimate
dors_corrs[34,5]<-fill$p.value

# 20 vs. 25
fill<-cor.test(asin(FTP_FLSdors$pupate[FTP_FLSdors$Temp==25]),
               asin(FTP_FLSdors$pupate[FTP_FLSdors$Temp==20]))
# sample size
length(asin(FTP_FLSdors$pupate[FTP_FLSdors$Temp==25]))

dors_corrs[35,4]<-fill$estimate
dors_corrs[35,5]<-fill$p.value

dors_corrs[33:35,]
# none significantly different from 0


# dorsalis correlations plot (Fig. 7B)---------------------------------------
# to display them properly, we to label them individually with a single label
# (rather than in subset combinations of trait x comparison)
dors_corrs[,6]<-paste(dors_corrs[,2],dors_corrs[,3])
colnames(dors_corrs)[6]<-"Trait_comparison"
dors_corrs$Trait_comparison<-as.factor(dors_corrs$Trait_comparison)
dors_corrs$Trait_comparison<-factor(dors_corrs$Trait_comparison,levels=rev(dors_corrs$Trait_comparison))

p_dors_corrs<-ggplot(data=dors_corrs, aes(x=Trait_comparison, y=r))+
  geom_col(position="dodge") + theme_classic2() + coord_flip() +
  theme(axis.text.y= element_text(size=7), 
        axis.title.y=element_text(angle=0, vjust=0.5), 
        panel.border = element_rect(colour = "black", fill=NA)) +
  geom_hline(yintercept=0) +
  expand_limits(y=c(-1,1)) +
  scale_x_discrete(expand=c(0.02,0.1) 
                   # note: this next line makes the plot really hard to
                   # interpret, but in the final graph, the traits are going to
                   # be be indicated using brackets and sub-brackets. It's
                   # easier to take out the traits with one line of code here
                   # than by changing each by hand in Illustrator. But be sure
                   # to do the labels while looking at a reference plot that has
                   # the full labels!
                   , labels=rev(dors_corrs$Comparison)
  ) +
  scale_y_continuous(breaks=pretty_breaks(5)) + 
  labs(x= "Trait", y= "r", fill="Comparison", title= "dorsalis pairwise correlations", family="sans")

p_dors_corrs

dors_corrs

# placida trait correlations -----------------------------------------------------
# make a data frame for collecting the correlations and their P-values
# the dummy rows will make nice blank spaces between each trait on the plot 
plac_corrs<-as.data.frame(matrix(nrow = 9*4,ncol = 5))
colnames(plac_corrs)<-c("Species","Trait","Comparison","r","P")
plac_corrs$Species<-"C_placida"
plac_corrs$Trait<-c("Day 1 L","Day 1 L","Day 1 L","Day 1 L",
                    "Day 15 L","Day 15 L","Day 15 L","Day 15 L",
                    "Pupae L","Pupae L","Pupae L","Pupae L",
                    "Adults L","Adults L","Adults L","Adults L",
                    "Pupae M","Pupae M","Pupae M","Pupae M",
                    "Adults M","Adults M","Adults M","Adults M",
                    "Larvae D","Larvae D","Larvae D","Larvae D",
                    "Pupae D","Pupae D","Pupae D","Pupae D",
                    "Larval S","Larval S","Larval S","Larval S")
plac_corrs$Comparison<-c("15 vs 20","15 vs 25","20 vs 25","dummy")
# the final dummy row will create a weird, extra-large space at that end of the 
# graph, so we'll drop it now.
plac_corrs<-plac_corrs[1:35,]

# Day 1 length
FTP_FT1plac<-FT_prep(larval_set1,"C_placida","Tamano_dia_1_mm")

any(FTP_FT1plac$N<5)
# No.


# what do the family x temperature correlations look like?
# 15 vs. 20
fill<-cor.test(log(FTP_FT1plac$Tamano_dia_1_mm[FTP_FT1plac$Temp==15]),
               log(FTP_FT1plac$Tamano_dia_1_mm[FTP_FT1plac$Temp==20]))
# sample size
length(log(FTP_FT1plac$Tamano_dia_1_mm[FTP_FT1plac$Temp==15]))

plac_corrs[1,4]<-fill$estimate
plac_corrs[1,5]<-fill$p.value

# 15 vs. 25
fill<-cor.test(log(FTP_FT1plac$Tamano_dia_1_mm[FTP_FT1plac$Temp==15]),
               log(FTP_FT1plac$Tamano_dia_1_mm[FTP_FT1plac$Temp==25]))
# sample size
length(log(FTP_FT1plac$Tamano_dia_1_mm[FTP_FT1plac$Temp==15]))

plac_corrs[2,4]<-fill$estimate
plac_corrs[2,5]<-fill$p.value

# 20 vs. 25
fill<-cor.test(log(FTP_FT1plac$Tamano_dia_1_mm[FTP_FT1plac$Temp==25]),
               log(FTP_FT1plac$Tamano_dia_1_mm[FTP_FT1plac$Temp==20]))
# sample size
length(log(FTP_FT1plac$Tamano_dia_1_mm[FTP_FT1plac$Temp==25]))

plac_corrs[3,4]<-fill$estimate
plac_corrs[3,5]<-fill$p.value

plac_corrs[1:4,]
# all significantly different from 0

# Day 16 length
FTP_FT15plac<-FT_prep(larval_set1,"C_placida","Tamano_dia_15_mm")

any(FTP_FT15plac$N<5)
# No.

# what do the family x temperature correlations look like?
# 15 vs. 20
fill<-cor.test(log(FTP_FT15plac$Tamano_dia_15_mm[FTP_FT15plac$Temp==15]),
               log(FTP_FT15plac$Tamano_dia_15_mm[FTP_FT15plac$Temp==20]))
# sample size
length(log(FTP_FT15plac$Tamano_dia_15_mm[FTP_FT15plac$Temp==15]))

plac_corrs[5,4]<-fill$estimate
plac_corrs[5,5]<-fill$p.value

# 15 vs. 25
fill<-cor.test(log(FTP_FT15plac$Tamano_dia_15_mm[FTP_FT15plac$Temp==15]),
               log(FTP_FT15plac$Tamano_dia_15_mm[FTP_FT15plac$Temp==25 ]))
# sample size
length(log(FTP_FT15plac$Tamano_dia_15_mm[FTP_FT15plac$Temp==15]))

plac_corrs[6,4]<-fill$estimate
plac_corrs[6,5]<-fill$p.value

# 20 vs. 25
fill<-cor.test(log(FTP_FT15plac$Tamano_dia_15_mm[FTP_FT15plac$Temp==25]),
               log(FTP_FT15plac$Tamano_dia_15_mm[FTP_FT15plac$Temp==20]))
# sample size
length(log(FTP_FT15plac$Tamano_dia_15_mm[FTP_FT15plac$Temp==25]))

plac_corrs[7,4]<-fill$estimate
plac_corrs[7,5]<-fill$p.value

plac_corrs[5:8,]
# all significantly different from 0: 15-20 is positive and 15-25 and 20-25 are
# NEGATIVE


# Pupal length
FTP_FTPplac<-FT_prep(larval_set1,"C_placida","Tamano_pupa_mm")

any(FTP_FTPplac$N<5)
# Yes.

# what do the family x temperature correlations look like?
# 15 vs. 20
fill<-cor.test(log(FTP_FTPplac$Tamano_pupa_mm[FTP_FTPplac$Temp==15 &
                                                FTP_FTPplac$N>=5 &
                                                as.numeric(as.character(FTP_FTPplac$Pareja))!=1]),
               log(FTP_FTPplac$Tamano_pupa_mm[FTP_FTPplac$Temp==20 &
                                                as.numeric(as.character(FTP_FTPplac$Pareja))!=1 &
                                                as.numeric(as.character(FTP_FTPplac$Pareja))!=8]))
# sample size
length(log(FTP_FTPplac$Tamano_pupa_mm[FTP_FTPplac$Temp==15 &
                                        FTP_FTPplac$N>=5 &
                                        as.numeric(as.character(FTP_FTPplac$Pareja))!=1]))

plac_corrs[9,4]<-fill$estimate
plac_corrs[9,5]<-fill$p.value

# 15 vs. 25
fill<-cor.test(log(FTP_FTPplac$Tamano_pupa_mm[FTP_FTPplac$Temp==15 &
                                                FTP_FTPplac$N>=5 &
                                                as.numeric(as.character(FTP_FTPplac$Pareja))!=1]),
               log(FTP_FTPplac$Tamano_pupa_mm[FTP_FTPplac$Temp==25 &
                                                as.numeric(as.character(FTP_FTPplac$Pareja))!=1 &
                                                as.numeric(as.character(FTP_FTPplac$Pareja))!=8]))

# sample size
length(log(FTP_FTPplac$Tamano_pupa_mm[FTP_FTPplac$Temp==15 &
                                        FTP_FTPplac$N>=5 &
                                        as.numeric(as.character(FTP_FTPplac$Pareja))!=1]))

plac_corrs[10,4]<-fill$estimate
plac_corrs[10,5]<-fill$p.value

# 20 vs. 25
fill<-cor.test(log(FTP_FTPplac$Tamano_pupa_mm[FTP_FTPplac$Temp==25]),
               log(FTP_FTPplac$Tamano_pupa_mm[FTP_FTPplac$Temp==20]))
# sample size
length(log(FTP_FTPplac$Tamano_pupa_mm[FTP_FTPplac$Temp==25]))

plac_corrs[11,4]<-fill$estimate
plac_corrs[11,5]<-fill$p.value

plac_corrs[9:12,]
# not significantly different from 0

# Adult length
FTP_FTAplac<-FT_prep(larval_set1,"C_placida","Tamano_adulto_mm")

any(FTP_FTAplac$N<5)
# No.
# Note (Jan 3, 2019): There are three families that have N>=5 for 15C, but
# that's really not enough points to do a decent correlation, in my opinion.

# what do the family x temperature correlations look like?
# there are not enough 15 degree families for adults, so there are no 
# 15v20 and 15v25 comparisons
plac_corrs[13,4]<-NA
plac_corrs[13,5]<-NA

plac_corrs[14,4]<-NA
plac_corrs[14,5]<-NA

# 20 vs. 25
fill<-cor.test(log(FTP_FTAplac$Tamano_adulto_mm[FTP_FTAplac$Temp==25]),
               log(FTP_FTAplac$Tamano_adulto_mm[FTP_FTAplac$Temp==20]))
# sample size
length(log(FTP_FTAplac$Tamano_adulto_mm[FTP_FTAplac$Temp==25]))

plac_corrs[15,4]<-fill$estimate
plac_corrs[15,5]<-fill$p.value

plac_corrs[13:16,]
# correlation not significantly different from 0

# Pupal mass
FTP_FPPplac<-FT_prep(larval_set1,"C_placida","Peso_pupa_g")
# need to omit Pareja 1 from 15C comparisons 

any(FTP_FPPplac$N<5)
FTP_FPPplac
# Yes. Pareja 8 only has 2 individuals in 15C and needs to be omitted from 15C
# correlations.
# the list of families to include in the 15C correlations is:
as.character(FTP_FPPplac$Pareja[FTP_FPPplac$Temp==15 & FTP_FPPplac$N>=5])

# what do the family x temperature correlations look like?
# 15 vs. 20
fill<-cor.test(log(FTP_FPPplac$Peso_pupa_g[FTP_FPPplac$Temp==15 
                                           & as.character(FTP_FPPplac$Pareja) %in% 
                                             as.character(FTP_FPPplac$Pareja[FTP_FPPplac$Temp==15 & FTP_FPPplac$N>=5])]),
               log(FTP_FPPplac$Peso_pupa_g[FTP_FPPplac$Temp==20 &
                                             as.character(FTP_FPPplac$Pareja) %in% 
                                             as.character(FTP_FPPplac$Pareja[FTP_FPPplac$Temp==15 & FTP_FPPplac$N>=5])]))
# sample size
length(log(FTP_FPPplac$Peso_pupa_g[FTP_FPPplac$Temp==15 
                                   & as.character(FTP_FPPplac$Pareja) %in% 
                                     as.character(FTP_FPPplac$Pareja[FTP_FPPplac$Temp==15 & FTP_FPPplac$N>=5])]))

plac_corrs[17,4]<-fill$estimate
plac_corrs[17,5]<-fill$p.value

# 15 vs. 25
fill<-cor.test(log(FTP_FPPplac$Peso_pupa_g[FTP_FPPplac$Temp==15
                                           & as.character(FTP_FPPplac$Pareja) %in% 
                                             as.character(FTP_FPPplac$Pareja[FTP_FPPplac$Temp==15 & FTP_FPPplac$N>=5])]),
               log(FTP_FPPplac$Peso_pupa_g[FTP_FPPplac$Temp==25  
                                           & as.character(FTP_FPPplac$Pareja) %in% 
                                             as.character(FTP_FPPplac$Pareja[FTP_FPPplac$Temp==15 & FTP_FPPplac$N>=5])]))
# sample size
length(log(FTP_FPPplac$Peso_pupa_g[FTP_FPPplac$Temp==15
                                   & as.character(FTP_FPPplac$Pareja) %in% 
                                     as.character(FTP_FPPplac$Pareja[FTP_FPPplac$Temp==15 & FTP_FPPplac$N>=5])]))

plac_corrs[18,4]<-fill$estimate
plac_corrs[18,5]<-fill$p.value

# 20 vs. 25
fill<-cor.test(log(FTP_FPPplac$Peso_pupa_g[FTP_FPPplac$Temp==25]),
               log(FTP_FPPplac$Peso_pupa_g[FTP_FPPplac$Temp==20]))
# sample size
length(log(FTP_FPPplac$Peso_pupa_g[FTP_FPPplac$Temp==25]))

plac_corrs[19,4]<-fill$estimate
plac_corrs[19,5]<-fill$p.value

plac_corrs[17:20,]
# all correlations not significantly different from 0

# Adult mass
FTP_FPAplac<-FT_prep(larval_set1,"C_placida","Peso_adulto_g")

any(FTP_FPAplac$N<5)
# No.

# what do the family x temperature correlations look like?
# there are not enough suitable 15 degree families for adults, so there are no 
# 15v20 and 15v25 comparisons
plac_corrs[21,4]<-NA
plac_corrs[21,5]<-NA

plac_corrs[22,4]<-NA
plac_corrs[22,5]<-NA

# 20 vs. 25
fill<-cor.test(log(FTP_FPAplac$Peso_adulto_g[FTP_FPAplac$Temp==25]),
               log(FTP_FPAplac$Peso_adulto_g[FTP_FPAplac$Temp==20]))
# sample size
length(log(FTP_FPAplac$Peso_adulto_g[FTP_FPAplac$Temp==25]))

plac_corrs[23,4]<-fill$estimate
plac_corrs[23,5]<-fill$p.value

plac_corrs[21:24,]
# the 20-25 correlation is not significantly different from 0

# larval development time
FTP_FLDTplac<-FT_prep(larval_set1[larval_set1$pupate==1,],
                      "C_placida","ldays")
# note that Pareja 1 has no 15C group

any(FTP_FLDTplac$N<5)
# Yes.Pareja #8 only has two individuals in 15C. Pareja #1 must also be omitted
# from 15C correlations.
# The list of families to include:
as.character(FTP_FLDTplac$Pareja[FTP_FLDTplac$Temp==15 & FTP_FLDTplac$N>=5])

# what do the family x temperature correlations look like?
# 15 vs. 20
fill<-cor.test(FTP_FLDTplac$ldays[FTP_FLDTplac$Temp==15
                                  & as.character(FTP_FLDTplac$Pareja) %in% 
                                    as.character(FTP_FLDTplac$Pareja[FTP_FLDTplac$Temp==15 & FTP_FLDTplac$N>=5])],
               FTP_FLDTplac$ldays[FTP_FLDTplac$Temp==20
                                  & as.character(FTP_FLDTplac$Pareja) %in% 
                                    as.character(FTP_FLDTplac$Pareja[FTP_FLDTplac$Temp==15 & FTP_FLDTplac$N>=5])])
# sample size
length(FTP_FLDTplac$ldays[FTP_FLDTplac$Temp==15
                          & as.character(FTP_FLDTplac$Pareja) %in% 
                            as.character(FTP_FLDTplac$Pareja[FTP_FLDTplac$Temp==15 & FTP_FLDTplac$N>=5])])

plac_corrs[25,4]<-fill$estimate
plac_corrs[25,5]<-fill$p.value

# 15 vs. 25
fill<-cor.test(FTP_FLDTplac$ldays[FTP_FLDTplac$Temp==15
                                  & as.character(FTP_FLDTplac$Pareja) %in% 
                                    as.character(FTP_FLDTplac$Pareja[FTP_FLDTplac$Temp==15 & FTP_FLDTplac$N>=5])],
               FTP_FLDTplac$ldays[FTP_FLDTplac$Temp==25
                                  & as.character(FTP_FLDTplac$Pareja) %in% 
                                    as.character(FTP_FLDTplac$Pareja[FTP_FLDTplac$Temp==15 & FTP_FLDTplac$N>=5])])
# sample size
length(FTP_FLDTplac$ldays[FTP_FLDTplac$Temp==15
                          & as.character(FTP_FLDTplac$Pareja) %in% 
                            as.character(FTP_FLDTplac$Pareja[FTP_FLDTplac$Temp==15 & FTP_FLDTplac$N>=5])])

plac_corrs[26,4]<-fill$estimate
plac_corrs[26,5]<-fill$p.value

# 20 vs. 25
fill<-cor.test(FTP_FLDTplac$ldays[FTP_FLDTplac$Temp==25],
               FTP_FLDTplac$ldays[FTP_FLDTplac$Temp==20])
# sample size
length(FTP_FLDTplac$ldays[FTP_FLDTplac$Temp==25])

plac_corrs[27,4]<-fill$estimate
plac_corrs[27,5]<-fill$p.value

plac_corrs[25:28,]
# 15v20 is significantly different from 0

# pupal development time
FTP_FPDTplac<-FT_prep(larval_set1[larval_set1$emerge==1,],
                      "C_placida","pdays")
# again, Pareja 1 has no 15C data

any(FTP_FPDTplac$N<5)
# Yes. Again, we need to omit families 1 and 8 from the 15C analysis, along with
# families 4 & 5:
as.character(FTP_FPDTplac$Pareja[FTP_FPDTplac$Temp==15 & FTP_FPDTplac$N>=5])

# what do the family x temperature correlations look like?
# 15 vs. 20
fill<-cor.test(FTP_FPDTplac$pdays[FTP_FPDTplac$Temp==15
                                  & as.character(FTP_FPDTplac$Pareja) %in% 
                                    as.character(FTP_FPDTplac$Pareja[FTP_FPDTplac$Temp==15 & FTP_FPDTplac$N>=5])],
               FTP_FPDTplac$pdays[FTP_FPDTplac$Temp==20
                                  & as.character(FTP_FPDTplac$Pareja) %in% 
                                    as.character(FTP_FPDTplac$Pareja[FTP_FPDTplac$Temp==15 & FTP_FPDTplac$N>=5])])
# sample size
length(FTP_FPDTplac$pdays[FTP_FPDTplac$Temp==15
                          & as.character(FTP_FPDTplac$Pareja) %in% 
                            as.character(FTP_FPDTplac$Pareja[FTP_FPDTplac$Temp==15 & FTP_FPDTplac$N>=5])])

plac_corrs[29,4]<-fill$estimate
plac_corrs[29,5]<-fill$p.value

# 15 vs. 25
fill<-cor.test(FTP_FPDTplac$pdays[FTP_FPDTplac$Temp==15
                                  & as.character(FTP_FPDTplac$Pareja) %in% 
                                    as.character(FTP_FPDTplac$Pareja[FTP_FPDTplac$Temp==15 & FTP_FPDTplac$N>=5])],
               FTP_FPDTplac$pdays[FTP_FPDTplac$Temp==25
                                  & as.character(FTP_FPDTplac$Pareja) %in% 
                                    as.character(FTP_FPDTplac$Pareja[FTP_FPDTplac$Temp==15 & FTP_FPDTplac$N>=5])])
# sample size
length(FTP_FPDTplac$pdays[FTP_FPDTplac$Temp==15
                          & as.character(FTP_FPDTplac$Pareja) %in% 
                            as.character(FTP_FPDTplac$Pareja[FTP_FPDTplac$Temp==15 & FTP_FPDTplac$N>=5])])

plac_corrs[30,4]<-fill$estimate
plac_corrs[30,5]<-fill$p.value

# 20 vs. 25
fill<-cor.test(FTP_FPDTplac$pdays[FTP_FPDTplac$Temp==25],
               FTP_FPDTplac$pdays[FTP_FPDTplac$Temp==20])
# sample size
length(FTP_FPDTplac$pdays[FTP_FPDTplac$Temp==25])

plac_corrs[31,4]<-fill$estimate
plac_corrs[31,5]<-fill$p.value

plac_corrs[29:32,]
# none significantly different from 0

# larval survival
# to calculate this, we'll just take the average of "pupate" for all the larvae,
# since pupation (survival) = 1 and death equals 0
FTP_FLSplac<-FT_prep(larval_set1,
                     "C_placida","pupate")

any(FTP_FLSplac$N<5)
# No.

# what do the family x temperature correlations look like?
# 15 vs. 20
fill<-cor.test(asin(FTP_FLSplac$pupate[FTP_FLSplac$Temp==15]),
               asin(FTP_FLSplac$pupate[FTP_FLSplac$Temp==20]))
# sample size
length(asin(FTP_FLSplac$pupate[FTP_FLSplac$Temp==15]))

plac_corrs[33,4]<-fill$estimate
plac_corrs[33,5]<-fill$p.value

# 15 vs. 25
fill<-cor.test(asin(FTP_FLSplac$pupate[FTP_FLSplac$Temp==15]),
               asin(FTP_FLSplac$pupate[FTP_FLSplac$Temp==25]))
# sample size
length(asin(FTP_FLSplac$pupate[FTP_FLSplac$Temp==15]))

plac_corrs[34,4]<-fill$estimate
plac_corrs[34,5]<-fill$p.value

# 20 vs. 25
fill<-cor.test(asin(FTP_FLSplac$pupate[FTP_FLSplac$Temp==25]),
               asin(FTP_FLSplac$pupate[FTP_FLSplac$Temp==20]))
# sample size
length(asin(FTP_FLSplac$pupate[FTP_FLSplac$Temp==25]))

plac_corrs[35,4]<-fill$estimate
plac_corrs[35,5]<-fill$p.value

plac_corrs[33:35,]
# none significantly different from 0


# placida correlations plot (Fig. 7C)----------------------------------------
# to display them properly, we to label them individually with a single label
# (rather than in subset combinations of trait x comparison)
plac_corrs[,7]<-paste(plac_corrs[,2],plac_corrs[,3])
colnames(plac_corrs)[7]<-"Trait_comparison"
plac_corrs$Trait_comparison<-as.factor(plac_corrs$Trait_comparison)
plac_corrs$Trait_comparison<-factor(plac_corrs$Trait_comparison,levels=rev(plac_corrs$Trait_comparison))

p_plac_corrs<-ggplot(data=plac_corrs, aes(x=Trait_comparison, y=r))+
  geom_col(position="dodge") + theme_classic2() + coord_flip() +
  theme(axis.text.y= element_text(size=7), 
        axis.title.y=element_text(angle=0, vjust=0.5), 
        panel.border = element_rect(colour = "black", fill=NA)) +
  geom_hline(yintercept=0) +
  expand_limits(y=c(-1,1)) +
  scale_x_discrete(expand=c(0.02,0.1) 
                   # note: this next line makes the plot really hard to
                   # interpret, but in the final graph, the traits are going to
                   # be be indicated using brackets and sub-brackets. It's
                   # easier to take out the traits with one line of code here
                   # than by changing each by hand in Illustrator. But be sure
                   # to do the labels while looking at a reference plot that has
                   # the full labels!
                   , labels=rev(plac_corrs$Comparison)
  ) +
  scale_y_continuous(breaks=pretty_breaks(5)) + 
  labs(x= "Trait", y= "r", fill="Comparison", title= "placida pairwise correlations", family="sans")

p_plac_corrs

plac_corrs

# Save correlation plots ------------------------------------------------
# ggsave("Fig7_dila.pdf",plot=p_dila_corrs,device="pdf")
# ggsave("Fig7_dors.pdf",plot=p_dors_corrs,device="pdf")
# ggsave("Fig7_plac.pdf",plot=p_plac_corrs,device="pdf")

# Tables S1-3. Broad-sense heritability estimates -----------------------------
# Goal: To calculate the range of possible broad-sense heritability estimates
# for each trait, depending on whether the measured individuals are full or
# half siblings.

# This code is based on Carlos Garcia-Robledo's original code for
# Garcia-Robledo and Horovitz (2012).
# dilaticollis broad-sense heritabilities (Table S1)--------------------------
# Day 1 length
any(FTP_FT1dila$N<5)
# no too-small families

# at 15C
# calculate dilaticollis heritability of day 1 length at 15C
dila_T1_BSH15<-lmer(Tamano_dia_1_mm~1+(1|Pareja),
                    data=dila[dila$Tratamiento_temperatura_C==15 &
                                as.numeric(as.character(dila$Pareja)) %in% 
                                as.numeric(as.character(FTP_FT1dila$Pareja[FTP_FT1dila$Temp==15])),])

summary(dila_T1_BSH15)
# this gives us the variance between families at given temperature (Pareja) &
# the variance between members of the same family (Residual)
# rather than try to reference these (not possible?), I'll hand-enter them
vars_dila_T1_BSH15<-c(0.002933,0.009657)

# calculate the total variance contributed by both levels
VarTot_dila_T1_BSH15<-sum(vars_dila_T1_BSH15)

#Here are the BSH estimates for half sibs and full sibs
# halfsib heritability
4*vars_dila_T1_BSH15[1]/VarTot_dila_T1_BSH15

# fullsib heritability
2*vars_dila_T1_BSH15[1]/VarTot_dila_T1_BSH15

# in the table, we will report the observational variance component 
# (var(Family)/var(Total)) and let readers do the half-sib/full-sib math.
vars_dila_T1_BSH15[1]/VarTot_dila_T1_BSH15

# calculate dilaticollis heritability of day 1 length at 20C
dila_T1_BSH20<-lmer(Tamano_dia_1_mm~1+(1|Pareja),
                    data=dila[dila$Tratamiento_temperatura_C==20 &
                                as.numeric(as.character(dila$Pareja)) %in% 
                                as.numeric(as.character(FTP_FT1dila$Pareja[FTP_FT1dila$Temp==20])),])

summary(dila_T1_BSH20)
# this gives us the variance between families at given temperature (Pareja) &
# the variance between members of the same family (Residual)
# rather than try to reference these (not possible?), I'll hand-enter them
vars_dila_T1_BSH20<-c(0.003604,0.020738)

# calculate the total variance contributed by both levels
VarTot_dila_T1_BSH20<-sum(vars_dila_T1_BSH20)

#Here are the BSH estimates for half sibs and full sibs
# halfsib heritability
4*vars_dila_T1_BSH20[1]/VarTot_dila_T1_BSH20

# fullsib heritability
2*vars_dila_T1_BSH20[1]/VarTot_dila_T1_BSH20

# in the table, we will report the observational variance component 
# (var(Family)/var(Total)) and let readers do the half-sib/full-sib math.
vars_dila_T1_BSH20[1]/VarTot_dila_T1_BSH20


# calculate dilaticollis heritability of day 1 length at 25C
dila_T1_BSH25<-lmer(Tamano_dia_1_mm~1+(1|Pareja),
                    data=dila[dila$Tratamiento_temperatura_C==25 &
                                as.numeric(as.character(dila$Pareja)) %in% 
                                as.numeric(as.character(FTP_FT1dila$Pareja[FTP_FT1dila$Temp==25])),])

summary(dila_T1_BSH25)
# this gives us the variance between families at given temperature (Pareja) &
# the variance between members of the same family (Residual)
# rather than try to reference these (not possible?), I'll hand-enter them
vars_dila_T1_BSH25<-c(0.003115,0.035165)

# calculate the total variance contributed by both levels
VarTot_dila_T1_BSH25<-sum(vars_dila_T1_BSH25)

#Here are the BSH estimates for half sibs and full sibs
# halfsib heritability
4*vars_dila_T1_BSH25[1]/VarTot_dila_T1_BSH25

# fullsib heritability
2*vars_dila_T1_BSH25[1]/VarTot_dila_T1_BSH25

# in the table, we will report the observational variance component 
# (var(Family)/var(Total)) and let readers do the half-sib/full-sib math.
vars_dila_T1_BSH25[1]/VarTot_dila_T1_BSH25

# Day 16 length
any(FTP_FT15dila$N<5)
# there are some too-small families at 15C that will need to be filtered out in
# the model data specification

# calculate dilaticollis heritability of day 16 length at 15C
dila_T15_BSH15<-lmer(Tamano_dia_15_mm~1+(1|Pareja),
                    data=dila[dila$Tratamiento_temperatura_C==15 &
                                as.numeric(as.character(dila$Pareja)) %in% 
                                as.numeric(as.character(
                                  FTP_FT15dila$Pareja[FTP_FT15dila$Temp==15 & FTP_FT15dila$N>=5])),])

summary(dila_T15_BSH15)
# this gives us the variance between families at given temperature (Pareja) &
# the variance between members of the same family (Residual)
# rather than try to reference these (not possible?), I'll hand-enter them
vars_dila_T15_BSH15<-c(0.0383,0.1243)

# calculate the total variance contributed by both levels
VarTot_dila_T15_BSH15<-sum(vars_dila_T15_BSH15)

#Here are the BSH estimates for half sibs and full sibs
# halfsib heritability
4*vars_dila_T15_BSH15[1]/VarTot_dila_T15_BSH15

# fullsib heritability
2*vars_dila_T15_BSH15[1]/VarTot_dila_T15_BSH15

# in the table, we will report the observational variance component 
# (var(Family)/var(Total)) and let readers do the half-sib/full-sib math.
vars_dila_T15_BSH15[1]/VarTot_dila_T15_BSH15


# calculate dilaticollis heritability of day 16 length at 20C
dila_T15_BSH20<-lmer(Tamano_dia_15_mm~1+(1|Pareja),
                     data=dila[dila$Tratamiento_temperatura_C==20 &
                                 as.numeric(as.character(dila$Pareja)) %in% 
                                 as.numeric(as.character(
                                   FTP_FT15dila$Pareja[FTP_FT15dila$Temp==20])),])

summary(dila_T15_BSH20)
# this gives us the variance between families at given temperature (Pareja) &
# the variance between members of the same family (Residual)
# rather than try to reference these (not possible?), I'll hand-enter them
vars_dila_T15_BSH20<-c(0.002891,0.419551)

# calculate the total variance contributed by both levels
VarTot_dila_T15_BSH20<-sum(vars_dila_T15_BSH20)

#Here are the BSH estimates for half sibs and full sibs
# halfsib heritability
4*vars_dila_T15_BSH20[1]/VarTot_dila_T15_BSH20

# fullsib heritability
2*vars_dila_T15_BSH20[1]/VarTot_dila_T15_BSH20

# in the table, we will report the observational variance component 
# (var(Family)/var(Total)) and let readers do the half-sib/full-sib math.
vars_dila_T15_BSH20[1]/VarTot_dila_T15_BSH20


# calculate dilaticollis heritability of day 16 length at 25C
dila_T15_BSH25<-lmer(Tamano_dia_15_mm~1+(1|Pareja),
                     data=dila[dila$Tratamiento_temperatura_C==25 &
                                 as.numeric(as.character(dila$Pareja)) %in% 
                                 as.numeric(as.character(
                                   FTP_FT15dila$Pareja[FTP_FT15dila$Temp==25])),])

summary(dila_T15_BSH25)
# this gives us the variance between families at given temperature (Pareja) &
# the variance between members of the same family (Residual)
# rather than try to reference these (not possible?), I'll hand-enter them
vars_dila_T15_BSH25<-c(0.0000,0.8521)

# calculate the total variance contributed by both levels
VarTot_dila_T15_BSH25<-sum(vars_dila_T15_BSH25)

#Here are the BSH estimates for half sibs and full sibs
# halfsib heritability
4*vars_dila_T15_BSH25[1]/VarTot_dila_T15_BSH25

# fullsib heritability
2*vars_dila_T15_BSH25[1]/VarTot_dila_T15_BSH25

# in the table, we will report the observational variance component 
# (var(Family)/var(Total)) and let readers do the half-sib/full-sib math.
vars_dila_T15_BSH25[1]/VarTot_dila_T15_BSH25


# Pupal length
# no one made it to pupation at 15C. This means that none of the remaining traits,
# (for pupae, adults, and development/survival times) need to be calculated.

# calculate dilaticollis heritability of pupal length at 20C
dila_TP_BSH20<-lmer(Tamano_pupa_mm~1+(1|Pareja),
                     data=dila[dila$Tratamiento_temperatura_C==20 &
                                 as.numeric(as.character(dila$Pareja)) %in% 
                                 as.numeric(as.character(
                                   FTP_FTPdila$Pareja[FTP_FTPdila$Temp==20])),])

summary(dila_TP_BSH20)
# this gives us the variance between families at given temperature (Pareja) &
# the variance between members of the same family (Residual)
# rather than try to reference these (not possible?), I'll hand-enter them
vars_dila_TP_BSH20<-c(0.0561,0.12020)

# calculate the total variance contributed by both levels
VarTot_dila_TP_BSH20<-sum(vars_dila_TP_BSH20)

#Here are the BSH estimates for half sibs and full sibs
# halfsib heritability
4*vars_dila_TP_BSH20[1]/VarTot_dila_TP_BSH20

# fullsib heritability
2*vars_dila_TP_BSH20[1]/VarTot_dila_TP_BSH20

# in the table, we will report the observational variance component 
# (var(Family)/var(Total)) and let readers do the half-sib/full-sib math.
vars_dila_TP_BSH20[1]/VarTot_dila_TP_BSH20


# calculate dilaticollis heritability of pupal length at 25C
dila_TP_BSH25<-lmer(Tamano_pupa_mm~1+(1|Pareja),
                     data=dila[dila$Tratamiento_temperatura_C==25 &
                                 as.numeric(as.character(dila$Pareja)) %in% 
                                 as.numeric(as.character(
                                   FTP_FTPdila$Pareja[FTP_FTPdila$Temp==25])),])

summary(dila_TP_BSH25)
# this gives us the variance between families at given temperature (Pareja) &
# the variance between members of the same family (Residual)
# rather than try to reference these (not possible?), I'll hand-enter them
vars_dila_TP_BSH25<-c(0.0008936,0.108870)

# calculate the total variance contributed by both levels
VarTot_dila_TP_BSH25<-sum(vars_dila_TP_BSH25)

#Here are the BSH estimates for half sibs and full sibs
# halfsib heritability
4*vars_dila_TP_BSH25[1]/VarTot_dila_TP_BSH25

# fullsib heritability
2*vars_dila_TP_BSH25[1]/VarTot_dila_TP_BSH25

# in the table, we will report the observational variance component 
# (var(Family)/var(Total)) and let readers do the half-sib/full-sib math.
vars_dila_TP_BSH25[1]/VarTot_dila_TP_BSH25



# Adult length
# no one made it to adulthood at 15C. 

# calculate dilaticollis heritability of adult length at 20C
dila_TA_BSH20<-lmer(Tamano_adulto_mm~1+(1|Pareja),
                    data=dila[dila$Tratamiento_temperatura_C==20 &
                                as.numeric(as.character(dila$Pareja)) %in% 
                                as.numeric(as.character(
                                  FTP_FTAdila$Pareja[FTP_FTAdila$Temp==20])),])

summary(dila_TA_BSH20)
# this gives us the variance between families at given temperature (Pareja) &
# the variance between members of the same family (Residual)
# rather than try to reference these (not possible?), I'll hand-enter them
vars_dila_TA_BSH20<-c(0.03627,0.04599)

# calculate the total variance contributed by both levels
VarTot_dila_TA_BSH20<-sum(vars_dila_TA_BSH20)

#Here are the BSH estimates for half sibs and full sibs
# halfsib heritability
4*vars_dila_TA_BSH20[1]/VarTot_dila_TA_BSH20

# fullsib heritability
2*vars_dila_TA_BSH20[1]/VarTot_dila_TA_BSH20

# in the table, we will report the observational variance component 
# (var(Family)/var(Total)) and let readers do the half-sib/full-sib math.
vars_dila_TA_BSH20[1]/VarTot_dila_TA_BSH20


# calculate dilaticollis heritability of adult length at 25C
dila_TA_BSH25<-lmer(Tamano_adulto_mm~1+(1|Pareja),
                    data=dila[dila$Tratamiento_temperatura_C==25 &
                                as.numeric(as.character(dila$Pareja)) %in% 
                                as.numeric(as.character(
                                  FTP_FTAdila$Pareja[FTP_FTAdila$Temp==25])),])

summary(dila_TA_BSH25)
# this gives us the variance between families at given temperature (Pareja) &
# the variance between members of the same family (Residual)
# rather than try to reference these (not possible?), I'll hand-enter them
vars_dila_TA_BSH25<-c(5.849*10^-17,0.01040)

# calculate the total variance contributed by both levels
VarTot_dila_TA_BSH25<-sum(vars_dila_TA_BSH25)

#Here are the BSH estimates for half sibs and full sibs
# halfsib heritability
4*vars_dila_TA_BSH25[1]/VarTot_dila_TA_BSH25

# fullsib heritability
2*vars_dila_TA_BSH25[1]/VarTot_dila_TA_BSH25

# in the table, we will report the observational variance component 
# (var(Family)/var(Total)) and let readers do the half-sib/full-sib math.
vars_dila_TA_BSH25[1]/VarTot_dila_TA_BSH25



# Pupal mass
# no one made it to pupation at 15C. 

# calculate dilaticollis heritability of pupal mass at 20C
dila_PP_BSH20<-lmer(Peso_pupa_g~1+(1|Pareja),
                    data=dila[dila$Tratamiento_temperatura_C==20 &
                                as.numeric(as.character(dila$Pareja)) %in% 
                                as.numeric(as.character(
                                  FTP_FPPdila$Pareja[FTP_FPPdila$Temp==20])),])

summary(dila_PP_BSH20)
# this gives us the variance between families at given temperature (Pareja) &
# the variance between members of the same family (Residual)
# rather than try to reference these (not possible?), I'll hand-enter them
vars_dila_PP_BSH20<-c(2.451e-07,1.029e-06)

# calculate the total variance contributed by both levels
VarTot_dila_PP_BSH20<-sum(vars_dila_PP_BSH20)

#Here are the BSH estimates for half sibs and full sibs
# halfsib heritability
4*vars_dila_PP_BSH20[1]/VarTot_dila_PP_BSH20

# fullsib heritability
2*vars_dila_PP_BSH20[1]/VarTot_dila_PP_BSH20

# in the table, we will report the observational variance component 
# (var(Family)/var(Total)) and let readers do the half-sib/full-sib math.
vars_dila_PP_BSH20[1]/VarTot_dila_PP_BSH20


# calculate dilaticollis heritability of pupal mass at 25C
dila_PP_BSH25<-lmer(Peso_pupa_g~1+(1|Pareja),
                    data=dila[dila$Tratamiento_temperatura_C==25 &
                                as.numeric(as.character(dila$Pareja)) %in% 
                                as.numeric(as.character(
                                  FTP_FPPdila$Pareja[FTP_FPPdila$Temp==25])),])

summary(dila_PP_BSH25)
# this gives us the variance between families at given temperature (Pareja) &
# the variance between members of the same family (Residual)
# rather than try to reference these (not possible?), I'll hand-enter them
vars_dila_PP_BSH25<-c(1.617e-07,1.878e-06)

# calculate the total variance contributed by both levels
VarTot_dila_PP_BSH25<-sum(vars_dila_PP_BSH25)

#Here are the BSH estimates for half sibs and full sibs
# halfsib heritability
4*vars_dila_PP_BSH25[1]/VarTot_dila_PP_BSH25

# fullsib heritability
2*vars_dila_PP_BSH25[1]/VarTot_dila_PP_BSH25

# in the table, we will report the observational variance component 
# (var(Family)/var(Total)) and let readers do the half-sib/full-sib math.
vars_dila_PP_BSH25[1]/VarTot_dila_PP_BSH25


# Adult mass
# no one made it to adulthood at 15C. 

# calculate dilaticollis heritability of adult length at 20C
dila_PA_BSH20<-lmer(Peso_adulto_g~1+(1|Pareja),
                    data=dila[dila$Tratamiento_temperatura_C==20 &
                                as.numeric(as.character(dila$Pareja)) %in% 
                                as.numeric(as.character(
                                  FTP_FPAdila$Pareja[FTP_FPAdila$Temp==20])),])

summary(dila_PA_BSH20)
# this gives us the variance between families at given temperature (Pareja) &
# the variance between members of the same family (Residual)
# rather than try to reference these (not possible?), I'll hand-enter them
vars_dila_PA_BSH20<-c(2.433e-07,6.709e-07)

# calculate the total variance contributed by both levels
VarTot_dila_PA_BSH20<-sum(vars_dila_PA_BSH20)

#Here are the BSH estimates for half sibs and full sibs
# halfsib heritability
4*vars_dila_PA_BSH20[1]/VarTot_dila_PA_BSH20

# fullsib heritability
2*vars_dila_PA_BSH20[1]/VarTot_dila_PA_BSH20

# in the table, we will report the observational variance component 
# (var(Family)/var(Total)) and let readers do the half-sib/full-sib math.
vars_dila_PA_BSH20[1]/VarTot_dila_PA_BSH20


# calculate dilaticollis heritability of adult mass at 25C
dila_PA_BSH25<-lmer(Peso_adulto_g~1+(1|Pareja),
                    data=dila[dila$Tratamiento_temperatura_C==25 &
                                as.numeric(as.character(dila$Pareja)) %in% 
                                as.numeric(as.character(
                                  FTP_FPAdila$Pareja[FTP_FPAdila$Temp==25])),])

summary(dila_PA_BSH25)
# this gives us the variance between families at given temperature (Pareja) &
# the variance between members of the same family (Residual)
# rather than try to reference these (not possible?), I'll hand-enter them
vars_dila_PA_BSH25<-c(0.000e+00,7.037e-07)

# calculate the total variance contributed by both levels
VarTot_dila_PA_BSH25<-sum(vars_dila_PA_BSH25)

#Here are the BSH estimates for half sibs and full sibs
# halfsib heritability
4*vars_dila_PA_BSH25[1]/VarTot_dila_PA_BSH25

# fullsib heritability
2*vars_dila_PA_BSH25[1]/VarTot_dila_PA_BSH25

# in the table, we will report the observational variance component 
# (var(Family)/var(Total)) and let readers do the half-sib/full-sib math.
vars_dila_PA_BSH25[1]/VarTot_dila_PA_BSH25


# larval development time
# no one made it to pupation at 15C. 

# calculate dilaticollis heritability of larval development time at 20C
dila_LDT_BSH20<-lmer(ldays~1+(1|Pareja),
                    data=dila[dila$Tratamiento_temperatura_C==20 &
                                as.numeric(as.character(dila$Pareja)) %in% 
                                as.numeric(as.character(
                                  FTP_FLDTdila$Pareja[FTP_FLDTdila$Temp==20])),])

summary(dila_LDT_BSH20)
# this gives us the variance between families at given temperature (Pareja) &
# the variance between members of the same family (Residual)
# rather than try to reference these (not possible?), I'll hand-enter them
vars_dila_LDT_BSH20<-c(1.868e-17,4.410e+02)

# calculate the total variance contributed by both levels
VarTot_dila_LDT_BSH20<-sum(vars_dila_LDT_BSH20)

#Here are the BSH estimates for half sibs and full sibs
# halfsib heritability
4*vars_dila_LDT_BSH20[1]/VarTot_dila_LDT_BSH20

# fullsib heritability
2*vars_dila_LDT_BSH20[1]/VarTot_dila_LDT_BSH20

# in the table, we will report the observational variance component 
# (var(Family)/var(Total)) and let readers do the half-sib/full-sib math.
vars_dila_LDT_BSH20[1]/VarTot_dila_LDT_BSH20


# calculate dilaticollis heritability of larval development time at 25C
dila_LDT_BSH25<-lmer(ldays~1+(1|Pareja),
                    data=dila[dila$Tratamiento_temperatura_C==25 &
                                as.numeric(as.character(dila$Pareja)) %in% 
                                as.numeric(as.character(
                                  FTP_FLDTdila$Pareja[FTP_FLDTdila$Temp==25])),])

summary(dila_LDT_BSH25)
# this gives us the variance between families at given temperature (Pareja) &
# the variance between members of the same family (Residual)
# rather than try to reference these (not possible?), I'll hand-enter them
vars_dila_LDT_BSH25<-c(1.779,209.149)

# calculate the total variance contributed by both levels
VarTot_dila_LDT_BSH25<-sum(vars_dila_LDT_BSH25)

#Here are the BSH estimates for half sibs and full sibs
# halfsib heritability
4*vars_dila_LDT_BSH25[1]/VarTot_dila_LDT_BSH25

# fullsib heritability
2*vars_dila_LDT_BSH25[1]/VarTot_dila_LDT_BSH25

# in the table, we will report the observational variance component 
# (var(Family)/var(Total)) and let readers do the half-sib/full-sib math.
vars_dila_LDT_BSH25[1]/VarTot_dila_LDT_BSH25

# pupal development time
# no one made it to pupation at 15C. 

# calculate dilaticollis heritability of pupal development time at 20C
dila_PDT_BSH20<-lmer(pdays~1+(1|Pareja),
                     data=dila[dila$Tratamiento_temperatura_C==20 &
                                 as.numeric(as.character(dila$Pareja)) %in% 
                                 as.numeric(as.character(
                                   FTP_FPDTdila$Pareja[FTP_FPDTdila$Temp==20])),])

summary(dila_PDT_BSH20)
# this gives us the variance between families at given temperature (Pareja) &
# the variance between members of the same family (Residual)
# rather than try to reference these (not possible?), I'll hand-enter them
vars_dila_PDT_BSH20<-c(1.682e-13,8.121e+01)

# calculate the total variance contributed by both levels
VarTot_dila_PDT_BSH20<-sum(vars_dila_PDT_BSH20)

#Here are the BSH estimates for half sibs and full sibs
# halfsib heritability
4*vars_dila_PDT_BSH20[1]/VarTot_dila_PDT_BSH20

# fullsib heritability
2*vars_dila_PDT_BSH20[1]/VarTot_dila_PDT_BSH20

# in the table, we will report the observational variance component 
# (var(Family)/var(Total)) and let readers do the half-sib/full-sib math.
vars_dila_PDT_BSH20[1]/VarTot_dila_PDT_BSH20


# calculate dilaticollis heritability of pupal development time at 25C
dila_PDT_BSH25<-lmer(pdays~1+(1|Pareja),
                     data=dila[dila$Tratamiento_temperatura_C==25 &
                                 as.numeric(as.character(dila$Pareja)) %in% 
                                 as.numeric(as.character(
                                   FTP_FPDTdila$Pareja[FTP_FPDTdila$Temp==25])),])

summary(dila_PDT_BSH25)
# this gives us the variance between families at given temperature (Pareja) &
# the variance between members of the same family (Residual)
# rather than try to reference these (not possible?), I'll hand-enter them
vars_dila_PDT_BSH25<-c(0.00,59.22)

# calculate the total variance contributed by both levels
VarTot_dila_PDT_BSH25<-sum(vars_dila_PDT_BSH25)

#Here are the BSH estimates for half sibs and full sibs
# halfsib heritability
4*vars_dila_PDT_BSH25[1]/VarTot_dila_PDT_BSH25

# fullsib heritability
2*vars_dila_PDT_BSH25[1]/VarTot_dila_PDT_BSH25

# in the table, we will report the observational variance component 
# (var(Family)/var(Total)) and let readers do the half-sib/full-sib math.
vars_dila_PDT_BSH25[1]/VarTot_dila_PDT_BSH25


# larval survival
# no one made it to pupation at 15C. 

# calculate dilaticollis heritability of larval survival at 20C
dila_LS_BSH20<-lmer(pupate~1+(1|Pareja),
                     data=dila[dila$Tratamiento_temperatura_C==20 &
                                 as.numeric(as.character(dila$Pareja)) %in% 
                                 as.numeric(as.character(
                                   FTP_FLSdila$Pareja[FTP_FLSdila$Temp==20])),])

summary(dila_LS_BSH20)
# this gives us the variance between families at given temperature (Pareja) &
# the variance between members of the same family (Residual)
# rather than try to reference these (not possible?), I'll hand-enter them
vars_dila_LS_BSH20<-c(0.0006114,0.2050348)

# calculate the total variance contributed by both levels
VarTot_dila_LS_BSH20<-sum(vars_dila_LS_BSH20)

#Here are the BSH estimates for half sibs and full sibs
# halfsib heritability
4*vars_dila_LS_BSH20[1]/VarTot_dila_LS_BSH20

# fullsib heritability
2*vars_dila_LS_BSH20[1]/VarTot_dila_LS_BSH20

# in the table, we will report the observational variance component 
# (var(Family)/var(Total)) and let readers do the half-sib/full-sib math.
vars_dila_LS_BSH20[1]/VarTot_dila_LS_BSH20


# calculate dilaticollis heritability of larval survival at 25C
dila_LS_BSH25<-lmer(pupate~1+(1|Pareja),
                     data=dila[dila$Tratamiento_temperatura_C==25 &
                                 as.numeric(as.character(dila$Pareja)) %in% 
                                 as.numeric(as.character(
                                   FTP_FLSdila$Pareja[FTP_FLSdila$Temp==25])),])

summary(dila_LS_BSH25)
# this gives us the variance between families at given temperature (Pareja) &
# the variance between members of the same family (Residual)
# rather than try to reference these (not possible?), I'll hand-enter them
vars_dila_LS_BSH25<-c(5.908e-05,1.751e-01)

# calculate the total variance contributed by both levels
VarTot_dila_LS_BSH25<-sum(vars_dila_LS_BSH25)

#Here are the BSH estimates for half sibs and full sibs
# halfsib heritability
4*vars_dila_LS_BSH25[1]/VarTot_dila_LS_BSH25

# fullsib heritability
2*vars_dila_LS_BSH25[1]/VarTot_dila_LS_BSH25

# in the table, we will report the observational variance component 
# (var(Family)/var(Total)) and let readers do the half-sib/full-sib math.
vars_dila_LS_BSH25[1]/VarTot_dila_LS_BSH25


# dorsalis broad-sense heritabilities (Table S2)------------------------------
# Day 1 length
any(FTP_FT1dors$N<5)
# some too-small families at 15C will need to be filtered out in the data selection

# at 15C
# calculate dorsalis heritability of day 1 length at 15C
dors_T1_BSH15<-lmer(Tamano_dia_1_mm~1+(1|Pareja),
                    data=dors[dors$Tratamiento_temperatura_C==15 &
                                as.numeric(as.character(dors$Pareja)) %in% 
                                as.numeric(as.character(FTP_FT1dors$Pareja[FTP_FT1dors$Temp==15& FTP_FT1dors$N>=5])),])
summary(dors_T1_BSH15)
# this gives us the variance between families at given temperature (Pareja) &
# the variance between members of the same family (Residual)
# rather than try to reference these (not possible?), I'll hand-enter them
vars_dors_T1_BSH15<-c(0.002927,0.008219)

# calculate the total variance contributed by both levels
VarTot_dors_T1_BSH15<-sum(vars_dors_T1_BSH15)

#Here are the BSH estimates for half sibs and full sibs
# halfsib heritability
4*vars_dors_T1_BSH15[1]/VarTot_dors_T1_BSH15

# fullsib heritability
2*vars_dors_T1_BSH15[1]/VarTot_dors_T1_BSH15

# in the table, we will report the observational variance component 
# (var(Family)/var(Total)) and let readers do the half-sib/full-sib math.
vars_dors_T1_BSH15[1]/VarTot_dors_T1_BSH15


# calculate dorsalis heritability of day 1 length at 20C
dors_T1_BSH20<-lmer(Tamano_dia_1_mm~1+(1|Pareja),
                    data=dors[dors$Tratamiento_temperatura_C==20 &
                                as.numeric(as.character(dors$Pareja)) %in% 
                                as.numeric(as.character(FTP_FT1dors$Pareja[FTP_FT1dors$Temp==20])),])
summary(dors_T1_BSH20)
# this gives us the variance between families at given temperature (Pareja) &
# the variance between members of the same family (Residual)
# rather than try to reference these (not possible?), I'll hand-enter them
vars_dors_T1_BSH20<-c(0.001837,0.013659)

# calculate the total variance contributed by both levels
VarTot_dors_T1_BSH20<-sum(vars_dors_T1_BSH20)

#Here are the BSH estimates for half sibs and full sibs
# halfsib heritability
4*vars_dors_T1_BSH20[1]/VarTot_dors_T1_BSH20

# fullsib heritability
2*vars_dors_T1_BSH20[1]/VarTot_dors_T1_BSH20

# in the table, we will report the observational variance component 
# (var(Family)/var(Total)) and let readers do the half-sib/full-sib math.
vars_dors_T1_BSH20[1]/VarTot_dors_T1_BSH20


# calculate dorsalis heritability of day 1 length at 25C
dors_T1_BSH25<-lmer(Tamano_dia_1_mm~1+(1|Pareja),
                    data=dors[dors$Tratamiento_temperatura_C==25 &
                                as.numeric(as.character(dors$Pareja)) %in% 
                                as.numeric(as.character(FTP_FT1dors$Pareja[FTP_FT1dors$Temp==25])),])

summary(dors_T1_BSH25)
# this gives us the variance between families at given temperature (Pareja) &
# the variance between members of the same family (Residual)
# rather than try to reference these (not possible?), I'll hand-enter them
vars_dors_T1_BSH25<-c(0.001067,0.014799)

# calculate the total variance contributed by both levels
VarTot_dors_T1_BSH25<-sum(vars_dors_T1_BSH25)

#Here are the BSH estimates for half sibs and full sibs
# halfsib heritability
4*vars_dors_T1_BSH25[1]/VarTot_dors_T1_BSH25

# fullsib heritability
2*vars_dors_T1_BSH25[1]/VarTot_dors_T1_BSH25

# in the table, we will report the observational variance component 
# (var(Family)/var(Total)) and let readers do the half-sib/full-sib math.
vars_dors_T1_BSH25[1]/VarTot_dors_T1_BSH25


# Day 16 length
FTP_FT15dors
# not enough for 15C analyses

# calculate dorsalis heritability of day 16 length at 20C
dors_T15_BSH20<-lmer(Tamano_dia_15_mm~1+(1|Pareja),
                     data=dors[dors$Tratamiento_temperatura_C==20 &
                                 as.numeric(as.character(dors$Pareja)) %in% 
                                 as.numeric(as.character(
                                   FTP_FT15dors$Pareja[FTP_FT15dors$Temp==20])),])

summary(dors_T15_BSH20)
# this gives us the variance between families at given temperature (Pareja) &
# the variance between members of the same family (Residual)
# rather than try to reference these (not possible?), I'll hand-enter them
vars_dors_T15_BSH20<-c(0.009542,0.165717)

# calculate the total variance contributed by both levels
VarTot_dors_T15_BSH20<-sum(vars_dors_T15_BSH20)

#Here are the BSH estimates for half sibs and full sibs
# halfsib heritability
4*vars_dors_T15_BSH20[1]/VarTot_dors_T15_BSH20

# fullsib heritability
2*vars_dors_T15_BSH20[1]/VarTot_dors_T15_BSH20

# in the table, we will report the observational variance component 
# (var(Family)/var(Total)) and let readers do the half-sib/full-sib math.
vars_dors_T15_BSH20[1]/VarTot_dors_T15_BSH20


# calculate dorsalis heritability of day 16 length at 25C
dors_T15_BSH25<-lmer(Tamano_dia_15_mm~1+(1|Pareja),
                     data=dors[dors$Tratamiento_temperatura_C==25 &
                                 as.numeric(as.character(dors$Pareja)) %in% 
                                 as.numeric(as.character(
                                   FTP_FT15dors$Pareja[FTP_FT15dors$Temp==25])),])
summary(dors_T15_BSH25)
# this gives us the variance between families at given temperature (Pareja) &
# the variance between members of the same family (Residual)
# rather than try to reference these (not possible?), I'll hand-enter them
vars_dors_T15_BSH25<-c(0.0571,0.7634)

# calculate the total variance contributed by both levels
VarTot_dors_T15_BSH25<-sum(vars_dors_T15_BSH25)

#Here are the BSH estimates for half sibs and full sibs
# halfsib heritability
4*vars_dors_T15_BSH25[1]/VarTot_dors_T15_BSH25

# fullsib heritability
2*vars_dors_T15_BSH25[1]/VarTot_dors_T15_BSH25

# in the table, we will report the observational variance component 
# (var(Family)/var(Total)) and let readers do the half-sib/full-sib math.
vars_dors_T15_BSH25[1]/VarTot_dors_T15_BSH25


# Pupal length
# no one made it to pupation at 15C. This means that none of the remaining traits,
# (for pupae, adults, and development/survival times) need to be calculated.

# calculate dorsalis heritability of pupal length at 20C
dors_TP_BSH20<-lmer(Tamano_pupa_mm~1+(1|Pareja),
                    data=dors[dors$Tratamiento_temperatura_C==20 &
                                as.numeric(as.character(dors$Pareja)) %in% 
                                as.numeric(as.character(
                                  FTP_FTPdors$Pareja[FTP_FTPdors$Temp==20])),])

summary(dors_TP_BSH20)
# this gives us the variance between families at given temperature (Pareja) &
# the variance between members of the same family (Residual)
# rather than try to reference these (not possible?), I'll hand-enter them
vars_dors_TP_BSH20<-c(0.0000,0.1051)

# calculate the total variance contributed by both levels
VarTot_dors_TP_BSH20<-sum(vars_dors_TP_BSH20)

#Here are the BSH estimates for half sibs and full sibs
# halfsib heritability
4*vars_dors_TP_BSH20[1]/VarTot_dors_TP_BSH20

# fullsib heritability
2*vars_dors_TP_BSH20[1]/VarTot_dors_TP_BSH20

# in the table, we will report the observational variance component 
# (var(Family)/var(Total)) and let readers do the half-sib/full-sib math.
vars_dors_TP_BSH20[1]/VarTot_dors_TP_BSH20


# calculate dorsalis heritability of pupal length at 25C
dors_TP_BSH25<-lmer(Tamano_pupa_mm~1+(1|Pareja),
                    data=dors[dors$Tratamiento_temperatura_C==25 &
                                as.numeric(as.character(dors$Pareja)) %in% 
                                as.numeric(as.character(
                                  FTP_FTPdors$Pareja[FTP_FTPdors$Temp==25])),])

summary(dors_TP_BSH25)
# this gives us the variance between families at given temperature (Pareja) &
# the variance between members of the same family (Residual)
# rather than try to reference these (not possible?), I'll hand-enter them
vars_dors_TP_BSH25<-c(0.003808,0.106425)

# calculate the total variance contributed by both levels
VarTot_dors_TP_BSH25<-sum(vars_dors_TP_BSH25)

#Here are the BSH estimates for half sibs and full sibs
# halfsib heritability
4*vars_dors_TP_BSH25[1]/VarTot_dors_TP_BSH25

# fullsib heritability
2*vars_dors_TP_BSH25[1]/VarTot_dors_TP_BSH25

# in the table, we will report the observational variance component 
# (var(Family)/var(Total)) and let readers do the half-sib/full-sib math.
vars_dors_TP_BSH25[1]/VarTot_dors_TP_BSH25



# Adult length
# no one made it to adulthood at 15C. 

# calculate dorsalis heritability of adult length at 20C
dors_TA_BSH20<-lmer(Tamano_adulto_mm~1+(1|Pareja),
                    data=dors[dors$Tratamiento_temperatura_C==20 &
                                as.numeric(as.character(dors$Pareja)) %in% 
                                as.numeric(as.character(
                                  FTP_FTAdors$Pareja[FTP_FTAdors$Temp==20])),])

summary(dors_TA_BSH20)
# this gives us the variance between families at given temperature (Pareja) &
# the variance between members of the same family (Residual)
# rather than try to reference these (not possible?), I'll hand-enter them
vars_dors_TA_BSH20<-c(0.00000,0.06609)

# calculate the total variance contributed by both levels
VarTot_dors_TA_BSH20<-sum(vars_dors_TA_BSH20)

#Here are the BSH estimates for half sibs and full sibs
# halfsib heritability
4*vars_dors_TA_BSH20[1]/VarTot_dors_TA_BSH20

# fullsib heritability
2*vars_dors_TA_BSH20[1]/VarTot_dors_TA_BSH20

# in the table, we will report the observational variance component 
# (var(Family)/var(Total)) and let readers do the half-sib/full-sib math.
vars_dors_TA_BSH20[1]/VarTot_dors_TA_BSH20


# calculate dorsalis heritability of adult length at 25C
dors_TA_BSH25<-lmer(Tamano_adulto_mm~1+(1|Pareja),
                    data=dors[dors$Tratamiento_temperatura_C==25 &
                                as.numeric(as.character(dors$Pareja)) %in% 
                                as.numeric(as.character(
                                  FTP_FTAdors$Pareja[FTP_FTAdors$Temp==25])),])

summary(dors_TA_BSH25)
# this gives us the variance between families at given temperature (Pareja) &
# the variance between members of the same family (Residual)
# rather than try to reference these (not possible?), I'll hand-enter them
vars_dors_TA_BSH25<-c(0.00392,0.07339)

# calculate the total variance contributed by both levels
VarTot_dors_TA_BSH25<-sum(vars_dors_TA_BSH25)

#Here are the BSH estimates for half sibs and full sibs
# halfsib heritability
4*vars_dors_TA_BSH25[1]/VarTot_dors_TA_BSH25

# fullsib heritability
2*vars_dors_TA_BSH25[1]/VarTot_dors_TA_BSH25

# in the table, we will report the observational variance component 
# (var(Family)/var(Total)) and let readers do the half-sib/full-sib math.
vars_dors_TA_BSH25[1]/VarTot_dors_TA_BSH25



# Pupal mass
# no one made it to pupation at 15C. 

# calculate dorsalis heritability of pupal mass at 20C
dors_PP_BSH20<-lmer(Peso_pupa_g~1+(1|Pareja),
                    data=dors[dors$Tratamiento_temperatura_C==20 &
                                as.numeric(as.character(dors$Pareja)) %in% 
                                as.numeric(as.character(
                                  FTP_FPPdors$Pareja[FTP_FPPdors$Temp==20])),])

summary(dors_PP_BSH20)
# this gives us the variance between families at given temperature (Pareja) &
# the variance between members of the same family (Residual)
# rather than try to reference these (not possible?), I'll hand-enter them
vars_dors_PP_BSH20<-c(1.594e-07,2.513e-06)

# calculate the total variance contributed by both levels
VarTot_dors_PP_BSH20<-sum(vars_dors_PP_BSH20)

#Here are the BSH estimates for half sibs and full sibs
# halfsib heritability
4*vars_dors_PP_BSH20[1]/VarTot_dors_PP_BSH20

# fullsib heritability
2*vars_dors_PP_BSH20[1]/VarTot_dors_PP_BSH20

# in the table, we will report the observational variance component 
# (var(Family)/var(Total)) and let readers do the half-sib/full-sib math.
vars_dors_PP_BSH20[1]/VarTot_dors_PP_BSH20


# calculate dorsalis heritability of pupal mass at 25C
dors_PP_BSH25<-lmer(Peso_pupa_g~1+(1|Pareja),
                    data=dors[dors$Tratamiento_temperatura_C==25 &
                                as.numeric(as.character(dors$Pareja)) %in% 
                                as.numeric(as.character(
                                  FTP_FPPdors$Pareja[FTP_FPPdors$Temp==25])),])

summary(dors_PP_BSH25)
# this gives us the variance between families at given temperature (Pareja) &
# the variance between members of the same family (Residual)
# rather than try to reference these (not possible?), I'll hand-enter them
vars_dors_PP_BSH25<-c(3.647e-07,3.355e-06)

# calculate the total variance contributed by both levels
VarTot_dors_PP_BSH25<-sum(vars_dors_PP_BSH25)

#Here are the BSH estimates for half sibs and full sibs
# halfsib heritability
4*vars_dors_PP_BSH25[1]/VarTot_dors_PP_BSH25

# fullsib heritability
2*vars_dors_PP_BSH25[1]/VarTot_dors_PP_BSH25

# in the table, we will report the observational variance component 
# (var(Family)/var(Total)) and let readers do the half-sib/full-sib math.
vars_dors_PP_BSH25[1]/VarTot_dors_PP_BSH25


# Adult mass
# no one made it to adulthood at 15C. 
# calculate dorsalis heritability of adult mass at 20C
dors_PA_BSH20<-lmer(Peso_adulto_g~1+(1|Pareja),
                    data=dors[dors$Tratamiento_temperatura_C==20 &
                                as.numeric(as.character(dors$Pareja)) %in% 
                                as.numeric(as.character(
                                  FTP_FPAdors$Pareja[FTP_FPAdors$Temp==20])),])

summary(dors_PA_BSH20)
# this gives us the variance between families at given temperature (Pareja) &
# the variance between members of the same family (Residual)
# rather than try to reference these (not possible?), I'll hand-enter them
vars_dors_PA_BSH20<-c(0.000e+00,9.157e-07)

# calculate the total variance contributed by both levels
VarTot_dors_PA_BSH20<-sum(vars_dors_PA_BSH20)

#Here are the BSH estimates for half sibs and full sibs
# halfsib heritability
4*vars_dors_PA_BSH20[1]/VarTot_dors_PA_BSH20

# fullsib heritability
2*vars_dors_PA_BSH20[1]/VarTot_dors_PA_BSH20

# in the table, we will report the observational variance component 
# (var(Family)/var(Total)) and let readers do the half-sib/full-sib math.
vars_dors_PA_BSH20[1]/VarTot_dors_PA_BSH20


# calculate dorsalis heritability of adult mass at 25C
dors_PA_BSH25<-lmer(Peso_adulto_g~1+(1|Pareja),
                    data=dors[dors$Tratamiento_temperatura_C==25 &
                                as.numeric(as.character(dors$Pareja)) %in% 
                                as.numeric(as.character(
                                  FTP_FPAdors$Pareja[FTP_FPAdors$Temp==25])),])

summary(dors_PA_BSH25)
# this gives us the variance between families at given temperature (Pareja) &
# the variance between members of the same family (Residual)
# rather than try to reference these (not possible?), I'll hand-enter them
vars_dors_PA_BSH25<-c(9.361e-08,1.894e-06)

# calculate the total variance contributed by both levels
VarTot_dors_PA_BSH25<-sum(vars_dors_PA_BSH25)

#Here are the BSH estimates for half sibs and full sibs
# halfsib heritability
4*vars_dors_PA_BSH25[1]/VarTot_dors_PA_BSH25

# fullsib heritability
2*vars_dors_PA_BSH25[1]/VarTot_dors_PA_BSH25

# in the table, we will report the observational variance component 
# (var(Family)/var(Total)) and let readers do the half-sib/full-sib math.
vars_dors_PA_BSH25[1]/VarTot_dors_PA_BSH25


# larval development time
# no one made it to pupation at 15C. 
# calculate dorsalis heritability of larval development time at 20C
dors_LDT_BSH20<-lmer(ldays~1+(1|Pareja),
                     data=dors[dors$Tratamiento_temperatura_C==20 &
                                 as.numeric(as.character(dors$Pareja)) %in% 
                                 as.numeric(as.character(
                                   FTP_FLDTdors$Pareja[FTP_FLDTdors$Temp==20])),])

summary(dors_LDT_BSH20)
# this gives us the variance between families at given temperature (Pareja) &
# the variance between members of the same family (Residual)
# rather than try to reference these (not possible?), I'll hand-enter them
vars_dors_LDT_BSH20<-c(0.0,340.9)

# calculate the total variance contributed by both levels
VarTot_dors_LDT_BSH20<-sum(vars_dors_LDT_BSH20)

#Here are the BSH estimates for half sibs and full sibs
# halfsib heritability
4*vars_dors_LDT_BSH20[1]/VarTot_dors_LDT_BSH20

# fullsib heritability
2*vars_dors_LDT_BSH20[1]/VarTot_dors_LDT_BSH20

# in the table, we will report the observational variance component 
# (var(Family)/var(Total)) and let readers do the half-sib/full-sib math.
vars_dors_LDT_BSH20[1]/VarTot_dors_LDT_BSH20


# calculate dorsalis heritability of larval development time at 25C
dors_LDT_BSH25<-lmer(ldays~1+(1|Pareja),
                     data=dors[dors$Tratamiento_temperatura_C==25 &
                                 as.numeric(as.character(dors$Pareja)) %in% 
                                 as.numeric(as.character(
                                   FTP_FLDTdors$Pareja[FTP_FLDTdors$Temp==25])),])

summary(dors_LDT_BSH25)
# this gives us the variance between families at given temperature (Pareja) &
# the variance between members of the same family (Residual)
# rather than try to reference these (not possible?), I'll hand-enter them
vars_dors_LDT_BSH25<-c(0.6256,210.7220)

# calculate the total variance contributed by both levels
VarTot_dors_LDT_BSH25<-sum(vars_dors_LDT_BSH25)

#Here are the BSH estimates for half sibs and full sibs
# halfsib heritability
4*vars_dors_LDT_BSH25[1]/VarTot_dors_LDT_BSH25

# fullsib heritability
2*vars_dors_LDT_BSH25[1]/VarTot_dors_LDT_BSH25

# in the table, we will report the observational variance component 
# (var(Family)/var(Total)) and let readers do the half-sib/full-sib math.
vars_dors_LDT_BSH25[1]/VarTot_dors_LDT_BSH25

# pupal development time
# no one made it to pupation at 15C. 
# calculate dorsalis heritability of pupal development time at 20C
dors_PDT_BSH20<-lmer(pdays~1+(1|Pareja),
                     data=dors[dors$Tratamiento_temperatura_C==20 &
                                 as.numeric(as.character(dors$Pareja)) %in% 
                                 as.numeric(as.character(
                                   FTP_FPDTdors$Pareja[FTP_FPDTdors$Temp==20])),])

summary(dors_PDT_BSH20)
# this gives us the variance between families at given temperature (Pareja) &
# the variance between members of the same family (Residual)
# rather than try to reference these (not possible?), I'll hand-enter them
vars_dors_PDT_BSH20<-c(17.28,250.58)

# calculate the total variance contributed by both levels
VarTot_dors_PDT_BSH20<-sum(vars_dors_PDT_BSH20)

#Here are the BSH estimates for half sibs and full sibs
# halfsib heritability
4*vars_dors_PDT_BSH20[1]/VarTot_dors_PDT_BSH20

# fullsib heritability
2*vars_dors_PDT_BSH20[1]/VarTot_dors_PDT_BSH20

# in the table, we will report the observational variance component 
# (var(Family)/var(Total)) and let readers do the half-sib/full-sib math.
vars_dors_PDT_BSH20[1]/VarTot_dors_PDT_BSH20


# calculate dorsalis heritability of pupal development time at 25C
dors_PDT_BSH25<-lmer(pdays~1+(1|Pareja),
                     data=dors[dors$Tratamiento_temperatura_C==25 &
                                 as.numeric(as.character(dors$Pareja)) %in% 
                                 as.numeric(as.character(
                                   FTP_FPDTdors$Pareja[FTP_FPDTdors$Temp==25])),])

summary(dors_PDT_BSH25)
# this gives us the variance between families at given temperature (Pareja) &
# the variance between members of the same family (Residual)
# rather than try to reference these (not possible?), I'll hand-enter them
vars_dors_PDT_BSH25<-c(9.269,143.136)

# calculate the total variance contributed by both levels
VarTot_dors_PDT_BSH25<-sum(vars_dors_PDT_BSH25)

#Here are the BSH estimates for half sibs and full sibs
# halfsib heritability
4*vars_dors_PDT_BSH25[1]/VarTot_dors_PDT_BSH25

# fullsib heritability
2*vars_dors_PDT_BSH25[1]/VarTot_dors_PDT_BSH25

# in the table, we will report the observational variance component 
# (var(Family)/var(Total)) and let readers do the half-sib/full-sib math.
vars_dors_PDT_BSH25[1]/VarTot_dors_PDT_BSH25


# larval survival
# no one made it to pupation at 15C. 
# calculate dorsalis heritability of larval survival at 20C
dors_LS_BSH20<-lmer(pupate~1+(1|Pareja),
                    data=dors[dors$Tratamiento_temperatura_C==20 &
                                as.numeric(as.character(dors$Pareja)) %in% 
                                as.numeric(as.character(
                                  FTP_FLSdors$Pareja[FTP_FLSdors$Temp==20])),])

summary(dors_LS_BSH20)
# this gives us the variance between families at given temperature (Pareja) &
# the variance between members of the same family (Residual)
# rather than try to reference these (not possible?), I'll hand-enter them
vars_dors_LS_BSH20<-c(5.555e-17,1.751e-01)

# calculate the total variance contributed by both levels
VarTot_dors_LS_BSH20<-sum(vars_dors_LS_BSH20)

#Here are the BSH estimates for half sibs and full sibs
# halfsib heritability
4*vars_dors_LS_BSH20[1]/VarTot_dors_LS_BSH20

# fullsib heritability
2*vars_dors_LS_BSH20[1]/VarTot_dors_LS_BSH20

# in the table, we will report the observational variance component 
# (var(Family)/var(Total)) and let readers do the half-sib/full-sib math.
vars_dors_LS_BSH20[1]/VarTot_dors_LS_BSH20


# calculate dorsalis heritability of larval survival at 25C
dors_LS_BSH25<-lmer(pupate~1+(1|Pareja),
                    data=dors[dors$Tratamiento_temperatura_C==25 &
                                as.numeric(as.character(dors$Pareja)) %in% 
                                as.numeric(as.character(
                                  FTP_FLSdors$Pareja[FTP_FLSdors$Temp==25])),])

summary(dors_LS_BSH25)
# this gives us the variance between families at given temperature (Pareja) &
# the variance between members of the same family (Residual)
# rather than try to reference these (not possible?), I'll hand-enter them
vars_dors_LS_BSH25<-c(0.002589,0.193668)

# calculate the total variance contributed by both levels
VarTot_dors_LS_BSH25<-sum(vars_dors_LS_BSH25)

#Here are the BSH estimates for half sibs and full sibs
# halfsib heritability
4*vars_dors_LS_BSH25[1]/VarTot_dors_LS_BSH25

# fullsib heritability
2*vars_dors_LS_BSH25[1]/VarTot_dors_LS_BSH25

# in the table, we will report the observational variance component 
# (var(Family)/var(Total)) and let readers do the half-sib/full-sib math.
vars_dors_LS_BSH25[1]/VarTot_dors_LS_BSH25

# placida broad-sense heritabilities (Table S3)------------------------------
# Day 1 length
any(FTP_FT1plac$N<5)
# some too-small families at 15C will need to be filtered out in the data
# selection

# at 15C
# calculate placida heritability of day 1 length at 15C
plac_T1_BSH15<-lmer(Tamano_dia_1_mm~1+(1|Pareja),
                    data=plac[plac$Tratamiento_temperatura_C==15 &
                                as.numeric(as.character(plac$Pareja)) %in% 
                                as.numeric(as.character(FTP_FT1plac$Pareja[FTP_FT1plac$Temp==15& FTP_FT1plac$N>=5])),])



summary(plac_T1_BSH15)
# this gives us the variance between families at given temperature (Pareja) &
# the variance between members of the same family (Residual)
# rather than try to reference these (not possible?), I'll hand-enter them
vars_plac_T1_BSH15<-c(0.004303,0.015025)

# calculate the total variance contributed by both levels
VarTot_plac_T1_BSH15<-sum(vars_plac_T1_BSH15)

#Here are the BSH estimates for half sibs and full sibs
# halfsib heritability
4*vars_plac_T1_BSH15[1]/VarTot_plac_T1_BSH15

# fullsib heritability
2*vars_plac_T1_BSH15[1]/VarTot_plac_T1_BSH15

# in the table, we will report the observational variance component 
# (var(Family)/var(Total)) and let readers do the half-sib/full-sib math.
vars_plac_T1_BSH15[1]/VarTot_plac_T1_BSH15


# calculate placida heritability of day 1 length at 20C
plac_T1_BSH20<-lmer(Tamano_dia_1_mm~1+(1|Pareja),
                    data=plac[plac$Tratamiento_temperatura_C==20 &
                                as.numeric(as.character(plac$Pareja)) %in% 
                                as.numeric(as.character(FTP_FT1plac$Pareja[FTP_FT1plac$Temp==20])),])

summary(plac_T1_BSH20)
# this gives us the variance between families at given temperature (Pareja) &
# the variance between members of the same family (Residual)
# rather than try to reference these (not possible?), I'll hand-enter them
vars_plac_T1_BSH20<-c(1.117e-22,3.183e-02)

# calculate the total variance contributed by both levels
VarTot_plac_T1_BSH20<-sum(vars_plac_T1_BSH20)

#Here are the BSH estimates for half sibs and full sibs
# halfsib heritability
4*vars_plac_T1_BSH20[1]/VarTot_plac_T1_BSH20

# fullsib heritability
2*vars_plac_T1_BSH20[1]/VarTot_plac_T1_BSH20

# in the table, we will report the observational variance component 
# (var(Family)/var(Total)) and let readers do the half-sib/full-sib math.
vars_plac_T1_BSH20[1]/VarTot_plac_T1_BSH20


# calculate placida heritability of day 1 length at 25C
plac_T1_BSH25<-lmer(Tamano_dia_1_mm~1+(1|Pareja),
                    data=plac[plac$Tratamiento_temperatura_C==25 &
                                as.numeric(as.character(plac$Pareja)) %in% 
                                as.numeric(as.character(FTP_FT1plac$Pareja[FTP_FT1plac$Temp==25])),])

summary(plac_T1_BSH25)
# this gives us the variance between families at given temperature (Pareja) &
# the variance between members of the same family (Residual)
# rather than try to reference these (not possible?), I'll hand-enter them
vars_plac_T1_BSH25<-c(0.006495,0.011764)

# calculate the total variance contributed by both levels
VarTot_plac_T1_BSH25<-sum(vars_plac_T1_BSH25)

#Here are the BSH estimates for half sibs and full sibs
# halfsib heritability
4*vars_plac_T1_BSH25[1]/VarTot_plac_T1_BSH25

# fullsib heritability
2*vars_plac_T1_BSH25[1]/VarTot_plac_T1_BSH25

# in the table, we will report the observational variance component 
# (var(Family)/var(Total)) and let readers do the half-sib/full-sib math.
vars_plac_T1_BSH25[1]/VarTot_plac_T1_BSH25


# Day 16 length
FTP_FT15plac
# we actually have enough larvae at 15C to do the analysis!

# calculate placida heritability of day 16 length at 20C
plac_T15_BSH15<-lmer(Tamano_dia_15_mm~1+(1|Pareja),
                     data=plac[plac$Tratamiento_temperatura_C==15 &
                                 as.numeric(as.character(plac$Pareja)) %in% 
                                 as.numeric(as.character(
                                   FTP_FT15plac$Pareja[FTP_FT15plac$Temp==15])),])

summary(plac_T15_BSH15)
# this gives us the variance between families at given temperature (Pareja) &
# the variance between members of the same family (Residual)
# rather than try to reference these (not possible?), I'll hand-enter them
vars_plac_T15_BSH15<-c(0.006745,0.049186)

# calculate the total variance contributed by both levels
VarTot_plac_T15_BSH15<-sum(vars_plac_T15_BSH15)

#Here are the BSH estimates for half sibs and full sibs
# halfsib heritability
4*vars_plac_T15_BSH15[1]/VarTot_plac_T15_BSH15

# fullsib heritability
2*vars_plac_T15_BSH15[1]/VarTot_plac_T15_BSH15

# in the table, we will report the observational variance component 
# (var(Family)/var(Total)) and let readers do the half-sib/full-sib math.
vars_plac_T15_BSH15[1]/VarTot_plac_T15_BSH15



# calculate placida heritability of day 16 length at 20C
plac_T15_BSH20<-lmer(Tamano_dia_15_mm~1+(1|Pareja),
                     data=plac[plac$Tratamiento_temperatura_C==20 &
                                 as.numeric(as.character(plac$Pareja)) %in% 
                                 as.numeric(as.character(
                                   FTP_FT15plac$Pareja[FTP_FT15plac$Temp==20])),])

summary(plac_T15_BSH20)
# this gives us the variance between families at given temperature (Pareja) &
# the variance between members of the same family (Residual)
# rather than try to reference these (not possible?), I'll hand-enter them
vars_plac_T15_BSH20<-c(0.01862,0.14256)

# calculate the total variance contributed by both levels
VarTot_plac_T15_BSH20<-sum(vars_plac_T15_BSH20)

#Here are the BSH estimates for half sibs and full sibs
# halfsib heritability
4*vars_plac_T15_BSH20[1]/VarTot_plac_T15_BSH20

# fullsib heritability
2*vars_plac_T15_BSH20[1]/VarTot_plac_T15_BSH20

# in the table, we will report the observational variance component 
# (var(Family)/var(Total)) and let readers do the half-sib/full-sib math.
vars_plac_T15_BSH20[1]/VarTot_plac_T15_BSH20


# calculate placida heritability of day 16 length at 25C
plac_T15_BSH25<-lmer(Tamano_dia_15_mm~1+(1|Pareja),
                     data=plac[plac$Tratamiento_temperatura_C==25 &
                                 as.numeric(as.character(plac$Pareja)) %in% 
                                 as.numeric(as.character(
                                   FTP_FT15plac$Pareja[FTP_FT15plac$Temp==25])),])

summary(plac_T15_BSH25)
# this gives us the variance between families at given temperature (Pareja) &
# the variance between members of the same family (Residual)
# rather than try to reference these (not possible?), I'll hand-enter them
vars_plac_T15_BSH25<-c(0.0004829,0.1526041)

# calculate the total variance contributed by both levels
VarTot_plac_T15_BSH25<-sum(vars_plac_T15_BSH25)

#Here are the BSH estimates for half sibs and full sibs
# halfsib heritability
4*vars_plac_T15_BSH25[1]/VarTot_plac_T15_BSH25

# fullsib heritability
2*vars_plac_T15_BSH25[1]/VarTot_plac_T15_BSH25

# in the table, we will report the observational variance component 
# (var(Family)/var(Total)) and let readers do the half-sib/full-sib math.
vars_plac_T15_BSH25[1]/VarTot_plac_T15_BSH25


# Tamano_pupa_mm
FTP_FTPplac[FTP_FTPplac$Temp==15,]
# There are 6 families with >=5 pupae at 15C. We'll do the analysis.

# calculate placida heritability of pupal length at 15C
plac_TP_BSH15<-lmer(Tamano_pupa_mm~1+(1|Pareja),
                    data=plac[plac$Tratamiento_temperatura_C==15 &
                                as.numeric(as.character(plac$Pareja)) %in% 
                                as.numeric(as.character(
                                  FTP_FTPplac$Pareja[FTP_FTPplac$Temp==15 & 
                                                       FTP_FTPplac$N>=5])),])

summary(plac_TP_BSH15)
# this gives us the variance between families at given temperature (Pareja) &
# the variance between members of the same family (Residual)
# rather than try to reference these (not possible?), I'll hand-enter them
vars_plac_TP_BSH15<-c(0.02795,0.17360)

# calculate the total variance contributed by both levels
VarTot_plac_TP_BSH15<-sum(vars_plac_TP_BSH15)

#Here are the BSH estimates for half sibs and full sibs
# halfsib heritability
4*vars_plac_TP_BSH15[1]/VarTot_plac_TP_BSH15

# fullsib heritability
2*vars_plac_TP_BSH15[1]/VarTot_plac_TP_BSH15

# in the table, we will report the observational variance component 
# (var(Family)/var(Total)) and let readers do the half-sib/full-sib math.
vars_plac_TP_BSH15[1]/VarTot_plac_TP_BSH15


# calculate placida heritability of pupal length at 20C
plac_TP_BSH20<-lmer(Tamano_pupa_mm~1+(1|Pareja),
                    data=plac[plac$Tratamiento_temperatura_C==20 &
                                as.numeric(as.character(plac$Pareja)) %in% 
                                as.numeric(as.character(
                                  FTP_FTPplac$Pareja[FTP_FTPplac$Temp==20])),])

summary(plac_TP_BSH20)
# this gives us the variance between families at given temperature (Pareja) &
# the variance between members of the same family (Residual)
# rather than try to reference these (not possible?), I'll hand-enter them
vars_plac_TP_BSH20<-c(0.0000,0.1889)

# calculate the total variance contributed by both levels
VarTot_plac_TP_BSH20<-sum(vars_plac_TP_BSH20)

#Here are the BSH estimates for half sibs and full sibs
# halfsib heritability
4*vars_plac_TP_BSH20[1]/VarTot_plac_TP_BSH20

# fullsib heritability
2*vars_plac_TP_BSH20[1]/VarTot_plac_TP_BSH20

# in the table, we will report the observational variance component 
# (var(Family)/var(Total)) and let readers do the half-sib/full-sib math.
vars_plac_TP_BSH20[1]/VarTot_plac_TP_BSH20


# calculate placida heritability of pupal length at 25C
plac_TP_BSH25<-lmer(Tamano_pupa_mm~1+(1|Pareja),
                    data=plac[plac$Tratamiento_temperatura_C==25 &
                                as.numeric(as.character(plac$Pareja)) %in% 
                                as.numeric(as.character(
                                  FTP_FTPplac$Pareja[FTP_FTPplac$Temp==25])),])

summary(plac_TP_BSH25)
# this gives us the variance between families at given temperature (Pareja) &
# the variance between members of the same family (Residual)
# rather than try to reference these (not possible?), I'll hand-enter them
vars_plac_TP_BSH25<-c(0.01047,0.13317)

# calculate the total variance contributed by both levels
VarTot_plac_TP_BSH25<-sum(vars_plac_TP_BSH25)

#Here are the BSH estimates for half sibs and full sibs
# halfsib heritability
4*vars_plac_TP_BSH25[1]/VarTot_plac_TP_BSH25

# fullsib heritability
2*vars_plac_TP_BSH25[1]/VarTot_plac_TP_BSH25

# in the table, we will report the observational variance component 
# (var(Family)/var(Total)) and let readers do the half-sib/full-sib math.
vars_plac_TP_BSH25[1]/VarTot_plac_TP_BSH25



# Adult length
FTP_FTAplac[FTP_FTAplac$Temp==15,]
# no families have enough adults at 15C to analyze 

FTP_FTAplac[FTP_FTAplac$Temp==20,]

# calculate placida heritability of adult length at 20C
plac_TA_BSH20<-lmer(Tamano_adulto_mm~1+(1|Pareja),
                    data=plac[plac$Tratamiento_temperatura_C==20 &
                                as.numeric(as.character(plac$Pareja)) %in% 
                                as.numeric(as.character(
                                  FTP_FTAplac$Pareja[FTP_FTAplac$Temp==20])),])

summary(plac_TA_BSH20)
# this gives us the variance between families at given temperature (Pareja) &
# the variance between members of the same family (Residual)
# rather than try to reference these (not possible?), I'll hand-enter them
vars_plac_TA_BSH20<-c(0.00000,0.2228)

# calculate the total variance contributed by both levels
VarTot_plac_TA_BSH20<-sum(vars_plac_TA_BSH20)

#Here are the BSH estimates for half sibs and full sibs
# halfsib heritability
4*vars_plac_TA_BSH20[1]/VarTot_plac_TA_BSH20

# fullsib heritability
2*vars_plac_TA_BSH20[1]/VarTot_plac_TA_BSH20

# in the table, we will report the observational variance component 
# (var(Family)/var(Total)) and let readers do the half-sib/full-sib math.
vars_plac_TA_BSH20[1]/VarTot_plac_TA_BSH20


# calculate placida heritability of adult length at 25C
plac_TA_BSH25<-lmer(Tamano_adulto_mm~1+(1|Pareja),
                    data=plac[plac$Tratamiento_temperatura_C==25 &
                                as.numeric(as.character(plac$Pareja)) %in% 
                                as.numeric(as.character(
                                  FTP_FTAplac$Pareja[FTP_FTAplac$Temp==25])),])

summary(plac_TA_BSH25)
# this gives us the variance between families at given temperature (Pareja) &
# the variance between members of the same family (Residual)
# rather than try to reference these (not possible?), I'll hand-enter them
vars_plac_TA_BSH25<-c(1.915e-16,1.565e-01)

# calculate the total variance contributed by both levels
VarTot_plac_TA_BSH25<-sum(vars_plac_TA_BSH25)

#Here are the BSH estimates for half sibs and full sibs
# halfsib heritability
4*vars_plac_TA_BSH25[1]/VarTot_plac_TA_BSH25

# fullsib heritability
2*vars_plac_TA_BSH25[1]/VarTot_plac_TA_BSH25

# in the table, we will report the observational variance component 
# (var(Family)/var(Total)) and let readers do the half-sib/full-sib math.
vars_plac_TA_BSH25[1]/VarTot_plac_TA_BSH25



# Pupal mass
# we have enough 15C pupae for the analysis
# calculate placida heritability of pupal mass at 15C
plac_PP_BSH15<-lmer(Peso_pupa_g~1+(1|Pareja),
                    data=plac[plac$Tratamiento_temperatura_C==15 &
                                as.numeric(as.character(plac$Pareja)) %in% 
                                as.numeric(as.character(
                                  FTP_FPPplac$Pareja[FTP_FPPplac$Temp==15 &
                                                       FTP_FPPplac$N>=5])),])

summary(plac_PP_BSH15)
# this gives us the variance between families at given temperature (Pareja) &
# the variance between members of the same family (Residual)
# rather than try to reference these (not possible?), I'll hand-enter them
vars_plac_PP_BSH15<-c(1.185e-06,8.292e-06)

# calculate the total variance contributed by both levels
VarTot_plac_PP_BSH15<-sum(vars_plac_PP_BSH15)

#Here are the BSH estimates for half sibs and full sibs
# halfsib heritability
4*vars_plac_PP_BSH15[1]/VarTot_plac_PP_BSH15

# fullsib heritability
2*vars_plac_PP_BSH15[1]/VarTot_plac_PP_BSH15

# in the table, we will report the observational variance component 
# (var(Family)/var(Total)) and let readers do the half-sib/full-sib math.
vars_plac_PP_BSH15[1]/VarTot_plac_PP_BSH15


# calculate placida heritability of pupal mass at 20C
plac_PP_BSH20<-lmer(Peso_pupa_g~1+(1|Pareja),
                    data=plac[plac$Tratamiento_temperatura_C==20 &
                                as.numeric(as.character(plac$Pareja)) %in% 
                                as.numeric(as.character(
                                  FTP_FPPplac$Pareja[FTP_FPPplac$Temp==20])),])

summary(plac_PP_BSH20)
# this gives us the variance between families at given temperature (Pareja) &
# the variance between members of the same family (Residual)
# rather than try to reference these (not possible?), I'll hand-enter them
vars_plac_PP_BSH20<-c(0.000e+00,9.605e-06)

# calculate the total variance contributed by both levels
VarTot_plac_PP_BSH20<-sum(vars_plac_PP_BSH20)

#Here are the BSH estimates for half sibs and full sibs
# halfsib heritability
4*vars_plac_PP_BSH20[1]/VarTot_plac_PP_BSH20

# fullsib heritability
2*vars_plac_PP_BSH20[1]/VarTot_plac_PP_BSH20

# in the table, we will report the observational variance component 
# (var(Family)/var(Total)) and let readers do the half-sib/full-sib math.
vars_plac_PP_BSH20[1]/VarTot_plac_PP_BSH20


# calculate placida heritability of pupal mass at 25C
plac_PP_BSH25<-lmer(Peso_pupa_g~1+(1|Pareja),
                    data=plac[plac$Tratamiento_temperatura_C==25 &
                                as.numeric(as.character(plac$Pareja)) %in% 
                                as.numeric(as.character(
                                  FTP_FPPplac$Pareja[FTP_FPPplac$Temp==25])),])

summary(plac_PP_BSH25)
# this gives us the variance between families at given temperature (Pareja) &
# the variance between members of the same family (Residual)
# rather than try to reference these (not possible?), I'll hand-enter them
vars_plac_PP_BSH25<-c(0.000e+00,6.907e-06)

# calculate the total variance contributed by both levels
VarTot_plac_PP_BSH25<-sum(vars_plac_PP_BSH25)

#Here are the BSH estimates for half sibs and full sibs
# halfsib heritability
4*vars_plac_PP_BSH25[1]/VarTot_plac_PP_BSH25

# fullsib heritability
2*vars_plac_PP_BSH25[1]/VarTot_plac_PP_BSH25

# in the table, we will report the observational variance component 
# (var(Family)/var(Total)) and let readers do the half-sib/full-sib math.
vars_plac_PP_BSH25[1]/VarTot_plac_PP_BSH25


# Adult mass
# not enough made it to adulthood at 15C 
# calculate placida heritability of adult mass at 20C
plac_PA_BSH20<-lmer(Peso_adulto_g~1+(1|Pareja),
                    data=plac[plac$Tratamiento_temperatura_C==20 &
                                as.numeric(as.character(plac$Pareja)) %in% 
                                as.numeric(as.character(
                                  FTP_FPAplac$Pareja[FTP_FPAplac$Temp==20])),])

summary(plac_PA_BSH20)
# this gives us the variance between families at given temperature (Pareja) &
# the variance between members of the same family (Residual)
# rather than try to reference these (not possible?), I'll hand-enter them
vars_plac_PA_BSH20<-c(0.000e+00,5.229e-06)

# calculate the total variance contributed by both levels
VarTot_plac_PA_BSH20<-sum(vars_plac_PA_BSH20)

#Here are the BSH estimates for half sibs and full sibs
# halfsib heritability
4*vars_plac_PA_BSH20[1]/VarTot_plac_PA_BSH20

# fullsib heritability
2*vars_plac_PA_BSH20[1]/VarTot_plac_PA_BSH20

# in the table, we will report the observational variance component 
# (var(Family)/var(Total)) and let readers do the half-sib/full-sib math.
vars_plac_PA_BSH20[1]/VarTot_plac_PA_BSH20


# calculate placida heritability of adult mass at 25C
plac_PA_BSH25<-lmer(Peso_adulto_g~1+(1|Pareja),
                    data=plac[plac$Tratamiento_temperatura_C==25 &
                                as.numeric(as.character(plac$Pareja)) %in% 
                                as.numeric(as.character(
                                  FTP_FPAplac$Pareja[FTP_FPAplac$Temp==25])),])

summary(plac_PA_BSH25)
# this gives us the variance between families at given temperature (Pareja) &
# the variance between members of the same family (Residual)
# rather than try to reference these (not possible?), I'll hand-enter them
vars_plac_PA_BSH25<-c(0.000e+00,2.719e-06)

# calculate the total variance contributed by both levels
VarTot_plac_PA_BSH25<-sum(vars_plac_PA_BSH25)

#Here are the BSH estimates for half sibs and full sibs
# halfsib heritability
4*vars_plac_PA_BSH25[1]/VarTot_plac_PA_BSH25

# fullsib heritability
2*vars_plac_PA_BSH25[1]/VarTot_plac_PA_BSH25

# in the table, we will report the observational variance component 
# (var(Family)/var(Total)) and let readers do the half-sib/full-sib math.
vars_plac_PA_BSH25[1]/VarTot_plac_PA_BSH25


# larval development time
FTP_FLDTplac[FTP_FLDTplac$Temp==15,]
# we have enough to do the analysis at 15C

# calculate placida heritability of larval development time at 15C
plac_LDT_BSH15<-lmer(ldays~1+(1|Pareja),
                     data=plac[plac$Tratamiento_temperatura_C==15 &
                                 as.numeric(as.character(plac$Pareja)) %in% 
                                 as.numeric(as.character(
                                   FTP_FLDTplac$Pareja[FTP_FLDTplac$Temp==15 &
                                                         FTP_FLDTplac$N>=5])),])

summary(plac_LDT_BSH15)
# this gives us the variance between families at given temperature (Pareja) &
# the variance between members of the same family (Residual)
# rather than try to reference these (not possible?), I'll hand-enter them
vars_plac_LDT_BSH15<-c(5.393e-15,5.203e+03)

# calculate the total variance contributed by both levels
VarTot_plac_LDT_BSH15<-sum(vars_plac_LDT_BSH15)

#Here are the BSH estimates for half sibs and full sibs
# halfsib heritability
4*vars_plac_LDT_BSH15[1]/VarTot_plac_LDT_BSH15

# fullsib heritability
2*vars_plac_LDT_BSH15[1]/VarTot_plac_LDT_BSH15

# in the table, we will report the observational variance component 
# (var(Family)/var(Total)) and let readers do the half-sib/full-sib math.
vars_plac_LDT_BSH15[1]/VarTot_plac_LDT_BSH15


# calculate placida heritability of larval development time at 20C
plac_LDT_BSH20<-lmer(ldays~1+(1|Pareja),
                     data=plac[plac$Tratamiento_temperatura_C==20 &
                                 as.numeric(as.character(plac$Pareja)) %in% 
                                 as.numeric(as.character(
                                   FTP_FLDTplac$Pareja[FTP_FLDTplac$Temp==20])),])

summary(plac_LDT_BSH20)
# this gives us the variance between families at given temperature (Pareja) &
# the variance between members of the same family (Residual)
# rather than try to reference these (not possible?), I'll hand-enter them
vars_plac_LDT_BSH20<-c(1.128e-12,4.175e+02)

# calculate the total variance contributed by both levels
VarTot_plac_LDT_BSH20<-sum(vars_plac_LDT_BSH20)

#Here are the BSH estimates for half sibs and full sibs
# halfsib heritability
4*vars_plac_LDT_BSH20[1]/VarTot_plac_LDT_BSH20

# fullsib heritability
2*vars_plac_LDT_BSH20[1]/VarTot_plac_LDT_BSH20

# in the table, we will report the observational variance component 
# (var(Family)/var(Total)) and let readers do the half-sib/full-sib math.
vars_plac_LDT_BSH20[1]/VarTot_plac_LDT_BSH20


# calculate placida heritability of larval development time at 25C
plac_LDT_BSH25<-lmer(ldays~1+(1|Pareja),
                     data=plac[plac$Tratamiento_temperatura_C==25 &
                                 as.numeric(as.character(plac$Pareja)) %in% 
                                 as.numeric(as.character(
                                   FTP_FLDTplac$Pareja[FTP_FLDTplac$Temp==25])),])

summary(plac_LDT_BSH25)
# this gives us the variance between families at given temperature (Pareja) &
# the variance between members of the same family (Residual)
# rather than try to reference these (not possible?), I'll hand-enter them
vars_plac_LDT_BSH25<-c(8.936,208.890)

# calculate the total variance contributed by both levels
VarTot_plac_LDT_BSH25<-sum(vars_plac_LDT_BSH25)

#Here are the BSH estimates for half sibs and full sibs
# halfsib heritability
4*vars_plac_LDT_BSH25[1]/VarTot_plac_LDT_BSH25

# fullsib heritability
2*vars_plac_LDT_BSH25[1]/VarTot_plac_LDT_BSH25

# in the table, we will report the observational variance component 
# (var(Family)/var(Total)) and let readers do the half-sib/full-sib math.
vars_plac_LDT_BSH25[1]/VarTot_plac_LDT_BSH25

# pupal development time
# not enough made it to adulthood at 15C. 

# calculate placida heritability of pupal development time at 20C
plac_PDT_BSH20<-lmer(pdays~1+(1|Pareja),
                     data=plac[plac$Tratamiento_temperatura_C==20 &
                                 as.numeric(as.character(plac$Pareja)) %in% 
                                 as.numeric(as.character(
                                   FTP_FPDTplac$Pareja[FTP_FPDTplac$Temp==20])),])

summary(plac_PDT_BSH20)
# this gives us the variance between families at given temperature (Pareja) &
# the variance between members of the same family (Residual)
# rather than try to reference these (not possible?), I'll hand-enter them
vars_plac_PDT_BSH20<-c(0.00,38.28)

# calculate the total variance contributed by both levels
VarTot_plac_PDT_BSH20<-sum(vars_plac_PDT_BSH20)

#Here are the BSH estimates for half sibs and full sibs
# halfsib heritability
4*vars_plac_PDT_BSH20[1]/VarTot_plac_PDT_BSH20

# fullsib heritability
2*vars_plac_PDT_BSH20[1]/VarTot_plac_PDT_BSH20

# in the table, we will report the observational variance component 
# (var(Family)/var(Total)) and let readers do the half-sib/full-sib math.
vars_plac_PDT_BSH20[1]/VarTot_plac_PDT_BSH20


# calculate placida heritability of pupal development time at 25C
plac_PDT_BSH25<-lmer(pdays~1+(1|Pareja),
                     data=plac[plac$Tratamiento_temperatura_C==25 &
                                 as.numeric(as.character(plac$Pareja)) %in% 
                                 as.numeric(as.character(
                                   FTP_FPDTplac$Pareja[FTP_FPDTplac$Temp==25])),])

summary(plac_PDT_BSH25)
# this gives us the variance between families at given temperature (Pareja) &
# the variance between members of the same family (Residual)
# rather than try to reference these (not possible?), I'll hand-enter them
vars_plac_PDT_BSH25<-c(0.8819,16.7424)

# calculate the total variance contributed by both levels
VarTot_plac_PDT_BSH25<-sum(vars_plac_PDT_BSH25)

#Here are the BSH estimates for half sibs and full sibs
# halfsib heritability
4*vars_plac_PDT_BSH25[1]/VarTot_plac_PDT_BSH25

# fullsib heritability
2*vars_plac_PDT_BSH25[1]/VarTot_plac_PDT_BSH25

# in the table, we will report the observational variance component 
# (var(Family)/var(Total)) and let readers do the half-sib/full-sib math.
vars_plac_PDT_BSH25[1]/VarTot_plac_PDT_BSH25


# larval survival
FTP_FLSplac[FTP_FLSplac$Temp==15,] 
# enough families for calculating the heritability at 15C

# calculate placida heritability of larval survival at 15C
plac_LS_BSH15<-lmer(pupate~1+(1|Pareja),
                    data=plac[plac$Tratamiento_temperatura_C==15 &
                                as.numeric(as.character(plac$Pareja)) %in% 
                                as.numeric(as.character(
                                  FTP_FLSplac$Pareja[FTP_FLSplac$Temp==15])),])

summary(plac_LS_BSH15)
# this gives us the variance between families at given temperature (Pareja) &
# the variance between members of the same family (Residual)
# rather than try to reference these (not possible?), I'll hand-enter them
vars_plac_LS_BSH15<-c(0,0.2414)

# calculate the total variance contributed by both levels
VarTot_plac_LS_BSH15<-sum(vars_plac_LS_BSH15)

#Here are the BSH estimates for half sibs and full sibs
# halfsib heritability
4*vars_plac_LS_BSH15[1]/VarTot_plac_LS_BSH15

# fullsib heritability
2*vars_plac_LS_BSH15[1]/VarTot_plac_LS_BSH15

# in the table, we will report the observational variance component 
# (var(Family)/var(Total)) and let readers do the half-sib/full-sib math.
vars_plac_LS_BSH15[1]/VarTot_plac_LS_BSH15



# calculate placida heritability of larval survival at 20C
plac_LS_BSH20<-lmer(pupate~1+(1|Pareja),
                    data=plac[plac$Tratamiento_temperatura_C==20 &
                                as.numeric(as.character(plac$Pareja)) %in% 
                                as.numeric(as.character(
                                  FTP_FLSplac$Pareja[FTP_FLSplac$Temp==20])),])

summary(plac_LS_BSH20)
# this gives us the variance between families at given temperature (Pareja) &
# the variance between members of the same family (Residual)
# rather than try to reference these (not possible?), I'll hand-enter them
vars_plac_LS_BSH20<-c(0.0000,0.1152)

# calculate the total variance contributed by both levels
VarTot_plac_LS_BSH20<-sum(vars_plac_LS_BSH20)

#Here are the BSH estimates for half sibs and full sibs
# halfsib heritability
4*vars_plac_LS_BSH20[1]/VarTot_plac_LS_BSH20

# fullsib heritability
2*vars_plac_LS_BSH20[1]/VarTot_plac_LS_BSH20

# in the table, we will report the observational variance component 
# (var(Family)/var(Total)) and let readers do the half-sib/full-sib math.
vars_plac_LS_BSH20[1]/VarTot_plac_LS_BSH20


# calculate placida heritability of larval survival at 25C
plac_LS_BSH25<-lmer(pupate~1+(1|Pareja),
                    data=plac[plac$Tratamiento_temperatura_C==25 &
                                as.numeric(as.character(plac$Pareja)) %in% 
                                as.numeric(as.character(
                                  FTP_FLSplac$Pareja[FTP_FLSplac$Temp==25])),])

summary(plac_LS_BSH25)
# this gives us the variance between families at given temperature (Pareja) &
# the variance between members of the same family (Residual)
# rather than try to reference these (not possible?), I'll hand-enter them
vars_plac_LS_BSH25<-c(0.005599,0.118839)

# calculate the total variance contributed by both levels
VarTot_plac_LS_BSH25<-sum(vars_plac_LS_BSH25)

#Here are the BSH estimates for half sibs and full sibs
# halfsib heritability
4*vars_plac_LS_BSH25[1]/VarTot_plac_LS_BSH25

# fullsib heritability
2*vars_plac_LS_BSH25[1]/VarTot_plac_LS_BSH25

# in the table, we will report the observational variance component 
# (var(Family)/var(Total)) and let readers do the half-sib/full-sib math.
vars_plac_LS_BSH25[1]/VarTot_plac_LS_BSH25


# Literature Cited --------------------------------------------------------
# Crawley, M. J. 2013. The R book. John Wiley & Sons. 
# García-Robledo, C. & Horvitz, C. C. 2012. Jack of all trades masters novel
# host plants: positive genetic correlations in specialist and generalist insect
# herbivores expanding their diets to novel hosts. Journal of evolutionary
# biology 25: 38-53.
# Therneau, T. M. 2017. Contrasts, populations, and "type III" tests.