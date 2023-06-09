#### Isotopic Sex Differences Meta-Analysis
#### Joshua Bauld

#### Setting up dataframe and calculating mean differences and lnVR for isotopic sex differences.

remove(list = ls())


## clear console and plots

library(metafor)
library(broom)
library(ggtext)
library(ggbeeswarm)
library(ggpubr)
library(gt)
library(gtsummary)
library(flextable)
library(scales)
library(ape)
library(rotl)
library(tidyverse)

## set wd

# read in data 

df_all <- 


## Calculating Effect Sizes

df_meanN <-escalc(measure="MD", 
                 n1i=n_male, n2i=n_female, 
                 m1i=M_M_d15N, m2i=M_F_d15N, 
                 sd1i=SD_M_d15N, sd2i=SD_F_d15N,
                 data=df_all,
                 var.names=c("Mean_DiffN","Mean_DiffN.sv"), add.measure=FALSE,
                 append=TRUE)

str(df_meanN)

df_VRN <-escalc(measure="VR", 
               n1i=n_male, n2i=n_female, 
               m1i=M_M_d15N, m2i=M_F_d15N, 
               sd1i=SD_M_d15N, sd2i=SD_F_d15N,
               data=df_all,
               var.names=c("lnVRN","lnVRN.sv"), add.measure=FALSE,
               append=TRUE)


str(df_VRN)

## NA's created when n=1 for either sex.


df_meanC <-escalc(measure="MD", 
                  n1i=n_male, n2i=n_female, 
                  m1i=M_M_d13C, m2i=M_F_d13C, 
                  sd1i=SD_M_d13C, sd2i=SD_F_d13C,
                  data=df_all,
                  var.names=c("Mean_DiffC","Mean_DiffC.sv"), add.measure=FALSE,
                  append=TRUE)
str(df_meanC)

df_VRC <-escalc(measure="VR", 
                n1i=n_male, n2i=n_female, 
                m1i=M_M_d13C, m2i=M_F_d13C, 
                sd1i=SD_M_d13C, sd2i=SD_F_d13C,
                data=df_all,
                var.names=c("lnVRC","lnVRC.sv"), add.measure=FALSE,
                append=TRUE)


str(df_VRC)

## NA's created when n=1 for either sex.

## combining new effect sizes into one data frame

## first drop all columns, aside from paper number and effect sizes, from df_sd, df_meanC and df_sdC


df_VRN <- df_VRN %>%
  select(Paper_Number, lnVRN, lnVRN.sv)

df_meanC <- df_meanC %>%
  select(Paper_Number, Mean_DiffC, Mean_DiffC.sv)

df_VRC <- df_VRC %>%
  select(Paper_Number, lnVRC, lnVRC.sv)


## then add these to df_effects

df_effects <- bind_cols(df_meanN, df_VRN, df_meanC, df_VRC)

## drop unneeded paper number columns and rename first one

df_effects <- df_effects %>%
  select(-Paper_Number...26, -Paper_Number...29, -Paper_Number...32)

df_effects <- df_effects %>%
  mutate(Paper_Number = Paper_Number...1)

## calculate ssd, size 

df_effects <- df_effects %>% 
  mutate(SSD = if_else(Weight_Female.Kg. < Weight_Male.Kg., ((Weight_Male.Kg.-Weight_Female.Kg.)/Weight_Female.Kg.), -((Weight_Female.Kg.-Weight_Male.Kg.)/Weight_Male.Kg.))) %>%
  mutate(Size = ((Weight_Female.Kg.+Weight_Male.Kg.)/2))

## drop unneeded variables

df_effects <- df_effects %>%
  select(Paper_Number, Year, Subphylum, Genus, Species, Diet, Gape_Lim, SSD, Size, Tissue_type, n_female, n_male, Mean_DiffN, Mean_DiffN.sv, lnVRN, lnVRN.sv, Mean_DiffC, Mean_DiffC.sv, lnVRC, lnVRC.sv)

names(df_effects)


## drop na's - not those resulting from effect size caluclations, but those in which information was not gathered during data collection e.g. no body mass data for a particular species

df_effects <- df_effects %>%
  drop_na(SSD) %>%
  drop_na(Size) %>%
  drop_na(Diet) 

head(df_effects)

## combine genus and species to a single column

df_effects <- df_effects %>% 
  unite(col = Species, c(Genus, Species), sep = " ")



## write csv

write_csv(df_effects, "Isotopic_Sex_Differences_Effect_Sizes.csv")
