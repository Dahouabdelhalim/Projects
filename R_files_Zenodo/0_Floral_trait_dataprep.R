# This R script combines the floral trait data for 17 floral traits 
# across 52 species using data from greenhouse measurements, previously 
# published studies, and photographs. It then imputes 34 missing trait 
# values, and performs factor analysis on the imputed data. It uses the 
# input files "Costus_traits_greenhouse.csv" and "Costus_traits_Maas.csv"

library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(corrplot)
library(Hmisc)
library(PerformanceAnalytics)
library(FactoMineR)
library(factoextra) 
library(missMDA)

setwd("~/Dropbox/Costus grant/Dena and Kathleen/FloralEvolution/Dryad_Deposit")

#We have floral trait data from Maas (and other monographs), the greenhouse, the field, 
#and a few photo sources (e.g., Dave Skinner's website, our own field photos). 
#First we asked, for species that overlap, do the measures correspond? 
#This analysis is at the end of the script. Based on how well the measures agreed
#between monographs and our greenhouse data, we decided we could use both. We gave
#precedence to our greenhouse data because we measured all the traits, whereas 
#Maas was often missing some. The one discrepancy is that Maas measured stamen and 
#labellum length from above the corolla tube, whereas we measured from the ovary.
#To correct this, we added the length of the corolla tube to Maas' stamen and 
#labellum measures.

#first read in data from monographs, primarily authored by P.J.Maas
Maas <- read.csv("Costus_traits_Maas.csv")
#next read in date from traits measured in the greenhouse and field and on photos from field and other collections
Greenhouse <- read.csv("Costus_traits_greenhouse.csv")

#Maas sometimes reported max and min values for each trait. When he did, we took the midpoint. 
#When he only reported one value, we used that.
#subset the data for just the columns we want

Maas <- select(Maas, !contains("Max"))
Maas <- select(Maas, !contains("Min"))

#subset Maas data set getting rid of verbal descriptions of color and variables relating to extrafloral nectaries
Maas <- Maas %>% select(Species_Maas, sp_tip_label, Syndrome, Stamen_Length_mid,	Stamen_Width_mid,
                        Corolla_Length_mid,	Corolla_Tube_Length_mid,	Corolla_Lobe_Length_mid,
                        Labellum_Length_mid,	Labellum_Width_mid,	Anther_Length_mid, bract_color_K,
                        corolla_color_K,	lab_color_K,	yellow_labstripe,	red_labstripe,	
                        stamen_color_K,	stamen_tip_color,	nectar_guides, Style_length_mid)

#Maas did not report stamen exsertion directly, so we calculated it as the
#stamen length minus the labellum length
#add calculated exsertion column to Maas
Maas <- dplyr::mutate(Maas, Stamen_exsertion_mid=Stamen_Length_mid - Labellum_Length_mid)


#adjust his labellum and stamen columns by adding tube length (which he left out)
Maas <- dplyr::mutate(Maas, Labellum_Length_mid = Labellum_Length_mid + Corolla_Tube_Length_mid)
Maas <- dplyr::mutate(Maas, Stamen_Length_mid = Stamen_Length_mid + Corolla_Tube_Length_mid)

#make sure factors are factors
Maas$bract_color_K <- as.factor(Maas$bract_color_K)
Maas$corolla_color_K <- as.factor(Maas$corolla_color_K)
Maas$lab_color_K <- as.factor(Maas$lab_color_K)
Maas$yellow_labstripe <- as.factor(Maas$yellow_labstripe)
Maas$red_labstripe <- as.factor(Maas$red_labstripe)
Maas$stamen_color_K <- as.factor(Maas$stamen_color_K)
Maas$stamen_tip_color <- as.factor(Maas$stamen_tip_color)
Maas$nectar_guides <- as.factor(Maas$nectar_guides)

#remove rows for samples not in phylogeny
Maas$sp_tip_label[Maas$sp_tip_label == ""] <- NA
Maas <- filter(Maas, !is.na(sp_tip_label))


#subset greenhouse data set getting rid of unnecessary columns
Greenhouse <- Greenhouse %>% select(Greenhouse_ID, tip_label,	sp_tip_label,	Syndrome,	Corolla_Length,
                        Corolla_Lobe_Length,	Corolla_Tube_Length,	Stamen_exsertion,	Labellum_Length,
                        Labellum_Width,	Stamen_Length,
                        Stamen_Width,	Anther_Length, Style_length, bract_color_K,
                        corolla_color_K,	lab_color_K,	yellow_labstripe,	red_labstripe,	
                        stamen_color_K,	stamen_tip_color,	nectar_guides)

#make sure factors are factors
Greenhouse$bract_color_K <- as.factor(Greenhouse$bract_color_K)
Greenhouse$corolla_color_K <- as.factor(Greenhouse$corolla_color_K)
Greenhouse$lab_color_K <- as.factor(Greenhouse$lab_color_K)
Greenhouse$yellow_labstripe <- as.factor(Greenhouse$yellow_labstripe)
Greenhouse$red_labstripe <- as.factor(Greenhouse$red_labstripe)
Greenhouse$stamen_color_K <- as.factor(Greenhouse$stamen_color_K)
Greenhouse$stamen_tip_color <- as.factor(Greenhouse$stamen_tip_color)
Greenhouse$nectar_guides <- as.factor(Greenhouse$nectar_guides)

#summarize (average or mode, depending on whether the data is continuous or categorical)
#data within greenhouse and Maas data frames
#first take average/mode of individuals measured multiple times
#then take average/mode of different individuals from the same population
#then take average/mode of different populations of the same species

#make a function to calculate the mode of categorical variables
Mode <- function(x) {
  uniqx <- unique(x)
  uniqx[which.max(tabulate(match(x, uniqx)))]
}

#summarize continuous traits by individual by taking the mean
#also add "_mid" to the trait names to make them comparable to the Maas data
Greenhouse_by_ind <- Greenhouse %>% group_by(Greenhouse_ID) %>% summarise(
  tibble(
    across(where(is.double), mean, na.rm = TRUE, .names = "{.col}_mid"),
    across(where(is.factor), unique)
  )
)
Greenhouse_by_ind <- as.data.frame(Greenhouse_by_ind)

#now take the mode of categorical variables
Greenhouse_by_ind <- Greenhouse_by_ind %>% group_by(Greenhouse_ID) %>% mutate(bract_color_K=Mode(bract_color_K))
Greenhouse_by_ind <- Greenhouse_by_ind %>% group_by(Greenhouse_ID) %>% mutate(corolla_color_K=Mode(corolla_color_K))
Greenhouse_by_ind <- Greenhouse_by_ind %>% group_by(Greenhouse_ID) %>% mutate(lab_color_K=Mode(lab_color_K))
Greenhouse_by_ind <- Greenhouse_by_ind %>% group_by(Greenhouse_ID) %>% mutate(yellow_labstripe=Mode(yellow_labstripe))
Greenhouse_by_ind <- Greenhouse_by_ind %>% group_by(Greenhouse_ID) %>% mutate(red_labstripe=Mode(red_labstripe))
Greenhouse_by_ind <- Greenhouse_by_ind %>% group_by(Greenhouse_ID) %>% mutate(stamen_color_K=Mode(stamen_color_K))
Greenhouse_by_ind <- Greenhouse_by_ind %>% group_by(Greenhouse_ID) %>% mutate(stamen_tip_color=Mode(stamen_tip_color))
Greenhouse_by_ind <- Greenhouse_by_ind %>% group_by(Greenhouse_ID) %>% mutate(nectar_guides=Mode(nectar_guides))

#now remove duplicates by ind
Greenhouse_by_ind <- Greenhouse_by_ind %>% group_by(Greenhouse_ID) %>% distinct()

#add back in columns that were deleted
tip_labels <- dplyr::select(Greenhouse, Greenhouse_ID,tip_label, sp_tip_label, Syndrome)
tip_labels <- distinct(tip_labels)

Greenhouse_by_ind <- full_join(Greenhouse_by_ind, tip_labels, by = "Greenhouse_ID")

#summarize by population
#first take the mean of continuous traits
Greenhouse_by_pop <- Greenhouse_by_ind %>% group_by(tip_label) %>% summarise(
  tibble(
    across(where(is.double), mean, na.rm = TRUE),
    across(where(is.factor), unique)
    )
  )
Greenhouse_by_pop <- as.data.frame(Greenhouse_by_pop)

#now take the mode of categorical variables
Greenhouse_by_pop <- Greenhouse_by_pop %>% group_by(tip_label) %>% mutate(bract_color_K=Mode(bract_color_K))
Greenhouse_by_pop <- Greenhouse_by_pop %>% group_by(tip_label) %>% mutate(corolla_color_K=Mode(corolla_color_K))
Greenhouse_by_pop <- Greenhouse_by_pop %>% group_by(tip_label) %>% mutate(lab_color_K=Mode(lab_color_K))
Greenhouse_by_pop <- Greenhouse_by_pop %>% group_by(tip_label) %>% mutate(yellow_labstripe=Mode(yellow_labstripe))
Greenhouse_by_pop <- Greenhouse_by_pop %>% group_by(tip_label) %>% mutate(red_labstripe=Mode(red_labstripe))
Greenhouse_by_pop <- Greenhouse_by_pop %>% group_by(tip_label) %>% mutate(stamen_color_K=Mode(stamen_color_K))
Greenhouse_by_pop <- Greenhouse_by_pop %>% group_by(tip_label) %>% mutate(stamen_tip_color=Mode(stamen_tip_color))
Greenhouse_by_pop <- Greenhouse_by_pop %>% group_by(tip_label) %>% mutate(nectar_guides=Mode(nectar_guides))

#now remove duplicates by pop
Greenhouse_by_pop <- Greenhouse_by_pop %>% group_by(tip_label) %>% distinct()

#add back in columns that were deleted
tip_labels <- dplyr::select(Greenhouse,tip_label, sp_tip_label, Syndrome)
tip_labels <- distinct(tip_labels)

Greenhouse_by_pop <- full_join(Greenhouse_by_pop, tip_labels, by = "tip_label")

#now summarise by species tip label to get means for each tip in the species tree
Greenhouse_by_species <- Greenhouse_by_pop %>% group_by(sp_tip_label) %>% summarise(
  tibble(
    across(where(is.double), mean, na.rm=TRUE),
    across(where(is.factor), unique)
    )
)

#take mode for polymorphic species
Greenhouse_by_species <- Greenhouse_by_species %>% group_by(sp_tip_label) %>% mutate(bract_color_K=Mode(bract_color_K))
Greenhouse_by_species <- Greenhouse_by_species %>% group_by(sp_tip_label) %>% mutate(corolla_color_K=Mode(corolla_color_K))
Greenhouse_by_species <- Greenhouse_by_species %>% group_by(sp_tip_label) %>% mutate(lab_color_K=Mode(lab_color_K))
Greenhouse_by_species <- Greenhouse_by_species %>% group_by(sp_tip_label) %>% mutate(yellow_labstripe=Mode(yellow_labstripe))
Greenhouse_by_species <- Greenhouse_by_species %>% group_by(sp_tip_label) %>% mutate(red_labstripe=Mode(red_labstripe))
Greenhouse_by_species <- Greenhouse_by_species %>% group_by(sp_tip_label) %>% mutate(stamen_color_K=Mode(stamen_color_K))
Greenhouse_by_species <- Greenhouse_by_species %>% group_by(sp_tip_label) %>% mutate(stamen_tip_color=Mode(stamen_tip_color))
Greenhouse_by_species <- Greenhouse_by_species %>% group_by(sp_tip_label) %>% mutate(nectar_guides=Mode(nectar_guides))

#now remove duplicate rows
Greenhouse_by_species <- Greenhouse_by_species %>% group_by(sp_tip_label) %>% distinct()

#remove rows for samples not in phylogeny
Greenhouse_by_species <- filter(Greenhouse_by_species, !is.na(sp_tip_label))
#remove species with poor data
Greenhouse_by_species <- Greenhouse_by_species %>% filter(sp_tip_label != "Costus_sp_19207", sp_tip_label != "Costus_sp_19113", sp_tip_label != "Costus_lucanusianus_99103" )

#add the Syndrome back in
sp_tip_labels <- dplyr::select(Greenhouse,sp_tip_label, Syndrome)
sp_tip_labels <- distinct(sp_tip_labels)
Greenhouse_by_species <- inner_join(Greenhouse_by_species, sp_tip_labels, by = "sp_tip_label")

#join greenhouse and Maas data together by species tip label
#give precedence to my data
#first remove species for which I or Maas have data with multiple missing traits, then find sp in Maas that I didn't sample
Greenhouse_by_species <- Greenhouse_by_species %>% filter(sp_tip_label!="Costus_vinosus_19273", sp_tip_label!= "Costus_dubius_98074", sp_tip_label!="Costus_vargasii_19210")
Maas <- filter(Maas, sp_tip_label!="Costus_zingiberoides_98162")
Maas_new <- anti_join(Maas, Greenhouse_by_species, by = "sp_tip_label")
#drop rows not in phylogeny
Maas_new <- filter(Maas_new, !is.na(sp_tip_label))

#remove columns not present in both
Maas_new <- Maas_new %>% select(-Species_Maas)
traits <- union(Greenhouse_by_species,Maas_new)
traits$Syndrome <- as.factor(traits$Syndrome)
traits$sp_tip_label <- as.factor(traits$sp_tip_label)

#rename variables so they are easier to interpret
colnames(traits)<-gsub("_mid","",colnames(traits))
colnames(traits)<-gsub("_K","",colnames(traits))

## fix missing values -- NaN to NA
traits$Stamen_exsertion[is.nan(traits$Stamen_exsertion)]<-NA
traits$Labellum_Length[is.nan(traits$Labellum_Length)]<-NA
traits$Labellum_Width[is.nan(traits$Labellum_Width)]<-NA
traits$Stamen_Length[is.nan(traits$Stamen_Length)]<-NA
traits$Stamen_Width[is.nan(traits$Stamen_Width)]<-NA
traits$Style_length[is.nan(traits$Style_length)]<-NA

#remove nectar guide character, since it is already captured by 
#the yellow_labstripe and red_labstripe characters
traits <- traits %>% select(-nectar_guides)

###########################
#checking the dataset variable by variable
#not used in paper, but used to check data
pdf(file="corolla_length.pdf", width = 10, height = 10, pointsize = 6)
cor_len <- ggplot(data = traits, aes(x=Syndrome, y=Corolla_Length)) + 
  geom_dotplot(binaxis='y', stackdir='center') +
  geom_text_repel(aes(label = traits$sp_tip_label), box.padding = unit(0.5,"lines"))
cor_len
dev.off()

pdf(file="corolla_lobe.pdf", width = 10, height = 10, pointsize = 6)
cor_lobe <- ggplot(data = traits, aes(x=Syndrome, y=Corolla_Lobe_Length)) + 
  geom_dotplot(binaxis='y', stackdir='center') +
  geom_text_repel(aes(label = traits$sp_tip_label), box.padding = unit(0.5,"lines"))
cor_lobe
dev.off()

pdf(file="corolla_tube.pdf", width = 10, height = 10, pointsize = 5)
cor_tube <- ggplot(data = traits, aes(x=Syndrome, y=Corolla_Tube_Length)) + 
  geom_dotplot(binaxis='y', stackdir='center') +
  geom_text_repel(aes(label = traits$sp_tip_label), box.padding = unit(0.5,"lines"))
cor_tube
dev.off()

pdf(file="stamen_exsert.pdf", width = 10, height = 10, pointsize = 5)
stam_ex <- ggplot(data = traits, aes(x=Syndrome, y=Stamen_exsertion)) + 
  geom_dotplot(binaxis='y', stackdir='center') +
  geom_text_repel(aes(label = traits$sp_tip_label), box.padding = unit(0.5,"lines"))
stam_ex
dev.off()

pdf(file="labellum_length.pdf", width = 10, height = 10, pointsize = 5)
lab_length <- ggplot(data = traits, aes(x=Syndrome, y=Labellum_Length)) + 
  geom_dotplot(binaxis='y', stackdir='center') +
  geom_text_repel(aes(label = traits$sp_tip_label), box.padding = unit(0.5,"lines"))
lab_length
dev.off()

pdf(file="labellum_width.pdf", width = 10, height = 10, pointsize = 5)
lab_width <- ggplot(data = traits, aes(x=Syndrome, y=Labellum_Width)) + 
  geom_dotplot(binaxis='y', stackdir='center') +
  geom_text_repel(aes(label = traits$sp_tip_label), box.padding = unit(0.5,"lines"))
lab_width
dev.off()

pdf(file="stamen_length.pdf", width = 10, height = 10, pointsize = 5)
stam_len <- ggplot(data = traits, aes(x=Syndrome, y=Stamen_Length)) + 
  geom_dotplot(binaxis='y', stackdir='center') +
  geom_text_repel(aes(label = traits$sp_tip_label), box.padding = unit(0.5,"lines"))
stam_len
dev.off()

pdf(file="stamen_width.pdf", width = 10, height = 10, pointsize = 5)
stam_wid <- ggplot(data = traits, aes(x=Syndrome, y=Stamen_Width)) + 
  geom_dotplot(binaxis='y', stackdir='center') +
  geom_text_repel(aes(label = traits$sp_tip_label), box.padding = unit(0.5,"lines"))
stam_wid
dev.off()

pdf(file="anther_length.pdf", width = 10, height = 10, pointsize = 5)
anth_len <- ggplot(data = traits, aes(x=Syndrome, y=Anther_Length)) + 
  geom_dotplot(binaxis='y', stackdir='center') +
  geom_text_repel(aes(label = traits$sp_tip_label), box.padding = unit(0.5,"lines"))
anth_len
dev.off()

pdf(file="style_length.pdf", width = 10, height = 10, pointsize = 5)
styl_len <- ggplot(data = traits, aes(x=Syndrome, y=Style_length)) + 
  geom_dotplot(binaxis='y', stackdir='center') +
  geom_text_repel(aes(label = traits$sp_tip_label), box.padding = unit(1,"lines"))
styl_len
dev.off()

pdf(file="bract_color.pdf", width = 10, height = 10, pointsize = 5)
brac_col <- ggplot(data = traits, aes(x=Syndrome, y=bract_color)) + 
  geom_jitter() +
  geom_text_repel(aes(label = traits$sp_tip_label), box.padding = unit(0.25,"lines"), max.overlaps = 10)
brac_col
dev.off()

pdf(file="corolla_color.pdf", width = 10, height = 10, pointsize = 5)
cor_col <- ggplot(data = traits, aes(x=Syndrome, y=corolla_color)) + 
  geom_jitter() +
  geom_text_repel(aes(label = traits$sp_tip_label), box.padding = unit(0.5,"lines"), max.overlaps = 15)
cor_col
dev.off()

pdf(file="labellum_color.pdf", width = 10, height = 10, pointsize = 5)
lab_col <- ggplot(data = traits, aes(x=Syndrome, y=lab_color)) + 
  geom_jitter() +
  geom_text_repel(aes(label = traits$sp_tip_label), box.padding = unit(0.5,"lines"), max.overlaps = 15)
lab_col
dev.off()

pdf(file="yellow_stripe.pdf", width = 10, height = 10, pointsize = 5)
yel_stripe <- ggplot(data = traits, aes(x=Syndrome, y=yellow_labstripe)) + 
  geom_jitter() +
  geom_text_repel(aes(label = traits$sp_tip_label), box.padding = unit(0.5,"lines"), max.overlaps = 15)
yel_stripe
dev.off()

pdf(file="red_stripe.pdf", width = 10, height = 10, pointsize = 5)
red_stripe <- ggplot(data = traits, aes(x=Syndrome, y=red_labstripe)) + 
  geom_jitter() +
  geom_text_repel(aes(label = traits$sp_tip_label), box.padding = unit(0.5,"lines"), max.overlaps = 15)
red_stripe
dev.off()

pdf(file="stamen_color.pdf", width = 10, height = 10, pointsize = 5)
stam_col <- ggplot(data = traits, aes(x=Syndrome, y=stamen_color)) + 
  geom_jitter() +
  geom_text_repel(aes(label = traits$sp_tip_label), box.padding = unit(0.5,"lines"), max.overlaps = 15)
stam_col
dev.off()

pdf(file="stamen_tip_color.pdf", width = 10, height = 10, pointsize = 5)
sttip_col <- ggplot(data = traits, aes(x=Syndrome, y=stamen_tip_color)) + 
  geom_jitter() +
  geom_text_repel(aes(label = traits$sp_tip_label), box.padding = unit(0.5,"lines"), max.overlaps = 15)
sttip_col
dev.off()

##################################################
#now do analysis of mixed data

#first impute missing data
#how many missing values are there?
sum(is.na(traits))
#35 out of 884 measurements missing; all missing values are in continuous traits
#there are 10 continuous traits and 7 categorical traits

#find optimal number of components for imputation of missing data
nc <- estim_ncpFAMD(traits, method.cv = "Kfold", verbose = TRUE, sup.var = c(1,19))
nc
#5

#impute missing data
df <- imputeFAMD(traits, ncp = 5, sup.var = c(1,19))
write_csv(df$completeObs, "imputed.csv")

#now perform factor analysis of mixed data on the imputed dataset
res.famd <- FAMD(df$completeObs, 
                 sup.var = c(1,19),  ## Set the tip label and syndrome so it is not included in the analysis for now
                 graph = TRUE, 
                 ncp=10) #only output the first 10 components
#you will get graphs of the individuals, continuous variables, categorical variables, and all 
summary(res.famd)
inds <- as.data.frame(res.famd$ind$coord)
famd_output <- bind_cols(traits$sp_tip_label, traits$Syndrome, inds)
#output coordinates to be used as new multivariate trait values
write.csv(famd_output,"famd_coords.csv")
loadings <- bind_rows(as.data.frame(res.famd$quali.var$coord), as.data.frame(res.famd$quanti.var$coord))
#output loadings to better understand how traits contribute to dimensions
write.csv(loadings, "famd_loadings.csv")

#use different program to visualize the results more nicely
#visualize categorical variables
fviz_famd_var(res.famd)
#visualize data points with ellipses
pdf(file = "color_biplot.pdf", width = 10, height = 10, pointsize = 5)
fviz_famd_ind(res.famd, 
              axes = c(1, 2),
              geom = c("point", "text"),
              repel = TRUE,
              habillage = 19,
              addEllipses = TRUE) #this gives you ellipses based on expert prediction of syndrome
dev.off()


####################################################################

#all below is what I used to check correspondence between greenhouse and Maas data

#first we need to import the data, without alterations
#this copies the code above but does not adjust labellum and stamen lengths

setwd("~/Dropbox/Dena and Kathleen/FloralEvolution/Floral_traits")

Maas <- read_sheet(ss="1O5E5VY5EJToUudhQ9sqEO-RorOEeElT8Ocyg9jkldRE", sheet="Maas_traits")
Greenhouse <- read_sheet(ss="1P0NNeONfMfebcPupPFTaNXSRFt5iNs2GNtHzep4GB3A", sheet="Greenhouse_traits")

#Maas sometimes reported max and min values for each trait. When he did, we took the midpoint. 
#When he only reported one value, we used that.
#subset the data for just the columns we want
Maas <- select(Maas, !contains("Max"))
Maas <- select(Maas, !contains("Min"))

#subset Maas data set getting rid of stigma exsertion and anther exsertion
Maas <- Maas %>% select(Species_Maas, sp_tip_label, Syndrome, Stamen_Length_mid,	Stamen_Width_mid,
                        Corolla_Length_mid,	Corolla_Tube_Length_mid,	Corolla_Lobe_Length_mid,
                        Labellum_Length_mid,	Labellum_Width_mid,	Anther_Length_mid, bract_color_K,
                        corolla_color_K,	lab_color_K,	yellow_labstripe,	red_labstripe,	
                        stamen_color_K,	stamen_tip_color,	nectar_guides, Style_length_mid)

#make sure factors are factors
Maas$bract_color_K <- as.factor(Maas$bract_color_K)
Maas$corolla_color_K <- as.factor(Maas$corolla_color_K)
Maas$lab_color_K <- as.factor(Maas$lab_color_K)
Maas$yellow_labstripe <- as.factor(Maas$yellow_labstripe)
Maas$red_labstripe <- as.factor(Maas$red_labstripe)
Maas$stamen_color_K <- as.factor(Maas$stamen_color_K)
Maas$stamen_tip_color <- as.factor(Maas$stamen_tip_color)
Maas$nectar_guides <- as.factor(Maas$nectar_guides)

#subset greenhouse data set getting rid of unnecessary columns
Greenhouse <- Greenhouse %>% select(Greenhouse_ID, tip_label,	sp_tip_label,	Syndrome,	Corolla_Length,
                                    Corolla_Lobe_Length,	Corolla_Tube_Length,	Stamen_exsertion,	Labellum_Length,
                                    Labellum_Width,	Stamen_Length,
                                    Stamen_Width,	Anther_Length, Style_length, bract_color_K,
                                    corolla_color_K,	lab_color_K,	yellow_labstripe,	red_labstripe,	
                                    stamen_color_K,	stamen_tip_color,	nectar_guides)

#make sure factors are factors
Greenhouse$bract_color_K <- as.factor(Greenhouse$bract_color_K)
Greenhouse$corolla_color_K <- as.factor(Greenhouse$corolla_color_K)
Greenhouse$lab_color_K <- as.factor(Greenhouse$lab_color_K)
Greenhouse$yellow_labstripe <- as.factor(Greenhouse$yellow_labstripe)
Greenhouse$red_labstripe <- as.factor(Greenhouse$red_labstripe)
Greenhouse$stamen_color_K <- as.factor(Greenhouse$stamen_color_K)
Greenhouse$stamen_tip_color <- as.factor(Greenhouse$stamen_tip_color)
Greenhouse$nectar_guides <- as.factor(Greenhouse$nectar_guides)


#summarize (average or mode, depending on whether the data is continuous or categorical)
#data within greenhouse and Maas data frames
#first take average/mode of individuals measured multiple times
#then take average/mode of different individuals from the same population
#then take average/mode of different populations of the same species

#make a function to calculate the mode of categorical variables
Mode <- function(x) {
  uniqx <- unique(x)
  uniqx[which.max(tabulate(match(x, uniqx)))]
}

#summarize continuous traits by individual by taking the mean
#also add "_mid" to the trait names to make them comparable to the Maas data
Greenhouse_by_ind <- Greenhouse %>% group_by(Greenhouse_ID) %>% summarise(
  tibble(
    across(where(is.double), mean, na.rm = TRUE, .names = "{.col}_mid"),
    across(where(is.factor), unique)
  )
)
Greenhouse_by_ind <- as.data.frame(Greenhouse_by_ind)


#now take the mode of categorical variables
Greenhouse_by_ind <- Greenhouse_by_ind %>% group_by(Greenhouse_ID) %>% mutate(bract_color_K=Mode(bract_color_K))
Greenhouse_by_ind <- Greenhouse_by_ind %>% group_by(Greenhouse_ID) %>% mutate(corolla_color_K=Mode(corolla_color_K))
Greenhouse_by_ind <- Greenhouse_by_ind %>% group_by(Greenhouse_ID) %>% mutate(lab_color_K=Mode(lab_color_K))
Greenhouse_by_ind <- Greenhouse_by_ind %>% group_by(Greenhouse_ID) %>% mutate(yellow_labstripe=Mode(yellow_labstripe))
Greenhouse_by_ind <- Greenhouse_by_ind %>% group_by(Greenhouse_ID) %>% mutate(red_labstripe=Mode(red_labstripe))
Greenhouse_by_ind <- Greenhouse_by_ind %>% group_by(Greenhouse_ID) %>% mutate(stamen_color_K=Mode(stamen_color_K))
Greenhouse_by_ind <- Greenhouse_by_ind %>% group_by(Greenhouse_ID) %>% mutate(stamen_tip_color=Mode(stamen_tip_color))
Greenhouse_by_ind <- Greenhouse_by_ind %>% group_by(Greenhouse_ID) %>% mutate(nectar_guides=Mode(nectar_guides))

#now remove duplicates by ind
Greenhouse_by_ind <- Greenhouse_by_ind %>% group_by(Greenhouse_ID) %>% distinct()

#add back in columns that were deleted
tip_labels <- dplyr::select(Greenhouse, Greenhouse_ID,tip_label, sp_tip_label, Syndrome)
tip_labels <- distinct(tip_labels)

Greenhouse_by_ind <- full_join(Greenhouse_by_ind, tip_labels, by = "Greenhouse_ID")


#summarize by population
#first take the mean of continuous traits
Greenhouse_by_pop <- Greenhouse_by_ind %>% group_by(tip_label) %>% summarise(
  tibble(
    across(where(is.double), mean, na.rm = TRUE),
    across(where(is.factor), unique)
  )
)
Greenhouse_by_pop <- as.data.frame(Greenhouse_by_pop)

#now take the mode of categorical variables
Greenhouse_by_pop <- Greenhouse_by_pop %>% group_by(tip_label) %>% mutate(bract_color_K=Mode(bract_color_K))
Greenhouse_by_pop <- Greenhouse_by_pop %>% group_by(tip_label) %>% mutate(corolla_color_K=Mode(corolla_color_K))
Greenhouse_by_pop <- Greenhouse_by_pop %>% group_by(tip_label) %>% mutate(lab_color_K=Mode(lab_color_K))
Greenhouse_by_pop <- Greenhouse_by_pop %>% group_by(tip_label) %>% mutate(yellow_labstripe=Mode(yellow_labstripe))
Greenhouse_by_pop <- Greenhouse_by_pop %>% group_by(tip_label) %>% mutate(red_labstripe=Mode(red_labstripe))
Greenhouse_by_pop <- Greenhouse_by_pop %>% group_by(tip_label) %>% mutate(stamen_color_K=Mode(stamen_color_K))
Greenhouse_by_pop <- Greenhouse_by_pop %>% group_by(tip_label) %>% mutate(stamen_tip_color=Mode(stamen_tip_color))
Greenhouse_by_pop <- Greenhouse_by_pop %>% group_by(tip_label) %>% mutate(nectar_guides=Mode(nectar_guides))

#now remove duplicates by pop
Greenhouse_by_pop <- Greenhouse_by_pop %>% group_by(tip_label) %>% distinct()

#add back in columns that were deleted
tip_labels <- dplyr::select(Greenhouse,tip_label, sp_tip_label, Syndrome)
tip_labels <- distinct(tip_labels)

Greenhouse_by_pop <- full_join(Greenhouse_by_pop, tip_labels, by = "tip_label")

#now summarise by species tip label to get means for each tip in the species tree

Greenhouse_by_species <- Greenhouse_by_pop %>% group_by(sp_tip_label) %>% summarise(
  tibble(
    across(where(is.double), mean, na.rm=TRUE),
    across(where(is.factor), unique)
  )
)

#take mode for polymorphic species
Greenhouse_by_species <- Greenhouse_by_species %>% group_by(sp_tip_label) %>% mutate(bract_color_K=Mode(bract_color_K))
Greenhouse_by_species <- Greenhouse_by_species %>% group_by(sp_tip_label) %>% mutate(corolla_color_K=Mode(corolla_color_K))
Greenhouse_by_species <- Greenhouse_by_species %>% group_by(sp_tip_label) %>% mutate(lab_color_K=Mode(lab_color_K))
Greenhouse_by_species <- Greenhouse_by_species %>% group_by(sp_tip_label) %>% mutate(yellow_labstripe=Mode(yellow_labstripe))
Greenhouse_by_species <- Greenhouse_by_species %>% group_by(sp_tip_label) %>% mutate(red_labstripe=Mode(red_labstripe))
Greenhouse_by_species <- Greenhouse_by_species %>% group_by(sp_tip_label) %>% mutate(stamen_color_K=Mode(stamen_color_K))
Greenhouse_by_species <- Greenhouse_by_species %>% group_by(sp_tip_label) %>% mutate(stamen_tip_color=Mode(stamen_tip_color))
Greenhouse_by_species <- Greenhouse_by_species %>% group_by(sp_tip_label) %>% mutate(nectar_guides=Mode(nectar_guides))

#now remove duplicate rows
Greenhouse_by_species <- Greenhouse_by_species %>% group_by(sp_tip_label) %>% distinct()

#remove rows for samples not in phylogeny
Greenhouse_by_species <- filter(Greenhouse_by_species, !is.na(sp_tip_label))
#remove species with poor data
Greenhouse_by_species <- Greenhouse_by_species %>% filter(sp_tip_label != "Costus_sp_19207", sp_tip_label != "Costus_sp_19113" )

#add the Syndrome back in
sp_tip_labels <- dplyr::select(Greenhouse,sp_tip_label, Syndrome)
sp_tip_labels <- distinct(sp_tip_labels)
Greenhouse_by_species <- inner_join(Greenhouse_by_species, sp_tip_labels, by = "sp_tip_label")

#use inner join to isolate species where we have overlapping data

overlap <- inner_join(Greenhouse_by_species, Maas, by = "sp_tip_label")

overlap_numbers <- select(overlap, where(is.double))

#pairwise correlations
df<-overlap_numbers

#make a dataframe of correlations
correlations <- data.frame(trait = character(), estimate = double(), p.value=double(), stringsAsFactors=FALSE)

corolla_length <- cor.test(df$Corolla_Length_mid.x,df$Corolla_Length_mid.y, use= "complete.obs")
plot(df$Corolla_Length_mid.x~df$Corolla_Length_mid.y)
abline(0,1)
abline(lm(df$Corolla_Length_mid.x~df$Corolla_Length_mid.y), col = "red")  
#corolla_length looks good r = 0.811 falls right on 1:1 line
correlations[nrow(correlations) + 1,] = 
  c("corolla_length", corolla_length$estimate, corolla_length$p.value)

corolla_lobe_length <- cor.test(df$Corolla_Lobe_Length_mid.x,df$Corolla_Lobe_Length_mid.y, use= "complete.obs")
plot(df$Corolla_Lobe_Length_mid.x~df$Corolla_Lobe_Length_mid.y)
abline(0,1)
abline(lm(df$Corolla_Lobe_Length_mid.x~df$Corolla_Lobe_Length_mid.y), col = "red")  
#corolla_lobe_length looks good r = 0.843 falls higher than 1:1 line
correlations[nrow(correlations) + 1,] = 
  c("corolla_lobe_length", corolla_lobe_length$estimate, corolla_lobe_length$p.value)

corolla_tube_length <- cor.test(df$Corolla_Tube_Length_mid.x,df$Corolla_Tube_Length_mid.y, use= "complete.obs")
plot(df$Corolla_Tube_Length_mid.x~df$Corolla_Tube_Length_mid.y)
abline(0,1)
abline(lm(df$Corolla_Tube_Length_mid.x~df$Corolla_Tube_Length_mid.y), col = "red")  
#corolla_tube_length r = 0.658 right on 1:1 line
correlations[nrow(correlations) + 1,] = 
  c("corolla_tube_length", corolla_tube_length$estimate, corolla_tube_length$p.value)

Labellum_length <- cor.test(df$Labellum_Length_mid.x,df$Labellum_Length_mid.y, use= "complete.obs")
plot(df$Labellum_Length_mid.x~df$Labellum_Length_mid.y)
abline(0,1)
abline(lm(df$Labellum_Length_mid.x~df$Labellum_Length_mid.y), col = "red")  
#labellum length is r = 0.881 and parallels the 1:1 line, but my measures are consistently higher (by nearly 20mm)
correlations[nrow(correlations) + 1,] = 
  c("Labellum_length", Labellum_length$estimate, Labellum_length$p.value)

Labellum_Width <- cor.test(df$Labellum_Width_mid.x,df$Labellum_Width_mid.y, use= "complete.obs")
plot(df$Labellum_Width_mid.x~df$Labellum_Width_mid.y)
abline(0,1)
abline(lm(df$Labellum_Width_mid.x~df$Labellum_Width_mid.y), col = "red")  
#labellum width is r = 0.883 and close to 1:1, although at high values, Maas' measures are larger
correlations[nrow(correlations) + 1,] = 
  c("Labellum_Width", Labellum_Width$estimate, Labellum_Width$p.value)

Stamen_length <- cor.test(df$Stamen_Length_mid.x,df$Stamen_Length_mid.y, use= "complete.obs")
plot(df$Stamen_Length_mid.x~df$Stamen_Length_mid.y)
abline(0,1)
abline(lm(df$Stamen_Length_mid.x~df$Stamen_Length_mid.y), col = "red")  
#stamen length has high correlation 0.813, but my measures are consistently higher by 10-20 mm
correlations[nrow(correlations) + 1,] = 
  c("Stamen_length", Stamen_length$estimate, Stamen_length$p.value)

Stamen_Width <- cor.test(df$Stamen_Width_mid.x,df$Stamen_Width_mid.y, use= "complete.obs")
plot(df$Stamen_Width_mid.x~df$Stamen_Width_mid.y)
abline(0,1)
abline(lm(df$Stamen_Width_mid.x~df$Stamen_Width_mid.y), col = "red")  
#stamen width high correlation 0.761 and close to 1:1 line
correlations[nrow(correlations) + 1,] = 
  c("Stamen_Width", Stamen_Width$estimate, Stamen_Width$p.value)

Anther_length <- cor.test(df$Anther_Length_mid.x,df$Anther_Length_mid.y, use= "complete.obs")
plot(df$Anther_Length_mid.x~df$Anther_Length_mid.y)
abline(0,1)
abline(lm(df$Anther_Length_mid.x~df$Anther_Length_mid.y), col = "red")  
#anther length r = 0.776 close to 1:1
correlations[nrow(correlations) + 1,] = 
  c("Anther_length", Anther_length$estimate, Anther_length$p.value)

write_csv(correlations, "correlations.csv")
#all categorical traits agree


