#Ryan Stephens; Finalized May 8, 2021
rm(list=ls()) # clears workspace
setwd("C:/Users/ryans/Documents/Isotopic_Routing_Small_Mammals/TDF_Lit_Review")

##################################
#loads files into R
##################################
TDF<- read.csv("Literature_Review_TDF_data.csv",header=T)#Literature review data
##################################


################################################################################################
#Get general data set ready
################################################################################################
library(dplyr)
#Diet SD (if not SD, take SE and calculate SD using sample size of diet samples analyzed)
TDF<- TDF %>%  mutate(SD_Diet_13C_Final=ifelse(Diet_d13C_SD_SE== "SD",Diet_d13C_Var,Diet_d13C_Var*sqrt(TDF$Diet_n)))#d13C
TDF<- TDF %>%  mutate(SD_Diet_15N_Final=ifelse(Diet_d15N_SD_SE== "SD",Diet_d15N_Var,Diet_d15N_Var*sqrt(TDF$Diet_n)))#d15N

#TDF SD (if not SD, take SE and calculate SD using sample size of animals)
TDF<- TDF %>%  mutate(SD_TDF_C13_Final=ifelse(TDF_d13C_SD_SE== "SD",TDF_13C_Var,TDF_13C_Var*sqrt(TDF$Animal_n)))#d13C
TDF<- TDF %>%  mutate(SD_TDF_15N_Final=ifelse(TDF_d15N_SD_SE== "SD",TDF_15N_Var,TDF_15N_Var*sqrt(TDF$Animal_n)))#d15N

TDF<- mutate(TDF,Genus_species = paste(Genus,Species,sep = '_'))#makes new column for genus-species
TDF$Genus_species<-as.factor(TDF$Genus_species)

#Remove data from field study that was used in second part of paper to test TDFs
TDF<-TDF %>% filter(Reference != "Stephens et al. (2021)")#TDF used to test the efficacy of field studies

#Remove these blood measures that are not used
TDF<-TDF %>% filter(Tissue != "Blood Cells" & 
                  Tissue != "Blood Plasma" &
                  Tissue != "Blood Serum" & 
                  Tissue != "Red Blood Cells")

library(forcats)#use forcats to rename levels
TDF<-mutate(TDF, Tissue = fct_recode(Tissue, "Blood" = "Whole Blood"))#rename "Whole Blood" to "Blood"

#combine reported and estimated C/N to make one C/N column
TDF<- TDF %>% mutate(Diet_C.N = ifelse(is.na(Reported_diet_C.N), Estimated_diet_C.N, Reported_diet_C.N))

#total number of studies
TDF%>%summarize(Total_studies = n_distinct(Reference))
TDF <- droplevels(TDF)#drop unused levels
################################################################################################




################################################################################################
#Carbon
################################################################################################
#If studies didn't use hair, then TDFs from another tissue type was converted to hair
Carbon<-TDF %>% filter(TDF_13C != "NA")#remove studies that did not record d13C
Carbon<-Carbon %>% filter(Lipid_corrected_diet == "N")#Limit to diets that were not lipid extracted
Carbon<- mutate(Carbon,Diet_Unique = paste(Reference,Diet,Genus_species,sep = '-'))#unique ID for each diet
Carbon_data<-select(Carbon, Reference, Diet_Unique, Study, Genus_species, Diet, Diet_Source,#selected columns 
                    Diet_isotope = Diet_d13C, Diet_C.N,Consumer_type = Class, TDF = TDF_13C, Tissue)
head(Carbon_data) 

library(tidyr)
C_Tissue_wide <- spread(Carbon_data, Tissue, TDF)#put tissue type into different columns  
head(C_Tissue_wide)

library(GGally)
ggpairs(C_Tissue_wide, columns = c(10,11,12,13,15))#test plot for correlations

################################################
#Hair
################################################
C_Hair<-filter(C_Tissue_wide,Hair != "NA")#filter out samples with TDF values for hair
C_Hair<-select(C_Hair, Reference:Consumer_type, TDF = Hair)#Select hair and make it the TDF column
C_Hair$Original_Tissue <- "Hair"#add column for the original tissue type
head(C_Hair)
str(C_Hair)
################################################


################################################
#Muscle to hair TDF conversion
################################################
#get offset
C_Hair_Muscle<-filter(C_Tissue_wide,Hair != "NA" & Muscle != "NA")#filter only diets with both hair and muscle
C_Hair_Muscle<-mutate(C_Hair_Muscle, difference = Hair - Muscle)#offset between hair and muscle
C_mean_Hair_muscle_Offset<-round(mean(C_Hair_Muscle$difference), digits = 2)#mean difference 
print(C_mean_Hair_muscle_Offset)
#adjust values
C_Muscle<-filter(C_Tissue_wide,is.na(Hair) & Muscle != "NA")#select muscle TDF values when hair has none
C_Muscle<-mutate(C_Muscle, TDF = Muscle + C_mean_Hair_muscle_Offset)#add offset to adjust it to hair
C_Muscle$Original_Tissue <- "Muscle"#add "original tissue" column and fill with "muscle"
C_Muscle<-select(C_Muscle, Reference:Consumer_type, TDF:Original_Tissue)#select needed columns
head(C_Muscle)
str(C_Muscle)
################################################


################################################
#Collagen to hair TDF conversion
################################################
#get offset
C_Hair_Collagen<-filter(C_Tissue_wide,Hair != "NA" & Collagen != "NA")#filter only diets with both hair and collagen
C_Hair_Collagen<-mutate(C_Hair_Collagen, difference = Hair - Collagen)#offset between hair and collagen
C_mean_Hair_Collagen_Offset<-round(mean(C_Hair_Collagen$difference), digits = 2)#mean difference 
print(C_mean_Hair_Collagen_Offset)
#adjust values
C_Collagen<-filter(C_Tissue_wide,is.na(Hair) & is.na(Muscle) & 
                   Collagen != "NA")#select collagen TDF values when hair and muscle have none
C_Collagen<-mutate(C_Collagen, TDF = Collagen + C_mean_Hair_Collagen_Offset)#add offset to adjust it to hair
C_Collagen$Original_Tissue <- "Collagen"#add "original tissue" column and fill with "collagen"
C_Collagen<-select(C_Collagen, Reference:Consumer_type, TDF:Original_Tissue)#select needed columns
head(C_Collagen)
str(C_Collagen)
################################################


################################################
#Blood to hair TDF conversion
################################################
#get offset
C_Hair_Blood<-filter(C_Tissue_wide,Hair != "NA" & Blood != "NA")#filter only diets with both hair and blood
C_Hair_Blood<-mutate(C_Hair_Blood, difference = Hair - Blood)#offset between hair and blood
C_mean_Hair_Blood_Offset<-round(mean(C_Hair_Blood$difference), digits = 2)#mean difference 
print(C_mean_Hair_Blood_Offset)
#adjust values
C_Blood<-filter(C_Tissue_wide,is.na(Hair) & is.na(Muscle) & is.na(Collagen) &
                     Blood != "NA")#select blood TDF values when hair, muscle, and collagen have none
C_Blood<-mutate(C_Blood, TDF = Blood + C_mean_Hair_Blood_Offset)#add offset to adjust it to hair
C_Blood$Original_Tissue <- "Blood"#add "original tissue" column and fill with "blood"
C_Blood<-select(C_Blood, Reference:Consumer_type, TDF:Original_Tissue)#select needed columns
head(C_Blood)
str(C_Blood)
################################################


################################################
#Skin to hair TDF conversion
################################################
#skin of marine mammals has similar isotopic values to muscle and we used the difference between hair 
#and muscle as the offset between skin and hair
#adjust values
C_Skin<-filter(C_Tissue_wide,is.na(Hair) & is.na(Muscle) & is.na(Collagen) & is.na(Blood) &
                  Skin != "NA")#select skin TDF values when hair, muscle, collagen, and blood have none
C_Skin<-mutate(C_Skin, TDF = Skin + C_mean_Hair_muscle_Offset)#add offset to adjust it to hair (muscle)
C_Skin$Original_Tissue <- "Skin"#add "original tissue" column and fill with "skin"
C_Skin<-select(C_Skin, Reference:Consumer_type, TDF:Original_Tissue)#select needed columns
head(C_Skin)
str(C_Skin)
################################################



################################################
#Liver to hair TDF conversion
################################################
#get offset
C_Hair_Liver<-filter(C_Tissue_wide,Hair != "NA" & Liver != "NA")#filter only diets with both hair and liver
C_Hair_Liver<-mutate(C_Hair_Liver, difference = Hair - Liver)#offset between hair and liver
C_mean_Hair_Liver_Offset<-round(mean(C_Hair_Liver$difference), digits = 2)#mean difference 
print(C_mean_Hair_Liver_Offset)

#adjust values
C_Liver<-filter(C_Tissue_wide,is.na(Hair) & is.na(Muscle) & is.na(Collagen) & is.na(Blood) & is.na(Skin) &
                 Liver != "NA")#select liver TDF values when hair, muscle, collagen, blood, and skin have none
str(C_Liver)#none
################################################


Carbon_adjusted<-rbind(C_Hair, C_Muscle, C_Collagen, C_Blood, C_Skin)#join hair TDF to those estimated from other tissues
Carbon_adjusted<- mutate(Carbon_adjusted,Diet_Unique_tissue = paste(Diet_Unique, Original_Tissue ,sep = '-'))#unique ID for each diet and tissue
head(Carbon_adjusted)
str(Carbon_adjusted)

#Add in SD and sample size
Carbon<- mutate(Carbon,Diet_Unique_tissue = paste(Diet_Unique, Tissue ,sep = '-'))#unique ID for each diet and tissue
SD<-select(Carbon, Diet_Unique_tissue, SD_Diet_13C_Final, SD_TDF_C13_Final, Animal_n)#SD data
Carbon_adjusted_sd<-left_join(Carbon_adjusted,SD,by="Diet_Unique_tissue")#join SD data to tdf data
Carbon_adjusted_sd<-rename(Carbon_adjusted_sd, SD_Diet_isotope = SD_Diet_13C_Final, SD_TDF = SD_TDF_C13_Final)
Carbon_adjusted_sd$Isotope <- "d13C"
head(Carbon_adjusted_sd)
str(Carbon_adjusted_sd)
################################################################################################





################################################################################################
#Nitrogen
################################################################################################
Nitrogen<-TDF %>% filter(TDF_15N != "NA")#remove studies that did not record d13C
#Nitrogen<-Nitrogen %>% filter(Lipid_corrected_diet == "N")#Limit to diets that were not lipid extracted
Nitrogen<- mutate(Nitrogen,Diet_Unique = paste(Reference,Diet, Genus_species,sep = '-'))#unique ID for each diet
Nitrogen_data<-select(Nitrogen, Reference, Diet_Unique, Study, Genus_species, Diet, Diet_Source,#selected columns
                      Diet_isotope = Diet_d15N, Diet_C.N,Consumer_type = Class, TDF = TDF_15N, Tissue)
head(Nitrogen_data) 

library(tidyr)
N_Tissue_wide <- spread(Nitrogen_data, Tissue, TDF)#put tissue type into different columns  
head(N_Tissue_wide)

library(GGally)
ggpairs(N_Tissue_wide, columns = c(10,11,12,13,15))#hair and blood are not significantly correlated

################################################
#Hair
################################################
N_Hair<-filter(N_Tissue_wide,Hair != "NA")#filter out samples with TDF values for hair
N_Hair<-select(N_Hair, Reference:Consumer_type, TDF = Hair)#Select hair and make it the TDF column
N_Hair$Original_Tissue <- "Hair"#add column for the original tissue type
head(N_Hair)
str(N_Hair)
################################################


################################################
#Muscle to hair TDF conversion
################################################
#adjust values
N_Muscle<-filter(N_Tissue_wide,is.na(Hair) & Muscle != "NA")#select muscle TDF values when hair has none
N_Muscle<-mutate(N_Muscle, TDF = Muscle)#use muscle TDF for hair TDF
N_Muscle$Original_Tissue <- "Muscle"#add "original tissue" column and fill with "muscle"
N_Muscle<-select(N_Muscle, Reference:Consumer_type, TDF:Original_Tissue)#select needed columns
head(N_Muscle)
str(N_Muscle)
################################################


################################################
#Collagen to hair TDF conversion
################################################
#adjust values
N_Collagen<-filter(N_Tissue_wide,is.na(Hair) & is.na(Muscle) & 
                     Collagen != "NA")#select collagen TDF values when hair and muscle have none
str(N_Collagen)#none
################################################


################################################
#Skin to hair TDF conversion
################################################
#skin of marine mammals has similar isotopic values to muscle and we used the difference between hair 
#and muscle as the offset between skin and hair (which is 0 for d15N)
#adjust values
N_Skin<-filter(N_Tissue_wide,is.na(Hair) & is.na(Muscle) &
                 Skin != "NA")#select skin TDF values when hair, muscle, and collagen have none
N_Skin<-mutate(N_Skin, TDF = Skin)#use skin TDF for hair TDF
N_Skin$Original_Tissue <- "Skin"#add "original tissue" column and fill with "skin"
N_Skin<-select(N_Skin, Reference:Consumer_type, TDF:Original_Tissue)#select needed columns
head(N_Skin)
str(N_Skin)
################################################


################################################
#Liver to hair TDF conversion
################################################
#adjust values
N_Liver<-filter(N_Tissue_wide,is.na(Hair) & is.na(Muscle) & is.na(Skin) &
                  Liver != "NA")#select liver TDF values when hair, muscle, collagen, blood, and skin have none
N_Liver<-mutate(N_Liver, TDF = Liver)#use liver TDF for hair TDF
N_Liver$Original_Tissue <- "Liver"#add "original tissue" column and fill with "liver"
N_Liver<-select(N_Liver, Reference:Consumer_type, TDF:Original_Tissue)#select needed columns
head(N_Liver)
str(N_Liver)
################################################


Nitrogen_adjusted<-rbind(N_Hair, N_Muscle, N_Skin, N_Liver)#join hair TDF to hair TDF estimated from other tissues
Nitrogen_adjusted<- mutate(Nitrogen_adjusted,Diet_Unique_tissue = paste(Diet_Unique, Original_Tissue ,sep = '-'))#unique ID for each diet and tissue
head(Nitrogen_adjusted)
str(Nitrogen_adjusted)
#Add in SD and sample size
Nitrogen<- mutate(Nitrogen,Diet_Unique_tissue = paste(Diet_Unique, Tissue ,sep = '-'))#unique ID for each diet and tissue
SD<-select(Nitrogen, Diet_Unique_tissue, SD_Diet_15N_Final, SD_TDF_15N_Final, Animal_n)#SD data

Nitrogen_adjusted_sd<-left_join(Nitrogen_adjusted,SD,by="Diet_Unique_tissue")#join SD data to TDF data

Nitrogen_adjusted_sd<-rename(Nitrogen_adjusted_sd, SD_Diet_isotope = SD_Diet_15N_Final, SD_TDF = SD_TDF_15N_Final)
Nitrogen_adjusted_sd$Isotope <- "d15N"
head(Nitrogen_adjusted_sd)


Nitrogen_adjusted_sd<-Nitrogen_adjusted_sd %>% filter(Diet_Unique != "Codron et al. (2012)-Natural-Damaliscus_pygargus" &#Remove field study where there may have been a diet switch
                                                        Diet_Unique != "Codron et al. (2012)-Natural-Antidorcas_marsupialis")#Remove field study where TDF of d15N TDF was very large (11.6)
str(Nitrogen_adjusted_sd)
################################################################################################




################################################################################################
#Join datasets 
################################################################################################
TDF_C_N_data<-rbind(Carbon_adjusted_sd, Nitrogen_adjusted_sd)
head(TDF_C_N_data)
#Summary
TDF_C_N_data%>%
  group_by(Isotope)%>%
  summarise(Studies = n_distinct(Reference), Species = n_distinct(Genus_species), 
            Mean_TDF = mean(TDF), sd_TDF = sd(TDF), min_TDF = min(TDF), max_TDF = max(TDF), TDF = n())

TDF_C_N_data%>%
  group_by(Isotope, Original_Tissue)%>%
  summarise(TDF = n())

TDF_C_N_data<- select(TDF_C_N_data, -Diet_Unique, -Diet_Unique_tissue)

write.csv(TDF_C_N_data,"MetaAnalysis_TDF_data.csv",row.names = F)#write file to csv
################################################################################################




################################################################################################
#Carbon TDF tissue comparison table
################################################################################################
#Make table with mean offsets between tissue types, along with their correlation values and t-tests
Carbon<-TDF %>% filter(TDF_13C != "NA")#remove studies that did not record d13C
Carbon<-Carbon %>% filter(Lipid_corrected_diet == "N")#Limit to diets that were not lipid extracted
Carbon<- mutate(Carbon,Diet_Unique = paste(Reference,Diet_d13C,Genus_species,sep = '-'))#unique ID for each diet
Carbon_data<-select(Carbon, Reference, Diet_Unique, Study, Genus_species, Diet_Source, Diet_isotope = Diet_d13C, Diet_C.N,#selected columns
                    Consumer_type = Class, TDF = TDF_13C, Original_Tissue = Tissue)
library(tidyr)
Carbon_Tissue_Wide <- spread(Carbon_data, Original_Tissue, TDF)#put tissues into different columns  

C_Tissues<-select(Carbon_Tissue_Wide, Diet_Source, Consumer_type, Collagen, Hair, Liver, Muscle, Blood)
head(C_Tissues)

library(GGally)
ggpairs(C_Tissues, columns = 3:7)#test graph

#Hair & Collagen
C_Hair_Collagen<-filter(C_Tissues,Hair != "NA" & Collagen != "NA")#filter only diets with both hair and Collagen
C_Hair_Collagen<-select(C_Hair_Collagen, Diet_Source, Consumer_type, Tissue_1 = Hair, Tissue_2 = Collagen)
C_Hair_Collagen$Comparison<-"Hair - Collagen"
#Hair & Muscle
C_Hair_Muscle<-filter(C_Tissues,Hair != "NA" & Muscle != "NA")#filter only diets with both hair and Muscle
C_Hair_Muscle<-select(C_Hair_Muscle, Diet_Source, Consumer_type, Tissue_1 = Hair, Tissue_2 = Muscle)
C_Hair_Muscle$Comparison<-"Hair - Muscle"
#Hair & Liver
C_Hair_Liver<-filter(C_Tissues,Hair != "NA" & Liver != "NA")#filter only diets with both hair and Liver
C_Hair_Liver<-select(C_Hair_Liver, Diet_Source, Consumer_type, Tissue_1 = Hair, Tissue_2 = Liver)
C_Hair_Liver$Comparison<-"Hair - Liver"
#Hair & Blood
C_Hair_Blood<-filter(C_Tissues,Hair != "NA" & Blood != "NaN")#filter only diets with both hair and Blood
C_Hair_Blood<-select(C_Hair_Blood, Diet_Source, Consumer_type, Tissue_1 = Hair, Tissue_2 = Blood)
C_Hair_Blood$Comparison<-"Hair - Blood"
#Muscle & Liver
C_Muscle_Liver<-filter(C_Tissues,Muscle != "NA" & Liver != "NA")#filter only diets with both Muscle and Liver
C_Muscle_Liver<-select(C_Muscle_Liver, Diet_Source, Consumer_type, Tissue_1 = Muscle, Tissue_2 = Liver)
C_Muscle_Liver$Comparison<-"Muscle - Liver"
#Muscle & Blood
C_Muscle_Blood<-filter(C_Tissues,Muscle != "NA" & Blood != "NaN")#filter only diets with both Muscle and Blood
C_Muscle_Blood<-select(C_Muscle_Blood, Diet_Source, Consumer_type, Tissue_1 = Muscle, Tissue_2 = Blood)
C_Muscle_Blood$Comparison<-"Muscle - Blood"
#Muscle & Collagen
C_Muscle_Collagen<-filter(C_Tissues,Muscle != "NA" & Collagen != "NaN")#filter only diets with both muscle and collagen
C_Muscle_Collagen<-select(C_Muscle_Collagen, Diet_Source, Consumer_type, Tissue_1 = Muscle, Tissue_2 = Collagen)
C_Muscle_Collagen$Comparison<-"Muscle - Collagen"
#Liver & Blood
C_Liver_Blood<-filter(C_Tissues,Liver != "NA" & Blood != "NaN")#filter only diets with both Liver and Blood
C_Liver_Blood<-select(C_Liver_Blood, Diet_Source, Consumer_type, Tissue_1 = Liver, Tissue_2 = Blood)
C_Liver_Blood$Comparison<-"Liver - Blood"
#Liver & Collagen
C_Liver_Collagen<-filter(C_Tissues,Liver != "NA" & Collagen != "NaN")#filter only diets with both Liver and collagen
C_Liver_Collagen<-select(C_Liver_Collagen, Diet_Source, Consumer_type, Tissue_1 = Liver, Tissue_2 = Collagen)
C_Liver_Collagen$Comparison<-"Liver - Collagen"


C_All_Tissue_Comparisons<-rbind(C_Hair_Collagen, C_Hair_Muscle, C_Hair_Liver,#bind together
                                C_Muscle_Liver, C_Hair_Blood, C_Muscle_Blood, C_Muscle_Collagen, C_Liver_Blood, C_Liver_Collagen)
C_All_Tissue_Comparisons<- mutate(C_All_Tissue_Comparisons, Diff = Tissue_1 - Tissue_2)#difference between tissue types
head(C_All_Tissue_Comparisons)


#############################
#Mean difference and sample size
#############################
C_Mean<-C_All_Tissue_Comparisons %>% 
  group_by (Comparison) %>%#Grouping variables
  summarise(Offset = round(mean(Diff),2), sd = round(sd(Diff),2), n =n())
print(C_Mean)
#############################


#############################
#Correlations
#############################
library(broom)
C_Correlation<-C_All_Tissue_Comparisons %>% group_by(Comparison) %>% do(tidy(cor.test(.$Tissue_1, .$Tissue_2)))#correlation for comparisons
C_Correlation<-as.data.frame(C_Correlation)
C_Correlation$estimate<-round(C_Correlation$estimate, 2)#round to 2 sig digits and keep trailing zero
C_Correlation<-mutate(C_Correlation, p.value = ifelse(p.value < .0001, sprintf('"%s"',"<0.0001"),#if less than 0.00001 than 0.00001
                                                      ifelse(p.value < .001, sprintf('"%s"',"<0.001"),#if less than 0.0001 than 0.0001    
                                                             sprintf('"%.3f"',round(p.value,digits=3)))))#if greater than 0.0001 than use three digits
C_Correlation<-select(C_Correlation, Comparison, Corr = estimate, P_corr = p.value)#select needed variables
print(C_Correlation)
#############################


#############################
#t-test
#############################
library(broom)
C_t_test<-C_All_Tissue_Comparisons %>% group_by(Comparison) %>% do(tidy(t.test(.$Tissue_1, .$Tissue_2,paired=TRUE)))     
C_t_test<-as.data.frame(C_t_test)
C_t_test$statistic<-round(C_t_test$statistic, 2)#round to 2 sig digits and keep trailing zero
C_t_test<-mutate(C_t_test, p.value = ifelse(p.value < .0001, sprintf('"%s"',"<0.0001"),#if less than 0.00001 than 0.00001
                                                      ifelse(p.value < .001, sprintf('"%s"',"<0.001"),#if less than 0.0001 than 0.0001    
                                                             sprintf('"%.3f"',round(p.value,digits=3)))))#if greater than 0.0001 than use three digits
C_t_test<-select(C_t_test, Comparison, t_stat = statistic, P_t = p.value)#select needed variables
print(C_t_test)
#############################

C_Offset_corr<- left_join(C_Mean, C_Correlation,  by = "Comparison")#join mean offset and correlation
C_Offset_corr_ttest<- left_join(C_Offset_corr, C_t_test,  by = "Comparison")#join t-test
C_Offset_corr_ttest$Isotope<- "Carbon"
################################################################################################



################################################################################################
#Nitrogen TDF tissue comparison table
################################################################################################
#Make table with mean offsets between tissue types, along with their correlation values and t-tests
Nitrogen<-TDF %>% filter(TDF_13C != "NA")#remove studies that did not record d13C
#Nitrogen<-Nitrogen %>% filter(Lipid_corrected_diet == "N")#Don't run to Keep diets that were lipid extracted
Nitrogen<- mutate(Nitrogen,Diet_Unique = paste(Reference,Diet,Genus_species,sep = '-'))#unique ID for each diet
Nitrogen_data<-select(Nitrogen, Reference, Diet_Unique, Study, Genus_species, Diet_Source, Diet_isotope = Diet_d15N, Diet_C.N,#selected columns
                    Consumer_type = Class, TDF = TDF_15N, Original_Tissue = Tissue)
head(Nitrogen_data)  
str(Nitrogen_data)
Nitrogen_data<-Nitrogen_data %>% filter(Diet_Unique != "Codron et al. (2012)-Natural-Damaliscus_pygargus" &#Remove field study where there may have been a diet switch
                                                        Diet_Unique != "Codron et al. (2012)-Natural-Antidorcas_marsupialis")#Remove field study where sd of d15N TDF was very large (11.6)

Nitrogen_Tissue_Wide <- spread(Nitrogen_data, Original_Tissue, TDF)#put tissues into different columns  
head(Nitrogen_Tissue_Wide)

N_Tissues<-select(Nitrogen_Tissue_Wide, Diet_Source, Consumer_type, Collagen, Hair, Liver, Muscle, Blood)
head(N_Tissues)

library(GGally)
ggpairs(N_Tissues, columns = 4:7)#test graphs
ggpairs(N_Tissues, columns = 3:4)

#Hair & Collagen
N_Hair_Collagen<-filter(N_Tissues,Hair != "NA" & Collagen != "NA")#filter only diets with both hair and Collagen
N_Hair_Collagen<-select(N_Hair_Collagen, Diet_Source, Consumer_type, Tissue_1 = Hair, Tissue_2 = Collagen)
N_Hair_Collagen$Comparison<-"Hair - Collagen"
#Hair & Muscle
N_Hair_Muscle<-filter(N_Tissues,Hair != "NA" & Muscle != "NA")#filter only diets with both hair and Muscle
N_Hair_Muscle<-select(N_Hair_Muscle, Diet_Source, Consumer_type, Tissue_1 = Hair, Tissue_2 = Muscle)
N_Hair_Muscle$Comparison<-"Hair - Muscle"
#Hair & Liver
N_Hair_Liver<-filter(N_Tissues,Hair != "NA" & Liver != "NA")#filter only diets with both hair and Liver
N_Hair_Liver<-select(N_Hair_Liver, Diet_Source, Consumer_type, Tissue_1 = Hair, Tissue_2 = Liver)
N_Hair_Liver$Comparison<-"Hair - Liver"
#Hair & Blood
N_Hair_Blood<-filter(N_Tissues,Hair != "NA" & Blood != "NaN")#filter only diets with both hair and Blood
N_Hair_Blood<-select(N_Hair_Blood, Diet_Source, Consumer_type, Tissue_1 = Hair, Tissue_2 = Blood)
N_Hair_Blood$Comparison<-"Hair - Blood"
#Muscle & Liver
N_Muscle_Liver<-filter(N_Tissues,Muscle != "NA" & Liver != "NA")#filter only diets with both Muscle and Liver
N_Muscle_Liver<-select(N_Muscle_Liver, Diet_Source, Consumer_type, Tissue_1 = Muscle, Tissue_2 = Liver)
N_Muscle_Liver$Comparison<-"Muscle - Liver"
#Muscle & Blood
N_Muscle_Blood<-filter(N_Tissues,Muscle != "NA" & Blood != "NaN")#filter only diets with both Muscle and Blood
N_Muscle_Blood<-select(N_Muscle_Blood, Diet_Source, Consumer_type, Tissue_1 = Muscle, Tissue_2 = Blood)
N_Muscle_Blood$Comparison<-"Muscle - Blood"
#Muscle & Collagen
N_Muscle_Collagen<-filter(N_Tissues,Muscle != "NA" & Collagen != "NaN")#filter only diets with both muscle and collagen
N_Muscle_Collagen<-select(N_Muscle_Collagen, Diet_Source, Consumer_type, Tissue_1 = Muscle, Tissue_2 = Collagen)
N_Muscle_Collagen$Comparison<-"Muscle - Collagen"
#Liver & Blood
N_Liver_Blood<-filter(N_Tissues,Liver != "NA" & Blood != "NaN")#filter only diets with both Liver and Blood
N_Liver_Blood<-select(N_Liver_Blood, Diet_Source, Consumer_type, Tissue_1 = Liver, Tissue_2 = Blood)
N_Liver_Blood$Comparison<-"Liver - Blood"
#Liver & Collagen
N_Liver_Collagen<-filter(N_Tissues,Liver != "NA" & Collagen != "NaN")#filter only diets with both Liver and collagen
#only two observations


N_All_Tissue_Comparisons<-rbind(N_Hair_Collagen, N_Hair_Muscle, N_Hair_Liver,#bind together
                                N_Muscle_Liver, N_Hair_Blood, N_Muscle_Blood, N_Muscle_Collagen, N_Liver_Blood)
N_All_Tissue_Comparisons<- mutate(N_All_Tissue_Comparisons, Diff = Tissue_1 - Tissue_2)#difference between tissue types
head(N_All_Tissue_Comparisons)


#############################
#Mean difference and sample size
#############################
N_Mean<-N_All_Tissue_Comparisons %>% 
  group_by (Comparison) %>%#Grouping variables
  summarise(Offset = round(mean(Diff),2), sd = round(sd(Diff),2), n =n())
print(N_Mean)
#############################


#############################
#Correlations
#############################
library(broom)
N_Correlation<-N_All_Tissue_Comparisons %>% group_by(Comparison) %>% do(tidy(cor.test(.$Tissue_1, .$Tissue_2)))#correlation for comparisons
N_Correlation<-as.data.frame(N_Correlation)
N_Correlation$estimate<-round(N_Correlation$estimate, 2)#round to 2 sig digits and keep trailing zero
N_Correlation<-mutate(N_Correlation, p.value = ifelse(p.value < .0001, sprintf('"%s"',"<0.0001"),#if less than 0.00001 than 0.00001
                                            ifelse(p.value < .001, sprintf('"%s"',"<0.001"),#if less than 0.0001 than 0.0001    
                                                   sprintf('"%.3f"',round(p.value,digits=3)))))#if greater than 0.0001 than use three digits
N_Correlation<-select(N_Correlation, Comparison, Corr = estimate, P_corr = p.value)#select needed variables
print(N_Correlation)
#############################


#############################
#t-test
#############################
library(broom)
N_t_test<-N_All_Tissue_Comparisons %>% group_by(Comparison) %>% do(tidy(t.test(.$Tissue_1, .$Tissue_2,paired=TRUE)))     
N_t_test<-as.data.frame(N_t_test)
N_t_test$statistic<-round(N_t_test$statistic, 2)#round to 2 sig digits and keep trailing zero
N_t_test<-mutate(N_t_test, p.value = ifelse(p.value < .0001, sprintf('"%s"',"<0.0001"),#if less than 0.00001 than 0.00001
                                                      ifelse(p.value < .001, sprintf('"%s"',"<0.001"),#if less than 0.0001 than 0.0001    
                                                             sprintf('"%.3f"',round(p.value,digits=3)))))#if greater than 0.0001 than use three digits
N_t_test<-select(N_t_test, Comparison, t_stat = statistic, P_t = p.value)#select needed variables
print(N_t_test)
#############################

N_Offset_corr<- left_join(N_Mean, N_Correlation,  by = "Comparison")#join mean offset and correlation
N_Offset_corr_ttest<- left_join(N_Offset_corr, N_t_test,  by = "Comparison")#join t-test
N_Offset_corr_ttest$Isotope<- "Nitrogen"
################################################################################################




################################################################################################
#Merge tables carbon and nitrogen tables and save
################################################################################################
C_N_offsets<-rbind(C_Offset_corr_ttest, N_Offset_corr_ttest)
C_N_offsets<-mutate(C_N_offsets, sd = sprintf('%.2f',sd),#keep trailing zeros
                             t_stat = sprintf('%.2f', t_stat))
C_N_offsets
head(C_N_offsets)
str(C_N_offsets)

C_N_offsets_arranged<-C_N_offsets %>%#arrange table
          arrange(Isotope, match(Comparison, c("Hair - Collagen",
                                               "Hair - Muscle", 
                                               "Hair - Liver",
                                               "Hair - Blood",
                                               "Muscle - Liver",
                                               "Muscle - Blood",
                                               "Liver - Blood")))

C_N_offsets_Table<-select(C_N_offsets_arranged, TDF = Isotope, 'Tissue comparison' = Comparison, Offset,#rename headers
                              sd, n, t = t_stat, "P-value" = P_t, r = Corr, "P-value2" = P_corr)

write.csv(C_N_offsets_Table, "Tissue_TDF_Comparisons.csv", quote=FALSE, row.names = FALSE)#save table as csv file
################################################################################################







################################################################################################
#Carbon tissue comparison graph
################################################################################################
library(tidyr)
C_graph_data<- C_All_Tissue_Comparisons%>%#dataset from above
  mutate(Tissue_names = Comparison)%>%#duplicate column of names 
  separate(Tissue_names, c("Tissue_1_name","Tissue_2_name"), sep = "( - )")#split column into two with names
C_graph_data<-mutate(C_graph_data, Diet_Source = fct_recode(Diet_Source, "C[3]" = "C3", "C[4]" = "C4"))#rename diet sources
#Reorder factor levels
C_graph_data<-mutate(C_graph_data, Consumer_type  =fct_relevel(Consumer_type , "Herbivore", "Omnivore", "Carnivore"))
C_graph_data<-mutate(C_graph_data, Tissue_1_name =fct_relevel(Tissue_1_name, "Hair", "Muscle", "Liver"))
C_graph_data<-mutate(C_graph_data, Tissue_2_name =fct_relevel(Tissue_2_name, "Muscle", "Liver", "Blood", "Collagen"))
head(C_graph_data)

#Correlations
library(broom)
C_Correlation<-C_graph_data %>% group_by(Comparison) %>% do(tidy(cor.test(.$Tissue_1, .$Tissue_2)))#correlation for comparisons
C_Correlation<-as.data.frame(C_Correlation)
C_Correlation$estimate<-round(C_Correlation$estimate, 2)#round to 2 sig digits and keep trailing zero
C_Correlation<-mutate(C_Correlation, R_text = paste("italic(r) ==~", estimate))#text for graphing
C_Correlation<-mutate(C_Correlation, p.value = ifelse(p.value < .0001, sprintf('"%s"',"0.0001"),#if less than 0.00001 than 0.00001
                               ifelse(p.value < .001, sprintf('"%s"',"0.001"),#if less than 0.0001 than 0.0001    
                                      sprintf('"%.3f"',round(p.value,digits=3)))))#if greater than 0.0001 than use three digits
C_Correlation<-mutate(C_Correlation, p.value_graph = ifelse(p.value=="\\"0.0001\\"" | p.value=="\\"0.001\\"",
                                     paste("italic(P)<", p.value), paste("italic(P)==", p.value)))#P-value for graphing
C_Correlation<-select(C_Correlation, Comparison, R_text, p.value_graph)#select needed variables
print(C_Correlation)

C_graph_data_cor<-left_join(C_graph_data, C_Correlation, by = "Comparison")#join correlation data
head(C_graph_data_cor)

#Sample size
C_Sample_Size<-C_graph_data%>%
  group_by(Comparison)%>%
  summarise(n = n())%>%
  mutate(Sample_size_text = paste("italic(n) ==~", n))#text for graphing

C_graph_data_cor_n<-left_join(C_graph_data_cor, C_Sample_Size, by = "Comparison")#join sample size data
head(C_graph_data_cor_n)

label_parse <- function(breaks) {#parse text for legend
  parse(text = breaks)
}

library(lemon)
library(ggplot2)
C_Graph<- ggplot(C_graph_data_cor_n, aes(x =Tissue_1 , y = Tissue_2)) +
#one-to-one line 
  geom_abline(aes(intercept =0, slope = 1), color="gray50", size = 1, linetype = "dotted")+
#Data points
  geom_point(size=3,aes(fill = Diet_Source, color = Diet_Source, shape = Consumer_type), alpha=.65)+
  scale_shape_manual(values =c(21,22,24))+#costum shapes
  scale_color_manual(values =c("chartreuse4","darkorange3","dodgerblue3", "mediumorchid3"),label = label_parse)+#custom colors
  scale_fill_manual(values =c("chartreuse4","darkorange3","dodgerblue3", "mediumorchid3"))+#custom colors
#Regression line
  geom_smooth(method=lm,color="black", alpha=.25)+#regression line
#Facets
  facet_rep_grid(Tissue_2_name~Tissue_1_name) + coord_capped_cart(bottom='both', left='both')+#use package "lemon" to keep tick marks
  theme_bw() + theme(panel.border=element_blank(), axis.line=element_line(size=.3))+ #remove lines of graph and background color
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ #removes gridlines, but keeps border
  theme(panel.spacing = unit(-0.5, "lines"))+#add space between panels
#r-squared
  geom_text(data=C_graph_data_cor_n,aes(x=-Inf, y= Inf, hjust = -.25, vjust = 2, label=(R_text)), parse=T, #location (C3-C, Mixed-C, C3-N, Mixed-N)
            size=3.5,inherit.aes=F, color = "gray40", alpha = .65)+ #size and colors
#P-value
  geom_text(data=C_graph_data_cor_n,aes(x=-Inf, y= Inf, hjust = -.1, vjust = 4, label=(p.value_graph)), parse=T, #location (C3-C, Mixed-C, C3-N, Mixed-N)
            size=3.5,inherit.aes=F, color = "gray40", alpha = .65)+ #size and colors
#sample size
  geom_text(data=C_graph_data_cor_n,aes(x=-Inf, y= Inf, hjust = -.25, vjust = 6, label=(Sample_size_text)), parse=T, #location (C3-C, Mixed-C, C3-N, Mixed-N)
            size=3.5,inherit.aes=F, color = "gray40", alpha = .65)+ #size and colors
#Legend 
  guides(fill=FALSE,#remove fill behind color legend for Food source
         color = guide_legend(title = "Diet source"),#rename Diet_Source
         shape = guide_legend(title = "Consumer type"))+#rename Consumer_type
  theme(#see label_parsel function above along with in scale_color_manual to have subscripts
    legend.position = c(0.81, 0.265),#position of legend
    legend.title = element_text(size=11),#text size of legend title
    legend.text=element_text(size=9.5),#text size of legend labels
    legend.text.align = 0,#make text align on left for system
    legend.box.background = element_rect(colour = "gray90"),#add rectangle
    legend.background=element_blank(),#background of legend
    legend.margin = margin(.1,.1,.1,.1, unit="cm"),#spacing of box around legend
    legend.spacing.y = unit(0, "cm"))+#vertical spacing
#Axes
  scale_x_continuous(expand = c(0, 0),limits = c(-6, 8),breaks=c(-6, -4, -2, 0, 2, 4, 6))+#specify limits and breaks
  scale_y_continuous(expand = c(0, 0),limits = c(-6, 8),breaks=c(-6, -4, -2, 0, 2, 4, 6))+#specify limits and breaks
  ylab(expression("TDF " * delta^13 * "C"))+
  xlab(expression("TDF " * delta^13 * "C"))+
  theme(axis.title.y = element_text(size=11),
        axis.text.y  = element_text(angle=0, vjust=0.5, size=9.5, colour="black"))+
  theme(axis.title.x = element_text(size=11),
        axis.text.x  = element_text(angle=0, vjust=0.5, size=9.5, colour="black"))

C_Graph

# extract legend
library(cowplot)
library(rlang)
Legend <- get_legend(C_Graph)

ggdraw(C_Graph+theme(legend.position="none")) + 
  
  draw_line(x = c(0.53,.53),y = c(.72, .955),color = "white", size = 58)+#hacky way to cover unused panels
  draw_line(x = c(0.8,.8),y = c(.51, .955),color = "white", size = 52)+#hacky way to cover unused panels
  draw_plot(Legend, x = .4, y = .73, width = .5, height = .1)#placement of legend over covered area
#Save figure
ggsave("C_Tissue_TDF_Comparison.tiff",
       plot = last_plot(), width =5.5, height = 7, units = "in",
       dpi = 600) 
################################################################################################




################################################################################################
#Nitrogen tissue comparison graph
################################################################################################
library(tidyr)
N_graph_data<- N_All_Tissue_Comparisons%>%#dataset from above
  mutate(Tissue_names = Comparison)%>% #duplicate column of names 
  separate(Tissue_names, c("Tissue_1_name","Tissue_2_name"), sep = "( - )")#split column into two with names
N_graph_data<-mutate(N_graph_data, Diet_Source = fct_recode(Diet_Source, "C[3]" = "C3", "C[4]" = "C4"))#rename diet sources
#Reorder factor levels
N_graph_data<-mutate(N_graph_data, Consumer_type  =fct_relevel(Consumer_type , "Herbivore", "Omnivore", "Carnivore"))
N_graph_data<-mutate(N_graph_data, Tissue_1_name =fct_relevel(Tissue_1_name, "Hair", "Muscle", "Liver"))
N_graph_data<-mutate(N_graph_data, Tissue_2_name =fct_relevel(Tissue_2_name, "Muscle", "Liver", "Blood", "Collagen"))
head(N_graph_data)

#Correlations
library(broom)
N_Correlation<-N_graph_data %>% group_by(Comparison) %>% do(tidy(cor.test(.$Tissue_1, .$Tissue_2)))#correlation for comparisons
N_Correlation<-as.data.frame(N_Correlation)
N_Correlation$estimate<-round(N_Correlation$estimate, 2)#round to 2 sig digits and keep trailing zero
N_Correlation<-mutate(N_Correlation, R_text = paste("italic(r) ==~", estimate))#text for graphing
N_Correlation<-mutate(N_Correlation, p.value = ifelse(p.value < .0001, sprintf('"%s"',"0.0001"),#if less than 0.00001 than 0.00001
                                                      ifelse(p.value < .001, sprintf('"%s"',"0.001"),#if less than 0.0001 than 0.0001    
                                                             sprintf('"%.3f"',round(p.value,digits=3)))))#if greater than 0.0001 than use three digits
N_Correlation<-mutate(N_Correlation, p.value_graph = ifelse(p.value=="\\"0.0001\\"" | p.value=="\\"0.001\\"",
                                                            paste("italic(P)<", p.value), paste("italic(P)==", p.value)))#P-value for graphing
N_Correlation<-select(N_Correlation, Comparison, R_text, p.value_graph)#select needed variables
print(N_Correlation)

N_graph_data_cor<-left_join(N_graph_data, N_Correlation, by = "Comparison")#join correlation data
head(N_graph_data_cor)

#Sample size
N_Sample_Size<-N_graph_data%>%
  group_by(Comparison)%>%
  summarise(n = n())%>%
  mutate(Sample_size_text = paste("italic(n) ==~", n))#text for graphing

N_graph_data_cor_n<-left_join(N_graph_data_cor, N_Sample_Size, by = "Comparison")#join sample size data
head(N_graph_data_cor_n)

#make data frame for labeling facets with less than 10 samples
Tissue_1_name<-as.factor(c("Liver"))
Tissue_2_name<-as.factor(c("Collagen"))
Tissue_1<-c(3)#TDF value for the label
Tissue_2<-c(3)#TDF value for the label
No_data_label<-tibble(Tissue_1_name, Tissue_2_name, Tissue_1, Tissue_2)
No_data_label$Label<- "Only\\n2 pairwise\\nsamples"

label_parse <- function(breaks) {#parse text for legend
  parse(text = breaks)
}

library(lemon)
library(ggplot2)
N_Graph<- ggplot(N_graph_data_cor_n, aes(x =Tissue_1 , y = Tissue_2)) +
#one-to-one line 
  geom_abline(aes(intercept =0, slope = 1), color="gray50", size = 1, linetype = "dotted")+
#Data points
  geom_point(size=3,aes(fill = Diet_Source, color = Diet_Source, shape = Consumer_type), alpha=.65)+
  scale_shape_manual(values =c(21,22,24))+#costum shapes
  scale_color_manual(values =c("chartreuse4","darkorange3","dodgerblue3", "mediumorchid3"),label = label_parse)+#custom colors
  scale_fill_manual(values =c("chartreuse4","darkorange3","dodgerblue3", "mediumorchid3"))+#custom colors
#Regression line
  geom_smooth(method=lm,color="black", alpha=.25)+#regression line
#Facets
  facet_rep_grid(Tissue_2_name~Tissue_1_name) + coord_capped_cart(bottom='both', left='both')+#use package "lemon" to keep tick marks
  theme_bw() + theme(panel.border=element_blank(), axis.line=element_line(size=.3))+ #remove lines of graph and background color
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ #removes gridlines, but keeps border
  theme(panel.spacing = unit(-0.5, "lines"))+#add space between panels
#r-squared
  geom_text(data=N_graph_data_cor_n,aes(x=-Inf, y= Inf, hjust = -.25, vjust = 2, label=(R_text)), parse=T, #location (C3-C, Mixed-C, C3-N, Mixed-N)
            size=3.5,inherit.aes=F, color = "gray40", alpha = .65)+ #size and colors
#P-value
  geom_text(data=N_graph_data_cor_n,aes(x=-Inf, y= Inf, hjust = -.1, vjust = 4, label=(p.value_graph)), parse=T, #location (C3-C, Mixed-C, C3-N, Mixed-N)
            size=3.5,inherit.aes=F, color = "gray40", alpha = .65)+ #size and colors
#sample size
  geom_text(data=N_graph_data_cor_n,aes(x=-Inf, y= Inf, hjust = -.25, vjust = 6, label=(Sample_size_text)), parse=T, #location (C3-C, Mixed-C, C3-N, Mixed-N)
            size=3.5,inherit.aes=F, color = "gray40", alpha = .65)+ #size and colors
#Text for less than 10 comparisons
  geom_text(data = No_data_label, aes(x = Tissue_1, y = Tissue_2, label = Label),size = 3.5, color = "gray65")+
#Legend 
  guides(fill=FALSE,#remove fill behind color legend for Food source
         color = guide_legend(title = "Diet source"),#rename Diet_Source
         shape = guide_legend(title = "Consumer type"))+#rename Consumer_type
  theme(#see label_parsel function above along with in scale_color_manual to have subscripts
    legend.position = c(0.81, 0.265),#position of legend
    legend.title = element_text(size=11),#text size of legend title
    legend.text=element_text(size=9.5),#text size of legend labels
    legend.text.align = 0,#make text align on left for system
    legend.box.background = element_rect(colour = "gray90"),#add rectangle
    legend.background=element_blank(),#background of legend
    legend.margin = margin(.1,.1,.1,.1, unit="cm"),#spacing of box around legend
    legend.spacing.y = unit(0, "cm"))+#vertical spacing
#Axes
  scale_x_continuous(expand = c(0, 0),limits = c(0, 7.5),breaks=c(0, 2, 4, 6))+#specify limits and breaks
  scale_y_continuous(expand = c(0, 0),limits = c(0, 7.5),breaks=c(0, 2, 4, 6))+#specify limits and breaks
  ylab(expression("TDF " * delta^15 * "N"))+
  xlab(expression("TDF " * delta^15 * "N"))+
  theme(axis.title.y = element_text(size=11),
        axis.text.y  = element_text(angle=0, vjust=0.5, size=9.5, colour="black"))+
  theme(axis.title.x = element_text(size=11),
        axis.text.x  = element_text(angle=0, vjust=0.5, size=9.5, colour="black"))

N_Graph

# extract legend
library(cowplot)
library(rlang)
Legend <- get_legend(N_Graph)

ggdraw(N_Graph+theme(legend.position="none")) + 
  
  draw_line(x = c(0.52,.52),y = c(.72, .955),color = "white", size = 60)+#hacky way to cover unused panels
  draw_line(x = c(0.78,.78),y = c(.51, .955),color = "white", size = 52)+#hacky way to cover unused panels
  draw_plot(Legend, x = .4, y = .73, width = .5, height = .1)#placement of legend over covered area
#Save figure
ggsave("N_Tissue_TDF_Comparison.tiff",
       plot = last_plot(), width =5.5, height = 7, units = "in",
       dpi = 600) 
################################################################################################

