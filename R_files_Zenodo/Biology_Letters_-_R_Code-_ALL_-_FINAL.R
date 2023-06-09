
##############################################
########## Packages and working directory#####
##############################################

library(dplyr)
library(ggplot2)
library(tidyr)
library(vegan)
library(gridExtra)
library(cowplot)

##############################################
########## Data input             ############
##############################################



snail_inf_full <- read.csv("Survey_Total_Prevalence.csv", header=TRUE)
attach(snail_inf_full)

##remove dead snails. Trematode != "NA" should catch all dead/empty snails, 

snail_inf_full %>% filter(Trematode != "NA")  -> snail_inf_nodead 

##remove superfluous or qualitative data columns 
snail_inf_nodead %>% dplyr::select(-c(photo_ID, notes)) -> snail_inf 

##############################################
#####   Assessing species prevalences    #####
##############################################
### Convert long data into wide data so that 
### its a site by species matrix
### mutate - add a new column that has all 1s, 
### this allows you to give a value for each trematode species once you spread

forcats::fct_explicit_na(snail_inf$Sex, na_level = "NA") -> snail_inf$Sex

snail_inf %>% group_by(Estuary, site, Region, Species, Variant, date, number, Sex) %>% 
  dplyr::select(Estuary, site, Species, Variant, Region, date, number, Sex, Trematode) %>% 
  mutate(presence = 1) %>% spread(Trematode, presence) -> parasites_spread



### change all NAs to 0, since they 
### are snails where that species is absent. 
parasites_spread[10:33][is.na(parasites_spread[10:33])] <- 0 

### remove column '0'
parasites_spread[c(1:8,10:33)] -> parasites_spread_clean

##############################################################
########Calculate Parasite Prevalence by Site ################
##############################################################


#remove number column, date, sex column so that it does not get summarized. 
##Assess prevalence of each species for each site
parasites_spread_clean %>% ungroup() %>% 
  group_by(Estuary, Region, site, Species) %>% 
  dplyr::select(-c(number, Variant, date, Sex)) %>% 
  summarize_all(funs(mean))  -> parasites_by_site

### check that is looks fine
head(parasites_by_site)

### determine the number of individuals sampled per site
parasites_spread_clean %>% ungroup() %>% group_by(Estuary, site, Region, Species) %>% 
  dplyr::select(-c(number, date, Sex)) %>% dplyr::summarise(n=n()) -> n_by_site

### check n_by_site
n_by_site

### merge the two dataframe so you have one with the prevalence 
### of each species and the total number os snails per site
merge(n_by_site, parasites_by_site) -> parasites_summary_by_site

### check to make sure each variable in parasites_summary_by_site is categorized correctly
head(parasites_summary_by_site)

############################################
####    Prevalence/Abundance per site   ####
############################################

### calculate the total PREVALENCE of trematode infection by summing together
### all the individual relative abundances of each parasite. 
### actual ABUNDANCE of parasites is calculated by multiplying the 
### PREVALENCE by the number of snails sampled
parasites_summary_by_site %>% ungroup() %>% group_by(Estuary, Region, site, Species, n) %>% 
  mutate(prevalence_sum = (ACAN + AUST + CATA + CLOA + EUHA + HIMA + HIMB + INF + LGXI +
                             MESO + PARO + PHOC + PROB + PYGI + REBU + RECE + REMA + RENB + RENC + 
                             RENIC + REPO + SMCY + SMMI + STIC)) %>% 
  mutate(abundance = n*prevalence_sum) -> parasites_by_site_prevalence


############################################
### PRESENCE/ABSENCE & Richness per site ###
############################################

### To calculate Richness, convert the relative abundance of each species into a binary - 
### either the species is present at that site or it is not.
parasites_by_site_prevalence %>% ungroup() %>% group_by(Estuary, site, Species, Region, n, 
                                          prevalence_sum, abundance) %>% 
  mutate_all(funs(ifelse(. > 0, 1, 0))) -> presence_by_site 

### Calculate richness from presence_by_site. Exclude INF because they are unlikely to be a unique species
presence_by_site %>% mutate(richness = (ACAN + AUST + CATA + CLOA + EUHA + HIMA + HIMB + LGXI + 
                                          MESO + PARO + PHOC + PROB + PYGI + REBU + RECE + REMA + RENB + 
                                          RENC + RENIC + REPO + SMCY + SMMI + STIC)) -> richness_by_site 

### Update species names
parasites_by_site_prevalence %>% mutate(Acha = STIC, Acsp = ACAN, Ausp = AUST, Cajo = CATA, Clmi = CLOA, 
                                        Euca = EUHA, Hirh = HIMA, Hisb = HIMB, INF = INF, Lgxi = LGXI, 
                                        Meap = MESO, Pasp = PARO, Phov = PHOC, Pruc = PROB, Pysp = PYGI, 
                                        Rebu = (REBU + RENB), Rece = (RECE + RENC + RENIC), Rema = REMA, 
                                        Repo = REPO, Smcy = SMCY, Smmi = SMMI) -> parasite_by_site_newnames

write.csv(parasite_by_site_newnames, "parasite_by_site_newnames_FINAL.csv")


##############################################
########### Prevalence per Estuary ###########
##############################################

### remove number column, date, sex column so that it does not get summarized. 
### Assess prevalence of each species for each site
parasites_spread_clean  %>% ungroup() %>% group_by(Estuary, Region, Species) %>% 
  dplyr::select(-c(site, date, number, Sex, Variant)) %>%
  summarise_all(funs(mean)) -> parasites_by_Estuary

###determine the number of individuals sampled per Estuary
parasites_spread_clean %>% ungroup() %>% group_by(Estuary, Region, Species) %>% 
  dplyr::select(-c(site, date, number, Sex)) %>% 
  dplyr::summarise(n = n()) -> n_by_Estuary

### merge the two dataframe so you have one with the prevalence 
### of each species and the total number os snails per Estuary
merge(n_by_Estuary, parasites_by_Estuary) -> parasites_summary_by_Estuary


#####################################################
####         Abundance/Prevalence per estuary    ####
#####################################################

### calculate the total PREVALENCE of trematode infection by summing together
### all the individual relative abundances of each parasite. 
### actual ABUNDANCE of parasites is calculated by multiplying 
### the PREVALENCE by the number of snails sampled
parasites_summary_by_Estuary %>% 
  mutate(prevalence_sum = (ACAN + AUST + CATA + CLOA + EUHA + HIMA + HIMB + INF + LGXI +
                             MESO + PARO + PHOC + PROB + PYGI + REBU + RECE + REMA + RENB + RENC + 
                             RENIC + REPO + SMCY + SMMI + STIC)) %>% 
  mutate(abundance = n*prevalence_sum) -> parasites_by_Estuary_prevalence

#################################################
###  PRESENCE/ABSENCE & Richness per Estuary  ###
#################################################

### To calculate Richness, convert the relative abundance of each species into a binary - 
### either the species is present at that Estuary or it is not.
parasites_by_Estuary_prevalence %>% ungroup() %>% group_by(Estuary, Species, Region, n, prevalence_sum, abundance) %>% 
  mutate_all(funs(ifelse(. > 0, 1, 0))) -> presence_by_Estuary


### Calculate richness from presence_by_site. 
### Exclude INF because they are unlikely to be a unique species
presence_by_Estuary %>% mutate(richness = (ACAN + AUST + CATA + CLOA + EUHA + HIMA + HIMB + LGXI +
                                             MESO + PARO + PHOC + PROB + PYGI + REBU + RECE + REMA + RENB + RENC + 
                                             RENIC + REPO + SMCY + SMMI + STIC)) -> richness_by_Estuary 

### add richness column to parasites_by_site_prevalence dataset
parasites_by_Estuary_prevalence$richness <- richness_by_Estuary$richness

### CHange names to new names
parasites_by_Estuary_prevalence %>% mutate(Acha = STIC, Acsp = ACAN, Ausp = AUST, Cajo = CATA, Clmi = CLOA, 
       Euca = EUHA, Hirh = HIMA, Hisb = HIMB, INF = INF, Lgxi = LGXI, 
       Meap = MESO, Pasp = PARO, Phov = PHOC, Pruc = PROB, Pysp = PYGI, 
       Rebu = (REBU + RENB), Rece = (RECE + RENC + RENIC), Rema = REMA, 
       Repo = REPO, Smcy = SMCY, Smmi = SMMI) -> parasites_by_Estuary_newnames

### Export parasites_by_Estuary_complete as .csv

write.csv(parasites_by_Estuary_newnames, "parasite_by_Estuary_newnames_FINAL.csv")

#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#
#_#_# CASTE RATIO CASTE RATIO CASTE RATIO #_#_#
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#

###########################################
############# Caste Ratio Data ############
###########################################
setwd("~/Documents/Texas/Research/GRIP Trematodes/Data/Emlyn's Data/Data for Final Analyses")
caste <- read.csv("Survey_Total_caste.csv", header=TRUE)
attach(caste)

### Check that all the variables are correctly assigned and that there are no extraneous factors
head(caste)
summary(caste)

caste$date <- as.Date(caste$date, "%m/%d/%y")
caste$site <- as.factor(caste$site)
caste$Estuary <- as.factor(caste$Estuary)

### Check what variables are in snail_inf. 
### All we want is the size and sex of each snail that we dissected
head(snail_inf)
snail_inf$date <- as.Date(snail_inf$date, "%m/%d/%y")
snail_inf$site <- as.factor(snail_inf$site)

### Join caste and snail inf to get extra information about snails, such as size and sex
inner_join(caste, snail_inf, by=c("Estuary", "site", "Species", "Variant", "date", "number")) %>% 
  dplyr::select(-c(infection_no, Trematode)) -> caste_long

### check caste_long
head(caste_long)

### add a suffix to each of the estuary-specific variables
head(parasites_by_Estuary_newnames)

### add Estuary tag to each variable so we know that these are estuary level values
parasite_by_Estuary_complete_tagged <- parasites_by_Estuary_newnames
colnames(parasite_by_Estuary_complete_tagged) <- paste(colnames(parasite_by_Estuary_complete_tagged), "Est", sep = "_")

###Rename first 4 columns so that they can merge with caste_long
parasite_by_Estuary_complete_tagged <- plyr::rename(parasite_by_Estuary_complete_tagged, 
                                              c(Region_Est = "Region", Estuary_Est = "Estuary", 
                                                 Species_Est = "Species", n_Est = "n"))

### left join allows you to add Estuary-level prevalence to each dissected snail
left_join(caste_long, parasite_by_Estuary_complete_tagged, by = c("Estuary", "Region","Species")) -> 
  caste_Estuary_prevalence

### Check caste_Estuary_prevalence
head(caste_Estuary_prevalence)

caste_Estuary_prevalence$site <- as.factor(caste_Estuary_prevalence$site)

### REPEAT WITH SITE LEVEL
### add SITE tag to each variable so we know that these are SITE level values
parasite_by_site_complete_tagged <- parasite_by_site_newnames
colnames(parasite_by_site_complete_tagged) <- paste(colnames(parasite_by_site_complete_tagged), "Site", sep = "_")

### Check parasites_by_site_prevalence_tagged
parasite_by_site_complete_tagged

###Rename first 3 columns so that they can merge with caste_Estuary_prevalence
parasite_by_site_complete_tagged <- plyr::rename(parasite_by_site_complete_tagged, c(site_Site = "site", 
                                              Estuary_Site = "Estuary", Region_Site = "Region", 
                                              Species_Site = "Species", 
                                              n_Site = "n"))

caste_Estuary_prevalence$site <- as.factor(caste_Estuary_prevalence$site)
parasite_by_site_complete_tagged$site <- droplevels(parasite_by_site_complete_tagged$site)
parasite_by_site_complete_tagged$site <- as.factor(parasite_by_site_complete_tagged$site)


### left join allows you to add site-level prevalence to each dissected snail
left_join(caste_Estuary_prevalence, parasite_by_site_complete_tagged, 
          by = c("Estuary", "site", "Region", "Species")) -> 
  caste_TOTAL_minus_biomass



###Clean up data, remove early infections, remove coinfections
caste_TOTAL_minus_biomass %>% 
  filter(infection_number == 1 & species != "STIC") %>% 
  filter(Dissection.notes != "coinfection, cryptic" &
           Dissection.notes != "lots of intermediates, coinfected with REBU" &
           Dissection.notes != "MAYBE TOO EARLY" &
           Dissection.notes != "MAYBE TOO EARLY?" &
           Dissection.notes != "Might be too early? Lots of cercaria but not many reproductives" &
           Dissection.notes != "potentially coinfected with SMMI?" &
           Dissection.notes != "probably too young" &
           Dissection.notes != "probably too young" &
           Dissection.notes != "very early, maybe don't count. ") -> caste_TOTAL_clean

######### Add competitive heirarchy location
head(caste_TOTAL_clean)

### Make a dataframe that has all the species and their rank

rank <- data.frame(species=c("PARO", "HIMA", "HIMB", "CLOA", "ACAN", "EUHA"), species_rank = c(1,2,3,4,4,5))
rank$species <- as.factor(rank$species)
rank$species_rank <- as.numeric(rank$species_rank)
caste_TOTAL_clean_rank <- left_join(caste_TOTAL_clean, rank)
caste_TOTAL_clean_rank %>% dplyr::select(species,species_rank)


### Add Longitude
  Estuary <- levels(caste_TOTAL_clean_rank$Estuary)
  latitude <- c(8.893, 37.90, 34.34,32.7518,34.418, 37.6775, 
                 32.7912, 35.321, 9.362, 9.2378, 8.895, 
                 33.012, 7.635, 32.759)
  latitude.DF <- data.frame(Estuary, latitude)
  
  ##merge longitude data with caste data.table
  merge(caste_TOTAL_clean_rank, latitude.DF) -> caste_with_lat
  
  caste_with_lat$species <- as.factor(caste_with_lat$species)
  
  summary(caste_with_lat$Variant)

  ##################################################
  ########## FIND OUTLIERS FOR DATA ################
  ##################################################
  
  ###Find outliers using boxplots from total reproductives, total soldiers, and OVM
  
  ##Outliers based on Total Repro ####
  
  caste_with_lat$snail_no <- 1:nrow(caste_with_lat)

     ggplot(subset(caste_with_lat, caste_with_lat$species=="ACAN"), aes(x=species, y=total_repro)) + geom_boxplot() + labs(x="") +
    geom_text(aes(label=ifelse(total_repro <quantile(total_repro, 0.25) - 2.5 * IQR(total_repro) |
                                 total_repro > quantile(total_repro, 0.75) + 2.5 * IQR(total_repro),
                               as.character(snail_no),""), hjust=1.1))

    ggplot(subset(caste_with_lat, caste_with_lat$species=="CLOA"), aes(x=species, y=total_repro)) + geom_boxplot() + labs(x="") +
    geom_text(aes(label=ifelse(total_repro <quantile(total_repro, 0.25) - 2.5 * IQR(total_repro) |
                                 total_repro > quantile(total_repro, 0.75) + 2.5 * IQR(total_repro),
                               as.character(snail_no),""), hjust=1.1))
    
    ggplot(subset(caste_with_lat, caste_with_lat$species=="EUHA"), aes(x=species, y=total_repro)) + geom_boxplot() + labs(x="") +
      geom_text(aes(label=ifelse(total_repro <quantile(total_repro, 0.25) - 2.5 * IQR(total_repro) |
                                   total_repro > quantile(total_repro, 0.75) + 2.5 * IQR(total_repro),
                                 as.character(snail_no),""), hjust=1.1))
    
    ### HIMA outlier = 20
    ggplot(subset(caste_with_lat, caste_with_lat$species=="HIMA"), aes(x=species, y=total_repro)) + geom_boxplot() + labs(x="") +
      geom_text(aes(label=ifelse(total_repro <quantile(total_repro, 0.25) - 2.5 * IQR(total_repro) |
                                   total_repro > quantile(total_repro, 0.75) + 2.5 * IQR(total_repro),
                                 as.character(snail_no),""), hjust=1.1))
    
    ### HIMB outliers = 41, 55, 196
    FigA = ggplot(subset(caste_with_lat, caste_with_lat$species=="HIMB"), aes(x=species, y=total_repro)) + geom_boxplot() + labs(x="") +
      geom_text(aes(label=ifelse(total_repro <quantile(total_repro, 0.25) - 2.5 * IQR(total_repro) |
                                   total_repro > quantile(total_repro, 0.75) + 2.5 * IQR(total_repro),
                                 as.character(snail_no),""), hjust=1.1))
    
    ### Paro outliers = 111
    FigB = ggplot(subset(caste_with_lat, caste_with_lat$species=="PARO"), aes(x=species, y=total_repro)) + geom_boxplot() + labs(x="") +
      geom_text(aes(label=ifelse(total_repro <quantile(total_repro, 0.25) - 2.5 * IQR(total_repro) |
                                   total_repro > quantile(total_repro, 0.75) + 2.5 * IQR(total_repro),
                                 as.character(snail_no),""), hjust=1.1))
    
### Outliers based on total_sold ####
    ggplot(subset(caste_with_lat, caste_with_lat$species=="ACAN"), aes(x=species, y=total_sold)) + geom_boxplot() + labs(x="") +
      geom_text(aes(label=ifelse(total_sold <quantile(total_sold, 0.25) - 2.5 * IQR(total_sold) |
                                   total_sold > quantile(total_sold, 0.75) + 2.5 * IQR(total_sold),
                                 as.character(snail_no),""), hjust=1.1))
    
    ### ClOA outliers = 164
    FigC = ggplot(subset(caste_with_lat, caste_with_lat$species=="CLOA"), aes(x=species, y=total_sold)) + geom_boxplot() + labs(x="") +
      geom_text(aes(label=ifelse(total_sold <quantile(total_sold, 0.25) - 2.5 * IQR(total_sold) |
                                   total_sold > quantile(total_sold, 0.75) + 2.5 * IQR(total_sold),
                                 as.character(snail_no),""), hjust=1.1))
    

    ggplot(subset(caste_with_lat, caste_with_lat$species=="EUHA"), aes(x=species, y=total_sold)) + geom_boxplot() + labs(x="") +
      geom_text(aes(label=ifelse(total_sold <quantile(total_sold, 0.25) - 2.5 * IQR(total_sold) |
                                   total_sold > quantile(total_sold, 0.75) + 2.5 * IQR(total_sold),
                                 as.character(snail_no),""), hjust=1.1))
    

    ggplot(subset(caste_with_lat, caste_with_lat$species=="HIMA"), aes(x=species, y=total_sold)) + geom_boxplot() + labs(x="") +
      geom_text(aes(label=ifelse(total_sold <quantile(total_sold, 0.25) - 2.5 * IQR(total_sold) |
                                   total_sold > quantile(total_sold, 0.75) + 2.5 * IQR(total_sold),
                                 as.character(snail_no),""), hjust=1.1))
    
    ### HIMB outliers = 55
    FigD = ggplot(subset(caste_with_lat, caste_with_lat$species=="HIMB"), aes(x=species, y=total_sold)) + geom_boxplot() + labs(x="") +
      geom_text(aes(label=ifelse(total_sold <quantile(total_sold, 0.25) - 2.5 * IQR(total_sold) |
                                   total_sold > quantile(total_sold, 0.75) + 2.5 * IQR(total_sold),
                                 as.character(snail_no),""), hjust=1.1))
    
    ### Paro outliers = 166
    FigE = ggplot(subset(caste_with_lat, caste_with_lat$species=="PARO"), aes(x=species, y=total_sold)) + geom_boxplot() + labs(x="") +
      geom_text(aes(label=ifelse(total_sold <quantile(total_sold, 0.25) - 2.5 * IQR(total_sold) |
                                   total_sold > quantile(total_sold, 0.75) + 2.5 * IQR(total_sold),
                                 as.character(snail_no),""), hjust=1.1))
    
    
### Outliers based on OVM
    
    caste_with_lat %>% mutate(OVM = ((mantle_sold_total+mid_sold_total)/(total_sold))) -> caste_with_OVM
  
    ggplot(subset(caste_with_OVM, caste_with_OVM$species=="ACAN"), aes(x=species, y=OVM)) + geom_boxplot() + labs(x="") +
      geom_text(aes(label=ifelse(OVM <quantile(OVM, 0.25) - 2.5 * IQR(OVM) |
                                   OVM > quantile(OVM, 0.75) + 2.5 * IQR(OVM),
                                 as.character(snail_no),""), hjust=1.1))
    
    ### ClOA outliers = 162
    FigF = ggplot(subset(caste_with_OVM, caste_with_OVM$species=="CLOA"), aes(x=species, y=OVM)) + geom_boxplot() + labs(x="") +
      geom_text(aes(label=ifelse(OVM <quantile(OVM, 0.25) - 2.5 * IQR(OVM) |
                                   OVM > quantile(OVM, 0.75) + 2.5 * IQR(OVM),
                                 as.character(snail_no),""), hjust=1.1))
    
      ggplot(subset(caste_with_OVM, caste_with_OVM$species=="EUHA"), aes(x=species, y=OVM)) + geom_boxplot() + labs(x="") +
      geom_text(aes(label=ifelse(OVM <quantile(OVM, 0.25) - 2.5 * IQR(OVM) |
                                   OVM > quantile(OVM, 0.75) + 2.5 * IQR(OVM),
                                 as.character(snail_no),""), hjust=1.1))
    
    ggplot(subset(caste_with_OVM, caste_with_OVM$species=="HIMA"), aes(x=species, y=OVM)) + geom_boxplot() + labs(x="") +
      geom_text(aes(label=ifelse(OVM <quantile(OVM, 0.25) - 2.5 * IQR(OVM) |
                                   OVM > quantile(OVM, 0.75) + 2.5 * IQR(OVM),
                                 as.character(snail_no),""), hjust=1.1))
    

    ggplot(subset(caste_with_OVM, caste_with_OVM$species=="HIMB"), aes(x=species, y=OVM)) + geom_boxplot() + labs(x="") +
      geom_text(aes(label=ifelse(OVM <quantile(OVM, 0.25) - 2.5 * IQR(OVM) |
                                   OVM > quantile(OVM, 0.75) + 2.5 * IQR(OVM),
                                 as.character(snail_no),""), hjust=1.1))
    

    ggplot(subset(caste_with_OVM, caste_with_OVM$species=="PARO"), aes(x=species, y=OVM)) + geom_boxplot() + labs(x="") +
      geom_text(aes(label=ifelse(OVM <quantile(OVM, 0.25) - 2.5 * IQR(OVM) |
                                   OVM > quantile(OVM, 0.75) + 2.5 * IQR(OVM),
                                 as.character(snail_no),""), hjust=1.1))
  ###plot figures 
  
  setwd("~/Documents/Texas/Research/GRIP Trematodes/Data/Emlyn's Data/Data for Final Analyses/For Paper 1 - All analyses/Full Data for Ryan and Mark")
  pdf("Survey_outliers.pdf", height = 8, width = 8)
  
  ggdraw() +
    draw_plot(FigA, x = 0.0, y = 0.5, width = 0.33, height = .5) +
    draw_plot(FigB, x = 0.33, y = 0.5, width = 0.33, height = .5) +
    draw_plot(FigC, x = 0.66, y = 0.5, width = 0.33, height = .5) +
    draw_plot(FigD, x = 0, y = 0, width = 0.33, height = .5) +
    draw_plot(FigE, x = 0.33, y = 0, width = 0.33, height = .5) +
    draw_plot(FigF, x = 0.66, y = 0, width = 0.33, height = .5) +
   
    draw_plot_label(c("A", "B","C", "D", "E", "F"), 
                    c(0, .33, 0.66, 0, 0.33, 0.66), 
                    c(1, 1, 1, .5, .5, .5), size = 14)
  
  dev.off()
  
  
  
  summary(caste_with_lat$Species) 
  
  ### remove 7 snails that are outliers (2.5 IQR above or below 50%....)
  ###What are the outliers?
  
  caste_with_lat %>% filter(snail_no == "41" | snail_no == "55" | snail_no == "196" |
                              snail_no == "166" | snail_no == "111" | snail_no == "164" |
                              snail_no == "162") 
  
  caste_with_lat %>% filter(snail_no != "41" & snail_no != "55" & snail_no != "196" &
                     snail_no != "166" & snail_no != "111" & snail_no != "164" & 
                       snail_no != "162") -> caste_no_outliers
  
  caste_no_outliers %>% dplyr::select(-c(snail_no)) -> caste_final
  
  summary(caste_final$species)

  setwd("~/Documents/Texas/Research/GRIP Trematodes/Data/Emlyn's Data/Data for Final Analyses")
write.csv(caste_final, "caste_TOTAL_FINAL.csv")



###########################################################
#_#_#_#_#_#   GLMM Models for Survey Data     #_#_#_#_#_#_#
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#

###---------Packages and working directory--------### ####

setwd("~/Documents/Texas/Research/GRIP Trematodes")
library(dplyr)
library(ggplot2)
library(tidyr)
library(wiqid) #for model selection
library(lmerTest)
library(scales)
library(lme4)
library(car)
library(MASS)
library(vcd)
library(DHARMa)
library(visreg)
library(cowplot)

###----------Data input and manipulations---------### ####
### try original data from MS submission 
caste <- read.csv("caste_TOTAL_FINAL.csv", header=TRUE)
attach(caste)
head(caste)



### Select only the variables of interest. Filter only snail species C. californica ####
caste %>% dplyr::select(Estuary, site, Species, Region, number, species, ratio, 
                        visc_sold_total, mid_sold_total, mantle_sold_total, total_sold, 
                        total_repro, size, prevalence_sum_Est, prevalence_sum_Site,
                        latitude, species_rank) %>% 
  filter(Species == "C. californica") -> caste_for_models


summary(caste$species)

### manipulate dependent variables of interest  ####
caste_for_models %>% mutate(sold_log = log(total_sold), 
                            prop_mantle = (mantle_sold_total/total_sold),
                            prop_mantle_log = log(mantle_sold_total/total_sold), 
                            prop_OVM = (mantle_sold_total+mid_sold_total)/total_sold) %>%
  mutate(prop_OVM_log = log(prop_OVM)) %>%
  remove_missing() %>% droplevels()-> caste_for_models_2

summary(caste_for_models$species)


### Rescale so that the variables are not on very different scales ####
###Standardize total_repro, prevalence, latitude (p 122)
caste_for_models_2$total_repro_s <- (caste_for_models_2$total_repro - mean(caste_for_models_2$total_repro))/
  sd(caste_for_models_2$total_repro)
caste_for_models_2$total_sold_s <- (caste_for_models_2$total_sold - mean(caste_for_models_2$total_sold))/
  sd(caste_for_models_2$total_sold)
caste_for_models_2$prevalence_sum_Est_s <- (caste_for_models_2$prevalence_sum_Est - mean(caste_for_models_2$prevalence_sum_Est))/
  sd(caste_for_models_2$prevalence_sum_Est)
caste_for_models_2$latitude_s <- (caste_for_models_2$latitude - mean(caste_for_models_2$latitude))/
  sd(caste_for_models_2$latitude)
caste_for_models_2$prevalence_sum_Site_s <- (caste_for_models_2$prevalence_sum_Site - mean(caste_for_models_2$prevalence_sum_Site))/
  sd(caste_for_models_2$prevalence_sum_Site)
caste_for_models_2$size_s <- (caste_for_models_2$size - mean(caste_for_models_2$size))/
  sd(caste_for_models_2$size)

### Format variables
caste_for_models_2$species <- factor(caste_for_models_2$species, levels = c("PARO", "HIMA", "HIMB", "CLOA", "ACAN", "EUHA"))
caste_for_models_2$mantle_sold_total <- as.integer(caste_for_models_2$mantle_sold_total)
caste_for_models_2$visc_sold_total <- as.integer(caste_for_models_2$visc_sold_total)
caste_for_models_2$ID <- seq.int(nrow(caste_for_models_2))

### Summary statistics - How many colonies, estuaries, and sites for each species?

caste_for_models_2 %>% filter(species == "ACAN") %>% 
  dplyr::select(Estuary, site) %>% droplevels() -> ACAN_sum

caste_for_models_2 %>% filter(species == "CLOA") %>% 
  dplyr::select(Estuary, site) %>% droplevels() -> CLOA_sum

caste_for_models_2 %>% filter(species == "EUHA") %>% 
  dplyr::select(Estuary, site) %>% droplevels() -> EUHA_sum

caste_for_models_2 %>% filter(species == "HIMA") %>% 
  dplyr::select(Estuary, site) %>% droplevels() -> HIMA_sum

caste_for_models_2 %>% filter(species == "HIMB") %>% 
  dplyr::select(Estuary, site) %>% droplevels() -> HIMB_sum

caste_for_models_2 %>% filter(species == "PARO") %>% 
  dplyr::select(Estuary, site) %>% droplevels() ->PARO_sum






###------------------Total Soldiers---------------### ####
###Trying things out to figure out distribution, overdispersion, etc. 
    ###total soldiers fits a poisson distribution
    distplot(caste_for_models_2$total_sold, type="poisson")
    histogram(caste_for_models_2$total_sold)
    ### First try Poisson GLMM with 1|Estuary/Site with Estuary level prevalence. #### 
    ### including (1|Est/Site) was singular because Estuary-level variance is ~0. So using (1/site) instead
    fit_poisson_sold_Est_old <- glmer(as.integer(total_sold) ~ total_repro_s + size_s +  
                                        species + prevalence_sum_Est_s + (1|site), family=poisson(), data = caste_for_models_2)
    
    # summary(fit_poisson_sold_Est)
    # Anova(fit_poisson_sold_Est)
    
    ###Is the high significance due to overdispersion? ####
    overdisp_fun <- function(model) {
      rdf <- df.residual(model)
      rp <- residuals(model,type="pearson")
      Pearson.chisq <- sum(rp^2)
      prat <- Pearson.chisq/rdf
      pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
      c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
    }
    
    overdisp_fun(fit_poisson_sold_Est_old) ###YES
    
    ###Add a random effect of colony (1|ID) to account for overdispersion
    
    ### CHECK VIF Values for snail size and total repro to make sure they aren't collinear. ####
    ###Since I can't do that for a glmer, just do it for glm
    fit_sold_VIF <- vif(glm(as.integer(total_sold) ~ total_repro_s + species + size_s +  
                              prevalence_sum_Est + prevalence_sum_Site  + latitude_s, 
                            family=poisson(), data = caste_for_models_2))
    
    fit_sold_VIF
    

###FULL FUNCTIONING MODEL ####
###This model has all three of our "spatial" variables - Estuary-level prevalence, site-level prevalence, and latitude. 
fit_poisson1_sold_FULL <- glmer(as.integer(total_sold) ~ total_repro_s:species + total_repro_s + size_s +  
                                  species + prevalence_sum_Est + prevalence_sum_Site  + latitude_s + 
                                  (1|Estuary/site) + (1|ID), glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)),
                                family=poisson(), data = caste_for_models_2)

overdisp_fun(fit_poisson1_sold_FULL)
Anova(fit_poisson1_sold_FULL, type = "III") 
summary(fit_poisson1_sold_FULL)


### CHECK ASSUMPTIONS ####
### How well does model fit? ##NORMAL RESIDUAL PLOTS DONT APPLY FOR POISSON. use DHARMa package
resDHARMa = simulateResiduals(fit_poisson1_sold_FULL)
plot(resDHARMa)

par(mfrow=c(1,3))
plot(fitted(fit_poisson1_sold_FULL),residuals(fit_poisson1_sold_FULL))
hist(residuals(fit_poisson1_sold_FULL))
qqnorm(residuals(final_glmer_OVM_FULL))
qqline(residuals(final_glmer_OVM_FULL))

### Simple model ####
### The simple model included NONE of the spatial variables    
fit_poisson1_sold_simple <- glmer(as.integer(total_sold) ~ total_repro_s:species + total_repro_s + size_s +  
                                    species + (1|Estuary/site) + (1|ID), 
                                  glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)),
                                  family=poisson(), data = caste_for_models_2)

overdisp_fun(fit_poisson1_sold_simple)
Anova(fit_poisson1_sold_simple, type = "III")
summary(fit_poisson1_sold_simple)  

### CHECK ASSUMPTIONS ####
### How well does model fit? ##NORMAL RESIDUAL PLOTS DONT APPLY FOR POISSON. use DHARMa package
resDHARMa = simulateResiduals(fit_poisson1_sold_simple)
plot(resDHARMa)  

### Estuary-level prevalence ####
### Only include the spatial variables Estuary-level prevalence
fit_poisson1_sold_Est <- glmer(as.integer(total_sold) ~ total_repro_s:species + total_repro_s + size_s +  
                                 species + prevalence_sum_Est + (1|Estuary/site) + (1|ID), 
                               glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)),
                               family=poisson(), data = caste_for_models_2)

overdisp_fun(fit_poisson1_sold_Est)
Anova(fit_poisson1_sold_Est, type="III")
summary(fit_poisson1_sold_Est)

### CHECK ASSUMPTIONS ####
### How well does model fit? ##NORMAL RESIDUAL PLOTS DONT APPLY FOR POISSON. use DHARMa package
resDHARMa = simulateResiduals(fit_poisson1_sold_Est)
plot(resDHARMa)

### Site-level prevalence ####
fit_poisson1_sold_Site <- glmer(as.integer(total_sold) ~ total_repro_s:species + total_repro_s + size_s +  
                                  species + prevalence_sum_Site + (1|Estuary/site) + (1|ID), 
                                glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)),
                                family=poisson(), data = caste_for_models_2)

overdisp_fun(fit_poisson1_sold_Site )
Anova(fit_poisson1_sold_Site, type = "III")
summary(fit_poisson1_sold_Site )


### CHECK ASSUMPTIONS ####
### How well does model fit? ##NORMAL RESIDUAL PLOTS DONT APPLY FOR POISSON. use DHARMa package
resDHARMa = simulateResiduals(fit_poisson1_sold_Site)
plot(resDHARMa)


### Latitude ####
fit_poisson1_sold_latitude <- glmer(as.integer(total_sold) ~ total_repro_s:species + total_repro_s + size_s +  
                                      species + latitude_s + (1|Estuary/site) + (1|ID), 
                                    glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)),
                                    family=poisson(), data = caste_for_models_2)

overdisp_fun(fit_poisson1_sold_latitude)
Anova(fit_poisson1_sold_latitude, type = "III")
summary(fit_poisson1_sold_latitude)

### CHECK ASSUMPTIONS ####
### How well does model fit? ##NORMAL RESIDUAL PLOTS DONT APPLY FOR POISSON. use DHARMa package
resDHARMa = simulateResiduals(fit_poisson1_sold_latitude)
plot(resDHARMa)

### Model Comparison 1 - Compare estuary-level prevalence, site-level prevalence, and latitude and simple model ####
### Site-level prevalence is included in best fit model - use site-level prevalence going forward
Total_Sold_AICtable1 <- AICtable(AICc(fit_poisson1_sold_FULL, fit_poisson1_sold_Site, 
                                      fit_poisson1_sold_Est, fit_poisson1_sold_latitude, fit_poisson1_sold_simple))
Total_Sold_AICtable1

Anova(fit_poisson1_sold_Site, type = "III")


fit_poisson1_sold_norepro <- glmer(as.integer(total_sold) ~ total_repro_s:species + size_s +  
                                         species + prevalence_sum_Site + (1|Estuary/site) + (1|ID), 
                                       glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)),
                                       family=poisson(), data = caste_for_models_2)

Anova(fit_poisson1_sold_norepro, type = "III")
summary(fit_poisson1_sold_norepro)
### ADD INTERACTION: Site-level prevalence:species ####
### Because of issues of overfitting, only at this point do we investigate 
### an interaction between prevalence and species identity
fit_poisson1_sold_Sitespecies <- glmer(as.integer(total_sold) ~ total_repro_s:species + total_repro_s + size_s +  
                                         species + prevalence_sum_Site:species + prevalence_sum_Site +  (1|Estuary/site) + (1|ID), 
                                       glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)),
                                       family=poisson(), data = caste_for_models_2)

overdisp_fun(fit_poisson1_sold_Sitespecies )
Anova(fit_poisson1_sold_Sitespecies, type = "III")
summary(fit_poisson1_sold_Sitespecies)


### Backward Selection ####
###Remove snail size ####
fit_poisson1_sold_nosize <- glmer(as.integer(total_sold) ~ total_repro_s:species + total_repro_s + prevalence_sum_Site +
                                    species + (1|Estuary/site) + (1|ID), 
                                  glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)),
                                  family=poisson(), data = caste_for_models_2)

overdisp_fun(fit_poisson1_sold_nosize)
Anova(fit_poisson1_sold_nosize, type = "III")
summary(fit_poisson1_sold_nosize)

### Remove site-level prevalence ####
fit_poisson1_sold_nosize_noprev <- glmer(as.integer(total_sold) ~ total_repro_s:species + total_repro_s +
                                           species + (1|Estuary/site) + (1|ID), 
                                         glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)),
                                         family=poisson(), data = caste_for_models_2)

overdisp_fun(fit_poisson1_sold_nosize_noprev)
Anova(fit_poisson1_sold_nosize_noprev, type = "III")
summary(fit_poisson1_sold_nosize_noprev)

### Remove total reproduction:species interaction
fit_poisson1_sold_nosize_noprev_noint <- glmer(as.integer(total_sold) ~  species + total_repro_s + (1|Estuary/site) + (1|ID), 
                                               glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)),
                                               family=poisson(), data = caste_for_models_2)

overdisp_fun(fit_poisson1_sold_nosize_noprev_noint)
Anova(fit_poisson1_sold_nosize_noprev_noint, type = "II")
summary(fit_poisson1_sold_nosize_noprev_noint)

### Remove total reproduction
fit_poisson1_sold_nosize_noprev_noint_norepro <- glmer(as.integer(total_sold) ~  species + (1|Estuary/site) + (1|ID), 
                                                       glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)),
                                                       family=poisson(), data = caste_for_models_2)

overdisp_fun(fit_poisson1_sold_nosize_noprev_noint_norepro)
Anova(fit_poisson1_sold_nosize_noprev_noint_norepro, type = "II")
summary(fit_poisson1_sold_nosize_noprev_noint_norepro)


### Model Comparison 2 - Compare all models ####
Total_Sold_AICtable <- AICtable(AICc(fit_poisson1_sold_Sitespecies, fit_poisson1_sold_FULL, fit_poisson1_sold_Site, 
                                     fit_poisson1_sold_Est, fit_poisson1_sold_latitude, fit_poisson1_sold_simple,  
                                     fit_poisson1_sold_nosize, fit_poisson1_sold_nosize_noprev, fit_poisson1_sold_nosize_noprev_noint, 
                                     fit_poisson1_sold_nosize_noprev_noint_norepro))

#setwd("~/Documents/Texas/Research/GRIP Trematodes/Data/Emlyn's Data/Data for Final Analyses/For Ch 2 - All analyses/Full Data for Ryan and Mark/Chapter 2 -  Biology Letters")
write.csv(Total_Sold_AICtable, "AIC_sold_12_12_19.csv")

Anova_sold <-  Anova(fit_poisson1_sold_Site, type = "III") 
summary(fit_poisson1_sold_Site)
write.csv(Anova_sold, "Anova_Sold_12_12_19.csv")



###-------------- NEW TEST - binomial GLM  --------### ####
###----------------Distribution: Mantle------------### ####

### Model selection  ####
### Is soldier deployment normally distributed? No - needs to be log transformed ####

hist(caste_for_models_2$prop_mantle, breaks=20)

### CHECK VIF Values for snail size and total repro to make sure they aren't collinear. ####
###Since I can't do that for a glmer, just do it for glm

final_glmer_mantle_FULL <- glmer(prop_mantle ~ size_s + prevalence_sum_Est + prevalence_sum_Site + latitude_s + 
                                   total_repro_s + species + (1|Estuary/site), weights = total_sold,
                                 data = caste_for_models_2, family = "binomial")

cor(data.frame(caste_for_models_2$total_repro_s,  caste_for_models_2$size_s))

Anova(final_glmer_mantle_FULL, type = "II")
summary(final_glmer_mantle_FULL)

### Model Assumptions - FIT! ####
par(mfrow=c(1,3))
plot(fitted(final_glmer_mantle_FULL),residuals(final_glmer_mantle_FULL))
hist(residuals(final_glmer_mantle_FULL, type = "pearson"))
qqplot(final_glmer_mantle_FULL)
qqnorm(residuals(final_glmer_mantle_FULL, type = "pearson"))
qqline(residuals(final_glmer_mantle_FULL, type = "pearson"))

### Model selection  ####

### CHECK VIF Values for snail size and total repro to make sure they aren't collinear. ####
###Since I can't do that for a glmer, just do it for glm
fit_mantle_VIF <- vif(glmer(prop_mantle ~ size_s + prevalence_sum_Est + prevalence_sum_Site + latitude_s + 
                              total_repro_s + species + (1|Estuary/site), weights = total_sold,
                            data = caste_for_models_2, family = "binomial"))

cor(data.frame(caste_for_models_2$total_repro_s,  caste_for_models_2$size_s))

### FULL MODEL ####
### This includes all three "spatial" parameters - site-level prevalence, estuary-level prevalence, and latitude
final_glmer_mantle_FULL <- glmer(prop_mantle ~ size_s + prevalence_sum_Est + 
                                   prevalence_sum_Site + latitude_s + total_repro_s + 
                                   species + (1|Estuary/site), weights = total_sold,
                                 family = "binomial", data = caste_for_models_2)

Anova(final_glmer_mantle_FULL, type = "II")

### Model Assumptions - FIT! ####
par(mfrow=c(1,3))
plot(fitted(final_glmer_mantle_FULL),residuals(final_glmer_mantle_FULL))
hist(residuals(final_glmer_mantle_FULL))
qqplot(final_glmer_mantle_FULL)
qqnorm(residuals(final_glmer_mantle_FULL))
qqline(residuals(final_glmer_mantle_FULL))

### Estuary-level prevalence ####
final_glmer_mantle_Est <- glmer(prop_mantle ~ size_s + prevalence_sum_Est + 
                                  total_repro_s + species + (1|Estuary/site), weights = total_sold,
                                family = "binomial", data = caste_for_models_2)

Anova(final_glmer_mantle_Est, type = "II")

### Model Assumptions - FIT! ####
par(mfrow=c(1,3))
plot(fitted(final_glmer_mantle_Est),residuals(final_glmer_mantle_Est))
hist(residuals(final_glmer_mantle_Est))
qqnorm(residuals(final_glmer_mantle_Est))
qqline(residuals(final_glmer_mantle_Est))

### Site-Level Prevalence ####
final_glmer_mantle_Site <-  glmer(prop_mantle ~ size_s + prevalence_sum_Site + 
                                    total_repro_s + species + (1|Estuary/site), weights = total_sold,
                                  family = "binomial", data = caste_for_models_2)
Anova(final_glmer_mantle_Site, type = "II")

### Model Assumptions - FIT! ####
par(mfrow=c(1,3))
plot(fitted(final_glmer_mantle_Site),residuals(final_glmer_mantle_Site))
hist(residuals(final_glmer_mantle_Site))
qqnorm(residuals(final_glmer_mantle_Site))
qqline(residuals(final_glmer_mantle_Site))

### latitude Prevalence ####

final_glmer_mantle_latitude <-  glmer(prop_mantle ~ size_s + latitude_s + 
                                        total_repro_s + species + (1|Estuary/site), weights = total_sold,
                                      family = "binomial", data = caste_for_models_2)

Anova(final_glmer_mantle_latitude, type = "II")

### Model Assumptions - FIT! ####
par(mfrow=c(1,3))
plot(fitted(final_glmer_mantle_latitude),residuals(final_glmer_mantle_latitude))
hist(residuals(final_glmer_mantle_latitude))
qqnorm(residuals(final_glmer_mantle_latitude))
qqline(residuals(final_glmer_mantle_latitude))

### Simple ####

?glmer

final_glmer_mantle_simple <-  glmer(prop_mantle ~ size_s + 
                                      total_repro_s + species + (1|Estuary/site), weights = total_sold,
                                    family = "binomial", data = caste_for_models_2)
Anova(final_glmer_mantle_simple, type = "II")

### Model Assumptions - FIT! ####
par(mfrow=c(1,3))
plot(fitted(final_glmer_mantle_simple),residuals(final_glmer_mantle_simple))
hist(residuals(final_glmer_mantle_simple))
qqnorm(residuals(final_glmer_mantle_simple))
qqline(residuals(final_glmer_mantle_simple))

### Model Comparison 1 - Prev_EST vs prev_Site vs Latitude and best model ####
### Simple model is the best fit!
Total_mantle_AICtable1 <- AICtable(AICc(final_glmer_mantle_simple, final_glmer_mantle_latitude, 
                                        final_glmer_mantle_Site, final_glmer_mantle_Est,
                                        final_glmer_mantle_FULL))
Total_mantle_AICtable1 

### Backward Selection starting with 
### Site-level prevalence model (which is worse than simple, but only by <2) ####
### Remove total Repro ####
final_glmer_mantle_norepro <- glmer(prop_mantle ~ size_s + prevalence_sum_Site + species + (1|Estuary/site), weights = total_sold,
                                    family = "binomial", data = caste_for_models_2)
Anova(final_glmer_mantle_norepro, type = "II")

### Model Assumptions - FIT! ####
par(mfrow=c(1,3))
plot(fitted(final_glmer_mantle_norepro),residuals(final_glmer_mantle_norepro))
hist(residuals(final_glmer_mantle_norepro))
qqnorm(residuals(final_glmer_mantle_norepro))
qqline(residuals(final_glmer_mantle_norepro))

### Remove site-level prevalence ####
final_glmer_mantle_norepro_noprev <- glmer(prop_mantle ~ size_s + species + (1|Estuary/site), weights = total_sold,
                                           family = "binomial", data = caste_for_models_2)
Anova(final_glmer_mantle_norepro_noprev, type = "II")

### Model Assumptions - FIT! ####
par(mfrow=c(1,3))
plot(fitted(final_glmer_mantle_norepro_noprev),residuals(final_glmer_mantle_norepro_noprev))
hist(residuals(final_glmer_mantle_norepro_noprev))
qqnorm(residuals(final_glmer_mantle_norepro_noprev))
qqline(residuals(final_glmer_mantle_norepro_noprev))


### Model Comparison 1 - Prev_EST vs prev_Site vs Latitude and best model ####
Total_mantle_AICtable <- AICtable(AICc(final_glmer_mantle_simple, final_glmer_mantle_latitude, 
                                       final_glmer_mantle_Site, final_glmer_mantle_Est,
                                       final_glmer_mantle_FULL, final_glmer_mantle_norepro, 
                                       final_glmer_mantle_norepro_noprev))
write.csv(Total_mantle_AICtable, "AIC_mantle_12_12_19.csv")

summary(final_glmer_mantle_norepro)
Anova_mantle <-  Anova(final_glmer_mantle_norepro, type = "II") 
write.csv(Anova_mantle, "Anova_mantle_12_12_19.csv")


###-------------- Distribution: OVM------------### ####
### What type of distribution - Binomial!####

hist(caste_for_models_2$prop_OVM, breaks = 20)

### GLMER! Make sure that the model will run. Yup!
final_glmer_OVM_FULL <- glmer(prop_OVM ~ size_s + prevalence_sum_Est + 
                                prevalence_sum_Site + latitude_s + 
                                species + total_repro_s + (1|Estuary/site), 
                              weights = total_sold,
                              family = "binomial", data = caste_for_models_2)
summary(final_glmer_OVM_Est)
Anova(final_glmer_OVM_Est, type = "II")


### CHECK VIF Values for snail size and total repro to make sure they aren't collinear. ####
###Since I can't do that for a glmer, just do it for glm
fit_OVM_VIF <- vif(glmer(prop_OVM ~ size_s + prevalence_sum_Est + 
                           prevalence_sum_Site + latitude_s + 
                           species + total_repro_s + (1|Estuary/site), 
                         weights = total_sold,
                         family = "binomial", data = caste_for_models_2))

### FULL MODEL ####
final_glmer_OVM_Full <- glmer(prop_OVM ~ size_s + prevalence_sum_Est + 
                                prevalence_sum_Site + latitude_s + 
                                species + total_repro_s + (1|Estuary/site), 
                              weights = total_sold,
                              family = "binomial", data = caste_for_models_2)
Anova(final_glmer_OVM_Full)

### Model Assumptions - FIT! ####
par(mfrow=c(1,3))
plot(fitted(final_glmer_OVM_Full),residuals(final_glmer_OVM_Full,  type = "pearson"))
hist(residuals(final_glmer_OVM_Full,  type = "pearson"))  
qqnorm(residuals(final_glmer_OVM_Full,  type = "pearson"))
qqline(residuals(final_glmer_OVM_Full,  type = "pearson"))

### Estuary level prevalence ####
final_glmer_OVM_Estuary <- glmer(prop_OVM ~ size_s + prevalence_sum_Est + 
                                   species + total_repro_s + (1|Estuary/site), 
                                 weights = total_sold,
                                 family = "binomial", data = caste_for_models_2)
Anova(final_glmer_OVM_Estuary)

### Model Assumptions - FIT! ####
par(mfrow=c(1,3))
plot(fitted(final_glmer_OVM_Estuary),residuals(final_glmer_OVM_Estuary))
hist(residuals(final_glmer_OVM_Estuary))
qqnorm(residuals(final_glmer_OVM_Estuary))
qqline(residuals(final_glmer_OVM_Estuary))

### Site-Level Prevalence ####
final_glmer_OVM_Site <- glmer(prop_OVM ~ size_s + prevalence_sum_Site + 
                                species + total_repro_s + (1|Estuary/site), 
                              weights = total_sold,
                              family = "binomial", data = caste_for_models_2)

Anova(final_glmer_OVM_Site, type = "II")

### Model Assumptions - FIT! ####
par(mfrow=c(1,3))
plot(fitted(final_glmer_OVM_Site),residuals(final_glmer_OVM_Site))
hist(residuals(final_glmer_OVM_Site))
qqnorm(residuals(final_glmer_OVM_Site))
qqline(residuals(final_glmer_OVM_Site))

### latitude  ####
final_glmer_OVM_latitude <- glmer(prop_OVM ~ size_s + latitude_s + 
                                    species + total_repro_s + (1|Estuary/site), 
                                  weights = total_sold,
                                  family = "binomial", data = caste_for_models_2)
Anova(final_glmer_OVM_latitude, type = "II")

### Model Assumptions - FIT! ####
par(mfrow=c(1,3))
plot(fitted(final_glmer_OVM_latitude),residuals(final_glmer_OVM_latitude))
hist(residuals(final_glmer_OVM_latitude))
qqnorm(residuals(final_glmer_OVM_latitude))
qqline(residuals(final_glmer_OVM_latitude))

### Simple ####
final_glmer_OVM_simple <- glmer(prop_OVM ~ size_s + 
                                  species + total_repro_s + (1|Estuary/site), 
                                weights = total_sold,
                                family = "binomial", data = caste_for_models_2)

Anova(final_glmer_OVM_simple, type = "II")

### Model Assumptions - FIT! ####
par(mfrow=c(1,3))
plot(fitted(final_glmer_OVM_simple),residuals(final_glmer_OVM_simple))
hist(residuals(final_glmer_OVM_simple))
qqnorm(residuals(final_glmer_OVM_simple))
qqline(residuals(final_glmer_OVM_simple))

### Model comparison between simple, full and prev/lat -> USE LATITUDE
Total_OVM_AICtable1 <- AICtable(AICc(final_glmer_OVM_simple, final_glmer_OVM_latitude, 
                                     final_glmer_OVM_Site, final_glmer_OVM_Estuary,
                                     final_glmer_OVM_Full))
Total_OVM_AICtable1
### Backward Selection starting with latitude
### Remove repro####
final_glmer_OVM_latitude_norepro <- glmer(prop_OVM ~ size_s + latitude_s + 
                                            species + (1|Estuary/site), 
                                          weights = total_sold,
                                          family = "binomial", data = caste_for_models_2)
Anova(final_glmer_OVM_latitude_norepro, type = "II")

### Model Assumptions - FIT! ####
par(mfrow=c(1,3))
plot(fitted(final_glmer_OVM_latitude_norepro),residuals(final_glmer_OVM_latitude_norepro))
hist(residuals(final_glmer_OVM_latitude_norepro))
qqnorm(residuals(final_glmer_OVM_latitude_norepro))
qqline(residuals(final_glmer_OVM_latitude_norepro))

### Remove latitude ####
final_glmer_OVM_latitude_norepro_nolat <- glmer(prop_OVM ~ size_s + 
                                                  species + (1|Estuary/site), 
                                                weights = total_sold,
                                                family = "binomial", data = caste_for_models_2)
Anova(final_glmer_OVM_latitude_norepro_nolat, type = "II")

### Model Assumptions - FIT! ####
par(mfrow=c(1,3))
plot(fitted(final_glmer_OVM_latitude_norepro_nolat),residuals(final_glmer_OVM_latitude_norepro_nolat))
hist(residuals(final_glmer_OVM_latitude_norepro_nolat, ))
qqnorm(residuals(final_glmer_OVM_latitude_norepro_nolat))
qqline(residuals(final_glmer_OVM_latitude_norepro_nolat))

?residuals()
### Model Comparison 1 - Prev_EST vs prev_Site vs Latitude and best model ####
Total_OVM_AICtable <- AICtable(AICc(final_glmer_OVM_latitude_norepro_nolat, 
                                    final_glmer_OVM_latitude_norepro,
                                    final_glmer_OVM_simple, final_glmer_OVM_latitude, 
                                    final_glmer_OVM_Site, final_glmer_OVM_Estuary,
                                    final_glmer_OVM_Full))
write.csv(Total_OVM_AICtable, "AIC_OVM_12_12_19.csv")

Anova_OVM <-  Anova(final_glmer_OVM_latitude_norepro_nolat, type = "II") 
summary(final_glmer_OVM_latitude_norepro_nolat)
write.csv(Anova_OVM, "Anova_OVM_12_12_19.csv")




###----------------- RESIDUAL PLOTS---------------### ####
library(visreg)
library(effects)

### PLOT PREVALENCE  ####
blank_data <- data.frame(species = c("PARO", "PARO", "HIMA", "HIMA", "HIMB", "HIMB",
                                     "CLOA", "CLOA", "ACAN", "ACAN", "EUHA", "EUHA"),
                         x = 0, y = c(7.4, 8.4, 4.2, 5.2, 5.7, 6.7, 6.7, 7.7, 6.2, 7.2, 5,6))

plot_prev_species <- visreg(fit_poisson1_sold_Site, "prevalence_sum_Site", by="species", 
                            gg=TRUE, band=TRUE, 
                            line=list(col="red"),
                            points=list(size=4, pch=21, col="black", bg="black", alpha=0.4)) + 
  theme_classic(base_size = 14) + 
  theme(strip.background = element_rect(colour="white", fill="white"),
        strip.text.x = element_text(size=12, face="bold"), 
        axis.text = element_text(colour = "black"),
        legend.position = "none") + 
  facet_wrap(~species, ncol=2, scale="free") + 
  geom_point(alpha = 0.4, cex=4) + 
  geom_blank(data = blank_data, aes(x = x, y = y)) + 
  labs(x="Site-Level Prevalence ", y= "Log(Total Soldier) Residuals") 

plot_prev_delta <- visreg(fit_poisson1_sold_Site, "prevalence_sum_Site", type="contrast",
                          gg=TRUE, band=TRUE, line=list(col="red"),
                          fill=list(fill="grey80"),
                          points=list(size=4, pch=21, col="black", bg="black", alpha=0.4)) + 
  theme_classic(base_size = 14) + geom_hline(yintercept=0, lty=2) +
  theme(axis.text = element_text(colour = "black")) + 
  labs(x="Site-Level Prevalence", y= expression(Delta*"Log(Total Soldiers) Residuals"))

### Combine plots ####
setwd("~/Documents/Texas/Research/GRIP Trematodes/Data/Emlyn's Data/Data for Final Analyses/For Ch 2 - All analyses/Full Data for Ryan and Mark")
pdf("Partial Residuals - Prevalence.pdf", height = 6, width = 10)

ggdraw() +
  draw_plot(plot_prev_species, x = 0.0, y = 0, width = .45, height = 1) +
  draw_plot(plot_prev_delta, x = 0.45, y = 0, width = 0.55, height = 1) +
  draw_plot_label(c("A", "B"), c(0, .45), c(1, 1), size = 14)
dev.off()

#### 08/13/19 - Trying to change axes of plot_prev_species ####

### Because our total repro variable is standardized, 
### I need to determine a good "standardized" value 
### figure out the equation to change standardized to unstandardized
caste_for_models_2 %>% dplyr::select(total_repro_s, total_repro) %>% 
  arrange(total_repro) -> compare_repros

lm(compare_repros$total_repro_s~compare_repros$total_repro) 

### Check that this equation gets the transformed data - CLOSE ENOUGH
compare_repros %>% mutate(check = (-0.7641664 + 0.0002466*total_repro))

### What is the total_repro_s value equivalent of 1000 reproductives?

repro_100 <- (-0.7641664 + 0.0002466*100)
repro_1000 <- (-0.7641664 + 0.0002466*1000)
repro_10000 <- (-0.7641664 + 0.0002466*10000)


plot_prev_species2 <- 
  visreg(fit_poisson1_sold_Site, "prevalence_sum_Site", by="species", 
         gg=TRUE, band=TRUE, 
         line=list(col="red"), cond=list(total_repro_s = repro_1000),
         points=list(size=4, pch=21, col="black", bg="black", alpha=0.4)) + 
  theme_classic(base_size = 14) + 
  theme(strip.background = element_rect(colour="white", fill="white"),
        strip.text.x = element_text(size=12, face="bold"), 
        axis.text = element_text(colour = "black"),
        legend.position = "none") + 
  facet_wrap(~species, ncol=2, scale="free") + 
  geom_point(alpha = 0.4, cex=4) + 
  geom_blank(data = blank_data, aes(x = x, y = y)) + 
  labs(x="Site-Level Prevalence ", y= "Log(Total Soldier) Residuals") + 
  scale_y_continuous(breaks = c(log(40),log(60),log(80),log(100), log(120), log(140), log(160), log(200), log(250), log(300), log(350), log(400), log(600), 
                                log(800), log(1000), log(1400),log(1800), log(2200), log(2600), log(3000), 
                                log(3600), log(4200), log(4800)), labels=c("40","60","80","100","120", "140","160", "200", "250","300",
                                                                           "350","400", "600", "800", "1000", "1400", "1800",
                                                                           "2200", "2600", "3000", "3600", "4200", "4800"))


### Combine plots ####
setwd("~/Documents/Texas/Research/GRIP Trematodes/Data/Emlyn's Data/Data for Final Analyses/For Ch 2 - All analyses/Full Data for Ryan and Mark")
pdf("Partial Residuals - Prevalence2.pdf", height = 6, width = 10)

ggdraw() +
  draw_plot(plot_prev_species2, x = 0.0, y = 0, width = .45, height = 1) +
  draw_plot(plot_prev_delta, x = 0.45, y = 0, width = 0.55, height = 1) +
  draw_plot_label(c("A", "B"), c(0, .45), c(1, 1), size = 14)
dev.off()


### Separate between species - best fit model
visreg_species <- visreg(fit_poisson1_sold_Site, "prevalence_sum_Site", by="species", plot=FALSE)


total_repro_visreg_PARO <- subset(visreg_species, species %in% c("PARO")) 
total_repro_visreg_HIMA <- subset(visreg_species, species %in% c("HIMA"))
total_repro_visreg_HIMB <- subset(visreg_species, species %in% c("HIMB"))
total_repro_visreg_CLOA <- subset(visreg_species, species %in% c("CLOA"))
total_repro_visreg_ACAN <- subset(visreg_species, species %in% c("ACAN"))
total_repro_visreg_EUHA <- subset(visreg_species, species %in% c("EUHA"))

figPARO <-  plot(total_repro_visreg_PARO, gg=TRUE, overlay=FALSE, line=list(col="red"),
                 cond=list(total_repro_s = repro_1000)) + 
  theme_classic(base_size = 14) +
  theme(strip.background = element_rect(colour="white", fill="white"),
        strip.text.x = element_text(size=12, face="bold"), 
        axis.text = element_text(colour = "black"),
        legend.position = "none", 
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  geom_point(col="black",alpha = 0.4, cex=4)  +
  labs(y= "") + 
  scale_y_continuous(limits = c(7.2,8.25), breaks = c(log(1000), log(1400),log(1800), log(2200), log(2600), log(3000), 
                                                      log(3400), log(3800), log(4200)), labels=c("1000", "1400", "1800",
                                                                                                 "2200", "2600", "3000", "3400", "3800", "4200"))

figHIMA <- plot(total_repro_visreg_HIMA, gg=TRUE, overlay=FALSE, line=list(col="red"),
                cond=list(total_repro_s = repro_1000)) + 
  theme_classic(base_size = 14) +
  theme(strip.background = element_rect(colour="white", fill="white"),
        strip.text.x = element_text(size=12, face="bold"), 
        axis.text = element_text(colour = "black"),
        legend.position = "none", 
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  geom_point(col="black",alpha = 0.4, cex=4)  +
  labs(y= "") + 
  scale_y_continuous(limits = c(4.2,5.25), 
                     breaks = c(log(20), log(40), log(60), log(80), log(100), 
                                log(120), log(140), log(160), log(180)),
                     labels=c("  20", " 40", " 60"," 80", " 100"," 120", " 140", " 160", " 180"))

figHIMB <- plot(total_repro_visreg_HIMB, gg=TRUE, overlay=FALSE, line=list(col="red"),
                cond=list(total_repro_s = repro_1000)) + 
  theme_classic(base_size = 14) +
  theme(strip.background = element_rect(colour="white", fill="white"),
        strip.text.x = element_text(size=12, face="bold"), 
        axis.text = element_text(colour = "black"),
        legend.position = "none", 
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  geom_point(col="black",alpha = 0.4, cex=4)  +
  labs( y= "Total Soldier # (residuals)") + 
  scale_y_continuous(limits = c(5.6,6.65), 
                     breaks = c(log(300), log(400), log(500), log(600), log(700), log(800)),
                     labels = c(" 300", " 400"," 500", " 600", " 700", " 800"))

figCLOA <- plot(total_repro_visreg_CLOA, gg=TRUE, overlay=FALSE, line=list(col="red"),
                cond=list(total_repro_s = repro_1000)) + 
  theme_classic(base_size = 14) +
  theme(strip.background = element_rect(colour="white", fill="white"),
        strip.text.x = element_text(size=12, face="bold"), 
        axis.text = element_text(colour = "black"),
        legend.position = "none", 
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  geom_point(col="black",alpha = 0.4, cex=4)  +
  labs( y= "") + 
  scale_y_continuous(limits = c(6.5,7.55), 
                     breaks = c(log(800), log(1000), log(1200), log(1400), log(1600), log(1800)),
                     labels = c("800", "1000","1200", "1400", "1600", "1800"))


figACAN <- plot(total_repro_visreg_ACAN, gg=TRUE, overlay=FALSE, line=list(col="red"),
                cond=list(total_repro_s = repro_1000)) + 
  theme_classic(base_size = 14) +
  theme(strip.background = element_rect(colour="white", fill="white"),
        strip.text.x = element_text(size=12, face="bold"), 
        axis.text = element_text(colour = "black"),
        legend.position = "none", 
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  geom_point(col="black",alpha = 0.4, cex=4)  +
  labs(y= "") + 
  scale_y_continuous(limits = c(6.2,7.2), 
                     breaks = c(log(400), log(600), log(800), log(1000), log(1200), log(1400)),
                     labels = c("400", "600", "800", "1000","1200", "1400"))


figEUHA <- plot(total_repro_visreg_EUHA, gg=TRUE, overlay=FALSE, line=list(col="red"),
                cond=list(total_repro_s = repro_1000)) + 
  theme_classic(base_size = 14) +
  theme(strip.background = element_rect(colour="white", fill="white"),
        strip.text.x = element_text(size=12, face="bold"), 
        axis.text = element_text(colour = "black"),
        legend.position = "none", 
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  geom_point(col="black",alpha = 0.4, cex=4) +
  labs( y= "") + 
  scale_y_continuous(limits = c(5.0,6.05),breaks = c(log(150), log(200), log(250), log(300), log(350), log(400)), 
                     labels=c(" 150", " 200", " 250"," 300", " 350"," 400"))

### Combine plots ####
setwd("~/Documents/Texas/Research/GRIP Trematodes/Data/Emlyn's Data/Data for Final Analyses/For Ch 2 - All analyses/Full Data for Ryan and Mark")
pdf("Partial Residuals - Prevalence no species interaction.pdf", height = 6, width = 10)

ggdraw() +
  draw_plot(figPARO, x = 0.04, y = .66, width = .205, height = .35) +
  draw_plot(figHIMA, x = 0.245, y = .66, width = .205, height = .35) +
  draw_plot(figHIMB, x = 0.04, y = .35, width = .205, height = .35) +
  draw_plot(figCLOA, x = 0.245, y = .35, width = .205, height = .35) +
  draw_plot(figACAN, x = 0.04, y = 0.04, width = .205, height = .35) +
  draw_plot(figEUHA, x = 0.245, y = 0.04, width = .205, height = .35) +
  draw_plot(plot_prev_delta, x = 0.45, y = 0, width = 0.55, height = 1) +
  draw_plot_label(c("A", "B"), c(0, .45), c(1, 1), size = 14)
dev.off()


### PLOT model with species interaction ####
### Separate between species
### Separate plot into each species so I can manipualte axes  
visreg_species_int <- visreg(fit_poisson1_sold_Sitespecies, "prevalence_sum_Site", by="species", plot=FALSE)


total_repro_visreg_PARO_int <- subset(visreg_species_int, species %in% c("PARO")) 
total_repro_visreg_HIMA_int <- subset(visreg_species_int, species %in% c("HIMA"))
total_repro_visreg_HIMB_int <- subset(visreg_species_int, species %in% c("HIMB"))
total_repro_visreg_CLOA_int <- subset(visreg_species_int, species %in% c("CLOA"))
total_repro_visreg_ACAN_int <- subset(visreg_species_int, species %in% c("ACAN"))
total_repro_visreg_EUHA_int <- subset(visreg_species_int, species %in% c("EUHA"))

### PARO ####
figPARO_int <-  plot(total_repro_visreg_PARO_int, gg=TRUE, overlay=FALSE, line=list(col="blue"),
                     cond=list(total_repro_s = repro_1000)) + 
  theme_classic(base_size = 14) +
  theme(strip.background = element_rect(colour="white", fill="white"),
        strip.text.x = element_text(size=12, face="bold"), 
        axis.text = element_text(colour = "black"),
        legend.position = "none", 
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  geom_point(col="black",alpha = 0.4, cex=4)  +
  labs(y= "") +
  scale_y_continuous(limits = c(7.2,8.25), 
                     breaks = c(log(1000), log(1400),log(1800), log(2200), log(2600), log(3000), 
                                log(3400), log(3800), log(4200)), labels=c("1000", "1400", "1800",
                                                                           "2200", "2600", "3000", 
                                                                           "3400", "3800", "4200")) +
  geom_line(data = total_repro_visreg_PARO$fit, aes(x=prevalence_sum_Site, y = visregFit), col="red", cex=1)


### HIMA ####
figHIMA_int <- plot(total_repro_visreg_HIMA_int, gg=TRUE, overlay=FALSE, line=list(col="blue"),
                    cond=list(total_repro_s = repro_1000)) + 
  theme_classic(base_size = 14) +
  theme(strip.background = element_rect(colour="white", fill="white"),
        strip.text.x = element_text(size=12, face="bold"), 
        axis.text = element_text(colour = "black"),
        legend.position = "none", 
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  geom_point(col="black",alpha = 0.4, cex=4)  +
  labs(y= "") + 
  scale_y_continuous(limits = c(4.2,5.25), 
                     breaks = c(log(20), log(40), log(60), log(80), log(100), 
                                log(120), log(140), log(160), log(180)),
                     labels=c("  20", " 40", " 60"," 80", " 100"," 120", " 140", " 160", " 180")) +
  geom_line(data = total_repro_visreg_HIMA$fit, aes(x=prevalence_sum_Site, y = visregFit), col="red", cex=1)



figHIMB_int <- plot(total_repro_visreg_HIMB_int, gg=TRUE, overlay=FALSE, line=list(col="blue"),
                    cond=list(total_repro_s = repro_1000)) + 
  theme_classic(base_size = 14) +
  theme(strip.background = element_rect(colour="white", fill="white"),
        strip.text.x = element_text(size=12, face="bold"), 
        axis.text = element_text(colour = "black"),
        legend.position = "none", 
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  geom_point(col="black",alpha = 0.4, cex=4)  +
  labs( y= "Total Soldier # (residuals)") + 
  scale_y_continuous(limits = c(5.6,6.65), 
                     breaks = c(log(300), log(400), log(500), log(600), log(700), log(800)),
                     labels = c(" 300", " 400"," 500", " 600", " 700", " 800")) +
  geom_line(data = total_repro_visreg_HIMB$fit, aes(x=prevalence_sum_Site, y = visregFit), col="red", cex=1)



figCLOA_int <- plot(total_repro_visreg_CLOA_int, gg=TRUE, overlay=FALSE, line=list(col="blue"),
                    cond=list(total_repro_s = repro_1000)) + 
  theme_classic(base_size = 14) +
  theme(strip.background = element_rect(colour="white", fill="white"),
        strip.text.x = element_text(size=12, face="bold"), 
        axis.text = element_text(colour = "black"),
        legend.position = "none", 
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  geom_point(col="black",alpha = 0.4, cex=4)  +
  labs( y= "") + 
  scale_y_continuous(limits = c(6.5,7.55), 
                     breaks = c(log(800), log(1000), log(1200), log(1400), log(1600), log(1800)),
                     labels = c("800", "1000","1200", "1400", "1600", "1800")) +
  geom_line(data = total_repro_visreg_CLOA$fit, aes(x=prevalence_sum_Site, y = visregFit), col="red", cex=1)




figACAN_int <- plot(total_repro_visreg_ACAN_int, gg=TRUE, overlay=FALSE, line=list(col="blue"),
                    cond=list(total_repro_s = repro_1000)) + 
  theme_classic(base_size = 14) +
  theme(strip.background = element_rect(colour="white", fill="white"),
        strip.text.x = element_text(size=12, face="bold"), 
        axis.text = element_text(colour = "black"),
        legend.position = "none", 
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  geom_point(col="black",alpha = 0.4, cex=4)  +
  labs(y= "") + 
  scale_y_continuous(limits = c(6.2,7.2), 
                     breaks = c(log(400), log(600), log(800), log(1000), log(1200), log(1400)),
                     labels = c("400", "600", "800", "1000","1200", "1400")) +
  geom_line(data = total_repro_visreg_ACAN$fit, aes(x=prevalence_sum_Site, y = visregFit), col="red", cex=1)




figEUHA_int <- plot(total_repro_visreg_EUHA_int, gg=TRUE, overlay=FALSE, line=list(col="blue"),
                    cond=list(total_repro_s = repro_1000)) + 
  theme_classic(base_size = 14) +
  theme(strip.background = element_rect(colour="white", fill="white"),
        strip.text.x = element_text(size=12, face="bold"), 
        axis.text = element_text(colour = "black"),
        legend.position = "none", 
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  geom_point(col="black",alpha = 0.4, cex=4) +
  labs( y= "") + 
  scale_y_continuous(limits = c(5.0,6.05),breaks = c(log(150), log(200), log(250), log(300), log(350), log(400)), 
                     labels=c(" 150", " 200", " 250"," 300", " 350"," 400")) +
  geom_line(data = total_repro_visreg_EUHA$fit, aes(x=prevalence_sum_Site, y = visregFit), col="red", cex=1)


### Combine plots ####
setwd("~/Documents/Texas/Research/GRIP Trematodes/Data/Emlyn's Data/Data for Final Analyses/For Ch 2 - All analyses/Full Data for Ryan and Mark")
pdf("Partial Residuals - Prevalence - Interaction and non-interaction models.pdf", height = 6, width = 10)

ggdraw() +
  draw_plot(figPARO_int, x = 0.04, y = .66, width = .205, height = .35) +
  draw_plot(figHIMA_int, x = 0.245, y = .66, width = .205, height = .35) +
  draw_plot(figHIMB_int, x = 0.04, y = .35, width = .205, height = .35) +
  draw_plot(figCLOA_int, x = 0.245, y = .35, width = .205, height = .35) +
  draw_plot(figACAN_int, x = 0.04, y = 0.04, width = .205, height = .35) +
  draw_plot(figEUHA_int, x = 0.245, y = 0.04, width = .205, height = .35) +
  draw_plot(plot_prev_delta, x = 0.45, y = 0, width = 0.55, height = 1) +
  draw_plot_label(c("A", "B"), c(0, .45), c(1, 1), size = 14)
dev.off()


###Plot for CEID presentation
setwd("~/Documents/Texas/Research/GRIP Trematodes/Data/Emlyn's Data/Data for Final Analyses/For Ch 2 - All analyses/Full Data for Ryan and Mark")
pdf("CEID TALK Partial Residuals - Prevalence - Interaction and non-interaction models.pdf", height = 6, width = 8)

ggdraw() +
  draw_plot(figPARO_int, x = 0, y = .5, width = .33, height = .5) +
  draw_plot(figHIMA_int, x = 0.33, y = .5, width = .33, height = .5) +
  draw_plot(figHIMB_int, x = 0.67, y = .5, width = .33, height = .5) +
  draw_plot(figCLOA_int, x = 0.0, y = .0, width = .33, height = .5) +
  draw_plot(figACAN_int, x = 0.33, y = 0.0, width = .33, height = .5) +
  draw_plot(figEUHA_int, x = 0.67, y = 0.0, width = .33, height = .5)

dev.off()

