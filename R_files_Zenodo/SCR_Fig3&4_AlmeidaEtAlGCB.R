#'@title Spatially-clustered rarefaction analysis for Almeida et al. GCB
#'@author Ryan J. Almeida, Kevin G. Smith
#'@description Upload data, run rarefaction analyses, calculate species richness, effective number of species,
#'and Pe values as described in Almeida et al. GCB. Note that the analyses in the manuscript use the data
#'from the "Post2" survey as post-disturbance data. Post1 and Post2 are qualitatively similar and do not change the
#'overall interpretations of our findings.

##First step is to download file GoldenrodData_Master.xlsx and put in wd (check using #getwd(), set using #setwd())


#Use readxl package to import data from "GoldenrodData_Master.xlsx"
library(readxl)
#Use Cairo to export high quality figures
library(Cairo)
#Use tidyverse for data manipulation and plotting
library(tidyverse)
#Use vegan for easy calculations of richness from species x site matrices
library(vegan)

#################### DATA CLEANING FUNCTIONS ####################
#'
#'@param comm, a species x site matrix
#'@param spp_list, a vector of species names. Each observation should correspond with a row in
#'comm in order (i.e. spp_list[n] is the species in comm[n,])
#'@return a species x site matrix that only includes species eligible for the analyses outlined
#'in Almeida et al. GCB
#'
#'@description for each site supplied to dataClean, the species that were not eligible for the 
#'analysis will be removed. For complete details, see Methods section of Almeida et al. GCB
dataclean <- function(comm, spp_list){
  
  rownames(comm) <- spp_list[,1]
  comm_original <- comm
  
  #Combine Acanaloniids
  comm["Acanalonia_conica",] <- comm_original["Acanalonia_conica",] +
    comm_original["Acanaloniidae_nymph",]
  comm <- comm[!rownames(comm) %in% "Acanaloniidae_nymph",]
  
  #Anisopodidae -> Asilidae
  pos <- which(rownames(comm)== "Anisopodidae")
  rownames(comm)[pos] <- "Asilidae"
  
  #Chrysodeixis includens -> Ctenoplusia oxygramma
  pos <- which(rownames(comm)== "Chrysodeixis_includens")
  rownames(comm)[pos] <- "Ctenoplusia_oxygramma"
  
  #Euschistus servus -> Pentatomidae
  pos <- which(rownames(comm)== "Euschistus_servus")
  rownames(comm)[pos] <- "Pentatomidae"
  
  #Combine crab spiders
  comm["Thomisidae",] <- comm_original["Misumenia_vatia",] + 
    comm_original["Green_crab_sider_brown_back/eyes",] +
    comm_original["Thomisidae",]
  comm <- comm[!rownames(comm) %in% "Misumenia_vatia",]
  comm <- comm[!rownames(comm) %in% "Green_crab_sider_brown_back/eyes",]
  
  #Walking cone -> Coleophora
  pos <- which(rownames(comm)== "Walking_cone")
  rownames(comm)[pos] <- "Coleophora"
  
  #Upside down moth -> Arta statlis
  pos <- which(rownames(comm)== "Upside_down_moth")
  rownames(comm)[pos] <- "Arta_statlis"
  
  #Brown spider white striped legs -> Yunohamella lyrica
  pos <- which(rownames(comm)== "Brown_Spider_White_Striped_Legs")
  rownames(comm)[pos] <- "Yunohamella_lyrica"
  
  #Draeculacephala -> Scaphytopis spp
  pos <- which(rownames(comm)== "Draeculacephala_spp.")
  rownames(comm)[pos] <- "Scaphytopius_spp"
  
  #Brown horned moth -> Acrolophus popaenella
  pos <- which(rownames(comm)== "Brown_horned_moth")
  rownames(comm)[pos] <- "Acrolophus_popaenella"
  
  #Light brown patterned moth -> Argyrotaenia velutinana
  comm["Light/dark_brown_argyle_moth",] <- comm_original["Light_brown_moth_dark_stripes",] +
    comm_original["Light/dark_brown_argyle_moth",]
  comm <- comm[!rownames(comm) %in% "Light_brown_moth_dark_stripes",]
  pos <- which(rownames(comm)== "Light/dark_brown_argyle_moth")
  rownames(comm)[pos] <- "Argyrotaenia_velutinana"
  
  #Combine Harmonia axyridis
  comm["Harmonia_axyridis",] <- comm_original["Harmonia_axyridis",] +
    comm_original["H._axyridis_larvae",]
  comm <- comm[!rownames(comm) %in% "H._axyridis_larvae",]
  
  #Combine Halyomorphys halys
  comm["Halyomorpha_halys",] <- comm_original["Halyomorpha_halys",] +
    comm_original["Halyomorpha_halys_nymph",]
  comm <- comm[!rownames(comm) %in% "Halyomorpha_halys_nymph",]
  
  #Combine Graphocephala
  comm["Graphocephala_spp.",] <- comm_original["Graphocephala_spp.",] +
    comm_original["Hopper_nymph",]
  comm <- comm[!rownames(comm) %in% "Hopper_nymph",]
  
  #Combine Chrysopa
  comm["Chrysopa_spp.",] <- comm_original["Chrysopa_spp.",] +
    comm_original["Lacewing_larvae",]
  comm <- comm[!rownames(comm) %in% "Lacewing_larvae",]
  
  #Combine Melanoplus
  comm["Melanoplus_differentialis",] <- comm_original["Melanoplus_differentialis",] +
    comm_original["Melanoplus_femurrubrum",] +
    comm_original["Melanoplus_nymph",]
  comm <- comm[!rownames(comm) %in% "Melanoplus_femurrubrum",]
  comm <- comm[!rownames(comm) %in% "Melanoplus_nymph",]
  pos <- which(rownames(comm)== "Melanoplus_differentialis")
  rownames(comm)[pos] <- "Melanoplus_spp"
  
  #Combine Zelus
  comm["Zelus_sp.",] <- comm_original["Zelus_sp.",] +
    comm_original["Zelus_nymph",]
  comm <- comm[!rownames(comm) %in% "Zelus_nymph",]
  
  #Remove species due to phenology
  comm <- comm[!rownames(comm) %in% "Brachypnoea_sp.",]
  comm <- comm[!rownames(comm) %in% "Entylia_sp.",]
  comm <- comm[!rownames(comm) %in% "Eutreta_sp.",]
  comm <- comm[!rownames(comm) %in% "Green_Araneus_sp",]
  comm <- comm[!rownames(comm) %in% "Limonius_sp",]
  comm <- comm[!rownames(comm) %in% "Neoconocephalus_spp.",]
  comm <- comm[!rownames(comm) %in% "Neoscona_spp",]
  comm <- comm[!rownames(comm) %in% "Pale_Jumper_Black_Eyes",]
  comm <- comm[!rownames(comm) %in% "Xenox_tigrinus",]
  
  #Remove transient species
  comm <- comm[!rownames(comm) %in% "Augochlora_spp.",]
  comm <- comm[!rownames(comm) %in% "Bombus_sp",]
  comm <- comm[!rownames(comm) %in% "Calliphoridae",]
  comm <- comm[!rownames(comm) %in% "Chrysotoxum_sp.",]
  comm <- comm[!rownames(comm) %in% "Drosophilidae",]
  comm <- comm[!rownames(comm) %in% "Halictidae",]
  comm <- comm[!rownames(comm) %in% "Hesperiidae",]
  comm <- comm[!rownames(comm) %in% "Sarcophagidae",]
  comm <- comm[!rownames(comm) %in% "Sciaridae",]
  comm <- comm[!rownames(comm) %in% "Tachinidae",]
  comm <- comm[!rownames(comm) %in% "Toxomerus_spp.",]
  comm <- comm[!rownames(comm) %in% "Tritoxa_incurva",]
  comm <- comm[!rownames(comm) %in% "Trupanea_sp.",]
  
  #Remove milkweed species
  comm <- comm[!rownames(comm) %in% "Ocyptamus_spp.",]
  comm <- comm[!rownames(comm) %in% "Syrphid_larvae",]
  comm <- comm[!rownames(comm) %in% "Labidomera_clivicollis",]
  comm <- comm[!rownames(comm) %in% "Labidomera_clivicollis_larvae",]
  comm <- comm[!rownames(comm) %in% "Aphis_nerii",]
  
  #Sort rows alphabetically
  comm <- comm[order(rownames(comm)),]
  
  return(comm)
}

#'
#'@param comm, a species x site matrix
#'@param spp_list, a vector of species names. Each observation should correspond with a row in
#'comm in order (i.e. spp_list[n] is the species in comm[n,])
#'@return a species x site matrix that only includes species eligible for the analyses outlined
#'in Almeida et al. GCB for the Post3 survey
#'
#'@description for each site supplied to datacleanPost3, the species that were not eligible for the 
#'analysis will be removed. For complete details, see Methods section of Almeida et al. GCB
datacleanPost3 <- function(comm, spp_list){
  rownames(comm) <- spp_list[,1]
  
  comm <- comm[!rownames(comm) %in% "Calliphoridae",]
  comm <- comm[!rownames(comm) %in% "Drosophila",]
  comm <- comm[!rownames(comm) %in% "Halictidae",]
  comm <- comm[!rownames(comm) %in% "Tachinidae",]
  comm <- comm[!rownames(comm) %in% "Toxomerus_spp.",]
  comm <- comm[!rownames(comm) %in% "Trupanea_sp.",]
  comm <- comm[!rownames(comm) %in% "Ocypteamus_spp.",]
  comm <- comm[!rownames(comm) %in% "Labidomera_clivicollis",]
  comm <- comm[!rownames(comm) %in% "TFF_spp._(Drosophila?)",]
  comm <- comm[!rownames(comm) %in% "metallic_gold_sweat_bee_fly",]
  comm <- comm[!rownames(comm) %in% "brown_fly_short_yelllow_antennae",]
  
  
  #Sort rows alphabetically
  comm <- comm[order(rownames(comm)),]
}

##Upload the list of species
spp_list <- read_xlsx("GoldenrodData_Master.xlsx", sheet = "spp_list", col_names = F)
##Upload list of species for Post3
spp_listPost3 <- read_xlsx("GoldenrodData_Master.xlsx", sheet = "spp_list2018", col_names = F)

#################### DATA UPLOAD AND CLEANING ####################
#Each sampling period (Pre, Post1, Post2) consists of 3 days of sampling (i.e. Pre1-3, Post1.1-1.3, Post2.1-2.3, Post3.1-3.3)
### Pre 1
Bag1_Pre1 <- read_xlsx("GoldenrodData_Master.xlsx",
                       sheet = "Pre1_Bag1")
Bag1_Pre1_Clean <- dataclean(as.data.frame(Bag1_Pre1), as.data.frame(spp_list))

Bag2_Pre1 <- read_xlsx("GoldenrodData_Master.xlsx",
                       sheet = "Pre1_Bag2")
Bag2_Pre1_Clean <- dataclean(as.data.frame(Bag2_Pre1), as.data.frame(spp_list))

Bag3_Pre1 <- read_xlsx("GoldenrodData_Master.xlsx",
                       sheet = "Pre1_Bag3")
Bag3_Pre1_Clean <- dataclean(as.data.frame(Bag3_Pre1), as.data.frame(spp_list))

Bag4_Pre1 <- read_xlsx("GoldenrodData_Master.xlsx",
                       sheet = "Pre1_Bag4")
Bag4_Pre1_Clean <- dataclean(as.data.frame(Bag4_Pre1), as.data.frame(spp_list))

Bag6_Pre1 <- read_xlsx("GoldenrodData_Master.xlsx",
                       sheet = "Pre1_Bag6")
Bag6_Pre1_Clean <- dataclean(as.data.frame(Bag6_Pre1), as.data.frame(spp_list))

Bag7_Pre1 <- read_xlsx("GoldenrodData_Master.xlsx",
                       sheet = "Pre1_Bag7")
Bag7_Pre1_Clean <- dataclean(as.data.frame(Bag7_Pre1), as.data.frame(spp_list))

Bag8_Pre1 <- read_xlsx("GoldenrodData_Master.xlsx",
                       sheet = "Pre1_Bag8")
Bag8_Pre1_Clean <- dataclean(as.data.frame(Bag8_Pre1), as.data.frame(spp_list))

Bag9_Pre1 <- read_xlsx("GoldenrodData_Master.xlsx",
                       sheet = "Pre1_Bag9")
Bag9_Pre1_Clean <- dataclean(as.data.frame(Bag9_Pre1), as.data.frame(spp_list))

Bag10_Pre1 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Pre1_Bag10")
Bag10_Pre1_Clean <- dataclean(as.data.frame(Bag10_Pre1), as.data.frame(spp_list))

Bag11_Pre1 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Pre1_Bag11")
Bag11_Pre1_Clean <- dataclean(as.data.frame(Bag11_Pre1), as.data.frame(spp_list))

Bag12_Pre1 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Pre1_Bag12")
Bag12_Pre1_Clean <- dataclean(as.data.frame(Bag12_Pre1), as.data.frame(spp_list))

Bag13_Pre1 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Pre1_Bag13")
Bag13_Pre1_Clean <- dataclean(as.data.frame(Bag13_Pre1), as.data.frame(spp_list))

Bag14_Pre1 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Pre1_Bag14")
Bag14_Pre1_Clean <- dataclean(as.data.frame(Bag14_Pre1), as.data.frame(spp_list))

Bag15_Pre1 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Pre1_Bag15")
Bag15_Pre1_Clean <- dataclean(as.data.frame(Bag15_Pre1), as.data.frame(spp_list))

Bag16_Pre1 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Pre1_Bag16")
Bag16_Pre1_Clean <- dataclean(as.data.frame(Bag16_Pre1), as.data.frame(spp_list))

Bag17_Pre1 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Pre1_Bag17")
Bag17_Pre1_Clean <- dataclean(as.data.frame(Bag17_Pre1), as.data.frame(spp_list))

Bag18_Pre1 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Pre1_Bag18")
Bag18_Pre1_Clean <- dataclean(as.data.frame(Bag18_Pre1), as.data.frame(spp_list))

Bag19_Pre1 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Pre1_Bag19")
Bag19_Pre1_Clean <- dataclean(as.data.frame(Bag19_Pre1), as.data.frame(spp_list))

Bag20_Pre1 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Pre1_Bag20")
Bag20_Pre1_Clean <- dataclean(as.data.frame(Bag20_Pre1), as.data.frame(spp_list))

### Pre 2 
Bag1_Pre2 <- read_xlsx("GoldenrodData_Master.xlsx",
                       sheet = "Pre2_Bag1")
Bag1_Pre2_Clean <- dataclean(as.data.frame(Bag1_Pre2), as.data.frame(spp_list))

Bag2_Pre2 <- read_xlsx("GoldenrodData_Master.xlsx",
                       sheet = "Pre2_Bag2")
Bag2_Pre2_Clean <- dataclean(as.data.frame(Bag2_Pre2), as.data.frame(spp_list))

Bag3_Pre2 <- read_xlsx("GoldenrodData_Master.xlsx",
                       sheet = "Pre2_Bag3")
Bag3_Pre2_Clean <- dataclean(as.data.frame(Bag3_Pre2), as.data.frame(spp_list))

Bag4_Pre2 <- read_xlsx("GoldenrodData_Master.xlsx",
                       sheet = "Pre2_Bag4")
Bag4_Pre2_Clean <- dataclean(as.data.frame(Bag4_Pre2), as.data.frame(spp_list))

Bag6_Pre2 <- read_xlsx("GoldenrodData_Master.xlsx",
                       sheet = "Pre2_Bag6")
Bag6_Pre2_Clean <- dataclean(as.data.frame(Bag6_Pre2), as.data.frame(spp_list))

Bag7_Pre2 <- read_xlsx("GoldenrodData_Master.xlsx",
                       sheet = "Pre2_Bag7")
Bag7_Pre2_Clean <- dataclean(as.data.frame(Bag7_Pre2), as.data.frame(spp_list))

Bag8_Pre2 <- read_xlsx("GoldenrodData_Master.xlsx",
                       sheet = "Pre2_Bag8")
Bag8_Pre2_Clean <- dataclean(as.data.frame(Bag8_Pre2), as.data.frame(spp_list))

Bag9_Pre2 <- read_xlsx("GoldenrodData_Master.xlsx",
                       sheet = "Pre2_Bag9")
Bag9_Pre2_Clean <- dataclean(as.data.frame(Bag9_Pre2), as.data.frame(spp_list))

Bag10_Pre2 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Pre2_Bag10")
Bag10_Pre2_Clean <- dataclean(as.data.frame(Bag10_Pre2), as.data.frame(spp_list))

Bag11_Pre2 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Pre2_Bag11")
Bag11_Pre2_Clean <- dataclean(as.data.frame(Bag11_Pre2), as.data.frame(spp_list))

Bag12_Pre2 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Pre2_Bag12")
Bag12_Pre2_Clean <- dataclean(as.data.frame(Bag12_Pre2), as.data.frame(spp_list))

Bag13_Pre2 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Pre2_Bag13")
Bag13_Pre2_Clean <- dataclean(as.data.frame(Bag13_Pre2), as.data.frame(spp_list))

Bag14_Pre2 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Pre2_Bag14")
Bag14_Pre2_Clean <- dataclean(as.data.frame(Bag14_Pre2), as.data.frame(spp_list))

Bag15_Pre2 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Pre2_Bag15")
Bag15_Pre2_Clean <- dataclean(as.data.frame(Bag15_Pre2), as.data.frame(spp_list))

Bag16_Pre2 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Pre2_Bag16")
Bag16_Pre2_Clean <- dataclean(as.data.frame(Bag16_Pre2), as.data.frame(spp_list))

Bag17_Pre2 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Pre2_Bag17")
Bag17_Pre2_Clean <- dataclean(as.data.frame(Bag17_Pre2), as.data.frame(spp_list))

Bag18_Pre2 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Pre2_Bag18")
Bag18_Pre2_Clean <- dataclean(as.data.frame(Bag18_Pre2), as.data.frame(spp_list))

Bag19_Pre2 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Pre2_Bag19")
Bag19_Pre2_Clean <- dataclean(as.data.frame(Bag19_Pre2), as.data.frame(spp_list))

Bag20_Pre2 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Pre2_Bag20")
Bag20_Pre2_Clean <- dataclean(as.data.frame(Bag20_Pre2), as.data.frame(spp_list))

### Pre 3 
Bag1_Pre3 <- read_xlsx("GoldenrodData_Master.xlsx",
                       sheet = "Pre3_Bag1")
Bag1_Pre3_Clean <- dataclean(as.data.frame(Bag1_Pre3), as.data.frame(spp_list))

Bag2_Pre3 <- read_xlsx("GoldenrodData_Master.xlsx",
                       sheet = "Pre3_Bag2")
Bag2_Pre3_Clean <- dataclean(as.data.frame(Bag2_Pre3), as.data.frame(spp_list))

Bag3_Pre3 <- read_xlsx("GoldenrodData_Master.xlsx",
                       sheet = "Pre3_Bag3")
Bag3_Pre3_Clean <- dataclean(as.data.frame(Bag3_Pre3), as.data.frame(spp_list))

Bag4_Pre3 <- read_xlsx("GoldenrodData_Master.xlsx",
                       sheet = "Pre3_Bag4")
Bag4_Pre3_Clean <- dataclean(as.data.frame(Bag4_Pre3), as.data.frame(spp_list))

Bag6_Pre3 <- read_xlsx("GoldenrodData_Master.xlsx",
                       sheet = "Pre3_Bag6")
Bag6_Pre3_Clean <- dataclean(as.data.frame(Bag6_Pre3), as.data.frame(spp_list))

Bag7_Pre3 <- read_xlsx("GoldenrodData_Master.xlsx",
                       sheet = "Pre3_Bag7")
Bag7_Pre3_Clean <- dataclean(as.data.frame(Bag7_Pre3), as.data.frame(spp_list))

Bag8_Pre3 <- read_xlsx("GoldenrodData_Master.xlsx",
                       sheet = "Pre3_Bag8")
Bag8_Pre3_Clean <- dataclean(as.data.frame(Bag8_Pre3), as.data.frame(spp_list))

Bag9_Pre3 <- read_xlsx("GoldenrodData_Master.xlsx",
                       sheet = "Pre3_Bag9")
Bag9_Pre3_Clean <- dataclean(as.data.frame(Bag9_Pre3), as.data.frame(spp_list))

Bag10_Pre3 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Pre3_Bag10")
Bag10_Pre3_Clean <- dataclean(as.data.frame(Bag10_Pre3), as.data.frame(spp_list))

Bag11_Pre3 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Pre3_Bag11")
Bag11_Pre3_Clean <- dataclean(as.data.frame(Bag11_Pre3), as.data.frame(spp_list))

Bag12_Pre3 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Pre3_Bag12")
Bag12_Pre3_Clean <- dataclean(as.data.frame(Bag12_Pre3), as.data.frame(spp_list))

Bag13_Pre3 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Pre3_Bag13")
Bag13_Pre3_Clean <- dataclean(as.data.frame(Bag13_Pre3), as.data.frame(spp_list))

Bag14_Pre3 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Pre3_Bag14")
Bag14_Pre3_Clean <- dataclean(as.data.frame(Bag14_Pre3), as.data.frame(spp_list))

Bag15_Pre3 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Pre3_Bag15")
Bag15_Pre3_Clean <- dataclean(as.data.frame(Bag15_Pre3), as.data.frame(spp_list))

Bag16_Pre3 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Pre3_Bag16")
Bag16_Pre3_Clean <- dataclean(as.data.frame(Bag16_Pre3), as.data.frame(spp_list))

Bag17_Pre3 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Pre3_Bag17")
Bag17_Pre3_Clean <- dataclean(as.data.frame(Bag17_Pre3), as.data.frame(spp_list))

Bag18_Pre3 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Pre3_Bag18")
Bag18_Pre3_Clean <- dataclean(as.data.frame(Bag18_Pre3), as.data.frame(spp_list))

Bag19_Pre3 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Pre3_Bag19")
Bag19_Pre3_Clean <- dataclean(as.data.frame(Bag19_Pre3), as.data.frame(spp_list))

Bag20_Pre3 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Pre3_Bag20")
Bag20_Pre3_Clean <- dataclean(as.data.frame(Bag20_Pre3), as.data.frame(spp_list))

### Post 1.1 
Bag1_Post1 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Post1_Bag1")
Bag1_Post1_Clean <- dataclean(as.data.frame(Bag1_Post1), as.data.frame(spp_list))

Bag2_Post1 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Post1_Bag2")
Bag2_Post1_Clean <- dataclean(as.data.frame(Bag2_Post1), as.data.frame(spp_list))

Bag3_Post1 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Post1_Bag3")
Bag3_Post1_Clean <- dataclean(as.data.frame(Bag3_Post1), as.data.frame(spp_list))

Bag4_Post1 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Post1_Bag4")
Bag4_Post1_Clean <- dataclean(as.data.frame(Bag4_Post1), as.data.frame(spp_list))

Bag6_Post1 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Post1_Bag6")
Bag6_Post1_Clean <- dataclean(as.data.frame(Bag6_Post1), as.data.frame(spp_list))

Bag7_Post1 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Post1_Bag7")
Bag7_Post1_Clean <- dataclean(as.data.frame(Bag7_Post1), as.data.frame(spp_list))

Bag8_Post1 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Post1_Bag8")
Bag8_Post1_Clean <- dataclean(as.data.frame(Bag8_Post1), as.data.frame(spp_list))

Bag9_Post1 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Post1_Bag9")
Bag9_Post1_Clean <- dataclean(as.data.frame(Bag9_Post1), as.data.frame(spp_list))

Bag10_Post1 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post1_Bag10")
Bag10_Post1_Clean <- dataclean(as.data.frame(Bag10_Post1), as.data.frame(spp_list))

Bag11_Post1 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post1_Bag11")
Bag11_Post1_Clean <- dataclean(as.data.frame(Bag11_Post1), as.data.frame(spp_list))

Bag12_Post1 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post1_Bag12")
Bag12_Post1_Clean <- dataclean(as.data.frame(Bag12_Post1), as.data.frame(spp_list))

Bag13_Post1 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post1_Bag13")
Bag13_Post1_Clean <- dataclean(as.data.frame(Bag13_Post1), as.data.frame(spp_list))

Bag14_Post1 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post1_Bag14")
Bag14_Post1_Clean <- dataclean(as.data.frame(Bag14_Post1), as.data.frame(spp_list))

Bag15_Post1 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post1_Bag15")
Bag15_Post1_Clean <- dataclean(as.data.frame(Bag15_Post1), as.data.frame(spp_list))

Bag16_Post1 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post1_Bag16")
Bag16_Post1_Clean <- dataclean(as.data.frame(Bag16_Post1), as.data.frame(spp_list))

Bag17_Post1 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post1_Bag17")
Bag17_Post1_Clean <- dataclean(as.data.frame(Bag17_Post1), as.data.frame(spp_list))

Bag18_Post1 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post1_Bag18")
Bag18_Post1_Clean <- dataclean(as.data.frame(Bag18_Post1), as.data.frame(spp_list))

Bag19_Post1 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post1_Bag19")
Bag19_Post1_Clean <- dataclean(as.data.frame(Bag19_Post1), as.data.frame(spp_list))

Bag20_Post1 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post1_Bag20")
Bag20_Post1_Clean <- dataclean(as.data.frame(Bag20_Post1), as.data.frame(spp_list))

### Post 1.2 
Bag1_Post2 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Post2_Bag1")
Bag1_Post2_Clean <- dataclean(as.data.frame(Bag1_Post2), as.data.frame(spp_list))

Bag2_Post2 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Post2_Bag2")
Bag2_Post2_Clean <- dataclean(as.data.frame(Bag2_Post2), as.data.frame(spp_list))

Bag3_Post2 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Post2_Bag3")
Bag3_Post2_Clean <- dataclean(as.data.frame(Bag3_Post2), as.data.frame(spp_list))

Bag4_Post2 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Post2_Bag4")
Bag4_Post2_Clean <- dataclean(as.data.frame(Bag4_Post2), as.data.frame(spp_list))

Bag6_Post2 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Post2_Bag6")
Bag6_Post2_Clean <- dataclean(as.data.frame(Bag6_Post2), as.data.frame(spp_list))

Bag7_Post2 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Post2_Bag7")
Bag7_Post2_Clean <- dataclean(as.data.frame(Bag7_Post2), as.data.frame(spp_list))

Bag8_Post2 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Post2_Bag8")
Bag8_Post2_Clean <- dataclean(as.data.frame(Bag8_Post2), as.data.frame(spp_list))

Bag9_Post2 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Post2_Bag9")
Bag9_Post2_Clean <- dataclean(as.data.frame(Bag9_Post2), as.data.frame(spp_list))

Bag10_Post2 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post2_Bag10")
Bag10_Post2_Clean <- dataclean(as.data.frame(Bag10_Post2), as.data.frame(spp_list))

Bag11_Post2 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post2_Bag11")
Bag11_Post2_Clean <- dataclean(as.data.frame(Bag11_Post2), as.data.frame(spp_list))

Bag12_Post2 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post2_Bag12")
Bag12_Post2_Clean <- dataclean(as.data.frame(Bag12_Post2), as.data.frame(spp_list))

Bag13_Post2 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post2_Bag13")
Bag13_Post2_Clean <- dataclean(as.data.frame(Bag13_Post2), as.data.frame(spp_list))

Bag14_Post2 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post2_Bag14")
Bag14_Post2_Clean <- dataclean(as.data.frame(Bag14_Post2), as.data.frame(spp_list))

Bag15_Post2 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post2_Bag15")
Bag15_Post2_Clean <- dataclean(as.data.frame(Bag15_Post2), as.data.frame(spp_list))

Bag16_Post2 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post2_Bag16")
Bag16_Post2_Clean <- dataclean(as.data.frame(Bag16_Post2), as.data.frame(spp_list))

Bag17_Post2 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post2_Bag17")
Bag17_Post2_Clean <- dataclean(as.data.frame(Bag17_Post2), as.data.frame(spp_list))

Bag18_Post2 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post2_Bag18")
Bag18_Post2_Clean <- dataclean(as.data.frame(Bag18_Post2), as.data.frame(spp_list))

Bag19_Post2 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post2_Bag19")
Bag19_Post2_Clean <- dataclean(as.data.frame(Bag19_Post2), as.data.frame(spp_list))

Bag20_Post2 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post2_Bag20")
Bag20_Post2_Clean <- dataclean(as.data.frame(Bag20_Post2), as.data.frame(spp_list))

### Post 1.3 
Bag1_Post3 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Post3_Bag1")
Bag1_Post3_Clean <- dataclean(as.data.frame(Bag1_Post3), as.data.frame(spp_list))

Bag2_Post3 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Post3_Bag2")
Bag2_Post3_Clean <- dataclean(as.data.frame(Bag2_Post3), as.data.frame(spp_list))

Bag3_Post3 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Post3_Bag3")
Bag3_Post3_Clean <- dataclean(as.data.frame(Bag3_Post3), as.data.frame(spp_list))

Bag4_Post3 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Post3_Bag4")
Bag4_Post3_Clean <- dataclean(as.data.frame(Bag4_Post3), as.data.frame(spp_list))

Bag6_Post3 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Post3_Bag6")
Bag6_Post3_Clean <- dataclean(as.data.frame(Bag6_Post3), as.data.frame(spp_list))

Bag7_Post3 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Post3_Bag7")
Bag7_Post3_Clean <- dataclean(as.data.frame(Bag7_Post3), as.data.frame(spp_list))

Bag8_Post3 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Post3_Bag8")
Bag8_Post3_Clean <- dataclean(as.data.frame(Bag8_Post3), as.data.frame(spp_list))

Bag9_Post3 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Post3_Bag9")
Bag9_Post3_Clean <- dataclean(as.data.frame(Bag9_Post3), as.data.frame(spp_list))

Bag10_Post3 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post3_Bag10")
Bag10_Post3_Clean <- dataclean(as.data.frame(Bag10_Post3), as.data.frame(spp_list))

Bag11_Post3 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post3_Bag11")
Bag11_Post3_Clean <- dataclean(as.data.frame(Bag11_Post3), as.data.frame(spp_list))

Bag12_Post3 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post3_Bag12")
Bag12_Post3_Clean <- dataclean(as.data.frame(Bag12_Post3), as.data.frame(spp_list))

Bag13_Post3 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post3_Bag13")
Bag13_Post3_Clean <- dataclean(as.data.frame(Bag13_Post3), as.data.frame(spp_list))

Bag14_Post3 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post3_Bag14")
Bag14_Post3_Clean <- dataclean(as.data.frame(Bag14_Post3), as.data.frame(spp_list))

Bag15_Post3 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post3_Bag15")
Bag15_Post3_Clean <- dataclean(as.data.frame(Bag15_Post3), as.data.frame(spp_list))

Bag16_Post3 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post3_Bag16")
Bag16_Post3_Clean <- dataclean(as.data.frame(Bag16_Post3), as.data.frame(spp_list))

Bag17_Post3 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post3_Bag17")
Bag17_Post3_Clean <- dataclean(as.data.frame(Bag17_Post3), as.data.frame(spp_list))

Bag18_Post3 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post3_Bag18")
Bag18_Post3_Clean <- dataclean(as.data.frame(Bag18_Post3), as.data.frame(spp_list))

Bag19_Post3 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post3_Bag19")
Bag19_Post3_Clean <- dataclean(as.data.frame(Bag19_Post3), as.data.frame(spp_list))

Bag20_Post3 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post3_Bag20")
Bag20_Post3_Clean <- dataclean(as.data.frame(Bag20_Post3), as.data.frame(spp_list))

### Post 2.1 
Bag1_Post4 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Post4_Bag1")
Bag1_Post4_Clean <- dataclean(as.data.frame(Bag1_Post4), as.data.frame(spp_list))

Bag2_Post4 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Post4_Bag2")
Bag2_Post4_Clean <- dataclean(as.data.frame(Bag2_Post4), as.data.frame(spp_list))

Bag3_Post4 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Post4_Bag3")
Bag3_Post4_Clean <- dataclean(as.data.frame(Bag3_Post4), as.data.frame(spp_list))

Bag4_Post4 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Post4_Bag4")
Bag4_Post4_Clean <- dataclean(as.data.frame(Bag4_Post4), as.data.frame(spp_list))

Bag6_Post4 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Post4_Bag6")
Bag6_Post4_Clean <- dataclean(as.data.frame(Bag6_Post4), as.data.frame(spp_list))

Bag7_Post4 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Post4_Bag7")
Bag7_Post4_Clean <- dataclean(as.data.frame(Bag7_Post4), as.data.frame(spp_list))

Bag8_Post4 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Post4_Bag8")
Bag8_Post4_Clean <- dataclean(as.data.frame(Bag8_Post4), as.data.frame(spp_list))

Bag9_Post4 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Post4_Bag9")
Bag9_Post4_Clean <- dataclean(as.data.frame(Bag9_Post4), as.data.frame(spp_list))

Bag10_Post4 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post4_Bag10")
Bag10_Post4_Clean <- dataclean(as.data.frame(Bag10_Post4), as.data.frame(spp_list))

Bag11_Post4 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post4_Bag11")
Bag11_Post4_Clean <- dataclean(as.data.frame(Bag11_Post4), as.data.frame(spp_list))

Bag12_Post4 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post4_Bag12")
Bag12_Post4_Clean <- dataclean(as.data.frame(Bag12_Post4), as.data.frame(spp_list))

Bag13_Post4 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post4_Bag13")
Bag13_Post4_Clean <- dataclean(as.data.frame(Bag13_Post4), as.data.frame(spp_list))

Bag14_Post4 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post4_Bag14")
Bag14_Post4_Clean <- dataclean(as.data.frame(Bag14_Post4), as.data.frame(spp_list))

Bag15_Post4 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post4_Bag15")
Bag15_Post4_Clean <- dataclean(as.data.frame(Bag15_Post4), as.data.frame(spp_list))

Bag16_Post4 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post4_Bag16")
Bag16_Post4_Clean <- dataclean(as.data.frame(Bag16_Post4), as.data.frame(spp_list))

Bag17_Post4 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post4_Bag17")
Bag17_Post4_Clean <- dataclean(as.data.frame(Bag17_Post4), as.data.frame(spp_list))

Bag18_Post4 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post4_Bag18")
Bag18_Post4_Clean <- dataclean(as.data.frame(Bag18_Post4), as.data.frame(spp_list))

Bag19_Post4 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post4_Bag19")
Bag19_Post4_Clean <- dataclean(as.data.frame(Bag19_Post4), as.data.frame(spp_list))

Bag20_Post4 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post4_Bag20")
Bag20_Post4_Clean <- dataclean(as.data.frame(Bag20_Post4), as.data.frame(spp_list))

### Post 2.2 
Bag1_Post5 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Post5_Bag1")
Bag1_Post5_Clean <- dataclean(as.data.frame(Bag1_Post5), as.data.frame(spp_list))

Bag2_Post5 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Post5_Bag2")
Bag2_Post5_Clean <- dataclean(as.data.frame(Bag2_Post5), as.data.frame(spp_list))

Bag3_Post5 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Post5_Bag3")
Bag3_Post5_Clean <- dataclean(as.data.frame(Bag3_Post5), as.data.frame(spp_list))

Bag4_Post5 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Post5_Bag4")
Bag4_Post5_Clean <- dataclean(as.data.frame(Bag4_Post5), as.data.frame(spp_list))

Bag6_Post5 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Post5_Bag6")
Bag6_Post5_Clean <- dataclean(as.data.frame(Bag6_Post5), as.data.frame(spp_list))

Bag7_Post5 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Post5_Bag7")
Bag7_Post5_Clean <- dataclean(as.data.frame(Bag7_Post5), as.data.frame(spp_list))

Bag8_Post5 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Post5_Bag8")
Bag8_Post5_Clean <- dataclean(as.data.frame(Bag8_Post5), as.data.frame(spp_list))

Bag9_Post5 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Post5_Bag9")
Bag9_Post5_Clean <- dataclean(as.data.frame(Bag9_Post5), as.data.frame(spp_list))

Bag10_Post5 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post5_Bag10")
Bag10_Post5_Clean <- dataclean(as.data.frame(Bag10_Post5), as.data.frame(spp_list))

Bag11_Post5 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post5_Bag11")
Bag11_Post5_Clean <- dataclean(as.data.frame(Bag11_Post5), as.data.frame(spp_list))

Bag12_Post5 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post5_Bag12")
Bag12_Post5_Clean <- dataclean(as.data.frame(Bag12_Post5), as.data.frame(spp_list))

Bag13_Post5 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post5_Bag13")
Bag13_Post5_Clean <- dataclean(as.data.frame(Bag13_Post5), as.data.frame(spp_list))

Bag14_Post5 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post5_Bag14")
Bag14_Post5_Clean <- dataclean(as.data.frame(Bag14_Post5), as.data.frame(spp_list))

Bag15_Post5 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post5_Bag15")
Bag15_Post5_Clean <- dataclean(as.data.frame(Bag15_Post5), as.data.frame(spp_list))

Bag16_Post5 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post5_Bag16")
Bag16_Post5_Clean <- dataclean(as.data.frame(Bag16_Post5), as.data.frame(spp_list))

Bag17_Post5 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post5_Bag17")
Bag17_Post5_Clean <- dataclean(as.data.frame(Bag17_Post5), as.data.frame(spp_list))

Bag18_Post5 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post5_Bag18")
Bag18_Post5_Clean <- dataclean(as.data.frame(Bag18_Post5), as.data.frame(spp_list))

Bag19_Post5 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post5_Bag19")
Bag19_Post5_Clean <- dataclean(as.data.frame(Bag19_Post5), as.data.frame(spp_list))

Bag20_Post5 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post5_Bag20")
Bag20_Post5_Clean <- dataclean(as.data.frame(Bag20_Post5), as.data.frame(spp_list))

### Post 2.3 
Bag1_Post6 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Post6_Bag1")
Bag1_Post6_Clean <- dataclean(as.data.frame(Bag1_Post6), as.data.frame(spp_list))

Bag2_Post6 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Post6_Bag2")
Bag2_Post6_Clean <- dataclean(as.data.frame(Bag2_Post6), as.data.frame(spp_list))

Bag3_Post6 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Post6_Bag3")
Bag3_Post6_Clean <- dataclean(as.data.frame(Bag3_Post6), as.data.frame(spp_list))

Bag4_Post6 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Post6_Bag4")
Bag4_Post6_Clean <- dataclean(as.data.frame(Bag4_Post6), as.data.frame(spp_list))

Bag6_Post6 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Post6_Bag6")
Bag6_Post6_Clean <- dataclean(as.data.frame(Bag6_Post6), as.data.frame(spp_list))

Bag7_Post6 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Post6_Bag7")
Bag7_Post6_Clean <- dataclean(as.data.frame(Bag7_Post6), as.data.frame(spp_list))

Bag8_Post6 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Post6_Bag8")
Bag8_Post6_Clean <- dataclean(as.data.frame(Bag8_Post6), as.data.frame(spp_list))

Bag9_Post6 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Post6_Bag9")
Bag9_Post6_Clean <- dataclean(as.data.frame(Bag9_Post6), as.data.frame(spp_list))

Bag10_Post6 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post6_Bag10")
Bag10_Post6_Clean <- dataclean(as.data.frame(Bag10_Post6), as.data.frame(spp_list))

Bag11_Post6 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post6_Bag11")
Bag11_Post6_Clean <- dataclean(as.data.frame(Bag11_Post6), as.data.frame(spp_list))

Bag12_Post6 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post6_Bag12")
Bag12_Post6_Clean <- dataclean(as.data.frame(Bag12_Post6), as.data.frame(spp_list))

Bag13_Post6 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post6_Bag13")
Bag13_Post6_Clean <- dataclean(as.data.frame(Bag13_Post6), as.data.frame(spp_list))

Bag14_Post6 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post6_Bag14")
Bag14_Post6_Clean <- dataclean(as.data.frame(Bag14_Post6), as.data.frame(spp_list))

Bag15_Post6 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post6_Bag15")
Bag15_Post6_Clean <- dataclean(as.data.frame(Bag15_Post6), as.data.frame(spp_list))

Bag16_Post6 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post6_Bag16")
Bag16_Post6_Clean <- dataclean(as.data.frame(Bag16_Post6), as.data.frame(spp_list))

Bag17_Post6 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post6_Bag17")
Bag17_Post6_Clean <- dataclean(as.data.frame(Bag17_Post6), as.data.frame(spp_list))

Bag18_Post6 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post6_Bag18")
Bag18_Post6_Clean <- dataclean(as.data.frame(Bag18_Post6), as.data.frame(spp_list))

Bag19_Post6 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post6_Bag19")
Bag19_Post6_Clean <- dataclean(as.data.frame(Bag19_Post6), as.data.frame(spp_list))

Bag20_Post6 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post6_Bag20")
Bag20_Post6_Clean <- dataclean(as.data.frame(Bag20_Post6), as.data.frame(spp_list))

### Post 3.1 
Bag1_Post7 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Post7_Bag1")
Bag1_Post7_Clean <- datacleanPost3(as.data.frame(Bag1_Post7), as.data.frame(spp_listPost3))

Bag2_Post7 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Post7_Bag2")
Bag2_Post7_Clean <- datacleanPost3(as.data.frame(Bag2_Post7), as.data.frame(spp_listPost3))

Bag3_Post7 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Post7_Bag3")
Bag3_Post7_Clean <- datacleanPost3(as.data.frame(Bag3_Post7), as.data.frame(spp_listPost3))

Bag4_Post7 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Post7_Bag4")
Bag4_Post7_Clean <- datacleanPost3(as.data.frame(Bag4_Post7), as.data.frame(spp_listPost3))

Bag6_Post7 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Post7_Bag6")
Bag6_Post7_Clean <- datacleanPost3(as.data.frame(Bag6_Post7), as.data.frame(spp_listPost3))

Bag7_Post7 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Post7_Bag7")
Bag7_Post7_Clean <- datacleanPost3(as.data.frame(Bag7_Post7), as.data.frame(spp_listPost3))

Bag8_Post7 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Post7_Bag8")
Bag8_Post7_Clean <- datacleanPost3(as.data.frame(Bag8_Post7), as.data.frame(spp_listPost3))

Bag9_Post7 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Post7_Bag9")
Bag9_Post7_Clean <- datacleanPost3(as.data.frame(Bag9_Post7), as.data.frame(spp_listPost3))

Bag10_Post7 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post7_Bag10")
Bag10_Post7_Clean <- datacleanPost3(as.data.frame(Bag10_Post7), as.data.frame(spp_listPost3))

Bag11_Post7 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post7_Bag11")
Bag11_Post7_Clean <- datacleanPost3(as.data.frame(Bag11_Post7), as.data.frame(spp_listPost3))

Bag12_Post7 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post7_Bag12")
Bag12_Post7_Clean <- datacleanPost3(as.data.frame(Bag12_Post7), as.data.frame(spp_listPost3))

Bag13_Post7 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post7_Bag13")
Bag13_Post7_Clean <- datacleanPost3(as.data.frame(Bag13_Post7), as.data.frame(spp_listPost3))

Bag14_Post7 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post7_Bag14")
Bag14_Post7_Clean <- datacleanPost3(as.data.frame(Bag14_Post7), as.data.frame(spp_listPost3))

Bag15_Post7 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post7_Bag15")
Bag15_Post7_Clean <- datacleanPost3(as.data.frame(Bag15_Post7), as.data.frame(spp_listPost3))

Bag16_Post7 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post7_Bag16")
Bag16_Post7_Clean <- datacleanPost3(as.data.frame(Bag16_Post7), as.data.frame(spp_listPost3))

Bag17_Post7 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post7_Bag17")
Bag17_Post7_Clean <- datacleanPost3(as.data.frame(Bag17_Post7), as.data.frame(spp_listPost3))

Bag18_Post7 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post7_Bag18")
Bag18_Post7_Clean <- datacleanPost3(as.data.frame(Bag18_Post7), as.data.frame(spp_listPost3))

Bag19_Post7 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post7_Bag19")
Bag19_Post7_Clean <- datacleanPost3(as.data.frame(Bag19_Post7), as.data.frame(spp_listPost3))

Bag20_Post7 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post7_Bag20")
Bag20_Post7_Clean <- datacleanPost3(as.data.frame(Bag20_Post7), as.data.frame(spp_listPost3))

### Post 3.2 
Bag1_Post8 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Post8_Bag1")
Bag1_Post8_Clean <- datacleanPost3(as.data.frame(Bag1_Post8), as.data.frame(spp_listPost3))

Bag2_Post8 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Post8_Bag2")
Bag2_Post8_Clean <- datacleanPost3(as.data.frame(Bag2_Post8), as.data.frame(spp_listPost3))

Bag3_Post8 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Post8_Bag3")
Bag3_Post8_Clean <- datacleanPost3(as.data.frame(Bag3_Post8), as.data.frame(spp_listPost3))

Bag4_Post8 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Post8_Bag4")
Bag4_Post8_Clean <- datacleanPost3(as.data.frame(Bag4_Post8), as.data.frame(spp_listPost3))

Bag6_Post8 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Post8_Bag6")
Bag6_Post8_Clean <- datacleanPost3(as.data.frame(Bag6_Post8), as.data.frame(spp_listPost3))

Bag7_Post8 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Post8_Bag7")
Bag7_Post8_Clean <- datacleanPost3(as.data.frame(Bag7_Post8), as.data.frame(spp_listPost3))

Bag8_Post8 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Post8_Bag8")
Bag8_Post8_Clean <- datacleanPost3(as.data.frame(Bag8_Post8), as.data.frame(spp_listPost3))

Bag9_Post8 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Post8_Bag9")
Bag9_Post8_Clean <- datacleanPost3(as.data.frame(Bag9_Post8), as.data.frame(spp_listPost3))

Bag10_Post8 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post8_Bag10")
Bag10_Post8_Clean <- datacleanPost3(as.data.frame(Bag10_Post8), as.data.frame(spp_listPost3))

Bag11_Post8 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post8_Bag11")
Bag11_Post8_Clean <- datacleanPost3(as.data.frame(Bag11_Post8), as.data.frame(spp_listPost3))

Bag12_Post8 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post8_Bag12")
Bag12_Post8_Clean <- datacleanPost3(as.data.frame(Bag12_Post8), as.data.frame(spp_listPost3))

Bag13_Post8 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post8_Bag13")
Bag13_Post8_Clean <- datacleanPost3(as.data.frame(Bag13_Post8), as.data.frame(spp_listPost3))

Bag14_Post8 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post8_Bag14")
Bag14_Post8_Clean <- datacleanPost3(as.data.frame(Bag14_Post8), as.data.frame(spp_listPost3))

Bag15_Post8 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post8_Bag15")
Bag15_Post8_Clean <- datacleanPost3(as.data.frame(Bag15_Post8), as.data.frame(spp_listPost3))

Bag16_Post8 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post8_Bag16")
Bag16_Post8_Clean <- datacleanPost3(as.data.frame(Bag16_Post8), as.data.frame(spp_listPost3))

Bag17_Post8 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post8_Bag17")
Bag17_Post8_Clean <- datacleanPost3(as.data.frame(Bag17_Post8), as.data.frame(spp_listPost3))

Bag18_Post8 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post8_Bag18")
Bag18_Post8_Clean <- datacleanPost3(as.data.frame(Bag18_Post8), as.data.frame(spp_listPost3))

Bag19_Post8 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post8_Bag19")
Bag19_Post8_Clean <- datacleanPost3(as.data.frame(Bag19_Post8), as.data.frame(spp_listPost3))

Bag20_Post8 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post8_Bag20")
Bag20_Post8_Clean <- datacleanPost3(as.data.frame(Bag20_Post8), as.data.frame(spp_listPost3))

### Post 3.3 
Bag1_Post9 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Post9_Bag1")
Bag1_Post9_Clean <- datacleanPost3(as.data.frame(Bag1_Post9), as.data.frame(spp_listPost3))

Bag2_Post9 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Post9_Bag2")
Bag2_Post9_Clean <- datacleanPost3(as.data.frame(Bag2_Post9), as.data.frame(spp_listPost3))

Bag3_Post9 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Post9_Bag3")
Bag3_Post9_Clean <- datacleanPost3(as.data.frame(Bag3_Post9), as.data.frame(spp_listPost3))

Bag4_Post9 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Post9_Bag4")
Bag4_Post9_Clean <- datacleanPost3(as.data.frame(Bag4_Post9), as.data.frame(spp_listPost3))

Bag6_Post9 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Post9_Bag6")
Bag6_Post9_Clean <- datacleanPost3(as.data.frame(Bag6_Post9), as.data.frame(spp_listPost3))

Bag7_Post9 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Post9_Bag7")
Bag7_Post9_Clean <- datacleanPost3(as.data.frame(Bag7_Post9), as.data.frame(spp_listPost3))

Bag8_Post9 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Post9_Bag8")
Bag8_Post9_Clean <- datacleanPost3(as.data.frame(Bag8_Post9), as.data.frame(spp_listPost3))

Bag9_Post9 <- read_xlsx("GoldenrodData_Master.xlsx",
                        sheet = "Post9_Bag9")
Bag9_Post9_Clean <- datacleanPost3(as.data.frame(Bag9_Post9), as.data.frame(spp_listPost3))

Bag10_Post9 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post9_Bag10")
Bag10_Post9_Clean <- datacleanPost3(as.data.frame(Bag10_Post9), as.data.frame(spp_listPost3))

Bag11_Post9 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post9_Bag11")
Bag11_Post9_Clean <- datacleanPost3(as.data.frame(Bag11_Post9), as.data.frame(spp_listPost3))

Bag12_Post9 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post9_Bag12")
Bag12_Post9_Clean <- datacleanPost3(as.data.frame(Bag12_Post9), as.data.frame(spp_listPost3))

Bag13_Post9 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post9_Bag13")
Bag13_Post9_Clean <- datacleanPost3(as.data.frame(Bag13_Post9), as.data.frame(spp_listPost3))

Bag14_Post9 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post9_Bag14")
Bag14_Post9_Clean <- datacleanPost3(as.data.frame(Bag14_Post9), as.data.frame(spp_listPost3))

Bag15_Post9 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post9_Bag15")
Bag15_Post9_Clean <- datacleanPost3(as.data.frame(Bag15_Post9), as.data.frame(spp_listPost3))

Bag16_Post9 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post9_Bag16")
Bag16_Post9_Clean <- datacleanPost3(as.data.frame(Bag16_Post9), as.data.frame(spp_listPost3))

Bag17_Post9 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post9_Bag17")
Bag17_Post9_Clean <- datacleanPost3(as.data.frame(Bag17_Post9), as.data.frame(spp_listPost3))

Bag18_Post9 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post9_Bag18")
Bag18_Post9_Clean <- datacleanPost3(as.data.frame(Bag18_Post9), as.data.frame(spp_listPost3))

Bag19_Post9 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post9_Bag19")
Bag19_Post9_Clean <- datacleanPost3(as.data.frame(Bag19_Post9), as.data.frame(spp_listPost3))

Bag20_Post9 <- read_xlsx("GoldenrodData_Master.xlsx",
                         sheet = "Post9_Bag20")
Bag20_Post9_Clean <- datacleanPost3(as.data.frame(Bag20_Post9), as.data.frame(spp_listPost3))

##Sum across three days of sampling to get one sp x site matrix for each survey

### Pre Totals ----

#For each site, sum across three days of pre-treatment sampling to get total for pre-treatment sampling period

Bag1_PreTotal_Clean <- Bag1_Pre1_Clean +
  Bag1_Pre2_Clean + Bag1_Pre3_Clean

#Bag1 has 3 extra plants. Let's randomly remove 3. Columns 3, 9, 17 were chosen.
Bag1_PreTotal_Clean <- select(Bag1_PreTotal_Clean, -c(3,9,17))

Bag2_PreTotal_Clean <- Bag2_Pre1_Clean +
  Bag2_Pre2_Clean + Bag2_Pre3_Clean

#Bag1 has 3 extra plants. Let's randomly remove one that was not removed during treatment. Column 12 was chosen.
Bag2_PreTotal_Clean <- select(Bag2_PreTotal_Clean, -12)

Bag3_PreTotal_Clean <- Bag3_Pre1_Clean +
  Bag3_Pre2_Clean + Bag3_Pre3_Clean

Bag4_PreTotal_Clean <- Bag4_Pre1_Clean +
  Bag4_Pre2_Clean + Bag4_Pre3_Clean

Bag6_PreTotal_Clean <- Bag6_Pre1_Clean +
  Bag6_Pre2_Clean + Bag6_Pre3_Clean

Bag7_PreTotal_Clean <- Bag7_Pre1_Clean +
  Bag7_Pre2_Clean + Bag7_Pre3_Clean

Bag8_PreTotal_Clean <- Bag8_Pre1_Clean +
  Bag8_Pre2_Clean + Bag8_Pre3_Clean

Bag9_PreTotal_Clean <- Bag9_Pre1_Clean +
  Bag9_Pre2_Clean + Bag9_Pre3_Clean

Bag10_PreTotal_Clean <- Bag10_Pre1_Clean +
  Bag10_Pre2_Clean + Bag10_Pre3_Clean

Bag11_PreTotal_Clean <- Bag11_Pre1_Clean +
  Bag11_Pre2_Clean + Bag11_Pre3_Clean

Bag12_PreTotal_Clean <- Bag12_Pre1_Clean +
  Bag12_Pre2_Clean + Bag12_Pre3_Clean

Bag13_PreTotal_Clean <- Bag13_Pre1_Clean +
  Bag13_Pre2_Clean + Bag13_Pre3_Clean

Bag14_PreTotal_Clean <- Bag14_Pre1_Clean +
  Bag14_Pre2_Clean + Bag14_Pre3_Clean

Bag15_PreTotal_Clean <- Bag15_Pre1_Clean +
  Bag15_Pre2_Clean + Bag15_Pre3_Clean

Bag16_PreTotal_Clean <- Bag16_Pre1_Clean +
  Bag16_Pre2_Clean + Bag16_Pre3_Clean

Bag17_PreTotal_Clean <- Bag17_Pre1_Clean +
  Bag17_Pre2_Clean + Bag17_Pre3_Clean

#Bag17 has one extra plant. Let's remove one randomly. Column 2 was chosen.
Bag17_PreTotal_Clean <- select(Bag17_PreTotal_Clean, -2)

Bag18_PreTotal_Clean <- Bag18_Pre1_Clean +
  Bag18_Pre2_Clean + Bag18_Pre3_Clean

Bag19_PreTotal_Clean <- Bag19_Pre1_Clean +
  Bag19_Pre2_Clean + Bag19_Pre3_Clean

Bag20_PreTotal_Clean <- Bag20_Pre1_Clean +
  Bag20_Pre2_Clean + Bag20_Pre3_Clean

##Group sites by treatment for Pre survey
PreTotal_Treatments <- list(Bag2_PreTotal_Clean, Bag4_PreTotal_Clean, Bag7_PreTotal_Clean,
                            Bag9_PreTotal_Clean, Bag10_PreTotal_Clean, Bag12_PreTotal_Clean,
                            Bag13_PreTotal_Clean, Bag16_PreTotal_Clean, Bag18_PreTotal_Clean,
                            Bag20_PreTotal_Clean)

#Group sites by control for Pre survey
PreTotal_Controls <- list(Bag1_PreTotal_Clean, Bag3_PreTotal_Clean, Bag6_PreTotal_Clean,
                          Bag8_PreTotal_Clean, Bag11_PreTotal_Clean, Bag14_PreTotal_Clean,
                          Bag15_PreTotal_Clean, Bag17_PreTotal_Clean, Bag19_PreTotal_Clean)

##Create regional treatment spp x site matrix, where each column is a pre-destruction treatment site
Pre_Treatment_Region <- matrix(nrow = nrow(Bag2_PreTotal_Clean), ncol = length(PreTotal_Treatments))
for (i in 1:length(PreTotal_Treatments)){
  Pre_Treatment_Region[,i] <- rowSums(PreTotal_Treatments[[i]])
}
rownames(Pre_Treatment_Region) <- rownames(Bag2_PreTotal_Clean) #add species names
colnames(Pre_Treatment_Region) <- c("Bag2", "Bag4", "Bag7", "Bag9", "Bag10", #name columns with corresponding sites
                                    "Bag12", "Bag13", "Bag16", "Bag18", "Bag20")

##Create regional control spp x site matrix, where each column is a pre-destruction control site
Pre_Controls_Region <- matrix(nrow = nrow(Bag1_PreTotal_Clean), ncol = length(PreTotal_Controls))
for (i in 1:length(PreTotal_Controls)){
  Pre_Controls_Region[,i] <- rowSums(PreTotal_Controls[[i]])
}
rownames(Pre_Controls_Region) <- rownames(Bag1_PreTotal_Clean)#add species names
colnames(Pre_Controls_Region) <- c("Bag1", "Bag3", "Bag6", "Bag8", "Bag11", #name columns with corresponding sites
                                   "Bag14", "Bag15", "Bag17", "Bag19")

### Post1 Totals ----

#For each site, sum across three days of immediate post-treatment sampling to get total for Post1 sampling period

Bag1_Post1Total_Clean <- Bag1_Post1_Clean +
  Bag1_Post2_Clean + Bag1_Post3_Clean

#Remove extra plants from Bag1
Bag1_Post1Total_Clean <- select(Bag1_Post1Total_Clean, -c(3,9,17))

Bag2_Post1Total_Clean <- Bag2_Post1_Clean +
  Bag2_Post2_Clean + Bag2_Post3_Clean

Bag3_Post1Total_Clean <- Bag3_Post1_Clean +
  Bag3_Post2_Clean + Bag3_Post3_Clean

Bag4_Post1Total_Clean <- Bag4_Post1_Clean +
  Bag4_Post2_Clean + Bag4_Post3_Clean

Bag6_Post1Total_Clean <- Bag6_Post1_Clean +
  Bag6_Post2_Clean + Bag6_Post3_Clean

Bag7_Post1Total_Clean <- Bag7_Post1_Clean +
  Bag7_Post2_Clean + Bag7_Post3_Clean

Bag8_Post1Total_Clean <- Bag8_Post1_Clean +
  Bag8_Post2_Clean + Bag8_Post3_Clean

Bag9_Post1Total_Clean <- Bag9_Post1_Clean +
  Bag9_Post2_Clean + Bag9_Post3_Clean

Bag10_Post1Total_Clean <- Bag10_Post1_Clean +
  Bag10_Post2_Clean + Bag10_Post3_Clean

Bag11_Post1Total_Clean <- Bag11_Post1_Clean +
  Bag11_Post2_Clean + Bag11_Post3_Clean

Bag12_Post1Total_Clean <- Bag12_Post1_Clean +
  Bag12_Post2_Clean + Bag12_Post3_Clean

Bag13_Post1Total_Clean <- Bag13_Post1_Clean +
  Bag13_Post2_Clean + Bag13_Post3_Clean

Bag14_Post1Total_Clean <- Bag14_Post1_Clean +
  Bag14_Post2_Clean + Bag14_Post3_Clean

Bag15_Post1Total_Clean <- Bag15_Post1_Clean +
  Bag15_Post2_Clean + Bag15_Post3_Clean

Bag16_Post1Total_Clean <- Bag16_Post1_Clean +
  Bag16_Post2_Clean + Bag16_Post3_Clean

Bag17_Post1Total_Clean <- Bag17_Post1_Clean +
  Bag17_Post2_Clean + Bag17_Post3_Clean

#Remove extra plant from Bag17
Bag17_Post1Total_Clean <- select(Bag17_Post1Total_Clean, -2)

Bag18_Post1Total_Clean <- Bag18_Post1_Clean +
  Bag18_Post2_Clean + Bag18_Post3_Clean

Bag19_Post1Total_Clean <- Bag19_Post1_Clean +
  Bag19_Post2_Clean + Bag19_Post3_Clean

Bag20_Post1Total_Clean <- Bag20_Post1_Clean +
  Bag20_Post2_Clean + Bag20_Post3_Clean

Post1Total_Treatments <- list(Bag2_Post1Total_Clean, Bag4_Post1Total_Clean, Bag7_Post1Total_Clean,
                              Bag9_Post1Total_Clean, Bag10_Post1Total_Clean, Bag12_Post1Total_Clean,
                              Bag13_Post1Total_Clean, Bag16_Post1Total_Clean, Bag18_Post1Total_Clean,
                              Bag20_Post1Total_Clean)

Post1Total_Controls <- list(Bag1_Post1Total_Clean, Bag3_Post1Total_Clean, Bag6_Post1Total_Clean,
                            Bag8_Post1Total_Clean, Bag11_Post1Total_Clean, Bag14_Post1Total_Clean,
                            Bag15_Post1Total_Clean, Bag17_Post1Total_Clean, Bag19_Post1Total_Clean)

Post1_Treatment_Region <- matrix(nrow = nrow(Bag2_Post1Total_Clean), ncol = length(Post1Total_Treatments))
for (i in 1:length(Post1Total_Treatments)){
  Post1_Treatment_Region[,i] <- rowSums(Post1Total_Treatments[[i]])
}
rownames(Post1_Treatment_Region) <- rownames(Bag2_Post1Total_Clean)
colnames(Post1_Treatment_Region) <- c("Bag2", "Bag4", "Bag7", "Bag9", "Bag10",
                                      "Bag12", "Bag13", "Bag16", "Bag18", "Bag20")

Post1_Controls_Region <- matrix(nrow = nrow(Bag1_Post1Total_Clean), ncol = length(Post1Total_Controls))
for (i in 1:length(Post1Total_Controls)){
  Post1_Controls_Region[,i] <- rowSums(Post1Total_Controls[[i]])
}
rownames(Post1_Controls_Region) <- rownames(Bag1_Post1Total_Clean)
colnames(Post1_Controls_Region) <- c("Bag1", "Bag3", "Bag6", "Bag8", "Bag11",
                                     "Bag14", "Bag15", "Bag17", "Bag19")

### Post2 Totals ----

#For each site, sum across three days of second post-treatment sampling to get total for Post2 sampling period

Bag1_Post2Total_Clean <- Bag1_Post4_Clean +
  Bag1_Post5_Clean + Bag1_Post6_Clean

#Remove extra plants from Bag1
Bag1_Post2Total_Clean <- select(Bag1_Post2Total_Clean, -c(3,9,17))

Bag2_Post2Total_Clean <- Bag2_Post4_Clean +
  Bag2_Post5_Clean + Bag2_Post6_Clean

Bag3_Post2Total_Clean <- Bag3_Post4_Clean +
  Bag3_Post5_Clean + Bag3_Post6_Clean

Bag4_Post2Total_Clean <- Bag4_Post4_Clean +
  Bag4_Post5_Clean + Bag4_Post6_Clean

Bag6_Post2Total_Clean <- Bag6_Post4_Clean +
  Bag6_Post5_Clean + Bag6_Post6_Clean

Bag7_Post2Total_Clean <- Bag7_Post4_Clean +
  Bag7_Post5_Clean + Bag7_Post6_Clean

Bag8_Post2Total_Clean <- Bag8_Post4_Clean +
  Bag8_Post5_Clean + Bag8_Post6_Clean

Bag9_Post2Total_Clean <- Bag9_Post4_Clean +
  Bag9_Post5_Clean + Bag9_Post6_Clean

Bag10_Post2Total_Clean <- Bag10_Post4_Clean +
  Bag10_Post5_Clean + Bag10_Post6_Clean

Bag11_Post2Total_Clean <- Bag11_Post4_Clean +
  Bag11_Post5_Clean + Bag11_Post6_Clean

Bag12_Post2Total_Clean <- Bag12_Post4_Clean +
  Bag12_Post5_Clean + Bag12_Post6_Clean

Bag13_Post2Total_Clean <- Bag13_Post4_Clean +
  Bag13_Post5_Clean + Bag13_Post6_Clean

Bag14_Post2Total_Clean <- Bag14_Post4_Clean +
  Bag14_Post5_Clean + Bag14_Post6_Clean

Bag15_Post2Total_Clean <- Bag15_Post4_Clean +
  Bag15_Post5_Clean + Bag15_Post6_Clean

Bag16_Post2Total_Clean <- Bag16_Post4_Clean +
  Bag16_Post5_Clean + Bag16_Post6_Clean

Bag17_Post2Total_Clean <- Bag17_Post4_Clean +
  Bag17_Post5_Clean + Bag17_Post6_Clean

#Remove extra plant from Bag17
Bag17_Post2Total_Clean <- select(Bag17_Post2Total_Clean, -2)

Bag18_Post2Total_Clean <- Bag18_Post4_Clean +
  Bag18_Post5_Clean + Bag18_Post6_Clean

Bag19_Post2Total_Clean <- Bag19_Post4_Clean +
  Bag19_Post5_Clean + Bag19_Post6_Clean

Bag20_Post2Total_Clean <- Bag20_Post4_Clean +
  Bag20_Post5_Clean + Bag20_Post6_Clean

Post2Total_Treatments <- list(Bag2_Post2Total_Clean, Bag4_Post2Total_Clean, Bag7_Post2Total_Clean,
                              Bag9_Post2Total_Clean, Bag10_Post2Total_Clean, Bag12_Post2Total_Clean,
                              Bag13_Post2Total_Clean, Bag16_Post2Total_Clean, Bag18_Post2Total_Clean,
                              Bag20_Post2Total_Clean)

Post2Total_Controls <- list(Bag1_Post2Total_Clean, Bag3_Post2Total_Clean, Bag6_Post2Total_Clean,
                            Bag8_Post2Total_Clean, Bag11_Post2Total_Clean, Bag14_Post2Total_Clean,
                            Bag15_Post2Total_Clean, Bag17_Post2Total_Clean, Bag19_Post2Total_Clean)

Post2_Treatment_Region <- matrix(nrow = nrow(Bag2_Post2Total_Clean), ncol = length(Post2Total_Treatments))
for (i in 1:length(Post2Total_Treatments)){
  Post2_Treatment_Region[,i] <- rowSums(Post2Total_Treatments[[i]])
}
rownames(Post2_Treatment_Region) <- rownames(Bag2_Post2Total_Clean)
colnames(Post2_Treatment_Region) <- c("Bag2", "Bag4", "Bag7", "Bag9", "Bag10",
                                      "Bag12", "Bag13", "Bag16", "Bag18", "Bag20")

Post2_Controls_Region <- matrix(nrow = nrow(Bag1_Post2Total_Clean), ncol = length(Post2Total_Controls))
for (i in 1:length(Post2Total_Controls)){
  Post2_Controls_Region[,i] <- rowSums(Post2Total_Controls[[i]])
}
rownames(Post2_Controls_Region) <- rownames(Bag1_Post2Total_Clean)
colnames(Post2_Controls_Region) <- c("Bag1", "Bag3", "Bag6", "Bag8", "Bag11",
                                     "Bag14", "Bag15", "Bag17", "Bag19")



### Remove colonists----

#The null model cannot account for species that colonized the experience after the pre-disturbance survey
#As a result, analyzing the data without these species may be important

#Remove colonist species for Post2 Treatment sites

Post2Total_Treatments_NoCol <- Post2Total_Treatments #a list where each object is a spp x site matrix
Post2_Treatment_Region_NoCol <- matrix(nrow = nrow(Bag2_Post1Total_Clean), ncol = length(Post1Total_Treatments)) # a single spp x site matrix where each site is a treatment site
colnames(Post2_Treatment_Region_NoCol) <- c("Bag2", "Bag4", "Bag7", "Bag9", "Bag10", 
                                            "Bag12", "Bag13", "Bag16", "Bag18", "Bag20")
rownames(Post2_Treatment_Region_NoCol) <- rownames(Post2_Treatment_Region)

for (i in 1:length(PreTotal_Treatments)){ #for each site
  site <- Post2Total_Treatments[[i]] #assign the current site to an object
  for (j in 1:nrow(PreTotal_Treatments[[i]])){ #for each species in a site
    if (sum(Post2Total_Treatments[[i]][j,]) > 0){ #if that species was observed in Post2
      if (sum(PreTotal_Treatments[[i]][j,]) == 0){ #but that species was not observed in Pre-Survey
        site[j,] <- 0 #Remove species from that site in Post2 survey
      }
    }
  }
  Post2Total_Treatments_NoCol[[i]] <- site #Store site without colonists into new list
  Post2_Treatment_Region_NoCol[,i] <- rowSums(site) #Sum across all plants in site to get site total and add to region
}

#Remove colonist species for Post2 Control sites

Post2Total_Controls_NoCol <- Post2Total_Controls

for (i in 1:length(PreTotal_Controls)){
  site <- Post2Total_Controls[[i]]
  for (j in 1:nrow(PreTotal_Controls[[i]])){
    if (sum(Post2Total_Controls[[i]][j,]) > 0){
      if (sum(PreTotal_Controls[[i]][j,]) == 0){
        site[j,] <- 0
      }
    }
  }
  Post2Total_Controls_NoCol[[i]] <- site
}

#Remove colonist species for Post1 Treatment sites

Post1Total_Treatments_NoCol <- Post1Total_Treatments
Post1_Treatment_Region_NoCol <- matrix(nrow = nrow(Bag2_Post1Total_Clean), ncol = length(Post1Total_Treatments))
colnames(Post1_Treatment_Region_NoCol) <- c("Bag2", "Bag4", "Bag7", "Bag9", "Bag10", 
                                            "Bag12", "Bag13", "Bag16", "Bag18", "Bag20")
rownames(Post1_Treatment_Region_NoCol) <- rownames(Post1_Treatment_Region)

for (i in 1:length(PreTotal_Treatments)){
  site <- Post1Total_Treatments[[i]]
  for (j in 1:nrow(PreTotal_Treatments[[i]])){
    if (sum(Post1Total_Treatments[[i]][j,]) > 0){
      if (sum(PreTotal_Treatments[[i]][j,]) == 0){
        site[j,] <- 0
      }
    }
  }
  Post1Total_Treatments_NoCol[[i]] <- site
  Post1_Treatment_Region_NoCol[,i] <- rowSums(site)
}

#Remove colonist species for Post1 Control sites

Post1Total_Controls_NoCol <- Post1Total_Controls

for (i in 1:length(PreTotal_Controls)){
  site <- Post1Total_Controls[[i]]
  for (j in 1:nrow(PreTotal_Controls[[i]])){
    if (sum(Post1Total_Controls[[i]][j,]) > 0){
      if (sum(PreTotal_Controls[[i]][j,]) == 0){
        site[j,] <- 0
      }
    }
  }
  Post1Total_Controls_NoCol[[i]] <- site
}

### Post 3 Totals ----
Bag1_Post3Total <- Bag1_Post7 +
  Bag1_Post8 + Bag1_Post9

Bag2_Post3Total <- Bag2_Post7 +
  Bag2_Post8 + Bag2_Post9

Bag3_Post3Total <- Bag3_Post7 +
  Bag3_Post8 + Bag3_Post9

Bag4_Post3Total <- Bag4_Post7 +
  Bag4_Post8 + Bag4_Post9

Bag6_Post3Total <- Bag6_Post7 +
  Bag6_Post8 + Bag6_Post9

Bag7_Post3Total <- Bag7_Post7 +
  Bag7_Post8 + Bag7_Post9

Bag8_Post3Total <- Bag8_Post7 +
  Bag8_Post8 + Bag8_Post9

Bag9_Post3Total <- Bag9_Post7 +
  Bag9_Post8 + Bag9_Post9

Bag10_Post3Total <- Bag10_Post7 +
  Bag10_Post8 + Bag10_Post9

Bag11_Post3Total <- Bag11_Post7 +
  Bag11_Post8 + Bag11_Post9

Bag12_Post3Total <- Bag12_Post7 +
  Bag12_Post8 + Bag12_Post9

Bag13_Post3Total <- Bag13_Post7 +
  Bag13_Post8 + Bag13_Post9

Bag14_Post3Total <- Bag14_Post7 +
  Bag14_Post8 + Bag14_Post9

Bag15_Post3Total <- Bag15_Post7 +
  Bag15_Post8 + Bag15_Post9

Bag16_Post3Total <- Bag16_Post7 +
  Bag16_Post8 + Bag16_Post9

Bag17_Post3Total <- Bag17_Post7 +
  Bag17_Post8 + Bag17_Post9

Bag18_Post3Total <- Bag18_Post7 +
  Bag18_Post8 + Bag18_Post9

Bag19_Post3Total <- Bag19_Post7 +
  Bag19_Post8 + Bag19_Post9

Bag20_Post3Total <- Bag20_Post7 +
  Bag20_Post8 + Bag20_Post9

Post3Total_Treatments <- list(Bag2_Post3Total, Bag4_Post3Total, Bag7_Post3Total,
                              Bag9_Post3Total, Bag10_Post3Total, Bag12_Post3Total,
                              Bag13_Post3Total, Bag16_Post3Total, Bag18_Post3Total,
                              Bag20_Post3Total)

for (i in 1:length(Post3Total_Treatments)){
  site <- Post3Total_Treatments[[i]]
  rownames(site) <- t(spp_listPost3[,1])
  Post3Total_Treatments[[i]] <- site
}

Post3Total_Controls <- list(Bag1_Post3Total, Bag3_Post3Total, Bag6_Post3Total,
                            Bag8_Post3Total, Bag11_Post3Total, Bag14_Post3Total,
                            Bag15_Post3Total, Bag17_Post3Total, Bag19_Post3Total)

for (i in 1:length(Post3Total_Controls)){
  site <- Post3Total_Controls[[i]]
  rownames(site) <- t(spp_listPost3[,1])
  Post3Total_Controls[[i]] <- site
}

Post3_Treatment_Region <- matrix(nrow = nrow(Bag2_Post3Total), ncol = length(Post3Total_Treatments))
for (i in 1:length(Post3Total_Treatments)){
  Post3_Treatment_Region[,i] <- rowSums(Post3Total_Treatments[[i]])
}
rownames(Post3_Treatment_Region) <- rownames(Bag2_Post3Total)
colnames(Post3_Treatment_Region) <- c("Bag2", "Bag4", "Bag7", "Bag9", "Bag10",
                                      "Bag12", "Bag13", "Bag16", "Bag18", "Bag20")

Post3_Controls_Region <- matrix(nrow = nrow(Bag1_Post3Total), ncol = length(Post3Total_Controls))
for (i in 1:length(Post3Total_Controls)){
  Post3_Controls_Region[,i] <- rowSums(Post3Total_Controls[[i]])
}
rownames(Post3_Controls_Region) <- rownames(Bag1_Post3Total)
colnames(Post3_Controls_Region) <- c("Bag1", "Bag3", "Bag6", "Bag8", "Bag11",
                                     "Bag14", "Bag15", "Bag17", "Bag19")


### Define variables----

#number of species in experiment
nSpp <- nrow(Bag1_Post1_Clean)

#number of treatment sites
nTreatments <- length(Post1Total_Treatments)

#number of control sites
nControls <- length(Post1Total_Controls)

#Treatment site names
treatmentSiteNames <- c("Bag2", "Bag4", "Bag7", "Bag9", "Bag10",
                        "Bag12", "Bag13", "Bag16", "Bag18", "Bag20")

#Control site names
controlSiteNames <- c("Bag1", "Bag3", "Bag6", "Bag11", "Bag10",
                      "Bag14", "Bag15", "Bag17", "Bag19")

#Cleaned species list
sppNames <- rownames(Bag1_Post1_Clean)
sppNamesPost3 <- rownames(Bag1_Post3Total)

#number of plants to rarefy down to
nRarefy <- 8




#################### NULL MODEL FUNCTIONS ####################

#'Conduct a spatially clustered rarefaction
#'
#'@param comm a site by species matrix 
#'@param samplesize an integer representing subsample size
#'@return a matrix made up of 1000 rarefied site by species matrices
#'
#'@description a matrix with 10 sites and 50 spp rarefied to sample size 5
#' will return a matrix with 10 rows and 50,000 columns, 
#' where each interval of 50 columns has sequences of 5 rows randomly changed to zeros 
spatially_clustered_rarefaction <- function(comm, samplesize){
  
  comm_samples <- 1:nrow(comm)
  comm_SBR_Full <- matrix(nrow = nrow(comm), ncol = 1)
  comm_SBR_Sim <- comm
  
  
  for (i in 1:1000){
    x <- rep(seq(1:nrow(comm)), times = 2)
    start <- sample(x, size = 1)
    comm_post_samples <- x[seq(start, start+7)]
    for (j in comm_post_samples){
      comm_SBR_Sim[j,] <- 0}
    comm_SBR_Full <- cbind(comm_SBR_Full, comm_SBR_Sim)
    comm_SBR_Sim <- comm
  }
  comm_SBR_Full <- comm_SBR_Full[,-1]
  return(comm_SBR_Full)
}

# Note that mechanistically, the spatially clustered rarefaction works almost identically to
# sample based rarefaction with the same output data structure; 
# the difference is in this case, our "sample" is just a sequence of plants. 
# For this reason, the following functions are written for sample-based rarefaction

#'Determine expected local richness based on sample based rarefaction
#'
#'@param SBR_comm a rarefied community matrix,
#'i.e. the output of spatially_clustered_rarefaction function
#'@return a matrix containing 1000 rarefied richness values
#'
#'@discription mean of output matrix is expected richness following
#'rarefaction
SBR_richness <- function(SBR_comm){
  
  SBR_comm_original <- SBR_comm
  SBR_comm_SBR_richnesses <- matrix(nrow = 1, ncol = 1000)
  SBR_sp_counter <- 0
  
  for (i in 1:1000){
    SBR_comm_SBR_Sim <- SBR_comm[,1:(ncol(SBR_comm_original)/1000)]
    for (j in 1:ncol(SBR_comm_SBR_Sim)){
      if (sum(SBR_comm_SBR_Sim[,j] > 0)){
        SBR_sp_counter <- SBR_sp_counter + 1}}
    SBR_comm_SBR_richnesses[,i] <- SBR_sp_counter
    SBR_sp_counter <- 0
    SBR_comm <- SBR_comm[,-(1:(ncol(SBR_comm_original)/1000))]
  }
  return(SBR_comm_SBR_richnesses)
}


#'Determine expected regional richness based on sample based rarefaction
#'
#'@param SBR_region a list containing the rarefied community matrix for
#'each site in a region (i.e. the output of spatially_clustered_rarefaction for
#'each site)
#'@return a matrix containing 1000 rarefied richness values
#'
#'@discription mean of output matrix is expected richness following
#'rarefaction
SBR_Regional_richnesses <- function(SBR_region){
  sp_counter <- 0
  RegSp <- 0
  RegionalSBR_Test_Matrix <- matrix(nrow = 1, ncol = 1000)
  
  for (j in 1:1000){
    for (i in 1:(ncol(SBR_region[[1]])/1000)){
      for (x in SBR_region){
        if (sum(x[,i + ((j-1)*(ncol(x)/1000))]) >0 ){
          sp_counter <- sp_counter + 1}}
      if (sp_counter > 0){
        RegSp <- RegSp + 1
        sp_counter <- 0}
      sp_counter <- 0}
    RegionalSBR_Test_Matrix[,j] <- RegSp
    RegSp <- 0}
  
  return(RegionalSBR_Test_Matrix)
}

#' Determine local Pe values based on sample based rarefaction
#' 
#' @param SBR_comm a rarefied community matrix,
#'the output of spatially_clustered_rarefaction function
#'@return a matrix containing local Pe value for each species
#'
#'@description Pe values are the proportion of times a species was
#'simulated to go extinct across 1000 rarefaction simulations
SBR_LocalPe <- function(SBR_comm){
  SBR_locEx <- matrix(nrow = 1, ncol = (ncol(SBR_comm)/1000))
  ex_num  <- 0
  
  for (i in 1:ncol(SBR_locEx)){
    for (j in 1:1000){
      if (sum(SBR_comm[,(i+((j-1)*(ncol(SBR_comm)/1000)))])==0){
        ex_num <- ex_num + 1}}
    SBR_locEx[1,i] <- ex_num/1000
    ex_num <- 0
    j <- 0}
  
  colnames(SBR_locEx) <- colnames(SBR_comm[,1:(ncol(SBR_comm)/1000)])
  return(SBR_locEx)
}

#'Determine regional Pe values based on sample based rarefaction
#'
#'@param SBR_region a list containing the rarefied community matrices for
#'each site in a region (i.e. the output of spatially_clustered_rarefaction for
#'each site)
#'@return a matrix containing regional Pe value for each species
#'
#'@discription Pe values are the proportion of times a species was
#'simulated to go extinct across 1000 rarefaction simulations
SBR_RegionalPe <- function(region){
  local_ex_counter <- 0
  regional_ex_counter<- 0
  RegionalPe_Matrix <- matrix(nrow = 1, ncol = ncol(region[[1]])/1000)
  
  for (i in 1:(ncol(region[[1]])/1000)){
    for (j in 1:1000){
      for (x in region){
        if (sum(x[,i + ((j-1)*(ncol(x)/1000))]) == 0 ){
          local_ex_counter <- local_ex_counter + 1}}
      if (local_ex_counter == length(region)){
        regional_ex_counter <- regional_ex_counter + 1
        local_ex_counter <- 0}
      local_ex_counter <- 0}
    RegionalPe_Matrix[,i] <- regional_ex_counter/1000
    regional_ex_counter <- 0}
  
  return(RegionalPe_Matrix)
}

#'Species-area relationship
#'
#'@param comm a species by site matrix
#'@param aNew the percentage of habitat remaining post-disturbance
#'@return postS, a prediction of post-disturbance species richness
#'
#'@description calculate species richness following a loss of area following
#'Brooks et al., 2002
SAR <- function(comm, aNew){
  preS <- specnumber(comm)
  postS <- preS*((1-aNew)^0.25)
  return(postS)
}
#################### RAREFACTION SIMULATIONS ##########

#Create list where each element will be a rarefied site
Treatment_Rarefactions_Region <- vector("list", nTreatments)

#Conduct rarefaction simulation
for (i in 1:nTreatments){ #for each site
  Treatment_Rarefactions_Region[[i]] <- spatially_clustered_rarefaction(t(PreTotal_Treatments[[i]]), nRarefy) #rarefy each site to 50% of remaining plants
}


#################### SPECIES RICHNESS ANALYSES ##########

#' @description Calculate SAR-predicted, null-predicted, and observed richnesses at both local and regional
#' scale for Pre, Post1, and Post2 sampling periods. For treatment sites, do this both including and excluding
#' colonist species. Run code sequentially to get final output dataframe.
#' @output TotalRichnesses, a tidy dataframe where each observation represents a richness with
#' eight variables: site, richness, survey, colonists, treatment, null, scale, SAR

### Observed Local Species richness ----

#Calculate local pre-treatment richness

PreTreatment_Richnesses<- data.frame("site" = 1:nTreatments, "richness" = 1:nTreatments, "survey" = 1:nTreatments) #create a df where each observation is a site
PreTreatment_Richnesses$survey <- "Pre"
PreTreatment_Richnesses$colonists <- NA
PreTreatment_Richnesses$treatment <- "Treatment"
PreTreatment_Richnesses$null <- 0

for (i in 1:nTreatments){ #for each site
  PreTreatment_Richnesses$richness[i] <- specnumber(Pre_Treatment_Region[,i]) #calculate richness at each site
}
PreTreatment_Richnesses$site <- treatmentSiteNames #assign site names to rows

#Post 2 Treatment
Post2Treatment_Richnesses<- data.frame("site" = 1:nTreatments, "richness" = 1:nTreatments, "survey" = 1:nTreatments) #create a df where each observation is a site
Post2Treatment_Richnesses$survey <- "Post2"
Post2Treatment_Richnesses$colonists <- 1
Post2Treatment_Richnesses$treatment <- "Treatment"
Post2Treatment_Richnesses$null <- 0

for (i in 1:nTreatments){ #for each site
  Post2Treatment_Richnesses$richness[i] <- specnumber(Post2_Treatment_Region[,i]) #calculate richness at each site
}
Post2Treatment_Richnesses$site <- treatmentSiteNames #assign site names to rows


#Post 1 Treatment
Post1Treatment_Richnesses<- data.frame("site" = 1:nTreatments, "richness" = 1:nTreatments, "survey" = 1:nTreatments) #create a df where each observation is a site
Post1Treatment_Richnesses$survey <- "Post1"
Post1Treatment_Richnesses$colonists <- 1
Post1Treatment_Richnesses$treatment <- "Treatment"
Post1Treatment_Richnesses$null <- 0

for (i in 1:nTreatments){ #for each site
  Post1Treatment_Richnesses$richness[i] <- specnumber(Post1_Treatment_Region[,i]) #calculate richness at each site
}
Post1Treatment_Richnesses$site <- treatmentSiteNames #assign site names to rows


#Pre Control
PreControls_Richnesses<- data.frame("site" = 1:nControls, "richness" = 1:nControls, "survey" = 1:nControls) #create a df where each observation is a site
PreControls_Richnesses$survey <- "Pre"
PreControls_Richnesses$colonists <- NA
PreControls_Richnesses$treatment <- "Control"
PreControls_Richnesses$null <- 0

for (i in 1:nControls){ #for each site
  PreControls_Richnesses$richness[i] <- specnumber(Pre_Controls_Region[,i]) #calculate richness at each site
}
PreControls_Richnesses$site <- controlSiteNames #assign site names to rows


#Post1 Control
Post1Controls_Richnesses<- data.frame("site" = 1:nControls, "richness" = 1:nControls, "survey" = 1:nControls) #create a df where each observation is a site
Post1Controls_Richnesses$survey <- "Post1"
Post1Controls_Richnesses$colonists <- 1
Post1Controls_Richnesses$treatment <- "Control"
Post1Controls_Richnesses$null <- 0

for (i in 1:nControls){ #for each site
  Post1Controls_Richnesses$richness[i] <- specnumber(Post1_Controls_Region[,i]) #calculate richness at each site
}
Post1Controls_Richnesses$site <- controlSiteNames #assign site names to rows


#Post2 Control
Post2Controls_Richnesses<- data.frame("site" = 1:nControls, "richness" = 1:nControls, "survey" = 1:nControls) #create a df where each observation is a site
Post2Controls_Richnesses$survey <- "Post2"
Post2Controls_Richnesses$colonists <- 1
Post2Controls_Richnesses$treatment <- "Control"
Post2Controls_Richnesses$null <- 0

for (i in 1:nControls){ #for each site
  Post2Controls_Richnesses$richness[i] <- specnumber(Post2_Controls_Region[,i]) #calculate richness at each site
}
Post2Controls_Richnesses$site <- controlSiteNames #assign site names to rows

#Post1 Treatment No colonists
Post1Treatment_Richnesses_NoCol<- data.frame("site" = 1:nTreatments, "richness" = 1:nTreatments, "survey" = 1:nTreatments) #create a df where each observation is a site
Post1Treatment_Richnesses_NoCol$survey <- "Post1"
Post1Treatment_Richnesses_NoCol$colonists <- 0
Post1Treatment_Richnesses_NoCol$treatment <- "Treatment"
Post1Treatment_Richnesses_NoCol$null <- 0

for (i in 1:nTreatments){ #for each site
  Post1Treatment_Richnesses_NoCol$richness[i] <- specnumber(Post1_Treatment_Region_NoCol[,i]) #calculate richness at each site
}
Post1Treatment_Richnesses_NoCol$site <- treatmentSiteNames #assign site names to rows


#Post2 Treatment No colonists
Post2Treatment_Richnesses_NoCol<- data.frame("site" = 1:nTreatments, "richness" = 1:nTreatments, "survey" = 1:nTreatments) #create a df where each observation is a site
Post2Treatment_Richnesses_NoCol$survey <- "Post2"
Post2Treatment_Richnesses_NoCol$colonists <- 0
Post2Treatment_Richnesses_NoCol$treatment <- "Treatment"
Post2Treatment_Richnesses_NoCol$null <- 0

for (i in 1:nTreatments){ #for each site
  Post2Treatment_Richnesses_NoCol$richness[i] <- specnumber(Post2_Treatment_Region_NoCol[,i]) #calculate richness at each site
}
Post2Treatment_Richnesses_NoCol$site <- treatmentSiteNames #assign site names to rows


### Local Null_Expected Richnesses ----

#Null treatment richnesses Post1 (same as Post2)
NullTreatmentPost1_Richnesses <- data.frame("site" = 1:nTreatments, "richness" = 1:nTreatments, "survey" = 1:nTreatments)
NullTreatmentPost1_Richnesses$survey <- "Post1"
NullTreatmentPost1_Richnesses$colonists <- 0
NullTreatmentPost1_Richnesses$treatment <- "Treatment"
NullTreatmentPost1_Richnesses$null <- 1

for (i in 1:nTreatments){
  NullTreatmentPost1_Richnesses$richness[i] <- mean(SBR_richness(Treatment_Rarefactions_Region[[i]]))
}

NullTreatmentPost1_Richnesses$site <- treatmentSiteNames

#Null treatment richnesses Post2 (same as Post1)
NullTreatmentPost2_Richnesses <- NullTreatmentPost1_Richnesses
NullTreatmentPost2_Richnesses$survey <- "Post2"

#Add colonists Post2
NullTreatment_Richnesses_colonistsPost2 <- NullTreatmentPost2_Richnesses
NullTreatment_Richnesses_colonistsPost2$colonists <- 1

for (i in 1:nTreatments){
  colonists <- 0
  for (j in 1:nSpp){
    if (sum(Post2Total_Treatments[[i]][j,]) > 0){
      if (sum(PreTotal_Treatments[[i]][j,]) == 0){
        colonists <- colonists + 1
      }
    }
  }
  NullTreatment_Richnesses_colonistsPost2$richness[i] <- NullTreatmentPost2_Richnesses$richness[i] + colonists
}

#Add colonists Post1
NullTreatment_Richnesses_colonistsPost1 <- NullTreatmentPost1_Richnesses
NullTreatment_Richnesses_colonistsPost1$colonists <- 1

for (i in 1:nTreatments){
  colonists <- 0
  for (j in 1:nSpp){
    if (sum(Post1Total_Treatments[[i]][j,]) > 0){
      if (sum(PreTotal_Treatments[[i]][j,]) == 0){
        colonists <- colonists + 1
      }
    }
  }
  NullTreatment_Richnesses_colonistsPost1$richness[i] <- NullTreatmentPost1_Richnesses$richness[i] + colonists
}

### Combine local richnesses (Observed to Null) ----
TotalRichnesses <- rbind(PreTreatment_Richnesses, Post1Treatment_Richnesses, Post2Treatment_Richnesses,
                         PreControls_Richnesses, Post1Controls_Richnesses, Post2Controls_Richnesses,
                         Post1Treatment_Richnesses_NoCol, Post2Treatment_Richnesses_NoCol,
                         NullTreatmentPost1_Richnesses, NullTreatmentPost2_Richnesses, 
                         NullTreatment_Richnesses_colonistsPost1, NullTreatment_Richnesses_colonistsPost2)

#Since we need to add regional richness soon, let's label all of these as local 
TotalRichnesses$scale <- "Local"
#Next we'll do species-area predictions. Let's make a column to denote whether or not these are SAR predictions.
TotalRichnesses$SAR <- 0

### Add local SAR richnesses to total richnesses----

#Add local SAR Post1(same as post 2)
SARTreatmentPost1_Richnesses <- data.frame("site" = 1:nTreatments, "richness" = 1:nTreatments, "survey" = 1:nTreatments)
SARTreatmentPost1_Richnesses$survey <- "Post1"
SARTreatmentPost1_Richnesses$colonists <- 0
SARTreatmentPost1_Richnesses$treatment <- "Treatment"
SARTreatmentPost1_Richnesses$null <- 0
SARTreatmentPost1_Richnesses$scale <- "Local"
SARTreatmentPost1_Richnesses$SAR <- 1

for (i in 1:nTreatments){
  SARTreatmentPost1_Richnesses$richness[i] <- SAR(Pre_Treatment_Region[,i], 0.5)
}

SARTreatmentPost1_Richnesses$site <- treatmentSiteNames
TotalRichnesses <- rbind(TotalRichnesses, SARTreatmentPost1_Richnesses)

#Add local SAR Post 2 (same as Post1)
SARTreatmentPost2_Richnesses <- SARTreatmentPost1_Richnesses
SARTreatmentPost2_Richnesses$survey <- "Post2"
TotalRichnesses <- rbind(TotalRichnesses, SARTreatmentPost2_Richnesses)

#Add local SAR Post1 + Colonists
SARTreatmentPost1_colonists_Richnesses <- SARTreatmentPost1_Richnesses
SARTreatmentPost1_colonists_Richnesses$colonists <- 1

for (i in 1:nTreatments){
  colonists <- 0
  for (j in 1:nSpp){
    if (sum(Post1Total_Treatments[[i]][j,]) > 0){
      if (sum(PreTotal_Treatments[[i]][j,]) == 0){
        colonists <- colonists + 1
      }
    }
  }
  SARTreatmentPost1_colonists_Richnesses$richness[i] <- SARTreatmentPost1_colonists_Richnesses$richness[i] + colonists
}
TotalRichnesses <- rbind(TotalRichnesses, SARTreatmentPost1_colonists_Richnesses)

#Add local SAR Post2 + Colonists
SARTreatmentPost2_colonists_Richnesses <- SARTreatmentPost2_Richnesses
SARTreatmentPost2_colonists_Richnesses$colonists <- 1

for (i in 1:nTreatments){
  colonists <- 0
  for (j in 1:nSpp){
    if (sum(Post2Total_Treatments[[i]][j,]) > 0){
      if (sum(PreTotal_Treatments[[i]][j,]) == 0){
        colonists <- colonists + 1
      }
    }
  }
  SARTreatmentPost2_colonists_Richnesses$richness[i] <- SARTreatmentPost2_colonists_Richnesses$richness[i] + colonists
}
TotalRichnesses <- rbind(TotalRichnesses, SARTreatmentPost2_colonists_Richnesses)


### Regional Species Richness Treatments----

#Pre
PreRegional_Abund <- vector(length = nSpp)
PreRegionalRichness <- data.frame("site" = "Region", "richness" = 0, "survey" = "Pre", 
                                  "colonists" = NA, "treatment" = "Treatment", 
                                  "null" = 0, "scale" = "Regional", "SAR" = 0)
for (site in PreTotal_Treatments){
  PreRegional_Abund <- PreRegional_Abund + rowSums(site)
}

PreRegionalRichness$richness <- specnumber(PreRegional_Abund)
TotalRichnesses <- rbind(TotalRichnesses, PreRegionalRichness)

#Post1
Post1Regional_Abund <- vector(length = nSpp)
Post1RegionalRichness <- data.frame("site" = "Region", "richness" = 0, "survey" = "Post1", 
                                    "colonists" = 1, "treatment" = "Treatment", 
                                    "null" = 0, "scale" = "Regional", "SAR" = 0)
for (site in Post1Total_Treatments){
  Post1Regional_Abund <- Post1Regional_Abund + rowSums(site)
}

Post1RegionalRichness$richness <- specnumber(Post1Regional_Abund)
TotalRichnesses <- rbind(TotalRichnesses, Post1RegionalRichness)

#Post2
Post2Regional_Abund <- vector(length = nSpp)
Post2RegionalRichness <- data.frame("site" = "Region", "richness" = 0, "survey" = "Post2", 
                                    "colonists" = 1, "treatment" = "Treatment", 
                                    "null" = 0, "scale" = "Regional", "SAR" = 0)
for (site in Post2Total_Treatments){
  Post2Regional_Abund <- Post2Regional_Abund + rowSums(site)
}

Post2RegionalRichness$richness <- specnumber(Post2Regional_Abund)
TotalRichnesses <- rbind(TotalRichnesses, Post2RegionalRichness)

#Null Post1
NullPost1Regional_Abund <- vector(length = nSpp)
NullPost1RegionalRichness <- data.frame("site" = "Region", "richness" = 0, "survey" = "Post1", 
                                        "colonists" = 0, "treatment" = "Treatment", 
                                        "null" = 1, "scale" = "Regional", "SAR" = 0)

NullPost1RegionalRichness$richness <- mean(SBR_Regional_richnesses(Treatment_Rarefactions_Region))

TotalRichnesses <- rbind(TotalRichnesses, NullPost1RegionalRichness)

#Null Post2
NullPost2Regional_Abund <- vector(length = nSpp)
NullPost2RegionalRichness <- data.frame("site" = "Region", "richness" = 0, "survey" = "Post2", 
                                        "colonists" = 0, "treatment" = "Treatment", 
                                        "null" = 1, "scale" = "Regional", "SAR" = 0)

NullPost2RegionalRichness$richness <- mean(SBR_Regional_richnesses(Treatment_Rarefactions_Region))

TotalRichnesses <- rbind(TotalRichnesses, NullPost2RegionalRichness)

#Null + Colonists Post1

NullPost1RegionalRichness_colonists <- NullPost1RegionalRichness
NullPost1RegionalRichness_colonists$colonists <- 1

Post1_Regional_colonists <- 0
for (i in 1:length(PreRegional_Abund)){
  if (Post1Regional_Abund[i] > 0){
    if (PreRegional_Abund[i] == 0){
      Post1_Regional_colonists <- Post1_Regional_colonists + 1
    }
  }
}

NullPost1RegionalRichness_colonists$richness <- NullPost1RegionalRichness$richness + colonists
TotalRichnesses <- rbind(TotalRichnesses, NullPost1RegionalRichness_colonists)

#Null + Colonists Post2

NullPost2RegionalRichness_colonists <- NullPost2RegionalRichness
NullPost2RegionalRichness_colonists$colonists <- 1

Post2_Regional_colonists <- 0
for (i in 1:length(PreRegional_Abund)){
  if (Post2Regional_Abund[i] > 0){
    if (PreRegional_Abund[i] == 0){
      Post2_Regional_colonists <- Post2_Regional_colonists + 1
    }
  }
}

NullPost2RegionalRichness_colonists$richness <- NullPost2RegionalRichness$richness + colonists
TotalRichnesses <- rbind(TotalRichnesses, NullPost2RegionalRichness_colonists)

#Post1 No Colonists
Post1RegionalRichness_NoCol <- Post1RegionalRichness
Post1RegionalRichness_NoCol$colonists <- 0
Post1RegionalRichness_NoCol$richness <- Post1RegionalRichness$richness - Post1_Regional_colonists
TotalRichnesses <- rbind(TotalRichnesses, Post1RegionalRichness_NoCol)


#Post2 No colonists
Post2RegionalRichness_NoCol <- Post2RegionalRichness
Post2RegionalRichness_NoCol$colonists <- 0
Post2RegionalRichness_NoCol$richness <- Post2RegionalRichness$richness - Post2_Regional_colonists
TotalRichnesses <- rbind(TotalRichnesses, Post2RegionalRichness_NoCol)

### Regional SAR richness----

#Post1 Regional SAR
SARPost1RegionalRichness <- data.frame("site" = "Region", "richness" = 0, "survey" = "Post1", 
                                       "colonists" = 0, "treatment" = "Treatment", 
                                       "null" = 0, "scale" = "Regional", "SAR" = 1)

SARPost1RegionalRichness$richness <- SAR(PreRegional_Abund, 0.5)
TotalRichnesses <- rbind(TotalRichnesses, SARPost1RegionalRichness)

#Post2 Regional SAR
SARPost2RegionalRichness <- data.frame("site" = "Region", "richness" = 0, "survey" = "Post2", 
                                       "colonists" = 0, "treatment" = "Treatment", 
                                       "null" = 0, "scale" = "Regional", "SAR" = 1)

SARPost2RegionalRichness$richness <- SAR(PreRegional_Abund, 0.5)
TotalRichnesses <- rbind(TotalRichnesses, SARPost2RegionalRichness)

#Post1 Regional SAR + Colonists
SARPost1RegionalRichness_colonists <- SARPost1RegionalRichness
SARPost1RegionalRichness_colonists$colonists <- 1
SARPost1RegionalRichness_colonists$richness <- SARPost1RegionalRichness_colonists$richness + Post1_Regional_colonists
TotalRichnesses <- rbind(TotalRichnesses, SARPost1RegionalRichness_colonists)

#Post2 Regional SAR + Colonists
SARPost2RegionalRichness_colonists <- SARPost2RegionalRichness
SARPost2RegionalRichness_colonists$colonists <- 1
SARPost2RegionalRichness_colonists$richness <- SARPost2RegionalRichness_colonists$richness + Post2_Regional_colonists
TotalRichnesses <- rbind(TotalRichnesses, SARPost2RegionalRichness_colonists)


### Regional Species Richness Controls----

#Pre
PreRegionalControls_Abund <- vector(length = nSpp)
PreRegionalControlsRichness <- data.frame("site" = "Region", "richness" = 0, "survey" = "Pre", 
                                          "colonists" = NA, "treatment" = "Control", 
                                          "null" = 0, "scale" = "Regional", "SAR" = 0)

for (site in PreTotal_Controls){
  PreRegionalControls_Abund <- PreRegionalControls_Abund + rowSums(site)
}

PreRegionalControlsRichness$richness <- specnumber(PreRegionalControls_Abund)
TotalRichnesses <- rbind(TotalRichnesses, PreRegionalControlsRichness)

#Post1
Post1RegionalControls_Abund <- vector(length = nSpp)
Post1RegionalControlsRichness <- data.frame("site" = "Region", "richness" = 0, "survey" = "Post1", 
                                            "colonists" = 1, "treatment" = "Control", 
                                            "null" = 0, "scale" = "Regional", "SAR" = 0)

for (site in Post1Total_Controls){
  Post1RegionalControls_Abund <- Post1RegionalControls_Abund + rowSums(site)
}

Post1RegionalControlsRichness$richness <- specnumber(Post1RegionalControls_Abund)
TotalRichnesses <- rbind(TotalRichnesses, Post1RegionalControlsRichness)

#Post2
Post2RegionalControls_Abund <- vector(length = nSpp)
Post2RegionalControlsRichness <- data.frame("site" = "Region", "richness" = 0, "survey" = "Post2", 
                                            "colonists" = 1, "treatment" = "Control", 
                                            "null" = 0, "scale" = "Regional", "SAR" = 0)

for (site in Post2Total_Controls){
  Post2RegionalControls_Abund <- Post2RegionalControls_Abund + rowSums(site)
}

Post2RegionalControlsRichness$richness <- specnumber(Post2RegionalControls_Abund)
TotalRichnesses <- rbind(TotalRichnesses, Post2RegionalControlsRichness)


#################### EXTINCTION ANALYSIS ##########

#' @description Calculate SAR-predicted, null-predicted, and observed extinctions at both local and regional
#' scale for all Pre, Post1, and Post2 sampling periods. Run code sequentially to get final output dataframe.
#' @output TotalExtinctions, a tidy dataframe where each observation represents extinctions with
#' eight variables: site, richness, survey, colonists, treatment, null, scale, SAR

### Local Extinctions----
#Filter data to get richness for each group without colonists
PreLocalS <- filter(TotalRichnesses, survey == "Pre", treatment == "Treatment", scale == "Local")
Post1LocalS <- filter(TotalRichnesses, survey == "Post1", treatment == "Treatment", scale == "Local", 
                      null == 0, colonists == 0, SAR == 0)
Post2LocalS <- filter(TotalRichnesses, survey == "Post1", treatment == "Treatment", scale == "Local", 
                      null == 0, colonists == 0, SAR == 0)
NullPost1LocalS <- filter(TotalRichnesses, survey == "Post1", treatment == "Treatment", scale == "Local", 
                          null == 1, colonists == 0, SAR == 0)
NullPost2LocalS <- filter(TotalRichnesses, survey == "Post2", treatment == "Treatment", scale == "Local", 
                          null == 1, colonists == 0, SAR == 0)
SARPost1LocalS <- filter(TotalRichnesses, survey == "Post1", treatment == "Treatment", scale == "Local", 
                         null == 0, colonists == 0, SAR == 1)
SARPost2LocalS <- filter(TotalRichnesses, survey == "Post2", treatment == "Treatment", scale == "Local", 
                         null == 0, colonists == 0, SAR == 1)

#Calculate Extinctions
#Post1
Post1LocalExtinctions <- Post1LocalS
Post1LocalExtinctions <- rename(Post1LocalExtinctions, extinctions = richness)
Post1LocalExtinctions$extinctions <- PreLocalS$richness - Post1LocalS$richness

#Post2
Post2LocalExtinctions <- Post2LocalS
Post2LocalExtinctions <- rename(Post2LocalExtinctions, extinctions = richness)
Post2LocalExtinctions$extinctions <- PreLocalS$richness - Post2LocalS$richness

#NullPost1
NullPost1LocalExtinctions <- NullPost1LocalS
NullPost1LocalExtinctions <- rename(NullPost1LocalExtinctions, extinctions = richness)
NullPost1LocalExtinctions$extinctions <- PreLocalS$richness - NullPost1LocalS$richness

#NullPost2
NullPost2LocalExtinctions <- NullPost2LocalS
NullPost2LocalExtinctions <- rename(NullPost2LocalExtinctions, extinctions = richness)
NullPost2LocalExtinctions$extinctions <- PreLocalS$richness - NullPost2LocalS$richness

#SARPost1
SARPost1LocalExtinctions <- SARPost1LocalS
SARPost1LocalExtinctions <- rename(SARPost1LocalExtinctions, extinctions = richness)
SARPost1LocalExtinctions$extinctions <- PreLocalS$richness - SARPost1LocalS$richness

#SARPost2
SARPost2LocalExtinctions <- SARPost2LocalS
SARPost2LocalExtinctions <- rename(SARPost2LocalExtinctions, extinctions = richness)
SARPost2LocalExtinctions$extinctions <- PreLocalS$richness - SARPost2LocalS$richness

LocalExtinctions <- rbind(Post1LocalExtinctions, Post2LocalExtinctions, NullPost1LocalExtinctions,
                          NullPost2LocalExtinctions, SARPost1LocalExtinctions, SARPost2LocalExtinctions)

### Regional Extinctions ----
PreRegionalS <- filter(TotalRichnesses, survey == "Pre", treatment == "Treatment", scale == "Regional")
Post1RegionalS <- filter(TotalRichnesses, survey == "Post1", treatment == "Treatment", scale == "Regional", 
                         null == 0, colonists == 0, SAR == 0)
Post2RegionalS <- filter(TotalRichnesses, survey == "Post2", treatment == "Treatment", scale == "Regional", 
                         null == 0, colonists == 0, SAR == 0)
NullPost1RegionalS <- filter(TotalRichnesses, survey == "Post1", treatment == "Treatment", scale == "Regional", 
                             null == 1, colonists == 0, SAR == 0)
NullPost2RegionalS <- filter(TotalRichnesses, survey == "Post2", treatment == "Treatment", scale == "Regional", 
                             null == 1, colonists == 0, SAR == 0)
SARPost1RegionalS <- filter(TotalRichnesses, survey == "Post1", treatment == "Treatment", scale == "Regional", 
                            null == 0, colonists == 0, SAR == 1)
SARPost2RegionalS <- filter(TotalRichnesses, survey == "Post2", treatment == "Treatment", scale == "Regional", 
                            null == 0, colonists == 0, SAR == 1)

#Calculate Extinctions
#Post1
Post1RegionalExtinctions <- Post1RegionalS
Post1RegionalExtinctions <- rename(Post1RegionalExtinctions, extinctions = richness)
Post1RegionalExtinctions$extinctions <- PreRegionalS$richness - Post1RegionalS$richness

#Post2
Post2RegionalExtinctions <- Post2RegionalS
Post2RegionalExtinctions <- rename(Post2RegionalExtinctions, extinctions = richness)
Post2RegionalExtinctions$extinctions <- PreRegionalS$richness - Post2RegionalS$richness

#NullPost1
NullPost1RegionalExtinctions <- NullPost1RegionalS
NullPost1RegionalExtinctions <- rename(NullPost1RegionalExtinctions, extinctions = richness)
NullPost1RegionalExtinctions$extinctions <- PreRegionalS$richness - NullPost1RegionalS$richness

#NullPost2
NullPost2RegionalExtinctions <- NullPost2RegionalS
NullPost2RegionalExtinctions <- rename(NullPost2RegionalExtinctions, extinctions = richness)
NullPost2RegionalExtinctions$extinctions <- PreRegionalS$richness - NullPost2RegionalS$richness

#SARPost1
SARPost1RegionalExtinctions <- SARPost1RegionalS
SARPost1RegionalExtinctions <- rename(SARPost1RegionalExtinctions, extinctions = richness)
SARPost1RegionalExtinctions$extinctions <- PreRegionalS$richness - SARPost1RegionalS$richness

#SARPost2
SARPost2RegionalExtinctions <- SARPost2RegionalS
SARPost2RegionalExtinctions <- rename(SARPost2RegionalExtinctions, extinctions = richness)
SARPost2RegionalExtinctions$extinctions <- PreRegionalS$richness - SARPost2RegionalS$richness

RegionalExtinctions <- rbind(Post1RegionalExtinctions, Post2RegionalExtinctions, NullPost1RegionalExtinctions,
                             NullPost2RegionalExtinctions, SARPost1RegionalExtinctions, SARPost2RegionalExtinctions)

### Total Extinctions ----
TotalExtinctions <- rbind(LocalExtinctions, RegionalExtinctions)

#################### EXTINCTION PROBABILITIES ##########

#' @description Calculate local and regional simulated extinction probabilities (Pe)
#' @output TotalPe, a tidy data frame where each observation is a Pe value with five
#' variables: Sp, Value, Extinct, Survey, and scale

### Local Pe ----

TreatmentPostLocalPe <- matrix(nrow = nSpp, ncol = nTreatments)

for (i in 1:length(Treatment_Rarefactions_Region)){
  TreatmentPostLocalPe[,i] <- t(SBR_LocalPe(Treatment_Rarefactions_Region[[i]]))
}
rownames(TreatmentPostLocalPe) <- sppNames

#PeCategories Post1
LocalPeValuesPost1 <- as.data.frame(matrix(nrow = nSpp*nTreatments, ncol = 4))
colnames(LocalPeValuesPost1) <- c("Sp","Value", "Extinct", "Survey")

pos <- 0
for (i in 1:nSpp){
  for (j in 1:nTreatments){
    pos <- pos + 1
    LocalPeValuesPost1[pos,"Sp"] <- rownames(TreatmentPostLocalPe)[i]
    pe <- TreatmentPostLocalPe[i,j]
    LocalPeValuesPost1[pos,"Value"] <- pe
    if (Post1_Treatment_Region[i,j] == 0){
      LocalPeValuesPost1[pos,"Extinct"] <- 1
    } else{
      LocalPeValuesPost1[pos,"Extinct"] <- 0
    }
  }
}
LocalPeValuesPost1_Clean <- filter(LocalPeValuesPost1, LocalPeValuesPost1[,"Value"] != 1)
LocalPeValuesPost1_Clean$Survey <- "Post1"

#PeCategories Post2
LocalPeValuesPost2 <- as.data.frame(matrix(nrow = nSpp*nTreatments, ncol = 4))
colnames(LocalPeValuesPost2) <- c("Sp","Value", "Extinct", "Survey")

pos <- 0
for (i in 1:nSpp){
  for (j in 1:nTreatments){
    pos <- pos + 1
    LocalPeValuesPost2[pos,"Sp"] <- rownames(TreatmentPostLocalPe)[i]
    pe <- TreatmentPostLocalPe[i,j]
    LocalPeValuesPost2[pos,"Value"] <- pe
    if (Post2_Treatment_Region[i,j] == 0){
      LocalPeValuesPost2[pos,"Extinct"] <- 1
    } else{
      LocalPeValuesPost2[pos,"Extinct"] <- 0
    }
  }
}
LocalPeValuesPost2_Clean <- filter(LocalPeValuesPost2, LocalPeValuesPost2[,"Value"] != 1)
LocalPeValuesPost2_Clean$Survey <- "Post2"

#Combine Pe Values
TotalLocalPe <- rbind(LocalPeValuesPost1_Clean, LocalPeValuesPost2_Clean)
TotalLocalPe$scale <- "Local"

### Regional Pe ----
TreatmentRegionalPe <- t(SBR_RegionalPe(Treatment_Rarefactions_Region))
TreatmentRegionalPe <- as.data.frame(TreatmentRegionalPe)
colnames(TreatmentRegionalPe) <- "Value"
TreatmentRegionalPe$Sp <- sppNames
TreatmentRegionalPe$Post1<- 0
TreatmentRegionalPe$Post2 <- 0

for (i in 1:nSpp){
  if (TreatmentRegionalPe$Value[i] < 1){
    if (sum(Post1_Treatment_Region[i,]) == 0){
      TreatmentRegionalPe$Post1[i] <- 1
    }
    else{
      TreatmentRegionalPe$Post1[i]  <- 0
    }
    if (sum(Post2_Treatment_Region[i,]) == 0){
      TreatmentRegionalPe$Post2[i] <- 1
    }
    else{
      TreatmentRegionalPe$Post2[i]  <- 0
    }
  }
  
}
TreatmentRegionalPe_complete <- filter(TreatmentRegionalPe, Value < 1)
TotalRegionalPe <- gather(TreatmentRegionalPe_complete, Survey, Extinct, 3, 4)
TotalRegionalPe$scale <- "Regional"

#Combine all Pe Values
TotalPe <- rbind(TotalLocalPe, TotalRegionalPe)

### Rarefaction 95% Confidence Intervals ----
Total95CI <- data.frame("site" = 1:nTreatments, "survey" = 0, "lower" = 0, "upper" = 0, "colonists" = 0)

##Local Post1
LocalPost1CI <- data.frame("site" = 1:nTreatments, "survey" = "Post1", "lower" = 0, "upper" = 0, "colonists" = 0)

for (i in 1:nTreatments){ # for each treatment site
  site <- Treatment_Rarefactions_Region[[i]] #extract rarefied site simualtions
  vals <- SBR_richness(site) #extract null richness estimates for each site
  upper <- quantile(vals, 0.975) #determine upper 97.5% richness estimate
  lower <- quantile(vals, 0.025) #determine lower 2.5% richness estimate
  LocalPost1CI$site[i] <- treatmentSiteNames[i] #add site name to df
  LocalPost1CI$lower[i] <- lower #add lower 95%CI to df
  LocalPost1CI$upper[i] <- upper #add upper 95%CI to df
}

##Local Post1 + Colonists
LocalPost1CI_Colonists <- LocalPost1CI #initialize df to store CIs
LocalPost1CI_Colonists$colonists <- 1 #indicate that these include colonists


LocalPost1Richnesses <- NullTreatment_Richnesses_colonistsPost1$richness - 
  NullTreatmentPost1_Richnesses$richness #Determine number of colonist species for Post1

LocalPost1CI_Colonists$lower <- LocalPost1CI_Colonists$lower + LocalPost1Richnesses #add colonists to lower estimate
LocalPost1CI_Colonists$upper <- LocalPost1CI_Colonists$upper + LocalPost1Richnesses #add colonists to lower estimate

##Local Post2
LocalPost2CI <- LocalPost1CI
LocalPost2CI$survey <- "Post2" 

##Local Post2 + Colonists
LocalPost2CI_Colonists <- LocalPost2CI #initialize df to store CIs
LocalPost2CI_Colonists$colonists <- 1 #indicate that these include colonists


LocalPost2Richnesses <- NullTreatment_Richnesses_colonistsPost2$richness - 
  NullTreatmentPost2_Richnesses$richness #Determine number of colonist species for Post2

LocalPost2CI_Colonists$lower <- LocalPost2CI_Colonists$lower + LocalPost2Richnesses #add colonists to lower estimate
LocalPost2CI_Colonists$upper <- LocalPost2CI_Colonists$upper + LocalPost2Richnesses #add colonists to lower estimate


##Regional Post1 
RegionalPost1CI <- data.frame("site" = "Region", "survey" = "Post1", "lower" = 0, "upper" = 0, "colonists" = 0) #initialize df
RegionalPost1Richness <- SBR_Regional_richnesses(Treatment_Rarefactions_Region) #calclate regional rarefied richness estimates
RegionalPost1CI$lower <- quantile(RegionalPost1Richness,0.025) #extract lower 2.5% richness estimate
RegionalPost1CI$upper <- quantile(RegionalPost1Richness,0.975) #extract upper 97.5% richness estimate

##Regional Post1 + Colonists
RegionalPost1CI_Colonists <- RegionalPost1CI #intitialize df
RegionalPost1CI_Colonists$colonists <- 1 #indicate that this includes colonists
RegionalPost1CI_Colonists$lower <- RegionalPost1CI$lower + Post1_Regional_colonists #add colonists to lower estimate
RegionalPost1CI_Colonists$upper <- RegionalPost1CI$upper + Post1_Regional_colonists #add colonists to upper estimate

##Regional Post2 
RegionalPost2CI <-RegionalPost1CI
RegionalPost2CI$survey <- "Post2"

##Regional Post2 + Colonists
RegionalPost2CI_Colonists <- RegionalPost2CI #intitialize df
RegionalPost2CI_Colonists$colonists <- 1 #indicate that this includes colonists
RegionalPost2CI_Colonists$lower <- RegionalPost2CI$lower + Post2_Regional_colonists #add colonists to lower estimate
RegionalPost2CI_Colonists$upper <- RegionalPost2CI$upper + Post2_Regional_colonists #add colonists to upper estimate

#Combine all richness CI's to single df
TotalNull95CI <- rbind(LocalPost1CI, LocalPost1CI_Colonists, LocalPost2CI, LocalPost2CI_Colonists,
                       RegionalPost1CI, RegionalPost1CI_Colonists, RegionalPost2CI, RegionalPost2CI_Colonists)
#Finally, we need to add ENS CI's. Let's make a column for that.

TotalNull95CI$ENS <- 0

#Null ENS CI Post1
ENSPost1CI <- data.frame("site" = 1:nTreatments, "survey" = "Post1", "lower" = 0, 
                         "upper" = 0, "colonists" = 0, "ENS" = 1) #initialize df

for (i in 1:nTreatments){ #for each site
  site <- Treatment_Rarefactions_Region[[i]] #assign current site to object
  vals <- SBR_ENS(site) #calculate ENS estimates for sites
  ENSPost1CI$site[i] <- treatmentSiteNames[i] #assign site name
  ENSPost1CI$lower[i] <- quantile(vals, 0.025) #add lower 2.5% estimate to df
  ENSPost1CI$upper[i] <- quantile(vals, 0.975) #add upper 2.5% estimate to df
}

#ENS CI Post2
ENSPost2CI <- ENSPost1CI
ENSPost2CI$survey <- "Post2"

#Combine all CI
Total95CI <- rbind(TotalNull95CI, ENSPost1CI, ENSPost2CI)


#################### FINAL DATAFRAMES----
#check out data frame containing richness values
head(TotalRichnesses)
#check out data frame containing extinction values
head(TotalExtinctions)
#Check out data frame contianing Pe values
head(TotalPe)
#Check out data frame containing 95% confidence intervals for rarefaction estimates
head(TotalNull95CI)

#It may be helpful to have these dataframes on hand, so that this script doesn't need to be constantly run or so
#we can look at it in Excel. Let's write them to the working directory.
write.csv(TotalRichnesses, "SCR_TotalRichness_AlmeidaEtAlGCB.csv")
write.csv(TotalExtinctions, "SCR_TotalExtinction_AlmeidaEtAlGCB.csv")
write.csv(TotalPe, "SCR_TotalPe_AlmeidaEtAlGCB.csv")
write.csv(TotalNull95CI, "SCR_TotalNull95CI_AlmeidaEtAlGCB.csv")

#################### LOGISTIC REGRESSION----

###LOCAL POST1 LOGISTIC REGRESSION###

#Extract Local Post1 Pe values
LogRegLocalPe_Post1 <- filter(TotalPe, Survey == "Post1", scale == "Local")

#For ease of interpretation of the odds ratio, let's make Pe values scale from 0-100
LogRegLocalPe_Post1$Value <- 100*LogRegLocalPe_Post1$Value

#Run logistic regression
localMod <- glm(Extinct ~ Value, data = LogRegLocalPe_Post1, family = binomial(link = "logit"))

#check results
summary(localMod)

#New fitted values for full range of variation
Local_fitted_range <-  data.frame(Value = seq(0, 100, by = 0.5))
Local_fitted_range$fitted = predict(localMod, newdata = Local_fitted_range, type = "response")
Local_fitted_range$Pe <- Local_fitted_range$Value/100
Local_fitted_range$scale <- "Local"

###REGIONAL POST1 LOGREG###

#Extract regional Post1 Pe values
LogRegRegionalPe_Post1 <- filter(TotalPe, Survey == "Post1", scale == "Regional")

#For ease of interpretation of the odds ratio, let's make Pe values scale from 0-100
LogRegRegionalPe_Post1$Value <- 100*LogRegRegionalPe_Post1$Value

#Run logistic regression
regionalMod <- glm(Extinct ~ Value, data = LogRegRegionalPe_Post1, family = binomial(link = "logit"))

#New fitted values for full range of variation
Regional_fitted_range <-  data.frame(Value = seq(0, 100, by = 0.5))
Regional_fitted_range$fitted = predict(regionalMod, newdata = Regional_fitted_range, type = "response")
Regional_fitted_range$Pe <- Regional_fitted_range$Value/100
Regional_fitted_range$scale <- "Regional"

#Total fitted values
logRegFitted <- rbind(Local_fitted_range, Regional_fitted_range)

#check out results
summary(regionalMod)
#################### CREATE FIGURE 3A ----

#Extract pre-survey richness data 
PreS <- filter(TotalRichnesses, survey == "Pre", treatment == "Treatment")
PreS_mean <- summarise(group_by(PreS, scale), 
                       mean = mean(richness), lower = mean(richness) - 1.96*(sd(richness)/sqrt(nTreatments)), upper = mean(richness) + 1.96*(sd(richness)/sqrt(nTreatments)))
PreS_mean$group <- "Pre-disturbance"

#Extract null richness data without colonists
NullS <- filter(TotalRichnesses, survey == "Post1", treatment == "Treatment", colonists == 0, null == 1)
NullS_mean <- summarise(group_by(NullS, scale), 
                        mean = mean(richness), lower = mean(richness) - 1.96*(sd(richness)/sqrt(nTreatments)), upper = mean(richness) + 1.96*(sd(richness)/sqrt(nTreatments)))
NullS_mean$group <- "Spatially-clustered rarefaction"
NullCIs <- filter(TotalNull95CI, site == "Region", colonists == 0, survey  == "Post1")
NullS_mean[2,]$lower <- NullCIs$lower
NullS_mean[2,]$upper <- NullCIs$upper

#Extract SAR richness data without colonists
SARS <- filter(TotalRichnesses, survey == "Post1", treatment == "Treatment", colonists == 0, SAR == 1)
SARS_mean <- summarise(group_by(SARS, scale), 
                       mean = mean(richness), lower = mean(richness) - 1.96*(sd(richness)/sqrt(nTreatments)), upper = mean(richness) + 1.96*(sd(richness)/sqrt(nTreatments)))
SARS_mean$group <- "Species-Area Relationship"

#Extract Post1 SAR richness data with colonists
Post1S <- filter(TotalRichnesses, survey == "Post1", treatment == "Treatment", colonists == 1, null == 0, SAR == 0)
Post1S_mean <- summarise(group_by(Post1S, scale), 
                         mean = mean(richness), lower = mean(richness) - 1.96*(sd(richness)/sqrt(nTreatments)), upper = mean(richness) + 1.96*(sd(richness)/sqrt(nTreatments)))
Post1S_mean$group <- "Post-disturbance"

#Total data for fig2
Fig3adf <- rbind(PreS_mean, NullS_mean, SARS_mean, Post1S_mean)

#Re-order groups
Fig3adf$group <- factor(Fig3adf$group, 
                       levels = c("Pre-disturbance", "Species-Area Relationship", "Spatially-clustered rarefaction", "Post-disturbance"))

#Colorblind palette with black:
cbbPalette <- c("#000000", "#E69F00", "#619CFF", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#Create figure 3a
Figure3a <- ggplot(Fig3adf, aes(x = scale, y = mean, colour = group, group = group, shape = group)) + #plot mean richness value at local and regional scale for each group
  geom_errorbar(aes(ymin = lower, ymax = upper), width =0.05) + #add error bars to estimates
  geom_line(aes(linetype = group), size = 1.1) + #connect grouped values with line
  geom_point(size = 2.7) + #add points for clarity
  scale_shape_manual(values = c(17,16,15,18)) + #change point shapes to distinguish between groups
  scale_color_manual(values = cbbPalette) +
  scale_linetype_manual(values = c("solid", "longdash", "dashed", "dotted")) +
  ylab("Species richness") + #label y axis
  xlab(element_blank())+ #remove label from x axis
  theme_bw() + #remove greyscale background
  theme(legend.title = element_blank()) + #remove legend title
  theme(panel.grid = element_blank()) + #remove gridlines
  ylim(c(0,70))

#Check out Figure 3a
Figure3a

#Export high quality image using Cairo
# tiff(filename="Figure3a_AlmeidaEtAlGCB.tiff",
#      type="cairo",
#      units = 'mm',
#      width=150,
#      height=100,
#      pointsize=12,
#      res = 600)
# 
# Figure3a
# 
# dev.off()


#################### CREATE FIGURE 3B ----

#Extract null extinctions data without colonists
NullExtinctions <- filter(TotalExtinctions, survey == "Post1", treatment == "Treatment", null == 1)
NullExtinctions_mean <- summarise(group_by(NullExtinctions, scale), 
                                  mean = mean(extinctions), lower = mean(extinctions) - 1.96*(sd(extinctions)/sqrt(nTreatments)), upper = mean(extinctions) + 1.96*(sd(extinctions)/sqrt(nTreatments)))
NullExtinctions_mean$group <- "Spatially-clustered rarefaction"
NullCIs <- filter(TotalNull95CI, site == "Region", colonists == 0, survey  == "Post1")
NullExtinctions_mean[2,]$lower <- PreS_mean$mean[2] - NullCIs$upper
NullExtinctions_mean[2,]$upper <- PreS_mean$mean[2] -NullCIs$lower

#Extract SAR extinctions data without colonists
SARExtinctions <- filter(TotalExtinctions, survey == "Post1", treatment == "Treatment",  SAR == 1)
SARExtinctions_mean <- summarise(group_by(SARExtinctions, scale), 
                                 mean = mean(extinctions), lower = mean(extinctions) - 1.96*(sd(extinctions)/sqrt(nTreatments)), upper = mean(extinctions) + 1.96*(sd(extinctions)/sqrt(nTreatments)))
SARExtinctions_mean$group <- "Species-Area Relationship"

#Extract Post1 SAR extinctions data with colonists
Post1Extinctions <- filter(TotalExtinctions, survey == "Post1", treatment == "Treatment", colonists == 0, null == 0, SAR == 0)
Post1Extinctions_mean <- summarise(group_by(Post1Extinctions, scale), 
                                   mean = mean(extinctions), lower = mean(extinctions) - 1.96*(sd(extinctions)/sqrt(nTreatments)), upper = mean(extinctions) + 1.96*(sd(extinctions)/sqrt(nTreatments)))
Post1Extinctions_mean$group <- "Post-disturbance"

#Total data for fig3
Fig3bdf <- rbind(NullExtinctions_mean, SARExtinctions_mean, Post1Extinctions_mean)

#Re-order groups
Fig3bdf$group <- factor(Fig3bdf$group, levels = c("Species-Area Relationship", "Spatially-clustered rarefaction", "Post-disturbance"))

# Colorblind palette with black:
cbbPalette <- c("#000000", "#E69F00", "#619CFF", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#Create figure 3b
Figure3b <- ggplot(Fig3bdf, aes(x = scale, y = mean, colour = group, group = group, shape = group)) + #plot mean richness value at local and regional scale for each group
  geom_errorbar(aes(ymin = lower, ymax = upper), width =0.05) + #add error bars to estimates
  geom_line(aes(linetype = group), size = 1.1) + #connect grouped values with line
  geom_point(size = 2.7) + #add points for clarity
  scale_shape_manual(values = c(16,15,18)) + #change point shapes to distinguish between groups +
  scale_color_manual(values = cbbPalette[-1])+
  scale_linetype_manual(values = c("longdash", "dashed", "dotted")) +
  ylab("Species extinctions") + #label y axis
  xlab(element_blank())+ #remove label from x axis
  theme_bw() + #remove greyscale background
  theme(legend.title = element_blank()) + #remove legend title
  theme(panel.grid = element_blank()) #remove gridlines

#Check out Figure3b
Figure3b

#Export high quality image using Cairo
# tiff(filename="Figure3B_AlmeidaEtAlGCB.tiff",
#      type="cairo",
#      units = 'mm',
#      width=150,
#      height=100,
#      pointsize=12,
#      res = 600)
# 
# Figure3b
# 
# dev.off()

#################### CREATE FIGURE 4 ----

Post1Pe <- filter(TotalPe, Survey == "Post1")
Post1Pe$Extinct[Post1Pe$Extinct == 1] <- "Extinct"
Post1Pe$Extinct[Post1Pe$Extinct == 0] <- "Extant"
Post1Pe$fitted <- c(fitted(localMod), fitted(regionalMod))*85

Figure4 <- ggplot(Post1Pe, aes(x = Value, fill = Extinct)) +
  geom_line(aes(y = fitted), size  = 1.25) +
  geom_histogram( binwidth = 0.05, colour = "black") +
  facet_wrap(~scale) +
  scale_fill_manual(values = c("light blue", "salmon")) +
  scale_y_continuous(name = "Number of Species", sec.axis = sec_axis(~./85, name = "Regression-Predicted\\nExtinction Probability")) +
  theme_bw() + #remove greyscale background
  theme(legend.title = element_blank()) + #remove legend title
  theme(panel.grid = element_blank()) + #remove gridlines
  xlab("Null-Simulated Extinction Probability")

#check out Figure4
Figure4
  

#Export high quality image using Cairo
png(filename="Figure4_AlmeidaEtAlGCB.png",
     type="cairo",
     units = 'mm',
     width=300,
     height=160,
     pointsize=12,
     res = 600)

Figure4

dev.off()
