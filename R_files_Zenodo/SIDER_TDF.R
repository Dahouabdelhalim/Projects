#Ryan Stephens; Finalized Dec 15, 2020

rm(list=ls()) # clears workspace
setwd("~/Isotopic_Routing_Small_Mammals/SIDER")

#########################################################################################
#Data from literature review of TDFs
#########################################################################################
TDF_Review<- read.csv("Literature_Review_TDF_data.csv",header=TRUE)#data includes all tissue types, not just hair
str(TDF_Review)
head(TDF_Review)
library(dplyr)
library(tidyr)


##############################
#Clean up data 
##############################
TDF_Review<-TDF_Review %>% mutate(TDF_13C = ifelse(#If diet was lipid corrected, than remove TDF for d13C
                                  Lipid_corrected_diet == "Y" & TDF_13C != "NA", NA, TDF_13C))

TDF_Review<-TDF_Review %>% mutate(TDF_15N = ifelse(#Remove field study where sd of d15N TDF was very large (11.6)
                                  Reference == "Codron et al. (2012)" & Common_name == "Springbok", NA, TDF_15N))   

TDF_Review<-TDF_Review %>% mutate(TDF_15N = ifelse(#Remove field study where TDF for d15N may have been influenced by a diet switch
                                  Reference == "Codron et al. (2012)" & Common_name == "Blesbok", NA, TDF_15N))   

#Remove data from field study since these are getting tested (generally this would otherwise be retained)                               
TDF_Review<-TDF_Review %>% filter(Reference != "Stephens et al. (2021)")#TDF used to test the efficacy of field studies

TDF_Review<-TDF_Review %>% filter(Tissue != "Skin")#Remove skin since it isn't included in SIDER

TDF_Review %>% count(!is.na(TDF_13C))#number of carbon TDFs
TDF_Review %>% count(!is.na(TDF_15N))#number of nitrogen TDFs
##############################
                                  

##############################
#Put data in SIDER format
##############################
#remove duplicate data that are also in SIDER
TDF_Review<-TDF_Review %>% replace_na(list(In_SIDER_database = ""))#change NA to blank
TDF_Review<-TDF_Review %>% filter(In_SIDER_database != "x") 

TDF_Review<- mutate(TDF_Review,species = paste(Genus,Species,sep = '_'))#new column for species
TDF_Review$taxonomic.class<-"mammalia"#add taxonomic.class

#select needed columns, order and rename to match SIDER
Review_isotope_data<-select(TDF_Review, species, habitat=Habitat, taxonomic.class, tissue=Tissue, diet.type=Class,
             source.iso.13C=Diet_d13C, source.iso.15N=Diet_d15N,
             delta13C=TDF_13C, delta15N=TDF_15N, citation=Citation )
head(Review_isotope_data)

library(forcats)
#rename levels of diet.type to match SIDER
levels(Review_isotope_data$diet.type)
Review_isotope_data<-mutate(Review_isotope_data, diet.type = fct_recode(diet.type,"carnivore" = "Carnivore",
                                                                                  "herbivore" = "Herbivore",
                                                                                  "Omnivore" = "Omnivore"))
#rename levels of tissue to match SIDER (collapse all blood types into "blood")
levels(Review_isotope_data$tissue)
Review_isotope_data<-mutate(Review_isotope_data, tissue = fct_recode(tissue,"blood" = "Blood Cells",
                                                                            "blood" = "Blood Plasma",
                                                                            "blood" = "Blood Serum",
                                                                            "blood" = "Red Blood Cells",
                                                                            "blood" = "Whole Blood",
                                                                            "collagen" = "Collagen",
                                                                            "hair" = "Hair",
                                                                            "liver" = "Liver",
                                                                            "muscle" = "Muscle"))
str(Review_isotope_data)
#Change factor to character to match SIDER
Review_isotope_data$species<-as.character(Review_isotope_data$species)
Review_isotope_data$habitat<-as.character(Review_isotope_data$habitat)
Review_isotope_data$tissue<-as.character(Review_isotope_data$tissue)
Review_isotope_data$diet.type<-as.character(Review_isotope_data$diet.type)
Review_isotope_data$citation<-as.character(Review_isotope_data$citation)
##############################

head(Review_isotope_data)
str(Review_isotope_data)
#########################################################################################




#########################################################################################
#SIDER models
#########################################################################################
#install.packages("remotes")
#remotes::install_github("healyke/SIDER")
library(SIDER)
#code modified from:
#https://github.com/healyke/SIDER/blob/master/vignettes/Introduction-to-SIDER.Rmd
#https://github.com/healyke/SIDER/blob/master/vignettes/Add-new-data-for-SIDER-run.Rmd


###########################################
#General model settings
###########################################

#Isotope data
isotope_data_SIDER <- scrumpSider(iso.data = "all")#save all isotope dataset as an object
str(isotope_data_SIDER)
isotope_data_SIDER %>% count(!is.na(delta13C), taxonomic.class)#number of carbon TDFs for each taxonomic group
isotope_data_SIDER %>% count(!is.na(delta15N), taxonomic.class)#number of nitrogen TDFs for each taxonomic group


isotope_data_SIDER_Updated <- rbind(isotope_data_SIDER, Review_isotope_data)#SIDER data updated with review data
str(isotope_data_SIDER_Updated)
isotope_data_SIDER_Updated %>% count(!is.na(delta13C), taxonomic.class)#number of carbon TDFs for each taxonomic group
isotope_data_SIDER_Updated %>% count(!is.na(delta15N), taxonomic.class)#number of nitrogen TDFs for each taxonomic group

#Tree data
combined_trees <- scrumpSider(tree = "all")

#Fixed effects for C & N models
formula_c <- delta13C ~ diet.type + habitat
formula_n <- delta15N ~ diet.type + habitat

#Random structure shared by all models
random_terms <- ( ~ animal + species + tissue)
prior <- list(R = list(V = 1, nu=0.002), 
              G = list(G1=list(V = 1, nu=0.002),
                       G2=list(V = 1, nu=0.002), 
                       G3=list(V = 1, nu=0.002)))

# model run settings
nitt   <- c(1200000)
burnin <- c(200000)
thin   <- c(500)
parameters <- c(nitt, thin, burnin)
n_chains <- c(2)

# convergence settings
convergence =  c(1.1)
ESS = c(1000)
###########################################


#########################################################################################
#Myodes gapperi
#########################################################################################
#########################################################################################
#Myodes gapperi - terrestrial herbivore
MYGA <- recipeSider(species = "Myodes_gapperi", 
                              habitat = "terrestrial", 
                              taxonomic.class = "mammalia", 
                              tissue = "hair", 
                              diet.type = "herbivore", 
                              tree = combined_trees)



#Carbon TDF SIDER model 
MYGA_tdf_SIDER_data_c <- prepareSider(MYGA,
                           isotope_data_SIDER,#SIDER data only 
                           combined_trees, 
                           "carbon")

MYGA_C_SIDER <- imputeSider(mulTree.data = MYGA_tdf_SIDER_data_c, 
                             formula = formula_c, 
                             random.terms = random_terms,
                             prior = prior, 
                             output = "MYGA_C_SIDER",
                             parameters = parameters,
                             chains = n_chains, 
                             convergence =  convergence, 
                             ESS = ESS)

#Nitrogen TDF SIDER model 
MYGA_tdf_SIDER_data_n <- prepareSider(MYGA,
                           isotope_data_SIDER,#SIDER data only 
                           combined_trees, 
                           "nitrogen")

MYGA_N_SIDER <- imputeSider(mulTree.data = MYGA_tdf_SIDER_data_n, 
                             formula = formula_n, 
                             random.terms = random_terms,
                             prior = prior, 
                             output = "MYGA_N_SIDER",
                             parameters = parameters,
                             chains = n_chains, 
                             convergence =  convergence, 
                             ESS = ESS)


#Carbon TDF SIDER updated model 
MYGA_tdf_SIDER_updated_data_c <- prepareSider(MYGA,
                                 isotope_data_SIDER_Updated,#SIDER updated data
                                 combined_trees, 
                                 "carbon")

MYGA_C_SIDER_Updated <- imputeSider(mulTree.data = MYGA_tdf_SIDER_updated_data_c, 
                                 formula = formula_c, 
                                 random.terms = random_terms,
                                 prior = prior, 
                                 output = "MYGA_C_SIDER_Updated",
                                 parameters = parameters,
                                 chains = n_chains, 
                                 convergence =  convergence, 
                                 ESS = ESS)

#Nitrogen TDF SIDER updated model 
MYGA_tdf_SIDER_updated_data_n <- prepareSider(MYGA,
                                 isotope_data_SIDER_Updated,#SIDER updated data 
                                 combined_trees, 
                                 "nitrogen")

MYGA_N_SIDER_Updated <- imputeSider(mulTree.data = MYGA_tdf_SIDER_updated_data_n, 
                                 formula = formula_n, 
                                 random.terms = random_terms,
                                 prior = prior, 
                                 output = "MYGA_N_SIDER_Updated",
                                 parameters = parameters,
                                 chains = n_chains, 
                                 convergence =  convergence, 
                                 ESS = ESS)

#save.image(file='MYGA.RData')#Save model runs
load('MYGA.RData')#Load in model runs

#bundle the global estimates of TDFs from each of the 4 models into a wide format data.frame
MYGA_wide <- data.frame(MYGA_C_SIDER = as.numeric(MYGA_C_SIDER$tdf_global),
                        MYGA_C_Update = as.numeric(MYGA_C_SIDER_Updated$tdf_global), 
                        MYGA_N_SIDER  = as.numeric(MYGA_N_SIDER$tdf_global),
                        MYGA_N_Update = as.numeric(MYGA_N_SIDER_Updated$tdf_global))

summary(MYGA_wide)# get some summary statistics of these estimates
MYGA_wide_long <- tidyr::gather(MYGA_wide, model, TDF)# tidy to long format for nice plotting

MYGA_TDF_Summary <- MYGA_wide_long %>%# calculate means and standard deviations 
  group_by(model) %>% 
  summarise(mean = round(mean(TDF),2), sd = round(sd(TDF),2))

MYGA_TDF_Summary$Species <- "M. gapperi"#add species
MYGA_TDF_Summary$Speces_Abr <- "MYGA"#add species abreviation
#########################################################################################



#########################################################################################
#Peromyscus maniculatus
#########################################################################################
#########################################################################################
#Peromyscus maniculatus- terrestrial omnivore
PEMA <- recipeSider(species = "Peromyscus_maniculatus", 
                         habitat = "terrestrial", 
                         taxonomic.class = "mammalia", 
                         tissue = "hair", 
                         diet.type = "omnivore", 
                         tree = combined_trees)



#Carbon TDF SIDER model 
PEMA_tdf_SIDER_data_c <- prepareSider(PEMA,
                                 isotope_data_SIDER,#SIDER data only 
                                 combined_trees, 
                                 "carbon")

PEMA_C_SIDER <- imputeSider(mulTree.data = PEMA_tdf_SIDER_data_c, 
                                 formula = formula_c, 
                                 random.terms = random_terms,
                                 prior = prior, 
                                 output = "PEMA_C_SIDER",
                                 parameters = parameters,
                                 chains = n_chains, 
                                 convergence =  convergence, 
                                 ESS = ESS)

#Nitrogen TDF SIDER model 
PEMA_tdf_SIDER_data_n <- prepareSider(PEMA,
                                 isotope_data_SIDER,#SIDER data only 
                                 combined_trees, 
                                 "nitrogen")

PEMA_N_SIDER <- imputeSider(mulTree.data = PEMA_tdf_SIDER_data_n, 
                                 formula = formula_n, 
                                 random.terms = random_terms,
                                 prior = prior, 
                                 output = "PEMA_N_SIDER",
                                 parameters = parameters,
                                 chains = n_chains, 
                                 convergence =  convergence, 
                                 ESS = ESS)


#Carbon TDF SIDER updated model 
PEMA_tdf_SIDER_updated_data_c <- prepareSider(PEMA,
                                         isotope_data_SIDER_Updated,#SIDER updated data
                                         combined_trees, 
                                         "carbon")

PEMA_C_SIDER_Updated <- imputeSider(mulTree.data = PEMA_tdf_SIDER_updated_data_c, 
                                         formula = formula_c, 
                                         random.terms = random_terms,
                                         prior = prior, 
                                         output = "PEMA_C_SIDER_Updated",
                                         parameters = parameters,
                                         chains = n_chains, 
                                         convergence =  convergence, 
                                         ESS = ESS)

#Nitrogen TDF SIDER updated model 
PEMA_tdf_SIDER_updated_data_n <- prepareSider(PEMA,
                                         isotope_data_SIDER_Updated,#SIDER updated data 
                                         combined_trees, 
                                         "nitrogen")

PEMA_N_SIDER_Updated <- imputeSider(mulTree.data = PEMA_tdf_SIDER_updated_data_n, 
                                         formula = formula_n, 
                                         random.terms = random_terms,
                                         prior = prior, 
                                         output = "PEMA_N_SIDER_Updated",
                                         parameters = parameters,
                                         chains = n_chains, 
                                         convergence =  convergence, 
                                         ESS = ESS)

#save.image(file='PEMA.RData')#Save model runs
load('PEMA.RData')#Load in model runs

#bundle the global estimates of TDFs from each of the 4 models into a wide format data.frame
PEMA_wide <- data.frame(PEMA_C_SIDER = as.numeric(PEMA_C_SIDER$tdf_global),
                        PEMA_C_Update = as.numeric(PEMA_C_SIDER_Updated$tdf_global), 
                        PEMA_N_SIDER  = as.numeric(PEMA_N_SIDER$tdf_global),
                        PEMA_N_Update = as.numeric(PEMA_N_SIDER_Updated$tdf_global))

summary(PEMA_wide)# get some summary statistics of these estimates
PEMA_wide_long <- tidyr::gather(PEMA_wide, model, TDF)# tidy to long format for nice plotting

PEMA_TDF_Summary <- PEMA_wide_long %>%# calculate means and standard deviations 
  group_by(model) %>% 
  summarise(mean = round(mean(TDF),2), sd = round(sd(TDF),2))

PEMA_TDF_Summary$Species <- "P. maniculatus"#add species
PEMA_TDF_Summary$Speces_Abr <- "PEMA"#add species abreviation
#########################################################################################




#########################################################################################
#Blarina brevicauda
#########################################################################################
#########################################################################################
#Blarina brevicauda- terrestrial carnivore
BLBR <- recipeSider(species = "Blarina_brevicauda", 
                    habitat = "terrestrial", 
                    taxonomic.class = "mammalia", 
                    tissue = "hair", 
                    diet.type = "carnivore", 
                    tree = combined_trees)



#Carbon TDF SIDER model 
BLBR_tdf_SIDER_data_c <- prepareSider(BLBR,
                                      isotope_data_SIDER,#SIDER data only 
                                      combined_trees, 
                                      "carbon")

BLBR_C_SIDER <- imputeSider(mulTree.data = BLBR_tdf_SIDER_data_c, 
                            formula = formula_c, 
                            random.terms = random_terms,
                            prior = prior, 
                            output = "BLBR_C_SIDER",
                            parameters = parameters,
                            chains = n_chains, 
                            convergence =  convergence, 
                            ESS = ESS)

#Nitrogen TDF SIDER model 
BLBR_tdf_SIDER_data_n <- prepareSider(BLBR,
                                      isotope_data_SIDER,#SIDER data only 
                                      combined_trees, 
                                      "nitrogen")

BLBR_N_SIDER <- imputeSider(mulTree.data = BLBR_tdf_SIDER_data_n, 
                            formula = formula_n, 
                            random.terms = random_terms,
                            prior = prior, 
                            output = "BLBR_N_SIDER",
                            parameters = parameters,
                            chains = n_chains, 
                            convergence =  convergence, 
                            ESS = ESS)


#Carbon TDF SIDER updated model 
BLBR_tdf_SIDER_updated_data_c <- prepareSider(BLBR,
                                              isotope_data_SIDER_Updated,#SIDER updated data
                                              combined_trees, 
                                              "carbon")

BLBR_C_SIDER_Updated <- imputeSider(mulTree.data = BLBR_tdf_SIDER_updated_data_c, 
                                    formula = formula_c, 
                                    random.terms = random_terms,
                                    prior = prior, 
                                    output = "BLBR_C_SIDER_Updated",
                                    parameters = parameters,
                                    chains = n_chains, 
                                    convergence =  convergence, 
                                    ESS = ESS)

#Nitrogen TDF SIDER updated model 
BLBR_tdf_SIDER_updated_data_n <- prepareSider(BLBR,
                                              isotope_data_SIDER_Updated,#SIDER updated data 
                                              combined_trees, 
                                              "nitrogen")

BLBR_N_SIDER_Updated <- imputeSider(mulTree.data = BLBR_tdf_SIDER_updated_data_n, 
                                    formula = formula_n, 
                                    random.terms = random_terms,
                                    prior = prior, 
                                    output = "BLBR_N_SIDER_Updated",
                                    parameters = parameters,
                                    chains = n_chains, 
                                    convergence =  convergence, 
                                    ESS = ESS)

#save.image(file='BLBR.RData')#Save model runs
load('BLBR.RData')#Load in model runs

#bundle the global estimates of TDFs from each of the 4 models into a wide format data.frame
BLBR_wide <- data.frame(BLBR_C_SIDER = as.numeric(BLBR_C_SIDER$tdf_global),
                        BLBR_C_Update = as.numeric(BLBR_C_SIDER_Updated$tdf_global), 
                        BLBR_N_SIDER  = as.numeric(BLBR_N_SIDER$tdf_global),
                        BLBR_N_Update = as.numeric(BLBR_N_SIDER_Updated$tdf_global))

summary(BLBR_wide)# get some summary statistics of these estimates
BLBR_wide_long <- tidyr::gather(BLBR_wide, model, TDF)# tidy to long format for nice plotting

BLBR_TDF_Summary <- BLBR_wide_long %>%# calculate means and standard deviations 
  group_by(model) %>% 
  summarise(mean = round(mean(TDF),2), sd = round(sd(TDF),2))

BLBR_TDF_Summary$Species <- "B. brevicauda"#add species
BLBR_TDF_Summary$Speces_Abr <- "BLBR"#add species abreviation
#########################################################################################




#########################################################################################
#Napaeozapus insignis
#########################################################################################
#########################################################################################
#Napaeozapus insignis - terrestrial omnivore
NAIN <- recipeSider(species = "Napaeozapus_insignis", 
                    habitat = "terrestrial", 
                    taxonomic.class = "mammalia", 
                    tissue = "hair", 
                    diet.type = "omnivore", 
                    tree = combined_trees)



#Carbon TDF SIDER model 
NAIN_tdf_SIDER_data_c <- prepareSider(NAIN,
                                      isotope_data_SIDER,#SIDER data only 
                                      combined_trees, 
                                      "carbon")

NAIN_C_SIDER <- imputeSider(mulTree.data = NAIN_tdf_SIDER_data_c, 
                            formula = formula_c, 
                            random.terms = random_terms,
                            prior = prior, 
                            output = "NAIN_C_SIDER",
                            parameters = parameters,
                            chains = n_chains, 
                            convergence =  convergence, 
                            ESS = ESS)

#Nitrogen TDF SIDER model 
NAIN_tdf_SIDER_data_n <- prepareSider(NAIN,
                                      isotope_data_SIDER,#SIDER data only 
                                      combined_trees, 
                                      "nitrogen")

NAIN_N_SIDER <- imputeSider(mulTree.data = NAIN_tdf_SIDER_data_n, 
                            formula = formula_n, 
                            random.terms = random_terms,
                            prior = prior, 
                            output = "NAIN_N_SIDER",
                            parameters = parameters,
                            chains = n_chains, 
                            convergence =  convergence, 
                            ESS = ESS)


#Carbon TDF SIDER updated model 
NAIN_tdf_SIDER_updated_data_c <- prepareSider(NAIN,
                                              isotope_data_SIDER_Updated,#SIDER updated data
                                              combined_trees, 
                                              "carbon")

NAIN_C_SIDER_Updated <- imputeSider(mulTree.data = NAIN_tdf_SIDER_updated_data_c, 
                                    formula = formula_c, 
                                    random.terms = random_terms,
                                    prior = prior, 
                                    output = "NAIN_C_SIDER_Updated",
                                    parameters = parameters,
                                    chains = n_chains, 
                                    convergence =  convergence, 
                                    ESS = ESS)

#Nitrogen TDF SIDER updated model 
NAIN_tdf_SIDER_updated_data_n <- prepareSider(NAIN,
                                              isotope_data_SIDER_Updated,#SIDER updated data 
                                              combined_trees, 
                                              "nitrogen")

NAIN_N_SIDER_Updated <- imputeSider(mulTree.data = NAIN_tdf_SIDER_updated_data_n, 
                                    formula = formula_n, 
                                    random.terms = random_terms,
                                    prior = prior, 
                                    output = "NAIN_N_SIDER_Updated",
                                    parameters = parameters,
                                    chains = n_chains, 
                                    convergence =  convergence, 
                                    ESS = ESS)

#save.image(file='NAIN.RData')#Save model runs
load('NAIN.RData')#Load in model runs

#bundle the global estimates of TDFs from each of the 4 models into a wide format data.frame
NAIN_wide <- data.frame(NAIN_C_SIDER = as.numeric(NAIN_C_SIDER$tdf_global),
                        NAIN_C_Update = as.numeric(NAIN_C_SIDER_Updated$tdf_global), 
                        NAIN_N_SIDER  = as.numeric(NAIN_N_SIDER$tdf_global),
                        NAIN_N_Update = as.numeric(NAIN_N_SIDER_Updated$tdf_global))

summary(NAIN_wide)# get some summary statistics of these estimates
NAIN_wide_long <- tidyr::gather(NAIN_wide, model, TDF)# tidy to long format for nice plotting

NAIN_TDF_Summary <- NAIN_wide_long %>%# calculate means and standard deviations 
  group_by(model) %>% 
  summarise(mean = round(mean(TDF),2), sd = round(sd(TDF),2))

NAIN_TDF_Summary$Species <- "N. insignis"#add species
NAIN_TDF_Summary$Speces_Abr <- "NAIN"#add species abreviation
#########################################################################################





#########################################################################################
#All small mammals
#########################################################################################

All<-rbind(MYGA_TDF_Summary, PEMA_TDF_Summary, BLBR_TDF_Summary, NAIN_TDF_Summary)
head(All)

All<-mutate(All, Isotope = ifelse(#Add isotope - if string has "_C_" then put d13C, if not put d15N
  grepl("_C_", model),"d13C", "d15N"))
All<-mutate(All, Method = ifelse(#Add method - if string has "SIDER" then put d13C, if not put d15N
  grepl("SIDER", model),"SIDER", "SIDER_update"))

SIDER_TDF<-select(All, Species, Speces_Abr, Method, Isotope, TDF_Mean = mean, TDF_SD = sd)
head(SIDER_TDF)

write.csv(SIDER_TDF, "SIDER_TDF.csv", row.names = F)





