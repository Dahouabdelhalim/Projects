#Code to run mixing models for each small mammal species (scripts have been merged here)
#Mydoes gapperi: lines 8 - 965
#Napaeozapus insignis: lines 975 - 1933
#Peromyscus maniculatus: lines 1944 - 2901
#Blarina bevicauda: lines 2912 - 3869


################################################################################
#Mydoes gapperi
################################################################################
rm(list=ls())# clears workspace
setwd("~/Isotopic_Routing_Small_Mammals/Mixing_Models/MYGA")#set working directory


################################################################################
#Diet
################################################################################

#############################
#load isotopic data
#############################
#Take stomach samples, convert to diet, and select species
library(dplyr)
Isotopes<- read.csv("Field_study_data.csv",header=T)
Diet<-filter(Isotopes, Type =="Stomach")#select stomach data
Diet<-Diet %>% mutate(d15N=ifelse(Species == "MYGA",d15N-1.85,d15N-1))#subtract 1 to account digestive enzymes of stomach (1.85 for MGYA)
Diet<-mutate(Diet, Type = "Diet")#code 'type' as diet
MYGA_Diet<- filter(Diet,Species =="MYGA") #Filter out to desired species

write.csv(MYGA_Diet, "MYGA_Diet.csv")

library(MixSIAR)
mix.filename.MYGA_Diet <-("~/Isotopic_Routing_Small_Mammals/Mixing_Models/MYGA/MYGA_Diet.csv") #load isotopic data for MixSIAR

mix.MYGA.Diet <- load_mix_data(filename=mix.filename.MYGA_Diet, #name of the CSV file with mix/consumer data
                               iso_names=c("d13C","d15N"), #tracers/isotopes
                               factors=c(NULL),            #factors
                               fac_random=c(NULL),         #random effects
                               fac_nested=c(NULL),         #Nested effects
                               cont_effects=NULL)          #continuous effects
#############################


#############################
# Load food source data
#############################
source.filename <- ("~/Isotopic_Routing_Small_Mammals/Mixing_Models/MYGA/Bartlett_diet_items.csv")

source_MYGA_Diet <- load_source_data(filename=source.filename, #name of the CSV file with food source data
                                     source_factors=NULL,     #column headings of random/fixed effects you have source data by
                                     conc_dep=TRUE,           #TRUE or FALSE for concentration dependence data 
                                     data_type="raw",         #"means" = means + sd or "raw" = raw data
                                     mix.MYGA.Diet)

# Load discrimination data
#No TDF for diet samples
Type<- c("AM Fungi","Arthropods", "Red Maple", "EM Fungi", "Berries")
TDF_MYGA_Field <- as.data.frame(Type)
TDF_MYGA_Field$Meand13C<-0	
TDF_MYGA_Field$SDd13C<-0		
TDF_MYGA_Field$Meand15N<-0		
TDF_MYGA_Field$SDd15N<-0	
write.csv(TDF_MYGA_Field, "Discrimination_Diet.csv", row.names=FALSE)

Discrimination_Diet <- ("~/Isotopic_Routing_Small_Mammals/Mixing_Models/MYGA/Discrimination_Diet.csv")
discr_Diet <- load_discr_data(filename=Discrimination_Diet, mix.MYGA.Diet)
#############################


#############################
# Make isospace plot
#############################
plot_data(filename="MYGA_Diet_isospace", #name to save the isospace plot as
          plot_save_pdf=FALSE,       #TRUE or FALSE, should mix.SIAR save the plot as a .pdf?
          plot_save_png=TRUE,      #TRUE or FALSE, should mix.SIAR save the plot as a .png?
          mix.MYGA.Diet,source_MYGA_Diet,discr_Diet)

# If 2 isotopes/tracers, calculate normalized surface area of the convex hull polygon(s)
#   *Note 1: discrimination SD is added to the source SD (see calc_area.r for details)
if(mix.MYGA.Diet$n.iso==2) calc_area(source=source_MYGA_Diet,mix=mix.MYGA.Diet,discr=discr_Diet)
#############################


#############################
# Define prior
#############################
kw.alpha <- c(1,0,1,1,1)

kw.alpha <- kw.alpha*length(kw.alpha)/sum(kw.alpha)# Generate alpha hyperparameters scaling sum(alpha)=n.sources

kw.alpha[which(kw.alpha==0)] <- 0.01# the Dirichlet hyperparameters for the alpha.prior cannot be 0 (but can set = .01)

#plot using "plot_prior"
plot_prior(alpha.prior=kw.alpha,
           source=source_MYGA_Diet,
           plot_save_pdf=FALSE,
           plot_save_png=TRUE,
           filename="MYGA_Diet_prior_plot_kw_inf")
#   RED = your prior
#   DARK GREY = "uninformative"/generalist (alpha = 1)
#   LIGHT GREY = "uninformative" Jeffrey's prior (alpha = 1/n.sources)
#############################


#############################
# Write JAGS model file (define model structure)
##############################
# There are 3 error term options available:
#   1. Residual * Process (resid_err = TRUE, process_err = TRUE)
#   2. Residual only (resid_err = TRUE, process_err = FALSE)
#   3. Process only (resid_err = FALSE, process_err = TRUE)
MYGA_Diet_model <- "MYGA_Diet_model.txt"
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(MYGA_Diet_model, resid_err, process_err, mix.MYGA.Diet, source_MYGA_Diet)
##############################


##############################
# Run MYGA Diet model
##############################
# MCMC options:
#run <- list(chainLength=30000, burn=15000, thin=25, chains=3, calcDIC=TRUE)
costum_run <- list(chainLength=200000, burn=50000, thin=50, chains=3, calcDIC=TRUE)

#Model
MYGA_Diet<- run_model(run=costum_run,mix.MYGA.Diet,source_MYGA_Diet,discr_Diet,MYGA_Diet_model,
                      alpha.prior=kw.alpha, resid_err, process_err)
##############################


##############################
# Process JAGS output
##############################
# Choose output options (see ?output_options for details)
output_options <- list(summary_save = TRUE,                 
                       summary_name = "MYGA_Diet_summary_statistics", 
                       sup_post = FALSE,                    
                       plot_post_save_pdf = TRUE,           
                       plot_post_name = "MYGA_Diet_posterior_density",
                       sup_pairs = FALSE,             
                       plot_pairs_save_pdf = TRUE,    
                       plot_pairs_name = "MYGA_Diet_pairs_plot",
                       sup_xy = TRUE,           
                       plot_xy_save_pdf = TRUE,
                       plot_xy_name = "MYGA_Diet_xy_plot",
                       gelman = TRUE,
                       heidel = FALSE,  
                       geweke = TRUE,   
                       diag_save = TRUE,
                       diag_name = "MYGA_Diet_diagnostics",
                       indiv_effect = FALSE,       
                       plot_post_save_png = TRUE, 
                       plot_pairs_save_png = TRUE,
                       plot_xy_save_png = TRUE)
##############################


##############################
# Create diagnostics, summary statistics, and posterior plots
##############################
output_JAGS(MYGA_Diet, mix.MYGA.Diet, source_MYGA_Diet, output_options)
print(MYGA_Diet)
getwd()
save.image("MYGA_Diet_Jags.RData") 
##############################


##################################################
#Put data into a dataframe
##################################################
setwd("~/Isotopic_Routing_Small_Mammals/Mixing_Models/MYGA")

require(R2jags) # Load packages 

# Load data 
load("~/Isotopic_Routing_Small_Mammals/Mixing_Models/MYGA/MYGA_Diet_Jags.RData")
attach.jags(MYGA_Diet) # call model file
print(MYGA_Diet)

#assign names to each food source https://github.com/brianstock/MixSIAR/issues/30
MYGA_Diet_Post <- data.frame(Source = "Diet", AM_Fungi = p.global[,1],Arthropods = p.global[,2],
                             Berries = p.global[,3], EM_Fungi = p.global[,4], Red_Maple = p.global[,5], 
                             Fungi = p.global[,1]+p.global[,4]) 

MYGA_Diet_Post<-data.frame(Number=rownames(MYGA_Diet_Post), MYGA_Diet_Post)

require(tidyr)
MYGA_Diet_Posterior <- MYGA_Diet_Post %>% gather(Source,value,AM_Fungi:Fungi)#Use tidyr to gather posterior for each food source
MYGA_Diet_Posterior["Model"] <- "Diet"#Add model type
MYGA_Diet_Posterior["Species"] <- "MYGA"#Add Species
#save raw posterior data
write.csv(MYGA_Diet_Posterior, file = "MYGA_Diet_Posterior.csv", row.names=F)#save raw posterior data
write.csv(MYGA_Diet_Posterior, file = "~/Isotopic_Routing_Small_Mammals/Mixing_Models/Mixing_Model_Figure/MYGA_Diet_Posterior.csv", row.names=F)
################################################################################











################################################################################
#Field
################################################################################
setwd("~/Isotopic_Routing_Small_Mammals/Mixing_Models/MYGA")#set working directory
#############################
#load isotopic data
#############################
#Take stomach samples, convert to Field, and select species
Hair<- read.csv("Bartlett_isotope_hair_data.csv",header=T)
library(dplyr)
MYGA_Field<-filter(Hair, Species=="MYGA")
write.csv(MYGA_Field, "MYGA_Field.csv")

library(MixSIAR)
mix.filename.MYGA_Field <-("~/Isotopic_Routing_Small_Mammals/Mixing_Models/MYGA/MYGA_Field.csv") #load isotopic data for MixSIAR

mix.MYGA.Field <- load_mix_data(filename=mix.filename.MYGA_Field, #name of the CSV file with mix/consumer data
                                iso_names=c("d13C","d15N"), #tracers/isotopes
                                factors=c(NULL),            #factors
                                fac_random=c(NULL),         #random effects
                                fac_nested=c(NULL),         #Nested effects
                                cont_effects=NULL)          #continuous effects
#############################


#############################
# Load food source data
#############################
source.filename <- ("~/Isotopic_Routing_Small_Mammals/Mixing_Models/MYGA/Bartlett_diet_items.csv")

source_MYGA_Field <- load_source_data(filename=source.filename, #name of the CSV file with food source data
                                      source_factors=NULL,     #column headings of random/fixed effects you have source data by
                                      conc_dep=TRUE,           #TRUE or FALSE for concentration dependence data 
                                      data_type="raw",         #"means" = means + sd or "raw" = raw data
                                      mix.MYGA.Field)

# Load discrimination data
#Select species and method of estimated TDF
Est_DFT<- read.csv("Estimated_TDF.csv",header=T)#discrimination factors
MYGA_Field_d13C_TDF<-filter(Est_DFT,Species_Abr == "MYGA" & Method =="Field" & Isotope == "d13C")
MYGA_Field_d15N_TDF<-filter(Est_DFT,Species_Abr == "MYGA" & Method =="Field" & Isotope == "d15N")


Type<- c("AM Fungi","Arthropods","Red Maple", "EM Fungi", "Berries")
TDF_MYGA_Field <- as.data.frame(Type)
TDF_MYGA_Field$Meand13C<-MYGA_Field_d13C_TDF$TDF_Mean	
TDF_MYGA_Field$SDd13C<-MYGA_Field_d13C_TDF$TDF_SD			
TDF_MYGA_Field$Meand15N<-MYGA_Field_d15N_TDF$TDF_Mean			
TDF_MYGA_Field$SDd15N<-MYGA_Field_d15N_TDF$TDF_SD		
write.csv(TDF_MYGA_Field, "MYGA_Discrimination_Field.csv", row.names=FALSE)

MYGA_Discrimination_Field <- ("~/Isotopic_Routing_Small_Mammals/Mixing_Models/MYGA/MYGA_Discrimination_Field.csv")
MYGA_Discrimination_Field <- load_discr_data(filename=MYGA_Discrimination_Field, mix.MYGA.Field)
#############################


#############################
# Make isospace plot
#############################
plot_data(filename="MYGA_Field_isospace", #name to save the isospace plot as
          plot_save_pdf=FALSE,       #TRUE or FALSE, should mix.SIAR save the plot as a .pdf?
          plot_save_png=TRUE,      #TRUE or FALSE, should mix.SIAR save the plot as a .png?
          mix.MYGA.Field,source_MYGA_Field,MYGA_Discrimination_Field)

# If 2 isotopes/tracers, calculate normalized surface area of the convex hull polygon(s)
#   *Note 1: discrimination SD is added to the source SD (see calc_area.r for details)
if(mix.MYGA.Field$n.iso==2) calc_area(source=source_MYGA_Field,mix=mix.MYGA.Field,discr=MYGA_Discrimination_Field)
#############################


#############################
# Define prior
#############################
kw.alpha <- c(1,0,1,1,1)

kw.alpha <- kw.alpha*length(kw.alpha)/sum(kw.alpha)# Generate alpha hyperparameters scaling sum(alpha)=n.sources

kw.alpha[which(kw.alpha==0)] <- 0.01# the Dirichlet hyperparameters for the alpha.prior cannot be 0 (but can set = .01)

#plot using "plot_prior"
plot_prior(alpha.prior=kw.alpha,
           source=source_MYGA_Field,
           plot_save_pdf=FALSE,
           plot_save_png=TRUE,
           filename="MYGA_Field_prior_plot_kw_inf")
#   RED = your prior
#   DARK GREY = "uninformative"/generalist (alpha = 1)
#   LIGHT GREY = "uninformative" Jeffrey's prior (alpha = 1/n.sources)
#############################


#############################
# Write JAGS model file (define model structure)
##############################
# There are 3 error term options available:
#   1. Residual * Process (resid_err = TRUE, process_err = TRUE)
#   2. Residual only (resid_err = TRUE, process_err = FALSE)
#   3. Process only (resid_err = FALSE, process_err = TRUE)
MYGA_Field_model <- "MYGA_Field_model.txt"
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(MYGA_Field_model, resid_err, process_err, mix.MYGA.Field, source_MYGA_Field)
##############################


##############################
# Run MYGA Field model
##############################
# MCMC options:
#run <- list(chainLength=30000, burn=15000, thin=25, chains=3, calcDIC=TRUE)
costum_run <- list(chainLength=200000, burn=50000, thin=50, chains=3, calcDIC=TRUE)

#Model
MYGA_Field<- run_model(run=costum_run,mix.MYGA.Field,source_MYGA_Field,MYGA_Discrimination_Field,MYGA_Field_model,
                       alpha.prior=kw.alpha, resid_err, process_err)
##############################


##############################
# Process JAGS output
##############################
# Choose output options (see ?output_options for details)
output_options <- list(summary_save = TRUE,                 
                       summary_name = "MYGA_Field_summary_statistics", 
                       sup_post = FALSE,                    
                       plot_post_save_pdf = TRUE,           
                       plot_post_name = "MYGA_Field_posterior_density",
                       sup_pairs = FALSE,             
                       plot_pairs_save_pdf = TRUE,    
                       plot_pairs_name = "MYGA_Field_pairs_plot",
                       sup_xy = TRUE,           
                       plot_xy_save_pdf = FALSE,
                       plot_xy_name = "MYGA_Field_xy_plot",
                       gelman = TRUE,
                       heidel = FALSE,  
                       geweke = TRUE,   
                       diag_save = TRUE,
                       diag_name = "MYGA_Field_diagnostics",
                       indiv_effect = FALSE,       
                       plot_post_save_png = FALSE, 
                       plot_pairs_save_png = FALSE,
                       plot_xy_save_png = FALSE)
##############################


##############################
# Create diagnostics, summary statistics, and posterior plots
##############################
output_JAGS(MYGA_Field, mix.MYGA.Field, source_MYGA_Field, output_options)
print(MYGA_Field)
getwd()
save.image("MYGA_Field_Jags.RData") 
##############################


##################################################
#Put data into a dataframe
##################################################
setwd("~/Isotopic_Routing_Small_Mammals/Mixing_Models/MYGA")

require(R2jags) # Load packages 

# Load data 
load("~/Isotopic_Routing_Small_Mammals/Mixing_Models/MYGA/MYGA_Field_Jags.RData")
attach.jags(MYGA_Field) # call model file
print(MYGA_Field)

#assign names to each food source https://github.com/brianstock/MixSIAR/issues/30
MYGA_Field_Post <- data.frame(Source = "Field", AM_Fungi = p.global[,1],Arthropods = p.global[,2],
                              Berries = p.global[,3], EM_Fungi = p.global[,4], Red_Maple = p.global[,5], 
                              Fungi = p.global[,1]+p.global[,4]) 

MYGA_Field_Post<-data.frame(Number=rownames(MYGA_Field_Post), MYGA_Field_Post)#add posterior number
require(tidyr)
MYGA_Field_Posterior <- MYGA_Field_Post %>% gather(Source,value,AM_Fungi:Fungi)#Use tidyr to gather posterior for each food source
MYGA_Field_Posterior["Model"] <- "Field"#Add model type
MYGA_Field_Posterior["Species"] <- "MYGA"#Add Species

#save raw posterior data
write.csv(MYGA_Field_Posterior, file = "MYGA_Field_Posterior.csv", row.names=FALSE)#save raw posterior data
write.csv(MYGA_Field_Posterior, file = "~/Isotopic_Routing_Small_Mammals/Mixing_Models/Mixing_Model_Figure/MYGA_Field_Posterior.csv", row.names=FALSE)
################################################################################










################################################################################
#SIDER
################################################################################
setwd("~/Isotopic_Routing_Small_Mammals/Mixing_Models/MYGA")#set working directory
#############################
#load isotopic data
#############################
#Take stomach samples, convert to SIDER, and select species
Hair<- read.csv("Bartlett_isotope_hair_data.csv",header=T)
library(dplyr)
MYGA_SIDER<-filter(Hair, Species=="MYGA")
write.csv(MYGA_SIDER, "MYGA_SIDER.csv")

library(MixSIAR)
mix.filename.MYGA_SIDER <-("~/Isotopic_Routing_Small_Mammals/Mixing_Models/MYGA/MYGA_SIDER.csv") #load isotopic data for MixSIAR

mix.MYGA.SIDER <- load_mix_data(filename=mix.filename.MYGA_SIDER, #name of the CSV file with mix/consumer data
                                iso_names=c("d13C","d15N"), #tracers/isotopes
                                factors=c(NULL),            #factors
                                fac_random=c(NULL),         #random effects
                                fac_nested=c(NULL),         #Nested effects
                                cont_effects=NULL)          #continuous effects
#############################


#############################
# Load food source data
#############################
source.filename <- ("~/Isotopic_Routing_Small_Mammals/Mixing_Models/MYGA/Bartlett_diet_items.csv")

source_MYGA_SIDER <- load_source_data(filename=source.filename, #name of the CSV file with food source data
                                      source_factors=NULL,     #column headings of random/fixed effects you have source data by
                                      conc_dep=TRUE,           #TRUE or FALSE for concentration dependence data 
                                      data_type="raw",         #"means" = means + sd or "raw" = raw data
                                      mix.MYGA.SIDER)

# Load discrimination data
#Select species and method of estimated TDF
Est_DFT<- read.csv("Estimated_TDF.csv",header=T)#discrimination factors
MYGA_SIDER_d13C_TDF<-filter(Est_DFT,Species_Abr == "MYGA" & Method =="SIDER" & Isotope == "d13C")
MYGA_SIDER_d15N_TDF<-filter(Est_DFT,Species_Abr == "MYGA" & Method =="SIDER" & Isotope == "d15N")


Type<- c("AM Fungi","Arthropods", "Red Maple", "EM Fungi", "Berries")
TDF_MYGA_SIDER <- as.data.frame(Type)
TDF_MYGA_SIDER$Meand13C<-MYGA_SIDER_d13C_TDF$TDF_Mean	
TDF_MYGA_SIDER$SDd13C<-MYGA_SIDER_d13C_TDF$TDF_SD			
TDF_MYGA_SIDER$Meand15N<-MYGA_SIDER_d15N_TDF$TDF_Mean			
TDF_MYGA_SIDER$SDd15N<-MYGA_SIDER_d15N_TDF$TDF_SD		
write.csv(TDF_MYGA_SIDER, "MYGA_Discrimination_SIDER.csv", row.names=FALSE)

MYGA_Discrimination_SIDER <- ("~/Isotopic_Routing_Small_Mammals/Mixing_Models/MYGA/MYGA_Discrimination_SIDER.csv")
MYGA_Discrimination_SIDER <- load_discr_data(filename=MYGA_Discrimination_SIDER, mix.MYGA.SIDER)
#############################


#############################
# Make isospace plot
#############################
plot_data(filename="MYGA_SIDER_isospace", #name to save the isospace plot as
          plot_save_pdf=FALSE,       #TRUE or FALSE, should mix.SIAR save the plot as a .pdf?
          plot_save_png=TRUE,      #TRUE or FALSE, should mix.SIAR save the plot as a .png?
          mix.MYGA.SIDER,source_MYGA_SIDER,MYGA_Discrimination_SIDER)

# If 2 isotopes/tracers, calculate normalized surface area of the convex hull polygon(s)
#   *Note 1: discrimination SD is added to the source SD (see calc_area.r for details)
if(mix.MYGA.SIDER$n.iso==2) calc_area(source=source_MYGA_SIDER,mix=mix.MYGA.SIDER,discr=MYGA_Discrimination_SIDER)
#############################


#############################
# Define prior
#############################
kw.alpha <- c(1,0,1,1,1)

kw.alpha <- kw.alpha*length(kw.alpha)/sum(kw.alpha)# Generate alpha hyperparameters scaling sum(alpha)=n.sources

kw.alpha[which(kw.alpha==0)] <- 0.01# the Dirichlet hyperparameters for the alpha.prior cannot be 0 (but can set = .01)

#plot using "plot_prior"
plot_prior(alpha.prior=kw.alpha,
           source=source_MYGA_SIDER,
           plot_save_pdf=FALSE,
           plot_save_png=TRUE,
           filename="MYGA_SIDER_prior_plot_kw_inf")
#   RED = your prior
#   DARK GREY = "uninformative"/generalist (alpha = 1)
#   LIGHT GREY = "uninformative" Jeffrey's prior (alpha = 1/n.sources)
#############################


#############################
# Write JAGS model file (define model structure)
##############################
# There are 3 error term options available:
#   1. Residual * Process (resid_err = TRUE, process_err = TRUE)
#   2. Residual only (resid_err = TRUE, process_err = FALSE)
#   3. Process only (resid_err = FALSE, process_err = TRUE)
MYGA_SIDER_model <- "MYGA_SIDER_model.txt"
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(MYGA_SIDER_model, resid_err, process_err, mix.MYGA.SIDER, source_MYGA_SIDER)
##############################


##############################
# Run MYGA SIDER model
##############################
# MCMC options:
run <- list(chainLength=30000, burn=15000, thin=25, chains=3, calcDIC=TRUE)
costum_run <- list(chainLength=200000, burn=50000, thin=50, chains=3, calcDIC=TRUE)

#Model
MYGA_SIDER<- run_model(run=costum_run,mix.MYGA.SIDER,source_MYGA_SIDER,MYGA_Discrimination_SIDER,MYGA_SIDER_model,
                       alpha.prior=kw.alpha, resid_err, process_err)
##############################


##############################
# Process JAGS output
##############################
# Choose output options (see ?output_options for details)
output_options <- list(summary_save = TRUE,                 
                       summary_name = "MYGA_SIDER_summary_statistics", 
                       sup_post = FALSE,                    
                       plot_post_save_pdf = TRUE,           
                       plot_post_name = "MYGA_SIDER_posterior_density",
                       sup_pairs = FALSE,             
                       plot_pairs_save_pdf = TRUE,    
                       plot_pairs_name = "MYGA_SIDER_pairs_plot",
                       sup_xy = TRUE,           
                       plot_xy_save_pdf = FALSE,
                       plot_xy_name = "MYGA_SIDER_xy_plot",
                       gelman = TRUE,
                       heidel = FALSE,  
                       geweke = TRUE,   
                       diag_save = TRUE,
                       diag_name = "MYGA_SIDER_diagnostics",
                       indiv_effect = FALSE,       
                       plot_post_save_png = FALSE, 
                       plot_pairs_save_png = FALSE,
                       plot_xy_save_png = FALSE)
##############################


##############################
# Create diagnostics, summary statistics, and posterior plots
##############################
output_JAGS(MYGA_SIDER, mix.MYGA.SIDER, source_MYGA_SIDER, output_options)
print(MYGA_SIDER)
getwd()
save.image("MYGA_SIDER_Jags.RData") 
##############################


##################################################
#Put data into a dataframe
##################################################
setwd("~/Isotopic_Routing_Small_Mammals/Mixing_Models/MYGA")

require(R2jags) # Load packages 

# Load data 
load("~/Isotopic_Routing_Small_Mammals/Mixing_Models/MYGA/MYGA_SIDER_Jags.RData")
attach.jags(MYGA_SIDER) # call model file
print(MYGA_SIDER)

#assign names to each food source https://github.com/brianstock/MixSIAR/issues/30
MYGA_SIDER_Post <- data.frame(Source = "SIDER", AM_Fungi = p.global[,1],Arthropods = p.global[,2],
                              Berries = p.global[,3], EM_Fungi = p.global[,4], Red_Maple = p.global[,5], 
                              Fungi = p.global[,1]+p.global[,4]) 

MYGA_SIDER_Post<-data.frame(Number=rownames(MYGA_SIDER_Post), MYGA_SIDER_Post)#add posterior number
require(tidyr)
MYGA_SIDER_Posterior <- MYGA_SIDER_Post %>% gather(Source,value,AM_Fungi:Fungi)#Use tidyr to gather posterior for each food source
MYGA_SIDER_Posterior["Model"] <- "SIDER"#Add model type
MYGA_SIDER_Posterior["Species"] <- "MYGA"#Add Species

#save raw posterior data
write.csv(MYGA_SIDER_Posterior, file = "MYGA_SIDER_Posterior.csv", row.names=FALSE)#save raw posterior data
write.csv(MYGA_SIDER_Posterior, file = "~/Isotopic_Routing_Small_Mammals/Mixing_Models/Mixing_Model_Figure/MYGA_SIDER_Posterior.csv", row.names=FALSE)
################################################################################











################################################################################
#SIDER_update
################################################################################
setwd("~/Isotopic_Routing_Small_Mammals/Mixing_Models/MYGA")#set working directory
#############################
#load isotopic data
#############################
#Take stomach samples, convert to SIDER_update, and select species
Hair<- read.csv("Bartlett_isotope_hair_data.csv",header=T)
library(dplyr)
MYGA_SIDER_update<-filter(Hair, Species=="MYGA")
write.csv(MYGA_SIDER_update, "MYGA_SIDER_update.csv")

library(MixSIAR)
mix.filename.MYGA_SIDER_update <-("~/Isotopic_Routing_Small_Mammals/Mixing_Models/MYGA/MYGA_SIDER_update.csv") #load isotopic data for MixSIAR

mix.MYGA.SIDER_update <- load_mix_data(filename=mix.filename.MYGA_SIDER_update, #name of the CSV file with mix/consumer data
                                       iso_names=c("d13C","d15N"), #tracers/isotopes
                                       factors=c(NULL),            #factors
                                       fac_random=c(NULL),         #random effects
                                       fac_nested=c(NULL),         #Nested effects
                                       cont_effects=NULL)          #continuous effects
#############################


#############################
# Load food source data
#############################
source.filename <- ("~/Isotopic_Routing_Small_Mammals/Mixing_Models/MYGA/Bartlett_diet_items.csv")

source_MYGA_SIDER_update <- load_source_data(filename=source.filename, #name of the CSV file with food source data
                                             source_factors=NULL,     #column headings of random/fixed effects you have source data by
                                             conc_dep=TRUE,           #TRUE or FALSE for concentration dependence data 
                                             data_type="raw",         #"means" = means + sd or "raw" = raw data
                                             mix.MYGA.SIDER_update)

# Load discrimination data
#Select species and method of estimated TDF
Est_DFT<- read.csv("Estimated_TDF.csv",header=T)#discrimination factors
MYGA_SIDER_update_d13C_TDF<-filter(Est_DFT,Species_Abr == "MYGA" & Method =="SIDER_update" & Isotope == "d13C")
MYGA_SIDER_update_d15N_TDF<-filter(Est_DFT,Species_Abr == "MYGA" & Method =="SIDER_update" & Isotope == "d15N")


Type<- c("AM Fungi","Arthropods", "Red Maple", "EM Fungi", "Berries")
TDF_MYGA_SIDER_update <- as.data.frame(Type)
TDF_MYGA_SIDER_update$Meand13C<-MYGA_SIDER_update_d13C_TDF$TDF_Mean	
TDF_MYGA_SIDER_update$SDd13C<-MYGA_SIDER_update_d13C_TDF$TDF_SD			
TDF_MYGA_SIDER_update$Meand15N<-MYGA_SIDER_update_d15N_TDF$TDF_Mean			
TDF_MYGA_SIDER_update$SDd15N<-MYGA_SIDER_update_d15N_TDF$TDF_SD		
write.csv(TDF_MYGA_SIDER_update, "MYGA_Discrimination_SIDER_update.csv", row.names=FALSE)

MYGA_Discrimination_SIDER_update <- ("~/Isotopic_Routing_Small_Mammals/Mixing_Models/MYGA/MYGA_Discrimination_SIDER_update.csv")
MYGA_Discrimination_SIDER_update <- load_discr_data(filename=MYGA_Discrimination_SIDER_update, mix.MYGA.SIDER_update)
#############################


#############################
# Make isospace plot
#############################
plot_data(filename="MYGA_SIDER_update_isospace", #name to save the isospace plot as
          plot_save_pdf=FALSE,       #TRUE or FALSE, should mix.SIAR save the plot as a .pdf?
          plot_save_png=TRUE,      #TRUE or FALSE, should mix.SIAR save the plot as a .png?
          mix.MYGA.SIDER_update,source_MYGA_SIDER_update,MYGA_Discrimination_SIDER_update)

# If 2 isotopes/tracers, calculate normalized surface area of the convex hull polygon(s)
#   *Note 1: discrimination SD is added to the source SD (see calc_area.r for details)
if(mix.MYGA.SIDER_update$n.iso==2) calc_area(source=source_MYGA_SIDER_update,mix=mix.MYGA.SIDER_update,discr=MYGA_Discrimination_SIDER_update)
#############################


#############################
# Define prior
#############################
kw.alpha <- c(1,0,1,1,1)

kw.alpha <- kw.alpha*length(kw.alpha)/sum(kw.alpha)# Generate alpha hyperparameters scaling sum(alpha)=n.sources

kw.alpha[which(kw.alpha==0)] <- 0.01# the Dirichlet hyperparameters for the alpha.prior cannot be 0 (but can set = .01)

#plot using "plot_prior"
plot_prior(alpha.prior=kw.alpha,
           source=source_MYGA_SIDER_update,
           plot_save_pdf=FALSE,
           plot_save_png=TRUE,
           filename="MYGA_SIDER_update_prior_plot_kw_inf")
#   RED = your prior
#   DARK GREY = "uninformative"/generalist (alpha = 1)
#   LIGHT GREY = "uninformative" Jeffrey's prior (alpha = 1/n.sources)
#############################


#############################
# Write JAGS model file (define model structure)
##############################
# There are 3 error term options available:
#   1. Residual * Process (resid_err = TRUE, process_err = TRUE)
#   2. Residual only (resid_err = TRUE, process_err = FALSE)
#   3. Process only (resid_err = FALSE, process_err = TRUE)
MYGA_SIDER_update_model <- "MYGA_SIDER_update_model.txt"
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(MYGA_SIDER_update_model, resid_err, process_err, mix.MYGA.SIDER_update, source_MYGA_SIDER_update)
##############################


##############################
# Run MYGA SIDER_update model
##############################
# MCMC options:
#run <- list(chainLength=30000, burn=15000, thin=25, chains=3, calcDIC=TRUE)
costum_run <- list(chainLength=200000, burn=50000, thin=50, chains=3, calcDIC=TRUE)

#Model
MYGA_SIDER_update<- run_model(run=costum_run,mix.MYGA.SIDER_update,source_MYGA_SIDER_update,MYGA_Discrimination_SIDER_update,MYGA_SIDER_update_model,
                              alpha.prior=kw.alpha, resid_err, process_err)
##############################


##############################
# Process JAGS output
##############################
# Choose output options (see ?output_options for details)
output_options <- list(summary_save = TRUE,                 
                       summary_name = "MYGA_SIDER_update_summary_statistics", 
                       sup_post = FALSE,                    
                       plot_post_save_pdf = TRUE,           
                       plot_post_name = "MYGA_SIDER_update_posterior_density",
                       sup_pairs = FALSE,             
                       plot_pairs_save_pdf = TRUE,    
                       plot_pairs_name = "MYGA_SIDER_update_pairs_plot",
                       sup_xy = TRUE,           
                       plot_xy_save_pdf = FALSE,
                       plot_xy_name = "MYGA_SIDER_update_xy_plot",
                       gelman = TRUE,
                       heidel = FALSE,  
                       geweke = TRUE,   
                       diag_save = TRUE,
                       diag_name = "MYGA_SIDER_update_diagnostics",
                       indiv_effect = FALSE,       
                       plot_post_save_png = FALSE, 
                       plot_pairs_save_png = FALSE,
                       plot_xy_save_png = FALSE)
##############################


##############################
# Create diagnostics, summary statistics, and posterior plots
##############################
output_JAGS(MYGA_SIDER_update, mix.MYGA.SIDER_update, source_MYGA_SIDER_update, output_options)
print(MYGA_SIDER_update)
getwd()
save.image("MYGA_SIDER_update_Jags.RData") 
##############################


##################################################
#Put data into a dataframe
##################################################
setwd("~/Isotopic_Routing_Small_Mammals/Mixing_Models/MYGA")

require(R2jags) # Load packages 

# Load data 
load("~/Isotopic_Routing_Small_Mammals/Mixing_Models/MYGA/MYGA_SIDER_update_Jags.RData")
attach.jags(MYGA_SIDER_update) # call model file
print(MYGA_SIDER_update)

#assign names to each food source https://github.com/brianstock/MixSIAR/issues/30
MYGA_SIDER_update_Post <- data.frame(Source = "SIDER_update", AM_Fungi = p.global[,1],Arthropods = p.global[,2],
                                     Berries = p.global[,3], EM_Fungi = p.global[,4], Red_Maple = p.global[,5], 
                                     Fungi = p.global[,1]+p.global[,4]) 

MYGA_SIDER_update_Post<-data.frame(Number=rownames(MYGA_SIDER_update_Post), MYGA_SIDER_update_Post)#add posterior number
require(tidyr)
MYGA_SIDER_update_Posterior <- MYGA_SIDER_update_Post %>% gather(Source,value,AM_Fungi:Fungi)#Use tidyr to gather posterior for each food source
MYGA_SIDER_update_Posterior["Model"] <- "SIDER_update"#Add model type
MYGA_SIDER_update_Posterior["Species"] <- "MYGA"#Add Species

#save raw posterior data
write.csv(MYGA_SIDER_update_Posterior, file = "MYGA_SIDER_update_Posterior.csv", row.names=FALSE)#save raw posterior data
write.csv(MYGA_SIDER_update_Posterior, file = "~/Isotopic_Routing_Small_Mammals/Mixing_Models/Mixing_Model_Figure/MYGA_SIDER_update_Posterior.csv", row.names=FALSE)
################################################################################









################################################################################
#Meta_analysis
################################################################################
setwd("~/Isotopic_Routing_Small_Mammals/Mixing_Models/MYGA")#set working directory
#############################
#load isotopic data
#############################
#Take stomach samples, convert to Meta_analysis, and select species
Hair<- read.csv("Bartlett_isotope_hair_data.csv",header=T)
library(dplyr)
MYGA_Meta_analysis<-filter(Hair, Species=="MYGA")
write.csv(MYGA_Meta_analysis, "MYGA_Meta_analysis.csv")

library(MixSIAR)
mix.filename.MYGA_Meta_analysis <-("~/Isotopic_Routing_Small_Mammals/Mixing_Models/MYGA/MYGA_Meta_analysis.csv") #load isotopic data for MixSIAR

mix.MYGA.Meta_analysis <- load_mix_data(filename=mix.filename.MYGA_Meta_analysis, #name of the CSV file with mix/consumer data
                                        iso_names=c("d13C","d15N"), #tracers/isotopes
                                        factors=c(NULL),            #factors
                                        fac_random=c(NULL),         #random effects
                                        fac_nested=c(NULL),         #Nested effects
                                        cont_effects=NULL)          #continuous effects
#############################


#############################
# Load food source data
#############################
source.filename <- ("~/Isotopic_Routing_Small_Mammals/Mixing_Models/MYGA/Bartlett_diet_items.csv")

source_MYGA_Meta_analysis <- load_source_data(filename=source.filename, #name of the CSV file with food source data
                                              source_factors=NULL,     #column headings of random/fixed effects you have source data by
                                              conc_dep=TRUE,           #TRUE or FALSE for concentration dependence data 
                                              data_type="raw",         #"means" = means + sd or "raw" = raw data
                                              mix.MYGA.Meta_analysis)

# Load discrimination data
#Select species and method of estimated TDF
Est_DFT<- read.csv("Estimated_TDF.csv",header=T)#discrimination factors
MYGA_Meta_analysis_d13C_TDF<-filter(Est_DFT,Species_Abr == "MYGA" & Method =="Meta_analysis" & Isotope == "d13C")
MYGA_Meta_analysis_d15N_TDF<-filter(Est_DFT,Species_Abr == "MYGA" & Method =="Meta_analysis" & Isotope == "d15N")


Type<- c("AM Fungi","Arthropods", "Red Maple", "EM Fungi", "Berries")
TDF_MYGA_Meta_analysis <- as.data.frame(Type)
TDF_MYGA_Meta_analysis$Meand13C<-MYGA_Meta_analysis_d13C_TDF$TDF_Mean	
TDF_MYGA_Meta_analysis$SDd13C<-MYGA_Meta_analysis_d13C_TDF$TDF_SD			
TDF_MYGA_Meta_analysis$Meand15N<-MYGA_Meta_analysis_d15N_TDF$TDF_Mean			
TDF_MYGA_Meta_analysis$SDd15N<-MYGA_Meta_analysis_d15N_TDF$TDF_SD		
write.csv(TDF_MYGA_Meta_analysis, "MYGA_Discrimination_Meta_analysis.csv", row.names=FALSE)

MYGA_Discrimination_Meta_analysis <- ("~/Isotopic_Routing_Small_Mammals/Mixing_Models/MYGA/MYGA_Discrimination_Meta_analysis.csv")
MYGA_Discrimination_Meta_analysis <- load_discr_data(filename=MYGA_Discrimination_Meta_analysis, mix.MYGA.Meta_analysis)
#############################


#############################
# Make isospace plot
#############################
plot_data(filename="MYGA_Meta_analysis_isospace", #name to save the isospace plot as
          plot_save_pdf=FALSE,       #TRUE or FALSE, should mix.SIAR save the plot as a .pdf?
          plot_save_png=TRUE,      #TRUE or FALSE, should mix.SIAR save the plot as a .png?
          mix.MYGA.Meta_analysis,source_MYGA_Meta_analysis,MYGA_Discrimination_Meta_analysis)

# If 2 isotopes/tracers, calculate normalized surface area of the convex hull polygon(s)
#   *Note 1: discrimination SD is added to the source SD (see calc_area.r for details)
if(mix.MYGA.Meta_analysis$n.iso==2) calc_area(source=source_MYGA_Meta_analysis,mix=mix.MYGA.Meta_analysis,discr=MYGA_Discrimination_Meta_analysis)
#############################


#############################
# Define prior
#############################
kw.alpha <- c(1,0,1,1,1)

kw.alpha <- kw.alpha*length(kw.alpha)/sum(kw.alpha)# Generate alpha hyperparameters scaling sum(alpha)=n.sources

kw.alpha[which(kw.alpha==0)] <- 0.01# the Dirichlet hyperparameters for the alpha.prior cannot be 0 (but can set = .01)

#plot using "plot_prior"
plot_prior(alpha.prior=kw.alpha,
           source=source_MYGA_Meta_analysis,
           plot_save_pdf=FALSE,
           plot_save_png=TRUE,
           filename="MYGA_Meta_analysis_prior_plot_kw_inf")
#   RED = your prior
#   DARK GREY = "uninformative"/generalist (alpha = 1)
#   LIGHT GREY = "uninformative" Jeffrey's prior (alpha = 1/n.sources)
#############################


#############################
# Write JAGS model file (define model structure)
##############################
# There are 3 error term options available:
#   1. Residual * Process (resid_err = TRUE, process_err = TRUE)
#   2. Residual only (resid_err = TRUE, process_err = FALSE)
#   3. Process only (resid_err = FALSE, process_err = TRUE)
MYGA_Meta_analysis_model <- "MYGA_Meta_analysis_model.txt"
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(MYGA_Meta_analysis_model, resid_err, process_err, mix.MYGA.Meta_analysis, source_MYGA_Meta_analysis)
##############################


##############################
# Run MYGA Meta_analysis model
##############################
# MCMC options:
run <- list(chainLength=30000, burn=15000, thin=25, chains=3, calcDIC=TRUE)
costum_run <- list(chainLength=200000, burn=50000, thin=50, chains=3, calcDIC=TRUE)

#Model
MYGA_Meta_analysis<- run_model(run=costum_run,mix.MYGA.Meta_analysis,source_MYGA_Meta_analysis,MYGA_Discrimination_Meta_analysis,MYGA_Meta_analysis_model,
                               alpha.prior=kw.alpha, resid_err, process_err)
##############################


##############################
# Process JAGS output
##############################
# Choose output options (see ?output_options for details)
output_options <- list(summary_save = TRUE,                 
                       summary_name = "MYGA_Meta_analysis_summary_statistics", 
                       sup_post = FALSE,                    
                       plot_post_save_pdf = TRUE,           
                       plot_post_name = "MYGA_Meta_analysis_posterior_density",
                       sup_pairs = FALSE,             
                       plot_pairs_save_pdf = TRUE,    
                       plot_pairs_name = "MYGA_Meta_analysis_pairs_plot",
                       sup_xy = TRUE,           
                       plot_xy_save_pdf = FALSE,
                       plot_xy_name = "MYGA_Meta_analysis_xy_plot",
                       gelman = TRUE,
                       heidel = FALSE,  
                       geweke = TRUE,   
                       diag_save = TRUE,
                       diag_name = "MYGA_Meta_analysis_diagnostics",
                       indiv_effect = FALSE,       
                       plot_post_save_png = FALSE, 
                       plot_pairs_save_png = FALSE,
                       plot_xy_save_png = FALSE)
##############################


##############################
# Create diagnostics, summary statistics, and posterior plots
##############################
output_JAGS(MYGA_Meta_analysis, mix.MYGA.Meta_analysis, source_MYGA_Meta_analysis, output_options)
print(MYGA_Meta_analysis)
getwd()
save.image("MYGA_Meta_analysis_Jags.RData") 
##############################


##################################################
#Put data into a dataframe
##################################################
setwd("~/Isotopic_Routing_Small_Mammals/Mixing_Models/MYGA")

require(R2jags) # Load packages 

# Load data 
load("~/Isotopic_Routing_Small_Mammals/Mixing_Models/MYGA/MYGA_Meta_analysis_Jags.RData")
attach.jags(MYGA_Meta_analysis) # call model file
print(MYGA_Meta_analysis)

#assign names to each food source https://github.com/brianstock/MixSIAR/issues/30
MYGA_Meta_analysis_Post <- data.frame(Source = "Meta_analysis", AM_Fungi = p.global[,1],Arthropods = p.global[,2],
                                      Berries = p.global[,3], EM_Fungi = p.global[,4], Red_Maple = p.global[,5], 
                                      Fungi = p.global[,1]+p.global[,4]) 

MYGA_Meta_analysis_Post<-data.frame(Number=rownames(MYGA_Meta_analysis_Post), MYGA_Meta_analysis_Post)#add posterior number

require(tidyr)
MYGA_Meta_analysis_Posterior <- MYGA_Meta_analysis_Post %>% gather(Source,value,AM_Fungi:Fungi)#Use tidyr to gather posterior for each food source
MYGA_Meta_analysis_Posterior["Model"] <- "Meta_analysis"#Add model type
MYGA_Meta_analysis_Posterior["Species"] <- "MYGA"#Add Species

#save raw posterior data
write.csv(MYGA_Meta_analysis_Posterior, file = "MYGA_Meta_analysis_Posterior.csv", row.names=FALSE)#save raw posterior data
write.csv(MYGA_Meta_analysis_Posterior, file = "~/Isotopic_Routing_Small_Mammals/Mixing_Models/Mixing_Model_Figure/MYGA_Meta_analysis_Posterior.csv", row.names=FALSE)
################################################################################









################################################################################
#Napaeozapus insignis
################################################################################
rm(list=ls())# clears workspace
setwd("~/Isotopic_Routing_Small_Mammals/Mixing_Models/NAIN")#set working directory


################################################################################
#Diet
################################################################################

#############################
#load isotopic data
#############################
#Take stomach samples, convert to diet, and select species
library(dplyr)
Isotopes<- read.csv("Field_study_data.csv",header=T)
Diet<-filter(Isotopes, Type =="Stomach")#select stomach data
Diet<-Diet %>% mutate(d15N=ifelse(Species == "MYGA",d15N-1.85,d15N-1))#subtract 1 to account digestive enzymes of stomach (1.85 for MGYA)
Diet<-mutate(Diet, Type = "Diet")#code 'type' as diet
NAIN_Diet<- filter(Diet,Species =="NAIN") #Filter out to desired species

write.csv(NAIN_Diet, "NAIN_Diet.csv")

library(MixSIAR)
mix.filename.NAIN_Diet <-("~/Isotopic_Routing_Small_Mammals/Mixing_Models/NAIN/NAIN_Diet.csv") #load isotopic data for MixSIAR

mix.NAIN.Diet <- load_mix_data(filename=mix.filename.NAIN_Diet, #name of the CSV file with mix/consumer data
                               iso_names=c("d13C","d15N"), #tracers/isotopes
                               factors=c(NULL),            #factors
                               fac_random=c(NULL),         #random effects
                               fac_nested=c(NULL),         #Nested effects
                               cont_effects=NULL)          #continuous effects
#############################


#############################
# Load food source data
#############################
source.filename <- ("~/Isotopic_Routing_Small_Mammals/Mixing_Models/NAIN/Bartlett_diet_items.csv")

source_NAIN_Diet <- load_source_data(filename=source.filename, #name of the CSV file with food source data
                                     source_factors=NULL,     #column headings of random/fixed effects you have source data by
                                     conc_dep=TRUE,           #TRUE or FALSE for concentration dependence data 
                                     data_type="raw",         #"means" = means + sd or "raw" = raw data
                                     mix.NAIN.Diet)

# Load discrimination data
#No TDF for diet samples
Type<- c("AM Fungi","Arthropods","Red Maple", "EM Fungi", "Berries")
TDF_NAIN_Field <- as.data.frame(Type)
TDF_NAIN_Field$Meand13C<-0	
TDF_NAIN_Field$SDd13C<-0		
TDF_NAIN_Field$Meand15N<-0		
TDF_NAIN_Field$SDd15N<-0	
write.csv(TDF_NAIN_Field, "Discrimination_Diet.csv", row.names=FALSE)

Discrimination_Diet <- ("~/Isotopic_Routing_Small_Mammals/Mixing_Models/NAIN/Discrimination_Diet.csv")
discr_Diet <- load_discr_data(filename=Discrimination_Diet, mix.NAIN.Diet)
#############################


#############################
# Make isospace plot
#############################
plot_data(filename="NAIN_Diet_isospace", #name to save the isospace plot as
          plot_save_pdf=FALSE,       #TRUE or FALSE, should mix.SIAR save the plot as a .pdf?
          plot_save_png=TRUE,      #TRUE or FALSE, should mix.SIAR save the plot as a .png?
          mix.NAIN.Diet,source_NAIN_Diet,discr_Diet)

# If 2 isotopes/tracers, calculate normalized surface area of the convex hull polygon(s)
#   *Note 1: discrimination SD is added to the source SD (see calc_area.r for details)
if(mix.NAIN.Diet$n.iso==2) calc_area(source=source_NAIN_Diet,mix=mix.NAIN.Diet,discr=discr_Diet)
#############################


#############################
# Define prior
#############################
kw.alpha <- c(1,1,1,1,1)

kw.alpha <- kw.alpha*length(kw.alpha)/sum(kw.alpha)# Generate alpha hyperparameters scaling sum(alpha)=n.sources

kw.alpha[which(kw.alpha==0)] <- 0.01# the Dirichlet hyperparameters for the alpha.prior cannot be 0 (but can set = .01)

#plot using "plot_prior"
plot_prior(alpha.prior=kw.alpha,
           source=source_NAIN_Diet,
           plot_save_pdf=FALSE,
           plot_save_png=TRUE,
           filename="NAIN_Diet_prior_plot_kw_inf")
#   RED = your prior
#   DARK GREY = "uninformative"/generalist (alpha = 1)
#   LIGHT GREY = "uninformative" Jeffrey's prior (alpha = 1/n.sources)
#############################


#############################
# Write JAGS model file (define model structure)
##############################
# There are 3 error term options available:
#   1. Residual * Process (resid_err = TRUE, process_err = TRUE)
#   2. Residual only (resid_err = TRUE, process_err = FALSE)
#   3. Process only (resid_err = FALSE, process_err = TRUE)
NAIN_Diet_model <- "NAIN_Diet_model.txt"
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(NAIN_Diet_model, resid_err, process_err, mix.NAIN.Diet, source_NAIN_Diet)
##############################


##############################
# Run NAIN Diet model
##############################
# MCMC options:
#run <- list(chainLength=30000, burn=15000, thin=25, chains=3, calcDIC=TRUE)
costum_run <- list(chainLength=200000, burn=50000, thin=50, chains=3, calcDIC=TRUE)

#Model
NAIN_Diet<- run_model(run=costum_run,mix.NAIN.Diet,source_NAIN_Diet,discr_Diet,NAIN_Diet_model,
                      alpha.prior=kw.alpha, resid_err, process_err)
##############################


##############################
# Process JAGS output
##############################
# Choose output options (see ?output_options for details)
output_options <- list(summary_save = TRUE,                 
                       summary_name = "NAIN_Diet_summary_statistics", 
                       sup_post = FALSE,                    
                       plot_post_save_pdf = TRUE,           
                       plot_post_name = "NAIN_Diet_posterior_density",
                       sup_pairs = FALSE,             
                       plot_pairs_save_pdf = TRUE,    
                       plot_pairs_name = "NAIN_Diet_pairs_plot",
                       sup_xy = TRUE,           
                       plot_xy_save_pdf = FALSE,
                       plot_xy_name = "NAIN_Diet_xy_plot",
                       gelman = TRUE,
                       heidel = FALSE,  
                       geweke = TRUE,   
                       diag_save = TRUE,
                       diag_name = "NAIN_Diet_diagnostics",
                       indiv_effect = FALSE,       
                       plot_post_save_png = FALSE, 
                       plot_pairs_save_png = FALSE,
                       plot_xy_save_png = FALSE)
##############################


##############################
# Create diagnostics, summary statistics, and posterior plots
##############################
output_JAGS(NAIN_Diet, mix.NAIN.Diet, source_NAIN_Diet, output_options)
print(NAIN_Diet)
getwd()
save.image("NAIN_Diet_Jags.RData") 
##############################


##################################################
#Put data into a dataframe
##################################################
setwd("~/Isotopic_Routing_Small_Mammals/Mixing_Models/NAIN")

require(R2jags) # Load packages 

# Load data 
load("~/Isotopic_Routing_Small_Mammals/Mixing_Models/NAIN/NAIN_Diet_Jags.RData")
attach.jags(NAIN_Diet) # call model file
print(NAIN_Diet)

#assign names to each food source https://github.com/brianstock/MixSIAR/issues/30
NAIN_Diet_Post <- data.frame(Source = "Diet", AM_Fungi = p.global[,1], Arthropods = p.global[,2],
                             Berries = p.global[,3], EM_Fungi = p.global[,4], Red_Maple = p.global[,5], 
                             Fungi = p.global[,1]+p.global[,4])

NAIN_Diet_Post<-data.frame(Number=rownames(NAIN_Diet_Post), NAIN_Diet_Post)

require(tidyr)
NAIN_Diet_Posterior <- NAIN_Diet_Post %>% gather(Source,value,AM_Fungi:Fungi)#Use tidyr to gather posterior for each food source
NAIN_Diet_Posterior["Model"] <- "Diet"#Add model type
NAIN_Diet_Posterior["Species"] <- "NAIN"#Add Species
#save raw posterior data
write.csv(NAIN_Diet_Posterior, file = "NAIN_Diet_Posterior.csv", row.names=F)#save raw posterior data
write.csv(NAIN_Diet_Posterior, file = "~/Isotopic_Routing_Small_Mammals/Mixing_Models/Mixing_Model_Figure/NAIN_Diet_Posterior.csv", row.names=F)
################################################################################











################################################################################
#Field
################################################################################
setwd("~/Isotopic_Routing_Small_Mammals/Mixing_Models/NAIN")#set working directory
#############################
#load isotopic data
#############################
#Take stomach samples, convert to Field, and select species
Hair<- read.csv("Bartlett_isotope_hair_data.csv",header=T)
library(dplyr)
NAIN_Field<-filter(Hair, Species=="NAIN")
write.csv(NAIN_Field, "NAIN_Field.csv")

library(MixSIAR)
mix.filename.NAIN_Field <-("~/Isotopic_Routing_Small_Mammals/Mixing_Models/NAIN/NAIN_Field.csv") #load isotopic data for MixSIAR

mix.NAIN.Field <- load_mix_data(filename=mix.filename.NAIN_Field, #name of the CSV file with mix/consumer data
                                iso_names=c("d13C","d15N"), #tracers/isotopes
                                factors=c(NULL),            #factors
                                fac_random=c(NULL),         #random effects
                                fac_nested=c(NULL),         #Nested effects
                                cont_effects=NULL)          #continuous effects
#############################


#############################
# Load food source data
#############################
source.filename <- ("~/Isotopic_Routing_Small_Mammals/Mixing_Models/NAIN/Bartlett_diet_items.csv")

source_NAIN_Field <- load_source_data(filename=source.filename, #name of the CSV file with food source data
                                      source_factors=NULL,     #column headings of random/fixed effects you have source data by
                                      conc_dep=TRUE,           #TRUE or FALSE for concentration dependence data 
                                      data_type="raw",         #"means" = means + sd or "raw" = raw data
                                      mix.NAIN.Field)

# Load discrimination data
#Select species and method of estimated TDF
Est_DFT<- read.csv("Estimated_TDF.csv",header=T)#discrimination factors
NAIN_Field_d13C_TDF<-filter(Est_DFT,Species_Abr == "NAIN" & Method =="Field" & Isotope == "d13C")
NAIN_Field_d15N_TDF<-filter(Est_DFT,Species_Abr == "NAIN" & Method =="Field" & Isotope == "d15N")


Type<- c("AM Fungi","Arthropods","Red Maple", "EM Fungi", "Berries")
TDF_NAIN_Field <- as.data.frame(Type)
TDF_NAIN_Field$Meand13C<-NAIN_Field_d13C_TDF$TDF_Mean	
TDF_NAIN_Field$SDd13C<-NAIN_Field_d13C_TDF$TDF_SD			
TDF_NAIN_Field$Meand15N<-NAIN_Field_d15N_TDF$TDF_Mean			
TDF_NAIN_Field$SDd15N<-NAIN_Field_d15N_TDF$TDF_SD		
write.csv(TDF_NAIN_Field, "NAIN_Discrimination_Field.csv", row.names=FALSE)

NAIN_Discrimination_Field <- ("~/Isotopic_Routing_Small_Mammals/Mixing_Models/NAIN/NAIN_Discrimination_Field.csv")
NAIN_Discrimination_Field <- load_discr_data(filename=NAIN_Discrimination_Field, mix.NAIN.Field)
#############################


#############################
# Make isospace plot
#############################
plot_data(filename="NAIN_Field_isospace", #name to save the isospace plot as
          plot_save_pdf=FALSE,       #TRUE or FALSE, should mix.SIAR save the plot as a .pdf?
          plot_save_png=TRUE,      #TRUE or FALSE, should mix.SIAR save the plot as a .png?
          mix.NAIN.Field,source_NAIN_Field,NAIN_Discrimination_Field)

# If 2 isotopes/tracers, calculate normalized surface area of the convex hull polygon(s)
#   *Note 1: discrimination SD is added to the source SD (see calc_area.r for details)
if(mix.NAIN.Field$n.iso==2) calc_area(source=source_NAIN_Field,mix=mix.NAIN.Field,discr=NAIN_Discrimination_Field)
#############################


#############################
# Define prior
#############################
kw.alpha <- c(1,1,1,1,1)

kw.alpha <- kw.alpha*length(kw.alpha)/sum(kw.alpha)# Generate alpha hyperparameters scaling sum(alpha)=n.sources

kw.alpha[which(kw.alpha==0)] <- 0.01# the Dirichlet hyperparameters for the alpha.prior cannot be 0 (but can set = .01)

#plot using "plot_prior"
plot_prior(alpha.prior=kw.alpha,
           source=source_NAIN_Field,
           plot_save_pdf=FALSE,
           plot_save_png=TRUE,
           filename="NAIN_Field_prior_plot_kw_inf")
#   RED = your prior
#   DARK GREY = "uninformative"/generalist (alpha = 1)
#   LIGHT GREY = "uninformative" Jeffrey's prior (alpha = 1/n.sources)
#############################


#############################
# Write JAGS model file (define model structure)
##############################
# There are 3 error term options available:
#   1. Residual * Process (resid_err = TRUE, process_err = TRUE)
#   2. Residual only (resid_err = TRUE, process_err = FALSE)
#   3. Process only (resid_err = FALSE, process_err = TRUE)
NAIN_Field_model <- "NAIN_Field_model.txt"
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(NAIN_Field_model, resid_err, process_err, mix.NAIN.Field, source_NAIN_Field)
##############################


##############################
# Run NAIN Field model
##############################
# MCMC options:
#run <- list(chainLength=30000, burn=15000, thin=25, chains=3, calcDIC=TRUE)
costum_run <- list(chainLength=200000, burn=50000, thin=50, chains=3, calcDIC=TRUE)

#Model
NAIN_Field<- run_model(run=costum_run,mix.NAIN.Field,source_NAIN_Field,NAIN_Discrimination_Field,NAIN_Field_model,
                       alpha.prior=kw.alpha, resid_err, process_err)
##############################


##############################
# Process JAGS output
##############################
# Choose output options (see ?output_options for details)
output_options <- list(summary_save = TRUE,                 
                       summary_name = "NAIN_Field_summary_statistics", 
                       sup_post = FALSE,                    
                       plot_post_save_pdf = TRUE,           
                       plot_post_name = "NAIN_Field_posterior_density",
                       sup_pairs = FALSE,             
                       plot_pairs_save_pdf = TRUE,    
                       plot_pairs_name = "NAIN_Field_pairs_plot",
                       sup_xy = TRUE,           
                       plot_xy_save_pdf = FALSE,
                       plot_xy_name = "NAIN_Field_xy_plot",
                       gelman = TRUE,
                       heidel = FALSE,  
                       geweke = TRUE,   
                       diag_save = TRUE,
                       diag_name = "NAIN_Field_diagnostics",
                       indiv_effect = FALSE,       
                       plot_post_save_png = FALSE, 
                       plot_pairs_save_png = FALSE,
                       plot_xy_save_png = FALSE)
##############################


##############################
# Create diagnostics, summary statistics, and posterior plots
##############################
output_JAGS(NAIN_Field, mix.NAIN.Field, source_NAIN_Field, output_options)
print(NAIN_Field)
getwd()
save.image("NAIN_Field_Jags.RData") 
##############################


##################################################
#Put data into a dataframe
##################################################
setwd("~/Isotopic_Routing_Small_Mammals/Mixing_Models/NAIN")

require(R2jags) # Load packages 

# Load data 
load("~/Isotopic_Routing_Small_Mammals/Mixing_Models/NAIN/NAIN_Field_Jags.RData")
attach.jags(NAIN_Field) # call model file
print(NAIN_Field)

#assign names to each food source https://github.com/brianstock/MixSIAR/issues/30
NAIN_Field_Post <- data.frame(Source = "Field", AM_Fungi = p.global[,1], Arthropods = p.global[,2],
                              Berries = p.global[,3], EM_Fungi = p.global[,4], Red_Maple = p.global[,5], 
                              Fungi = p.global[,1]+p.global[,4]) 

NAIN_Field_Post<-data.frame(Number=rownames(NAIN_Field_Post), NAIN_Field_Post)#add posterior number
require(tidyr)
NAIN_Field_Posterior <- NAIN_Field_Post %>% gather(Source,value,AM_Fungi:Fungi)#Use tidyr to gather posterior for each food source
NAIN_Field_Posterior["Model"] <- "Field"#Add model type
NAIN_Field_Posterior["Species"] <- "NAIN"#Add Species

#save raw posterior data
write.csv(NAIN_Field_Posterior, file = "NAIN_Field_Posterior.csv", row.names=FALSE)#save raw posterior data
write.csv(NAIN_Field_Posterior, file = "~/Isotopic_Routing_Small_Mammals/Mixing_Models/Mixing_Model_Figure/NAIN_Field_Posterior.csv", row.names=FALSE)
################################################################################










################################################################################
#SIDER
################################################################################
setwd("~/Isotopic_Routing_Small_Mammals/Mixing_Models/NAIN")#set working directory
#############################
#load isotopic data
#############################
#Take stomach samples, convert to SIDER, and select species
Hair<- read.csv("Bartlett_isotope_hair_data.csv",header=T)
library(dplyr)
NAIN_SIDER<-filter(Hair, Species=="NAIN")
write.csv(NAIN_SIDER, "NAIN_SIDER.csv")

library(MixSIAR)
mix.filename.NAIN_SIDER <-("~/Isotopic_Routing_Small_Mammals/Mixing_Models/NAIN/NAIN_SIDER.csv") #load isotopic data for MixSIAR

mix.NAIN.SIDER <- load_mix_data(filename=mix.filename.NAIN_SIDER, #name of the CSV file with mix/consumer data
                                iso_names=c("d13C","d15N"), #tracers/isotopes
                                factors=c(NULL),            #factors
                                fac_random=c(NULL),         #random effects
                                fac_nested=c(NULL),         #Nested effects
                                cont_effects=NULL)          #continuous effects
#############################


#############################
# Load food source data
#############################
source.filename <- ("~/Isotopic_Routing_Small_Mammals/Mixing_Models/NAIN/Bartlett_diet_items.csv")

source_NAIN_SIDER <- load_source_data(filename=source.filename, #name of the CSV file with food source data
                                      source_factors=NULL,     #column headings of random/fixed effects you have source data by
                                      conc_dep=TRUE,           #TRUE or FALSE for concentration dependence data 
                                      data_type="raw",         #"means" = means + sd or "raw" = raw data
                                      mix.NAIN.SIDER)

# Load discrimination data
#Select species and method of estimated TDF
Est_DFT<- read.csv("Estimated_TDF.csv",header=T)#discrimination factors
NAIN_SIDER_d13C_TDF<-filter(Est_DFT,Species_Abr == "NAIN" & Method =="SIDER" & Isotope == "d13C")
NAIN_SIDER_d15N_TDF<-filter(Est_DFT,Species_Abr == "NAIN" & Method =="SIDER" & Isotope == "d15N")


Type<- c("AM Fungi","Arthropods","Red Maple", "EM Fungi", "Berries")
TDF_NAIN_SIDER <- as.data.frame(Type)
TDF_NAIN_SIDER$Meand13C<-NAIN_SIDER_d13C_TDF$TDF_Mean	
TDF_NAIN_SIDER$SDd13C<-NAIN_SIDER_d13C_TDF$TDF_SD			
TDF_NAIN_SIDER$Meand15N<-NAIN_SIDER_d15N_TDF$TDF_Mean			
TDF_NAIN_SIDER$SDd15N<-NAIN_SIDER_d15N_TDF$TDF_SD		
write.csv(TDF_NAIN_SIDER, "NAIN_Discrimination_SIDER.csv", row.names=FALSE)

NAIN_Discrimination_SIDER <- ("~/Isotopic_Routing_Small_Mammals/Mixing_Models/NAIN/NAIN_Discrimination_SIDER.csv")
NAIN_Discrimination_SIDER <- load_discr_data(filename=NAIN_Discrimination_SIDER, mix.NAIN.SIDER)
#############################


#############################
# Make isospace plot
#############################
plot_data(filename="NAIN_SIDER_isospace", #name to save the isospace plot as
          plot_save_pdf=FALSE,       #TRUE or FALSE, should mix.SIAR save the plot as a .pdf?
          plot_save_png=TRUE,      #TRUE or FALSE, should mix.SIAR save the plot as a .png?
          mix.NAIN.SIDER,source_NAIN_SIDER,NAIN_Discrimination_SIDER)

# If 2 isotopes/tracers, calculate normalized surface area of the convex hull polygon(s)
#   *Note 1: discrimination SD is added to the source SD (see calc_area.r for details)
if(mix.NAIN.SIDER$n.iso==2) calc_area(source=source_NAIN_SIDER,mix=mix.NAIN.SIDER,discr=NAIN_Discrimination_SIDER)
#############################


#############################
# Define prior
#############################
kw.alpha <- c(1,1,1,1,1)

kw.alpha <- kw.alpha*length(kw.alpha)/sum(kw.alpha)# Generate alpha hyperparameters scaling sum(alpha)=n.sources

kw.alpha[which(kw.alpha==0)] <- 0.01# the Dirichlet hyperparameters for the alpha.prior cannot be 0 (but can set = .01)

#plot using "plot_prior"
plot_prior(alpha.prior=kw.alpha,
           source=source_NAIN_SIDER,
           plot_save_pdf=FALSE,
           plot_save_png=TRUE,
           filename="NAIN_SIDER_prior_plot_kw_inf")
#   RED = your prior
#   DARK GREY = "uninformative"/generalist (alpha = 1)
#   LIGHT GREY = "uninformative" Jeffrey's prior (alpha = 1/n.sources)
#############################


#############################
# Write JAGS model file (define model structure)
##############################
# There are 3 error term options available:
#   1. Residual * Process (resid_err = TRUE, process_err = TRUE)
#   2. Residual only (resid_err = TRUE, process_err = FALSE)
#   3. Process only (resid_err = FALSE, process_err = TRUE)
NAIN_SIDER_model <- "NAIN_SIDER_model.txt"
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(NAIN_SIDER_model, resid_err, process_err, mix.NAIN.SIDER, source_NAIN_SIDER)
##############################


##############################
# Run NAIN SIDER model
##############################
# MCMC options:
#run <- list(chainLength=30000, burn=15000, thin=25, chains=3, calcDIC=TRUE)
costum_run <- list(chainLength=200000, burn=50000, thin=50, chains=3, calcDIC=TRUE)

#Model
NAIN_SIDER<- run_model(run=costum_run,mix.NAIN.SIDER,source_NAIN_SIDER,NAIN_Discrimination_SIDER,NAIN_SIDER_model,
                       alpha.prior=kw.alpha, resid_err, process_err)
##############################


##############################
# Process JAGS output
##############################
# Choose output options (see ?output_options for details)
output_options <- list(summary_save = TRUE,                 
                       summary_name = "NAIN_SIDER_summary_statistics", 
                       sup_post = FALSE,                    
                       plot_post_save_pdf = TRUE,           
                       plot_post_name = "NAIN_SIDER_posterior_density",
                       sup_pairs = FALSE,             
                       plot_pairs_save_pdf = TRUE,    
                       plot_pairs_name = "NAIN_SIDER_pairs_plot",
                       sup_xy = TRUE,           
                       plot_xy_save_pdf = FALSE,
                       plot_xy_name = "NAIN_SIDER_xy_plot",
                       gelman = TRUE,
                       heidel = FALSE,  
                       geweke = TRUE,   
                       diag_save = TRUE,
                       diag_name = "NAIN_SIDER_diagnostics",
                       indiv_effect = FALSE,       
                       plot_post_save_png = FALSE, 
                       plot_pairs_save_png = FALSE,
                       plot_xy_save_png = FALSE)
##############################


##############################
# Create diagnostics, summary statistics, and posterior plots
##############################
output_JAGS(NAIN_SIDER, mix.NAIN.SIDER, source_NAIN_SIDER, output_options)
print(NAIN_SIDER)
getwd()
save.image("NAIN_SIDER_Jags.RData") 
##############################


##################################################
#Put data into a dataframe
##################################################
setwd("~/Isotopic_Routing_Small_Mammals/Mixing_Models/NAIN")

require(R2jags) # Load packages 

# Load data 
load("~/Isotopic_Routing_Small_Mammals/Mixing_Models/NAIN/NAIN_SIDER_Jags.RData")
attach.jags(NAIN_SIDER) # call model file
print(NAIN_SIDER)

#assign names to each food source https://github.com/brianstock/MixSIAR/issues/30
NAIN_SIDER_Post <- data.frame(Source = "SIDER", AM_Fungi = p.global[,1], Arthropods = p.global[,2],
                              Berries = p.global[,3], EM_Fungi = p.global[,4], Red_Maple = p.global[,5], 
                              Fungi = p.global[,1]+p.global[,4]) 

NAIN_SIDER_Post<-data.frame(Number=rownames(NAIN_SIDER_Post), NAIN_SIDER_Post)#add posterior number
require(tidyr)
NAIN_SIDER_Posterior <- NAIN_SIDER_Post %>% gather(Source,value,AM_Fungi:Fungi)#Use tidyr to gather posterior for each food source
NAIN_SIDER_Posterior["Model"] <- "SIDER"#Add model type
NAIN_SIDER_Posterior["Species"] <- "NAIN"#Add Species

#save raw posterior data
write.csv(NAIN_SIDER_Posterior, file = "NAIN_SIDER_Posterior.csv", row.names=FALSE)#save raw posterior data
write.csv(NAIN_SIDER_Posterior, file = "~/Isotopic_Routing_Small_Mammals/Mixing_Models/Mixing_Model_Figure/NAIN_SIDER_Posterior.csv", row.names=FALSE)
################################################################################











################################################################################
#SIDER_update
################################################################################
setwd("~/Isotopic_Routing_Small_Mammals/Mixing_Models/NAIN")#set working directory
#############################
#load isotopic data
#############################
#Take stomach samples, convert to SIDER_update, and select species
Hair<- read.csv("Bartlett_isotope_hair_data.csv",header=T)
library(dplyr)
NAIN_SIDER_update<-filter(Hair, Species=="NAIN")
write.csv(NAIN_SIDER_update, "NAIN_SIDER_update.csv")

library(MixSIAR)
mix.filename.NAIN_SIDER_update <-("~/Isotopic_Routing_Small_Mammals/Mixing_Models/NAIN/NAIN_SIDER_update.csv") #load isotopic data for MixSIAR

mix.NAIN.SIDER_update <- load_mix_data(filename=mix.filename.NAIN_SIDER_update, #name of the CSV file with mix/consumer data
                                       iso_names=c("d13C","d15N"), #tracers/isotopes
                                       factors=c(NULL),            #factors
                                       fac_random=c(NULL),         #random effects
                                       fac_nested=c(NULL),         #Nested effects
                                       cont_effects=NULL)          #continuous effects
#############################


#############################
# Load food source data
#############################
source.filename <- ("~/Isotopic_Routing_Small_Mammals/Mixing_Models/NAIN/Bartlett_diet_items.csv")

source_NAIN_SIDER_update <- load_source_data(filename=source.filename, #name of the CSV file with food source data
                                             source_factors=NULL,     #column headings of random/fixed effects you have source data by
                                             conc_dep=TRUE,           #TRUE or FALSE for concentration dependence data 
                                             data_type="raw",         #"means" = means + sd or "raw" = raw data
                                             mix.NAIN.SIDER_update)

# Load discrimination data
#Select species and method of estimated TDF
Est_DFT<- read.csv("Estimated_TDF.csv",header=T)#discrimination factors
NAIN_SIDER_update_d13C_TDF<-filter(Est_DFT,Species_Abr == "NAIN" & Method =="SIDER_update" & Isotope == "d13C")
NAIN_SIDER_update_d15N_TDF<-filter(Est_DFT,Species_Abr == "NAIN" & Method =="SIDER_update" & Isotope == "d15N")


Type<- c("AM Fungi","Arthropods","Red Maple", "EM Fungi", "Berries")
TDF_NAIN_SIDER_update <- as.data.frame(Type)
TDF_NAIN_SIDER_update$Meand13C<-NAIN_SIDER_update_d13C_TDF$TDF_Mean	
TDF_NAIN_SIDER_update$SDd13C<-NAIN_SIDER_update_d13C_TDF$TDF_SD			
TDF_NAIN_SIDER_update$Meand15N<-NAIN_SIDER_update_d15N_TDF$TDF_Mean			
TDF_NAIN_SIDER_update$SDd15N<-NAIN_SIDER_update_d15N_TDF$TDF_SD		
write.csv(TDF_NAIN_SIDER_update, "NAIN_Discrimination_SIDER_update.csv", row.names=FALSE)

NAIN_Discrimination_SIDER_update <- ("~/Isotopic_Routing_Small_Mammals/Mixing_Models/NAIN/NAIN_Discrimination_SIDER_update.csv")
NAIN_Discrimination_SIDER_update <- load_discr_data(filename=NAIN_Discrimination_SIDER_update, mix.NAIN.SIDER_update)
#############################


#############################
# Make isospace plot
#############################
plot_data(filename="NAIN_SIDER_update_isospace", #name to save the isospace plot as
          plot_save_pdf=FALSE,       #TRUE or FALSE, should mix.SIAR save the plot as a .pdf?
          plot_save_png=TRUE,      #TRUE or FALSE, should mix.SIAR save the plot as a .png?
          mix.NAIN.SIDER_update,source_NAIN_SIDER_update,NAIN_Discrimination_SIDER_update)

# If 2 isotopes/tracers, calculate normalized surface area of the convex hull polygon(s)
#   *Note 1: discrimination SD is added to the source SD (see calc_area.r for details)
if(mix.NAIN.SIDER_update$n.iso==2) calc_area(source=source_NAIN_SIDER_update,mix=mix.NAIN.SIDER_update,discr=NAIN_Discrimination_SIDER_update)
#############################


#############################
# Define prior
#############################
kw.alpha <- c(1,1,1,1,1)

kw.alpha <- kw.alpha*length(kw.alpha)/sum(kw.alpha)# Generate alpha hyperparameters scaling sum(alpha)=n.sources

kw.alpha[which(kw.alpha==0)] <- 0.01# the Dirichlet hyperparameters for the alpha.prior cannot be 0 (but can set = .01)

#plot using "plot_prior"
plot_prior(alpha.prior=kw.alpha,
           source=source_NAIN_SIDER_update,
           plot_save_pdf=FALSE,
           plot_save_png=TRUE,
           filename="NAIN_SIDER_update_prior_plot_kw_inf")
#   RED = your prior
#   DARK GREY = "uninformative"/generalist (alpha = 1)
#   LIGHT GREY = "uninformative" Jeffrey's prior (alpha = 1/n.sources)
#############################


#############################
# Write JAGS model file (define model structure)
##############################
# There are 3 error term options available:
#   1. Residual * Process (resid_err = TRUE, process_err = TRUE)
#   2. Residual only (resid_err = TRUE, process_err = FALSE)
#   3. Process only (resid_err = FALSE, process_err = TRUE)
NAIN_SIDER_update_model <- "NAIN_SIDER_update_model.txt"
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(NAIN_SIDER_update_model, resid_err, process_err, mix.NAIN.SIDER_update, source_NAIN_SIDER_update)
##############################


##############################
# Run NAIN SIDER_update model
##############################
# MCMC options:
#run <- list(chainLength=30000, burn=15000, thin=25, chains=3, calcDIC=TRUE)

costum_run <- list(chainLength=200000, burn=50000, thin=50, chains=3, calcDIC=TRUE)

#Model
NAIN_SIDER_update<- run_model(run=costum_run,mix.NAIN.SIDER_update,source_NAIN_SIDER_update,NAIN_Discrimination_SIDER_update,NAIN_SIDER_update_model,
                              alpha.prior=kw.alpha, resid_err, process_err)
##############################


##############################
# Process JAGS output
##############################
# Choose output options (see ?output_options for details)
output_options <- list(summary_save = TRUE,                 
                       summary_name = "NAIN_SIDER_update_summary_statistics", 
                       sup_post = FALSE,                    
                       plot_post_save_pdf = TRUE,           
                       plot_post_name = "NAIN_SIDER_update_posterior_density",
                       sup_pairs = FALSE,             
                       plot_pairs_save_pdf = TRUE,    
                       plot_pairs_name = "NAIN_SIDER_update_pairs_plot",
                       sup_xy = TRUE,           
                       plot_xy_save_pdf = FALSE,
                       plot_xy_name = "NAIN_SIDER_update_xy_plot",
                       gelman = TRUE,
                       heidel = FALSE,  
                       geweke = TRUE,   
                       diag_save = TRUE,
                       diag_name = "NAIN_SIDER_update_diagnostics",
                       indiv_effect = FALSE,       
                       plot_post_save_png = FALSE, 
                       plot_pairs_save_png = FALSE,
                       plot_xy_save_png = FALSE)
##############################


##############################
# Create diagnostics, summary statistics, and posterior plots
##############################
output_JAGS(NAIN_SIDER_update, mix.NAIN.SIDER_update, source_NAIN_SIDER_update, output_options)
print(NAIN_SIDER_update)
getwd()
save.image("NAIN_SIDER_update_Jags.RData") 
##############################


##################################################
#Put data into a dataframe
##################################################
setwd("~/Isotopic_Routing_Small_Mammals/Mixing_Models/NAIN")

require(R2jags) # Load packages 

# Load data 
load("~/Isotopic_Routing_Small_Mammals/Mixing_Models/NAIN/NAIN_SIDER_update_Jags.RData")
attach.jags(NAIN_SIDER_update) # call model file
print(NAIN_SIDER_update)

#assign names to each food source https://github.com/brianstock/MixSIAR/issues/30
NAIN_SIDER_update_Post <- data.frame(Source = "SIDER_update", AM_Fungi = p.global[,1], Arthropods = p.global[,2],
                                     Berries = p.global[,3], EM_Fungi = p.global[,4], Red_Maple = p.global[,5], 
                                     Fungi = p.global[,1]+p.global[,4]) 

NAIN_SIDER_update_Post<-data.frame(Number=rownames(NAIN_SIDER_update_Post), NAIN_SIDER_update_Post)#add posterior number
require(tidyr)
NAIN_SIDER_update_Posterior <- NAIN_SIDER_update_Post %>% gather(Source,value,AM_Fungi:Fungi)#Use tidyr to gather posterior for each food source
NAIN_SIDER_update_Posterior["Model"] <- "SIDER_update"#Add model type
NAIN_SIDER_update_Posterior["Species"] <- "NAIN"#Add Species

#save raw posterior data
write.csv(NAIN_SIDER_update_Posterior, file = "NAIN_SIDER_update_Posterior.csv", row.names=FALSE)#save raw posterior data
write.csv(NAIN_SIDER_update_Posterior, file = "~/Isotopic_Routing_Small_Mammals/Mixing_Models/Mixing_Model_Figure/NAIN_SIDER_update_Posterior.csv", row.names=FALSE)
################################################################################









################################################################################
#Meta_analysis
################################################################################
setwd("~/Isotopic_Routing_Small_Mammals/Mixing_Models/NAIN")#set working directory
#############################
#load isotopic data
#############################
#Take stomach samples, convert to Meta_analysis, and select species
Hair<- read.csv("Bartlett_isotope_hair_data.csv",header=T)
library(dplyr)
NAIN_Meta_analysis<-filter(Hair, Species=="NAIN")
write.csv(NAIN_Meta_analysis, "NAIN_Meta_analysis.csv")

library(MixSIAR)
mix.filename.NAIN_Meta_analysis <-("~/Isotopic_Routing_Small_Mammals/Mixing_Models/NAIN/NAIN_Meta_analysis.csv") #load isotopic data for MixSIAR

mix.NAIN.Meta_analysis <- load_mix_data(filename=mix.filename.NAIN_Meta_analysis, #name of the CSV file with mix/consumer data
                                        iso_names=c("d13C","d15N"), #tracers/isotopes
                                        factors=c(NULL),            #factors
                                        fac_random=c(NULL),         #random effects
                                        fac_nested=c(NULL),         #Nested effects
                                        cont_effects=NULL)          #continuous effects
#############################


#############################
# Load food source data
#############################
source.filename <- ("~/Isotopic_Routing_Small_Mammals/Mixing_Models/NAIN/Bartlett_diet_items.csv")

source_NAIN_Meta_analysis <- load_source_data(filename=source.filename, #name of the CSV file with food source data
                                              source_factors=NULL,     #column headings of random/fixed effects you have source data by
                                              conc_dep=TRUE,           #TRUE or FALSE for concentration dependence data 
                                              data_type="raw",         #"means" = means + sd or "raw" = raw data
                                              mix.NAIN.Meta_analysis)

# Load discrimination data
#Select species and method of estimated TDF
Est_DFT<- read.csv("Estimated_TDF.csv",header=T)#discrimination factors
NAIN_Meta_analysis_d13C_TDF<-filter(Est_DFT,Species_Abr == "NAIN" & Method =="Meta_analysis" & Isotope == "d13C")
NAIN_Meta_analysis_d15N_TDF<-filter(Est_DFT,Species_Abr == "NAIN" & Method =="Meta_analysis" & Isotope == "d15N")


Type<- c("AM Fungi","Arthropods","Red Maple", "EM Fungi", "Berries")
TDF_NAIN_Meta_analysis <- as.data.frame(Type)
TDF_NAIN_Meta_analysis$Meand13C<-NAIN_Meta_analysis_d13C_TDF$TDF_Mean	
TDF_NAIN_Meta_analysis$SDd13C<-NAIN_Meta_analysis_d13C_TDF$TDF_SD			
TDF_NAIN_Meta_analysis$Meand15N<-NAIN_Meta_analysis_d15N_TDF$TDF_Mean			
TDF_NAIN_Meta_analysis$SDd15N<-NAIN_Meta_analysis_d15N_TDF$TDF_SD		
write.csv(TDF_NAIN_Meta_analysis, "NAIN_Discrimination_Meta_analysis.csv", row.names=FALSE)

NAIN_Discrimination_Meta_analysis <- ("~/Isotopic_Routing_Small_Mammals/Mixing_Models/NAIN/NAIN_Discrimination_Meta_analysis.csv")
NAIN_Discrimination_Meta_analysis <- load_discr_data(filename=NAIN_Discrimination_Meta_analysis, mix.NAIN.Meta_analysis)
#############################


#############################
# Make isospace plot
#############################
plot_data(filename="NAIN_Meta_analysis_isospace", #name to save the isospace plot as
          plot_save_pdf=FALSE,       #TRUE or FALSE, should mix.SIAR save the plot as a .pdf?
          plot_save_png=TRUE,      #TRUE or FALSE, should mix.SIAR save the plot as a .png?
          mix.NAIN.Meta_analysis,source_NAIN_Meta_analysis,NAIN_Discrimination_Meta_analysis)

# If 2 isotopes/tracers, calculate normalized surface area of the convex hull polygon(s)
#   *Note 1: discrimination SD is added to the source SD (see calc_area.r for details)
if(mix.NAIN.Meta_analysis$n.iso==2) calc_area(source=source_NAIN_Meta_analysis,mix=mix.NAIN.Meta_analysis,discr=NAIN_Discrimination_Meta_analysis)
#############################


#############################
# Define prior
#############################
kw.alpha <- c(1,1,1,1,1)

kw.alpha <- kw.alpha*length(kw.alpha)/sum(kw.alpha)# Generate alpha hyperparameters scaling sum(alpha)=n.sources

kw.alpha[which(kw.alpha==0)] <- 0.01# the Dirichlet hyperparameters for the alpha.prior cannot be 0 (but can set = .01)

#plot using "plot_prior"
plot_prior(alpha.prior=kw.alpha,
           source=source_NAIN_Meta_analysis,
           plot_save_pdf=FALSE,
           plot_save_png=TRUE,
           filename="NAIN_Meta_analysis_prior_plot_kw_inf")
#   RED = your prior
#   DARK GREY = "uninformative"/generalist (alpha = 1)
#   LIGHT GREY = "uninformative" Jeffrey's prior (alpha = 1/n.sources)
#############################


#############################
# Write JAGS model file (define model structure)
##############################
# There are 3 error term options available:
#   1. Residual * Process (resid_err = TRUE, process_err = TRUE)
#   2. Residual only (resid_err = TRUE, process_err = FALSE)
#   3. Process only (resid_err = FALSE, process_err = TRUE)
NAIN_Meta_analysis_model <- "NAIN_Meta_analysis_model.txt"
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(NAIN_Meta_analysis_model, resid_err, process_err, mix.NAIN.Meta_analysis, source_NAIN_Meta_analysis)
##############################


##############################
# Run NAIN Meta_analysis model
##############################
# MCMC options:
#run <- list(chainLength=30000, burn=15000, thin=25, chains=3, calcDIC=TRUE)
costum_run <- list(chainLength=200000, burn=50000, thin=50, chains=3, calcDIC=TRUE)

#Model
NAIN_Meta_analysis<- run_model(run=costum_run,mix.NAIN.Meta_analysis,source_NAIN_Meta_analysis,NAIN_Discrimination_Meta_analysis,NAIN_Meta_analysis_model,
                               alpha.prior=kw.alpha, resid_err, process_err)
##############################


##############################
# Process JAGS output
##############################
# Choose output options (see ?output_options for details)
output_options <- list(summary_save = TRUE,                 
                       summary_name = "NAIN_Meta_analysis_summary_statistics", 
                       sup_post = FALSE,                    
                       plot_post_save_pdf = TRUE,           
                       plot_post_name = "NAIN_Meta_analysis_posterior_density",
                       sup_pairs = FALSE,             
                       plot_pairs_save_pdf = TRUE,    
                       plot_pairs_name = "NAIN_Meta_analysis_pairs_plot",
                       sup_xy = TRUE,           
                       plot_xy_save_pdf = FALSE,
                       plot_xy_name = "NAIN_Meta_analysis_xy_plot",
                       gelman = TRUE,
                       heidel = FALSE,  
                       geweke = TRUE,   
                       diag_save = TRUE,
                       diag_name = "NAIN_Meta_analysis_diagnostics",
                       indiv_effect = FALSE,       
                       plot_post_save_png = FALSE, 
                       plot_pairs_save_png = FALSE,
                       plot_xy_save_png = FALSE)
##############################


##############################
# Create diagnostics, summary statistics, and posterior plots
##############################
output_JAGS(NAIN_Meta_analysis, mix.NAIN.Meta_analysis, source_NAIN_Meta_analysis, output_options)
print(NAIN_Meta_analysis)
getwd()
save.image("NAIN_Meta_analysis_Jags.RData") 
##############################


##################################################
#Put data into a dataframe
##################################################
setwd("~/Isotopic_Routing_Small_Mammals/Mixing_Models/NAIN")

require(R2jags) # Load packages 

# Load data 
load("~/Isotopic_Routing_Small_Mammals/Mixing_Models/NAIN/NAIN_Meta_analysis_Jags.RData")
attach.jags(NAIN_Meta_analysis) # call model file
print(NAIN_Meta_analysis)

#assign names to each food source https://github.com/brianstock/MixSIAR/issues/30
NAIN_Meta_analysis_Post <- data.frame(Source = "Meta_analysis", AM_Fungi = p.global[,1], Arthropods = p.global[,2],
                                      Berries = p.global[,3], EM_Fungi = p.global[,4], Red_Maple = p.global[,5], 
                                      Fungi = p.global[,1]+p.global[,4]) 

NAIN_Meta_analysis_Post<-data.frame(Number=rownames(NAIN_Meta_analysis_Post), NAIN_Meta_analysis_Post)#add posterior number

require(tidyr)
NAIN_Meta_analysis_Posterior <- NAIN_Meta_analysis_Post %>% gather(Source,value,AM_Fungi:Fungi)#Use tidyr to gather posterior for each food source
NAIN_Meta_analysis_Posterior["Model"] <- "Meta_analysis"#Add model type
NAIN_Meta_analysis_Posterior["Species"] <- "NAIN"#Add Species

#save raw posterior data
write.csv(NAIN_Meta_analysis_Posterior, file = "NAIN_Meta_analysis_Posterior.csv", row.names=FALSE)#save raw posterior data
write.csv(NAIN_Meta_analysis_Posterior, file = "~/Isotopic_Routing_Small_Mammals/Mixing_Models/Mixing_Model_Figure/NAIN_Meta_analysis_Posterior.csv", row.names=FALSE)
################################################################################










################################################################################
#Peromyscus maniculatus
################################################################################
rm(list=ls())# clears workspace
setwd("~/Isotopic_Routing_Small_Mammals/Mixing_Models/PEMA")#set working directory


################################################################################
#Diet
################################################################################

#############################
#load isotopic data
#############################
#Take stomach samples, convert to diet, and select species
library(dplyr)
Isotopes<- read.csv("Field_study_data.csv",header=T)
Diet<-filter(Isotopes, Type =="Stomach")#select stomach data
Diet<-Diet %>% mutate(d15N=ifelse(Species == "MYGA",d15N-1.85,d15N-1))#subtract 1 to account digestive enzymes of stomach (1.85 for MGYA)
Diet<-mutate(Diet, Type = "Diet")#code 'type' as diet
PEMA_Diet<- filter(Diet,Species =="PEMA") #Filter out to desired species

write.csv(PEMA_Diet, "PEMA_Diet.csv")

library(MixSIAR)
mix.filename.PEMA_Diet <-("~/Isotopic_Routing_Small_Mammals/Mixing_Models/PEMA/PEMA_Diet.csv") #load isotopic data for MixSIAR

mix.PEMA.Diet <- load_mix_data(filename=mix.filename.PEMA_Diet, #name of the CSV file with mix/consumer data
                               iso_names=c("d13C","d15N"), #tracers/isotopes
                               factors=c(NULL),            #factors
                               fac_random=c(NULL),         #random effects
                               fac_nested=c(NULL),         #Nested effects
                               cont_effects=NULL)          #continuous effects
#############################


#############################
# Load food source data
#############################
source.filename <- ("~/Isotopic_Routing_Small_Mammals/Mixing_Models/PEMA/Bartlett_diet_items.csv")

source_PEMA_Diet <- load_source_data(filename=source.filename, #name of the CSV file with food source data
                                     source_factors=NULL,     #column headings of random/fixed effects you have source data by
                                     conc_dep=TRUE,           #TRUE or FALSE for concentration dependence data 
                                     data_type="raw",         #"means" = means + sd or "raw" = raw data
                                     mix.PEMA.Diet)

# Load discrimination data
#No TDF for diet samples
Type<- c("AM Fungi","Arthropods","Red Maple", "EM Fungi", "Berries")
TDF_PEMA_Field <- as.data.frame(Type)
TDF_PEMA_Field$Meand13C<-0	
TDF_PEMA_Field$SDd13C<-0		
TDF_PEMA_Field$Meand15N<-0		
TDF_PEMA_Field$SDd15N<-0	
write.csv(TDF_PEMA_Field, "Discrimination_Diet.csv", row.names=FALSE)

Discrimination_Diet <- ("~/Isotopic_Routing_Small_Mammals/Mixing_Models/PEMA/Discrimination_Diet.csv")
discr_Diet <- load_discr_data(filename=Discrimination_Diet, mix.PEMA.Diet)
#############################


#############################
# Make isospace plot
#############################
plot_data(filename="PEMA_Diet_isospace", #name to save the isospace plot as
          plot_save_pdf=FALSE,       #TRUE or FALSE, should mix.SIAR save the plot as a .pdf?
          plot_save_png=TRUE,      #TRUE or FALSE, should mix.SIAR save the plot as a .png?
          mix.PEMA.Diet,source_PEMA_Diet,discr_Diet)

# If 2 isotopes/tracers, calculate normalized surface area of the convex hull polygon(s)
#   *Note 1: discrimination SD is added to the source SD (see calc_area.r for details)
if(mix.PEMA.Diet$n.iso==2) calc_area(source=source_PEMA_Diet,mix=mix.PEMA.Diet,discr=discr_Diet)
#############################


#############################
# Define prior
#############################
kw.alpha <- c(1,1,1,1,1)

kw.alpha <- kw.alpha*length(kw.alpha)/sum(kw.alpha)# Generate alpha hyperparameters scaling sum(alpha)=n.sources

kw.alpha[which(kw.alpha==0)] <- 0.01# the Dirichlet hyperparameters for the alpha.prior cannot be 0 (but can set = .01)

#plot using "plot_prior"
plot_prior(alpha.prior=kw.alpha,
           source=source_PEMA_Diet,
           plot_save_pdf=FALSE,
           plot_save_png=TRUE,
           filename="PEMA_Diet_prior_plot_kw_inf")
#   RED = your prior
#   DARK GREY = "uninformative"/generalist (alpha = 1)
#   LIGHT GREY = "uninformative" Jeffrey's prior (alpha = 1/n.sources)
#############################


#############################
# Write JAGS model file (define model structure)
##############################
# There are 3 error term options available:
#   1. Residual * Process (resid_err = TRUE, process_err = TRUE)
#   2. Residual only (resid_err = TRUE, process_err = FALSE)
#   3. Process only (resid_err = FALSE, process_err = TRUE)
PEMA_Diet_model <- "PEMA_Diet_model.txt"
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(PEMA_Diet_model, resid_err, process_err, mix.PEMA.Diet, source_PEMA_Diet)
##############################


##############################
# Run PEMA Diet model
##############################
# MCMC options:
#run <- list(chainLength=30000, burn=15000, thin=25, chains=3, calcDIC=TRUE)
costum_run <- list(chainLength=200000, burn=50000, thin=50, chains=3, calcDIC=TRUE)

#Model
PEMA_Diet<- run_model(run=costum_run,mix.PEMA.Diet,source_PEMA_Diet,discr_Diet,PEMA_Diet_model,
                      alpha.prior=kw.alpha, resid_err, process_err)
##############################


##############################
# Process JAGS output
##############################
# Choose output options (see ?output_options for details)
output_options <- list(summary_save = TRUE,                 
                       summary_name = "PEMA_Diet_summary_statistics", 
                       sup_post = FALSE,                    
                       plot_post_save_pdf = TRUE,           
                       plot_post_name = "PEMA_Diet_posterior_density",
                       sup_pairs = FALSE,             
                       plot_pairs_save_pdf = TRUE,    
                       plot_pairs_name = "PEMA_Diet_pairs_plot",
                       sup_xy = TRUE,           
                       plot_xy_save_pdf = FALSE,
                       plot_xy_name = "PEMA_Diet_xy_plot",
                       gelman = TRUE,
                       heidel = FALSE,  
                       geweke = TRUE,   
                       diag_save = TRUE,
                       diag_name = "PEMA_Diet_diagnostics",
                       indiv_effect = FALSE,       
                       plot_post_save_png = FALSE, 
                       plot_pairs_save_png = FALSE,
                       plot_xy_save_png = FALSE)
##############################


##############################
# Create diagnostics, summary statistics, and posterior plots
##############################
output_JAGS(PEMA_Diet, mix.PEMA.Diet, source_PEMA_Diet, output_options)
print(PEMA_Diet)
getwd()
save.image("PEMA_Diet_Jags.RData") 
##############################


##################################################
#Put data into a dataframe
##################################################
setwd("~/Isotopic_Routing_Small_Mammals/Mixing_Models/PEMA")

require(R2jags) # Load packages 

# Load data 
load("~/Isotopic_Routing_Small_Mammals/Mixing_Models/PEMA/PEMA_Diet_Jags.RData")
attach.jags(PEMA_Diet) # call model file
print(PEMA_Diet)

#assign names to each food source https://github.com/brianstock/MixSIAR/issues/30
PEMA_Diet_Post <- data.frame(Source = "Diet", AM_Fungi = p.global[,1], Arthropods = p.global[,2],
                             Berries = p.global[,3], EM_Fungi = p.global[,4], Red_Maple = p.global[,5], 
                             Fungi = p.global[,1]+p.global[,4])

PEMA_Diet_Post<-data.frame(Number=rownames(PEMA_Diet_Post), PEMA_Diet_Post)

require(tidyr)
PEMA_Diet_Posterior <- PEMA_Diet_Post %>% gather(Source,value,AM_Fungi:Fungi)#Use tidyr to gather posterior for each food source
PEMA_Diet_Posterior["Model"] <- "Diet"#Add model type
PEMA_Diet_Posterior["Species"] <- "PEMA"#Add Species
#save raw posterior data
write.csv(PEMA_Diet_Posterior, file = "PEMA_Diet_Posterior.csv", row.names=F)#save raw posterior data
write.csv(PEMA_Diet_Posterior, file = "~/Isotopic_Routing_Small_Mammals/Mixing_Models/Mixing_Model_Figure/PEMA_Diet_Posterior.csv", row.names=F)
################################################################################











################################################################################
#Field
################################################################################
setwd("~/Isotopic_Routing_Small_Mammals/Mixing_Models/PEMA")#set working directory
#############################
#load isotopic data
#############################
#Take stomach samples, convert to Field, and select species
Hair<- read.csv("Bartlett_isotope_hair_data.csv",header=T)
library(dplyr)
PEMA_Field<-filter(Hair, Species=="PEMA")
write.csv(PEMA_Field, "PEMA_Field.csv")

library(MixSIAR)
mix.filename.PEMA_Field <-("~/Isotopic_Routing_Small_Mammals/Mixing_Models/PEMA/PEMA_Field.csv") #load isotopic data for MixSIAR

mix.PEMA.Field <- load_mix_data(filename=mix.filename.PEMA_Field, #name of the CSV file with mix/consumer data
                                iso_names=c("d13C","d15N"), #tracers/isotopes
                                factors=c(NULL),            #factors
                                fac_random=c(NULL),         #random effects
                                fac_nested=c(NULL),         #Nested effects
                                cont_effects=NULL)          #continuous effects
#############################


#############################
# Load food source data
#############################
source.filename <- ("~/Isotopic_Routing_Small_Mammals/Mixing_Models/PEMA/Bartlett_diet_items.csv")

source_PEMA_Field <- load_source_data(filename=source.filename, #name of the CSV file with food source data
                                      source_factors=NULL,     #column headings of random/fixed effects you have source data by
                                      conc_dep=TRUE,           #TRUE or FALSE for concentration dependence data 
                                      data_type="raw",         #"means" = means + sd or "raw" = raw data
                                      mix.PEMA.Field)

# Load discrimination data
#Select species and method of estimated TDF
Est_DFT<- read.csv("Estimated_TDF.csv",header=T)#discrimination factors
PEMA_Field_d13C_TDF<-filter(Est_DFT,Species_Abr == "PEMA" & Method =="Field" & Isotope == "d13C")
PEMA_Field_d15N_TDF<-filter(Est_DFT,Species_Abr == "PEMA" & Method =="Field" & Isotope == "d15N")


Type<- c("AM Fungi","Arthropods","Red Maple", "EM Fungi", "Berries")
TDF_PEMA_Field <- as.data.frame(Type)
TDF_PEMA_Field$Meand13C<-PEMA_Field_d13C_TDF$TDF_Mean	
TDF_PEMA_Field$SDd13C<-PEMA_Field_d13C_TDF$TDF_SD			
TDF_PEMA_Field$Meand15N<-PEMA_Field_d15N_TDF$TDF_Mean			
TDF_PEMA_Field$SDd15N<-PEMA_Field_d15N_TDF$TDF_SD		
write.csv(TDF_PEMA_Field, "PEMA_Discrimination_Field.csv", row.names=FALSE)

PEMA_Discrimination_Field <- ("~/Isotopic_Routing_Small_Mammals/Mixing_Models/PEMA/PEMA_Discrimination_Field.csv")
PEMA_Discrimination_Field <- load_discr_data(filename=PEMA_Discrimination_Field, mix.PEMA.Field)
#############################


#############################
# Make isospace plot
#############################
plot_data(filename="PEMA_Field_isospace", #name to save the isospace plot as
          plot_save_pdf=FALSE,       #TRUE or FALSE, should mix.SIAR save the plot as a .pdf?
          plot_save_png=TRUE,      #TRUE or FALSE, should mix.SIAR save the plot as a .png?
          mix.PEMA.Field,source_PEMA_Field,PEMA_Discrimination_Field)

# If 2 isotopes/tracers, calculate normalized surface area of the convex hull polygon(s)
#   *Note 1: discrimination SD is added to the source SD (see calc_area.r for details)
if(mix.PEMA.Field$n.iso==2) calc_area(source=source_PEMA_Field,mix=mix.PEMA.Field,discr=PEMA_Discrimination_Field)
#############################


#############################
# Define prior
#############################
kw.alpha <- c(1,1,1,1,1)

kw.alpha <- kw.alpha*length(kw.alpha)/sum(kw.alpha)# Generate alpha hyperparameters scaling sum(alpha)=n.sources

kw.alpha[which(kw.alpha==0)] <- 0.01# the Dirichlet hyperparameters for the alpha.prior cannot be 0 (but can set = .01)

#plot using "plot_prior"
plot_prior(alpha.prior=kw.alpha,
           source=source_PEMA_Field,
           plot_save_pdf=FALSE,
           plot_save_png=TRUE,
           filename="PEMA_Field_prior_plot_kw_inf")
#   RED = your prior
#   DARK GREY = "uninformative"/generalist (alpha = 1)
#   LIGHT GREY = "uninformative" Jeffrey's prior (alpha = 1/n.sources)
#############################


#############################
# Write JAGS model file (define model structure)
##############################
# There are 3 error term options available:
#   1. Residual * Process (resid_err = TRUE, process_err = TRUE)
#   2. Residual only (resid_err = TRUE, process_err = FALSE)
#   3. Process only (resid_err = FALSE, process_err = TRUE)
PEMA_Field_model <- "PEMA_Field_model.txt"
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(PEMA_Field_model, resid_err, process_err, mix.PEMA.Field, source_PEMA_Field)
##############################


##############################
# Run PEMA Field model
##############################
# MCMC options:
#run <- list(chainLength=30000, burn=15000, thin=25, chains=3, calcDIC=TRUE)
costum_run <- list(chainLength=200000, burn=50000, thin=50, chains=3, calcDIC=TRUE)

#Model
PEMA_Field<- run_model(run=costum_run,mix.PEMA.Field,source_PEMA_Field,PEMA_Discrimination_Field,PEMA_Field_model,
                       alpha.prior=kw.alpha, resid_err, process_err)
##############################


##############################
# Process JAGS output
##############################
# Choose output options (see ?output_options for details)
output_options <- list(summary_save = TRUE,                 
                       summary_name = "PEMA_Field_summary_statistics", 
                       sup_post = FALSE,                    
                       plot_post_save_pdf = TRUE,           
                       plot_post_name = "PEMA_Field_posterior_density",
                       sup_pairs = FALSE,             
                       plot_pairs_save_pdf = TRUE,    
                       plot_pairs_name = "PEMA_Field_pairs_plot",
                       sup_xy = TRUE,           
                       plot_xy_save_pdf = FALSE,
                       plot_xy_name = "PEMA_Field_xy_plot",
                       gelman = TRUE,
                       heidel = FALSE,  
                       geweke = TRUE,   
                       diag_save = TRUE,
                       diag_name = "PEMA_Field_diagnostics",
                       indiv_effect = FALSE,       
                       plot_post_save_png = FALSE, 
                       plot_pairs_save_png = FALSE,
                       plot_xy_save_png = FALSE)
##############################


##############################
# Create diagnostics, summary statistics, and posterior plots
##############################
output_JAGS(PEMA_Field, mix.PEMA.Field, source_PEMA_Field, output_options)
print(PEMA_Field)
getwd()
save.image("PEMA_Field_Jags.RData") 
##############################


##################################################
#Put data into a dataframe
##################################################
setwd("~/Isotopic_Routing_Small_Mammals/Mixing_Models/PEMA")

require(R2jags) # Load packages 

# Load data 
load("~/Isotopic_Routing_Small_Mammals/Mixing_Models/PEMA/PEMA_Field_Jags.RData")
attach.jags(PEMA_Field) # call model file
print(PEMA_Field)

#assign names to each food source https://github.com/brianstock/MixSIAR/issues/30
PEMA_Field_Post <- data.frame(Source = "Field", AM_Fungi = p.global[,1], Arthropods = p.global[,2],
                              Berries = p.global[,3], EM_Fungi = p.global[,4], Red_Maple = p.global[,5], 
                              Fungi = p.global[,1]+p.global[,4]) 

PEMA_Field_Post<-data.frame(Number=rownames(PEMA_Field_Post), PEMA_Field_Post)#add posterior number
require(tidyr)
PEMA_Field_Posterior <- PEMA_Field_Post %>% gather(Source,value,AM_Fungi:Fungi)#Use tidyr to gather posterior for each food source
PEMA_Field_Posterior["Model"] <- "Field"#Add model type
PEMA_Field_Posterior["Species"] <- "PEMA"#Add Species

#save raw posterior data
write.csv(PEMA_Field_Posterior, file = "PEMA_Field_Posterior.csv", row.names=FALSE)#save raw posterior data
write.csv(PEMA_Field_Posterior, file = "~/Isotopic_Routing_Small_Mammals/Mixing_Models/Mixing_Model_Figure/PEMA_Field_Posterior.csv", row.names=FALSE)
################################################################################










################################################################################
#SIDER
################################################################################
setwd("~/Isotopic_Routing_Small_Mammals/Mixing_Models/PEMA")#set working directory
#############################
#load isotopic data
#############################
#Take stomach samples, convert to SIDER, and select species
Hair<- read.csv("Bartlett_isotope_hair_data.csv",header=T)
library(dplyr)
PEMA_SIDER<-filter(Hair, Species=="PEMA")
write.csv(PEMA_SIDER, "PEMA_SIDER.csv")

library(MixSIAR)
mix.filename.PEMA_SIDER <-("~/Isotopic_Routing_Small_Mammals/Mixing_Models/PEMA/PEMA_SIDER.csv") #load isotopic data for MixSIAR

mix.PEMA.SIDER <- load_mix_data(filename=mix.filename.PEMA_SIDER, #name of the CSV file with mix/consumer data
                                iso_names=c("d13C","d15N"), #tracers/isotopes
                                factors=c(NULL),            #factors
                                fac_random=c(NULL),         #random effects
                                fac_nested=c(NULL),         #Nested effects
                                cont_effects=NULL)          #continuous effects
#############################


#############################
# Load food source data
#############################
source.filename <- ("~/Isotopic_Routing_Small_Mammals/Mixing_Models/PEMA/Bartlett_diet_items.csv")

source_PEMA_SIDER <- load_source_data(filename=source.filename, #name of the CSV file with food source data
                                      source_factors=NULL,     #column headings of random/fixed effects you have source data by
                                      conc_dep=TRUE,           #TRUE or FALSE for concentration dependence data 
                                      data_type="raw",         #"means" = means + sd or "raw" = raw data
                                      mix.PEMA.SIDER)

# Load discrimination data
#Select species and method of estimated TDF
Est_DFT<- read.csv("Estimated_TDF.csv",header=T)#discrimination factors
PEMA_SIDER_d13C_TDF<-filter(Est_DFT,Species_Abr == "PEMA" & Method =="SIDER" & Isotope == "d13C")
PEMA_SIDER_d15N_TDF<-filter(Est_DFT,Species_Abr == "PEMA" & Method =="SIDER" & Isotope == "d15N")


Type<- c("AM Fungi","Arthropods","Red Maple", "EM Fungi", "Berries")
TDF_PEMA_SIDER <- as.data.frame(Type)
TDF_PEMA_SIDER$Meand13C<-PEMA_SIDER_d13C_TDF$TDF_Mean	
TDF_PEMA_SIDER$SDd13C<-PEMA_SIDER_d13C_TDF$TDF_SD			
TDF_PEMA_SIDER$Meand15N<-PEMA_SIDER_d15N_TDF$TDF_Mean			
TDF_PEMA_SIDER$SDd15N<-PEMA_SIDER_d15N_TDF$TDF_SD		
write.csv(TDF_PEMA_SIDER, "PEMA_Discrimination_SIDER.csv", row.names=FALSE)

PEMA_Discrimination_SIDER <- ("~/Isotopic_Routing_Small_Mammals/Mixing_Models/PEMA/PEMA_Discrimination_SIDER.csv")
PEMA_Discrimination_SIDER <- load_discr_data(filename=PEMA_Discrimination_SIDER, mix.PEMA.SIDER)
#############################


#############################
# Make isospace plot
#############################
plot_data(filename="PEMA_SIDER_isospace", #name to save the isospace plot as
          plot_save_pdf=FALSE,       #TRUE or FALSE, should mix.SIAR save the plot as a .pdf?
          plot_save_png=TRUE,      #TRUE or FALSE, should mix.SIAR save the plot as a .png?
          mix.PEMA.SIDER,source_PEMA_SIDER,PEMA_Discrimination_SIDER)

# If 2 isotopes/tracers, calculate normalized surface area of the convex hull polygon(s)
#   *Note 1: discrimination SD is added to the source SD (see calc_area.r for details)
if(mix.PEMA.SIDER$n.iso==2) calc_area(source=source_PEMA_SIDER,mix=mix.PEMA.SIDER,discr=PEMA_Discrimination_SIDER)
#############################


#############################
# Define prior
#############################
kw.alpha <- c(1,1,1,1,1)

kw.alpha <- kw.alpha*length(kw.alpha)/sum(kw.alpha)# Generate alpha hyperparameters scaling sum(alpha)=n.sources

kw.alpha[which(kw.alpha==0)] <- 0.01# the Dirichlet hyperparameters for the alpha.prior cannot be 0 (but can set = .01)

#plot using "plot_prior"
plot_prior(alpha.prior=kw.alpha,
           source=source_PEMA_SIDER,
           plot_save_pdf=FALSE,
           plot_save_png=TRUE,
           filename="PEMA_SIDER_prior_plot_kw_inf")
#   RED = your prior
#   DARK GREY = "uninformative"/generalist (alpha = 1)
#   LIGHT GREY = "uninformative" Jeffrey's prior (alpha = 1/n.sources)
#############################


#############################
# Write JAGS model file (define model structure)
##############################
# There are 3 error term options available:
#   1. Residual * Process (resid_err = TRUE, process_err = TRUE)
#   2. Residual only (resid_err = TRUE, process_err = FALSE)
#   3. Process only (resid_err = FALSE, process_err = TRUE)
PEMA_SIDER_model <- "PEMA_SIDER_model.txt"
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(PEMA_SIDER_model, resid_err, process_err, mix.PEMA.SIDER, source_PEMA_SIDER)
##############################


##############################
# Run PEMA SIDER model
##############################
# MCMC options:
#run <- list(chainLength=30000, burn=15000, thin=25, chains=3, calcDIC=TRUE)
costum_run <- list(chainLength=200000, burn=50000, thin=50, chains=3, calcDIC=TRUE)

#Model
PEMA_SIDER<- run_model(run=costum_run,mix.PEMA.SIDER,source_PEMA_SIDER,PEMA_Discrimination_SIDER,PEMA_SIDER_model,
                       alpha.prior=kw.alpha, resid_err, process_err)
##############################


##############################
# Process JAGS output
##############################
# Choose output options (see ?output_options for details)
output_options <- list(summary_save = TRUE,                 
                       summary_name = "PEMA_SIDER_summary_statistics", 
                       sup_post = FALSE,                    
                       plot_post_save_pdf = TRUE,           
                       plot_post_name = "PEMA_SIDER_posterior_density",
                       sup_pairs = FALSE,             
                       plot_pairs_save_pdf = TRUE,    
                       plot_pairs_name = "PEMA_SIDER_pairs_plot",
                       sup_xy = TRUE,           
                       plot_xy_save_pdf = FALSE,
                       plot_xy_name = "PEMA_SIDER_xy_plot",
                       gelman = TRUE,
                       heidel = FALSE,  
                       geweke = TRUE,   
                       diag_save = TRUE,
                       diag_name = "PEMA_SIDER_diagnostics",
                       indiv_effect = FALSE,       
                       plot_post_save_png = FALSE, 
                       plot_pairs_save_png = FALSE,
                       plot_xy_save_png = FALSE)
##############################


##############################
# Create diagnostics, summary statistics, and posterior plots
##############################
output_JAGS(PEMA_SIDER, mix.PEMA.SIDER, source_PEMA_SIDER, output_options)
print(PEMA_SIDER)
getwd()
save.image("PEMA_SIDER_Jags.RData") 
##############################


##################################################
#Put data into a dataframe
##################################################
setwd("~/Isotopic_Routing_Small_Mammals/Mixing_Models/PEMA")

require(R2jags) # Load packages 

# Load data 
load("~/Isotopic_Routing_Small_Mammals/Mixing_Models/PEMA/PEMA_SIDER_Jags.RData")
attach.jags(PEMA_SIDER) # call model file
print(PEMA_SIDER)

#assign names to each food source https://github.com/brianstock/MixSIAR/issues/30
PEMA_SIDER_Post <- data.frame(Source = "SIDER", AM_Fungi = p.global[,1], Arthropods = p.global[,2],
                              Berries = p.global[,3], EM_Fungi = p.global[,4], Red_Maple = p.global[,5], 
                              Fungi = p.global[,1]+p.global[,4]) 

PEMA_SIDER_Post<-data.frame(Number=rownames(PEMA_SIDER_Post), PEMA_SIDER_Post)#add posterior number
require(tidyr)
PEMA_SIDER_Posterior <- PEMA_SIDER_Post %>% gather(Source,value,AM_Fungi:Fungi)#Use tidyr to gather posterior for each food source
PEMA_SIDER_Posterior["Model"] <- "SIDER"#Add model type
PEMA_SIDER_Posterior["Species"] <- "PEMA"#Add Species

#save raw posterior data
write.csv(PEMA_SIDER_Posterior, file = "PEMA_SIDER_Posterior.csv", row.names=FALSE)#save raw posterior data
write.csv(PEMA_SIDER_Posterior, file = "~/Isotopic_Routing_Small_Mammals/Mixing_Models/Mixing_Model_Figure/PEMA_SIDER_Posterior.csv", row.names=FALSE)
################################################################################











################################################################################
#SIDER_update
################################################################################
setwd("~/Isotopic_Routing_Small_Mammals/Mixing_Models/PEMA")#set working directory
#############################
#load isotopic data
#############################
#Take stomach samples, convert to SIDER_update, and select species
Hair<- read.csv("Bartlett_isotope_hair_data.csv",header=T)
library(dplyr)
PEMA_SIDER_update<-filter(Hair, Species=="PEMA")
write.csv(PEMA_SIDER_update, "PEMA_SIDER_update.csv")

library(MixSIAR)
mix.filename.PEMA_SIDER_update <-("~/Isotopic_Routing_Small_Mammals/Mixing_Models/PEMA/PEMA_SIDER_update.csv") #load isotopic data for MixSIAR

mix.PEMA.SIDER_update <- load_mix_data(filename=mix.filename.PEMA_SIDER_update, #name of the CSV file with mix/consumer data
                                       iso_names=c("d13C","d15N"), #tracers/isotopes
                                       factors=c(NULL),            #factors
                                       fac_random=c(NULL),         #random effects
                                       fac_nested=c(NULL),         #Nested effects
                                       cont_effects=NULL)          #continuous effects
#############################


#############################
# Load food source data
#############################
source.filename <- ("~/Isotopic_Routing_Small_Mammals/Mixing_Models/PEMA/Bartlett_diet_items.csv")

source_PEMA_SIDER_update <- load_source_data(filename=source.filename, #name of the CSV file with food source data
                                             source_factors=NULL,     #column headings of random/fixed effects you have source data by
                                             conc_dep=TRUE,           #TRUE or FALSE for concentration dependence data 
                                             data_type="raw",         #"means" = means + sd or "raw" = raw data
                                             mix.PEMA.SIDER_update)

# Load discrimination data
#Select species and method of estimated TDF
Est_DFT<- read.csv("Estimated_TDF.csv",header=T)#discrimination factors
PEMA_SIDER_update_d13C_TDF<-filter(Est_DFT,Species_Abr == "PEMA" & Method =="SIDER_update" & Isotope == "d13C")
PEMA_SIDER_update_d15N_TDF<-filter(Est_DFT,Species_Abr == "PEMA" & Method =="SIDER_update" & Isotope == "d15N")


Type<- c("AM Fungi","Arthropods","Red Maple", "EM Fungi", "Berries")
TDF_PEMA_SIDER_update <- as.data.frame(Type)
TDF_PEMA_SIDER_update$Meand13C<-PEMA_SIDER_update_d13C_TDF$TDF_Mean	
TDF_PEMA_SIDER_update$SDd13C<-PEMA_SIDER_update_d13C_TDF$TDF_SD			
TDF_PEMA_SIDER_update$Meand15N<-PEMA_SIDER_update_d15N_TDF$TDF_Mean			
TDF_PEMA_SIDER_update$SDd15N<-PEMA_SIDER_update_d15N_TDF$TDF_SD		
write.csv(TDF_PEMA_SIDER_update, "PEMA_Discrimination_SIDER_update.csv", row.names=FALSE)

PEMA_Discrimination_SIDER_update <- ("~/Isotopic_Routing_Small_Mammals/Mixing_Models/PEMA/PEMA_Discrimination_SIDER_update.csv")
PEMA_Discrimination_SIDER_update <- load_discr_data(filename=PEMA_Discrimination_SIDER_update, mix.PEMA.SIDER_update)
#############################


#############################
# Make isospace plot
#############################
plot_data(filename="PEMA_SIDER_update_isospace", #name to save the isospace plot as
          plot_save_pdf=FALSE,       #TRUE or FALSE, should mix.SIAR save the plot as a .pdf?
          plot_save_png=TRUE,      #TRUE or FALSE, should mix.SIAR save the plot as a .png?
          mix.PEMA.SIDER_update,source_PEMA_SIDER_update,PEMA_Discrimination_SIDER_update)

# If 2 isotopes/tracers, calculate normalized surface area of the convex hull polygon(s)
#   *Note 1: discrimination SD is added to the source SD (see calc_area.r for details)
if(mix.PEMA.SIDER_update$n.iso==2) calc_area(source=source_PEMA_SIDER_update,mix=mix.PEMA.SIDER_update,discr=PEMA_Discrimination_SIDER_update)
#############################


#############################
# Define prior
#############################
kw.alpha <- c(1,1,1,1,1)

kw.alpha <- kw.alpha*length(kw.alpha)/sum(kw.alpha)# Generate alpha hyperparameters scaling sum(alpha)=n.sources

kw.alpha[which(kw.alpha==0)] <- 0.01# the Dirichlet hyperparameters for the alpha.prior cannot be 0 (but can set = .01)

#plot using "plot_prior"
plot_prior(alpha.prior=kw.alpha,
           source=source_PEMA_SIDER_update,
           plot_save_pdf=FALSE,
           plot_save_png=TRUE,
           filename="PEMA_SIDER_update_prior_plot_kw_inf")
#   RED = your prior
#   DARK GREY = "uninformative"/generalist (alpha = 1)
#   LIGHT GREY = "uninformative" Jeffrey's prior (alpha = 1/n.sources)
#############################


#############################
# Write JAGS model file (define model structure)
##############################
# There are 3 error term options available:
#   1. Residual * Process (resid_err = TRUE, process_err = TRUE)
#   2. Residual only (resid_err = TRUE, process_err = FALSE)
#   3. Process only (resid_err = FALSE, process_err = TRUE)
PEMA_SIDER_update_model <- "PEMA_SIDER_update_model.txt"
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(PEMA_SIDER_update_model, resid_err, process_err, mix.PEMA.SIDER_update, source_PEMA_SIDER_update)
##############################


##############################
# Run PEMA SIDER_update model
##############################
# MCMC options:
#run <- list(chainLength=30000, burn=15000, thin=25, chains=3, calcDIC=TRUE)
costum_run <- list(chainLength=200000, burn=50000, thin=50, chains=3, calcDIC=TRUE)

#Model
PEMA_SIDER_update<- run_model(run=costum_run,mix.PEMA.SIDER_update,source_PEMA_SIDER_update,PEMA_Discrimination_SIDER_update,PEMA_SIDER_update_model,
                              alpha.prior=kw.alpha, resid_err, process_err)
##############################


##############################
# Process JAGS output
##############################
# Choose output options (see ?output_options for details)
output_options <- list(summary_save = TRUE,                 
                       summary_name = "PEMA_SIDER_update_summary_statistics", 
                       sup_post = FALSE,                    
                       plot_post_save_pdf = TRUE,           
                       plot_post_name = "PEMA_SIDER_update_posterior_density",
                       sup_pairs = FALSE,             
                       plot_pairs_save_pdf = TRUE,    
                       plot_pairs_name = "PEMA_SIDER_update_pairs_plot",
                       sup_xy = TRUE,           
                       plot_xy_save_pdf = FALSE,
                       plot_xy_name = "PEMA_SIDER_update_xy_plot",
                       gelman = TRUE,
                       heidel = FALSE,  
                       geweke = TRUE,   
                       diag_save = TRUE,
                       diag_name = "PEMA_SIDER_update_diagnostics",
                       indiv_effect = FALSE,       
                       plot_post_save_png = FALSE, 
                       plot_pairs_save_png = FALSE,
                       plot_xy_save_png = FALSE)
##############################


##############################
# Create diagnostics, summary statistics, and posterior plots
##############################
output_JAGS(PEMA_SIDER_update, mix.PEMA.SIDER_update, source_PEMA_SIDER_update, output_options)
print(PEMA_SIDER_update)
getwd()
save.image("PEMA_SIDER_update_Jags.RData") 
##############################


##################################################
#Put data into a dataframe
##################################################
setwd("~/Isotopic_Routing_Small_Mammals/Mixing_Models/PEMA")

require(R2jags) # Load packages 

# Load data 
load("~/Isotopic_Routing_Small_Mammals/Mixing_Models/PEMA/PEMA_SIDER_update_Jags.RData")
attach.jags(PEMA_SIDER_update) # call model file
print(PEMA_SIDER_update)

#assign names to each food source https://github.com/brianstock/MixSIAR/issues/30
PEMA_SIDER_update_Post <- data.frame(Source = "SIDER_update", AM_Fungi = p.global[,1], Arthropods = p.global[,2],
                                     Berries = p.global[,3], EM_Fungi = p.global[,4], Red_Maple = p.global[,5], 
                                     Fungi = p.global[,1]+p.global[,4]) 

PEMA_SIDER_update_Post<-data.frame(Number=rownames(PEMA_SIDER_update_Post), PEMA_SIDER_update_Post)#add posterior number
require(tidyr)
PEMA_SIDER_update_Posterior <- PEMA_SIDER_update_Post %>% gather(Source,value,AM_Fungi:Fungi)#Use tidyr to gather posterior for each food source
PEMA_SIDER_update_Posterior["Model"] <- "SIDER_update"#Add model type
PEMA_SIDER_update_Posterior["Species"] <- "PEMA"#Add Species

#save raw posterior data
write.csv(PEMA_SIDER_update_Posterior, file = "PEMA_SIDER_update_Posterior.csv", row.names=FALSE)#save raw posterior data
write.csv(PEMA_SIDER_update_Posterior, file = "~/Isotopic_Routing_Small_Mammals/Mixing_Models/Mixing_Model_Figure/PEMA_SIDER_update_Posterior.csv", row.names=FALSE)
################################################################################









################################################################################
#Meta_analysis
################################################################################
setwd("~/Isotopic_Routing_Small_Mammals/Mixing_Models/PEMA")#set working directory
#############################
#load isotopic data
#############################
#Take stomach samples, convert to Meta_analysis, and select species
Hair<- read.csv("Bartlett_isotope_hair_data.csv",header=T)
library(dplyr)
PEMA_Meta_analysis<-filter(Hair, Species=="PEMA")
write.csv(PEMA_Meta_analysis, "PEMA_Meta_analysis.csv")

library(MixSIAR)
mix.filename.PEMA_Meta_analysis <-("~/Isotopic_Routing_Small_Mammals/Mixing_Models/PEMA/PEMA_Meta_analysis.csv") #load isotopic data for MixSIAR

mix.PEMA.Meta_analysis <- load_mix_data(filename=mix.filename.PEMA_Meta_analysis, #name of the CSV file with mix/consumer data
                                        iso_names=c("d13C","d15N"), #tracers/isotopes
                                        factors=c(NULL),            #factors
                                        fac_random=c(NULL),         #random effects
                                        fac_nested=c(NULL),         #Nested effects
                                        cont_effects=NULL)          #continuous effects
#############################


#############################
# Load food source data
#############################
source.filename <- ("~/Isotopic_Routing_Small_Mammals/Mixing_Models/PEMA/Bartlett_diet_items.csv")

source_PEMA_Meta_analysis <- load_source_data(filename=source.filename, #name of the CSV file with food source data
                                              source_factors=NULL,     #column headings of random/fixed effects you have source data by
                                              conc_dep=TRUE,           #TRUE or FALSE for concentration dependence data 
                                              data_type="raw",         #"means" = means + sd or "raw" = raw data
                                              mix.PEMA.Meta_analysis)

# Load discrimination data
#Select species and method of estimated TDF
Est_DFT<- read.csv("Estimated_TDF.csv",header=T)#discrimination factors
PEMA_Meta_analysis_d13C_TDF<-filter(Est_DFT,Species_Abr == "PEMA" & Method =="Meta_analysis" & Isotope == "d13C")
PEMA_Meta_analysis_d15N_TDF<-filter(Est_DFT,Species_Abr == "PEMA" & Method =="Meta_analysis" & Isotope == "d15N")


Type<- c("AM Fungi","Arthropods","Red Maple", "EM Fungi", "Berries")
TDF_PEMA_Meta_analysis <- as.data.frame(Type)
TDF_PEMA_Meta_analysis$Meand13C<-PEMA_Meta_analysis_d13C_TDF$TDF_Mean	
TDF_PEMA_Meta_analysis$SDd13C<-PEMA_Meta_analysis_d13C_TDF$TDF_SD			
TDF_PEMA_Meta_analysis$Meand15N<-PEMA_Meta_analysis_d15N_TDF$TDF_Mean			
TDF_PEMA_Meta_analysis$SDd15N<-PEMA_Meta_analysis_d15N_TDF$TDF_SD		
write.csv(TDF_PEMA_Meta_analysis, "PEMA_Discrimination_Meta_analysis.csv", row.names=FALSE)

PEMA_Discrimination_Meta_analysis <- ("~/Isotopic_Routing_Small_Mammals/Mixing_Models/PEMA/PEMA_Discrimination_Meta_analysis.csv")
PEMA_Discrimination_Meta_analysis <- load_discr_data(filename=PEMA_Discrimination_Meta_analysis, mix.PEMA.Meta_analysis)
#############################


#############################
# Make isospace plot
#############################
plot_data(filename="PEMA_Meta_analysis_isospace", #name to save the isospace plot as
          plot_save_pdf=FALSE,       #TRUE or FALSE, should mix.SIAR save the plot as a .pdf?
          plot_save_png=TRUE,      #TRUE or FALSE, should mix.SIAR save the plot as a .png?
          mix.PEMA.Meta_analysis,source_PEMA_Meta_analysis,PEMA_Discrimination_Meta_analysis)

# If 2 isotopes/tracers, calculate normalized surface area of the convex hull polygon(s)
#   *Note 1: discrimination SD is added to the source SD (see calc_area.r for details)
if(mix.PEMA.Meta_analysis$n.iso==2) calc_area(source=source_PEMA_Meta_analysis,mix=mix.PEMA.Meta_analysis,discr=PEMA_Discrimination_Meta_analysis)
#############################


#############################
# Define prior
#############################
kw.alpha <- c(1,1,1,1,1)

kw.alpha <- kw.alpha*length(kw.alpha)/sum(kw.alpha)# Generate alpha hyperparameters scaling sum(alpha)=n.sources

kw.alpha[which(kw.alpha==0)] <- 0.01# the Dirichlet hyperparameters for the alpha.prior cannot be 0 (but can set = .01)

#plot using "plot_prior"
plot_prior(alpha.prior=kw.alpha,
           source=source_PEMA_Meta_analysis,
           plot_save_pdf=FALSE,
           plot_save_png=TRUE,
           filename="PEMA_Meta_analysis_prior_plot_kw_inf")
#   RED = your prior
#   DARK GREY = "uninformative"/generalist (alpha = 1)
#   LIGHT GREY = "uninformative" Jeffrey's prior (alpha = 1/n.sources)
#############################


#############################
# Write JAGS model file (define model structure)
##############################
# There are 3 error term options available:
#   1. Residual * Process (resid_err = TRUE, process_err = TRUE)
#   2. Residual only (resid_err = TRUE, process_err = FALSE)
#   3. Process only (resid_err = FALSE, process_err = TRUE)
PEMA_Meta_analysis_model <- "PEMA_Meta_analysis_model.txt"
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(PEMA_Meta_analysis_model, resid_err, process_err, mix.PEMA.Meta_analysis, source_PEMA_Meta_analysis)
##############################


##############################
# Run PEMA Meta_analysis model
##############################
# MCMC options:
#run <- list(chainLength=30000, burn=15000, thin=25, chains=3, calcDIC=TRUE)
costum_run <- list(chainLength=200000, burn=50000, thin=50, chains=3, calcDIC=TRUE)

#Model
PEMA_Meta_analysis<- run_model(run=costum_run,mix.PEMA.Meta_analysis,source_PEMA_Meta_analysis,PEMA_Discrimination_Meta_analysis,PEMA_Meta_analysis_model,
                               alpha.prior=kw.alpha, resid_err, process_err)
##############################


##############################
# Process JAGS output
##############################
# Choose output options (see ?output_options for details)
output_options <- list(summary_save = TRUE,                 
                       summary_name = "PEMA_Meta_analysis_summary_statistics", 
                       sup_post = FALSE,                    
                       plot_post_save_pdf = TRUE,           
                       plot_post_name = "PEMA_Meta_analysis_posterior_density",
                       sup_pairs = FALSE,             
                       plot_pairs_save_pdf = TRUE,    
                       plot_pairs_name = "PEMA_Meta_analysis_pairs_plot",
                       sup_xy = TRUE,           
                       plot_xy_save_pdf = FALSE,
                       plot_xy_name = "PEMA_Meta_analysis_xy_plot",
                       gelman = TRUE,
                       heidel = FALSE,  
                       geweke = TRUE,   
                       diag_save = TRUE,
                       diag_name = "PEMA_Meta_analysis_diagnostics",
                       indiv_effect = FALSE,       
                       plot_post_save_png = FALSE, 
                       plot_pairs_save_png = FALSE,
                       plot_xy_save_png = FALSE)
##############################


##############################
# Create diagnostics, summary statistics, and posterior plots
##############################
output_JAGS(PEMA_Meta_analysis, mix.PEMA.Meta_analysis, source_PEMA_Meta_analysis, output_options)
print(PEMA_Meta_analysis)
getwd()
save.image("PEMA_Meta_analysis_Jags.RData") 
##############################


##################################################
#Put data into a dataframe
##################################################
setwd("~/Isotopic_Routing_Small_Mammals/Mixing_Models/PEMA")

require(R2jags) # Load packages 

# Load data 
load("~/Isotopic_Routing_Small_Mammals/Mixing_Models/PEMA/PEMA_Meta_analysis_Jags.RData")
attach.jags(PEMA_Meta_analysis) # call model file
print(PEMA_Meta_analysis)

#assign names to each food source https://github.com/brianstock/MixSIAR/issues/30
PEMA_Meta_analysis_Post <- data.frame(Source = "Meta_analysis", AM_Fungi = p.global[,1], Arthropods = p.global[,2],
                                      Berries = p.global[,3], EM_Fungi = p.global[,4], Red_Maple = p.global[,5], 
                                      Fungi = p.global[,1]+p.global[,4]) 

PEMA_Meta_analysis_Post<-data.frame(Number=rownames(PEMA_Meta_analysis_Post), PEMA_Meta_analysis_Post)#add posterior number

require(tidyr)
PEMA_Meta_analysis_Posterior <- PEMA_Meta_analysis_Post %>% gather(Source,value,AM_Fungi:Fungi)#Use tidyr to gather posterior for each food source
PEMA_Meta_analysis_Posterior["Model"] <- "Meta_analysis"#Add model type
PEMA_Meta_analysis_Posterior["Species"] <- "PEMA"#Add Species

#save raw posterior data
write.csv(PEMA_Meta_analysis_Posterior, file = "PEMA_Meta_analysis_Posterior.csv", row.names=FALSE)#save raw posterior data
write.csv(PEMA_Meta_analysis_Posterior, file = "~/Isotopic_Routing_Small_Mammals/Mixing_Models/Mixing_Model_Figure/PEMA_Meta_analysis_Posterior.csv", row.names=FALSE)
################################################################################










################################################################################
#Blarina bevicauda
################################################################################
rm(list=ls())# clears workspace
setwd("~/Isotopic_Routing_Small_Mammals/Mixing_Models/BLBR")#set working directory


################################################################################
#Diet
################################################################################

#############################
#load isotopic data
#############################
#Take stomach samples, convert to diet, and select species
library(dplyr)
Isotopes<- read.csv("Field_study_data.csv",header=T)
Diet<-filter(Isotopes, Type =="Stomach")#select stomach data
Diet<-Diet %>% mutate(d15N=ifelse(Species == "MYGA",d15N-1.85,d15N-1))#subtract 1 to account digestive enzymes of stomach (1.85 for MGYA)
Diet<-mutate(Diet, Type = "Diet")#code 'type' as diet
BLBR_Diet<- filter(Diet,Species =="BLBR") #Filter out to desired species

write.csv(BLBR_Diet, "BLBR_Diet.csv")

library(MixSIAR)
mix.filename.BLBR_Diet <-("~/Isotopic_Routing_Small_Mammals/Mixing_Models/BLBR/BLBR_Diet.csv") #load isotopic data for MixSIAR

mix.BLBR.Diet <- load_mix_data(filename=mix.filename.BLBR_Diet, #name of the CSV file with mix/consumer data
                               iso_names=c("d13C","d15N"), #tracers/isotopes
                               factors=c(NULL),            #factors
                               fac_random=c(NULL),         #random effects
                               fac_nested=c(NULL),         #Nested effects
                               cont_effects=NULL)          #continuous effects
#############################


#############################
# Load food source data
#############################
source.filename <- ("~/Isotopic_Routing_Small_Mammals/Mixing_Models/BLBR/Bartlett_diet_items.csv")

source_BLBR_Diet <- load_source_data(filename=source.filename, #name of the CSV file with food source data
                                     source_factors=NULL,     #column headings of random/fixed effects you have source data by
                                     conc_dep=TRUE,           #TRUE or FALSE for concentration dependence data 
                                     data_type="raw",         #"means" = means + sd or "raw" = raw data
                                     mix.BLBR.Diet)

# Load discrimination data
#No TDF for diet samples
Type<- c("AM Fungi","Arthropods", "Red Maple", "EM Fungi", "Berries")
TDF_BLBR_Field <- as.data.frame(Type)
TDF_BLBR_Field$Meand13C<-0	
TDF_BLBR_Field$SDd13C<-0		
TDF_BLBR_Field$Meand15N<-0		
TDF_BLBR_Field$SDd15N<-0	
write.csv(TDF_BLBR_Field, "Discrimination_Diet.csv", row.names=FALSE)

Discrimination_Diet <- ("~/Isotopic_Routing_Small_Mammals/Mixing_Models/BLBR/Discrimination_Diet.csv")
discr_Diet <- load_discr_data(filename=Discrimination_Diet, mix.BLBR.Diet)
#############################


#############################
# Make isospace plot
#############################
plot_data(filename="BLBR_Diet_isospace", #name to save the isospace plot as
          plot_save_pdf=FALSE,       #TRUE or FALSE, should mix.SIAR save the plot as a .pdf?
          plot_save_png=TRUE,      #TRUE or FALSE, should mix.SIAR save the plot as a .png?
          mix.BLBR.Diet,source_BLBR_Diet,discr_Diet)

# If 2 isotopes/tracers, calculate normalized surface area of the convex hull polygon(s)
#   *Note 1: discrimination SD is added to the source SD (see calc_area.r for details)
if(mix.BLBR.Diet$n.iso==2) calc_area(source=source_BLBR_Diet,mix=mix.BLBR.Diet,discr=discr_Diet)
#############################


#############################
# Define prior
#############################
kw.alpha <- c(1,1,1,0,1)

kw.alpha <- kw.alpha*length(kw.alpha)/sum(kw.alpha)# Generate alpha hyperparameters scaling sum(alpha)=n.sources

kw.alpha[which(kw.alpha==0)] <- 0.01# the Dirichlet hyperparameters for the alpha.prior cannot be 0 (but can set = .01)

#plot using "plot_prior"
plot_prior(alpha.prior=kw.alpha,
           source=source_BLBR_Diet,
           plot_save_pdf=FALSE,
           plot_save_png=TRUE,
           filename="BLBR_Diet_prior_plot_kw_inf")
#   RED = your prior
#   DARK GREY = "uninformative"/generalist (alpha = 1)
#   LIGHT GREY = "uninformative" Jeffrey's prior (alpha = 1/n.sources)
#############################


#############################
# Write JAGS model file (define model structure)
##############################
# There are 3 error term options available:
#   1. Residual * Process (resid_err = TRUE, process_err = TRUE)
#   2. Residual only (resid_err = TRUE, process_err = FALSE)
#   3. Process only (resid_err = FALSE, process_err = TRUE)
BLBR_Diet_model <- "BLBR_Diet_model.txt"
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(BLBR_Diet_model, resid_err, process_err, mix.BLBR.Diet, source_BLBR_Diet)
##############################


##############################
# Run BLBR Diet model
##############################
# MCMC options:
#run <- list(chainLength=30000, burn=15000, thin=25, chains=3, calcDIC=TRUE)
costum_run <- list(chainLength=200000, burn=50000, thin=50, chains=3, calcDIC=TRUE)

#Model
BLBR_Diet<- run_model(run=costum_run,mix.BLBR.Diet,source_BLBR_Diet,discr_Diet,BLBR_Diet_model,
                      alpha.prior=kw.alpha, resid_err, process_err)
##############################


##############################
# Process JAGS output
##############################
# Choose output options (see ?output_options for details)
output_options <- list(summary_save = TRUE,                 
                       summary_name = "BLBR_Diet_summary_statistics", 
                       sup_post = FALSE,                    
                       plot_post_save_pdf = TRUE,           
                       plot_post_name = "BLBR_Diet_posterior_density",
                       sup_pairs = FALSE,             
                       plot_pairs_save_pdf = TRUE,    
                       plot_pairs_name = "BLBR_Diet_pairs_plot",
                       sup_xy = TRUE,           
                       plot_xy_save_pdf = TRUE,
                       plot_xy_name = "BLBR_Diet_xy_plot",
                       gelman = TRUE,
                       heidel = FALSE,  
                       geweke = TRUE,   
                       diag_save = TRUE,
                       diag_name = "BLBR_Diet_diagnostics",
                       indiv_effect = FALSE,       
                       plot_post_save_png = TRUE, 
                       plot_pairs_save_png = TRUE,
                       plot_xy_save_png = TRUE)
##############################


##############################
# Create diagnostics, summary statistics, and posterior plots
##############################
output_JAGS(BLBR_Diet, mix.BLBR.Diet, source_BLBR_Diet, output_options)
print(BLBR_Diet)
getwd()
save.image("BLBR_Diet_Jags.RData") 
##############################


##################################################
#Put data into a dataframe
##################################################
setwd("~/Isotopic_Routing_Small_Mammals/Mixing_Models/BLBR")

require(R2jags) # Load packages 

# Load data 
load("~/Isotopic_Routing_Small_Mammals/Mixing_Models/BLBR/BLBR_Diet_Jags.RData")
attach.jags(BLBR_Diet) # call model file
print(BLBR_Diet)

#assign names to each food source https://github.com/brianstock/MixSIAR/issues/30
BLBR_Diet_Post <- data.frame(Source = "Diet", AM_Fungi = p.global[,1],Arthropods = p.global[,2],
                             Berries = p.global[,3], EM_Fungi = p.global[,4], Red_Maple = p.global[,5], 
                             Fungi = p.global[,1]+p.global[,4]) 

BLBR_Diet_Post<-data.frame(Number=rownames(BLBR_Diet_Post), BLBR_Diet_Post)

require(tidyr)
BLBR_Diet_Posterior <- BLBR_Diet_Post %>% gather(Source,value,AM_Fungi:Fungi)#Use tidyr to gather posterior for each food source
BLBR_Diet_Posterior["Model"] <- "Diet"#Add model type
BLBR_Diet_Posterior["Species"] <- "BLBR"#Add Species
#save raw posterior data
write.csv(BLBR_Diet_Posterior, file = "BLBR_Diet_Posterior.csv", row.names=F)#save raw posterior data
write.csv(BLBR_Diet_Posterior, file = "~/Isotopic_Routing_Small_Mammals/Mixing_Models/Mixing_Model_Figure/BLBR_Diet_Posterior.csv", row.names=F)
################################################################################











################################################################################
#Field
################################################################################
setwd("~/Isotopic_Routing_Small_Mammals/Mixing_Models/BLBR")#set working directory
#############################
#load isotopic data
#############################
#Take stomach samples, convert to Field, and select species
Hair<- read.csv("Bartlett_isotope_hair_data.csv",header=T)
library(dplyr)
BLBR_Field<-filter(Hair, Species=="BLBR")
write.csv(BLBR_Field, "BLBR_Field.csv")

library(MixSIAR)
mix.filename.BLBR_Field <-("~/Isotopic_Routing_Small_Mammals/Mixing_Models/BLBR/BLBR_Field.csv") #load isotopic data for MixSIAR

mix.BLBR.Field <- load_mix_data(filename=mix.filename.BLBR_Field, #name of the CSV file with mix/consumer data
                                iso_names=c("d13C","d15N"), #tracers/isotopes
                                factors=c(NULL),            #factors
                                fac_random=c(NULL),         #random effects
                                fac_nested=c(NULL),         #Nested effects
                                cont_effects=NULL)          #continuous effects
#############################


#############################
# Load food source data
#############################
source.filename <- ("~/Isotopic_Routing_Small_Mammals/Mixing_Models/BLBR/Bartlett_diet_items.csv")

source_BLBR_Field <- load_source_data(filename=source.filename, #name of the CSV file with food source data
                                      source_factors=NULL,     #column headings of random/fixed effects you have source data by
                                      conc_dep=TRUE,           #TRUE or FALSE for concentration dependence data 
                                      data_type="raw",         #"means" = means + sd or "raw" = raw data
                                      mix.BLBR.Field)

# Load discrimination data
#Select species and method of estimated TDF
Est_DFT<- read.csv("Estimated_TDF.csv",header=T)#discrimination factors
BLBR_Field_d13C_TDF<-filter(Est_DFT,Species_Abr == "BLBR" & Method =="Field" & Isotope == "d13C")
BLBR_Field_d15N_TDF<-filter(Est_DFT,Species_Abr == "BLBR" & Method =="Field" & Isotope == "d15N")


Type<- c("AM Fungi","Arthropods","Red Maple", "EM Fungi", "Berries")
TDF_BLBR_Field <- as.data.frame(Type)
TDF_BLBR_Field$Meand13C<-BLBR_Field_d13C_TDF$TDF_Mean	
TDF_BLBR_Field$SDd13C<-BLBR_Field_d13C_TDF$TDF_SD			
TDF_BLBR_Field$Meand15N<-BLBR_Field_d15N_TDF$TDF_Mean			
TDF_BLBR_Field$SDd15N<-BLBR_Field_d15N_TDF$TDF_SD		
write.csv(TDF_BLBR_Field, "BLBR_Discrimination_Field.csv", row.names=FALSE)

BLBR_Discrimination_Field <- ("~/Isotopic_Routing_Small_Mammals/Mixing_Models/BLBR/BLBR_Discrimination_Field.csv")
BLBR_Discrimination_Field <- load_discr_data(filename=BLBR_Discrimination_Field, mix.BLBR.Field)
#############################


#############################
# Make isospace plot
#############################
plot_data(filename="BLBR_Field_isospace", #name to save the isospace plot as
          plot_save_pdf=FALSE,       #TRUE or FALSE, should mix.SIAR save the plot as a .pdf?
          plot_save_png=TRUE,      #TRUE or FALSE, should mix.SIAR save the plot as a .png?
          mix.BLBR.Field,source_BLBR_Field,BLBR_Discrimination_Field)

# If 2 isotopes/tracers, calculate normalized surface area of the convex hull polygon(s)
#   *Note 1: discrimination SD is added to the source SD (see calc_area.r for details)
if(mix.BLBR.Field$n.iso==2) calc_area(source=source_BLBR_Field,mix=mix.BLBR.Field,discr=BLBR_Discrimination_Field)
#############################


#############################
# Define prior
#############################
kw.alpha <- c(1,1,1,0,1)

kw.alpha <- kw.alpha*length(kw.alpha)/sum(kw.alpha)# Generate alpha hyperparameters scaling sum(alpha)=n.sources

kw.alpha[which(kw.alpha==0)] <- 0.01# the Dirichlet hyperparameters for the alpha.prior cannot be 0 (but can set = .01)

#plot using "plot_prior"
plot_prior(alpha.prior=kw.alpha,
           source=source_BLBR_Field,
           plot_save_pdf=FALSE,
           plot_save_png=TRUE,
           filename="BLBR_Field_prior_plot_kw_inf")
#   RED = your prior
#   DARK GREY = "uninformative"/generalist (alpha = 1)
#   LIGHT GREY = "uninformative" Jeffrey's prior (alpha = 1/n.sources)
#############################


#############################
# Write JAGS model file (define model structure)
##############################
# There are 3 error term options available:
#   1. Residual * Process (resid_err = TRUE, process_err = TRUE)
#   2. Residual only (resid_err = TRUE, process_err = FALSE)
#   3. Process only (resid_err = FALSE, process_err = TRUE)
BLBR_Field_model <- "BLBR_Field_model.txt"
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(BLBR_Field_model, resid_err, process_err, mix.BLBR.Field, source_BLBR_Field)
##############################


##############################
# Run BLBR Field model
##############################
# MCMC options:
#run <- list(chainLength=30000, burn=15000, thin=25, chains=3, calcDIC=TRUE)
costum_run <- list(chainLength=200000, burn=50000, thin=50, chains=3, calcDIC=TRUE)

#Model
BLBR_Field<- run_model(run=costum_run,mix.BLBR.Field,source_BLBR_Field,BLBR_Discrimination_Field,BLBR_Field_model,
                       alpha.prior=kw.alpha, resid_err, process_err)
##############################


##############################
# Process JAGS output
##############################
# Choose output options (see ?output_options for details)
output_options <- list(summary_save = TRUE,                 
                       summary_name = "BLBR_Field_summary_statistics", 
                       sup_post = FALSE,                    
                       plot_post_save_pdf = TRUE,           
                       plot_post_name = "BLBR_Field_posterior_density",
                       sup_pairs = FALSE,             
                       plot_pairs_save_pdf = TRUE,    
                       plot_pairs_name = "BLBR_Field_pairs_plot",
                       sup_xy = TRUE,           
                       plot_xy_save_pdf = FALSE,
                       plot_xy_name = "BLBR_Field_xy_plot",
                       gelman = TRUE,
                       heidel = FALSE,  
                       geweke = TRUE,   
                       diag_save = TRUE,
                       diag_name = "BLBR_Field_diagnostics",
                       indiv_effect = FALSE,       
                       plot_post_save_png = FALSE, 
                       plot_pairs_save_png = FALSE,
                       plot_xy_save_png = FALSE)
##############################


##############################
# Create diagnostics, summary statistics, and posterior plots
##############################
output_JAGS(BLBR_Field, mix.BLBR.Field, source_BLBR_Field, output_options)
print(BLBR_Field)
getwd()
save.image("BLBR_Field_Jags.RData") 
##############################


##################################################
#Put data into a dataframe
##################################################
setwd("~/Isotopic_Routing_Small_Mammals/Mixing_Models/BLBR")

require(R2jags) # Load packages 

# Load data 
load("~/Isotopic_Routing_Small_Mammals/Mixing_Models/BLBR/BLBR_Field_Jags.RData")
attach.jags(BLBR_Field) # call model file
print(BLBR_Field)

#assign names to each food source https://github.com/brianstock/MixSIAR/issues/30
BLBR_Field_Post <- data.frame(Source = "Field", AM_Fungi = p.global[,1],Arthropods = p.global[,2],
                              Berries = p.global[,3], EM_Fungi = p.global[,4], Red_Maple = p.global[,5], 
                              Fungi = p.global[,1]+p.global[,4]) 

BLBR_Field_Post<-data.frame(Number=rownames(BLBR_Field_Post), BLBR_Field_Post)#add posterior number
require(tidyr)
BLBR_Field_Posterior <- BLBR_Field_Post %>% gather(Source,value,AM_Fungi:Fungi)#Use tidyr to gather posterior for each food source
BLBR_Field_Posterior["Model"] <- "Field"#Add model type
BLBR_Field_Posterior["Species"] <- "BLBR"#Add Species

#save raw posterior data
write.csv(BLBR_Field_Posterior, file = "BLBR_Field_Posterior.csv", row.names=FALSE)#save raw posterior data
write.csv(BLBR_Field_Posterior, file = "~/Isotopic_Routing_Small_Mammals/Mixing_Models/Mixing_Model_Figure/BLBR_Field_Posterior.csv", row.names=FALSE)
################################################################################










################################################################################
#SIDER
################################################################################
setwd("~/Isotopic_Routing_Small_Mammals/Mixing_Models/BLBR")#set working directory
#############################
#load isotopic data
#############################
#Take stomach samples, convert to SIDER, and select species
Hair<- read.csv("Bartlett_isotope_hair_data.csv",header=T)
library(dplyr)
BLBR_SIDER<-filter(Hair, Species=="BLBR")
write.csv(BLBR_SIDER, "BLBR_SIDER.csv")

library(MixSIAR)
mix.filename.BLBR_SIDER <-("~/Isotopic_Routing_Small_Mammals/Mixing_Models/BLBR/BLBR_SIDER.csv") #load isotopic data for MixSIAR

mix.BLBR.SIDER <- load_mix_data(filename=mix.filename.BLBR_SIDER, #name of the CSV file with mix/consumer data
                                iso_names=c("d13C","d15N"), #tracers/isotopes
                                factors=c(NULL),            #factors
                                fac_random=c(NULL),         #random effects
                                fac_nested=c(NULL),         #Nested effects
                                cont_effects=NULL)          #continuous effects
#############################


#############################
# Load food source data
#############################
source.filename <- ("~/Isotopic_Routing_Small_Mammals/Mixing_Models/BLBR/Bartlett_diet_items.csv")

source_BLBR_SIDER <- load_source_data(filename=source.filename, #name of the CSV file with food source data
                                      source_factors=NULL,     #column headings of random/fixed effects you have source data by
                                      conc_dep=TRUE,           #TRUE or FALSE for concentration dependence data 
                                      data_type="raw",         #"means" = means + sd or "raw" = raw data
                                      mix.BLBR.SIDER)

# Load discrimination data
#Select species and method of estimated TDF
Est_DFT<- read.csv("Estimated_TDF.csv",header=T)#discrimination factors
BLBR_SIDER_d13C_TDF<-filter(Est_DFT,Species_Abr == "BLBR" & Method =="SIDER" & Isotope == "d13C")
BLBR_SIDER_d15N_TDF<-filter(Est_DFT,Species_Abr == "BLBR" & Method =="SIDER" & Isotope == "d15N")


Type<- c("AM Fungi","Arthropods", "Red Maple", "EM Fungi", "Berries")
TDF_BLBR_SIDER <- as.data.frame(Type)
TDF_BLBR_SIDER$Meand13C<-BLBR_SIDER_d13C_TDF$TDF_Mean	
TDF_BLBR_SIDER$SDd13C<-BLBR_SIDER_d13C_TDF$TDF_SD			
TDF_BLBR_SIDER$Meand15N<-BLBR_SIDER_d15N_TDF$TDF_Mean			
TDF_BLBR_SIDER$SDd15N<-BLBR_SIDER_d15N_TDF$TDF_SD		
write.csv(TDF_BLBR_SIDER, "BLBR_Discrimination_SIDER.csv", row.names=FALSE)

BLBR_Discrimination_SIDER <- ("~/Isotopic_Routing_Small_Mammals/Mixing_Models/BLBR/BLBR_Discrimination_SIDER.csv")
BLBR_Discrimination_SIDER <- load_discr_data(filename=BLBR_Discrimination_SIDER, mix.BLBR.SIDER)
#############################


#############################
# Make isospace plot
#############################
plot_data(filename="BLBR_SIDER_isospace", #name to save the isospace plot as
          plot_save_pdf=FALSE,       #TRUE or FALSE, should mix.SIAR save the plot as a .pdf?
          plot_save_png=TRUE,      #TRUE or FALSE, should mix.SIAR save the plot as a .png?
          mix.BLBR.SIDER,source_BLBR_SIDER,BLBR_Discrimination_SIDER)

# If 2 isotopes/tracers, calculate normalized surface area of the convex hull polygon(s)
#   *Note 1: discrimination SD is added to the source SD (see calc_area.r for details)
if(mix.BLBR.SIDER$n.iso==2) calc_area(source=source_BLBR_SIDER,mix=mix.BLBR.SIDER,discr=BLBR_Discrimination_SIDER)
#############################


#############################
# Define prior
#############################
kw.alpha <- c(1,1,1,0,1)

kw.alpha <- kw.alpha*length(kw.alpha)/sum(kw.alpha)# Generate alpha hyperparameters scaling sum(alpha)=n.sources

kw.alpha[which(kw.alpha==0)] <- 0.01# the Dirichlet hyperparameters for the alpha.prior cannot be 0 (but can set = .01)

#plot using "plot_prior"
plot_prior(alpha.prior=kw.alpha,
           source=source_BLBR_SIDER,
           plot_save_pdf=FALSE,
           plot_save_png=TRUE,
           filename="BLBR_SIDER_prior_plot_kw_inf")
#   RED = your prior
#   DARK GREY = "uninformative"/generalist (alpha = 1)
#   LIGHT GREY = "uninformative" Jeffrey's prior (alpha = 1/n.sources)
#############################


#############################
# Write JAGS model file (define model structure)
##############################
# There are 3 error term options available:
#   1. Residual * Process (resid_err = TRUE, process_err = TRUE)
#   2. Residual only (resid_err = TRUE, process_err = FALSE)
#   3. Process only (resid_err = FALSE, process_err = TRUE)
BLBR_SIDER_model <- "BLBR_SIDER_model.txt"
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(BLBR_SIDER_model, resid_err, process_err, mix.BLBR.SIDER, source_BLBR_SIDER)
##############################


##############################
# Run BLBR SIDER model
##############################
# MCMC options:
run <- list(chainLength=30000, burn=15000, thin=25, chains=3, calcDIC=TRUE)
costum_run <- list(chainLength=200000, burn=50000, thin=50, chains=3, calcDIC=TRUE)

#Model
BLBR_SIDER<- run_model(run=costum_run,mix.BLBR.SIDER,source_BLBR_SIDER,BLBR_Discrimination_SIDER,BLBR_SIDER_model,
                       alpha.prior=kw.alpha, resid_err, process_err)
##############################


##############################
# Process JAGS output
##############################
# Choose output options (see ?output_options for details)
output_options <- list(summary_save = TRUE,                 
                       summary_name = "BLBR_SIDER_summary_statistics", 
                       sup_post = FALSE,                    
                       plot_post_save_pdf = TRUE,           
                       plot_post_name = "BLBR_SIDER_posterior_density",
                       sup_pairs = FALSE,             
                       plot_pairs_save_pdf = TRUE,    
                       plot_pairs_name = "BLBR_SIDER_pairs_plot",
                       sup_xy = TRUE,           
                       plot_xy_save_pdf = FALSE,
                       plot_xy_name = "BLBR_SIDER_xy_plot",
                       gelman = TRUE,
                       heidel = FALSE,  
                       geweke = TRUE,   
                       diag_save = TRUE,
                       diag_name = "BLBR_SIDER_diagnostics",
                       indiv_effect = FALSE,       
                       plot_post_save_png = FALSE, 
                       plot_pairs_save_png = FALSE,
                       plot_xy_save_png = FALSE)
##############################


##############################
# Create diagnostics, summary statistics, and posterior plots
##############################
output_JAGS(BLBR_SIDER, mix.BLBR.SIDER, source_BLBR_SIDER, output_options)
print(BLBR_SIDER)
getwd()
save.image("BLBR_SIDER_Jags.RData") 
##############################


##################################################
#Put data into a dataframe
##################################################
setwd("~/Isotopic_Routing_Small_Mammals/Mixing_Models/BLBR")

require(R2jags) # Load packages 

# Load data 
load("~/Isotopic_Routing_Small_Mammals/Mixing_Models/BLBR/BLBR_SIDER_Jags.RData")
attach.jags(BLBR_SIDER) # call model file
print(BLBR_SIDER)

#assign names to each food source https://github.com/brianstock/MixSIAR/issues/30
BLBR_SIDER_Post <- data.frame(Source = "SIDER", AM_Fungi = p.global[,1],Arthropods = p.global[,2],
                              Berries = p.global[,3], EM_Fungi = p.global[,4], Red_Maple = p.global[,5], 
                              Fungi = p.global[,1]+p.global[,4]) 

BLBR_SIDER_Post<-data.frame(Number=rownames(BLBR_SIDER_Post), BLBR_SIDER_Post)#add posterior number
require(tidyr)
BLBR_SIDER_Posterior <- BLBR_SIDER_Post %>% gather(Source,value,AM_Fungi:Fungi)#Use tidyr to gather posterior for each food source
BLBR_SIDER_Posterior["Model"] <- "SIDER"#Add model type
BLBR_SIDER_Posterior["Species"] <- "BLBR"#Add Species

#save raw posterior data
write.csv(BLBR_SIDER_Posterior, file = "BLBR_SIDER_Posterior.csv", row.names=FALSE)#save raw posterior data
write.csv(BLBR_SIDER_Posterior, file = "~/Isotopic_Routing_Small_Mammals/Mixing_Models/Mixing_Model_Figure/BLBR_SIDER_Posterior.csv", row.names=FALSE)
################################################################################











################################################################################
#SIDER_update
################################################################################
setwd("~/Isotopic_Routing_Small_Mammals/Mixing_Models/BLBR")#set working directory
#############################
#load isotopic data
#############################
#Take stomach samples, convert to SIDER_update, and select species
Hair<- read.csv("Bartlett_isotope_hair_data.csv",header=T)
library(dplyr)
BLBR_SIDER_update<-filter(Hair, Species=="BLBR")
write.csv(BLBR_SIDER_update, "BLBR_SIDER_update.csv")

library(MixSIAR)
mix.filename.BLBR_SIDER_update <-("~/Isotopic_Routing_Small_Mammals/Mixing_Models/BLBR/BLBR_SIDER_update.csv") #load isotopic data for MixSIAR

mix.BLBR.SIDER_update <- load_mix_data(filename=mix.filename.BLBR_SIDER_update, #name of the CSV file with mix/consumer data
                                       iso_names=c("d13C","d15N"), #tracers/isotopes
                                       factors=c(NULL),            #factors
                                       fac_random=c(NULL),         #random effects
                                       fac_nested=c(NULL),         #Nested effects
                                       cont_effects=NULL)          #continuous effects
#############################


#############################
# Load food source data
#############################
source.filename <- ("~/Isotopic_Routing_Small_Mammals/Mixing_Models/BLBR/Bartlett_diet_items.csv")

source_BLBR_SIDER_update <- load_source_data(filename=source.filename, #name of the CSV file with food source data
                                             source_factors=NULL,     #column headings of random/fixed effects you have source data by
                                             conc_dep=TRUE,           #TRUE or FALSE for concentration dependence data 
                                             data_type="raw",         #"means" = means + sd or "raw" = raw data
                                             mix.BLBR.SIDER_update)

# Load discrimination data
#Select species and method of estimated TDF
Est_DFT<- read.csv("Estimated_TDF.csv",header=T)#discrimination factors
BLBR_SIDER_update_d13C_TDF<-filter(Est_DFT,Species_Abr == "BLBR" & Method =="SIDER_update" & Isotope == "d13C")
BLBR_SIDER_update_d15N_TDF<-filter(Est_DFT,Species_Abr == "BLBR" & Method =="SIDER_update" & Isotope == "d15N")


Type<- c("AM Fungi","Arthropods", "Red Maple", "EM Fungi", "Berries")
TDF_BLBR_SIDER_update <- as.data.frame(Type)
TDF_BLBR_SIDER_update$Meand13C<-BLBR_SIDER_update_d13C_TDF$TDF_Mean	
TDF_BLBR_SIDER_update$SDd13C<-BLBR_SIDER_update_d13C_TDF$TDF_SD			
TDF_BLBR_SIDER_update$Meand15N<-BLBR_SIDER_update_d15N_TDF$TDF_Mean			
TDF_BLBR_SIDER_update$SDd15N<-BLBR_SIDER_update_d15N_TDF$TDF_SD		
write.csv(TDF_BLBR_SIDER_update, "BLBR_Discrimination_SIDER_update.csv", row.names=FALSE)

BLBR_Discrimination_SIDER_update <- ("~/Isotopic_Routing_Small_Mammals/Mixing_Models/BLBR/BLBR_Discrimination_SIDER_update.csv")
BLBR_Discrimination_SIDER_update <- load_discr_data(filename=BLBR_Discrimination_SIDER_update, mix.BLBR.SIDER_update)
#############################


#############################
# Make isospace plot
#############################
plot_data(filename="BLBR_SIDER_update_isospace", #name to save the isospace plot as
          plot_save_pdf=FALSE,       #TRUE or FALSE, should mix.SIAR save the plot as a .pdf?
          plot_save_png=TRUE,      #TRUE or FALSE, should mix.SIAR save the plot as a .png?
          mix.BLBR.SIDER_update,source_BLBR_SIDER_update,BLBR_Discrimination_SIDER_update)

# If 2 isotopes/tracers, calculate normalized surface area of the convex hull polygon(s)
#   *Note 1: discrimination SD is added to the source SD (see calc_area.r for details)
if(mix.BLBR.SIDER_update$n.iso==2) calc_area(source=source_BLBR_SIDER_update,mix=mix.BLBR.SIDER_update,discr=BLBR_Discrimination_SIDER_update)
#############################


#############################
# Define prior
#############################
kw.alpha <- c(1,1,1,0,1)

kw.alpha <- kw.alpha*length(kw.alpha)/sum(kw.alpha)# Generate alpha hyperparameters scaling sum(alpha)=n.sources

kw.alpha[which(kw.alpha==0)] <- 0.01# the Dirichlet hyperparameters for the alpha.prior cannot be 0 (but can set = .01)

#plot using "plot_prior"
plot_prior(alpha.prior=kw.alpha,
           source=source_BLBR_SIDER_update,
           plot_save_pdf=FALSE,
           plot_save_png=TRUE,
           filename="BLBR_SIDER_update_prior_plot_kw_inf")
#   RED = your prior
#   DARK GREY = "uninformative"/generalist (alpha = 1)
#   LIGHT GREY = "uninformative" Jeffrey's prior (alpha = 1/n.sources)
#############################


#############################
# Write JAGS model file (define model structure)
##############################
# There are 3 error term options available:
#   1. Residual * Process (resid_err = TRUE, process_err = TRUE)
#   2. Residual only (resid_err = TRUE, process_err = FALSE)
#   3. Process only (resid_err = FALSE, process_err = TRUE)
BLBR_SIDER_update_model <- "BLBR_SIDER_update_model.txt"
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(BLBR_SIDER_update_model, resid_err, process_err, mix.BLBR.SIDER_update, source_BLBR_SIDER_update)
##############################


##############################
# Run BLBR SIDER_update model
##############################
# MCMC options:
#run <- list(chainLength=30000, burn=15000, thin=25, chains=3, calcDIC=TRUE)
costum_run <- list(chainLength=200000, burn=50000, thin=50, chains=3, calcDIC=TRUE)

#Model
BLBR_SIDER_update<- run_model(run=costum_run,mix.BLBR.SIDER_update,source_BLBR_SIDER_update,BLBR_Discrimination_SIDER_update,BLBR_SIDER_update_model,
                              alpha.prior=kw.alpha, resid_err, process_err)
##############################


##############################
# Process JAGS output
##############################
# Choose output options (see ?output_options for details)
output_options <- list(summary_save = TRUE,                 
                       summary_name = "BLBR_SIDER_update_summary_statistics", 
                       sup_post = FALSE,                    
                       plot_post_save_pdf = TRUE,           
                       plot_post_name = "BLBR_SIDER_update_posterior_density",
                       sup_pairs = FALSE,             
                       plot_pairs_save_pdf = TRUE,    
                       plot_pairs_name = "BLBR_SIDER_update_pairs_plot",
                       sup_xy = TRUE,           
                       plot_xy_save_pdf = FALSE,
                       plot_xy_name = "BLBR_SIDER_update_xy_plot",
                       gelman = TRUE,
                       heidel = FALSE,  
                       geweke = TRUE,   
                       diag_save = TRUE,
                       diag_name = "BLBR_SIDER_update_diagnostics",
                       indiv_effect = FALSE,       
                       plot_post_save_png = FALSE, 
                       plot_pairs_save_png = FALSE,
                       plot_xy_save_png = FALSE)
##############################


##############################
# Create diagnostics, summary statistics, and posterior plots
##############################
output_JAGS(BLBR_SIDER_update, mix.BLBR.SIDER_update, source_BLBR_SIDER_update, output_options)
print(BLBR_SIDER_update)
getwd()
save.image("BLBR_SIDER_update_Jags.RData") 
##############################


##################################################
#Put data into a dataframe
##################################################
setwd("~/Isotopic_Routing_Small_Mammals/Mixing_Models/BLBR")

require(R2jags) # Load packages 

# Load data 
load("~/Isotopic_Routing_Small_Mammals/Mixing_Models/BLBR/BLBR_SIDER_update_Jags.RData")
attach.jags(BLBR_SIDER_update) # call model file
print(BLBR_SIDER_update)

#assign names to each food source https://github.com/brianstock/MixSIAR/issues/30
BLBR_SIDER_update_Post <- data.frame(Source = "SIDER_update", AM_Fungi = p.global[,1],Arthropods = p.global[,2],
                                     Berries = p.global[,3], EM_Fungi = p.global[,4], Red_Maple = p.global[,5], 
                                     Fungi = p.global[,1]+p.global[,4]) 

BLBR_SIDER_update_Post<-data.frame(Number=rownames(BLBR_SIDER_update_Post), BLBR_SIDER_update_Post)#add posterior number
require(tidyr)
BLBR_SIDER_update_Posterior <- BLBR_SIDER_update_Post %>% gather(Source,value,AM_Fungi:Fungi)#Use tidyr to gather posterior for each food source
BLBR_SIDER_update_Posterior["Model"] <- "SIDER_update"#Add model type
BLBR_SIDER_update_Posterior["Species"] <- "BLBR"#Add Species

#save raw posterior data
write.csv(BLBR_SIDER_update_Posterior, file = "BLBR_SIDER_update_Posterior.csv", row.names=FALSE)#save raw posterior data
write.csv(BLBR_SIDER_update_Posterior, file = "~/Isotopic_Routing_Small_Mammals/Mixing_Models/Mixing_Model_Figure/BLBR_SIDER_update_Posterior.csv", row.names=FALSE)
################################################################################









################################################################################
#Meta_analysis
################################################################################
setwd("~/Isotopic_Routing_Small_Mammals/Mixing_Models/BLBR")#set working directory
#############################
#load isotopic data
#############################
#Take stomach samples, convert to Meta_analysis, and select species
Hair<- read.csv("Bartlett_isotope_hair_data.csv",header=T)
library(dplyr)
BLBR_Meta_analysis<-filter(Hair, Species=="BLBR")
write.csv(BLBR_Meta_analysis, "BLBR_Meta_analysis.csv")

library(MixSIAR)
mix.filename.BLBR_Meta_analysis <-("~/Isotopic_Routing_Small_Mammals/Mixing_Models/BLBR/BLBR_Meta_analysis.csv") #load isotopic data for MixSIAR

mix.BLBR.Meta_analysis <- load_mix_data(filename=mix.filename.BLBR_Meta_analysis, #name of the CSV file with mix/consumer data
                                        iso_names=c("d13C","d15N"), #tracers/isotopes
                                        factors=c(NULL),            #factors
                                        fac_random=c(NULL),         #random effects
                                        fac_nested=c(NULL),         #Nested effects
                                        cont_effects=NULL)          #continuous effects
#############################


#############################
# Load food source data
#############################
source.filename <- ("~/Isotopic_Routing_Small_Mammals/Mixing_Models/BLBR/Bartlett_diet_items.csv")

source_BLBR_Meta_analysis <- load_source_data(filename=source.filename, #name of the CSV file with food source data
                                              source_factors=NULL,     #column headings of random/fixed effects you have source data by
                                              conc_dep=TRUE,           #TRUE or FALSE for concentration dependence data 
                                              data_type="raw",         #"means" = means + sd or "raw" = raw data
                                              mix.BLBR.Meta_analysis)

# Load discrimination data
#Select species and method of estimated TDF
Est_DFT<- read.csv("Estimated_TDF.csv",header=T)#discrimination factors
BLBR_Meta_analysis_d13C_TDF<-filter(Est_DFT,Species_Abr == "BLBR" & Method =="Meta_analysis" & Isotope == "d13C")
BLBR_Meta_analysis_d15N_TDF<-filter(Est_DFT,Species_Abr == "BLBR" & Method =="Meta_analysis" & Isotope == "d15N")


Type<- c("AM Fungi","Arthropods", "Red Maple", "EM Fungi", "Berries")
TDF_BLBR_Meta_analysis <- as.data.frame(Type)
TDF_BLBR_Meta_analysis$Meand13C<-BLBR_Meta_analysis_d13C_TDF$TDF_Mean	
TDF_BLBR_Meta_analysis$SDd13C<-BLBR_Meta_analysis_d13C_TDF$TDF_SD			
TDF_BLBR_Meta_analysis$Meand15N<-BLBR_Meta_analysis_d15N_TDF$TDF_Mean			
TDF_BLBR_Meta_analysis$SDd15N<-BLBR_Meta_analysis_d15N_TDF$TDF_SD		
write.csv(TDF_BLBR_Meta_analysis, "BLBR_Discrimination_Meta_analysis.csv", row.names=FALSE)

BLBR_Discrimination_Meta_analysis <- ("~/Isotopic_Routing_Small_Mammals/Mixing_Models/BLBR/BLBR_Discrimination_Meta_analysis.csv")
BLBR_Discrimination_Meta_analysis <- load_discr_data(filename=BLBR_Discrimination_Meta_analysis, mix.BLBR.Meta_analysis)
#############################


#############################
# Make isospace plot
#############################
plot_data(filename="BLBR_Meta_analysis_isospace", #name to save the isospace plot as
          plot_save_pdf=FALSE,       #TRUE or FALSE, should mix.SIAR save the plot as a .pdf?
          plot_save_png=TRUE,      #TRUE or FALSE, should mix.SIAR save the plot as a .png?
          mix.BLBR.Meta_analysis,source_BLBR_Meta_analysis,BLBR_Discrimination_Meta_analysis)

# If 2 isotopes/tracers, calculate normalized surface area of the convex hull polygon(s)
#   *Note 1: discrimination SD is added to the source SD (see calc_area.r for details)
if(mix.BLBR.Meta_analysis$n.iso==2) calc_area(source=source_BLBR_Meta_analysis,mix=mix.BLBR.Meta_analysis,discr=BLBR_Discrimination_Meta_analysis)
#############################


#############################
# Define prior
#############################
kw.alpha <- c(1,1,1,0,1)

kw.alpha <- kw.alpha*length(kw.alpha)/sum(kw.alpha)# Generate alpha hyperparameters scaling sum(alpha)=n.sources

kw.alpha[which(kw.alpha==0)] <- 0.01# the Dirichlet hyperparameters for the alpha.prior cannot be 0 (but can set = .01)

#plot using "plot_prior"
plot_prior(alpha.prior=kw.alpha,
           source=source_BLBR_Meta_analysis,
           plot_save_pdf=FALSE,
           plot_save_png=TRUE,
           filename="BLBR_Meta_analysis_prior_plot_kw_inf")
#   RED = your prior
#   DARK GREY = "uninformative"/generalist (alpha = 1)
#   LIGHT GREY = "uninformative" Jeffrey's prior (alpha = 1/n.sources)
#############################


#############################
# Write JAGS model file (define model structure)
##############################
# There are 3 error term options available:
#   1. Residual * Process (resid_err = TRUE, process_err = TRUE)
#   2. Residual only (resid_err = TRUE, process_err = FALSE)
#   3. Process only (resid_err = FALSE, process_err = TRUE)
BLBR_Meta_analysis_model <- "BLBR_Meta_analysis_model.txt"
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(BLBR_Meta_analysis_model, resid_err, process_err, mix.BLBR.Meta_analysis, source_BLBR_Meta_analysis)
##############################


##############################
# Run BLBR Meta_analysis model
##############################
# MCMC options:
run <- list(chainLength=30000, burn=15000, thin=25, chains=3, calcDIC=TRUE)
costum_run <- list(chainLength=200000, burn=50000, thin=50, chains=3, calcDIC=TRUE)

#Model
BLBR_Meta_analysis<- run_model(run=costum_run,mix.BLBR.Meta_analysis,source_BLBR_Meta_analysis,BLBR_Discrimination_Meta_analysis,BLBR_Meta_analysis_model,
                               alpha.prior=kw.alpha, resid_err, process_err)
##############################


##############################
# Process JAGS output
##############################
# Choose output options (see ?output_options for details)
output_options <- list(summary_save = TRUE,                 
                       summary_name = "BLBR_Meta_analysis_summary_statistics", 
                       sup_post = FALSE,                    
                       plot_post_save_pdf = TRUE,           
                       plot_post_name = "BLBR_Meta_analysis_posterior_density",
                       sup_pairs = FALSE,             
                       plot_pairs_save_pdf = TRUE,    
                       plot_pairs_name = "BLBR_Meta_analysis_pairs_plot",
                       sup_xy = TRUE,           
                       plot_xy_save_pdf = FALSE,
                       plot_xy_name = "BLBR_Meta_analysis_xy_plot",
                       gelman = TRUE,
                       heidel = FALSE,  
                       geweke = TRUE,   
                       diag_save = TRUE,
                       diag_name = "BLBR_Meta_analysis_diagnostics",
                       indiv_effect = FALSE,       
                       plot_post_save_png = FALSE, 
                       plot_pairs_save_png = FALSE,
                       plot_xy_save_png = FALSE)
##############################


##############################
# Create diagnostics, summary statistics, and posterior plots
##############################
output_JAGS(BLBR_Meta_analysis, mix.BLBR.Meta_analysis, source_BLBR_Meta_analysis, output_options)
print(BLBR_Meta_analysis)
getwd()
save.image("BLBR_Meta_analysis_Jags.RData") 
##############################


##################################################
#Put data into a dataframe
##################################################
setwd("~/Isotopic_Routing_Small_Mammals/Mixing_Models/BLBR")

require(R2jags) # Load packages 

# Load data 
load("~/Isotopic_Routing_Small_Mammals/Mixing_Models/BLBR/BLBR_Meta_analysis_Jags.RData")
attach.jags(BLBR_Meta_analysis) # call model file
print(BLBR_Meta_analysis)

#assign names to each food source https://github.com/brianstock/MixSIAR/issues/30
BLBR_Meta_analysis_Post <- data.frame(Source = "Meta_analysis", AM_Fungi = p.global[,1],Arthropods = p.global[,2],
                                      Berries = p.global[,3], EM_Fungi = p.global[,4], Red_Maple = p.global[,5], 
                                      Fungi = p.global[,1]+p.global[,4]) 

BLBR_Meta_analysis_Post<-data.frame(Number=rownames(BLBR_Meta_analysis_Post), BLBR_Meta_analysis_Post)#add posterior number

require(tidyr)
BLBR_Meta_analysis_Posterior <- BLBR_Meta_analysis_Post %>% gather(Source,value,AM_Fungi:Fungi)#Use tidyr to gather posterior for each food source
BLBR_Meta_analysis_Posterior["Model"] <- "Meta_analysis"#Add model type
BLBR_Meta_analysis_Posterior["Species"] <- "BLBR"#Add Species

#save raw posterior data
write.csv(BLBR_Meta_analysis_Posterior, file = "BLBR_Meta_analysis_Posterior.csv", row.names=FALSE)#save raw posterior data
write.csv(BLBR_Meta_analysis_Posterior, file = "~/Isotopic_Routing_Small_Mammals/Mixing_Models/Mixing_Model_Figure/BLBR_Meta_analysis_Posterior.csv", row.names=FALSE)
################################################################################

















