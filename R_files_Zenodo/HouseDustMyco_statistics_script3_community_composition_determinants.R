
#########################Filename: HouseDustMyco_statistics_script3_community_composition_determinants.R

#Publication: Analyzing indoor mycobiomes through a large-scale citizen science study in Norway
#Authors: Pedro M. Martin-Sanchez, Eva Lena F. Estensmo, Luis N. Morgado, Sundy Maurice, Ingeborg B. Engh, Inger Skrede and Håvard Kauserud

##Related to the sections: 
        #Materials and Methods- 2.5. Statistical analyses
        #Results- 3.2. Determinants of the community composition
        #Supplemental Information

## Goals: 
#(1) To compare the fungal community composition in the house-dust samples - NMDS ordination plots 
#(2) To calculate the contribution of the different variables to explain the compositional variation of the mycobiomes   

#Inputs: 
        #The rarefied fungal OTU table ("OTU_table_rar2000x10.xlsx")
        #The metadata file for dust samples ("metadata_file_dust_samples.xlsx")
        #The taxonomy file for the OTUs in the rarefied table ("taxonomy_for_rarTable_funguild.xlsx")

#Outputs: 
        #Supplemental file 2 - Table S4 (Top 200 most abundant OTUs)
        #Plots for the Figure 3 and S8 (NMDS ordinations)
        #Data for the Tables 1 and S3 (PERMANOVA -Adonis2) 
        #Plots for the Figure 4 (VPA) 

#########################
#PREPARATIONS
#########################

#Load packages
library(openxlsx)
library(tidyverse)
library(vegan)
library(metagenomeSeq)
library(ggrepel)
library(ggordiplots)


#Set your path
setwd("C:/WP1_data_analyses/script_preparation/")

#Load the rarefied OTU table (6,632 OTUs from 807 dust samples)
OTU_tab_rar=read.xlsx("OTU_table_rar2000x10.xlsx")

#Different format after moving the column "sample" to the row names
OTU_tab_rar_for_analyses=OTU_tab_rar %>% remove_rownames() %>% column_to_rownames("sample")

#Load the metadata
metadata=read.xlsx("metadata_file_dust_samples.xlsx")
#Remark: Some geographic variables(latitude, longitude, municipality, county and region) as well as the construction year of the study houses
# were excluded from this file to protect personal data of volunteers   

#Load the taxonomy file for the 6,632 OTUs in the rarefied table, including the FUNGuild annotations
taxonomy=read.xlsx("taxonomy_for_rarTable_funguild.xlsx")

#########################
# Master table for analyses
#########################

#Merge all metadata and taxonomic assigments with the OTU table
master_table=OTU_tab_rar %>% 
        gather(OTU_id, rarRead, -sample) %>% 
        left_join(taxonomy) %>% 
        left_join(metadata)


#########################
# Some initial calculations on abundance and distribution of OTUs
#########################

#Relative Abundances of OTUs in the complete dataset  
RA_OTU_all= OTU_tab_rar_for_analyses %>% 
        colSums() %>% 
        as.data.frame() %>% 
        rename(sum_rarRead =".") %>% 
        tibble::rownames_to_column("OTU_id")%>%
        mutate (RA = sum_rarRead /sum(sum_rarRead))

#OTU distribution (number of samples where is present) in the complete dataset
OTU_distribution_all= OTU_tab_rar_for_analyses %>% 
        decostand("pa") %>%
        colSums() %>% 
        as.data.frame() %>% 
        rename(number_of_samples =".") %>% 
        rownames_to_column("OTU_id")

#Number of houses where the OTUs are present
number_of_houses_OTUS_all<-master_table %>% 
        mutate(PA=sign(rarRead)) %>% 
        group_by(OTU_id,house) %>% 
        summarise(location_per_house=sum(PA)) %>% 
        filter(location_per_house>0)%>%
        mutate(PA_houses=sign(house)) %>%
        group_by(OTU_id) %>% 
        summarise(number_of_houses=sum(PA_houses))

##### Prepare the Excel table for the Top 200 most abundant OTUS (basis for the Supplemental file 2 - Table S4)
RA_OTU_top_200=RA_OTU_all %>%
        top_n(200)

Top200_OTUS_table = RA_OTU_top_200 %>%  left_join(OTU_distribution_all)%>%
        left_join(number_of_houses_OTUS_all) %>%left_join(taxonomy)

#Save excel with Top-200 most abundant OTUS
write.xlsx(Top200_OTUS_table, "Top200_OTUS_table.xlsx")


#########################
# Transformation of the rarefied OTU table
#########################
#We initially tried three different transformations: Logarithmic (Log), Hellinger (Hel) and Cumulative Sum Scaling (CSS)
#After a preliminary evaluation on the complete dataset (ordinations for ALL dust samples), we selected the Hellinger transformation for further analyses 

#Logarithmic
otutable_log = OTU_tab_rar_for_analyses %>% decostand("log")

#Helinger
otutable_hel = OTU_tab_rar_for_analyses %>% decostand("hellinger")

##CSS transformation
NORM_otutable = OTU_tab_rar_for_analyses %>%
        t() %>% 
        as.data.frame() %>% 
        newMRexperiment()

p_no_outliers=cumNormStatFast(NORM_otutable)
NORM_otutable2<-cumNorm(NORM_otutable,p_no_outliers)

otutable_CSS <-MRcounts (NORM_otutable2, norm = TRUE, log = T) %>%
        t() %>%
        as.data.frame()


##################################################
# 1. Non-metric Multidimensional Scaling (NMDS) analyses 
#########################
#To compare fungal community composition of the house-dust samples
#They allowed us to obtain the ordination plots for both samples and species (OTUs)
#NMDS analyses were done for 3 datasets:
        #ALL dust samples
        #INDOOR dust samples
        #OUTDOOR dust samples


#########################
# 1.1. NMDS - ALL SAMPLES 
#########################

# Prepare metadata
all_metadata=master_table %>% 
        select(-c(2:19)) %>%
        distinct(sample, .keep_all = T)

#Ordinations for the complete rarefied table comparing: (i) without any transformation, (ii) Log, (iii) Hel and (iv) CSS. 
all_ord <-metaMDS(OTU_tab_rar_for_analyses, k=2, autotransform = F, maxit=200, smin = 1e-7, sfgrmin = 1e-7, try = 200) # 14: no. of iterations >= maxit; 186: stress ratio > sratmax
all_ord_log <-metaMDS(otutable_log, k=2, autotransform = F, maxit=200, smin = 1e-7, sfgrmin = 1e-7, try = 200) # 85: no. of iterations >= maxit; 115: stress ratio > sratmax
all_ord_hel <-metaMDS(otutable_hel, k=2, autotransform = F, maxit=200, smin = 1e-7, sfgrmin = 1e-7, try = 200)   # 87: no. of iterations >= maxit; 113: stress ratio > sratmax
all_ord_CSS <-metaMDS(otutable_CSS, k=2, autotransform = F, maxit=200, smin = 1e-7, sfgrmin = 1e-7, try = 200) # 61: no. of iterations >= maxit; 139: stress ratio > sratmax

#Check information of the MDS objects generated. No convergent solution in any of the cases.
all_ord_hel
#Call:
#        metaMDS(comm = otutable_hel, k = 2, try = 200, autotransform = F,      maxit = 200, smin = 1e-07, sfgrmin = 1e-07) 
#global Multidimensional Scaling using monoMDS
#Data:     otutable_hel 
#Distance: bray 
#Dimensions: 2 
#Stress:     0.2131583 
#Stress type 1, weak ties
#No convergent solutions - best solution after 200 tries
#Scaling: centring, PC rotation, halfchange scaling 
#Species: expanded scores based on 'otutable_hel

all_ord
#Stress:     0.2283401 

all_ord_log
#Stress:     0.2138411

all_ord_CSS
#Stress:     0.2176118 

## Save ordinations (metaMDS objects) / read (load) them
#saveRDS(all_ord_hel, file ="all_ord_hel.rds")
all_ord_hel_final <- readRDS("C:/WP1_data_analyses/October/all_ord_hel.rds")

##After the comparison of the different ordinations, we decided to use the rarefied Hellinger-transformed OTU table  for further analyses, 
# as the ordination from the Hel option reported the lower stress value    -> Why???

#Extract ordination sites for SAMPLES and merge with all metadata
all_ord_sites_metadata_hel=all_ord_hel$points %>% as.data.frame() %>% rownames_to_column("sample") %>% left_join(all_metadata)

#Ordination plot for ALL dust samples
all_ord_sites_metadata_hel %>% ggplot(aes(x=MDS1,y=MDS2))+
        geom_point (aes(fill = location), shape=21, color="black",size=4) +
        scale_fill_manual(values = c("orange","blue","darkgrey"))+
        labs (title = "ALL samples", subtitle="House compartment")+
        theme(legend.title = element_text(colour="black", size=10, 
                                          face="bold"))+
        theme_bw()+
        guides(fill = guide_legend(reverse=TRUE))
####Remark: This plot was used in Figure 3a


#Plotting the 200 most abundant OTUS adding the labels (species names) to the top 20

RA_OTU_top_20=RA_OTU_all %>% 
        top_n(20) %>%
        left_join(taxonomy) %>% 
        select(OTU_id, species) %>%
        mutate (OTU_label= OTU_id)%>%
        unite(species_id,c("OTU_label","species"))

#Extract the ordination sites for SPECIES ( Top 200 OTUs)
all_ord_species_tax_RA_top200=all_ord_hel$species %>% as.data.frame() %>% rownames_to_column("OTU_id") %>% 
        left_join(taxonomy) %>% right_join(RA_OTU_top_200) %>% left_join(RA_OTU_top_20)

#Plot coloring the different Phylla for the Top 200 OTUs
all_ord_species_tax_RA_top200 %>% 
        filter(!is.na(phyllum)) %>% 
        ggplot(aes(x=MDS1,y=MDS2))+
        geom_point(aes(size=RA,fill=phyllum), shape=21, col= "black", stroke = 0.2, alpha= 1/1.5)+
        scale_fill_manual(values = c("grey40","red", "yellow"))+
        scale_size_area(max_size = 30)+
        labs (title = "ALL SAMPLES", subtitle="top 200 most abundant OTUs")+
        xlim(-2.1, 2.1)+ ylim(-1.8, 2.1)+
        scale_color_manual(values = c("black","red", "yellow"))+
        theme_bw()+
        labs(col="Phyllum", size="Relative abundance")+
        theme (panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
        geom_label_repel(aes(label=species_id, col=phyllum),size= 3, na.rm = T, box.padding = 3, label.size = 0.2, segment.size = 0.2)
####Remark: This plot was used in Figure 3c


#########################
# 1.2. NMDS - INDOOR samples 
#########################

#Filtering the data matrix to get the OTU table for indoor samples
data_for_ord_indoor=master_table %>% 
        filter(indoor_outdoor=="indoor")

indoor_otutable=data_for_ord_indoor %>% 
        select(sample, OTU_id, rarRead) %>% 
        group_by(OTU_id) %>% 
        filter(sum(rarRead)>0) %>% 
        ungroup() %>% 
        spread(OTU_id, rarRead) %>% 
        column_to_rownames('sample')

#Hellinger transformation
indoor_otutable_hel = indoor_otutable %>% decostand("hel")

#Select metadata for indoor samples
indoor_metadata=data_for_ord_indoor %>% 
        select(-c(2:19)) %>%
        distinct(sample, .keep_all = T)

#Ordination for the Hellinger-transformed table
indoor_ord_hel <-metaMDS(indoor_otutable_hel, k=2, autotransform = F, maxit=200, smin = 1e-7, sfgrmin = 1e-7, try = 200) 

#Check the resulting ordination
indoor_ord_hel
# Remarks: No convergent solutions - best solution after 200 tries,
#Stress:     0.2441975 

#Extract ordination sites for SAMPLES and merge with all metadata
in_ord_sites_metadata_hel=indoor_ord_hel$points %>% as.data.frame() %>% rownames_to_column("sample") %>% left_join(indoor_metadata)

#Reorder the categories of the variable region for a better plot

#Remark: NOTE THAT THE VARIABLE "REGION" HAS BEEN EXCLUDED OF THE METADATA FILE
# In order to protect the personal data of the volunteers
# Therefore, THIS ANALYSIS (to get the Figure S8c) CANNOT BE RUN!!

in_ord_sites_metadata_hel$region <- factor(in_ord_sites_metadata_hel$region,
                                            levels = c('West','East','South', 'Mid', 'North','Svalbard'), ordered = TRUE)

#Set a color blind friendly palette for the regions
cbf_3 <- c("#999999", "#E69F00", "#56B4E9","#0072B2" , 
           "#F0E442", "#009E73", "#D55E00", "#CC79A7")

#Ordination plot for INDOOR dust samples
in_ord_sites_metadata_hel %>% ggplot(aes(x=MDS1,y=MDS2,col=region))+
        geom_point (aes(fill = region), shape=21, color="black",size=5,  stroke = 0.2) +
        scale_fill_manual(values = cbf_3)+
        theme_bw()+
        theme (panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
        labs (title = "INDOOR samples", subtitle="Geographical regions")
####Remark: This plot was used in Figure S8c 


#Plotting the 200 most abundant OTUS adding the labels (species names) to the top 20

#Relative AbundanceS for the indoor OTUS 
RA_OTU_indoor= indoor_otutable %>% 
        colSums() %>% 
        as.data.frame() %>% 
        rename(sum_rarRead =".") %>% 
        tibble::rownames_to_column("OTU_id")%>%
        mutate (RA = sum_rarRead /sum(sum_rarRead))

#Filter the top 200
RA_OTU_top_200_indoor=RA_OTU_indoor %>% 
        top_n(200)

#Filter the top 20 for adding the labels (species names)
RA_OTU_top_20_indoor=RA_OTU_indoor %>% # to put labels of the Top20 species
        top_n(20) %>% 
        left_join(taxonomy) %>% 
        select(OTU_id, species) %>%
        mutate (OTU_label= OTU_id)%>%
        unite(species_id,c("OTU_label","species"))

#Extract the ordination sites for SPECIES (Top 200 OTUs)
indoor_ord_species_tax_RA_top200=indoor_ord_hel$species %>% as.data.frame() %>% rownames_to_column("OTU_id") %>% 
        left_join(taxonomy) %>% right_join(RA_OTU_top_200_indoor) %>% left_join(RA_OTU_top_20_indoor)


#Plot coloring the different Phylla for the Top 200 OTUs
indoor_ord_species_tax_RA_top200 %>% 
        filter(!is.na(phyllum)) %>% 
        ggplot(aes(x=MDS1,y=MDS2))+
        geom_point(aes(size=RA,fill=phyllum), shape=21, col= "black", stroke = 0.2, alpha= 1/1.5)+
        scale_fill_manual(values = c("grey40","red", "yellow", "light blue", "blue"))+
        scale_size_area(max_size = 30)+
        labs (title = "INDOOR", subtitle="top 200 most abundant OTUs")+
        xlim(-1, 1.8)+ ylim(-1, 1.7)+
        scale_color_manual(values = c("black","red", "yellow", "light blue", "blue"))+
        theme_bw()+
        labs(col="Phyllum", size="Relative abundance")+
        geom_label_repel(aes(label=species_id, col=phyllum),size= 3, na.rm = T, box.padding = 3, label.size = 0.2, segment.size = 0.2)
####Remark: This plot was used in Figure S8d


#########################
# 1.3. NMDS - OUTDOOR samples 
#########################

#Filtering the data matrix to get the OTU table for outdoor samples
data_for_ord_outdoor=master_table %>% 
        filter(indoor_outdoor=="outdoor")

outdoor_otutable=data_for_ord_outdoor %>% 
        select(sample, OTU_id, rarRead) %>% 
        group_by(OTU_id) %>% 
        filter(sum(rarRead)>0) %>% 
        ungroup() %>% 
        spread(OTU_id, rarRead) %>% 
        column_to_rownames('sample')

#Hellinger transformation
outdoor_otutable_hel = outdoor_otutable %>% decostand("hel")

#Select metadata for outdoor samples
outdoor_metadata=data_for_ord_outdoor %>% 
        select(-c(2:19)) %>%
        distinct(sample, .keep_all = T)

#Ordination for the Hellinger-transformed table
outdoor_ord_hel <-metaMDS(outdoor_otutable_hel, k=2, autotransform = F, maxit=200, smin = 1e-7, sfgrmin = 1e-7, try = 200) # 15: no. of iterations >= maxit; 185: stress ratio > sratmax

#Check the resulting ordination
outdoor_ord_hel
# Remarks: No convergent solutions - best solution after 200 tries,
#Stress:     0.256098 

#Extract ordination sites for SAMPLES and merge with all metadata 
out_ord_sites_metadata_hel=outdoor_ord_hel$points %>% as.data.frame() %>% rownames_to_column("sample") %>% left_join(outdoor_metadata)

#Reorder the categories of the variable region for better plotting

#Remark: NOTE THAT THE VARIABLE "REGION" HAS BEEN EXCLUDED OF THE METADATA FILE
# In order to protect the personal data of the volunteers
# Therefore, THIS ANALYSIS (to get the Figure S8a) CANNOT BE RUN!!

out_ord_sites_metadata_hel$region <- factor(out_ord_sites_metadata_hel$region,
                                           levels = c('West','East','South', 'Mid', 'North','Svalbard'), ordered = TRUE)

#Ordination plot for OUTDOOR dust samples
out_ord_sites_metadata_hel %>% ggplot(aes(x=MDS1,y=MDS2,col=region))+
        geom_point (aes(fill = region), shape=21, color="black",size=5,  stroke = 0.2) +
        scale_fill_manual(values = cbf_3)+
        theme_bw()+
        theme (panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
        labs (title = "OUTDOOR samples", subtitle="Geographical regions")
####Remark: This plot was used in Figure S8a 


#Plotting the 200 most abundant OTUS adding the labels (species names) to the top 20

#Relative AbundanceS for the outdoor OTUS 
RA_OTU_outdoor= outdoor_otutable %>% 
        colSums() %>% 
        as.data.frame() %>% 
        rename(sum_rarRead =".") %>% 
        tibble::rownames_to_column("OTU_id")%>%
        mutate (RA = sum_rarRead /sum(sum_rarRead))

#Filter the top 200
RA_OTU_top_200_outdoor=RA_OTU_outdoor %>% # to plot the Top200 OTUs
        top_n(200)

#Filter the top 20 for adding the labels (species names)
RA_OTU_top_20_outdoor=RA_OTU_outdoor %>%
        top_n(20) %>% 
        left_join(taxonomy) %>% 
        select(OTU_id, species) %>%
        mutate (OTU_label= OTU_id)%>%
        unite(species_id,c("OTU_label","species"))

#Extract the ordination sites for SPECIES (Top 200 OTUs)
outdoor_ord_species_tax_RA_top200=outdoor_ord_hel$species %>% as.data.frame() %>% rownames_to_column("OTU_id") %>% 
        left_join(taxonomy) %>% right_join(RA_OTU_top_200_outdoor) %>% left_join(RA_OTU_top_20_outdoor)

#Plot coloring the different Phylla for the Top 200 OTUs
outdoor_ord_species_tax_RA_top200 %>% 
        filter(!is.na(phyllum)) %>% 
        ggplot(aes(x=MDS1,y=MDS2))+
        geom_point(aes(size=RA,fill=phyllum), shape=21, col= "black", stroke = 0.2, alpha= 1/1.5)+
        scale_fill_manual(values = c("grey40","red", "yellow"))+
        scale_size_area(max_size = 30)+
        labs (title = "OUTDOOR", subtitle="top 200 most abundant OTUs")+
        scale_color_manual(values = c("black","red", "yellow"))+
        theme_bw()+
        labs(col="Phyllum", size="Relative abundance")+
        geom_label_repel(aes(label=species_id, col=phyllum),size= 3, na.rm = T, box.padding = 3, label.size = 0.2, segment.size = 0.2)


##################################################
# 2. Linear-regression of continuous variables with the NMDS plots
#########################
#We used the function gg-envfit from the package ggordiplots
#On the previous NMDS ordinations from the rarefied Hellinger-transformed tables 
#Similar analyses were done for the 3 datasets:
#ALL dust samples
#INDOOR dust samples
#OUTDOOR dust samples


#########################
# 2.1. Linear-regression of continuous variables  - ALL SAMPLES 
#########################

# Collect all continuous variables for the analysis
all_ord_sites_metadata_hel_all_var<-all_ord_sites_metadata_hel %>% 
        select(BIO1, BIO4, BIO9, BIO10, BIO11, BIO12, growing_season_length, swe_4, sca_2,
               solar_radiation, #longitude, latitude,
               people, females, children, dust_coverage)
#Remark: Note that geographic coordinates (latitude and longitude) have been intentionally excluded in this example

#gg_envfit
all_data_all_var_envfit<-gg_envfit(all_ord_hel,all_ord_sites_metadata_hel_all_var, alpha = 1)

#Extract the arrows filtering the significant variables
all_data_all_var_envfit_arrows<-all_data_all_var_envfit$df_arrows  %>% 
        filter(p.val<0.05)
#Plot
all_ord_sites_metadata_hel %>% ggplot(aes(x=MDS1,y=MDS2,col=location))+
        geom_point (aes(fill = location), shape=21, color="black",size=5,  stroke = 0.2) +
        scale_fill_manual(values = c("orange","blue","darkgrey"))+
        geom_segment(data = all_data_all_var_envfit_arrows,
                     aes(x = 0, xend =x, y = 0, yend =y),
                     arrow = arrow(length = unit(0.25, "cm")), colour = "black")+
        geom_label_repel(data = all_data_all_var_envfit_arrows, 
                         aes(x = x, y = y, label=var),
                         size = 4,
                         hjust = 0, col="black")+
        theme_bw()+
        xlim(-2.1, 2.1)+ ylim(-1.8, 2.1)+
        theme (panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
        guides(fill = guide_legend(reverse=TRUE))
labs (title = "ALL SAMPLES_house_compartments", subtitle="All significant continuous variables")

#Get R2 and p-values associated with the vectors
all_gg_envfit_results<-cbind(as.data.frame(all_data_all_var_envfit$plot$plot_env$fit$vectors[1]),as.data.frame(all_data_all_var_envfit$plot$plot_env$fit$vectors[2]),as.data.frame(all_data_all_var_envfit$plot$plot_env$fit$vectors[4])) %>%
        rownames_to_column("EnvVariables")

all_gg_envfit_results2 =all_gg_envfit_results %>% rename(arrows.NMDS1_ALL=arrows.NMDS1, arrows.NMDS2_ALL=arrows.NMDS2, r_ALL=r, pvals_ALL = pvals)


#########################
# 2.2. Linear-regression of continuous variables  - INDOOR SAMPLES 
#########################

# Collect all continuous variables for the analysis
indoor_sites_metadata_hel_for_envfit<-in_ord_sites_metadata_hel %>% 
        select(BIO1, BIO4, BIO9, BIO10, BIO11, BIO12, growing_season_length, swe_4, sca_2,
               solar_radiation, #longitude, latitude,
               people, females, children, dust_coverage)
#gg_envfit
indoor_data_all_var_envfit<-gg_envfit(indoor_ord_hel,indoor_sites_metadata_hel_for_envfit, alpha = 1)

#Extract the arrows filtering the significant variables
indoor_data_all_var_envfit_arrows<-indoor_data_all_var_envfit$df_arrows %>% 
        filter(p.val<0.05)

#reorder the categories in the variable region for better plotting

#Remark: NOTE THAT THE VARIABLE "REGION" HAS BEEN EXCLUDED OF THE METADATA FILE
# In order to protect the personal data of the volunteers
# Therefore, THIS ANALYSIS (to get the Figure S8c) CANNOT BE RUN!!

in_ord_sites_metadata_hel$region <- factor(in_ord_sites_metadata_hel$region,
                                           levels = c('West','East','South', 'Mid', 'North','Svalbard'), ordered = TRUE)

#Plot
in_ord_sites_metadata_hel %>% ggplot(aes(x=MDS1,y=MDS2,col=region))+
        geom_point (aes(fill = region), shape=21, color="black",size=5,  stroke = 0.2) +
        scale_fill_manual(values = cbf_3)+
        geom_segment(data = indoor_data_all_var_envfit_arrows,
                     aes(x = 0, xend =x, y = 0, yend =y),
                     arrow = arrow(length = unit(0.25, "cm")), colour = "black")+
        geom_label_repel(data = indoor_data_all_var_envfit_arrows, 
                         aes(x = x, y = y, label=var),
                         size = 4,
                         hjust = 0, col="black")+
        theme_bw()+
        theme (panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
        labs (title = "INDOOR_region", subtitle="All significant continuous variables")
#####Remark: Plot used in the Figure S8c

#Get R2 and p-values associated with the vectors
indoor_gg_envfit_results<-cbind(as.data.frame(indoor_data_all_var_envfit$plot$plot_env$fit$vectors[1]),as.data.frame(indoor_data_all_var_envfit$plot$plot_env$fit$vectors[2]),as.data.frame(indoor_data_all_var_envfit$plot$plot_env$fit$vectors[4])) %>%
        rownames_to_column("EnvVariables")
indoor_gg_envfit_results2 = indoor_gg_envfit_results %>% rename(arrows.NMDS1_IN=arrows.NMDS1, arrows.NMDS2_IN=arrows.NMDS2, r_IN=r, pvals_IN = pvals)


#########################
# 2.3. Linear-regression of continuous variables  - OUTDOOR SAMPLES 
#########################


# Collect all continuous variables for the analysis
outdoor_ord_sites_metadata_hel_for_envfit<-out_ord_sites_metadata_hel %>%
        select(BIO1, BIO4, BIO9, BIO10, BIO11, BIO12, growing_season_length, swe_4, sca_2,
               solar_radiation, #longitude, latitude,
               people, females, children, dust_coverage)

#gg_envfit        
outdoor_data_all_var_envfit<-gg_envfit(outdoor_ord_hel,outdoor_ord_sites_metadata_hel_for_envfit, alpha = 1)

#Extract the arrows filtering the significant variables
outdoor_data_all_var_envfit_arrows<-outdoor_data_all_var_envfit$df_arrows %>% 
        filter(p.val<0.05)

#reorder the categories in the variable region for better plotting 

#Remark: NOTE THAT THE VARIABLE "REGION" HAS BEEN EXCLUDED OF THE METADATA FILE
# In order to protect the personal data of the volunteers
# Therefore, THIS ANALYSIS (to get the Figure S8a) CANNOT BE RUN!!

out_ord_sites_metadata_hel$region <- factor(out_ord_sites_metadata_hel$region,
                                            levels = c('West','East','South', 'Mid', 'North','Svalbard'), ordered = TRUE)
#Plot
out_ord_sites_metadata_hel %>% ggplot(aes(x=MDS1,y=MDS2,col=region))+
        geom_point (aes(fill = region), shape=21, color="black",size=5,  stroke = 0.2) +
        scale_fill_manual(values = cbf_3)+
        geom_segment(data = outdoor_data_all_var_envfit_arrows,
                     aes(x = 0, xend =x, y = 0, yend =y),
                     arrow = arrow(length = unit(0.25, "cm")), colour = "black")+
        geom_label_repel(data = outdoor_data_all_var_envfit_arrows, 
                         aes(x = x, y = y, label=var),
                         size = 4,
                         hjust = 0, col="black")+
        theme_bw()+
        theme (panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
        labs (title = "OUTDOOR_regions", subtitle="All significant continuous variables")
#####Remark: Plot used in the Figure S8a

#Get R2 and p-values associated with the vectors
outdoor_gg_envfit_results<-cbind(as.data.frame(outdoor_data_all_var_envfit$plot$plot_env$fit$vectors[1]),as.data.frame(outdoor_data_all_var_envfit$plot$plot_env$fit$vectors[2]),as.data.frame(outdoor_data_all_var_envfit$plot$plot_env$fit$vectors[4])) %>%
        rownames_to_column("EnvVariables")
outdoor_gg_envfit_results2 =outdoor_gg_envfit_results %>% rename(arrows.NMDS1_OUT=arrows.NMDS1, arrows.NMDS2_OUT=arrows.NMDS2, r_OUT=r, pvals_OUT = pvals)

#Merge all R2 and p-values in a table (used for the supplemental Table S3)
gg_envfit_results = all_gg_envfit_results2 %>%  left_join(outdoor_gg_envfit_results2) %>% left_join(indoor_gg_envfit_results2)


##################################################
# 3. PERMANOVA analysis using ADONIS
#########################
#We used the function adonis2 from the vegan package
#On the rarefied Hellinger-transformed tables 
#Similar analyses were done for the 3 datasets:
        #ALL dust samples
        #INDOOR dust samples
        #OUTDOOR dust samples


#########################
# 3.1. PERMANOVA analysis using ADONIS  - ALL SAMPLES 
#########################

#Select the metadata 
all_metadata_for_adonis=all_metadata %>% 
        select(-house,-seq_depth)
#Remark: Note that "latitude", "longitude" and "construction_year" were already excluded to protect the personal data of volunteers
# Therefore, this example shows slightly different results compared to those detailed in the article considering all these variables

#Preparation of the table with results
mylist=colnames(all_metadata_for_adonis)[2:30] 
mylist

results_all=data.frame(matrix(nrow=length(mylist),ncol=3))
colnames(results_all)=c("Group","R2","p")

#Adonis test in loop
for (i in (1:length(mylist))) { 
        sub=data.frame(all_metadata_for_adonis[,colnames(all_metadata_for_adonis)==mylist[i]])
        colnames(sub)=c("variable")
        model_data<-adonis2(otutable_hel~variable, permutations = 999, data=sub)
        print(mylist[i])
        print(model_data)
        results_all$Group[i]=mylist[i]
        results_all$p[i]=model_data$`Pr(>F)`[1]
        results_all$R2[i] =model_data$R2[1]
}

#Results used for the Table 1
#Correction of p-values using the Bonferroni method
results_Adonis_dust_all_samp<-results_all %>%
        mutate (p_bonferroni= (p.adjust (results_all$p, method="bonferroni"))) %>% 
        arrange(desc(R2))%>%
        mutate(dataset="all_samples")

#Plotting results as a bar chart
results_Adonis_dust_all_samp %>%
        ggplot(aes(x=reorder(Group,R2),y=R2,fill=p_bonferroni))+
        geom_bar(stat="identity")+
        coord_flip()+
        scale_fill_gradient(low = "red", high = "black")+
        xlab("") +
        labs(title= "Adonis test", subtitle = "ALL SAMPLES")


#########################
# 3.2. PERMANOVA analysis using ADONIS  - INDOOR SAMPLES 
#########################

#Select the metadata 
indoor_metadata_for_adonis=indoor_metadata %>% 
        select(-indoor_outdoor, -house, -seq_depth)
#Remark: Note that latitude and longitude were already excluded to protect the personal data of volunteers

#Preparation of the table with results
mylist=colnames(indoor_metadata_for_adonis)[2:30] 
mylist

results_indoor=data.frame(matrix(nrow=length(mylist),ncol=3))
colnames(results_indoor)=c("Group","R2","p")

#Adonis test in loop
for (i in (1:length(mylist))) { 
        sub=data.frame(indoor_metadata_for_adonis[,colnames(indoor_metadata_for_adonis)==mylist[i]])
        colnames(sub)=c("variable")
        model_data<-adonis2(indoor_otutable_hel~variable, permutations = 999, data=sub)
        print(mylist[i])
        print(model_data)
        results_indoor$Group[i]=mylist[i]
        results_indoor$p[i]=model_data$`Pr(>F)`[1]
        results_indoor$R2[i] =model_data$R2[1]
}


#Results used for the Table 1
#Correction of p-values using the Bonferroni method
results_Adonis_dust_indoor<-results_indoor %>%
        mutate (p_bonferroni= (p.adjust (results_indoor$p, method="bonferroni"))) %>% 
        arrange(desc(R2))%>%
        mutate(dataset="indoor")

#Plotting results as a bar chart
results_Adonis_dust_indoor %>% 
        ggplot(aes(x=reorder(Group,R2),y=R2,fill=p_bonferroni))+
        geom_bar(stat="identity")+
        coord_flip()+
        scale_fill_gradient(low = "red", high = "black")+
        xlab("")+
        labs(title= "Adonis test", subtitle= "INDOOR")


#########################
# 3.3. PERMANOVA analysis using ADONIS  - OUTDOOR SAMPLES 
#########################

#Select the metadata 
outdoor_metadata_for_adonis=outdoor_metadata %>%
        select(-indoor_outdoor, -location, -house,-seq_depth)
#Remark: Note that latitude and longitude were already excluded to protect the personal data of volunteers

#Preparation of the table with results
mylist=colnames(outdoor_metadata_for_adonis)[2:29] 
mylist

results_outdoor=data.frame(matrix(nrow=length(mylist),ncol=3))
colnames(results_outdoor)=c("Group","R2","p")

#Adonis test in loop
for (i in (1:length(mylist))) { 
        sub=data.frame(outdoor_metadata_for_adonis[,colnames(outdoor_metadata_for_adonis)==mylist[i]])
        colnames(sub)=c("variable")
        model_data<-adonis2(outdoor_otutable_hel~variable, permutations = 999, data=sub)
        print(mylist[i])
        print(model_data)
        results_outdoor$Group[i]=mylist[i]
        results_outdoor$p[i]=model_data$`Pr(>F)`[1]
        results_outdoor$R2[i] =model_data$R2[1]
}


#Results used for the Table 1
#Correction of p-values using the Bonferroni method
results_Adonis_dust_outdoor<-results_outdoor %>%
        mutate (p_bonferroni= (p.adjust (results_outdoor$p, method="bonferroni"))) %>% 
        arrange(desc(R2))%>%
        mutate(dataset="outdoor")

#Plotting results as a bar chart
results_Adonis_dust_outdoor %>%
        ggplot(aes(x=reorder(Group,R2),y=R2,fill=p_bonferroni))+
        geom_bar(stat="identity")+
        coord_flip()+
        scale_fill_gradient(low = "red", high = "black")+
        xlab("")+
        labs(title= "Adonis test", subtitle= "OUTDOOR_hel")


##################################################
# 4. Variation Partitioning Analysis (VPA)
#########################
#On the rarefied Hellinger-transformed tables 
#Similar analyses were done for the 3 datasets:
        #ALL dust samples
        #INDOOR dust samples
        #OUTDOOR dust samples


#########################
# 4.1. VPA  - ALL SAMPLES 
#########################

#Groups of variables for the analysis
building_vars=all_metadata %>% 
        select( building_type, building_material, ventilation_type, water_damage, moisture_problem, odor_problem, pest_type, construction_year, dust_coverage)

occupant_vars=all_metadata %>% 
        select(people, females, children, asthma, allergy_type, pet_type)

climate_vars=all_metadata %>% 
        select(BIO1, BIO4, BIO9, BIO10, BIO11, BIO12, growing_season_length, swe_4, sca_2)

#Adding the variable house compartment (location) for comparison 
house_compartment=all_metadata %>% 
        select(location)

#VPA
varpart_all_dataset.bc <- varpart(vegdist(otutable_hel), building_vars, occupant_vars, climate_vars, house_compartment)

#Check results
varpart_all_dataset.bc

# VPA 
plot(varpart_all_dataset.bc)
plot(varpart_all_dataset.bc, digits=4, bg=2:5, cex=0.7)
plot(varpart_all_dataset.bc, digits=4, bg=2:5, title(main = "ALL SAMPLES_vegdist; X1=Building, X2=Occupant, X3=Climate, X4=House Compartment"))
####Remark: this plot used for the Figure 4
 

#########################
# 4.1. VPA  - INDOOR SAMPLES 
#########################

#Groups of variables for the analysis
building_vars_indoor=indoor_metadata %>% 
        select( building_type, building_material, ventilation_type, water_damage, moisture_problem, odor_problem, pest_type, construction_year, dust_coverage)

occupant_vars_indoor=indoor_metadata %>% 
        select(people, females, children, asthma, allergy_type, pet_type)

climate_vars_indoor=indoor_metadata %>% 
        select(BIO1, BIO4, BIO9, BIO10, BIO11, BIO12, growing_season_length, swe_4, sca_2)

#Adding the variable house compartment (location) for comparison 
house_compartment_indoor=indoor_metadata %>% 
        select(location)

#VPA
varpart_indoor_dataset.bc <- varpart(vegdist(indoor_otutable_hel),building_vars_indoor, occupant_vars_indoor, climate_vars_indoor, house_compartment_indoor)
varpart_indoor_dataset.bc

#VPA plot
plot(varpart_indoor_dataset.bc)
plot(varpart_indoor_dataset.bc, digits=4, bg=2:5, cex=0.7)
plot(varpart_indoor_dataset.bc, digits=4, bg=2:5, title(main = "INDOOR_vegdist; X1=Building, X2=Occupant, X3=Climate, X4=Room"))
####Remark: this plot used for the Figure 4


#########################
# 4.1. VPA  - OUTDOOR SAMPLES 
#########################

#Groups of variables for the analysis
building_vars_outdoor=outdoor_metadata %>% 
        select( building_type, building_material, ventilation_type, water_damage, moisture_problem, odor_problem, pest_type, construction_year, dust_coverage)

occupant_vars_outdoor=outdoor_metadata %>% 
        select(people, females, children, asthma, allergy_type, pet_type)

climate_vars_outdoor=outdoor_metadata %>% 
        select(BIO1, BIO4, BIO9, BIO10, BIO11, BIO12, growing_season_length, swe_4, sca_2)

#VPA of only 3 groups, excluding the variable house compartment (no sense for the outdoor dataset)
varpart_outdoor_dataset.bc <- varpart(vegdist(outdoor_otutable_hel),building_vars_outdoor, occupant_vars_outdoor, climate_vars_outdoor)
varpart_outdoor_dataset.bc

#Plot
plot(varpart_outdoor_dataset.bc)
plot(varpart_outdoor_dataset.bc, digits=4, bg=2:5, cex=0.7)
plot(varpart_outdoor_dataset.bc, digits=4, bg=2:5, title(main = "OUTDOOR_vegdist; X1=Building, X2=Occupant, X3=Climate"))
####Remark: this plot used for the Figure 4

#An additional VPA adding a new group of variables related to the land cover 
landcover_vars_outdoor=outdoor_metadata %>% 
        select(ar50,geonorge_bedrock_nutrient, solar_radiation, area)

#VPA
varpart_outdoor_dataset2.bc <- varpart(vegdist(outdoor_otutable_hel),building_vars_outdoor, occupant_vars_outdoor, climate_vars_outdoor, landcover_vars_outdoor)
varpart_outdoor_dataset2.bc

#Plot
plot(varpart_outdoor_dataset2.bc)
plot(varpart_outdoor_dataset2.bc, digits=4, bg=2:5, cex=0.7)
plot(varpart_outdoor_dataset2.bc, digits=4, bg=2:5, title(main = "OUTDOOR_vegdist; X1=Building, X2=Occupant, X3=Climate, X4=Land cover"))


#################### END OF THE SCRIPT ##########################
