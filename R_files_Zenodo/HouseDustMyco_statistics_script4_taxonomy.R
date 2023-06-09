
#########################Filename: HouseDustMyco_statistics_script4_taxonomy.R

#Publication: Analyzing indoor mycobiomes through a large-scale citizen science study in Norway
#Authors: Pedro M. Martin-Sanchez, Eva Lena F. Estensmo, Luis N. Morgado, Sundy Maurice, Ingeborg B. Engh, Inger Skrede and Håvard Kauserud

##Related to the sections: 
        #Materials and Methods- 2.5. Statistical analyses
        #Results- 3.3. Dominant fungi in house dust
        #Supplemental Information

## Goals: 
#(1) To evaluate the taxonomic distribution of the dominant fungi in  house-dust samples - Relative abundances of taxa
#(2) To evaluate the FUNGuild distribution of the dominant fungi - Relative abundances of trophic modes and guilds

#Inputs: 
        #The rarefied fungal OTU table ("OTU_table_rar2000x10.xlsx")
        #The metadata file for dust samples ("metadata_file_dust_samples.xlsx")
        #The taxonomy file for the OTUs in the rarefied table ("taxonomy_for_rarTable_funguild.xlsx")

#Outputs: 
        #Plots for the Figure 5 (Taxonomic distribution)
        #Plots for the Figure S9 (FUNGuild distribution) 


#########################
#PREPARATIONS
#########################

#Load packages
library(openxlsx)
library(tidyverse)
library(RColorBrewer)


#Set your path
setwd("C:/WP1_data_analyses/script_preparation/")

#Load the rarefied OTU table (6,632 OTUs from 807 dust samples)
OTU_tab_rar=read.xlsx("OTU_table_rar2000x10.xlsx")

#Different format after moving the column "sample" to the row names
#OTU_tab_rar_for_analyses=OTU_tab_rar %>% remove_rownames() %>% column_to_rownames("sample")

#Load the metadata
metadata=read.xlsx("metadata_file_dust_samples.xlsx")
#Remark: Some geographic variables(latitude, longitude, municipality, county and region) as well as the construction year of the study houses
# were excluded from this file to protect personal data of volunteers   

#Load the taxonomy file for the 6,632 OTUs in the rarefied table, including the FUNGuild annotations
taxonomy=read.xlsx("taxonomy_for_rarTable_funguild.xlsx")


### Merge OTU table with taxonomy
OTU_tab_rar_with_tax<-OTU_tab_rar %>% 
        remove_rownames() %>% 
        column_to_rownames("sample") %>% 
        t() %>% 
        as.data.frame() %>% 
        rownames_to_column("OTU_id") %>% 
        left_join(taxonomy)

### Preparation of data:
#Select and order columns (excluding FUNGuild annotation), calculate read abundances per sample, and merge metadata
OTU_tab_rar_with_tax_and_metadata<-OTU_tab_rar_with_tax %>% 
        select(OTU_id,c(809:824) ,c(2:808)) %>% 
        gather(sample, read_abundance, -c(1:17)) %>% 
        left_join(metadata)


##################################################
# 1. TAXONOMIC DISTRIBUTION
#########################
#To show the most abundant taxa by house compartments (location)
#Relative abundances of taxa are mean values per sample calculated based on the rarefied matrix
#For four taxonomic levels: 
        #Phyllum
        #Order
        #Genus
        #Species

#Select the metadata of interest:sample and house compartment (location)
OTU_tab_rar_with_tax_and_metadata_location=OTU_tab_rar_with_tax_and_metadata %>%  
        select(sample, location) %>% 
        distinct(sample,.keep_all = T)

#########################
# 1.1. TAXONOMIC DISTRIBUTION - PHYLLUM
#########################

#Calculation of data: mean relative abundances and standard deviations
tax_phyllum_per_sample<-OTU_tab_rar_with_tax_and_metadata %>% 
        group_by(sample, phyllum) %>% 
        summarise(sum_phyllum=sum(read_abundance)) %>% 
        left_join(OTU_tab_rar_with_tax_and_metadata_location) %>% 
        group_by(sample) %>% 
        mutate(Seq_per_sample=sum(sum_phyllum)) %>% 
        mutate(RSA_per_sample_per_phyllum=sum_phyllum/Seq_per_sample) %>% 
        group_by(location, phyllum) %>% 
        summarise(sd_rsa=sd(RSA_per_sample_per_phyllum),avg_relative_abundance=mean(RSA_per_sample_per_phyllum)) %>% 
        filter(!is.na(phyllum))%>%
        #top_n(3) %>% #To filter the top 3 phyla
        mutate(new_std_min=avg_relative_abundance-sd_rsa) %>% # get the min of sd, and correction when min<0 (min2)
        mutate(new_std_min_2=ifelse(new_std_min<0,0,new_std_min))
#Remark: Not assigned (NA) phyla were excluded: 1.9% outside, 0.9% living room (central), 1.1% bathroom

#Plot
tax_phyllum_plot=tax_phyllum_per_sample %>% 
        ggplot(aes(x=reorder(phyllum,avg_relative_abundance), y=avg_relative_abundance, fill=location))+
        geom_bar(stat = "identity", position=position_dodge(0.9))+
        geom_errorbar(aes(ymin = new_std_min_2, ymax = avg_relative_abundance+sd_rsa), 
                      width = 0.1,position =position_dodge(0.9))+
        scale_fill_manual(values= c("orange", "blue", "darkgrey"))+
        ylab("Relative abundance")+
        xlab("")+
        labs(fill="House compartment")+
        labs(title= "Phylla")+
        theme_bw()+
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
              legend.position = c(0.8,0.2))+
        guides(fill = guide_legend(reverse=TRUE))+
        coord_flip()
tax_phyllum_plot
#Remark: This plot was used for Figure 5a


#########################
# 1.2. TAXONOMIC DISTRIBUTION - ORDER
#########################

#Calculation of data for the top 15 orders: mean relative abundances and standard deviations
tax_order_per_sample<-OTU_tab_rar_with_tax_and_metadata %>% 
        group_by(sample, order) %>% 
        summarise(sum_order=sum(read_abundance)) %>% 
        left_join(OTU_tab_rar_with_tax_and_metadata_location) %>% 
        group_by(sample) %>% 
        mutate(Seq_per_sample=sum(sum_order)) %>% 
        mutate(RSA_per_sample_per_order=sum_order/Seq_per_sample) %>% 
        group_by(location, order) %>% 
        summarise(sd_rsa=sd(RSA_per_sample_per_order),avg_relative_abundance=mean(RSA_per_sample_per_order)) %>% 
        filter(!is.na(order)) %>%
        filter(order=="Eurotiales" |order=="Capnodiales" |order=="Pucciniales" |order=="Saccharomycetales"|
                       order=="Agaricales" |order=="Pleosporales" |order=="Lecanorales" | order=="Helotiales"| 
                       order=="Dothideales" |order=="Polyporales" |order=="Malasseziales" |order=="Chaetothyriales" |
                       order=="Mucorales" |order=="Tremellales" |order=="Filobasidiales")%>%
        mutate(new_std_min=avg_relative_abundance-sd_rsa) %>% 
        mutate(new_std_min_2=ifelse(new_std_min<0,0,new_std_min))
#Remark: Not assigned (NA) orders were excluded: 5.5% outside, 1.9% living room (central), 2.5% bathroom

#Plot
tax_order_plot=tax_order_per_sample %>%
        ggplot(aes(x=reorder(order,avg_relative_abundance), y=avg_relative_abundance, fill=location))+
        geom_bar(stat = "identity", position=position_dodge(0.9))+
        geom_errorbar(aes(ymin = new_std_min_2, ymax = avg_relative_abundance+sd_rsa), 
                      width = 0.1,position =position_dodge(0.9))+
        scale_fill_manual(values= c("orange", "blue", "darkgrey"))+
        ylab("Relative abundance")+
        xlab("")+
        labs(title= "Top Orders")+
        theme_bw()+
        labs(fill="House compartment")+
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
              legend.position = c(0.8,0.2))+
        guides(fill = guide_legend(reverse=TRUE))+
        coord_flip()
tax_order_plot
#Remark: This plot was used for Figure 5b

#########################
# 1.3. TAXONOMIC DISTRIBUTION - GENUS
#########################

#Calculation of data for the top 20 orders: mean relative abundances and standard deviations
tax_genera_per_sample<-OTU_tab_rar_with_tax_and_metadata %>% 
        group_by(sample, genus) %>% 
        summarise(sum_genus=sum(read_abundance)) %>% 
        left_join(OTU_tab_rar_with_tax_and_metadata_location) %>% 
        group_by(sample) %>% 
        mutate(Seq_per_sample=sum(sum_genus)) %>% 
        mutate(RSA_per_sample_per_genus=sum_genus/Seq_per_sample) %>% 
        group_by(location, genus) %>% 
        summarise(sd_rsa=sd(RSA_per_sample_per_genus),avg_relative_abundance=mean(RSA_per_sample_per_genus)) %>% 
        filter(!is.na(genus)) %>% 
        filter(genus=="Cladosporium" | genus=="Penicillium" | genus=="Aspergillus"| #Select the top 20 to avoid the bias of using "top_()"
                       genus=="Saccharomyces"| genus=="Thekopsora"| genus=="Verrucocladosporium"| genus=="Malassezia"|
                       genus=="Botrytis"| genus=="Aureobasidium"| genus=="Epicoccum"|genus=="Scoliciosporum"|genus=="Lycoperdon"|
                       genus=="Hypogymnia"|genus=="Strobilurus"|genus=="Fomitopsis"|genus=="Debaryomyces"|genus=="Cylindrobasidium"|
                       genus=="Melampsora"|genus=="Fomes"|genus=="Naganishia")%>%
        mutate(new_std_min=avg_relative_abundance-sd_rsa) %>% 
        mutate(new_std_min_2=ifelse(new_std_min<0,0,new_std_min))
#Remark: Not assigned (NA) genera were excluded: 16.3% outside, 7.3% living room (central), 8% bathroom

#Plot
tax_genera_plot=tax_genera_per_sample %>%
        ggplot(aes(x=reorder(genus,avg_relative_abundance), y=avg_relative_abundance, fill=location))+
        geom_bar(stat = "identity", position=position_dodge(0.9))+
        geom_errorbar(aes(ymin = new_std_min_2, ymax = avg_relative_abundance+sd_rsa), 
                      width = 0.1,position =position_dodge(0.9))+
        scale_fill_manual(values= c("orange", "blue", "darkgrey"))+
        ylab("Relative abundance")+
        xlab("")+
        labs(title= "Top 20 Genera")+
        theme_bw()+
        labs(fill="House compartment")+
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
              legend.position = c(0.8,0.2))+
        guides(fill = guide_legend(reverse=TRUE))+
        coord_flip()
tax_genera_plot
#Remark: This plot was used for Figure 5c


#########################
# 1.4. TAXONOMIC DISTRIBUTION - SPECIES
#########################

#Calculation of data for species with >1% reads: mean relative abundances and standard deviations
tax_species_per_sample<-OTU_tab_rar_with_tax_and_metadata %>% 
        group_by(sample, species) %>% 
        summarise(sum_species=sum(read_abundance)) %>% 
        left_join(OTU_tab_rar_with_tax_and_metadata_location) %>% 
        group_by(sample) %>% 
        mutate(Seq_per_sample=sum(sum_species)) %>% 
        mutate(RSA_per_sample_per_species=sum_species/Seq_per_sample) %>% 
        group_by(location, species) %>% 
        summarise(sd_rsa=sd(RSA_per_sample_per_species),avg_relative_abundance=mean(RSA_per_sample_per_species)) %>% 
        filter(!is.na(species)) %>% 
        filter(avg_relative_abundance>0.01)%>% #21 species
        mutate(new_std_min=avg_relative_abundance-sd_rsa) %>% 
        mutate(new_std_min_2=ifelse(new_std_min<0,0,new_std_min))
#Remark: Not assigned (NA) species were excluded: 26.8% outside, 25.6% living room (central), 26.3% bathroom

#Plot
tax_species_plot=tax_species_per_sample %>%
        ggplot(aes(x=reorder(species,avg_relative_abundance), y=avg_relative_abundance, fill=location))+
        geom_bar(stat = "identity", position=position_dodge(0.9))+
        geom_errorbar(aes(ymin = new_std_min_2, ymax = avg_relative_abundance+sd_rsa), 
                      width = 0.1,position =position_dodge(0.9))+
        scale_fill_manual(values= c("orange", "blue", "darkgrey"))+
        ylab("Relative abundance")+
        xlab("")+
        labs(title= "Species RA>1%")+
        theme_bw()+
        labs(fill="House compartment")+
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
              legend.position = c(0.8,0.2))+
        guides(fill = guide_legend(reverse=TRUE))+
        coord_flip()
tax_species_plot


##################################################
# 2. FUNGUILD DISTRIBUTION 
#########################
#To show the most abundant trophic modes and guilds by house compartments (location)
#Relative abundances are  calculated on the rarefied matrix


#########################
# 2.1. FUNGUILD DISTRIBUTION - TROPHIC MODES
#########################

#Summary of relative abundances (% of rarefied reads) of the 7 trophic modes. 9.4% reads were Not Assigned (NA; "-")
tax_RA_trophic<-OTU_tab_rar_with_tax_and_metadata %>% 
        group_by(Trophic.Mode) %>% 
        summarise(sum_tropic=sum(read_abundance)) %>% 
        #filter(Trophic.Mode!="-")  %>% # to exclude the NA OTUs 
        mutate(Relative_abundance=sum_tropic/sum(sum_tropic))

##Comparing indoor vs. outdoor
trophic_indoor_outdoor<-OTU_tab_rar_with_tax_and_metadata %>% 
        group_by(indoor_outdoor, Trophic.Mode) %>% 
        summarise(sum_trophic =sum(read_abundance)) %>% 
        #filter(Trophic.Mode!="-") %>% 
        group_by(indoor_outdoor) %>% 
        mutate(Relative_abundance=sum_trophic/sum(sum_trophic))

#Number of levels in Trophic.Mode
colourCount1 <- length(unique(trophic_indoor_outdoor$Trophic.Mode)) 
colourCount1 #8

#Plot
trophic_indoor_outdoor %>% 
        ggplot(aes(x=indoor_outdoor, y=Relative_abundance, fill=Trophic.Mode))+
        geom_bar(stat = "identity")+
        #coord_flip()+
        scale_fill_manual(values = colorRampPalette(brewer.pal(12, 
                                                               "Spectral"))(colourCount1))+
        ylab("Relative abundance")+
        xlab("Trophic modes")+
        labs(fill="")+
        labs(title= "Trophic modes", subtitle = "IN vs. OUT")+
        theme_bw()

##Comparing house compartments (location)
trophic_location<-OTU_tab_rar_with_tax_and_metadata %>% 
        group_by(location, Trophic.Mode) %>% 
        summarise(sum_trophic =sum(read_abundance)) %>% 
        #filter(Trophic.Mode!="-") %>% 
        group_by(location) %>% 
        mutate(Relative_abundance=sum_trophic/sum(sum_trophic))

#Number of levels in Trophic.Mode
colourCount2 <- length(unique(trophic_location$Trophic.Mode)) # number of levels
colourCount2 #8

#Plot
trophic_location %>% 
        ggplot(aes(x=location, y=Relative_abundance, fill=Trophic.Mode))+
        geom_bar(stat = "identity")+
        scale_fill_manual(values = colorRampPalette(brewer.pal(12, 
                                                               "Spectral"))(colourCount2))+
        ylab("Relative abundance")+
        xlab("")+
        labs(fill="Trophic modes")+
        labs(title= "Trophic modes", subtitle = "House compartments")+
        theme_bw()
#Remark: This plot was used for Figure S9a


#########################
# 2.2. FUNGUILD DISTRIBUTION - GUILDS
#########################

#Summary of relative abundances (% of rarefied reads) of the 131 Guilds. 9.4% reads correspond to NA
tax_RA_guild<-OTU_tab_rar_with_tax_and_metadata %>% 
        group_by(Guild) %>% 
        summarise(sum_guild=sum(read_abundance)) %>% 
        #filter(Guild!="-")  %>% # to exclude the NA OTUs 
        mutate(Relative_abundance=sum_guild/sum(sum_guild))

##Comparing indoor vs. outdoor
guild_indoor_outdoor<-OTU_tab_rar_with_tax_and_metadata %>% 
        group_by(indoor_outdoor, Guild) %>% 
        summarise(sum_guild =sum(read_abundance)) %>%
        #filter(Guild!="-") %>% 
        group_by(indoor_outdoor) %>% 
        mutate(Relative_abundance=sum_guild/sum(sum_guild))%>%
        filter(Relative_abundance>=0.01)

#Number of levels in Guild
colourCount3 <- length(unique(guild_indoor_outdoor$Guild)) # number of levels
colourCount3 #13

#Plot
guild_indoor_outdoor %>% 
        ggplot(aes(x=indoor_outdoor, y=Relative_abundance, fill=Guild))+
        geom_bar(stat = "identity")+
        scale_fill_manual(values = colorRampPalette(brewer.pal(12, 
                                                               "Spectral"))(colourCount3))+
        ylab("Relative abundance")+
        xlab("")+
        labs(fill="Guilds")+
        labs(title= "Guilds; RA>=1%", subtitle = "IN vs. OUT")+
        theme_bw()

##Comparing house compartments (location)
guild_location<-OTU_tab_rar_with_tax_and_metadata %>% 
        group_by(location, Guild) %>% 
        summarise(sum_guild =sum(read_abundance)) %>% 
        #filter(Guild!="-") %>% 
        group_by(location) %>% 
        mutate(Relative_abundance=sum_guild/sum(sum_guild))%>%
        filter(Relative_abundance>=0.01)

#Number of levels in Guild
colourCount4 <- length(unique(guild_location$Guild))
colourCount4 #13

#Plot
guild_location %>% 
        ggplot(aes(x=location, y=Relative_abundance, fill=Guild))+
        geom_bar(stat = "identity")+
        #coord_flip()+
        scale_fill_manual(values = colorRampPalette(brewer.pal(12, 
                                                               "Spectral"))(colourCount4))+
        ylab("Relative abundance")+
        xlab("")+
        ylim(0, 1)+
        labs(fill="Guilds")+
        labs(title= "Guilds; RA>=1%", subtitle = "House compartments")+
        theme_bw()
#Remark: This plot was used for Figure S9b


#################### END OF THE SCRIPT ##########################
