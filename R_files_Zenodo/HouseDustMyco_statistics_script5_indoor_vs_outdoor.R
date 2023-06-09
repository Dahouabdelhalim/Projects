
#########################Filename: HouseDustMyco_statistics_script5_indoor_vs_outdoor.R

#Publication: Analyzing indoor mycobiomes through a large-scale citizen science study in Norway
#Authors: Pedro M. Martin-Sanchez, Eva Lena F. Estensmo, Luis N. Morgado, Sundy Maurice, Ingeborg B. Engh, Inger Skrede and Håvard Kauserud

##Related to the sections: 
        #Materials and Methods- 2.5. Statistical analyses
        #Results- 3.3. Indoor versus outdoor mycobiomes
        #Supplemental Information

## Goals: 
#(1) To evaluate the distribution of OTUs by house compartments - Venn diagrams
#(2) To identify indicator species of indoor and outdoor environments

#Inputs: 
        #The rarefied fungal OTU table ("OTU_table_rar2000x10.xlsx")
        #The metadata file for dust samples ("metadata_file_dust_samples.xlsx")
        #The taxonomy file for the OTUs in the rarefied table ("taxonomy_for_rarTable_funguild.xlsx")

#Outputs: 
        #Plots for the Figure 6 (Venn diagrams showing the distribution of OTUs by house compartments) 
        #Data for the Table 2 (Indicator species)         
        #Plots for the Figure S10 (Taxonomic affiliation at order level of indicator species) 


#########################
#PREPARATIONS
#########################

#Load packages
library(openxlsx)
library(tidyverse)
library(vegan)
library(VennDiagram)
library(indicspecies)
library(RColorBrewer)


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


##################################################
# 1. Venn diagrams showing the distribution of dust mycobiomes across the three house compartments
#########################
#Calculate number and proportion of OTUs by three ways:
        #(1.1) across overall data 
        #(1.2) after removing low-abundance OTUs (< 10 reads per sample)
        #(1.3) comparing at a house-by-house basis 
                #(1.4) GLM model


#########################
# 1.1. Venn diagram by house compartment - Overall data
#########################

#Summarise the number of houses where the OTUs are present for each house compartment
data_for_shared_OTUs_analyses<-master_table %>% 
        mutate(PA=sign(rarRead)) %>% 
        group_by(OTU_id, location) %>% 
        summarise(number_of_houses=sum(PA)) %>% 
        filter(number_of_houses>0)

#Table with lists of OTUs for each house compartment (including NAs)
lists_OTUS_location<-data_for_shared_OTUs_analyses %>% 
        rownames_to_column("keys") %>% 
        spread(location,OTU_id) %>% 
        select(-keys, -number_of_houses)

#Number of OTUs per house compartments, after removing NAs
OTUS_central <-na.omit(lists_OTUS_location$central)
OTUS_bathroom<-na.omit(lists_OTUS_location$bathroom)
OTUS_outside <-na.omit(lists_OTUS_location$outside)

length(OTUS_central) #4224
length(OTUS_bathroom)#4643
length(OTUS_outside) #4170

#Numbers of OTUs for intersections
length(intersect(OTUS_central,OTUS_bathroom))#3150
length(intersect(OTUS_bathroom,OTUS_outside))#2956
length(intersect(OTUS_central,OTUS_outside)) #2807
length(intersect(intersect(OTUS_central,OTUS_bathroom),OTUS_outside))#2408

#Plotting the venn diagram
grid.newpage(); 
venn.plot <- draw.triple.venn(
        area1 = 4324,
        area2 = 4643,
        area3 = 4170,
        n12 = 3150,
        n23 = 2956,
        n13 = 2807,
        n123 = 2408,
        category = c("Living room", "Bathroom", "Outside"),
        fill = c("darkblue", "orange", "grey"),
        lty = "blank",
        cex = 1,
        cat.cex = 4,
        cat.col = c("darkblue", "orange", "grey"),
        print.mode = "percent", #Change to "raw" to plot the number of OTUs instead percentages 
);
grid.draw(venn.plot);

#Remarks: This plot corresponds to Figure 6a-left
#Results show a relatively low overlap between living room-bathroom or indoor-outdoor;
#perhaps this is related to the rare fungi (OTUS with low number of reads)


#########################
# 1.2. Venn diagram by house compartment - after removing low-abundance OTUs
#########################

#Preparation of data filtering OTUs with >=10 reads per sample
data_for_shared_OTUs_more_10reads<-master_table %>% filter(rarRead>=10)%>%
        mutate(PA=sign(rarRead)) %>% 
        group_by(OTU_id, location) %>% 
        summarise(number_of_houses=sum(PA)) %>% 
        filter(number_of_houses>0)

#Table with lists of OTUs for each house compartment (including NAs)
lists_OTUS_location_more_10reads<-data_for_shared_OTUs_more_10reads %>% 
        rownames_to_column("keys") %>% 
        spread(location,OTU_id) %>% 
        select(-keys, -number_of_houses)

#Number of OTUs per house compartments, after removing NAs
OTUS_central_more_10reads <-na.omit(lists_OTUS_location_more_10reads$central)
OTUS_bathroom_more_10reads<-na.omit(lists_OTUS_location_more_10reads$bathroom)
OTUS_outside_more_10reads <-na.omit(lists_OTUS_location_more_10reads$outside)

length(OTUS_central_more_10reads) #830
length(OTUS_bathroom_more_10reads)#1008
length(OTUS_outside_more_10reads) #1104

#Numbers of OTUs for intersections
length(intersect(OTUS_central_more_10reads,OTUS_bathroom_more_10reads))#505
length(intersect(OTUS_bathroom_more_10reads,OTUS_outside_more_10reads))#450
length(intersect(OTUS_central_more_10reads,OTUS_outside_more_10reads)) #400
length(intersect(intersect(OTUS_central_more_10reads,OTUS_bathroom_more_10reads),OTUS_outside_more_10reads))#326

#Plotting the venn diagram
grid.newpage(); 
venn.plot_more_10reads <- draw.triple.venn(
        area1 = 830,
        area2 = 1008,
        area3 = 1104,
        n12 = 505,
        n23 = 450,
        n13 = 400,
        n123 = 326,
        category = c("Living room", "Bathroom", "Outside"),
        fill = c("darkblue", "orange", "grey"),
        lty = "blank",
        cex = 1,
        cat.cex = 4,
        cat.col = c("darkblue", "orange", "grey"),
        print.mode = "percent", #Change to "raw" to plot the number of OTUs instead percentages 
);
grid.draw(venn.plot_more_10reads);

#Remarks: This plot corresponds to Figure 6a-right
#Basically, removal/reduction of rare fungi: (i) decreases the overlaps (LR-BR or IN-OUT),
#and (ii) increases the proportion of OTUS that are exclusively present outside 


#########################
# 1.3. Venn diagram by house compartment - house-by-house basis 
#########################
##Calculation of percents of shared OTUS for each sector of the Venn diagram

#Preparation of data
data_for_shared_OTUs_analyses_per_house<-master_table %>% 
        mutate(PA=sign(rarRead)) %>% 
        group_by(OTU_id, location, house) %>% 
        summarise(OTUs_per_house=sum(PA)) %>% 
        filter(OTUs_per_house>0)

#Number of OTUs (richness) per house
richness_per_house = master_table %>%
        mutate(PA=sign(rarRead)) %>% 
        group_by(house, OTU_id) %>% 
        summarise(test=sum(PA)) %>% 
        filter(test>0) %>% 
        mutate(PA2=1) %>% 
        group_by(house) %>% 
        summarise(rich_per_house=sum(PA2))

##Separate calculation for each sector:

#Sector living room - bathroom - outside
data_for_shared_OTUs_bath_cent_out=data_for_shared_OTUs_analyses_per_house %>% 
        group_by(house, OTU_id) %>% 
        summarise(OTUs_per_location_per_house=sum(OTUs_per_house)) %>% 
        filter(OTUs_per_location_per_house==3)

Per_house_bath_cent_out =data_for_shared_OTUs_bath_cent_out %>% 
        group_by(house) %>% 
        tally() %>% 
        left_join(richness_per_house) %>% 
        mutate(out_bath_cent_perc=100*(n/rich_per_house)) %>% 
        rename(out_bath_cent_n=n)

#Sector living room - bathroom
data_for_shared_OTUs_bath_cent=data_for_shared_OTUs_analyses_per_house %>% 
        group_by(house, OTU_id) %>%
        filter(location=="bathroom"| location=="central") %>% 
        summarise(OTUs_per_location_per_house=sum(OTUs_per_house)) %>% 
        filter(OTUs_per_location_per_house==2)

Per_house_bath_cent=data_for_shared_OTUs_bath_cent %>% 
        group_by(house) %>% 
        tally() %>% 
        full_join(Per_house_bath_cent_out) %>% #to subtract the number of OTUS in out_bath_cent
        mutate(bath_cent_n=n-out_bath_cent_n) %>%
        mutate(bath_cent_perc=100*(bath_cent_n/rich_per_house))

Per_house_bath_cent= Per_house_bath_cent[c(1, 6,7)]

#Sector bathroom - outside
data_for_shared_OTUs_bath_out=data_for_shared_OTUs_analyses_per_house %>% 
        group_by(house, OTU_id) %>% 
        filter(location=="bathroom"| location=="outside") %>% 
        summarise(OTUs_per_location_per_house=sum(OTUs_per_house)) %>% 
        filter(OTUs_per_location_per_house==2)

Per_house_bath_out=data_for_shared_OTUs_bath_out %>% 
        group_by(house) %>% 
        tally() %>% 
        full_join(Per_house_bath_cent_out) %>% # to subtract the number of OTUS in out_bath_cent
        mutate(bath_out_n=n-out_bath_cent_n) %>%
        mutate(bath_out_perc=100*(bath_out_n/rich_per_house))

Per_house_bath_out= Per_house_bath_out[c(1, 6,7)]

#Sector living room - outside
data_for_shared_OTUs_cent_out=data_for_shared_OTUs_analyses_per_house %>% 
        group_by(house, OTU_id) %>% 
        filter(location=="central"| location=="outside") %>% 
        summarise(OTUs_per_location_per_house=sum(OTUs_per_house)) %>% 
        filter(OTUs_per_location_per_house==2)

Per_house_cent_out=data_for_shared_OTUs_cent_out %>% 
        group_by(house) %>% 
        tally() %>%
        full_join(Per_house_bath_cent_out) %>% # to subtract the number of OTUS in out_bath_cent
        mutate(cent_out_n=n-out_bath_cent_n) %>%
        mutate(cent_out_perc=100*(cent_out_n/rich_per_house))

Per_house_cent_out= Per_house_cent_out[c(1, 6,7)]

#Sector outside
data_for_shared_OTUs_out=data_for_shared_OTUs_analyses_per_house %>% 
        group_by(house, OTU_id) %>% 
        filter(location=="outside") %>% 
        summarise(OTUs_per_location_per_house=sum(OTUs_per_house))

Per_house_out=data_for_shared_OTUs_out %>% 
        group_by(house) %>% 
        tally() %>% 
        full_join(Per_house_bath_cent_out) %>% # to subtract the number of OTUS in other sectors
        full_join(Per_house_bath_out) %>% 
        full_join(Per_house_cent_out) %>% 
        mutate(out_n = n-(out_bath_cent_n+bath_out_n+cent_out_n)) %>%
        mutate(out_perc=100*(out_n/rich_per_house))

Per_house_out= Per_house_out[c(1, 10,11)]

#Sector bathroom
data_for_shared_OTUs_bath=data_for_shared_OTUs_analyses_per_house %>% 
        group_by(house, OTU_id) %>% 
        filter(location=="bathroom") %>% 
        summarise(OTUs_per_location_per_house=sum(OTUs_per_house))

Per_house_bath=data_for_shared_OTUs_bath %>% 
        group_by(house) %>% 
        tally() %>%
        full_join(Per_house_bath_cent_out) %>% #to subtract the number of OTUS in other sectors
        full_join(Per_house_bath_out) %>% 
        full_join(Per_house_bath_cent) %>% 
        mutate(bath_n = n-(out_bath_cent_n+bath_out_n+bath_cent_n)) %>%
        mutate(bath_perc=100*(bath_n/rich_per_house))

Per_house_bath= Per_house_bath[c(1, 10,11)]

#Sector living room
data_for_shared_OTUs_cent=data_for_shared_OTUs_analyses_per_house %>% 
        group_by(house, OTU_id) %>% 
        filter(location=="central") %>% 
        summarise(OTUs_per_location_per_house=sum(OTUs_per_house))

Per_house_cent=data_for_shared_OTUs_cent %>% 
        group_by(house) %>% 
        tally() %>%
        full_join(Per_house_bath_cent_out) %>% #to subtract the number of OTUS in other sectors
        full_join(Per_house_cent_out) %>% 
        full_join(Per_house_bath_cent) %>% 
        mutate(cent_n = n-(out_bath_cent_n+cent_out_n+bath_cent_n)) %>%
        mutate(cent_perc=100*(cent_n/rich_per_house))

Per_house_cent= Per_house_cent[c(1, 10,11)]

#Compile all data 
table_shared_otus_per_house_location = Per_house_bath_cent_out %>%
        full_join(Per_house_bath_cent) %>%
        full_join(Per_house_bath_out) %>%
        full_join(Per_house_cent_out) %>%
        full_join(Per_house_bath) %>%
        full_join(Per_house_cent)%>%
        full_join(Per_house_out)

##There are some houses where one house compartment (location) is missing:
#"outside" in houses 44, 109, 126, 246, 269 
#"central" in the house 88

# Thus, we delete these houses
table_shared_otus_per_house_location = table_shared_otus_per_house_location %>%
        filter(house !="44", house != "88", house !="109", house != "126", house !="246", house != "269")

#calculate averages of percentages for each sector and their standard deviations  (total number of houses =265)
average_percentages_shared_OTUs_per_house=table_shared_otus_per_house_location %>% 
        select(house, contains("_perc")) %>% 
        gather(shared_locations, perc_shared, -house) %>%
        group_by(shared_locations) %>% 
        summarise(avg_shared=mean(perc_shared), sd_shared=sd(perc_shared))

#Remark: They are the values included in the Venn diagram (Figure 6b) 


#########################
# 1.4. GLM model - shared OTUs (indoor-outdoor) 
#########################
##To check the effect of the selected environmental variables on the shared OTUS indoor and outdoor (mean percents per house)

#Calculate the sum of percentages of the shared OTUs between indoor and outdoor (in_out_perc_sum)
#Select the data for the GLM model and merge the metadata
shared_otus_per_house_indoor_outdoor_metadata = table_shared_otus_per_house_location%>%
        mutate(in_out_perc_sum=out_bath_cent_perc+bath_out_perc+cent_out_perc)%>%
        select(house, in_out_perc_sum)%>%
        left_join(metadata)%>%
        distinct(house, .keep_all = T)

#Selection of variables. 
data_for_glm_analyses2<-shared_otus_per_house_indoor_outdoor_metadata %>%
        select(#region, county, municipality,latitude, longitude,
               BIO1, BIO4, BIO9, BIO10, BIO11, BIO12, growing_season_length, swe_4, sca_2,
               ar50,geonorge_bedrock_nutrient, solar_radiation, area,
               building_type, building_material, ventilation_type, water_damage, moisture_problem, odor_problem, pest_type,
               people, females, children, asthma, allergy_type, pet_type)

#Remark: NOTE THAT THE VARIABLES "LATITUDE", "LONGITUDE", "MUNICIPALITY", COUNTY" AND "REGION" HAVE BEEN EXCLUDED OF THE METADATA FILE
# In order to protect the personal data (addresses) of the volunteers
# Therefore, this example shows slightly different results compared to those detailed in the article considering al these variables 
  
# Construct the null model used as reference in the first screening (below) 
shared_OTU_per_house<-shared_otus_per_house_indoor_outdoor_metadata %>% select(in_out_perc_sum) 

shared_OTU_null_model=glm(shared_OTU_per_house$in_out_perc_sum~1, data =data_for_glm_analyses2) # The null model is used as reference for the first screening; now there is no offset 
AIC(shared_OTU_null_model) #### AIC  1840.472
BIC(shared_OTU_null_model) #### BIC  1847.632

# First screening of selected variable in a loop
mylist=colnames(data_for_glm_analyses2)[1:26]
results=data.frame(matrix(nrow=length(mylist),ncol=4))
colnames(results)=c("Variable","p", "AIC", "BIC")

for (i in (1:length(mylist))) {
        sub=data.frame(data_for_glm_analyses2[,colnames(data_for_glm_analyses2)==mylist[i]], data_for_glm_analyses2)
        colnames(sub)=c("variable")
        model_data<-glm(shared_OTU_per_house$in_out_perc_sum~variable, data=sub)
        print(mylist[i])
        print(model_data)
        plot(resid(model_data))
        results$Variable[i]=mylist[i]
        results$p[i]=summary(model_data)$coefficients[8]
        results$AIC[i] =AIC(model_data)
        results$BIC[i] =BIC(model_data)
}


Results_shared_OTU_GLM<-results %>% arrange(AIC) # Get and see the results
Results_shared_OTU_GLM
# Variables wih AIC < 1840.472 that was for null model (those vars can improve the explanatory power of the GLM model) 
#                   Variable           p      AIC      BIC
#1              building_type -2.7763779 1836.286 1850.605
#2                       BIO1  0.1170602 1839.992 1850.732
#3                      BIO10  0.1221478 1840.060 1850.799
#4          building_material -1.8620185 1840.097 1854.416

#BUT ALL OF THEM ARE NOT SIGNIFICANT (P>0.05)
#No variable wih BIC<1847.632


##################################################
# 2. Indicator Species analysis
#########################
        #(2.1) Indoor vs. Outdoor; complete dataset 
        #(2.2) Living room vs. Bathroom;  indoor dataset


#########################
# 2.1. Indicator Species analysis - Indoor vs. Outdoor
#########################

#Prepare metadata in the right order
all_metadata=master_table %>% 
        select(-c(2:19)) %>%
        distinct(sample, .keep_all = T)

#Indicator species
test_indoor_outdoor<-multipatt(OTU_tab_rar_for_analyses, all_metadata$indoor_outdoor)

#Summary
Indc_ind_species_indoor_outdoor_summary<-summary.multipatt(test_indoor_outdoor)

#Multilevel pattern analysis
#---------------------------
#        Association function: IndVal.g
#Significance level (alpha): 0.05

#Total number of species: 6632
#Selected number of species: 791 
#Number of species associated to 1 group: 791 

#List of species associated to each combination: 
#       Group indoor  #sps.  241
#       Group outdoor  #sps.  550 


#Prepare results and add complementary information of OTUs (taxonomy, FUNGuild, Relative abundance, distribution as #of samples and #of houses)
Indc_ind_species_indoor_outdoor_summary2<-test_indoor_outdoor$sign %>% #index1=Indoor; Index2=outdoor
        rownames_to_column("OTU") %>% rename(OTU_id="OTU")%>%
        filter(p.value<=0.05) %>%
        mutate(Ind_sp_location=ifelse(index==1,"indoor",
                                      ifelse(index==2,"outdoor",
                                             "other")))%>%
        left_join(taxonomy)%>%
        left_join(RA_OTU_all)%>%
        left_join(OTU_distribution_all)%>%
        left_join(number_of_houses_OTUS_all)

######Remark: These data were used for the Table 2 

#Plot of indicator species with IndVal>0,5
Indc_ind_species_indoor_outdoor_summary2 %>%  
        unite(OTU_sp,c("OTU_id","genus","species"))%>% 
        filter(Ind_sp_location!="other",stat>0.5) %>% 
        ggplot(aes(x = Ind_sp_location,y = reorder(OTU_sp,stat),fill = stat)) +
        geom_tile() +
        scale_x_discrete(position = "top")+
        scale_fill_viridis_c(option = "inferno")+
        theme(axis.text.x = element_text(size = 12),
              axis.text.y = element_text(size = 12),
              axis.title.x = element_blank(),
              axis.title.y = element_blank())+
        labs(title = "Indicator species - Indoor_Outdoor; IndVal>0.5")


## Figure with taxonomic assignment (major Orders) of the indicator species (Figure S10)

#Data preparation
indic_table_order<-Indc_ind_species_indoor_outdoor_summary2 %>% 
        select(OTU_id,c(7:14))%>%
        mutate(count=1)%>%
        group_by(Ind_sp_location, order) %>% 
        summarise(sum_order=sum(count)) %>% 
        filter(!is.na(order)) %>% 
        group_by(Ind_sp_location) %>% 
        mutate(Relative_abundance=sum_order/sum(sum_order)) %>% 
        mutate(new_order=ifelse(Relative_abundance<0.025, "z_others", order)) 
#Remark: Removing the NA filter we can check the percents of unidentified OTUs (at order level)
#Non-identified: 9.1% and 22.7% for indoor and outdoor indicators, respectively

# Number of levels for color palette
colourCount <- length(unique(indic_table_order$new_order)) 

#Plot
indic_table_order %>% 
        ggplot(aes(x=Ind_sp_location, y=Relative_abundance, fill=new_order))+
        geom_bar(stat = "identity")+
        coord_flip()+
        scale_fill_manual(values = colorRampPalette(brewer.pal(12, 
                                                               "Paired"))(colourCount))+
        guides(fill = guide_legend(reverse=TRUE))+ # Reverse the order of legend elements
        ylab("Relative abundance")+
        xlab("Environment")+
        labs(fill="orders >2.5%")+
        labs(title= "order", subtitle = "Environment")+
        theme_bw()
####Remark: This plot was used for Figure S10


#########################
# 2.2. Indicator Species analysis - Living room vs. Bathroom (INDOOR dataset)
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

#Select metadata for indoor samples
indoor_metadata=data_for_ord_indoor %>% 
        select(-c(2:19)) %>%
        distinct(sample, .keep_all = T)

#Indicator species
test_location_ind<-multipatt(indoor_otutable, indoor_metadata$location) ### same order of samples in otutable and metadata!!

#Summary
Indc_species_location_ind_summary<-summary.multipatt(test_location_ind)

#Multilevel pattern analysis
#---------------------------
#Association function: IndVal.g
#Significance level (alpha): 0.05

#Total number of species: 5817
#Selected number of species: 98 
#Number of species associated to 1 group: 98 
#List of species associated to each combination: 
#       Group bathroom  #sps.  53
#       Group central  #sps.  45 

#Prepare results and add complementary information of OTUs (taxonomy, FUNGuild, Relative abundance, distribution as #of samples and #of houses)
Indc_species_location_ind_summary2<-test_location_ind$sign %>% 
        rownames_to_column("OTU") %>% rename(OTU_id="OTU")%>%
        filter(p.value<=0.05) %>% #filter significant values
        mutate(Ind_sp_location=ifelse(index==1,"Bathroom",
                                      ifelse(index==2,"Central",
                                             "Other"))) %>% 
        left_join(taxonomy)%>%
        left_join(RA_OTU_all)%>%
        left_join(OTU_distribution_all)%>%
        left_join(number_of_houses_OTUS_all)

#Plot of indicator species with IndVal>0.2
Indc_species_location_ind_summary2 %>% 
        unite(OTU_sp,c("OTU_id","genus","species"))%>% 
        filter(stat>0.2) %>% 
        ggplot(aes(x = Ind_sp_location,y =reorder(OTU_sp,stat),fill = stat)) +
        geom_tile() +
        scale_x_discrete(position = "top")+
        scale_fill_viridis_c(option = "inferno")+
        theme(axis.text.x = element_text(size = 12, angle = 30, hjust = 0),
              axis.text.y = element_text(size = 12),
              axis.title.x = element_blank(),
              axis.title.y = element_blank())+
        labs(title = "Indicator species - House compartments (INDOOR); IndVal>0.2")

#Remark: Only 2 OTUs associated with bathrooms with IndVal>0.5



#################### END OF THE SCRIPT ##########################
