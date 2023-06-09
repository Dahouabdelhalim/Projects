
#########################Filename: bioinformatics_script_filtering_otutable.R

#Publication: Analyzing indoor mycobiomes through a large-scale citizen science study in Norway
#Authors: Pedro M. Martin-Sanchez, Eva Lena F. Estensmo, Luis N. Morgado, Sundy Maurice, Ingeborg B. Engh, Inger Skrede and Håvard Kauserud

##Related to the sections: 
        #Materials and Methods- 2.4. Bioinformatics pipeline
        #Results- 3.1. Data features and overall fungal diversity
        #Supplemental Information

## Goals: 
#(1) Quality filtering of OTU table: 
        #Removal of OTUs with  <10 reads,
        #Removal of OTUs with <70% identity percent in taxonomic assignment
        #Filtering OTUs assigned to the kingdom Fungi
#(2) Assessing positive controls (mock community samples) and technical replicates (dust samples)
#(3) Rarefying the final fungal OTU table at 2,000 reads per sample

#Inputs: Complete OTU table after bioinformatics analyses including post-clustering curation by LULU ("edited_lulu_curated_otutable.xlsx")
#Outputs: 
        #Filtered OTU table for all organisms ("filtered_otutable_all.xlsx")
        #Filtered OTU table for fungi ("filtered_otutable_fungi.xlsx")
        #Final clean and rarefied OTU table ("OTU_table_rar2000x10.xlsx")
        #Plot for the Figure S4 (NMDS ordination plot for the 17 technical replicates)
        #Plot for the FigureS5 (OTUs vs. reads plot)

#########################
#PREPARATIONS
#########################

#Load packages

library(openxlsx)
library(tidyverse)
library(vegan)

#Set your path
setwd("C:/WP1_data_analyses/script_preparation/")

#Load the complete OTU table (including all organisms): 11,625 OTUs (rows) and 857 study samples (columns).
OTU_tab<-read.xlsx("edited_lulu_curated_otutable.xlsx")
#Remark: The table includes other data: sizes of OTUs (OTU_abund), lenght of their representative sequences and data related to the taxonomic assignment.   

# Transform the column "ident_perc" as numeric because is chr
OTU_tab$ident_perc = as.numeric (OTU_tab$ident_perc) 
is.numeric (OTU_tab$ident_perc) #double check

# Delete 5 samples that should be initially excluded of the dataset and delete the OTUs without reads afterwards 
OTU_tab1 <- OTU_tab %>%
  select(-c(contains ("cont")))%>%
  gather(sample, read_counts, -c(1:13))%>%
  filter(read_counts>0) %>% 
  spread(sample,read_counts)
OTU_tab1[is.na(OTU_tab1)]<-0 #transform all NA to 0

#Check numbers of rows (11,614 OTUs) and columns (865 = 13+ 852 samples including including dust samples, controls (extraction blanks, PCR negatives and mock samples) and technical replicates
dim(OTU_tab1)


#########################
#Filter 1- Removal of OTUs with <10 reads
#########################

#Filter OTUs with >= 10 reads per OTU (8,534 OTUsleft; 3,080 OTUS removed)
OTU_tab2 <- OTU_tab1 %>%  
  filter(OTU_abund >=10) 


#########################
#Filter 2- Removal of OTUs with <70% identity percent in the taxonomic assignment
#########################

# Filter OTUs with identity percent (tax) >= 70%  (8,058 OTUs left; 476 OTUS removed) 
OTU_tab3 <- OTU_tab2 %>%  
  filter(ident_perc>=70) #%>%
  #summarize(how_many_samples=n()) #To check different cutoffs


#########################
#Checking members of mock community (samples Mock1-Mock9)
#########################
# Mock samples includes three OTUs:  OTU20 Mycena bellarium, OTU50 Pycnoporellus fulgens and OTU77 Inonotus hispidus

# Transforming the table to allow the next tests
OTU_tab3_transformed <- OTU_tab3 %>% 
  gather(sample, read_counts, -c(1:13)) 

# Check the distribution of OTU20 - only 10 samples: 9 mock samples + 50B (23 reads)
OTU20_samples_reads <- OTU_tab3_transformed %>%  
  filter(OTU_id=="OTU20", read_counts >0)%>%
  select(sample,read_counts)
OTU20_samples_reads

# Check the distribution of OTU50 - only 10 samples: 9 mock samples + 127B (2 reads)
OTU50_samples_reads <- OTU_tab3_transformed %>%  
  filter(OTU_id=="OTU50", read_counts >0)%>%
  select(sample,read_counts)
OTU50_samples_reads

# Check the distribution of OTU77 - only in the 9 mock samples
OTU77_samples_reads <- OTU_tab3_transformed %>%  
  filter(OTU_id=="OTU77", read_counts >0)%>%
  select(sample,read_counts)
OTU77_samples_reads

#Remark: THERE IS NO SIGNIFICANT TAG SWITCHING. So no additional filters are needed to correct this problem. 
#Remark: The three mock OTUs should be deleted later (see below).

#########################
#Assessing microbial community composition between technical replicates - including all organisms 
#########################

#Select the OTU table for the 17 technical replicates
OTU_tab_replicates=OTU_tab3 %>%
  select(OTU_id, contains("_1"), contains("_2"), contains("269S"))  %>%
  gather(sample, reads,-OTU_id) %>% # create the otutable for the replicates, removing all first 13 columns except OTU_id
  group_by(OTU_id) %>% 
  filter(reads>0) %>% 
  spread(sample,reads)
OTU_tab_replicates[is.na(OTU_tab_replicates)]<-0 #transform all NA to 0

#Transpose the table to get the proper format
OTU_tab_replicates_t=OTU_tab_replicates %>% 
  column_to_rownames("OTU_id") %>% 
  t() %>% 
  data.frame()

#Hellinger transformation
replicates_Hell<-OTU_tab_replicates_t %>% 
  decostand(method = "hellinger")

#Nonmetric multidimensional scaling plot (NMDS)
NMDS_test=metaMDS(replicates_Hell, distance = "bray", k = 2, try = 20, trymax = 20, 
                  engine = c("monoMDS"))

nmds1=NMDS_test$points[,1]
nmds2=NMDS_test$points[,2]

labels=rownames(replicates_Hell)
nmds_axis_and_samples=cbind(labels,nmds1,nmds2) %>% as.data.frame() %>% 
  separate(labels,c("plot"),"_") %>% 
  mutate(nmds1=as.numeric(as.character(nmds1))) %>% 
  mutate(nmds2=as.numeric(as.character(nmds2)))

plot(NMDS_test,c("sites"))
text(NMDS_test,c("sites"),labels)

#NMDS plot
nmds_axis_and_samples %>% ggplot(aes(x=nmds1, y=nmds2, label=plot, fill=plot))+
  geom_point()+
  geom_label()
#Remark: This is the plot for the Supplemental Figure S4

#Remark:Based on that plot, for the 17 dust samples both duplicates showed similar microbial community composition.
#We decided to keep those replicates with higher number of reads for further analyses and discard the other ones.


#########################
#Assessing Negative controls (Extraction blanks and PCR negatives) - Manually as detailed in the Supplemental file 1   
#########################
# six extraction blanks and three PCR negatives contained a relatively low number of reads, representing an average of 4.1±2.6 OTUs per negative control.
# We decided to delete two rare OTUs (< 10 reads in two samples): OTU310 and OTU2405.
# The remaining 22 OTUs were kept because they were widely distributed in the dataset and correspond to ubiquitous fungi in the built environment.


#########################
#Delete samples (controls and replicates) and the mentioned OTUs  
#########################

#Delete mock samples (9), replicates (20) and negatives (12). 
#Delete 5 OTUs (3 from mock community and 2 from negatives)
#Rename the kept replicates.
filtered_otutable_all <- OTU_tab3 %>%
  gather(sample, reads,-c(1:13)) %>% 
  group_by(OTU_id) %>% 
  filter(sample != "Mock1", sample != "Mock2", sample != "Mock3", sample != "Mock4", sample != "Mock5", sample != "Mock6", sample != "Mock7", sample != "Mock8", sample != "Mock9") %>%
  filter(sample != "269S2",sample != "342B_1",sample != "342S_1",sample != "342U_2",sample != "343B_1",sample != "343S_1", sample != "343U_2",sample != "344B_2", sample != "344S_1", sample != "358U_1",sample != "359B_1", sample !="359S_1", sample != "359U_1", sample != "88B_1") %>%
  filter(sample != "3S_oct_1", sample != "3S_oct_2", sample != "3B_oct_1",sample != "3B_oct_2", sample != "3U_oct_1", sample != "3U_oct_2") %>%
  filter(sample != "NegExt3",sample != "NegExt5",sample != "NegExt6",sample != "NegExt7", sample != "NegExt8",sample != "NegExt9",sample != "NegPCR9") %>% #negative controls
  filter(sample!="NegPCR1",sample!="NegExt1",sample!="NegExt4",sample!="NegPCR7",sample!="NegPCR8")%>%
  filter(OTU_id != "OTU20", OTU_id != "OTU50", OTU_id != "OTU77") %>% #3 mock OTUs
  filter(OTU_id != "OTU310", OTU_id != "OTU2405") %>% # 2 rare OTUs in Negatives
  filter(reads>0) %>% 
  spread(sample,reads)%>%
  rename("269S"="269S1", "342B" = "342B_2", "342S" = "342S_2","342U" = "342U_1","343B" = "343B_2", "343S"="343S_2","343U"="343U_1", "344B"="344B_1","344S"="344S_2") %>%
  rename ("358U"="358U_2","359B"="359B_2", "359S"="359S_2", "359U"="359U_2","88B"="88B_2")
filtered_otutable_all[is.na(filtered_otutable_all)]<-0 #transform all NA to 0

#Checks of sample names, #of samples and # of OTUs
colnames(filtered_otutable_all)
dim(filtered_otutable_all) #[1] 8033 OTUs (16 less) and  824 (13+ 811 samples)

# Save the filtered OTU table including all organisms  
write.xlsx(filtered_otutable_all,"filtered_otutable_all.xlsx")


#########################
#Filter 3- OTUs assigned to the kingdom Fungi  
#########################

#Filter Fungi  (7110 OTUs left (Fungi 88.5% ); 923 removed- plants and others)
filtered_otutable_fungi <- filtered_otutable_all %>%  
  filter(Domain == "d:Fungi")

#save the filtered fungal OTU table
write.xlsx(filtered_otutable_fungi,"filtered_otutable_fungi.xlsx")


#########################
#Plot OTUs vs reads  
#########################

#Preparation of the data 
OTU_vs_reads_fungi <- filtered_otutable_fungi %>%
  gather(sample, read_counts, -c(1:13)) %>% 
  mutate(PA=sign(read_counts)) %>%  
  group_by(sample) %>%
  summarise(OTUs_per_samp=sum(PA), reads_per_samp=sum(read_counts))

#Check the total number of reads
Total_number_reads = OTU_vs_reads_fungi %>% summarise (sum(reads_per_samp))
Total_number_reads #22,622,815 reads

#Plot
OTU_vs_reads_fungi %>% ggplot(aes(x=reads_per_samp,y=OTUs_per_samp))+
  geom_point()+
  geom_smooth(method="loess")
#Remark: This is the plot for the Supplemental Figure S5

#Test if they are fitting a linear model
test_lm=lm(OTUs_per_samp~reads_per_samp,OTU_vs_reads_fungi)
summary(test_lm)


#########################
#Rarefying the final fungal OTU table at 2,000 reads   
#########################

#Filter samples with >=2000 reads (removing 4 samples: 44U, 109U, 126U and 246U)
filtered_otutable_fungi_2000reads <- filtered_otutable_fungi %>% 
  gather(sample, reads, -c(1:13)) %>% 
  group_by(sample) %>% 
  mutate(seq_depth=sum(reads)) %>% 
  filter(seq_depth>=2000) %>% 
  ungroup() %>% 
  group_by(OTU_id) %>%
  filter(reads>0) %>% 
  select(-seq_depth)%>%
  spread(sample, reads)
filtered_otutable_fungi_2000reads[is.na(filtered_otutable_fungi_2000reads)]<-0 #transform all NA to 0

#Check the numbers of OTUs and samples
dim(filtered_otutable_fungi_2000reads) #[1] 7,108 OTUs (2 less) and  820 (13+ 807 samples)

#Prepare the OTU table for the rarefaction
filtered_otutable_fungi_2000reads_rarefaction = filtered_otutable_fungi_2000reads %>% 
  remove_rownames() %>% 
  select(OTU_id,c(14:820)) %>% # discard the unnecessary columns
  column_to_rownames("OTU_id") %>% 
  t() %>%
  as.data.frame() %>% 
  rownames_to_column("sample")

#Edit sample names
filtered_otutable_fungi_2000reads_rarefaction2<-filtered_otutable_fungi_2000reads_rarefaction %>%
  mutate(addingName="rar") %>% # create a new column 
  unite("newSampleName","sample","addingName") %>% 
  column_to_rownames("newSampleName") %>% 
  as.data.frame()

#Check the right format:columns are OTUs and rows are the new sample names - like "XXX_rar"
colnames(filtered_otutable_fungi_2000reads_rarefaction2)
rownames(filtered_otutable_fungi_2000reads_rarefaction2)

# Rarefying at 2,000 reads 10 times and merge the resulting tables
res <- lapply(as.list(1:10), function(x) as.data.frame(rrarefy(filtered_otutable_fungi_2000reads_rarefaction2, sample=2000)))  ### change to 10 in ...as.list(1:2)
X <- lapply(res, rownames_to_column, "sample")
X <- lapply(res, function(x) t(x))
test.list_merge10<-do.call("cbind",X)

test.list_merge_melting<-test.list_merge10 %>% as.data.frame() %>% 
  rownames_to_column("OTU_id") %>% 
  gather(sample,rar_read_abund, -OTU_id) %>% 
  separate(sample,c("rarSample","replicate"),"_") %>% 
  select(-replicate) %>%   
  group_by(OTU_id,rarSample) %>% 
  summarize(newRar_reads=median(rar_read_abund))

test.list_merge_melting_wide<-test.list_merge_melting %>% 
  spread(OTU_id,newRar_reads) %>% 
  dplyr::rename(sample=rarSample)

# Delete all OTUS with 0 rarRead from the rarefied table (6632 instead 7108 OTUS)
clean_rar_table_fungi<-test.list_merge_melting_wide %>% 
  gather(OTU_id, rarRead, -sample) %>%
  filter(rarRead>0)%>%
  spread(OTU_id, rarRead)
clean_rar_table_fungi[is.na(clean_rar_table_fungi)]<-0 #transform all NA to 0

#Save the rarefied table
write.xlsx(clean_rar_table_fungi,"OTU_table_rar2000x10.xlsx")  


############ END OF THE SCRIPT #################