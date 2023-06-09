########## packges load ###################################################################################
library(vegan)
library(car)
library (DESeq2); packageVersion("DESeq2")
library(Rfast) #to take row maximums
library (ggthemes)
library(data.table); packageVersion("data.table")
library(adaptiveGPCA)
library(ggrepel)
library(phyloseq)
library(plyr)
############ load clean and rerified data for analysis ###############################################################
#data_in <- readRDS(file="Crane_CleanData_rarefied_1028.rds")
data_in <- readRDS(file="Crane_CleanData_rarefied_8000.rds")

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  #====================================================================================================================
###########################  Look at the data ##########################################################################
ps4<-data_in
taxa_names(ps4) <- paste0("SV", seq(ntaxa(ps4))) # rename to short names
OTU1 = as(otu_table(ps4), "matrix")
if(taxa_are_rows(ps4)){OTU1 <- t(OTU1)}
OTUdf = as.data.frame(OTU1)
write.table(OTUdf,file="data.nameNew.csv",sep=",") # creates table with abundances per sample of each OTU
# create a taxanomic table 
kl<-tax_table(ps4)
KLdf = as.data.frame(kl)
write.table(KLdf,file="TaxonomicTableShortNamesNew.csv",row.names = FALSE,sep=",") 

# sequences 
T<-taxa_names(data_in)
TT<-as.data.frame(T)
write.table(TT,file="SequencesNew.csv",sep=",")
# If you wabt the whole sequence use the "data_in" object
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++ Add field per Control Not control +++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Collect the relevant data
Week <- data_in@sam_data[["WeekSinceBreeding"]]
Loc <- data_in@sam_data[["Location"]]
FeedType <- data_in@sam_data[["Feeding_type"]]

#+++++++++++++++++ Add field per Control Not control Per week of sampling +++++++++++++++++++++++++++++++++
ControlGroupWeek <- Week
ControlGroupWeek[Week >2 & Week <7] <- 0 # Pre-Migration (20-Aug to 20 Sep)
ControlGroupWeek[Week >10 & Week <16] <- 1 # Fall (15-Oct to 16-Nov) => includes hula and alternative

ControlGroupWeek[FeedType == "fields" & Week >20 & Week <29] <- 2 # Winter (27-Dec to 10-Feb) => Fields only
ControlGroupWeek[FeedType == "feeding_station" & Week >20 & Week <29] <- 3 # Winter (27-Dec to 10-Feb) => feeding satation Only

ControlGroupWeekF <- factor(ControlGroupWeek, labels = c("Pre_migration", "Fall",
                                                         "Winter_fields","Winter_feeding_station" ))

data_in@sam_data[["ControlGroupWeek"]]<-ControlGroupWeekF

##--- Add column with discription including control---------------------------------------
ControlGroupWeek<- data_in@sam_data[["ControlGroupWeek"]]
Week<- data_in@sam_data[["WeekSinceArrival"]]
Loc <- data_in@sam_data[["Location"]]

GroupWithControl <- data_in@sam_data[["Week"]]
GroupWithControl <-0
GroupWithControl[ControlGroupWeek== "Pre_migration"] <- 0
GroupWithControl[ControlGroupWeek== "Fall" & Week != 6] <- 1
GroupWithControl[Loc == "Hula" & Week == 6] <- 2
GroupWithControl[Loc == "Emek_Israel" & Week == 6] <-3
GroupWithControl[ControlGroupWeek== "Winter_fields"] <- 4
GroupWithControl[ControlGroupWeek== "Winter_feeding_station"] <- 5

GroupWithControl <- factor(GroupWithControl, labels = c("Pre_migration","Fall",
                                                        "Control_Hula", "Control_Other",
                                                        "Winter_fields","Winter_feeding_station"))

data_in@sam_data[["GroupWithControl"]]<-GroupWithControl
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
temp<-data.frame(sample_data(data_in)) # crate a data frame so we can look that is is correct
# Sample size
tapply(temp$WeekSinceBreeding, temp$ControlGroupWeek, length)

#Pre_migration                   Fall          Winter_fields   Winter_feeding_station 
#--    36                     42 (33+9)            32                     38 

# save
saveRDS(data_in, file="data_in_rarefied_8000.rds")
# if we want to write a table
a1<-data.frame(data_in@sam_data)
write.csv(a1, file = "MetaToSeeTemp.csv")