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
data_in <- readRDS(file="data_in_rarefied_8000.rds")

# Sample size
temp<-data.frame(sample_data(data_in))
tapply(temp$WeekSinceBreeding, temp$ControlGroupWeek, length)
# Pre_migration   Fall   Winter_fields   Spring_fields 
#     36            42            32            38 

tapply(temp$WeekSinceBreeding, temp$GroupWithControl, length)

###########################  Obundance table  ##########################################################################

# Most abundant Phylum
glom <- tax_glom(data_in, taxrank = 'Phylum')
ps2_rel = transform_sample_counts(glom, function(x){x / sum(x)}) # convert to relative obundance
# create dataframe from phyloseq object
dat <- psmelt(ps2_rel)

# convert Phylum to a character vector from a factor
dat$Phylum <- as.character(dat$Phylum)
# convert to relative abundance - per season
aa_Phyl<-plyr::ddply(dat, ~Phylum*ControlGroupWeek, summarize, 
                     median.abundance = round(median(Abundance), digits = 3),
                     mean.abundance = round(mean(Abundance), digits = 3)*100,
                     sd.abundance=round(sd(Abundance), digits = 3)*100,
                     se.abundance = round(sd(Abundance)/(length(Abundance))^0.5, digits = 3)*100,
                     N = length(Abundance))

PhylaMain_season <- aa_Phyl[(aa_Phyl$Phylum=="Firmicutes" | aa_Phyl$Phylum=="Actinobacteria" |
                               aa_Phyl$Phylum=="Epsilonbacteraeota" |
                               aa_Phyl$Phylum=="Fusobacteria"| aa_Phyl$Phylum=="Proteobacteria" |
                               aa_Phyl$Phylum=="Tenericutes" | aa_Phyl$Phylum=="Cyanobacteria" |
                               aa_Phyl$Phylum=="Bacteroidetes"),]

# all
aa_Phyl<-plyr::ddply(dat, ~Phylum, summarize, 
                     median.abundance = round(median(Abundance), digits = 3),
                     mean.abundance = round(mean(Abundance), digits = 3)*100,
                     sd.abundance=round(sd(Abundance), digits = 3)*100,
                     se.abundance = round(sd(Abundance)/(length(Abundance))^0.5, digits = 3)*100)

PhylaMain_all <- aa_Phyl[(aa_Phyl$Phylum=="Firmicutes" | aa_Phyl$Phylum=="Actinobacteria" |
                            aa_Phyl$Phylum=="Epsilonbacteraeota" |
                            aa_Phyl$Phylum=="Fusobacteria"| aa_Phyl$Phylum=="Proteobacteria" |
                            aa_Phyl$Phylum=="Tenericutes" | aa_Phyl$Phylum=="Cyanobacteria" |
                            aa_Phyl$Phylum=="Bacteroidetes"),]

# reshape
library(reshape2)
PhylaMain_season1<-PhylaMain_season[,c(1,2,4,6)]
PhylaMain_seasoneshape_se<-dcast(PhylaMain_season1,Phylum~ControlGroupWeek,value.var="se.abundance")
PhylaMain_seasoneshape_mean<-dcast(PhylaMain_season1,Phylum~ControlGroupWeek,value.var="mean.abundance")

PhylaMain_all1_mean<-PhylaMain_all[,c(1,3)]
PhylaMain_all1_se<-PhylaMain_all[,c(1,5)]
All_se<-cbind(PhylaMain_seasoneshape_se,PhylaMain_all1_se)
All_mean<-cbind(PhylaMain_seasoneshape_mean,PhylaMain_all1_mean)

write.csv(All_mean , "PhylaMainFinal.csv")
write.csv(All_se , "PhylaMainFinal_se.csv")

# ===================================================================================================
############################ Core microbiome  ##################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# for all cranes together
#--- create the prevelnce and abundance table ------------------------------------

#--(1) abundance with SE
ps2_rel = transform_sample_counts(data_in, function(x){x / sum(x)}) # convert to relative obundance
dat <- psmelt(ps2_rel)
dat$OTUID <- as.character(dat$OTUID)
aa_OTU<-plyr::ddply(dat, ~OTUID, summarize, 
                  median.abundance = round(median(Abundance), digits = 3),
                  mean.abundance = round(mean(Abundance), digits = 3)*100,
                  sd.abundance=round(sd(Abundance), digits = 3)*100,
                  se.abundance = round(sd(Abundance)/(length(Abundance))^0.5, digits = 3)*100,
                  N = length(Abundance))
aa_OTU$OTUID <-as.numeric(aa_OTU$OTUID)


prev0 = apply(X = otu_table(data_in),
              MARGIN = ifelse(taxa_are_rows(data_in), yes = 1, no = 2),
              FUN = function(x){sum(x > 0)})
prevdf = data.frame(Prevalence = prev0,
                    TotalAbundance = taxa_sums(data_in),
                    tax_table(data_in))

prevdf$OTUID <-as.numeric(prevdf$OTUID)
NumberOfSamples<-nrow(sample_data(data_in))
prevdf$RelAbund<-(prevdf$TotalAbundance/sum(prevdf$TotalAbundance))*100
prevdf$RelPrev<-(prevdf$Prevalence/NumberOfSamples)*100

## Joim together
Abund_Prev <- left_join(aa_OTU, prevdf, 
                           by = c("OTUID" = "OTUID"))
Abund_Prev <- Abund_Prev[order(-Abund_Prev$RelPrev),]
write.csv(Abund_Prev, file = "Abundance_SE&Prevalence.csv")
# number of ASVs present in less than 5% of the birds
nrow(prevdf[prevdf$RelPrev<0.05,])/nrow(prevdf)


# calcuale percent per class
prevdf_Firmicutes<-prevdf[prevdf$Phylum=="Firmicutes",]
PerClass_Firmicutes <- plyr::ddply(prevdf_Firmicutes, "Class", function(df1){cbind(mean(df1$TotalAbundance),sum(df1$TotalAbundance))}) 
PerClass_Firmicutes$percent <- (100*PerClass_Firmicutes$'2')/sum(as.numeric(PerClass_Firmicutes$"2"))

prevdf_Actinobacteria<-prevdf[prevdf$Phylum=="Actinobacteria",]
PerClass_Actinobacteria <- plyr::ddply(prevdf_Actinobacteria, "Class", function(df1){cbind(mean(df1$TotalAbundance),sum(df1$TotalAbundance))}) 
PerClass_Actinobacteria$percent <- (100*PerClass_Actinobacteria$'2')/sum(as.numeric(PerClass_Actinobacteria$"2"))

prevdf_Proteobacteria<-prevdf[prevdf$Phylum=="Proteobacteria",]
PerClass_Proteobacteria <- plyr::ddply(prevdf_Proteobacteria, "Class", function(df1){cbind(mean(df1$TotalAbundance),sum(df1$TotalAbundance))}) 
PerClass_Proteobacteria$percent <- (100*PerClass_Proteobacteria$'2')/sum(as.numeric(PerClass_Proteobacteria$"2"))
#-----------------------------------------------------------
# Most abundant Genus
data_in_genera = tax_glom(data_in, "Genus", NArm = TRUE) # collapse to genus level
ps2_rel = transform_sample_counts(data_in_genera, function(x){x / sum(x)}) # convert to relative obundance
dat <- psmelt(ps2_rel)
dat$Genus <- as.character(dat$Genus)
aa_Genus<-plyr::ddply(dat, ~Genus, summarize, 
                    median.abundance = round(median(Abundance), digits = 3),
                    mean.abundance = round(mean(Abundance), digits = 3)*100,
                    sd.abundance=round(sd(Abundance), digits = 3)*100,
                    se.abundance = round(sd(Abundance)/(length(Abundance))^0.5, digits = 3)*100,
                    N = length(Abundance))

# Prevelance per genus--------------------------------------------------------------------------
prev0_genera  = apply(X = otu_table(data_in_genera),
                      MARGIN = ifelse(taxa_are_rows(data_in_genera), yes = 1, no = 2),
                      FUN = function(x){sum(x > 0)})
Tax<-tax_table(data_in_genera)
prevdf_genera  = data.frame(Prevalence = prev0_genera ,
                            Tax@.Data)

NumberOfSamples<-nrow(sample_data(data_in_genera))
prevdf_genera$Relative.Prevelence  = prevdf_genera$Prevalence/NumberOfSamples
# Two genera with lowabundance but high prevalence
aa_Gen$mean.abundance[aa_Gen$Genus=="Campylobacter"]*100 # 95% pevelance, 2%  abundance
aa_Gen$mean.abundance[aa_Gen$Genus=="Turicibacter"]*100 # 91% pevelance, 4.9%  abundance
aa_Gen$mean.abundance[aa_Gen$Genus=="Enterococcus"]*100 # 91% pevelance, 0.6%  abundance
aa_Gen$mean.abundance[aa_Gen$Genus=="Fusobacterium"]*100 # 85% pevelance, 1.6%  abundance