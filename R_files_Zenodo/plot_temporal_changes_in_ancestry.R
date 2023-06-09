#load the required library 
library(ggplot2)

#read in data files for heterozygosity rate, ancestry-informative markers (AIM), and the index of admixture
#before plotting, we summarized for each population the mean and standard error for all three metrics

#estimates of heterozygosity rate for samples with temporal replicates were obtained using the "per_sample_heterozygosity.sh" script, and are listed in the heterozygosity_temporal.txt output file; additional information for each sample was added for plotting ("Heterozygosity_percent", "Population", "Year","Population_Year"); this annotated file was saved as a .csv
Heterozygosity_data<-read.csv(file="heterozygosity_temporal.csv",header = TRUE, sep = ,)

#estimates of allele frequency for Western-Cuba ancestry informative markers (AIMs) were obtained using the "identify_WC_AIMs.sh" script, and are listed in the WC_AF.txt output file; additional information for each sample was added for plotting ("Population", "Year","Population_Year"); this annotated file was saved as a .csv
WC_AF_data<-read.csv(file="WC_AF.csv",header = TRUE, sep = ,)

#the index of admixture (HA) was calculated following Keller & Taylor (2010), in excel; a .csv with the separate components of the HA metric is provided as part of this data archive (HA.csv)
#the individual STRUCTURE membership coefficients used in these calculations were obtained using a STRUCTURE analysis at K=6, with prior population information
HA_data<-read.csv(file="HA.csv",header = TRUE, sep = ,)

#order populations and sampling times
WC_AF_data$Population_Year = factor(WC_AF_data$Population_Year, levels=c("MON2003", "MON2018", "MIA2003", "MIA2018", "ALA2003", "ALA2018", "ORA2003", "ORA2018", "TAM2003", "TAM2018", "STP2003", "STP2018"))
HA_data$Population_Year = factor(HA_data$Population_Year, levels=c("MON2003", "MON2018", "MIA2003", "MIA2018", "ALA2003", "ALA2018", "ORA2003", "ORA2018", "TAM2003", "TAM2018", "STP2003", "STP2018"))
Heterozygosity_data$Population_Year = factor(Heterozygosity_data$Population_Year, levels=c("MON2003", "MON2018", "MIA2003", "MIA2018", "ALA2003", "ALA2018", "ORA2003", "ORA2018", "TAM2003", "TAM2018", "STP2003", "STP2018"))

#plot AIM allele frequency per population and year
ggplot(WC_AF_data, aes(Population_Year, AF_WC_AIMs))+
  geom_jitter(aes(colour=as.factor(Year)), size=2, width=0.2)+
  geom_errorbar(aes(ymin=Mean_AF_WC_AIMs-SE_AF_WC_AIMs, ymax=Mean_AF_WC_AIMs+SE_AF_WC_AIMs), width=0)+
  geom_point(aes(x=Population_Year, y=Mean_AF_WC_AIMs, colour=as.factor(Year)), size=5)+
  geom_point(aes(x=Population_Year, y=Mean_AF_WC_AIMs), colour="black",pch=21, size=5)+
  theme_bw()+
  theme(axis.title=element_text(size=12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  scale_colour_manual(values = c("grey70","#2d7a4bfe"))

#plot the index of admixture per population and year
ggplot(HA_data, aes(Population_Year, HA))+
  geom_jitter(aes(colour=as.factor(Year)), size=2, width=0.2)+
  geom_errorbar(aes(ymin=Mean_HA-SE_HA, ymax=Mean_HA+SE_HA), width=0)+
  geom_point(aes(x=Population_Year, y=Mean_HA, colour=as.factor(Year)), size=5)+
  geom_point(aes(x=Population_Year, y=Mean_HA), colour="black",pch=21, size=5)+
  theme_bw()+
  theme(axis.title=element_text(size=12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  scale_colour_manual(values = c("grey70","#2d7a4bfe"))

#plot heterozygosity rate per population and year
ggplot(Heterozygosity_data, aes(Population_Year, Heterozygosity_percent))+
  geom_jitter(aes(colour=as.factor(Year)), size=2, width=0.2)+
  geom_errorbar(aes(ymin=Mean_heterozygosity-SE_heterozygosity, ymax=Mean_heterozygosity+SE_heterozygosity), width=0)+
  geom_point(aes(x=Population_Year, y=Mean_heterozygosity, colour=as.factor(Year)), size=5)+
  geom_point(aes(x=Population_Year, y=Mean_heterozygosity), colour="black",pch=21, size=5)+
  theme_bw()+
  theme(axis.title=element_text(size=12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  scale_colour_manual(values = c("grey70","#2d7a4bfe"))+
  ylim(5,12)

#build linear models with AIM allele frequency, index of admixture, or heterozygosity rate as the response, and population IDs and sampling year (2003 / 2018) as the predictor.
lm_HA<-lm(HA_data$HA~HA_data$Population+as.factor(HA_data$Year))
anova(lm_HA)

lm_het<-lm(Heterozygosity_data$Heterozygosity~Heterozygosity_data$Population+as.factor(Heterozygosity_data$Year))
anova(lm_het)

lm_AIM<-lm(WC_AF_data$AF_WC_AIMs~WC_AF_data$Population+as.factor(WC_AF_data$Year))
anova(lm_AIM)


#plot allele frequency histograms for markers that were polymorphic in 2003 and became fixed in 2018 (<population>.fixed_2018_sites.txt), as well as for those that were polymorphic in 2003, irrespective of 2018 allele frequencies (<population>.all_2003_polymorphic_sites.txt)
#these allele frequency input files were obtained using the "AF_2013_2018.sh" script
#read in input files for each population
ALA_fixed_2018_sites<-read.csv(file="ALA.fixed_2018_sites.txt",header = FALSE, sep = '\\t')
ALA_all_2003_polymorphic_sites<-read.csv(file="ALA.all_2003_polymorphic_sites.txt",header = FALSE, sep = '\\t')

ORA_fixed_2018_sites<-read.csv(file="ORA.fixed_2018_sites.txt",header = FALSE, sep = '\\t')
ORA_all_2003_polymorphic_sites<-read.csv(file="ORA.all_2003_polymorphic_sites.txt",header = FALSE, sep = '\\t')

STP_fixed_2018_sites<-read.csv(file="STP.fixed_2018_sites.txt",header = FALSE, sep = '\\t')
STP_all_2003_polymorphic_sites<-read.csv(file="STP.all_2003_polymorphic_sites.txt",header = FALSE, sep = '\\t')

TAM_fixed_2018_sites<-read.csv(file="TAM.fixed_2018_sites.txt",header = FALSE, sep = '\\t')
TAM_all_2003_polymorphic_sites<-read.csv(file="TAM.all_2003_polymorphic_sites.txt",header = FALSE, sep = '\\t')

MIA_fixed_2018_sites<-read.csv(file="MIA.fixed_2018_sites.txt",header = FALSE, sep = '\\t')
MIA_all_2003_polymorphic_sites<-read.csv(file="MIA.all_2003_polymorphic_sites.txt",header = FALSE, sep = '\\t')

MON_fixed_2018_sites<-read.csv(file="MON.fixed_2018_sites.txt",header = FALSE, sep = '\\t')
MON_all_2003_polymorphic_sites<-read.csv(file="MON.all_2003_polymorphic_sites.txt",header = FALSE, sep = '\\t')


#make allele frequency histograms for each marker category (example here given for the MON population)
ggplot(MON_fixed_2018_sites, aes(x=V4))+
  geom_histogram(position="identity", fill="black")+
  theme_bw()+
  theme(axis.title=element_text(size=12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position="none")

ggplot(MON_all_2003_polymorphic_sites, aes(x=V4))+
  geom_histogram(position="identity", fill="black")+
  theme_bw()+
  theme(axis.title=element_text(size=12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position="none")

