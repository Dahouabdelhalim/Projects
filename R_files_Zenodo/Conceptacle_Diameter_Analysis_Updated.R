############################################################################################

### CONCEPTACLE DIAMETER ###

setwd("/Users/billymoore/Documents/Conceptacles Paper/Publication/Proc B Submission/Manuscript Revision_28.03.21/Data")
ConcDiameter <- read.csv(file = "Conceptacle Diameter_Updated.csv", header = TRUE)
attach(ConcDiameter)

library(ggpubr)
library(Rmisc)
library(ggplot2)
library(glmmTMB)

############################################################################################

### FIGURES ###

# Convert label names from abbreviations (C,F,T,S) to words for site and variability #
Site.labs <- c("Shell Island", "Tallon Island")
names(Site.labs) <- c("S", "T")
Variability.labs <- c("Low Variability pH Treatment", "High Variability pH Treatment")
names(Variability.labs) <- c("C", "F")

#-------------------------------------------------------------------------------#

### Supplementary Figure 8: Conceptacle Length by Generation and Mean Treatment pH ###

(Len_Number_Treat_Gen <- ggplot(ConcDiameter, aes(x=factor(Generation), y=Mean.Conceptacle.Diameter.um, fill=MeanTreatmentpH)) + 
    geom_boxplot()+
    theme_classic()+
    xlab("Generation") +
    ylab("Conceptacle Diameter (μm)")+
    labs(fill="Mean Treatment pH")+
    theme(legend.position = c(0.85, 0.9), legend.title = element_text(size=12), 
          legend.text = element_text(size=12), axis.text=element_text(size=12),
          axis.title=element_text (size=13)) +
    scale_fill_manual(values=c('#CCCCCC','#E69F00'),
                      labels=c("Present Day", "Ocean Acidification")))

#-------------------------------------------------------------------------------#
 
### Supplementary Figure 9: Conceptacle Length by Generation, Mean Treatment pH, Site and Variability ###
  
(Len_split_plot_Site_Var <- ggplot(aes(factor(Generation), Mean.Conceptacle.Diameter.um, fill = MeanTreatmentpH), data = ConcDiameter) + 
   geom_boxplot() + 
   facet_grid(Site.of.origin ~ pHVariability) + 
   labs(fill="Mean Treatment pH:") +
   xlab("Generation") + 
   ylab("Conceptacle Diameter (μm)") +
   theme_bw()+
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
         legend.position = "top", legend.title = element_text(size=12), 
         legend.text = element_text(size=12), axis.text=element_text(size=10),
         axis.title=element_text (size=13),strip.text.x = element_text(size=10), 
         strip.text.y = element_text(size=10)) +
    theme(strip.background = element_rect(fill="white", size=0, color="black")) +
    theme(strip.text = element_text(face="bold",size=10, color="black")) +
    theme(panel.spacing = unit(1.2, "lines"))+
   facet_grid(Site.of.origin ~ pHVariability, labeller = labeller(Site.of.origin = Site.labs, pHVariability = Variability.labs))+
   scale_fill_manual(values=c('#CCCCCC','#E69F00'),
                     labels=c("Present Day", "Ocean Acidification")))

############################################################################################

### MODEL ###

LenQPMod <- glmmTMB(Mean.Conceptacle.Diameter.um ~ Generation * MeanTreatmentpH * Site.of.origin * pHVariability + (1|Header) + (1|Waterbath), family = ziGamma(link = "log"))   
summary(LenQPMod)
summary(aov(LenQPMod))

############################################################################################



