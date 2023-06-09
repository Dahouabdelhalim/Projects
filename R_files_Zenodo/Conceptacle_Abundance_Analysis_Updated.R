############################################################################################

### CONCEPTACLE ABUNDANCE  ###

setwd("/Users/billymoore/Documents/Conceptacles Paper/Publication/Proc B Submission/Manuscript Revision_28.03.21/Data")
ConcNumb <- read.csv(file = "Conceptacle Abundance_Updated.csv", header = TRUE)

library(ggpubr)
library(Rmisc)
library(ggplot2)
library(glmmTMB)

############################################################################################

# Normalise Data, by the area of the photographs (mm2)

ConcAbundance = ConcNumb %>%
   mutate(Conceptacle.number = Conceptacle.number / 16)

attach(ConcAbundance)
############################################################################################

### FIGURES ###

# Convert label names from abbreviations (C,F,T,S) to words for site and variability #
Site.labs <- c("Shell Island", "Tallon Island")
names(Site.labs) <- c("S", "T")
Variability.labs <- c("Low Variability pH Treatment", "High Variability pH Treatment")
names(Variability.labs) <- c("C", "F")

#-------------------------------------------------------------------------------#

### Supplementary Figure 5: Conceptacle Abundance by Generation and Mean Treatment pH ###

(Number_Treat_Gen <- ggplot(ConcAbundance, aes(x=factor(Generation), y=Conceptacle.number, fill=MeanTreatmentpH)) + 
   geom_boxplot() +
   theme_classic()+
   xlab("Generation") +
   ylab(bquote('Conceptacle Abundance (' ~ conceptacles~ mm^2*')')) +
   labs(fill="Mean Treatment pH")+
   scale_y_continuous(breaks = seq(0, 5, by = 1)) +
   theme(legend.position = c(0.15, 0.9), legend.title = element_text(size=12), 
         legend.text = element_text(size=12), axis.text=element_text(size=12),
         axis.title=element_text (size=13)) +
   scale_fill_manual(values=c('#CCCCCC','#E69F00'),
                     labels=c("Present Day", "Ocean Acidification")))

#-------------------------------------------------------------------------------#
  
### Figure 2: Conceptacle Abundance by Generation, Mean Treatment pH, Site and Variability ###
  
(split_plot_Site_Var <- ggplot(aes(factor(Generation), Conceptacle.number, fill = MeanTreatmentpH), data = ConcAbundance) + 
   geom_boxplot() + 
   facet_grid(Site.of.origin ~ pHVariability) + 
   labs(fill="Mean Treatment pH:") +
   xlab("Generation") + 
   ylab(bquote('Conceptacle Abundance (' ~ conceptacles~ mm^2*')')) +
   scale_y_continuous(breaks = seq(0, 5, by = 1)) +
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

### GLM MODEL ###

NumQPMod <- glmmTMB(Conceptacle.number ~ Generation * MeanTreatmentpH * Site.of.origin * pHVariability + (1|Header) + (1|Waterbath), family = compois)
summary(NumQPMod)
NumQPModTab <- summary(aov(NumQPMod));NumQPModTab

############################################################################################





