########################################################################################################

### RECIPROCAL TRANSPLANT EXPERIMENT ###

setwd("/Users/billymoore/Documents/Conceptacles Paper/Publication/Proc B Submission/Manuscript Revision_28.03.21/Data")
RTransplant <- read.csv(file = "Reciprocal Transplant_Updated.csv", header = TRUE)
attach(RTransplant)

library(ggpubr)
library(Rmisc)
library(ggplot2)
library(plyr)
library(dplyr)
library(multcomp)
library(glmmTMB)

# Changes the A, OA treatment abbreviations #
RTransplant <- RTransplant %>%
mutate(Transplant.Mean.pH = revalue(Transplant.Mean.pH, c("A"="Present Day", "OA"="Ocean Acidification")))

# Normalise Data, by the area of the photograph (mm2)

RTransplant <- RTransplant %>%
mutate(Conceptacle.number = Conceptacle.number / 16)

############################################################################################

# Convert label names from abbreviations (C,F,T,S) to words for site and variability #
Site.labs <- c("Shell Island", "Tallon Island")
names(Site.labs) <- c("S", "T")
Variability.labs <- c("Low Variability pH Treatment", "High Variability pH Treatment")
names(Variability.labs) <- c("C", "F")

############################################################################################

### CONCEPTACLE ABUNDANCE ### 

### Manuscript Figure 3: Conceptacle Abundance by Generation 2-6 Treatment and Transplant Treatment ###

(Num_Transplant <- ggplot(RTransplant, aes(x=factor(Transplant.Mean.pH), y=Conceptacle.number, fill=Gen.2.6.MeanTreatmentpH)) + 
    geom_boxplot()+
    theme_classic()+
    xlab("Transplant Treatment") +
     ylab(bquote('Conceptacle Abundance (' ~ conceptacles~ mm^2*')')) +
     scale_y_continuous(breaks = seq(0, 6, by = 1)) +
    labs(fill="Generation 2-6 Treatment")+
        theme(legend.position = c(0.80, 0.88), legend.title = element_text(size=12), 
        legend.text = element_text(size=12), axis.text=element_text(size=12),
        axis.title=element_text (size=13)) +
    scale_fill_manual(values=c('#CCCCCC','#E69F00'),
                 labels=c("Present Day", "Ocean Acidification")))

 
### Conceptacle Abundance Reciprocal Transplant Model ### 

TransplantNumPModQ <- glmmTMB(Conceptacle.number ~ Transplant.Mean.pH * Gen.2.6.MeanTreatmentpH * Site.of.origin * Variability + (1|Header) + (1|Waterbath), family = compois)
summary(TransplantNumPModQ)
TransplantNumPModQTab <- summary(aov(TransplantNumPModQ));TransplantNumPModQTab

#-------------------------------------------------------------------------------#

### Wilcoxon Rank Sum Pairwise Comparisons for Conceptacle Abundance Reciprocal Transplant ###

# OA - OA vs A - OA numbers
OAOA <- RTransplant %>%
    filter(Transplant.Mean.pH == "Ocean Acidification", Gen.2.6.MeanTreatmentpH == "OA" ) %>%
    pull(Conceptacle.number)

AOA <- RTransplant %>%
    filter(Transplant.Mean.pH == "Ocean Acidification", Gen.2.6.MeanTreatmentpH == "A" ) %>%
    pull(Conceptacle.number)

# A - A vs OA - A numbers
AA <- RTransplant %>%
    filter(Transplant.Mean.pH == "Present Day", Gen.2.6.MeanTreatmentpH == "A" ) %>%
    pull(Conceptacle.number)

OAA <- RTransplant %>%
    filter(Transplant.Mean.pH == "Present Day", Gen.2.6.MeanTreatmentpH == "OA" ) %>%
    pull(Conceptacle.number)

wilcox.test(AOA, OAOA, alternative = "two.sided")

wilcox.test(AA, OAA, alternative = "two.sided")

wilcox.test(AA, AOA, alternative = "two.sided")


#################################################################################

### CONCEPTACLE SIZE ### 

### Supplementary Figure 10: Conceptacle Diameter by Generation 2-6 Treatment and Transplant Treatment ###

(Len_Transplant <- ggplot(RTransplant, aes(x=factor(Transplant.Mean.pH), y=Conceptacle.Diameter, fill=Gen.2.6.MeanTreatmentpH)) + 
     geom_boxplot()+
     theme_classic()+
     xlab("Transplant Treatment") +
     ylab("Conceptacle Diameter (Î¼m)")+
     labs(fill="Generation 2-6 Treatment")+
     theme(legend.position = c(0.7, 0.9), legend.title = element_text(size=12), 
           legend.text = element_text(size=12), axis.text=element_text(size=12),
           axis.title=element_text (size=13)) +
     scale_fill_manual(values=c('#CCCCCC','#E69F00'),
                       labels=c("Present Day", "Ocean Acidification")))

#-------------------------------------------------------------------------------#

### Conceptacle diameter reciprocal transplant Quassipoisson Model ###

# Model:
TransplantLenPModQ <- glmmTMB(Conceptacle.Diameter ~ Transplant.Mean.pH * Gen.2.6.MeanTreatmentpH * Site.of.origin * Variability + (1|Header) + (1|Waterbath), family = ziGamma(link = "log"))
summary(TransplantLenPModQ)
TransplantLenPModQTab <- summary(aov(TransplantLenPModQ));TransplantLenPModQTab






