####################################################################################################

setwd("/Users/billymoore/Documents/Conceptacles Paper/Publication/Proc B Submission/Final Docs for Submission/Data")
AllData <- read.csv(file = "Growth and Recruitment Comparison.csv", header = TRUE)

library(ggpubr)
library(Rmisc)
library(ggplot2)
library(car)
library(dplyr)

####################################################################################################

AllData = AllData %>%
  mutate(Conceptacle.Number = Conceptacle.Number / 16)

attach(AllData)
####################################################################################################

### FIGURES ###

#--------------------------------------------------------------------------------------------#

# Create Figure Labels 
Generation.labs <- c("Generation 2", "Generation 3", "Generation 4", "Generation 5","Generation 6")
names(Generation.labs) <- c("2","3","4","5","6")

#--------------------------------------------------------------------------------------------#

# Supplementary Figure 7:  Total Recruit Area by Conceptacle Abundance and Generation
(Recruit_Number_Generation <- ggplot(AllData, aes(x=Conceptacle.Number, y=Recruit.Area, colour = Treatment)) +
  geom_point()+
  facet_wrap(~Generation, labeller = labeller(Generation = Generation.labs), scales = "free", as.table=T)+
  geom_smooth(method=glm, se=FALSE)+
  ylab(bquote("Total Recruit Area ("~cm^2*")")) +
  xlab(bquote('Conceptacle Abundance (' ~ conceptacles~ mm^2*')')) +
  labs(colour="Mean Treatment pH")+
  ylim(-50,400)+
  xlim(0,6)+
  theme_classic()+
  theme(legend.position = c(0.81, 0.43), legend.title = element_text(size=11), 
        legend.text = element_text(size=11), axis.text=element_text(size=11),
        panel.grid.major = element_blank(), axis.title=element_text (size=12),
        panel.grid.minor = element_blank()) + 
  theme(strip.background = element_rect(fill="white", size=0, color="black")) +
  theme(strip.text = element_text(face="bold",size=11, color="black"))+
  theme(panel.spacing = unit(1.4, "lines"))+
  scale_colour_manual(values=c('#999999','#E69F00'),
                      labels=c("Present Day", "Ocean Acidification")))

#--------------------------------------------------------------------------------------------#

# Manuscript Figure 6:  Conceptacle Abundance by Growth and Generation
(Number_Growth_Generation <- ggplot(AllData, aes(x=Growth, y=Conceptacle.Number, colour = Treatment)) +
  geom_point()+
  facet_wrap(~Generation, labeller = labeller(Generation = Generation.labs), scales = "free", as.table=T)+
  geom_smooth(method=glm, se=FALSE)+
  xlab(bquote('Growth (' ~ mm^2~ d^-1*')')) +
  ylab(bquote('Conceptacle Abundance (' ~ conceptacles~ mm^2*')')) +
  labs(colour="Mean Treatment pH")+
  ylim(-0.5,6)+
  xlim(0,25)+
  theme_classic()+
  theme(axis.title.x = element_text(face="bold"))+
  theme(legend.position = c(0.81, 0.43), legend.title = element_text(size=11), 
        legend.text = element_text(size=11), axis.text=element_text(size=11,),
        panel.grid.major = element_blank(), axis.title=element_text (size=12),
        panel.grid.minor = element_blank()) + 
  theme(strip.background = element_rect(fill="white", size=0, color="black")) +
  theme(strip.text = element_text(face="bold",size=11, color="black"))+
  theme(panel.spacing = unit(1.4, "lines"))+
  scale_colour_manual(values=c('#999999','#E69F00'),
                      labels=c("Present Day", "Ocean Acidification")))


####################################################################################################

### Models ###

--------------------------------------------------------------------------------------------

# Conceptacle Abundance by Mean Treatment pH and Generation  
Pmodel.1 = glm(Conceptacle.Number ~ Growth * Treatment * Generation, family = quasipoisson,
               data = AllData)
summary(Pmodel.1)
summary(aov(Pmodel.1))

--------------------------------------------------------------------------------------------

# Total Recruit Area by Conceptacle Abundance and Generation  
Pmodel.2 = glm(Recruit.Area ~ Conceptacle.Number * Treatment * Generation, family = quasipoisson,
               data = AllData)
summary(Pmodel.2)
summary(aov(Pmodel.2))


