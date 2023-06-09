#AA figure
#Fraction of essential amino acids (EAA) and non-essential amino acids (NEAA) that come from dietary
#protein based on equations developed by Hobbie (2017, Rapid Communications in Mass Spectrometry, 31:639-648)

#Ryan Stephens - Dec. 8, 2021


########################################################################################
#Generate data
########################################################################################
#Create sequence of protein content in diet that ranges from 0.05 to 0.9 with increments of 0.01
Protein<-seq(0.05, 0.9, by = 0.01)
Protein_data<-as.data.frame(Protein)#Make into dataframe
head(Protein_data)

library(dplyr)
#Collagen EAA
Collagen<-mutate(Protein_data, Type = "Collagen")#Make column for tissue type
Collagen_AA<-mutate(Collagen, AA = "EAA")#Make column for amino acid type
Collagen_AA_Percent<-mutate(Collagen_AA, Percent = 0.215)#constant proportion of EAA from diet
head(Collagen_AA_Percent)

#Collagen NEAA from protein in diet
Collagen<-mutate(Protein_data, Type = "Collagen")#Make column for tissue type
Collagen_NEAA<-mutate(Collagen, AA = "NEAA")#Make column for amino acid type
#Equation from Hobbie (2017, Rapid Communications in Mass Spectrometry, 31:639-648)
#Fraction of non-essential in Collagen from dietary protein= 0.785 + 0.160 × log(fraction of protein in diet)
Collagen_NEAA_Percent<-mutate(Collagen_NEAA, Percent = 0.785 + 0.160 * log(Protein))
head(Collagen_NEAA_Percent)

#Muscle EAA
Muscle<-mutate(Protein_data, Type = "Muscle")#Make column for tissue type
Muscle_AA<-mutate(Muscle, AA = "EAA")#Make column for amino acid type
Muscle_AA_Percent<-mutate(Muscle_AA, Percent = 0.39)#constant proportion of EAA from diet
head(Muscle_AA_Percent)

#Muscle NEAA from protein in diet
Muscle<-mutate(Protein_data, Type = "Muscle")#Make column for tissue type
Muscle_NEAA<-mutate(Muscle, AA = "NEAA")#Make column for amino acid type
#Equation from Hobbie (2017, Rapid Communications in Mass Spectrometry, 31:639-648)
#Fraction of non-essential in muscle from dietary protein= 0.61 + 0.082 × log(fraction of protein in diet)
Muscle_NEAA_Percent<-mutate(Muscle_NEAA, Percent = 0.61 + 0.082 * log(Protein))
head(Muscle_NEAA_Percent)


AA_data<-rbind(Collagen_AA_Percent, Collagen_NEAA_Percent, Muscle_AA_Percent, Muscle_NEAA_Percent)
AA_data$AA <- factor(AA_data$AA, levels = c("NEAA","EAA"))#change order of consumers
head(AA_data)
str(AA_data)
########################################################################################




########################################################################################
#AA figure
########################################################################################
#Make labels for figure
dat_text <- data.frame(
  label = c("bolditalic('de novo')~bold(NEAA)","bolditalic('de novo')~bold(NEAA)",
            "bold(NEAA~from~protein)","bold(NEAA~from~protein)",
            "bold(EAA~from~protein)","bold(EAA~from~protein)"),
  Type   = c("Collagen", "Muscle",
             "Collagen", "Muscle",
             "Collagen", "Muscle"),
  AA = c("EAA", "EAA","EAA", "EAA","EAA", "EAA"),
  x     = c(.07, .07, .07, .07, .07, .07),
  y     = c(.92, .96,
            .38, .58,
            .11, .19)
)

library(ggplot2)
library(lemon)

ggplot(AA_data, aes(x=Protein, y=Percent, fill=AA)) + 
#Fill area (note that de novo NEAA is the background)
  geom_area( size=.5, colour="whitesmoke")+
  scale_fill_manual(values =c("deeppink4","dodgerblue3"))+#custom color for fill
  theme(panel.background = element_rect(fill = "orange"))+#background color
#Facets and background
  facet_rep_grid(~ Type ) + coord_capped_cart(bottom='both', left='both')+#use package "lemon" to keep tick marks
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1), axis.line=element_line(size=.3))+#remove lines of graph and background color
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+#removes gridlines, but keeps border
  theme(strip.background = element_blank())+#remove facet background
  theme(strip.text.x = element_text(size=12))+#text size of facets 
#Text labels
  geom_text(data = dat_text, mapping = aes(x = x, y = y, label = label),
            parse = T, size = 3.5, hjust = 0, alpha = 0.65, fontface = "bold", color = "white")+
#Axes
  scale_x_continuous(expand = c(0, 0),limits = c(.05, .9),breaks=c(0.05, .25,.5,.75))+#custom scale for x-axis
  scale_y_continuous(expand = c(0, 0),limits = c(0, 1),breaks=c(0,.2,.4,.6,.8,1))+#custom scale for y-axis
  theme(legend.position="none")+#remove legend
  xlab("Protein fraction in diet")+#x axis label
  ylab("Dietary source contribution")+#y axis label
  theme(axis.title.x= element_text(size=12), axis.ticks.length.x=unit(-1, "mm"),
        axis.text.x  = element_text(angle=0, hjust = .5,size=9, colour="black",
                                    margin = unit(c(t = 3, r = 0, b = 0, l = 0), "mm")))+
  theme(axis.title.y = element_text(size=12), axis.ticks.length.y=unit(-1, "mm"),
        axis.text.y  = element_text(angle=0, vjust=0.5, size=9, colour="black",
                                    margin = unit(c(t = 0, r = 3, b = 0, l = 0), "mm")))

ggsave("AA_Figure.jpeg",#save figure
       plot = last_plot(), width =6, height = 3, units = "in",
       dpi = 600) 
########################################################################################