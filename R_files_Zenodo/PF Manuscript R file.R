#####Generalised Linear Models for Parrotfish Community Manuscript#####

####Model for Parrotfish Density####
library(glmmTMB)
library(lme4)
library(dplyr)
library(ggplot2)
setwd("C:/Users/wenze/OneDrive/Desktop/PF Data Dryad")
Data=read.csv("PF Size.csv")
Data$Size.Class <- factor(Data$Size.Class,levels=c("Small", "Medium", "Large"))
g0 =glmmTMB(Total.Abundance~Depth.Class*Aspect*Size.Class +Island +(1|Island:Location), data=Data, family="nbinom2")
summary(g0)
confint(g0)

####Model for Parrotfish Biomass####
Data=read.csv("PF Biomass Transect Data.csv")
Data$New.Biomass=Data$New.Biomass+100
library(lme4)
g4=lmer(log(New.Biomass)~Aspect*Depth.Class +Island +(1|Island:Location), data=Data)
summary(g4)
confint(g4)
qqnorm(residuals(g4))
qqline(residuals(g4))


####Model for Parrotfish Bioerosion####
Data=read.csv("PF Biomass Transect Data.csv")
#Adding a small constant#
Data$Total.Bioerosion=Data$Total.Bioerosion+0.001
g4=lmer(log(Total.Bioerosion)~Aspect*Depth.Class +Island +(1|Island:Location), data=Data)
confint(g4)
summary(g4)


####Model for Parrotfish Species Richness####
Data =read.csv("Species Richness Transect Data.csv")
g1 =glmmTMB(No.of.Species~Depth.Class*Aspect +Island+ (1|Island:Location), data=Data, family="nbinom2")
g1a <- update(g1, control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
?glmmTMBControl
summary(g1a)
confint(g1a)
qqnorm(residuals(g1a))
qqline(residuals(g1))




####COde for Graph Construction####
##Graph for Parrotfish Density##
library(ggplot2)
library(dplyr)
se <- function(x) {sd(x) / sqrt(length(x-1))}
Data =read.csv("PF Abundance Biomass Data.xlsx.csv")
Data$Exposure <- paste(Data$Aspect,Data$Depth.Class)
Data=Data %>% group_by(Island,Aspect,Location,Depth.Class,Exposure) %>% summarise("Abundance"=mean(Abundance))
Data=Data %>% group_by(Exposure,Depth.Class,Aspect) %>% summarise("S.E"=se(Abundance), "nabundance"=mean(Abundance))
Data$Exposure <- factor(Data$Exposure,levels = c("East Deep","East Shallow", "West Deep", "West Shallow"))

newpalette=c("dodgerblue3","firebrick1")
newgreypal=c("grey60","grey30")
a = ggplot(data=Data, aes(Exposure,nabundance,fill=Depth.Class)) + geom_bar(stat="identity", position=position_dodge()) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs( x="", y="Density per 250m2") + geom_errorbar(aes(ymin=nabundance-S.E, ymax=nabundance+S.E), width=.2, position=position_dodge(.9))+scale_fill_manual(values=newgreypal)+theme_bw() +theme(axis.text=element_text(size=18),axis.title=element_text(size=18))+ theme(legend.position = "none")+ggtitle("")+ xlab("") + ylab("")
a

##Graph for Parrotfish Biomass##
setwd("C:/Users/wenze/OneDrive/Desktop/PF Data Dryad")
Data =read.csv("PF Abundance Biomass Data.xlsx.csv")
#Getting biomass in kg
Data$Biomass = Data$Biomass/1000
Data$Exposure <- paste(Data$Aspect,Data$Depth.Class)
Data=Data %>% group_by(Island,Aspect,Location,Depth.Class,Exposure) %>% summarise("Abundance"=mean(Biomass))
Data=Data %>% group_by(Exposure,Depth.Class,Aspect) %>% summarise("S.E"=se(Abundance), "nabundance"=mean(Abundance))
Data$Exposure <- factor(Data$Exposure,levels = c("East Deep", "East Shallow","West Deep", "West Shallow"))

newpalette=c("dodgerblue3","firebrick1")


b = ggplot(data=Data, aes(Exposure,nabundance,fill=Depth.Class)) + geom_bar(stat="identity", position=position_dodge()) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(title="Parrotfish Biomass across Exposure Classes", x="", y="Biomass in kg per 250m2") + geom_errorbar(aes(ymin=nabundance-S.E, ymax=nabundance+S.E), width=.2, position=position_dodge(.9)) +scale_fill_manual(values=newgreypal)+theme_bw()  +theme(axis.text=element_text(size=18),axis.title=element_text(size=18))+ theme(legend.position = "none") +ggtitle("")+ xlab("") + ylab("")
b

##Graph for Parrotfish Bioerosion##
Data=read.csv("PF Biomass Transect Data.csv")
Data=Data %>% group_by(Island,Aspect,Depth.Class,Location) %>% summarise("Bio"=mean(Total.Bioerosion))
Data$Depth.Class = ifelse(Data$Depth.Class=="Intermediate","Deep","Shallow")
Data$Exposure <- paste(Data$Aspect,Data$Depth.Class)
Data=Data %>% group_by(Exposure,Depth.Class,Aspect) %>% summarise("S.E"=se(Bio), "nabundance"=mean(Bio))
Data$Exposure <- factor(Data$Exposure,levels = c("East Deep","East Shallow", "West Deep", "West Shallow"))

d = ggplot(data=Data, aes(Exposure,nabundance,fill=Depth.Class)) + geom_bar(stat="identity", position=position_dodge()) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(title="Parrotfish Bioerosion across Exposure Classes", x="Exposure", y="Bioerosion in kg m-2 yr-1") + geom_errorbar(aes(ymin=nabundance-S.E, ymax=nabundance+S.E), width=.2, position=position_dodge(.9)) +scale_fill_manual(values=newgreypal)+theme_bw()  +theme(axis.text=element_text(size=18),axis.title=element_text(size=18)) + theme(legend.position = "none")+ggtitle("")+ xlab("") + ylab("")
d


##Species Richness##
setwd("C:/Users/wenze/OneDrive/Desktop/PF Data Dryad")
Data =read.csv("Raw fish transect data.csv")
actual =subset(Data, Phase!="")
newset =actual %>% group_by(Island,Aspect,Depth.Class,Location) %>% summarise("No.of.Species"=n_distinct(Species.Name))
Data=newset %>% group_by(Depth.Class,Aspect) %>% summarise("S.E"=se(No.of.Species), "nabundance"=mean(No.of.Species))
Data$Exposure <- paste(Data$Aspect,Data$Depth.Class)
Data$Exposure <- factor(Data$Exposure,levels = c("East Deep","East Shallow", "West Deep",  "West Shallow"))

c = ggplot(data=Data, aes(Exposure,nabundance,fill=Depth.Class)) + geom_bar(stat="identity", position=position_dodge()) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(title="Parrotfish Species Richness across Exposure Classes", x="Exposure", y="Species richness per 250m2") + geom_errorbar(aes(ymin=nabundance-S.E, ymax=nabundance+S.E), width=.2, position=position_dodge(.9)) +scale_fill_manual(values=newgreypal)+theme_bw() +theme(axis.text=element_text(size=18),axis.title=element_text(size=18))+ theme(legend.position = "none")+ggtitle("")+ xlab("") + ylab("")
c

##Variation in Parrotfish Density as a function of Exposure and mediated by body shape and size##
Data =read.csv("Size and shape data.csv")
Data=Data[Data$Abundance>0,]
Data=Data[Data$Species.Name!="Chlorurus capistratoides",]
Data=Data[Data$Species.Name!="Chlorurus enneacanthus",]
Data=Data[Data$Species.Name!="Cetoscarus ocellatus",]
Data=Data[Data$Species.Name!="Scarus flavipectoralis",]
Data=Data[Data$Species.Name!="Scarus tricolor",]
Data$Species.Name <- factor(Data$Species.Name,levels = c("Scarus scaber","Chlorurus sordidus","Scarus rubroviolaceus","Scarus psittacus", "Scarus frenatus","Hipposcarus harid","Scarus viridifucatus","Scarus caudofasciatus","Scarus ghobban","Scarus russelli", "Scarus niger", "Scarus prasiognathos", "Chlorurus strongylocephalus"))
Data$Disc.Size.Class <- factor(Data$Disc.Size.Class,levels =c("Small","Large"))
Data$Exp <- factor(Data$Exp,levels =c("Low","Intermediate","High"))

Disc.gr= ggplot(data=Data, aes(Exp, Species.Name, size=Ab.point.within.size)) + geom_point(aes(colour=Species.Name)) + scale_color_grey(start = 0.8, end = 0.2) +facet_wrap(~ Disc.Size.Class)+theme_bw() +theme(axis.text=element_text(size=16),axis.title=element_text(size=16, face = "bold", vjust = 10)) +theme(strip.text.x = element_text(size = 16)) +labs(title="Body size",y ="", x ="Wave exposure")+scale_size_continuous(range = c(6,14)) + theme(plot.title = element_text(size = 16, hjust=0.5, face = "bold"))
Disc.gr


###Alluvial Plot###
library("easyalluvial")
library("ggalluvial")
Data =read.csv("New Alluvial dataset.csv")


Data$Category = factor(Data$Category,levels=c("Density", "Biomass", "Bioerosion"))


library(ggplot2)
ggplot(data = Data,
       aes(x =Category, y = Percentage, alluvium = Species)) +
  geom_alluvium(aes(fill = Species, colour = Species),
                alpha = .75, decreasing = FALSE) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(type = "qual", palette = "Dark2") +
  scale_color_brewer(type = "qual", palette = "Dark2") +
  ggtitle("") +
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16), legend.key.size = unit(0.8, "cm"),legend.key.width = unit(0.8,"cm"), legend.text =element_text(size=14, face="italic")) + theme(axis.title.x = element_blank()) + theme(legend.title = element_text(size = 16)) + theme(axis.text.x = element_text(size=14, face="bold"))

