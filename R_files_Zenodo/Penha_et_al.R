## This script is from the following manuscript
##Haemosporidian parasites and incubation period influence plumage coloration in tanagers (Passeriformes: Thraupidae)

## Victor Aguiar de Souza Penha1* (https://orcid.org/0000-0002-9036-3862), Fabricius Maia Chaves Bicalho Domingos2 
## (https://orcid.org/0000-0003-2069-9317), Alan Fecchio3 (https://orcid.org/0000-0002-7319-0234), Jeffrey A. Bell4 
## (https://orcid.org/0000-0003-2069-9317), Jason D. Weckstein5 (https://orcid.org/0000-0001-7941-5724), Robert E. Ricklefs6 
## (https://orcid.org/0000-0001-7649-8800), Erika Martins Braga7 (https://orcid.org/0000-0001-5550-7157), Patrícia de Abreu Moreira8 
## (https://orcid.org/0000-0002-6020-449X), Leticia Soares9 (https://orcid.org/0000-0002-6933-8048), Steven Latta10, 
## Graziela Tolesano-Pascoli11 (https://orcid.org/0000-0001-8219-191X), Renata Duarte Alquezar12 (https://orcid.org/0000-0001-8294-722X), 
## Kleber Del-Claro13 (https://orcid.org/0000-0001-8886-9568), Lilian Tonelli Manica2 (https://orcid.org/0000-0001-6005-7103)

##1Graduate program in Ecology and Conservation, Federal University of Paraná, Curitiba, Paraná, Brazil. 
##2Zoology Department, Federal University of Paraná, Curitiba, Paraná, Brazil.
##3Centro de Investigación Esquel de Montaña y Estepa Patagónica (CIEMEP), CONICET – Universidad Nacional de la Patagonia San Juan Bosco, Esquel, Chubut, Argentina 4Department of Biology, University of North Dakota, Grand Forks, United States.  
##5Academy of Natural Sciences of Drexel University and Department of Biodiversity, Earth, and Environmental Science, Drexel University, Philadelphia, Pennsylvania, United States.
##6Department of Biology, University of Missouri-Saint Louis, Saint Louis, Missouri, United States. 
##7Malaria Laboratory, Federal University of Minas Gerais, Belo Horizonte, Minas Gerais, Brazil.
##8Federal University of Ouro Preto, Ouro Preto, Minas Gerais, Brazil. 
##9Department of Biology, Western Ontario University, London, Ontario, Canada. 
##10Conservation and Field Research, National Aviary, Pittsburgh, PA, United States.
##11Zoology Department, Institute of Biological Sciences, University of Brasilia, Brasilia, Distrito Federal, Brazil.
##12Animal Behavior Laboratory, Graduate Program in Ecology, University of Brasilia, Brasilia, Distrito Federal, Brazil.
##13Behavioral Ecology and Interactions Laboratory, Graduate Program in Ecology and Conservation of Natural Resources, Federal University of Uberlândia, Uberlândia, Minas Gerais, Brazil. 
##*Corresponding author: victoraspenha@gmail.com 
##12Animal Behavior Laboratory, Graduate Program in Ecology, Universidade de Brasília, Brasília, Distrito Federal, Brazil.
##13Behavioral Ecology and Interactions Laboratory, Graduate Program in Ecology and Conservation of Natural Resources, Universidade Federal de Uberlândia, Uberlândia, Minas Gerais, Brazil. 
##*Author for correspondence: victoraspenha@gmail.com 


## Remove leftovers
rm(list=ls())

##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################

## Required packages
library(arm)
library(corrplot)
require(usdm)
require(tidyverse)
require(rcompanion)
library(lsr)
library(ggpubr)
library(ggplot2)
library(MuMIn)
library(mnormt)
require(ape)
require(geiger)
require(caper)
require(phytools)
require(adephylo)
library(phylobase)
library(picante)
library(phylosignal)
library(devtools)
library(MASS)
library(vegan)
require(nlme)
require(mvMORPH)
require(OUwie)
require(tidyr)
require(plyr)
library(lme4)
library(regclass)
library(car)
library(caper)
library(AICcmodavg)
library(ape)
library(raster)
library(sp)
library(phangorn)
library(psych)
library(mnormt)
require(ape)
require(geiger)
require(caper)
require(phytools)
require(adephylo)
library(phylobase)
library(picante)
library(phylosignal)
library(devtools)
library(MASS)
library(vegan)
require(nlme)
require(mvMORPH)
require(OUwie)
require(tidyr)
require(plyr)
library(dplyr)
library(mgcv)
library(BBmisc)
library(heatmaply)
library(geiger)

## Data
## change directory
fulldata <- read.csv("/data.csv")
dim(fulldata)
fulldata <- fulldata[order(fulldata$Species),]

## Shults data: you can get it in here 
## https://onlinelibrary.wiley.com/doi/epdf/10.1111/evo.13196 
## change directory
color <- read.csv("/evo13196-sup-0003-table1.csv")
color$Species <- gsub(" ","_",color$Species)
color <- color[which(color$Species %in% fulldata$Species),]
color <- color[order(color$Species),]


## Getting data from their database
fulldata$Species == color$Species
fulldata$AvgSexDich3 <- color$AvgSexDich
fulldata$MalePCA1 <- color$MalePhyloPC1
fulldata$FemalePCA1 <- color$FemalePhyloPC1

### Checking variables
colnames(fulldata)
fulldata$AvgSexDich3 <- as.numeric(fulldata$AvgSexDich3)
fulldata$MalePCA1 <- as.numeric(fulldata$MalePCA1)
fulldata$FemalePCA1 <- as.numeric(fulldata$FemalePCA1)
fulldata$LinRich2 <- fulldata$LinRich / fulldata$Captures

## Multicolinearity models
mod <- lm(log(AvgSexDich3)~sqrt(Prevalence)+sqrt(LinRich2)+Clutch_Size+scale(Incubation)+
            sqrt(Body_Size), fulldata)
VIF(mod)
## Nothing

## Phylogeny: https://www.sciencedirect.com/science/article/abs/pii/S1055790314000578
## change directory
host <- read.nexus("/TanagerMCCTreeTranslated.nex")
recorte <- which(!host$tip.label %in% fulldata$Species)
host2 <- drop.tip(host, recorte)
host2$tip.label


## Models
model_BM <- gls(normalize(AvgSexDich3)~scale(Prevalence)+sqrt(LinRich2)+as.factor(Clutch_Size)+scale(Incubation)+
                  scale(Body_Size), fulldata,
                method = "ML", correlation = corBrownian(phy = host2, form = ~Species))
model_OU <- gls(normalize(AvgSexDich3)~scale(Prevalence)+sqrt(LinRich2)+as.factor(Clutch_Size)+scale(Incubation)+
                  scale(Body_Size), fulldata,
                method = "ML", correlation = corPagel(1,phy = host2, form = ~Species,fixed = F))
AIC(model_BM, model_OU)

## Model selection
selecao <- dredge(model_OU)
## change directory
write.csv(selecao, "/Model1.csv")
selecao2 <- model.avg(selecao)
summary(selecao2)
confint(selecao2)

## check model's fit
binnedplot(predict(model_OU), resid(model_OU), xlab="Expected values", ylab="Average residual", 
           main="Binned residual plot", col.pts="blue", col.int="red")
plot(model_OU)
## Plotting models
## change directory
jpeg("/Dichromatism.jpg", width = 8, height = 4, units = 'in', res = 800)
a <- ggplot(data = fulldata,aes(x = Prevalence, y= normalize(AvgSexDich3)))+
  geom_blank()+geom_point()+
  xlab("Haemosporidian parasite prevalence")+ 
  ylab("Dichromatism")+ theme_bw() + 
  geom_smooth(method = "glm", color = "black", se = T)
b <- ggplot(data = fulldata,aes(x = Incubation, y= normalize(AvgSexDich3)))+
  geom_blank()+geom_point()+
  xlab("Incubation period (number of days)")+ theme_bw() + 
  ylab(NULL)+ theme_bw(base_size = 12, base_family = "Times New Roman")+
  geom_smooth(method = "glm", color = "black", se = T)
ggarrange(a,b, ncol=2)
dev.off()

## Plotting non-significant values
## change directory
jpeg("/Dichromatism3.jpg", width = 4, height = 4, units = 'in', res = 800)
ggplot(data = fulldata,aes(x = sqrt(LinRich2), y= normalize(AvgSexDich3)))+
  geom_blank()+geom_point()+
  xlab("Haemosporidian parasite lineage richness")+ 
  ylab("Dichromatism")+ theme_bw()
dev.off()

## change directory
jpeg("/Dichromatism4.jpg", width = 4, height = 4, units = 'in', res = 800)
ggplot(fulldata, aes(x = as.factor(Clutch_Size), y = normalize(AvgSexDich3))) + geom_boxplot() + 
  geom_jitter(position = position_jitter(0.2)) + theme_bw() + 
  xlab("Clutch size (number of eggs)")+ 
  ylab("Dichromatism")+ theme_bw()
dev.off()

## change directory
jpeg("/Dichromatism5.jpg", width = 4, height = 4, units = 'in', res = 800)
ggplot(data = fulldata,aes(x = scale(Body_Size), y= normalize(AvgSexDich3)))+
  geom_blank()+geom_point()+
  xlab("Body size")+ 
  ylab("Dichromatism")+ theme_bw()
dev.off()

##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################
## Males 
## Models
model_BM_males <- gls(normalize(MalePCA1)~scale(Prevalence)+sqrt(LinRich2)+as.factor(Clutch_Size)+scale(Incubation)+
                        scale(Body_Size), fulldata,
                method = "ML", correlation = corBrownian(phy = host2, form = ~Species))
model_OU_males<- gls(normalize(MalePCA1)~scale(Prevalence)+sqrt(LinRich2)+factor(Clutch_Size)+scale(Incubation)+
                       scale(Body_Size), fulldata,
                method = "ML", correlation = corPagel(1,phy = host2, form = ~Species,fixed = F))
AIC(model_BM_males, model_OU_males)

## Model selection
selecao_male <- dredge(model_OU_males)
## change directory
write.csv(selecao_male, "/Model2.csv")
selecao2_male <- model.avg(selecao_male)
summary(selecao2_male)
confint(selecao2_male)

## Model
binnedplot(predict(model_OU_males), resid(model_OU_males), xlab="Expected values", ylab="Average residual", 
           main="Binned residual plot", col.pts="blue", col.int="red")

## Ploting non-sginificant results
## change directory
jpeg("/Male.jpg", width = 4, height = 4, units = 'in', res = 800)
ggplot(data = fulldata,aes(x = scale(Prevalence), y= normalize(MalePC1)))+
  geom_blank()+geom_point()+
  xlab("Haemosporidian parasite prevalence")+ 
  ylab("Male plumage coloration complexity")+ theme_bw()
dev.off()

## change directory
jpeg("/Male2.jpg", width = 4, height = 4, units = 'in', res = 800)
ggplot(data = fulldata,aes(x = sqrt(LinRich2), y= normalize(MalePC1)))+
  geom_blank()+geom_point()+
  xlab("Haemosporidian parasite lineage richness")+ 
  ylab("Male plumage coloration complexity")+ theme_bw()
dev.off()

## change directory
jpeg("/Male3.jpg", width = 4, height = 4, units = 'in', res = 800)
ggplot(data = fulldata,aes(x = Incubation, y= normalize(MalePC1)))+
  geom_blank()+geom_point()+
  xlab("Incubation period (number of days)")+ theme_bw() + 
  ylab("Male plumage coloration complexity")+ theme_bw(base_size = 12, base_family = "Times New Roman")
dev.off()

## change directory
jpeg("/Male4.jpg", width = 4, height = 4, units = 'in', res = 800)
ggplot(fulldata, aes(x = as.factor(Clutch_Size), y = normalize(MalePC1))) + geom_boxplot() + 
  geom_jitter(position = position_jitter(0.2)) + theme_bw() + 
  xlab("Clutch size (number of eggs)")+ 
  ylab("Male plumage coloration complexity")+ theme_bw()
dev.off()

## change directory
jpeg("/Male5.jpg", width = 4, height = 4, units = 'in', res = 800)
ggplot(data = fulldata,aes(x = scale(Body_Size), y= normalize(MalePC1)))+
  geom_blank()+geom_point()+
  xlab("Body size")+ 
  ylab("Male plumage coloration complexity")+ theme_bw()
dev.off()

## Females 
## Models
model_BM_females <- gls(normalize(FemalePCA1)~scale(Prevalence)+sqrt(LinRich2)+factor(Clutch_Size)+scale(Incubation)+
                          scale(Body_Size), fulldata,
                      method = "ML", correlation = corBrownian(phy = host2, form = ~Species))
model_OU_females<- gls(normalize(FemalePCA1)~scale(Prevalence)+sqrt(LinRich2)+factor(Clutch_Size)+scale(Incubation)+
                         scale(Body_Size), fulldata,
                     method = "ML", correlation = corPagel(1,phy = host2, form = ~Species,fixed = F))
AIC(model_BM_females, model_OU_females)

## Model selection
selecao_female <- dredge(model_BM_females)
## change directory
write.csv(selecao_female, "/Model3.csv")
selecao2_female <- model.avg(selecao_female)
summary(selecao2_female)
confint(selecao2_female)

## Mode
binnedplot(predict(model_BM_females), resid(model_BM_females), xlab="Expected values", ylab="Average residual", 
           main="Binned residual plot", col.pts="blue", col.int="red")

## Plotting models
## change directory
jpeg("/Female.jpg", width = 4, height = 4, units = 'in', res = 800)
ggplot(data = fulldata,aes(x = Incubation, y= FemalePCA1))+
  geom_blank()+geom_point()+
  xlab("Incubation period (number of days)")+ 
  ylab("Female Plumage complexity")+ theme_bw() +
  geom_smooth(method = "glm", color = "black", se = T)
dev.off()

## Plotting non-significant results
## change directory
jpeg("/Female.jpg", width = 4, height = 4, units = 'in', res = 800)
ggplot(data = fulldata,aes(x = scale(Prevalence), y= normalize(FemalePCA1)))+
  geom_blank()+geom_point()+
  xlab("Haemosporidian parasite prevalence")+ 
  ylab("Female plumage coloration complexity")+ theme_bw()
dev.off()

## change directory
jpeg("/Female2.jpg", width = 4, height = 4, units = 'in', res = 800)
ggplot(data = fulldata,aes(x = sqrt(LinRich2), y= normalize(FemalePCA1)))+
  geom_blank()+geom_point()+
  xlab("Haemosporidian parasite lineage richness")+ 
  ylab("Female plumage coloration complexity")+ theme_bw()
dev.off()

## change directory
jpeg("/Female3.jpg", width = 4, height = 4, units = 'in', res = 800)
ggplot(fulldata, aes(x = as.factor(Clutch_Size), y = normalize(FemalePCA1))) + geom_boxplot() + 
  geom_jitter(position = position_jitter(0.2)) + theme_bw() + 
  xlab("Clutch size (number of eggs)")+ 
  ylab("Female plumage coloration complexity")+ theme_bw()
dev.off()

## change directory
jpeg("/Female4.jpg", width = 4, height = 4, units = 'in', res = 800)
ggplot(data = fulldata,aes(x = scale(Body_Size), y= normalize(FemalePCA1)))+
  geom_blank()+geom_point()+
  xlab("Body size")+ 
  ylab("Female plumage coloration complexity")+ theme_bw()
dev.off()
## End of code