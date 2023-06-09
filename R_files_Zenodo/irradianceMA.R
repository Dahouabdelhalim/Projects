### SCRIPT FOR METANALYSIS Irradiance  ####
### Helmut Hillebrand 08.03.2022 ######

#Clear all existing data and close graphics
rm(list=ls())
graphics.off()

#open library
library(metafor)
library(ggplot2)
library(plyr)
library(car)
library(psych)
require(gplots)
require(gridExtra)
library(ggExtra)
library(multcomp)
library(ggrepel)
library(cowplot)
library(grid)
library(egg)
library(dplyr)
library(broom)


#set directory
setwd("~/R/irradiance")

metadata <- read.csv("codedata/metalight.csv")
#check for successful upload
names(metadata)
summary(metadata)

# the file aready contains Hedges'd and its sampling variance
# additionally produce LRR
metadata$LRR<-log(metadata$BM.TRT/metadata$BM.CON)
metadata$LRR[metadata$LRR=="-Inf"]<-NA

#check distributions
hist(metadata$d)
hist(metadata$LRR)
plot(LRR~d,metadata)

#create sampling variance for LRR
metadata$var.LRR<-(metadata$BM.TRT.SD**2/(metadata$BM.TRT**2*metadata$BM.TRT.N))+(metadata$BM.CON.SD**2/(metadata$BM.CON**2*metadata$BM.CON.N))
metadata$var.LRR[metadata$var.LRR=="Inf"]<-NA

##########
## LRR ###
##########

# get rid of NA in effect size
metadata<-metadata[!is.na(metadata$LRR),]
# studies have LRR but no sampling variance 
metadata<-metadata[!is.na(metadata$var.LRR),]
hist((1/sqrt(metadata$var.LRR)))

summary(metadata)

# get rid of studies which need to be excluded for reasons stated in $EXCLUDE
metadata<-metadata[is.na(metadata$EXCLUDE),]

#description of variables see meta-information
summary(metadata)

#summarize some categories, delete some variables not needed for the analysis, correct mistakes and typoes

#typoes in characterizing Lab and Field
metadata$LAB.FIELD[metadata$LAB.FIELD=="field"]<-"Field"
metadata$LAB.FIELD[metadata$LAB.FIELD=="lab"]<-"Lab"
metadata$LAB.FIELD<-as.factor(metadata$LAB.FIELD)

#get habitats into simple categories
unique(metadata$HABITAT)
metadata$HABITAT[metadata$HABITAT=="coastal "]<-"coastal"
metadata$HABITAT[metadata$HABITAT=="Lab"]<-"culture"
metadata$HABITAT[metadata$HABITAT=="seawater"]<-"culture"
metadata$HABITAT[metadata$HABITAT=="atlantic"]<-"offshore"
metadata$HABITAT[metadata$HABITAT=="South Pacific"]<-"offshore"
metadata$HABITAT[metadata$HABITAT=="Antarctica"]<-"offshore"
metadata$HABITAT[metadata$HABITAT=="antarctica"]<-"offshore"
metadata$HABITAT[metadata$HABITAT=="ocean"]<-"offshore"
metadata$HABITAT[metadata$HABITAT=="North Sea"]<-"coastal"
metadata$HABITAT[metadata$HABITAT=="estuary"]<-"coastal"
metadata$HABITAT<-as.factor(metadata$HABITAT)

#correct experiment time
unique(metadata$EXP.TIME)
metadata$EXP.TIME[metadata$EXP.TIME=="spring-summer"]<-NA
metadata$EXP.TIME[metadata$EXP.TIME=="annual"]<-NA
metadata$EXP.TIME<-as.character(metadata$EXP.TIME)
metadata$EXP.TIME[is.na(metadata$EXP.TIME)]<-"undefined"
metadata$EXP.TIME<-as.factor(metadata$EXP.TIME)


#create binary categories for system and organism type
unique(metadata$ORGANISM)
metadata$SYSTEM<-"benthos"
metadata$SYSTEM[metadata$ORGANISM=="phytoplankton"]<-"plankton"
metadata$SYSTEM<-as.factor(metadata$SYSTEM)
metadata$ORG.TYPE<-"macrophyte"
metadata$ORG.TYPE[metadata$ORGANISM=="phytoplankton"]<-"microalgae"
metadata$ORG.TYPE[metadata$ORGANISM=="microalgae"]<-"microalgae"
metadata$ORG.TYPE[metadata$ORGANISM=="phytobenthos"]<-"microalgae"
metadata$ORG.TYPE<-as.factor(metadata$ORG.TYPE)

#correct double entries in experiment unit type
unique(metadata$EXP.UNIT.TYPE)
metadata$EXP.UNIT.TYPE[metadata$EXP.UNIT.TYPE=="in situ incubation"]<-"incubation"

#refine NA in other treatments
unique(metadata$OTHER.TRTS)
metadata$OTHER.TRTS[is.na(metadata$OTHER.TRTS)]<-"none"

#get unified categories for light treatment type
unique(metadata$LIGHT.TRT.TYPE)
metadata$LIGHT.TRT.TYPE<-as.character(metadata$LIGHT.TRT.TYPE)
metadata$LIGHT.TRT.TYPE[metadata$LIGHT.TRT.TYPE=="reduced light"] <-"light reduced"
metadata$LIGHT.TRT.TYPE[metadata$LIGHT.TRT.TYPE=="light gradient"] <-"gradient"
metadata$LIGHT.TRT.TYPE[metadata$LIGHT.TRT.TYPE=="shading"] <-"shading screen"
metadata$LIGHT.TRT.TYPE[metadata$LIGHT.TRT.TYPE=="Crassostrea virginica oysters glued to lines"] <-"other"
metadata$LIGHT.TRT.TYPE[metadata$LIGHT.TRT.TYPE=="incubation at different depths"] <-"other"
metadata$LIGHT.TRT.TYPE[metadata$LIGHT.TRT.TYPE=="color_light reduced"] <-"other"
metadata$LIGHT.TRT.TYPE[metadata$LIGHT.TRT.TYPE=="inhibition"]<- "shading screen"
metadata$LIGHT.TRT.TYPE[is.na(metadata$LIGHT.TRT.TYPE)]<-"other"
metadata$LIGHT.TRT.TYPE<-as.factor(metadata$LIGHT.TRT.TYPE)

#correct typo in day length
metadata$DAYLENGTH[metadata$DAYLENGTH=="Natural"]<-"natural"


#create simpler and consistent response categories
unique(metadata$RESP.CAT)
metadata$RESP.CAT[metadata$RESP.CAT=="biom"]<-"biomass"
metadata$RESP.CAT[metadata$RESP.CAT=="PP"& metadata$BM.MEASURE=="Primary Production (µg CL-1d-1)"]<-"area-specific production"
metadata$RESP.CAT[metadata$RESP.CAT=="PP"& metadata$BM.MEASURE=="production at saturating irradiance"]<-"mass-specific production"
metadata$RESP.CAT[metadata$RESP.CAT=="PP"& metadata$BM.MEASURE=="PP"]<-"area-specific production"
metadata$RESP.CAT[metadata$RESP.CAT=="phys"& metadata$BM.MEASURE=="Fv/Fm"]<-"quantum yield"
metadata$RESP.CAT[metadata$RESP.CAT=="phys"& metadata$BM.MEASURE=="Pmax"]<-"mass-specific production"
metadata$RESP.CAT[metadata$RESP.CAT=="phys"& metadata$BM.MEASURE=="Chl a (pg cell-1)"]<-"cellular content"
metadata$RESP.CAT[metadata$RESP.CAT=="productivity"& metadata$BM.MEASURE=="Fv/Fm"]<-"quantum yield"
metadata$RESP.CAT[metadata$RESP.CAT=="productivity"& metadata$BM.MEASURE=="net photosynthesis (cell)"]<-"mass-specific production"
metadata$RESP.CAT[metadata$RESP.CAT=="productivity"& metadata$BM.MEASURE=="net PS"]<-"mass-specific production"
metadata$RESP.CAT[metadata$RESP.CAT=="productivity"& metadata$BM.MEASURE=="Photosynthetic oxygen evolution rate curves"]<-"mass-specific production"
metadata$RESP.CAT[metadata$RESP.CAT=="productivity"& metadata$BM.MEASURE=="Pmax"]<-"mass-specific production"
metadata$RESP.CAT[metadata$RESP.CAT=="productivity"& metadata$BM.MEASURE=="PS capacity"]<-"mass-specific production"
metadata$RESP.CAT[metadata$RESP.CAT=="productivity"& metadata$BM.MEASURE=="photosynthesis"]<-"mass-specific production"
metadata$RESP.CAT[metadata$RESP.CAT=="productivity"& metadata$BM.MEASURE=="Photosynthetic electron transport rate curves"]<-"quantum yield"
metadata$RESP.CAT[metadata$RESP.CAT=="productivity"& metadata$BM.MEASURE=="effective absorption cross section of PS II (?PSII) "]<-"quantum yield"
metadata$RESP.CAT[metadata$RESP.CAT=="productivity"& metadata$BM.MEASURE=="effective absorption cross section of PS II (?PSII)"]<-"quantum yield"
metadata$RESP.CAT[metadata$RESP.CAT=="productivity"& metadata$BM.MEASURE=="maximum electron transport rate"]<-"quantum yield"
metadata$RESP.CAT[metadata$RESP.CAT=="productivity"& metadata$BM.MEASURE=="ETRmax (e? PSII?1 s?1)" ]<-"quantum yield"
metadata$RESP.CAT[metadata$RESP.CAT=="productivity"& metadata$BM.MEASURE=="maximum electron transport rate "]<-"quantum yield"
metadata$RESP.CAT[metadata$RESP.CAT=="productivity"& metadata$BM.MEASURE=="ETRmax (e? PSII?1 s?1) "]<-"quantum yield"
metadata$RESP.CAT[metadata$RESP.CAT=="productivity"& metadata$BM.MEASURE=="rETRmax "]<-"quantum yield"
metadata$RESP.CAT[metadata$RESP.CAT=="productivity"& metadata$BM.MEASURE=="Primary Production (µg CL-1d-1)"]<-"area-specific production"
metadata$RESP.CAT[metadata$RESP.CAT=="productivity"& metadata$BM.MEASURE=="production"]<-"area-specific production"
metadata$RESP.CAT[metadata$RESP.CAT=="productivity"& metadata$BM.MEASURE=="calcification rate"]<-"mass-specific production"

metadata$RESP.CAT<-as.character(metadata$RESP.CAT)
metadata$RESP.CAT[metadata$RESP.CAT=="cellular content"& metadata$BM.MEASURE=="Chl a (pg/cell)"]<-"pigment content"
metadata$RESP.CAT[metadata$RESP.CAT=="cellular content"& metadata$BM.MEASURE=="Chl a (µg/g)"]<-"pigment content"
metadata$RESP.CAT[metadata$RESP.CAT=="cellular content"& metadata$BM.MEASURE=="chl a (mg g FW)"]<-"pigment content"
metadata$RESP.CAT[metadata$RESP.CAT=="cellular content"& metadata$BM.MEASURE=="Chl a (pg cell)"]<-"pigment content"
metadata$RESP.CAT[metadata$RESP.CAT=="cellular content"& metadata$BM.MEASURE=="Chl a (ng per cell)"]<-"pigment content"
metadata$RESP.CAT[metadata$RESP.CAT=="cellular content"& metadata$BM.MEASURE=="Chl a (µg/g)"]<-"pigment content"
metadata$RESP.CAT[metadata$RESP.CAT=="cellular content"& metadata$BM.MEASURE=="cellular Chl a (pg/cell)"]<-"pigment content"
metadata$RESP.CAT[metadata$RESP.CAT=="cellular content"& metadata$BM.MEASURE=="chl a (mg g FW)"]<-"pigment content"
metadata$RESP.CAT[metadata$RESP.CAT=="cellular content"& metadata$BM.MEASURE=="Chl a (pg cell-1) "]<-"pigment content"
metadata$RESP.CAT[metadata$RESP.CAT=="cellular content"& metadata$BM.MEASURE=="Chl a pg cell-1"]<-"pigment content"
metadata$RESP.CAT[metadata$RESP.CAT=="cellular content"& metadata$BM.MEASURE=="Chl a (pg cell-1"]<-"pigment content"
metadata$RESP.CAT[metadata$RESP.CAT=="cellular content"& metadata$BM.MEASURE=="Chla (pg cell)"]<-"pigment content"
metadata$RESP.CAT[metadata$RESP.CAT=="cellular content"& metadata$BM.MEASURE=="chl a (pg cell)"]<-"pigment content"
metadata$RESP.CAT[metadata$RESP.CAT=="cellular content"& metadata$BM.MEASURE=="chl a (mg g dw)"]<-"pigment content"
metadata$RESP.CAT[metadata$RESP.CAT=="cellular content"& metadata$BM.MEASURE=="Cellular cholorphyll a"]<-"pigment content"
metadata$RESP.CAT[metadata$RESP.CAT=="cellular content"& metadata$BM.MEASURE=="Chla  (pg cell)"]<-"pigment content"
metadata$RESP.CAT[metadata$RESP.CAT=="cellular content"& metadata$BM.MEASURE=="Chl a (pg cell-1)"]<-"pigment content"

metadata$RESP.CAT[metadata$RESP.CAT=="cellular content"]<- "cont.carb"
metadata$RESP.CAT[metadata$RESP.CAT=="pigment content"]<- "cont.pigm"
metadata$RESP.CAT[metadata$RESP.CAT=="area-specific production"]<- "prod.abs"
metadata$RESP.CAT[metadata$RESP.CAT=="mass-specific production"]<- "prod.spec"
metadata$RESP.CAT[metadata$RESP.CAT=="quantum yield"]<- "quant.etr"
metadata$RESP.CAT[metadata$RESP.CAT=="growth rate"]<- "growth"
metadata$RESP.CAT<-as.factor(metadata$RESP.CAT)

#get binary categories for response type
metadata$RESPONSE.TYPE[metadata$RESPONSE.TYPE=="abundance"]<-"biom"
metadata$RESPONSE.TYPE[metadata$RESPONSE.TYPE=="dQ"]<-"phys"
metadata$RESPONSE.TYPE[metadata$RESPONSE.TYPE=="dSS"]<-"biom"
metadata$RESPONSE.TYPE[metadata$RESPONSE.TYPE=="phys"&metadata$RESP.CAT=="prod.abs"]<-"biom"

summary(metadata)

#create a confined data frame with the variables needed
metadata<-metadata[,c("WOS","AUTHORS","YEAR","EXP.TIME","LAB.FIELD","HABITAT","LATITUDE","SYSTEM","ORG.TYPE",
                      "DURATION","EXP.UNIT.TYPE", "SIZE.VOL","SIZE.AREA","LIGHT.TRT.TYPE","DAYLENGTH","TEMP",
                      "SRP","DIN","TP","TN","Light.CON.s","Light.CON.d", "RESP.CAT","RESPONSE.TYPE",
                      "Light.REG.TRT","OTHER.TRTS","LRR","var.LRR")]

summary(metadata)


#combine moderators measured in different units

#size
#measured as area (m2) or volume (L), but not both
#assuming that the benthic areas are 1 m in height, we derive volumes
metadata$SIZE<-metadata$SIZE.VOL
metadata$SIZE[!is.na(metadata$SIZE.AREA)]<-metadata$SIZE.AREA[!is.na(metadata$SIZE.AREA)]*1000
hist(log(metadata$SIZE+1))
plot(SIZE~SIZE.AREA,metadata)
metadata$SIZE.VOL<-NULL
metadata$SIZE.AREA<-NULL

#control irradiance in micromol s-1 and mol d-1, neither or
metadata$CON.IRR<-metadata$Light.CON.s
#assuming 12 hours light 
metadata$CON.IRR[!is.na(metadata$Light.CON.d)]<-metadata$Light.CON.d[!is.na(metadata$Light.CON.d)]/12/3600*1000000
hist(log(metadata$CON.IRR))
metadata$Light.CON.s<-NULL
metadata$Light.CON.d<-NULL


summary(metadata)
metadata$DURATION[metadata$DURATION==0]<-NA

# check for log transformation
pairs.panels(metadata[,c(6,9,13,15,24:25)]) 
metadata$DURATION<-log(metadata$DURATION)
metadata$SIZE<-log(metadata$SIZE)
metadata$LATITUDE<-abs(metadata$LATITUDE)
metadata$DIN<-log(metadata$DIN)
metadata$CON.IRR<-log(metadata$CON.IRR)
pairs.panels(metadata[,c(6,9,13,15,24:25)]) 

#get rid of 1 strudy with a var.LRR = 0
metadata<-metadata[metadata$var.LRR>0,]

#### create a size variable for plotting
### radius of points will be proportional to the inverse standard errors
### hence the area of the points will be proportional to inverse variances
hist(log(metadata$var.LRR))
metadata$pointsize <- 1/log(metadata$var.LRR+1.1)
metadata$pointsize <- metadata$pointsize / 5+1
metadata$WOS<-as.factor(metadata$WOS)
summary(metadata)

#address the multiple ES per study problem
testsum <-
  metadata %>%
  dplyr::group_by(WOS, AUTHORS,YEAR,EXP.TIME, LAB.FIELD, SYSTEM, ORG.TYPE,RESP.CAT,RESPONSE.TYPE,HABITAT) %>%
  dplyr::summarize(N= length(LRR),
                   timesteps = length(unique(DURATION)),
                   levels = length(unique(Light.REG.TRT)),
                   sites = length(unique(as.character(LATITUDE))),
                   treatm = length(unique(as.character(OTHER.TRTS))),
                   mean = mean(LRR, na.rm=T),
                   sd = sd(LRR,na.rm=T), 
                   cv=abs(sd/mean),
                   contexts=levels*sites*treatm*timesteps)

plot(timesteps~N, testsum)
plot(treatm~N, testsum)
plot(levels~N, testsum)
plot(sites~N, testsum)
plot(contexts~N, testsum)

summary(testsum)
hist(testsum$N)


#### plotting data with single predictors
names(metadata)

#create a graphical variable to separate shapes of symbols for lab.field and org.type
metadata$USI<-with(metadata, interaction(LAB.FIELD,ORG.TYPE,drop=TRUE))


dur.LRR<-ggplot(metadata, aes(DURATION, LRR, shape=USI, colour=RESPONSE.TYPE,size=pointsize))+
  geom_point(alpha=.25)+
  ylab("Log response ratio")+
  xlab("Duration (ln days)")+
  theme_bw()+
  scale_color_manual(values = c("darkgreen", "darkorange"))+
  scale_shape_manual(values=c(0,15,1,16),labels=c("F macro", "L macro", "F micro", "L micro")) + 
  #scale_fill_manual(breaks=c("Lab", "Field"), values=c("Lab"="white", "Field"="lightgrey"))+
  geom_hline(aes(yintercept=0), color = "darkgrey")+
  theme(axis.title.y=element_text(size=18, face="plain", colour="black",vjust=0.3),axis.text.y=element_text(size=12,face="bold",colour="black",angle=0,hjust=0.4))+
  theme(axis.title.x=element_text(size=18,face="plain",colour="black",vjust=0),axis.text.x=element_text(size=12,face="bold",colour="black"))+
  theme(axis.ticks=element_line(colour="black",size=1),axis.ticks.length=unit(0.3,"cm"))+
  theme(panel.border=element_rect(colour="black",size=1.5))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm"))+
  theme(legend.justification=c(1,1), legend.position=c(.95,0.95))+
  theme(legend.text = element_text(size=18))+
  theme(legend.title = element_blank())+
  guides(color = "none", size = "none",shape = guide_legend(override.aes = list(size=5)))

dur.LRR

lat.LRR<-ggplot(metadata, aes(LATITUDE, LRR, shape=USI, colour=RESPONSE.TYPE,size=pointsize))+
  geom_point(alpha=.25)+
  ylab("Log response ratio")+
  xlab("Absolute latitude (°N or °S)")+
  theme_bw()+
  scale_color_manual(values = c("darkgreen", "darkorange"))+
  scale_shape_manual(values=c(0,15,1,16),labels=c("F macro", "L macro", "F micro", "L micro")) + 
  geom_hline(aes(yintercept=0), color = "darkgrey")+
  theme(axis.title.y=element_text(size=18, face="plain", colour="black",vjust=0.3),axis.text.y=element_text(size=12,face="bold",colour="black",angle=0,hjust=0.4))+
  theme(axis.title.x=element_text(size=18,face="plain",colour="black",vjust=0),axis.text.x=element_text(size=12,face="bold",colour="black"))+
  theme(axis.ticks=element_line(colour="black",size=1),axis.ticks.length=unit(0.3,"cm"))+
  theme(panel.border=element_rect(colour="black",size=1.5))+
  theme(legend.position="none")+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm"))+
  guides(color = "none", size = "none",shape = guide_legend(override.aes = list(size=5)))

lat.LRR


light.LRR<-ggplot(metadata, aes(Light.REG.TRT, LRR, shape=USI, colour=RESPONSE.TYPE,size=pointsize))+
  geom_point(alpha=.25)+
  ylab("Log response ratio")+
  xlab("Remaining light (% control)")+
  theme_bw()+
  scale_color_manual(values = c("darkgreen", "darkorange"))+
  scale_shape_manual(values=c(0,15,1,16),labels=c("F macro", "L macro", "F micro", "L micro")) + 
  geom_hline(aes(yintercept=0), color = "darkgrey")+
  theme(axis.title.y=element_text(size=18, face="plain", colour="black",vjust=0.3),axis.text.y=element_text(size=12,face="bold",colour="black",angle=0,hjust=0.4))+
  theme(axis.title.x=element_text(size=18,face="plain",colour="black",vjust=0),axis.text.x=element_text(size=12,face="bold",colour="black"))+
  theme(axis.ticks=element_line(colour="black",size=1),axis.ticks.length=unit(0.3,"cm"))+
  theme(panel.border=element_rect(colour="black",size=1.5))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm"))  +
  theme(legend.justification=c(1,1), legend.position=c(.95,0.95))+
  theme(legend.text = element_text(size=18))+
  theme(legend.title = element_blank())+
  guides(color = guide_legend(override.aes = list(size=5,alpha=1)), size = "none",shape = "none")
light.LRR

irr.LRR<-ggplot(metadata, aes(CON.IRR, LRR, shape=USI, colour=RESPONSE.TYPE,size=pointsize))+
  geom_point(alpha=.25)+
  ylab("Log response ratio")+
  xlab("Incoming light (micromole m-2 s-1)")+
  theme_bw()+
  scale_color_manual(values = c("darkgreen", "darkorange"))+
  scale_shape_manual(values=c(0,15,1,16),labels=c("F macro", "L macro", "F micro", "L micro")) + 
  geom_hline(aes(yintercept=0), color = "darkgrey")+
  theme(axis.title.y=element_text(size=18, face="plain", colour="black",vjust=0.3),axis.text.y=element_text(size=12,face="bold",colour="black",angle=0,hjust=0.4))+
  theme(axis.title.x=element_text(size=18,face="plain",colour="black",vjust=0),axis.text.x=element_text(size=12,face="bold",colour="black"))+
  theme(axis.ticks=element_line(colour="black",size=1),axis.ticks.length=unit(0.3,"cm"))+
  theme(panel.border=element_rect(colour="black",size=1.5))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm"))+
  theme(legend.justification=c(1,1), legend.position=c(.25,0.95))+
  theme(legend.text = element_text(size=18))+
  theme(legend.title = element_blank())+
  guides(color = guide_legend(override.aes = list(size=5,alpha=1)), size = "none",shape = "none")


irr.LRR

size.LRR<-ggplot(metadata, aes(SIZE, LRR, shape=USI, colour=RESPONSE.TYPE,size=pointsize))+
  geom_point(alpha=.25)+
  ylab("Log response ratio")+
  xlab("Experiment Size")+
  theme_bw()+
  scale_color_manual(values = c("darkgreen", "darkorange"))+
  scale_shape_manual(values=c(0,15,1,16),labels=c("F macro", "L macro", "F micro", "L micro")) + 
  geom_hline(aes(yintercept=0), color = "darkgrey")+
  theme(axis.title.y=element_text(size=18, face="plain", colour="black",vjust=0.3),axis.text.y=element_text(size=12,face="bold",colour="black",angle=0,hjust=0.4))+
  theme(axis.title.x=element_text(size=18,face="plain",colour="black",vjust=0),axis.text.x=element_text(size=12,face="bold",colour="black"))+
  theme(axis.ticks=element_line(colour="black",size=1),axis.ticks.length=unit(0.3,"cm"))+
  theme(panel.border=element_rect(colour="black",size=1.5))+
  theme(legend.position="none")+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm"))+
  guides(color = "none", size = "none",shape = guide_legend(override.aes = list(size=5)))

size.LRR

temp.LRR<-ggplot(metadata, aes(TEMP, LRR,  shape=USI, colour=RESPONSE.TYPE,size=pointsize))+
  geom_point(alpha=.25)+
  ylab("Log response ratio")+
  xlab("Temperature (°C)")+
  theme_bw()+
  scale_color_manual(values = c("darkgreen", "darkorange"))+
  scale_shape_manual(values=c(0,15,1,16),labels=c("F macro", "L macro", "F micro", "L micro")) + 
  geom_hline(aes(yintercept=0), color = "darkgrey")+
  theme(axis.title.y=element_text(size=18, face="plain", colour="black",vjust=0.3),axis.text.y=element_text(size=12,face="bold",colour="black",angle=0,hjust=0.4))+
  theme(axis.title.x=element_text(size=18,face="plain",colour="black",vjust=0),axis.text.x=element_text(size=12,face="bold",colour="black"))+
  theme(axis.ticks=element_line(colour="black",size=1),axis.ticks.length=unit(0.3,"cm"))+
  theme(panel.border=element_rect(colour="black",size=1.5))+
  theme(legend.position="none")+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm"))+
  guides(color = "none", size = "none",shape = guide_legend(override.aes = list(size=5)))

temp.LRR


year.LRR<-ggplot(metadata, aes(YEAR, LRR,  shape=USI, colour=RESPONSE.TYPE,size=pointsize))+
  geom_point(alpha=.25)+
  ylab("Log response ratio")+
  xlab("Publication year")+
  theme_bw()+
  scale_color_manual(values = c("darkgreen", "darkorange"))+
  scale_shape_manual(values=c(0,15,1,16),labels=c("F macro", "L macro", "F micro", "L micro")) + 
  geom_hline(aes(yintercept=0), color = "darkgrey")+
  theme(axis.title.y=element_text(size=18, face="plain", colour="black",vjust=0.3),axis.text.y=element_text(size=12,face="bold",colour="black",angle=0,hjust=0.4))+
  theme(axis.title.x=element_text(size=18,face="plain",colour="black",vjust=0),axis.text.x=element_text(size=12,face="bold",colour="black"))+
  theme(axis.ticks=element_line(colour="black",size=1),axis.ticks.length=unit(0.3,"cm"))+
  theme(panel.border=element_rect(colour="black",size=1.5))+
  theme(legend.position="none")+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm"))+
  guides(color = "none", size = "none",shape = guide_legend(override.aes = list(size=5)))

year.LRR

TN<-ggplot(metadata, aes(DIN, LRR,  shape=USI, colour=RESPONSE.TYPE,size=pointsize))+
  geom_point(alpha=.25)+
  ylab("Log response ratio")+
  xlab("Total nitrogen (ln µM)")+
  theme_bw()+
  scale_color_manual(values = c("darkgreen", "darkorange"))+
  scale_shape_manual(values=c(0,15,1,16),labels=c("F macro", "L macro", "F micro", "L micro")) + 
  geom_hline(aes(yintercept=0), color = "darkgrey")+
  theme(axis.title.y=element_text(size=18, face="plain", colour="black",vjust=0.3),axis.text.y=element_text(size=12,face="bold",colour="black",angle=0,hjust=0.4))+
  theme(axis.title.x=element_text(size=18,face="plain",colour="black",vjust=0),axis.text.x=element_text(size=12,face="bold",colour="black"))+
  theme(axis.ticks=element_line(colour="black",size=1),axis.ticks.length=unit(0.3,"cm"))+
  theme(panel.border=element_rect(colour="black",size=1.5))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm"))+
  theme(legend.justification=c(1,1), legend.position=c(.95,0.95))+
  theme(legend.text = element_text(size=18))+
  theme(legend.title = element_blank())+
  guides(color = "none", size = "none",shape = guide_legend(override.aes = list(size=5)))


TN



##########################################
####      Multilevel Meta-analyses  ######
##########################################

#in this section we run the main analyses. 
#as we have interdependent effect sizes (multiple from each experiment) we use a 
#multi-level meta-analysis that accounts for this by random effects and creates reliable
#estimated and confidence intervals


summary(metadata)

#to correct for non-independence, I create two unique identifiers
# for papers AUTHORS & YEAR
metadata$UPI<-with(metadata, interaction(WOS,YEAR,drop=TRUE))

# for experiments within papers 
names(metadata)
metadata$UEI<-with(metadata, interaction(EXP.TIME,LAB.FIELD,SYSTEM,RESPONSE.TYPE,ORG.TYPE,HABITAT,RESP.CAT,OTHER.TRTS,drop=TRUE))
unique(metadata$UEI)
#as the character string is not handy, I use numbers
metadata$UEI<-as.numeric(metadata$UEI)

summary(metadata)

#multi-level MA WITHOUT moderators
mod.test<-rma.mv(LRR,var.LRR,random = ~ 1 | UPI/UEI, 
                 method="ML",data=metadata[metadata$var.LRR>0,])
summary(mod.test)

nomoderator<-tidy(summary(mod.test))
nomoderator$k<-mod.test$k
nomoderator$AIC<-mod.test$fit.stats$REML[3]
nomoderator$model<-"no.mod"

#publication bias for model without moderators
ranktest(mod.test)

regtest(LRR,var.LRR,data=metadata)
tes(metadata$LRR,metadata$var.LRR)
trimfill(mod.test)
trimfill(rma(LRR,var.LRR,data=metadata))
fsn(LRR,var.LRR,data=metadata, type="Rosenberg")


#model with a subset of moderators
test.full<-rma.mv(LRR,var.LRR,mods = ~RESPONSE.TYPE+RESP.CAT+Light.REG.TRT+LAB.FIELD+ORG.TYPE+SYSTEM+HABITAT+TEMP+DURATION,
                           random = ~ 1 | UPI/UEI, 
                 method="REML",data=metadata)
summary(test.full)
allmoderator<-tidy(summary(test.full))
allmoderator$k<-test.full$k
allmoderator$AIC<-test.full$fit.stats$REML[3]
allmoderator$model<-"full"

output<-rbind(nomoderator,allmoderator)


#model with all moderators
test.complete<-rma.mv(LRR,var.LRR,mods = ~RESPONSE.TYPE+RESP.CAT+Light.REG.TRT+LAB.FIELD+ORG.TYPE+SYSTEM+HABITAT+TEMP+DURATION+
                        YEAR+EXP.UNIT.TYPE+LIGHT.TRT.TYPE,
                  random = ~ 1 | UPI/UEI, 
                  method="REML",data=metadata)
summary(test.complete)
complete.mod<-tidy(summary(test.complete))
complete.mod$k<-test.complete$k
complete.mod$AIC<-test.complete$fit.stats$REML[3]
complete.mod$model<-"complete"
output<-rbind(output,complete.mod)

#I2
X<-model.matrix(test.complete)
W <- diag(1/test.complete$vi)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * test.complete$sigma2 / (sum(test.complete$sigma2) + (test.complete$k-test.complete$p)/sum(diag(P)))

#publication bias for model without moderators

ranktest(test.complete)
regtest(rma(LRR,var.LRR,mods = ~RESPONSE.TYPE+RESP.CAT+Light.REG.TRT+LAB.FIELD+ORG.TYPE+SYSTEM+HABITAT+TEMP+DURATION+
                            YEAR+EXP.UNIT.TYPE+LIGHT.TRT.TYPE,
                              method="REML",data=metadata))

tes(rma(LRR,var.LRR,mods = ~RESPONSE.TYPE+RESP.CAT+Light.REG.TRT+LAB.FIELD+ORG.TYPE+SYSTEM+HABITAT+TEMP+DURATION+
          YEAR+EXP.UNIT.TYPE+LIGHT.TRT.TYPE,
        method="REML",data=metadata), test="chi2")
(mod.test$sigma2 - test.complete$sigma2) / mod.test$sigma2
(sum(mod.test$sigma2) - sum(test.complete$sigma2)) / sum(mod.test$sigma2)



# models for those moderators that had too many NA, separate for biomass and phasiology as well as complete

TN.biom<-rma.mv(LRR,var.LRR,mods = ~TN,
                          random = ~ 1 | UPI/UEI, 
                          method="REML",data=metadata[metadata$RESPONSE.TYPE=="biom",])

summary(TN.biom)
TN.phys<-rma.mv(LRR,var.LRR,mods = ~TN,
                random = ~ 1 | UPI/UEI, 
                method="REML",data=metadata[metadata$RESPONSE.TYPE=="biom",])

summary(TN.phys)write.csv(output2,file="output2.csv")
TN.all<-rma.mv(LRR,var.LRR,mods = ~TN,
                random = ~ 1 | UPI/UEI, 
                method="REML",data=metadata)

summary(TN.all)


IRR.biom<-rma.mv(LRR,var.LRR,mods = ~CON.IRR,
                random = ~ 1 | UPI/UEI, 
                method="REML",data=metadata[metadata$RESPONSE.TYPE=="biom",])

summary(IRR.biom)
IRR.phys<-rma.mv(LRR,var.LRR,mods = ~CON.IRR,
                random = ~ 1 | UPI/UEI, 
                method="REML",data=metadata[metadata$RESPONSE.TYPE=="biom",])

summary(IRR.phys)write.csv(output2,file="output2.csv")
IRR.all<-rma.mv(LRR,var.LRR,mods = ~CON.IRR,
               random = ~ 1 | UPI/UEI, 
               method="REML",data=metadata)

summary(IRR.all)


LAT.biom<-rma.mv(LRR,var.LRR,mods = ~abs(LATITUDE),
                random = ~ 1 | UPI/UEI, 
                method="REML",data=metadata[metadata$RESPONSE.TYPE=="biom",])

summary(LAT.biom)
LAT.phys<-rma.mv(LRR,var.LRR,mods = ~abs(LATITUDE),
                random = ~ 1 | UPI/UEI, 
                method="REML",data=metadata[metadata$RESPONSE.TYPE=="biom",])

summary(LAT.phys)write.csv(output2,file="output2.csv")
LAT.all<-rma.mv(LRR,var.LRR,mods = ~abs(LATITUDE),
               random = ~ 1 | UPI/UEI, 
               method="REML",data=metadata)

summary(LAT.all)


size.biom<-rma.mv(LRR,var.LRR,mods = ~SIZE,
                random = ~ 1 | UPI/UEI, 
                method="REML",data=metadata[metadata$RESPONSE.TYPE=="biom",])

summary(size.biom)
size.phys<-rma.mv(LRR,var.LRR,mods = ~SIZE,
                random = ~ 1 | UPI/UEI, 
                method="REML",data=metadata[metadata$RESPONSE.TYPE=="biom",])

summary(size.phys)write.csv(output2,file="output2.csv")
size.all<-rma.mv(LRR,var.LRR,mods = ~SIZE,
               random = ~ 1 | UPI/UEI, 
               method="REML",data=metadata)

summary(size.all)


##########################################
#### Getting data for effect sizes  ######
##########################################

#the following loop creates rma.mv estimates for mean effect sizes
#and their CI for each category per variable split for biomass/physiology

USI<-unique(metadata$RESPONSE.TYPE)
USI.SYSTEM<-unique(metadata$SYSTEM)
USI.HABITAT<-unique(metadata$HABITAT)
USI.ORGANISM<-unique(metadata$ORG.TYPE)
USI.LF<-unique(metadata$LAB.FIELD)
USI.R2<-unique(metadata$RESP.CAT)
USI.EUT<-unique(metadata$EXP.UNIT.TYPE)
USI.LTT<-unique(metadata$LIGHT.TRT.TYPE)


output3<-data.frame()#opens empty data frame

for(i in 1:length(USI)){
  temp<-metadata[metadata$RESPONSE.TYPE==USI[i], ]#creates a temp data set for each experimental unit
  if(dim(temp)[1]>1){ # checks whether more than 1 date has been measure
    mod1<-rma.mv(LRR,var.LRR,method="REML", random = ~ 1 | UPI/UEI, data=temp)
    ES.ML<-mod1$b
    CIL.ML<-mod1$ci.lb
    CIU.ML<-mod1$ci.ub
    k.ML<-mod1$k
    cat1<-"All"
    cat2<-temp$RESPONSE.TYPE[1]
    cat3<-"All"
    #gets the regression coefficients
    output3<-rbind(output3,data.frame(cat1,cat2,cat3,ES.ML,CIL.ML,CIU.ML,k.ML))#binds this in a table together with the first 6 columns of the original file (to get treatment, mesocosm id, etc, to be changed if needed)
  }
    for (i in 1:length (USI.HABITAT)){
    d2<-temp[temp$HABITAT==USI.HABITAT[i], ]
    if(dim(d2)[1]>1){ # checks whether more than 1 date has been measure
    mod1<-rma.mv(LRR,var.LRR,method="REML", random = ~ 1 | UPI/UEI, data=d2)
    ES.ML<-mod1$b
    CIL.ML<-mod1$ci.lb
    CIU.ML<-mod1$ci.ub
    k.ML<-mod1$k
    cat1<-"Habitat"
    cat2<-d2$RESPONSE.TYPE[1]
    cat3<-d2$HABITAT[1]
    #gets the regression coefficients
    output3<-rbind(output3,data.frame(cat1,cat2,cat3,ES.ML,CIL.ML,CIU.ML,k.ML))#binds this in a table together with the first 6 columns of the original file (to get treatment, mesocosm id, etc, to be changed if needed)
    rm(d2)
    }
  }
  for (i in 1:length (USI.ORGANISM)){
    d2<-temp[temp$ORG.TYPE==USI.ORGANISM[i], ]
    if(dim(d2)[1]>1){ # checks whether more than 1 date has been measure
      mod1<-rma.mv(LRR,var.LRR,method="REML", random = ~ 1 | UPI/UEI, data=d2)
      ES.ML<-mod1$b
      CIL.ML<-mod1$ci.lb
      CIU.ML<-mod1$ci.ub
      k.ML<-mod1$k
      cat1<-"Organism group"
      cat2<-d2$RESPONSE.TYPE[1]
      cat3<-d2$ORG.TYPE[1]
      #gets the regression coefficients
      output3<-rbind(output3,data.frame(cat1,cat2,cat3,ES.ML,CIL.ML,CIU.ML,k.ML))#binds this in a table together with the first 6 columns of the original file (to get treatment, mesocosm id, etc, to be changed if needed)
      rm(d2)
    }
  }
  for (i in 1:length (USI.EUT)){
    d2<-temp[temp$EXP.UNIT.TYPE==USI.EUT[i], ]
    if(dim(d2)[1]>1){ # checks whether more than 1 date has been measure
      mod1<-rma.mv(LRR,var.LRR,method="REML", random = ~ 1 | UPI/UEI, data=d2)
      ES.ML<-mod1$b
      CIL.ML<-mod1$ci.lb
      CIU.ML<-mod1$ci.ub
      k.ML<-mod1$k
      cat1<-"Exp Unit Type"
      cat2<-d2$RESPONSE.TYPE[1]
      cat3<-d2$EXP.UNIT.TYPE[1]
      #gets the regression coefficients
      output3<-rbind(output3,data.frame(cat1,cat2,cat3,ES.ML,CIL.ML,CIU.ML,k.ML))#binds this in a table together with the first 6 columns of the original file (to get treatment, mesocosm id, etc, to be changed if needed)
      rm(d2)
    }
  }
  for (i in 1:length (USI.LTT)){
    d2<-temp[temp$LIGHT.TRT.TYPE==USI.LTT[i], ]
    if(dim(d2)[1]>1){ # checks whether more than 1 date has been measure
      mod1<-rma.mv(LRR,var.LRR,method="REML", random = ~ 1 | UPI/UEI, data=d2)
      ES.ML<-mod1$b
      CIL.ML<-mod1$ci.lb
      CIU.ML<-mod1$ci.ub
      k.ML<-mod1$k
      cat1<-"Light Treat Type"
      cat2<-d2$RESPONSE.TYPE[1]
      cat3<-d2$LIGHT.TRT.TYPE[1]
      #gets the regression coefficients
      output3<-rbind(output3,data.frame(cat1,cat2,cat3,ES.ML,CIL.ML,CIU.ML,k.ML))#binds this in a table together with the first 6 columns of the original file (to get treatment, mesocosm id, etc, to be changed if needed)
      rm(d2)
    }
  }
  for (i in 1:length (USI.SYSTEM)){
    d2<-temp[temp$SYSTEM==USI.SYSTEM[i], ]
    if(dim(d2)[1]>1){ # checks whether more than 1 date has been measure
      mod1<-rma.mv(LRR,var.LRR,method="REML", random = ~ 1 | UPI/UEI, data=d2)
      ES.ML<-mod1$b
      CIL.ML<-mod1$ci.lb
      CIU.ML<-mod1$ci.ub
      k.ML<-mod1$k
      cat1<-"System"
      cat2<-d2$RESPONSE.TYPE[1]
      cat3<-d2$SYSTEM[1]
      #gets the regression coefficients
      output3<-rbind(output3,data.frame(cat1,cat2,cat3,ES.ML,CIL.ML,CIU.ML,k.ML))#binds this in a table together with the first 6 columns of the original file (to get treatment, mesocosm id, etc, to be changed if needed)
      rm(d2)
    }
  }
  for (i in 1:length (USI.LF)){
    d2<-temp[temp$LAB.FIELD==USI.LF[i], ]
    if(dim(d2)[1]>1){ # checks whether more than 1 date has been measure
      mod1<-rma.mv(LRR,var.LRR,method="REML", random = ~ 1 | UPI/UEI, data=d2)
      ES.ML<-mod1$b
      CIL.ML<-mod1$ci.lb
      CIU.ML<-mod1$ci.ub
      k.ML<-mod1$k
      cat1<-"Lab.FIELD"
      cat2<-d2$RESPONSE.TYPE[1]
      cat3<-d2$LAB.FIELD[1]
      #gets the regression coefficients
      output3<-rbind(output3,data.frame(cat1,cat2,cat3,ES.ML,CIL.ML,CIU.ML,k.ML))#binds this in a table together with the first 6 columns of the original file (to get treatment, mesocosm id, etc, to be changed if needed)
      rm(d2)
    }
  }
  for (i in 1:length (USI.R2)){
    d2<-temp[temp$RESP.CAT==USI.R2[i], ]
    if(dim(d2)[1]>1){ # checks whether more than 1 date has been measure
      mod1<-rma.mv(LRR,var.LRR,method="REML", random = ~ 1 | UPI/UEI, data=d2)
       ES.ML<-mod1$b
      CIL.ML<-mod1$ci.lb
      CIU.ML<-mod1$ci.ub
      k.ML<-mod1$k
      cat1<-"Response"
      cat2<-d2$RESPONSE.TYPE[1]
      cat3<-d2$RESP.CAT[1]
      #gets the regression coefficients
      output3<-rbind(output3,data.frame(cat1,cat2,cat3,ES.ML,CIL.ML,CIU.ML,k.ML))#binds this in a table together with the first 6 columns of the original file (to get treatment, mesocosm id, etc, to be changed if needed)
      rm(d2)
      
    }
  }
}

names(output3)
write.csv(output3,file="output3.csv")



# for plotting, create certain orders on the x-axis
summary(output3)
output3$order<-NA
output3$order[output3$cat1=="All"]<-1
output3$order[output3$cat1=="Response"]<-2
output3$order[output3$cat1=="Habitat"]<-3
output3$order[output3$cat1=="Organism group"]<-4
output3$order[output3$cat1=="Lab.FIELD"]<-5
output3$order[output3$cat1=="System"]<-6

output3.appendix<-output3[is.na(output3$order),]
output3<-output3[!is.na(output3$order),]
output3$cat1<-with(output3, reorder(cat1,order))

output3$cat3<-as.factor(output3$cat3)
summary(output3$cat3)
output3$cat3<-as.character(output3$cat3)
output3$cat3[output3$cat3=="abundance"]<-"Abu"
output3$cat3[output3$cat3=="biomass"]<-"Biom"
output3$cat3[output3$cat3=="growth"]<-"Gr.r"
output3$cat3[output3$cat3=="prod.abs"]<-"P.abs"
output3$cat3[output3$cat3=="prod.spec"]<-"P.sp"
output3$cat3[output3$cat3=="cont.carb"]<-"C.org"
output3$cat3[output3$cat3=="cont.pigm"]<-"Pigm"
output3$cat3[output3$cat3=="quant.etr"]<-"Quant"
output3$cat3[output3$cat3=="coastal"]<-"Coast"
output3$cat3[output3$cat3=="offshore"]<-"Offsh"
output3$cat3[output3$cat3=="culture"]<-"Cult"
output3$cat3[output3$cat3=="macrophyte"]<-"Macrophyte"
output3$cat3[output3$cat3=="microalgae"]<-"Microalgae"
output3$cat3[output3$cat3=="benthos"]<-"Benthos"
output3$cat3[output3$cat3=="plankton"]<-"Plankton"


output3$cat3<-as.factor(output3$cat3)
output3$order2<-output3$order
output3$order2[output3$cat3=="All"]<-1
output3$order2[output3$cat3=="Abu"]<-2.1
output3$order2[output3$cat3=="Biom"]<-2.2
output3$order2[output3$cat3=="Gr.r"]<-2.3
output3$order2[output3$cat3=="P.abs"]<-2.4
output3$order2[output3$cat3=="P.sp"]<-2.5
output3$order2[output3$cat3=="C.org"]<-2.6
output3$order2[output3$cat3=="Pigm"]<-2.7
output3$order2[output3$cat3=="Quant"]<-2.8
output3$order2[output3$cat3=="Coast"]<-3.1
output3$order2[output3$cat3=="Offsh"]<-3.2
output3$order2[output3$cat3=="Cult"]<-3.3
output3$order2[output3$cat3=="Macrophyte"]<-4.1
output3$order2[output3$cat3=="Microalgae"]<-4.2
output3$order2[output3$cat3=="Lab"]<-5.1
output3$order2[output3$cat3=="Field"]<-5.2
output3$order2[output3$cat3=="Benthos"]<-6.1
output3$order2[output3$cat3=="Plankton"]<-6.2
output3$cat3<-with(output3, reorder(cat3,order2))


pd<-position_dodge(.5)
final.LRR.cat<-ggplot(output3, aes(cat3, ES.ML, col=cat2))+
  scale_color_manual(values = c("darkgreen","darkorange"))+
  scale_shape_manual(values=c(15,17)) +
  geom_point(size=5,position=pd)+
  geom_errorbar(aes(ymin=CIL.ML,ymax=CIU.ML),position=pd,width=0,size=1.6,alpha=.5)+
  #geom_line(position=pd)+
  ylab("Log response ratio")+
  theme_bw()+
  xlab("Category")+
  theme(axis.title.y=element_text(size=18, colour="black",vjust=0.7),axis.text.y=element_text(size=12,face="bold",colour="black",angle=0,hjust=0.9))+
  theme(axis.title.x=element_text(size=18,colour="black",vjust=0.2),axis.text.x=element_text(size=12,colour="black",face="bold"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(plot.margin=unit(c(0.2,0.2,.2,0.2),"cm"))+
  theme(panel.border=element_rect(colour="black",size=1.5))+
  theme(axis.ticks=element_line(colour="black",size=1),axis.ticks.length=unit(0.3,"cm"))+
  geom_hline(yintercept=0,colour="grey")+theme(legend.title = element_blank())+
  geom_text(aes(cat3,CIL.ML-0.05,label = k.ML),position=pd, size = 3)+
  facet_wrap(~cat1,scales = "free")
final.LRR.cat



pd<-position_dodge(.5)
summary(output3.appendix)
final.LRR.append3<-ggplot(output3.appendix, aes(cat3, ES.ML, col=cat2))+
  scale_color_manual(values = c("darkgreen","darkorange"))+
  scale_shape_manual(values=c(15,17)) +
  geom_point(size=5,position=pd)+
  geom_errorbar(aes(ymin=CIL.ML,ymax=CIU.ML),position=pd,width=0,size=1.6,alpha=.5)+
  #geom_line(position=pd)+
  ylab("Log response ratio")+
  theme_bw()+
  xlab("Category")+
  theme(axis.title.y=element_text(size=18, colour="black",vjust=0.7),axis.text.y=element_text(size=12,face="bold",colour="black",angle=0,hjust=0.9))+
  theme(axis.title.x=element_text(size=18,colour="black",vjust=0.2),axis.text.x=element_text(size=12,colour="black",face="bold"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(plot.margin=unit(c(0.2,0.2,.2,0.2),"cm"))+
  theme(panel.border=element_rect(colour="black",size=1.5))+
  theme(axis.ticks=element_line(colour="black",size=1),axis.ticks.length=unit(0.3,"cm"))+
  geom_hline(yintercept=0,colour="grey")+theme(legend.title = element_blank())+
  geom_text(aes(cat3,CIL.ML-0.05,label = k.ML),position=pd, size = 3)+
  facet_wrap(~cat1,scales = "free")
final.LRR.append3



# create publication ready figures and table

tiff(file = "figures/fig1.tiff", width = 4200, height = 2800, units = "px", res = 300)
grid.arrange(tag_facet(final.LRR.cat, 
                       hjust = -0.5,
                       open = "", close = "",
                       tag_pool = LETTERS,size=6))
dev.off() 


tiff(file = "figures/fig2.tiff", width = 4500, height = 3600, units = "px", res = 300)
plot_grid(light.LRR,dur.LRR,temp.LRR,year.LRR, ncol=2,labels="AUTO",align="v", label_size = 18)
dev.off() 

tiff(file = "figures/appD_figS2.tiff", width = 3600, height = 3600, units = "px", res = 300)
plot_grid(irr.LRR,TN,size.LRR,lat.LRR, ncol=2,labels="AUTO",align="v", label_size = 18)
dev.off() 

tiff(file = "figures/appD_figS3.tiff", width = 3200, height = 2400, units = "px", res = 300)
grid.arrange(tag_facet(final.LRR.append3, 
                       hjust = -0.5,
                       open = "", close = "",
                       tag_pool = LETTERS,size=6))

dev.off() 


tiff(file = "figures/appC_FigS1.tiff", width = 2400, height = 2400, units = "px", res = 300)
par(mfrow=c(1,2))
funnel(mod.test,yaxis = "wi", main="A: no moderators")
funnel(test.complete,yaxis = "wi", main="B: all moderators")
dev.off()


library(kableExtra)
kbl(complete.mod[,c(1,3,4,6)], caption = "", digits=4) %>%
  kable_styling(bootstrap_options = "striped",
                full_width = F) 


######################################
######## single moderator analyses ###
######################################

#Creating the univariate ML-meta-analyses


"Year"
aoh1<-rma.mv(LRR,var.LRR,method="REML", random = ~ 1 | UPI/UEI, data=metadata,
         mods=~YEAR)
summary(aoh1)
year<-tidy(summary(aoh1))
year$k<-aoh1$k
year$AIC<-aoh1$fit.stats$REML[3]
year$model<-"year"

"Lab.Field"
aoh1<-rma.mv(LRR,var.LRR,method="REML", random = ~ 1 | UPI/UEI, data=metadata,
         mods=~LAB.FIELD)
summary(aoh1)
Lab.field<-tidy(summary(aoh1))
Lab.field$k<-aoh1$k
Lab.field$AIC<-aoh1$fit.stats$REML[3]
Lab.field$model<-"Lab.field"

"Experiment time"
aoh1<-rma.mv(LRR,var.LRR,method="REML", random = ~ 1 | UPI/UEI, data=metadata,
             mods=~EXP.TIME)
summary(aoh1)
exp.time<-tidy(summary(aoh1))
exp.time$k<-aoh1$k
exp.time$AIC<-aoh1$fit.stats$REML[3]
exp.time$model<-"exp.time"


"Habitat"
aoh1<-rma.mv(LRR,var.LRR,method="REML", random = ~ 1 | UPI/UEI, data=metadata,
         mods=~HABITAT)
summary(aoh1)
habitat<-tidy(summary(aoh1))
habitat$k<-aoh1$k
habitat$AIC<-aoh1$fit.stats$REML[3]
habitat$model<-"habitat"

"Latitude"
aoh1<-rma.mv(LRR,var.LRR,method="REML", random = ~ 1 | UPI/UEI, data=metadata,
         mods=~abs(LATITUDE))
summary(aoh1)
lat<-tidy(summary(aoh1))
lat$k<-aoh1$k
lat$AIC<-aoh1$fit.stats$REML[3]
lat$model<-"lat"

"Organism"
aoh1<-rma.mv(LRR,var.LRR,method="REML", random = ~ 1 | UPI/UEI, data=metadata,
         mods=~ORG.TYPE)
summary(aoh1)
organism<-tidy(summary(aoh1))
organism$k<-aoh1$k
organism$AIC<-aoh1$fit.stats$REML[3]
organism$model<-"organism"

"System"
aoh1<-rma.mv(LRR,var.LRR,method="REML", random = ~ 1 | UPI/UEI, data=metadata,
             mods=~SYSTEM)
summary(aoh1)
system<-tidy(summary(aoh1))
system$k<-aoh1$k
system$AIC<-aoh1$fit.stats$REML[3]
system$model<-"system"

"Duration"
aoh1<-rma.mv(LRR,var.LRR,method="REML", random = ~ 1 | UPI/UEI, data=metadata,
         mods=~DURATION)
summary(aoh1)
duration<-tidy(summary(aoh1))
duration$k<-aoh1$k
duration$AIC<-aoh1$fit.stats$REML[3]
duration$model<-"duration"

"Exp.unit.type"
aoh1<-rma.mv(LRR,var.LRR,method="REML", random = ~ 1 | UPI/UEI, data=metadata,
         mods=~EXP.UNIT.TYPE)
summary(aoh1)
exp.unit<-tidy(summary(aoh1))
exp.unit$k<-aoh1$k
exp.unit$AIC<-aoh1$fit.stats$REML[3]
exp.unit$model<-"exp.unit"

"Size"
aoh1<-rma.mv(LRR,var.LRR,method="REML", random = ~ 1 | UPI/UEI, data=metadata,
         mods=~SIZE)
summary(aoh1)
size<-tidy(summary(aoh1))
size$k<-aoh1$k
size$AIC<-aoh1$fit.stats$REML[3]
size$model<-"size"


"Temperature"
aoh1<-rma.mv(LRR,var.LRR,method="REML", random = ~ 1 | UPI/UEI, data=metadata,
         mods=~TEMP)
summary(aoh1)
temperature<-tidy(summary(aoh1))
temperature$k<-aoh1$k
temperature$AIC<-aoh1$fit.stats$REML[3]
temperature$model<-"temperature"

"Remaining light"
aoh1<-rma.mv(LRR,var.LRR,method="REML", random = ~ 1 | UPI/UEI, data=metadata,
         mods=~Light.REG.TRT)
summary(aoh1)
rem.light<-tidy(summary(aoh1))
rem.light$k<-aoh1$k
rem.light$AIC<-aoh1$fit.stats$REML[3]
rem.light$model<-"rem.light"



"IRR"
aoh1<-rma.mv(LRR,var.LRR,method="REML", random = ~ 1 | UPI/UEI, data=metadata,
         mods=~CON.IRR)
summary(aoh1)
control.irr<-tidy(summary(aoh1))
control.irr$k<-aoh1$k
control.irr$AIC<-aoh1$fit.stats$REML[3]
control.irr$model<-"control.irr"

"Response.Type"
aoh1<-rma.mv(LRR,var.LRR,method="REML", random = ~ 1 | UPI/UEI, data=metadata,
         mods=~RESPONSE.TYPE)
summary(aoh1)
resptype<-tidy(summary(aoh1))
resptype$k<-aoh1$k
resptype$AIC<-aoh1$fit.stats$REML[3]
resptype$model<-"resptype"

"Response.Type2 = Category"
aoh1<-rma.mv(LRR,var.LRR,method="REML", random = ~ 1 | UPI/UEI, data=metadata,
         mods=~RESP.CAT)
summary(aoh1)
respcat<-tidy(summary(aoh1))
respcat$k<-aoh1$k
respcat$AIC<-aoh1$fit.stats$REML[3]
respcat$model<-"respcat"


"Light treatment type"
aoh1<-rma.mv(LRR,var.LRR,method="REML", random = ~ 1 | UPI/UEI, data=metadata,
         mods=~LIGHT.TRT.TYPE)
summary(aoh1)
Light.treat<-tidy(summary(aoh1))
Light.treat$k<-aoh1$k
Light.treat$AIC<-aoh1$fit.stats$REML[3]
Light.treat$model<-"Light.treat"

"TN"
aoh1<-rma.mv(LRR,var.LRR,method="REML", random = ~ 1 | UPI/UEI, data=metadata,
         mods=~TN)
summary(aoh1)
TN<-tidy(summary(aoh1))
TN$k<-aoh1$k
TN$AIC<-aoh1$fit.stats$REML[3]
TN$model<-"TN"


output4<-rbind(year,habitat,organism, system, duration, 
               temperature,rem.light,control.irr,respcat,resptype,
               exp.unit, size, Lab.field,lat, TN,Light.treat)
write.csv(output4, file="single.csv")



##########################
### plot predicted  ######
##########################


### fit mixed-effects model with Light.REG.TRT and response type as single predictors
res<-rma.mv(LRR,var.LRR,method="REML", 
            random = ~ 1 | UPI/UEI, 
            data=metadata,
             mods=~Light.REG.TRT+RESPONSE.TYPE)
summary(res)

### do the prediction
preds <- predict(test.complete,addx=TRUE)

predicted<-as.data.frame((preds$X))
predicted$LRR<-preds$pred
predicted$cilb<-preds$ci.lb
predicted$ciub<-preds$ci.ub

summary(predicted)



predicted.LRR<-ggplot(predicted, aes(Light.REG.TRT, (exp(LRR)*100),  colour=as.factor(RESPONSE.TYPEbiom),fill=as.factor(RESPONSE.TYPEbiom)))+
  geom_point(alpha=.25,shape=16)+
  geom_smooth(method="loess",se=TRUE)+
  #geom_smooth(aes(Light.REG.TRT, (exp(cilb)*100)),method="lm",se=FALSE,linetype="dashed")+
  #geom_smooth(aes(Light.REG.TRT, (exp(ciub)*100)),method="lm",se=FALSE,linetype="dashed")+
  ylab("Predicted remaining performance (%)")+
  xlab("Remaining Light (%)")+
  theme_bw()+
  scale_color_manual(values = c("darkorange", "darkgreen"))+
  scale_fill_manual(values = c("orange", "lightgreen"))+
  #scale_shape_manual(values=c(0,15,1,16),labels=c("F macro", "L macro", "F micro", "L micro")) + 
  geom_hline(aes(yintercept=100), color = "darkgrey")+
  geom_hline(aes(yintercept=50), color = "darkgrey",alpha=.2)+
  geom_hline(aes(yintercept=75), color = "darkgrey",alpha=.2)+
  theme(axis.title.y=element_text(size=18, face="plain", colour="black",vjust=0.3),axis.text.y=element_text(size=12,face="bold",colour="black",angle=0,hjust=0.4))+
  theme(axis.title.x=element_text(size=18,face="plain",colour="black",vjust=0),axis.text.x=element_text(size=12,face="bold",colour="black"))+
  theme(axis.ticks=element_line(colour="black",size=1),axis.ticks.length=unit(0.3,"cm"))+
  theme(panel.border=element_rect(colour="black",size=1.5))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm"))+
  theme(legend.justification=c(1,1), legend.position=c(.95,0.95))+
  theme(legend.text = element_text(size=18))+
  theme(legend.title = element_blank())+
  guides(color = "none", fill = "none",shape = guide_legend(override.aes = list(size=5)))


predicted.LRR

### calculate prediction for the model with only that
pred.gen <- predict(res,addx=TRUE)

predicted2<-as.data.frame((pred.gen$X))
predicted2$LRR<-pred.gen$pred
predicted2$cilb<-pred.gen$ci.lb
predicted2$ciub<-pred.gen$ci.ub
summary(predicted2)
plot(LRR~Light.REG.TRT,predicted2)


predicted.LRR2<-predicted.LRR+
  geom_smooth(data=predicted2,aes( Light.REG.TRT, (exp(LRR)*100),colour=as.factor(RESPONSE.TYPEbiom)),
              method="loess",se=FALSE,linetype="dashed")
predicted.LRR2

tiff(file = "figures/fig3.tiff", width = 2400, height = 1800, units = "px", res = 300)
predicted.LRR2
dev.off() 

#####################
### creating a map###
#####################


latlong <- read.csv("data/latlong.csv", sep=";")
names(latlong)
summary(latlong)



library(maps)
library(mapdata)


worldmap <- map_data("worldHires")


largemap<-ggplot()+
  geom_map(data = worldmap, aes(x = long, y = lat, map_id = region),
           map = worldmap, colour = NA, fill = "grey60")+
  geom_point(data = latlong,
            aes(x = LONGITUDE, y = LATITUDE),
            colour = "blue", size=2) +
  xlab(expression("Longitude"))+ 
  ylab(expression("Latitude"))
largemap   


tiff(file = "figures/appD_figS1.tiff", width = 4800, height = 2400, units = "px", res = 400)
plot_grid(largemap)
dev.off() 
