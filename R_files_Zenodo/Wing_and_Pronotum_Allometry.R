library(ggplot2)
#install.packages("ggpubr")
library(ggpubr)

alodata.leg<-read.csv("AllometryAnalysisDataTransformed_final.csv", header=TRUE)
wing.pronotum.dat<-read.csv("AllometryAnalysisData_final.csv", header = TRUE)
is.data.frame(alodata)
is.data.frame(wing.pronotum.dat)
summary(alodata)
Genus<-alodata$Genus
colors.aloplot<-c("darkslateblue", "lightskyblue3", "plum4", "goldenrod2", "sienna3", "maroon")
Species.aloplot<-c("Aetalion reticulatum", "Entylia carinata","Ennya chrysura", "Membracis mexicana","Metheisa lucillodes", "Polyglypta costata")
names(Species.aloplot) <- c("Aet","Eca", "Ech", "Mem", "Met", "Pol")
genus<-factor(Genus, levels =c("Aet", "Ech", "Eca", "Mem", "Met", "Pol"))
stage<-SAdata$Stage
stage<-as.factor(stage)

Species2<-c("Aetalion reticulatum", "Entylia carinata","Ennya chrysura", "Membracis mexicana","Metheisa lucillodes", "Polyglypta costata")
names(Species2) <- c("Aet", "Eca", "Ech", "Mem", "Met", "Pol")

alodata$log.Head.Width<-as.numeric(alodata$log.Head.Width)
alodata$log.Size.measurement<-as.numeric(alodata$log.Size.measurement)

alodata<-subset(alodata.leg, Tissue !="Leg")
## Allometry Plot (final)
tiff("Pronotum vs Wings Allometry Plot w Equation and R2 Final.tiff", units="in", width=10, height=5, res=200)
allometry.plot<-ggplot(data=alodata, aes(x=log.Head.Width, y=log.Size.measurement, colour=Genus)) +geom_point() +geom_smooth(method = "lm", alpha = .15)
allometry.plot +  
  facet_grid(Tissue~Genus, labeller = labeller(Genus=Species.aloplot))+ theme_bw() + scale_color_manual(values=colors.aloplot)+ 
  labs(y="log(Trait Size, mm)", x="log(Body Size, mm)") + 
  stat_regline_equation(label.x=-0.2, label.y=0.7) + 
  stat_regline_equation(label.x = -0.2, label.y = 0.6, aes(label= ..rr.label..)) + 
  theme(strip.text.x = element_text(face = "italic"), legend.position = "none")
dev.off()


## Allometry Plot (final, with Leg)
tiff("Pronotum vs Wings Allometry Plot w Equation and R2 Final with Leg.tiff", units="in", width=10, height=5, res=200)
allometry.plot<-ggplot(data=alodata.leg, aes(x=log.Head.Width, y=log.Size.measurement, colour=Genus)) +geom_point() +geom_smooth(method = "lm", alpha = .15)
allometry.plot +  
  facet_grid(Tissue~Genus, labeller = labeller(Genus=Species.aloplot))+ theme_bw() + scale_color_manual(values=colors.aloplot)+ 
  labs(y="log(Trait Size, mm)", x="log(Body Size, mm)") + 
  stat_regline_equation(label.x=-0.2, label.y=0.7) + 
  stat_regline_equation(label.x = -0.2, label.y = 0.6, aes(label= ..rr.label..)) + 
  theme(strip.text.x = element_text(face = "italic"), legend.position = "none")
dev.off()


## Plot of the raw Surface Area data of the Wings and the Pronotum 
tiff("Pronotum vs Wings Surface Area Plot CHECK.tiff", units="in", width=10, height=5, res=200)
wing.pronotum.plot<-ggplot(data = wing.pronotum.dat, aes(y = Pronotum.SA, x = Wing.SA, colour = Genus)) + geom_point(aes(shape=Stage), size=2) +geom_smooth(method = "lm", alpha = .15) 
wing.pronotum.plot + theme_classic() + scale_color_manual(values=colors.aloplot, labels = Species2, guide =guide_legend(label.theme = element_text(angle = 0, face = "italic"))) + labs(x= expression('Wing Surface Area (mm'^2*')'), y=expression('Pronotum Surface Area (mm'^2*')'), colour="Species", shape="Stages")
dev.off()


## Pronotum v Wing - Slopes
AetalionDat<-subset(wing.pronotum.dat, Genus=="Aet")
EntyliaDat<-subset(wing.pronotum.dat, Genus=="Eca")
EnnyaDat<-subset(wing.pronotum.dat, Genus=="Ech")
MembracisDat<-subset(wing.pronotum.dat, Genus=="Mem")
MetheisaDat<-subset(wing.pronotum.dat, Genus=="Met")
PolyglyptaDat<-subset(wing.pronotum.dat, Genus=="Pol")

modAet<-lm(Pronotum.SA ~ Wing.SA, data= AetalionDat)
summary(modAet)
modEca<-lm(Pronotum.SA ~ Wing.SA, data= EntyliaDat)
summary(modEca)
modEch<-lm(Pronotum.SA ~ Wing.SA, data= EnnyaDat)
summary(modEch)
modMet<-lm(Pronotum.SA ~ Wing.SA, data= MetheisaDat)
summary(modMet)
modMem<-lm(Pronotum.SA ~ Wing.SA, data= MembracisDat)
summary(modMem)
modPol<-lm(Pronotum.SA ~ Wing.SA, data= PolyglyptaDat)
summary(modPol)


## Significance Testing -- Linear Models with head size
AetalionData<-subset(alodata, Genus=="Aet")
EntyliaData<-subset(alodata, Genus=="Eca")
EnnyaData<-subset(alodata, Genus=="Ech")
MembracisData<-subset(alodata, Genus=="Mem")
MetheisaData<-subset(alodata, Genus=="Met")
PolyglyptaData<-subset(alodata, Genus=="Pol")
WingData<-subset(alodata, Tissue =="Wing")
PronotumData<-subset(alodata, Tissue=="Pronotum")


mod1<-lm(Size.sqrt ~ Head.Width + AetalionData$Tissue, data= AetalionData)
anova(mod1)
mod1.2<-lm(Ratio ~ Tissue, data=AetalionData)
anova(mod1.2)

mod2<-lm(Size.sqrt ~ Head.Width + Tissue, data= EntyliaData)
anova(mod2)
mod2.2<-lm(Ratio ~ Tissue, data=EntyliaData)
anova(mod2.2)

mod3<-lm(Size.sqrt ~ Head.Width + Tissue, data= EnnyaData)
anova(mod3)
mod3.2<-lm(Ratio ~ Tissue, data= EnnyaData)
anova(mod3.2)

mod4<-lm(Size.sqrt ~ Head.Width + Tissue, data= MembracisData)
anova(mod4)
mod4.2<-lm(Ratio ~  Tissue, data= MembracisData)
anova(mod4.2)

mod5<-lm(Size.sqrt ~ Head.Width + Tissue, data= MetheisaData)
anova(mod5)
mod5.2<-lm(Ratio ~ Tissue, data= MetheisaData)
anova(mod5.2)

mod6<-lm(Size.sqrt ~ Head.Width + Tissue, data= PolyglyptaData)
anova(mod6)
mod6.2<-lm(Ratio ~ Tissue, data= PolyglyptaData)
anova(mod6.2)

mod.genusW<-lm(Ratio ~ Genus, data=WingData)
anova(mod.genusW)

mod.genusP<-lm(Ratio ~ Genus, data=PronotumData)
anova(mod.genusP)

## Means of the ratios
Aet.ratio<-subset(alodata, Genus=="Aet")
Aet.ratio<-alodata$Ratio


## Barplot of the Ratios
tiff("Pronotum vs Wings Ratio Plot.tiff", units="in", width=10, height=5, res=200)
ratio.plot<-ggplot(data=alodata, aes(x=factor(Stage), y=Ratio, colour=Genus)) + geom_boxplot()
ratio.plot +  
  facet_grid(Tissue~Genus)+ theme_bw() + scale_color_manual(values=colors.aloplot) + labs(y="Ratio of Trait Size to Body Size", x="Developmental Stage") + 
  stat_regline_equation(label.x=-0.1, label.y=0.7) + theme(legend.position = "none") 
dev.off()

