library(lmerTest)
library(interactions)
library(sjPlot)
library(ggplot2)
library(sciplot)
library(MuMIn)
library(multcomp)
library(emmeans)
library(effects)
library(knitr)

fulldata<-read.csv(file.choose(), header = TRUE) #analysis_final.csv
names(fulldata)
fulldata$Harvest <- as.factor(fulldata$Harvest)
fulldata$Light <- as.factor(fulldata$Light)
fulldata$Fungi <- as.factor(fulldata$Fungi)
fulldata$Tree <- as.factor(fulldata$Tree)
fulldata$Shelf <- as.factor(fulldata$Shelf)
levels(fulldata$Fungi)

options(scipen = 999) 

####################################################################################
#Question 1:
#HOW DOES LIGHT AVAILABILITY AFFECT RESOURCE EXCHANGE?

lightNC <- lmer(N.C ~ Harvest + Light + Fungi+
                  Harvest:Light + 
                  Light:Fungi +
                  Harvest:Fungi+
                 # Harvest:Light:Fungi+
                  (1|Shelf) + (1|Tree), data=fulldata)

NCtable <- anova(lightNC) # 3-way interaction highly non-significant, so dropped
capture.output(NCtable,file="S2_NC.doc")

# ' Simple effect' contrasts to test effect of Fungi for the 3 Harvests
emmeans(lightNC, pairwise~Fungi|Harvest,adjustSigma=FALSE)

# Graph: Harvest x Fungi interaction (Fig. 2)
nc.out <- emmeans(lightNC, pairwise~Fungi|Harvest,adjustSigma=FALSE)
nc <- data.frame(head(nc.out,6))
ggplot(nc, aes(emmeans.Harvest, emmeans.emmean, fill=emmeans.Fungi)) + geom_bar(stat="identity", position="dodge2", size=1.0)+ 
  scale_fill_manual("Fungi", values = c("Pisolithus" = "gray21", "Rhizopogon" = "gray85"))+
  geom_errorbar(aes(ymin=emmeans.emmean-emmeans.SE, ymax=emmeans.emmean+emmeans.SE), width=0.4, 
                position=position_dodge(width=0.9)) + theme_minimal(base_size=16)+
  xlab("Harvest (weeks)")+ ylab("N:C transfer (µmol/µmol)") +ylim(0, 0.0017)+
  theme_bw()+theme(panel.grid.major=element_blank(), 
                   panel.border=element_blank(), axis.line=element_line(colour="black"))
# to create a high-res version:
dev.copy(png, file='Figure_2.png', width = 6, 
         height = 5, units = 'in', res = 300)
dev.off()

lightPC <- lmer(P.C ~ Harvest + Light + Fungi+
                  Harvest:Light + 
                  Light:Fungi+
                  Harvest:Fungi+
                  #  Harvest:Light:Fungi+
                  (1|Shelf) + (1|Tree), data=fulldata)
capture.output(anova(lightPC), file="S3_PC.doc")

# Graph: Main effect of Harvest for P:C (Supp Fig. S4)
emmeans(lightPC, pairwise~Harvest)
pc.out <- emmeans(lightPC, pairwise~Harvest)
pc <- data.frame(head(pc.out))
ggplot(pc, aes(emmeans.Harvest, emmeans.emmean)) + geom_bar(stat="identity", position="dodge")+ 
  geom_errorbar(aes(ymin=emmeans.emmean-emmeans.SE, ymax=emmeans.emmean+emmeans.SE), width=0.4, position=position_dodge(width=0.9)) + theme_minimal(base_size=16)+
  xlab("Harvest (weeks)")+ ylab("P:C transfer (µmol/µmol)")+ ylim(0,0.000035)+
  theme_bw()+theme(panel.grid.major=element_blank(), 
                   panel.border=element_blank(), axis.line=element_line(colour="black"))

lightC <- lmer(total.C_umol ~ Harvest + Light + Fungi+
                 Harvest:Light + 
                 Light:Fungi +
                 Harvest:Fungi+
                 Harvest:Light:Fungi+
                 (1|Shelf) + (1|Tree), data=fulldata)
capture.output(anova(lightC), file="S4_Ctot.doc")

# Graph: total C: Light x Harvest x Fungi (Fig. 3)
emmeans(lightC, pairwise~Light|Fungi:Harvest,adjustSigma=FALSE)
ctot.out <- emmeans(lightC, pairwise~Harvest|Fungi:Light,adjustSigma=FALSE)
ctot <- data.frame(head(ctot.out))
ggplot(ctot, aes(emmeans.Harvest, emmeans.emmean, fill=emmeans.Fungi, color=emmeans.Light )) + geom_bar(stat="identity", position="dodge2", size=1.5)+ 
  scale_color_manual("Light", values = c("High" = "orange", "Low" = "steelblue3"))+
  scale_fill_manual("Fungi", values = c("Pisolithus" = "black", "Rhizopogon" = "gray85"))+
  geom_errorbar(aes(ymin=emmeans.emmean-emmeans.SE, ymax=emmeans.emmean+emmeans.SE), width=0.4, position=position_dodge(width=0.9)) + theme_minimal(base_size=16)+
  xlab("Harvest (weeks)")+ ylab("Total C transfer to fungi (µmol)") + ylim(-13000,1500050)+
  theme_bw()+theme(panel.grid.major=element_blank(), 
                   panel.border=element_blank(), axis.line=element_line(colour="black"))

dev.copy(png, file='Figure_3.png', width = 6, 
         height = 5, units = 'in', res = 300)
dev.off()


C.resp <- lmer(respired.C_umol ~ Harvest + Light + Fungi+
                 Harvest:Light + 
                 Light:Fungi +
                 Harvest:Light:Fungi+
                 (1|Shelf) + (1|Tree), data=fulldata)
anova(C.resp)

# Graph: fungal respired C: Light x Harvest x Fungi (very similar to C total, so don't use)
cresp.out <- emmeans(C.resp, pairwise~Harvest|Fungi:Light,adjustSigma=FALSE)
cresp <- data.frame(head(cresp.out))
ggplot(cresp, aes(emmeans.Harvest, emmeans.emmean, fill=emmeans.Fungi, color=emmeans.Light )) + geom_bar(stat="identity", position="dodge2", size=1.5)+ 
  scale_color_manual("Light", values = c("High" = "orange", "Low" = "steelblue3"))+
  scale_fill_manual("Fungi", values = c("Pisolithus" = "black", "Rhizopogon" = "gray85"))+
  geom_errorbar(aes(ymin=emmeans.emmean-emmeans.SE, ymax=emmeans.emmean+emmeans.SE), width=0.4, position=position_dodge(width=0.9)) + theme_minimal(base_size=16)+
  xlab("Harvest (weeks)")+ ylab("Total C respired by fungi (µmol)") + ylim(-13000,1500050)+
  theme_bw()+theme(panel.grid.major=element_blank(), 
                   panel.border=element_blank(), axis.line=element_line(colour="black"))

# calculate mean and SE for proportion of total C that was respired, by Fungi:
fulldata$propn_resp <- fulldata$respired.C_umol/fulldata$total.C_umol
rhizo <- subset(fulldata, Fungi=='Rhizopogon')
piso <- subset(fulldata, Fungi=="Pisolithus")
mean(rhizo$propn_resp)
sd(rhizo$propn_resp)/sqrt(length(rhizo$propn_resp))
mean(piso$propn_resp)
sd(piso$propn_resp)/sqrt(length(piso$propn_resp))

C.bio <- lmer(biomass.C_umol ~ Harvest + Light + Fungi+
                  Harvest:Light + 
                  Light:Fungi+
                  Harvest:Fungi+
                 # Harvest:Light:Fungi+
                  (1|Shelf) + (1|Tree), data=fulldata)
capture.output(anova(C.bio), file="S5_Cbio.doc")

# Graph: fungal biomass C: Harvest x Fungi interaction (Supp. Fig. S5)
cbio.out <- emmeans(C.bio, pairwise~Fungi|Harvest,adjustSigma=FALSE)
cbio <- data.frame(head(cbio.out,6))
ggplot(cbio, aes(emmeans.Harvest, emmeans.emmean, fill=emmeans.Fungi)) + geom_bar(stat="identity", position="dodge2", size=1.0)+ 
  scale_fill_manual("Fungi", values = c("Pisolithus" = "gray21", "Rhizopogon" = "gray85"))+
  geom_errorbar(aes(ymin=emmeans.emmean-emmeans.SE, ymax=emmeans.emmean+emmeans.SE), width=0.4, 
                position=position_dodge(width=0.9)) + theme_minimal(base_size=16)+
  xlab("Harvest (weeks)")+ ylab("C in fungal biomass (µmol)") +ylim(0, 5000)+
  theme_bw()+theme(panel.grid.major=element_blank(), 
                   panel.border=element_blank(), axis.line=element_line(colour="black"))

# calculate and analyze fungal carbon  use efficiency (CUE):
fulldata$CUE_fungi <- fulldata$biomass.C_umol/fulldata$total.C_umol

lightCUE <- lmer(CUE_fungi ~ Harvest + Light + Fungi+
                 Harvest:Light + 
                 Light:Fungi +
                 Harvest:Fungi+
              #    Harvest:Light:Fungi+
                 (1|Shelf) + (1|Tree), data=fulldata)
capture.output(anova(lightCUE), file="S6_CUE.doc")

# Graph: fungal carbon use efficiency: Harvest x Fungi interaction (Supp. Fig. S6)
cue.out <- emmeans(lightCUE, pairwise~Fungi|Harvest,adjustSigma=FALSE)
cue <- data.frame(head(cue.out,6))
ggplot(cue, aes(emmeans.Harvest, emmeans.emmean, fill=emmeans.Fungi)) + geom_bar(stat="identity", position="dodge2", size=1.0)+ 
  scale_fill_manual("Fungi", values = c("Pisolithus" = "gray21", "Rhizopogon" = "gray85"))+
  geom_errorbar(aes(ymin=emmeans.emmean-emmeans.SE, ymax=emmeans.emmean+emmeans.SE), width=0.4, 
                position=position_dodge(width=0.9)) + theme_minimal(base_size=16)+
  xlab("Harvest (weeks)")+ ylab("fungal CUE (biomass C/total C)") +ylim(-0.01,0.1)+
  theme_bw()+theme(panel.grid.major=element_blank(), 
                   panel.border=element_blank(), axis.line=element_line(colour="black"))

lightN <- lmer(total.N ~ Harvest + Light + Fungi+
                 Harvest:Light + 
                 Light:Fungi +
                 Harvest:Fungi+
                # Harvest:Light:Fungi+
                 (1|Shelf) + (1|Tree), data=fulldata)
capture.output(anova(lightN), file="S7_N.doc")

# Graph: total N: Light x Fungi interaction (Fig. 4)
ntot.out <- emmeans(lightN, pairwise~Light|Fungi,adjustSigma=FALSE)
ntot <- data.frame(head(ntot.out,6))
ggplot(ntot, aes(emmeans.Fungi, emmeans.emmean, fill=emmeans.Light )) + geom_bar(stat="identity", position="dodge")+ 
  scale_fill_manual("Light", values = c("High" = "orange", "Low" = "steelblue3"))+
  geom_errorbar(aes(ymin=emmeans.emmean-emmeans.SE, ymax=emmeans.emmean+emmeans.SE), width=0.4, position=position_dodge(width=0.9)) + theme_minimal(base_size=16)+
  xlab("EM Fungal species")+ ylab("N transfer to seedling (µmol)") + ylim(0,140)+
  theme_bw()+theme(panel.grid.major=element_blank(), 
                   panel.border=element_blank(), axis.line=element_line(colour="black"))

dev.copy(png, file='Figure_4.png', width = 6, 
         height = 4, units = 'in', res = 300)
dev.off()

# Graph: total N: main effect of Harvest (Supp. Fig. S7)
emmeans(lightN, pairwise~Harvest)
ntot.out2 <- emmeans(lightN, pairwise~Harvest)
ntot2 <- data.frame(head(ntot.out2))
ggplot(ntot2, aes(emmeans.Harvest, emmeans.emmean)) + geom_bar(stat="identity", position="dodge")+ 
  geom_errorbar(aes(ymin=emmeans.emmean-emmeans.SE, ymax=emmeans.emmean+emmeans.SE), width=0.4, position=position_dodge(width=0.9)) + theme_minimal(base_size=16)+
  xlab("Harvest (weeks)")+ ylab("N transfer to plant (µmol)")+ ylim(0,160)+
  theme_bw()+theme(panel.grid.major=element_blank(), 
                   panel.border=element_blank(), axis.line=element_line(colour="black"))


lightP <- lmer(total.P ~ Harvest + Light + Fungi+
                 Harvest:Light + 
                 Light:Fungi+
                 Harvest:Fungi+
               #  Harvest:Light:Fungi+
                 (1|Shelf) + (1|Tree), data=fulldata)
capture.output(anova(lightP), file="S8_P.doc")
emmeans(lightP, pairwise~Light)
emmeans(lightP, pairwise~Fungi)
emmeans(lightP, pairwise~Harvest)

# Graph: Total P transfer: main effect of Harvest (Supp. Fig. S8)
ptot.out <- emmeans(lightP, pairwise~Harvest)
ptot <- data.frame(head(ptot.out))
ggplot(ptot, aes(emmeans.Harvest, emmeans.emmean)) + geom_bar(stat="identity", position="dodge")+ 
  geom_errorbar(aes(ymin=emmeans.emmean-emmeans.SE, ymax=emmeans.emmean+emmeans.SE), width=0.4, position=position_dodge(width=0.9)) + theme_minimal(base_size=16)+
  xlab("Harvest (weeks)")+ ylab("P transfer to plant (µmol)")+ ylim(0,6.0)+
  theme_bw()+theme(panel.grid.major=element_blank(), 
                   panel.border=element_blank(), axis.line=element_line(colour="black"))

light.fungalmass <- lmer(fungalmass_ug ~ Harvest + Light + Fungi+
                          Harvest:Light + 
                          Light:Fungi+
                          Harvest:Fungi+
                      #    Harvest:Light:Fungi+
                          (1|Shelf) + (1|Tree), data=fulldata)
capture.output(anova(light.fungalmass), file="S9_fungalmass.doc")

# Graph: fungal biomass: Harvest x Fungi interaction (Supp. Fig. S9)
fungalmass.out <- emmeans(light.fungalmass, pairwise~Fungi|Harvest,adjustSigma=FALSE)
fungalmass <- data.frame(head(fungalmass.out,6))
ggplot(fungalmass, aes(emmeans.Harvest, emmeans.emmean, fill=emmeans.Fungi)) + geom_bar(stat="identity", position="dodge2", size=1.0)+ 
  scale_fill_manual("Fungi", values = c("Pisolithus" = "gray21", "Rhizopogon" = "gray85"))+
  geom_errorbar(aes(ymin=emmeans.emmean-emmeans.SE, ymax=emmeans.emmean+emmeans.SE), width=0.4, 
                position=position_dodge(width=0.9)) + theme_minimal(base_size=16)+
  xlab("Harvest (weeks)")+ ylab("fungal biomass (µg)") +ylim(0,100000)+
  theme_bw()+theme(panel.grid.major=element_blank(), 
                   panel.border=element_blank(), axis.line=element_line(colour="black"))

light.plantmass <- lmer(plantmass_g ~ Harvest + Light + Fungi+
                 Harvest:Light + 
                 Light:Fungi+
                 Harvest:Fungi+
#                 Harvest:Light:Fungi+
                 (1|Shelf) + (1|Tree), data=fulldata)
capture.output(anova(light.plantmass), file="S10_plantmass.doc")

# Graph: Plant biomass ~ Light x Fungi interaction (Supp. Fig. S10)
plantmass.out2 <- emmeans(light.plantmass, pairwise~Light|Fungi,adjustSigma=FALSE)
plantmass2 <- data.frame(head(plantmass.out2,6))
ggplot(plantmass2, aes(emmeans.Fungi, emmeans.emmean, fill=emmeans.Light )) + geom_bar(stat="identity", position="dodge")+ 
  scale_fill_manual("Light", values = c("High" = "orange", "Low" = "steelblue3"))+
  geom_errorbar(aes(ymin=emmeans.emmean-emmeans.SE, ymax=emmeans.emmean+emmeans.SE), width=0.4, position=position_dodge(width=0.9)) + theme_minimal(base_size=16)+
  xlab("EM Fungal species")+ ylab("Plant biomass (g)") + ylim(0,0.5)+
  theme_bw()+theme(panel.grid.major=element_blank(), 
                   panel.border=element_blank(), axis.line=element_line(colour="black"))

# Figure: plant biomass ~ Harvest (Supp. Fig. S11)
plantmass.out1 <- emmeans(light.plantmass, pairwise~Harvest)
plantmass1 <- data.frame(head(plantmass.out1))
ggplot(plantmass1, aes(emmeans.Harvest, emmeans.emmean)) + geom_bar(stat="identity", position="dodge")+ 
  geom_errorbar(aes(ymin=emmeans.emmean-emmeans.SE, ymax=emmeans.emmean+emmeans.SE), width=0.4, position=position_dodge(width=0.9)) + theme_minimal(base_size=16)+
  xlab("Harvest (weeks)")+ ylab("Plant biomass (g)")+ ylim(0,0.4)+
  theme_bw()+theme(panel.grid.major=element_blank(), 
                   panel.border=element_blank(), axis.line=element_line(colour="black"))

# check normality assumption:
hist(resid(lightNC)) 
hist(resid(lightPC)) 
hist(resid(lightC))
hist(resid(C.resp))
hist(resid(C.bio))
hist(resid(lightN))
hist(resid(lightP))

##########################################################################
# Question 2:
#DO RESOURCE EXCHANGE RATIOS OR ABSOLUTE FLUXES AFFECT PINE SEEDLING GROWTH AND FUNGAL BIOMASS ACCUMULATION?
#Use model selection on mixed linear models to compare effects of N:C and P:C ratios, 
# absolute N/P, and absolute C on fungal and plant biomass.
#Comparing models using AICc, from model fit using ML

##########################################################################
#DO RESOURCE EXCHANGE RATIOS AFFECT FUNGAL BIOMASS ACUMULATION?

fung.mass_model <- lmer(fungalmass_ug ~ N.C + P.C + total.N + total.P + total.C_umol +
                            (1|Shelf) + (1|Tree)
                          , data=fulldata, REML=T, na.action="na.fail")

d2 <- dredge(fung.mass_model, rank=AICc, REML=F)
importance(d2)
subset(d2, delta < 4)
summary(model.avg(d2, cumsum(weight) <= .95, rank.args = list(REML = TRUE)))


##########################################################################
# DO RESOURCE EXCHANGE RATIOS AFFECT PLANT BIOMASS ACUMULATION?

plant.mass_model <- lmer(plantmass_g ~ N.C + P.C + total.N + total.P + total.C_umol +
                           (1|Shelf) + (1|Tree)
                         , data=fulldata, REML = T, na.action="na.fail")

d4 <- dredge(plant.mass_model, rank=AICc, REML=F)
importance(d4)
subset(d4, delta < 4) 
summary(model.avg(d4, cumsum(weight) <= .95, rank.args = list(REML = TRUE)))





