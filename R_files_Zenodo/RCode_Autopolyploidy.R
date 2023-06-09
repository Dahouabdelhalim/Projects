#RCode for "Autopolyploidy alters nodule-level interactions in the legume-rhizobium mutualism"

#ANOVA of Nodule Number
a1 <- aov(c$NodNo ~ c$Ploidy * c$Strain)
summary(a1)

#-------------------------
#Correlation Tests of Internal Nodule Traits
cor.test(c$NodArea, c$ZoneArea)
cor.test(c$NodArea, c$SymArea)
cor.test(c$ZoneArea, c$SymArea) 

#------------------------------
#MANOVA of Internal Nodule Traits
man <- manova(cbind(NodArea, ZoneArea, SymArea) ~ Ploidy*Strain, data = c)
summary(man)
summary.aov(man)

#------------------------------
#Individual ANOVA of Nodule Area
a5 <- aov(c$NodArea ~ c$Ploidy * c$Strain)
summary(a5)
lsmeans(a5, pairwise ~ Ploidy | Strain)

#Calculate Cohen's F for Nod Area
cohens_f(a5)

#Calculate power within my data set - Nod Area
pwr.anova.test(k=2, n=16,f= 0.2485414, sig.level=0.176)

#Number of samples needed to get power of 0.9
pwr.anova.test(k=2,f= 0.2485414, p= 0.9, sig.level=0.05)

#-------------------------
#Individual ANOVA of Zone Area
a6 <- aov(c$ZoneArea ~ c$Ploidy * c$Strain)
summary(a6)
lsmeans(a6, pairwise ~ Ploidy | Strain)

#Calculate Cohen's F for Zone Area
cohens_f(a6)

#Power within my data set - Zone Area
pwr.anova.test(k=2, n=16, f= 0.3150130, sig.level=0.0893)

#Number of samples needed to get power of 0.90
pwr.anova.test(k=2,f= 0.3150130, p= 0.9, sig.level=0.05)

#-----------------
#ANOVA of symbiosome area
a7 <- aov(c$SymArea ~ c$Ploidy * c$Strain)
summary(a7)
lsmeans(a7, pairwise ~ Ploidy | Strain)
