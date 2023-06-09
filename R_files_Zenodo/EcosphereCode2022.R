###################################### LOAD PACKAGES ##################################### 

setwd("~/Desktop/R") # Set the working directory
library(plyr) # Necessary for Rmisc
library(Rmisc) # Necessary for summarySE
library(dplyr) # Necessary for group_by
library(reshape) # Necessary for melt & cast
library(ggplot2) # Necessary for ggplots
library(MASS) # Necessary for glm.nb
library(binr) # Necessary for bins
library(scales) # Necessary to add commas to ggplot continuous variables.
library(vegan) # Needed for diversity indeces
library(cluster) # Necessary to make silhouettes.
source("Codes/hcoplot.R")
source("Codes/cleanplot.pca.R")
source("Codes/evplot.R") # Necessary for evplot function
library(gclus) # Necessary for hcoplots.
library(dendextend) # Necessary for dendrogram
library(grid) # Necessary for pushViewPort
library(pairwiseAdonis) # Necessary for pairwise PERMANOVA
library(chisq.posthoc.test) # Necessary for post-hoc Chi-Square
library(ggpubr) # Necessary for ggarrange
library(ape) # Necessary for pcoa
library(hilldiv) # Necessary for diversity comparisons
library(RColorBrewer) # Necessary for color stuff later on
library(scales)   
library(cowplot) # Alternative multiple figure display

# The errorbars overlapped, so use position_dodge to move them horizontally
pd <- position_dodge(0.1) # move them .05 to the left and right

##################################### DATAFRAME SETUP #################################### 

Basal <- read.csv("Data/EcosphereData2022.csv")
Basal$Plot <- as.factor(Basal$Plot)
Basal$Subplot <- as.factor(Basal$Subplot)
Basal$Treatment <- as.factor(Basal$Treatment)
Basal$Family <- as.factor(Basal$Family)

BasalW <- droplevels(subset(Basal, Basal$Origin == "Recruited")) # Exclude Planted trees
BasalP <- droplevels(subset(Basal, Basal$Origin == "Planted")) # Exclude Recruits

Census <- read.csv("Data/LateCensus.csv") 
Slope <- read.csv("Data/TuxtlasSlope.csv") # Load this up in order to add the treatment.

Census$Treatment <- with(Slope, Treatment[match(Census$Plot, Plot)]) # Match Treatment by plot. 
Census$Treatment <- relevel(Census$Treatment, ref = "Wind") 
Census$SuperPlot <- paste(Census$Treatment, Census$Plot)

##################################### DATAFRAME SETUP CENSUS 2018 #################################### 

CensusPlot <- aggregate(cbind(Ocotea.uxpanapana, Trophis.mexicana, Virola.guatemalensis, Psychotria.limonensis, Faramea.occidentalis, Nectandra.ambigens) ~ Plot, Census, sum)
CensusPlot$Treatment <- with(Slope, Treatment[match(CensusPlot$Plot, Plot)]) # Match Treatment by plot. 
CensusPlot$SuperPlot <- paste(CensusPlot$Treatment, CensusPlot$Plot)
CensusMelt <- melt(CensusPlot, id=c("Plot", "Treatment", "SuperPlot"))
colnames(CensusMelt) <- c("Plot", "Treatment", "SuperPlot", "Species", "Frequency") # Rename columns
CensusMelt <- droplevels(subset(CensusMelt, CensusMelt$Frequency != "0")) # Subset without zeroes.
CensusLS <- CensusMelt[rep(row.names(CensusMelt), CensusMelt$Frequency), 1:5] # This is to multiply rows by the value they hold.
CensusLS <- CensusLS[-5]

################################## SPECIES X TREATMENT ################################### 

summary(Census)
Table1 <- ddply(Census,~ Treatment,summarise,OCUX=sum(Ocotea.uxpanapana),VIGU=sum(Virola.guatemalensis), TRME=sum(Trophis.mexicana), FAOC=sum(Faramea.occidentalis), PSLI=sum(Psychotria.limonensis), NEAM=sum(Nectandra.ambigens), Total=sum(Total))
TableP <- ddply(Census,~ Plot,summarise,OCUX=sum(Ocotea.uxpanapana),VIGU=sum(Virola.guatemalensis), TRME=sum(Trophis.mexicana), FAOC=sum(Faramea.occidentalis), PSLI=sum(Psychotria.limonensis), NEAM=sum(Nectandra.ambigens), Total=sum(Total))

TableP <- subset(TableP, select = -c(Total) )
TableMeltX <- melt(TableP, id.vars='Plot') # Melt the data for plotting
TableMeltX$Treatment <- with(Census, Treatment[match(TableMeltX$Plot, Plot)])
colnames(TableMeltX) <- c("Plot", "Species", "Number", "Treatment") # Rename columns

TableMelt <- melt(Table1, id.vars='Treatment') # Melt the data for plotting
colnames(TableMelt) <- c("Treatment", "Species", "Number") # Rename columns
TableMelt2 <- subset(TableMelt, Species != "Total") # Without Total
Positions <- c("PSLI", "OCUX", "VIGU", "NEAM", "TRME", "FAOC")
TableMelt2$Treatment <- as.factor(TableMelt2$Treatment)

Table2

Fig1 <- ggplot(TableMelt2, aes(x=Species, y=Number, fill=Treatment)) + 
  geom_bar(stat = "identity", color = "black") +
  xlab("Species") + 
  ylab("N") +  
  scale_fill_manual(values=c("#9999CC", "#CC6666", "#66CC99")) +
  scale_x_discrete(limits = Positions) +
  theme_minimal()

TableMeltX$Treatment <- as.factor(TableMeltX$Treatment)
TableMeltX$Treatment <- relevel(TableMeltX$Treatment, ref = "Wind") 
TableMeltX$Treatment <- relevel(TableMeltX$Treatment, ref = "Animal") 

TableMeltX$Species <- relevel(TableMeltX$Species, ref = "TRME") 
TableMeltX$Species <- relevel(TableMeltX$Species, ref = "NEAM") 
TableMeltX$Species <- relevel(TableMeltX$Species, ref = "VIGU") 
TableMeltX$Species <- relevel(TableMeltX$Species, ref = "OCUX") 
TableMeltX$Species <- relevel(TableMeltX$Species, ref = "PSLI") 

tgc <- summarySE(TableMeltX, measurevar="Number", groupvars=c("Species","Treatment"), na.rm = TRUE)
tgc$Treatment <- revalue(tgc$Treatment, c("Control"="Nat. Succ."))

ggplot(tgc, aes(x=Species, y=Number, fill=Treatment)) + 
  geom_bar(stat="identity", position="dodge", color = "black", width = 0.5) +
  geom_errorbar(aes(ymin=Number-se, ymax=Number+se), width = .08, position = position_dodge(0.5)) +
  xlab("Species") + 
  ylab("N") +  
  geom_text(label = c("a", "b", "b", #PSLI
                      "a", "a", "b", #OCUX
                      "a", " b·", " c·", #VIGU
                      "a", "b", "b", #NEAM
                      "a", "a", "a", #TRME
                      "a", "a", "a"),#FAOC
            aes(y = Number+se, x = Species),vjust = -0.5, size = 5, position = position_dodge(0.5)) +
  scale_fill_manual(values=c("#9999CC", "#66CC99", "#CC6666")) +
  theme_bw(base_size = 15) 

# Models
Table3 <- Table2
Table3$Treatment <- with(Census, Treatment[match(Table3$Plot, Plot)])
Table3$Treatment <- as.factor(Table3$Treatment)

lmPSLI <- lm(PSLI ~ Treatment, data = Table3)
summary(lmPSLI)
# F-statistic: 2.337 on 2 and 21 DF,  p-value: 0.1213
summary(m1 <- glm(PSLI ~ Treatment, family="poisson", data=Table3))
# Animal x Wind    -0.84157    0.16917  -4.975 6.54e-07 ***
# Animal x Control -1.09003    0.18510  -5.889 3.89e-09 ***
# Wind x Control  -0.2485     0.2136  -1.163    0.245    

lmOCUX <- lm(OCUX ~ Treatment, data = Table3)
summary(lmOCUX)
# F-statistic: 2.972 on 2 and 21 DF,  p-value: 0.073

summary(m1 <- glm(OCUX ~ Treatment, family="poisson", data=Table3))
# Animal x Wind    -0.1591     0.1787   -0.89    0.374 
# Animal x Control -1.9169     0.3387   -5.66 1.51e-08 ***
# Wind x Control -1.7579     0.3424  -5.134 2.84e-07 ***  

lmVIGU <- lm(VIGU ~ Treatment, data = Table3)
summary(lmVIGU)
# F-statistic: 4.631 on 2 and 21 DF,  p-value: 0.02157
summary(m1 <- glm(VIGU ~ Treatment, family="poisson", data=Table3))
# Animal x Wind   -1.0561     0.2902  -3.639 0.000274 ***
# Animal x Control -1.8827     0.4057  -4.641 3.47e-06 ***
# Wind x Control  -0.8267     0.4532  -1.824 0.068113 .  

lmNEAM <- lm(NEAM ~ Treatment, data = Table3)
summary(lmNEAM)
# F-statistic: 2.935 on 2 and 21 DF,  p-value: 0.07515
summary(m1 <- glm(NEAM ~ Treatment, family="poisson", data=Table3))
# Animal x Wind   -2.6856     0.5967  -4.501 6.77e-06 ***
# Animal x Control -3.7842     1.0113  -3.742 0.000183 ***
# Wind x Control  -1.0986     1.1547  -0.951   0.3414  

lmTRME <- lm(TRME ~ Treatment, data = Table3)
summary(lmTRME)
# F-statistic: 0.2414 on 2 and 21 DF,  p-value: 0.7877
summary(m1 <- glm(TRME ~ Treatment, family="poisson", data=Table3))
# Animal x Wind   -2.877e-01  5.401e-01  -0.533    0.594
# Animal x Control -6.931e-01  6.124e-01  -1.132    0.258
# Wind x Control  -0.4055     0.6455  -0.628    0.530

lmFAOC <- lm(FAOC ~ Treatment, data = Table3)
summary(lmFAOC)
# F-statistic: 0.1489 on 2 and 21 DF,  p-value: 0.8625
summary(m1 <- glm(FAOC ~ Treatment, family="poisson", data=Table3))
# Animal x Wind   0.6931     1.2247   0.566   0.5714 
# Animal x Control 0.6931     1.2247   0.566   0.5714 
# Wind x Control  8.901e-12  1.000e+00   0.000   1.0000 

aggregate(CensusMelt$Frequency, by=list(Category=CensusMelt$Species), FUN=sum, na.rm="TRUE")

136+205+69+48+5+18
56
##################################### RECRUITMENT X TREATMENT #################################### 

# First with planted trees. 
Table1 <- t(table(Basal$Species, Basal$Plot))
Table2 <- as.data.frame(cbind(Table1, Total = rowSums(Table1)))
Table2[,1] <- rownames(Table2)
Table2 <- Table2[-c(2:144)]
names(Table2)[1] <- "Plot"
names(Table2)[2] <- "Recruits"

Table2$Treatment <- with(Basal, Treatment[match(Table2$Plot, Plot)]) # Match Treatment by Plot

Table2$Treatment <- relevel(Table2$Treatment, ref = "Wind") 
Table2$Treatment <- relevel(Table2$Treatment, ref = "Animal") 
Table2$Treatment <- relevel(Table2$Treatment, ref = "Forest") 

Table1 <- aggregate(Table2$Recruits, by=list(Category=Table2$Treatment), FUN=sum, na.rm="TRUE")

Table2$Treatment <- relevel(Table2$Treatment, ref = "Wind") 
Table2$Treatment <- relevel(Table2$Treatment, ref = "Animal") 
Table2$Treatment <- relevel(Table2$Treatment, ref = "Forest") 

tgc <- summarySE(Table2, measurevar="Recruits", groupvars=c("Treatment"), na.rm = TRUE)

Sup1 <- ggplot(tgc, aes(x=Treatment, y=Recruits, fill=Treatment)) + 
  geom_bar(stat = "identity", color = "black") +
  geom_errorbar(aes(ymin=Recruits-se, ymax=Recruits+se), width=.1) +
  xlab("Habitat") + 
  ylab("Plant Abundance") + 
  scale_fill_manual(values=c("peachpuff4", "#9999CC", "#66CC99", "#CC6666")) +
  geom_text(position=pd, label = c("a", "a", "b", "a"), aes(y = Recruits+se, x = Treatment),vjust = -0.5, size = 6) +
  guides(fill = FALSE) +
  theme_bw(base_size = 20) +
  ylim(0,390) 
Sup1

# Treatment N Recruits          sd          se          ci
# 1       Wind 8  326.000 88.63407922 31.33687923 74.09994460
# 2     Animal 8  214.125 44.57237293 15.75871358 37.26343630
# 3     Forest 8  189.000 43.79171481 15.48270925 36.61078977
# 4 Nat. Succ. 8  172.125 48.51932605 17.15417223 40.56317168

# First Model
Table2$Treatment <- relevel(Table2$Treatment, ref = "Wind") 
lm1 <- lm(Recruits ~ Treatment, data = Table2)
#par(mfrow=c(2,2)); plot(lm1) # Looks okay... 

#ggplot(Table2, aes(x=Recruits)) +
#geom_histogram(aes(y=..density..),binwidth=10,colour="black", fill="white") +
#geom_density(alpha=.2, fill="#FF6666") # A little skew to the right.

# THIS INCLUDES ALL PLANTS
summary(lm1) # lm(formula = Recruits ~ Treatment, data = Table2)
# F-statistic: 10.89044 on 3 and 28 DF,  p-value: 6.499047e-05

# Forest x Animal : t = 0.84593,    0.40477
# Forest x Wind : t = 4.61262, 7.9934e-05 ***
# Forest x Control : t = -0.56816,    0.57446  
# Animal x Wind : t = 3.76669, 0.00078266 ***
# Animal x Control : t = -1.41409, 0.16835998
# Wind x Control : t = -5.18078, 1.6920e-05 ***

# Recruits Only 
BasalX <- droplevels(subset(Basal, Basal$Origin == "Recruited")) # Exclude Planted trees
Table3 <- t(table(BasalX$Species, BasalX$Plot))
Table4 <- as.data.frame(cbind(Table3, Total = rowSums(Table3)))
Table4[,1] <- rownames(Table4)
Table4 <- Table4[-c(2:141)]
names(Table4)[1] <- "Plot"
names(Table4)[2] <- "Recruits"

# Does skewness of species abundance distributions of non-planted species differ by treatment?

Table4DF <- as.data.frame(Table4)

ggplot(Table4DF, aes(x = Recruits)) +
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666") +
  xlim(0,500) +
  scale_color_manual(values = c("peachpuff4", "#9999CC", "#66CC99", "#CC6666")) +
  scale_fill_manual(values = c("peachpuff4", "#9999CC", "#66CC99", "#CC6666")) +
  facet_grid(~Treatment, scales = "free") 

Table4DF.F <- droplevels(subset(Table4DF, Table4DF$Treatment == "Forest"))
Table4DF.A <- droplevels(subset(Table4DF, Table4DF$Treatment == "Animal"))
Table4DF.W <- droplevels(subset(Table4DF, Table4DF$Treatment == "Wind"))
Table4DF.C <- droplevels(subset(Table4DF, Table4DF$Treatment == "Nat. Succ."))

shapiro.test(Table4DF.F$Recruits) # W = 0.82034, p-value = 0.04705
shapiro.test(Table4DF.A$Recruits) # W = 0.959, p-value = 0.8005
shapiro.test(Table4DF.W$Recruits) # W = 0.87251, p-value = 0.1595
shapiro.test(Table4DF.C$Recruits) # W = 0.89176, p-value = 0.243

Table4$Treatment <- with(Basal, Treatment[match(Table4$Plot, Plot)]) # Match Treatment by Plot

Table4$Treatment <- relevel(Table4$Treatment, ref = "Wind") 
Table4$Treatment <- relevel(Table4$Treatment, ref = "Animal") 
Table4$Treatment <- relevel(Table4$Treatment, ref = "Forest") 

Table5 <- aggregate(Table4$Recruits, by=list(Category=Table4$Treatment), FUN=sum, na.rm="TRUE")

tgc <- summarySE(Table4, measurevar="Recruits", groupvars=c("Treatment"), na.rm = TRUE)

Sup1R <- ggplot(tgc, aes(x=Treatment, y=Recruits, fill=Treatment)) + 
  geom_bar(stat = "identity", color = "black") +
  geom_errorbar(aes(ymin=Recruits-se, ymax=Recruits+se), width=.1) +
  xlab("Habitat") + 
  ylab("Plant Abundance") + 
  ylim(0,390) + #
  geom_text(position=pd, label = c("a", "", "", "a"), aes(y = Recruits+se, x = Treatment),vjust = -0.5, size = 6) +
  scale_fill_manual(values=c("peachpuff4", "#9999CC", "#66CC99", "#CC6666")) +
  guides(fill = FALSE) +
  theme_bw(base_size = 20) +
  ylim(0,390) 

Sup1R

# Treatment N Recruits          sd          se          ci
#1     Forest 8  189.000 43.79171481 15.48270925 36.61078977
#2     Animal 8  155.250 40.12391521 14.18594626 33.54443257
#3       Wind 8  264.250 76.61545909 27.08765533 64.05212672
#4 Nat. Succ. 8  172.125 48.51932605 17.15417223 40.56317168

# Second Model
par(mfrow=c(2,2)); plot(lm1) # Looks okay... 

ggplot(Table2, aes(x=Recruits)) +
  geom_histogram(aes(y=..density..),binwidth=10,colour="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666") # A little skew to the right.

Table4$Treatment <- relevel(Table4$Treatment, ref = "Wind") 
Table4$Treatment <- relevel(Table4$Treatment, ref = "Animal") 
Table4$Treatment <- relevel(Table4$Treatment, ref = "Forest") 

# ONLY RECRUITS

lm2 <- lm(Recruits ~ Treatment, data = Table4)
summary(lm2) # lm(formula = Recruits ~ Treatment, data = Table4)
# F-statistic: 6.294506 on 3 and 28 DF,  p-value: 0.002115681
# Forest x Animal : t = -1.24533  0.2233328
# Forest x Wind : t =  2.77662  0.0096847 ** 
# Forest x Control : t = -0.62266  0.5385442 
# Animal x Wind: t = 4.02195 0.00039618 ***
# Animal x Control: t = 0.62266 0.53854416  
# Wind x Control:  t = -3.39928 0.00204631 ** 

# How many trees were planted? 
Table6 <- table(Basal$Treatment,Basal$Origin)
Table7 <- as.data.frame(Table6)

colnames(Table7) <- c("Treatment", "Origin", "Frequency") 
Table7$Origin <- relevel(Table7$Origin, ref = "Recruited") 

ggplot(Table7, aes(x= Treatment, y= Frequency, group = Origin, fill = Origin)) + 
  geom_bar(stat="identity", color="black")+
  xlab("Habitat") + ylab("Trees") + 
  theme_bw(base_size = 20) 

# Richness total with all plants. 
Table1.3 <- table(Basal$Species, Basal$Plot)
Table1.4 <- as.data.frame(Table1.3) 

colnames(Table1.4) <- c("Species", "Plot", "Frequency")
Table1.4$Frequency2[Table1.4$Frequency > 0] <- 1
Table1.4[is.na(Table1.4)] <- 0
Table1.5 <- aggregate(Table1.4$Frequency2, by=list(Category=Table1.4$Plot), FUN=sum, na.rm="TRUE")
colnames(Table1.5) <- c("Plot", "Richness")
Table1.5$Treatment <- with(Basal, Treatment[match(Table1.5$Plot, Plot)]) # Match Treatment by Plot

Table1.5$Treatment <- as.factor(Table1.5$Treatment)

Table1.5$Treatment <- relevel(Table1.5$Treatment, ref = "Wind") 
Table1.5$Treatment <- relevel(Table1.5$Treatment, ref = "Animal") 
Table1.5$Treatment <- relevel(Table1.5$Treatment, ref = "Forest") 

tgc <- summarySE(Table1.5, measurevar="Richness", groupvars=c("Treatment"))

SupRich <- ggplot(tgc, aes(x=Treatment, y=Richness, fill=Treatment)) + 
  geom_bar(stat="identity", color="black")+
  geom_errorbar(aes(ymin=Richness-se, ymax=Richness+se), width=.1) +
  xlab("Habitat") + 
  ylab("Species Density") + 
  ylim(0,42) +
  scale_fill_manual(values=c("peachpuff4", "#9999CC", "#66CC99", "#CC6666")) +
  geom_text(position=pd, label = c("ab", "a", "b", "c"), aes(y = Richness+se, x = Treatment),vjust = -0.5, size = 6) +
  theme_bw(base_size = 20) +
  theme(legend.position="none", axis.title.y = element_text(margin = margin(t = 0, r = 13, b = 0, l = 0)), axis.title.x = element_blank()) 
SupRich

#    Treatment N Richness          sd          se          ci
# 1     Forest 8   34.125 6.243568119 2.207434678 5.219753573
# 2     Animal 8   32.000 3.927922024 1.388730150 3.283824991
# 3       Wind 8   37.125 4.389516732 1.551928524 3.669727824
# 4 Nat.Succ. 8   23.625 2.559994420 0.905094707 2.140208894

Table1.5$Treatment <- relevel(Table1.5$Treatment, ref = "Wind") 
lmx <- lm(data = Table1.5, Richness ~ Treatment)
summary(lmx)

# Richness, All 
# F-statistic: 13.37525 on 3 and 28 DF,  p-value: 1.335289e-05
# Forest x Animal: t = -0.94895    0.35076 
# Forest x Wind: t = 1.33970    0.19112  
# Forest x Control:t = -4.68894 6.4901e-05 ***
# Animal x Wind: t = 2.28865 0.02985023 *
# Animal x Control: t = -3.73999 0.00083996 ***
# Wind x Control: t = -6.02864 1.6944e-06 ***

# Richness recruit only
Table1.6 <- table(BasalW$Species, BasalW$Plot)
Table1.7 <- as.data.frame(Table1.6) 

colnames(Table1.7) <- c("Species", "Plot", "Frequency")
Table1.7$Frequency2[Table1.7$Frequency > 0] <- 1
Table1.7[is.na(Table1.7)] <- 0
Table1.8 <- aggregate(Table1.7$Frequency2, by=list(Category=Table1.7$Plot), FUN=sum, na.rm="TRUE")
colnames(Table1.8) <- c("Plot", "Richness")
Table1.8$Treatment <- with(Basal, Treatment[match(Table1.8$Plot, Plot)]) # Match Treatment by Plot

Table1.8$Treatment <- as.factor(Table1.8$Treatment)

Table1.8$Treatment <- relevel(Table1.8$Treatment, ref = "Wind") 
Table1.8$Treatment <- relevel(Table1.8$Treatment, ref = "Animal") 
Table1.8$Treatment <- relevel(Table1.8$Treatment, ref = "Forest") 

tgc <- summarySE(Table1.8, measurevar="Richness", groupvars=c("Treatment"))

SupRichR <- ggplot(tgc, aes(x=Treatment, y=Richness, fill=Treatment)) + 
  geom_bar(stat="identity", color="black")+
  geom_errorbar(aes(ymin=Richness-se, ymax=Richness+se), width=.1) +
  xlab("") + 
  ylab("Species Density") + 
  ylim(0,42) +
  scale_fill_manual(values=c("peachpuff4","#9999CC","#66CC99", "#CC6666")) +
  geom_text(position=pd, label = c("ab", "", "", "c"), aes(y = Richness+se, x = Treatment),vjust = -0.5, size = 6) +
  theme_bw(base_size = 20) +
  theme(legend.position="none", axis.title.y = element_text(margin = margin(t = 0, r = 13, b = 0, l = 0)), axis.title.x = element_blank()) 
SupRichR

Table1.8$Treatment <- relevel(Table1.8$Treatment, ref = "Wind") 
lmx <- lm(data = Table1.8, Richness ~ Treatment)
summary(lmx)

# Richness, Recruits 
# F-statistic: 8.962066 on 3 and 28 DF,  p-value: 0.000254
# Forest x Animal: t = -4.25190 0.00021314 ***
# Forest x Wind: t = -2.04522 0.05033414 .
# Forest x Control:t = -4.52100 0.00010261 ***
# Animal x Wind: t = 2.20668 0.03570548 * 
# Animal x Control: t = -0.26911 0.78982016   
# Wind x Control: t = -2.47579 0.019609 * 

#pushViewport(viewport(layout=grid.layout(nrow=2,ncol=1))) # All plants.
#print(Sup1, vp=viewport(layout.pos.row=2, layout.pos.col=1))
#print(SupRich, vp=viewport(layout.pos.row=1, layout.pos.col=1)) # 800 x 1200

#pushViewport(viewport(layout=grid.layout(nrow=2,ncol=1))) # Recruits only.
#print(Sup1R, vp=viewport(layout.pos.row=2, layout.pos.col=1))
#print(SupRichR, vp=viewport(layout.pos.row=1, layout.pos.col=1)) # 800 x 1200

# With Plant Abundance as covariate

sum(Table1.8[,2]) # Recruits Species
sum(Table4[,2]) # Recruited Abundance

# Recruited Scale

RecTable <- Table1.8
RecTable$Recruits <- with(Table4, Recruits[match(RecTable$Plot, Plot)]) # Match Treatment by Plot
RecTable$Treatment <- relevel(RecTable$Treatment, ref = "Forest") 

lmx <- lm(data = RecTable, Richness ~ Treatment + Recruits)
summary(lmx)
anova(lmx)

lmx2 <- lm(data = Table1.8, Richness ~ Treatment)
summary(lmx2)

sum(Table1.5[,2]) # Planted Species
sum(Table2[,2]) # Planted Abundance

lmx <- lm(data = Tabla1.1, Diversity ~ Treatment)
Tabla1.1$Treatment <- with(Basal, Treatment[match(Tabla1.1$Plot, Plot)]) # Match Treatment by Plot

# All 

AllTable <- Table1.5

AllTable$Recruits <- with(Table2, Recruits[match(AllTable$Plot, Plot)]) # Match Treatment by Plot
AllTable$Treatment <- relevel(AllTable$Treatment, ref = "Wind") 

lmx <- lm(data = AllTable, Richness ~ Treatment + Recruits)
summary(lmx)


Table1.5
Table2

##################################### DIVERSITY #################################### 

Tabla1 <- diversity(t(Table1.3), "shannon") # SW Diversity Index. All
Tabla2 <- diversity(t(Table1.6), "shannon") # SW Diversity Index. Recruits Only

# Using Shannon. All
Tabla1.1 <- as.data.frame(Tabla1)
Tabla1.1$Plot <- row.names(Tabla1.1)
colnames(Tabla1.1) <- c("Diversity", "Plot")
Tabla1.1$Treatment <- with(Basal, Treatment[match(Tabla1.1$Plot, Plot)]) # Match Treatment by Plot

Tabla1.1$Treatment <- relevel(Tabla1.1$Treatment, ref = "Wind") 
Tabla1.1$Treatment <- relevel(Tabla1.1$Treatment, ref = "Animal") 
Tabla1.1$Treatment <- relevel(Tabla1.1$Treatment, ref = "Forest") 

tgc <- summarySE(Tabla1.1, measurevar="Diversity", groupvars=c("Treatment"))

SupDivdAll <- ggplot(tgc, aes(x=Treatment, y=Diversity, fill=Treatment)) + 
  geom_bar(stat="identity", color="black")+
  geom_errorbar(aes(ymin=Diversity-se, ymax=Diversity+se), width=.1) +
  xlab("") + 
  ylab("Shannon Index (H)") + 
  scale_y_continuous(breaks=seq(0,3.25,0.5), limits = c(0,3.25)) +
  scale_fill_manual(values=c("peachpuff4","#9999CC","#66CC99", "#CC6666")) +
  geom_text(position=pd, label = c("b", "a", "a", "b"), aes(y = Diversity+se, x = Treatment),vjust = -0.5, size = 6) +
  theme_bw(base_size = 20) +
  theme(legend.position="none", axis.title.y = element_text(margin = margin(t = 0, r = 13, b = 0, l = 0)))

# # Using Shannon. Recruits Only

Tabla2.1 <- as.data.frame(Tabla2)
Tabla2.1$Plot <- row.names(Tabla2.1)
colnames(Tabla2.1) <- c("Diversity", "Plot")
Tabla2.1$Treatment <- with(Basal, Treatment[match(Tabla2.1$Plot, Plot)]) # Match Treatment by Plot

Tabla2.1$Treatment <- relevel(Tabla2.1$Treatment, ref = "Wind") 
Tabla2.1$Treatment <- relevel(Tabla2.1$Treatment, ref = "Animal") 
Tabla2.1$Treatment <- relevel(Tabla2.1$Treatment, ref = "Forest") 

tgc <- summarySE(Tabla2.1, measurevar="Diversity", groupvars=c("Treatment"))

SupDivdR <- ggplot(tgc, aes(x=Treatment, y=Diversity, fill=Treatment)) + 
  geom_bar(stat="identity", color="black")+
  geom_errorbar(aes(ymin=Diversity-se, ymax=Diversity+se), width=.1) +
  xlab("") + 
  ylab("Shannon Index (H)") + 
  scale_y_continuous(breaks=seq(0,3.25,0.5), limits = c(0,3.25)) +
  scale_fill_manual(values=c("peachpuff4","#9999CC","#66CC99", "#CC6666")) +
  geom_text(position=pd, label = c("b", "", "", "b"), aes(y = Diversity+se, x = Treatment),vjust = -0.5, size = 6) +
  theme_bw(base_size = 20) +
  theme(legend.position="none", axis.title.y = element_text(margin = margin(t = 0, r = 13, b = 0, l = 0))) 
SupDivdR

Tabla1.1$Treatment <- relevel(Tabla1.1$Treatment, ref = "Wind") 
lmx <- lm(data = Tabla1.1, Diversity ~ Treatment)
summary(lmx)

# Diversity All
# F-statistic: 8.14 on 3 and 28 DF,  p-value: 0.000473
# Forest x Animal: t = 3.88  0.00058 ***
# Forest x Wind: t = 3.86  0.00062 ***
# Forest x Control:t = 0.86  0.39916 
# Animal x Wind: t = -0.03  0.97934
# Animal x Control: t =  -3.03  0.00527 **
# Wind x Control: t = -3.00  0.00562 **  

Tabla2.1$Treatment <- relevel(Tabla2.1$Treatment, ref = "Wind") 
lmx <- lm(data = Tabla2.1, Diversity ~ Treatment)
summary(lmx)

# Diversity only Recruits 
# F-statistic: 0.927 on 3 and 28 DF,  p-value: 0.441
# Forest x Animal: t = 1.27     0.22  
# Forest x Wind: t = 1.57     0.13 
# Forest x Control:t = 0.86     0.40  
# Animal x Wind: t = 0.30     0.76 
# Animal x Control: t = -0.40     0.69  
# Wind x Control: t = -0.71     0.49  

##################################### DESCRIPTIVES #################################### 
BasalJ <- droplevels(subset(BasalW, BasalW$Treatment != "Forest"))
BasalJ$Species  <-  as.factor(BasalJ$Species)

summary(BasalJ$Species) # 72 excluding NAs. 76
unique(BasalJ$Species) # 70 excluding NAs 
summary(BasalJ$Family) # 32 excluding NAs
unique(BasalJ$Family) # 33 total

##################################### DISPERSAL MODE x TREATMENT (Recruits Only) #################################### 

# BasalW includes only recruited species. 
BasalAW <- droplevels(subset(BasalW, BasalW$Dispersal.Mode == "Biotic" | BasalW$Dispersal.Mode == "Abiotic")) # This removes the NA. 
BasalAW1 <- droplevels(subset(BasalAW, BasalAW$Dispersal.Mode == "Biotic")) # Only Biotic
BasalAW2 <- droplevels(subset(BasalAW, BasalAW$Dispersal.Mode == "Abiotic")) # Only Biotic

Table3.3 <- table(BasalAW$Species, BasalAW$Plot)
Table3.4 <- as.data.frame(Table3.3) # Table 3.4 includes only recruits.

Table1.3 <- table(Basal$Species, Basal$Plot)
colnames(Table3.4) <- c("Species", "Plot", "Frequency")
Table3.4$Frequency2[Table3.4$Frequency > 0] <- 1
Table3.4[is.na(Table3.4)] <- 0
Table3.4$Dispersal.Mode <- with(Basal, Dispersal.Mode[match(Table3.4$Species, Species)]) # Match Dispersal Mode by Species
Table3.5 <- aggregate(Table3.4$Frequency2, by=list(Category=Table3.4$Plot, Table3.4$Dispersal.Mode), FUN=sum, na.rm="TRUE")
colnames(Table3.5) <- c("Plot", "Dispersal.Mode", "Richness")
Table3.5$Treatment <- with(Basal, Treatment[match(Table3.5$Plot, Plot)]) # Match Treatment by Plot
Table3.5$Treatment <- as.factor(Table3.5$Treatment)
Table3.5$Treatment <- relevel(Table3.5$Treatment, ref = "Wind") 
Table3.5$Treatment <- relevel(Table3.5$Treatment, ref = "Animal") 
Table3.5$Treatment <- relevel(Table3.5$Treatment, ref = "Forest") 

tgc <- summarySE(Table3.5, measurevar="Richness", groupvars=c("Treatment","Dispersal.Mode"))
tgc$Dispersal.Mode <- revalue(tgc$Dispersal.Mode, c("Abiotic"="Abiotic Dispersal", "Biotic"="Biotic Dispersal"))

DPRichRecruits <- ggplot(tgc, aes(x=Treatment, y=Richness, fill=Treatment)) + 
  geom_bar(stat="identity", color="black")+
  geom_errorbar(aes(ymin=Richness-se, ymax=Richness+se), width=.1) +
  xlab("Habitat") + 
  ylab("Species Density") + 
  ylim(0,35) +
  scale_fill_manual(values=c("peachpuff4", "#9999CC","#66CC99", "#CC6666")) +
  facet_grid(~Dispersal.Mode, scales = "free") +
  theme_bw(base_size = 20) +
  theme(legend.position="none", axis.title.y = element_text(margin = margin(t = 0, r = 13, b = 0, l = 0)), axis.title.x = element_blank()) +
  geom_text(position=pd, label = c("c", "b", "", "b","a", "", "c", "d"), aes(y = Richness+se, x = Treatment),vjust = -0.5, size = 6) 

DPRichRecruits

Table2.5 <- droplevels(subset(Table3.5, Table3.5$Dispersal.Mode != ""))
Table2.5A <- droplevels(subset(Table2.5, Table2.5$Dispersal.Mode == "Biotic")) # Only Animal-Dispersed
Table2.5A$Treatment <- relevel(Table2.5A$Treatment, ref = "Wind") 
lmx <- lm(data = Table2.5A, Richness ~ Treatment)
summary(lmx)

# Biotically-Dispersed, Richness, Recruits Only
# F-statistic: 19.06878 on 3 and 28 DF,  p-value: 6.178753e-07
# Forest x Animal: t = -6.27748 8.7024e-07 ***
# Forest x Wind: t = -4.48392 0.00011352 ***
# Forest x Control:t = -6.78993 2.2471e-07 ***
# Animal x Wind: t = 1.79357    0.08369 .   
# Animal x Control: t = -0.51245    0.61236  
# Wind x Control: t =-2.30601 0.02872753 *   

Table2.5W <- droplevels(subset(Table2.5, Table2.5$Dispersal.Mode == "Abiotic"))
Table2.5W$Treatment <- relevel(Table2.5W$Treatment, ref = "Wind") 
lmx <- lm(data = Table2.5W, Richness ~ Treatment)
summary(lmx)

# Abiotically-Dispersed, Richness, Recruits Only
# F-statistic: 7.685315 on 3 and 28 DF,  p-value: 0.0006736998

# Forest x Animal: t = 3.218 0.003250 **
# Forest x Wind: t =  4.478 0.000115 ***
# Forest x Control: t =  3.638 0.001099 ** 
# Animal x Wind: t = 1.259  0.21829
# Animal x Control: t =  0.420  0.67784 
# Wind x Control:  t =  -0.840 0.408255 

# Individuals by Dispersal Mode
Table3.7 <- aggregate(Table3.4$Frequency, by=list(Category=Table3.4$Plot, Table3.4$Dispersal.Mode), FUN=sum, na.rm="TRUE")
colnames(Table3.7) <- c("Plot", "Dispersal.Mode", "Recruits")

Table3.7$Treatment <- with(Basal, Treatment[match(Table3.7$Plot, Plot)]) # Match Treatment by Plot
Table3.7$Treatment <- relevel(Table3.7$Treatment, ref = "Wind") 
Table3.7$Treatment <- relevel(Table3.7$Treatment, ref = "Animal") 
Table3.7$Treatment <- relevel(Table3.7$Treatment, ref = "Forest") 

tgc <- summarySE(Table3.7, measurevar="Recruits", groupvars=c("Treatment","Dispersal.Mode"))
tgc$Dispersal.Mode <- revalue(tgc$Dispersal.Mode, c("Abiotic"="Abiotic Dispersal", "Biotic"="Biotic Dispersal"))

DPindiRec <- ggplot(tgc, aes(x=Treatment, y=Recruits, fill=Treatment)) + 
  geom_bar(stat="identity", color="black")+
  geom_errorbar(aes(ymin=Recruits-se, ymax=Recruits+se), width=.1) +
  xlab("Habitat") + 
  ylab("Plant Abundance") + 
  ylim(0,280) +
  scale_fill_manual(values=c("peachpuff4", "#9999CC","#66CC99", "#CC6666")) +
  facet_grid(~Dispersal.Mode, scales = "free") +
  theme_bw(base_size = 20) +
  theme(legend.position="none") +
  geom_text(position=pd, label = c("b", "bc", "", "c","bc·", "", "a·", "·c"), aes(y = Recruits+se, x = Treatment),vjust = -0.5, size = 6) 

DPindiRec

Table3.7A <- droplevels(subset(Table3.7, Table3.7$Dispersal.Mode == "Biotic"))
Table3.7A$Treatment <- as.factor(Table3.7A$Treatment)
Table3.7A$Treatment <- relevel(Table3.7A$Treatment, ref = "Forest") 
lmx <- lm(data = Table3.7A, Recruits ~ Treatment)
summary(lmx)

# Abundance, Biotic, Recruits Only

# F-statistic: 5.193144 on 3 and 28 DF,  p-value: 0.005595676
# Forest x Animal: t = -1.73579  0.093595 . 
# Forest x Wind: t = 1.85281  0.074477 . 
# Forest x Control: t = -1.30184  0.203579 
# Animal x Wind: t = 3.58860  0.0012511 ** 
# Animal x Control: t = 0.43395  0.6676519 
# Wind x Control: t = -3.15465  0.0038181 ** 

Table3.7W <- droplevels(subset(Table3.7, Table3.7$Dispersal.Mode == "Abiotic"))
Table3.7W$Treatment <- relevel(Table3.7W$Treatment, ref = "Forest") 
lmx <- lm(data = Table3.7W, Recruits ~ Treatment)
summary(lmx)

# Abundance, Abiotic, Recruits Only

# F-statistic: 8.239 on 3 and 28 DF,  p-value: 0.0004373
# Forest x Animal: t =  2.000  0.05530 .
# Forest x Control: t =  3.000  0.00562 **
# Forest x Wind: t = 4.869 3.97e-05 ***
# Animal x Wind: 2.869 0.007740 **
# Animal x Control: t = 1.000 0.325907
# Wind x Control: t = -1.869  0.07205 .

# Diversity by Dispersal Mode
BasalAW1.1 <- table(BasalAW1$Species, BasalAW1$Plot) # Biotic = 1
BasalAW2.1 <- table(BasalAW2$Species, BasalAW2$Plot) # Abiotic = 2

BasalAW1.2 <- diversity(t(BasalAW1.1), "shannon") # Inverse-S  Diversity Index. Recruits Only
BasalAW2.2 <- diversity(t(BasalAW2.1), "shannon") # Inverse-S  Diversity Index. Recruits Only

# Shannon-Diversity Index. Only Recruits Biotic. 
BasalAW1.3 <- as.data.frame(BasalAW1.2)
BasalAW1.3$Plot <- row.names(BasalAW1.3)
colnames(BasalAW1.3) <- c("Diversity", "Plot")
BasalAW1.3$Treatment <- with(Basal, Treatment[match(BasalAW1.3$Plot, Plot)]) # Match Treatment by Plot

BasalAW1.3$Treatment <- relevel(BasalAW1.3$Treatment, ref = "Wind") 
BasalAW1.3$Treatment <- relevel(BasalAW1.3$Treatment, ref = "Animal") 
BasalAW1.3$Treatment <- relevel(BasalAW1.3$Treatment, ref = "Forest") 

# Shannon-Diversity Index. Only Recruits Abiotic. 
BasalAW2.3 <- as.data.frame(BasalAW2.2)
BasalAW2.3$Plot <- row.names(BasalAW2.3)
colnames(BasalAW2.3) <- c("Diversity", "Plot")
BasalAW2.3$Treatment <- with(Basal, Treatment[match(BasalAW2.3$Plot, Plot)]) # Match Treatment by Plot

BasalAW2.3$Treatment <- relevel(BasalAW2.3$Treatment, ref = "Wind") 
BasalAW2.3$Treatment <- relevel(BasalAW2.3$Treatment, ref = "Animal") 
BasalAW2.3$Treatment <- relevel(BasalAW2.3$Treatment, ref = "Forest") 

# Now rejoining Biotic and Abiotic

tgc1 <- summarySE(BasalAW1.3, measurevar="Diversity", groupvars=c("Treatment"))
tgc2 <- summarySE(BasalAW2.3, measurevar="Diversity", groupvars=c("Treatment"))

tgc1$Dispersal.Mode <- "Biotic"
tgc2$Dispersal.Mode <- "Abiotic"

tgc <- rbind(tgc1,tgc2)

tgc$Dispersal.Mode <- revalue(tgc$Dispersal.Mode, c("Abiotic"="Abiotic Dispersal", "Biotic"="Biotic Dispersal"))

DPDivRecruits <- ggplot(tgc, aes(x=Treatment, y=Diversity, fill=Treatment)) + 
  geom_bar(stat="identity", color="black")+
  geom_errorbar(aes(ymin=Diversity-se, ymax=Diversity+se), width=.1) +
  xlab("Habitat") + 
  ylab("Shannon Index (H)") + 
  scale_y_continuous(breaks=seq(0,3,0.5), limits = c(0,3)) +
  scale_fill_manual(values=c("peachpuff4", "#9999CC","#66CC99", "#CC6666")) +
  facet_grid(~Dispersal.Mode, scales = "free") +
  theme_bw(base_size = 20) +
  theme(legend.position="none", axis.title.y = element_text(margin = margin(t = 0, r = 13, b = 0, l = 0)), axis.title.x = element_blank()) + 
  geom_text(position=pd, label = c("c", "b", "", "b","b", "", "b", "b"), aes(y = Diversity+se, x = Treatment),vjust = -0.5, size = 6) 

DPDivRecruits

BasalAW1.3$Treatment <- relevel(BasalAW1.3$Treatment, ref = "Forest") 
lmx <- lm(data = BasalAW1.3, Diversity ~ Treatment) # Recruits Biotic
summary(lmx)
# F-statistic: 0.228 on 3 and 28 DF,  p-value: 0.876

# Diversity, Abiotic, Recruits Only
# Forest x Animal: t = 0.31     0.76   
# Forest x Wind: t = 0.40     0.69 
# Forest x Control: t = -0.34     0.73 
# Animal x Wind: t = 0.09     0.93 
# Animal x Control: t = -0.66     0.52 
# Wind x Control: t = -0.74     0.46 

BasalAW2.3$Treatment <- relevel(BasalAW2.3$Treatment, ref = "Wind") 
lmx <- lm(data = BasalAW2.3, Diversity ~ Treatment) # Recruits Abiotic
summary(lmx)
# F-statistic: 12.8 on 3 and 28 DF,  p-value: 1.87e-05

# Diversity, Abiotic, Recruits Only
# Forest x Animal: t = 4.75  5.4e-05 
# Forest x Control: t = 5.28  1.3e-05 
# Forest x Wind: t = 4.93  3.4e-05 
# Animal x Wind: t = 0.32     0.75 
# Animal x Control: t = -0.03     0.97 
# Wind x Control: t = -0.35     0.73    

##################################### DISPERSAL MODE x TREATMENT (ALL) #################################### 

BasalA <- droplevels(subset(Basal, Basal$Dispersal.Mode == "Biotic" | Basal$Dispersal.Mode == "Abiotic"))
BasalA1 <- droplevels(subset(BasalA, BasalA$Dispersal.Mode == "Biotic"))
BasalA2 <- droplevels(subset(BasalA, BasalA$Dispersal.Mode == "Abiotic"))

Table4.1 <- table(BasalA$Species, BasalA$Plot)
Table4.2 <- as.data.frame(Table4.1)

colnames(Table4.2) <- c("Species", "Plot", "Frequency")
Table4.2$Frequency2[Table4.2$Frequency > 0] <- 1
Table4.2[is.na(Table4.2)] <- 0
Table4.2$Dispersal.Mode <- with(Basal, Dispersal.Mode[match(Table4.2$Species, Species)]) # Match Dispersal Mode by Species
Table4.3 <- aggregate(Table4.2$Frequency2, by=list(Category=Table4.2$Plot, Table4.2$Dispersal.Mode), FUN=sum, na.rm="TRUE")
colnames(Table4.3) <- c("Plot", "Dispersal.Mode", "Richness")
Table4.3$Treatment <- with(Basal, Treatment[match(Table4.3$Plot, Plot)]) # Match Treatment by Plot

Table4.3$Treatment <- relevel(Table4.3$Treatment, ref = "Wind") 
Table4.3$Treatment <- relevel(Table4.3$Treatment, ref = "Animal") 
Table4.3$Treatment <- relevel(Table4.3$Treatment, ref = "Forest") 

tgc <- summarySE(Table4.3, measurevar="Richness", groupvars=c("Treatment","Dispersal.Mode"))
tgc$Dispersal.Mode <- revalue(tgc$Dispersal.Mode, c("Abiotic"="Abiotic Dispersal", "Biotic"="Biotic Dispersal"))

DPRichAll <- ggplot(tgc, aes(x=Treatment, y=Richness, fill=Treatment)) + 
  geom_bar(stat="identity", color="black")+
  geom_errorbar(aes(ymin=Richness-se, ymax=Richness+se), width=.1) +
  xlab("Habitat") + 
  ylab("Species Density") + 
  scale_fill_manual(values=c("peachpuff4", "#9999CC","#66CC99", "#CC6666")) +
  geom_text(position=pd, label = c("c", "b", "a", "b", "a", "b", "c", "d"), aes(y = Richness+se, x = Treatment),vjust = -0.5, size = 6) +
  ylim(0,35) +
  facet_grid(~Dispersal.Mode, scales = "free") +
  theme_bw(base_size = 20) +
  theme(legend.position="none", axis.title.y = element_text(margin = margin(t = 0, r = 13, b = 0, l = 0)), axis.title.x = element_blank())

DPRichAll

Table4.3X <- droplevels(subset(Table4.3, Table4.3$Dispersal.Mode != ""))
Table4.3A <- droplevels(subset(Table4.3X, Table4.3X$Dispersal.Mode == "Biotic"))
Table4.3A$Treatment <- relevel(Table4.3A$Treatment, ref = "Forest") 
lmx <- lm(data = Table4.3A, Richness ~ Treatment)
summary(lmx)

# Animal-Dispersed, Richness, All
# F-statistic: 16.08157 on 3 and 28 DF,  p-value: 2.850738e-06
# Forest x Animal: t = -2.24624 0.03276134 *  
# Forest x Wind: t =  -4.36769 0.00015571 ***
# Forest x Control: t = -6.61394 3.5667e-07 ***
# Animal x Wind: t = -2.12145 0.04286925 *   
# Animal x Control: t = 4.36769 0.00015571 ***
# Wind x Control: t =  -2.24624 0.03276134 *

Table4.3W <- droplevels(subset(Table4.3X, Table4.3X$Dispersal.Mode == "Abiotic"))
Table4.3W$Treatment <- relevel(Table4.3W$Treatment, ref = "Animal") 
lmx <- lm(data = Table4.3W, Richness ~ Treatment)
summary(lmx)

# Wind-Dispersed, Richness, All 
# F-statistic: 134.1 on 3 and 28 DF,  p-value: < 2.2e-16
# Forest x Animal: t = 4.647 7.28e-05 ***
# Forest x Control: 5.253 1.39e-05 ***
# Forest x Wind: t = 18.991  < 2e-16 ***
# Animal x Wind: t = 14.344 1.99e-14 ***
# Animal x Control: t = 0.606    0.549 
# Wind x Control: t = -13.74 5.75e-14 ***

# By Individuals, All

Table4.4 <- aggregate(Table4.2$Frequency, by=list(Category=Table4.2$Plot, Table4.2$Dispersal.Mode), FUN=sum, na.rm="TRUE")

colnames(Table4.4) <- c("Plot", "Dispersal.Mode", "Recruits")
Table4.4$Treatment <- with(Basal, Treatment[match(Table4.4$Plot, Plot)]) # Match Treatment by Plot

Table4.4$Treatment <- relevel(Table4.4$Treatment, ref = "Wind") 
Table4.4$Treatment <- relevel(Table4.4$Treatment, ref = "Animal") 
Table4.4$Treatment <- relevel(Table4.4$Treatment, ref = "Forest") 

tgc <- summarySE(Table4.4, measurevar="Recruits", groupvars=c("Treatment","Dispersal.Mode"))
tgc$Dispersal.Mode <- revalue(tgc$Dispersal.Mode, c("Abiotic"="Abiotic Dispersal", "Biotic"="Biotic Dispersal"))

DPindiAll <- ggplot(tgc, aes(x=Treatment, y=Recruits, fill=Treatment)) + 
  geom_bar(stat="identity", color="black")+
  geom_errorbar(aes(ymin=Recruits-se, ymax=Recruits+se), width=.1) +
  xlab("Habitat") + 
  ylab(expression(Plant~Abundance)) + 
  ylim(0,280) +
  scale_fill_manual(values=c("peachpuff4", "#9999CC","#66CC99", "#CC6666")) +
  facet_grid(~Dispersal.Mode, scales = "free") +
  geom_text(position=pd, label = c("b", "bc", "a", "c","bc·", "·ab", "a·", "·c"), aes(y = Recruits+se, x = Treatment),vjust = -0.5, size = 6) +
  theme_bw(base_size = 20) +
  theme(legend.position="none")
DPindiAll

Table4.4A <- droplevels(subset(Table4.4, Table4.4$Dispersal.Mode == "Biotic"))
Table4.4A$Treatment <- relevel(Table4.4A$Treatment, ref = "Forest") 
lmx <- lm(data = Table4.4A, Recruits ~ Treatment)
summary(lmx)

# Animal-Dispersed, Abundance, All  
# F-statistic: 3.313561 on 3 and 28 DF,  p-value: 0.03428174
# Forest x Animal: t = 0.55176   0.585493 
# Forest x Wind: t = 1.82319   0.078968 . 
# Forest x Control: t = -1.28103   0.210691 
# Animal x Wind: t =  1.27144   0.214034 
# Animal x Control: t = -1.83279   0.077489 . 
# Wind x Control: t = -3.10422  0.0043337 ** 

Table4.4W <- droplevels(subset(Table4.4, Table4.4$Dispersal.Mode == "Abiotic"))
Table4.4W$Treatment <- relevel(Table4.4W$Treatment, ref = "Wind") 
lmx <- lm(data = Table4.4W, Recruits ~ Treatment)
summary(lmx)

# Wind-Dispersed, Individuals, All 
# F-statistic:  53.4449 on 3 and 28 DF,  p-value: 1.036635e-11

# Forest x Animal: t = 1.45992   0.155441 
# Forest x Wind: t =   11.39375 5.0103e-12 ***
# Forest x Control: t =  2.18989   0.037026 *
# Animal x Wind: t =  9.93383 1.1142e-10 ***
# Animal x Control: t = 0.72996  0.4714770
# Wind x Control: t = -9.20386 5.8144e-10 ***

pushViewport(viewport(layout=grid.layout(nrow=2,ncol=1))) # This sets the multiplot.
print(DPRichAll, vp=viewport(layout.pos.row=1, layout.pos.col=1)) # 1000 x 800
print(DPindiAll, vp=viewport(layout.pos.row=2, layout.pos.col=1))

pushViewport(viewport(layout=grid.layout(nrow=2,ncol=1))) # This sets the multiplot.
print(DPRichRecruits, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(DPindiRec, vp=viewport(layout.pos.row=2, layout.pos.col=1))

# Diversity by Dispersal Mode
BasalA1.4 <- table(BasalA1$Species, BasalA1$Plot) # Biotic = 1
BasalA2.4 <- table(BasalA2$Species, BasalA2$Plot) # Abiotic = 2

BasalA1.5 <- diversity(t(BasalA1.4), "shannon") # Inverse-S  Diversity Index. Recruits Only
BasalA2.5 <- diversity(t(BasalA2.4), "shannon") # Inverse-S  Diversity Index. Recruits Only

# Shannon-Diversity Index. Only Recruits Biotic. 
BasalA1.6 <- as.data.frame(BasalA1.5)
BasalA1.6$Plot <- row.names(BasalA1.6)
colnames(BasalA1.6) <- c("Diversity", "Plot")
BasalA1.6$Treatment <- with(Basal, Treatment[match(BasalA1.6$Plot, Plot)]) # Match Treatment by Plot

BasalA1.6$Treatment <- relevel(BasalA1.6$Treatment, ref = "Wind") 
BasalA1.6$Treatment <- relevel(BasalA1.6$Treatment, ref = "Animal") 
BasalA1.6$Treatment <- relevel(BasalA1.6$Treatment, ref = "Forest") 

# Shannon-Diversity Index. Only Recruits Abiotic. 
BasalA2.6 <- as.data.frame(BasalA2.5)
BasalA2.6$Plot <- row.names(BasalA2.6)
colnames(BasalA2.6) <- c("Diversity", "Plot")
BasalA2.6$Treatment <- with(Basal, Treatment[match(BasalA2.6$Plot, Plot)]) # Match Treatment by Plot

BasalA2.6$Treatment <- relevel(BasalA2.6$Treatment, ref = "Wind") 
BasalA2.6$Treatment <- relevel(BasalA2.6$Treatment, ref = "Animal") 
BasalA2.6$Treatment <- relevel(BasalA2.6$Treatment, ref = "Forest") 

# Now rejoining Biotic and Abiotic

tgc1 <- summarySE(BasalA1.6, measurevar="Diversity", groupvars=c("Treatment"))
tgc2 <- summarySE(BasalA2.6, measurevar="Diversity", groupvars=c("Treatment"))

tgc1$Dispersal.Mode <- "Biotic"
tgc2$Dispersal.Mode <- "Abiotic"

tgc <- rbind(tgc1,tgc2)

tgc$Dispersal.Mode <- revalue(tgc$Dispersal.Mode, c("Abiotic"="Abiotic Dispersal", "Biotic"="Biotic Dispersal"))

DPDivAll <- ggplot(tgc, aes(x=Treatment, y=Diversity, fill=Treatment)) + 
  geom_bar(stat="identity", color="black")+
  geom_errorbar(aes(ymin=Diversity-se, ymax=Diversity+se), width=.1) +
  xlab("Habitat") + 
  ylab("Shannon Index (H)") + 
  scale_y_continuous(breaks=seq(0,3,0.5), limits = c(0,3)) +
  scale_fill_manual(values=c("peachpuff4", "#9999CC","#66CC99", "#CC6666")) +
  facet_grid(~Dispersal.Mode, scales = "free") +
  theme_bw(base_size = 20) +
  theme(legend.position="none", axis.title.y = element_text(margin = margin(t = 0, r = 13, b = 0, l = 0)), axis.title.x = element_blank()) + 
  geom_text(position=pd, label = c("c", "b", "a", "b","b", "a", "b", "b"), aes(y = Diversity+se, x = Treatment),vjust = -0.5, size = 6) 

DPDivAll

BasalA1.6$Treatment <- relevel(BasalA1.6$Treatment, ref = "Animal") 
lmx <- lm(data = BasalA1.6, Diversity ~ Treatment) # All Biotic
summary(lmx)
# F-statistic: 5.73 on 3 and 28 DF,  p-value: 0.00345

# Diversity, Abiotic, Recruits Only
# Forest x Animal: t = 3.35   0.0023 ** 
# Forest x Control: t =  -0.34   0.7347   
# Forest x Wind: t =  0.40   0.6909 
# Animal x Wind: t = -2.95  0.00638 ** 
# Animal x Control: t = -3.69  0.00095 ***
# Wind x Control: t = -0.74   0.4630 

BasalA2.6$Treatment <- relevel(BasalA2.6$Treatment, ref = "Wind") 
lmx <- lm(data = BasalA2.6, Diversity ~ Treatment) # Recruits Abiotic
summary(lmx)
# F-statistic: 58.3 on 3 and 28 DF,  p-value: 3.7e-12

# Diversity, Abiotic, Recruits Only
# Forest x Animal: t = 5.82  2.9e-06 ***
# Forest x Wind: t = 13.17  1.6e-13 ***
# Forest x Control: t =  5.79  3.3e-06 ***
# Animal x Wind: t = 7.35  5.3e-08 ***
# Animal x Control: t = -0.04     0.97 
# Wind x Control: t = -7.39  4.8e-08 ***

# Covariate Model for Dispersal Mode

AllTable2A <- Table4.3A
AllTable2W <- Table4.3W

AllTable2A$Recruits <- with(Table4.4A, Recruits[match(AllTable2A$Plot, Plot)]) 
AllTable2W$Recruits <- with(Table4.4W, Recruits[match(AllTable2W$Plot, Plot)]) 

AllTable2A$Treatment <- relevel(AllTable2A$Treatment, ref = "Wind") 
AllTable2A$Treatment <- relevel(AllTable2A$Treatment, ref = "Forest") 

lmx <- lm(data = AllTable2A, Richness ~ Treatment + Recruits)
summary(lmx)

AllTable2W$Treatment <- relevel(AllTable2W$Treatment, ref = "Wind") 
lmx <- lm(data = AllTable2W, Richness ~ Treatment + Recruits)
summary(lmx)

sum(Table2.5[,3]) # Recruited Species
sum(Table3.7[,3]) # Recruited Abundance

# Covariate Model for All Dispersal Mode
Table4.4A
Table4.4W

Table4.3A
Table4.3W
##################################### LIFE HISTORY x TREATMENT (Recruits Only) #################################### 

BasalPLS <- droplevels(subset(BasalW, BasalW$Life.History == "Pioneer" | BasalW$Life.History == "Non-Pioneer"))
BasalPLS1 <- droplevels(subset(BasalPLS, BasalPLS$Life.History == "Pioneer")) # Only Pioneer
BasalPLS2 <- droplevels(subset(BasalPLS, BasalPLS$Life.History == "Non-Pioneer")) # Only Non-Pioneer

Table2.1 <- table(BasalPLS$Species, BasalPLS$Plot)
Table2.1 <- as.data.frame(Table2.1) # Includes only recruits.

colnames(Table2.1) <- c("Species", "Plot", "Frequency")
Table2.1$Frequency2[Table2.1$Frequency > 0] <- 1
Table2.1[is.na(Table2.1)] <- 0
Table2.1$Life.History <- with(Basal, Life.History[match(Table2.1$Species, Species)]) # Match Dispersal Mode by Species
Table2.2 <- aggregate(Table2.1$Frequency2, by=list(Category=Table2.1$Plot, Table2.1$Life.History), FUN=sum, na.rm="TRUE")
colnames(Table2.2) <- c("Plot", "Life.History", "Richness")
Table2.2$Treatment <- with(Basal, Treatment[match(Table2.2$Plot, Plot)]) # Match Treatment by Plot

Table2.2$Treatment <- relevel(Table2.2$Treatment, ref = "Wind") 
Table2.2$Treatment <- relevel(Table2.2$Treatment, ref = "Animal") 
Table2.2$Treatment <- relevel(Table2.2$Treatment, ref = "Forest") 

tgc <- summarySE(Table2.2, measurevar="Richness", groupvars=c("Treatment","Life.History"))

LFSpeciesR <- ggplot(tgc, aes(x=Treatment, y=Richness, fill=Treatment)) + 
  geom_bar(stat="identity", color="black")+
  geom_errorbar(aes(ymin=Richness-se, ymax=Richness+se), width=.1) +
  xlab("Habitat") + 
  ylab("Species Density") + 
  ylim(0,33) +
  scale_fill_manual(values=c("peachpuff4", "#9999CC","#66CC99", "#CC6666")) +
  facet_grid(~Life.History, scales = "free") +
  geom_text(position=pd, label = c("a", "", "", "c","c", "", "", "b"), aes(y = Richness+se, x = Treatment),vjust = -0.5, size = 6) +
  theme_bw(base_size = 20) +
  theme(legend.position="none", axis.title.y = element_text(margin = margin(t = 0, r = 13, b = 0, l = 0)), axis.title.x = element_blank()) 
LFSpeciesR

Table3.9 <- droplevels(subset(Table2.2, Table2.2$Life.History != ""))
Table3.9LS <- droplevels(subset(Table3.9, Table3.9$Life.History == "Non-Pioneer"))
Table3.9LS$Treatment <- relevel(Table3.9LS$Treatment, ref = "Wind") 
lmx <- lm(data = Table3.9LS, Richness ~ Treatment)
summary(lmx)

# Non-Pioneer, Species Richness, Recruits Only
# F-statistic: 55.24199 on 3 and 28 DF,  p-value: 6.998297e-12

# Forest x Animal: t = -10.34264 4.5513e-11 ***
# Forest x Wind: t =  -8.89335 1.1993e-09 ***
# Forest x Control: t =  -11.59430 3.3386e-12 ***
# Animal x Wind: t = 1.44929    0.15837 
# Animal x Control: t =  -1.25166    0.22105
# Wind x Control: t = 2.70094   0.011602 * 

Table3.9P <- droplevels(subset(Table3.9, Table3.9$Life.History == "Pioneer"))
Table3.9P$Treatment <- relevel(Table3.9P$Treatment, ref = "Wind") 
lmx <- lm(data = Table3.9P, Richness ~ Treatment)
summary(lmx)

# Pioneer, Species Richness, Recruits Only
# F-statistic: 37.60621 on 3 and 28 DF,  p-value: 5.905602e-10

# Forest x Animal: t = 7.67280 2.3347e-08 ***
# Forest x Wind: t =  9.07636 7.8155e-10 ***
# Forest x Control: t = 8.98279 9.7232e-10 ***
# Animal x Wind: t = 1.40356    0.17144 
# Animal x Control: t = 1.30999    0.20085 
# Wind x Control: t =-0.09357    0.92612 

# By individuals.
Table3.8 <- aggregate(Table2.1$Frequency, by=list(Category=Table2.1$Plot, Table2.1$Life.History), FUN=sum, na.rm="TRUE")
colnames(Table3.8) <- c("Plot", "Life.History", "Recruits")
Table3.8$Treatment <- with(Basal, Treatment[match(Table3.8$Plot, Plot)]) # Match Treatment by Plot

Table3.8$Treatment <- relevel(Table3.8$Treatment, ref = "Wind") 
Table3.8$Treatment <- relevel(Table3.8$Treatment, ref = "Animal") 
Table3.8$Treatment <- relevel(Table3.8$Treatment, ref = "Forest") 

tgc2 <- summarySE(Table3.8, measurevar="Recruits", groupvars=c("Treatment","Life.History"))

LFAbundanceR <- ggplot(tgc2, aes(x=Treatment, y=Recruits, fill=Treatment)) + 
  geom_bar(stat="identity", color="black")+
  geom_errorbar(aes(ymin=Recruits-se, ymax=Recruits+se), width=.1) +
  xlab("Habitat") + 
  ylab(expression(Plant~Abundance)) + 
  ylim(0,260) + 
  scale_fill_manual(values=c("peachpuff4", "#9999CC","#66CC99", "#CC6666")) +
  facet_grid(~Life.History, scales = "free") +
  theme_bw(base_size = 20) +
  theme(legend.position="none") + 
  geom_text(position=pd, label = c("a", "", "", "b","c", "", "", "b"), aes(y = Recruits+se, x = Treatment),vjust = -0.5, size = 6) 

LFAbundanceR

Table5.1 <- droplevels(subset(Table3.8, Table3.8$Life.History != ""))
Table5.1LS <- droplevels(subset(Table5.1, Table5.1$Life.History == "Non-Pioneer"))
Table5.1LS$Treatment <- relevel(Table5.1LS$Treatment, ref = "Wind") 
lmx <- lm(data = Table5.1LS, Recruits ~ Treatment)
summary(lmx)

# Non-Pioneer, Individuals, Recruits Only

# F-statistic: 49.87761 on 3 and 28 DF,  p-value: 2.339947e-11
# Forest x Animal: t = -10.43500 3.7288e-11 ***
# Forest x Wind: t = -6.40256 6.2385e-07 ***
# Forest x Control: t = -10.70311 2.1032e-11 ***
# Animal x Wind: t = 4.03243 0.00038518 ***
# Animal x Control: t = -0.26811 0.79057696
# Wind x Control: t = -4.30055 0.00018682 ***

Table5.1P <- droplevels(subset(Table5.1, Table5.1$Life.History == "Pioneer"))

Table5.1P$Treatment <- relevel(Table5.1P$Treatment, ref = "Wind") 
Table5.1P$Treatment <- relevel(Table5.1P$Treatment, ref = "Animal") 
Table5.1P$Treatment <- relevel(Table5.1P$Treatment, ref = "Forest") 

lmx <- lm(data = Table5.1P, Recruits ~ Treatment)
summary(lmx)

# Pioneer, Individuals, Recruits Only

# F-statistic: 22.49976 on 3 and 28 DF,  p-value: 1.281247e-07
# Forest x Animal: t = 4.69534 6.3777e-05 ***
# Forest x Wind: t = 7.95361 1.1592e-08 ***
# Forest x Control: t = 5.75494 3.5459e-06 ***
# Animal x Wind: t = 3.25826  0.0029371 ** 
# Animal x Control: t = 1.05960  0.2983831  
# Wind x Control: t = -2.19867  0.0363302 *  

# By Species
NP <- droplevels(subset(tgc, tgc$Life.History == "Non-Pioneer"))
P <- droplevels(subset(tgc, tgc$Life.History == "Pioneer"))
NP$Ratio <- NP$Richness / P$Richness
# Forest 4.8478261
# Animal  0.6347826 
# Wind 0.5328467 
# Control 0.3309859

# By individuals
NP <- droplevels(subset(tgc2, tgc2$Life.History == "Non-Pioneer"))
P <- droplevels(subset(tgc2, tgc2$Life.History == "Pioneer"))
P$Ratio <- NP$Recruits / P$Recruits
# Forest 9.6056338 
# Animal 0.4642452  
# Wind 0.5292793 
# Control 0.3620178

# Diversity by Life History
BasalPLS1.1 <- table(BasalPLS1$Species, BasalPLS1$Plot) # Pioneer = 1
BasalPLS2.1 <- table(BasalPLS2$Species, BasalPLS2$Plot) # Non-Pioner = 2

BasalPLS1.2 <- diversity(t(BasalPLS1.1), "shannon") # Inverse-S  Diversity Index. Recruits Only
BasalPLS2.2 <- diversity(t(BasalPLS2.1), "shannon") # Inverse-S  Diversity Index. Recruits Only

# Shannon-Diversity Index. Only Recruits Pioneer 
BasalPLS1.3 <- as.data.frame(BasalPLS1.2)
BasalPLS1.3$Plot <- row.names(BasalPLS1.3)
colnames(BasalPLS1.3) <- c("Diversity", "Plot")
BasalPLS1.3$Treatment <- with(Basal, Treatment[match(BasalPLS1.3$Plot, Plot)]) # Match Treatment by Plot

BasalPLS1.3$Treatment <- relevel(BasalPLS1.3$Treatment, ref = "Wind") 
BasalPLS1.3$Treatment <- relevel(BasalPLS1.3$Treatment, ref = "Animal") 
BasalPLS1.3$Treatment <- relevel(BasalPLS1.3$Treatment, ref = "Forest") 

# Shannon-Diversity Index. Only Recruits Non-Pioneer 
BasalPLS2.3 <- as.data.frame(BasalPLS2.2)
BasalPLS2.3$Plot <- row.names(BasalPLS2.3)
colnames(BasalPLS2.3) <- c("Diversity", "Plot")
BasalPLS2.3$Treatment <- with(Basal, Treatment[match(BasalA2.6$Plot, Plot)]) # Match Treatment by Plot

BasalPLS2.3$Treatment <- relevel(BasalPLS2.3$Treatment, ref = "Wind") 
BasalPLS2.3$Treatment <- relevel(BasalPLS2.3$Treatment, ref = "Animal") 
BasalPLS2.3$Treatment <- relevel(BasalPLS2.3$Treatment, ref = "Forest") 

# Now rejoining Pioneer and Non-Pioneer
tgc1 <- summarySE(BasalPLS1.3, measurevar="Diversity", groupvars=c("Treatment"))
tgc2 <- summarySE(BasalPLS2.3, measurevar="Diversity", groupvars=c("Treatment"))

tgc1$Life.History <- "Pioneer"
tgc2$Life.History <- "Non-Pioneer"

tgc <- rbind(tgc1,tgc2)

LFDivRecruits <- ggplot(tgc, aes(x=Treatment, y=Diversity, fill=Treatment)) + 
  geom_bar(stat="identity", color="black")+
  geom_errorbar(aes(ymin=Diversity-se, ymax=Diversity+se), width=.1) +
  xlab("Habitat") + 
  ylab("Shannon Index (H)") + 
  scale_y_continuous(breaks=seq(0,3,0.5), limits = c(0,2.8)) +
  scale_fill_manual(values=c("peachpuff4", "#9999CC","#66CC99", "#CC6666")) +
  facet_grid(~Life.History, scales = "free") +
  theme_bw(base_size = 20) +
  theme(legend.position="none", axis.title.y = element_text(margin = margin(t = 0, r = 13, b = 0, l = 0)), axis.title.x = element_blank()) +
  geom_text(position=pd, label = c("a", "", "", "b","b", "", "", "a"), aes(y = Diversity+se, x = Treatment),vjust = -0.5, size = 6) 

LFDivRecruits

BasalPLS1.3$Treatment <- relevel(BasalPLS1.3$Treatment, ref = "Wind") 
lmx <- lm(data = BasalPLS1.3, Diversity ~ Treatment) # Pioneer Recruits
summary(lmx) 

# F-statistic: 12.9 on 3 and 28 DF,  p-value: 1.74e-05
# Forest x Animal: t = 4.96  3.1e-05 ***
# Forest x Control: t = 5.38  9.8e-06 ***
# Forest x Wind: t = 4.86  4.0e-05 ***
# Animal x Wind: t = -0.09     0.93  
# Animal x Control: t = 0.42     0.67 
# Wind x Control: t = 0.52     0.61    

BasalPLS2.3$Treatment <- relevel(BasalPLS2.3$Treatment, ref = "Wind") 
lmx <- lm(data = BasalPLS2.3, Diversity ~ Treatment) # Non-Pioneer Recruits
summary(lmx) 

# F-statistic: 12.1 on 3 and 28 DF,  p-value: 2.97e-05
# Forest x Animal: t = -4.00  0.00043 ***
# Forest x Control: t = -5.89  2.4e-06 ***
# Forest x Wind: t = -3.12  0.00421 ** 
# Animal x Wind: t = 0.88  0.38626 
# Animal x Control: t = -1.90  0.06793 . 
# Wind x Control: t = -2.78   0.0096 **   

##################################### LIFE HISTORY x TREATMENT (ALL) #################################### 
BasalPLS2 <- droplevels(subset(Basal, Basal$Life.History == "Pioneer" | Basal$Life.History == "Non-Pioneer"))
BasalPL1 <- droplevels(subset(BasalPLS2, BasalPLS2$Life.History == "Pioneer" ))
BasalPL2 <- droplevels(subset(BasalPLS2, BasalPLS2$Life.History == "Non-Pioneer"))

Table2.3 <- table(BasalPLS2$Species, BasalPLS2$Plot)
Table2.3 <- as.data.frame(Table2.3) # Includes all.

colnames(Table2.3) <- c("Species", "Plot", "Frequency")
Table2.3$Frequency2[Table2.3$Frequency > 0] <- 1
Table2.3[is.na(Table2.3)] <- 0
Table2.3$Life.History <- with(Basal, Life.History[match(Table2.3$Species, Species)]) # Match Life History by Species
Table2.4 <- aggregate(Table2.3$Frequency2, by=list(Category=Table2.3$Plot, Table2.3$Life.History), FUN=sum, na.rm="TRUE")
colnames(Table2.4) <- c("Plot", "Life.History", "Richness")
Table2.4$Treatment <- with(Basal, Treatment[match(Table2.4$Plot, Plot)]) # Match Treatment by Plot

Table2.4$Treatment <- relevel(Table2.4$Treatment, ref = "Wind") 
Table2.4$Treatment <- relevel(Table2.4$Treatment, ref = "Animal") 
Table2.4$Treatment <- relevel(Table2.4$Treatment, ref = "Forest") 

tgc <- summarySE(Table2.4, measurevar="Richness", groupvars=c("Treatment","Life.History"))

LFSpeciesAll <- ggplot(tgc, aes(x=Treatment, y=Richness, fill=Treatment)) + 
  geom_bar(stat="identity", color="black")+
  geom_errorbar(aes(ymin=Richness-se, ymax=Richness+se), width=.1) +
  xlab("Habitat") + 
  ylab("Species Density") + 
  ylim(0,33) +
  scale_fill_manual(values=c("peachpuff4", "#9999CC", "#66CC99", "#CC6666")) +
  geom_text(position=pd, label = c("a", "b", "b", "c","c", "b", "a", "b"), aes(y = Richness+se, x = Treatment),vjust = -0.5, size = 6) +
  facet_grid(~Life.History, scales = "free") +
  theme_bw(base_size = 20) +
  theme(legend.position="none", axis.title.y = element_text(margin = margin(t = 0, r = 13, b = 0, l = 0)), axis.title.x = element_blank())
LFSpeciesAll

Table2.4LS <- droplevels(subset(Table2.4, Table2.4$Life.History == "Non-Pioneer"))
Table2.4LS$Treatment <- relevel(Table2.4LS$Treatment, ref = "Wind") 
lmx <- lm(data = Table2.4LS, Richness ~ Treatment)
summary(lmx)

# Non-Pioneers, Richness, All
# F-statistic: 43.99361 on 3 and 28 DF,  p-value: 1.003242e-10

# Forest x Animal: t =  -6.92954 1.5613e-07 ***
# Forest x Wind: t = -6.28192 8.6000e-07 ***
# Forest x Control: t = -11.39813 4.9658e-12 ***
# Animal x Wind: t =  0.64762 0.52250526  
# Animal x Control: t = -4.455 0.000123 *
# Wind x Control: t = -5.11621 2.0185e-05 ***

Table2.4P <- droplevels(subset(Table2.4, Table2.4$Life.History == "Pioneer"))
Table2.4P$Treatment <- relevel(Table2.4P$Treatment, ref = "Wind") 
lmx <- lm(data = Table2.4P, Richness ~ Treatment)
summary(lmx)

# Pioneers, Richness, All
# F-statistic: 59.87365 on 3 and 28 DF,  p-value: 2.667357e-12

# Forest x Animal: t =  -12.42836 6.4833e-13 *** 
# Forest x Wind: t = 12.42836 6.4833e-13 *** 
# Forest x Control: t =   9.86052 1.3112e-10 *** 
# Animal x Wind: t = 2.77327  0.0097628 ** 
# Animal x Control: t = 0.20543  0.8387244     
# Wind x Control: t = -2.56784  0.0158606 *

# Life History, By individuals, All
Table4.7 <- aggregate(Table2.3$Frequency, by=list(Category=Table2.3$Plot, Table2.3$Life.History), FUN=sum, na.rm="TRUE")
colnames(Table4.7) <- c("Plot", "Life.History", "Frequency")
Table4.7$Treatment <- with(Basal, Treatment[match(Table4.7$Plot, Plot)]) # Match Treatment by Plot

Table4.7$Treatment <- relevel(Table4.7$Treatment, ref = "Wind") 
Table4.7$Treatment <- relevel(Table4.7$Treatment, ref = "Animal") 
Table4.7$Treatment <- relevel(Table4.7$Treatment, ref = "Forest") 

tgc <- summarySE(Table4.7, measurevar="Frequency", groupvars=c("Treatment","Life.History"))

LFAllAbundance <- ggplot(tgc, aes(x=Treatment, y=Frequency, fill=Treatment)) + 
  geom_bar(stat="identity", color="black")+
  geom_errorbar(aes(ymin=Frequency-se, ymax=Frequency+se), width=.1) +
  xlab("Habitat") + 
  ylab(expression(Plant~Abundance)) + 
  ylim(0,260) + 
  scale_fill_manual(values=c("peachpuff4", "#9999CC","#66CC99", "#CC6666")) +
  geom_text(position=pd, label = c("a", "b·", "a·", "b","a", "b", "c", "b"), aes(y = Frequency+se, x = Treatment),vjust = -0.5, size = 6) +
  facet_grid(~Life.History, scales = "free") +
  theme_bw(base_size = 20) +
  theme(legend.position="none")
LFAllAbundance

Table4.8<- aggregate(Table2.3$Frequency, by=list(Category=Table2.3$Plot, Table2.3$Species), FUN=sum, na.rm="TRUE")
colnames(Table4.8) <- c("Plot", "Species", "Recruits")
Table4.8$Treatment <- with(Basal, Treatment[match(Table4.8$Plot, Plot)]) # Match Treatment by Plot

Table4.8$Treatment <- relevel(Table4.8$Treatment, ref = "Wind") 
Table4.8$Treatment <- relevel(Table4.8$Treatment, ref = "Animal") 
Table4.8$Treatment <- relevel(Table4.8$Treatment, ref = "Forest") 

Table4.8$Life.History <- with(Basal, Life.History[match(Table4.8$Species, Species)]) # Match Treatment by Plot

Table4.21 <- aggregate(Table4.8$Recruits, by=list(Category=Table4.8$Plot, Table4.8$Life.History), FUN=sum, na.rm="TRUE")

colnames(Table4.21) <- c("Plot", "Life.History", "Recruits")
Table4.21$Treatment <- with(Basal, Treatment[match(Table4.21$Plot, Plot)]) # Match Treatment by Plot

Table4.21$Treatment <- relevel(Table4.21$Treatment, ref = "Wind") 
Table4.21$Treatment <- relevel(Table4.21$Treatment, ref = "Animal") 
Table4.21$Treatment <- relevel(Table4.21$Treatment, ref = "Forest") 

tgc <- summarySE(Table4.21, measurevar="Recruits", groupvars=c("Treatment","Life.History")) 

LFAllAbundance <- ggplot(tgc, aes(x=Treatment, y=Recruits, fill=Treatment)) + 
  geom_bar(stat="identity", color="black")+
  geom_errorbar(aes(ymin=Recruits-se, ymax=Recruits+se), width=.1) +
  xlab("Habitat") + 
  ylab("Plant Abundance") + 
  ylim(0,260) + 
  scale_fill_manual(values=c("peachpuff4", "#9999CC","#66CC99", "#CC6666")) +
  facet_grid(~Life.History, scales = "free") +
  geom_text(position=pd, label = c("a", "b·", "a·", "b","c", "b", "a", "b"), aes(y = Recruits+se, x = Treatment),vjust = -0.5, size = 6) +
  theme_bw(base_size = 20) +
  theme(legend.position="none")

LFAllAbundance

Table4.9LS <- droplevels(subset(Table4.8, Table4.8$Life.History == "Non-Pioneer")) 
Table4.9LS$Treatment <- relevel(Table4.9LS$Treatment, ref = "Wind") 
lmx <- lm(data = Table4.9LS, Recruits ~ Treatment)
summary(lmx)

# Non-Pioneers, Individuals, All Plants

# F-statistic:  7.81567 on 3 and 2972 DF,  p-value: 3.399092e-05
# Forest x Animal: t = -3.27821  0.0010567 **  
# Forest x Wind: t = -1.36295  0.1730015  
# Forest x Control: t = -4.44518 9.1034e-06 ***
# Animal x Wind: t =  1.91526  0.0555557 . 
# Animal x Control: t = -1.16697  0.2433160
# Wind x Control: t = -3.08223  0.0020734 ** 

Table4.9P <- droplevels(subset(Table4.8, Table4.8$Life.History == "Pioneer")) 
Table4.9P$Treatment <- relevel(Table4.9P$Treatment, ref = "Wind") 
lmx <- lm(data = Table4.9P, Recruits ~ Treatment)
summary(lmx)

# Pioneers, Individuals, All plants

# F-statistic: 19.44634 on 3 and 1404 DF,  p-value: 2.352505e-12

# Forest x Animal: t =  5.01840 5.8767e-07 ***
# Forest x Wind: t = 7.47173 1.3856e-13 ***
# Forest x Control: t = 4.62459 4.0988e-06 ***
# Animal x Wind: t = 2 2.45332   0.014275 * 
# Animal x Control: t = -0.39381   0.693782 
# Wind x Control: t = -2.84713  0.0044757 **

pushViewport(viewport(layout=grid.layout(nrow=2,ncol=1))) # This sets the multiplot.
print(LFSpeciesAll, vp=viewport(layout.pos.row=1, layout.pos.col=1)) # 1000 x 800
print(LFAllAbundance, vp=viewport(layout.pos.row=2, layout.pos.col=1))

pushViewport(viewport(layout=grid.layout(nrow=2,ncol=1))) # This sets the multiplot.
print(LFSpeciesR, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(LFAbundanceR, vp=viewport(layout.pos.row=2, layout.pos.col=1))

# Diversity by Life History
BasalPLS1.4 <- table(BasalPL1$Species, BasalPL1$Plot) # Pioneer = 1
BasalPLS2.4 <- table(BasalPL2$Species, BasalPL2$Plot) # Non-Pioneer = 2

BasalPLS1.5 <- diversity(t(BasalPLS1.4), "shannon") # Inverse-S  Diversity Index. 
BasalPLS2.5 <- diversity(t(BasalPLS2.4), "shannon") # Inverse-S  Diversity Index. 

# Shannon-Diversity Index. Recruits Pioneer 
BasalPLS1.6 <- as.data.frame(BasalPLS1.5)
BasalPLS1.6$Plot <- row.names(BasalPLS1.6)
colnames(BasalPLS1.6) <- c("Diversity", "Plot")
BasalPLS1.6$Treatment <- with(Basal, Treatment[match(BasalPLS1.6$Plot, Plot)]) # Match Treatment by Plot

BasalPLS1.6$Treatment <- relevel(BasalPLS1.6$Treatment, ref = "Wind") 
BasalPLS1.6$Treatment <- relevel(BasalPLS1.6$Treatment, ref = "Animal") 
BasalPLS1.6$Treatment <- relevel(BasalPLS1.6$Treatment, ref = "Forest") 

# Shannon-Diversity Index. Recruits NonPioneer 
BasalPLS2.6 <- as.data.frame(BasalPLS2.5)
BasalPLS2.6$Plot <- row.names(BasalPLS2.6)
colnames(BasalPLS2.6) <- c("Diversity", "Plot")
BasalPLS2.6$Treatment <- with(Basal, Treatment[match(BasalPLS2.6$Plot, Plot)]) # Match Treatment by Plot

BasalPLS2.6$Treatment <- relevel(BasalPLS2.6$Treatment, ref = "Wind") 
BasalPLS2.6$Treatment <- relevel(BasalPLS2.6$Treatment, ref = "Animal") 
BasalPLS2.6$Treatment <- relevel(BasalPLS2.6$Treatment, ref = "Forest") 

# Now rejoining Pioneer and NonPioneer
tgc1 <- summarySE(BasalPLS1.6, measurevar="Diversity", groupvars=c("Treatment"))
tgc2 <- summarySE(BasalPLS2.6, measurevar="Diversity", groupvars=c("Treatment"))

tgc1$Life.History <- "Pioneer"
tgc2$Life.History <- "Non-Pioneer"

tgc <- rbind(tgc1,tgc2)

LFDivAll <- ggplot(tgc, aes(x=Treatment, y=Diversity, fill=Treatment)) + 
  geom_bar(stat="identity", color="black")+
  geom_errorbar(aes(ymin=Diversity-se, ymax=Diversity+se), width=.1) +
  xlab("Habitat") + 
  ylab("Shannon Index (H)") + 
  scale_y_continuous(breaks=seq(0,3,0.5), limits = c(0,2.8)) +
  scale_fill_manual(values=c("peachpuff4", "#9999CC","#66CC99", "#CC6666")) +
  facet_grid(~Life.History, scales = "free") +
  theme_bw(base_size = 20) +
  theme(legend.position="none", axis.title.y = element_text(margin = margin(t = 0, r = 13, b = 0, l = 0)), axis.title.x = element_blank()) + 
  geom_text(position=pd, label = c("a", "a", "a", "b","b", "a", "a", "a"), aes(y = Diversity+se, x = Treatment),vjust = -0.5, size = 6) 

LFDivAll

BasalPLS1.6$Treatment <- relevel(BasalPLS1.6$Treatment, ref = "Wind") 
lmx <- lm(data = BasalPLS1.6, Diversity ~ Treatment) # Pioneer All
summary(lmx) 

# F-statistic: 20.7 on 3 and 28 DF,  p-value: 2.84e-07
# Forest x Animal: t = 6.71  2.8e-07 ***
# Forest x Control: t = 5.80  3.1e-06 ***
# Forest x Wind: t = 6.65  3.3e-07 ***
# Animal x Wind: t = -0.06     0.95
# Animal x Control: t = -0.91     0.37
# Wind x Control: t = -0.85     0.40   

BasalPLS2.6$Treatment <- relevel(BasalPLS2.6$Treatment, ref = "Wind") 
lmx <- lm(data = BasalPLS2.6, Diversity ~ Treatment) # Non-Pioneer All
summary(lmx) 

# F-statistic: 15.4 on 3 and 28 DF,  p-value: 4.18e-06

# Forest x Animal: t = -0.48     0.64 
# Forest x Control: t = -5.85  2.7e-06 ***
# Forest x Wind: t = -0.49     0.62  
# Animal x Wind: t = -0.02     0.99
# Animal x Control: t = -5.37  1.0e-05 ***
# Wind x Control: t = -5.36  1.0e-05 *** 

# Covariate Model for Life History

RecTable2L <- Table3.9LS
RecTable2P <- Table3.9P

RecTable2L$Recruits <- with(Table5.1LS, Recruits[match(RecTable2L$Plot, Plot)]) 
RecTable2P$Recruits <- with(Table5.1P, Recruits[match(RecTable2P$Plot, Plot)]) 

RecTable2L$Treatment <- relevel(RecTable2L$Treatment, ref = "Wind") 
lmx <- lm(data = RecTable2L, Richness ~ Treatment + Recruits)
summary(lmx) # Late-Succ

RecTable2P$Treatment <- relevel(RecTable2P$Treatment, ref = "Wind") 
lmx <- lm(data = RecTable2P, Richness ~ Treatment + Recruits)
summary(lmx) # Pioneer

sum(Table2.2[,3]) # Recruited Species
sum(Table3.8[,3]) # Recruited Abundance

Table3.7

# Covariate model all life history

AllTable2L <- Table2.4LS
AllTable2P <- Table2.4P

AllTable2L$Recruits <- with(Table4.9LS, Recruits[match(AllTable2L$Plot, Plot)]) 
AllTable2P$Recruits <- with(Table4.9P, Recruits[match(AllTable2P$Plot, Plot)]) 

AllTable2L$Treatment <- relevel(AllTable2L$Treatment, ref = "Wind") 
lmx <- lm(data = AllTable2L, Richness ~ Treatment + Recruits)
summary(lmx) # Late-Succ

AllTable2P$Treatment <- relevel(AllTable2P$Treatment, ref = "Wind") 
lmx <- lm(data = AllTable2P, Richness ~ Treatment + Recruits)
summary(lmx) # Pioneer

sum(Table2.4[,3]) # All Species
sum(Table4.8[,3]) # All Abundance


Table4.9P

##################################### Simpson's Diversity Index #################################### 
# Shannon-Weiner is embedded in previous sections. 
Tablax1 <- diversity(t(Table1.3), "invsimpson") # Simpsons Diversity Index. All
Tablax2 <- diversity(t(Table1.6), "invsimpson") # Simpsons Diversity Index. Recruits Only

# Using Simpsons All
Tablax1.1 <- as.data.frame(Tablax1)
Tablax1.1$Plot <- row.names(Tablax1.1)
colnames(Tablax1.1) <- c("Diversity", "Plot")
Tablax1.1$Treatment <- with(Basal, Treatment[match(Tablax1.1$Plot, Plot)]) # Match Treatment by Plot

Tablax1.1$Treatment <- relevel(Tablax1.1$Treatment, ref = "Wind") 
Tablax1.1$Treatment <- relevel(Tablax1.1$Treatment, ref = "Animal") 
Tablax1.1$Treatment <- relevel(Tablax1.1$Treatment, ref = "Forest") 

tgc <- summarySE(Tablax1.1, measurevar="Diversity", groupvars=c("Treatment"))

SupDivdAll2 <- ggplot(tgc, aes(x=Treatment, y=Diversity, fill=Treatment)) + 
  geom_bar(stat="identity", color="black")+
  geom_errorbar(aes(ymin=Diversity-se, ymax=Diversity+se), width=.1) +
  xlab("") + 
  ylab("Simpson's Index (1/D)") + 
  scale_y_continuous(breaks=seq(0,17,3), limits = c(0,17)) +
  scale_fill_manual(values=c("peachpuff4","#9999CC","#66CC99", "#CC6666")) +
  geom_text(position=pd, label = c("b", "a", "a", "b"), aes(y = Diversity+se, x = Treatment),vjust = -0.5, size = 6) +
  theme_bw(base_size = 20) +
  theme(legend.position="none", axis.title.y = element_text(margin = margin(t = 0, r = 13, b = 0, l = 0))) 
SupDivdAll2

# Using Simpson's Recruits Only

Tablax2.1 <- as.data.frame(Tablax2)
Tablax2.1$Plot <- row.names(Tablax2.1)
colnames(Tablax2.1) <- c("Diversity", "Plot")
Tablax2.1$Treatment <- with(Basal, Treatment[match(Tablax2.1$Plot, Plot)]) # Match Treatment by Plot

Tablax2.1$Treatment <- relevel(Tablax2.1$Treatment, ref = "Wind") 
Tablax2.1$Treatment <- relevel(Tablax2.1$Treatment, ref = "Animal") 
Tablax2.1$Treatment <- relevel(Tablax2.1$Treatment, ref = "Forest") 

tgc <- summarySE(Tablax2.1, measurevar="Diversity", groupvars=c("Treatment"))

SupDivdR2 <- ggplot(tgc, aes(x=Treatment, y=Diversity, fill=Treatment)) + 
  geom_bar(stat="identity", color="black")+
  geom_errorbar(aes(ymin=Diversity-se, ymax=Diversity+se), width=.1) +
  xlab("") + 
  ylab("Simpson's Index (1/D)") + 
  scale_y_continuous(breaks=seq(0,17,3), limits = c(0,17)) +
  scale_fill_manual(values=c("peachpuff4","#9999CC","#66CC99", "#CC6666")) +
  geom_text(position=pd, label = c("b", "", "", "b"), aes(y = Diversity+se, x = Treatment),vjust = -0.5, size = 6) +
  theme_bw(base_size = 20) +
  theme(legend.position="none", axis.title.y = element_text(margin = margin(t = 0, r = 13, b = 0, l = 0))) 
SupDivdR2

Tablax1.1$Treatment <- relevel(Tablax1.1$Treatment, ref = "Wind") 
lmx <- lm(data = Tablax1.1, Diversity ~ Treatment)
summary(lmx)

# Diversity All
# F-statistic: 13.3 on 3 and 28 DF,  p-value: 1.42e-05

# Forest x Animal: t = 5.51  6.9e-06 ***
# Forest x Control: t = 1.65     0.11 
# Forest x Wind: t = 4.67  6.9e-05 ***
# Animal x Wind: t = -0.84  0.40642  
# Animal x Control: t =  -3.86  0.00061 ***
# Wind x Control: t = -3.02   0.0054 ** 

Tablax2.1$Treatment <- relevel(Tablax2.1$Treatment, ref = "Animal") 
lmx <- lm(data = Tablax2.1, Diversity ~ Treatment)
summary(lmx)

# Diversity only Recruits 
# F-statistic: 2.67 on 3 and 28 DF,  p-value: 0.0668
# Forest x Animal: t = 2.34    0.027 *
# Forest x Control:t = 1.92    0.066 .  
# Forest x Wind: t = 2.51    0.018 *
# Animal x Wind: t = 00.17    0.864  
# Animal x Control: t = -0.42    0.676   
# Wind x Control: t = 0.17    0.864 

# DISPERSAL MODE x TREATMENT (Recruits) # 
BasalAW1.2X <- diversity(t(BasalAW1.1), "invsimpson") # Inverse-S  Diversity Index. Recruits Only
BasalAW2.2X <- diversity(t(BasalAW2.1), "invsimpson") # Inverse-S  Diversity Index. Recruits Only

# invsimpson-Diversity Index. Only Recruits Biotic. 
BasalAW1.3X <- as.data.frame(BasalAW1.2X)
BasalAW1.3X$Plot <- row.names(BasalAW1.3X)
colnames(BasalAW1.3X) <- c("Diversity", "Plot")
BasalAW1.3X$Treatment <- with(Basal, Treatment[match(BasalAW1.3X$Plot, Plot)]) # Match Treatment by Plot

BasalAW1.3X$Treatment <- relevel(BasalAW1.3X$Treatment, ref = "Wind") 
BasalAW1.3X$Treatment <- relevel(BasalAW1.3X$Treatment, ref = "Animal") 
BasalAW1.3X$Treatment <- relevel(BasalAW1.3X$Treatment, ref = "Forest") 

# invsimpson-Diversity Index. Only Recruits Abiotic. 
BasalAW2.3X <- as.data.frame(BasalAW2.2X)
BasalAW2.3X$Plot <- row.names(BasalAW2.3X)
colnames(BasalAW2.3X) <- c("Diversity", "Plot")
BasalAW2.3X$Treatment <- with(Basal, Treatment[match(BasalAW2.3X$Plot, Plot)]) # Match Treatment by Plot

BasalAW2.3X$Treatment <- relevel(BasalAW2.3X$Treatment, ref = "Wind") 
BasalAW2.3X$Treatment <- relevel(BasalAW2.3X$Treatment, ref = "Animal") 
BasalAW2.3X$Treatment <- relevel(BasalAW2.3X$Treatment, ref = "Forest") 

# Now rejoining Biotic and Abiotic

tgc1 <- summarySE(BasalAW1.3X, measurevar="Diversity", groupvars=c("Treatment"))
tgc2 <- summarySE(BasalAW2.3X, measurevar="Diversity", groupvars=c("Treatment"))

tgc1$Dispersal.Mode <- "Biotic"
tgc2$Dispersal.Mode <- "Abiotic"

tgc <- rbind(tgc1,tgc2)

tgc$Dispersal.Mode <- revalue(tgc$Dispersal.Mode, c("Abiotic"="Abiotic Dispersal", "Biotic"="Biotic Dispersal"))

DPDivRecruits2 <- ggplot(tgc, aes(x=Treatment, y=Diversity, fill=Treatment)) + 
  geom_bar(stat="identity", color="black")+
  geom_errorbar(aes(ymin=Diversity-se, ymax=Diversity+se), width=.1) +
  xlab("Habitat") + 
  ylab("Simpson's Index (1/D)") + 
  scale_y_continuous(breaks=seq(0,15,3), limits = c(0,15)) +
  scale_fill_manual(values=c("peachpuff4", "#9999CC","#66CC99", "#CC6666")) +
  facet_grid(~Dispersal.Mode, scales = "free") +
  theme_bw(base_size = 20) +
  theme(legend.position="none", axis.title.y = element_text(margin = margin(t = 0, r = 13, b = 0, l = 0)), axis.title.x = element_blank()) +
  geom_text(position=pd, label = c("b", "c", "", "c","b", "", "b", "b"), aes(y = Diversity+se, x = Treatment),vjust = -0.5, size = 6) 

DPDivRecruits2

BasalAW1.3X$Treatment <- relevel(BasalAW1.3X$Treatment, ref = "Wind") 
lmx <- lm(data = BasalAW1.3X, Diversity ~ Treatment) # Recruits Biotic
summary(lmx)
# F-statistic: 1.15 on 3 and 28 DF,  p-value: 0.345

# Diversity, Biotic, Recruits Only
# Forest x Animal: t = 1.64     0.11 
# Forest x Wind: t = 1.58     0.13 
# Forest x Control: t = 0.97     0.34
# Animal x Wind: t = -0.06     0.95  
# Animal x Control: t = -0.67     0.51
# Wind x Control: t = -0.60     0.55 

BasalAW2.3X$Treatment <- relevel(BasalAW2.3X$Treatment, ref = "Forest") 
lmx <- lm(data = BasalAW2.3X, Diversity ~ Treatment) # Recruits Abiotic
summary(lmx)
# -statistic:  5.3 on 3 and 28 DF,  p-value: 0.00506

# Diversity, Abiotic, Recruits Only
# Forest x Animal: t = 2.96  0.00614 ** 
# Forest x Wind: t = 43.72  0.00088 ***
# Forest x Control: t = 2.77  0.00983 ** 
# Animal x Wind: t = 0.32     0.75 
# Animal x Control: t = -0.03     0.97 
# Wind x Control: t = -0.35     0.73    

# DISPERSAL MODE x TREATMENT (ALL) 
BasalA1.5X <- diversity(t(BasalA1.4), "invsimpson") # Inverse-S  Diversity Index. Recruits Only
BasalA2.5X <- diversity(t(BasalA2.4), "invsimpson") # Inverse-S  Diversity Index. Recruits Only

# invsimpson-Diversity Index. Only Recruits Biotic. 
BasalA1.6X <- as.data.frame(BasalA1.5X)
BasalA1.6X$Plot <- row.names(BasalA1.6X)
colnames(BasalA1.6X) <- c("Diversity", "Plot")
BasalA1.6X$Treatment <- with(Basal, Treatment[match(BasalA1.6X$Plot, Plot)]) # Match Treatment by Plot

BasalA1.6X$Treatment <- relevel(BasalA1.6X$Treatment, ref = "Wind") 
BasalA1.6X$Treatment <- relevel(BasalA1.6X$Treatment, ref = "Animal") 
BasalA1.6X$Treatment <- relevel(BasalA1.6X$Treatment, ref = "Forest") 

# invsimpson-Diversity Index. Only Recruits Abiotic. 
BasalA2.6X <- as.data.frame(BasalA2.5X)
BasalA2.6X$Plot <- row.names(BasalA2.6X)
colnames(BasalA2.6X) <- c("Diversity", "Plot")
BasalA2.6X$Treatment <- with(Basal, Treatment[match(BasalA2.6X$Plot, Plot)]) # Match Treatment by Plot

BasalA2.6X$Treatment <- relevel(BasalA2.6X$Treatment, ref = "Wind") 
BasalA2.6X$Treatment <- relevel(BasalA2.6X$Treatment, ref = "Animal") 
BasalA2.6X$Treatment <- relevel(BasalA2.6X$Treatment, ref = "Forest") 

# Now rejoining Biotic and Abiotic

tgc1 <- summarySE(BasalA1.6X, measurevar="Diversity", groupvars=c("Treatment"))
tgc2 <- summarySE(BasalA2.6X, measurevar="Diversity", groupvars=c("Treatment"))

tgc1$Dispersal.Mode <- "Biotic"
tgc2$Dispersal.Mode <- "Abiotic"

tgc <- rbind(tgc1,tgc2)

tgc$Dispersal.Mode <- revalue(tgc$Dispersal.Mode, c("Abiotic"="Abiotic Dispersal", "Biotic"="Biotic Dispersal"))

DPDivAll2 <- ggplot(tgc, aes(x=Treatment, y=Diversity, fill=Treatment)) + 
  geom_bar(stat="identity", color="black")+
  geom_errorbar(aes(ymin=Diversity-se, ymax=Diversity+se), width=.1) +
  xlab("Habitat") + 
  ylab("Simpson's Index (1/D)") + 
  scale_y_continuous(breaks=seq(0,15,3), limits = c(0,15)) +
  scale_fill_manual(values=c("peachpuff4", "#9999CC","#66CC99", "#CC6666")) +
  facet_grid(~Dispersal.Mode, scales = "free") +
  theme_bw(base_size = 20) +
  theme(legend.position="none", axis.title.y = element_text(margin = margin(t = 0, r = 13, b = 0, l = 0)), axis.title.x = element_blank()) + 
  geom_text(position=pd, label = c("b", "c", "a", "c","b", "a", "b", "b"), aes(y = Diversity+se, x = Treatment),vjust = -0.5, size = 6) 
DPDivAll2

BasalA1.6X$Treatment <- relevel(BasalA1.6X$Treatment, ref = "Wind") 
lmx <- lm(data = BasalA1.6X, Diversity ~ Treatment) # All Biotic
summary(lmx)

# F-statistic: 12.7 on 3 and 28 DF,  p-value: 1.99e-05
# Diversity, Biotic, All 
# Forest x Animal: t = 5.69  4.3e-06 ***
# Forest x Wind: t =  1.47     0.15  
# Forest x Control: t =  0.90     0.37  
# Animal x Wind: t = -4.22  0.00023 ***
# Animal x Control: t = -4.78  5.0e-05 ***
# Wind x Control: t = -0.56  0.57966    

BasalA2.6X$Treatment <- relevel(BasalA2.6X$Treatment, ref = "Wind") 
lmx <- lm(data = BasalA2.6X, Diversity ~ Treatment) # Recruits Abiotic
summary(lmx)
# F-statistic: 45.5 on 3 and 28 DF,  p-value: 6.87e-11

# Diversity, Abiotic, All Only
# Forest x Animal: t = 2.20   0.0360 * 
# Forest x Wind: t = 10.74  1.9e-11 ***
# Forest x Control: t =  2.06   0.0489 *
# Animal x Wind: t = 8.54  2.8e-09 ***
# Animal x Control: t = -0.14    0.887 
# Wind x Control: t = -8.68  2.0e-09 ***

# LIFE HISTORY x TREATMENT (Recruits Only) 
BasalPLS1.2X <- diversity(t(BasalPLS1.1), "invsimpson") # Inverse-S  Diversity Index. Recruits Only
BasalPLS2.2X <- diversity(t(BasalPLS2.1), "invsimpson") # Inverse-S  Diversity Index. Recruits Only

# invsimpson-Diversity Index. Only Recruits Pioneer 
BasalPLS1.3X <- as.data.frame(BasalPLS1.2X)
BasalPLS1.3X$Plot <- row.names(BasalPLS1.3X)
colnames(BasalPLS1.3X) <- c("Diversity", "Plot")
BasalPLS1.3X$Treatment <- with(Basal, Treatment[match(BasalPLS1.3X$Plot, Plot)]) # Match Treatment by Plot

BasalPLS1.3X$Treatment <- relevel(BasalPLS1.3X$Treatment, ref = "Wind") 
BasalPLS1.3X$Treatment <- relevel(BasalPLS1.3X$Treatment, ref = "Animal") 
BasalPLS1.3X$Treatment <- relevel(BasalPLS1.3X$Treatment, ref = "Forest") 

# invsimpson-Diversity Index. Only Recruits Non-Pioneer 
BasalPLS2.3X <- as.data.frame(BasalPLS2.2X)
BasalPLS2.3X$Plot <- row.names(BasalPLS2.3X)
colnames(BasalPLS2.3X) <- c("Diversity", "Plot")
BasalPLS2.3X$Treatment <- with(Basal, Treatment[match(BasalPLS2.3X$Plot, Plot)]) # Match Treatment by Plot

BasalPLS2.3X$Treatment <- relevel(BasalPLS2.3X$Treatment, ref = "Wind") 
BasalPLS2.3X$Treatment <- relevel(BasalPLS2.3X$Treatment, ref = "Animal") 
BasalPLS2.3X$Treatment <- relevel(BasalPLS2.3X$Treatment, ref = "Forest") 

# Now rejoining Pioneer and Non-Pioneer
tgc1 <- summarySE(BasalPLS1.3X, measurevar="Diversity", groupvars=c("Treatment"))
tgc2 <- summarySE(BasalPLS2.3X, measurevar="Diversity", groupvars=c("Treatment"))

tgc1$Life.History <- "Pioneer"
tgc2$Life.History <- "Non-Pioneer"

tgc <- rbind(tgc1,tgc2)

LFDivRecruits2 <- ggplot(tgc, aes(x=Treatment, y=Diversity, fill=Treatment)) + 
  geom_bar(stat="identity", color="black")+
  geom_errorbar(aes(ymin=Diversity-se, ymax=Diversity+se), width=.1) +
  xlab("Habitat") + 
  ylab("Simpson's Index (1/D)") + 
  scale_y_continuous(breaks=seq(0,11,2), limits = c(0,11)) +
  scale_fill_manual(values=c("peachpuff4", "#9999CC","#66CC99", "#CC6666")) +
  facet_grid(~Life.History, scales = "free") +
  theme_bw(base_size = 20) +
  theme(legend.position="none", axis.title.y = element_text(margin = margin(t = 0, r = 13, b = 0, l = 0)), axis.title.x = element_blank()) +
  geom_text(position=pd, label = c("a", "", "", "b","b", "", "", "c"), aes(y = Diversity+se, x = Treatment),vjust = -0.5, size = 6) 

LFDivRecruits2

BasalPLS1.3X$Treatment <- relevel(BasalPLS1.3X$Treatment, ref = "Wind") 
lmx <- lm(data = BasalPLS1.3X, Diversity ~ Treatment) # Pioneer Recruits
summary(lmx) 

# F-statistic: 9.02 on 3 and 28 DF,  p-value: 0.000243
# Forest x Animal: t =  4.30  0.00019 ***
# Forest x Wind: t = 3.46  0.00175 ** 
# Forest x Control: t = 4.63  7.5e-05 ***
# Animal x Wind: t = -0.84  0.40784   
# Animal x Control: t =  0.33  0.74023 
# Wind x Control: t = 1.18   0.2498  

BasalPLS2.3X$Treatment <- relevel(BasalPLS2.3X$Treatment, ref = "Forest") 
lmx <- lm(data = BasalPLS2.3X, Diversity ~ Treatment) # Non-Pioneer Recruits
summary(lmx) 

# F-statistic: 3.91 on 3 and 28 DF,  p-value: 0.0189

# Forest x Animal: t = -2.17   0.0383 * 
# Forest x Wind: t = -1.78   0.0857 . 
# Forest x Control: t = -3.38   0.0022 ** 
# Animal x Wind: t =  0.39    0.697  
# Animal x Control: t =  0.39    0.697 
# Wind x Control: t = -1.60    0.122      

# LIFE HISTORY x TREATMENT (ALL) 
BasalPLS1.5X <- diversity(t(BasalPLS1.4), "invsimpson") # Inverse-S  Diversity Index. ALL
BasalPLS2.5X <- diversity(t(BasalPLS2.4), "invsimpson") # Inverse-S  Diversity Index. ALL

# invsimpson-Diversity Index. Recruits Pioneer 
BasalPLS1.6X <- as.data.frame(BasalPLS1.5X)
BasalPLS1.6X$Plot <- row.names(BasalPLS1.6X)
colnames(BasalPLS1.6X) <- c("Diversity", "Plot")
BasalPLS1.6X$Treatment <- with(Basal, Treatment[match(BasalPLS1.6X$Plot, Plot)]) # Match Treatment by Plot

BasalPLS1.6X$Treatment <- relevel(BasalPLS1.6X$Treatment, ref = "Wind") 
BasalPLS1.6X$Treatment <- relevel(BasalPLS1.6X$Treatment, ref = "Animal") 
BasalPLS1.6X$Treatment <- relevel(BasalPLS1.6X$Treatment, ref = "Forest") 

# invsimpson-Diversity Index. Recruits NonPioneer 
BasalPLS2.6X <- as.data.frame(BasalPLS2.5X)
BasalPLS2.6X$Plot <- row.names(BasalPLS2.6X)
colnames(BasalPLS2.6X) <- c("Diversity", "Plot")
BasalPLS2.6X$Treatment <- with(Basal, Treatment[match(BasalPLS2.6X$Plot, Plot)]) # Match Treatment by Plot

BasalPLS2.6X$Treatment <- relevel(BasalPLS2.6X$Treatment, ref = "Wind") 
BasalPLS2.6X$Treatment <- relevel(BasalPLS2.6X$Treatment, ref = "Animal") 
BasalPLS2.6X$Treatment <- relevel(BasalPLS2.6X$Treatment, ref = "Forest") 

# Now rejoining Pioneer and NonPioneer
tgc1 <- summarySE(BasalPLS1.6X, measurevar="Diversity", groupvars=c("Treatment"))
tgc2 <- summarySE(BasalPLS2.6X, measurevar="Diversity", groupvars=c("Treatment"))

tgc1$Life.History <- "Pioneer"
tgc2$Life.History <- "Non-Pioneer"

tgc <- rbind(tgc1,tgc2)

LFDivAll2 <- ggplot(tgc, aes(x=Treatment, y=Diversity, fill=Treatment)) + 
  geom_bar(stat="identity", color="black")+
  geom_errorbar(aes(ymin=Diversity-se, ymax=Diversity+se), width=.1) +
  xlab("Habitat") + 
  ylab("Simpson's Index (1/D)") + 
  scale_y_continuous(breaks=seq(0,11,2), limits = c(0,11)) +
  scale_fill_manual(values=c("peachpuff4", "#9999CC","#66CC99", "#CC6666")) +
  facet_grid(~Life.History, scales = "free") +
  theme_bw(base_size = 20) +
  theme(legend.position="none", axis.title.y = element_text(margin = margin(t = 0, r = 13, b = 0, l = 0)), axis.title.x = element_blank()) +
  geom_text(position=pd, label = c("a", "a", "a", "b","b", "a", "ac", "c"), aes(y = Diversity+se, x = Treatment),vjust = -0.5, size = 6) 

LFDivAll2

BasalPLS1.6X$Treatment <- relevel(BasalPLS1.6X$Treatment, ref = "Forest") 
lmx <- lm(data = BasalPLS1.6X, Diversity ~ Treatment) # Pioneer All
summary(lmx) 

# F-statistic: 16.8 on 3 and 28 DF,  p-value: 1.93e-06

# Forest x Animal: t = 6.76  2.4e-07 ***
# Forest x Wind: t = 5.22  1.5e-05 ***
# Forest x Control: t = 4.40  0.00014 ***
# Animal x Wind: t = -1.55    0.134
# Animal x Control: t = -2.36    0.025 *
# Wind x Control: t = -0.82     0.42  

BasalPLS2.6X$Treatment <- relevel(BasalPLS2.6X$Treatment, ref = "Forest") 
lmx <- lm(data = BasalPLS2.6X, Diversity ~ Treatment) # Non-Pioneer All
summary(lmx) 

# F-statistic: 5.43 on 3 and 28 DF,  p-value: 0.00453

# Forest x Animal: t = 1.06    0.298 
# Forest x Wind: t = 0.53    0.603   
# Forest x Control: t = -2.65    0.013 *
# Animal x Wind: t = -0.54  0.59670
# Animal x Control: t = -3.71  0.00091 ***
# Wind x Control: t = -3.18   0.0036 ** 

################################## MEGA FIGURES ################################# 

# Life History 
pushViewport(viewport(layout=grid.layout(nrow=3,ncol=1))) # This sets the multiplot.
print(LFSpeciesAll, vp=viewport(layout.pos.row=1, layout.pos.col=1)) # 800 x 1200 
print(LFAllAbundance, vp=viewport(layout.pos.row=3, layout.pos.col=1))
print(LFDivAll, vp=viewport(layout.pos.row=2, layout.pos.col=1))
# print(LFDivAll2, vp=viewport(layout.pos.row=3, layout.pos.col=1))

pushViewport(viewport(layout=grid.layout(nrow=3,ncol=1))) # This sets the multiplot.
print(LFSpeciesR, vp=viewport(layout.pos.row=1, layout.pos.col=1)) # 800 x 1200
print(LFAbundanceR, vp=viewport(layout.pos.row=3, layout.pos.col=1))
print(LFDivRecruits, vp=viewport(layout.pos.row=2, layout.pos.col=1))
# print(LFDivRecruits2, vp=viewport(layout.pos.row=3, layout.pos.col=1))

# Dispersal Mode

pushViewport(viewport(layout=grid.layout(nrow=3,ncol=1))) # This sets the multiplot.
print(DPRichAll, vp=viewport(layout.pos.row=1, layout.pos.col=1)) # 800 x 1200
print(DPindiAll, vp=viewport(layout.pos.row=3, layout.pos.col=1))
print(DPDivAll, vp=viewport(layout.pos.row=2, layout.pos.col=1))
# print(DPDivAll2, vp=viewport(layout.pos.row=3, layout.pos.col=1))

pushViewport(viewport(layout=grid.layout(nrow=3,ncol=1))) # This sets the multiplot.
print(DPRichRecruits, vp=viewport(layout.pos.row=1, layout.pos.col=1)) # 800 x 1200
print(DPindiRec, vp=viewport(layout.pos.row=3, layout.pos.col=1))
print(DPDivRecruits, vp=viewport(layout.pos.row=2, layout.pos.col=1))
# print(DPDivRecruits2, vp=viewport(layout.pos.row=3, layout.pos.col=1))

# Supplemental
pushViewport(viewport(layout=grid.layout(nrow=3,ncol=1))) # This sets the multiplot.
print(SupRich, vp=viewport(layout.pos.row=1, layout.pos.col=1)) # 600 x 1200
print(Sup1, vp=viewport(layout.pos.row=3, layout.pos.col=1))
print(SupDivdAll, vp=viewport(layout.pos.row=2, layout.pos.col=1))
#print(SupDivdAll2, vp=viewport(layout.pos.row=2, layout.pos.col=1))

pushViewport(viewport(layout=grid.layout(nrow=3,ncol=1))) # This sets the multiplot.
print(SupRichR, vp=viewport(layout.pos.row=1, layout.pos.col=1)) # 600 x 1200
print(Sup1R, vp=viewport(layout.pos.row=3, layout.pos.col=1))
print(SupDivdR, vp=viewport(layout.pos.row=2, layout.pos.col=1))
#print(SupDivdR2, vp=viewport(layout.pos.row=2, layout.pos.col=1))

################################## BASAL AREA X TREATMENT ################################# 

Basal$Treatment <- relevel(Basal$Treatment, ref = "Wind") 
Basal$Treatment <- relevel(Basal$Treatment, ref = "Animal") 
Basal$Treatment <- relevel(Basal$Treatment, ref = "Forest") 

T3 <- table(Basal$Total.Area/100, Basal$Treatment)
T4 <- aggregate(Basal$Total.Area/100, by=list(Category=Basal$Treatment, Basal$Origin), FUN=sum, na.rm="TRUE")
colnames(T4) <- c("Treatment", "Origin", "Total.Area")

ggplot(data = T4, aes(y = Total.Area, x = Treatment, fill = Origin)) +
  geom_bar(position = "stack", stat="identity", color = "Black") +
  xlab("Habitat") +
  ylab(expression(Total~Basal~Area~(m^2))) + 
  theme(axis.title.x = element_text(size = 20, color = "white")) + 
  theme(axis.text.x = element_text(size = 5.5, color = "black")) + 
  theme(axis.text.y = element_text(size = 7.5, color = "black")) +
  theme(legend.position="none") + theme_bw(base_size = 20)

aggregate(Basal$Total.Area/100, by=list(Category = Basal$Treatment, Category = Basal$Origin), FUN = sum)

# Treatment    Origin Total.Area
# 1     Animal   Planted   2613.056
# 2       Wind   Planted   2509.206
# 3     Forest Recruited   3160.816
# 4     Animal Recruited    619.711
# 5       Wind Recruited    792.475
# 6 Nat. Succ. Recruited   1364.562
# 6   Control Recruited   1364.562
#7 AnimalTotal 2613.056 + 619.711 = 3232.767
#8 WindTotal 2509.206 +  792.475 = 3301.681

2613.056/2564.766 # Treatments are actually fairly equal on area taken up by planted trees.

# Do treatments differ in total basal area per plot?

T10 <- aggregate(Basal$Total.Area/100, by=list(Category=Basal$Treatment, Basal$Plot), FUN=sum, na.rm="TRUE")
colnames(T10) <- c("Treatment", "Plot", "Total.Area")

tgc <- summarySE(T10, measurevar="Total.Area", groupvars=c("Treatment")) 
Fig3C <- ggplot(tgc, aes(x=Treatment, y=Total.Area, fill=Treatment)) + 
  geom_bar(stat="identity", color="black")+
  geom_errorbar(aes(ymin=Total.Area-se, ymax=Total.Area+se), width=.1) +
  xlab("Habitat") + 
  ylab(expression(Plot~Basal~Area~(m^2))) + 
  ylim(0,500) +
  scale_fill_manual(values=c("peachpuff4", "#9999CC","#66CC99", "#CC6666")) +
  geom_text(position=pd, label = c("a", "a", "a", "b"), aes(y = Total.Area+se, x = Treatment),vjust = -0.5, size = 6) +
  theme_bw(base_size = 15) +
  theme(legend.position="none")
Fig3C

T10$Treatment <- relevel(T10$Treatment, ref = "Wind") 
T10$Treatment <- relevel(T10$Treatment, ref = "Animal") 
T10$Treatment <- relevel(T10$Treatment, ref = "Forest") 

lmx <- lm(Total.Area ~ Treatment, data = T10)
summary(lmx)

# Total Basal Area All 
# F-statistic:  7.06 on 3 and 28 DF,  p-value: 0.001115
# Forest x Animal: t = 0.145  0.88612  
# Forest x Wind: t = 0.283  0.77930  
# Forest Control: t = -3.608  0.00119 *
# Animal x Wind: t = 0.138 0.890898
# Animal x Control: t = -3.752 0.000813 *
# Wind x Control: t = -3.891 0.000563 *

# Recruits Only

T11 <- aggregate(BasalW$Total.Area/100, by=list(Category=BasalW$Treatment, BasalW$Plot), FUN=sum, na.rm="TRUE")
colnames(T11) <- c("Treatment", "Plot", "Total.Area")


T11$Treatment <- relevel(T11$Treatment, ref = "Wind") 
T11$Treatment <- relevel(T11$Treatment, ref = "Animal") 
T11$Treatment <- relevel(T11$Treatment, ref = "Forest") 

tgc <- summarySE(T11, measurevar="Total.Area", groupvars=c("Treatment")) 

Fig3CR <- ggplot(tgc, aes(x=Treatment, y=Total.Area, fill=Treatment)) + 
  geom_bar(stat="identity", color="black")+
  geom_errorbar(aes(ymin=Total.Area-se, ymax=Total.Area+se), width=.1) +
  xlab("Habitat") + 
  ylab(expression(Plot~Basal~Area~(m^2))) + 
  ylim(0,500) +
  scale_fill_manual(values=c("peachpuff4", "#9999CC","#66CC99", "#CC6666")) +
  geom_text(position=pd, label = c("a", "", "", "b"), aes(y = Total.Area+se, x = Treatment),vjust = -0.5, size = 6) +
  theme_bw(base_size = 15) +
  theme(legend.position="none")
Fig3CR

T11$Treatment <- relevel(T11$Treatment, ref = "Wind") 
T11$Treatment <- relevel(T11$Treatment, ref = "Animal") 
T11$Treatment <- relevel(T11$Treatment, ref = "Forest") 

lmx <- lm(Total.Area ~ Treatment, data = T11)
summary(lmx)


# Total Basal Area Recruits Only
# F-statistic: 16.09 on 3 and 28 DF,  p-value: 2.841e-06
# Forest x Animal: -6.202 1.06e-06 ***
# Forest x Wind: -5.780 3.31e-06 ***
# Forest x Control:  -4.384 0.000149 ***
# Animal x Wind: 0.422   0.6765  
# Animal x Control: 1.818   0.0798 . 
# Wind x Control:  1.396   0.1736    

# Planted Trees Only 

summary(BasalP$Species)

T1x1 <- aggregate(BasalP$Total.Area/100, by=list(Category=BasalP$Treatment, BasalP$Plot), FUN=sum, na.rm="TRUE")
colnames(T1x1) <- c("Treatment", "Plot", "Total.Area")

tgc <- summarySE(T1x1, measurevar="Total.Area", groupvars=c("Treatment")) 

ggplot(tgc, aes(x=Treatment, y=Total.Area, fill=Treatment)) + 
  geom_bar(stat="identity", color="black")+
  geom_errorbar(aes(ymin=Total.Area-se, ymax=Total.Area+se), width=.1) +
  xlab("Habitat") + 
  ylab(expression(Plot~Basal~Area~(m^2))) + 
  #ylim(0,500) +
  scale_fill_manual(values=c("#9999CC","#66CC99")) +
  geom_text(position=pd, label = c("a", "a"), aes(y = Total.Area+se, x = Treatment),vjust = -0.5, size = 6) +
  theme_bw(base_size = 15) +
  theme(legend.position="none")

lmx <- lm(Total.Area ~ Treatment, data = T1x1)
summary(lmx)

BasalP$AreaHa <- BasalP$Total.AreaMeters /900*10000

T1x2 <- aggregate(BasalP$AreaHa, by=list(Category=BasalP$Treatment, BasalP$Plot), FUN=mean, na.rm="TRUE")
colnames(T1x2) <- c("Treatment", "Plot", "AreaHa")

tgc <- summarySE(T1x2, measurevar="AreaHa", groupvars=c("Treatment")) 

ggplot(tgc, aes(x=Treatment, y=AreaHa, fill=Treatment)) + 
  geom_bar(stat="identity", color="black")+
  geom_errorbar(aes(ymin=AreaHa-se, ymax=AreaHa+se), width=.1) +
  xlab("Habitat") + 
  ylab(expression(Tree~Basal~Area~(m^2/ha))) + 
  #ylim(0,500) +
  scale_fill_manual(values=c("#9999CC","#66CC99")) +
  geom_text(position=pd, label = c("a", "a"), aes(y = AreaHa+se, x = Treatment),vjust = -0.5, size = 6) +
  theme_bw(base_size = 15) +
  theme(legend.position="none")

T1x1$Treatment <- relevel(T1x1$Treatment, ref = "Wind") 
T1x1$Treatment <- relevel(T1x1$Treatment, ref = "Animal") 
T1x1$Treatment <- relevel(T1x1$Treatment, ref = "Forest") 

lmx <- lm(AreaHa ~ Treatment, data = T1x2)
summary(lmx)

table(BasalW$Total.Area, BasalW$Treatment)
T12 <- aggregate(BasalW$Total.Area, by=list(Category=BasalW$Treatment), FUN=max, na.rm="TRUE")

# Distribution of Basal Area (BoxPlots)

ggplot(Basal, aes(x=Plot, y=log(Total.Area/100), fill=Treatment)) + 
  geom_boxplot(show.legend = FALSE) +
  xlab("Plot") +
  ylab(expression(Log~Total~Basal~Area~(m^2))) + 
  theme(axis.title.x = element_text(size = 10, color = "white")) + 
  theme(axis.text.x = element_text(size = 5, color = "black")) + 
  theme(axis.text.y = element_text(size = 7.5, color = "black")) +
  scale_fill_manual(values=c("peachpuff4", "#9999CC","#66CC99", "#CC6666")) +
  facet_grid(~Treatment, scales = "free") +
  theme(legend.position="none") + theme_bw(base_size = 20) 

ggplot(Basal, aes(x=Plot, y=log(Total.Area/100), fill=Treatment)) + 
  geom_boxplot(show.legend = FALSE) +
  xlab("Plot") +
  ylab(expression(Log~Total~Basal~Area~(m^2))) + 
  theme(axis.title.x = element_text(size = 10, color = "white")) + 
  theme(axis.text.x = element_text(size = 5, color = "black")) + 
  theme(axis.text.y = element_text(size = 7.5, color = "black")) +
  scale_fill_manual(values=c("snow4","snow4", "snow4", "snow4")) +
  facet_grid(~Treatment, scales = "free") +
  theme(legend.position="none") + theme_bw(base_size = 20) 

# To calculate totals per category.
TX3 <- TX[-3]
TX4 <- cast(TX3, Category ~ Group.2, value = "Total.Area")
TX4[is.na(TX4)] <- 0

Basal3 <- na.omit(Basal)
#aggregate(Species~Total.Diameter, Basal, FUN=max)
aggregate(Basal$Total.Diameter, by=list(Category=Basal$Species), FUN=max)

BasalTJX <- droplevels(subset(Basal, Basal$Origin == "Planted")) # Exclude Forest. Only Planted.
bintable1 <- as.data.frame(aggregate(BasalTJX$Total.Area, by=list(Category=BasalTJX$Species), FUN=max))
bintable2 <- as.data.frame(aggregate(BasalTJX$Virtual.Diameter, by=list(Category=BasalTJX$Species), FUN=max))
bintable3 <- as.data.frame(aggregate(BasalTJX$Total.Diameter, by=list(Category=BasalTJX$Species), FUN=max))

bintable1$Virtual.Diameter <- bintable2$x
bintable1$Total.Diameter <- bintable3$x

max(bintable1$x) # 7783.1
max(bintable2$x) # 99.5

Basal3 <- Basal[!is.na(Basal$Total.Area),] # Exclude the NA
# Astrocaryum mexicanum  V. Diameter = 10.5 | Area = 50.3
# Cecropia obtusifolia  V. Diameter = 58.6 | Area = 1916.7
# Heliocarpus appendiculatus V. Diameter = 99.5 | Area = 7783.1

max(Basal3$Total.Area)
max(Basal3$Virtual.Diameter)

Basal3$Bins4 <- cut(x = na.omit(Basal3$Total.Area), breaks = c(0,50.3,1916.7,7783.2,28562.2)) # Edited category.
Basal3$Bins5 <- cut(x = na.omit(Basal3$Virtual.Diameter), breaks = c(0,10.5,58.6,99.5,190.7)) # Edited category.

TX5 <- aggregate(Basal3$Total.Area, by=list(Category=Basal3$Bins4, Basal3$Treatment), FUN=sum)
TX5$Total.Area <- as.numeric(TX5$x)
TX6 <- aggregate(TX5$Total.Area, by=list(Category=TX5$Category), FUN=max)

# Diameter = D
TX5D <- aggregate(Basal3$Virtual.Diameter, by=list(Category=Basal3$Bins5, Basal3$Treatment), FUN=sum)
TX5D$Virtual.Diameter <- as.numeric(TX5D$x)
TX6D <- aggregate(TX5D$Virtual.Diameter, by=list(Category=TX5D$Category), FUN=max)

# First find total sum of area per treatment 
TX7 <- aggregate(TX5$Total.Area, by=list(Category=TX5$Group.2), FUN=sum)
TX7D <- aggregate(TX5D$Virtual.Diameter, by=list(Category=TX5D$Group.2), FUN=sum)

# Forest Diameter 

tgc <- summarySE(Basal3, measurevar="Total.Area", groupvars=c("Treatment"))
tgc <- summarySE(Basal3, measurevar="Virtual.Diameter", groupvars=c("Treatment"))

# Then make new total column

TX5$TotalB <- as.numeric(c("316081.6","316081.6","316081.6","316081.6",
                           "323276.7","323276.7","323276.7",
                           "330168.1","330168.1","330168.1",
                           "136456.2","136456.2","136456.2","136456.2"))

TX5D$TotalB <- as.numeric(c("12604.6","12604.6","12604.6","12604.6",
                            "16325.8","16325.8","16325.8",
                            "20742.7","20742.7","20742.7",
                            "10462.8","10462.8","10462.8","10462.8"))

# Then add percentage column
TX5$Percent <- (TX5$Total.Area/TX5$TotalB)*100
TX5$Percent <- sprintf(TX5$Percent, fmt='%#.2f')

TX5 <- ddply(TX5, .(Group.2), # This code will add the position of the text.
             transform, pos = TotalB - cumsum(Total.Area) + (0.5 * Total.Area))

ggplot(TX5, aes(fill=Category, y=Total.Area, x=Group.2)) + 
  geom_bar(position="stack", stat="identity", color = "Black") +
  xlab("Habitat") +
  scale_fill_brewer(palette="Set3", labels = c("Understory", "Midstory", "Planting Canopy", "Forest Canopy")) +
  labs(fill = "Occupied Strata", y = 'Total Basal Area' ~(cm^2)) + 
  geom_text(data=TX5, aes(x = Group.2, y = pos, label = paste0(Percent,"%")), size=4) +
  scale_y_continuous(label=comma) +
  theme_bw(base_size = 20)

# By diameter

TX5D$Percent <- (TX5D$Virtual.Diameter/TX5$TotalB)*100
TX5D$Percent <- sprintf(TX5D$Percent, fmt='%#.2f')

TX5D <- ddply(TX5D, .(Group.2), # This code will add the position of the text.
              transform, pos = TotalB - cumsum(Virtual.Diameter) + (0.5 * Virtual.Diameter))

ggplot(TX5D, aes(fill=Category, y=Virtual.Diameter, x=Group.2)) + 
  geom_bar(position="stack", stat="identity", color = "Black") +
  xlab("Habitat") +
  scale_fill_brewer(palette="Set3", labels = c("Understory", "Midstory", "Planting Canopy", "Forest Canopy")) +
  labs(fill = "Occupied Strata", y = 'Total Diameter' ~(cm^2)) + 
  geom_text(data=TX5D, aes(x = Group.2, y = pos, label = paste0(Percent,"%")), size=4) +
  scale_y_continuous(label=comma) +
  theme_bw(base_size = 20)


TX5.5<-TX5
TX5.5$Total.Aream <- TX5.5$Total.Area/100

TX5.5$TotalBm <- TX5.5$TotalB/100
TX5.5$Percent <- (TX5.5$Total.Area/TX5.5$TotalBm)
TX5.5$Percent <- sprintf(TX5.5$Percent, fmt='%#.2f')

TX5.5 <- ddply(TX5.5, .(Group.2), # This code will add the position of the text.
               transform, pos = (TotalBm) - cumsum(Total.Aream) + (0.5 * Total.Aream))

# In meters scale
Fig4B <- ggplot(TX5.5, aes(fill=Category, y=Total.Aream, x=Group.2)) + 
  geom_bar(position="stack", stat="identity", color = "Black") +
  xlab("Habitat") +
  scale_fill_brewer(palette="Set3", labels = c("Understory", "Midstory", "Planting Canopy", "Forest Canopy")) +
  labs(fill = "Occupied Strata", y = 'Total Habitat Basal Area' ~(m^2)) + 
  geom_text(data=TX5.5, aes(x = Group.2, y = pos, label = paste0(Percent,"%")), size=4) +
  scale_y_continuous(label=comma) +
  theme_bw(base_size = 15)
Fig4B

table(Basal$Origin, Basal$Treatment)

1242 + 2114    +   1377

ggplot(TX5, aes(fill=Category, y=Total.Area, x=Group.2)) + 
  geom_bar(position="stack", stat="identity", color = "Black") +
  xlab("Habitat") +
  scale_fill_brewer(palette="Set3", labels = c("0 - 73.9","73.9 - 1,916.7", "1,916.7 - 7,783.2","7,783.2 - 28,562.2")) +
  labs(fill = 'Basal Area' ~(cm^2), y = 'Total Basal Area' ~(cm^2)) + 
  # geom_text(data=T3, aes(x = Treatment, y = pos, label = paste0(Percent,"%")), size=4) +
  scale_y_continuous(label=comma) +
  theme_bw(base_size = 15)


TX5

M1 <- as.table(rbind(c(21597.2, 12688.2, 25692.0,14219.4 ), 
                     c(93813.6, 193526.7, 231540.7, 92510.0),
                     c(59713.3, 117061.8, 72935.4, 18921.1),
                     c(140957.5, 0, 0, 10805.7))) # Updated.

dimnames(M1) <- list(Category = c("Understory","Midstory","Planting Canopy", "Forest Canopy"),
                     Treatment = c("Forest","Animal","Wind","Control"))                     

chisq.test(M1)
# X-squared = 404885, df = 9, p-value < 2.2e-16

chisq.posthoc.test(M1, method = "bonferroni")

# Dimension     Value      Forest     Animal       Wind   Control
# 1      Understory Residuals    3.300461  -75.20758   29.42036  58.53796
# 2      Understory  p values    0.015444    0.00000    0.00000   0.00000
# 3        Midstory Residuals -342.527569   62.30877  204.87534  99.30315
# 4        Midstory  p values    0.000000    0.00000    0.00000   0.00000
# 5 Planting Canopy Residuals  -83.728509  187.90125  -35.17378 -95.89483
# 6 Planting Canopy  p values    0.000000    0.00000    0.00000   0.00000
# 7   Forest Canopy Residuals  596.912816 -269.53872 -273.60364 -66.54277
# 8   Forest Canopy  p values    0.000000    0.00000    0.00000   0.00000

# Post-Hoc with Bonferoni Correction (Joined Canopy Categories)
# 59713.3 + 140957.5 = 200670.8 # Forest
# 18921.1 + 10805.7 = 29726.8 # Control 

M2 <- as.table(rbind(c(21597.2, 12688.2, 25692.0,14219.4), 
                     c(93813.6, 193526.7, 231540.7, 92510.0),
                     c(200670.8,117061.8,72935.4, 29726.8))) # Updated

dimnames(M2) <- list(Category = c("Understory","Midstory","Canopy"),
                     Treatment = c("Forest","Animal","Wind","Control"))

chisq.test(M2)
# X-squared = 148558, df = 6, p-value < 2.2e-16

chisq.posthoc.test(M2, method = "bonferroni")

# Dimension     Value      Forest    Animal       Wind    Control
# 1 Understory Residuals    3.300461 -75.20758   29.42036   58.53796
# 2 Understory  p values    0.011583   0.00000    0.00000    0.00000
# 3   Midstory Residuals -342.527569  62.30877  204.87534   99.30315
# 4   Midstory  p values    0.000000   0.00000    0.00000    0.00000
# 5     Canopy Residuals  349.146370 -25.06177 -225.01430 -131.88443
# 6     Canopy  p values    0.000000   0.00000    0.00000    0.00000


# To calculate percentages of the parts: 
TX4.5$Category2 <- c("Understory", "Midstory", "Planting Canopy", "Forest Canopy")
sum(TX4.5$Forest) # 316081.6
sum(TX4.5$Animal) # 323276.7
sum(TX4.5$Wind) # 330168.1
sum(TX4.5$`Nat. Succ.`) # 136456.2

# How much of the total did canopy trees make up?
# Forest
(59713.3+140957.5)/316081.6*100 # 63.487%
# Animal
(117061.8)/323276.7*100 # 36.211%
# Wind
(72935.4)/330168.1*100 # 22.090
# Natural Succession
(18921.1+10805.7)/136456.2*100 # 21.785

# How much of the total did midstory trees make up?
# Forest
(91347.1)/316081.6*100 # 28.89985
# Animal
(189978.8 )/323276.7*100 # 58.76662
# Wind
(223745.0)/330168.1*100 # 67.767
# Natural Succession
(88362.3)/136456.2*100 # 64.75506

# How much of the total did understory trees make up?
# Forest
(24063.7)/316081.6*100 # 7.613129
# Animal
(16236.1 )/323276.7*100 # 5.022354
# Wind
(33487.7)/330168.1*100 # 10.14262
# Natural Succession
(18367.1)/136456.2*100 # 13.46007

# Lets look at it a different way:
chisq.test(M2, correct = FALSE)$stdres # These are the residuals. 
sig.lvl <- 0.05

# Let's make the adjustment manually:
adj.sig <- sig.lvl/(nrow(M2)*ncol(M2))
adj.sig # 0.0042

# The critical z-value:
qnorm(adj.sig/2) # two-tailed so divided by two.
# -2.86526

# None of our residuals are within those values, so therefore significant differences all around.

# To determine how much area is taken up by the ten largest trees vs. the rest. 
Basal$Treatment <- relevel(Basal$Treatment, ref = "Wind") 
Basal$Treatment <- relevel(Basal$Treatment, ref = "Animal") 
Basal$Treatment <- relevel(Basal$Treatment, ref = "Forest") 

T1 <- aggregate(Basal$Total.Area, list(Basal$Treatment), FUN=sum, na.rm = TRUE)
T1

DB1 <- Basal[with(Basal,order(-Total.Area)),]
DB1 <- droplevels(subset(DB1, DB1$Treatment == "Forest"))
DB1 <- DB1[1:10,-9:-59]
F10 <- sum(DB1$Total.Area)

DB2 <- Basal[with(Basal,order(-Total.Area)),]
DB2 <- droplevels(subset(DB2, DB2$Treatment == "Animal"))
DB2 <- DB2[1:10,-9:-59]
A10 <- sum(DB2$Total.Area)

DB3 <- Basal[with(Basal,order(-Total.Area)),]
DB3 <- droplevels(subset(DB3, DB3$Treatment == "Wind"))
DB3 <- DB3[1:10,-9:-59]
W10 <- sum(DB3$Total.Area)

DB4 <- Basal[with(Basal,order(-Total.Area)),]
DB4 <- droplevels(subset(DB4, DB4$Treatment == "Nat. Succ."))
DB4 <- DB4[1:10,-9:-59]
C10 <- sum(DB4$Total.Area)

T1$Top10 <- c(F10,A10,W10,C10)
T1$Percent <- (T1$Top10/T1$x)*100
T1$Below10 <- T1$x - T1$Top10

T2<-T1
T2$Top10 <- T2$Below10
T3 <- rbind(T1,T2)
T3$TB <- c("Top", "Top", "Top","Top","Below", "Below", "Below", "Below")
names(T3)[1]<-"Treatment"
names(T3)[2]<-"TotalBasalArea"
names(T3)[3]<-"BasalArea"
names(T3)[6]<-"Top10orNot"
T3$Percent <- (T3$BasalArea/T3$TotalBasalArea)*100
T3$Percent <- sprintf(T3$Percent, fmt='%#.2f')

T3 <- ddply(T3, .(Treatment), # This code will add the position of the text.
            transform, pos = cumsum(BasalArea/100) - (0.5 * BasalArea/100))

T3$Treatment <- relevel(T3$Treatment, ref = "Wind") 
T3$Treatment <- relevel(T3$Treatment, ref = "Animal") 
T3$Treatment <- relevel(T3$Treatment, ref = "Forest") 

Fig4A <- ggplot(data = T3, aes(y = BasalArea/100 , x = Treatment, fill = Top10orNot)) + 
  geom_bar(colour="black",stat="identity") +
  geom_text(data=T3, aes(x = Treatment, y = pos, label = paste0(Percent,"%")), size=4) +
  xlab("Habitat") +
  ylab(expression(Total~Plot~Basal~Area~(m^2))) + 
  theme(axis.title.x = element_text(size = 20, color = "white")) + 
  theme(axis.text.x = element_text(size = 7.5, color = "black")) + 
  theme(axis.text.y = element_text(size = 7.5, color = "black")) +
  scale_fill_brewer(name = "Trees", labels = c("Remaining", "10 Largest"), palette = "Greens") +
  theme(legend.position="none") + theme_bw(base_size = 20)

ggplot(data = T3, aes(y = BasalArea/100 , x = Treatment, fill = Top10orNot)) + 
  geom_bar(colour="black",stat="identity") +
  geom_text(data=T3, aes(x = Treatment, y = pos, label = paste0(Percent,"%")), size=4) +
  xlab("Habitat") +
  ylab(expression(Total~Basal~Area~(m^2))) + 
  theme(axis.title.x = element_text(size = 20, color = "white")) + 
  theme(axis.text.x = element_text(size = 7.5, color = "black")) + 
  theme(axis.text.y = element_text(size = 7.5, color = "black")) +
  scale_fill_brewer(name = "Trees", labels = c("Remaining", "10 Largest"), palette = "Greys") +
  theme(legend.position="none") + theme_bw(base_size = 20)

# Basal Area per tree x Treatment
Basal$Treatment <- relevel(Basal$Treatment, ref = "Forest") 
Basal$Total.AreaMeters <- Basal$Total.Area/100
Basal$AreaHa <- Basal$Total.AreaMeters /900*10000
BasalW$Treatment <- as.factor(BasalW$Treatment)
BasalW$Treatment <- relevel(BasalW$Treatment, ref = "Wind") 
BasalW$Treatment <- relevel(BasalW$Treatment, ref = "Animal") 
BasalW$Treatment <- relevel(BasalW$Treatment, ref = "Forest") 
BasalW$Total.AreaMeters <- BasalW$Total.Area/100
BasalW$AreaHa <- BasalW$Total.AreaMeters /900*10000

tgc2 <- summarySE(Basal, measurevar="Total.AreaMeters", groupvars=c("Treatment"), na.rm = TRUE)
ggplot(tgc2, aes(x=Treatment, y=Total.AreaMeters, fill=Treatment)) + 
  geom_bar(stat="identity", color="black")+
  geom_errorbar(aes(ymin=Total.AreaMeters-se, ymax=Total.AreaMeters+se), width=.1) +
  xlab("Habitat") + 
  ylab(expression(Average~Basal~Area~per~Tree~(m^2))) + 
  ylim(0,2.7) +
  scale_fill_manual(values=c("peachpuff4", "#9999CC","#66CC99", "#CC6666")) +
  theme_bw(base_size = 15) +
  #geom_text(label = c("a", "a", "b", "b"), aes(y = Total.AreaMeters+se, x = Treatment),vjust = -0.5, size = 6) +
  theme(legend.position="none")

ggplot(tgc2, aes(x=Treatment, y=Total.AreaMeters, fill=Treatment)) + 
  geom_bar(stat="identity", color="black")+
  geom_errorbar(aes(ymin=Total.AreaMeters-se, ymax=Total.AreaMeters+se), width=.1) +
  xlab("Habitat") + 
  ylab(expression(Average~Basal~Area~per~Tree~(m^2))) + 
  ylim(0,2.7) +
  scale_fill_manual(values=c("snow4","snow4", "snow4", "snow4")) +
  theme_bw(base_size = 15) +
  geom_text(label = c("a", "a", "b", "b"), aes(y = Total.AreaMeters+se, x = Treatment),vjust = -0.5, size = 6) +
  theme(legend.position="none")

Basal$Treatment <- as.factor(Basal$Treatment)
Basal$Treatment <- relevel(Basal$Treatment, ref = "Wind") 
Basal$Treatment <- relevel(Basal$Treatment, ref = "Animal") 
Basal$Treatment <- relevel(Basal$Treatment, ref = "Forest") 

pushViewport(viewport(layout=grid.layout(nrow=1,ncol=2))) # This sets the multiplot.
print(Fig4A, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(Fig4B, vp=viewport(layout.pos.row=1, layout.pos.col=2))

################################## Basal Area (m2/ha) X TREATMENT ################################# 
Basal$Treatment <- relevel(Basal$Treatment, ref = "Wind") 
Basal$Treatment <- relevel(Basal$Treatment, ref = "Animal") 
Basal$Treatment <- relevel(Basal$Treatment, ref = "Forest") 

tgc3 <- summarySE(Basal, measurevar="AreaHa", groupvars=c("Treatment"), na.rm = TRUE)

Fig3A <- ggplot(tgc3, aes(x=Treatment, y=AreaHa, fill=Treatment)) + 
  geom_bar(stat="identity", color="black")+
  geom_errorbar(aes(ymin=AreaHa-se, ymax=AreaHa+se), width=.1) +
  xlab("Habitat") + 
  ylim(0,30) +
  ylab(expression(Individual~Tree~Basal~Area~(m^2/ha))) + 
  geom_text(label = c("a", "a", "b", "b"), aes(y = AreaHa+se, x = Treatment),vjust = -0.5, size = 6) +
  scale_fill_manual(values=c("peachpuff4", "#9999CC","#66CC99", "#CC6666")) +
  theme_bw(base_size = 15) +
  theme(legend.position="none")
Fig3A

Basal$Treatment <- relevel(Basal$Treatment, ref = "Wind") 

lm2 <- lm(AreaHa ~ Treatment, data = Basal) 
summary(lm2)

# Tree Basal Area All Plants
# F-statistic: 7.107 on 3 and 7205 DF,  p-value: 9.161e-05
# Forest x Animal: t = -0.749 0.453735  
# Forest x Wind: t = -3.300 0.000970 ***
# Forest x Control: t = -3.818 0.000136 ***
# Animal x Wind: t = -2.581  0.00988 ** 
# Animal x Control: t = -3.199  0.00138 **
# Wind x Control: t = -1.067  0.28618 

# Only Recruits
BasalW$Treatment <- relevel(BasalW$Treatment, ref = "Wind") 
BasalW$Treatment <- relevel(BasalW$Treatment, ref = "Animal") 
BasalW$Treatment <- relevel(BasalW$Treatment, ref = "Forest") 

tgc4 <- summarySE(BasalW, measurevar="AreaHa", groupvars=c("Treatment"), na.rm = TRUE)

Fig3AR <- ggplot(tgc4, aes(x=Treatment, y=AreaHa, fill=Treatment)) + 
  geom_bar(stat="identity", color="black")+
  geom_errorbar(aes(ymin=AreaHa-se, ymax=AreaHa+se), width=.1) +
  xlab("Habitat") + 
  ylim(0,30) +
  ylab(expression(Individual~Tree~Basal~Area~(m^2/ha))) + 
  geom_text(label = c("a", "", "", "b"), aes(y = AreaHa+se, x = Treatment),vjust = -0.5, size = 6) +
  scale_fill_manual(values=c("peachpuff4", "#9999CC","#66CC99", "#CC6666")) +
  theme_bw(base_size = 15) +
  theme(legend.position="none")

Fig3AR

BasalW$Treatment <- relevel(BasalW$Treatment, ref = "Wind") 

lm2 <- lm(AreaHa ~ Treatment, data = BasalW) 
summary(lm2)

# Tree Basal Area Recruits only
# F-statistic: 17.59 on 3 and 6240 DF,  p-value: 2.282e-11
# Forest x Animal: -5.587 2.41e-08 ***
# Forest x Wind: -6.847 8.28e-12 ***
# Forest X Control: -3.969 7.28e-05 ***
# Animal x Wind: -0.466   0.6410 
# Animal x Control: 1.689   0.0913 .
# Wind x Control: 2.390   0.0169 *  

pushViewport(viewport(layout=grid.layout(nrow=1,ncol=2))) # This sets the multiplot. 1000 x 500
print(Fig3A, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(Fig3C, vp=viewport(layout.pos.row=1, layout.pos.col=2))

pushViewport(viewport(layout=grid.layout(nrow=1,ncol=2))) # This sets the multiplot.
print(Fig3AR, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(Fig3CR, vp=viewport(layout.pos.row=1, layout.pos.col=2))

# Stacked BarPlot for Area Covered by Planted Trees Only
BasalTJ <- droplevels(subset(Basal, Basal$Treatment != "Forest")) # Exclude Forest
BasalTJ <- droplevels(subset(BasalTJ, BasalTJ$Treatment != "Nat. Succ.")) # Exclude Forest
BasalTJ 
BasalTJ$Origin <- as.factor(BasalTJ$Origin)

# Stacked
TJ <- aggregate(BasalTJ$Total.AreaMeters, by = list(Category=BasalTJ$Origin,BasalTJ$Treatment), FUN=sum, na.rm=TRUE)
colnames(TJ) <- c("Origin", "Treatment", "Total.AreaMeters") # Rename columns
TJ <- droplevels(subset(TJ, TJ$Treatment != "Nat. Succ.")) 
TJ <- droplevels(subset(TJ, TJ$Treatment != "Forest")) 
TJ$Origin <- as.factor(TJ$Origin)

TJ$Origin <- relevel(TJ$Origin, ref = "Recruited") 
TJ

Tax <- as.data.frame(table(BasalTJ$Origin, BasalTJ$Treatment))
colnames(Tax) <- c("Origin", "Treatment", "Trees") # Rename columns
# Animal Wind
# Planted      471  494
# Recruited   1242 2114
Tax$Total.AreaMeters <- TJ$Total.AreaMeters

Tax <- ddply(Tax, .(Treatment), # This code will add the position of the text.
             transform, pos = TotalB - Total.AreaMeters + (0.5 * Total.AreaMeters))

Tax$pos2 <- as.numeric(c("1500", "2900", "1700", "2950"))

Tax$TotalB <- as.numeric(c("3232.767", "3232.767", "3301.681", "3301.681"))

ggplot(TJ, aes(fill=Origin, y=Total.AreaMeters, x=Treatment)) + 
  geom_bar(position="stack", stat="identity", color = "black") + 
  xlab("Planting Treatment") +
  labs(y = 'Total Basal Area' ~(m^2)) + 
  geom_text(data=Tax, aes(x = Treatment, y = pos2, label = Trees), size=5) +
  scale_y_continuous(label=comma) +
  scale_fill_brewer(palette="Set3") +
  theme_bw(base_size = 15) 

# How many trees?

ggplot(TX5, aes(fill=Category, y=Total.Area, x=Group.2)) + 
  geom_bar(position="stack", stat="identity", color = "Black") +
  xlab("Habitat") +
  scale_fill_brewer(palette="Set3", labels = c("Understory", "Midstory", "Planting Canopy", "Forest Canopy")) +
  labs(fill = "Occupied Strata", y = 'Total Basal Area' ~(cm^2)) + 
  geom_text(data=TX5, aes(x = Group.2, y = pos2, label = paste0(Percent,"%")), size=4) +
  scale_y_continuous(label=comma) +
  theme_bw(base_size = 20)


BasalPalmLess2 <- droplevels(subset(Basal, Basal$Species != "Astrocaryum mexicanum")) # Exclude Palms
tgc5 <- summarySE(BasalPalmLess2, measurevar="AreaHa", groupvars=c("Treatment"), na.rm = TRUE)


BasalPalmLess <- droplevels(subset(Basal, Basal$Family != "Arecaceae")) # Exclude Palms

tgc5 <- summarySE(BasalPalmLess, measurevar="AreaHa", groupvars=c("Treatment"), na.rm = TRUE)

ggplot(tgc5, aes(x=Treatment, y=AreaHa, fill=Treatment)) + 
  geom_bar(stat="identity", color="black")+
  geom_errorbar(aes(ymin=AreaHa-se, ymax=AreaHa+se), width=.1) +
  xlab("Habitat") + 
  ylab(expression(Basal~Area~(m^2/ha))) + 
  #ylim(0,2.7) +
  scale_fill_manual(values=c("peachpuff4", "#9999CC","#66CC99", "#CC6666")) +
  theme_bw(base_size = 20) +
  theme(legend.position="none")

BasalPalmLess$Treatment <- relevel(BasalPalmLess$Treatment, ref = "Forest") 

lm2 <- lm(AreaHa ~ Treatment, data = BasalPalmLess) 
summary(lm2)

# F-statistic: 23.77 on 3 and 6473 DF,  p-value: 2.687e-15

# Forest x Animal: -5.467 4.75e-08 *
# Forest x Wind:-7.658 2.17e-14 *
# Forest X Control:-7.740 1.15e-14 *
# Animal x Wind:-2.440  0.01472 *
# Animal x Control:-3.047  0.00232 *
# Wind x Control: -1.030   0.3029

# Without Palms and Only Recruits

BasalWPalmLess <- droplevels(subset(BasalW, BasalW$Family != "Arecaceae")) # Exclude Palms

tgc6 <- summarySE(BasalWPalmLess, measurevar="AreaHa", groupvars=c("Treatment"), na.rm = TRUE)

ggplot(tgc6, aes(x=Treatment, y=AreaHa, fill=Treatment)) + 
  geom_bar(stat="identity", color="black")+
  geom_errorbar(aes(ymin=AreaHa-se, ymax=AreaHa+se), width=.1) +
  xlab("Habitat") + 
  ylab(expression(Basal~Area~(m^2/ha))) + 
  #ylim(0,2.7) +
  scale_fill_manual(values=c("peachpuff4", "#9999CC","#66CC99", "#CC6666")) +
  theme_bw(base_size = 20) +
  theme(legend.position="none")

BasalWPalmLess$Treatment <- relevel(BasalWPalmLess$Treatment, ref = "Wind") 

lm2 <- lm(AreaHa ~ Treatment, data = BasalWPalmLess) 
summary(lm2)

# F-statistic: 37.48 on 3 and 5483 DF,  p-value: < 2.2e-16

# Forest x Animal: -8.537  < 2e-16 
# Forest x Wind: -10.360  < 2e-16
# Forest X Control: -7.946 2.32e-15
# Animal x Wind: -1.368 0.171443 
# Animal x Control: 0.799 0.424508
# Wind x Control:  2.281   0.0226

# Basal Area per plot x Treatment
T6 <- aggregate(Basal$Total.AreaMeters, by = list(Category=Basal$Plot,Basal$Treatment), FUN=sum, na.rm=TRUE)
colnames(T6) <- c("Plot", "Treatment", "Total.AreaMeters")

tgc2 <- summarySE(T6, measurevar="Total.AreaMeters", groupvars=c("Treatment"), na.rm = TRUE)

ggplot(tgc2, aes(x=Treatment, y=Total.AreaMeters, fill=Treatment)) + 
  geom_bar(stat="identity", color="black")+
  geom_errorbar(aes(ymin=Total.AreaMeters-se, ymax=Total.AreaMeters+se), width=.1) +
  xlab("Habitat") + 
  ylab(expression(Basal~Area~(m^2))) + 
  scale_fill_manual(values=c("peachpuff4", "#9999CC","#66CC99", "#CC6666")) +
  theme_bw(base_size = 20) +
  # geom_text(label = c("a", "a", "b", "b"), aes(y = Total.AreaMeters+se, x = Treatment),vjust = -0.5, size = 6) +
  theme(legend.position="none")

BasalX$Total.AreaMeters <- BasalX$Total.Area/100
T6.2 <- aggregate(BasalX$Total.AreaMeters, by = list(Category=BasalX$Plot,BasalX$Treatment), FUN=sum, na.rm=TRUE)
colnames(T6.2) <- c("Plot", "Treatment", "Total.AreaMeters")

tgc3 <- summarySE(T6.2, measurevar="Total.AreaMeters", groupvars=c("Treatment"), na.rm = TRUE)

ggplot(tgc3, aes(x=Treatment, y=Total.AreaMeters, fill=Treatment)) + 
  geom_bar(stat="identity", color="black")+
  geom_errorbar(aes(ymin=Total.AreaMeters-se, ymax=Total.AreaMeters+se), width=.1) +
  xlab("Habitat") + 
  ylab(expression(Basal~Area~(m^2))) + 
  scale_fill_manual(values=c("peachpuff4", "#9999CC","#66CC99", "#CC6666")) +
  theme_bw(base_size = 20) +
  # geom_text(label = c("a", "a", "b", "b"), aes(y = Total.AreaMeters+se, x = Treatment),vjust = -0.5, size = 6) +
  theme(legend.position="none")

ggplot(tgc2, aes(x=Treatment, y=Total.AreaMeters, fill=Treatment)) + 
  geom_bar(stat="identity", color="black")+
  geom_errorbar(aes(ymin=Total.AreaMeters-se, ymax=Total.AreaMeters+se), width=.1) +
  xlab("Habitat") + 
  ylab(expression(Average~Basal~Area~per~Plot~(m^2))) + 
  ylim(0,2.7) +
  scale_fill_manual(values=c("snow4","snow4", "snow4", "snow4")) +
  theme_bw(base_size = 20) +
  geom_text(label = c("a", "a", "b", "b"), aes(y = Total.AreaMeters+se, x = Treatment),vjust = -0.5, size = 6) +
  theme(legend.position="none")

T6$Treatment <- relevel(T6$Treatment, ref = "Wind") 
lmx <- lm(data = T6, Total.AreaMeters ~ Treatment)
summary(lmx)

# F-statistic:  7.06 on 3 and 28 DF,  p-value: 0.001115

# Forest x Animal     8.994     62.232   0.145  0.88612    
# Forest x Wind      17.608     62.232   0.283  0.77930    
# Forest x Control -224.532     62.232  -3.608  0.00119 ** 
# Animal x Wind       8.614     62.232   0.138 0.890898    
# Animal x Control -233.526     62.232  -3.752 0.000813 ***
# Wind x Control -242.140     62.232  -3.891 0.000563 ***

T6.2$Treatment <- relevel(T6.2$Treatment, ref = "Animal") 
lmx <- lm(data = T6.2, Total.AreaMeters ~ Treatment)
summary(lmx) # Recruits Only

# F-statistic: 14.05 on 3 and 28 DF,  p-value: 8.919e-06

# Forest x Animal    -274.7       51.9  -5.293 1.25e-05 ***
# Forest x Wind      -303.0       51.9  -5.838 2.83e-06 ***
# Forest x Control   -224.5       51.9  -4.326 0.000174 ***
# Animal x Wind     28.30      51.90   0.545   0.5899    
# Animal x Control    50.16      51.90   0.966  0.34209    
# Wind x Control    78.46      51.90   1.512   0.1418  

######################################### NMDS ######################################### 
######################################### NMDS With All by Species ######################################### 
Basal$Treatment <- relevel(Basal$Treatment, ref = "Wind") 
Basal$Treatment <- relevel(Basal$Treatment, ref = "Animal") 
Basal$Treatment <- relevel(Basal$Treatment, ref = "Forest") 

T5 <- as.data.frame.matrix(table(Basal$Species,Basal$Plot))
T5 <- t(as.data.frame.matrix(table(Basal$Species,Basal$Plot)))
# T5.2 <- t(T5[c(25,26,27,28,29,30,31,32,2,4,7,12,14,16,20,22,3,5,9,11,13,18,21,23,1,6,8,10,15,17,19,24)]) # Keeps Species ID & reorganizes for later separation of treatment with NMDS.
# T5.2 allows us to order the plots by treatment, which is useful for illustrating the NDMS. 
# We want to use T5 for the PERMANOVA though. 

LSNMDS<-metaMDS(T5, distance="bray", k=2, trymax=35, autotransform=TRUE) # k is the number of dimensions
LSNMDS # Stress:  0.06380049 
# Stress: similarity of observed distance to ordination distance. < 0.15 to indidates acceptable fit.

# stressplot(LSNMDS)

tempo <- as.data.frame(as.factor(1:32))
colnames(tempo) <- c("Plot")
tempo$Treatment <- with(Basal, Treatment[match(tempo$Plot, Plot)]) # Match Treatment by Plot

# Pull points from NMDS
NMDS1 <- LSNMDS$points[,1] ##also found using scores(birdMDS)
NMDS2 <- LSNMDS$points[,2]
NMDS<-cbind(tempo, NMDS1, NMDS2)
NMDS$Habitat <- NMDS$Treatment

NMDS$Treatment <- relevel(NMDS$Treatment, ref = "Wind") 
NMDS$Treatment <- relevel(NMDS$Treatment, ref = "Animal") 
NMDS$Treatment <- relevel(NMDS$Treatment, ref = "Forest") 

# Plot NMDS (Treatment)
ggplot(NMDS, aes(NMDS1, NMDS2, color=Habitat))+
  geom_point(position=position_jitter(.1), shape=4, size = 2)+##separates overlapping points
  stat_ellipse(type='t',size =.5)+ ##draws 95% confidence interval ellipses
  #stat_ellipse(aes(group = Treatment), color = "black", type='t',size =1) +
  #geom_text(data=NMDS,aes(x=NMDS1,y=NMDS2,label=Plot),size=3,vjust=0) +  # add the site labels
  scale_color_manual(values=c("peachpuff4","#9999CC", "#66CC99", "#CC6666")) +
  theme_bw(base_size = 15) 

vare.pca <- prcomp(T5)
ewf <- as.data.frame(scores(vare.pca, choices=c(1,2), display = "species"))
ewf

NMDS$points

spe.bray <- vegdist(T5) 

PERMTREAT <- adonis2(spe.bray~Treatment, data=tempo, permutations = 999, method="bray", strata="Treatment")
PERMTREAT # F = 17.352 p = 0.001 

adonis(spe.bray~Treatment, data=tempo)

pairwise.adonis2(spe.bray~Treatment, data=tempo)
# Animal vs Wind: F = 8.316 0.37265  0.001 ***
# Animal vs Control: F = 5.5008 0.28208  0.002 **
# Control vs Wind: F =  5.42 0.27909  0.001 ***
# Forest vs Animal: F = 27.432 0.6621  0.001 ***
# Forest vs Wind: F = 8.316 0.37265  0.001 ***
# Forest vs Control: F = 25.868 0.64884  0.001 ***

# Adjust p-values
k <- c(0.001,
       0.001,
       0.001,
       0.001,
       0.001,
       0.001)

p.adjust(k, method = "fdr", n = length(k)) # Ajuste a los p-values igual.

# Multivariate dispersion
dispersion<-betadisper(spe.bray, group=tempo$Treatment)
permutest(dispersion) # F = 3.1899, Df = 3,28, p = 0.036
anova(dispersion)
plot(dispersion, hull=FALSE, ellipse=TRUE) # Standard Deviation Ellipse?

# Pairwise Dispersion
(mod.HSD <- TukeyHSD(dispersion))
plot(mod.HSD)

######################################### NMDS Without Planted Trees by Species ######################################### 
BasalW$Species <- as.factor(BasalW$Species)

# Without Planted Trees (Species Scale)
T7 <- t(as.data.frame.matrix(table(BasalW$Species,BasalW$Plot)))
T7.2 <- t(T7[c(25,26,27,28,29,30,31,32,2,4,7,12,14,16,20,22,3,5,9,11,13,18,21,23,1,6,8,10,15,17,19,24)]) # Keeps Species ID & reorganizes for later separation of treatment with NMDS.

LSNMDS<-metaMDS(T7, distance="bray", k=2, trymax=35, autotransform=TRUE,weakties = TRUE) # k is the number of dimensions
LSNMDS # stress is (nearly) zero: you may have insufficient data

#LSNMDS<-metaMDS(T7, distance="euc", k=2, trymax=35, autotransform=TRUE) # k is the number of dimensions
# LSNMDS # 0.1276526. Doesn't work with bray-curtis.

# Pull points from NMDS
NMDS1 <- LSNMDS$points[,1] ##also found using scores(birdMDS)
NMDS2 <- LSNMDS$points[,2]
NMDS<-cbind(tempo, NMDS1, NMDS2)

NMDS$Treatment <- relevel(NMDS$Treatment, ref = "Wind") 
NMDS$Treatment <- relevel(NMDS$Treatment, ref = "Animal") 
NMDS$Treatment <- relevel(NMDS$Treatment, ref = "Forest") 

# Plot NMDS (Treatment)
ggplot(NMDS, aes(NMDS1/100, NMDS2, color=Treatment))+
  geom_point(position=position_jitter(.1), shape=4)+##separates overlapping points
  stat_ellipse(type='t',size =.5)+ ##draws 95% confidence interval ellipses
  #geom_text(data=NMDS,aes(x=NMDS1,y=NMDS2,label=Plot),size=3,vjust=0) +  # add the site labels
  scale_color_manual(values=c("peachpuff4","#9999CC", "#66CC99", "#CC6666")) +
  theme_bw(base_size = 20) 

scores(NMDS, choices=c(1,2))
vare.pca <- prcomp(T7)
scores(vare.pca, choices=c(1,2), display = "species")


spe.bray2 <- vegdist(T7) 
PERMTREAT2 <- adonis2(spe.bray2~Treatment, data=tempo, permutations = 999, method="bray")
# F = 13.701, p = 0.001 ***

adonis(spe.bray2~Treatment, data=tempo)

pairwise.adonis2(spe.bray2~Treatment, data=tempo)
# Animal vs Wind: F = 2.4035 0.14652   0.01 **
# Animal vs Control: F = 2.2467 0.13829  0.011 *
# Control vs Wind: F =  2.5092 0.15199  0.013 *
# Forest vs Animal: F = 23.343 0.6251  0.001 ***
# Forest vs Wind: F = 24.265 0.63413  0.001 ***
# Forest vs Control: F = 25.868 0.64884  0.001 ***

# Adjust p-values
k <- c(0.01,
       0.011,
       0.013,
       0.001,
       0.001,
       0.001)

p.adjust(k, method = "fdr", n = length(k)) # Ajuste a los p-values igual.
# [1] 0.013 0.013 0.013 0.002 0.002 0.002

# Multivariate dispersion
dispersion<-betadisper(spe.bray2, group=tempo$Treatment)
permutest(dispersion) # F = 0.0761, Df = 3,28, p = 0.977 
anova(dispersion)
plot(dispersion, hull=FALSE, ellipse = TRUE, label = FALSE, segments = FALSE, conf = 0.95) 
# Useful, but uses different method to NMDS. 

# TuxtMatrix <- as.data.frame(T5)
#write.csv(TuxtMatrix, "DataPerma.csv")

# TuxtMatrix2 <- as.data.frame(T7)
#write.csv(TuxtMatrix2, "DataPermaSinPlantados.csv")

y <- metaMDS(T5.2, k=3) # The number of reduced dimensions. 

Basal$Treatment <- relevel(Basal$Treatment, ref = "Wind") 
Basal$Treatment <- relevel(Basal$Treatment, ref = "Animal") 
Basal$Treatment <- relevel(Basal$Treatment, ref = "Forest") 

T6 <- as.data.frame.matrix(table(Basal$Family,Basal$Plot))
T6 <- t(T6[c(25,26,27,28,29,30,31,32,2,4,7,12,14,16,20,22,3,5,9,11,13,18,21,23,1,6,8,10,15,17,19,24)]) # Keeps Species ID & reorganizes for later separation of treatment with NMDS.
y <- metaMDS(T6, k=2) # The number of reduced dimensions. 
# y <- metaMDS(T7, k=2) # The number of reduced dimensions. 

tablespec <- table(Basal$Plot, Basal$Treatment)
tablespec$Plot <- tempo
tempo <- as.data.frame(as.factor(1:32))
colnames(tempo) <- c("Plot")
tempo$Treatment <- with(Basal, Treatment[match(tempo$Plot, Plot)]) # Match Treatment by Plot

######################################### MDS Without Planted Trees ######################################### 
tempo <- as.data.frame(as.factor(1:32))
colnames(tempo) <- c("Plot")
tempo$Treatment <- with(Basal, Treatment[match(tempo$Plot, Plot)]) # Match Treatment by Plot

tempo$Treatment <- relevel(tempo$Treatment, ref = "Wind") 
tempo$Treatment <- relevel(tempo$Treatment, ref = "Animal") 
tempo$Treatment <- relevel(tempo$Treatment, ref = "Forest") 

# Plot MDS
par(mfrow=c(1,1))
T6 <- t(as.data.frame.matrix(table(BasalW$Species,BasalW$Plot))) # Excludes planted trees! 

JaccardNMDS2<-metaMDS(T6, distance="jaccard", k=2, trymax=35, autotransform=TRUE, binary = FALSE, weakties = FALSE) # k is the number of dimensions
JaccardNMDS2 # Stress:     7.125718e-05 

new.pcoa <- pcoa(vegdist(T6,"jaccard", binary = FALSE))
biplot(new.pcoa, plot.axes=c(1,2), xlab = "", main = "")

JaccardVegX <- vegdist(T6, method = "jaccard", binary = FALSE) # This one. 

pcoa <- cmdscale(JaccardVegX, eig = TRUE)

col_vector <- tempo$Treatment
col_vector

col_palette <- c("#CC6666","#9999CC", "#66CC99", "peachpuff4")

ordiplot(pcoa, display = 'sites', text = 'sites')
ordiellipse(pcoa, groups = tempo$Treatment, draw = "polygon", lty = 1, col = c("peachpuff4","#9999CC", "#66CC99", "#CC6666"), pch = 16, conf = 0.95, kind = "se")
ordiellipse(pcoa, groups = tempo$Treatment, draw = "polygon", lty = 1, col = c("peachpuff4","#9999CC", "#66CC99", "#CC6666"), conf = 0.95, kind = "se", legend(0.4,-0.17, legend=unique(col_vector), col=unique(col_palette), pch = 16))

# PERMANOVA, Recruit-Only
T6 <- t(as.data.frame.matrix(table(BasalW$Species,BasalW$Plot))) # Excludes planted trees.
JaccardVeg <- vegdist(T6, method = "jaccard", binary = FALSE) # This one. 
adonis2(JaccardVeg~Treatment, data=tempo, permutations = 999)

# PERMANOVA, Planted-Trees Included
T5 <- t(as.data.frame.matrix(table(Basal$Species,Basal$Plot)))
JaccardVeg2 <- vegdist(T5, method = "jaccard", binary = FALSE) # This one. 
adonis2(JaccardVeg2~Treatment, data=tempo, permutations = 999)

# Pairwise Comparisons, Recruit-Only

HillQ <- as.data.frame.matrix(t(table(BasalW$Plot, BasalW$Species)))

div_test(HillQ,qvalue=0,hierarchy=tempo,posthoc=TRUE)
div_test(HillQ,qvalue=1,hierarchy=tempo,posthoc=TRUE)
div_test(HillQ,qvalue=2,hierarchy=tempo,posthoc=TRUE)

# Pairwise Comparisons, Planted Trees Included
HillQ2 <- as.data.frame.matrix(t(table(Basal$Plot, Basal$Species)))

div_test(HillQ2,qvalue=0,hierarchy=tempo,posthoc=TRUE)
div_test(HillQ2,qvalue=1,hierarchy=tempo,posthoc=TRUE)
div_test(HillQ2,qvalue=2,hierarchy=tempo,posthoc=TRUE)

# Cqn =  Sørensen-type overlap
# Uqn = Jaccard-type overlap
# Vqn = Sørensen-type turnover-complement
# SqN = Jaccard-type turnover-complement

######################################### GEOM-AREA  ###################################### 
# With planted species.
BasalAnimal <- droplevels(subset(Basal, Basal$Treatment == "Animal"))
BasalWind <- droplevels(subset(Basal, Basal$Treatment == "Wind"))
BasalControl <- droplevels(subset(Basal, Basal$Treatment == "Nat. Succ."))
BasalForest <- droplevels(subset(Basal, Basal$Treatment == "Forest"))

BasalFoX <- droplevels(subset(BasalForest, BasalForest$Species == "Astrocaryum mexicanum"))
summary(BasalFoX$Total.Area) # Actual Max is 50.30000

oop <- as.data.frame(summary(BasalAnimal$Family))
oop$Family <- row.names(oop)
names(oop)[1]<- "Frequency"
oop$Frequency <- as.numeric(oop$Frequency)
newdata <- oop[order(oop$Frequency),]
newdata

Top10Anim <- list("Fabaceae","Malvaceae","Piperaceae","Euphorbiaceae","Buseraceae","Apocynaceae","Moraceae","Melastomataceae","Urticaceae","Solanaceae")
Top10Wind <- list("Malvaceae","Apocynaceae","Burseraceae","Euphorbiaceae","Fabaceae","Piperaceae","Melastomataceae","Urticaceae","Boraginaceae","Vochysiaceae")
Top10Cont <- list("Burseraceae","Melastomataceae","Malvaceae","Piperaceae","Solanceae","Apocynaceae","Euphorbiaceae","Fabaceae","Urticaceae","Boraginaceae")
Top10Fore <- list("Arecaceae","Moraceae","Rubiaceae","Piperaceae","Euphorbiaceae","Sapotaceae","Lauraceae","Meliaceae","Clusiaceae","Urticaceae")

BasalAnimal$BinsD <- cut(x = na.omit(log(BasalAnimal$Total.Diameter)), breaks = c(-1,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5))
BasalAnimal$BinsDV <- cut(x = na.omit(log(BasalAnimal$Virtual.Diameter)), breaks = c(-0.5,0.05,0.6,1.15,1.7,2.25,2.8,3.35,3.9,4.45,5))
BasalAnimal$BinsD <- revalue(BasalAnimal$BinsD, c("(-1,0.5]"="I","(0.5,1]"="I", "(1,1.5]"="II", "(1.5,2]"="III", "(2,2.5]"="IV","(2.5,3]"="V", "(3,3.5]"="VI", "(3.5,4]"="VII",  "(4,4.5]"="VIII",  "(4.5,5]"="IX",  "(5,5.5]"="X"))
BasalAnimal$BinsDV <- revalue(BasalAnimal$BinsDV, c("(-0.5,0.05]"="I","(0.05,0.6]"="II", "(0.6,1.15]"="III", "(1.15,1.7]"="IV","(1.7,2.25]"="V", "(2.25,2.8]"="VI", "(2.8,3.35]"="VII",  "(3.35,3.9]"="VIII",  "(3.9,4.45]"="IX",  "(4.45,5]"="X"))

BasalWind$BinsD <- cut(x = na.omit(log(BasalWind$Total.Diameter)), breaks = c(-1,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5))
BasalWind$BinsDV <- cut(x = na.omit(log(BasalWind$Virtual.Diameter)), breaks = c(-0.5,0.05,0.6,1.15,1.7,2.25,2.8,3.35,3.9,4.45,5))
BasalWind$BinsD <- revalue(BasalWind$BinsD, c("(-1,0.5]"="I","(0.5,1]"="I", "(1,1.5]"="II", "(1.5,2]"="III", "(2,2.5]"="IV","(2.5,3]"="V", "(3,3.5]"="VI", "(3.5,4]"="VII",  "(4,4.5]"="VIII",  "(4.5,5]"="IX",  "(5,5.5]"="X"))
BasalWind$BinsDV <- revalue(BasalWind$BinsDV, c("(-0.5,0.05]"="I","(0.05,0.6]"="II", "(0.6,1.15]"="III", "(1.15,1.7]"="IV","(1.7,2.25]"="V", "(2.25,2.8]"="VI", "(2.8,3.35]"="VII",  "(3.35,3.9]"="VIII",  "(3.9,4.45]"="IX",  "(4.45,5]"="X"))

BasalControl$BinsD <- cut(x = na.omit(log(BasalControl$Total.Diameter)), breaks = c(-1,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5))
BasalControl$BinsDV <- cut(x = na.omit(log(BasalControl$Virtual.Diameter)), breaks = c(-0.5,0.05,0.6,1.15,1.7,2.25,2.8,3.35,3.9,4.45,5))
BasalControl$BinsDV <- revalue(BasalControl$BinsDV, c("(-0.5,0.05]"="I","(0.05,0.6]"="II", "(0.6,1.15]"="III", "(1.15,1.7]"="IV","(1.7,2.25]"="V", "(2.25,2.8]"="VI", "(2.8,3.35]"="VII",  "(3.35,3.9]"="VIII",  "(3.9,4.45]"="IX",  "(4.45,5]"="X"))
BasalControl$BinsDV <- revalue(BasalControl$BinsDV, c("(-1,0.5]"="I","(0.5,1]"="I", "(1,1.5]"="II", "(1.5,2]"="III", "(2,2.5]"="IV","(2.5,3]"="V", "(3,3.5]"="VI", "(3.5,4]"="VII",  "(4,4.5]"="VIII",  "(4.5,5]"="IX",  "(5,5.5]"="X"))

BasalForest$BinsD <- cut(x = na.omit(log(BasalForest$Total.Diameter)), breaks = c(-1,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5))
BasalForest$BinsDV <- cut(x = na.omit(log(BasalForest$Virtual.Diameter)), breaks = c(-0.5,0.1,0.7,1.3,1.9,2.5,3.1,3.7,4.3,4.9,5.5))
BasalForest$BinsD <- revalue(BasalForest$BinsD, c("(-1,0.5]"="I","(0.5,1]"="I", "(1,1.5]"="II", "(1.5,2]"="III", "(2,2.5]"="IV","(2.5,3]"="V", "(3,3.5]"="VI", "(3.5,4]"="VII",  "(4,4.5]"="VIII",  "(4.5,5]"="IX",  "(5,5.5]"="X"))
BasalForest$BinsDV <- revalue(BasalForest$BinsDV, c("(-0.5,0.1]"="I","(0.1,0.7]"="II", "(0.7,1.3]"="III", "(1.3,1.9]"="IV","(1.9,2.5]"="V", "(2.5,3.1]"="VI", "(3.1,3.7]"="VII",  "(3.7,4.3]"="VIII",  "(4.3,4.9]"="IX",  "(4.9,5.5]"="X"))

BasalAnimal$Top10Family <- BasalAnimal$Family
BasalWind$Top10Family <- BasalWind$Family
BasalControl$Top10Family <- BasalControl$Family
BasalForest$Top10Family <- BasalForest$Family

BasalAnimal$Top10Family <- recode(BasalAnimal$Top10Family,Fabaceae ="Fabaceae",Malvaceae="Malvaceae",Piperaceae="Piperaceae",Euphorbiaceae="Euphorbiaceae",Burseraceae="Burseraceae",Apocynaceae="Apocynaceae",Moraceae="Moraceae",Melastomataceae="Melastomataceae",Urticaceae="Urticaceae",Solanaceae="Solanaceae", .default = "Other")
BasalWind$Top10Family <- recode(BasalWind$Top10Family,Malvaceae="Malvaceae",Apocynaceae="Apocynaceae",Burseraceae="Burseraceae",Euphorbiaceae="Euphorbiaceae",Fabaceae="Fabaceae",Piperaceae="Piperaceae",Melastomataceae="Melastomataceae",Urticaceae="Urticaceae",Boraginaceae="Boraginaceae",Vochysiaceae="Vochysiaceae", .default = "Other")
BasalControl$Top10Family <- recode(BasalControl$Top10Family,Burseraceae="Burseraceae",Melastomataceae="Melastomataceae",Malvaceae="Malvaceae",Piperaceae="Piperaceae",Solanaceae="Solanaceae",Apocynaceae="Apocynaceae",Euphorbiaceae="Euphorbiaceae",Fabaceae="Fabaceae",Urticaceae="Urticaceae",Boraginaceae="Boraginaceae", .default = "Other")
BasalForest$Top10Family <- recode(BasalForest$Top10Family, Arecaceae="Arecaceae",Moraceae="Moraceae",Rubiaceae="Rubiaceae",Piperaceae="Piperaceae",Euphorbiaceae="Euphorbiaceae",Sapotaceae="Sapotaceae",Lauraceae="Lauraceae",Meliaceae="Meliaceae",Clusiaceae="Clusiaceae",Urticaceae="Urticaceae", .default = "Other")

summary(BasalForest$Top10Family)
summary(BasalControl$Top10Family)
summary(BasalWind$Top10Family)
summary(BasalAnimal$Top10Family)

# Virtual
AnimTableV <- table(BasalAnimal$Top10Family, BasalAnimal$BinsDV)
ContTableV <- table(BasalControl$Top10Family, BasalControl$BinsDV)
WindTableV <- table(BasalWind$Top10Family, BasalWind$BinsDV)
ForeTableV <- table(BasalForest$Top10Family, BasalForest$BinsDV)

AnimFamiV <- melt(AnimTableV, id.vars = c("Rank"), 
                  measure.vars = c("I","II","III","IV", "V",
                                   "VI", "VII", "VIII", "IX", "X"))
WindFamiV <- melt(WindTableV, id.vars = c("Rank"), 
                  measure.vars = c("I","II","III","IV", "V",
                                   "VI", "VII", "VIII", "IX", "X"))
ContFamiV <- melt(ContTableV, id.vars = c("Rank"), 
                  measure.vars = c("I","II","III","IV", "V",
                                   "VI", "VII", "VIII", "IX", "X"))
ForeFamiV <- melt(ForeTableV, id.vars = c("Rank"), 
                  measure.vars = c("I","II","III","IV", "V",
                                   "VI", "VII", "VIII", "IX", "X"))

Anim_PorcV <-  ddply(AnimFamiV, "Var.2", transform, Percent = value / sum(value) * 100)
WindFamiV <-  ddply(WindFamiV, "Var.2", transform, Percent = value / sum(value) * 100)
ContFamiV <-  ddply(ContFamiV, "Var.2", transform, Percent = value / sum(value) * 100)
ForeFamiV <-  ddply(ForeFamiV, "Var.2", transform, Percent = value / sum(value) * 100)

colnames(Anim_PorcV) <- c("Family", "DBHRank", "Frequency", "Percent") # Rename columns
colnames(WindFamiV) <- c("Family", "DBHRank", "Frequency", "Percent") # Rename columns
colnames(ContFamiV) <- c("Family", "DBHRank", "Frequency", "Percent") # Rename columns
colnames(ForeFamiV) <- c("Family", "DBHRank", "Frequency", "Percent") # Rename columns

# Reorder
Anim_PorcV$DBHRank <- relevel(Anim_PorcV$DBHRank, ref = "X") ; Anim_PorcV$DBHRank <- relevel(Anim_PorcV$DBHRank, ref = "IX") ; Anim_PorcV$DBHRank <- relevel(Anim_PorcV$DBHRank, ref = "VIII") ;Anim_PorcV$DBHRank <- relevel(Anim_PorcV$DBHRank, ref = "VII") ;Anim_PorcV$DBHRank <- relevel(Anim_PorcV$DBHRank, ref = "VI") ;Anim_PorcV$DBHRank <- relevel(Anim_PorcV$DBHRank, ref = "V") ;Anim_PorcV$DBHRank <- relevel(Anim_PorcV$DBHRank, ref = "IV") ;Anim_PorcV$DBHRank <- relevel(Anim_PorcV$DBHRank, ref = "III") ;Anim_PorcV$DBHRank <- relevel(Anim_PorcV$DBHRank, ref = "II") ;Anim_PorcV$DBHRank <- relevel(Anim_PorcV$DBHRank, ref = "I") 
WindFamiV$DBHRank <- relevel(WindFamiV$DBHRank, ref = "X") ; WindFamiV$DBHRank <- relevel(WindFamiV$DBHRank, ref = "IX") ; WindFamiV$DBHRank <- relevel(WindFamiV$DBHRank, ref = "VIII") ;WindFamiV$DBHRank <- relevel(WindFamiV$DBHRank, ref = "VII") ;WindFamiV$DBHRank <- relevel(WindFamiV$DBHRank, ref = "VI") ;WindFamiV$DBHRank <- relevel(WindFamiV$DBHRank, ref = "V") ;WindFamiV$DBHRank <- relevel(WindFamiV$DBHRank, ref = "IV") ;WindFamiV$DBHRank <- relevel(WindFamiV$DBHRank, ref = "III") ;WindFamiV$DBHRank <- relevel(WindFamiV$DBHRank, ref = "II") ;WindFamiV$DBHRank <- relevel(WindFamiV$DBHRank, ref = "I") 
ContFamiV$DBHRank <- relevel(ContFamiV$DBHRank, ref = "X") ; ContFamiV$DBHRank <- relevel(ContFamiV$DBHRank, ref = "IX") ; ContFamiV$DBHRank <- relevel(ContFamiV$DBHRank, ref = "VIII") ;ContFamiV$DBHRank <- relevel(ContFamiV$DBHRank, ref = "VII") ;ContFamiV$DBHRank <- relevel(ContFamiV$DBHRank, ref = "VI") ;ContFamiV$DBHRank <- relevel(ContFamiV$DBHRank, ref = "V") ;ContFamiV$DBHRank <- relevel(ContFamiV$DBHRank, ref = "IV") ;ContFamiV$DBHRank <- relevel(ContFamiV$DBHRank, ref = "III") ;ContFamiV$DBHRank <- relevel(ContFamiV$DBHRank, ref = "II") ;ContFamiV$DBHRank <- relevel(ContFamiV$DBHRank, ref = "I") 
ForeFamiV$DBHRank <- relevel(ForeFamiV$DBHRank, ref = "X") ; ForeFamiV$DBHRank <- relevel(ForeFamiV$DBHRank, ref = "IX") ; ForeFamiV$DBHRank <- relevel(ForeFamiV$DBHRank, ref = "VIII") ;ForeFamiV$DBHRank <- relevel(ForeFamiV$DBHRank, ref = "VII") ;ForeFamiV$DBHRank <- relevel(ForeFamiV$DBHRank, ref = "VI") ;ForeFamiV$DBHRank <- relevel(ForeFamiV$DBHRank, ref = "V") ;ForeFamiV$DBHRank <- relevel(ForeFamiV$DBHRank, ref = "IV") ;ForeFamiV$DBHRank <- relevel(ForeFamiV$DBHRank, ref = "III") ;ForeFamiV$DBHRank <- relevel(ForeFamiV$DBHRank, ref = "II") ;ForeFamiV$DBHRank <- relevel(ForeFamiV$DBHRank, ref = "I") 

Anim_PorcV$Family <- relevel(Anim_PorcV$Family, ref = "Moraceae"); 
Anim_PorcV$Family <- relevel(Anim_PorcV$Family, ref = "Melastomataceae");
Anim_PorcV$Family <- relevel(Anim_PorcV$Family, ref = "Malvaceae");
Anim_PorcV$Family <- relevel(Anim_PorcV$Family, ref = "Fabaceae");
Anim_PorcV$Family <- relevel(Anim_PorcV$Family, ref = "Euphorbiaceae");
Anim_PorcV$Family <- relevel(Anim_PorcV$Family, ref = "Burseraceae");
Anim_PorcV$Family <- relevel(Anim_PorcV$Family, ref = "Apocynaceae");

WindFamiV$Family <- relevel(WindFamiV$Family, ref = "Other"); 
WindFamiV$Family <- relevel(WindFamiV$Family, ref = "Melastomataceae"); 
WindFamiV$Family <- relevel(WindFamiV$Family, ref = "Malvaceae");
WindFamiV$Family <- relevel(WindFamiV$Family, ref = "Fabaceae");
WindFamiV$Family <- relevel(WindFamiV$Family, ref = "Euphorbiaceae");
WindFamiV$Family <- relevel(WindFamiV$Family, ref = "Burseraceae");
WindFamiV$Family <- relevel(WindFamiV$Family, ref = "Boraginaceae");
WindFamiV$Family <- relevel(WindFamiV$Family, ref = "Apocynaceae");

ContFamiV$Family <- relevel(ContFamiV$Family, ref = "Other"); 
ContFamiV$Family <- relevel(ContFamiV$Family, ref = "Melastomataceae"); 
ContFamiV$Family <- relevel(ContFamiV$Family, ref = "Malvaceae");
ContFamiV$Family <- relevel(ContFamiV$Family, ref = "Fabaceae");
ContFamiV$Family <- relevel(ContFamiV$Family, ref = "Euphorbiaceae");
ContFamiV$Family <- relevel(ContFamiV$Family, ref = "Burseraceae");
ContFamiV$Family <- relevel(ContFamiV$Family, ref = "Boraginaceae");
ContFamiV$Family <- relevel(ContFamiV$Family, ref = "Apocynaceae");

ForeFamiV$Family <- relevel(ForeFamiV$Family, ref = "Other");
ForeFamiV$Family <- relevel(ForeFamiV$Family, ref = "Moraceae");
ForeFamiV$Family <- relevel(ForeFamiV$Family, ref = "Meliaceae");
ForeFamiV$Family <- relevel(ForeFamiV$Family, ref = "Lauraceae");
ForeFamiV$Family <- relevel(ForeFamiV$Family, ref = "Euphorbiaceae");
ForeFamiV$Family <- relevel(ForeFamiV$Family, ref = "Clusiaceae");
ForeFamiV$Family <- relevel(ForeFamiV$Family, ref = "Arecaceae");

AnimArea <- ggplot(Anim_PorcV, aes(x= DBHRank, y= Percent, group = Family, fill = Family)) + xlab("Diameter Rank") + ylab("Proportional Abundance") + ggtitle("(B) Animal-Dispersed Plantings") + geom_area(position = "fill") +coord_flip() + theme_bw(base_size = 14)+ 
  scale_fill_manual(values=c("brown4","#D53E4F","#F46D43","#FDAE61","#FEE08B","#FFFFBF","#E6F598","#ABDDA4","#66C2A5","#3288BD","#5E4FA2")) 

WindArea <- ggplot(WindFamiV, aes(x= DBHRank, y= Percent, group = Family, fill = Family))+ xlab("Diameter Rank") + ylab("Proportional Abundance")+ ggtitle("(C) Wind-Dispersed Plantings")  + geom_area(position = "fill") +coord_flip() +  theme_bw(base_size = 14)+
  scale_fill_manual(values=c("brown4","#9E0142","#D53E4F","#F46D43","#FDAE61","#FEE08B","#FFFFBF","#ABDDA4","#66C2A5","#5E4FA2","slateblue1")) 

ContArea <- ggplot(ContFamiV, aes(x= DBHRank, y= Percent, group = Family, fill = Family))+ xlab("Diameter Rank") + ylab("Proportional Abundance") + ggtitle("(D) Natural Succession") + geom_area(position = "fill") +coord_flip() + theme_bw(base_size = 14)+
  scale_fill_manual(values=c("brown4","#9E0142","#D53E4F","#F46D43","#FDAE61","#FEE08B","#FFFFBF","#ABDDA4","#66C2A5","#3288BD","#5E4FA2")) 

ForeArea <- ggplot(ForeFamiV, aes(x= DBHRank, y= Percent, group = Family, fill = Family))+ xlab("Diameter Rank") + ylab("Proportional Abundance") + ggtitle("(A) Forest") + geom_area(position = "fill") +coord_flip() + theme_bw()+  theme_bw(base_size = 14)+
  scale_fill_manual(values=c("coral4","orange3","#F46D43","darkolivegreen4","olivedrab3","#E6F598","#ABDDA4","#66C2A5","pink1","plum3","#5E4FA2")) 

pushViewport(viewport(layout=grid.layout(nrow=2,ncol=2))) # This sets the multiplot.
print(ForeArea, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(AnimArea, vp=viewport(layout.pos.row=1, layout.pos.col=2))
print(WindArea, vp=viewport(layout.pos.row=2, layout.pos.col=1))
print(ContArea, vp=viewport(layout.pos.row=2, layout.pos.col=2))

ggplot(Basal, aes(x=Plot, y=Total.Diameter, fill=Treatment)) + 
  geom_boxplot(show.legend = FALSE) +
  xlab("Plot") +
  ylab(expression(Total~Diameter~(cm^2))) + 
  theme(axis.title.x = element_text(size = 10, color = "white")) + 
  theme(axis.text.x = element_text(size = 5, color = "black")) + 
  theme(axis.text.y = element_text(size = 7.5, color = "black")) +
  scale_fill_manual(values=c("peachpuff4", "#9999CC","#66CC99", "#CC6666")) +
  facet_grid(~Treatment, scales = "free") +
  theme(legend.position="none") + theme_bw(base_size = 20) 

ggplot(Basal, aes(x=Plot, y=log(Total.Diameter), fill=Treatment)) + 
  geom_boxplot(show.legend = FALSE) +
  xlab("Plot") +
  ylab(expression(Log~Total~Diameter~(cm^2))) + 
  theme(axis.title.x = element_text(size = 10, color = "white")) + 
  theme(axis.text.x = element_text(size = 5, color = "black")) + 
  theme(axis.text.y = element_text(size = 7.5, color = "black")) +
  scale_fill_manual(values=c("peachpuff4", "#9999CC","#66CC99", "#CC6666")) +
  facet_grid(~Treatment, scales = "free") +
  theme(legend.position="none") + theme_bw(base_size = 20) 

############################################# Extra Stuff ################################

Moraceae <- droplevels(subset(BasalAnimal, BasalAnimal$Family == "Moraceae")) 
table(Moraceae$Species, Moraceae$Origin)

tgc <- summarySE(BasalAnimal, measurevar="Total.AreaMeters", groupvars=c("Species"), na.rm = TRUE)
20.13942466 / 77.11512 * 100
sum(tgc$Total.AreaMeters) 

614/1511*100 # 40.64% of all individuals.

# How many trees of each planted species

BasalP <- droplevels(subset(Basal, Basal$Origin == "Planted")) # Exclude Recruits
BasalP$Species <- as.factor(BasalP$Species)
Exp <- as.table(summary(BasalP$Species))

table(BasalP$Species, BasalP$Total.Area)

Exp2 <- aggregate(BasalP$Total.Area/100, by=list(Category=BasalP$Species, BasalP$Treatment), FUN=sum, na.rm="TRUE")

Exp3 <- droplevels(subset(Exp2, Exp2$Group.2 == "Animal")) # Exclude Planted trees
Exp4 <- droplevels(subset(Exp2, Exp2$Group.2 == "Wind")) # Exclude Planted trees

# Define the number of colors you want
nb.cols <- 12
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)
# Create a ggplot with 18 colors 
# Use scale_fill_manual

n1 <- nrow(Exp3)                                        # Amount of default colors
hex_codes1 <- hue_pal()(n1)                             # Identify hex codes
hex_codes1
show_col(hex_codes1)                                    # Plot hex codes/colors

Exp3$Category <- revalue(Exp3$Category, c("Guarea grandifolia"="Guarea glabra"))

SPA <- ggplot(Exp3, aes(fill=Category, y=x, x=Group.2)) + 
  geom_bar(position="stack", stat="identity", color = "Black", width = 1) +
  scale_fill_manual(values=c("#66C2A5", "#C5A07A", "#DD927E", "#979EC1", "#BE93C6",
                             "#DB98AE", "#B1C968", "#CED843", "#FCD738", "#ECC978", "#D2BD9F",
                             "#B3B3B3")) +
  scale_y_continuous(label=comma, limits = c(0, 2650), breaks = seq(0, 2650, by = 500)) +
  labs(fill = "Planted Species", y = 'Total Basal Area' ~(cm^2)) + 
  theme_bw(base_size = 15) +
  theme(axis.title.x=element_blank(), legend.text=element_text(size=10))
SPA

mycolors <- colorRampPalette(brewer.pal(8, "Set3"))(nb.cols)

SPW <- ggplot(Exp4, aes(fill=Category, y=x, x=Group.2)) + 
  geom_bar(position="stack", stat="identity", color = "Black", width = 1) +
  scale_fill_manual(values=c("#8DD3C7","#D5EFBA","#EDECBD","#C3C0D6","#DF9AA1","#E48883","#96A8C1","#B8B29F","#F6B762","#8a914d","#CDD796","#FCCDE5")) +
  scale_y_continuous(label=comma, limits = c(0, 2650), breaks = seq(0, 2650, by = 500)) +
  labs(fill = "Planted Species", y = 'Total Basal Area' ~(cm^2)) + 
  theme_bw(base_size = 15) +
  theme(axis.title.x=element_blank(), legend.text=element_text(size=10))
SPW

pushViewport(viewport(layout=grid.layout(nrow=1,ncol=2))) # This sets the multiplot.
print(SPA, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(SPW, vp=viewport(layout.pos.row=1, layout.pos.col=2))

Exp5 <- aggregate(BasalControl$Total.Area/100, by=list(Category=BasalControl$Species), FUN=sum, na.rm="TRUE")

Exp5<- top_n(Exp5, 12, x)  

Exp5$Group.2 <- "Nat. Succ."
mycolors <- colorRampPalette(brewer.pal(8, "Accent"))(nb.cols)

SPC <- ggplot(Exp5, aes(fill=Category, y=x , x=Group.2)) + 
  geom_bar(position="stack", stat="identity", color = "Black", width = 1) +
  scale_fill_manual(values=c("#8DD3C7","#A7B7B5","#CFB2BE","#979EC1","#FEE290","#DAE49D","#5C86AB","#96659e","#d96ca0","#e67063","#9E5F33","#666666")) +
  labs(fill = "12 Largest Species", y = 'Total Basal Area' ~(cm^2)) + 
  scale_y_continuous(label=comma, limits = c(0, 2650), breaks = seq(0, 2650, by = 500)) +
  theme_bw(base_size = 15) +
  theme(axis.title.x=element_blank(), legend.text=element_text(size=10))
SPC

Exp6 <- aggregate(BasalForest$Total.Area/100, by=list(Category=BasalForest$Species), FUN=sum, na.rm="TRUE")
Exp6<- top_n(Exp6, 12, x)  
Exp6$Group.2 <- "Forest"

mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols)

SPF <- ggplot(Exp6, aes(fill=Category, y=x , x=Group.2)) + 
  geom_bar(position="stack", stat="identity", color = "Black", width = 1) +
  scale_fill_manual(values=c("#ed6b6c","#DD927E","#3D8B9A","#4AAA54","#B1C968","#AA5685","#EC761D","#c7a861","#ECC978","#BE842A","#C3655E","#F781BF")) +
  scale_y_continuous(label=comma, limits = c(0, 2650), breaks = seq(0, 2650, by = 500)) +
  labs(fill = "12 Largest Species", y = 'Total Basal Area' ~(cm^2), x = 'Habitat') +  
  theme_bw(base_size = 15) +
  theme(axis.title.x=element_blank(), legend.text=element_text(size=10))
SPF

plot_grid(SPA, SPW, SPC, SPF, ncol=2, nrow=2, align="v", labels=LETTERS[1:4])

