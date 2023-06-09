# Data Analysis Feb 9, 2019

# ========  Housekeeping  ===========
library(ggplot2)
library(data.table)
library(scales)
library(plyr)
library(dplyr)
library(reshape2)
library(lme4)
library(lmerTest)
library(cowplot)
library(car)
library(effects)
library(OpenStreetMap)
library(maps)
library(ggrepel)
library(sjPlot)
library(glmmTMB)
library(DHARMa)
library(emmeans)

powerTransform <- function(y, lambda1, lambda2 = NULL, method = "boxcox") {
  
  boxcoxTrans <- function(x, lam1, lam2 = NULL) {
    
    # if we set lambda2 to zero, it becomes the one parameter transformation
    lam2 <- ifelse(is.null(lam2), 0, lam2)
    
    if (lam1 == 0L) {
      log(y + lam2)
    } else {
      (((y + lam2)^lam1) - 1) / lam1
    }
  }
  
  switch(method
         , boxcox = boxcoxTrans(y, lambda1, lambda2)
         , tukey = y^lambda1
  )
}

latData <- read.csv("~/Downloads/longLat.csv")

# ===================================================
# ========  Preference Data ===================
WL <- read.csv("WinLose21Jan2019.csv")
WL$Latitude <- latData[match(WL$Genotype, latData$Genotype),"Latitude"]
WL$Longitude <- latData[match(WL$Genotype, latData$Genotype),"Longitude"]
WL$Population <- latData[match(WL$Genotype, latData$Genotype),"Population"]
WL$Population <- as.factor(WL$Population)
# =========== WL models =======
mod8 <- glmer(Preferred ~  Range+(1|Population)+(1 |Genotype) + (1|Pair), family = binomial(link = "logit"), data = WL)

# ====== Model Diagnostics
sr1 <- simulateResiduals(mod8)
testResiduals(sr1)
plot(sr1)
leveneTest(residuals(mod8) ~ WL$Range)
drop1(mod8, test = "Chisq")
# qqplot plus KS test indicate residuals are uniform, not overdispersed, and there are no significant outliers

# ==========================================================
# =========== No Choice Data ===========
NoChoice <- read.csv("NoChoice1Feb2019.csv")

# Remove obvious outliers
NoChoice <- subset(NoChoice, comsumed > 0)
NoChoice <- subset(NoChoice, comsumed < 1)

levels(NoChoice$Genotype)[levels(NoChoice$Genotype)=="Mt.Wilson_3"] <- "Mt.Wilson_2"
levels(NoChoice$Genotype)[levels(NoChoice$Genotype)=="ChilePI368930"] <- "ChilePI368950"

#NoChoice[!complete.cases(NoChoice),] # Check to fix naming problem


NoChoice$Latitude <- latData[match(NoChoice$Genotype, latData$Genotype),"Latitude"]
NoChoice$Longitude <- latData[match(NoChoice$Genotype, latData$Genotype),"Longitude"]
NoChoice$Population <- latData[match(NoChoice$Genotype, latData$Genotype),"Population"]
NoChoice$Population <- as.factor(NoChoice$Population)

# ======== No Choice model ======
NCbetamod <- glmmTMB(comsumed ~ Range + (1|Genotype) + (1| Population), family = beta_family(), NoChoice)

sr2 <- simulateResiduals(NCbetamod)
testResiduals(sr2)
plot(sr2)
leveneTest(residuals(NCbetamod) ~ NoChoice$Range) # heteroskedasticity!

#according to levene's test, there is heterosckedascity, to include in model according to https://stats.stackexchange.com/questions/213591/how-can-i-test-for-differences-in-variation-between-groups-in-a-mixed-model-lme/214007#214007
NoChoice <- transform(NoChoice,obs=factor(1:nrow(NoChoice)), Rangedummy=as.numeric(Range=="Native"))


NCbetamodb <- glmmTMB(comsumed~Range + (1|Genotype) + (1| Population) + (Rangedummy-1|obs), family = beta_family(),  NoChoice)

sr3 <- simulateResiduals(NCbetamodb)
testResiduals(sr3)
plot(sr3)
leveneTest(residuals(NCbetamodb) ~ NoChoice$Range)

anova(NCbetamod, NCbetamodb)

Anova(NCbetamodb)
drop1(NCbetamodb, test = "Chisq")


# Redoing
pairComp <- read.csv("Results/ProcessedData/compAlonePairedZ.csv")
# Remove obvious outliers
pairComp <- subset(pairComp, proprotionTotEate > 0)
pairComp <- subset(pairComp, proprotionTotEate < 1)

# summarySE(pairComp, measurevar = "proprotionTotEate", groupvars = c("Type", "group"))
# pairComp$Population <- latData[match(pairComp$Genotype, latData$Genotype),"Population"]
# pairComp$Population <- as.factor(pairComp$Population)
# 
# NCbetamodB <- glmmTMB(proprotionTotEate ~ Type * group + (1|Genotype) + (1|Population), family = beta_family(), pairComp)
# 
# leveneTest(residuals(NCbetamodB) ~ pairComp$Type)
# 
# pairComp <- transform(pairComp,obs=factor(1:nrow(pairComp)), Typedummy=as.numeric(Type=="Native"))
# 
# NCbetamodB2 <- glmmTMB(proprotionTotEate ~ Type * group + (1|Genotype) + (1|Population)+ (Typedummy-1|obs), family = beta_family(), pairComp)
# anova(NCbetamodB, NCbetamodB2)
# 
# Anova(NCbetamodB2)
# (ae <- allEffects(NCbetamodB))
# plot(ae)
# 
# NCbetamodC <- glmmTMB(propStd ~ Type * group + (1|Genotype), family = beta_family(), pairComp)
# Anova(NCbetamodC)
# summary(C.glmer)
# (ae <- allEffects(NCbetamodC))
# plot(ae)
# 
# joint_tests(NCbetamodB, by = "group")
# joint_tests(NCbetamodC, by = "group")


pairComp2 <- subset(pairComp, Genotype != "GILNS13_6")

NCbetamodC <- glmmTMB(proprotionTotEate ~ Type * group + (1|Genotype) + (1|Population), family = beta_family(), pairComp2)

leveneTest(residuals(NCbetamodC) ~ pairComp2$Type)

pairComp2 <- transform(pairComp2,obs=factor(1:nrow(pairComp2)), Typedummy=as.numeric(Type=="Native"))

NCbetamodC2 <- glmmTMB(proprotionTotEate ~ Type * group + (1|Genotype) + (1|Population)+ (Typedummy-1|obs), family = beta_family(), pairComp2)
anova(NCbetamodC, NCbetamodC2)

Anova(NCbetamodC2)



joint_tests(NCbetamodC, by = "group")
# ============= Insect Performance Data =========
InsectData <- read.csv("InsectsHerbData.csv")
InsectData$removeV <- InsectData$estFinDry/InsectData$estInitDry
levels(InsectData$Genotype)[levels(InsectData$Genotype)=="Mt.Wilson_3"] <- "Mt.Wilson_2"
InsectData$Latitude <- latData[match(InsectData$Genotype, latData$Genotype),"Latitude"]
InsectData$Longitude <- latData[match(InsectData$Genotype, latData$Genotype),"Longitude"]
InsectData$Population <- latData[match(InsectData$Genotype, latData$Genotype),"Population"]
InsectData$Death <- InsectData$InitialLooperNumber - InsectData$FinalLooperNumber

InsectData$Population <- as.factor(InsectData$Population)

Sd <- sd(InsectData$removeV)
Mn <- mean(InsectData$removeV)
IGsubset <- subset(InsectData, removeV > (Mn-2*Sd))  #(From 138 to 135)
# ============== Growth model =========
growM3 <- lmer(GrowthRate ~ Type + (1 | Genotype) + (1|Population), IGsubset)

sr3 <- simulateResiduals(growM3)
testResiduals(sr3)
plot(sr3)
leveneTest(residuals(growM3) ~ IGsubset$Type)

Anova(growM3)

# ============== Death models =========
deathmod7 <- glmer.nb(Death ~ Type + (1|Genotype) + (1|Population), IGsubset, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000000)))

sr4 <- simulateResiduals(deathmod7)
testResiduals(sr4)
plot(sr4) # model looks to be overdispersed
testZeroInflation(sr4)

deathmod6b <- glmmTMB(Death ~ Type + (1|Genotype) + (1|Population), family = poisson(), IGsubset) 
Anova(deathmod6b)

deathmod7c <- glmer(Death ~ Type + (1|Genotype) + (1|Population), family = "poisson", IGsubset)
Anova(deathmod7)

# April 10, Mort data redone =========
# using current experiment data:
deathmod8 <- glmmTMB(FinalLooperNumber/InitialLooperNumber ~ Type + (1|Genotype) + (1|Population), weights = InitialLooperNumber, family = "binomial", IGsubset)
Anova(deathmod8)

# using data from April 24-30, 2015
AprDeath <- read.csv("Results/ProcessedData/MortDatafromApr2015.csv")
deathmod9 <- glmmTMB(WeekAlive/WeekCases ~ Range + (1|Genotype) + (1|Population), weights = WeekCases, family = "binomial", AprDeath)
Anova(deathmod9)
ggplot(AprDeath, aes(x = Range, y = 1-WeekAlive/WeekCases)) + stat_summary(fun.data = "mean_se")


# Figures

# Figure 2: map
world_map <- map_data("world")
world_map$soybean <- ifelse(world_map$region == "Canada" | world_map$region == "USA"| world_map$region == "Mexico"| world_map$region == "Belize"| world_map$region == "El Salvador"| world_map$region == "Guatemala"|world_map$region == "Honduras"|world_map$region == "Nicaragua"|world_map$region == "Panama"|world_map$region == "Argentina"|world_map$region == "Bolivia"|world_map$region == "Brazil"|world_map$region == "Chile"|world_map$region == "Colombia"|world_map$region == "Ecuador"|world_map$region == "Guyana"|world_map$region == "Paraguay"|world_map$region == "Suriname"|world_map$region == "Uruguay"|world_map$region == "Venezuela"|world_map$region == "Australia"|world_map$region == "Peru", "yes", "no")

myColors3 <- c("white", "snow3")
class(myColors3)
names(myColors3) <- levels(world_map$soybean)
colScaleA <- scale_fill_manual(name = "soybean", values = myColors3)

p <- ggplot() + coord_equal() + xlab("") + ylab("")
base_world_messy <- p + geom_polygon(data=world_map, aes(x=long, y=lat, fill= soybean, group=group), colour="black") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = 'white', colour = 'white'), axis.line = element_line(colour = "white"), legend.position="none", axis.ticks=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank()) + colScaleA

base_world_messy

myColors4 <- c("red", "blue")
class(myColors4)
names(myColors4) <- levels(latData$PopGroup)
colScaleB <- scale_colour_manual(name = "PopGroup", values = myColors4)
base_world_messy + colScaleB +geom_label_repel(data = latData, aes(x = Longitude, y = Latitude, label = MaoID, colour = PopGroup), fontface = "bold", nudge_x = c(1, -1.5, 2, 2, -1), nudge_y = c(0.25, -0.25, 0.5, 0.5, -0.5), segment.colour = "black", min.segment.length = 0, segment.size = 0.75, force = 1) + geom_point(data = latData, aes(x=Longitude, y=Latitude, colour = PopGroup)) 

ggsave("MapGenotypes30Apr2019.jpeg")

#  ===== Figure 3
Fig3A <- ggplot(WL, aes(Range, Preferred)) + stat_summary(fun.data = "mean_se") + ylim(0,1) + theme(axis.title.x = element_blank()) + labs(y = "Preferrence") #+ theme_cowplot()

Fig3B <-ggplot(pairComp, aes(Type, proprotionTotEate, fill = group))   + theme(axis.title.x = element_blank(), legend.position = c(.7,.8)) + labs(y = "Proportion \\nconsumed") + geom_boxplot() + scale_fill_grey(start =.3, end = 0.8, name = "Experiment", labels = c("No Choice", "Choice")) 

Fig3C <- ggplot(IGsubset, aes(Type, GrowthRate)) + geom_boxplot()  + theme(axis.title.x = element_blank()) + labs(y = "RGR") #+ theme_cowplot()

Fig3D <- ggplot(IGsubset, aes(Type, DR)) + geom_boxplot()  + theme(axis.title.x = element_blank()) + labs(y = "Mortality rate")#+ theme_cowplot()

Fig3 <- cowplot::plot_grid(Fig3A, Fig3C, Fig3D, Fig3B, labels = "AUTO", label_size = 18, nrow = 2, align = "hv", axis = "l")

cowplot::save_plot("Fig3Herb24Apr2019.pdf", Fig3, base_width = 10, base_height = 6)

ggsave("Fig3Herb30Apr2019.jpeg", Fig3)
# Responding to reviewer comments again (April 5, 2019)

# mortality data- 
str(pairComp)
unique(pairComp$Genotype)
# Just looking over utilization worksheet
Udat <- read.csv("~/Documents/Projects/RapidEvolution/Results/ProcessedData/Book4.csv")
str(Udat)

# Growth Rate
Anova(lme4::lmer(Growth.Rate ~ Type + (1| Genotype), data = Udat))
modu1 <- lme4::lmer(Growth.Rate ~ Type + (1| Genotype), data = Udat)

sr1 <- simulateResiduals(modu1)
testResiduals(sr1)
plot(sr1)

leveneTest(residuals(modu1) ~ modu1@frame[["Type"]])
ggplot(Udat, aes(x = Type, y = Growth.Rate)) + geom_boxplot()


# Efficiency.Conversion.Ingested.Food
Anova(lme4::lmer(Consumption.Indices ~ Type + (1| Genotype), data = Udat))
modu1 <- lme4::lmer(Efficiency.Conversion.Digested.Food ~ Type + (1| Genotype), data = Udat)

sr1 <- simulateResiduals(modu1)
testResiduals(sr1)
plot(sr1)

min(Udat$Efficiency.Conversion.Ingested.Food)
CRL.bc <- boxcox(Efficiency.Conversion.Ingested.Food + 1.83 ~ Type, data=Udat)
(lambda <- CRL.bc$x[which.max(CRL.bc$y)])
Mod2.CRL <-  lme4::lmer(powerTransform((Efficiency.Conversion.Ingested.Food + 1.83), lambda) ~ Type + (1| Genotype), data = Udat)

sr1 <- simulateResiduals(Mod2.CRL)
testResiduals(sr1)
plot(sr1)

Anova(Mod2.CRL)
ggplot(Udat, aes(x = Type, y = Efficiency.Conversion.Ingested.Food)) + geom_boxplot()
ddply(Udat, "Type", summarise, MeanGro = mean(Growth.Rate))
ddply(Udat, "Type", summarise, MeanGro = mean(Efficiency.Conversion.Ingested.Food))
