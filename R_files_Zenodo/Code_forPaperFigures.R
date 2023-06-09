# Libraries required
library("lmerTest")
library("lme4")
library("lmerTest")
library("DHARMa")
library("glmmTMB")
library("vegan")
library("ggplot2")

#Import Data

#Convert .xlsx data to .csv file from Dryad, for file "Herbarium stigma sample pollen counts"
ATSso.h <- read.csv() #Convert .xlsx data to .csv file from Dryad, for file "Herbarium stigma sample pollen counts" and upload here

## Run models for conspecific pollen, heterospecific pollen richness, and heterospecific pollen proportion

#### negative binomial zero inflated model

#CP#

CP.count <- glmmTMB(Conspecific ~ TimeCat + (1|Spp), data=ATSsp.h, ziformula=~1,family="nbinom2") 
null.CP <-  glmmTMB(Conspecific ~ (1|Spp), data=ATSsp.hx2, ziformula=~1,family="nbinom2") 

summary(CP.count)
anova(CP.count, null.CP) #no difference between null and real model; it's not significant

#Does it fit assumptions of GLMM TMB Zero-Inflated Negative Binomial model
simulationOutput1 <- simulateResiduals(fittedModel = CP.count, n = 1000)
plotSimulatedResiduals(simulationOutput = simulationOutput1) #Plot of QQ plot looks good: residual vs. predicted quantile lines won't work
testUniformity(simulationOutput = simulationOutput1) #D = 0.0835, p-value = 0.2697
testDispersion(simulationOutput1) #ratioObsSim = 0.95424, p-value = 0.68
testZeroInflation(simulationOutput1) #Zero inflation test yields ratio ObsSim = 0.91848, p-value = 0.864

#HP Richness#

HP.rich <- glmmTMB(HPRich ~ TimeCat + (1|Spp), data=ATSsp.hx2, ziformula=~1,family="nbinom2") 
null.HPr <-  glmmTMB(HPRich ~ (1|Spp), data=ATSsp.hx2, ziformula=~1,family="nbinom2") 

summary(HP.rich)
anova(CP.count, null.CP) #no difference between null and real model; it's not significant

#Does it fit assumptions of GLMM TMB Zero-Inflated Negative Binomial model
simulationOutput1 <- simulateResiduals(fittedModel = HP.rich, n = 1000)
plotSimulatedResiduals(simulationOutput = simulationOutput1) #Plot of QQ plot looks good: residual vs. predicted quantile lines won't work
testUniformity(simulationOutput = simulationOutput1) #D = 0.045667, p-value = 0.9248
testDispersion(simulationOutput1) #ratioObsSim = 0.88676, p-value = 0.49
testZeroInflation(simulationOutput1) #Zero inflation test yields ratio ObsSim = 0.99, p-value = 1


#HP Proportion#

HP.prop <- glmmTMB(HPProp ~ TimeCat + (1|Spp), data=ATSsp.hx2, ziformula=~1,family="tweedie")
null.HPp <-  glmmTMB(HPProp ~ (1|Spp), data=ATSsp.hx2, ziformula=~1,family="tweedie") 

summary(HP.prop)
anova(HP.prop, null.HPp) #no difference between null and real model; it's not significant

#Does it fit assumptions of GLMM TMB Zero-Inflated Negative Binomial model
simulationOutput1 <- simulateResiduals(fittedModel = HP.prop, n = 1000)
plotSimulatedResiduals(simulationOutput = simulationOutput1) #Plot of QQ plot looks good: residual vs. predicted quantile lines won't work
testUniformity(simulationOutput = simulationOutput1) #D = 0.045667, p-value = 0.9248
testDispersion(simulationOutput1) #ratioObsSim = 0.88676, p-value = 0.49
testZeroInflation(simulationOutput1) #Zero inflation test yields ratio ObsSim = 0.99, p-value = 1

### Make the plots for the paper, in Figure 2

hprich <- ggplot(ATSsp.h, aes(TimeCat, HPRich)) + geom_boxplot() + 
  geom_dotplot(aes(color=Spp, fill=Spp), shape=Spp, binaxis="y", stackdir="center",dotsize=0.6, position=position_dodge(0.6)) + 
  ylab("Heterospecific pollen richness") + xlab("Time category") + theme_classic() + theme(text=element_text(size=15))

CP <- ggplot(ATSsp.h, aes(TimeCat, Conspecific)) + 
  geom_boxplot() + geom_jitter(aes(color=Spp, pch=Spp)) + ylim(0, 500) + 
  ylab("Number of conspecific pollen grains") + xlab("Time category") + 
  theme_classic() + theme(text=element_text(size=15))

hpprop <- ggplot(ATSsp.h, aes(TimeCat, HPProp)) + geom_boxplot() + geom_jitter(aes(color=Spp, pch=Spp)) + 
  ylab("Heterospecific pollen proportion") + xlab("Time category") + theme_classic() + 
  theme(text=element_text(size=15))

hprich <- ggplot(ATSsp.h, aes(TimeCat, HPRich)) + geom_boxplot() + 
  geom_jitter(aes(color=Spp, pch=Spp), height=0) +  ylab("Heterospecific pollen richness") +
  xlab("Time category") + theme_classic() + theme(text=element_text(size=15))

# Conduct compositional ordination in Figure 3
SPList <- c("DOVI","MEPO","MYSA","VARE")
Comp2 <- ATSsp.h[ATSsp.h$Code %in% SPList,]#select out only the four species with large enough sample size

dist <- vegdist(Comp2[,11:41], method="bray")
pcoa <- pcoa(dist)
plot(pcoa$vectors[,1], pcoa$vectors[,2], type="n")
text(pcoa$vectors[,1], pcoa$vectors[,2])

pcoa2 <- cmdscale(dist, eig=TRUE)

x <- pcoa2$points[,1]
y <- pcoa2$points[,2]
plot(x,y, xlab="PCoA axis 1", ylab="PCoA axis 2", type="n")
text(x, y, labels=Comp$Spp, cex=0.7)

pcoa <- as.data.frame(pcoa2$points)

var1 <- betadisper(vegdist(Comp2[,11:41], method="bray"), group=Comp2$SppTime)
permutest(var1) #no diff

var2 <- betadisper(vegdist(Comp2[,11:41], method="bray"), group=Comp2$Spp)
permutest(var2) #no diff


ad <- adonis(Comp2[,11:41]~Comp2$TimeCat, strata=Comp2$Spp)

#Create Figure 3
centroids <- as.data.frame(var1$centroids[,1:2])
centroids$Spp <- rep(c("DOVI","MEPO","MYSA","VARE"), each=3)
centroids$Time <- rep(c("Post1950","Pre1950"), 4)

ord2 <- ggplot(pcoa)

ord2 + geom_point(aes(V1, V2, colour=Comp2$Spp), alpha=0.5) + stat_ellipse(aes(V1,V2, colour=Comp2$Spp), type="norm", linetype=2) + geom_point(data=centroids, aes(PCoA1,PCoA2, colour=Spp, shape=Time),size=3) + theme_bw() + labs(x="PCoA1 28% of variation", y="PCoA2 30% of variation", fill="Sample time period")


