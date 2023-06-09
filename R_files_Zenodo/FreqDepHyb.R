#set path to directory containing data files
#setwd('~/Downloads/')

version #R version 4.1.2 (2021-11-01)
library(ggplot2) #packageVersion("ggplot2") ‘3.3.5’
library(patchwork) #packageVersion("patchwork") ‘1.1.1’
library(glmmTMB) #packageVersion("glmmTMB") ‘1.1.2.3’
library(DHARMa) #packageVersion("DHARMa") ‘0.4.4’
library(ggeffects) #packageVersion("ggeffects") ‘1.1.1’
library(emmeans) #packageVersion("emmeans") ‘1.7.2’
library(multcomp) #packageVersion("multcomp") ‘1.4.18’

####Overdispersion function from 
#http://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#overdispersion
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

##Observational Transect##
#####
transectlong<-read.csv('Transect2019.csv')

cor.test(transectlong[transectlong$Species=="guttatus",]$Freq, transectlong[transectlong$Species=="guttatus",]$Prop, method="pearson")
cor.test(transectlong[transectlong$Species=="nudatus",]$Freq, transectlong[transectlong$Species=="nudatus",]$Prop, method="pearson")

p1<-ggplot(transectlong, aes(Position,Freq),fill=Species)+
  geom_area(aes(fill=Species),position="stack",alpha=0.7)+
  scale_fill_manual(values = c("#ca0020","#0571b0"))+
  geom_point(data=transectlong[transectlong$Species=="guttatus",],
             aes(Position,Prop),size=4,shape=21,fill="#ca0020",colour="white")+
  xlab(NULL)+guides(size = "none")+
  theme_classic()+ylab(NULL)

p2<-ggplot(transectlong, aes(Position,Freq),fill=Species)+
  geom_area(aes(fill=Species),position="stack",alpha=0.7)+
  scale_fill_manual(values = c("#ca0020","#0571b0"))+
  geom_point(data=transectlong[transectlong$Species=="nudatus",],
             aes(Position,Prop),size=4,shape=24,fill="#0571b0",colour="white")+
  xlab("Position Along Transect (cm)")+guides(size = "none")+
  theme_classic()+ylab(NULL)

result<-p1/p2+ plot_layout(guides = 'collect')
gt <- patchwork::patchworkGrob(result)
gridExtra::grid.arrange(gt, left = "Seed Viability (points)\\nFruit Proportion (background)")
#saved as 5x5pdf; Figure2-transect.pdf
#####

#2018 transplant experiment
#####
#2018 data import fruit separated
freqhyb2 <- read.csv('FreqDepHyb2018.csv')
freqhyb2$Total <- freqhyb2$Normal+freqhyb2$Flat
freqhyb2$PropViable <- freqhyb2$Normal/freqhyb2$Total

#summing seeds counted per fruit
freqhyb <- aggregate (cbind(Normal=freqhyb2$Normal, Flat=freqhyb2$Flat),
                      by= list (Block=freqhyb2$Block, Family=freqhyb2$Family, Frequency=freqhyb2$Frequency,
                                Flower=freqhyb2$Flower, Fruit=freqhyb2$Fruit,
                                Neighbors=freqhyb2$Neighbors),
                      FUN=sum, na.rm= TRUE )

freqhyb$Total <- freqhyb$Normal+freqhyb$Flat
freqhyb$HybRate <- freqhyb$Flat/freqhyb$Total

#Analyze flowers per immigrant focal M. guttatus plant in 2018

#Flowers
#fit poisson model
Fpoi18<-glmmTMB(Flower~Neighbors+(1|Family),freqhyb,family=poisson)
overdisp_fun(Fpoi18) #p=3.252629e-50 
#DHARMa model diagnostics:
simulationOutput <- simulateResiduals(fittedModel = Fpoi18, plot = T)
plotQQunif(simulationOutput) #KS, dispersion, and outlier tests sig
testZeroInflation(simulationOutput) #zero inflated

#fit negative binomial model
Fnb218<-glmmTMB(Flower~Neighbors+(1|Family),freqhyb,family=nbinom2)
#DHARMa model diagnostics:
simulationOutput <- simulateResiduals(fittedModel = Fnb218, plot = T)
plotQQunif(simulationOutput) #all tests NS
testZeroInflation(simulationOutput) #not zero inflated

#Wald type II chi sq:
car::Anova(Fnb218) #Chisq=0.1146,df=1,p= 0.7349

# Seeds Per Flower
#fit poisson model
Tpoi18<-glmmTMB(Total~Neighbors+(1|Block/Family),freqhyb2,family=poisson)
overdisp_fun(Tpoi18) #p=0.0000
#DHARMa model diagnostics:
simulationOutput <- simulateResiduals(fittedModel = Tpoi18, plot = T)
plotQQunif(simulationOutput) #KS and dispersion tests sig
testZeroInflation(simulationOutput) # not zero inflated

#fit negative binomial model
Tnb218<-glmmTMB(Total~Neighbors+(1|Block/Family),freqhyb2,family=nbinom2) #model convergence problem
#DHARMa model diagnostics:
simulationOutput <- simulateResiduals(fittedModel = Tnb218, plot = T)#quantile deviations detected
plotQQunif(simulationOutput) #KS test p=0
testZeroInflation(simulationOutput) # zero inflated

#fit zero-inflated negative binomial model
ziTnb218<-glmmTMB(Total~Neighbors+(1|Block/Family),ziformula = ~1,freqhyb2,family=nbinom2)
#DHARMa model diagnostics:
simulationOutput <- simulateResiduals(fittedModel = ziTnb218, plot = T)#quantile deviations detected
plotQQunif(simulationOutput) #KS test p=0.02
testZeroInflation(simulationOutput) # not zero inflated

#Wald type II chi sq:
car::Anova(ziTnb218) #Chisq=0.5168,df=1,p= 0.4722

#Proportion Viable Seeds (=1-hybridization rate), fit binomial model
PropVia18<-glmmTMB(cbind(Normal,Flat)~Neighbors+(1|Family),freqhyb,family=binomial)
simulationOutput <- simulateResiduals(fittedModel = PropVia18, plot = T) #quantile deviations detected
plotQQunif(simulationOutput) #all tests sig

#fit betabinomal model
PropVia18<-glmmTMB(cbind(Normal,Flat)~Neighbors+(1|Family),freqhyb,family=betabinomial(link="logit"))
simulationOutput <- simulateResiduals(fittedModel = PropVia18, plot = T)
plotQQunif(simulationOutput) #dispersion test significant, but model diagnostics better

#Wald type II chi sq:
car::Anova(PropVia18) #Chisq=7.4921,df=1,p= 0.006197

summary(PropVia18) #est= 0.02, SE= 0.009, z-value= 2.74, p=0.006

#Plot fecundity data, hybridization & hybridization predictions

PropViadf<-ggpredict(PropVia18,terms="Neighbors")

Flow18<-ggplot(data=freqhyb,aes(x=Neighbors,y=Flower)) + ylim(0,40) +
  geom_point(position = "jitter",size=3,shape=21,fill="#ca0020",colour="white",alpha=0.7) + 
  xlab(NULL)+ylab("Flowers Per Plant")+theme_classic()+
  annotate("text",y=40,x=100,label="NS",vjust=1,hjust=1)

ToFl18<-ggplot(data=freqhyb2,aes(x=Neighbors,y=Total)) + ylim(0,400) +
  geom_point(position = "jitter",size=3,shape=21,fill="#ca0020",colour="white",alpha=0.7) + 
  xlab(NULL)+ylab("Seeds Per Flower")+theme_classic()+
  annotate("text",y=400,x=100,label="NS",vjust=1,hjust=1)

PropVia18<-ggplot(PropViadf, aes(x, predicted)) + ylim(0,1) +
  geom_point(data=freqhyb,aes(x=Neighbors,y=(Normal/Total)),position = "jitter",
             size=3,shape=21,fill="#ca0020",colour="white",alpha=0.7) + 
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  xlab("Conspecific Neighbors")+ylab("Seed Viability")+theme_classic()


#Figure 3 - saved as 3x5pdf#
Flow18/ToFl18/PropVia18+ plot_annotation(tag_levels = c('A'))
#####

#Supplemental sensitivity analysis of M. guttatus immigration experiment
#Analysis with log transformed neighbors
#Flowers
#fit poisson model
Fpoi18<-glmmTMB(Flower~log(Neighbors+1)+(1|Family),freqhyb,family=poisson)
overdisp_fun(Fpoi18) #p=1.945682e-51 
#DHARMa model diagnostics:
simulationOutput <- simulateResiduals(fittedModel = Fpoi18, plot = T)
plotQQunif(simulationOutput) #KS, disperpsion, and outlier tests sig
testZeroInflation(simulationOutput) #zero inflated

#fit negative binomial model
Fnb218<-glmmTMB(Flower~log(Neighbors+1)+(1|Family),freqhyb,family=nbinom2)
#DHARMa model diagnostics:
simulationOutput <- simulateResiduals(fittedModel = Fnb218, plot = T)
plotQQunif(simulationOutput) #all tests NS
testZeroInflation(simulationOutput) #not zero inflated

#Wald type II chi sq:
car::Anova(Fnb218) #Chisq=0.0196,df=1,p= 0.8886

# Seeds Per Flower
#fit poisson model
Tpoi18<-glmmTMB(Total~log(Neighbors+1)+(1|Block/Family),freqhyb2,family=poisson)
overdisp_fun(Tpoi18) #p=0.0000
#DHARMa model diagnostics:
simulationOutput <- simulateResiduals(fittedModel = Tpoi18, plot = T)
plotQQunif(simulationOutput) #KS test sig
testZeroInflation(simulationOutput) # not zero inflated

#fit negative binomial model
Tnb218<-glmmTMB(Total~log(Neighbors+1)+(1|Block/Family),freqhyb2,family=nbinom2) 
#DHARMa model diagnostics:
simulationOutput <- simulateResiduals(fittedModel = Tnb218, plot = T)#quantile deviations detected
plotQQunif(simulationOutput) #KS test=0
testZeroInflation(simulationOutput) # zero inflated

#fit negative binomial model
ziTnb218<-glmmTMB(Total~log(Neighbors+1)+(1|Block/Family),ziformula = ~1,freqhyb2,family=nbinom2)
#DHARMa model diagnostics:
simulationOutput <- simulateResiduals(fittedModel = ziTnb218, plot = T)#quantile deviations detected
plotQQunif(simulationOutput) #KS test p=0.0225
testZeroInflation(simulationOutput) # not zero inflated

#Wald type II chi sq:
car::Anova(ziTnb218) #Chisq=0.3444,df=1,p= 0.5573

#Proportion Viable Seeds (=1-hybridization rate)
PropVia18<-glmmTMB(cbind(Normal,Flat)~log(Neighbors+1)+(1|Family),freqhyb,family=binomial)
simulationOutput <- simulateResiduals(fittedModel = PropVia18, plot = T) #quantile deviations detected
plotQQunif(simulationOutput) #all tests sig

logNeighbors<-log(freqhyb$Neighbors+1)
PropVia18<-glmmTMB(cbind(Normal,Flat)~logNeighbors+(1|Family),freqhyb,family=betabinomial(link="logit"))
simulationOutput <- simulateResiduals(fittedModel = PropVia18, plot = T)
plotQQunif(simulationOutput) #dispersion test significant, but model diagnostics better

#Wald type II chi sq
car::Anova(PropVia18) #Chisq=7.9541,df=1,p= 0.004798

summary(PropVia18) #est= 0.28, SE= 0.099, z-value= 2.82, p=0.005

#Plot fecundity data, hybridization & hybridization predictions

PropViadf<-ggpredict(PropVia18,terms="logNeighbors")

Flow18<-ggplot(data=freqhyb,aes(x=logNeighbors,y=Flower)) + ylim(0,40) +
  geom_point(position = "jitter",size=3,shape=21,fill="#ca0020",colour="white",alpha=0.7) + 
  xlab(NULL)+ylab("Flowers Per Plant")+theme_classic()+
  annotate("text",y=40,x=5,label="NS",vjust=1,hjust=1)

ToFl18<-ggplot(data=freqhyb2,aes(x=log(Neighbors+1),y=Total)) + ylim(0,400) +
  geom_point(position = "jitter",size=3,shape=21,fill="#ca0020",colour="white",alpha=0.7) + 
  xlab(NULL)+ylab("Seeds Per Flower")+theme_classic()+
  annotate("text",y=400,x=5,label="NS",vjust=1,hjust=1)

PropVia18<-ggplot(PropViadf, aes(x, predicted)) + ylim(0,1) +
  geom_point(data=freqhyb,aes(x=logNeighbors,y=(Normal/Total)),position = "jitter",
             size=3,shape=21,fill="#ca0020",colour="white",alpha=0.7) + 
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  xlab("log(Conspecific Neighbors + 1)")+ylab("Seed Viability")+theme_classic()


#Figure S1 - saved as 3x5pdf#
Flow18/ToFl18/PropVia18+ plot_annotation(tag_levels = c('A'))

#Analysis without  outlier focal plant
no_outlier_freqhyb2<-subset(freqhyb2,freqhyb2$Neighbors<103)
#summing seeds counted per fruit
no_outlier_freqhyb <- aggregate (cbind(Normal=no_outlier_freqhyb2$Normal, Flat=no_outlier_freqhyb2$Flat),
                                 by= list (Block=no_outlier_freqhyb2$Block, Family=no_outlier_freqhyb2$Family, Frequency=no_outlier_freqhyb2$Frequency,
                                           Flower=no_outlier_freqhyb2$Flower, Fruit=no_outlier_freqhyb2$Fruit,
                                           Neighbors=no_outlier_freqhyb2$Neighbors),
                                 FUN=sum, na.rm= TRUE )

no_outlier_freqhyb$Total <- no_outlier_freqhyb$Normal+no_outlier_freqhyb$Flat
no_outlier_freqhyb$HybRate <- no_outlier_freqhyb$Flat/no_outlier_freqhyb$Total

#Flowers
#fit poisson model
Fpoi18<-glmmTMB(Flower~Neighbors+(1|Family),no_outlier_freqhyb,family=poisson)
overdisp_fun(Fpoi18) #p=4.063470e-49  
#DHARMa model diagnostics:
simulationOutput <- simulateResiduals(fittedModel = Fpoi18, plot = T)
plotQQunif(simulationOutput) #KS, and dispersion  tests sig
testZeroInflation(simulationOutput) #zero inflated

#fit negative binomial model
Fnb218<-glmmTMB(Flower~Neighbors+(1|Family),no_outlier_freqhyb,family=nbinom2)
#DHARMa model diagnostics:
simulationOutput <- simulateResiduals(fittedModel = Fnb218, plot = T)
plotQQunif(simulationOutput) #all tests NS
testZeroInflation(simulationOutput) #not zero inflated

#Wald type II chi sq:
car::Anova(Fnb218) #Chisq=0.9273,df=1,p= 0.3356

# Seeds Per Flower
#fit poisson model
Tpoi18<-glmmTMB(Total~Neighbors+(1|Block/Family),no_outlier_freqhyb2,family=poisson)
overdisp_fun(Tpoi18) #p=0.0000
#DHARMa model diagnostics:
simulationOutput <- simulateResiduals(fittedModel = Tpoi18, plot = T)
plotQQunif(simulationOutput) #KS and dispersion tests sig
testZeroInflation(simulationOutput) # not zero inflated

#fit negative binomial model
Tnb218<-glmmTMB(Total~Neighbors+(1|Block/Family),no_outlier_freqhyb2,family=nbinom2) #model convergence problem
#DHARMa model diagnostics:
simulationOutput <- simulateResiduals(fittedModel = Tnb218, plot = T)#quantile deviations detected
plotQQunif(simulationOutput) #KS test=0
testZeroInflation(simulationOutput) # zero inflated

#fit negative binomial model
ziTnb218<-glmmTMB(Total~Neighbors+(1|Block/Family),ziformula = ~1,no_outlier_freqhyb2,family=nbinom2)
#DHARMa model diagnostics:
simulationOutput <- simulateResiduals(fittedModel = ziTnb218, plot = T)#quantile deviations detected
plotQQunif(simulationOutput) #KS test p=0.03
testZeroInflation(simulationOutput) # not zero inflated

#Wald type II chi sq:
car::Anova(ziTnb218) #Chisq=0.3183,df=1,p= 0.5727

#Proportion Viable Seeds (=1-hybridization rate)
PropVia18<-glmmTMB(cbind(Normal,Flat)~Neighbors+(1|Family),no_outlier_freqhyb,family=binomial)
simulationOutput <- simulateResiduals(fittedModel = PropVia18, plot = T) #quantile deviations detected
plotQQunif(simulationOutput) #all tests sig

PropVia18<-glmmTMB(cbind(Normal,Flat)~Neighbors+(1|Family),no_outlier_freqhyb,family=betabinomial(link="logit"))
simulationOutput <- simulateResiduals(fittedModel = PropVia18, plot = T)
plotQQunif(simulationOutput) #dispersion test significant, but model diagnostics better

#Wald type II chi sq:
car::Anova(PropVia18) #Chisq=4.2758,df=1,p= 0.03866

summary(PropVia18) #est= 0.026, SE= 0.013, z-value= 2.07, p=0.039

#Plot fecundity data, hybridization & hybridization predictions

PropViadf<-ggpredict(PropVia18,terms="Neighbors")

Flow18<-ggplot(data=no_outlier_freqhyb,aes(x=Neighbors,y=Flower)) + ylim(0,40) +
  geom_point(position = "jitter",size=3,shape=21,fill="#ca0020",colour="white",alpha=0.7) + 
  xlab(NULL)+ylab("Flowers Per Plant")+theme_classic()+
  annotate("text",y=40,x=40,label="NS",vjust=1,hjust=1)

ToFl18<-ggplot(data=no_outlier_freqhyb2,aes(x=Neighbors,y=Total)) + ylim(0,400) +
  geom_point(position = "jitter",size=3,shape=21,fill="#ca0020",colour="white",alpha=0.7) + 
  xlab(NULL)+ylab("Seeds Per Flower")+theme_classic()+
  annotate("text",y=400,x=40,label="NS",vjust=1,hjust=1)

PropVia18<-ggplot(PropViadf, aes(x, predicted)) + ylim(0,1) +
  geom_point(data=no_outlier_freqhyb,aes(x=Neighbors,y=(Normal/Total)),position = "jitter",
             size=3,shape=21,fill="#ca0020",colour="white",alpha=0.7) + 
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
  xlab("Conspecific Neighbors (no outlier)")+ylab("Seed Viability")+theme_classic()


#Figure S2 - saved as 3x5pdf#
Flow18/ToFl18/PropVia18+ plot_annotation(tag_levels = c('A'))


#2019 frequency experiment

#load data!
data<-read.csv('FreqDepHyb2019.csv')
str(data)
data$Species <- as.factor(data$Species)
data$Block <- as.factor(data$Block)
data$Plant <- as.factor(data$Plant)

data$ImmigrantFreq <- data$ImmigrantFreq/100
data$HybRate <- data$Flat/(data$Normal+data$Flat)
data$Total <- data$Normal+data$Flat
data$fImmigrantFreq <- as.factor(data$ImmigrantFreq)

#Is whole plant total seed production freq dep?#
#sum plant flower counts#
data$flowercount<-1
agg_sum <- aggregate (cbind(Flat=data$Flat,Normal=data$Normal,
                            flowercount=data$flowercount),
                      by= list (Species=data$Species, Site=data$Site,
                                ImmigrantFreq=data$ImmigrantFreq,
                                Plant=data$Plant, Block=data$Block),
                      FUN=sum, na.rm= TRUE )

agg_sum$HybRate <- agg_sum$Flat/(agg_sum$Normal+agg_sum$Flat)
agg_sum$Total <- agg_sum$Normal+agg_sum$Flat
agg_sum$fImmigrantFreq <- as.factor(agg_sum$ImmigrantFreq)

#to test for frequency dependent home species advantage!
sumCS<-agg_sum[agg_sum$Site == 'CS',]
sumQV<-agg_sum[agg_sum$Site == 'QV',]
#subset blocks where both species are present
sumCS0<-sumCS[sumCS$ImmigrantFreq > 0,]
sumQV0<-sumQV[sumQV$ImmigrantFreq > 0,]

#to test for frequency dependent fecundity within species, separated by site
sumGuttatusCS<-agg_sum[agg_sum$Species == 'guttatus' & agg_sum$Site == 'CS',]
sumNudatusCS<-agg_sum[agg_sum$Species == 'nudatus' & agg_sum$Site == 'CS',]
sumGuttatusQV<-agg_sum[agg_sum$Species == 'guttatus' & agg_sum$Site == 'QV',]
sumNudatusQV<-agg_sum[agg_sum$Species == 'nudatus' & agg_sum$Site == 'QV',]

#*SUBSET FOR TOTAL SEED PER FLOWER ANALYSIS*#
GuttatusCS<-data[data$Species == 'guttatus' & data$Site == 'CS',]
NudatusCS<-data[data$Species == 'nudatus' & data$Site == 'CS',]
GuttatusQV<-data[data$Species == 'guttatus' & data$Site == 'QV',]
NudatusQV<-data[data$Species == 'nudatus' & data$Site == 'QV',]

#Lifetime fecundity analysis
#SEEP SITE (M. guttatus native)#
#####  
#Analyze viable seeds per plant in the seep - resident M. guttatus & immigrant M. nudatus

#Fit poisson model:
CSVipoi19<-glmmTMB(Normal~fImmigrantFreq*Species+(1|Block),sumCS0,family=poisson)
#Test for overdispersion:
overdisp_fun(CSVipoi19) #p=0.00; poisson model is significantly overdispersed
simulationOutput <- simulateResiduals(fittedModel = CSVipoi19, plot = T)
plotQQunif(simulationOutput) #KS and outlier tests sig

#Fit negative binomial model:
CSVinb219<-glmmTMB(Normal~fImmigrantFreq*Species+(1|Block),sumCS0,family=nbinom2)
simulationOutput <- simulateResiduals(fittedModel = CSVinb219, plot = T)
plotQQunif(simulationOutput) #KS test sig
testZeroInflation(simulationOutput) #zero inflated

#Fit zero-inflated negative binomial model:
ziCSVinb219<-glmmTMB(Normal~fImmigrantFreq*Species+(1|Block),ziformula = ~1,sumCS0,family=nbinom2)
simulationOutput <- simulateResiduals(fittedModel = ziCSVinb219, plot = T)
plotQQunif(simulationOutput) #all tests NS

#Wald type III Wald-Chi square tests :
car::Anova(ziCSVinb219, type=3) #frequency, species, and interaction significant
#Tukey posthoc contrasts:
pairs(emmeans(ziCSVinb219, ~c(Species,fImmigrantFreq)))

#cld tukey groupings:
cld(emmeans(ziCSVinb219, ~c(fImmigrantFreq,Species)))
#g0.05b, n0.05 a; 0.25g b, 0.25n b; 0.5g b, 0.5n a

#Estimate marginal means and bootstrap confidence intervals:
ziCSVi19df<- ggpredict(ziCSVinb219,terms=c("fImmigrantFreq","Species"),nsim=500)
plot(ziCSVi19df)

CSVinewx2 <- as.data.frame(seq(0,5, by = 1))
CSVinewx2$fImmigrantFreq<-as.factor(c(0.05,0.05,0.25,0.25,0.5,0.5))
CSVinewx2$Site= "CS"
CSVinewx2$Species= c("guttatus","nudatus","guttatus","nudatus","guttatus","nudatus")
CSVinewx2$Normal<-ziCSVi19df$predicted
CSVinewx2$se<-ziCSVi19df$std.error
CSVinewx2$conf.low <- ziCSVi19df$conf.low
CSVinewx2$conf.high <- ziCSVi19df$conf.high

CSVinewx2$Species <- relevel(as.factor(CSVinewx2$Species), ref = "nudatus")
sumCS0$Species <- relevel(as.factor(sumCS0$Species), ref = "nudatus")
anno_seep<-data.frame(xstar = c(1,2,1,2,1,2), ystar = c(47,299,341,242,42,190),
                      lab = c("A","B","B","B","A","B"),
                      fImmigrantFreq = c(0.05,0.05,0.25,0.25,0.5,0.5),
                      Species=c("nudatus","guttatus","nudatus","guttatus","nudatus","guttatus"))
CSVilabels=c("0.05" = "1N:18G\\n5%", 
             "0.25" = "4N:15G\\n25%","0.5" = "9N:10G\\n50%")

CSViPlot<-ggplot(CSVinewx2, aes(x=Species, y=Normal,fill=Species)) +   theme_classic()+
  facet_wrap(.~fImmigrantFreq, strip.position = "bottom", scales = "free_x",labeller = as_labeller(CSVilabels)) +  
  theme(panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),
        strip.placement = "outside",
        axis.text.x = element_text(face = "italic"),
        legend.position = "none") +
  geom_point(data=sumCS0,aes(x=Species,y=Normal,shape=Species),position = position_jitterdodge(.9),alpha=0.5,colour="white",size=3)+
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high), position = position_dodge2(width=0.5),size=1,width=0.5)+
  geom_point(aes(shape=Species,fill=Species),position = position_dodge2(width = 0.9),color="black",size=4,stroke=1) +
  scale_fill_manual(values = c("#0571b0","#ca0020"))+scale_colour_manual(values = c("#0571b0","#ca0020"))+
  scale_shape_manual(values=c(24,21))+ylab(NULL)+
  xlab("Seep Immigrant Frequency")+scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),breaks=c(0,4,40,400,4000),limits=c(0,4000))+
  geom_label(data = anno_seep, aes(x = xstar,  y = ystar, label = lab,fontface=2),nudge_y = 0.25,label.size = NA,fill="white",alpha=0.75) 
#####

#Lifetime fecundity analysis
#WASH SITE (M. nudatus native)#
#####
#Analyze viable seeds per plant in the wash - resident M. nudatus & immigrant M. guttatus

#Fit poisson model:
QVVipoi19<-glmmTMB(Normal~fImmigrantFreq*Species+(1|Block),sumQV0,family=poisson)
#Test for overdispersion:
overdisp_fun(QVVipoi19) #p=0.00; poisson model is significantly overdispersed
simulationOutput <- simulateResiduals(fittedModel = QVVipoi19, plot = T)
plotQQunif(simulationOutput) #all tests sig
testZeroInflation(simulationOutput) #Very zero inflated

#Fit negative binomial model:
QVVinb219<-glmmTMB(Normal~fImmigrantFreq*Species+(1|Block),sumQV0,family=nbinom2)
simulationOutput <- simulateResiduals(fittedModel = QVVinb219, plot = T)
plotQQunif(simulationOutput) #KS test sig
testZeroInflation(simulationOutput) #still zero inflated

ziQVVinb219<-glmmTMB(Normal~fImmigrantFreq*Species+(1|Block),ziformula = ~1,sumQV0,family=nbinom2)
simulationOutput <- simulateResiduals(fittedModel = ziQVVinb219, plot = T)
plotQQunif(simulationOutput) #all tests NS
testZeroInflation(simulationOutput) #no longer zero inflated

#Wald typeIII Wald-Chi square tests :
car::Anova(ziQVVinb219, type=3) # all terms significant

#Tukey post-hoc contrasts:
test(pairs(emmeans(ziQVVinb219, ~c(Species,fImmigrantFreq)))) 

cld(emmeans(ziQVVinb219, ~c(fImmigrantFreq,Species)))
#g0.05 a, n0.05 b; 0.25g ab, 0.25n ab; 0.5g ab, 0.5n ab

ziQVVi19df<-ggpredict(ziQVVinb219,terms=c("fImmigrantFreq","Species"),nsim=500)
plot(ziQVVi19df)

QVVinewx2 <- as.data.frame(seq(0,5, by = 1))
QVVinewx2$fImmigrantFreq<-as.factor(c(0.05,0.05,0.25,0.25,0.5,0.5))
QVVinewx2$Site= "QV"
QVVinewx2$Species= c("guttatus","nudatus","guttatus","nudatus","guttatus","nudatus")
QVVinewx2$Normal<-ziQVVi19df$predicted
QVVinewx2$se<-ziQVVi19df$std.error
QVVinewx2$conf.low <- ziQVVi19df$conf.low
QVVinewx2$conf.high <- ziQVVi19df$conf.high

anno_wash<-data.frame(xstar = c(1,2,1,2,1,2), ystar = c(172,569,602,298,483,386),
                      lab = c("A","B","AB","AB","AB","AB"),
                      fImmigrantFreq = c(0.05,0.05,0.25,0.25,0.5,0.5),
                      Species=c("guttatus","nudatus","guttatus","nudatus","guttatus","nudatus"))

QVVilabels=c("0.05" = "1G:18N\\n5%","0.25" = "4G:15N\\n25%","0.5" = "9G:10N\\n50%")

QVViPlot<-ggplot(QVVinewx2, aes(x=Species, y=Normal,fill=Species)) + theme_classic()+
  facet_wrap(.~fImmigrantFreq, strip.position = "bottom", scales = "free_x",labeller = as_labeller(QVVilabels)) +
  theme(panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),
        strip.placement = "outside",
        axis.text.x = element_text(face = "italic"),
        legend.position = "none") +
  geom_point(data=sumQV0,aes(x=Species,y=Normal,shape=Species),position = position_jitterdodge(.9),alpha=0.5,colour="white",size=3)+
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high), position = position_dodge2(width=0.5),size=1,width=0.5)+
  geom_point(aes(shape=Species,fill=Species),position = position_dodge2(width = 0.9),color="black",size=4,stroke=1) +
  scale_fill_manual(values = c("#ca0020","#0571b0"))+scale_colour_manual(values = c("#ca0020","#0571b0"))+
  scale_shape_manual(values=c(21,24))+ 
  ylab("Viable Seeds Per Plant")+theme(legend.position="none")+
  xlab("Wash Immigrant Frequency")+scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),breaks=c(0,4,40,400,4000),limits=c(0,4000))+
  geom_label(data = anno_wash, aes(x = xstar,  y = ystar, label = lab,fontface=2),nudge_y = 0.25,label.size = NA,fill="white",alpha=0.75)

#Plotting marginal means from models (fecundity=species*frequency) together for a panel
QVViPlot+CSViPlot+ plot_layout(guides = 'collect')
#saved as 8x4 PDF - Figure6-LifetimeFecundity.pdf
#####

#Fecundity component analysis#

#2019#
#Flowers per plant
#SEEP SITE (M. guttatus native)#
#guttatus, seep
#Fit poisson model:
GUCSFlpoi19<-glmmTMB(flowercount~fImmigrantFreq+(1|Block),sumGuttatusCS,family=poisson)
#Test for overdispersion:
overdisp_fun(GUCSFlpoi19) #p=8.382036e-76; poisson model is significantly overdispersed
#DHARMa model diagnostics:
simulationOutput <- simulateResiduals(fittedModel = GUCSFlpoi19, plot = T)
plotQQunif(simulationOutput) #all tests p<0.05
testZeroInflation(simulationOutput) #not zero inflated

#Fit negative binomial model:
GUCSFlnb219<-glmmTMB(flowercount~fImmigrantFreq+(1|Block),sumGuttatusCS,family=nbinom2)
#DHARMa model diagnostics:
simulationOutput <- simulateResiduals(fittedModel = GUCSFlnb219, plot = T)
plotQQunif(simulationOutput) #KS test p<0.05
testZeroInflation(simulationOutput) #not zero inflated

#Wald typeII Wald-Chi square tests:
car::Anova(GUCSFlnb219) #chisq=3.1432,df=3,p=0.3701

#Estimate marginal means and bootstrap confidence intervals:
GuCSFl19df<-ggpredict(GUCSFlnb219,terms="fImmigrantFreq",nsim=500)
plot(GuCSFl19df)
pairs((emmeans(GUCSFlnb219, ~fImmigrantFreq))) #all NS
cld((emmeans(GUCSFlnb219, ~fImmigrantFreq))) #all A

#nudatus, seep
#Fit poisson model:
NuCSFlpoi19<-glmmTMB(flowercount~fImmigrantFreq+(1|Block),sumNudatusCS,family=poisson)
#Test for overdispersion:
overdisp_fun(NuCSFlpoi19) #p=5.365664e-04; poisson model is significantly overdispersed
#DHARMa model diagnostics:
simulationOutput <- simulateResiduals(fittedModel = NuCSFlpoi19, plot = T)
plotQQunif(simulationOutput) #dispersion test p=0.088; marginal but since overdisp_fun says sig will fit nb
testZeroInflation(simulationOutput) #not zero inflated

#Fit negative binomial model:
NuCSFlnb219<-glmmTMB(flowercount~fImmigrantFreq+(1|Block),sumNudatusCS,family=nbinom2)
overdisp_fun(NuCSFlnb219) #p=0.8566308
#DHARMa model diagnostics:
simulationOutput <- simulateResiduals(fittedModel = NuCSFlpoi19, plot = T)
plotQQunif(simulationOutput) #dispersion test still p=0.088
testZeroInflation(simulationOutput) #not zero inflated
#Wald typeII Wald-Chi square tests:
car::Anova(NuCSFlnb219) #chisq=11.738,df=2,p=0.002825
#Tukey Post-hoc:
pairs((emmeans(NuCSFlnb219, ~fImmigrantFreq)))

#To get compact letter display for figures:
cld(emmeans(NuCSFlnb219, ~fImmigrantFreq)) #0.05 ab, 0.25a, 0.5b
NuCSFl19df<-ggpredict(NuCSFlnb219,terms="fImmigrantFreq",nsim=500)
plot(NuCSFl19df)

#plot predictions from models (fecundity=immigrant frequency + (1|block))

#FLOWER PREDS#
GuCSFlnewx <- as.data.frame(seq(0,3, by = 1))
GuCSFlnewx$fImmigrantFreq<-GuCSFl19df$x
GuCSFlnewx$Site= "CS"
GuCSFlnewx$Species= '"Resident"~italic("M. guttatus")'
GuCSFlnewx$flowercount<-GuCSFl19df$predicted
GuCSFlnewx$se<-GuCSFl19df$std.error
GuCSFlnewx$conf.low <- GuCSFl19df$conf.low
GuCSFlnewx$conf.high <- GuCSFl19df$conf.high
GuCSFlnewx<-subset(GuCSFlnewx, select = -c(`seq(0, 3, by = 1)`))

anno_guflseep<-data.frame(xstar = c(1,2,3,4), ystar = c(6.951914,5.560195,4.803992,5.335045),
                          lab = c("A","A","A","A"),
                          fImmigrantFreq = c(0,0.05,0.25,0.5))

GUCSFlPlot<-ggplot(GuCSFlnewx, aes(x=fImmigrantFreq, y=flowercount)) + theme_classic()+
  theme(legend.position = "none") +
  geom_jitter(data=sumGuttatusCS,aes(x=fImmigrantFreq,y=flowercount),alpha=0.5,fill="#ca0020",colour="white",size=3,shape=21)+
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high), position = position_dodge2(width=0.5),size=1,width=0.5)+
  geom_point(position = position_dodge2(width = 0.9),color="black",fill="#ca0020",size=4,stroke=1,shape=21) +
  ylab(NULL)+theme(legend.position="none")+
  xlab(NULL)+annotate("text",y=Inf,x=2.5,label="Resident"~italic("M. guttatus"),vjust=1,hjust=0.5)+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),limits=c(0,70))+
  scale_x_discrete(labels=c("0" = "0N:19G\\n0%", "0.05" = "1N:18G\\n5%","0.25" = "4N:15G\\n25%","0.5" = "9N:10G\\n50%"))+
  geom_label(data = anno_guflseep, aes(x = xstar,  y = ystar, label = lab,fontface=2,alpha=0.75),fill="white",nudge_y = 0.25,label.size = NA)


NuCSFlnewx <- as.data.frame(seq(0,2, by = 1))
NuCSFlnewx$fImmigrantFreq<-NuCSFl19df$x
NuCSFlnewx$Site= "CS"
NuCSFlnewx$Species= '"Immigrant"~italic("M. nudatus")'
NuCSFlnewx$flowercount<-NuCSFl19df$predicted
NuCSFlnewx$se<-NuCSFl19df$std.error
NuCSFlnewx$conf.low <- NuCSFl19df$conf.low
NuCSFlnewx$conf.high <- NuCSFl19df$conf.high
NuCSFlnewx<-subset(NuCSFlnewx, select = -c(`seq(0, 2, by = 1)`))

anno_nuflseep<-data.frame(xstar = c(1,2,3), ystar = c(9.607141,9.381954,2.760228),
                          lab = c("AB","A","B"),
                          fImmigrantFreq = c(0.05,0.25,0.5))

NUCSFlPlot<-ggplot(NuCSFlnewx, aes(x=fImmigrantFreq, y=flowercount)) + theme_classic()+
  theme(legend.position = "none") +
  geom_jitter(data=sumNudatusCS,aes(x=fImmigrantFreq,y=flowercount),alpha=0.5,fill="#0571b0",colour="white",size=3,shape=24)+
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high), position = position_dodge2(width=0.5),size=1,width=0.5)+
  geom_point(position = position_dodge2(width = 0.9),color="black",fill="#0571b0",size=4,stroke=1,shape=24) +
  ylab("Flowers Per Plant")+theme(legend.position="none")+
  xlab(NULL)+annotate("text",y=Inf,x=2,label="Immigrant"~italic("M. nudatus"),vjust=1,hjust=0.5)+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),limits=c(0,70))+
  scale_x_discrete(labels=c("0.05" = "1N:18G\\n5%","0.25" = "4N:15G\\n25%","0.5" = "9N:10G\\n50%"))+
  geom_label(data = anno_nuflseep, aes(x = xstar,  y = ystar, label = lab,fontface=2,alpha=0.75),fill="white",nudge_y = 0.25,label.size = NA)


#####
#nudatus, wash
#Analyze flowers per plant in the wash - resident M. nudatus only
#Fit poisson model:
NuQVFlpoi19<-glmmTMB(flowercount~fImmigrantFreq+(1|Block),sumNudatusQV,family=poisson)
#Test for overdispersion:
overdisp_fun(NuQVFlpoi19) #p=1.036546e-295; poisson model is significantly overdispersed
#DHARMa model diagnostics:
simulationOutput <- simulateResiduals(fittedModel = NuQVFlpoi19, plot = T)
plotQQunif(simulationOutput) #ks and dispersion tests sig
testZeroInflation(simulationOutput) #not zero inflated

#Fit negative binomial model:
NuQVFlnb219<-glmmTMB(flowercount~fImmigrantFreq+(1|Block),sumNudatusQV,family=nbinom2)
#Test for overdispersion:
overdisp_fun(NuQVFlnb219) #p=0.5755684
#DHARMa model diagnostics:
simulationOutput <- simulateResiduals(fittedModel = NuQVFlnb219, plot = T)
plotQQunif(simulationOutput) #all tests NS
testZeroInflation(simulationOutput) #not zero inflated

#Wald typeII Wald-Chi square tests for graphing:
car::Anova(NuQVFlnb219) #Chisq=3.1755,df=3,p=0.3653
#Tukey Post-hoc:
pairs(emmeans(NuQVFlnb219, ~fImmigrantFreq)) #all p>0.35
#Estimate marginal means and bootstrap confidence intervals:
NuQVFl19df<-ggpredict(NuQVFlnb219,terms="fImmigrantFreq",nsim=500)
plot(NuQVFl19df)

cld((emmeans(NuQVFlnb219, ~fImmigrantFreq))) #all A

#####
#guttatus, wash
#Fit poisson model:
GuQVFlpoi19<-glmmTMB(flowercount~fImmigrantFreq+(1|Block),sumGuttatusQV,family=poisson)
#Test for overdispersion:
overdisp_fun(GuQVFlpoi19) #p=3.305883e-46; poisson model is significantly overdispersed
#DHARMa model diagnostics:
simulationOutput <- simulateResiduals(fittedModel = GuQVFlpoi19, plot = T)
plotQQunif(simulationOutput) #KS test p=0.01
testZeroInflation(simulationOutput) #not zero inflated

#Fit negative binomial model:
GuQVFlnb219<-glmmTMB(flowercount~fImmigrantFreq+(1|Block),sumGuttatusQV,family=nbinom2)
#Test for overdispersion:
overdisp_fun(GuQVFlnb219) #p=0.5432238
#DHARMa model diagnostics:
simulationOutput <- simulateResiduals(fittedModel = GuQVFlnb219, plot = T)
plotQQunif(simulationOutput) #all tests NS
testZeroInflation(simulationOutput) #not zero inflated

#Wald typeII Wald-Chi square tests:
car::Anova(GuQVFlnb219) #Chisq=1.7457,df=2,p=0.4178
#Tukey Post-hoc:
test(pairs(emmeans(GuQVFlnb219, ~fImmigrantFreq))) #all p>0.42
cld((emmeans(GuQVFlnb219, ~fImmigrantFreq))) #all A

#Estimate marginal means and bootstrap confidence intervals:
GuQVFl19df<-ggpredict(GuQVFlnb219,terms="fImmigrantFreq",nsim=500)
plot(GuQVFl19df)

#plot predictions from models (fecundity=immigrant frequency + (1|block))

#FLOWER PREDS#
GuQVFlnewx <- as.data.frame(seq(0,2, by = 1))
GuQVFlnewx$fImmigrantFreq<-GuQVFl19df$x
GuQVFlnewx$Site="QV"
GuQVFlnewx$Species= '"Immigrant"~italic("M. guttatus")'
GuQVFlnewx$flowercount<-GuQVFl19df$predicted
GuQVFlnewx$conf.low <- GuQVFl19df$conf.low
GuQVFlnewx$conf.high <- GuQVFl19df$conf.high
GuQVFlnewx<-subset(GuQVFlnewx, select = -c(`seq(0, 2, by = 1)`))

anno_guflwash<-data.frame(xstar = c(1,2,3), ystar = c(21.64654,14.66915,17.00406),
                          lab = c("A","A","A"),
                          fImmigrantFreq = c(0.05,0.25,0.5))
GUQVFlPlot<-ggplot(GuQVFlnewx, aes(x=fImmigrantFreq, y=flowercount)) + theme_classic()+
  theme(legend.position = "none") +
  geom_jitter(data=sumGuttatusQV,aes(x=fImmigrantFreq,y=flowercount),alpha=0.5,fill="#ca0020",colour="white",size=3,shape=21)+
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high), position = position_dodge2(width=0.5),size=1,width=0.5)+
  geom_point(position = position_dodge2(width = 0.9),color="black",fill="#ca0020",size=4,stroke=1,shape=21) +
  ylab("Flowers Per Plant")+theme(legend.position="none")+
  xlab(NULL)+annotate("text",y=Inf,x=2,label="Immigrant"~italic("M. guttatus"),vjust=1,hjust=0.5)+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),limits=c(0,70))+
  scale_x_discrete(labels=c("0.05" = "1G:18N\\n5%","0.25" = "4G:15N\\n25%","0.5" = "9G:10N\\n50%"))+
  geom_label(data = anno_guflwash, aes(x = xstar,  y = ystar, label = lab,fontface=2,alpha=0.75),fill="white",nudge_y = 0.25,label.size = NA)

NuQVFlnewx <- as.data.frame(seq(0,3, by = 1))
NuQVFlnewx$fImmigrantFreq<-NuQVFl19df$x
NuQVFlnewx$Site="QV"
NuQVFlnewx$Species= '"Resident"~italic("M. nudatus")'
NuQVFlnewx$flowercount<-NuQVFl19df$predicted
NuQVFlnewx$conf.low <- NuQVFl19df$conf.low
NuQVFlnewx$conf.high <- NuQVFl19df$conf.high
NuQVFlnewx<-subset(NuQVFlnewx, select = -c(`seq(0, 3, by = 1)`))

anno_nuflwash<-data.frame(xstar = c(1,2,3,4), ystar = c(15.05084,12.86196,16.31740,10.28323),
                          lab = c("A","A","A","A"),
                          fImmigrantFreq = c(0,0.05,0.25,0.5))
NUQVFlPlot<-ggplot(NuQVFlnewx, aes(x=fImmigrantFreq, y=flowercount)) + theme_classic()+
  theme(legend.position = "none") +
  geom_jitter(data=sumNudatusQV,aes(x=fImmigrantFreq,y=flowercount),alpha=0.5,fill="#0571b0",colour="white",size=3,shape=24)+
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high), position = position_dodge2(width=0.5),size=1,width=0.5)+
  geom_point(position = position_dodge2(width = 0.9),color="black",fill="#0571b0",size=4,stroke=1,shape=24) +
  ylab(NULL)+theme(legend.position="none")+
  xlab(NULL)+annotate("text",y=Inf,x=2.5,label="Resident"~italic("M. nudatus"),vjust=1,hjust=0.5)+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),limits=c(0,70))+
  scale_x_discrete(labels=c("0" = "0G:19N\\n0%","0.05" = "1G:18N\\n5%","0.25" = "4G:15N\\n25%","0.5" = "9G:10N\\n50%"))+
  geom_label(data = anno_nuflwash, aes(x = xstar,  y = ystar, label = lab,fontface=2,alpha=0.75),fill="white",nudge_y = 0.25,label.size = NA)

#M. guttatus, seep
#Analyze total seeds per flower in the seep - resident M. guttatus only
#Fit poisson model:
GuCSToFlpoi19<-glmmTMB(Total~fImmigrantFreq+(1|Block/Plant),GuttatusCS,family=poisson)
#Test for overdispersion:
overdisp_fun(GuCSToFlpoi19) #p=0.0000; poisson model is significantly overdispersed
#DHARMa model diagnostics:
simulationOutput <- simulateResiduals(fittedModel = GuCSToFlpoi19, plot = T)
plotQQunif(simulationOutput) #KS test and dispersion test p=0
testZeroInflation(simulationOutput) #zero inflated

#Fit negative binomial model:
GuCSToFlnb219<-glmmTMB(Total~fImmigrantFreq+(1|Block/Plant),GuttatusCS,family=nbinom2)
#DHARMa model diagnostics:
simulationOutput <- simulateResiduals(fittedModel = GuCSToFlnb219, plot = T)
plotQQunif(simulationOutput) #KS test and dispersion p=0
testZeroInflation(simulationOutput) #zero inflated

#Fit zero-inflated negative binomial model:
ziGuCSToFlnb219<-glmmTMB(Total~fImmigrantFreq+(1|Block/Plant),ziformula = ~1,GuttatusCS,family=nbinom2)
#DHARMa model diagnostics:
simulationOutput <- simulateResiduals(fittedModel = ziGuCSToFlnb219, plot = T)
plotQQunif(simulationOutput) #not overdispersed
testZeroInflation(simulationOutput) #not zero inflated
#does model fit better w/ zi term?
anova(ziGuCSToFlnb219,GuCSToFlnb219) #yes p<0.001
#Wald typeII Wald-Chi square tests
car::Anova(ziGuCSToFlnb219)#Chisq=5.7957,df=3,p=0.122

#Estimate marginal means and bootstrap confidence intervals:
ziGuCSToFl19df<-ggpredict(ziGuCSToFlnb219,terms="fImmigrantFreq",nsim=500)
plot(ziGuCSToFl19df)
pairs((emmeans(ziGuCSToFlnb219, ~fImmigrantFreq))) #all NS
cld((emmeans(ziGuCSToFlnb219, ~fImmigrantFreq))) #all A

#Nudatus, seep
#Analyze total seeds per flower in the seep - immigrant M. nudatus only
#Fit poisson model:
NuCSToFlpoi19<-glmmTMB(Total~fImmigrantFreq+(1|Block/Plant),NudatusCS,family=poisson)
#Test for overdispersion:
overdisp_fun(NuCSToFlpoi19) #p=0.00000; poisson model is significantly overdispersed
#DHARMa model diagnostics:
simulationOutput <- simulateResiduals(fittedModel = NuCSToFlpoi19, plot = T)
plotQQunif(simulationOutput) #KS test p=0.02, disp test=0.056
testZeroInflation(simulationOutput) #not zero inflated

#Fit negative binomial model:
NuCSToFlnb219<-glmmTMB(Total~fImmigrantFreq+(1|Block/Plant),NudatusCS,family=nbinom2)
#DHARMa model diagnostics:
simulationOutput <- simulateResiduals(fittedModel = NuCSToFlnb219, plot = T)
plotQQunif(simulationOutput) #KS and dispersion test sig
testZeroInflation(simulationOutput) #not zero inflated

#Fit negative binomial model:
ziNuCSToFlnb219<-glmmTMB(Total~fImmigrantFreq+(1|Block/Plant),ziformula = ~1,NudatusCS,family=nbinom2)
#DHARMa model diagnostics:
simulationOutput <- simulateResiduals(fittedModel = ziNuCSToFlnb219, plot = T)
plotQQunif(simulationOutput) #all tests NS
testZeroInflation(simulationOutput) #not zero inflated

#does model fit better w/ zi term?
anova(ziNuCSToFlnb219,NuCSToFlnb219) #yes p<0.001
#Wald type 2 chisq:
car::Anova(ziNuCSToFlnb219) #chisq=12.209,df=2,p=0.002232
test(pairs(emmeans(ziNuCSToFlnb219, ~fImmigrantFreq)))

NuCSToFl19df<-ggpredict(ziNuCSToFlnb219,terms="fImmigrantFreq",nsim=500)
plot(NuCSToFl19df)
cld(emmeans(ziNuCSToFlnb219, ~fImmigrantFreq)) #0.05a, 0.25b,0.5a 

#TOTAL PREDS#
GuCSToFlnewx <- as.data.frame(seq(0,3, by = 1))
GuCSToFlnewx$fImmigrantFreq<-ziGuCSToFl19df$x
GuCSToFlnewx$Site= "CS"
GuCSToFlnewx$Species= '"Resident"~italic("M. guttatus")'
GuCSToFlnewx$Total<-ziGuCSToFl19df$predicted
GuCSToFlnewx$se<-ziGuCSToFl19df$std.error
GuCSToFlnewx$conf.low <- ziGuCSToFl19df$conf.low
GuCSToFlnewx$conf.high <- ziGuCSToFl19df$conf.high
GuCSToFlnewx<-subset(GuCSToFlnewx, select = -c(`seq(0, 3, by = 1)`))

dataCS<-data[data$Site == 'CS',]
dataGUCS<-dataCS[dataCS$Species == 'guttatus',]
anno_gutoflseep<-data.frame(xstar = c(1,2,3,4), ystar = c(77.27635,82.44504,79.94076,59.47241),
                            lab = c("A","A","A","A"),
                            fImmigrantFreq = c(0,0.05,0.25,0.5))

GUCSToFlPlot<-ggplot(GuCSToFlnewx, aes(x=fImmigrantFreq, y=Total)) + theme_classic()+
  theme(legend.position = "none") +
  geom_jitter(data=dataGUCS,aes(x=fImmigrantFreq,y=Total),alpha=0.5,fill="#ca0020",colour="white",size=3,shape=21)+
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high), position = position_dodge2(width=0.5),size=1,width=0.5)+
  geom_point(position = position_dodge2(width = 0.9),color="black",fill="#ca0020",size=4,stroke=1,shape=21) +
  ylab(NULL)+theme(legend.position="none")+xlab(NULL)+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),limits=c(0,450))+
  scale_x_discrete(labels=c("0" = "0N:19G\\n0%", "0.05" = "1N:18G\\n5%","0.25" = "4N:15G\\n25%","0.5" = "9N:10G\\n50%"))+
  geom_label(data = anno_gutoflseep, aes(x = xstar,  y = ystar, label = lab,fontface=2,alpha=0.75),fill="white",nudge_y = 0.25,label.size = NA)

NuCSToFlnewx <- as.data.frame(seq(0,2, by = 1))
NuCSToFlnewx$fImmigrantFreq<-NuCSToFl19df$x
NuCSToFlnewx$Site= "CS"
NuCSToFlnewx$Species= '"Immigrant"~italic("M. nudatus")'
NuCSToFlnewx$Total<-NuCSToFl19df$predicted
NuCSToFlnewx$se<-NuCSToFl19df$std.error
NuCSToFlnewx$conf.low <- NuCSToFl19df$conf.low
NuCSToFlnewx$conf.high <- NuCSToFl19df$conf.high
NuCSToFlnewx<-subset(NuCSToFlnewx, select = -c(`seq(0, 2, by = 1)`))

dataNUCS<-dataCS[dataCS$Species == 'nudatus',]
anno_nutoflseep<-data.frame(xstar = c(1,2,3), ystar = c(48.51213,86.68424,44.20205),
                            lab = c("A","B","A"),
                            fImmigrantFreq = c(0.05,0.25,0.5))
NUCSToFlPlot<-ggplot(NuCSToFlnewx, aes(x=fImmigrantFreq, y=Total)) + theme_classic()+
  theme(legend.position = "none") +
  geom_jitter(data=dataNUCS,aes(x=fImmigrantFreq,y=Total),alpha=0.5,fill="#0571b0",colour="white",size=3,shape=24)+
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high), position = position_dodge2(width=0.5),size=1,width=0.5)+
  geom_point(position = position_dodge2(width = 0.9),color="black",fill="#0571b0",size=4,stroke=1,shape=24) +
  ylab("Seeds Per Flower")+theme(legend.position="none")+xlab(NULL)+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),limits=c(0,450))+
  scale_x_discrete(labels=c("0.05" = "1N:18G\\n5%","0.25" = "4N:15G\\n25%","0.5" = "9N:10G\\n50%"))+
  geom_label(data = anno_nutoflseep, aes(x = xstar,  y = ystar, label = lab,fontface=2,alpha=0.75),fill="white",nudge_y = 0.25,label.size = NA)

#####
#nudatus, wash
#Analyze total seeds per flower in the wash - resident M. nudatu only
#Fit poisson model:
NuQVToFlpoi19<-glmmTMB(Total~fImmigrantFreq+(1|Block/Plant),NudatusQV,family=poisson)
#Test for overdispersion:
overdisp_fun(NuQVToFlpoi19) #p=0.0000; poisson model is significantly overdispersed
#DHARMa model diagnostics:
simulationOutput <- simulateResiduals(fittedModel = NuQVToFlpoi19, plot = T)
plotQQunif(simulationOutput) #KS & outlier tests p=0
testZeroInflation(simulationOutput) #zero inflated

#Fit negative binomial model:
NuQVToFlnb219<-glmmTMB(Total~fImmigrantFreq+(1|Block/Plant),NudatusQV,family=nbinom2)
#DHARMa model diagnostics:
simulationOutput <- simulateResiduals(fittedModel = NuQVToFlnb219, plot = T)
plotQQunif(simulationOutput) #KS & dispersion test p=0
testZeroInflation(simulationOutput) #zero inflated

#Fit zi negative binomial model:
ziNuQVToFlnb219<-glmmTMB(Total~fImmigrantFreq+(1|Block/Plant),ziformula = ~1,NudatusQV,family=nbinom2)
#DHARMa model diagnostics:
simulationOutput <- simulateResiduals(fittedModel = ziNuQVToFlnb219, plot = T)
plotQQunif(simulationOutput) #all tests NS
testZeroInflation(simulationOutput) #not zero inflated
#does model fit better w/ zi term?
anova(ziNuQVToFlnb219,NuQVToFlnb219) #yes p<0.001

#Wald typeII Wald-Chi square tests:
car::Anova(ziNuQVToFlnb219)#Chisq=6.4724,df=3,p=0.09076
#Estimate marginal means and bootstrap confidence intervals:
NuQVToFl19df<-ggpredict(ziNuQVToFlnb219,terms="fImmigrantFreq",nsim=500)
plot(NuQVToFl19df)
pairs((emmeans(ziNuQVToFlnb219, ~fImmigrantFreq)))
cld((emmeans(ziNuQVToFlnb219, ~fImmigrantFreq))) #All A

#####
#guttatus, wash
#Analyze total seeds per flower in the wash - immigrant M. guttatus only
#Fit poisson model:
GuQVToFlpoi19<-glmmTMB(Total~fImmigrantFreq+(1|Block/Plant),GuttatusQV,family=poisson)
#Test for overdispersion:
overdisp_fun(GuQVToFlpoi19) #p=0.0000; poisson model is significantly overdispersed
#DHARMa model diagnostics:
simulationOutput <- simulateResiduals(fittedModel = GuQVToFlpoi19, plot = T)
plotQQunif(simulationOutput) #outlier & KS test sig
testZeroInflation(simulationOutput) #=zero inflated

#Fit negative binomial model:
GuQVToFlnb219<-glmmTMB(Total~fImmigrantFreq+(1|Block/Plant),GuttatusQV,family=nbinom2)
#DHARMa model diagnostics:
simulationOutput <- simulateResiduals(fittedModel = GuQVToFlnb219, plot = T)
plotQQunif(simulationOutput) #KS & dispersion test sig
testZeroInflation(simulationOutput) #zero inflated

#Fit zi negative binomial model:
ziGuQVToFlnb219<-glmmTMB(Total~fImmigrantFreq+(1|Block/Plant),ziformula = ~1,GuttatusQV,family=nbinom2)
#DHARMa model diagnostics:
simulationOutput <- simulateResiduals(fittedModel = ziGuQVToFlnb219, plot = T)
plotQQunif(simulationOutput) #all tests NS
testZeroInflation(simulationOutput) #zero inflated
#does model fit better w/ zi term?
anova(ziGuQVToFlnb219,GuQVToFlnb219) #yes p<0.001

#Wald typeII Wald-Chi square tests:
car::Anova(ziGuQVToFlnb219)#Chisq=0.4386,df=2,p=0.8031
#Estimate marginal means and bootstrap confidence intervals:
GuQVToFl19df<-ggpredict(ziGuQVToFlnb219,terms="fImmigrantFreq",nsim=500)
plot(GuQVToFl19df)
pairs((emmeans(ziGuQVToFlnb219, ~fImmigrantFreq)))
cld((emmeans(ziGuQVToFlnb219, ~fImmigrantFreq))) #All A

#TOTAL PREDS#
GuQVToFlnewx <- as.data.frame(seq(0,2, by = 1))
GuQVToFlnewx$fImmigrantFreq<-GuQVToFl19df$x
GuQVToFlnewx$Site="QV"
GuQVToFlnewx$Species= '"Immigrant"~italic("M. guttatus")'
GuQVToFlnewx$Total<-GuQVToFl19df$predicted
GuQVToFlnewx$conf.low <- GuQVToFl19df$conf.low
GuQVToFlnewx$conf.high <- GuQVToFl19df$conf.high
GuQVToFlnewx<-subset(GuQVToFlnewx, select = -c(`seq(0, 2, by = 1)`))

dataQV<-data[data$Site == 'QV',]
dataGUQV<-dataQV[dataQV$Species == 'guttatus',]
anno_gutoflwash<-data.frame(xstar = c(1,2,3), ystar = c(107.56371,97.01046,73.16126),
                            lab = c("A","A","A"),
                            fImmigrantFreq = c(0.05,0.25,0.5))
GUQVToFlPlot<-ggplot(GuQVToFlnewx, aes(x=fImmigrantFreq, y=Total)) + theme_classic()+
  theme(legend.position = "none") +
  geom_jitter(data=dataGUQV,aes(x=fImmigrantFreq,y=Total),alpha=0.5,fill="#ca0020",colour="white",size=3,shape=21)+
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high), position = position_dodge2(width=0.5),size=1,width=0.5)+
  geom_point(position = position_dodge2(width = 0.9),color="black",fill="#ca0020",size=4,stroke=1,shape=21) +
  ylab("Seeds Per Flower")+theme(legend.position="none")+xlab(NULL)+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),limits=c(0,450))+
  scale_x_discrete(labels=c("0.05" = "1G:18N\\n5%","0.25" = "4G:15N\\n25%","0.5" = "9G:10N\\n50%"))+
  geom_label(data = anno_gutoflwash, aes(x = xstar,  y = ystar, label = lab,fontface=2,alpha=0.75),fill="white",nudge_y = 0.25,label.size = NA)

NuQVToFlnewx <- as.data.frame(seq(0,3, by = 1))
NuQVToFlnewx$fImmigrantFreq<-NuQVToFl19df$x
NuQVToFlnewx$Site="QV"
NuQVToFlnewx$Species= '"Resident"~italic("M. nudatus")'
NuQVToFlnewx$Total<-NuQVToFl19df$predicted
NuQVToFlnewx$conf.low <- NuQVToFl19df$conf.low
NuQVToFlnewx$conf.high <- NuQVToFl19df$conf.high
NuQVToFlnewx<-subset(NuQVToFlnewx, select = -c(`seq(0, 3, by = 1)`))
QVToFlnewx <- bind_rows(NuQVToFlnewx,GuQVToFlnewx)

dataNUQV<-dataQV[dataQV$Species == 'nudatus',]
anno_nutoflwash<-data.frame(xstar = c(1,2,3,4), ystar = c(65.93237,74.47355,59.51415,57.87583),
                            lab = c("A","A","A","A"),
                            fImmigrantFreq = c(0,0.05,0.25,0.5))
NUQVToFlPlot<-ggplot(NuQVToFlnewx, aes(x=fImmigrantFreq, y=Total)) + theme_classic()+
  theme(legend.position = "none") +
  geom_jitter(data=dataNUQV,aes(x=fImmigrantFreq,y=Total),alpha=0.5,fill="#0571b0",colour="white",size=3,shape=24)+
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high), position = position_dodge2(width=0.5),size=1,width=0.5)+
  geom_point(position = position_dodge2(width = 0.9),color="black",fill="#0571b0",size=4,stroke=1,shape=24) +
  ylab(NULL)+theme(legend.position="none")+xlab(NULL)+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),limits=c(0,450))+
  scale_x_discrete(labels=c("0" = "0G:19N\\n0%","0.05" = "1G:18N\\n5%","0.25" = "4G:15N\\n25%","0.5" = "9G:10N\\n50%"))+
  geom_label(data = anno_nutoflwash, aes(x = xstar,  y = ystar, label = lab,fontface=2,alpha=0.75),fill="white",nudge_y = 0.25,label.size = NA)


#Seed Viability#
#guttatus, seep
GuCSPropVia19<-glmmTMB(cbind(Normal,Flat)~fImmigrantFreq+(1|Block),sumGuttatusCS,family=binomial)
#DHARMa model diagnostics:
simulationOutput <- simulateResiduals(fittedModel = GuCSPropVia19, plot = T)
plotQQunif(simulationOutput) #KS & outlier tests: p=0

#betabinomial model
GuCSPropVia19<-glmmTMB(cbind(Normal,Flat)~fImmigrantFreq+(1|Block),sumGuttatusCS,family=betabinomial())
#DHARMa model diagnostics:
simulationOutput <- simulateResiduals(fittedModel = GuCSPropVia19, plot = T)
plotQQunif(simulationOutput) #all tests NS

#Wald typeII Wald-Chi square tests :
car::Anova(GuCSPropVia19)#Chisq=21.563,df=3,p=8.042e-05
#Tukey posthoc test
pairs(emmeans(GuCSPropVia19, ~fImmigrantFreq))
cld(emmeans(GuCSPropVia19, ~fImmigrantFreq))#0a,0.05a, 0.25ab,0.5b

#Estimate marginal means and bootstrap confidence intervals:
GuCSPropVia19df<-ggpredict(GuCSPropVia19,terms="fImmigrantFreq",nsim=500)
plot(GuCSPropVia19df)

#nudatus, seep
NuCSPropVia19<-glmmTMB(cbind(Normal,Flat)~fImmigrantFreq+(1|Block),sumNudatusCS,family=binomial)
#DHARMa model diagnostics:
simulationOutput <- simulateResiduals(fittedModel = NuCSPropVia19, plot = T)
plotQQunif(simulationOutput) #not overdispersed
#Wald typeII Wald-Chi square tests :
car::Anova(NuCSPropVia19)#Chisq=1.5171,df=2,p=0.4684
#Tukey posthoc test
pairs(emmeans(NuCSPropVia19, ~fImmigrantFreq))#all p>0.47
cld(emmeans(NuCSPropVia19, ~fImmigrantFreq))#all A

#Estimate marginal means and bootstrap confidence intervals:
NuCSPropVia19df<-ggpredict(NuCSPropVia19,terms="fImmigrantFreq",nsim=500)
plot(NuCSPropVia19df)

#PropVia PREDS#
GuCSVinewx <- as.data.frame(seq(0,3, by = 1))
GuCSVinewx$fImmigrantFreq<-GuCSPropVia19df$x
GuCSVinewx$Site= "CS"
GuCSVinewx$Species= '"Resident"~italic("M. guttatus")'
GuCSVinewx$Normal<-GuCSPropVia19df$predicted
GuCSVinewx$se<-GuCSPropVia19df$std.error
GuCSVinewx$conf.low <- GuCSPropVia19df$conf.low
GuCSVinewx$conf.high <- GuCSPropVia19df$conf.high
GuCSVinewx<-subset(GuCSVinewx, select = -c(`seq(0, 3, by = 1)`))

anno_guviseep<-data.frame(xstar = c(1,2,3,4), ystar = c(0.8800268,0.8586269,0.8039875,0.6990100),
                          lab = c("A","A","AB","B"),
                          fImmigrantFreq = c(0,0.05,0.25,0.5))

GUCSPropViaPlot<-ggplot(GuCSVinewx, aes(x=fImmigrantFreq, y=Normal)) + theme_classic()+ ylim(0,1)+
  theme(legend.position = "none") +
  geom_jitter(data=sumGuttatusCS,aes(x=fImmigrantFreq,y=(Normal/Total)),alpha=0.5,fill="#ca0020",colour="white",size=3,shape=21)+
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high), position = position_dodge2(width=0.5),size=1,width=0.5)+
  geom_point(position = position_dodge2(width = 0.9),color="black",fill="#ca0020",size=4,stroke=1,shape=21) +
  ylab(NULL)+theme(legend.position="none")+xlab("Seep Immigrant Frequency")+
  scale_x_discrete(labels=c("0" = "0N:19G\\n0%", "0.05" = "1N:18G\\n5%","0.25" = "4N:15G\\n25%","0.5" = "9N:10G\\n50%"))+
  geom_label(data = anno_guviseep, aes(x = xstar,  y = ystar, label = lab,fontface=2,alpha=0.75),fill="white",nudge_y = 0.1,label.size = NA)

NuCSVinewx <- as.data.frame(seq(0,2, by = 1))
NuCSVinewx$fImmigrantFreq<-NuCSPropVia19df$x
NuCSVinewx$Site= "CS"
NuCSVinewx$Species= '"Immigrant"~italic("M. nudatus")'
NuCSVinewx$Normal<-NuCSPropVia19df$predicted
NuCSVinewx$se<-NuCSPropVia19df$std.error
NuCSVinewx$conf.low <- NuCSPropVia19df$conf.low
NuCSVinewx$conf.high <- NuCSPropVia19df$conf.high
NuCSVinewx<-subset(NuCSVinewx, select = -c(`seq(0, 2, by = 1)`))
CSVinewx <- bind_rows(GuCSVinewx,NuCSVinewx)

anno_nuviseep<-data.frame(xstar = c(1,2,3), ystar = c(0.4421123,0.3978810,0.7084927),
                          lab = c("A","A","A"),
                          fImmigrantFreq = c(0.05,0.25,0.5))

NUCSPropViaPlot<-ggplot(NuCSVinewx, aes(x=fImmigrantFreq, y=Normal)) + theme_classic()+ ylim(0,1)+
  theme(legend.position = "none") +
  geom_jitter(data=sumNudatusCS,aes(x=fImmigrantFreq,y=(Normal/Total)),alpha=0.5,fill="#0571b0",colour="white",size=3,shape=24)+
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high), position = position_dodge2(width=0.5),size=1,width=0.5)+
  geom_point(position = position_dodge2(width = 0.9),color="black",fill="#0571b0",size=4,stroke=1,shape=24) +
  ylab("Seed Viability")+theme(legend.position="none")+xlab("Seep Immigrant Frequency")+
  scale_x_discrete(labels=c("0.05" = "1N:18G\\n5%","0.25" = "4N:15G\\n25%","0.5" = "9N:10G\\n50%"))+
  geom_label(data = anno_nuviseep, aes(x = xstar,  y = ystar, label = lab,fontface=2,alpha=0.75),fill="white",nudge_y = 0.1,label.size = NA)


#####
#####
#guttatus, wash
GuQVPropVia19<-glmmTMB(cbind(Normal,Flat)~fImmigrantFreq+(1|Block),sumGuttatusQV,family=binomial)
#DHARMa model diagnostics:
simulationOutput <- simulateResiduals(fittedModel = GuQVPropVia19, plot = T)
plotQQunif(simulationOutput) #KS test=0.04, dispersion & outlier tests NS

#Wald typeII Wald-Chi square tests :
car::Anova(GuQVPropVia19)#Chisq=10.602,df=2,p=0.004987
#Tukey posthoc test
pairs(emmeans(GuQVPropVia19, ~fImmigrantFreq))
cld(emmeans(GuQVPropVia19, ~fImmigrantFreq))#0.05a, 0.25ab,0.5b

#Estimate marginal means and bootstrap confidence intervals:
GuQVPropVia19df<-ggpredict(GuQVPropVia19,terms="fImmigrantFreq",nsim=500)
plot(GuQVPropVia19df)


#####
#nudatus, wash
NuQVPropVia19<-glmmTMB(cbind(Normal,Flat)~fImmigrantFreq+(1|Block),sumNudatusQV,family=binomial)
#DHARMa model diagnostics:
simulationOutput <- simulateResiduals(fittedModel = NuQVPropVia19, plot = T)
plotQQunif(simulationOutput) #all tests sig

#betabinomial model
NuQVPropVia19<-glmmTMB(cbind(Normal,Flat)~fImmigrantFreq+(1|Block),sumNudatusQV,family=betabinomial())
#DHARMa model diagnostics:
simulationOutput <- simulateResiduals(fittedModel = NuQVPropVia19, plot = T)
plotQQunif(simulationOutput) #all tests NS

#Wald typeII Wald-Chi square tests :
car::Anova(NuQVPropVia19)#Chisq=11.402,df=3,p=0.009737
#Tukey posthoc test
pairs(emmeans(NuQVPropVia19, ~fImmigrantFreq))
cld(emmeans(NuQVPropVia19, ~fImmigrantFreq))#0a,0.05ab, 0.25b,0.5ab

#Estimate marginal means and bootstrap confidence intervals:
NuQVPropVia19df<-ggpredict(NuQVPropVia19,terms="fImmigrantFreq",nsim=500)
plot(NuQVPropVia19df)

#Normal PREDS#
GuQVVinewx <- as.data.frame(seq(0,2, by = 1))
GuQVVinewx$fImmigrantFreq<-GuQVPropVia19df$x
GuQVVinewx$Site="QV"
GuQVVinewx$Species= '"Immigrant"~italic("M. guttatus")'
GuQVVinewx$Normal<-GuQVPropVia19df$predicted
GuQVVinewx$conf.low <- GuQVPropVia19df$conf.low
GuQVVinewx$conf.high <- GuQVPropVia19df$conf.high
GuQVVinewx<-subset(GuQVVinewx, select = -c(`seq(0, 2, by = 1)`))

anno_guviwash<-data.frame(xstar = c(1,2,3), ystar = c(0.2658431,0.3799397,0.6574370),
                          lab = c("A","AB","B"),
                          fImmigrantFreq = c(0.05,0.25,0.5))
GUQVPropViaPlot<-ggplot(GuQVVinewx, aes(x=fImmigrantFreq, y=Normal)) + theme_classic()+ ylim(0,1)+
  theme(legend.position = "none") +
  geom_jitter(data=sumGuttatusQV,aes(x=fImmigrantFreq,y=(Normal/Total)),alpha=0.5,fill="#ca0020",colour="white",size=3,shape=21)+
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high), position = position_dodge2(width=0.5),size=1,width=0.5)+
  geom_point(position = position_dodge2(width = 0.9),color="black",fill="#ca0020",size=4,stroke=1,shape=21) +
  ylab("Seed Viability")+theme(legend.position="none")+xlab("Wash Immigrant Frequency")+
  scale_x_discrete(labels=c("0.05" = "1G:18N\\n5%","0.25" = "4G:15N\\n25%","0.5" = "9G:10N\\n50%"))+
  geom_label(data = anno_guviwash, aes(x = xstar,  y = ystar, label = lab,fontface=2,alpha=0.75),fill="white",nudge_y = 0.1,label.size = NA)

NuQVVinewx <- as.data.frame(seq(0,3, by = 1))
NuQVVinewx$fImmigrantFreq<-NuQVPropVia19df$x
NuQVVinewx$Site="QV"
NuQVVinewx$Species= '"Resident"~italic("M. nudatus")'
NuQVVinewx$Normal<-NuQVPropVia19df$predicted
NuQVVinewx$conf.low <- NuQVPropVia19df$conf.low
NuQVVinewx$conf.high <- NuQVPropVia19df$conf.high
NuQVVinewx<-subset(NuQVVinewx, select = -c(`seq(0, 3, by = 1)`))

anno_nuviwash<-data.frame(xstar = c(1,2,3,4), ystar = c(0.8727001,0.8222331,0.7502716,0.8224417),
                          lab = c("A","AB","B","AB"),
                          fImmigrantFreq = c(0,0.05,0.25,0.5))

NUQVPropViaPlot<-ggplot(NuQVVinewx, aes(x=fImmigrantFreq, y=Normal)) + theme_classic()+ ylim(0,1)+
  theme(legend.position = "none") +
  geom_jitter(data=sumNudatusQV,aes(x=fImmigrantFreq,y=(Normal/Total)),alpha=0.5,fill="#0571b0",colour="white",size=3,shape=24)+
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high), position = position_dodge2(width=0.5),size=1,width=0.5)+
  geom_point(position = position_dodge2(width = 0.9),color="black",fill="#0571b0",size=4,stroke=1,shape=24) +
  ylab(NULL)+theme(legend.position="none")+xlab("Wash Immigrant Frequency")+
  scale_x_discrete(labels=c("0" = "0G:19N\\n5%","0.05" = "1G:18N\\n5%","0.25" = "4G:15N\\n25%","0.5" = "9G:10N\\n50%"))+
  geom_label(data = anno_nuviwash, aes(x = xstar,  y = ystar, label = lab,fontface=2,alpha=0.75),fill="white",nudge_y = 0.1,label.size = NA)

anno_viwash<-data.frame(xstar = c(1,2,3,1,2,3,4), ystar = c(0.2658431,0.3799397,0.6574370,0.8727001,0.8222331,0.7502716,0.8224417),
                        lab = c("A","AB","B","A","AB","B","AB"),
                        fImmigrantFreq = c(0.05,0.25,0.5,0,0.05,0.25,0.5),
                        Species=c('"Immigrant"~italic("M. guttatus")','"Immigrant"~italic("M. guttatus")','"Immigrant"~italic("M. guttatus")','"Resident"~italic("M. nudatus")','"Resident"~italic("M. nudatus")','"Resident"~italic("M. nudatus")','"Resident"~italic("M. nudatus")'))


Wash <- GUQVFlPlot+NUQVFlPlot+
  GUQVToFlPlot+NUQVToFlPlot+
  GUQVPropViaPlot+xlab(NULL)+NUQVPropViaPlot+xlab(NULL)+
  plot_layout(ncol=2)+plot_annotation(tag_levels = 'A')
wt <- patchwork::patchworkGrob(Wash)
gridExtra::grid.arrange(wt, bottom = "Wash Immigrant Frequency")
#Save as 6x7; Fig4-Wash.pdf

Seep<-NUCSFlPlot+GUCSFlPlot+
  NUCSToFlPlot+GUCSToFlPlot+
  NUCSPropViaPlot+xlab(NULL)+GUCSPropViaPlot+xlab(NULL)+
  plot_layout(ncol=2)+plot_annotation(tag_levels = 'A')
st <- patchwork::patchworkGrob(Seep)
gridExtra::grid.arrange(st, bottom = "Seep Immigrant Frequency")
#Save as 6x7; Fig5-Seep.pdf




