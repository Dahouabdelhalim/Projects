library(Hmisc)
library(lmtest)
library(lme4)
library(lmerTest)
library(emmeans)
library(car)
library(vegan)
library(MuMIn)
library(glmmTMB)
library(gridExtra)
library(RVAideMemoire)
library(DHARMa)
library(optimx)
library(patchwork)
library(tidyverse)

exp16 <-read_csv("Maron_etal_2020_JEcol_Data.csv")
exp16$Year <- factor(exp16$Year, levels=c("Seedlings","Adults_Y3"))  ## reorder factor levels

levels(exp16$Year) <- c("Seedlings","Adults")

### grass cover models, filter data set to just one species so extra lines are not included
gr1 <- lmer(sqrt(Graminoid)~Rodents*Competition+(1|SiteCode), data=exp16 %>% filter(Species=='ACMI'&Year=='Adults'))
anova(gr1)
emmeans(gr1, ~Competition, transform="response")

lit1 <- lmer(Litterdepth~Rodents*Competition+(1|SiteCode), data=exp16 %>% filter(Species=='ACMI'&Year=='Adults'))
anova(lit1)
emmeans(lit1, ~Competition)

##########################################################################################################################
## load species trait data
fec <- read.csv("Maron_etal_2020_JEcol_Trait.csv", header=T) 
fec$SeedSize_mg <- fec$SeedSize*1000
fec$lSS <- log10(fec$SeedSize_mg)
fec$Fec <- log10(fec$Fecundity)
fec$LMA <- fec$LMA*1000

## manually estimate effect sizes

#calculate rodent effect
effect1 <- exp16 %>% group_by(Species,Year,Rodents) %>% summarize(numb=mean(numb,na.rm=T), prec1=mean(prec1,na.rm=T)) %>% 
  pivot_wider(id_col=c("Species","Year"), names_from = "Rodents", values_from = c("numb","prec1"))

effect1$p_lrr <-  log(effect1$'prec1_+'/effect1$'prec1_-')
effect <- full_join(effect1,fec,by="Species")

#calculate competition effect
effect2 <- exp16 %>% group_by(Species,Year,Competition) %>% summarize(numb=mean(numb,na.rm=T), prec1=mean(prec1,na.rm=T)) %>% 
  pivot_wider(id_col=c("Species","Year"), names_from = "Competition", values_from = c("numb","prec1"))

effect2$p_lrrc <-  log(effect2$'prec1_+'/effect2$'prec1_-')

effectc <- full_join(effect2,fec,by="Species")

### construct and analyses models off competition effect size
lmcomp1 <- lmer(p_lrrc~lSS*Year+(1|Species), data=effectc)
anova(lmcomp1)
emtrends(lmcomp1, var="lSS", ~Year)

lmcomp2 <- lmer(p_lrrc~LMA*Year+(1|Species), data=effectc)
anova(lmcomp2)

lmcomp3 <- lmer(p_lrrc~Height*Year+(1|Species), data=effectc)
anova(lmcomp3)

lmrods <- lmer(p_lrr~lSS+Year+(1|Species), data=effect)
anova(lmrods)
summary(lmrods)
confint(lmrods)
emtrends(lmrods, var="lSS")


### Figure 1 panels ################################################################################
p1a <- ggplot(data=effectc, aes(x=lSS, y=p_lrrc, fill=Year, linetype=Year))+
  #geom_hline(aes(yintercept=0), linetype="longdash", color="grey", lwd=1.5)+
  geom_smooth(method="lm", se=T, color="black", fill="grey")+geom_point(size=3, stroke=1.5, pch=21)+
  scale_fill_manual(values=c("white","black"), label=c("Cumulative seedlings", "Adults"))+
  scale_linetype_manual(values=c("dashed","solid"), label=c("Cumulative seedlings",  "Adults"))+
  scale_x_continuous("log10(Seed size mg)")+
  scale_y_continuous("Competition effect size", limits=c(-3,1))+
  theme(legend.position=c(.6,.125),legend.key=element_blank(),legend.title = element_blank(), 
        legend.background = element_rect(fill="white", size=1, linetype="solid", color="black"))+
  theme(axis.line = element_line(size = 0, colour = "black"))+
  theme(panel.background = element_blank(), plot.margin=margin(.5,.5,.5,.7, "cm"))+
  theme(panel.border = element_rect(color="black", fill=NA, size=1.5)) +
  theme(text = element_text(size=16), axis.text.x=element_text(size=14))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(strip.background =  element_blank(), strip.text = element_text(size=16))+
  theme(axis.ticks.length=unit(0.2, "cm"),  
        axis.text.x = element_text(margin=margin(5,5,5,5,"pt"),colour="black"),
        axis.text.y = element_text(margin=margin(5,5,5,5,"pt"),colour="black"))+
  theme(legend.key.width = unit(2,"cm"))

p1b <- ggplot(data=effect)+
  #geom_hline(aes(yintercept=0), linetype="longdash", color="grey", lwd=1.5)+
  geom_smooth(data=effect, aes(x=lSS, y=p_lrr), method="lm", se=T, color="black")+
  geom_point(data=effect, aes(x=lSS, y=p_lrr, fill=Year), size=3, stroke=1.5, pch=21)+
  scale_fill_manual(values=c("white","black"))+
  scale_linetype_manual(values=c("dashed","solid"))+
  xlab("log10(Seed size mg)")+
  scale_y_continuous("Rodent effect size", limits=c(-3,1))+
  theme(legend.position="none",legend.key=element_blank(), 
        legend.background = element_rect(fill="white", size=1, linetype="solid", color="black"))+
  theme(axis.line = element_line(size = 0, colour = "black"))+
  theme(panel.background = element_blank(), plot.margin=margin(.5,.5,.5,.7, "cm"))+
  theme(panel.border = element_rect(color="black", fill=NA, size=1.5)) +
  theme(text = element_text(size=16), axis.text.x=element_text(size=14))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(strip.background =  element_blank(), strip.text = element_text(size=16))+
  theme(axis.ticks.length=unit(0.2, "cm"),  
        axis.text.x = element_text(margin=margin(5,5,5,5,"pt"),colour="black"),
        axis.text.y = element_text(margin=margin(5,5,5,5,"pt"),colour="black"))

# create .tiff of Figure 1
tiff("Fig1.tiff", width=10, height=5, units='in',
     res=600, compression='lzw')
p1a+p1b
dev.off()


#############################################################################################################
### FIGURE 3 ANALYSIS
### abundance across productivity gradient
exp16$obs <- 1:length(exp16$numb)

## construct model to estimate response of each species across productivity gradient
alm2 <- glmer(numb~prod_m2*Species*Year+(1|SiteCode/Species)+(1|obs), data=exp16 %>% filter(Trt=='++'), family="poisson",
              control = glmerControl(calc.derivs = FALSE,optimizer=c("bobyqa","Nelder_Mead"),optCtrl=list(maxfun=2e5)), nAGQ = 0)
overdisp.glmer(alm2) # overdispersion ratio calculator from RVAideMemoire

sim_alm2 <- simulateResiduals(fittedModel = alm2, n = 250)
plot(sim_alm2)

Anova(alm2)
prod1 <- as.data.frame(emtrends(alm2, ~Species:Year, var='prod_m2'))

pt1 <- full_join(fec,prod1,by="Species")

## analyze responses across productivity gradient depending on species' traits
plm1 <- lmer(prod_m2.trend~lSS*Year+(1|Species), data=pt1 )
anova(plm1)
summary(plm1)
confint(plm1)

plm2 <- lmer(prod_m2.trend~LMA*Year+(1|Species), data=pt1 )
anova(plm2)
summary(plm2)

plm3 <- lmer(prod_m2.trend~Height*Year+(1|Species), data=pt1 )
anova(plm3)
summary(plm3)
confint(plm3)

### create panels for productivity figure 3
p3a <- ggplot(data=pt1)+
  #geom_hline(aes(yintercept=0), linetype="longdash", color="grey", lwd=1.5)+
  geom_smooth(data=pt1, aes(x=lSS, y=prod_m2.trend*100), method="lm", se=T, color="black")+
  geom_point(data=pt1, aes(x=lSS, y=prod_m2.trend*100, fill=Year), size=3, stroke=1.5, pch=21)+
  scale_fill_manual(values=c("white","black"), labels=c("Cumulative seedlings","Adults"))+
  scale_linetype_manual(values=c("dashed","solid"))+
  xlab("log10(Seed size mg)")+
  scale_y_continuous("Response across productivity gradient", limits=c(-2.5,1.1))+
  theme(legend.position="none")+
  theme(axis.line = element_line(size = 0, colour = "black"))+
  theme(panel.background = element_blank(), plot.margin=margin(.5,.5,.5,.7, "cm"))+
  theme(panel.border = element_rect(color="black", fill=NA, size=1.5)) +
  theme(text = element_text(size=16), axis.text.x=element_text(size=14))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(strip.background =  element_blank(), strip.text = element_text(size=16))+
  theme(axis.ticks.length=unit(0.2, "cm"),  
        axis.text.x = element_text(margin=margin(5,5,5,5,"pt"),colour="black"),
        axis.text.y = element_text(margin=margin(5,5,5,5,"pt"),colour="black"))

p3b <- ggplot(data=pt1)+
  #geom_hline(aes(yintercept=0), linetype="longdash", color="grey", lwd=1.5)+
  geom_smooth(data=pt1, aes(x=Height, y=prod_m2.trend*100), method="lm", se=T, color="black")+
  geom_point(data=pt1, aes(x=Height, y=prod_m2.trend*100, fill=Year), size=3, stroke=1.5, pch=21)+
  scale_fill_manual(values=c("white","black"), labels=c("Cumulative seedlings","Adults"))+
  scale_linetype_manual(values=c("dashed","solid"))+
  xlab("Height (cm)")+
  scale_y_continuous("", limits=c(-2.5,1.1))+
  theme(legend.position=c(.7,.15),legend.key=element_blank(), legend.title = element_blank(), 
        legend.background = element_rect(fill="white", size=1, linetype="solid", color="black"))+
  theme(axis.line = element_line(size = 0, colour = "black"))+
  theme(panel.background = element_blank(), plot.margin=margin(.5,.5,.5,.7, "cm"))+
  theme(panel.border = element_rect(color="black", fill=NA, size=1.5)) +
  theme(text = element_text(size=16), axis.text.x=element_text(size=14))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(strip.background =  element_blank(), strip.text = element_text(size=16))+
  theme(axis.ticks.length=unit(0.2, "cm"),  
        axis.text.x = element_text(margin=margin(5,5,5,5,"pt"),colour="black"),
        axis.text.y = element_text(margin=margin(5,5,5,5,"pt"),colour="black"))

p3c <- ggplot(data=pt1)+
  #geom_hline(aes(yintercept=0), linetype="longdash", color="grey", lwd=1.5)+
  #geom_smooth(data=pt1, aes(x=LMA, y=prod_m2.trend*100), method="lm", se=T, color="black")+
  geom_point(data=pt1, aes(x=LMA, y=prod_m2.trend*100, fill=Year), size=3, stroke=1.5, pch=21)+
  scale_fill_manual(values=c("white","black"), labels=c("Cumulative seedlings","Adults"))+
  scale_linetype_manual(values=c("dashed","solid"))+
  xlab(expression('LMA g m'^-2))+
  scale_y_continuous("", limits=c(-2.5,1.1))+
  theme(legend.position="none")+
  theme(axis.line = element_line(size = 0, colour = "black"))+
  theme(panel.background = element_blank(), plot.margin=margin(.5,.5,.5,.7, "cm"))+
  theme(panel.border = element_rect(color="black", fill=NA, size=1.5)) +
  theme(text = element_text(size=16), axis.text.x=element_text(size=14))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(strip.background =  element_blank(), strip.text = element_text(size=16))+
  theme(axis.ticks.length=unit(0.2, "cm"),  
        axis.text.x = element_text(margin=margin(5,5,5,5,"pt"),colour="black"),
        axis.text.y = element_text(margin=margin(5,5,5,5,"pt"),colour="black"))

tiff("Fig3.tiff", width=15, height=5, units='in',
     res=600, compression='lzw')
p3a+p3b+p3c
dev.off()



#######################################################################################################
#####  FIGURE 5 ANALYSIS
##### adult abundance models
fec1 <- glmer(numb~log10(Fecundity)*Rodents*log10(Fecundity)*Competition+(1|SiteCode/Rodents/Competition)+(1|obs), 
              data=exp16 %>% filter(Year=='Adults'), family="poisson",
              control = glmerControl(calc.derivs = FALSE,optimizer="bobyqa",optCtrl=list(maxfun=2e5)), nAGQ = 1)
overdisp.glmer(fec1) # overdispersion ratio calculator from RVAideMemoire

sim_odpm1 <- simulateResiduals(fittedModel = fec1, n = 250)
plot(sim_odpm1)

summary(fec1)
Anova(fec1)
confint(emmeans(fec1, ~Rodents, type="response"))

emtrends(fec1, ~Competition, var='log10(Fecundity)')
confint(emtrends(fec1, ~Competition:Rodents, var='log10(Fecundity)'), level=.95)
r.squaredGLMM(fec1)

## summarize data for plotting Figure 5
exp16sum <- exp16 %>% filter(Year=="Adults") %>% group_by(Species,Competition,Rodents) %>%
  summarize(Fecundity=mean(Fecundity,na.rm=T),lSS=mean(log10(SeedSize_mg),na.rm=T),numb=mean(numb,na.rm=T)) %>% 
  mutate(Abund_SS=(numb*lSS),log_Abund_SS=log(numb)*lSS)
exp16sum$Comp1 <- factor(exp16sum$Competition, levels=c('-','+'))
complabs <- c("-Competition","+Competition") 
names(complabs) <- c("-","+")

p5 <- ggplot(data=exp16sum)+
  geom_point(data=exp16sum,aes(x=Fecundity, y=numb, shape=Rodents), size=3, stroke=1.5)+
  geom_smooth(data=exp16sum,aes(x=Fecundity, y=round(numb,0), linetype=Rodents), color="black", method="glm",
                           method.args = list(family = 'poisson'))+
  facet_grid(~Comp1, labeller = labeller(Comp1=complabs))+
  scale_linetype_manual(values=c("solid","dashed"))+
  scale_shape_manual(values=c(16,1))+
  scale_x_continuous(trans="log10",breaks = c(100,1000,10000),labels=c(100,1000,10000),limits=c(25,11000))+
  scale_y_continuous("Mean number of adult plants")+
  theme(legend.position=c(.75,.85),legend.key=element_rect(fill="white"), legend.box.background = element_rect(fill=NA), 
        legend.background = element_rect(fill="white", size=1, linetype="solid", color="black"))+
  theme(axis.line = element_line(size = 0, colour = "black"))+
  theme(panel.background = element_blank(), plot.margin=margin(.5,.5,.5,.7, "cm"))+
  theme(panel.border = element_rect(color="black", fill=NA, size=1.5)) +
  theme(text = element_text(size=16), axis.text.x=element_text(size=14))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(strip.background =  element_blank(), strip.text = element_text(size=16))+
  theme(axis.ticks.length=unit(0.2, "cm"),  
        axis.text.x = element_text(margin=margin(5,5,5,5,"pt"),colour="black"),
        axis.text.y = element_text(margin=margin(5,5,5,5,"pt"),colour="black"))+
  theme(legend.key.width = unit(1.75,"cm"))+
  guides(color=guide_legend(override.aes=list(color=NA,fill=NA)))

tiff("Fig5.tiff", width=8, height=5, units='in',
     res=600, compression='lzw')
p5
dev.off()


