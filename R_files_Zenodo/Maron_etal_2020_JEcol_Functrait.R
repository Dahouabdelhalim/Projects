library(lattice)
library(lmtest)
library(lme4)
library(lmerTest)
library(emmeans)
library(tidyverse)
library(FD)
library(car)
library(vegan)
library(MuMIn)
library(patchwork)
library(gridExtra)

## load data file
fd2 <-read_csv("Maron_etal_2020_JEcol_FuctTrait.csv")

## construct model and analysis for CWM LMA
mod1 <- lmer(CWM_LMA*1000~Rodents*Competition*Year+(1|SiteCode/Rodents/Competition), data=fd2,
             control = lmerControl(calc.derivs = FALSE,optCtrl=list(maxfun=2e5)))
anova(mod1)
summary(mod1)
r.squaredGLMM(mod1) 
confint(emmeans(mod1, pairwise~Competition|Year))

## construct model and analysis for CWM SeedSize
mod2 <- lmer(CWM_SeedSize~Rodents*Competition*Year+(1|SiteCode/Rodents/Competition), data=fd2,
             control = lmerControl(calc.derivs = FALSE,optCtrl=list(maxfun=2e5)))
anova(mod2)
summary(mod2)
r.squaredGLMM(mod2) 

## construct model and analysis for CWM SeedSize
mod3 <- lmer(CWM_Height~Rodents*Competition*Year+(1|SiteCode/Rodents/Competition), data=fd2,
             control = lmerControl(calc.derivs = FALSE,optCtrl=list(maxfun=2e5)))
anova(mod3)
summary(mod3)
r.squaredGLMM(mod3) 


#############################################################
### analyze CWM trait data across productivity gradient
mod1a <- lmer(CWM_SeedSize~scale(productivity_m2)*Year+(1|SiteCode), data=fd2 %>% filter(Competition=='+',Rodents=='+'))
anova(mod1a)
summary(mod1a)
emtrends(mod1a, ~Year, var='productivity_m2')
r.squaredGLMM(mod1a) 

mod1b <- lmer(CWM_LMA~scale(productivity_m2)*Year+(1|SiteCode), data=fd2 %>% filter(Competition=='+',Rodents=='+'))
anova(mod1b)
summary(mod1b)
emtrends(mod1b, ~Competition|Year, var='productivity_m2')
r.squaredGLMM(mod1b) 

mod1c <- lmer(CWM_Height~scale(productivity_m2)*Year+(1|SiteCode), data=fd2 %>% filter(Competition=='+',Rodents=='+'))
anova(mod1c)
summary(mod1c)
emtrends(mod1c, ~Year, var='productivity_m2')
r.squaredGLMM(mod1c) 

#############################################################################################
### code to create Fig 2 -- CWM plots
fd2$Year <- factor(fd2$Year, levels=c("Seedlings","Adults_Y3"))  ## reorder factor levels

fd2mean <- fd2 %>% group_by(Competition,Rodents,Year) %>%  summarize(CWM_sm=mean(CWM_SeedSize,na.rm=T),CWM_ssd=sd(CWM_SeedSize,na.rm=T),
                                                               CWM_lma=mean(CWM_LMA,na.rm=T),CWM_lsd=sd(CWM_LMA,na.rm=T),
                                                               CWM_hgt=mean(CWM_Height,na.rm=T),CWM_hsd=sd(CWM_Height,na.rm=T))
levels(fd2mean$Year) <- c("Seedlings","Adults")

p2a <- ggplot(fd2mean, aes(x=Competition, y=CWM_lma, group=Rodents, fill=Rodents)) + 
  geom_pointrange(position=position_dodge(width=.5), size=1.5, pch=21, stroke=2,
                  aes(ymin=CWM_lma-(CWM_lsd/sqrt(9)),ymax=CWM_lma+(CWM_lsd/sqrt(9))) )+
  #geom_text(x=.5,y=1.2,label="A", size=8)+
  scale_fill_manual(values=c("black","grey75"), labels=c("No","Yes"))+
  scale_x_discrete(name="Competition", labels=c("No","Yes"))+
  #ylab(bquote(~CWM[SeedSize]))+
  ylab(bquote(~CWM[LMA]))+ facet_wrap(~Year)+
  theme(legend.position=c(.85,.2),legend.key=element_blank(), 
        legend.background = element_rect(fill="white", size=1, linetype="solid", color="black"))+
  theme(axis.line = element_line(size = 0, colour = "black"))+
  theme(panel.background = element_blank(), plot.margin=margin(.5,.5,.5,.5, "cm"))+
  theme(panel.border = element_rect(color="black", fill=NA, size=1.5)) +
  theme(text = element_text(size=24), axis.text=element_text(size=24))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(strip.background =  element_blank())+
  theme(axis.ticks.length=unit(0.2, "cm"),  
        axis.text.x = element_text(margin=margin(5,5,5,5,"pt"),colour="black"),
        axis.text.y = element_text(margin=margin(5,5,5,5,"pt"),colour="black"))+ 
  guides(color=guide_legend(override.aes=list(fill=NA)))+
  guides(fill = guide_legend(override.aes = list(stroke = 2,size = 1)))

p2b <- ggplot(fd2mean, aes(x=Competition, y=CWM_sm, group=Rodents, fill=Rodents)) + 
  geom_pointrange(position=position_dodge(width=.5), size=1.5, pch=21, stroke=2,
                  aes(ymin=CWM_sm-(CWM_ssd/sqrt(9)),ymax=CWM_sm+(CWM_ssd/sqrt(9))) )+
  #geom_text(x=.5,y=1.2,label="A", size=8)+
  scale_fill_manual(values=c("black","grey75"), labels=c("No","Yes"))+
  scale_x_discrete(name="Competition", labels=c("No","Yes"))+
  ylab(bquote(~CWM[SeedSize]))+facet_wrap(~Year)+
  #ylab(bquote(~CWM[PCA]))+
  theme(legend.position=c(.85,.2),legend.key=element_blank(),
        legend.background = element_rect(fill="white", size=1, linetype="solid", color="black"))+
  theme(axis.line = element_line(size = 0, colour = "black"))+
  theme(panel.background = element_blank(), plot.margin=margin(.5,.5,.5,.5, "cm"))+
  theme(panel.border = element_rect(color="black", fill=NA, size=1.5)) +
  theme(text = element_text(size=24), axis.text=element_text(size=24))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(strip.background =  element_blank())+
  theme(axis.ticks.length=unit(0.2, "cm"),  
        axis.text.x = element_text(margin=margin(5,5,5,5,"pt"),colour="black"),
        axis.text.y = element_text(margin=margin(5,5,5,5,"pt"),colour="black"))+ 
  guides(color=guide_legend(override.aes=list(fill=NA)))+
  guides(fill = guide_legend(override.aes = list(stroke = 2,size = 1)))

p2c <- ggplot(fd2mean, aes(x=Competition, y=CWM_hgt, group=Rodents, fill=Rodents)) + 
  geom_pointrange(position=position_dodge(width=.5), size=1.5, pch=21, stroke=2,
                  aes(ymin=CWM_hgt-(CWM_hsd/sqrt(9)),ymax=CWM_hgt+(CWM_hsd/sqrt(9))) )+
  #geom_text(x=.5,y=1.2,label="A", size=8)+
  scale_fill_manual(values=c("black","grey75"), labels=c("No","Yes"))+
  scale_x_discrete(name="Competition", labels=c("No","Yes"))+
  #ylab(bquote(~CWM[SeedSize]))+
  ylab(bquote(~CWM[Height]))+ facet_wrap(~Year)+
  theme(legend.position=c(.85,.2),legend.key=element_blank(), 
        legend.background = element_rect(fill="white", size=1, linetype="solid", color="black"))+
  theme(axis.line = element_line(size = 0, colour = "black"))+
  theme(panel.background = element_blank(), plot.margin=margin(.5,.5,.5,.5, "cm"))+
  theme(panel.border = element_rect(color="black", fill=NA, size=1.5)) +
  theme(text = element_text(size=24), axis.text=element_text(size=24))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(strip.background =  element_blank())+
  theme(axis.ticks.length=unit(0.2, "cm"),  
        axis.text.x = element_text(margin=margin(5,5,5,5,"pt"),colour="black"),
        axis.text.y = element_text(margin=margin(5,5,5,5,"pt"),colour="black"))+ 
  guides(color=guide_legend(override.aes=list(fill=NA)))+
  guides(fill = guide_legend(override.aes = list(stroke = 2,size = 1)))

## create .tiff version of figure 2
tiff("CWMPlots.tiff", width=10, height=14, units='in',
     res=600, compression='lzw')
grid.arrange(p2b,p2c,p2a,nrow=3)
dev.off()

#############################################################################################
### code to create Fig 4 -- productivy comp plot
p2_a <- ggplot(fd2 %>% filter(Competition=='+',Rodents=='+'), 
       aes(x=productivity_m2, y=CWM_LMA, shape=Year, linetype=Year)) + 
  #geom_text(x=90,y=.9,label="B", size=8)+
  #stat_smooth(method = "lm", se=F, color="black") + 
  geom_point(size=5, stroke=1.5, color="black") + 
  scale_shape_manual(values=c(1,16), labels=c("Cumulative seedlings","Adults"))+
  scale_linetype_manual(values=c("dashed", "solid"), labels=c("Seedlings","Adults"))+
  ylab(bquote(~CWM[LMA]))+
  xlab(bquote('Productivity (g'~m^-2* 'dry mass)'))+
  theme(legend.position=c(.7,.15),legend.key=element_blank(),legend.title = element_blank(),
        legend.background = element_rect(fill="white", size=1, linetype="solid", color="black"))+
  theme(axis.line = element_line(size = 0, colour = "black"))+
  theme(panel.background = element_blank(), plot.margin=margin(.5,.5,.5,.5, "cm"))+
  theme(panel.border = element_rect(color="black", fill=NA, size=1.5)) +
  theme(text = element_text(size=16), axis.text = element_text(size=14))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(strip.background =  element_blank())+
  theme(axis.ticks.length=unit(0.2, "cm"),  
        axis.text.x = element_text(margin=margin(5,5,5,5,"pt"),colour="black"),
        axis.text.y = element_text(margin=margin(5,5,5,5,"pt"),colour="black"))

p2_b <- ggplot(fd2 %>% filter(Competition=='+',Rodents=='+'), 
               aes(x=productivity_m2, y=CWM_SeedSize, shape=Year, linetype=Year)) + 
  #geom_text(x=90,y=.9,label="B", size=8)+
  stat_smooth(method = "lm", se=F, color="black") + 
  geom_point(size=5, stroke=1.5, color="black") + 
  scale_shape_manual(values=c(1,16), labels=c("Seedlings","Adults"))+
  scale_linetype_manual(values=c("dashed", "solid"), labels=c("Seedlings","Adults"))+
  ylab(bquote(~CWM[SeedSize]))+
  xlab(bquote('Productivity (g'~m^-2* 'dry mass)'))+
  theme(legend.position="none",legend.key=element_blank(),
        legend.background = element_rect(fill="white", size=1, linetype="solid", color="black"))+
  theme(axis.line = element_line(size = 0, colour = "black"))+
  theme(panel.background = element_blank(), plot.margin=margin(.5,.5,.5,.5, "cm"))+
  theme(panel.border = element_rect(color="black", fill=NA, size=1.5)) +
  theme(text = element_text(size=16), axis.text = element_text(size=14))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(strip.background =  element_blank())+
  theme(axis.ticks.length=unit(0.2, "cm"),  
        axis.text.x = element_text(margin=margin(5,5,5,5,"pt"),colour="black"),
        axis.text.y = element_text(margin=margin(5,5,5,5,"pt"),colour="black"))

p2_c <- ggplot(fd2 %>% filter(Competition=='+',Rodents=='+'), 
               aes(x=productivity_m2, y=CWM_Height, shape=Year, linetype=Year)) + 
  #geom_text(x=90,y=.9,label="B", size=8)+
  stat_smooth(method = "lm", se=F, color="black") + 
  geom_point(size=5, stroke=1.5, color="black") + 
  scale_shape_manual(values=c(1,16))+
  scale_linetype_manual(values=c("dashed", "solid"))+
  ylab(bquote(~CWM[Height]))+
  xlab(bquote('Productivity (g'~m^-2* 'dry mass)'))+
  theme(legend.position="none",legend.key=element_blank(),
        legend.background = element_rect(fill="white", size=1, linetype="solid", color="black"))+
  theme(axis.line = element_line(size = 0, colour = "black"))+
  theme(panel.background = element_blank(), plot.margin=margin(.5,.5,.5,.5, "cm"))+
  theme(panel.border = element_rect(color="black", fill=NA, size=1.5)) +
  theme(text = element_text(size=16), axis.text = element_text(size=14))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(strip.background =  element_blank())+
  theme(axis.ticks.length=unit(0.2, "cm"),  
        axis.text.x = element_text(margin=margin(5,5,5,5,"pt"),colour="black"),
        axis.text.y = element_text(margin=margin(5,5,5,5,"pt"),colour="black"))


## create .tiff version of Figure 4
tiff("ProdPlots.tiff", width=15, height=4, units='in',
     res=600, compression='lzw')
grid.arrange(p2_b,p2_c,p2_a,nrow=1)
dev.off()




