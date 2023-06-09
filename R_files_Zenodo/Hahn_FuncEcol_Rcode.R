### Code for Hahn et al. (2021) Intraspecific correlations between growth and defense vary with resource availability and differ within-
##      and among populations. Functional Ecology.
## Contant Philip G Hahn, University of Florida, hahnp@ufl.edu

## load packages
library(lme4)
library(lmerTest)
library(emmeans)
library(car)
library(smatr)
library(MuMIn)
library(patchwork)
library(tidyverse)
library(glmmTMB)
library(DHARMa)


#############################################################################################################################################
#### SITE DATA AND ANALYSIS #################################################################################################################

site1 <- read_csv('Hahn_FuncEcol_SiteData_20210628.csv')

## env PCA
spc <- prcomp(site1[c("orgmat","nit","cec","pH","P","mat1","pwrq18")], scale=TRUE, center=TRUE)
summary(spc)
spc
plot(spc)
biplot(spc, pc.biplot=TRUE)

## extract pc axes
site1$pc1 <-  predict(spc)[,1]
site1$pc2 <-  predict(spc)[,2]

## extract pc components for plotting
pcrot <- as.data.frame(spc$rotation[,1:2], stringsAsFactors=TRUE)
pcrot$Envvar <- c("%OM","soil_N","CEC","pH","soil_P","temp","precip")
pcrot$PC1s <- pcrot$PC1 #scale loadings for plotting as vector arrows
pcrot$PC2s <- pcrot$PC2 #scale loadings for plotting as vector arrows


### PCA FIGURE ###
fig2B <- ggplot() + xlab("PC1 (45.1%)")+ylab("PC2 (27.9%)")+ 
  stat_ellipse(data=site1, geom="polygon" ,aes(x=pc1,y=pc2, color=Region, fill=Region), lwd=1, alpha=.1)+
  geom_point(data=site1 %>% filter(Site!='UWM',Site!='BEN', Site!='MS'),aes(x=pc1, y=pc2, color=Region), size=5,stroke=1.25) +
  geom_point(data=site1 %>% filter(Site=='MS'),aes(x=pc1, y=pc2), pch=22, size=5, color="darkgoldenrod1", fill="darkgoldenrod1") +
  geom_point(data=site1 %>% filter(Site=='BEN'),aes(x=pc1, y=pc2), pch=22, size=5, color="dodgerblue", fill="dodgerblue") +
  geom_point(data=site1 %>% filter(Site=='UWM'),aes(x=pc1, y=pc2), pch=24, size=5, color="dodgerblue", fill="dodgerblue") +
  geom_point(data=site1 %>% filter(Site=='UWM'),aes(x=pc1, y=pc2), pch=25, size=5, color="dodgerblue", fill="dodgerblue") +
  annotate(geom='text', x=-4.5,y=5, label="B", fontface='bold', size=8)+
  scale_fill_manual(values=c("darkgoldenrod1","dodgerblue"), labels=c("Low productivity","High productivity")) +
  scale_color_manual(values=c("darkgoldenrod1","dodgerblue"), labels=c("Low productivity","High productivity")) +
  geom_segment(data=pcrot, aes(x=0, y=0,xend=PC1*5,yend=PC2*5), arrow=arrow(length=unit(0.03, "npc")), color='grey', lwd=1)+
  geom_text(data=pcrot, aes(x=PC1*5.75, y=PC2*5.75, label=Envvar), fontface="bold")+
  theme(axis.line = element_line(size = 1, colour = "black"))+
  theme(panel.background = element_blank())+
  theme(panel.border = element_rect(color="black", fill=NA, size=1.5)) +
  theme(text = element_text(size=18))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.ticks.length=unit(0.2, "cm"),  
        axis.text.x = element_text(margin=margin(5,5,5,5,"pt"),colour="black"),
        axis.text.y = element_text(margin=margin(5,5,5,5,"pt"),colour="black"))+
  theme(plot.title = element_text(hjust = 0.5,size = 10))+
  theme(strip.background = element_rect(colour=NA, fill=NA))+
  theme(strip.text = element_text(size=18))+
  theme(legend.position = "none")+
  guides(lty=F, lty=guide_legend(override.aes = list(lty=NA)))  

## load productivity field data
prod1 <- read_csv("Hahn_WIMT_Productivity_07242018.csv")

prod2 <- full_join(site1,prod1,by="Site")
prod2$veghgt <- (prod2$veghgt1+prod2$veghgt2+prod2$veghgt3+prod2$veghgt4)/4;
prod2$light_reduce <- prod2$light_abv-prod2$light_blw;

head(prod2)
ggplot(prod2, aes(x=veghgt, y=light_reduce))+geom_point()+facet_wrap(~Site)

lgt1 <- lm(light_reduce ~ veghgt, data=prod2)
summary(lgt1)

prodm <- prod2 %>% group_by(Region,Site) %>% summarize(pc1=mean(pc1), pc2=mean(pc2), veghgt=mean(veghgt, na.rm=T), mat1=mean(mat1,na.rm=T), pwrq18=mean(pwrq18,na.rm=T),light_reduce=(mean(light_reduce,na.rm=T)))

## test relationship between veg height and pc1
veg1 <- lm(log(veghgt)~pc1, data=prodm)
anova(veg1)
v1 <- summary(veg1)

veg2 <- lm(log(veghgt)~pc2, data=prodm)
v2 <- summary(veg2)

## test differences in PC axes between regions
an1 <- lm(pc1~Region, data=site1)
anova(an1)
summary(an1)
an2 <- lm(pc2~Region, data=site1)
anova(an2)
summary(an2)

pc1x <- seq(min(prodm$pc1),max(prodm$pc1),.1)
vegy <- exp(predict(veg1,list(pc1=pc1x)))

fig2C <- ggplot() + xlab("PC1")+ylab("mean veg height (cm)")+ geom_line(data=, aes(x=pc1x, y=vegy), lwd=1.5) +
  geom_point(data=prodm %>% filter(Site!='MS',Site!='BEN'), aes(x=pc1, y=veghgt, color=Region), size=5)+
  geom_point(data=prodm %>% filter(Site=='MS'), aes(x=pc1, y=veghgt, color=Region), size=5, pch=22, color="darkgoldenrod1", fill="darkgoldenrod1")+
  geom_point(data=prodm %>% filter(Site=='BEN'), aes(x=pc1, y=veghgt, color=Region), size=5, pch=22, color="dodgerblue", fill="dodgerblue")+
  geom_point(data=prodm %>% filter(Site=='UWM'), aes(x=pc1, y=veghgt, color=Region), size=5, pch=24, color="dodgerblue", fill="dodgerblue")+
  geom_point(data=prodm %>% filter(Site=='UWM'), aes(x=pc1, y=veghgt, color=Region), size=5, pch=25, color="dodgerblue", fill="dodgerblue")+
  scale_color_manual(values=c("darkgoldenrod1","dodgerblue"), labels=c("Low resource","High resource")) +
  annotate(geom='text', x=-3,y=61, label="C", fontface='bold', size=8)+
  #annotate(geom='text', x=-1.5,y=50, label="italic(R)^2 == 0.84", parse=T, size=6)+
  theme(axis.line = element_line(size = 1, colour = "black"))+
  theme(panel.background = element_blank())+
  theme(panel.border = element_rect(color="black", fill=NA, size=1.5)) +
  theme(text = element_text(size=18))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.ticks.length=unit(0.2, "cm"),  
        axis.text.x = element_text(margin=margin(5,5,5,5,"pt"),colour="black"),
        axis.text.y = element_text(margin=margin(5,5,5,5,"pt"),colour="black"))+
  theme(plot.title = element_text(hjust = 0.5,size = 10))+
  theme(strip.background = element_rect(colour=NA, fill=NA))+
  theme(strip.text = element_text(size=18))+
  theme(legend.position = "top", legend.background = element_rect(fill="white", size=1.5, linetype="solid", color="black"))

tiff("Hahn_FuncEcol_Fig2BC.tiff", width=10, height=5, units='in',
     res=300, compression='lzw')
(fig2B+guides(color='none',fill='none'))+fig2C+plot_layout(guides="collect") & theme(legend.position = 'top')
dev.off()




#################################################################################################################################################
######### DATA AND ANALYSIS FOR ADDRESSING MAIN RESEARCH QUESTIONS #############################################################################
#### load data
s51 <- read_csv("Hahn_FuncEcol_FullTraits.csv")

## summarize data by Population
s5p <- s51 %>% group_by(Region,Site) %>% summarize(Biomass_g=mean(Biomass_g, na.rm=T),total=mean(total, na.rm=T),
                                                   thymol=mean(thymol, na.rm=T), carvacrol=mean(carvacrol, na.rm=T), totall=mean(totall,na.rm=T))



### SMA ANALYSIS - Q1 examining tradeoffs across populations from two regions ###################################################################
sma_pop1 <- sma(totall~Biomass_g, data=s5p)
summary(sma_pop1)
plot(sma_pop1)
plot(sma_pop1, which="residual")
plot(sma_pop1, which="qq")

## Q2 examining tradeoffs across populations within each region
sma_pop2 <- sma(totall~Biomass_g*Region, data=s5p, robust=T)
summary(sma_pop2)
plot(sma_pop2)
plot(sma_pop2, which="residual")
plot(sma_pop2, which="qq")

smap1 <- s5p
smap1$y1 <- (coef(sma_pop1)[1]+smap1$Biomass_g*coef(sma_pop1)[2]) 

smap1line <- with(s5p, line.cis(x=Biomass_g, y=totall, alpha=.66))

smap1fit <- fitted(sma_pop1, type="fitted", centered=F, interval="confidence")

smap1$y1L <- (smap1line[1,2]+smap1$Biomass_g*smap1line[2,2]) 
smap1$y1U <- (smap1line[1,3]+smap1$Biomass_g*smap1line[2,3]) 
smap1$yMT <- (coef(sma_pop2)[1,1]+smap1$Biomass_g*coef(sma_pop2)[1,2]) 
smap1$yWI <- (coef(sma_pop2)[2,1]+smap1$Biomass_g*coef(sma_pop2)[2,2]) 

p1b <- ggplot(smap1, aes(x=Biomass_g,y=exptotall)) + 
  #geom_ribbon(data=smap1, aes(ymin=y1L,ymax=y1U,x=Biomass_g), fill='grey')+
  #geom_smooth(data=s5p, aes(x=Biomass_g,y=totall), method="lm", se=T, color="black", size=1.5)+
  geom_line(data=smap1, aes(x=Biomass_g,y=y1), color="black", size=2)+
  #geom_line(data=smap1, aes(x=Biomass_g,y=y1L), color="black", size=2)+
  #geom_line(data=smap1, aes(x=Biomass_g,y=y1U), color="black", size=2)+
  geom_line(data=smap1 %>% filter(Region=='MT'), aes(x=Biomass_g,y=yMT), color="goldenrod", size=1.25, lty=3)+
  geom_line(data=smap1 %>% filter(Region=='WI'), aes(x=Biomass_g,y=yWI), color="dodgerblue", size=1.25, lty=1)+
  geom_point(data=smap1, aes(x=Biomass_g,y=totall,color=Region), size=3, pch=16) +
  geom_point(data=smap1 %>% filter(Site=='MS'), aes(x=Biomass_g,y=totall), size=3.25, color="goldenrod", pch=15) +
  geom_point(data=smap1 %>% filter(Site=='BEN'), aes(x=Biomass_g,y=totall), size=3.25, color="dodgerblue", pch=15) +
  annotate(geom='text', x=.3,y=4.05, label="A", fontface='bold', size=5)+
  scale_color_manual(values=c("goldenrod","dodgerblue"), labels=c("Low productivity","High productivity"))+
  scale_x_continuous("Biomass (mg)", limits = c(0.3,3.25), breaks=seq(.5,3,.5))+
  scale_y_continuous("Total terpenes log(mg/g)", limits = c(1.5,4.15), breaks=seq(1.5,4,.5))+
  #coord_cartesian(xlim = c(0.3,2.55), ylim = c(1.5,3.8))+
  theme(legend.key=element_blank(), legend.position = c(.7,.8), 
        legend.background = element_rect(fill="white", size=1, linetype="solid", color="black"))+  
  theme(panel.background = element_blank())+
  theme(panel.border = element_rect(color="black", fill=NA, size=2)) +
  theme(text = element_text(size=15))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.ticks.length=unit(0.3, "cm"),  
        axis.text.x = element_text(margin=margin(15,5,3,5,"pt"),colour="black"),
        axis.text.y = element_text(margin=margin(3,15,10,10,"pt"),colour="black"))

smap1 <- s5p
smap1$y1 <- exp(coef(sma_pop1)[1]+smap1$Biomass_g*coef(sma_pop1)[2]) 
smap1$yMT <- exp(coef(sma_pop2)[1,1]+smap1$Biomass_g*coef(sma_pop2)[1,2])
smap1$yWI <- exp(coef(sma_pop2)[2,1]+smap1$Biomass_g*coef(sma_pop2)[2,2]) 

p1b <- ggplot(smap1, aes(x=Biomass_g,y=total)) + 
  #geom_ribbon(data=smap1, aes(ymin=y1L,ymax=y1U,x=Biomass_g), fill='grey')+
  #geom_smooth(data=s5p, aes(x=Biomass_g,y=totall), method="lm", se=T, color="black", size=1.5)+
  geom_line(data=smap1, aes(x=Biomass_g,y=y1), color="black", size=2)+
  #geom_line(data=smap1, aes(x=Biomass_g,y=y1L), color="black", size=2)+
  #geom_line(data=smap1, aes(x=Biomass_g,y=y1U), color="black", size=2)+
  geom_line(data=smap1 %>% filter(Region=='MT'), aes(x=Biomass_g,y=yMT), color="goldenrod", size=1.25, lty=3)+
  geom_line(data=smap1 %>% filter(Region=='WI'), aes(x=Biomass_g,y=yWI), color="dodgerblue", size=1.25, lty=1)+
  geom_point(data=smap1, aes(x=Biomass_g,y=exp(totall),color=Region), size=3, pch=16) +
  geom_point(data=smap1 %>% filter(Site=='MS'), aes(x=Biomass_g,y=exp(totall)), size=3.25, color="goldenrod", pch=15) +
  geom_point(data=smap1 %>% filter(Site=='BEN'), aes(x=Biomass_g,y=exp(totall)), size=3.25, color="dodgerblue", pch=15) +
  annotate(geom='text', x=1,y=48, label="A) among populations", fontface='bold', size=5)+
  scale_color_manual(values=c("goldenrod","dodgerblue"), labels=c("Low resource","High resource"))+
  scale_x_continuous("Biomass (mg)", limits = c(0.3,2.6), breaks=seq(.5,2.5,.5))+
  scale_y_continuous("Total terpenes (mg/g)", limits = c(5,49), breaks=seq(0,40,10))+
  #coord_cartesian(xlim = c(0.3,2.55), ylim = c(1.5,3.8))+
  theme(legend.key=element_blank(), legend.position = c(.75,.6), 
        legend.background = element_rect(fill="white", size=1, linetype="solid", color="black"))+  
  theme(panel.background = element_blank())+
  theme(panel.border = element_rect(color="black", fill=NA, size=2)) +
  theme(text = element_text(size=15))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.ticks.length=unit(0.3, "cm"),  
        axis.text.x = element_text(margin=margin(15,5,3,5,"pt"),colour="black"),
        axis.text.y = element_text(margin=margin(3,15,10,10,"pt"),colour="black"))

######################################################################################################################################
##### SMAR for Q3 - tradeoffs among siblings across all 12 populations
pop1 <- s51 %>% filter(Site=='BEN'|Site=='MS') ## filter to populations with higher maternal plant-level replication

## summarize data by each maternal plant
pop1a <- pop1 %>% group_by(Region,Site,Plant) %>% summarize(Biomass_g=mean(Biomass_g, na.rm=T), total=mean(total, na.rm=T),
                                                            thymol=mean(thymol, na.rm=T), carvacrol=mean(carvacrol, na.rm=T), totall=mean(totall,na.rm=T))

## build model for Q3 - tradeoffs among siblings across all 12 populations
sma1 <- sma(totall~Biomass_g, data=pop1a )
summary(sma1)
plot(sma1)
plot(sma1, which="residual")
plot(sma1, which="qq")


sma2 <- sma(totall~Biomass_g*Region, data=pop1a, type="shift")
sma2
summary(sma2)
plot(sma2)
plot(sma2, which="residual")
plot(sma2, which="qq")

## extract model coefficients to make figure
smap2 <- pop1a
smap2$y1 <- exp(coef(sma1)[1]+smap2$Biomass_g*coef(sma1)[2]) 
smap2$yMT <- exp(coef(sma2)[1,1]+smap2$Biomass_g*coef(sma2)[1,2]) 
smap2$yWI <- exp(coef(sma2)[2,1]+smap2$Biomass_g*coef(sma2)[2,2]) 

## BIG POPS FIG
p2b <- ggplot(pop1a, aes(x=Biomass_g,y=total)) + 
  #geom_smooth(data=smap2, aes(x=Biomass_g,y=totall), method="lm", color="black", size=1.5, se=T)+
  geom_line(data=smap2 , aes(x=Biomass_g,y=y1), color="black", size=1.5)+
  geom_line(data=smap2 %>% filter(Region=='MT'), aes(x=Biomass_g,y=yMT), color="goldenrod", size=1.25, lty=2)+
  geom_line(data=smap2 %>% filter(Region=='WI'), aes(x=Biomass_g,y=yWI), color="dodgerblue", size=1.25, lty=2)+
  geom_point(data=pop1a, aes(x=Biomass_g,y=total,color=Region), size=2, pch=15) +
  annotate(geom='text', x=.95,y=48, label="B) within populations", fontface='bold', size=5)+
  scale_color_manual(values=c("goldenrod","dodgerblue"))+
  scale_x_continuous("Biomass (mg)", limits = c(0.3,2.6), breaks=seq(.5,2.5,.5))+
  scale_y_continuous("Total terpenes (mg/g)", limits = c(5,49), breaks=seq(0,43,10))+
  theme(legend.key=element_blank(), legend.position = 'none')+  
  theme(panel.background = element_blank())+
  theme(panel.border = element_rect(color="black", fill=NA, size=2)) +
  theme(text = element_text(size=15))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.ticks.length=unit(0.3, "cm"),  
        axis.text.x = element_text(margin=margin(15,5,3,5,"pt"),colour="black"),
        axis.text.y = element_text(margin=margin(3,15,10,10,"pt"),colour="black"))

pdf("FE-2021-00135_Fig3.pdf", width=5, height=8)
p1b/p2b
dev.off()


##############################################################################################################################
################# FIELD DATA #################################################################################################

dmg2 <- read_csv('Hahn_FuncEcol_FieldData_20210628.csv')

## smithson and Verkuilen 2006 transformation
dmg2$LeafDmg <- (dmg2$LfDmg*(length(dmg2$LfDmg)-1)+.5)/length(dmg2$LfDmg)
dmg2$SeedDmg <- (dmg2$PropDmg*(length(dmg2$PropDmg)-1)+.5)/length(dmg2$PropDmg)

## model for plant stem height
lmd3 <- lmer(Hgt~Region+(1|Site), data=dmg2)
anova(lmd3, ddf="Kenward-Roger")
lmd3m <- emmeans(lmd3, ~Region, transform="response") %>% as.data.frame()
summary(lmd3)
r.squaredGLMM(lmd3)

sim_lmd3 <- simulateResiduals(fittedModel = lmd3, n = 250)
plot(sim_lmd3) ## slight difference in variance between groups, but not major. Doesn't affect results.
plot(lmd3)
hist(resid(lmd3))

## model for leaf damage
lmd1 <- lmer(logit(LeafDmg)~Region+(1|Site), data=dmg2)
summary(lmd1)
anova(lmd1, ddf="Kenward-Roger")
r.squaredGLMM(lmd1)
emmeans(lmd1, ~1, transform="response")
lmd1m <- emmeans(lmd1, ~Region, transform="response") %>% as.data.frame()

sim_lmd1 <- simulateResiduals(fittedModel = lmd1, n = 250)
plot(sim_lmd1)
plot(lmd1)
hist(resid(lmd1))

## model for seed head damage
lmd2 <- glmer(cbind(NumbInsect,(Totheads-NumbInsect))~Region+(1|Site), data=dmg2, family='binomial')
Anova(lmd2, type=2)
summary(lmd2)
r.squaredGLMM(lmd2)
lmd2m <- emmeans(lmd2, ~Region, transform="response") %>% as.data.frame()

sim_lmd2 <- simulateResiduals(fittedModel = lmd2, n = 250)
plot(sim_lmd2)


## summarize data for plotting population means
dmg3 <- dmg2 %>% dplyr::select(Region,Site,Hgt,LfDmg,PropDmg,LeafDmg,SeedDmg) %>% group_by(Region,Site) %>% 
  summarize(Hgt=mean(Hgt, na.rm=T), LfDmg=mean(LfDmg, na.rm=T), PropDmg=mean(PropDmg, na.rm=T), 
            LeafDmg=mean(LeafDmg, na.rm=T), SeedDmg=mean(SeedDmg, na.rm=T))

set.seed(17)
phgt <- ggplot()+geom_jitter(data=dmg2, aes(x=Region, y=Hgt), height=0,width=.2, color='grey')+
  geom_jitter(data=dmg3, aes(x=Region, y=Hgt, color=Region), height=0,width=.3, size=3)+
  geom_point(data=lmd3m, aes(x=Region, y=emmean, color=Region), size=6, shape='square')+
  geom_errorbar(data=lmd3m, aes(x=Region, y=emmean, ymin=lower.CL, ymax=upper.CL, color=Region), width=0, lwd=2)+
  scale_color_manual(values=c("goldenrod","dodgerblue"), labels=c("Low resource", "High resource"), 
                     guide=guide_legend(direction="vertical"))+
  scale_y_continuous("Stem height (cm)")+
  theme(legend.key=element_blank(), legend.position = c(.35,.8), 
        legend.background = element_rect(fill=NA, size=1, linetype="solid", color="black"))+  
  theme(panel.background = element_blank())+
  theme(panel.border = element_rect(color="black", fill=NA, size=2)) +
  theme(text = element_text(size=18))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.ticks.length=unit(0.3, "cm"),  
        axis.text.x = element_text(margin=margin(15,5,3,5,"pt"),colour="black"),
        axis.text.y = element_text(margin=margin(3,15,10,10,"pt"),colour="black"))+
  annotate(geom='text', x=.5,y=105, label="A", fontface='bold', size=5)+
  annotate(geom='text', x=1.5,y=75, label="p=0.004", size=4)

set.seed(3)
pleaf <- ggplot()+geom_jitter(data=dmg2, aes(x=Region, y=LeafDmg), height=0,width=.2, color='grey')+
  geom_jitter(data=dmg3, aes(x=Region, y=LeafDmg, color=Region), height=0,width=.25, size=3)+
  geom_point(data=lmd1m, aes(x=Region, y=response, color=Region), size=6, shape='square')+
  geom_errorbar(data=lmd1m, aes(x=Region, y=response, ymin=lower.CL, ymax=upper.CL, color=Region), width=0, lwd=2)+
  scale_color_manual(values=c("goldenrod","dodgerblue"))+
  scale_y_continuous("Leaf damage", labels=scales::percent_format(accuracy=1), breaks=c(0,.05,.1,.15,.2,.25))+
  theme(legend.key=element_blank(), legend.position = "none")+  
  theme(panel.background = element_blank())+
  theme(panel.border = element_rect(color="black", fill=NA, size=2)) +
  theme(text = element_text(size=18))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.ticks.length=unit(0.3, "cm"),  
        axis.text.x = element_text(margin=margin(15,5,3,5,"pt"),colour="black"),
        axis.text.y = element_text(margin=margin(3,15,10,10,"pt"),colour="black"))+
  annotate(geom='text', x=.5,y=.26, label="B", fontface='bold', size=5)+
  annotate(geom='text', x=1.5,y=.2, label="p=0.27", size=4)
  
set.seed(1)
pseed <- ggplot()+geom_jitter(data=dmg2, aes(x=Region, y=SeedDmg), height=0,width=.2, color='grey')+
  geom_jitter(data=dmg3, aes(x=Region, y=SeedDmg, color=Region), height=0,width=.25, size=3)+
  geom_point(data=lmd2m, aes(x=Region, y=response, color=Region), size=6, shape='square')+
  geom_errorbar(data=lmd2m, aes(x=Region, y=response, ymin=asymp.LCL, ymax=asymp.UCL, color=Region), width=0, lwd=2)+
  scale_color_manual(values=c("goldenrod","dodgerblue"))+
  scale_y_continuous("Sead head damage", labels=scales::percent_format(accuracy=1))+
  theme(legend.key=element_blank(), legend.position = 'none')+  
  theme(panel.background = element_blank())+
  theme(panel.border = element_rect(color="black", fill=NA, size=2)) +
  theme(text = element_text(size=18))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.ticks.length=unit(0.3, "cm"),  
        axis.text.x = element_text(margin=margin(15,5,3,5,"pt"),colour="black"),
        axis.text.y = element_text(margin=margin(3,15,10,10,"pt"),colour="black"))+
  annotate(geom='text', x=.5,y=1, label="C", fontface='bold', size=5)+
  annotate(geom='text', x=1.5,y=.8, label="p<0.001", size=4)


pdf("FE-2021-00135_Fig4.pdf", width=5.5, height=12)
phgt/pleaf/pseed
dev.off()




