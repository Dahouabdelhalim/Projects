#open packages
library(dplyr)
library(ggplot2)
library(car)
library(lme4)
library(extrafont)
loadfonts(device="win")
#this "loadfonts()" function is not necessary on a mac

#clear workspace
rm(list=ls())

#import the beetle data
beetle_hole_data <- read.csv("beetle_hole_data.csv")
str(beetle_hole_data)

#aggregate data
aggregated_data<-beetle_hole_data%>%
  group_by(control_strike,strike,treetag)%>%
  summarise(Emergence.Holes.sum=sum(Emergence.Holes),
            Dieback.mean=mean(Dieback),
            distance.from.center.mean = mean(distance.from.center),
            diameter.mean = mean(DBH.mm))
str(aggregated_data)

#now create a variable for paired analyses (grouping paired strikes together in terms of random variation)
aggregated_data$strike_paired<-as.factor(gsub("N","",aggregated_data$strike))
aggregated_data$strike_paired

##################################################
#confirm that strike and control trees follow the same size distribution
ks.test(aggregated_data$diameter.mean[aggregated_data$control_strike=="Strike"],aggregated_data$diameter.mean[aggregated_data$control_strike=="Control"])

##################################################
#make figures

#create custom theme style
theme_basis<-theme(axis.title = element_text(family = "Arial", color="black", face="bold", size=14))+
  theme(axis.text.y=element_text(family = "Arial", 
                                 colour="black", face ="bold", size = 12)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),axis.ticks = element_line(colour="black"))+
  theme(axis.text.x = element_text(family = "Arial", colour="black", face ="bold", size = 12))+
  theme(plot.margin = unit(c(.5, .5, .5, .5), "cm"))

#create a figure with density plots for both groups of trees
density_plot<-ggplot(aggregated_data, aes(x = diameter.mean, linetype = control_strike, alpha = control_strike))+
  geom_density(fill = "black")+ 
  scale_alpha_manual(values = c(0.2,0.4))+
  scale_linetype_manual(values=c("dashed","solid"))+
  scale_y_continuous(name = "Probability Density", expand = c(0,0.00001),breaks = seq(0,0.003,0.0005), 
                     labels = c("0.000","","0.001","","0.002","","0.003"))+
  scale_x_continuous(name = "Diameter at Breast Height (mm)",breaks = c(seq(100,1000,100),1500), 
                labels = c("100", rep("",8),"1000","1500"), expand = c(.002,0))+
  theme(legend.position = c(.8,.8), legend.title = element_blank(), legend.background = element_blank(),
        legend.text = element_text(family = "Arial",colour="black", face ="bold",size=11))+
  theme(axis.text = element_text(color = "black"))+theme_basis

density_plot

ggsave("density_plot.png",density_plot,width = 3, height = 2,dpi = 600, scale = 1.9)
ggsave("density_plot.tiff",density_plot,width = 3, height = 2,dpi = 600, scale = 1.9, compression = "lzw")


#Diameter vs Holes Figure

diameterfigure<-ggplot(aggregated_data, aes(x = diameter.mean,y = Emergence.Holes.sum+1))+
  geom_point(size = 2,aes(shape = aggregated_data$control_strike))+
  scale_y_log10(name = "Emergence Holes (sum)",expand = c(0.02,0),
                breaks = c(seq(1,11,1),seq(21,101,10),seq(201,1001,100)),
                labels = c("0",rep("",9),"10",rep("",8),"100",rep("",8),"1000"))+
  scale_x_log10(name = "Diameter at Breast Height (mm)",breaks = c(seq(100,1000,100),1500), 
                labels = c("100", rep("",8),"1000","1500"), expand = c(.01,0))+
  geom_smooth(method = "lm",se=TRUE,
              aes(linetype = aggregated_data$control_strike),color = "black")+
  scale_linetype_manual(values = c("dashed","solid"))+
  scale_shape_manual(values = c(1,16))+
  theme(legend.position = c(.8,.13), legend.title = element_blank(), legend.background = element_blank(),
                        legend.text = element_text(family = "Arial",colour="black", face ="bold",size=11))+
  theme(axis.text = element_text(color = "black"))+theme_basis

diameterfigure
ggsave("Holes_Diameter.png",diameterfigure,dpi=600)
ggsave("Holes_Diameter.tiff",diameterfigure,dpi=600,compression = "lzw")
getwd()

#Distance vs Emergence Holes
distancefigure<-ggplot(aggregated_data, aes(x = distance.from.center.mean,y = Emergence.Holes.sum+1))+
  geom_point(size = 2,aes(shape = aggregated_data$control_strike))+
  scale_y_log10(name = "Emergence Holes (sum)",expand = c(0.02,0),
                breaks = c(seq(1,11,1),seq(21,101,10),seq(201,1001,100)),
                labels = c("0",rep("",9),"10",rep("",8),"100",rep("",8),"1000"))+
  scale_x_continuous(name = "Distance from Central Tree (m)",breaks = c(seq(0,40,5)), 
                     labels = c("0","","10","","20","","30","","40"),expand = c(.01,0))+
  geom_smooth(method = "lm",se=TRUE,
              aes(linetype = aggregated_data$control_strike),color = "black")+
  scale_linetype_manual(values = c("dashed","solid"))+
  scale_shape_manual(values = c(1,16))+
  theme(legend.position = c(.8,.13), legend.title = element_blank(), legend.background = element_blank(),
        legend.text = element_text(family = "Arial",colour="black", face ="bold",size=11))+
  theme(axis.text = element_text(color = "black"))+theme_basis

distancefigure
ggsave("Holes_Distance.png",distancefigure,dpi=600)
ggsave("Holes_Distance.tiff",distancefigure,dpi=600,compression = "lzw")
getwd()

###########################################

control_only<-aggregated_data[aggregated_data$control_strike=="Control",]
str(control_only)

colnames(control_only)[7]<-"Diameter (mm)"

#Distance with only control trees - plot the size of points based on tree DBH
distancefigure_control<-ggplot(control_only, 
                       aes(x = distance.from.center.mean,y = Emergence.Holes.sum+1))+
  geom_point(aes(size = `Diameter (mm)`),shape = 1)+
  scale_y_log10(name = "Emergence Holes (sum)",expand = c(0.02,0),
                breaks = c(seq(1,11,1),seq(21,101,10),seq(201,1001,100)),
                labels = c("0",rep("",9),"10",rep("",8),"100",rep("",8),"1000"))+
  scale_x_continuous(name = "Distance from Central Tree (m)",breaks = c(seq(0,40,5)), 
                     labels = c("0","","10","","20","","30","","40"),expand = c(.05,0))+
  geom_smooth(method = "lm",se=TRUE, color = "black", linetype = "dashed")+
  theme(legend.position = c(.85,.88), legend.background = element_blank(),legend.key=element_blank(),
        legend.text = element_text(family = "Arial",colour="black", face ="bold",size=11),
        legend.title = element_text(family = "Arial",colour="black", face ="bold",size=12),
        legend.key.height = unit(0.3,"cm"))+
  theme(axis.text = element_text(color = "black"))+theme_basis

distancefigure_control

ggsave("Holes_Distance_control.tiff",distancefigure_control,dpi=600, width = 5, height = 4,compression = "lzw")
getwd()

reg1_wfocal <- glmer(Emergence.Holes.sum~distance.from.center.mean + (1|strike), data = control_only, family = poisson)
reg2_wfocal <- glmer(Emergence.Holes.sum~1 + (1|strike), data = control_only, family = poisson)
summary(reg1_wfocal)
anova(reg1_wfocal,reg2_wfocal)

control_only_nofocal<-control_only[control_only$distance.from.center.mean!=0,]

reg1 <- glmer(Emergence.Holes.sum~distance.from.center.mean + (1|strike), data = control_only_nofocal, family = poisson)
reg2 <- glmer(Emergence.Holes.sum~1 + (1|strike), data = control_only_nofocal, family = poisson)
anova(reg1,reg2)


###########################################

#Dieback vs Emergence Holes
logit_dieback<-logit(aggregated_data$Dieback.mean)
str(logit_dieback)
breaks_x_logit <- logit(c(0,5,seq(10,100,10)))

Diebackfigure<-ggplot(aggregated_data, aes(x = logit_dieback,y = Emergence.Holes.sum+1))+
  geom_point(size = 2,aes(shape = aggregated_data$control_strike),
             position = position_jitter(width = 0.2))+
  scale_y_log10(name = "Emergence Holes (sum)",expand = c(0.02,0),
                breaks = c(seq(1,11,1),seq(21,101,10),seq(201,1001,100)),
                labels = c("0",rep("",9),"10",rep("",8),"100",rep("",8),"1000"))+
  scale_x_continuous(name = "Dieback (%)",breaks = breaks_x_logit,expand = c(.01,0),
                     labels = c("0","<5","10","20","30","40","50","60","70","80","90","100"))+
  geom_smooth(method = "lm",se=TRUE,
              aes(linetype = aggregated_data$control_strike),color = "black")+
  scale_linetype_manual(values = c("dashed","solid"))+
  scale_shape_manual(values = c(1,16))+
  theme(legend.position = c(.8,.13), legend.title = element_blank(), legend.background = element_blank(),
        legend.text = element_text(family = "Arial",colour="black", face ="bold",size=11))+
  theme(axis.text = element_text(color = "black"))+theme_basis

Diebackfigure

ggsave("Holes_Dieback.png",Diebackfigure,dpi=600)
ggsave("Holes_Dieback.tiff",Diebackfigure,dpi=600, compression = "lzw")


####################################
#test whether emergence holes are non-randomly associated with damaged crown sections
strike_only<-beetle_hole_data[beetle_hole_data$control_strike=="Strike",]#include only struck trees
strike_only$crown.dmg.fac<-as.factor(strike_only$Crown.Damage)

#create poisson glm to test whether being associated with damage increases strike likelihood
####this model controls for tree effects and accounts for total dieback
#make sure that this effect is not due to the effects of completlely dead trees
####to do this, remove dead trees
strike_only_alive<-strike_only[strike_only$Dieback<98,]
nrow(strike_only)-nrow(strike_only_alive)#172 trees removed

#create a plot demonstrating differences in beetle counts between crown damage/health across range of dieback values
#I will use both dead and living trees because the model works both ways
test_distributions<-strike_only%>%
  group_by(strike_tag)%>%
  summarise(total_damaged=sum(Crown.Damage),
            dieback = mean(Dieback),
            total_sections = length(Crown.Damage))

hist(test_distributions$total_damaged)
hist(test_distributions$dieback)
hist(test_distributions$total_sections)

#grab only the trees with total damaged sections >1 or <8 because these have no variation
damaged_only_selection<-test_distributions[test_distributions$total_damaged>1&test_distributions$total_damaged<8,]
nrow(damaged_only_selection)-nrow(test_distributions)#97 rows were removed
nrow(damaged_only_selection)#76 trees remain

#left join to grab the data that match these trees
damaged_only<-inner_join(damaged_only_selection,strike_only, by = "strike_tag")
str(damaged_only)#here are 608 rows - equaling 8 sections of 76 different trees

#create poisson glm to test whether being associated with damage increases strike likelihood
####this model controls for tree effects and accounts for total dieback
#now runt he model again with only these trees
mod_crown_damaged<-glmer(Emergence.Holes ~ crown.dmg.fac + Dieback + (1|strike_tag),data = damaged_only,
                       family = poisson)
plot(mod_crown_damaged)#great plot
qqnorm(resid(mod_crown_damaged));qqline(resid(mod_crown_damaged))#not perfect, but it is ok
summary(mod_crown_damaged)

mod_crown2_damaged<-glmer(Emergence.Holes ~ Dieback + (1|strike_tag),data = damaged_only,
                        family = poisson)
#the model fit is decent (not great), but barely fails to pass convergence test by 0.00000365
#lets check for convergence issues
tt <- getME(mod_crown2_damaged,"theta")
ll <- getME(mod_crown2_damaged,"lower")
min(tt[ll==0])
#this value is very far from zero, so we are all good

#double check the gradient calculation
derivs1 <- mod_crown2_damaged@optinfo$derivs
sc_grad1 <- with(derivs1,solve(Hessian,gradient))
max(abs(sc_grad1))
#calculating it by hand suggests that we actually are not having issues - the quick-and-dirty internal check was just flawed
#this lack-of-an-error is supported by the great model fit shown above

max(pmin(abs(sc_grad1),abs(derivs1$gradient)))

#we can then redo the caluculations with numerdiv
library(numDeriv)
dd <- update(mod_crown2_damaged,devFunOnly=TRUE)
pars <- unlist(getME(mod_crown2_damaged,c("theta","fixef")))
grad2 <- grad(dd,pars)
hess2 <- hessian(dd,pars)
sc_grad2 <- solve(hess2,grad2)
max(pmin(abs(sc_grad2),abs(grad2)))
#this also shows that the gradient is far below 0.001 (the threshold defined in lme4)

#try rerunning the model with these alterations
ss <- getME(mod_crown2_damaged,c("theta","fixef"))
mod_crown2_damaged_updated <- update(mod_crown2_damaged,start=ss,control=glmerControl(optCtrl=list(maxfun=2e4)))
#this also says it does not converge (which still seems overly conservative because the gradient is less than 0.002)
mod_crown2_damaged_bobyqa <- update(mod_crown2_damaged,start=ss,control=glmerControl(optimizer="bobyqa",
                                                 optCtrl=list(maxfun=2e5)))
#this doesn't give a warning (and the other version actually seemed to give a similarly close warning - so continue)
anova(mod_crown_damaged,mod_crown2_damaged_bobyqa,test = "LRT")

mod_crown3_damaged<-glmer(Emergence.Holes ~ crown.dmg.fac + (1|strike_tag),data = damaged_only,
                          family = poisson)
anova(mod_crown_damaged,mod_crown3_damaged,test = "LRT")
summary(mod_crown3_damaged)

#no effect of dieback!!
#but there is a strong effect of being directly under a damaged section of the crown
#this suggests that lightning may not follow wood grain in tropical forests, or just indicate a lower frequency of wood grain twist here

#calculate mean within each group of crown damage
mean(damaged_only$Emergence.Holes[damaged_only$crown.dmg.fac==1])
mean(damaged_only$Emergence.Holes[damaged_only$crown.dmg.fac==0])
damaged_only%>%
  group_by(crown.dmg.fac)%>%
  summarise(sd.holes = sd(Emergence.Holes),
            n.holes = length(Emergence.Holes))%>%
  mutate(se.holes = sd.holes/(n.holes^(1/2)))

#########################################################
#analyze per-tree beetle counts

#a generalized linear model is probably our best option
#rescale variables so that the generalized model will converge
aggregated_data$logit_dieback<-logit(aggregated_data$Dieback.mean)#logit transformation
aggregated_data$scaled_distance<-scale(aggregated_data$distance.from.center.mean)#rescaled (z-transformed)
aggregated_data$scaled_diameter<-scale(aggregated_data$diameter.mean)#log

#test for covariance
cor.test(aggregated_data$logit_dieback,aggregated_data$scaled_distance)
cor.test(aggregated_data$logit_dieback,aggregated_data$scaled_diameter)
cor.test(aggregated_data$scaled_distance,aggregated_data$scaled_diameter)

#calculate mean within each group of crown damage
mean(aggregated_data$Emergence.Holes.sum[aggregated_data$control_strike=="Strike"])
mean(aggregated_data$Emergence.Holes.sum[aggregated_data$control_strike=="Control"])
aggregated_data%>%
  group_by(control_strike)%>%
  summarise(sd.holes = sd(Emergence.Holes.sum),
            n.holes = length(Emergence.Holes.sum))%>%
  mutate(se.holes = sd.holes/(n.holes^(1/2)))

#now run the glm
#we will only include the individual contextual terms and their interactions with the control_strike (control versus strike)
gmod1<-glmer(Emergence.Holes.sum~control_strike + logit_dieback + scaled_distance +
             scaled_diameter + control_strike:scaled_diameter + control_strike:logit_dieback +
             control_strike:scaled_distance + (1|strike_paired), data = aggregated_data, 
             family = poisson, control = glmerControl(optCtrl=list(maxfun=2e4)))
plot(gmod1)
qqnorm(resid(gmod1));qqline(resid(gmod1))#there is one outlier in the residuals - removing this outlier makes the dieback:control/strike interaction non-significant, but that does not meaningfully change our interpretation of the results
########we will not remove the outlier because there is no biological reason to remove the outlier
summary(gmod1)
AIC(gmod1)
######THIS IS YOUR BEST FIT MODEL

#confirm that the random effect matters
gmod1_norand<-glm(Emergence.Holes.sum~control_strike + logit_dieback + scaled_distance +
               scaled_diameter + control_strike:scaled_diameter + control_strike:logit_dieback +
               control_strike:scaled_distance, data = aggregated_data, 
             family = poisson)
summary(gmod1_norand)
anova(gmod1,gmod1_norand)
#the random effect has a strong effect

#drop control_strike versus dieback and test for an effect
gmod2<-glmer(Emergence.Holes.sum~control_strike + logit_dieback + scaled_distance +
               scaled_diameter + control_strike:scaled_diameter +
               control_strike:scaled_distance + (1|strike_paired), data = aggregated_data, 
             family = poisson, control = glmerControl(optCtrl=list(maxfun=2e4)))
AIC(gmod1)-AIC(gmod2)
#if AIC decreases by 2 or more, then the term contributed significantly
#marginally significant effect of dieback versus control/strike

#drop control_strike versus distance and test for an effect
gmod3<-glmer(Emergence.Holes.sum~control_strike + logit_dieback + scaled_distance +
               scaled_diameter + control_strike:scaled_diameter + control_strike:logit_dieback +
               (1|strike_paired), data = aggregated_data, 
             family = poisson, control = glmerControl(optCtrl=list(maxfun=2e4)))
AIC(gmod1)-AIC(gmod3)
#significant effect of that interaction

#convergence warning is only marginally above tolerance, but worth investigating anyway
#lets check for convergence issues
tt <- getME(gmod3,"theta")
ll <- getME(gmod3,"lower")
min(tt[ll==0])
#this value is very far from zero, so we are all good

#double check the gradient calculation
derivs1 <- gmod3@optinfo$derivs
sc_grad1 <- with(derivs1,solve(Hessian,gradient))
max(abs(sc_grad1))
#calculating it by hand suggests that we actually are not having issues - the quick-and-dirty internal check was just flawed
#this lack-of-an-error is supported by the decent model fit shown above

max(pmin(abs(sc_grad1),abs(derivs1$gradient)))

#we can then redo the caluculations with numerdiv
library(numDeriv)
dd <- update(gmod3,devFunOnly=TRUE)
pars <- unlist(getME(gmod3,c("theta","fixef")))
grad2 <- grad(dd,pars)
hess2 <- hessian(dd,pars)
sc_grad2 <- solve(hess2,grad2)
max(pmin(abs(sc_grad2),abs(grad2)))
#this also shows that the gradient is far below 0.001 (the threshold defined in the package)

#try rerunning the model with these alterations
ss <- getME(gmod3,c("theta","fixef"))
gmod3_updated <- update(gmod3,start=ss,control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e4)))
#the model looked good with the original optimizer and now there is no convergence warning at all
#compare AIC values again
AIC(gmod1)-AIC(gmod3_updated)#significant effect of control_strike versus distance interaction

#drop control_strike versus diameter and test for an effect
gmod4<-glmer(Emergence.Holes.sum~control_strike + logit_dieback + scaled_distance +
               scaled_diameter +  control_strike:logit_dieback +
               control_strike:scaled_distance + (1|strike_paired), data = aggregated_data, 
             family = poisson, control = glmerControl(optCtrl=list(maxfun=2e4)))
AIC(gmod1)-AIC(gmod4)
#the control_strike versus diameter interaction is significant

#drop distance and test for an effect
gmod5<-glmer(Emergence.Holes.sum~ logit_dieback + scaled_distance +
               scaled_diameter + control_strike:scaled_diameter + control_strike:logit_dieback +
               control_strike:scaled_distance + (1|strike_paired), data = aggregated_data, 
             family = poisson, control = glmerControl(optCtrl=list(maxfun=2e4)))
AIC(gmod1)-AIC(gmod5)
#control_strike alone is exceptionally significant

