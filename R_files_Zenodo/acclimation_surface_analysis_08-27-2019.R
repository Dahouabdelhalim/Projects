

#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

#### Clean analysis and estimation of acclimation surface    ####

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# Load libraries
library(dplyr)
library(reshape2)
library(growthTools)
library(bbmle)
library(ggplot2)
library(gridExtra)
library(mleTools)
library(growthTools)
library(mgcv)

# If desired, re-run data processing & growth rate estimation:
#source("./code/acclimation_surface_data_processing_073119.R")

# equation for fitting droop model to acclimation surface, assuming replete nutrients
# see 'dynamic quota acclimation model_v3.nb'
droop.surf<-function(a,b,Topt,w,q0,q1,vmax,m,Ta,Tc){
  nbcurve(Ta,Topt,w,a,b)*(1-(q0+q1*Ta)/(vmax/nbcurve(Tc,Topt,w,a,b)+q0+q1*Tc))-m
}

# setting up custom color palette for surfaces
or2<-c('#fdd49e','#fdbb84','#fc8d59','#e34a33','#b30000')
bl2<-c('#d0d1e6','#a6bddb','#74a9cf','#2b8cbe','#045a8d')
or2.bl2.pal<-c(colorRampPalette(colors=rev(bl2))(100),colorRampPalette(colors=or2)(100))


#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

#### Load data    ####

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

dat<-read.csv("~/Acc_surface_growth_data.csv")

gr.cr <- subset(dat, Species=="C. reinhardtii")
gr.ma <- subset(dat, Species=="M. aeruginosa")


#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

#### C. reinhardtii analyses    ####

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

head(gr.cr)

# Desired analyses:
# (1) establish acclimated curve
# (2) plot acute curves?
# (3) characterize surface
#     (a) GAMs
#     (b) Droop model


#### (1) Establish acclimated curve ####

# Subset data to focus on acclimated curve
gr.cr.acc<-gr.cr[gr.cr$Acute.temp==gr.cr$Acclim.temp,]
gr.cr.acc$Pass<-as.factor(gr.cr.acc$Pass)

# Fit GAMs
gm.cr.0<-gam(mu~s(Acclim.temp,k=5),bs='ps',method='REML',data=gr.cr.acc)
gm.cr.1<-gam(mu~s(Acclim.temp,k=5)+Pass,bs='ps',method='REML',data=gr.cr.acc)
gm.cr.test<-gam(mu~s(Acclim.temp,k=5,by=Pass)+Pass,bs='ps',method='REML',data=gr.cr.acc)
summary(gm.cr.test)
summary(gm.cr.1)
AICtab(gm.cr.0,gm.cr.1,gm.cr.test)

# Plot results
g1.cr<-ggplot(gr.cr.acc,aes(x=Acclim.temp,y=mu))+
  stat_smooth(method='gam',formula=y~s(x,k=5),aes(fill=factor(Pass),colour=factor(Pass)))+
  geom_hline(yintercept=0)+
  geom_point(aes(colour=factor(Pass)),alpha=0.5)+
  scale_x_continuous('Acclimated temperature (C)')+
  scale_y_continuous('Growth rate (1/day)')+
  scale_colour_discrete('Pass')+
  scale_fill_discrete('Pass')+
  theme_bw()+
  ggtitle('Acclimated curve')
g1.cr

# suggests there are differences in elevation among passes, but no changes in shape across passes (interaction term is not significant)
lm1.acc.cr<-lm(mu~factor(Acclim.temp)*factor(Pass),data=gr.cr.acc)
anova(lm1.acc.cr)


#### (2) Plot acute curves (INCOMPLETE) ####
tmp<-gr.cr[gr.cr$Acute.temp==gr.cr$Acclim.temp,]
tmp<-tmp[,-which(names(tmp)=='Treatment.label')]
acclim.tab<-tmp[tmp$Acute.temp<42,-which(names(tmp)=='Pass')]

acuteCR<-gr.cr %>% filter(Treatment %in% c(1,3,4,6,7)) %>% filter(Acute.temp<41) %>%
  ggplot(aes(x=Acute.temp,y=mu))+
  stat_smooth(data=acclim.tab,method='gam',
              formula=y~s(x,k=4),colour=gray(0.5))+
  stat_smooth(method='gam',formula=y~s(x,k=4),
              colour='blue',fill='blue')+
  geom_point(colour='blue',alpha=0.4)+
  geom_hline(yintercept = 0)+
  geom_vline(aes(xintercept=Acclim.temp),size=0.2)+
  facet_grid(Pass~Treatment.label2)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_x_continuous('Acute temperature (°C)')+
  scale_y_continuous('Growth rate (1/day)')
acuteCR

#ggsave(filename = '~/Dropbox/Phyto acclimation/accData/results/poster/acute_vs_acclimated_trial4.pdf',acuteCR,width=10,height=3.5*1.8)


#### (3) Characterize acclimation surface ####

# data from Pass 1, excluding high temperatures
#gr.cr.P1<-gr.cr %>% filter(Pass==1 & Acclim.temp<40 & Acute.temp<43)
gr.cr.P1<-gr.cr %>% filter(Pass==1 & Acclim.temp<40)

tmpB<-gr.cr %>% filter(Pass==1) %>% group_by(Acclim.temp,Acute.temp) %>% summarise(mean.mu=mean(mu),var.mu=var(mu))
tmpB %>% filter(Acute.temp>40)
ggplot(tmpB,aes(x=Acute.temp,y=Acclim.temp))+
  geom_point(aes(colour=mean.mu),size=7)+
  scale_x_continuous('Acute temperature (°C)')+
  scale_y_continuous('Acclimated temperature (°C)')+
  scale_colour_gradientn('Growth \\n   rate \\n(1/day)',
                         colours=or2.bl2.pal,
                         breaks=seq(-2.0,2.0,0.5),limits=c(-2.1,2.1),
                         guide = guide_colourbar(ticks.colour = gray(gl),
                                                 ticks.linewidth = 1.,
                                                 frame.colour = "black"))+
  theme_bw()

tmpB %>% filter(Acute.temp>40)
ggplot(tmpB,aes(x=Acute.temp,y=Acclim.temp))+
  geom_point(aes(colour=var.mu),size=7)+
  scale_x_continuous('Acute temperature (°C)')+
  scale_y_continuous('Acclimated temperature (°C)')+
  scale_colour_distiller(palette = 'Spectral')+
  theme_bw()


# Excluding data where Acute.temp >= 43, because multiple negative growth rates occur already
# at ~41 C
# Excluding Acclim.temp >= 40 because these values aren't supported by a corresponding acute
# temperature treatment.
# table(gr.cr[,c('Acclim.temp','Acute.temp')])


####  ~ (a)  Fit GAM surface ####

# After some exploration, best approach appears to be:
# use method = REML - more resilient to overfitting given smaller data sets
# leave gamma = 1
# exclude models using isomorphic s() structure in favor of tensor product te(), 
# which better accommodates different levels of wiggliness in x vs. y dimensions
# model comparison favors p-splines

#mthd<-'GCV.Cp'
mthd<-'REML'
gm<-1
#gm<-1.4
#gm<-2
gm1.te.tp<-gam(mu~te(Acclim.temp,Acute.temp,k=5,bs='tp'),method = mthd,gamma=gm,data=gr.cr.P1)
gm1.te.cr<-gam(mu~te(Acclim.temp,Acute.temp,k=5,bs='cr'),method = mthd,gamma=gm,data=gr.cr.P1)
gm1.te.ps<-gam(mu~te(Acclim.temp,Acute.temp,k=5,bs='ps'),method = mthd,gamma=gm,data=gr.cr.P1)
gm1.te.ps.k6<-gam(mu~te(Acclim.temp,Acute.temp,k=6,bs='ps'),method = mthd,gamma=gm,data=gr.cr.P1)
gm1.te.ps.k7<-gam(mu~te(Acclim.temp,Acute.temp,k=7,bs='ps'),method = mthd,gamma=gm,data=gr.cr.P1)
gm1.te.ds<-gam(mu~te(Acclim.temp,Acute.temp,k=5,bs='ds'),method = mthd,gamma=gm,data=gr.cr.P1)
gm1.s.tp<-gam(mu~s(Acclim.temp,Acute.temp,bs='tp'),method = mthd,gamma=gm,data=gr.cr.P1)
gm1.s.tp.k20<-gam(mu~s(Acclim.temp,Acute.temp,bs='tp',k=20),method = mthd,gamma=gm,data=gr.cr.P1)
gm1.s.tp.k25<-gam(mu~s(Acclim.temp,Acute.temp,bs='tp',k=25),method = mthd,gamma=gm,data=gr.cr.P1)
gm1.s.tp.k30<-gam(mu~s(Acclim.temp,Acute.temp,bs='tp',k=30),method = mthd,gamma=gm,data=gr.cr.P1)
gm1.s.tp.k35<-gam(mu~s(Acclim.temp,Acute.temp,bs='tp',k=35),method = mthd,gamma=gm,data=gr.cr.P1)

# full model comparison table
AICtab(gm1.te.tp,gm1.te.cr,gm1.te.ps,gm1.te.ps.k6,gm1.te.ps.k7,gm1.te.ds,gm1.s.tp,gm1.s.tp.k25,gm1.s.tp.k30,gm1.s.tp.k35,gm1.s.tp.k20)

# with REML
# dAIC  df              
# gm1.s.tp.k35   0.0 34.2734845324714
# gm1.te.ps.k7   8.3 28.4755388977701
# gm1.s.tp      31.4 29.6428853633897
# gm1.s.tp.k30  31.4 29.6428853633897
# gm1.te.ps.k6  39.1 24.1172389248475
# gm1.s.tp.k25  85.6 24.656862314551 
# gm1.te.ps     98.0 16.7143579298919
# gm1.s.tp.k20  98.6 20.3337790483785
# gm1.te.tp    139.1 17.513787336352 
# gm1.te.ds    139.1 17.5137873363531
# gm1.te.cr    147.6 17.0648477332015

plot.cr.model<-function(gam.mod,title='gam surface'){
  
  # temperature grid:
  vals<-expand.grid(Acute.temp=seq(14,45,length.out = 100),Acclim.temp=seq(14,41,length.out = 100))
  acc.slices<-unique(gr.cr.P1$Acclim.temp)[c(-2,-5)]
  acute.vals<-expand.grid(Acute.temp=seq(14,45,length.out = 100),Acclim.temp=acc.slices)
  
  # predictions from GAM
  pds.gam<-predict(gam.mod,newdata=vals)
  pds.gam.slices<-predict(gam.mod,newdata=acute.vals,se=T)
  
  # organize predictions
  cr.surf<-data.frame(vals,gam=pds.gam)
  cr.surf<-melt(cr.surf,id.vars=c('Acute.temp','Acclim.temp'))
  cr.slice<-data.frame(acute.vals,mu=pds.gam.slices$fit,mu.se=pds.gam.slices$se.fit)
  
  a1.plot<-ggplot(cr.slice,aes(x=Acute.temp,y=mu))+
    geom_hline(yintercept = 0)+
    geom_ribbon(aes(fill=factor(Acclim.temp),ymin=mu-1.96*mu.se,ymax=mu+1.96*mu.se),alpha=0.2)+
    geom_line(aes(colour=factor(Acclim.temp)))+
    geom_point(data=filter(gr.cr.P1,Acclim.temp %in% acc.slices),aes(colour=factor(Acclim.temp)),alpha=0.4,size=0.5)+
    facet_wrap(~Acclim.temp,nrow=1)+
    coord_cartesian(ylim = c(-2.05,2))+
    scale_y_continuous('Growth rate')+
    theme_bw()+
    theme(legend.position = 'none')
  #a1.plot
  
  a2.plot<-ggplot(cr.surf,aes(x=Acute.temp,y=Acclim.temp,z=value))+
    geom_raster(aes(fill=value))+
    geom_abline(intercept=0,slope=1,linetype=2)+
    geom_contour(colour=gray(0.2),binwidth=0.1,size=0.1)+
    geom_contour(colour=gray(0.2),binwidth=0.5,size=0.4)+
    geom_contour(colour='red',breaks=c(0),size=1)+
    geom_point(data=emp.gr.cr,fill='black',colour='black',size=5.2,shape=21)+
    geom_point(data=emp.gr.cr,aes(fill=mu),colour='black',size=5,shape=21)+
    scale_x_continuous('Acute temperature (°C)',expand=c(0,0))+
    scale_y_continuous('Acclimated temperature (°C)',expand=c(0,0))+
    scale_fill_gradientn('Growth \\n   rate \\n(1/day)',
                         colours=or2.bl2.pal,
                         breaks=seq(-2.0,2.0,0.5),limits=c(-2.1,2.1),
                         guide = guide_colourbar(ticks.colour = gray(gl),
                                                 ticks.linewidth = 1.,
                                                 frame.colour = "black"))+
    coord_cartesian(ylim=c(15,38),xlim=c(15,45))+
    facet_wrap(~variable,nrow=2)+
    theme_bw()+
    ggtitle(title)
  #a2.plot
  
  grid.arrange(a2.plot,a1.plot,ncol=1,layout_matrix=rbind(c(1),c(1),c(2)))
}

# Model plots
plot.cr.model(gm1.te.ps.k7,'te.ps.k7')
#plot.cr.model(gm1.s.tp.k35,'s.tp.k35')

# Digression:
# something to perhaps remember - sample size slightly larger where acute=acclimated temperature... probably doesn't have a large effect though.
acc.slices<-unique(gr.cr.P1$Acclim.temp)[c(-2,-5)]
gr.cr.P1 %>% filter(Acclim.temp %in% acc.slices) %>%
  ggplot(aes(x=Acute.temp,y=mu))+
  geom_point(aes(shape=factor(TGB_lane),colour=factor(TGB_lane)))+
  facet_wrap(~Acclim.temp)+
  theme_bw()

# model comparison
AICtab(gm1.te.tp,gm1.te.cr,gm1.te.ps,gm1.te.ps.k6,gm1.te.ps.k7,gm1.te.ds)
# dAIC  df              
# gm1.te.ps.k7   0.0 28.4755388977701
# gm1.te.ps.k6  30.7 24.1172389248475
# gm1.te.ps     89.7 16.7143579298919
# gm1.te.tp    130.8 17.513787336352 
# gm1.te.ds    130.8 17.5137873363531
# gm1.te.cr    139.3 17.0648477332015

summary(gm1.te.ps.k7)
# 
# Family: gaussian 
# Link function: identity 
# 
# Formula:
#   mu ~ te(Acclim.temp, Acute.temp, k = 7, bs = "ps")
# 
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.58479    0.01345   43.49   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df     F p-value    
# te(Acclim.temp,Acute.temp) 24.46  27.83 122.9  <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.957   Deviance explained = 96.3%
# -REML = -8.6847  Scale est. = 0.02821   n = 156

gam.check(gm1.te.ps.k7,k.rep = 10000)
# Method: REML   Optimizer: outer newton
# full convergence after 5 iterations.
# Gradient range [-5.71043e-08,4.706706e-09]
# (score -8.684675 & scale 0.02820983).
# Hessian positive definite, eigenvalue range [1.602899,76.78971].
# Model rank =  49 / 49 
# 
# Basis dimension (k) checking results. Low p-value (k-index<1) may
# indicate that k is too low, especially if edf is close to k'.
# 
#                              k'  edf k-index p-value  
# te(Acclim.temp,Acute.temp) 48.0 24.5    0.85    0.01 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

surf.diag.cr<-plot.cr.model(gm1.te.ps.k7,'C. reinhardtii, model te.ps.k7')
ggsave(filename = './results/acclimation_surfaces/acclimation_surface_CR_model_selection.pdf',surf.diag.cr,width=6.75,height=7.75)

# Save preferred fit using short name:
gm1.cr<-gm1.te.ps.k7



####  ~ (b)  Fit Droop surface ####

# fit droop surface model

# only trimming out the two highest acclimated obs:
mD.cr<-mle2(mu~dnorm(mean=droop.surf(exp(a),exp(b),Topt,10*w,10*q0,q1,10*vmax,exp(m),Acute.temp,Acclim.temp),sd=exp(s)),start=list(a=-13.3787,b=-2.3797,Topt=30.81897,w=36.3326,q0=86.497,q1=-9.10588,vmax=128.017566,m=-9.74,s=log(1)),control=list(maxit=10000),data=gr.cr.P1)
summary(mD.cr)

pf.cr<-profile(mD.cr)
#plot(pf.cr)
# not very smooth profiles... (maybe b/c of how I'm cheating parscale? hmm) but fit seems stable

AICctab(mD.cr,gm1.cr,nobs=nrow(gr.cr.P1))
#        dAICc df              
# gm1.cr  0.0  28.4755388977701
# mD.cr  96.4  9 

# R2 - Droop
get.R2(predict(mD.cr),gr.cr.P1$mu)
# 0.9060624

# R2 - GAM
get.R2(predict(gm1.cr),gr.cr.P1$mu)
# 0.9634787

# Droop model somewhat worse predictive power, but with half as many 
# parameters as the GAM. NICE!

####  ~ (c)  Visualize surfaces ####

# temperature grid:
vals<-expand.grid(Acute.temp=seq(14,45,length.out = 100),Acclim.temp=seq(14,41,length.out = 100))

# predictions from GAM
pds.gam<-predict(gm1.cr,newdata=vals)
pds.drp<-predict(mD.cr,newdata=vals)

# combine data
cr.surf<-data.frame(vals,gam=pds.gam,droop=pds.drp)
cr.surf<-melt(cr.surf,id.vars=c('Acute.temp','Acclim.temp'))
head(cr.surf)

# sample locations
sample.locs<-unique(gr.cr[gr.cr$Acclim.temp<40,c('Acute.temp','Acclim.temp')])
sample.locs$value<-NA

cr.surf.plt<-ggplot(cr.surf,aes(x=Acute.temp,y=Acclim.temp,z=value))+
  geom_raster(aes(fill=value))+
  geom_abline(intercept=0,slope=1,linetype=2)+
  geom_contour(colour=gray(0.2),binwidth=0.25,size=0.2)+
  geom_contour(colour='red',breaks=c(0),size=1)+
  geom_point(data=sample.locs,alpha=0.5)+
  scale_x_continuous('Acute temperature (°C)',expand=c(0,0))+
  scale_y_continuous('Acclimated temperature (°C)',expand=c(0,0))+
  scale_fill_gradientn('Growth \\n   rate \\n(1/day)',
                         colours=or2.bl2.pal,
                         breaks=seq(-2.0,2.0,0.5),limits=c(-2.1,2.1),
                         guide = guide_colourbar(ticks.colour = gray(gl),
                                                 ticks.linewidth = 1.,
                                                 frame.colour = "black"))+
  coord_cartesian(xlim=c(15,45),ylim=c(15,38))+
  facet_wrap(~variable)+
  theme_bw()
cr.surf.plt

emp.gr.cr<-gr.cr.P1 %>% group_by(Acclim.temp,Acute.temp) %>% summarise(mu=median(mu))
#emp.gr.cr %>% filter(Acute.temp > 40)

emp.cr.surf.plt<-ggplot(emp.gr.cr,aes(x=Acute.temp,y=Acclim.temp))+
  geom_point(aes(colour=mu),size=7)+
  scale_x_continuous('Acute temperature (°C)')+
  scale_y_continuous('Acclimated temperature (°C)')+
  scale_colour_gradientn('Growth \\n   rate \\n(1/day)',
                         colours=or2.bl2.pal,
                         breaks=seq(-2.0,2.0,0.5),limits=c(-2.1,2.1),
                         guide = guide_colourbar(ticks.colour = gray(gl),
                                                 ticks.linewidth = 1.,
                                                 frame.colour = "black"))+
  theme_bw()
emp.cr.surf.plt

# combine panels
cr.surf.plt2<-cr.surf.plt+theme(legend.position = 'none')
grid.arrange(emp.cr.surf.plt,cr.surf.plt2,nrow=2,layout_matrix=rbind(c(1,1,NA),c(2,2,2)))
a1.cr<-arrangeGrob(emp.cr.surf.plt,cr.surf.plt2,nrow=2,layout_matrix=rbind(c(1,1,NA),c(2,2,2)))

ggsave(filename = '~/Dropbox/Phyto acclimation/accData/results/acclimation_surfaces/acclimation_surfaces_CR_P1_082719.pdf',a1.cr,width=8.25,height=8)

# Overlay empirical obs on estimated surfaces
emp.gr.cr$value<-NA

a2.cr<-ggplot(cr.surf,aes(x=Acute.temp,y=Acclim.temp,z=value))+
  geom_raster(aes(fill=value))+
  geom_abline(intercept=0,slope=1,linetype=2)+
  geom_contour(colour=gray(0.2),binwidth=0.1,size=0.1)+
  geom_contour(colour=gray(0.2),binwidth=0.5,size=0.4)+
  geom_contour(colour='red',breaks=c(0),size=1)+
  geom_point(data=emp.gr.cr,fill='black',colour='black',size=5.2,shape=21)+
  geom_point(data=emp.gr.cr,aes(fill=mu),colour='black',size=5,shape=21)+
  scale_x_continuous('Acute temperature (°C)',expand=c(0,0))+
  scale_y_continuous('Acclimated temperature (°C)',expand=c(0,0))+
  scale_fill_gradientn('Growth \\n   rate \\n(1/day)',
                         colours=or2.bl2.pal,
                         breaks=seq(-2.0,2.0,0.5),limits=c(-2.1,2.1),
                         guide = guide_colourbar(ticks.colour = gray(gl),
                                                 ticks.linewidth = 1.,
                                                 frame.colour = "black"))+
  coord_cartesian(ylim=c(15,38),xlim=c(15,45))+
  facet_wrap(~variable,nrow=2)+
  theme_bw()+
  ggtitle('C. reinhardtii')
a2.cr

# visualize implied acclimated curve:

# temperature grid:
vals2<-data.frame(Acute.temp=seq(14,45,length.out = 100),Acclim.temp=seq(14,45,length.out = 100))

# predictions from GAM
pds.gam.acc<-predict(gm1.cr,newdata=vals2,se=T)
vals2<-data.frame(vals2,mu=pds.gam.acc$fit,se=pds.gam.acc$se.fit)


plist<-list()

plist[[1]]<-ggplot(vals2,aes(x=Acclim.temp,y=mu))+
    geom_hline(yintercept = 0)+
    geom_ribbon(aes(ymin=mu-1.96*se,ymax=mu+1.96*se),fill='gray')+
    geom_line()+
    stat_smooth(method='gam',formula=y~s(x,k=4,bs='ps'))+
    geom_point(data=gr.cr.P1[gr.cr.P1$Acclim.temp==gr.cr.P1$Acute.temp,])+
    coord_cartesian(ylim = c(-1,1.5))+
    theme_bw()+
    ggtitle('k = 4')

plist[[2]]<-ggplot(vals2,aes(x=Acclim.temp,y=mu))+
  geom_hline(yintercept = 0)+
  geom_ribbon(aes(ymin=mu-1.96*se,ymax=mu+1.96*se),fill='gray')+
  geom_line()+
  stat_smooth(method='gam',formula=y~s(x,k=5,bs='ps'))+
  geom_point(data=gr.cr.P1[gr.cr.P1$Acclim.temp==gr.cr.P1$Acute.temp,])+
  coord_cartesian(ylim = c(-1,1.5))+
  theme_bw()+
  ggtitle('k = 5')

plist[[3]]<-ggplot(vals2,aes(x=Acclim.temp,y=mu))+
  geom_hline(yintercept = 0)+
  geom_ribbon(aes(ymin=mu-1.96*se,ymax=mu+1.96*se),fill='gray')+
  geom_line()+
  stat_smooth(method='gam',formula=y~s(x,k=6,bs='ps'))+
  geom_point(data=gr.cr.P1[gr.cr.P1$Acclim.temp==gr.cr.P1$Acute.temp,])+
  coord_cartesian(ylim = c(-1,1.5))+
  theme_bw()+
  ggtitle('k = 6')

plist[[4]]<-ggplot(vals2,aes(x=Acclim.temp,y=mu))+
  geom_hline(yintercept = 0)+
  geom_ribbon(aes(ymin=mu-1.96*se,ymax=mu+1.96*se),fill='gray')+
  geom_line()+
  stat_smooth(method='gam',formula=y~s(x,k=7,bs='ps'))+
  geom_point(data=gr.cr.P1[gr.cr.P1$Acclim.temp==gr.cr.P1$Acute.temp,])+
  coord_cartesian(ylim = c(-1,1.5))+
  theme_bw()+
  ggtitle('k = 7')

grid.arrange(plist[[1]],plist[[2]],plist[[3]],plist[[4]],top='C. reinhardtii acclimated curve')


ggplot(vals2,aes(x=Acclim.temp,y=mu))+
  geom_hline(yintercept = 0)+
  geom_ribbon(aes(ymin=mu-1.96*se,ymax=mu+1.96*se),fill='gray')+
  geom_line()+
  stat_smooth(method='gam',formula=y~s(x,k=9))+
  geom_point(data=gr.cr.P1[gr.cr.P1$Acclim.temp==gr.cr.P1$Acute.temp,])+
  coord_cartesian(ylim = c(-1,1.5))+
  theme_bw()

cr.acc.gm<-gam(mu~s(Acclim.temp,k=7,bs='ps'),data=gr.cr.P1[gr.cr.P1$Acclim.temp==gr.cr.P1$Acute.temp,])
summary(cr.acc.gm)

cr.acc.gm<-gam(mu~s(Acclim.temp,k=6,bs='ps'),data=gr.cr.P1[gr.cr.P1$Acclim.temp==gr.cr.P1$Acute.temp,])
summary(cr.acc.gm)


#### (4) Save surface for downstream analyses ####

#write.csv(x = gr.cr.P1,"./data/Acc_surface/derived_data/acc_surface_CR_t4_GAM_subset_082719.csv",row.names=F)
write.csv(x = cr.surf,"./data/Acc_surface/derived_data/acc_surface_CR_t4_GAM_subset_092419.csv",row.names=F)


#### (5) Examine residuals ####

head(emp.gr.cr)

pds.cr2<-predict(gm1.cr,newdata=emp.gr.cr)
emp.gr.cr2<-emp.gr.cr
emp.gr.cr2$predicted<-pds.cr2

ggplot(emp.gr.cr2,aes(x=Acute.temp))+
  geom_hline(yintercept = 0)+
  geom_point(aes(y=mu),colour='black')+
  geom_point(aes(y=predicted),colour='red',size=0.5)+
  geom_line(aes(y=predicted),colour='red')+
  facet_wrap(~Acclim.temp)+
  theme_bw()



#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

#### M. aeruginosa analyses                                  ####

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

#### (1) Establish acclimated curve ####

# Subset data to focus on acclimated curve
gr.ma.acc<-gr.ma[gr.ma$Acute.temp==gr.ma$Acclim.temp,]
gr.ma.acc$Pass<-as.factor(gr.ma.acc$Pass)

# Fit GAMs
gm.ma.0<-gam(mu~s(Acclim.temp,k=5),bs='ps',method='REML',data=gr.ma.acc)
gm.ma.1<-gam(mu~s(Acclim.temp,k=5)+Pass,bs='ps',method='REML',data=gr.ma.acc)
gm.ma.test<-gam(mu~s(Acclim.temp,k=5,by=Pass)+Pass,bs='ps',method='REML',data=gr.ma.acc)
summary(gm.ma.test)
summary(gm.ma.1)
AICtab(gm.ma.0,gm.ma.1,gm.ma.test)

# Plot results
g1.ma<-ggplot(gr.ma.acc,aes(x=Acclim.temp,y=mu))+
  stat_smooth(method='gam',formula=y~s(x,k=5),aes(fill=factor(Pass),colour=factor(Pass)))+
  geom_hline(yintercept=0)+
  geom_point(aes(colour=factor(Pass)),alpha=0.5)+
  scale_x_continuous('Acclimated temperature (C)')+
  scale_y_continuous('Growth rate (1/day)')+
  scale_colour_discrete('Pass')+
  scale_fill_discrete('Pass')+
  theme_bw()+
  ggtitle('Acclimated curve')
g1.ma

#ggsave(filename = '~/Dropbox/Phyto acclimation/accData/results/acclimation_surfaces/acclimated_curves_MA.pdf',g1.ma)

# suggests there are differences in elevation and shape among passes, even w/out hottest temperature...
lm1.acc.ma<-lm(mu~factor(Acclim.temp)*factor(Pass),data=gr.ma.acc[gr.ma.acc$Acclim.temp<43,])
anova(lm1.acc.ma)


#### (2) Plot acute curves (INCOMPLETE) ####

#### (3) Characterize acclimation surface ####

tmpB<-gr.ma %>% filter(Pass==1) %>% group_by(Acclim.temp,Acute.temp) %>% summarise(mean.mu=mean(mu),var.mu=var(mu))
tmpB %>% filter(Acute.temp>40)
ggplot(tmpB,aes(x=Acute.temp,y=Acclim.temp))+
  geom_point(aes(colour=mean.mu),size=7)+
  scale_x_continuous('Acute temperature (°C)')+
  scale_y_continuous('Acclimated temperature (°C)')+
  scale_colour_gradientn('Growth \\n   rate \\n(1/day)',
                         colours=or2.bl2.pal,
                         breaks=seq(-2.0,2.0,0.5),limits=c(-2.1,2.1),
                         guide = guide_colourbar(ticks.colour = gray(gl),
                                                 ticks.linewidth = 1.,
                                                 frame.colour = "black"))+
  theme_bw()

tmpB %>% filter(Acute.temp>40)
ggplot(tmpB,aes(x=Acute.temp,y=Acclim.temp))+
  geom_point(aes(colour=var.mu),size=7)+
  scale_x_continuous('Acute temperature (°C)')+
  scale_y_continuous('Acclimated temperature (°C)')+
  scale_colour_distiller(palette = 'Spectral')+
  theme_bw()


# data from Pass 1, excluding high temperatures
gr.ma.P1<-gr.ma %>% filter(Pass==1 & Acclim.temp<40 & Acute.temp<43)

# Notes: decided to exclude 
# Acute = 45 C (as low/negative growth rates occur already at ~41 C)
# Acclimated >37.6 C, b/c excluding acute 45 and because there's no acute cross at 41 C
# table(gr.ma[,c('Acclim.temp','Acute.temp')])

####  ~ (a)  Fit GAM surface ####

# Settled on method = 'REML' and gm = 1 (default)
#mthd<-'GCV.Cp'
mthd<-'REML'
gm<-1
#gm<-1.4
#gm<-2
gm2.te.tp<-gam(mu~te(Acclim.temp,Acute.temp,k=5,bs='tp'),method = mthd,gamma=gm,data=gr.ma.P1)
gm2.te.cr<-gam(mu~te(Acclim.temp,Acute.temp,k=5,bs='cr'),method = mthd,gamma=gm,data=gr.ma.P1)
gm2.te.ps<-gam(mu~te(Acclim.temp,Acute.temp,k=5,bs='ps'),method = mthd,gamma=gm,data=gr.ma.P1)
gm2.te.ps.k6<-gam(mu~te(Acclim.temp,Acute.temp,k=6,bs='ps'),method = mthd,gamma=gm,data=gr.ma.P1)
gm2.te.ps.k7<-gam(mu~te(Acclim.temp,Acute.temp,k=7,bs='ps'),method = mthd,gamma=gm,data=gr.ma.P1)
gm2.te.ds<-gam(mu~te(Acclim.temp,Acute.temp,k=5,bs='ds'),method = mthd,gamma=gm,data=gr.ma.P1)
gm2.s.tp<-gam(mu~s(Acclim.temp,Acute.temp,bs='tp'),method = mthd,gamma=gm,data=gr.ma.P1)
gm2.s.tp.k20<-gam(mu~s(Acclim.temp,Acute.temp,bs='tp',k=20),method = mthd,gamma=gm,data=gr.ma.P1)
gm2.s.tp.k25<-gam(mu~s(Acclim.temp,Acute.temp,bs='tp',k=25),method = mthd,gamma=gm,data=gr.ma.P1)
gm2.s.tp.k30<-gam(mu~s(Acclim.temp,Acute.temp,bs='tp',k=30),method = mthd,gamma=gm,data=gr.ma.P1)
gm2.s.tp.k35<-gam(mu~s(Acclim.temp,Acute.temp,bs='tp',k=35),method = mthd,gamma=gm,data=gr.ma.P1)

# model comparison
AICtab(gm2.te.tp,gm2.te.cr,gm2.te.ps,gm2.te.ps.k6,gm2.te.ps.k7,gm2.te.ds,gm2.s.tp,gm2.s.tp.k25,gm2.s.tp.k30,gm2.s.tp.k35,gm2.s.tp.k20)

# dAIC  df              
# gm2.te.ps.k7   0.0 35.713590947671 
# gm2.s.tp.k35   1.6 35.1523349041778
# gm2.s.tp       9.4 30.5091450910089
# gm2.s.tp.k30   9.4 30.5091450910089
# gm2.te.ps.k6  39.6 31.455435234698 
# gm2.s.tp.k25  50.7 25.6625805129856
# gm2.te.ps     76.9 25.245278036697 
# gm2.te.tp    105.8 24.8041615238945
# gm2.te.ds    105.8 24.8041615238946
# gm2.te.cr    107.9 24.7396952923601
# gm2.s.tp.k20 137.5 20.6923848400447

# stash empirical observations for plotting
emp.gr.ma<-gr.ma.P1 %>% group_by(Acclim.temp,Acute.temp) %>% summarise(mu=median(mu))
emp.gr.ma$value<-NA

# plotting function for visualizing fit
plot.ma.model<-function(gam.mod,title='gam surface'){
  
  # temperature grid:
  vals<-expand.grid(Acute.temp=seq(14,45,length.out = 100),Acclim.temp=seq(14,41,length.out = 100))
  acc.slices<-unique(gr.ma.P1$Acclim.temp)[c(-3,-5)]
  acute.vals<-expand.grid(Acute.temp=seq(14,45,length.out = 100),Acclim.temp=acc.slices)
  
  # predictions from GAM
  pds.gam<-predict(gam.mod,newdata=vals)
  pds.gam.slices<-predict(gam.mod,newdata=acute.vals,se=T)
  
  # organize predictions
  ma.surf<-data.frame(vals,gam=pds.gam)
  ma.surf<-melt(ma.surf,id.vars=c('Acute.temp','Acclim.temp'))
  ma.slice<-data.frame(acute.vals,mu=pds.gam.slices$fit,mu.se=pds.gam.slices$se.fit)
  
  a1.plot<-ggplot(ma.slice,aes(x=Acute.temp,y=mu))+
    geom_hline(yintercept = 0)+
    geom_ribbon(aes(fill=factor(Acclim.temp),ymin=mu-1.96*mu.se,ymax=mu+1.96*mu.se),alpha=0.2)+
    geom_line(aes(colour=factor(Acclim.temp)))+
    geom_point(data=filter(gr.ma.P1,Acclim.temp %in% acc.slices),aes(colour=factor(Acclim.temp)),alpha=0.4,size=0.5)+
    facet_wrap(~Acclim.temp,nrow=1)+
    coord_cartesian(ylim = c(-0.75,1.1))+
    scale_y_continuous('Growth rate')+
    theme_bw()+
    theme(legend.position = 'none')
  #a1.plot
  
  a2.plot<-ggplot(ma.surf,aes(x=Acute.temp,y=Acclim.temp,z=value))+
    geom_raster(aes(fill=value))+
    geom_abline(intercept=0,slope=1,linetype=2)+
    geom_contour(colour=gray(0.2),binwidth=0.1,size=0.1)+
    geom_contour(colour=gray(0.2),binwidth=0.5,size=0.4)+
    geom_contour(colour='red',breaks=c(0),size=1)+
    geom_point(data=emp.gr.ma,fill='black',colour='black',size=5.2,shape=21)+
    geom_point(data=emp.gr.ma,aes(fill=mu),colour='black',size=5,shape=21)+
    scale_x_continuous('Acute temperature (°C)',expand=c(0,0))+
    scale_y_continuous('Acclimated temperature (°C)',expand=c(0,0))+
    scale_fill_gradientn('Growth \\n   rate \\n(1/day)',
                         colours=or2.bl2.pal,
                         breaks=seq(-2.0,2.0,0.5),limits=c(-2.1,2.1),
                         guide = guide_colourbar(ticks.colour = gray(gl),
                                                 ticks.linewidth = 1.,
                                                 frame.colour = "black"))+
    coord_cartesian(ylim=c(15,38),xlim=c(15,45))+
    facet_wrap(~variable,nrow=2)+
    theme_bw()+
    ggtitle(title)
  #a2.plot
  
  grid.arrange(a2.plot,a1.plot,ncol=1,layout_matrix=rbind(c(1),c(1),c(2)))
}

# model comparison, excluding s() models:
AICtab(gm2.te.tp,gm2.te.cr,gm2.te.ps,gm2.te.ps.k6,gm2.te.ps.k7,gm2.te.ds)
# dAIC  df              
# gm2.te.ps.k7   0.0 35.713590947671 
# gm2.te.ps.k6  39.6 31.455435234698 
# gm2.te.ps     76.9 25.245278036697 
# gm2.te.tp    105.8 24.8041615238945
# gm2.te.ds    105.8 24.8041615238946
# gm2.te.cr    107.9 24.7396952923601

# visualize preferred fit
plot.ma.model(gm2.te.ps.k7)

# check adequacy of k=7:
gam.check(gm2.te.ps.k7,k.rep=10000)
# Method: REML   Optimizer: outer newton
# full convergence after 4 iterations.
# Gradient range [-5.984246e-06,2.139492e-06]
# (score -172.0421 & scale 0.001758535).
# Hessian positive definite, eigenvalue range [3.186829,70.23194].
# Model rank =  49 / 49 
# 
# Basis dimension (k) checking results. Low p-value (k-index<1) may
# indicate that k is too low, especially if edf is close to k'.
# 
#                              k'  edf k-index p-value
# te(Acclim.temp,Acute.temp) 48.0 32.9    1.06    0.74

summary(gm2.te.ps.k7)
# Family: gaussian 
# Link function: identity 
# 
# Formula:
#   mu ~ te(Acclim.temp, Acute.temp, k = 7, bs = "ps")
# 
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 0.481126   0.003532   136.2   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df   F p-value    
# te(Acclim.temp,Acute.temp) 32.87  34.78 383  <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =   0.99   Deviance explained = 99.2%
# -REML = -172.04  Scale est. = 0.0017585  n = 141

surf.diag.ma<-plot.ma.model(gm2.te.ps.k7,'M. aeruginosa, model te.ps.k7')
ggsave(filename = './results/acclimation_surfaces/acclimation_surface_MA_model_selection.pdf',surf.diag.ma,width=6.75,height=7.75)

# Save with short name
gm1.ma<-gm2.te.ps.k7
summary(gm1.ma)


####  ~ (b)  Fit Droop surface ####

# fit droop surface model
# gr.ma.P1<-gr.ma %>% filter(Pass==1 & Acclim.temp<40 & Acute.temp<43)
mD.ma<-mle2(mu~dnorm(mean=droop.surf(exp(a),exp(b),Topt,10*w,q0,q1/10,10*vmax,exp(m),Acute.temp,Acclim.temp),sd=exp(s)),start=list(a=-13.126,b=-3.06,Topt=24.61,w=12.28,q0=52.14,q1=-16,vmax=5.25,m=log(0.1),s=log(1)),control=list(maxit=10000),data=gr.ma.P1)
summary(mD.ma)

#pf<-profile(mD.ma)

AICctab(mD.ma,gm1.ma,nobs=nrow(gr.ma.P1))
#        dAICc df             
# gm1.ma   0.0 35.713590947671
# mD.ma  353.5 9 

# R2 - Droop
get.R2(predict(mD.ma),gr.ma.P1$mu)
# 0.8311315

# R2 - GAM
get.R2(predict(gm1.ma),gr.ma.P1$mu)
# 0.9920393

# Droop model achieves somewhat worse predictive power, but with half as many 
# parameters as the GAM, and is still damn good. NICE!

####  ~ (c)  Visualize surfaces ####

# temperature grid:
range(gr.ma.P1$Acclim.temp)
range(gr.ma.P1$Acute.temp)
vals<-expand.grid(Acute.temp=seq(15,45,length.out = 100),Acclim.temp=seq(15,40,length.out = 100))

# predictions from GAM
pds.gam<-predict(gm1.ma,newdata=vals)
pds.drp<-predict(mD.ma,newdata=vals)

# combine data
ma.surf<-data.frame(vals,gam=pds.gam,droop=pds.drp)
ma.surf<-melt(ma.surf,id.vars=c('Acute.temp','Acclim.temp'))
head(ma.surf)

# sample locations
sample.locs<-unique(gr.ma[gr.ma$Acclim.temp<43,c('Acute.temp','Acclim.temp')])
sample.locs$value<-NA

ma.surf.plt<-ggplot(ma.surf,aes(x=Acute.temp,y=Acclim.temp,z=value))+
  geom_raster(aes(fill=value))+
  geom_abline(intercept=0,slope=1,linetype=2)+
  geom_contour(colour=gray(0.2),binwidth=0.1,size=0.2)+
  geom_contour(colour='red',breaks=c(0),size=1)+
  geom_point(data=sample.locs,alpha=0.5)+
  scale_x_continuous('Acute temperature (°C)',limits = c(15,42),expand=c(0,0))+
  scale_y_continuous('Acclimated temperature (°C)',expand=c(0,0))+
  scale_fill_gradientn('Growth \\n   rate \\n(1/day)',
                         colours=or2.bl2.pal,
                         breaks=seq(-1.5,1.5,0.5),limits=c(-1.8,1.8),
                         guide = guide_colourbar(ticks.colour = gray(gl),
                                                 ticks.linewidth = 1.,
                                                 frame.colour = "black"))+
  coord_cartesian(ylim=c(15,38))+
  facet_wrap(~variable)+
  theme_bw()
ma.surf.plt

emp.gr.ma<-gr.ma.P1 %>% group_by(Acclim.temp,Acute.temp) %>% summarise(mu=mean(mu))

emp.ma.surf.plt<-ggplot(emp.gr.ma,aes(x=Acute.temp,y=Acclim.temp))+
  geom_point(aes(colour=mu),size=7)+
  scale_x_continuous('Acute temperature (°C)',limits = c(15,42))+
  scale_y_continuous('Acclimated temperature (°C)')+
  scale_colour_gradientn('Growth \\n   rate \\n(1/day)',
                         colours=or2.bl2.pal,
                         breaks=seq(-1.5,1.5,0.5),limits=c(-1.8,1.8),
                         guide = guide_colourbar(ticks.colour = gray(gl),
                                                 ticks.linewidth = 1.,
                                                 frame.colour = "black"))+
  theme_bw()
emp.ma.surf.plt

# combine panels
ma.surf.plt2<-ma.surf.plt+theme(legend.position = 'none')
grid.arrange(emp.ma.surf.plt,ma.surf.plt2,nrow=2,layout_matrix=rbind(c(1,1,NA),c(2,2,2)))
a1.ma<-arrangeGrob(emp.ma.surf.plt,ma.surf.plt2,nrow=2,layout_matrix=rbind(c(1,1,NA),c(2,2,2)))

ggsave(filename = '~/Dropbox/Phyto acclimation/accData/results/acclimation_surfaces/acclimation_surfaces_MA_P1_082719.pdf',a1.ma,width=8.25,height=8)

# Overlay empirical obs on estimated surfaces
emp.gr.ma$value<-NA

a2.ma<-ggplot(ma.surf,aes(x=Acute.temp,y=Acclim.temp,z=value))+
  geom_raster(aes(fill=value))+
  geom_abline(intercept=0,slope=1,linetype=2)+
  geom_contour(colour=gray(0.2),binwidth=0.1,size=0.1)+
  geom_contour(colour=gray(0.2),binwidth=0.5,size=0.4)+
  geom_contour(colour='red',breaks=c(0),size=1)+
  geom_point(data=emp.gr.ma,fill='black',colour='black',size=5.2,shape=21)+
  geom_point(data=emp.gr.ma,aes(fill=mu),colour='black',size=5,shape=21)+
  scale_x_continuous('Acute temperature (°C)',limits = c(15,42),expand=c(0,0))+
  scale_y_continuous('Acclimated temperature (°C)',expand=c(0,0))+
  scale_fill_gradientn('Growth \\n   rate \\n(1/day)',
                       colours=or2.bl2.pal,
                       breaks=seq(-1.5,1.5,0.5),limits=c(-1.8,1.8),
                       guide = guide_colourbar(ticks.colour = gray(gl),
                                               ticks.linewidth = 1.,
                                               frame.colour = "black"))+
  facet_wrap(~variable,nrow=2)+
  theme_bw()+
  ggtitle('M. aeruginosa')
a2.ma

#ggsave(filename = '~/Dropbox/Phyto acclimation/accData/results/acclimation_surfaces/acclimation_surfaces_wObs_MA_P1.pdf',a2.ma)

a2<-arrangeGrob(a2.cr,a2.ma,nrow=1)
ggsave(filename = '~/Dropbox/Phyto acclimation/accData/results/acclimation_surfaces/acclimation_surfaces_wObs_P1_082719.pdf',a2,width=8,height=6)


# visualize implied acclimated curve:

# temperature grid:
vals2<-data.frame(Acute.temp=seq(14,45,length.out = 100),Acclim.temp=seq(14,45,length.out = 100))

# predictions from GAM
pds.gam.acc<-predict(gm1.ma,newdata=vals2,se=T)
vals2<-data.frame(vals2,mu=pds.gam.acc$fit,se=pds.gam.acc$se.fit)
head(vals2)

predict(gm1.ma,newdata=data.frame(Acute.temp=35,Acclim.temp=35),se=T)
predict(gm1.ma,newdata=data.frame(Acute.temp=37,Acclim.temp=37),se=T)

plist.ma<-list()

plist.ma[[1]]<-ggplot(vals2,aes(x=Acclim.temp,y=mu))+
  geom_hline(yintercept = 0)+
  geom_ribbon(aes(ymin=mu-1.96*se,ymax=mu+1.96*se),fill='gray')+
  geom_line()+
  stat_smooth(method='gam',formula=y~s(x,k=4,bs='ps'))+
  geom_point(data=gr.ma.P1[gr.ma.P1$Acclim.temp==gr.ma.P1$Acute.temp,])+
  coord_cartesian(ylim = c(-2,1.5))+
  theme_bw()+
  ggtitle('k = 4')

plist.ma[[2]]<-ggplot(vals2,aes(x=Acclim.temp,y=mu))+
  geom_hline(yintercept = 0)+
  geom_ribbon(aes(ymin=mu-1.96*se,ymax=mu+1.96*se),fill='gray')+
  geom_line()+
  stat_smooth(method='gam',formula=y~s(x,k=5,bs='ps'))+
  geom_point(data=gr.ma.P1[gr.ma.P1$Acclim.temp==gr.ma.P1$Acute.temp,])+
  coord_cartesian(ylim = c(-2,1.5))+
  theme_bw()+
  ggtitle('k = 5')

plist.ma[[3]]<-ggplot(vals2,aes(x=Acclim.temp,y=mu))+
  geom_hline(yintercept = 0)+
  geom_ribbon(aes(ymin=mu-1.96*se,ymax=mu+1.96*se),fill='gray')+
  geom_line()+
  stat_smooth(method='gam',formula=y~s(x,k=6,bs='ps'))+
  geom_point(data=gr.ma.P1[gr.ma.P1$Acclim.temp==gr.ma.P1$Acute.temp,])+
  coord_cartesian(ylim = c(-2,1.5))+
  theme_bw()+
  ggtitle('k = 6')

plist.ma[[4]]<-ggplot(vals2,aes(x=Acclim.temp,y=mu))+
  geom_hline(yintercept = 0)+
  geom_ribbon(aes(ymin=mu-1.96*se,ymax=mu+1.96*se),fill='gray')+
  geom_line()+
  stat_smooth(method='gam',formula=y~s(x,k=7,bs='ps'))+
  geom_point(data=gr.ma.P1[gr.ma.P1$Acclim.temp==gr.ma.P1$Acute.temp,])+
  coord_cartesian(ylim = c(-2,1.5))+
  theme_bw()+
  ggtitle('k = 7')

grid.arrange(plist.ma[[1]],plist.ma[[2]],plist.ma[[3]],plist.ma[[4]],top='M. aeruginosa acclimated curve')



ggplot(vals2,aes(x=Acclim.temp,y=mu))+
  geom_hline(yintercept = 0)+
  geom_ribbon(aes(ymin=mu-1.96*se,ymax=mu+1.96*se),fill='gray')+
  geom_line()+
  stat_smooth(method='gam',formula=y~s(x,k=9))+
  geom_point(data=gr.cr.P1[gr.cr.P1$Acclim.temp==gr.cr.P1$Acute.temp,])+
  coord_cartesian(ylim = c(-1,1.5))+
  theme_bw()

cr.acc.gm<-gam(mu~s(Acclim.temp,k=7,bs='ps'),data=gr.cr.P1[gr.cr.P1$Acclim.temp==gr.cr.P1$Acute.temp,])
summary(cr.acc.gm)

cr.acc.gm<-gam(mu~s(Acclim.temp,k=6,bs='ps'),data=gr.cr.P1[gr.cr.P1$Acclim.temp==gr.cr.P1$Acute.temp,])
summary(cr.acc.gm)





#### (4) Save surface for downstream analyses ####

#write.csv(x = gr.ma.P1,"./data/Acc_surface/derived_data/acc_surface_MA_t4_GAM_subset_082719.csv",row.names=F)
write.csv(x = ma.surf,"./data/Acc_surface/derived_data/acc_surface_MA_t4_GAM_subset_092419.csv",row.names=F)



##### Combined CR and MA plot, GAMs only ####

## quick load data:
ma.surf<-read.csv("./data/Acc_surface/derived_data/acc_surface_MA_t4_GAM_subset_092419.csv")
cr.surf<-read.csv("./data/Acc_surface/derived_data/acc_surface_CR_t4_GAM_subset_092419.csv")

# pool data across species
ma.surf$Species<-'M. aeruginosa'
cr.surf$Species<-'C. reinhardtii'
sps.surf<-rbind(ma.surf,cr.surf)
emp.gr.ma$Species<-'M. aeruginosa'
emp.gr.cr$Species<-'C. reinhardtii'
emp.gr<-rbind(emp.gr.ma,emp.gr.cr)

# try out some plots
gl<-0.2

a2C<-sps.surf %>% filter(variable == 'gam') %>%
  ggplot(aes(x=Acute.temp,y=Acclim.temp,z=value))+
  geom_raster(aes(fill=value))+
  geom_abline(intercept=0,slope=1,linetype=2)+
  geom_contour(colour=gray(gl),binwidth=0.1,size=0.1)+
  geom_contour(colour=gray(gl),binwidth=0.5,size=0.4)+
  geom_contour(colour='black',breaks=c(0),size=1)+
  #geom_point(data=emp.gr,fill=gray(0.6),colour=gray(0.6),size=5.2,shape=21)+
  geom_point(data=emp.gr,aes(fill=mu),colour='black',size=5,shape=21)+
  scale_x_continuous('Acute temperature (°C)',expand=c(0,0))+
  scale_y_continuous('Acclimated temperature (°C)',expand=c(0,0))+
  #scale_fill_distiller('Growth \\n   rate \\n(1/day)',type='div',limits=c(-1.8,1.8),
  #                     palette = 'RdYlBu',breaks=seq(-1.5,1.5,0.5),
  #                     guide = guide_colourbar(ticks.colour = "black",
  #                                             ticks.linewidth = 1.,
  #                                             frame.colour = "black"))+
  scale_fill_gradientn('Growth \\n   rate \\n(1/day)',
                       colours=or2.bl2.pal,
                       breaks=seq(-1.5,1.5,0.5),limits=c(-1.8,1.8),
                       guide = guide_colourbar(ticks.colour = gray(gl),
                                               ticks.linewidth = 1.,
                                               frame.colour = "black"))+
  facet_wrap(~Species,nrow=1)+
  coord_cartesian(ylim=c(15,40),xlim=c(15,42))+
  theme_bw()+
  ggtitle('Acclimation surface GAMs')
a2C

ggsave('~/Dropbox/Phyto acclimation/accData/results/acclimation_surfaces/acclimation_surfaces_wObs_P1_gams_only_v4.pdf',a2C,width=8,height=4)


#### Separate z ranges:

gl<-0.0

a2C.cr<-sps.surf %>% filter(variable == 'gam' & Species=='C. reinhardtii') %>%
  ggplot(aes(x=Acute.temp,y=Acclim.temp,z=value))+
  geom_raster(aes(fill=value))+
  geom_abline(intercept=0,slope=1,linetype=2)+
  geom_contour(colour=gray(gl),binwidth=0.1,size=0.1)+
  geom_contour(colour=gray(gl),binwidth=0.5,size=0.4)+
  geom_contour(colour='black',breaks=c(0),size=1)+
  geom_point(data=emp.gr.cr,aes(fill=mu),colour='black',size=5,shape=21)+
  scale_x_continuous('Acute temperature (°C)',expand=c(0,0))+
  scale_y_continuous('Acclimated temperature (°C)',expand=c(0,0))+
  scale_fill_gradientn('Growth \\n   rate \\n(1/day)',
                       colours=or2.bl2.pal,
                       breaks=seq(-2,2,0.5),limits=c(-2.1,2.1),
                       guide = guide_colourbar(ticks.colour = gray(gl),
                                               ticks.linewidth = 1.,
                                               frame.colour = "black"))+
  facet_wrap(~Species,nrow=1)+
  coord_cartesian(ylim=c(15,39),xlim=c(15,45))+
  theme_bw()+
  ggtitle('Acclimation surface GAMs')
a2C.cr

a2C.ma<-sps.surf %>% filter(variable == 'gam' & Species=='M. aeruginosa') %>%
  ggplot(aes(x=Acute.temp,y=Acclim.temp,z=value))+
  geom_raster(aes(fill=value))+
  geom_abline(intercept=0,slope=1,linetype=2)+
  geom_contour(colour=gray(gl),binwidth=0.1,size=0.1)+
  geom_contour(colour=gray(gl),binwidth=0.5,size=0.4)+
  geom_contour(colour='black',breaks=c(0),size=1)+
  geom_point(data=emp.gr.ma,aes(fill=mu),colour='black',size=5,shape=21)+
  scale_x_continuous('Acute temperature (°C)',expand=c(0,0))+
  scale_y_continuous('Acclimated temperature (°C)',expand=c(0,0))+
  scale_fill_gradientn('Growth \\n   rate \\n(1/day)',
                       colours=or2.bl2.pal,
                       breaks=seq(-1.0,1.0,0.25),limits=c(-1.1,1.1),
                       guide = guide_colourbar(ticks.colour = gray(gl),
                                               ticks.linewidth = 1.,
                                               frame.colour = "black"))+
  facet_wrap(~Species,nrow=1)+
  coord_cartesian(ylim=c(15,39),xlim=c(15,42))+
  theme_bw()+
  ggtitle("    ")
a2C.ma

a2D<-arrangeGrob(a2C.cr,a2C.ma,nrow=1)
#grid.arrange(a2C.cr,a2C.ma)

ggsave('~/Dropbox/Phyto acclimation/accData/results/acclimation_surfaces/acclimation_surfaces_wObs_P1_gams_only_v5.pdf',a2D,width=10,height=4)

