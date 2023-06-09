rm(list=ls(all=TRUE))

library(car) 
library(ggplot2)
library(lme4)
library(plyr)
library(AICcmodavg)
library(ggplot2)
library(lattice)
library(doBy)
library(MASS)
library(pscl)
library(lmtest)
library(VGAM)  #for running truncated negbin models
  stats::resid #makes sure resid command from bases stats is used and not from VGAM (see p. 268 from Zuur 2009)
library(data.table)
library(vegan)
library(visreg)
library(nlme)
library(boot)
library(multcomp)

### start here ###
kpc15<-read.csv("KPC15data_EcolMono.csv")
  kpc15$BLOCK<-as.factor(kpc15$BLOCK)
  kpc15$GTYPE_1<-as.factor(kpc15$GTYPE_1)
      kpc15<-as.data.table(kpc15)

kpc15_exp<-kpc15[kpc15$EXP_MONO=="exp",]
kpc15_mono<-kpc15[kpc15$GENDIV=="1",][SEEDNUM %in% c("96","960")]  #create monos only dataset for comparing gtypes grown individually

#########################
# effectiveness of COMP treatments
summaryBy(BARECOVER13~COMP, FUN=c(mean,sd,length), data=kpc15_exp)
summaryBy(BARECOVER14~COMP, FUN=c(mean,sd,length), data=kpc15_exp)
summaryBy(BARECOVER15~COMP, FUN=c(mean,sd,length), data=kpc15_exp)
  
t.test(BARECOVER13~COMP, data=kpc15_exp)
t.test(BARECOVER14~COMP, data=kpc15_exp)
t.test(BARECOVER15~COMP, data=kpc15_exp)

#########################
# how did S_13 (gtype richness in year 1) respond to treatments?  lc=low comp, hc=high comp
s13_lc<-glm(S_13~SEEDNUM+GENDIV+PC2USE+I(SEEDNUM^2), family=poisson , data=subset(kpc15_exp,COMP=="low"))  
summary(s13_lc)  #sig effects of SEEDNUM, SEEDNUM^2 and GENDIV
    Anova(s13_lc,test.statistic="Wald",type="III")

s13_hc<-glm(S_13~SEEDNUM+GENDIV+PC2USE+I(SEEDNUM^2), family=poisson , data=subset(kpc15_exp,COMP=="control"))  
summary(s13_hc)  #marginally sig effects of SEEDNUM and SEEDNUM^2, plus GENDIV sig 
    Anova(s13_hc,test.statistic="Wald",type="III")

s13<-glm(S_13~COMP+SEEDNUM+GENDIV+PC2USE+I(SEEDNUM^2), family=poisson , data=kpc15_exp)  
summary(s13)  #sig effects of SEEDNUM, SEEDNUM^2 and GENDIV
    Anova(s13,test.statistic="Wald",type="III")
    

#how did S_15 respond to treatments?
s15_lc<-glm(S_15~SEEDNUM+GENDIV+PC2USE+I(SEEDNUM^2), family=poisson , data=subset(kpc15_exp,COMP=="low"))  
  summary(s15_lc)  # sig effects of SEEDNUM, SEEDNUM^2 and GD; marginally sig PC2
    Anova(s15_lc,test.statistic="Wald",type="III")

s15_hc<-glm(S_15~SEEDNUM+GENDIV+PC2USE, family=poisson , data=subset(kpc15_exp,COMP=="control"))  
  summary(s15_hc)  # sig effects of SEEDNUM only
    Anova(s15_hc,test.statistic="Wald",type="III")

s15<-glm(S_15~COMP+SEEDNUM+GENDIV+PC2USE+I(SEEDNUM^2), family=poisson , data=kpc15_exp)  
  summary(s15)  # 
   Anova(s15,test.statistic="Wald",type="III")


########################
##calculate genotype retention rates (S_13 or S_15 relative to GENDIV)
##  as an alternative to realized richness
kpc15_exp$S_13ret<-kpc15_exp$S_13/kpc15_exp$GENDIV
kpc15_exp$S_15ret<-kpc15_exp$S_15/kpc15_exp$GENDIV

#ret13_l<-glm(S_13ret~SEEDNUM+GENDIV+PC2USE+I(SEEDNUM^2), data=kpc15_exp[COMP=="low"])
#nonlinear regression based on michaelis-menton: Cmax=asymptote, K=value of SEEDNUM for half-asympote, xshift is amount of shift left to right for entire curve
  #parameters PC and GD allow for the curve to shift up or down in response to PC2USE and GENDIV
ret13_l_nlin<-nls(S_13ret~(GD*GENDIV)+(PC*PC2USE)+(Cmax*(SEEDNUM-xshift))/(K+SEEDNUM-xshift), data=kpc15_exp[COMP=="low"], 
                      start=list(GD=.05,PC=.05, Cmax=.15, K=10, xshift=30 ))
  summary(ret13_l_nlin)
    Anova(ret13_l_nlin,test.statistic="Wald",type="III")

#ret13_c<-glm(S_13ret~SEEDNUM+GENDIV+PC2USE+I(SEEDNUM^2), data=kpc15_exp[COMP=="control"])
#nonlinear regression based on michaelis-menton: Cmax=asymptote, K=value of SEEDNUM for half-asympote, xshift is amount of shift left to right for entire curve
  #parameters PC and GD allow for the curve to shift up or down in response to PC2USE and GENDIV
ret13_c_nlin<-nls(S_13ret~(GD*GENDIV)+(PC*PC2USE)+(Cmax*(SEEDNUM-xshift))/(K+SEEDNUM-xshift), data=kpc15_exp[COMP=="control"], 
                      start=list(GD=.05,PC=.05, Cmax=.15, K=10, xshift=30 ))
  summary(ret13_c_nlin)
    Anova(ret13_c_nlin,test.statistic="Wald",type="III")


ret15_l<-glm(S_15ret~SEEDNUM+GENDIV+PC2USE+I(SEEDNUM^2), data=kpc15_exp[COMP=="low"])  #note that this isn't driven by GD1 plots
 # AIC(ret15_l) #88.9
#nonlinear regression based on michaelis-menton: Cmax=asymptote, K=value of SEEDNUM for half-asympote, xshift is amount of shift left to right for entire curve
  #parameters PC and GD allow for the curve to shift up or down in response to PC2USE and GENDIV
ret15_l_nlin<-nls(S_15ret~(GD*GENDIV)+(PC*PC2USE)+(Cmax*(SEEDNUM-xshift))/(K+SEEDNUM-xshift), data=kpc15_exp[COMP=="low"], 
                      start=list(GD=.05,PC=.05, Cmax=.15, K=10, xshift=30 ))
#    AIC(ret15_l_nlin) #91.72 - poorer fit  (perhaps indicating that negative density dependence drove greater loss of gtypes)
  summary(ret15_l)
    Anova(ret15_l,test.statistic="Wald",type="III")

ret15_c<-glm(S_15ret~SEEDNUM+GENDIV+PC2USE, data=kpc15_exp[COMP=="control"])
  summary(ret15_c)
    Anova(ret15_c,test.statistic="Wald",type="III")


    
ret13<-glm(S_13ret~COMP+SEEDNUM+GENDIV+PC2USE+I(SEEDNUM^2), data=kpc15_exp)
  summary(ret13)
    Anova(ret13,test.statistic="Wald",type="III")
ret15<-glm(S_15ret~COMP+SEEDNUM+GENDIV+PC2USE+I(SEEDNUM^2), data=kpc15_exp)
  summary(ret15)
    Anova(ret15,test.statistic="Wald",type="III")


#########################
### 2015 abundance, with PC2 for nutrients as covariate  ###
## analysis with all plots gets to a model with many sig interactions, including COMP:SEEDNUM:PC2, COMP:GENDIV:PC2, COMP:SEEDNUM and COMP:GENDIV
## we take this as evidence to support running analyses separately by COMP treatment

## Low COMP plots only ##

#binomial models for presence/absence 
pc2_lowc_pres<-glm(PRES15~SEEDNUM*GENDIV*PC2USE, family=binomial, data=subset(kpc15_exp,COMP=="low"))
  pc2_lowc_pres2<-update(pc2_lowc_pres,.~.-SEEDNUM:GENDIV:PC2USE)
  pc2_lowc_pres3<-update(pc2_lowc_pres2,.~.-SEEDNUM:GENDIV)
  pc2_lowc_pres4<-update(pc2_lowc_pres3,.~.-GENDIV:PC2USE) #
  pc2_lowc_pres4<-glm(PRES15 ~ SEEDNUM + GENDIV + PC2USE + SEEDNUM:PC2USE,
                      family=binomial, data=subset(kpc15_exp,COMP=="low"))
  pc2_lowc_pres5<-update(pc2_lowc_pres4,.~.+I(SEEDNUM^2)) #
  pc2_lowc_pres6<-update(pc2_lowc_pres5,.~.-SEEDNUM:PC2USE) #final model
drop1(pc2_lowc_pres5, test="LRT")  
  lrtest(pc2_lowc_pres5,pc2_lowc_pres4)
  AIC(pc2_lowc_pres4, pc2_lowc_pres6)
summary(pc2_lowc_pres6)
  Anova(pc2_lowc_pres6, test.statistic = "Wald", type= "III")
  plot(pc2_lowc_pres4)

  visreg2d(pc2_lowc_pres4,"SEEDNUM" ,"PC2USE",plot.type="image")  #shows 2-way interaction b/w SEEDNUM and PC2
  visreg(pc2_lowc_pres4,"GENDIV",scale="response",ylab="Pr(survival to year 3)")  #shows positive effect of GENDIV
    visreg(pc2_lowc_pres4,"GENDIV",ylab="log odds(survival to year 3)")  #shows positive effect of GENDIV
  visreg(pc2_lowc_pres6,"SEEDNUM", scale="response" , ylab="Pr(survival to year 3)")  #

# binomial models excluding the highest nutrient plots
quantile(kpc15_exp[COMP=="low"]$PC2USE,c(.25,.5,.75))
pc2_low75<-glm(PRES15~SEEDNUM+GENDIV+PC2USE+I(SEEDNUM^2), family=binomial, data=subset(kpc15_exp,COMP=="low" & PC2USE <0.32940323))
  summary(pc2_low75) #lower 75th percentile of PC2 (GD p=.015, SN p=.001)
## 

#truncated negbin models; default link is log
pc2_lowc_trnb<-vglm(FL_SUM15~SEEDNUM*GENDIV*PC2USE, family=posnegbinomial,control=vglm.control(maxit=100), data=subset(kpc15_exp,PRES15==1 & COMP=="low"))  #use VGAM::vglm to run truncated NB GLM (see p. 268 of Zuur 2009)
  pc2_lowc_trnb2<-update(pc2_lowc_trnb,.~.-SEEDNUM:GENDIV:PC2USE)
  pc2_lowc_trnb3<-update(pc2_lowc_trnb2,.~.-SEEDNUM:PC2USE)
  pc2_lowc_trnb4<-update(pc2_lowc_trnb3,.~.-GENDIV:PC2USE)
  pc2_lowc_trnb5<-update(pc2_lowc_trnb4,.~.-SEEDNUM:GENDIV) #final model
  pc2_lowc_trnb5<-vglm(FL_SUM15~SEEDNUM+GENDIV+PC2USE, family=posnegbinomial,control=vglm.control(maxit=100), data=subset(kpc15_exp,PRES15==1 & COMP=="low" ))  #use VGAM::vglm to run truncated NB GLM (see p. 268 of Zuur 2009)
#  pc2_lowc_trnb5b<-vglm(FL_SUM15~SEEDNUM+GENDIV+PC2USE, family=posnegbinomial,control=vglm.control(maxit=100), data=subset(kpc15_exp,PRES15==1 & COMP=="low" & FL_SUM15<200))  #run without 2x 1GD outliers
summary(pc2_lowc_trnb5)  #note that the sig negative effect of GD goes away when 2 1GD plots with FL_SUM15>200 are dropped
    Anova(pc2_lowc_trnb5, test.statistic = "Chisq", type="III")
  plot(pc2_lowc_trnb5b)
  plot(pc2_lowc_trnb5,form=resid(., type="pearson")~fitted(.),abline=0)
  hatplot(pc2_lowc_trnb5b)




#####  
#####  
## High COMP plots only ##
#truncated negative binomial model
pc2_hc_trnb<-vglm(FL_SUM15~SEEDNUM+GENDIV+PC2USE+SEEDNUM:GENDIV+SEEDNUM:PC2USE+GENDIV:PC2USE+SEEDNUM:GENDIV:PC2USE, family=pospoisson,control=vglm.control(maxit=100), data=subset(kpc15_exp,PRES15==1 & COMP=="control"))  #use VGAM::vglm to run truncated NB GLM (see p. 268 of Zuur 2009)
  pc2_hc_trnb2<-vglm(FL_SUM15~SEEDNUM+GENDIV+PC2USE+SEEDNUM:GENDIV+SEEDNUM:PC2USE+GENDIV:PC2USE, family=pospoisson,control=vglm.control(maxit=100), data=subset(kpc15_exp,PRES15==1 & COMP=="control"))  #use VGAM::vglm to run truncated NB GLM (see p. 268 of Zuur 2009)
  pc2_hc_trnb3<-update(pc2_hc_trnb2,.~.-GENDIV:PC2USE)
  pc2_hc_trnb4<-update(pc2_hc_trnb3,.~.-SEEDNUM:GENDIV)
  pc2_hc_trnb5<-update(pc2_hc_trnb4,.~.-SEEDNUM:PC2USE)
  pc2_hc_trnb5<-vglm(FL_SUM15~SEEDNUM+GENDIV+PC2USE,family=posnegbinomial,control=vglm.control(maxit=100), data=subset(kpc15_exp,PRES15==1 & COMP=="control"  & FL_SUM15<30))
  pc2_hc_trnb5b<-vglm(FL_SUM15~SEEDNUM+GENDIV+PC2USE,family=pospoisson,control=vglm.control(maxit=100), data=subset(kpc15_exp,PRES15==1 & COMP=="control" & FL_SUM15<30 ))# & FL_SUM15<30))
summary(pc2_hc_trnb5) #note that one 8GD outlier with 33 individuals is omitted here to make it run without convergence errors (but the hurdle command above appears to be less finicky)
  AIC(pc2_hc_trnb5b)
plot(pc2_hc_trnb5,form=resid(., type="pearson")~fitted(.),abline=0)
  hatplot(pc2_hc_trnb5)
  
#binomial presence-absence model
pc2_hc_pres<-glm(PRES15~SEEDNUM*GENDIV*PC2USE, family=binomial, data=subset(kpc15_exp,COMP=="control"))
  pc2_hc_pres2<-update(pc2_hc_pres,.~. -SEEDNUM:GENDIV:PC2USE)
  pc2_hc_pres3<-update(pc2_hc_pres2,.~. -SEEDNUM:GENDIV)
  pc2_hc_pres4<-update(pc2_hc_pres3,.~. -SEEDNUM:PC2USE)
  pc2_hc_pres5<-update(pc2_hc_pres4,.~. -GENDIV:PC2USE) #final model, all interactions dropped
  pc2_hc_pres5<-glm(PRES15~SEEDNUM+GENDIV+PC2USE,family=binomial,data=subset(kpc15_exp,COMP=="control"))
drop1(pc2_hc_pres5, test="LRT")
  summary(pc2_hc_pres5)
    Anova(pc2_hc_pres5, test.statistic = "Wald", type="III")

###############################
###############################
###############################
# repeated measures analyses on abundances over time FL_SUM13, 14 and 15

#  first, create a simplified dataset with only the responses of interest
# and then convert from wide to long form with reshape()
kpc_rm_wide<-kpc15_exp[,.(V1,COMP,SEEDNUM,GENDIV,PC2USE,FL_SUM13,FL_SUM14,FL_SUM15)] #pull out a subset of variables
kpc_rm<-reshape(kpc_rm_wide,varying=c("FL_SUM13","FL_SUM14","FL_SUM15"),idvar="V1",direction="long",v.names="FL_SUM",sep="")
  kpc_rm$time<-as.integer(kpc_rm$time)
  kpc_rm$V1<-as.factor(kpc_rm$V1)
  
#then, run rm analysis
  #first, on low-comp plots....

#use this approach, with FL_SUM sqrt transformed  
flrm1<-lme(sqrt(FL_SUM)~SEEDNUM*GENDIV*PC2USE*time, random=~1|V1,method="ML", correlation=corAR1(form = ~time|V1,value=.640),data=subset(kpc_rm,COMP=="low")) 
  flrm2<-update(flrm1,.~. -SEEDNUM:GENDIV:PC2USE:time)
  flrm3<-update(flrm2,.~. -SEEDNUM:GENDIV:time)
  flrm4<-update(flrm3,.~. -SEEDNUM:PC2USE:time)
  flrm5<-update(flrm4,.~. -SEEDNUM:GENDIV:PC2USE)
  flrm6<-update(flrm5,.~. -SEEDNUM:GENDIV)  #best model
  flrm6<-lme(sqrt(FL_SUM)~SEEDNUM+GENDIV+PC2USE+time+SEEDNUM:PC2USE+GENDIV:PC2USE+SEEDNUM:time+GENDIV:time+PC2USE:time+GENDIV:PC2USE:time,random=~1|V1, method="ML", correlation=corAR1(form= ~time|V1, value=.64), data=subset(kpc_rm,COMP=="low" ))
#  flrm7<-update(flrm6,.~. -GENDIV:PC2USE:time)  # not supported by AIC
  flrm7<-update(flrm6,.~. +I(SEEDNUM^2))  #much better fit (22 AIC units lower); interaction with time not supported
    ACF(flrm1)
    summary(flrm7)  #results are robust to inclusion of plots 314 and 311 (outliers in 2015 with GD1 and AT04)
    drop1(flrm6,test="Chisq")
    Anova(flrm7,type="III",test.statistic="Wald")
      AIC(flrm6,flrm7)
    vif(flrm6)
    hist(resid(flrm6))

###
#then on high-comp plots...
flrm_hc1<-lme(sqrt(FL_SUM)~SEEDNUM*GENDIV*PC2USE*time, random=~1|V1,method="ML", correlation=corAR1(form = ~time|V1, value=-.5),data=subset(kpc_rm,COMP=="control")) 
  flrm_hc2<-update(flrm_hc1,.~. -SEEDNUM:GENDIV:PC2USE:time)
  flrm_hc3<-update(flrm_hc2,.~. -SEEDNUM:GENDIV:time)
  flrm_hc4<-update(flrm_hc3,.~. -SEEDNUM:GENDIV:PC2USE)
  flrm_hc5<-update(flrm_hc4,.~. -SEEDNUM:GENDIV)
  flrm_hc6<-update(flrm_hc5,.~. -GENDIV:PC2USE:time)  
  flrm_hc7<-update(flrm_hc6,.~. -GENDIV:time)  
  flrm_hc8<-update(flrm_hc7,.~. -GENDIV:PC2USE)  # nearly the best model, but AIC and parsimony support continuing to simplify...
  flrm_hc9<-update(flrm_hc8,.~. -SEEDNUM:PC2USE:time)  
  flrm_hc10<-update(flrm_hc9,.~. -PC2USE:time)  #final model
  flrm_hc10<-lme(sqrt(FL_SUM)~SEEDNUM+GENDIV+PC2USE+time+SEEDNUM:PC2USE+SEEDNUM:time, random=~1|V1,method="ML", correlation=corAR1(form = ~time|V1, value=-.5),data=subset(kpc_rm,COMP=="control")) 
  flrm_hc11<-update(flrm_hc10,.~. +I(SEEDNUM^2))  #
#  flrm_hc12<-update(flrm_hc11,.~. +I(SEEDNUM^2):time)  #
    ACF(flrm_hc1)
    drop1(flrm_hc10,test="Chisq")
    summary(flrm_hc11)
    Anova(flrm_hc11,type="III",test.statistic="F")
    vif(flrm_hc10)
    qqnorm(resid(flrm_hc10))
    plot(flrm_hc10)
    AIC(flrm_hc12,flrm_hc11)

###
#and just to verify that it's reasonable to split the analyses by comp
flrm_all<-lme(sqrt(FL_SUM)~COMP*SEEDNUM*GENDIV*PC2USE*time, random=~1|V1,method="ML", correlation=corAR1(form = ~time|V1),data=kpc_rm) 
  flrm_all1<-update(flrm_all, .~. -COMP:SEEDNUM:GENDIV:PC2USE:time)
  summary(flrm_all1)
  Anova(flrm_all1,type="III",test.statistic="F")
  
###############################
###############################
###############################
##   testing for 'good' versus 'poor' genotypes  
#  getting at gtype diffs using the monoculture-only data (96 & 960 seed plots, both comp treatments) in a similar way as above
mono_abund13<-glmer((FL_SUM13+.005)~GTYPE_1 + PC2USE + (1|SEEDNUM) + (1|COMP), family=Gamma(link="log"),data=kpc15_mono)  #accounts for PC2, not much change...
  summary(mono_abund13) 
  plot(mono_abund13)
    mono_ab13_mcp<-glht(mono_abund13, linfct=mcp(GTYPE_1="Tukey"))
    summary(mono_ab13_mcp) #also takes a long time to run

mono_abund15<-glmer((FL_SUM15+.005)~GTYPE_1 + PC2USE + (1|SEEDNUM) + (1|COMP), family=Gamma(link="log"),data=kpc15_mono)  #accounts for PC2, not much change...
  summary(mono_abund15) 
  plot(mono_abund15)
    mono_ab15_mcp<-glht(mono_abund15, linfct=mcp(GTYPE_1="Tukey"))
    summary(mono_ab15_mcp) #also takes a long time to run
    
###############################
###############################
###############################

#  comparing GD 4 or GD8 means with expectations based on monocultures 
###############################
###############################
###############################

  #0. First, create the necessary datasets (to be used for newdata in step 2)
kpc_indiv<-read.csv("KPC15data_indiv_EcolMono.csv")
  kpc_indiv<-as.data.table(kpc_indiv)
  setkey(kpc_indiv,PLOT)
  kpc_indiv$GTYPE_1<-as.factor(kpc_indiv$GTYPE_1)

  
#1. using all monocultures (kpc15_mono), run model on focal response (e.g., abund15, pres15) to get a predicted value for each genotype 
    #note this is done only for low COMP  (so 150 plots total)
  #yr 3 abundance, excluding zeroes, low comp
a15_no0_low<-vglm(FL_SUM15~PC2USE+SEEDNUM+GTYPE_1,family=pospoisson,data=subset(kpc15_mono,COMP=="low" & PRES15==1 )) #add FL_SUM15<200 as another condition for subsetting to plot exp minus outliers
  #yr 3 persistence, low comp
pres15_low<-glm(PRES15~PC2USE+SEEDNUM+GTYPE_1,family=binomial,data=subset(kpc15_mono,COMP=="low"))
  #yr 3 persistence, high comp
pres15_high<-glm(PRES15~PC2USE+SEEDNUM+GTYPE_1,family=binomial,data=subset(kpc15_mono,COMP=="control"))
#2. calculate the per-genotype expectations for all genotypes grown alone, extrapolating across all 5 levels of SEEDNUM
    #use the predict function (or predict.glm), predicting values for a 'newdata' dataset (with additional steps for vglm since AT01 isn't included in those models)
    #newdata needs to be a dataframe with one row per gtype & plot, with PC2USE, GTYPE (as GTYPE_1) and SEEDNUM (as SEEDNUM/GENDIV to give per gtype SEEDNUM)
        ##NOTE that predict is done on per-plot SEEDNUM, aggregating by averaging for abundance 
        ##  but for persistence, predict is based on per-gtype SEEDNUM and then aggregated by inclusion-exclusion
    #first for abundances (minus zeroes)
kpc_indiv_vglm<-subset(kpc_indiv,COMP=="low" & GTYPE_1!=1) #create subset without GTYPE AT01 b/c its absent from vglm analysis (no surviving pops)
  kpc_indiv_vglm$a15_no0_exp<-predictvglm(a15_no0_low,newdata=kpc_indiv_vglm,type="response") #create per-gtype predictions in new new df
  setkey(kpc_indiv_vglm,PLOT,GTYPE) #set key in the subset datatable to make the next line work correctly
  setDT(kpc_indiv,key=c("PLOT","GTYPE"))[kpc_indiv_vglm,a15_no0_exp := i.a15_no0_exp] #left join to kpc_indiv, keeping only the vglm predictions from the smaller df
  kpc_indiv[GTYPE_1==1 & COMP=="low",a15_no0_exp := 0]  #define the vglm prediction as 0 (replaces NA) for gtype AT01 to make the expected values by plot accurate
  
    #then for persistence: first, create a dataset with per-genotype seed numbers as SEEDNUM (so predictions are based on per-genotype seed numbers)
  kpc_indiv_copy<-kpc_indiv
    kpc_indiv_copy$SEEDNUM_PLOT<-kpc_indiv_copy$SEEDNUM #then create a new SEEDNUM variable specifically at plot scale
    kpc_indiv_copy$SEEDNUM<-kpc_indiv_copy$SEEDNUM_PLOT/kpc_indiv_copy$GENDIV #re-scale SEEDNUM for this dataset to be per-genotype in each plot
kpc_indiv_copy$pres15_exp<-ifelse(kpc_indiv_copy$COMP=="low",predict.glm(pres15_low,newdata=kpc_indiv_copy, type="response"),NA)
  setkey(kpc_indiv_copy,PLOT,GTYPE) #set key in the copy datatable to make the next line work correctly
  kpc_indiv[kpc_indiv_copy,pres15_exp :=i.pres15_exp]

kpc_indiv_copy$pres15_high_exp<-ifelse(kpc_indiv_copy$COMP=="control",predict.glm(pres15_high,newdata=kpc_indiv_copy, type="response"),NA)
  kpc_indiv[kpc_indiv_copy,pres15_high_exp :=i.pres15_high_exp]
     
  #3. calculate expectations for mixture plots from monoculture predictions (will also be done for 1GD plots from exp, as a check)
  # and add those predicted values to kpc15_exp
      #aggregating abundance values
expvals<-summaryBy(a15_no0_exp~PLOT, data=kpc_indiv, FUN=mean)
  setkey(expvals,PLOT)
  setkey(kpc15_exp,V1)
kpc15_exp<-kpc15_exp[expvals]
      
      #aggregating predicted persistence values - LOW COMP by default, but then re-run this code for CONTROL comp by selecting alt rows to run
gd1<-subset(kpc_indiv,GENDIV==1)
  gd1[,c("COMP","SEEDNUM","GENDIV","GTYPE","ABUND_13","ABUND_15","PC2USE","GTYPE_1","a15_no0_exp","pres15_high_exp"):=NULL] 
gd4<-subset(kpc_indiv,GENDIV==4)  
  gd4[,c("COMP","SEEDNUM","GENDIV","GTYPE","ABUND_13","ABUND_15","PC2USE","GTYPE_1","a15_no0_exp","pres15_high_exp"):=NULL] 
  gd4$G<-rep(c("G1","G2","G3","G4"),90)
  gd4<-reshape(gd4,idvar="PLOT", timevar="G",direction="wide")
  setnames(gd4,c("pres15_exp.G1","pres15_exp.G2","pres15_exp.G3","pres15_exp.G4"),c("G1","G2","G3","G4"))
  gd4$pres15_exp<-ifelse(is.na(gd4$G1),NA,
                  (gd4$G1+gd4$G2+gd4$G3+gd4$G4
                  -(gd4$G1*gd4$G2)-(gd4$G1*gd4$G3)-(gd4$G1*gd4$G4)-(gd4$G2*gd4$G3)-(gd4$G2*gd4$G4)-(gd4$G3*gd4$G4)
                  +(gd4$G1*gd4$G2*gd4$G3)+(gd4$G1*gd4$G2*gd4$G4)+(gd4$G1*gd4$G3*gd4$G4)+(gd4$G2*gd4$G3*gd4$G4)
                  -(gd4$G1*gd4$G2*gd4$G3*gd4$G4)))  #inclusion-exclusion works!!
  gd4[,c("G1","G2","G3","G4") :=NULL]
gd8<-subset(kpc_indiv,GENDIV==8)  
  gd8[,c("COMP","SEEDNUM","GENDIV","GTYPE","ABUND_13","ABUND_15","PC2USE","GTYPE_1","a15_no0_exp","pres15_high_exp"):=NULL] 
  gd8$G<-rep(c("G1","G2","G3","G4","G5","G6","G7","G8"),90)
  gd8<-reshape(gd8,idvar="PLOT", timevar="G",direction="wide")
  setnames(gd8,c("pres15_exp.G1","pres15_exp.G2","pres15_exp.G3","pres15_exp.G4","pres15_exp.G5","pres15_exp.G6","pres15_exp.G7","pres15_exp.G8"),c("G1","G2","G3","G4","G5","G6","G7","G8"))
  #next long batch of code aggregates persistence probabilities for 8gd treatments using inclusion-exclusion principle
  gd8$pres15_exp<-ifelse(is.na(gd8$G1),NA,
                  (gd8$G1+gd8$G2+gd8$G3+gd8$G4+gd8$G5+gd8$G6+gd8$G7+gd8$G8 
                  -(gd8$G1*gd8$G2)-(gd8$G1*gd8$G3)-(gd8$G1*gd8$G4)-(gd8$G1*gd8$G5)-(gd8$G1*gd8$G6)-(gd8$G1*gd8$G7)-(gd8$G1*gd8$G8)
                      -(gd8$G2*gd8$G3)-(gd8$G2*gd8$G4)-(gd8$G2*gd8$G5)-(gd8$G2*gd8$G6)-(gd8$G2*gd8$G7)-(gd8$G2*gd8$G8)
                      -(gd8$G3*gd8$G4)-(gd8$G3*gd8$G5)-(gd8$G3*gd8$G6)-(gd8$G3*gd8$G7)-(gd8$G3*gd8$G8)
                      -(gd8$G4*gd8$G5)-(gd8$G4*gd8$G6)-(gd8$G4*gd8$G7)-(gd8$G4*gd8$G8)
                      -(gd8$G5*gd8$G6)-(gd8$G5*gd8$G7)-(gd8$G5*gd8$G8)
                      -(gd8$G6*gd8$G7)-(gd8$G6*gd8$G8)
                      -(gd8$G7*gd8$G8)  
                  +(gd8$G1*gd8$G2*gd8$G3)+(gd8$G1*gd8$G2*gd8$G4)+(gd8$G1*gd8$G2*gd8$G5)+(gd8$G1*gd8$G2*gd8$G6)+(gd8$G1*gd8$G2*gd8$G7)+(gd8$G1*gd8$G2*gd8$G8)
                    +(gd8$G1*gd8$G3*gd8$G4)+(gd8$G1*gd8$G3*gd8$G5)+(gd8$G1*gd8$G3*gd8$G6)+(gd8$G1*gd8$G3*gd8$G7)+(gd8$G1*gd8$G3*gd8$G8)
                    +(gd8$G1*gd8$G4*gd8$G5)+(gd8$G1*gd8$G4*gd8$G6)+(gd8$G1*gd8$G4*gd8$G7)+(gd8$G1*gd8$G4*gd8$G8)
                    +(gd8$G1*gd8$G5*gd8$G6)+(gd8$G1*gd8$G5*gd8$G7)+(gd8$G1*gd8$G5*gd8$G8)
                    +(gd8$G1*gd8$G6*gd8$G7)+(gd8$G1*gd8$G6*gd8$G8)+(gd8$G1*gd8$G7*gd8$G8)
                    +(gd8$G2*gd8$G3*gd8$G4)+(gd8$G2*gd8$G3*gd8$G5)+(gd8$G2*gd8$G3*gd8$G6)+(gd8$G2*gd8$G3*gd8$G7)+(gd8$G2*gd8$G3*gd8$G8)
                    +(gd8$G2*gd8$G4*gd8$G5)+(gd8$G2*gd8$G4*gd8$G6)+(gd8$G2*gd8$G4*gd8$G7)+(gd8$G2*gd8$G4*gd8$G8)
                    +(gd8$G2*gd8$G5*gd8$G6)+(gd8$G2*gd8$G5*gd8$G7)+(gd8$G2*gd8$G5*gd8$G8)
                    +(gd8$G2*gd8$G6*gd8$G7)+(gd8$G2*gd8$G6*gd8$G8)+(gd8$G2*gd8$G7*gd8$G8)
                    +(gd8$G3*gd8$G4*gd8$G5)+(gd8$G3*gd8$G4*gd8$G6)+(gd8$G3*gd8$G4*gd8$G7)+(gd8$G3*gd8$G4*gd8$G8)
                    +(gd8$G3*gd8$G5*gd8$G6)+(gd8$G3*gd8$G5*gd8$G7)+(gd8$G3*gd8$G5*gd8$G8)
                    +(gd8$G3*gd8$G6*gd8$G7)+(gd8$G3*gd8$G6*gd8$G8)+(gd8$G3*gd8$G7*gd8$G8)
                    +(gd8$G4*gd8$G5*gd8$G6)+(gd8$G4*gd8$G5*gd8$G7)+(gd8$G4*gd8$G5*gd8$G8)
                    +(gd8$G4*gd8$G6*gd8$G7)+(gd8$G4*gd8$G6*gd8$G8)+(gd8$G4*gd8$G7*gd8$G8)
                    +(gd8$G5*gd8$G6*gd8$G7)+(gd8$G5*gd8$G6*gd8$G8)+(gd8$G5*gd8$G7*gd8$G8)
                    +(gd8$G6*gd8$G7*gd8$G8) 
                  -(gd8$G1*gd8$G2*gd8$G3*gd8$G4)-(gd8$G1*gd8$G2*gd8$G3*gd8$G5)-(gd8$G1*gd8$G2*gd8$G3*gd8$G6)-(gd8$G1*gd8$G2*gd8$G3*gd8$G7)-(gd8$G1*gd8$G2*gd8$G3*gd8$G8)
                    -(gd8$G1*gd8$G2*gd8$G4*gd8$G5)-(gd8$G1*gd8$G2*gd8$G4*gd8$G6)-(gd8$G1*gd8$G2*gd8$G4*gd8$G7)-(gd8$G1*gd8$G2*gd8$G4*gd8$G8)
                    -(gd8$G1*gd8$G2*gd8$G5*gd8$G6)-(gd8$G1*gd8$G2*gd8$G5*gd8$G7)-(gd8$G1*gd8$G2*gd8$G5*gd8$G8)
					          -(gd8$G1*gd8$G2*gd8$G6*gd8$G7)-(gd8$G1*gd8$G2*gd8$G6*gd8$G8)-(gd8$G1*gd8$G2*gd8$G7*gd8$G8)
                    -(gd8$G1*gd8$G3*gd8$G4*gd8$G5)-(gd8$G1*gd8$G3*gd8$G4*gd8$G6)-(gd8$G1*gd8$G3*gd8$G4*gd8$G7)-(gd8$G1*gd8$G3*gd8$G4*gd8$G8)
                    -(gd8$G1*gd8$G3*gd8$G5*gd8$G6)-(gd8$G1*gd8$G3*gd8$G5*gd8$G7)-(gd8$G1*gd8$G3*gd8$G5*gd8$G8)
				          	-(gd8$G1*gd8$G3*gd8$G6*gd8$G7)-(gd8$G1*gd8$G3*gd8$G6*gd8$G8)-(gd8$G1*gd8$G3*gd8$G7*gd8$G8)
                    -(gd8$G1*gd8$G4*gd8$G5*gd8$G6)-(gd8$G1*gd8$G4*gd8$G5*gd8$G7)-(gd8$G1*gd8$G4*gd8$G5*gd8$G8)
					          -(gd8$G1*gd8$G4*gd8$G6*gd8$G7)-(gd8$G1*gd8$G4*gd8$G6*gd8$G8)-(gd8$G1*gd8$G4*gd8$G7*gd8$G8)
					          -(gd8$G1*gd8$G5*gd8$G6*gd8$G7)-(gd8$G1*gd8$G5*gd8$G6*gd8$G8)-(gd8$G1*gd8$G5*gd8$G7*gd8$G8)
					          -(gd8$G1*gd8$G6*gd8$G7*gd8$G8)
                    -(gd8$G2*gd8$G3*gd8$G4*gd8$G5)-(gd8$G2*gd8$G3*gd8$G4*gd8$G6)-(gd8$G2*gd8$G3*gd8$G4*gd8$G7)-(gd8$G2*gd8$G3*gd8$G4*gd8$G8)
                    -(gd8$G2*gd8$G3*gd8$G5*gd8$G6)-(gd8$G2*gd8$G3*gd8$G5*gd8$G7)-(gd8$G2*gd8$G3*gd8$G5*gd8$G8)
					          -(gd8$G2*gd8$G3*gd8$G6*gd8$G7)-(gd8$G2*gd8$G3*gd8$G6*gd8$G8)-(gd8$G2*gd8$G3*gd8$G7*gd8$G8)
                    -(gd8$G2*gd8$G4*gd8$G5*gd8$G6)-(gd8$G2*gd8$G4*gd8$G5*gd8$G7)-(gd8$G2*gd8$G4*gd8$G5*gd8$G8)
					          -(gd8$G2*gd8$G4*gd8$G6*gd8$G7)-(gd8$G2*gd8$G4*gd8$G6*gd8$G8)-(gd8$G2*gd8$G4*gd8$G7*gd8$G8)
					          -(gd8$G2*gd8$G5*gd8$G6*gd8$G7)-(gd8$G2*gd8$G5*gd8$G6*gd8$G8)-(gd8$G2*gd8$G5*gd8$G7*gd8$G8)
					          -(gd8$G2*gd8$G6*gd8$G7*gd8$G8)
                    -(gd8$G3*gd8$G4*gd8$G5*gd8$G6)-(gd8$G3*gd8$G4*gd8$G5*gd8$G7)-(gd8$G3*gd8$G4*gd8$G5*gd8$G8)
					          -(gd8$G3*gd8$G4*gd8$G6*gd8$G7)-(gd8$G3*gd8$G4*gd8$G6*gd8$G8)-(gd8$G3*gd8$G4*gd8$G7*gd8$G8)
					          -(gd8$G3*gd8$G5*gd8$G6*gd8$G7)-(gd8$G3*gd8$G5*gd8$G6*gd8$G8)-(gd8$G3*gd8$G5*gd8$G7*gd8$G8)
					          -(gd8$G3*gd8$G6*gd8$G7*gd8$G8)					
                    -(gd8$G4*gd8$G5*gd8$G6*gd8$G7)-(gd8$G4*gd8$G5*gd8$G6*gd8$G8)-(gd8$G4*gd8$G5*gd8$G7*gd8$G8)-(gd8$G4*gd8$G6*gd8$G7*gd8$G8)
                    -(gd8$G5*gd8$G6*gd8$G7*gd8$G8) 
                  +(gd8$G1*gd8$G2*gd8$G3*gd8$G4*gd8$G5)+(gd8$G1*gd8$G2*gd8$G3*gd8$G4*gd8$G6)+(gd8$G1*gd8$G2*gd8$G3*gd8$G4*gd8$G7)+(gd8$G1*gd8$G2*gd8$G3*gd8$G4*gd8$G8)
					          +(gd8$G1*gd8$G2*gd8$G3*gd8$G5*gd8$G6)+(gd8$G1*gd8$G2*gd8$G3*gd8$G5*gd8$G7)+(gd8$G1*gd8$G2*gd8$G3*gd8$G5*gd8$G8)
				          	+(gd8$G1*gd8$G2*gd8$G3*gd8$G6*gd8$G7)+(gd8$G1*gd8$G2*gd8$G3*gd8$G6*gd8$G8)+(gd8$G1*gd8$G2*gd8$G3*gd8$G7*gd8$G8)
	        		  		+(gd8$G1*gd8$G2*gd8$G4*gd8$G5*gd8$G6)+(gd8$G1*gd8$G2*gd8$G4*gd8$G5*gd8$G7)+(gd8$G1*gd8$G2*gd8$G4*gd8$G5*gd8$G8)
	        	  			+(gd8$G1*gd8$G2*gd8$G4*gd8$G6*gd8$G7)+(gd8$G1*gd8$G2*gd8$G4*gd8$G6*gd8$G8)+(gd8$G1*gd8$G2*gd8$G4*gd8$G7*gd8$G8)
	          				+(gd8$G1*gd8$G2*gd8$G5*gd8$G6*gd8$G7)+(gd8$G1*gd8$G2*gd8$G5*gd8$G6*gd8$G8)+(gd8$G1*gd8$G2*gd8$G5*gd8$G7*gd8$G8)
		          			+(gd8$G1*gd8$G2*gd8$G6*gd8$G7*gd8$G8)
				          	+(gd8$G1*gd8$G3*gd8$G4*gd8$G5*gd8$G6)+(gd8$G1*gd8$G3*gd8$G4*gd8$G5*gd8$G7)+(gd8$G1*gd8$G3*gd8$G4*gd8$G5*gd8$G8)
				          	+(gd8$G1*gd8$G3*gd8$G4*gd8$G6*gd8$G7)+(gd8$G1*gd8$G3*gd8$G4*gd8$G6*gd8$G8)+(gd8$G1*gd8$G3*gd8$G4*gd8$G7*gd8$G8)
					          +(gd8$G1*gd8$G3*gd8$G5*gd8$G6*gd8$G7)+(gd8$G1*gd8$G3*gd8$G5*gd8$G6*gd8$G8)+(gd8$G1*gd8$G3*gd8$G5*gd8$G7*gd8$G8)
					          +(gd8$G1*gd8$G3*gd8$G6*gd8$G7*gd8$G8)
					          +(gd8$G1*gd8$G4*gd8$G5*gd8$G6*gd8$G7)+(gd8$G1*gd8$G4*gd8$G5*gd8$G6*gd8$G8)+(gd8$G1*gd8$G4*gd8$G5*gd8$G7*gd8$G8)
				          	+(gd8$G1*gd8$G4*gd8$G6*gd8$G7*gd8$G8)
				          	+(gd8$G1*gd8$G5*gd8$G6*gd8$G7*gd8$G8)
                    +(gd8$G2*gd8$G3*gd8$G4*gd8$G5*gd8$G6)+(gd8$G2*gd8$G3*gd8$G4*gd8$G5*gd8$G7)+(gd8$G2*gd8$G3*gd8$G4*gd8$G5*gd8$G8)
			          		+(gd8$G2*gd8$G3*gd8$G4*gd8$G6*gd8$G7)+(gd8$G2*gd8$G3*gd8$G4*gd8$G6*gd8$G8)+(gd8$G2*gd8$G3*gd8$G4*gd8$G7*gd8$G8)
			          		+(gd8$G2*gd8$G3*gd8$G5*gd8$G6*gd8$G7)+(gd8$G2*gd8$G3*gd8$G5*gd8$G6*gd8$G8)+(gd8$G2*gd8$G3*gd8$G5*gd8$G7*gd8$G8)
			          		+(gd8$G2*gd8$G3*gd8$G6*gd8$G7*gd8$G8)
			          		+(gd8$G2*gd8$G4*gd8$G5*gd8$G6*gd8$G7)+(gd8$G2*gd8$G4*gd8$G5*gd8$G6*gd8$G8)+(gd8$G2*gd8$G4*gd8$G5*gd8$G7*gd8$G8)
			          		+(gd8$G2*gd8$G4*gd8$G6*gd8$G7*gd8$G8)
			          		+(gd8$G2*gd8$G5*gd8$G6*gd8$G7*gd8$G8)
		          			+(gd8$G3*gd8$G4*gd8$G5*gd8$G6*gd8$G7)+(gd8$G3*gd8$G4*gd8$G5*gd8$G6*gd8$G8)+(gd8$G3*gd8$G4*gd8$G5*gd8$G7*gd8$G8)
		          			+(gd8$G3*gd8$G4*gd8$G6*gd8$G7*gd8$G8)
		          			+(gd8$G3*gd8$G5*gd8$G6*gd8$G7*gd8$G8)
                    +(gd8$G4*gd8$G5*gd8$G6*gd8$G7*gd8$G8) 
                  -(gd8$G1*gd8$G2*gd8$G3*gd8$G4*gd8$G5*gd8$G6)-(gd8$G1*gd8$G2*gd8$G3*gd8$G4*gd8$G5*gd8$G7)-(gd8$G1*gd8$G2*gd8$G3*gd8$G4*gd8$G5*gd8$G8)
		          			-(gd8$G1*gd8$G2*gd8$G3*gd8$G4*gd8$G6*gd8$G7)-(gd8$G1*gd8$G2*gd8$G3*gd8$G4*gd8$G6*gd8$G8)-(gd8$G1*gd8$G2*gd8$G3*gd8$G4*gd8$G7*gd8$G8)
		          			-(gd8$G1*gd8$G2*gd8$G3*gd8$G5*gd8$G6*gd8$G7)-(gd8$G1*gd8$G2*gd8$G3*gd8$G5*gd8$G6*gd8$G8)-(gd8$G1*gd8$G2*gd8$G3*gd8$G5*gd8$G7*gd8$G8)
		          			-(gd8$G1*gd8$G2*gd8$G3*gd8$G6*gd8$G7*gd8$G8)
		          			-(gd8$G1*gd8$G2*gd8$G4*gd8$G5*gd8$G6*gd8$G7)-(gd8$G1*gd8$G2*gd8$G4*gd8$G5*gd8$G6*gd8$G8)-(gd8$G1*gd8$G2*gd8$G4*gd8$G5*gd8$G7*gd8$G8)
		          			-(gd8$G1*gd8$G2*gd8$G4*gd8$G6*gd8$G7*gd8$G8)
		          			-(gd8$G1*gd8$G2*gd8$G5*gd8$G6*gd8$G7*gd8$G8)
		          			-(gd8$G1*gd8$G3*gd8$G4*gd8$G5*gd8$G6*gd8$G7)-(gd8$G1*gd8$G3*gd8$G4*gd8$G5*gd8$G6*gd8$G8)-(gd8$G1*gd8$G3*gd8$G4*gd8$G5*gd8$G7*gd8$G8)
		          			-(gd8$G1*gd8$G3*gd8$G4*gd8$G6*gd8$G7*gd8$G8)
		          			-(gd8$G1*gd8$G3*gd8$G5*gd8$G6*gd8$G7*gd8$G8)
		          			-(gd8$G1*gd8$G4*gd8$G5*gd8$G6*gd8$G7*gd8$G8) 
		          			-(gd8$G2*gd8$G3*gd8$G4*gd8$G5*gd8$G6*gd8$G7)-(gd8$G2*gd8$G3*gd8$G4*gd8$G5*gd8$G6*gd8$G8)-(gd8$G2*gd8$G3*gd8$G4*gd8$G5*gd8$G7*gd8$G8)
		          			-(gd8$G2*gd8$G3*gd8$G4*gd8$G6*gd8$G7*gd8$G8)
		          			-(gd8$G2*gd8$G3*gd8$G5*gd8$G6*gd8$G7*gd8$G8)
		          			-(gd8$G2*gd8$G4*gd8$G5*gd8$G6*gd8$G7*gd8$G8)
                    -(gd8$G3*gd8$G4*gd8$G5*gd8$G6*gd8$G7*gd8$G8)  
                  +(gd8$G1*gd8$G2*gd8$G3*gd8$G4*gd8$G5*gd8$G6*gd8$G7)+(gd8$G1*gd8$G2*gd8$G3*gd8$G4*gd8$G5*gd8$G6*gd8$G8)
                    +(gd8$G1*gd8$G2*gd8$G3*gd8$G4*gd8$G5*gd8$G7*gd8$G8)+(gd8$G1*gd8$G2*gd8$G3*gd8$G4*gd8$G6*gd8$G7*gd8$G8)
                    +(gd8$G1*gd8$G2*gd8$G3*gd8$G5*gd8$G6*gd8$G7*gd8$G8)+(gd8$G1*gd8$G2*gd8$G4*gd8$G5*gd8$G6*gd8$G7*gd8$G8)
                    +(gd8$G1*gd8$G3*gd8$G4*gd8$G5*gd8$G6*gd8$G7*gd8$G8)+(gd8$G2*gd8$G3*gd8$G4*gd8$G5*gd8$G6*gd8$G7*gd8$G8) #ok - all 8 7-prob combos done
                  -(gd8$G1*gd8$G2*gd8$G3*gd8$G4*gd8$G5*gd8$G6*gd8$G7*gd8$G8)
                    ))  
    gd8[,c("G1","G2","G3","G4","G5","G6","G7","G8") :=NULL]

expvals3<-rbind(gd1,gd4,gd8) # for low comp
setkey(expvals3,PLOT) # for low comp
kpc15_exp<-kpc15_exp[expvals3] # for low comp


gd1<-subset(kpc_indiv,GENDIV==1)
  gd1[,c("COMP","SEEDNUM","GENDIV","GTYPE","ABUND_13","ABUND_15","PC2USE","GTYPE_1","a15_no0_exp","pres15_exp"):=NULL]  #FOR HIGH COMP
gd4<-subset(kpc_indiv,GENDIV==4)  
  gd4[,c("COMP","SEEDNUM","GENDIV","GTYPE","ABUND_13","ABUND_15","PC2USE","GTYPE_1","a15_no0_exp","pres15_exp"):=NULL]   #FOR HIGH COMP
  gd4$G<-rep(c("G1","G2","G3","G4"),90)
  gd4<-reshape(gd4,idvar="PLOT", timevar="G",direction="wide")
  setnames(gd4,c("pres15_high_exp.G1","pres15_high_exp.G2","pres15_high_exp.G3","pres15_high_exp.G4"),c("G1","G2","G3","G4"))  #FOR HIGH COMP
  gd4$pres15_exp<-ifelse(is.na(gd4$G1),NA,
                  (gd4$G1+gd4$G2+gd4$G3+gd4$G4
                  -(gd4$G1*gd4$G2)-(gd4$G1*gd4$G3)-(gd4$G1*gd4$G4)-(gd4$G2*gd4$G3)-(gd4$G2*gd4$G4)-(gd4$G3*gd4$G4)
                  +(gd4$G1*gd4$G2*gd4$G3)+(gd4$G1*gd4$G2*gd4$G4)+(gd4$G1*gd4$G3*gd4$G4)+(gd4$G2*gd4$G3*gd4$G4)
                  -(gd4$G1*gd4$G2*gd4$G3*gd4$G4)))  #inclusion-exclusion works!!
  gd4[,c("G1","G2","G3","G4") :=NULL]
  setnames(gd4,"pres15_exp","pres15_high_exp")  #FOR HIGH COMP
gd8<-subset(kpc_indiv,GENDIV==8)  
  gd8[,c("COMP","SEEDNUM","GENDIV","GTYPE","ABUND_13","ABUND_15","PC2USE","GTYPE_1","a15_no0_exp","pres15_exp"):=NULL]   #FOR HIGH COMP
  gd8$G<-rep(c("G1","G2","G3","G4","G5","G6","G7","G8"),90)
  gd8<-reshape(gd8,idvar="PLOT", timevar="G",direction="wide")
  setnames(gd8,c("pres15_high_exp.G1","pres15_high_exp.G2","pres15_high_exp.G3","pres15_high_exp.G4","pres15_high_exp.G5","pres15_high_exp.G6","pres15_high_exp.G7","pres15_high_exp.G8"),c("G1","G2","G3","G4","G5","G6","G7","G8"))  #FOR HIGH COMP
  #next long batch of code aggregates persistence probabilities for 8gd treatments using inclusion-exclusion principle
  gd8$pres15_exp<-ifelse(is.na(gd8$G1),NA,
                  (gd8$G1+gd8$G2+gd8$G3+gd8$G4+gd8$G5+gd8$G6+gd8$G7+gd8$G8 
                  -(gd8$G1*gd8$G2)-(gd8$G1*gd8$G3)-(gd8$G1*gd8$G4)-(gd8$G1*gd8$G5)-(gd8$G1*gd8$G6)-(gd8$G1*gd8$G7)-(gd8$G1*gd8$G8)
                      -(gd8$G2*gd8$G3)-(gd8$G2*gd8$G4)-(gd8$G2*gd8$G5)-(gd8$G2*gd8$G6)-(gd8$G2*gd8$G7)-(gd8$G2*gd8$G8)
                      -(gd8$G3*gd8$G4)-(gd8$G3*gd8$G5)-(gd8$G3*gd8$G6)-(gd8$G3*gd8$G7)-(gd8$G3*gd8$G8)
                      -(gd8$G4*gd8$G5)-(gd8$G4*gd8$G6)-(gd8$G4*gd8$G7)-(gd8$G4*gd8$G8)
                      -(gd8$G5*gd8$G6)-(gd8$G5*gd8$G7)-(gd8$G5*gd8$G8)
                      -(gd8$G6*gd8$G7)-(gd8$G6*gd8$G8)
                      -(gd8$G7*gd8$G8)  
                  +(gd8$G1*gd8$G2*gd8$G3)+(gd8$G1*gd8$G2*gd8$G4)+(gd8$G1*gd8$G2*gd8$G5)+(gd8$G1*gd8$G2*gd8$G6)+(gd8$G1*gd8$G2*gd8$G7)+(gd8$G1*gd8$G2*gd8$G8)
                    +(gd8$G1*gd8$G3*gd8$G4)+(gd8$G1*gd8$G3*gd8$G5)+(gd8$G1*gd8$G3*gd8$G6)+(gd8$G1*gd8$G3*gd8$G7)+(gd8$G1*gd8$G3*gd8$G8)
                    +(gd8$G1*gd8$G4*gd8$G5)+(gd8$G1*gd8$G4*gd8$G6)+(gd8$G1*gd8$G4*gd8$G7)+(gd8$G1*gd8$G4*gd8$G8)
                    +(gd8$G1*gd8$G5*gd8$G6)+(gd8$G1*gd8$G5*gd8$G7)+(gd8$G1*gd8$G5*gd8$G8)
                    +(gd8$G1*gd8$G6*gd8$G7)+(gd8$G1*gd8$G6*gd8$G8)+(gd8$G1*gd8$G7*gd8$G8)
                    +(gd8$G2*gd8$G3*gd8$G4)+(gd8$G2*gd8$G3*gd8$G5)+(gd8$G2*gd8$G3*gd8$G6)+(gd8$G2*gd8$G3*gd8$G7)+(gd8$G2*gd8$G3*gd8$G8)
                    +(gd8$G2*gd8$G4*gd8$G5)+(gd8$G2*gd8$G4*gd8$G6)+(gd8$G2*gd8$G4*gd8$G7)+(gd8$G2*gd8$G4*gd8$G8)
                    +(gd8$G2*gd8$G5*gd8$G6)+(gd8$G2*gd8$G5*gd8$G7)+(gd8$G2*gd8$G5*gd8$G8)
                    +(gd8$G2*gd8$G6*gd8$G7)+(gd8$G2*gd8$G6*gd8$G8)+(gd8$G2*gd8$G7*gd8$G8)
                    +(gd8$G3*gd8$G4*gd8$G5)+(gd8$G3*gd8$G4*gd8$G6)+(gd8$G3*gd8$G4*gd8$G7)+(gd8$G3*gd8$G4*gd8$G8)
                    +(gd8$G3*gd8$G5*gd8$G6)+(gd8$G3*gd8$G5*gd8$G7)+(gd8$G3*gd8$G5*gd8$G8)
                    +(gd8$G3*gd8$G6*gd8$G7)+(gd8$G3*gd8$G6*gd8$G8)+(gd8$G3*gd8$G7*gd8$G8)
                    +(gd8$G4*gd8$G5*gd8$G6)+(gd8$G4*gd8$G5*gd8$G7)+(gd8$G4*gd8$G5*gd8$G8)
                    +(gd8$G4*gd8$G6*gd8$G7)+(gd8$G4*gd8$G6*gd8$G8)+(gd8$G4*gd8$G7*gd8$G8)
                    +(gd8$G5*gd8$G6*gd8$G7)+(gd8$G5*gd8$G6*gd8$G8)+(gd8$G5*gd8$G7*gd8$G8)
                    +(gd8$G6*gd8$G7*gd8$G8) 
                  -(gd8$G1*gd8$G2*gd8$G3*gd8$G4)-(gd8$G1*gd8$G2*gd8$G3*gd8$G5)-(gd8$G1*gd8$G2*gd8$G3*gd8$G6)-(gd8$G1*gd8$G2*gd8$G3*gd8$G7)-(gd8$G1*gd8$G2*gd8$G3*gd8$G8)
                    -(gd8$G1*gd8$G2*gd8$G4*gd8$G5)-(gd8$G1*gd8$G2*gd8$G4*gd8$G6)-(gd8$G1*gd8$G2*gd8$G4*gd8$G7)-(gd8$G1*gd8$G2*gd8$G4*gd8$G8)
                    -(gd8$G1*gd8$G2*gd8$G5*gd8$G6)-(gd8$G1*gd8$G2*gd8$G5*gd8$G7)-(gd8$G1*gd8$G2*gd8$G5*gd8$G8)
					          -(gd8$G1*gd8$G2*gd8$G6*gd8$G7)-(gd8$G1*gd8$G2*gd8$G6*gd8$G8)-(gd8$G1*gd8$G2*gd8$G7*gd8$G8)
                    -(gd8$G1*gd8$G3*gd8$G4*gd8$G5)-(gd8$G1*gd8$G3*gd8$G4*gd8$G6)-(gd8$G1*gd8$G3*gd8$G4*gd8$G7)-(gd8$G1*gd8$G3*gd8$G4*gd8$G8)
                    -(gd8$G1*gd8$G3*gd8$G5*gd8$G6)-(gd8$G1*gd8$G3*gd8$G5*gd8$G7)-(gd8$G1*gd8$G3*gd8$G5*gd8$G8)
				          	-(gd8$G1*gd8$G3*gd8$G6*gd8$G7)-(gd8$G1*gd8$G3*gd8$G6*gd8$G8)-(gd8$G1*gd8$G3*gd8$G7*gd8$G8)
                    -(gd8$G1*gd8$G4*gd8$G5*gd8$G6)-(gd8$G1*gd8$G4*gd8$G5*gd8$G7)-(gd8$G1*gd8$G4*gd8$G5*gd8$G8)
					          -(gd8$G1*gd8$G4*gd8$G6*gd8$G7)-(gd8$G1*gd8$G4*gd8$G6*gd8$G8)-(gd8$G1*gd8$G4*gd8$G7*gd8$G8)
					          -(gd8$G1*gd8$G5*gd8$G6*gd8$G7)-(gd8$G1*gd8$G5*gd8$G6*gd8$G8)-(gd8$G1*gd8$G5*gd8$G7*gd8$G8)
					          -(gd8$G1*gd8$G6*gd8$G7*gd8$G8)
                    -(gd8$G2*gd8$G3*gd8$G4*gd8$G5)-(gd8$G2*gd8$G3*gd8$G4*gd8$G6)-(gd8$G2*gd8$G3*gd8$G4*gd8$G7)-(gd8$G2*gd8$G3*gd8$G4*gd8$G8)
                    -(gd8$G2*gd8$G3*gd8$G5*gd8$G6)-(gd8$G2*gd8$G3*gd8$G5*gd8$G7)-(gd8$G2*gd8$G3*gd8$G5*gd8$G8)
					          -(gd8$G2*gd8$G3*gd8$G6*gd8$G7)-(gd8$G2*gd8$G3*gd8$G6*gd8$G8)-(gd8$G2*gd8$G3*gd8$G7*gd8$G8)
                    -(gd8$G2*gd8$G4*gd8$G5*gd8$G6)-(gd8$G2*gd8$G4*gd8$G5*gd8$G7)-(gd8$G2*gd8$G4*gd8$G5*gd8$G8)
					          -(gd8$G2*gd8$G4*gd8$G6*gd8$G7)-(gd8$G2*gd8$G4*gd8$G6*gd8$G8)-(gd8$G2*gd8$G4*gd8$G7*gd8$G8)
					          -(gd8$G2*gd8$G5*gd8$G6*gd8$G7)-(gd8$G2*gd8$G5*gd8$G6*gd8$G8)-(gd8$G2*gd8$G5*gd8$G7*gd8$G8)
					          -(gd8$G2*gd8$G6*gd8$G7*gd8$G8)
                    -(gd8$G3*gd8$G4*gd8$G5*gd8$G6)-(gd8$G3*gd8$G4*gd8$G5*gd8$G7)-(gd8$G3*gd8$G4*gd8$G5*gd8$G8)
					          -(gd8$G3*gd8$G4*gd8$G6*gd8$G7)-(gd8$G3*gd8$G4*gd8$G6*gd8$G8)-(gd8$G3*gd8$G4*gd8$G7*gd8$G8)
					          -(gd8$G3*gd8$G5*gd8$G6*gd8$G7)-(gd8$G3*gd8$G5*gd8$G6*gd8$G8)-(gd8$G3*gd8$G5*gd8$G7*gd8$G8)
					          -(gd8$G3*gd8$G6*gd8$G7*gd8$G8)					
                    -(gd8$G4*gd8$G5*gd8$G6*gd8$G7)-(gd8$G4*gd8$G5*gd8$G6*gd8$G8)-(gd8$G4*gd8$G5*gd8$G7*gd8$G8)-(gd8$G4*gd8$G6*gd8$G7*gd8$G8)
                    -(gd8$G5*gd8$G6*gd8$G7*gd8$G8) 
                  +(gd8$G1*gd8$G2*gd8$G3*gd8$G4*gd8$G5)+(gd8$G1*gd8$G2*gd8$G3*gd8$G4*gd8$G6)+(gd8$G1*gd8$G2*gd8$G3*gd8$G4*gd8$G7)+(gd8$G1*gd8$G2*gd8$G3*gd8$G4*gd8$G8)
					          +(gd8$G1*gd8$G2*gd8$G3*gd8$G5*gd8$G6)+(gd8$G1*gd8$G2*gd8$G3*gd8$G5*gd8$G7)+(gd8$G1*gd8$G2*gd8$G3*gd8$G5*gd8$G8)
				          	+(gd8$G1*gd8$G2*gd8$G3*gd8$G6*gd8$G7)+(gd8$G1*gd8$G2*gd8$G3*gd8$G6*gd8$G8)+(gd8$G1*gd8$G2*gd8$G3*gd8$G7*gd8$G8)
	        		  		+(gd8$G1*gd8$G2*gd8$G4*gd8$G5*gd8$G6)+(gd8$G1*gd8$G2*gd8$G4*gd8$G5*gd8$G7)+(gd8$G1*gd8$G2*gd8$G4*gd8$G5*gd8$G8)
	        	  			+(gd8$G1*gd8$G2*gd8$G4*gd8$G6*gd8$G7)+(gd8$G1*gd8$G2*gd8$G4*gd8$G6*gd8$G8)+(gd8$G1*gd8$G2*gd8$G4*gd8$G7*gd8$G8)
	          				+(gd8$G1*gd8$G2*gd8$G5*gd8$G6*gd8$G7)+(gd8$G1*gd8$G2*gd8$G5*gd8$G6*gd8$G8)+(gd8$G1*gd8$G2*gd8$G5*gd8$G7*gd8$G8)
		          			+(gd8$G1*gd8$G2*gd8$G6*gd8$G7*gd8$G8)
				          	+(gd8$G1*gd8$G3*gd8$G4*gd8$G5*gd8$G6)+(gd8$G1*gd8$G3*gd8$G4*gd8$G5*gd8$G7)+(gd8$G1*gd8$G3*gd8$G4*gd8$G5*gd8$G8)
				          	+(gd8$G1*gd8$G3*gd8$G4*gd8$G6*gd8$G7)+(gd8$G1*gd8$G3*gd8$G4*gd8$G6*gd8$G8)+(gd8$G1*gd8$G3*gd8$G4*gd8$G7*gd8$G8)
					          +(gd8$G1*gd8$G3*gd8$G5*gd8$G6*gd8$G7)+(gd8$G1*gd8$G3*gd8$G5*gd8$G6*gd8$G8)+(gd8$G1*gd8$G3*gd8$G5*gd8$G7*gd8$G8)
					          +(gd8$G1*gd8$G3*gd8$G6*gd8$G7*gd8$G8)
					          +(gd8$G1*gd8$G4*gd8$G5*gd8$G6*gd8$G7)+(gd8$G1*gd8$G4*gd8$G5*gd8$G6*gd8$G8)+(gd8$G1*gd8$G4*gd8$G5*gd8$G7*gd8$G8)
				          	+(gd8$G1*gd8$G4*gd8$G6*gd8$G7*gd8$G8)
				          	+(gd8$G1*gd8$G5*gd8$G6*gd8$G7*gd8$G8)
                    +(gd8$G2*gd8$G3*gd8$G4*gd8$G5*gd8$G6)+(gd8$G2*gd8$G3*gd8$G4*gd8$G5*gd8$G7)+(gd8$G2*gd8$G3*gd8$G4*gd8$G5*gd8$G8)
			          		+(gd8$G2*gd8$G3*gd8$G4*gd8$G6*gd8$G7)+(gd8$G2*gd8$G3*gd8$G4*gd8$G6*gd8$G8)+(gd8$G2*gd8$G3*gd8$G4*gd8$G7*gd8$G8)
			          		+(gd8$G2*gd8$G3*gd8$G5*gd8$G6*gd8$G7)+(gd8$G2*gd8$G3*gd8$G5*gd8$G6*gd8$G8)+(gd8$G2*gd8$G3*gd8$G5*gd8$G7*gd8$G8)
			          		+(gd8$G2*gd8$G3*gd8$G6*gd8$G7*gd8$G8)
			          		+(gd8$G2*gd8$G4*gd8$G5*gd8$G6*gd8$G7)+(gd8$G2*gd8$G4*gd8$G5*gd8$G6*gd8$G8)+(gd8$G2*gd8$G4*gd8$G5*gd8$G7*gd8$G8)
			          		+(gd8$G2*gd8$G4*gd8$G6*gd8$G7*gd8$G8)
			          		+(gd8$G2*gd8$G5*gd8$G6*gd8$G7*gd8$G8)
		          			+(gd8$G3*gd8$G4*gd8$G5*gd8$G6*gd8$G7)+(gd8$G3*gd8$G4*gd8$G5*gd8$G6*gd8$G8)+(gd8$G3*gd8$G4*gd8$G5*gd8$G7*gd8$G8)
		          			+(gd8$G3*gd8$G4*gd8$G6*gd8$G7*gd8$G8)
		          			+(gd8$G3*gd8$G5*gd8$G6*gd8$G7*gd8$G8)
                    +(gd8$G4*gd8$G5*gd8$G6*gd8$G7*gd8$G8) 
                  -(gd8$G1*gd8$G2*gd8$G3*gd8$G4*gd8$G5*gd8$G6)-(gd8$G1*gd8$G2*gd8$G3*gd8$G4*gd8$G5*gd8$G7)-(gd8$G1*gd8$G2*gd8$G3*gd8$G4*gd8$G5*gd8$G8)
		          			-(gd8$G1*gd8$G2*gd8$G3*gd8$G4*gd8$G6*gd8$G7)-(gd8$G1*gd8$G2*gd8$G3*gd8$G4*gd8$G6*gd8$G8)-(gd8$G1*gd8$G2*gd8$G3*gd8$G4*gd8$G7*gd8$G8)
		          			-(gd8$G1*gd8$G2*gd8$G3*gd8$G5*gd8$G6*gd8$G7)-(gd8$G1*gd8$G2*gd8$G3*gd8$G5*gd8$G6*gd8$G8)-(gd8$G1*gd8$G2*gd8$G3*gd8$G5*gd8$G7*gd8$G8)
		          			-(gd8$G1*gd8$G2*gd8$G3*gd8$G6*gd8$G7*gd8$G8)
		          			-(gd8$G1*gd8$G2*gd8$G4*gd8$G5*gd8$G6*gd8$G7)-(gd8$G1*gd8$G2*gd8$G4*gd8$G5*gd8$G6*gd8$G8)-(gd8$G1*gd8$G2*gd8$G4*gd8$G5*gd8$G7*gd8$G8)
		          			-(gd8$G1*gd8$G2*gd8$G4*gd8$G6*gd8$G7*gd8$G8)
		          			-(gd8$G1*gd8$G2*gd8$G5*gd8$G6*gd8$G7*gd8$G8)
		          			-(gd8$G1*gd8$G3*gd8$G4*gd8$G5*gd8$G6*gd8$G7)-(gd8$G1*gd8$G3*gd8$G4*gd8$G5*gd8$G6*gd8$G8)-(gd8$G1*gd8$G3*gd8$G4*gd8$G5*gd8$G7*gd8$G8)
		          			-(gd8$G1*gd8$G3*gd8$G4*gd8$G6*gd8$G7*gd8$G8)
		          			-(gd8$G1*gd8$G3*gd8$G5*gd8$G6*gd8$G7*gd8$G8)
		          			-(gd8$G1*gd8$G4*gd8$G5*gd8$G6*gd8$G7*gd8$G8) 
		          			-(gd8$G2*gd8$G3*gd8$G4*gd8$G5*gd8$G6*gd8$G7)-(gd8$G2*gd8$G3*gd8$G4*gd8$G5*gd8$G6*gd8$G8)-(gd8$G2*gd8$G3*gd8$G4*gd8$G5*gd8$G7*gd8$G8)
		          			-(gd8$G2*gd8$G3*gd8$G4*gd8$G6*gd8$G7*gd8$G8)
		          			-(gd8$G2*gd8$G3*gd8$G5*gd8$G6*gd8$G7*gd8$G8)
		          			-(gd8$G2*gd8$G4*gd8$G5*gd8$G6*gd8$G7*gd8$G8)
                    -(gd8$G3*gd8$G4*gd8$G5*gd8$G6*gd8$G7*gd8$G8)  
                  +(gd8$G1*gd8$G2*gd8$G3*gd8$G4*gd8$G5*gd8$G6*gd8$G7)+(gd8$G1*gd8$G2*gd8$G3*gd8$G4*gd8$G5*gd8$G6*gd8$G8)
                    +(gd8$G1*gd8$G2*gd8$G3*gd8$G4*gd8$G5*gd8$G7*gd8$G8)+(gd8$G1*gd8$G2*gd8$G3*gd8$G4*gd8$G6*gd8$G7*gd8$G8)
                    +(gd8$G1*gd8$G2*gd8$G3*gd8$G5*gd8$G6*gd8$G7*gd8$G8)+(gd8$G1*gd8$G2*gd8$G4*gd8$G5*gd8$G6*gd8$G7*gd8$G8)
                    +(gd8$G1*gd8$G3*gd8$G4*gd8$G5*gd8$G6*gd8$G7*gd8$G8)+(gd8$G2*gd8$G3*gd8$G4*gd8$G5*gd8$G6*gd8$G7*gd8$G8) #ok - all 8 7-prob combos done
                  -(gd8$G1*gd8$G2*gd8$G3*gd8$G4*gd8$G5*gd8$G6*gd8$G7*gd8$G8)
                    ))  
    gd8[,c("G1","G2","G3","G4","G5","G6","G7","G8") :=NULL]
    setnames(gd8,"pres15_exp","pres15_high_exp")   #FOR HIGH COMP

expvals4<-rbind(gd1,gd4,gd8)   #FOR HIGH COMP
  setkey(expvals4,PLOT)  #FOR HIGH COMP
kpc15_exp<-kpc15_exp[expvals4]  #FOR HIGH COMP

  #4. calculate the difference between observed response versus expected response
kpc15_exp$a15_no0_mean_diff<-kpc15_exp$FL_SUM15-kpc15_exp$a15_no0_exp.mean
kpc15_exp$pres15_diff<-kpc15_exp$PRES15-kpc15_exp$pres15_exp
  kpc15_exp$pres15_high_diff<-kpc15_exp$PRES15-kpc15_exp$pres15_high_exp
kpc15_exp$pres15_max_diff<-kpc15_exp$PRES15-kpc15_exp$pres15_exp.max
  kpc15_exp$pres15_high_max_diff<-kpc15_exp$PRES15-kpc15_exp$pres15_high_exp.max

  #5. run analysis on the difference (t-tests for obs vs exp).  Null hyp is DIFF=0 which would suggest sampling effects only are at play.
    #all diff calculated as obs-exp, with exp as the summed exp for each gtype in the plot, given the per-gtype seed number added 
    #5b. abundance diff, excluding zeroes (a15_no0_sum_diff); note that the distribution is really non-normal, so the k-w tests at the end are probably the best (and give the same results as anova)    
        #note that for abundance we're looking at data aggregated as means across gtypes, not sums
t.test(kpc15_exp[GENDIV==8]$FL_SUM15,kpc15_exp[GENDIV==8]$a15_no0_exp.mean,paired=T) #
  t.test(kpc15_exp[GENDIV==8 & SEEDNUM==960]$FL_SUM15,kpc15_exp[GENDIV==8 & SEEDNUM==960]$a15_no0_exp.mean,paired=T) #
  t.test(kpc15_exp[GENDIV==8 & SEEDNUM==320]$FL_SUM15,kpc15_exp[GENDIV==8 & SEEDNUM==320]$a15_no0_exp.mean,paired=T) #
  t.test(kpc15_exp[GENDIV==8 & SEEDNUM==160]$FL_SUM15,kpc15_exp[GENDIV==8 & SEEDNUM==160]$a15_no0_exp.mean,paired=T) #
  t.test(kpc15_exp[GENDIV==8 & SEEDNUM==96]$FL_SUM15,kpc15_exp[GENDIV==8 & SEEDNUM==96]$a15_no0_exp.mean,paired=T) #
  t.test(kpc15_exp[GENDIV==8 & SEEDNUM==32]$FL_SUM15,kpc15_exp[GENDIV==8 & SEEDNUM==32]$a15_no0_exp.mean,paired=T) #
t.test(kpc15_exp[GENDIV==4]$FL_SUM15,kpc15_exp[GENDIV==4]$a15_no0_exp.mean,paired=T) #
  t.test(kpc15_exp[GENDIV==4 & SEEDNUM==960]$FL_SUM15,kpc15_exp[GENDIV==4 & SEEDNUM==960]$a15_no0_exp.mean,paired=T) #
  t.test(kpc15_exp[GENDIV==4 & SEEDNUM==320]$FL_SUM15,kpc15_exp[GENDIV==4 & SEEDNUM==320]$a15_no0_exp.mean,paired=T) #
  t.test(kpc15_exp[GENDIV==4 & SEEDNUM==160]$FL_SUM15,kpc15_exp[GENDIV==4 & SEEDNUM==160]$a15_no0_exp.mean,paired=T) #
  t.test(kpc15_exp[GENDIV==4 & SEEDNUM==96]$FL_SUM15,kpc15_exp[GENDIV==4 & SEEDNUM==96]$a15_no0_exp.mean,paired=T) #
  t.test(kpc15_exp[GENDIV==4 & SEEDNUM==32]$FL_SUM15,kpc15_exp[GENDIV==4 & SEEDNUM==32]$a15_no0_exp.mean,paired=T) #
t.test(kpc15_exp[GENDIV==1]$FL_SUM15,kpc15_exp[GENDIV==1]$a15_no0_exp.mean,paired=T) # 
  t.test(kpc15_exp[GENDIV==1 & SEEDNUM==960]$FL_SUM15,kpc15_exp[GENDIV==1 & SEEDNUM==960]$a15_no0_exp.mean,paired=T) #
  t.test(kpc15_exp[GENDIV==1 & SEEDNUM==320]$FL_SUM15,kpc15_exp[GENDIV==1 & SEEDNUM==320]$a15_no0_exp.mean,paired=T) #
  t.test(kpc15_exp[GENDIV==1 & SEEDNUM==160]$FL_SUM15,kpc15_exp[GENDIV==1 & SEEDNUM==160]$a15_no0_exp.mean,paired=T) #
  t.test(kpc15_exp[GENDIV==1 & SEEDNUM==96]$FL_SUM15,kpc15_exp[GENDIV==1 & SEEDNUM==96]$a15_no0_exp.mean,paired=T) #
  t.test(kpc15_exp[GENDIV==1 & SEEDNUM==32]$FL_SUM15,kpc15_exp[GENDIV==1 & SEEDNUM==32]$a15_no0_exp.mean,paired=T) #


    #5c. persistence diff (pres15_diff)    
t.test(kpc15_exp[GENDIV==8]$PRES15,kpc15_exp[GENDIV==8]$pres15_exp,paired=T) #
  t.test(kpc15_exp[GENDIV==8 & SEEDNUM==960]$PRES15,kpc15_exp[GENDIV==8 & SEEDNUM==960]$pres15_exp,paired=T, alternative = "two") #
  t.test(kpc15_exp[GENDIV==8 & SEEDNUM==320]$PRES15,kpc15_exp[GENDIV==8 & SEEDNUM==320]$pres15_exp,paired=T,alternative = "two") #
  t.test(kpc15_exp[GENDIV==8 & SEEDNUM==160]$PRES15,kpc15_exp[GENDIV==8 & SEEDNUM==160]$pres15_exp,paired=T,alternative = "two") #
  t.test(kpc15_exp[GENDIV==8 & SEEDNUM==96]$PRES15,kpc15_exp[GENDIV==8 & SEEDNUM==96]$pres15_exp,paired=T,alternative = "two") #
  t.test(kpc15_exp[GENDIV==8 & SEEDNUM==32]$PRES15,kpc15_exp[GENDIV==8 & SEEDNUM==32]$pres15_exp,paired=T,alternative = "two") #
t.test(kpc15_exp[GENDIV==4]$PRES15,kpc15_exp[GENDIV==4]$pres15_exp,paired=T) #
  t.test(kpc15_exp[GENDIV==4 & SEEDNUM==960]$PRES15,kpc15_exp[GENDIV==4 & SEEDNUM==960]$pres15_exp,paired=T,alternative = "two") #
  t.test(kpc15_exp[GENDIV==4 & SEEDNUM==320]$PRES15,kpc15_exp[GENDIV==4 & SEEDNUM==320]$pres15_exp,paired=T,alternative = "two") #
  t.test(kpc15_exp[GENDIV==4 & SEEDNUM==160]$PRES15,kpc15_exp[GENDIV==4 & SEEDNUM==160]$pres15_exp,paired=T,alternative = "two") #
  t.test(kpc15_exp[GENDIV==4 & SEEDNUM==96]$PRES15,kpc15_exp[GENDIV==4 & SEEDNUM==96]$pres15_exp,paired=T,alternative = "two") #
  t.test(kpc15_exp[GENDIV==4 & SEEDNUM==32]$PRES15,kpc15_exp[GENDIV==4 & SEEDNUM==32]$pres15_exp,paired=T,alternative = "two") #
t.test(kpc15_exp[GENDIV==1]$PRES15,kpc15_exp[GENDIV==1]$pres15_exp,paired=T ,alternative = "less") 
  t.test(kpc15_exp[GENDIV==1 & SEEDNUM==960]$PRES15,kpc15_exp[GENDIV==1 & SEEDNUM==960]$pres15_exp,paired=T,alternative = "two") #
  t.test(kpc15_exp[GENDIV==1 & SEEDNUM==320]$PRES15,kpc15_exp[GENDIV==1 & SEEDNUM==320]$pres15_exp,paired=T,alternative = "two") #
  t.test(kpc15_exp[GENDIV==1 & SEEDNUM==160]$PRES15,kpc15_exp[GENDIV==1 & SEEDNUM==160]$pres15_exp,paired=T,alternative = "two") #
  t.test(kpc15_exp[GENDIV==1 & SEEDNUM==96]$PRES15,kpc15_exp[GENDIV==1 & SEEDNUM==96]$pres15_exp,paired=T,alternative = "two") #
  t.test(kpc15_exp[GENDIV==1 & SEEDNUM==32]$PRES15,kpc15_exp[GENDIV==1 & SEEDNUM==32]$pres15_exp,paired=T,alternative = "two") #

    #5d. persistence diff in high comp (pres15_diff)    
t.test(kpc15_exp[GENDIV==8]$PRES15,kpc15_exp[GENDIV==8]$pres15_high_exp,paired=T) # 
  t.test(kpc15_exp[GENDIV==8 & SEEDNUM==960]$PRES15,kpc15_exp[GENDIV==8 & SEEDNUM==960]$pres15_high_exp,paired=T) 
  t.test(kpc15_exp[GENDIV==8 & SEEDNUM==320]$PRES15,kpc15_exp[GENDIV==8 & SEEDNUM==320]$pres15_high_exp,paired=T) 
  t.test(kpc15_exp[GENDIV==8 & SEEDNUM==160]$PRES15,kpc15_exp[GENDIV==8 & SEEDNUM==160]$pres15_high_exp,paired=T) #
  t.test(kpc15_exp[GENDIV==8 & SEEDNUM==96]$PRES15,kpc15_exp[GENDIV==8 & SEEDNUM==96]$pres15_high_exp,paired=T) #
  t.test(kpc15_exp[GENDIV==8 & SEEDNUM==32]$PRES15,kpc15_exp[GENDIV==8 & SEEDNUM==32]$pres15_high_exp,paired=T) #
t.test(kpc15_exp[GENDIV==4]$PRES15,kpc15_exp[GENDIV==4]$pres15_high_exp,paired=T) #
  t.test(kpc15_exp[GENDIV==4 & SEEDNUM==960]$PRES15,kpc15_exp[GENDIV==4 & SEEDNUM==960]$pres15_high_exp,paired=T) #
  t.test(kpc15_exp[GENDIV==4 & SEEDNUM==320]$PRES15,kpc15_exp[GENDIV==4 & SEEDNUM==320]$pres15_high_exp,paired=T) #
  t.test(kpc15_exp[GENDIV==4 & SEEDNUM==160]$PRES15,kpc15_exp[GENDIV==4 & SEEDNUM==160]$pres15_high_exp,paired=T) #
  t.test(kpc15_exp[GENDIV==4 & SEEDNUM==96]$PRES15,kpc15_exp[GENDIV==4 & SEEDNUM==96]$pres15_high_exp,paired=T) #
  t.test(kpc15_exp[GENDIV==4 & SEEDNUM==32]$PRES15,kpc15_exp[GENDIV==4 & SEEDNUM==32]$pres15_high_exp,paired=T) #
t.test(kpc15_exp[GENDIV==1]$PRES15,kpc15_exp[GENDIV==1]$pres15_high_exp,paired=T ) #
  t.test(kpc15_exp[GENDIV==1 & SEEDNUM==960]$PRES15,kpc15_exp[GENDIV==1 & SEEDNUM==960]$pres15_high_exp,paired=T) #
  t.test(kpc15_exp[GENDIV==1 & SEEDNUM==320]$PRES15,kpc15_exp[GENDIV==1 & SEEDNUM==320]$pres15_high_exp,paired=T) #
  t.test(kpc15_exp[GENDIV==1 & SEEDNUM==160]$PRES15,kpc15_exp[GENDIV==1 & SEEDNUM==160]$pres15_high_exp,paired=T) #
  t.test(kpc15_exp[GENDIV==1 & SEEDNUM==96]$PRES15,kpc15_exp[GENDIV==1 & SEEDNUM==96]$pres15_high_exp,paired=T) #
  t.test(kpc15_exp[GENDIV==1 & SEEDNUM==32]$PRES15,kpc15_exp[GENDIV==1 & SEEDNUM==32]$pres15_high_exp,paired=T) #
  

###############################
###############################
###############################
