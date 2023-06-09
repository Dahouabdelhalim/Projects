### R-Script 'Cross -activity of honeybee queen mandibular pheromone in bumblebees provides evidence for sensory exploitation' 

# load required packages
library(afex)
library(gplots)
library(lme4)
library(lmerTest)
library(ordinal)
library(emmeans)
library(effects)
library(multcomp)
library(ggplot2)
library(ggthemes)
library(export) 
library(broom)
library(car)

###### A. ANALYSIS OF INDIVIDUAL WORKER DATA ####

workerdata=read.table("WorkerDataOvaryDevelopment.csv", header=TRUE, sep=",", dec=".")
workerdata$Treatment=factor(workerdata$Treatment, levels=c("Control", "n-C25", "QMP"))
workerdata$Body_size=-prcomp(cbind(workerdata$Head.width,workerdata$Mean.marginal.cell))$x[,1] # we summarize head width & mean marginal cell length into one summary measure body size based on a PCA
rot=-prcomp(cbind(workerdata$Head.width,workerdata$Mean.marginal.cell))$rotation[,1] # factor loadings
#as.vector(rot[1]*scale(workerdata$Head.width,center=T,scale=F)+rot[2]*scale(workerdata$Mean.marginal.cell,center=T,scale=F))
# PCA projection based measure of body size was given by
# Body_size = 0.7483498*(head_width_in_mm-mean_head_width)+0.6633042*(marginal_cell_length_in_mm-mean_marginal_cell_length)
#           = 0.7483498*(head_width_in_mm-4.159583)+0.6633042*(Mean_marginal_cell_length_in_mm-2.955417)

# A1. ANALYSIS OF PROPORTION OF WORKERS WITH ACTIVATED OVARIES (FIG. 1A) ####

set_treatment_contrasts() 

# here random intercept model with colony & nest-level overdispersion gave best AIC
# different models that were fitted:
fit1A=glmer(Developed_ovaries~Treatment+scale(Body_size)+(1|Colony/Nest)+scale(Total.workers.7days)+scale(Mortality), data=workerdata, family = binomial, control = glmerControl(check.conv.singular = .makeCC(action="warning", tol = 1e-10)))
fit1B=glmer(Developed_ovaries~Treatment+scale(Body_size)+(1+Treatment|Colony)+(1|Nest)+scale(Total.workers.7days)+scale(Mortality), data=workerdata, family = binomial, control = glmerControl(check.conv.grad = .makeCC("warning", tol = 1e-2)))
AIC(fit1A,fit1B) # best model (lowest AIC) is fit1A, used for FIG. 1B
bestfit1=fit1A
sum1=tidy(bestfit1)
sum1$p.value=sum1$p.value/2 # we use 1-sided p values
as.data.frame(sum1)

#                         term      estimate  std.error statistic      p.value  group
# 1                (Intercept)  4.954123e-01 0.20917800  2.368376 8.933173e-03  fixed
# 2             Treatmentn-C25 -3.487520e-01 0.22522801 -1.548440 6.075826e-02  fixed
# 3               TreatmentQMP -4.601221e-01 0.22717123 -2.025442 2.141102e-02  fixed
# 4           scale(Body.size)  2.845678e-01 0.09676099  2.940936 1.636112e-03  fixed
# 5 scale(Total.workers.7days) -6.093273e-01 0.11508123 -5.294758 5.958695e-08  fixed
# 6           scale(Mortality) -2.849052e-01 0.10913748 -2.610517 4.520277e-03  fixed
# 7        sd_(Intercept).Nest  4.289318e-05         NA        NA           NA   Nest
# 8      sd_(Intercept).Colony  3.805035e-01         NA        NA           NA Colony

# Effect plots of the significant additional covariates (see supplementary material)
# using effects package

plot(effect('scale(Body_size)', mod=bestfit1, confint=list(level=0.9)), ylim=c(0,1),
     type="response",ylab="Workers with activated ovaries (%)", main = "") 
graph2ppt(file="FigS1a.pptx",width=5,aspectr=sqrt(2))

plot(effect('Total.workers.7days', mod=bestfit1, confint=list(level=0.9)), ylim=c(0,1), 
     type="response",xlab="Nr of workers in the nest", ylab="Workers with activated ovaries (%)", main = "") 
graph2ppt(file="FigS1b.pptx",width=5,aspectr=sqrt(2))

plot(effect('Mortality', mod=bestfit1, confint=list(level=0.9)), ylim=c(0,1), 
     type="response",ylab="Workers with activated ovaries (%)", main = "") 
graph2ppt(file="FigS1c.pptx",width=5,aspectr=sqrt(2))

# Control-treatment posthoc comparisons with FDR corrected p values
emmeans1=summary(emmeans(bestfit1,"Treatment", contr = "trt.vs.ctrl", adjust="fdr"),type="response",level=0.90) # FDR posthoc comparisons, using 95% confidence bounds (=1-sided confidence intervals)
emmeans1$contrasts$p.value=emmeans1$contrasts$p.value/2 # we use 1 sided p values
emmeans1

# $emmeans
# Treatment  prob     SE  df asymp.LCL asymp.UCL
# Control   0.621 0.0492 Inf     0.538     0.698
# n-C25     0.537 0.0471 Inf     0.459     0.613
# QMP       0.509 0.0485 Inf     0.429     0.588
# 
# Confidence level used: 0.9 
# Intervals are back-transformed from the logit scale 
# 
# $contrasts
# contrast        odds.ratio    SE  df z.ratio p.value
# n-C25 / Control      0.706 0.159 Inf -1.548  0.0608 
# QMP / Control        0.631 0.143 Inf -2.025  0.0428 
# 
# P value adjustment: fdr method for 2 tests 
# Tests are performed on the log odds ratio scale

# Effect sizes in terms of odds ratios of workers not having activated ovaries compared to control with 95% confidence bounds
df_posthoc1=data.frame(summary(contrast(emmeans(bestfit1,"Treatment", type="response"),
                                        method="trt.vs.ctrl", group="Control", adjust="FDR", level=0.90)))
df_posthoc1$p.value=df_posthoc1$p.value/2 # we use 1-sided FDR corrected p values
df_posthoc1$odds.ratio.lower95CB=df_posthoc1$odds.ratio+1.96*df_posthoc1$SE
df_posthoc1$odds.ratio.upper95CB=df_posthoc1$odds.ratio-1.96*df_posthoc1$SE
df_posthoc1=df_posthoc1[,-c(3:6)]
df_posthoc1[,c(2:4)]=1/df_posthoc1[,c(2:4)] # to express odds ratio for workers to have developed ovaries in control relative to in each treatment, ie to reverse contrast compared to default
df_posthoc1$contrast=gsub(" / Control", "", df_posthoc1$contrast)
df_posthoc1

#   contrast odds.ratio    odds.ratio.lower95CB odds.ratio.upper95CB
# 1    n-C25   1.417298    0.9832465             2.537445
# 2      QMP   1.584267    1.0961849             2.855851 


# Post hoc effects plotted using ggplot2 (FIG. 1B)
df1=data.frame(emmeans1$emmeans)
df1$cols=rep(rgb(0,0,0),length(levels(workerdata$Treatment)))
df1$fills=c("grey75","steelblue","firebrick3")
ggplot(df1, aes(x=Treatment, y=prob*100, ymin=asymp.LCL*100, ymax=asymp.UCL*100, colour=I(cols), fill=I(fills))) + 
  geom_col(width=0.5, colour=NA) +
  geom_errorbar(aes(ymin=asymp.LCL*100, ymax=asymp.UCL*100), width=0, size=0.5) +
  theme_few(base_size=14) +
  labs(y = "Workers with activated ovaries (%)", x="Treatment") +
  theme(legend.position="none", axis.line=element_line(colour="black"),
        panel.background=element_rect(fill="white"), axis.text=element_text(colour="black")) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  coord_cartesian(ylim=c(0,100))
graph2ppt(file="Fig1B.pptx",width=5,aspectr=sqrt(2))

# We should note that the random slope model fit1B, where the treatment effects were also
# allowed to randomly vary across colonies or nests (subcolonies) produced qualitatively identical results
# (we did not present this in the main text though since this model had a worse AIC value than
# the random intercept model fit1A) :
#       df      AIC
# fit1A  8 775.1061
# fit1B 13 784.7355



###### B. ANALYSIS OF NEST-LEVEL WORKER DATA (EGG-LAYING & AGGRESSION) ####

nestdata=read.csv("WorkerDataEgglayingAggression.csv")
nestdata$Treatment=factor(nestdata$Treatment, levels=c("Control", "n-C25", "QMP"))
nestdata$Avg.body.size=as.vector(rot[1]*scale(nestdata$Nest.head.width,center=T,scale=F)+rot[2]*scale(nestdata$Nest.marginal.cel,center=T,scale=F)) # average body size measure based on first PCA component, using projection as above
nestdata$NDeadW=as.vector(round(nestdata$Mortality*7))
nestdata$NTotalW=as.vector(nestdata$NDeadW+nestdata$Total.workers.7days)

### B0. ANALYSIS OF QMP TOXICITY ####

fit0=glmer(NDeadW~Treatment+scale(NTotalW)+(1|Colony/Nest), family = poisson (link=log), data = nestdata, 
           control = glmerControl(check.conv.singular = .makeCC(action="warning", tol = 1e-10)))

sum0=tidy(fit0)
sum0$p.value=sum0$p.value/2
as.data.frame(sum0)

#                         term    estimate std.error  statistic      p.value       group
# 1                (Intercept)  0.82265289 0.2485068  3.3103835 0.0004658411       fixed
# 2             Treatmentn-C25 -0.08364980 0.2920648 -0.2864083 0.3872826940       fixed
# 3               TreatmentQMP -0.01774646 0.2918845 -0.0607996 0.4757594028       fixed
# 4             scale(NTotalW)  0.25434882 0.1449021  1.7553145 0.0396027615       fixed
# 5 sd_(Intercept).Nest:Colony  0.05480573        NA         NA           NA Nest:Colony
# 6      sd_(Intercept).Colony  0.37841614        NA         NA           NA      Colony

# Overall mortality was low, as there were just 2.2 dead workers on average
# per nest of on average 30 workers at the end of the 7 day experiment
summary(emmeans(fit0,~1), type="response", level=0.95)
# 1       rate    SE  df asymp.LCL asymp.UCL
# overall  2.2 0.405 Inf      1.54      3.16
mean(nestdata$NTotalW) # mean nr of initial workers + eclosed workers + dead workers = 29.97

emmeans0=summary(emmeans(fit0,"Treatment", contr = "trt.vs.ctrl", adjust="fdr"),
                 type="response", level=0.9) # FDR post hoc comparisons with 95% confidence bounds
emmeans0$contrasts$p.value=emmeans0$contrasts$p.value/2 # we used 1 sided p values
emmeans0

# $emmeans
# Treatment rate    SE  df asymp.LCL asymp.UCL
# Control   2.28 0.566 Inf      1.51      3.43
# n-C25     2.09 0.530 Inf      1.38      3.17
# QMP       2.24 0.551 Inf      1.49      3.35
# 
# Confidence level used: 0.9 
# Intervals are back-transformed from the log scale 
# 
# $contrasts
# contrast        ratio    SE  df z.ratio p.value
# n-C25 / Control 0.920 0.269 Inf -0.286  0.4758 
# QMP / Control   0.982 0.287 Inf -0.061  0.4758 
# 
# P value adjustment: fdr method for 2 tests 
# Tests are performed on the log scale 

# Effect sizes in terms of fold reduction in nr dead workers in treatment vs control groups with 95% confidence bounds
df_posthoc0=data.frame(summary(contrast(emmeans(fit0,"Treatment", type="response"),
                                        method="trt.vs.ctrl", group="Control", adjust="FDR", level=0.90))) 
df_posthoc0$p.value=df_posthoc0$p.value/2 # we use 1-sided FDR corrected p values
df_posthoc0$ratio.lower95CB=df_posthoc0$ratio+1.96*df_posthoc0$SE
df_posthoc0$ratio.upper95CB=df_posthoc0$ratio-1.96*df_posthoc0$SE
df_posthoc0=df_posthoc0[,-c(3:6)]
df_posthoc0[,c(2:4)]=1/df_posthoc0[,c(2:4)] # to express reverse contrast compared to default
df_posthoc0$contrast=gsub(" / Control", "", df_posthoc0$contrast)
df_posthoc0

#   contrast    ratio ratio.lower95CB ratio.upper95CB
# 1    n-C25 1.087248       0.6914370        2.542955
# 2      QMP 1.017905       0.6474836        2.378803



### B1. ANALYSIS OF NR OF EGGS LAID BY WORKERS (FIG. 1A) ####
# random slope or random intercept models with or without average worker body size included as extra covariate:
fit2A=glmer(Nr.worker.eggs~Treatment+(1|Colony/Nest)+scale(Total.workers.7days)+scale(Mortality), 
            family = poisson (link=log), data = nestdata, 
            control = glmerControl(check.conv.singular = .makeCC(action="warning", tol = 1e-10)))
fit2B=glmer(Nr.worker.eggs~Treatment+(1+Treatment|Colony)+(1|Nest)+scale(Total.workers.7days)+scale(Mortality), 
            family = poisson (link=log), data = nestdata, 
            control = glmerControl(check.conv.grad = .makeCC("warning", tol = 1e-2)))
fit2C=glmer(Nr.worker.eggs~Treatment+(1|Colony/Nest)+scale(Total.workers.7days)+scale(Mortality)+scale(Avg.body.size), 
            family = poisson (link=log), data = nestdata, 
            control = glmerControl(check.conv.singular = .makeCC(action="warning", tol = 1e-10)))
fit2D=glmer(Nr.worker.eggs~Treatment+(1+Treatment|Colony)+(1|Nest)+scale(Total.workers.7days)+scale(Mortality)+scale(Avg.body.size), 
            family = poisson (link=log), data = nestdata, 
            control = glmerControl(check.conv.grad = .makeCC("ignore", tol = 1e-2)))
AIC(fit2A,fit2B,fit2C,fit2D) # fit2A gives best AIC (random intercept for colony & nest and
# not including average worker body size as additional covariate)
# note: analysis of model 2C gives qualitatively identical results, cf. summary(fit2C)
bestfit2=fit2A

sum2=tidy(bestfit2)
sum2$p.value=sum2$p.value/2 # we use 1-sided p values
as.data.frame(sum2)

#                         term      estimate  std.error  statistic       p.value  group
# 1                (Intercept)  3.425631e+00 0.14488272 23.6441647 6.776446e-124  fixed
# 2             Treatmentn-C25 -3.732117e-01 0.20985790 -1.7784018  3.766894e-02  fixed
# 3               TreatmentQMP -6.196929e-01 0.20945881 -2.9585429  1.545486e-03  fixed
# 4 scale(Total.workers.7days) -3.291881e-02 0.08817686 -0.3733271  3.544525e-01  fixed
# 5           scale(Mortality)  4.055467e-02 0.08696546  0.4663307  3.204894e-01  fixed
# 6        sd_(Intercept).Nest  4.080398e-01         NA         NA            NA   Nest
# 7      sd_(Intercept).Colony  2.272369e-05         NA         NA            NA Colony

# Control-treatment posthoc comparisons with FDR corrected p values
emmeans2=summary(emmeans(bestfit2,"Treatment", contr = "trt.vs.ctrl", adjust="fdr"),
                 type="response", level=0.9) # FDR posthoc comparisons, using 95% confidence bounds (=1-sided confidence intervals)
emmeans2$contrasts$p.value=emmeans2$contrasts$p.value/2 # we used 1 sided p values
emmeans2

# $emmeans
# Treatment rate   SE  df asymp.LCL asymp.UCL
# Control   30.7 4.45 Inf      24.2      39.0
# n-C25     21.2 3.14 Inf      16.6      27.0
# QMP       16.5 2.50 Inf      12.9      21.2
# 
# Confidence level used: 0.9 
# Intervals are back-transformed from the log scale 
# 
# $contrasts
# contrast        ratio    SE  df z.ratio p.value
# n-C25 / Control 0.689 0.144 Inf -1.778  0.0377 
# QMP / Control   0.538 0.113 Inf -2.959  0.0031 
# 
# P value adjustment: fdr method for 2 tests 
# Tests are performed on the log scale

# Effect sizes in terms of fold reduction in nr of eggs laid by workers in treatment vs control groups with 95% confidence bounds
df_posthoc2=data.frame(summary(contrast(emmeans(bestfit2,"Treatment", type="response"),
                                        method="trt.vs.ctrl", group="Control", adjust="FDR", level=0.90))) 
df_posthoc2$p.value=df_posthoc2$p.value/2 # we use 1-sided FDR corrected p values
df_posthoc2$ratio.lower95CB=df_posthoc2$ratio+1.96*df_posthoc2$SE
df_posthoc2$ratio.upper95CB=df_posthoc2$ratio-1.96*df_posthoc2$SE
df_posthoc2=df_posthoc2[,-c(3:6)]
df_posthoc2[,c(2:4)]=1/df_posthoc2[,c(2:4)] # to express reverse contrast compared to default
df_posthoc2$contrast=gsub(" / Control", "", df_posthoc2$contrast)
df_posthoc2

#   contrast    ratio ratio.lower95CB ratio.upper95CB
# 1    n-C25 1.452392 1.029101        2.467207
# 2      QMP 1.858357 1.317480        3.152639


# effect plot calculated using emmeans & plotted with ggplot2 (FIG. 1A)
df2=data.frame(emmeans2$emmeans)
df2$cols=rep(rgb(0,0,0),length(levels(nestdata$Treatment)))
df2$fills=c("grey75","steelblue","firebrick3")
ggplot(df2, aes(x=Treatment, y=rate, ymin=asymp.LCL, ymax=asymp.UCL, colour=I(cols), fill=I(fills))) + 
  geom_col(width=0.5, colour=NA) +
  geom_errorbar(aes(ymin=asymp.LCL, ymax=asymp.UCL), width=0, size=0.5) +
  theme_few(base_size=14) +
  labs(y = "Eggs laid by workers (avg.)", x="Treatment") +
  theme(legend.position="none", axis.line=element_line(colour="black"),
        panel.background=element_rect(fill="white"), axis.text=element_text(colour="black")) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
graph2ppt(file="Fig1A.pptx",width=5,aspectr=sqrt(2))



### B2. ANALYSIS OF AGGRESSION DATA (FIG. 1C) ####

set_treatment_contrasts() 
# random slope or random intercept models
fit3A=glmer(cbind(Nr.aggressive,Nr.non.aggressive) ~ 
              Treatment+(1|Colony/Nest)+scale(Total.workers.7days)+scale(Mortality), 
            family=binomial(link=logit), data=nestdata)
fit3B=glmer(cbind(Nr.aggressive,Nr.non.aggressive) ~ 
              Treatment+(1+Treatment|Colony)+(1|Nest)+scale(Total.workers.7days)+scale(Mortality), 
            family=binomial(link=logit), data=nestdata,
            control = glmerControl(optimizer = "Nelder_Mead"))
AIC(fit3A,fit3B) # fit3A gives best AIC (random intercepts for colony & nest)
bestfit3=fit3A
sum3=tidy(bestfit3)
sum3$p.value=sum3$p.value/2 # we use 1-sided p values
sum3$p.value[3]=1-sum3$p.value[3] # since QMP effect was opposite to a priori hypothesis
as.data.frame(sum3)

#                         term    estimate  std.error   statistic      p.value  group
# 1                (Intercept) -1.70609580 0.16684203 -10.2258155 7.597181e-25  fixed
# 2             Treatmentn-C25 -0.45415923 0.19408664  -2.3399819 9.642337e-03  fixed
# 3               TreatmentQMP  0.04713116 0.18410529   0.2560011 6.010250e-01  fixed
# 4 scale(Total.workers.7days) -0.18186554 0.09670502  -1.8806216 3.001171e-02  fixed
# 5           scale(Mortality) -0.03310080 0.09040253  -0.3661491 3.571269e-01  fixed
# 6        sd_(Intercept).Nest  0.34201101         NA          NA           NA   Nest
# 7      sd_(Intercept).Colony  0.31631896         NA          NA           NA Colony

# Effect plot significant effect of worker number using effects package (ESM)
plot(effect('scale(Total.workers.7days)', mod=bestfit3, confint=list(level=0.9)),
     type="response",ylab="Aggressive behaviour (%)", main = "") 
graph2ppt(file="FigS2.pptx",width=5,aspectr=sqrt(2))

# Control-treatment posthoc comparisons with FDR corrected p values
emmeans3=summary(emmeans(bestfit3,"Treatment", contr = "trt.vs.ctrl", adjust="fdr"),
                 type="response", level=0.90) # FDR posthoc comparisons, using 95% confidence bounds (=1-sided confidence intervals)
emmeans3$contrasts$p.value=emmeans3$contrasts$p.value/2 # we use 1 sided p vals
emmeans3$contrasts$p.value[2]=1-emmeans3$contrasts$p.value[2]
emmeans3

# $emmeans
# Treatment  prob     SE  df asymp.LCL asymp.UCL
# Control   0.154 0.0217 Inf    0.1213     0.193
# n-C25     0.103 0.0158 Inf    0.0802     0.132
# QMP       0.160 0.0216 Inf    0.1274     0.199
# 
# Confidence level used: 0.9 
# Intervals are back-transformed from the logit scale 
# 
# $contrasts
# contrast        odds.ratio    SE  df z.ratio p.value
# n-C25 / Control      0.635 0.123 Inf -2.340  0.0193 
# QMP / Control        1.048 0.193 Inf  0.256  0.6010 
# 
# P value adjustment: fdr method for 2 tests 
# Tests are performed on the log odds ratio scale 

# Effect sizes in terms of fold reduction in aggression in treatment vs control groups with 95% confidence bounds
df_posthoc3=data.frame(summary(contrast(emmeans(bestfit3,"Treatment", type="response"),
                                        method="trt.vs.ctrl", group="Control", adjust="FDR", level=0.90))) 
df_posthoc3$odds.ratio.lower95CB=df_posthoc3$odds.ratio+1.96*df_posthoc3$SE
df_posthoc3$odds.ratio.upper95CB=df_posthoc3$odds.ratio-1.96*df_posthoc3$SE
df_posthoc3=df_posthoc3[,-c(3:6)]
df_posthoc3[,c(2:4)]=1/df_posthoc3[,c(2:4)] # to express reverse contrast compared to default
df_posthoc3$contrast=gsub(" / Control", "", df_posthoc3$contrast)
df_posthoc3

#   contrast odds.ratio  odds.ratio.lower95CB  odds.ratio.upper95CB
# 1    n-C25  1.5748487  1.1408559             2.541759
# 2      QMP  0.9539623  0.7010066             1.492540

# Effect plots
# using ggplot2 (FIG. 1C)
df3=data.frame(emmeans3$emmeans)
df3$cols=rep(rgb(0,0,0),length(levels(nestdata$Treatment)))
df3$fills=c("grey75","steelblue","firebrick3")
ggplot(df3, aes(x=Treatment, y=prob*100, ymin=asymp.LCL*100, ymax=asymp.UCL*100, colour=I(cols), fill=I(fills))) + 
  geom_col(width=0.5, colour=NA) +
  geom_errorbar(aes(ymin=asymp.LCL*100, ymax=asymp.UCL*100), width=0, size=0.5) +
  theme_few(base_size=14) +
  labs(y = "Aggressive behaviour (%)", x="Treatment") +
  theme(legend.position="none", axis.line=element_line(colour="black"),
        panel.background=element_rect(fill="white"), axis.text=element_text(colour="black")) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) 
graph2ppt(file="Fig1C.pptx",width=5,aspectr=sqrt(2))



##### C. ANALYSIS OF QUEEN DATA ####

queendata=read.csv("Queendata.csv")
queendata$Treatment=factor(queendata$Treatment, levels=c("Control", "n-C25", "QMP"))

# C1. ANALYSIS OF THE NR. OF VIABLE OOCYTES IN QUEENS (FIG. 2B) ####

# model with or without observation-level random effect included to deal with possible overdispersion
fit4A=glmer(Viable_oocytes~Treatment+(1|Queen_ID), data=queendata, family = poisson)
fit4B=glm(Viable_oocytes~Treatment, data=queendata, family = poisson)

AIC(fit4A,fit4B) # best model (lowest AIC) is fit4A, ie with overdispersion included
bestfit4=fit4A
sum4=tidy(bestfit4)
sum4$p.value=sum4$p.value/2 # we use 1-sided p values
as.data.frame(sum4)
#                      term    estimate std.error  statistic      p.value    group
# 1             (Intercept)  1.13090167 0.2532387  4.4657537 3.989367e-06    fixed
# 2          Treatmentn-C25 -0.08380108 0.3501800 -0.2393086 4.054331e-01    fixed
# 3            TreatmentQMP -1.02121265 0.4140656 -2.4663063 6.825726e-03    fixed
# 4 sd_(Intercept).Queen_ID  0.53813889        NA         NA           NA Queen_ID

# Control-treatment posthoc comparisons with FDR corrected p values
emmeans4=summary(emmeans(bestfit4,"Treatment", contr = "trt.vs.ctrl", adjust="fdr"),type="response",level=0.90) # FDR posthoc comparisons, 95% confidence bounds
emmeans4$contrasts$p.value=emmeans4$contrasts$p.value/2 # we use 1 sided p values
emmeans4

# $emmeans
# Treatment rate    SE  df asymp.LCL asymp.UCL
# Control   3.10 0.785 Inf      2.04      4.70
# n-C25     2.85 0.750 Inf      1.85      4.39
# QMP       1.12 0.388 Inf      0.63      1.98
# 
# Confidence level used: 0.9 
# Intervals are back-transformed from the log scale 
# 
# $contrasts
# contrast        ratio    SE  df z.ratio p.value
# n-C25 / Control  0.92 0.322 Inf -0.239  0.4054 
# QMP / Control    0.36 0.149 Inf -2.466  0.0137 
# 
# P value adjustment: fdr method for 2 tests 
# Tests are performed on the log scale 

# Effect sizes in terms of fold reduction in nr of viable oocytes visible in queen ovaries 
# in treatment vs control groups with 95% confidence bounds
df_posthoc4=data.frame(summary(contrast(emmeans(bestfit4,"Treatment", type="response"),
                                        method="trt.vs.ctrl", group="Control", adjust="FDR", level=0.90))) 
df_posthoc4$ratio.lower95CB=df_posthoc4$ratio+1.96*df_posthoc4$SE
df_posthoc4$ratio.upper95CB=df_posthoc4$ratio-1.96*df_posthoc4$SE
df_posthoc4=df_posthoc4[,-c(3:6)]
df_posthoc4[,c(2:4)]=1/df_posthoc4[,c(2:4)] # to express reverse contrast compared to default
df_posthoc4$contrast=gsub(" / Control", "", df_posthoc4$contrast)
df_posthoc4

#   contrast    ratio ratio.lower95CB ratio.upper95CB
# 1    n-C25 1.087413        0.644831        3.466992
# 2      QMP 2.776560        1.532683       14.735122


# effect plot (FIG. 2A)
df4=data.frame(emmeans4$emmeans)
df4$cols=rep(rgb(0,0,0),length(levels(queendata$Treatment)))
df4$fills=c("grey75","steelblue","firebrick3")
ggplot(df4, aes(x=Treatment, y=rate, ymin=asymp.LCL, ymax=asymp.UCL, colour=I(cols), fill=I(fills))) + 
  geom_col(width=0.5, colour=NA) +
  geom_errorbar(aes(ymin=asymp.LCL, ymax=asymp.UCL), width=0, size=0.5) +
  theme_few(base_size=14) +
  labs(y = "Viable mature oocytes (avg.)", x="Treatment") +
  theme(legend.position="none", axis.line=element_line(colour="black"),
        panel.background=element_rect(fill="white"), axis.text=element_text(colour="black")) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
graph2ppt(file="Fig2B.pptx",width=5,aspectr=sqrt(2))


### C2. ANALYSIS OF THE NR. OF EGGS LAID BY QUEENS (FIG. 2A) ####

# model with or without observation-level random effect included to deal with possible overdispersion
fit5A=glmer(Total_eggs_laid~Treatment+(1|Queen_ID), data=queendata, family = poisson)
fit5B=glm(Total_eggs_laid~Treatment, data=queendata, family = poisson)

AIC(fit5A,fit5B) # best model (lowest AIC) is fit5A
bestfit5=fit5A
sum5=tidy(bestfit5)
sum5$p.value=sum5$p.value/2 # we use 1-sided p values
as.data.frame(sum5)
#                      term   estimate std.error  statistic       p.value    group
# 1             (Intercept)  3.6428171 0.1471964 24.7480096 1.628549e-135    fixed
# 2          Treatmentn-C25 -0.1088504 0.2087735 -0.5213803  3.010509e-01    fixed
# 3            TreatmentQMP -0.4964961 0.2120406 -2.3415145  9.602839e-03    fixed
# 4 sd_(Intercept).Queen_ID  0.4346437        NA         NA            NA Queen_ID

# Control-treatment posthoc comparisons with FDR corrected p values
emmeans5=summary(emmeans(bestfit5,"Treatment", contr = "trt.vs.ctrl", adjust="fdr"),type="response",level=0.90) # FDR posthoc comparisons, 95% confidence bounds
emmeans5$contrasts$p.value=emmeans5$contrasts$p.value/2 # we use 1 sided p values
emmeans5
# $emmeans
# Treatment rate   SE  df asymp.LCL asymp.UCL
# Control   38.2 5.62 Inf      30.0      48.7
# n-C25     34.3 5.09 Inf      26.8      43.7
# QMP       23.3 3.56 Inf      18.1      29.9
# 
# Confidence level used: 0.9 
# Intervals are back-transformed from the log scale 
# 
# $contrasts
# contrast        ratio    SE  df z.ratio p.value
# n-C25 / Control 0.897 0.187 Inf -0.521  0.3011 
# QMP / Control   0.609 0.129 Inf -2.342  0.0192 
# 
# P value adjustment: fdr method for 2 tests 
# Tests are performed on the log scale

# Effect sizes in terms of fold reduction in nr of eggs laid by queen in control vs treatment groups with 95% confidence bounds
df_posthoc5=data.frame(summary(contrast(emmeans(bestfit5,"Treatment", type="response"),
                                        method="trt.vs.ctrl", group="Control", adjust="FDR", level=0.90))) 
df_posthoc5$ratio.lower95CB=df_posthoc5$ratio+1.96*df_posthoc5$SE
df_posthoc5$ratio.upper95CB=df_posthoc5$ratio-1.96*df_posthoc5$SE
df_posthoc5=df_posthoc5[,-c(3:6)]
df_posthoc5[,c(2:4)]=1/df_posthoc5[,c(2:4)] # to express reverse contrast compared to default
df_posthoc5$contrast=gsub(" / Control", "", df_posthoc5$contrast)
df_posthoc5

#   contrast    ratio ratio.lower95CB ratio.upper95CB
# 1    n-C25 1.114996       0.7912281        1.887251
# 2      QMP 1.642954       1.1606068        2.811350


# effect plot
df5=data.frame(emmeans5$emmeans)
df5$cols=rep(rgb(0,0,0),length(levels(queendata$Treatment)))
df5$fills=c("grey75","steelblue","firebrick3")
ggplot(df5, aes(x=Treatment, y=rate, ymin=asymp.LCL, ymax=asymp.UCL, colour=I(cols), fill=I(fills))) + 
  geom_col(width=0.5, colour=NA) +
  geom_errorbar(aes(ymin=asymp.LCL, ymax=asymp.UCL), width=0, size=0.5) +
  theme_few(base_size=14) +
  labs(y = "Eggs laid by queens (avg.)", x="Treatment") +
  theme(legend.position="none", axis.line=element_line(colour="black"),
        panel.background=element_rect(fill="white"), axis.text=element_text(colour="black")) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
graph2ppt(file="Fig2A.pptx",width=5,aspectr=sqrt(2))



# D. FISHER'S METHOD TO COMBINE p VALUES OF SEPARATE ANALYSES ####

# function to combine p values of multiple independent tests of the same hypothesis using Fisher's method:
pfishersMethod = function(pvals) pchisq(-2*sum(log(pvals)),df=2*length(pvals),lower.tail=F)

# D1. SIGNIFICANCE OF INHIBITION OF REPRODUCTION (EGG-LAYING) CAUSED BY QMP ACROSS BOTH CASTES ####

# One-sided FDR corrected p value for effect of QMP on worker egg laying and queen egg laying obtained above were 
# 0.0031 (cf. emmeans2 above) and 0.0192 (cf. emmeans5 above)
# combining these the overall one-sided p value for the inhibition in egg laying across both castes is
pfishersMethod(c(0.0031, 0.0192)) # p=0.00064, df = 4, or 2-sided this would be p=0.0013
-2*sum(log(c(0.0031, 0.0192))) # Chi squared=19.46

# D2. SIGNIFICANCE OF INHIBITION OF WORKERY OVARY ACTIVATION INDUCED BY PENTACOSANE ACROSS THIS & two EARLIER STUDIES ####

# For pentacosane, earlier studies reported 1-sided p values for the effect on worker ovary inhibition of
# p=0.001 (Van Oystaeyen et al. 2014, outcome variable=ovary regression) and p=0.001 (Holman 2014, outcome variable=nr of visible oocytes in worker ovaries), which together with the p=0.061 observed in our study (cf. emmeans1 above) results in
# a combined overall one-sided p value for pentacosane-induced worker ovary inhibition of
pfishersMethod(c(0.001, 0.001, 0.061)) # p=9*10^-6, df = 6, or 2-sided this would be p=2*10^-5
-2*sum(log(c(0.001, 0.001, 0.061))) # Chi squared=33.22
