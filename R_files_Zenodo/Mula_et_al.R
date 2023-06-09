# Mulá, Thorogood & Hämäläinen 2022, Behavioral Ecology
# Social information use about novel aposematic prey depends on the intensity of the observed cue

# packages needed
library(lme4)
library(DHARMa)
library(ggplot2)
library(Rmisc)
library(MASS)

# Save data as csv files and set a working directory where the data files are saved

#read the data file      
data<-read.csv("Mula et al. data 2022.csv",header=T, sep  = ";")
str(data) 

# change bird and demonstrator ids and experiment duration to be categorical variables
data$BirdID<-as.factor(data$BirdID)
data$DemID<-as.factor(data$DemID)
data$Duration.days<-as.factor(data$Duration.days)


########################

# FIRST PREY CHOICE IN THE EXPERIMENT

#create a dataset of first foraging trial only
data2<-subset(data, Trial=="1") 

table(data2$First.choice) # first choices
binom.test(31, 45, p = 0.5,	conf.level = 0.95) 

# number of successes = 31, number of trials = 45, p-value = 0.01609
# alternative hypothesis: true probability of success is not equal to 0.5
# 95 percent confidence interval:
#  0.5335090 0.8183412
# sample estimates:
# probability of success 
# 0.6888889 

# birds chose a square as a first prey significantly more often

table(data2$First.choice, data2$Treatment) # first choices in each treatment

# Chi-squared test to test whether first choices differ among treatments
chisq.test(data2$First.choice, data2$Treatment,correct=FALSE) 

#	Pearson's Chi-squared test
#data:  data2$First.choice and data2$Treatment
#X-squared = 0.20737, df = 2, p-value = 0.9015

# no evidence that treatment influences first choice

############################

# TIME BEFORE ATTACKING THE FIRST PREY 

# GLM with a negative binomial error distribution (response variable time before attack)

# including all observations, treatment as a fixed effect, demonstrator id as a random effect
time_model<-glmer.nb(Attack.time~Treatment+(1|DemID), data=data2)
summary(time_model)
# the model indicates that birds in the palatable treatment were slower to attack the prey

# the data includes two outliers (both in the palatable treatment) that attacked the prey after 923 and 2768s

# exclude these observations
data3<-subset(data, Attack.time<800)
hist(data3$Attack.time)
plot(data3$Attack.time~data3$Treatment)

max(data3$Attack.time)
min(data3$Attack.time)
mean(data3$Attack.time)
# after excluding two outliers, attack time varies from 5-325s, mean = 61s

# run the same model excluding outliers
time_model2<-glmer.nb(Attack.time~Treatment+(1|DemID), data=data3)
summary(time_model2)

#Random effects:
# Groups Name        Variance Std.Dev.
# DemID  (Intercept) 0.02517  0.1586  
# Number of obs: 43, groups:  DemID, 5

# Fixed effects:
#                  Estimate  Std. Error z value   Pr(>|z|)    
# (Intercept)      4.21348    0.26774   15.737    <2e-16 ***
# Treatmentstrong -0.33117    0.35501   -0.933     0.351
# Treatmentweak   -0.04987    0.35449   -0.141     0.888    

# by default, intercept is the palatable treatment, this can be changed with relevel function
data3$Treatment<-relevel(data3$Treatment, ref="weak") # intercept weak treatment 
data3$Treatment<-relevel(data3$Treatment, ref="strong") # intercept strong treatment
data3$Treatment<-relevel(data3$Treatment, ref="palatable") # intercept palatable treatment

# the model indicates that there are no differences in attack times among the treatments

# plot residuals to investigate model fit
simulationOutput <- simulateResiduals(fittedModel = time_model2, plot = F)
plot(simulationOutput)

# Figure 3: differences in attack times (excluding 2 outliers)

ggplot(data3, aes(x=Treatment, y=Attack.time)) + 
  geom_violin(trim=FALSE, fill="gray")+
  geom_boxplot(width=0.1)+
  scale_x_discrete(breaks=c("palatable", "weak", "strong"),labels=c("Palatable", "Weak", "Strong"))+
  theme_bw() + theme(plot.background = element_blank(),panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())+
  ylab(expression(paste("Time to attack the first prey item (sec)")))+
  xlab("Treatment")+
  theme(axis.title.y = element_text(size = 12))+
  theme(axis.title.x = element_text(size = 12))+
  theme(axis.text = element_text(size = 12))+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.title.x = element_text(margin = margin(t = 20, r = 10, b = 0, l = 0)))+
  theme(legend.position="none")


#######################

# PREY CHOICES IN THE FIRST TRIAL

# GLM with a binomial error distribution (response variable number of squares and crosses consumed)

# treatment and age as fixed effects and demonstrator id as a random effect 
model1<-glmer(cbind(Squares,Crosses)~Treatment+(1|DemID)+Age,family=binomial, data=data2)
summary(model1) 

# Random effects:
# Groups Name        Variance Std.Dev.
# DemID  (Intercept) 0.03377  0.1838  
# Number of obs: 45, groups:  DemID, 5

# Fixed effects:
#                  Estimate Std. Error  z value   Pr(>|z|)    
# (Intercept)       0.6232     0.1882    3.311   0.000929 ***
# Treatmentstrong  -0.5801     0.1860   -3.118   0.001819 ** 
# Treatmentweak    -0.1788     0.1869   -0.957   0.338736    
# AgeJuvenile      -0.3250     0.1692   -1.921   0.054719

# change intercept with relevel function to compare treatments
data2$Treatment<-relevel(data2$Treatment, ref="weak") # intercept weak treatment 
data2$Treatment<-relevel(data2$Treatment, ref="strong") # intercept strong treatment
data2$Treatment<-relevel(data2$Treatment, ref="palatable") # intercept palatable treatment

# the model indicates that birds in the strong treatment attacked fewer aposematic prey compared to palatable and weak treatments 
# there is no difference between weak and palatable treatments

# plot residuals to investigate model fit
simulationOutput <- simulateResiduals(fittedModel = model1, plot = F)
plot(simulationOutput)

############################

# LEARNING ACROSS TWO TRIALS

# a model with treatment * trial number interaction, age and experiment duration (1 or 2 days) as fixed effects and demonstrator and observer ids as random effects
trials_model<-glmer(cbind(Squares,Crosses)~Treatment*Trial+Duration.days+Age+(1|BirdID)+(1|DemID),family=binomial, data=data)
summary(trials_model) # interaction not significant, no difference in learning rates among treatments 

# a model with experiment duration * trial number interaction, age and treatment as fixed effects and demonstrator and observer ids as random effects
trials_model2<-glmer(cbind(Squares,Crosses)~Treatment+Duration.days*Trial+Age+(1|BirdID)+(1|DemID),family=binomial, data=data)
summary(trials_model2) # interaction p = 0.076, close to significant

# a model with experiment duration * treatment interaction, age and trial number as fixed effects and demonstrator and observer ids as random effects
trials_model3<-glmer(cbind(Squares,Crosses)~Treatment*Duration.days+Trial+Age+(1|BirdID)+(1|DemID),family=binomial, data=data)
summary(trials_model3) # interaction not significant

# a model with no interactions, treatment, duration and age as fixed effects and demonstrator and observer ids as random effects
trials_model4<-glmer(cbind(Squares,Crosses)~Treatment+Duration.days+Trial+Age+(1|BirdID)+(1|DemID),family=binomial, data=data)
summary(trials_model4) 

# comparing models with interactions to the model without interactions
anova(trials_model, trials_model4) # treatment*trial number, interaction term not significant
anova(trials_model2, trials_model4) # duration*trial number, interaction term close to significant (p = 0.074) but not at alpha level 0.05
anova(trials_model3, trials_model4) # duration*treatment, interaction term not significant

# final model
trials_model4<-glmer(cbind(Squares,Crosses)~Treatment+Duration.days+Trial+Age+(1|BirdID)+(1|DemID),family=binomial, data=data)
summary(trials_model4) 

#Random effects:
#  Groups Name        Variance Std.Dev.
#  BirdID (Intercept) 0.08271  0.28759 
#  DemID  (Intercept) 0.00119  0.03449 
#  Number of obs: 90, groups:  BirdID, 45; DemID, 5

# Fixed effects:
#                   Estimate Std. Error  z value  Pr(>|z|)    
#  (Intercept)       1.1263     0.2354   4.785   1.71e-06 ***
#  Treatmentstrong  -0.4346     0.1751  -2.482   0.013049 *  
#  Treatmentweak    -0.1900     0.1728  -1.100   0.271451    
#  Duration.days2    0.6884     0.1783   3.861   0.000113 ***
#  Trial            -0.7779     0.1111  -7.001   2.55e-12 ***
#  AgeJuvenile      -0.1748     0.1485  -1.177   0.239281  

# change the intercept using relevel function
data$Treatment<-relevel(data$Treatment, ref="strong")
data$Treatment<-relevel(data$Treatment, ref="weak")
data$Treatment<-relevel(data$Treatment, ref="palatable")

# birds reduced consumption of aposematic prey across the two trials
# slow birds (duration 2 days) consumed more aposematic prey than fast birds

# plot residuals to investigate model fit
simulationOutput <- simulateResiduals(fittedModel = trials_model4, plot = F)
plot(simulationOutput)

# Figure 4: learning across two trials

data$Trial<-as.factor(data$Trial) # for the graph, change trial to be a categorical variable

# Calculate relative predation risk for aposematic prey by dividing the number of squares consumed by 8 
# (because birds consumed 16 prey in each trial, they were expected to consume 8 of each prey type if they chose randomly)
data$predation_risk<-data$Squares/8

# This summarizes the predation risk (+ sd, se and ci) in each treatment and trial 
pred_risk_trials<- summarySE(data, measurevar="predation_risk", groupvars=c("Treatment","Trial"), na.rm=T)
pred_risk_trials

ggplot(data, aes(x=Trial, y=predation_risk, shape=factor(Treatment, labels=c("Palatable", "Weak", "Strong")),group=Treatment))+ # x axis = trial, y axis = predation risk, plot each treatment separately
geom_point(position=position_jitterdodge(dodge.width=0.7), color="gray47", show.legend=F) +  geom_point(data=pred_risk_trials, aes(x = Trial, y = predation_risk), size=3.5, position=position_dodge(width=0.6))+ # plot individual datapoints
  geom_errorbar(data=pred_risk_trials, aes(ymin=predation_risk-se, ymax=predation_risk+se),width=0.1, position=position_dodge(width=0.6))+ # add mean and error bars
  scale_shape_manual(values=c(16,2,8))+ # define symbols for each treatment 
  geom_line(data=pred_risk_trials, aes(linetype=Treatment),position=position_dodge(width=0.6))+ # add lines
  scale_linetype_manual(values=c("dashed", "solid", "dotted"), guide='none')+ # define linetypes for each treatment 
  theme_bw() + theme(plot.background = element_blank(),panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())+ 
  scale_y_continuous(limits=c(0,2.0), breaks=c(0,0.2,0.4,0.6,0.8,1.0,1.2,1.4, 1.6, 1.8))+   
  scale_x_discrete(breaks=c("1", "2"),labels=c("Trial 1", "Trial 2"))+ 
  labs(x="", 
       y="Relative predation risk for aposematic prey")+ 
  theme(axis.text = element_text(size =11))+ 
  theme(axis.text.y = element_text(size =11))+ 
  theme(axis.title.y = element_text(size =12))+ 
  theme(legend.position = c(0.8, 0.88), legend.title = element_blank())+ 
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))+ 
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))+ 
  geom_hline(yintercept = 1, linetype=3) # add a horizontal line at y = 1 to illustrate predation risk of 1 (when birds consume same number of both prey items)
  

####################

# VARIATION IN BEAK WIPING

# read the data file (second sheet in Excel file)
data4<-read.csv("Beak wiping variation.csv",header=T, sep  = ";")
str(data4) 

# Testing the effect of age on beak wiping behaviour
# GLM with a negative binomial error distribution (response variable number of beak wipes)
beak_wiping_model<-glm.nb(beak.wipes~Age, data=data4)
summary(beak_wiping_model)

# Coefficients:
#               Estimate  Std. Error  z value   Pr(>|z|)    
# (Intercept)   3.64388    0.18873    19.308   <2e-16 ***
# AgeJuvenile  -0.01954    0.33146    -0.059    0.953    

# the model indicates that age does not influence beak wiping behavior

# plot residuals to investigate model fit
simulationOutput <- simulateResiduals(fittedModel = beak_wiping_model, plot = F)
plot(simulationOutput)

# plot variation in beak wiping
hist(data4$beak.wipes, breaks=18, main="", col="gray",
     xlab="Number of beak wipes in 60s", xlim=c(0,180), las=1,
     ylab="Number of birds", ylim=c(0,10))

####################

# Supplementary analysis: 

# did variation in demonstrators' responses in the strong treatment influence observers' behaviour?

# include only birds from strong treatment and choices in the first trial
strong_treatment<-subset(data, Treatment=="strong" & Trial=="1")

# correlation between the number of demonstrators' beak wipes and the number of aposematic prey observers attacked
cor.test(strong_treatment$Squares,strong_treatment$Dem.wipes, method = "pearson")
# no significant correlation, p = 0.5572

# Supplementary Fig. S1

ggplot(data = strong_treatment) + 
  geom_point(mapping = aes(x = Dem.wipes, y = Squares), size=2)+
  theme_bw() + theme(plot.background = element_blank(),panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())+
  scale_y_continuous(limits=c(1,14), breaks=c(2,4,6,8,10,12,14))+ 
  scale_x_continuous(limits=c(60,90), breaks=c(60,65,70,75,80,85,90))+ 
  labs(x="Number of demonstrator beak wipes", #  
       y="Aposematic prey attacked in the first trial")+
  theme(axis.title.y = element_text(size =12))+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))+
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))



