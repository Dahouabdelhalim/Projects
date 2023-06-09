# Analyses for: Distinct developmental trajectories for risky and impulsive decision-making in chimpanzees
# Rosati, Emery Thompson, Atencia, & Buckholtz
# Risky Choice Task


########### LOAD PACKAGES ################################################################ 

library(lattice) # plots
library(ggplot2) #plots
library(effects) #plots of effects
library(gridExtra) #combining plots
library(lme4) #install GLMM package
library(emmeans) #multiple comparisons
library(optimx) #for mixed model optimization
library(plyr) # split data dpply function
library(lmerTest) #tests for linear mixed effects models
library(cowplot) #plot grid to integrate graphs
library(dplyr) #organize data
library(car) #type 3 anova for omnibus tests
library(bbmle) #AIC comparison
library(readxl) #to open data
library(tidyr) #long to wide format
library(arm) #binned residual plot
library(lsr) #cohensD


########### LOAD DATA ###############################################################

# set screen print width
options(width = 120)

## set working directory

# Open risk task data
data <- read_excel("Rosati_etal_ChimpanzeeAdolescence-Data-RiskyChoice.xlsx")

#rename sex factors
data$Sex <- as.character(data$Sex)
data$Sex[data$Sex=="F"]<-"Female"
data$Sex[data$Sex=="M"]<-"Male"

#age cohorts
data$Cohort <-"Younger"
data$Cohort[data$Age>=15]<-"Adult"
data$Cohort <- as.factor(data$Cohort)
data$Cohort <- relevel(data$Cohort, ref = "Younger")

#rename outcome factors
data$Outcome <- as.character(data$Outcome)
data$Outcome[data$Outcome=="1"]<-"Good"
data$Outcome[data$Outcome=="0"]<-"Bad"
data$Outcome[data$Outcome=="2"]<-"Safe"

#create binary outcome factors for some analyses
data$Outcome_Binary<-"Negative"
data$Outcome_Binary[data$Outcome=="Good"| data$Outcome=="Safe"]<-"Positive"
data$Outcome_Binary[data$Outcome=="Bad"]<-"Negative"

#rename previous outcome factors
data$Prev_Outcome <- as.character(data$Prev_Outcome)
data$Prev_Outcome[data$Prev_Outcome=="1"]<-"Good"
data$Prev_Outcome[data$Prev_Outcome=="0"]<-"Bad"
data$Prev_Outcome[data$Prev_Outcome=="2"]<-"Safe"

#rename trial type factors
data$Type <- as.character(data$Type)
data$Type[data$Type=="E"]<-"Exposure"
data$Type[data$Type=="C"]<-"Test"

#make sure switch is coded as a number
data$Switch <- as.numeric(data$Switch)

#specify factors on task data
framefacs <- c("Subject","Sex", "Type", "Outcome", "Outcome_Binary", "Prev_Outcome")
data[, framefacs] <- lapply(data[, framefacs], as.factor)

#check variables
lapply(list(data), str)

#subset choice (test trial) data
test_data <- data[ which(data$Type=='Test'), ]



######### CREATE SUBJECT DATA ###############################################

#individual subj means for risk by previous choice
risk_mean<-as.data.frame(ddply(test_data,c("Subject", "Sex", "Age", "Cohort"),summarise,
               	Pref_score	= mean(Pref_Score, na.rm=TRUE),
               	Mean_Risk	= mean(Risky, na.rm=TRUE)))
risk_mean


#subj means for emotional responses to good outcome trials - using all data (exposure and choice trials)
emotion_good<-data[which(data$Outcome == "Good"),]
mean_emotion_good<-as.data.frame(ddply(emotion_good,c("Subject"),summarise,
               	Mean_Score_Good	= mean(Score, na.rm=TRUE)))
mean_emotion_good

#subj means for emotional responses to bad outcome trials - using all data (exposure and choice trials)
emotion_bad<-data[which(data$Outcome == "Bad"),]
mean_emotion_bad<-as.data.frame(ddply(emotion_bad,c("Subject"),summarise,
               	Mean_Score_Bad	= mean(Score, na.rm=TRUE)))
mean_emotion_bad

#subj means for switch responses to good outcome trials - using all data (exposure and choice trial)
switch_good<-test_data[which(test_data$Outcome == "Good"),]
mean_switch_good<-as.data.frame(ddply(switch_good,c("Subject"),summarise,
               	Mean_Switch_Good	= mean(Switch, na.rm=TRUE)))
mean_switch_good

#subj means for switch responses to bad outcome trials - using all data (exposure and choice trial)
switch_bad<-test_data[which(test_data$Outcome == "Bad"),]
mean_switch_bad<-as.data.frame(ddply(switch_bad,c("Subject"),summarise,
               	Mean_Switch_Bad	= mean(Switch, na.rm=TRUE)))
mean_switch_bad

#create merged data
subj_data<-merge(risk_mean, mean_emotion_good, all = TRUE, by = c("Subject"))
subj_data<-merge(subj_data, mean_emotion_bad, all = TRUE, by = c("Subject"))
subj_data<-merge(subj_data, mean_switch_good, all = TRUE, by = c("Subject"))
subj_data<-merge(subj_data, mean_switch_bad, all = TRUE, by = c("Subject"))

#calculate an individual affect difference score - indexing difference in affect score to bad v. good risk outcomes
subj_data$Affect_DiffScore<-(subj_data$Mean_Score_Bad-subj_data$Mean_Score_Good)

#calculate an individual switching difference score - indexing difference in switching to bad v. good risk outcomes
subj_data$Switch_DiffScore<-(subj_data$Mean_Switch_Bad-subj_data$Mean_Switch_Good)

#check variables
lapply(list(subj_data), str)

#save for cross-task analysis

# check age ranges
min(subj_data$Age) #6
max(subj_data$Age) #25

# sample size by sex
length(which(subj_data$Sex == 'Female')) # 19
length(which(subj_data$Sex == 'Male')) # 21



######### RISKY CHOICE ANALYSIS ###############################################

#Histogram to check distributions of mean risky choice
ggplot(subj_data,aes(x=Mean_Risk, fill=Cohort))+ geom_histogram(position="dodge",binwidth=.1)
#ok

#overall means by cohort for risk
#subj means for risk by previous choice
risk_overall_mean<-as.data.frame(ddply(subj_data,c("Cohort"),summarise,
               	N		= length(!is.na(Mean_Risk)),
               	mean	= mean(Mean_Risk, na.rm=TRUE),
               	sd  	= sd(Mean_Risk, na.rm=TRUE),
               	se   	= sd / sqrt(N)))
risk_overall_mean

#t-test: initial comparison of adults and younger
t.test(subj_data$Mean_Risk ~ subj_data$Cohort, var.equal = TRUE)
#younger more risk prone    
#calculate effect size
cohensD(subj_data$Mean_Risk ~ subj_data$Cohort) 
# medium effect size 
      
#GLMMS of risky choice
#Base model of choice trials: account for sex, trial number (1-20 within choice trials), food pref score
Risk_model1 = glmer(Risky ~ (1|Subject) + Sex + Trial_type + Pref_Score, 
		family=binomial(link = "logit"),
		data=test_data)
summary(Risk_model1)
#positive effect of pref score in base model

#Add in prev outcome
Risk_model2 = glmer(Risky ~ (1|Subject) + Sex + Trial_type + Pref_Score + Prev_Outcome, 
		family=binomial(link = "logit"),
		data=test_data)
summary(Risk_model2)
anova(Risk_model1, Risk_model2)
#trend to improves fit 

emmeans(Risk_model2, pairwise ~ Prev_Outcome)
#trend in bad < good p = 0.065

#Add in age (continuous)
Risk_model3 = glmer(Risky ~ (1|Subject) + Sex + Trial_type + Pref_Score + Prev_Outcome + Age, 
		family=binomial(link = "logit"),
		data=test_data)
summary(Risk_model3)
anova(Risk_model2, Risk_model3)
#improves fit 

#Add in age (continuous) X outcome interaction
Risk_model4 = glmer(Risky ~ (1|Subject) + Sex + Trial_type + Pref_Score + Prev_Outcome*Age, 
		family=binomial(link = "logit"),
		data=test_data,
		control = glmerControl(optimizer = "optimx", optCtrl=list(method = "nlminb")))
summary(Risk_model4)
anova(Risk_model3, Risk_model4)
#no improvement

#Add in age (cohort): extra check for model 3
Risk_model3b = glmer(Risky ~ (1|Subject) + Sex + Trial_type + Pref_Score + Prev_Outcome + Cohort, 
		family=binomial(link = "logit"),
		data=test_data)
summary(Risk_model3b)
anova(Risk_model2, Risk_model3b)
#improves fit as when continuous age is used

#Add in age (cohort) interaction: extra check for model 4
Risk_model4b = glmer(Risky ~ (1|Subject) + Sex + Trial_type + Pref_Score + Prev_Outcome*Cohort  , 
		family=binomial(link = "logit"),
		data=test_data)
summary(Risk_model3b)
anova(Risk_model3b, Risk_model4b)
#no improvement as when continuous age is used

#omnibus test on full model
Anova(Risk_model4, type = "III")
#main effect of effect of age, trends for score and outcome

#full AIC comparison
AICtab(Risk_model1, Risk_model2, Risk_model3, Risk_model4, base = T, weights = T)
#3 is best by AIC  as well

#full AIC comparison using cohort
AICtab(Risk_model1, Risk_model2, Risk_model3b, Risk_model4b, base = T, weights = T)
#3b is best by AIC  as well

#plot effects
quartz.options(width=10, height=6)
plot(allEffects(Risk_model3))

#look at binned residual plot for logistic regression to assess assumptions
binnedplot(fitted(Risk_model3),residuals(Risk_model3, type="response"), cex.pts=1, col.int="black")
#ok

#save model 3 parameters
write.table(summary(Risk_model3)$coefficients, sep="\\t")



############### AFFECT SCORE ANALYSIS ###############################################

#subj means for emotion scores
subjmeans_score<-as.data.frame(ddply(data,c("Subject", "Sex", "Age", "Outcome"),summarise,
               	N		= length(!is.na(Score)),
               	mean	= mean(Score, na.rm=TRUE),
               	sd  	= sd(Score, na.rm=TRUE),
               	se   	= sd / sqrt(N)))
subjmeans_score
#check that matches summary from above

#overall means for scores
means_score<-as.data.frame(ddply(subjmeans_score,c("Outcome"),summarise,
               	N		= length(!is.na(mean)),
               	Mean	= mean(mean, na.rm=TRUE),
               	SD  	= sd(mean, na.rm=TRUE),
               	SE   	= SD / sqrt(N)))
means_score

#Histogram to check distributions of average scores by outcome
ggplot(subjmeans_score,aes(x=mean, fill=Outcome))+ geom_histogram(position="dodge",binwidth=.1)
#higher scores for bad outcomes
#use linear models as scores reflect intensity of response; cannot use alternative because of 0's

#LMMS
#Base model: account for trial number (including exposure trials, 1-30), sex, and trial type (exposure or choice)
Affect_model1 = lmer(Score ~ (1|Subject) + Sex +  Trial + Type, 
		data=data)
summary(Affect_model1)

#Add in outcome
Affect_model2 = lmer(Score ~ (1|Subject) + Sex + Trial +  Type + Outcome, 
		data=data)
summary(Affect_model2)
anova(Affect_model1, Affect_model2)
#improves fit

emmeans(Affect_model2, pairwise ~ Outcome) 
# Bad responses > good
# bad responses > safe

#Add in Age
Affect_model3 = lmer(Score ~ (1|Subject) + Sex + Trial + Type + Outcome + Age, 
		data=data)
summary(Affect_model3)
anova(Affect_model2, Affect_model3)
#no improvement

#Add in Age X outcome interaction
Affect_model4 = lmer(Score ~ (1|Subject) + Sex + Trial_type + Type + Outcome*Age, 
		data=data)
summary(Affect_model4)
anova(Affect_model2, Affect_model4)
#no improvement compared to 2nd model

#plot effects
quartz.options(width=10, height=6)
plot(allEffects(Affect_model3))

#look at residuals
qqnorm(resid(Affect_model2),main='Normality Check for Residuals')
qqline(resid(Affect_model2))
#some skew

#omnibus test
Anova(Affect_model3, type = "III")
#effect of outcome only, no age effect

#AIC comparison
AICtab(Affect_model1, Affect_model2, Affect_model3, Affect_model4, base = T, weights = T)
#2 is best by AIC as well

#save model 3 parameters to account for lack of age effect
write.table(summary(Affect_model3)$coefficients,  sep="\\t")




############## SWITCHING ANALYSIS ###############################################

#subj means for swtiching by outcome (good, bad, safe)
subjmeans_switching<-as.data.frame(ddply(test_data,c("Subject", "Sex", "Age", "Outcome"),summarise,
               	N		= length(!is.na(Switch)),
               	mean	= mean(Switch, na.rm=TRUE),
               	sd  	= sd(Switch, na.rm=TRUE),
               	se   	= sd / sqrt(N)))
subjmeans_switching
#check this matches summary above

#overall means for switching
means_switching<-as.data.frame(ddply(subjmeans_switching,c("Outcome"),summarise,
               	N		= length(!is.na(mean)),
               	Mean	= mean(mean, na.rm=TRUE),
               	SD  	= sd(mean, na.rm=TRUE),
               	SE   	= SD / sqrt(N)))
means_switching
#switching is rare: use binary outcomes collapsing good/safe for this analysis

#GLMMs looking at switching
#Base model: account for trial within type (1-20, only choice trials) and sex
Switch_model1 = glmer(Switch ~ (1|Subject) + Trial_type + Sex, 
		family=binomial(link = "logit"),	
		data=test_data)
summary(Switch_model1)

#Add in outcome: use binary outcome here
Switch_model2 = glmer(Switch ~ (1|Subject) + Trial + Sex + Outcome_Binary, 
		family=binomial(link = "logit"),	
		data=test_data)
summary(Switch_model2)
anova(Switch_model1, Switch_model2)
#improves fit

#Add in Age
Switch_model3 = glmer(Switch ~ (1|Subject) + Trial + Sex + Outcome_Binary + Age, 
		family=binomial(link = "logit"),	
		data=test_data)
summary(Switch_model3)
anova(Switch_model2, Switch_model3)
#No improvement

#Add in Age interaction
Switch_model4 = glmer(Switch ~ (1|Subject) + Trial + Sex + Outcome_Binary*Age, 
		family=binomial(link = "logit"),	
		data=test_data)
summary(Switch_model4)
anova(Switch_model2, Switch_model4)
#No improvement

#plot effects
quartz.options(width=10, height=6)
plot(allEffects(Switch_model2))

#look at binned residual plot 
binnedplot(fitted(Switch_model2),residuals(Risk_model3, type="response"), cex.pts=1, col.int="black")
#some values outside, but this is appropriate structure for this response outcome

#omnibus test
Anova(Switch_model3, type = "III")
#only effect of binary outcome, no age

#AIC comparison
AICtab(Switch_model1, Switch_model2, Switch_model3, Switch_model4, base = T, weights = T)
#2 is best by AIC  as well

#save model 3 parameters with age
write.table(summary(Switch_model3)$coefficients, sep="\\t")




########### RISK PLOTS ###############################################


#Convert parameter estimates for risky choice preferences
Risk_age_effects <- allEffects(Risk_model3, xlevels = list(Age = seq(6, 25, by = 1)))
Risk_effects_df <- as.data.frame(Risk_age_effects[["Age"]])

Risk_age_effects <-ggplot(Risk_effects_df, aes(x = Age, y = fit)) +
    geom_line(color = "royalblue4") +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, fill = "royalblue4") +
    geom_point(data = subj_data, aes(x = Age, y = Mean_Risk), position = "jitter", inherit.aes = FALSE, color = "grey50") +
    theme_bw() +
    ggtitle("(a) Risk Preferences by Age")+
    scale_x_continuous(name="Age (Years)", breaks=seq(6, 25, 2)) +
    expand_limits(y=c(0,1)) +
    theme(aspect.ratio=1/1) +
    scale_y_continuous(name="Predicted value (Choose Risk)") +
    theme(plot.title = element_text(hjust = 0.5),
    	axis.title.x = element_text(face="bold", size=12), 
    	axis.title.y = element_text(face="bold", size = 12), 				
    	axis.text.y = element_text(size=10), 
    	axis.text.x = element_text(size=10), 
    	legend.text = element_text(size=10),
    	legend.title = element_text(face="bold", size=12))
Risk_age_effects

#subj means for risk by previous choice
subjmeans_prevchoice_risk<-as.data.frame(ddply(test_data,c("Subject", "Sex", "Age", "Cohort", "Prev_Outcome"),summarise,
               	N		= length(!is.na(Risky)),
               	mean	= mean(Risky, na.rm=TRUE),
               	sd  	= sd(Risky, na.rm=TRUE),
               	se   	= sd / sqrt(N)))
subjmeans_prevchoice_risk

#Boxplot
subjmeans_prevchoice_risk$Prev_Outcome<-ordered(subjmeans_prevchoice_risk$Prev_Outcome, levels = c("Good", "Safe", "Bad"))
risk_plot <- ggplot(subjmeans_prevchoice_risk, aes(x=Cohort, y=mean, fill = Prev_Outcome)) + 
  		geom_boxplot() +
  		stat_summary(fun=mean, geom="point", position = position_dodge(width=0.75), shape=5, size=4) +  	
    		theme_bw() +
		theme(aspect.ratio=1/1) +
		ggtitle("(b) Prior Outcomes") +
		scale_fill_manual(values=c("lightskyblue1", "dodgerblue", "dodgerblue4"), name="Outcome") +    		scale_x_discrete(name="") +
    		scale_y_continuous(name="Mean Choices for Risk", breaks=seq(0,1,.25)) +   		
    	theme(  plot.title = element_text(hjust = 0.5),
	   	 		axis.text.x = element_text(face="bold", size = 12), 
   				axis.title.y = element_text(face="bold", size = 12), 
   				axis.text.y = element_text(size=12),
   				legend.text = element_text(size=10),
    			legend.title = element_text(face="bold", size=12))
risk_plot

#subj means for emotion scores
subjmeans_score<-as.data.frame(ddply(data,c("Subject", "Sex", "Age", "Cohort", "Outcome"),summarise,
               	N		= length(!is.na(Score)),
               	Mean	= mean(Score, na.rm=TRUE),
               	sd  	= sd(Score, na.rm=TRUE),
               	se   	= sd / sqrt(N)))
subjmeans_score

subjmeans_score$Outcome<-ordered(subjmeans_score$Outcome, levels = c("Good", "Safe", "Bad"))
#Bar graph of score
quartz.options(width=5, height=5)
score_plot <- ggplot(subjmeans_score, aes(x=Cohort, y=Mean, fill = Outcome)) + 
  		geom_boxplot() +
	stat_summary(fun=mean, geom="point", position = position_dodge(width=0.75), shape=5, size=4) +  	
	theme_bw() +
    	expand_limits(y=c(0, 1)) +
    	theme(aspect.ratio=1/1) +
		ggtitle("(c) Emotional Responses") +
		scale_fill_manual(values=c("lightskyblue1", "dodgerblue", "dodgerblue4"), name="Outcome") +
    	scale_x_discrete(name="") +
    	scale_y_continuous(name="Mean Affect Score", breaks=seq(0,1.5,.25)) +
   		theme(  plot.title = element_text(hjust = 0.5),
	   	 		axis.text.x = element_text(face="bold", size = 12), 
   				axis.title.y = element_text(face="bold", size = 12), 
   				axis.text.y = element_text(size=12),
   				legend.text = element_text(size=10),
    			legend.title = element_text(face="bold", size=12)) 
score_plot

#subj means for swtiching
subjmeans_switching<-as.data.frame(ddply(test_data,c("Subject", "Cohort", "Outcome"),summarise,
               	N		= length(!is.na(Switch)),
               	mean	= mean(Switch, na.rm=TRUE),
               	sd  	= sd(Switch, na.rm=TRUE),
               	se   	= sd / sqrt(N)))
subjmeans_switching

subjmeans_switching$Outcome<-ordered(subjmeans_switching$Outcome, levels = c("Good", "Safe", "Bad"))
#Bar graph of switching
quartz.options(width=5, height=5)
switching_plot <- ggplot(data=subjmeans_switching, aes(x=Outcome, y=mean, fill = Outcome)) + 
  		geom_boxplot() +
  		stat_summary(fun=mean, geom="point", position = position_dodge(width=0.75), shape=5, size=4) +  	
    		theme_bw() +
    	expand_limits(y=c(0,1)) +
    	theme(aspect.ratio=1/1) +
		ggtitle("(d) Choice Switching") +
		scale_fill_manual(values=c("lightskyblue1", "dodgerblue", "dodgerblue4"), name="Outcome") +    		scale_x_discrete(name="") +
    	scale_x_discrete(name="") +
    	scale_y_continuous(name="Mean Choice Switching", breaks=seq(0,1,.25)) +
   		theme(  plot.title = element_text(hjust = 0.5),
	   	 		axis.text.x = element_text(face="bold", size = 12), 
   				axis.title.y = element_text(face="bold", size = 12), 
   				axis.text.y = element_text(size=12),
   				legend.text = element_text(size=10),
    			legend.title = element_text(face="bold", size=12))
switching_plot

#function to save legend
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

#save plot legend
legend1 <- get_legend(risk_plot)
legend2 <- get_legend(switching_plot)

#combine graphs
quartz.options(width=10, height=10)
toprow<- grid.arrange(Risk_age_effects, 
			risk_plot + theme(legend.position="none"), 
			legend1, ncol = 3, widths=c(5,5,1))
		
bottomrow<- grid.arrange(score_plot + theme(legend.position="none"), 
			switching_plot+ theme(legend.position="none"), 
			legend2, ncol = 3, widths=c(5,5,1))

Combined_risk_plots <- grid.arrange(toprow, bottomrow, nrow = 2)
Combined_risk_plots




