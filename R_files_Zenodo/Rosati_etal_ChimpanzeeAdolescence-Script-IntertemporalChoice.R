# Analyses for: Distinct developmental trajectories for risky and impulsive decision-making in chimpanzees
# Rosati, Emery Thompson, Atencia, & Buckholtz
# Intertemporal Choice Task


######################## LOAD PACKAGES ######################

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
library(tidyr) #long format data
library(arm) #binned residual plot
library(lsr) #cohensD


######################## LOAD DATA ######################

# set screen print width
options(width = 120)

## set working directory

# Open discounting task data
data <- read_excel("Rosati_etal_ChimpanzeeAdolescence-Data-IntertemporalChoice.xlsx")

#rename sex factors
data$Sex <- as.character(data$Sex)
data$Sex[data$Sex=="F"]<-"Female"
data$Sex[data$Sex=="M"]<-"Male"

#set age cohorts
data$Cohort <-"Younger"
data$Cohort[data$Age>=15]<-"Adult"
data$Cohort <- as.factor(data$Cohort)
data$Cohort <- relevel(data$Cohort, ref = "Younger")

#rename trial type factors
data$Type <- as.character(data$Type)
data$Type[data$Type=="E"]<-"Exposure"
data$Type[data$Type=="C"]<-"Test"
data$Type[data$Type=="Num"]<-"NumberPretest"

#rename choice factors to use for some graphs
data$Choice <- as.character(data$Larger_reward)
data$Choice[data$Choice=="0"]<-"Immediate"
data$Choice[data$Choice=="1"]<-"Delayed"

#specify factors on task data
framefacs <- c("Subject","Sex", "Session","Type", "Cohort", "Choice")
data[, framefacs] <- lapply(data[, framefacs], as.factor)

#make sure affect score is coded as a number
data$Score <- as.numeric(data$Score)

#check variables
lapply(list(data), str)

#split out different subsets of data for analysis
test_data <- data[ which(data$Type=='Test'), ] #test trial data from test session
choice_data <- data[ which(data$Type=='Test' | data$Type=='NumberPretest'), ] #all choice trials both both sessions
affect_data <- data[ which(data$Type=='Test' | data$Type=='Exposure'), ] #all trials from test session to analyze affect scores



######### SUBJECT SUMMARY DATA AND CHECK ###############################################

#create summary subjectdata with subject means for discounting choices; affect score in response to larger, delayed reward; affect score in response to smaller, immediate reward
subj_data <- ddply(data[ which(data$Type=='Test'), ], c("Subject", "Sex", "Age", "Cohort"), summarise,  Discounting = mean(Larger_reward))
score3_mean <- ddply(data[ which(data$Choice=='Delayed'), ], c("Subject"), summarise, Mean3Score = mean(Score, na.rm=TRUE))
score1_mean <- ddply(data[ which(data$Choice=='Immediate'), ], c("Subject"), summarise, Mean1Score = mean(Score, na.rm=TRUE))
num_mean <- ddply(data[],c("Subject"), summarise, NumberPretest = mean(Num_Pretest, na.rm=TRUE))

#merge into a new dataframe
subj_data<-merge(subj_data, score3_mean, all = TRUE, by = c('Subject'))
subj_data<-merge(subj_data, score1_mean, all = TRUE, by = c('Subject'))
subj_data<-merge(subj_data, num_mean, all = TRUE, by = c('Subject'))

#check variables
lapply(list(subj_data), str)
#save subject data for cross-task analysis

# check age ranges
min(subj_data$Age) #6
max(subj_data$Age) #25

# sample size by sex
length(which(subj_data$Sex == 'Female')) # 19
length(which(subj_data$Sex == 'Male')) # 21

#split out subject data by cohort for some comparisons
subj_data_young <- subj_data[ which(subj_data$Cohort=='Younger'), ]
subj_data_adult <- subj_data[ which(subj_data$Cohort=='Adult'), ]

#make long format data as well
keycol <- "Condition"
valuecol <- "Mean_large"
gathercols <- c("Discounting","NumberPretest")
subj_data_long <- gather_(subj_data, keycol, valuecol, gathercols)
#make sure condition is a factor
subj_data_long$Condition <- as.factor(subj_data_long$Condition)
lapply(list(subj_data_long), str)



######### COMPARISON OF DISCOUNTING VERSUS NUMBER PRETEST SESSION CHOICES ######################################

#overall means for discounting versus number sessions
means_discount<-as.data.frame(ddply(subj_data_long,c("Condition"),summarise,
               	N		= length(!is.na(Mean_large)),
               	Mean	= mean(Mean_large, na.rm=TRUE),
               	sd  	= sd(Mean_large, na.rm=TRUE),
               	se   	= sd / sqrt(N)))
means_discount

#paired samples t test for discounting test trials versus number trials (all chimpanzees)
t.test(subj_data$Discounting,subj_data$Number, paired=TRUE)
#more likely to choose large in number pretest overall
cohensD(subj_data$Discounting,subj_data$Number, method = 'paired') 
# medium effect size  

#paired samples t test for discounting test trials versus number trials: younger only
t.test(subj_data_young$Discounting,subj_data_young$Number, paired=TRUE)
#more likely to choose large in number pretest in younger
cohensD(subj_data_young$Discounting,subj_data_young$Number, method = 'paired') 
# medium effect size  

#paired samples t test for discounting test trials versus number trials: adults only
t.test(subj_data_adult$Discounting,subj_data_adult$Number, paired=TRUE)
#more likely to choose large in number pretest in adults
cohensD(subj_data_adult$Discounting,subj_data_adult$Number, method = 'paired') 
# medium effect size  

#t-test: comparison of number preference by age
t.test(subj_data$NumberPretest ~ subj_data$Cohort, var.equal = TRUE)
#no difference   
#calculate effect size
cohensD(subj_data$NumberPretest ~ subj_data$Cohort) 
# medium effect size  



######### INTERTEMPORAL CHOICE ANALYSIS ###############################################

#Histogram to check distributions
ggplot(subj_data,aes(x=Discounting))+ geom_histogram(position="dodge",binwidth=.1)

#overall means by cohort
#subj means for larger reward
discount_overall_mean<-as.data.frame(ddply(subj_data,c("Cohort"),summarise,
               	N		= length(!is.na(Discounting)),
               	mean		= mean(Discounting, na.rm=TRUE),
               	sd  		= sd(Discounting, na.rm=TRUE),
               	se   	= sd / sqrt(N)))
discount_overall_mean

#t-test: comparison of adults and younger
t.test(subj_data$Discounting ~ subj_data$Cohort, var.equal = TRUE)
# no difference  cohensD(subj_data$Discounting ~ subj_data$Cohort) 
cohensD(subj_data$Discounting ~ subj_data$Cohort) 
# small effect size     


#GLMMS of inter-temporal choices
#Base model: account for sex, number pretest performance, and trial type within choice trials (1-14)
Discount_model1 = glmer(Larger_reward ~ (1|Subject) + Sex + Num_Pretest + Trial_type, 
		family=binomial(link = "logit"),
		data=test_data)
summary(Discount_model1)

#Add in age (continuous)
Discount_model2 = glmer(Larger_reward ~ (1|Subject) + Sex + Num_Pretest + Trial_type + Age, 
		family=binomial(link = "logit"),
		data=test_data)
summary(Discount_model2)
anova(Discount_model1, Discount_model2)
#no improvement

#Add in age cohort
Discount_model2b = glmer(Larger_reward ~ (1|Subject) + Sex + Num_Pretest + Trial_type + Cohort, 
		family=binomial(link = "logit"),
		data=test_data)
summary(Discount_model2b)
anova(Discount_model1, Discount_model2b)
# also no improvement

#omnibus test
Anova(Discount_model2, type = "III")

#AIC comparison
AICtab(Discount_model1, Discount_model2, base = T, weights = T)
# model 1 is best by AIC 

#plot effects
quartz.options(width=10, height=6)
plot(allEffects(Discount_model2))

#look at binned residual plot for logistic regression to assess assumptions
binnedplot(fitted(Discount_model2),residuals(Discount_model2, type="response"), cex.pts=1, col.int="black")
#looks ok

#save model 2 parameters
write.table(summary(Discount_model2)$coefficients, sep="\\t")


######### AFFECT SCORE ANALYSIS ###############################################


#Histogram to check distributions of mean affect scores
ggplot(subj_data,aes(x=Mean3Score))+ geom_histogram(position="dodge",binwidth=.1)
ggplot(subj_data,aes(x=Mean1Score))+ geom_histogram(position="dodge",binwidth=.1)
#ok

#subj means for emotion scores
subjmeans_score<-as.data.frame(ddply(affect_data,c("Subject","Cohort", "Choice"),summarise,
               	N		= length(!is.na(Score)),
               	mean		= mean(Score, na.rm=TRUE),
               	sd  		= sd(Score, na.rm=TRUE),
               	se   	= sd / sqrt(N)))
subjmeans_score
#check this matches summary above

#overall means for scores by cohort
means_score<-as.data.frame(ddply(subjmeans_score,c("Cohort","Choice"),summarise,
               	N		= length(!is.na(mean)),
               	Mean		= mean(mean, na.rm=TRUE),
               	SD  		= sd(mean, na.rm=TRUE),
               	SE   	= SD / sqrt(N)))
means_score

#overall means for scores across all individuals
means_score<-as.data.frame(ddply(subjmeans_score,c("Choice"),summarise,
               	N		= length(!is.na(mean)),
               	Mean		= mean(mean, na.rm=TRUE),
               	SD  		= sd(mean, na.rm=TRUE),
               	SE   	= SD / sqrt(N)))
means_score

#relevel so immediate choice is the baseline
affect_data$Choice<-relevel(affect_data$Choice, ref = "Immediate")

#LMMs to examine affective responses; includes exposure and test trials in test sessions
#Base model: account for trial type (test versus exposure), trial number within session (1-22) and sex
Affect_model1 = lmer(Score ~ (1|Subject) + Type + Trial + Sex, 
		data=affect_data)
summary(Affect_model1)
#effect of type (test > exposure) and trial

#Add in choice (immediate versus delayed)
Affect_model2 = lmer(Score ~ (1|Subject) + Type + Trial + Sex + Choice, 
		data=affect_data)
summary(Affect_model2)
anova(Affect_model1, Affect_model2)
#improves fit: lower scores in response to immediate rewards

#Add in age (continuous)
Affect_model3 = lmer(Score ~ (1|Subject) + Type + Trial + Sex + Choice + Age, 
		data=affect_data)
summary(Affect_model3)
anova(Affect_model2, Affect_model3)
#improves fit: lower scores with increasing age 

#Add in age (continuous) X choice
Affect_model4 = lmer(Score ~ (1|Subject) + Type + Trial + Sex + Choice*Age, 
		data=affect_data)
summary(Affect_model4)
anova(Affect_model3, Affect_model4)
#improves fit

#posthoc comparison of age trends
emt1<-emtrends(Affect_model4,  "Choice", var = "Age") 
pairs(emt1)
#different: more negative slopes for delayed

#Check for model 3 using age (cohort)
Affect_model3b = lmer(Score ~ (1|Subject) + Type + Trial + Sex + Choice + Cohort, 
		data=affect_data)
summary(Affect_model3b)
anova(Affect_model2, Affect_model3b)
#trend

#Check for model 4 using cohort X choice 
#relevel so delayed choice is the baseline to interpret interaction comparison
affect_data$Choice<-relevel(affect_data$Choice, ref = "Delayed")
Affect_model4b = lmer(Score ~ (1|Subject) + Type + Trial + Sex + Choice*Cohort, 
		data=affect_data)
summary(Affect_model4b)
anova(Affect_model3b, Affect_model4b)
anova(Affect_model2, Affect_model4b)
#improves fit, also in comparison to model 2

#examine interaction
emm1<-emmeans (Affect_model4b, ~ Choice | Cohort) #choice effect by age cohort
con1 <- contrast(emm1, interaction = "pairwise")
pairs(con1, by = NULL)
#choice effect is bigger in younger age cohort

#omnibus test
Anova(Affect_model4, type = "III")
#significant interaction between choice outcome and age

#AIC comparison
AICtab(Affect_model1, Affect_model2, Affect_model3, Affect_model4, base = T, weights = T)
# model 4 is best by AIC also
#AIC comparison with cohort
AICtab(Affect_model1, Affect_model2, Affect_model3b, Affect_model4b, base = T, weights = T)
# model 4b is best by AIC also when using cohort models

#plot effects
quartz.options(width=10, height=6)
plot(allEffects(Affect_model4))

#look at residuals
qqnorm(resid(Affect_model4),main='Normality Check for Residuals')
qqline(resid(Affect_model4))
#looks ok

#save model 4 parameters
write.table(summary(Affect_model4)$coefficients,  sep="\\t")


#################### INTERTEMPORAL CHOICE PLOTS ###############################################

#overall means for discounting
means_discount<-as.data.frame(ddply(subj_data_long,c("Subject","Condition"),summarise,
               	N		= length(!is.na(Mean_large)),
               	Mean	= mean(Mean_large, na.rm=TRUE),
               	sd  	= sd(Mean_large, na.rm=TRUE),
               	se   	= sd / sqrt(N)))
means_discount

levels (means_discount$Condition) <-c("Delay", "Number")
means_discount$Condition <- relevel(means_discount$Condition, ref = "Number")

#Bar graph of discounting preferences
quartz.options(width=5, height=5)
discount_plot <- ggplot(data=means_discount, aes(x=Condition, y=Mean, fill =  Condition)) + 
  		geom_boxplot() +
  		stat_summary(fun=mean, geom="point", position = position_dodge(width=0.75), shape=5, size=4) + 
  		theme_bw() +
		theme(aspect.ratio=1/1) +
		expand_limits(y=c(0,1))+
		#theme(legend.position="none") +
		ggtitle("(a) Number vs Delay Test") +
		scale_fill_manual(values=c("lightskyblue1", "dodgerblue4")) +
    		scale_x_discrete(name  ="") +
    		scale_y_continuous(name="Mean Choices for Larger Reward", breaks=seq(0,1,.2)) +
   		theme(  plot.title = element_text(hjust = 0.5),
	   	 		axis.text.x = element_text(face="bold", size = 12), 
   				axis.title.y = element_text(face="bold", size = 12), 
   				axis.text.y = element_text(size=12),
   				legend.text = element_text(size=10),
    			legend.title = element_text(face="bold", size=12))
discount_plot


#Convert parameter estimates for discount choice preferences
Discount_age_effects <- allEffects(Discount_model2, xlevels = list(Age = seq(6, 25, by = 1)))
Discount_effects_df <- as.data.frame(Discount_age_effects[["Age"]])

Discount_age_effects <-ggplot(Discount_effects_df, aes(x = Age, y = fit)) +
    geom_line(color = "royalblue4") +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, fill = "dodgerblue4") +
    geom_point(data = subj_data, aes(x = Age, y = Discounting), position = "jitter", inherit.aes = FALSE, color = "grey50") +
    theme_bw() +
    ggtitle("(b) Temporal Preferences by Age")+
    scale_x_continuous(name="Age (Years)", breaks=seq(6, 25, 2)) +
    expand_limits(y=c(0,1)) +
    theme(aspect.ratio=1/1) +
    scale_y_continuous(name="Predicted value (Choose Delayed)") +
    theme(plot.title = element_text(hjust = 0.5),
    	axis.title.x = element_text(face="bold", size=12), 
    	axis.title.y = element_text(face="bold", size = 12), 				
    	axis.text.y = element_text(size=10), 
    	axis.text.x = element_text(size=10), 
    	legend.text = element_text(size=10),
    	legend.title = element_text(face="bold", size=12))
Discount_age_effects


#Convert parameter estimates
Discount_affect_effects <- allEffects(Affect_model4, xlevels = list(Age = seq(6, 25, by = 1)))
Discount_affect_effects_df <- as.data.frame(Discount_affect_effects[["Choice:Age"]])

Discount_affect_effects <-ggplot(Discount_affect_effects_df, aes(x = Age, y = fit, color = Choice, fill = Choice)) +
    geom_line(values=c()) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
    #geom_point(data = subj_data, aes(x = Age, y = Mean1Score), position = "jitter", inherit.aes = FALSE, color = "lightskyblue3") +
    #geom_point(data = subj_data, aes(x = Age, y = Mean3Score), position = "jitter", inherit.aes = FALSE, color = "dodgerblue4") +
    theme_bw() +
    ggtitle("(c) Emotional Responses")+
	scale_fill_manual(values=c("dodgerblue4", "lightskyblue3")) +
    scale_colour_manual(values=c("dodgerblue4", "lightskyblue3")) +
    scale_x_continuous(name="Age (Years)", breaks=seq(6, 25, 2)) +
    expand_limits(y=c(0,2)) +
    theme(aspect.ratio=1/1) +
    scale_y_continuous(name="Predicted value (Affect Score)") +
    theme(plot.title = element_text(hjust = 0.5),
    	axis.title.x = element_text(face="bold", size=12), 
    	axis.title.y = element_text(face="bold", size = 12), 				
    	axis.text.y = element_text(size=10), 
    	axis.text.x = element_text(size=10), 
    	legend.text = element_text(size=10),
    	legend.title = element_text(face="bold", size=12))
Discount_affect_effects

#save plot legend
legend1 <- get_legend(Discount_affect_effects)

#combine graphs
quartz.options(width=10, height=5)
combined_discounting_plots <- grid.arrange(discount_plot + theme(legend.position="none"), 
	Discount_age_effects + theme(legend.position="none"), 
	Discount_affect_effects + theme(legend.position="none"), 
	legend1, 
	ncol = 4, widths=c(4, 4, 4, 1))



############### SUPPLEMENTARY FIGURE: AFFECT ######################################

#subj means for emotion scores
subjmeans_score<-as.data.frame(ddply(affect_data,c("Subject","Cohort", "Choice"),summarise,
               	N		= length(!is.na(Score)),
               	Mean	= mean(Score, na.rm=TRUE),
               	SD  	= sd(Score, na.rm=TRUE),
               	SE   	= SD / sqrt(N)))
subjmeans_score

subjmeans_score$Choice<-ordered(subjmeans_score$Choice, levels = c("Immediate", "Delayed"))

#Bar graph of score
quartz.options(width=5, height=5)
score_boxplot <- ggplot(subjmeans_score, aes(x=Cohort, y=Mean, fill = Choice)) + 
  	geom_boxplot() +
	stat_summary(fun=mean, geom="point", position = position_dodge(width=0.75), shape=5, size=4) +  	
	theme_bw() +
    expand_limits(y=c(0, 1)) +
    theme(aspect.ratio=1/1) +
	scale_fill_manual(values=c("lightskyblue1", "dodgerblue4"), name="Choice") +
    scale_x_discrete(name="") +
    scale_y_continuous(name="Mean Affect Score", breaks=seq(0,2,.5)) +
   	theme(  plot.title = element_text(hjust = 0.5),
	   	 	axis.text.x = element_text(face="bold", size = 12), 
   			axis.title.y = element_text(face="bold", size = 12), 
   			axis.text.y = element_text(size=12),
   			legend.text = element_text(size=10),
    		legend.title = element_text(face="bold", size=12)) 
score_boxplot


