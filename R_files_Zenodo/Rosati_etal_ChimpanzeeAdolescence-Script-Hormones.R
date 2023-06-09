# Analyses for: Distinct developmental trajectories for risky and impulsive decision-making in chimpanzees
# Rosati, Emery Thompson, Atencia, & Buckholtz
# Hormones

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
library(tidyr) #data conversions
library(car) #type 3 anova for omnibus tests
library(bbmle) #AIC comparison
library(readxl) #to open data


########### LOAD DATA ###############################################################

# set screen print width
options(width = 120)

# set working directory

# Open hormone data
data <- read_excel("Rosati_etal_ChimpanzeeAdolescence-Data-Hormones.xlsx")

#rename sex factors
data$Sex <- as.character(data$Sex)
data$Sex[data$Sex=="F"]<-"Female"
data$Sex[data$Sex=="M"]<-"Male"

#age cohorts
data$Cohort <-"Younger"
data$Cohort[data$Age>=15]<-"Adult"
data$Cohort <- as.factor(data$Cohort)
data$Cohort <- relevel(data$Cohort, ref = "Younger")

#log transformed hormone values
data$logValue <- log(data$Value)

#specify factors in hormone data
framefacs <- c("Subject", "Sex", "Hormone","Unit","Cohort")
data[, framefacs] <- lapply(data[, framefacs], as.factor)

#check variables
lapply(list(data), str)

#split
cortisol_data <- data[ which(data$Hormone=='Cortisol'), ]
testosterone_data <- data[ which(data$Hormone=='Testosterone'), ]


############## HORMONE SUMMARY DATA ##################################################


#Creates summaries of logged data by subject with info for some analyses
subj_logdata_long_info <- ddply(data, c("Subject", "Hormone", "Age", "Sex"), summarise, 
				MeanValue = mean(logValue))
subj_logdata_long_info

#organize long version
hormone_subj_logdata_info<-spread(subj_logdata_long_info, Hormone, MeanValue)  


#Creates summaries of logged data by subject to combine with other data in intertask relationaship analyses
subj_logdata_long <- ddply(data, c("Subject", "Age", "Sex","Hormone"), summarise, 
				MeanValue = mean(logValue))
subj_logdata_long

hormone_subj_logdata<-spread(subj_logdata_long, Hormone, MeanValue)  

# Histogram to check distributions on means of log transformation on cortisol
ggplot(hormone_subj_logdata,aes(x=Cortisol))+ geom_histogram(position="dodge",binwidth=.1)

# Histogram to check distributions on means of log transformation on testosterone
ggplot(hormone_subj_logdata,aes(x=Testosterone))+ geom_histogram(position="dodge",binwidth=.1)

#save averages of log transformed data for cross-task analysis

# check age ranges
min(hormone_subj_logdata$Age) #6
max(hormone_subj_logdata$Age) #25

# sample size by sex
length(which(hormone_subj_logdata$Sex == 'Female')) # 19
length(which(hormone_subj_logdata$Sex == 'Male')) # 21


######################## ANALYSES OF CORTISOL ##############################

#Cortisol averages across all samples
length(cortisol_data$Value) # 169
mean(cortisol_data$Value, na.rm=TRUE) #4.04 
min(cortisol_data$Value, na.rm=TRUE) #.05
max(cortisol_data$Value, na.rm=TRUE) #32.48

#Cortisol averages across cohorts
cort_means_cohort <- ddply(cortisol_data, c("Cohort"), summarise, 
				TotalN = length(Value),
				Mean_Cortisol = mean(Value, na.rm=TRUE))
cort_means_cohort

# Histogram to check distributions
ggplot(cortisol_data,aes(x=Value))+ geom_histogram(position="dodge",binwidth=1)
# cortisol is right skewed: use gamma log link to account for skew

#Look at relationship between average of subject's log transformed cortisol and age
cor.test(hormone_subj_logdata_info$Cortisol, hormone_subj_logdata_info$Age,  method = "pearson")
#increases with age

#GLMMs using gamma link to deal with data skew
#Base model: account for sex
Cort_model1 = glmer(Value ~ (1|Subject) + Sex, 
		family=Gamma (link = "log"),
		data=cortisol_data)
summary(Cort_model1)

#Add in age
Cort_model2 = glmer(Value ~ (1|Subject) + Sex + Age, 
		family=Gamma (link = "log"),
		data=cortisol_data)
summary(Cort_model2)
anova(Cort_model1, Cort_model2)
#improves fit

#Add in interaction
Cort_model3 = glmer(Value ~ (1|Subject) + Age*Sex, 
		family=Gamma (link = "log"),
		data=cortisol_data)
summary(Cort_model3)
anova(Cort_model2, Cort_model3)
#no improvement

#Check for model 2: Add in age as cohort
Cort_model2b = glmer(Value ~ (1|Subject) + Sex + Cohort, 
		family=Gamma (link = "log"),
		data=cortisol_data)
summary(Cort_model2b)
anova(Cort_model1, Cort_model2b)
#improves fit as above

#Check for model 3: Add in interaction with Cohort
Cort_model3b = glmer(Value ~ (1|Subject) + Cohort*Sex, 
		family=Gamma (link = "log"),
		data=cortisol_data)
summary(Cort_model3b)
anova(Cort_model2b, Cort_model3b)
#no improvement as above

#omnibus test
Anova(Cort_model3, type = "III")
#effect of age

#AIC comparison
AICtab(Cort_model1, Cort_model2, Cort_model3, base = T, weights = T)
# model 2 is best by AIC 

#plot effects
quartz.options(width=10, height=6)
plot(allEffects(Cort_model3))

#look at residuals: normality and homogeneity 
plot(fitted(Cort_model3),residuals(Cort_model3), xlab='Fitted Values', ylab='Residuals'); abline(h=0)
qqnorm(resid(Cort_model3),main='Normality Check for Residuals')
qqline(resid(Cort_model3))
#ok

#save model 3 parameters
write.table(summary(Cort_model3)$coefficients, sep="\\t")



########################### ANALYSES OF TESTOSTERONE ##############################


#Testosterone averages across all samples
length(testosterone_data$Value) #160
mean(testosterone_data$Value, na.rm=TRUE) #4364.824
min(testosterone_data$Value, na.rm=TRUE) #423.2
max(testosterone_data$Value, na.rm=TRUE) #22420

#testosterone averages across cohorts
testosterone_means_cohort <- ddply(testosterone_data, c("Cohort"), summarise, 
				TotalN = length(Value),
				Mean_Cortisol = mean(Value, na.rm=TRUE))
testosterone_means_cohort


# Histogram to check distributions
ggplot(testosterone_data,aes(x=Value))+ geom_histogram(position="dodge",binwidth=500)
# right skewed: use gamma log link

#Relationship between mean of log transformed testosterone and age
cor.test(hormone_subj_logdata_info$Testosterone, hormone_subj_logdata_info$Age,  method = "pearson")
#testosterone increases with age

#GLMMs using gamma link to account for skew
#Base model: account for sex
Test_model1 = glmer(Value ~ (1|Subject) + Sex, 
		family=Gamma (link = "log"),
		data=testosterone_data)
summary(Test_model1)

#Add in age
Test_model2 = glmer(Value ~ (1|Subject) + Sex + Age, 
		family=Gamma (link = "log"),
		data=testosterone_data)
summary(Test_model2)
anova(Test_model1, Test_model2)
#improves fit

#Add in interaction
Test_model3 = glmer(Value ~ (1|Subject) + Age*Sex, 
		family=Gamma (link = "log"),
		data=testosterone_data,
		control = glmerControl(optimizer = "optimx", optCtrl=list(method = "bobyqa")))
summary(Test_model3)
anova(Test_model2, Test_model3)
#trend

#posthoc comparison of age trends
emt1<-emtrends(Test_model3,  "Sex", var = "Age") 
pairs(emt1)
#different: more positive slopes for males

#Check for model 2: Add in age as cohort
Test_model2b = glmer(Value ~ (1|Subject) + Sex + Cohort, 
		family=Gamma (link = "log"),
		data=testosterone_data)
summary(Test_model2b)
anova(Test_model1, Test_model2b)
#improves fit

#Check for model 3: Add in interaction as cohort
Test_model3b = glmer(Value ~ (1|Subject) + Cohort*Sex, 
		family=Gamma (link = "log"),
		data=testosterone_data)
summary(Test_model3b)
anova(Test_model2b, Test_model3b)
#improves

#posthocs: pairwise
emmeans (Test_model3b, pairwise ~ Sex*Cohort) #pairwise
#increase in males, not females

#omnibus test
Anova(Test_model3, type = "III")
#Age X sex interaction

#AIC comparison
AICtab(Test_model1, Test_model2, Test_model3, base = T, weights = T)
# model 3 is best by AIC 

#plot effects
quartz.options(width=10, height=6)
plot(allEffects(Test_model3))

#look at residuals: normality and homogeneity 
plot(fitted(Test_model3),residuals(Test_model3), xlab='Fitted Values', ylab='Residuals'); abline(h=0)
qqnorm(resid(Test_model3),main='Normality Check for Residuals')
qqline(resid(Test_model3))
#ok

#save model 3 parameters
write.table(summary(Test_model3)$coefficients, sep="\\t")


######### HORMONE PLOTS ###############################################

#Creates averges of UNlogged data by subject
cortisol_subj_means <- ddply(cortisol_data, c("Subject", "Sex", "Cohort"), summarise, 
               	N_subj		= length(!is.na(Value)),
               	Mean_subj		= mean(Value, na.rm=TRUE),
               	SD_subj  		= sd(Value, na.rm=TRUE),
               	SE_subj   	= SD_subj/ sqrt(N_subj))
cortisol_subj_means

#Creates averges of UNlogged data
cortisol_all_means <- ddply(cortisol_subj_means, c("Cohort", "Sex"), summarise, 
               	N		= length(!is.na(Mean_subj)),
               	Nav		= mean(N_subj, na.rm=TRUE),
               	Mean		= mean(Mean_subj, na.rm=TRUE),
               	SD  		= sd(Mean_subj, na.rm=TRUE),
               	SE   	= SD/ sqrt(N))
cortisol_all_means


#Bar graph of cort 
quartz.options(width=5, height=5)
cort_plot <- ggplot(data=cortisol_all_means, aes(x=Cohort, y=Mean, group = Sex, colour = Sex)) + 
		geom_line(size = 1, position=position_dodge(.2)) + 
		geom_point(size = 3, position=position_dodge(.2)) +
    	geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=.1,position=position_dodge(.2)) +
    	theme_bw() +
		theme(aspect.ratio=1/1) +
		 ggtitle("(a) Cortisol") +
	scale_color_manual(values=c("indianred3", "deepskyblue4")) +
    	scale_x_discrete(name="") +
    	scale_y_continuous(name="Mean Cortisol (ng/ml)", breaks=seq(0,8, 1)) +
   		theme(  plot.title = element_text(hjust = 0.5),
   				legend.text = element_text(size = 14),
   				legend.title= element_blank(),
	   	 		axis.text.x = element_text(face="bold", size = 15), 
   				axis.title.y = element_text(face="bold", size = 10), 
   				axis.text.y = element_text(size=12))
cort_plot

#Creates averges of UNlogged data by subject
testosterone_subj_means <- ddply(testosterone_data, c("Subject", "Sex", "Cohort"), summarise, 
               	N_subj		= length(!is.na(Value)),
               	Mean_subj		= mean(Value, na.rm=TRUE),
               	SD_subj  		= sd(Value, na.rm=TRUE),
               	SE_subj   	= SD_subj/ sqrt(N_subj))
testosterone_subj_means

#Creates averges of UNlogged data 
testosterone_all_means <- ddply(testosterone_subj_means, c("Sex", "Cohort"), summarise, 
               	N		= length(!is.na(Mean_subj)),
               	Mean		= mean(Mean_subj, na.rm=TRUE),
               	SD  		= sd(Mean_subj, na.rm=TRUE),
               	SE   	= SD/ sqrt(N))
testosterone_all_means


#Bar graph of testosterone
quartz.options(width=5, height=5)
testosterone_plot <- ggplot(data=testosterone_all_means, aes(x=Cohort, y=Mean, group = Sex, colour = Sex)) + 
		geom_line(size = 1, position=position_dodge(.2)) + 
		geom_point(size = 3, position=position_dodge(.2)) +
    	geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=.1,position=position_dodge(.2)) +
    	theme_bw() +
		theme(aspect.ratio=1/1) +
		 ggtitle("(b) Testosterone") +
	scale_color_manual(values=c("indianred3", "deepskyblue4")) +
    	scale_x_discrete(name="") +
    	scale_y_continuous(name="Mean Testosterone (pg/ml)", breaks=seq(0,10000, 1000)) +
   		theme(  plot.title = element_text(hjust = 0.5),
   				legend.text = element_text(size = 14),
   				legend.title= element_blank(),
	   	 		axis.text.x = element_text(face="bold", size = 15), 
   				axis.title.y = element_text(face="bold", size = 10), 
   				axis.text.y = element_text(size=12))
testosterone_plot

#function to save legend
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

#save plot legend
legend <- get_legend(testosterone_plot)

quartz.options(width=10, height=5)
Hormone_plots <- grid.arrange(
	cort_plot + theme(legend.position="none"), 
	testosterone_plot + theme(legend.position="none"), 
	legend, ncol = 3, widths=c(2.1, 2.2, 0.5))
Hormone_plots


########################## SUPPLEMENTAL PLOTS #################################

#use wide format data: LOGGED values
head(hormone_subj_logdata_info) 

#Scatter plot: Cortisol difference score
Cortisol_scatter <- ggplot(hormone_subj_logdata_info, aes(x=Age, y=Cortisol, color=Sex)) +
	geom_point() +
    geom_smooth(method=lm, se=FALSE)+
    theme_bw()+
	ggtitle("(a) Cortisol") +
	scale_y_continuous(name="Log Transformed Cortisol (ng/ml)") +
	   scale_colour_manual(values=c("indianred3", "deepskyblue4")) +
	    	    	scale_x_continuous(name="Age (Years)")  +
   		theme(  plot.title = element_text(hjust = 0.5),
   				legend.text = element_text(size = 14),
	   	 		axis.title.x = element_text(face="bold", size = 12), 
	   	 		axis.text.x = element_text(size = 10), 
   				axis.title.y = element_text(face="bold", size = 12), 
   				axis.text.y = element_text(size=10))
Cortisol_scatter


#Scatter plot: Cortisol difference score
Testosterone_scatter <- ggplot(hormone_subj_logdata_info, aes(x=Age, y=Testosterone, color=Sex)) +
	geom_point() +
    geom_smooth(method=lm, se=FALSE)+
    theme_bw()+
	ggtitle("(b) Testosterone") +
	scale_y_continuous(name="Log Transformed Testosterone (pg/ml)") +
	   scale_colour_manual(values=c("indianred3", "deepskyblue4")) +
	    	    	scale_x_continuous(name="Age (Years)")  +
   		theme(  plot.title = element_text(hjust = 0.5),
   				legend.text = element_text(size = 14),
	   	 		axis.title.x = element_text(face="bold", size = 12), 
	   	 		axis.text.x = element_text(size = 10), 
   				axis.title.y = element_text(face="bold", size = 12), 
   				axis.text.y = element_text(size=10))
Testosterone_scatter


quartz.options(width=10, height=5)
Hormone_scatter_plots <- grid.arrange(
	Cortisol_scatter + theme(legend.position="none"), 
	Testosterone_scatter, 
	ncol = 2, widths=c(2.3, 3))
Hormone_scatter_plots
