# Analyses for: Distinct developmental trajectories for risky and impulsive decision-making in chimpanzees
# Rosati, Emery Thompson, Atencia, & Buckholtz
# Inter-task relationships


########### LOAD PACKAGES #############################################

library(lattice) # plots
library(ggplot2) #plots
library(effects) #plots of effects
library(gridExtra) #combining plots
library(lme4) #install GLMM package
library(emmeans) #multiple comparisons
library(optimx) #for mixed model optimization
library(plyr) # split data dpply function
library(tidyr) #data conversions
library(car) #type 3 anova for omnibus tests
library(bbmle) #AIC comparison
library(readxl) #to open data
library(lmerTest) #tests for linear mixed effects models
library(cowplot) #plot grid to integrate graphs
library(dplyr) #organize data
library(corrplot) #correlation plots
library(Hmisc) #graphing
library(lmtest) #comparing linear models
library(psych) #does bartlett's test on PCA
library(paran) #to perform parallel analysis for factor selection
library(gvlma) #check model assumptions
library(factoextra) #for PCA
library(ggcorrplot) #ggplot for correlations



########### LOAD TASK DATA #####################################################

# set screen print width
options(width = 120)

## set working directory

# Open combined task data
subj_data <- read_excel("Rosati_etal_ChimpanzeeAdolescence-Data-TaskRelations.xlsx")

# assign age cohorts
subj_data$Cohort <-"Younger"
subj_data$Cohort[subj_data$Age>=15]<-"Adult"
subj_data$Cohort <- as.factor(subj_data$Cohort)
subj_data$Cohort <- relevel(subj_data$Cohort, ref = "Younger")

#rename sex factors
subj_data$Sex <- as.character(subj_data$Sex)
subj_data$Sex[subj_data$Sex=="F"]<-"Female"
subj_data$Sex[subj_data$Sex=="M"]<-"Male"

#specify factors on task data
framefacs <- c("Subject","Sex", "Cohort")
subj_data[, framefacs] <- lapply(subj_data[, framefacs], as.factor)

#check variables
lapply(list(subj_data), str)


# check age ranges
min(subj_data$Age) #6
max(subj_data$Age) #25

# sample size by sex
length(which(subj_data$Sex == 'Female')) # 19
length(which(subj_data$Sex == 'Male')) # 21

#subset cohort data
adult_data<-subj_data[ which(subj_data$Cohort=='Adult'), ]
young_data<-subj_data[ which(subj_data$Cohort=='Younger'), ]




################### DATA FOR CORRELATONS #################################

#select columns that are numeric for all individuals
numeric_data<-as.data.frame(dplyr::select_if(subj_data, is.numeric))

#select variables for correlations
subset_data <- subset(numeric_data, select = c(Age, Risk, RiskAffectScore, Discounting, DiscountingAffectScore, Cortisol, Testosterone))

corr_data = rcorr(as.matrix(subset_data))
corr_data
corr_data.coeff = corr_data$r
#p values for plot
corr_data.p = corr_data$P

#save correlation coefficients
#save correlation p values



########################################################################
######### RELATIONSHIPS BETWEEN DISCOUNTING & RISK #####################
########################################################################

######## Look at pairwise correlations

#Relationship between discounting and risk
cor.test(subj_data$Discounting, subj_data$Risk,  method = "pearson")
#no sign correlation
#aligns with trial-by-trial results

#Risk and Age
cor.test(subj_data$Age, subj_data$Risk,  method = "pearson")
#support for correlation
#aligns with trial-by-trial results

#Discounting and Age
cor.test(subj_data$Age, subj_data$Discounting,  method = "pearson")
#no sign correlation


######## Does Age moderate relationship between discounting and risk?

#task comparison: predict discounting using risk, sex, and cohort
comp_model1 <- lm(Discounting ~ Risk + Cohort + Sex, data = subj_data)
summary(comp_model1)

#add in interaction between risk and cohort to test if age moderates effects
comp_model2 <- lm(Discounting ~ Risk*Cohort + Sex, data = subj_data)
summary(comp_model2)
lrtest(comp_model1, comp_model2)
#improves fit

#check assumptions of linear models
gvlma(comp_model2)
#ok

#posthoc comparison of trends by cohort
emt1<-emtrends(comp_model2,  "Cohort", var = "Risk") 
pairs(emt1)
#different: more positive slopes for adults


plot(allEffects(comp_model2))
# Moderates effect: YES
# Greater risk-seeking increases patience in adults
# downward trend in juveniles


###############################################################################
############ RELATIONSHIPS BETWEEN OTHER MEASURES and AGE #####################
###############################################################################

#Risk affect score and Age
cor.test(subj_data$Age, subj_data$RiskAffectScore,  method = "pearson")
#no correlation
#aligns with trial-by-trial results

#Discounting affect score and Age
cor.test(subj_data$Age, subj_data$DiscountingAffectScore,  method = "pearson")
# correlation: affect difference score decreases with age
#aligns with trial-by-trial results

#Cortisol and Age
cor.test(subj_data$Age, subj_data$Cortisol,  method = "pearson")
#correlation: cortisol increases with age

#Testosterone and Age
cor.test(subj_data$Age, subj_data$Testosterone,  method = "pearson")
#correlation: av. testosterone increases with age



###############################################################################
######### RELATIONSHIPS BETWEEN OTHER BEHAVIORAL MEASURES #####################
###############################################################################


#Relationship between Risk and risk affect score
cor.test(subj_data$RiskAffectScore, subj_data$Risk,  method = "pearson")
#no sign correlation

#Relationship between Risk and discounting affect score
cor.test(subj_data$DiscountingAffectScore, subj_data$Risk,  method = "pearson")
# sign correlation: higher discounting affect score associated with more risk-taking

#Relationship between Discounting and discounting affect score
cor.test(subj_data$DiscountingAffectScore, subj_data$Discounting,  method = "pearson")
#no sign correlation

#Relationship between Discounting and risk affect score
cor.test(subj_data$RiskAffectScore, subj_data$Discounting,  method = "pearson")
#no sign correlation

#Relationship between risk affect score and discounting affect score
cor.test(subj_data$DiscountingAffectScore, subj_data$RiskAffectScore,  method = "pearson")
#no sign correlation



##################################################################################
######### RELATIONSHIPS BEHAVIORAL MEASURES AND TESTOSTERONE #####################
##################################################################################


#Relationship between Risk and T
cor.test(subj_data$Testosterone, subj_data$Risk,  method = "pearson")
#no sign correlation

#Relationship between Discount and T
cor.test(subj_data$Testosterone, subj_data$Discounting,  method = "pearson")
#no sign correlation

#Relationship between discounting affect and T
cor.test(subj_data$Testosterone, subj_data$DiscountingAffectScore,  method = "pearson")
#no sign correlation

#Relationship between Risk affect and T
cor.test(subj_data$Testosterone, subj_data$RiskAffectScore,  method = "pearson")
#sign correlation: risk affect score increases with higher testosterone

#check if age is also included: is it due to general age co-relationships?
#include age in a base model
model_T_riskaffect1 <- lm (RiskAffectScore ~ Age, data = subj_data)
summary(model_T_riskaffect1)

#add in testosterone in a second model
model_T_riskaffect2 <- lm (RiskAffectScore ~ Age + Testosterone, data = subj_data)
summary(model_T_riskaffect2)
anova(model_T_riskaffect1, model_T_riskaffect2, test = 'LRT')
#including testosterone improves fit even when accounting for age


##################################################################################
############# RELATIONSHIPS BEHAVIORAL MEASURES AND CORTISOL #####################
##################################################################################


#Relationship between Risk and Cortisol
cor.test(subj_data$Cortisol, subj_data$Risk,  method = "pearson")
#trend

#check if age is also included: is it due to age?
#include age in base model
model_cort_risk1 <- lm (Risk ~ Age, data = subj_data)
summary(model_cort_risk1)

#add cortisol in a second model
model_cort_risk2 <- lm (Risk ~ Age + Cortisol, data = subj_data)
summary(model_cort_risk2)
anova(model_cort_risk1, model_cort_risk2, test = 'LRT')
#age accounts for variance, adding cortisol does not improve fit

#Relationship between Discount and Cortisol
cor.test(subj_data$Cortisol, subj_data$Discounting,  method = "pearson")
#no sign correlation

#Relationship between Risk affect and Cortisol
cor.test(subj_data$Cortisol, subj_data$RiskAffectScore,  method = "pearson")
#no sign correlation

#Relationship between discounting affect and Cortisol
cor.test(subj_data$Cortisol, subj_data$DiscountingAffectScore,  method = "pearson")
#sign correlation

#check if age is also included: is it due to age?
#Base model accounting for age
model_cort_discountaffect1 <- lm (DiscountingAffectScore ~ Age, data = subj_data)
summary(model_cort_discountaffect1)

#include cortisol in a second model
model_cort_discountaffect2 <- lm (DiscountingAffectScore ~ Age + Cortisol, data = subj_data)
summary(model_cort_discountaffect2)
anova(model_cort_discountaffect1, model_cort_discountaffect2, test = 'LRT')
#age accounts for variance, adding cortisol does not improve fit

#Relationship between Testosterone and Cortisol
cor.test(subj_data$Cortisol, subj_data$Testosterone,  method = "pearson")
# correlation

#check if age is also included: is it due to age?
#include age in a base model
model_cort_T1 <- lm (Cortisol ~ Age, data = subj_data)
summary(model_cort_T1)

#add in testosterone in a second model
model_cort_T2 <- lm (Cortisol ~ Age + Testosterone, data = subj_data)
summary(model_cort_T2)
anova(model_cort_T1, model_cort_T2, test = 'LRT')
#improves fit even accounting for age



#################################################################################
######################## PRINCIPAL COMPONENTS ANALYSIS ##########################
#################################################################################


PCA_data<-subset(subj_data, select = c(Risk, RiskAffectScore, Discounting, DiscountingAffectScore, Cortisol, Testosterone))

#Calculate the correlation matrix
cor_matrix <- cor(PCA_data)

#Assessing the quality of the corralation matrix (Sampling Adequacy)
#Bartlettâ€™s test
cortest.bartlett(cor_matrix, n = 40) 
#significant, can proceed

#conduct PCA
task.pca<-prcomp(na.omit(PCA_data), center = TRUE, scale = TRUE, retx = TRUE)
summary(task.pca)

paran(PCA_data, cfa = FALSE, graph = TRUE, color = TRUE, col = c("black", "red", "blue"))
#Keep component with eingenvalues > 1
#only first 2 components to retain

#select dimensions
task.pca$sdev ^ 2
#keep only the eigenvalues > 1
#retain PC1, PC2

fviz_eig(task.pca) #scree plot
#percent variance explained by PCA
round(task.pca$sdev^2/sum(task.pca$sdev^2)*100) # percent explained variance


###### graph of tasks: dimension 1 and 2 ######
fviz_pca_var(task.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

########### Plot dimensions #####

#loadings
task.pca$rotation
#save loadings

#PC1
load <- task.pca$rotation
sorted.loadings=load[order(load [, 1]), 1]
Main="Loadings plot for PC1"
xlabs= "Variables loading"
dotplot(sorted.loadings, main=Main, xlab=xlabs, cex=1.5, col="red")

#PC2
load <- task.pca$rotation
sorted.loadings=load[order(load [, 2]), 2]
Main="Loadings plot for PC2"
xlabs= "Variables loading"
dotplot(sorted.loadings, main=Main, xlab=xlabs, cex=1.5, col="red")


#### look at dimensions by age and sex 

#merge with subject demographics
subjs<-subset(subj_data, select = c(Subject, Sex, Age, Cohort))
combined_data <- cbind(subjs, task.pca$x)
lapply(list(combined_data), str)
head(combined_data)



############## COMPARISONS OF PCA DIMENSION 1 BY AGE AND SEX ##############

#simple correlations with age: PC1
cor.test(combined_data$Age, combined_data$PC1, method = "pearson") #PCA1 correlated
#correlation with age

#base model: no predictors
dim1_model0 <-lm(PC1 ~ 1, data = combined_data)

# Add in age
dim1_model1 <-lm(PC1 ~ Age, data = combined_data)
summary(dim1_model1)
lrtest(dim1_model0, dim1_model1)
#improved fit

#Add sex
dim1_model2 <-lm(PC1 ~ Age + Sex, data = combined_data)
summary(dim1_model2)
lrtest(dim1_model1, dim1_model2)
#no improvement

#Add sex X age interaction
dim1_model3 <-lm(PC1 ~ Age*Sex, data = combined_data)
summary(dim1_model3)
lrtest(dim1_model1, dim1_model3)
#no improvement

#check assumptions of linear model for model 1 (best fit)
gvlma(dim1_model1)
#ok except link function

plot(allEffects(dim1_model3))



############## COMPARISONS OF PCA DIMENSION 2 BY AGE AND SEX ##############

#simple correlations with age: PC2
cor.test(combined_data$Age, combined_data$PC2, method = "pearson") #PCA2 not correlated
# no correlations with age

#base model: no predictors
dim2_model0 <-lm(PC2 ~ 1, data = combined_data)

#Add age
dim2_model1 <-lm(PC2 ~ Age, data = combined_data)
summary(dim2_model1)
lrtest(dim2_model0, dim2_model1)
#no improvement

#Add sex
dim2_model2 <-lm(PC2 ~ Age + Sex, data = combined_data)
summary(dim2_model2)
lrtest(dim2_model1, dim2_model2)
#improvement

#Add sex X age interaction
dim2_model3 <-lm(PC2 ~ Age*Sex, data = combined_data)
summary(dim2_model3)
lrtest(dim2_model2, dim2_model3)
#no improvement

#check assumptions of linear model for model 2 (best fit)
gvlma(dim2_model2)
#ok

plot(allEffects(dim2_model3))

##########################################################################
############################ PLOTS #######################################
##########################################################################

#adjust variable lables for plots
plot_data <-subset_data
plot_data <- plot_data %>%
  rename('Risky Choice' = Risk,
  'Temporal Choice' = Discounting,
  'Risk Affect Score' = RiskAffectScore,
  'Temporal Affect Score' = DiscountingAffectScore)

#pairwise correlation
#plot
corr_viz2 = cor(plot_data, method = c("pearson"))
correl_plot<-ggcorrplot(corr_viz2, method = "circle",  type = "upper")

#PCA
#graph of tasks: dimension 1 and 2
PCA_data2<-subset(plot_data, select = c('Risky Choice', 'Risk Affect Score', 'Temporal Choice', 'Temporal Affect Score', 'Cortisol', 'Testosterone'))

#Calculate the correlation matrix
cor_matrix <- cor(PCA_data2)

#conduct PCA all
task.pca2<-prcomp(na.omit(PCA_data2), center = TRUE, scale = TRUE, retx = TRUE)
summary(task.pca2)


PCA_plot<-fviz_pca_var(task.pca2,
             title = "",
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

#combine graphs
quartz.options(width=10, height=5)
combined_corr_plots <- grid.arrange(correl_plot, PCA_plot, ncol = 2, widths=c(5,4.5))




