library(ggplot2) #needed for graphics
library(dplyr) #needed to summarize
library(lme4) #needed for mixed models
library(car) #needed for type 3 analyses on glmers
library(lmerTest) #needed for type 3 analyses with correct degrees of freedom on lmers
library(emmeans) #needed for post-hoc tests
library(visreg) # needed to visualize residuals
options(contrasts = c("contr.sum","contr.poly"))  #set the options for a type 3 analysis with proper intercept

#user sets the working directory
setwd()

#read in the data
census <- read.csv("Killifish_plasticity.csv")
#look at it
census

#create unique sire-dam info
census$dam_sire <- paste(census$sire,"-",census$dam,sep="")

#subset the data in tea and clear data sets
only_tea <- subset(census, census$water=="tea")
only_clear <- subset(census, census$water=="clear")

#rename the variable "propsolidblue" and "propanyblue" so it is specific to water treatment
colnames(only_tea)[24:29]<-c("propanyred_tea","propanyyellow_tea","propsolidred_tea","propsolidyellow_tea","propsolidblue_tea","propanyblue_tea")
colnames(only_clear)[24:29]<-c("propanyred_clear","propanyyellow_clear","propsolidred_clear","propsolidyellow_clear","propsolidblue_clear","propanyblue_clear")

#remove some descriptor columns to avoid duplicate columns
only_tea <- only_tea[,-c(11:23)] #remove unwanted columns to keep dataframe manageable
only_clear <- only_clear[,-c(1:23)] #remove unwanted columns to keep dataframe manageable

#merge the two datasets
census1b_merge <- merge(only_tea, only_clear, by='dam_sire') #merge by dam-sire combination

#create the solid blue plasticity variable
#this is the difference in the proportion of sons with solid blue anal fins 
#between each unique dam-sire family raised in tea-stained and clear water
census1b_merge$solidblueplasticity <- census1b_merge$propsolidblue_tea-census1b_merge$propsolidblue_clear

#create the solid red plasticity variable
census1b_merge$solidredplasticity <- census1b_merge$propsolidred_tea-census1b_merge$propsolidred_clear

#create the solid yellow plasticity variable
census1b_merge$solidyellowplasticity <- census1b_merge$propsolidyellow_tea-census1b_merge$propsolidyellow_clear

#create the any blue plasticity variable
census1b_merge$anyblueplasticity <- census1b_merge$propanyblue_tea-census1b_merge$propanyblue_clear

#create the any red plasticity variable
census1b_merge$anyredplasticity <- census1b_merge$propanyred_tea-census1b_merge$propanyred_clear

#create the any yellow plasticity variable
census1b_merge$anyyellowplasticity <- census1b_merge$propanyyellow_tea-census1b_merge$propanyyellow_clear

#make certain that sire is a factor
census1b_merge$sire <- factor(census1b_merge$sire)

#----------The plasticity data set is now set. --------------------
#full model with cross and male color pattern and solid blue plasticity
#this model considers all four color patterns and all four crosses
model1 <- lmer(solidblueplasticity~cross*malecolorpattern + (1|sire), data=census1b_merge)
summary(model1) #summary
anova(model1) #type 3 analysis with Sattherwaite degrees of freedom
model1_anova <- anova(model1) #assign to an object
shapiro.test(residuals(model1)) #test the residuals
outlierTest(model1) #test for outliers
Null_model <- lmer(solidblueplasticity~ (1|sire), data=census1b_merge) #stipulate the null model
whmodel1 <-anova(model1, Null_model) #test the full model

#table 3A is the following
anova(model1) #type 3 analysis with Sattherwaite degrees of freedom
anova(model1, Null_model) #test the full model over the null model

#because the residuals are non-normal, we also conducted a simulation where we
#randomized the data among the various treatment groups and re-analyzed the data
#we calculated the proportion of simulations that produced F-values greater than the
#ones that we observed

#-----------------------Calculate the p-value with a simulation ---------------------
#double-check p-values through simulation
cross_F<-vector() #create a vector to hold the F-values for cross
malecolor_F <- vector() #create a vector to hold the F-values for color pattern
cross_color_F <- vector() #create a vector to hold the F-values for the interaction
whole_model_Chi <- vector() #create a vector to hold the Chi-squared values

#run the simulation 10,000 times
for(i in 1:10000){
  tea_sample <- sample(census1b_merge$propsolidblue_tea) #randomly re-assign the tea-values across cross/color
  clear_sample <- sample(census1b_merge$propsolidblue_clear) #randomly re-assign the clear-values across cross/color
  plasticity <- tea_sample - clear_sample #calculate the new plasticity values
  holdmodel <- lmer(plasticity~census1b_merge$cross*census1b_merge$malecolorpattern + (1|census1b_merge$sire)) #conduct model on randomized data
  Null_model <-  lmer(plasticity~ (1|census1b_merge$sire)) #run the null model
   holdanova <- anova(holdmodel) #calculate F-values for simulated data
  model_anova <- anova(holdmodel, Null_model) #calculate Chi-squared values for the full model
  
  cross_F[i] <- holdanova$`F value`[1] #store F-values from the cross
  malecolor_F[i] <- holdanova$`F value`[2] #store F-values from the male color pattern
  cross_color_F[i] <- holdanova$`F value`[3] #store F-values from the interaction
  whole_model_Chi[i] <- model_anova$Chisq[2] #store Chi-values for the whole model
}

#calculate the proportion of F-values and Chi-squared values that were greater than that observed
P_sim_cross= sum(cross_F>model1_anova$`F value`[1])/10000 #proportion as extreme as the F-value we obtained
P_sim_color=sum(malecolor_F>model1_anova$`F value`[2])/10000 #proportion as extreme as the F-value we obtained
P_sim_crossxcolor=sum(cross_color_F>model1_anova$`F value`[3])/10000 #proportion as extreme as the F-value we obtained
P_sim_model=sum(whole_model_Chi>whmodel1$Chisq[2])/10000 #proportion as extreme as the F-value we obtained

#The simulated values below should be similar to those in table 3A.
c(P_sim_cross=P_sim_cross, P_sim_color=P_sim_color,P_sim_crossxcolor=P_sim_crossxcolor,P_sim_model=P_sim_model)

#---------Re-run the analysis for the condensed model with --------------------------------
#---------condensed color pattern and condensed cross -------------------------------------
#---------For this model we considered swamp crosses verses the other 3 crosses involving spring
#--------parents (swamp_other).  We also considered blue sires versus non-blue sires (Color_blue). ------------------------
model2 <- lmer(solidblueplasticity~swamp_other*Color_blue + (1|sire), data=census1b_merge)
summary(model2) #summary of model
model2_anova <-anova(model2) #anova with Sattherwaite degrees of freedom of terms
shapiro.test(residuals(model2)) #test the residuals
outlierTest(model2) #look for outliers
Null_model <- lmer(solidblueplasticity~ (1|sire), data=census1b_merge) #run the null model
whmodel2 <- anova(model2, Null_model) #compare model with the null model

#the table 3B has the following statistical analysis
anova(model2) #anova with Sattherwaite degrees of freedom of terms
anova(model2, Null_model) #compare model with the null model

#because the residuals are non-normal, we also conducted a simulation where we
#randomized the data among the various treatment groups and re-analyzed the data
#we calculated the proportion of simulations that produced F-values greater than the
#ones that we observed

#double-check p-values through simulation
cross_F<-vector() #create vector to hold F-values from cross
malecolor_F <- vector() #create vector to hold F-values from male color pattern
cross_color_F <- vector() #create vector to hold F-values from the interaction
whole_model_Chi <- vector() #create a vector to hold the Chi-squared values

#run the simulation 10,000 times
for(i in 1:10000){
  tea_sample <- sample(census1b_merge$propsolidblue_tea) #randomly sample tea
  clear_sample <- sample(census1b_merge$propsolidblue_clear) #randomly sample clear
  plasticity <- tea_sample - clear_sample #calculate plasticity
  holdmodel <- lmer(plasticity~census1b_merge$swamp_other*census1b_merge$Color_blue + (1|census1b_merge$sire), ) #run the model
  Null_model <-  lmer(plasticity~ (1|census1b_merge$sire)) #run the null model
  holdanova <- anova(holdmodel) #get F-values from Satherwaite for terms
  model_anova <- anova(holdmodel, Null_model) #get chisq values from full model test
  
  cross_F[i] <- holdanova$`F value`[1] #store F-values for cross
  malecolor_F[i] <- holdanova$`F value`[2] #store F-values for male color pattern
  cross_color_F[i] <- holdanova$`F value`[3] #store F-values from interaction
  whole_model_Chi[i] <- model_anova$Chisq[2] #store Chisq values from model
  
}

#calculate the proportion of F-values and Chi-squared values that were greater than that observed
P_sim_cross_cond = sum(cross_F>model2_anova$`F value`[1])/10000 #proportion as extreme as F for cross
P_sim_color_cond = sum(malecolor_F>model2_anova$`F value`[2])/10000 #proportion as extreme as F for color pattern
P_sim_crossxcolor_cond = sum(cross_color_F>model2_anova$`F value`[3])/10000 #proportion as extreme 
P_sim_model_cond=sum(whole_model_Chi>whmodel2$Chisq[2])/10000 #proportion as extreme as the F-value we obtained

#The simulated p-values below should be similar to those in table 3B.
c(P_sim_cross_cond =P_sim_cross_cond, P_sim_color_cond=P_sim_color_cond,
  P_sim_crossxcolor_cond =P_sim_crossxcolor_cond,P_sim_model_cond=P_sim_model_cond)

#----------------Make the graph for solid blue--------------------
#get means on sires*water to create the graph
balanced_sires <- census %>% group_by(sire,malecolorpattern,water,cross, Color_blue) %>%
              summarise(mean_propsolidblue = mean(propsolidblue),
                        mean_propanyblue=mean(propanyblue),
                        mean_propsolidred=mean(propsolidred),
                        mean_propsanyred = mean(propanyred),
                        mean_propsolidyellow=mean(propsolidyellow),
                        mean_propanyyellow = mean(propanyyellow))

#set the labels for the graph
labels <- c(spmxspf = "Spring", swmxspf="Swm x Spf", spmxswf = "Spm x Swf", swmxswf = "Swamp"  )

#set the order on the male color patterns
balanced_sires$malecolorpattern <- factor(balanced_sires$malecolorpattern, levels=c("rb","yb","rr","yy"), labels=c("R-B","Y-B","R-R","Y-Y"))

#set the order on the crosses
balanced_sires$cross <- factor(balanced_sires$cross, levels=c("spmxspf", "swmxspf","spmxswf", "swmxswf" ))

#create the graph for solid blue plasticity
Figure_3 <- ggplot(balanced_sires, aes(x=water, y=mean_propsolidblue, color=malecolorpattern, group=sire)) + 
  geom_point(color='black') + 
  geom_line() + 
  scale_color_manual(values=c("blue4","royalblue2","red","gold"))+
  facet_grid(.~cross, labeller=labeller(cross = labels)) + 
  labs(col="Color Pattern") +
  xlab("Water Treatment") +
  ylab("Proportion Solid Blue Sons") + theme_gray() + 
  theme_bw()+
  theme(
    legend.key = element_rect(colour = "white", fill = NA),
    legend.text=element_text(size=12), 
  )

#The graph blow should be the same 
#as that shown in figure 3. 
Figure_3

#compare the magnitude of plasticity between swamp versus non-swamp crosses
cross_means <- census1b_merge %>%
  group_by(swamp_other) %>%
  summarise(solid_blue_plasticity_cross = mean(solidblueplasticity))

#report numbers
cross_means

#magnitude of plasticity difference
#reference to this number can be found in paragraph 4 of the results
cross_means[2,2]/cross_means[1,2]


#--------------------------Now, re-run the whole thing for plasticity in 'any blue
#full model with cross and male color pattern and any blue plasticity
#this model considers all four color patterns and all four crosses
#original full model for any blue plasticity
anyblue_original_model <- lmer(anyblueplasticity~cross*malecolorpattern + (1|sire), data=census1b_merge) #the full model
Null_model <- lmer(anyblueplasticity~(1|sire), data=census1b_merge) #Null model
whmodel_test <- anova(anyblue_original_model, Null_model) #anova comparing full and null model
anyblue_original_model_anova <- anova(anyblue_original_model) #get analysis on terms with Sattherwaite
shapiro.test(residuals(anyblue_original_model)) #check the residuals
outlierTest(anyblue_original_model) #check the outliers

#table 3C should be the same as below
anova(anyblue_original_model)
anova(anyblue_original_model, Null_model)

#because the residuals are non-normal, we also conducted a simulation where we
#randomized the data among the various treatment groups and re-analyzed the data
#we calculated the proportion of simulations that produced F-values greater than the
#ones that we observed

#double-check p-values through simulation
cross_F<-vector() #create a vector to hold F-values on cross
malecolor_F <- vector() #create a vector to hold F-values on color pattern
cross_color_F <- vector() #create a vector to hold F-values on interaction
whole_model_Chi <- vector() #create a vector to hold the Chi-squared values

#run simulation
for(i in 1:10000){
  tea_sample <- sample(census1b_merge$propanyblue_tea) #randomize tea-stained values across treatments
  clear_sample <- sample(census1b_merge$propanyblue_clear) #randomize clear-water values across treatments
  plasticity <- tea_sample - clear_sample #calculate plasticity
  holdmodel <- lmer(plasticity~census1b_merge$cross*census1b_merge$malecolorpattern + (1|census1b_merge$sire)) #conduct the model
  Null_model <-  lmer(plasticity~ (1|census1b_merge$sire)) #make the null model
  holdanova <- anova(holdmodel) #hold the results of the type 3 analysis of terms
  model_anova <- anova(holdmodel, Null_model) #hold the results for the test of the entire model
  cross_F[i] <- holdanova$`F value`[1] #assign the F-value for cross
  malecolor_F[i] <- holdanova$`F value`[2] #assign the F-value for male color
  cross_color_F[i] <- holdanova$`F value`[3] #assign the F-value for the interactions
  whole_model_Chi[i] <- model_anova$Chisq[2] #assign the Chi-squared value for the whole model
}

#calculate the proportion of F-values and Chi-squared values that were greater than that observed
anyblue_psim_cross = sum(cross_F>anyblue_original_model_anova$`F value`[1])/10000
anyblue_psim_color = sum(malecolor_F>anyblue_original_model_anova$`F value`[2])/10000
anyblue_psim_crossxcolor = sum(cross_color_F>anyblue_original_model_anova$`F value`[3])/10000
anyblue_psim_model = sum(whole_model_Chi>whmodel_test$Chisq[2])/10000 #proportion as extreme as the F-value we obtained

#the simulated p-values in table 3C should be similar the ones listed below. 

c(anyblue_psim_cross =anyblue_psim_cross, anyblue_psim_color =anyblue_psim_color, 
  anyblue_psim_crossxcolor =anyblue_psim_crossxcolor, anyblue_psim_model =anyblue_psim_model)

#----------------Run the Analysis for Any blue plasticity condensed model------------------------------
anyblue_condense_model <- lmer(anyblueplasticity~swamp_other*Color_blue+ (1|sire), data=census1b_merge)
summary(anyblue_condense_model) #look at the summary
whmodel_test <- anova(anyblue_condense_model, Null_model) #test the full model
whmodel_test #report it
anyblue_condense_model_anova <- anova(anyblue_condense_model) #look at the test of model terms
anyblue_condense_model_anova #report it
shapiro.test(residuals(anyblue_condense_model)) #test the residuals
outlierTest(anyblue_condense_model) #look for outliers

#table 3D is the following
anova(anyblue_condense_model) #type 3 analysis
anova(anyblue_condense_model, Null_model) #whole model

#because the residuals are non-normal, we also conducted a simulation where we
#randomized the data among the various treatment groups and re-analyzed the data
#we calculated the proportion of simulations that produced F-values greater than the
#ones that we observed

#double-check p-values through simulation
cross_F<-vector()
malecolor_F <- vector()
cross_color_F <- vector()
whole_model_Chi <- vector() #create a vector to hold the Chi-squared values

#run the simulation 10,000 times

for(i in 1:10000){
  
  tea_sample <- sample(census1b_merge$propanyblue_tea) #randomize tea-sample
  clear_sample <- sample(census1b_merge$propanyblue_clear) #randomize clear-sample
  plasticity <- tea_sample - clear_sample #calculate plasticity
  holdmodel <- lmer(plasticity~census1b_merge$swamp_other*census1b_merge$Color_blue + (1|census1b_merge$sire)) #merge with descriptors and perform the model
  Null_model <-  lmer(plasticity~ (1|census1b_merge$sire)) #make the null model
  holdanova <- anova(holdmodel) #hold the results of the type 3 analysis of terms
  model_anova <- anova(holdmodel, Null_model) #hold the results of the overall model
  
  cross_F[i] <- holdanova$`F value`[1] #assign the F-value for cross
  malecolor_F[i] <- holdanova$`F value`[2] #assign the F-value for male color
  cross_color_F[i] <- holdanova$`F value`[3] #assign the F-value for the interaction
  whole_model_Chi[i] <- model_anova$Chisq[2] #assign the Chi-squred value for the whole model
  
}

#calculate the proportion of F-values and Chi-squared values that were greater than that observed
anyblue_psim_cross_cond = sum(cross_F>anyblue_condense_model_anova$`F value`[1])/10000 #proportion of simulations as extreme as F for cross
sum(malecolor_F>anyblue_condense_model_anova$`F value`[2])/10000 #proportion of simulations as extreme as F for color
sum(cross_color_F>anyblue_condense_model_anova$`F value`[3])/10000 #proportion of simulations as extreme as F for interaction
sum(whole_model_Chi>whmodel2$Chisq[2])/10000 #proportion as extreme as the F-value we obtained

#the simulated p-values below should be similar to those in table 3D
c(anyblue_psim_cross_cond = anyblue_psim_cross_cond, anyblue_psim_color_cond =anyblue_psim_color_cond, 
  anyblue_psim_crossxcolor_cond =anyblue_psim_crossxcolor_cond, anyblue_psim_model_cond =anyblue_psim_model_cond)


#make the graph
#this is the graph for supplemental figure 2
Supp_Fig2 <- ggplot(balanced_sires, aes(x=water, y=mean_propanyblue, color=malecolorpattern, group=sire)) + 
  geom_point(color='black') + 
  geom_line() + 
  scale_color_manual(values=c("blue4","royalblue2","red","gold"))+
  facet_grid(.~cross, labeller=labeller(cross = labels)) + 
  labs(col="Color Pattern") +
  xlab("Water Treatment") + 
  ylab("Proportion Any Blue Sons") + theme_gray() + 
  theme_bw()+
  theme(
    legend.key = element_rect(colour = "white", fill = NA),
    legend.text=element_text(size=12)
  )

#look at it
Supp_Fig2 

#the remaining text provides statistical support for the following statement
#There was no evidence of genetic variation in phenotypic plasticity in 
#the proportion of male offspring with any red, solid yellow, or any yellow 
#anal fins (Table 4, Table S2). #The data below are for table S2.

#------------Analyze Red plasticity-----------------------------
model3 <- lmer(solidredplasticity~cross*malecolorpattern + (1|sire), data=census1b_merge) #run the full model
summary(model3) #look at the summary
anova(model3) #look at the test of terms
Null_model <- lmer(solidredplasticity~(1|sire), data=census1b_merge) #run the null model
shapiro.test(residuals(model3)) #look at the residuals
outlierTest(model3) #test for outliers
leveneTest(residuals(model3), group=census1b_merge$cross) #look for heteroscedastic variances
leveneTest(residuals(model3), group=census1b_merge$malecolorpattern)#look for heteroscedastic variances
visreg(model3) #look at residuals
anova(model3, Null_model, test="F") #test the full model
emmeans(model3, pairwise~"cross")
emmeans(model3, pairwise~"cross*malecolorpattern")

#supplemental table 2A has the following data
anova(model3) #look at the test of terms
anova(model3, Null_model, test="F") #test the full model

#the residuals were normal and there was no evidence of heteroscedasticity
#no need to simulate p-value

#the following graph does not appear in the manuscript or supplemental material
#create the graph for solid red plasticity
ggplot(balanced_sires, aes(x=water, y=mean_propsolidred, color=malecolorpattern, group=sire)) + 
  geom_point(color='black') + 
  geom_line() + 
  scale_color_manual(values=c("blue4","red","royalblue2","gold"))+
  facet_grid(.~cross, labeller=labeller(cross = labels)) + 
  #guides(fill = guide_legend(title = 'male color pattern')) +
  labs(col="Color Pattern") +
  xlab("Water Treatment") +
  ylab("Proportion Solid Red Sons") + theme_gray() + 
  theme_bw()+
  theme(
    legend.key = element_rect(colour = "white", fill = NA),
    legend.text=element_text(size=12), 
  )


#--------------------------Run the analysis for any plasticity
model4 <- lmer(anyredplasticity~cross*malecolorpattern + (1|sire), data=census1b_merge) #make the full model
summary(model4 ) #look at the summary
anova(model4 ) #look at test of terms
Null_model <- lmer(anyredplasticity~(1|sire), data=census1b_merge) #null model
anova(model4, Null_model) #test the full model
visreg(model4 ) #look at the residuals
leveneTest(residuals(model4), group=census1b_merge$cross) #look for heteroscedastic variances
leveneTest(residuals(model4), group=census1b_merge$malecolorpattern)#look for heteroscedastic variances
outlierTest(model4) #check for outliers
shapiro.test(residuals(model4 )) #test that residuals are normal

#supplemental table 2B has the following data
anova(model4 ) #look at test of terms
anova(model4, Null_model) #test the full model

#no deviations - no need to simulate the p-value

#Make the graph 
#the following graph does not appear in the manuscript or supplemental materials
ggplot(balanced_sires, aes(x=water, y=mean_propsanyred, color=malecolorpattern, group=sire)) + 
  geom_point(color='black') + 
  geom_line() + 
  scale_color_manual(values=c("blue4","royalblue2","red","gold"))+
  facet_grid(.~cross, labeller=labeller(cross = labels)) + 
  #guides(fill = guide_legend(title = 'male color pattern')) +
  labs(col="Color Pattern") +
  xlab("Water Treatment") +
  ylab("Proportion Any Red Sons") + theme_gray() + 
  theme_bw()+
  theme(
    legend.key = element_rect(colour = "white", fill = NA),
    legend.text=element_text(size=12), 
  )


#run the model on solid yellow plasticity
model5 <- lmer(solidyellowplasticity~cross*malecolorpattern + (1|sire), data=census1b_merge)
summary(model5) #look at the summary
Null_model <- lmer(solidyellowplasticity~(1|sire), data=census1b_merge) #make the null model
anova(model5, Null_model) #test the full model
anova(model5) #test the term effects
visreg(model5) #look at the residuals
shapiro.test(residuals(model5)) #test the residuals for normality
outlierTest(model5) #look for outliers
leveneTest(residuals(model5), group=census1b_merge$cross) #look for heteroscedastic variances
leveneTest(residuals(model5), group=census1b_merge$malecolorpattern)#look for heteroscedastic variances

#supplemental table 2C
anova(model5) #test the term effects
anova(model5, Null_model) #test the full model

#no significant effects, no deviations from normality
#the following graph does not appear in the manuscript or supplemental materials
#make the graph
ggplot(balanced_sires, aes(x=water, y=mean_propsolidyellow, color=malecolorpattern, group=sire)) + 
  geom_point(color='black') + 
  geom_line() + 
  scale_color_manual(values=c("blue4","royalblue2","red","gold"))+
  facet_grid(.~cross, labeller=labeller(cross = labels)) + 
  #guides(fill = guide_legend(title = 'male color pattern')) +
  labs(col="Color Pattern") +
  xlab("Water Treatment") +
  ylab("Proportion Solid Yellow Sons") + theme_gray() + 
  theme_bw()+
  theme(
    legend.key = element_rect(colour = "white", fill = NA),
    legend.text=element_text(size=12), 
  )


model6 <- lmer(anyyellowplasticity~cross*malecolorpattern + (1|sire), data=census1b_merge)
summary(model6) #look at the summary
Null_model <- lmer(anyyellowplasticity~(1|sire), data=census1b_merge) #make the null model
anova(model6, Null_model) #test the full model
anova(model6) #test the term effects
shapiro.test(residuals(model6)) #test the residuals for normality
outlierTest(model6) #look for outliers

#there were no significant effects here
#we did run a simulation and found no evidence that an alternate approach would change
#the fixed effects of the model

#make the graph
#this graph does not appear in the manuscript or supplemental material
ggplot(balanced_sires, aes(x=water, y=mean_propanyyellow, color=malecolorpattern, group=sire)) + 
  geom_point(color='black') + 
  geom_line() + 
  scale_color_manual(values=c("blue4","royalblue2","red","gold"))+
  facet_grid(.~cross, labeller=labeller(cross = labels)) + 
  labs(col="Color Pattern") +
  xlab("Water Treatment") +
  ylab("Proportion Solid Yellow Sons") + theme_gray() + 
  theme_bw()+
  theme(
    legend.key = element_rect(colour = "white", fill = NA),
    legend.text=element_text(size=12), 
  )
