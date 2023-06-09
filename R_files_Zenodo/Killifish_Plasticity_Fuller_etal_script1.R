#this script summarizes the counts of the various color morphs across treatments and 
#analyzes the effect of treatments on the proportion of males with solid blue anal fins 
#and the proportion of males with any blue on their anal fins.  This corresponds to the 
#results found in the first three paragraph of the results. 
library(ggplot2) #make the graphs
library(dplyr) #summarize data
library(lme4) #perform mixed models
library(car) #perform type 3 models not available from lmerTest (i.e., glmer models)
library(lmerTest) #perform type 3 models using appropriate degrees of freedom for lmers
library(RVAideMemoire) #check for overdispersion
library(emmeans) #post-hoc treatment comparisons
library(visreg) #visually examine residuals
library(janitor) #add marginal totals to tables from dplyr

#set the contrasts to properly analyze interactions with unbalanced design for type 3 analyses
options(contrasts = c("contr.sum","contr.poly"))  

#user sets the working directory
setwd()

#read in the data
census <- read.csv("Killifish_plasticity.csv")
census
 
#make the sires as factors (not numeric data)
census$sire <- as.factor(census$sire)

#make each dam unique
census$newdam <- as.factor(paste(census$sire, census$dam,sep=""))

#-------------------make Table #1------------------------------------------
table1_counts <- census%>%
  group_by(cross, water) %>%
  summarise(
    solid_blue = sum(totalsolidblue),
    solid_red=sum(totalsolidred),
    solid_yellow=sum(totalsolidyellow),
    solid_orange=sum(totalsolidorange),
    red_blue_combo=sum(totalred_blue_combo),
    yellow_blue_combo=sum(totalyellow_blue_combo),
    orange_blue_combo=sum(totalorange_blue_combo)
  ) %>%
  adorn_totals("row")  %>%
  adorn_totals("col")

#look at the table
table1_counts
#--------------------------------------------------------
#the next bit of code uses dplyr to compare the proportion of sons
#expressing solid blue anal fins and any blue on their anal fins
#these results can be found in the second paragraph of the results
#----------the effects are examined at the level of sire by water treatment

sire_means <- census %>% group_by(Color_blue,water,swamp_other, sire) %>%
  summarise(
    mean_solidblue = mean(propsolidblue),
    mean_anyblue = mean(propanyblue)
  )


#get means for water effects
Water_effects <- sire_means %>% group_by(water) %>%
  summarise(mean_solidbuesire = mean(mean_solidblue))

#3.195 times more likely to be blue in tea-stained water
Water_effects[2,2]/Water_effects[1,2]  #3.2X more likely to be blue when raised in tea-stained

#3X more likely to be blue from blue dads
#get means for blue dads
Color_effects <- sire_means %>% group_by(Color_blue) %>%
  summarise(mean_solidbuesire = mean(mean_solidblue))

#3X more likely to be solid blue when dad is blue
Color_effects[1,2]/Color_effects[2,2]  #3X more likely to be blue when dad is blue

#get the effects for swamps
Swamp_effects <- sire_means %>% group_by(swamp_other) %>%
  summarise(mean_solidbuesire = mean(mean_solidblue))

Swamp_effects[2,2]/Swamp_effects[1,2] #swamp fish 7.4X more likely to be blue than others

#----------Interaction color * water ------------------------------
cross_water <- sire_means %>% group_by(swamp_other, water) %>%
  summarise(mean_solidbuesire = mean(mean_solidblue))

cross_water #look at it 
cross_water[4,3]/mean(0.0168,0.0185,0.0350) #swamp in tea 13.5X more likely to be blue than others
#----------Interaction color * water ------------------------------
color_water <- sire_means %>% group_by(Color_blue, water) %>%
  summarise(mean_solidbuesire = mean(mean_solidblue))

#males are 4X more likely to express blue when dad is blue and raised in tea-stained water
color_water[2,3]/(mean(0.0271,0.0158,0.0303))

#---------Interaction of cross*WATER------------------
cross_water <- sire_means %>% group_by(swamp_other, water) %>%
  summarise(mean_solidbuesire = mean(mean_solidblue))

#males 13.5X more likely to express blue when raised in tea-stained water and have
#swamp parents
cross_water[4,3]/mean(0.0168,0.0185,0.0350)

#------------------Third paragraph of the results-------------------


#-----------Analyze solid blue anal fins-------------------------------------

#Null Model
Null_model <- glmer(propsolidblue~ (1|sire/newdam), nAGQ = 0, family=binomial(link="logit"), weights = totalmaleswithcolor, data=census)

#original model with 4 cross types, 4 color patterns, and 2 water types
orig_model <- glmer(propsolidblue~ cross*malecolorpattern*water + (1|sire/newdam), nAGQ = 0, family=binomial(link="logit"), weights = totalmaleswithcolor, data=census)
#look at the summary
summary(orig_model)
#test the terms
Anova(orig_model, type=3)
#check for overdispersion
overdisp.glmer(orig_model)
#look at the residuals
plot(orig_model)

#condense the sire color pattern so that it is blue versus not blue for the sires
condense_color_model <- glmer(propsolidblue~ cross*Color_blue*water + (1|sire/newdam), nAGQ = 0, family=binomial(link="logit"), weights = totalmaleswithcolor, data=census)
#look at the summary
summary(condense_color_model)
#test the terms
Anova(condense_color_model, type=3)
#check overdispersion
overdisp.glmer(condense_color_model)
#plot residuals
plot(condense_color_model)

#condense cross
#consider populations as a function of swamp versus crosses with spring animals
condense_cross_model <- glmer(propsolidblue~ swamp_other*malecolorpattern*water + (1|sire/newdam), nAGQ = 0, family=binomial(link="logit"), weights = totalmaleswithcolor, data=census)
#look at summary
summary(condense_cross_model)
#test terms
Anova(condense_cross_model, type=3)
#check overdispersion
overdisp.glmer(condense_cross_model)
#plot residuals
plot(condense_cross_model)


#condense cross, condense color pattern, include all interactions
condense_both_full_model <- glmer(propsolidblue~ swamp_other*Color_blue*water + (1|sire/newdam), nAGQ = 0, family=binomial(link="logit"), weights = totalmaleswithcolor, data=census)
#look at summary
summary(condense_both_full_model)
#test the terms
Anova(condense_both_full_model, type=3)
#check overdispersion
overdisp.glmer(condense_both_full_model)
#look at residuals
plot(condense_both_full_model)
#test model versus null
anova(condense_both_full_model, Null_model)


#condense cross, condense color pattern
#remove non significant interactions - 3 way interaction and colorxwater
condense_both_removensint <- glmer(propsolidblue~ swamp_other+Color_blue+water+swamp_other:Color_blue + swamp_other*water  + (1|sire/newdam), nAGQ = 0, family=binomial(link="logit"), weights = totalmaleswithcolor, data=census)
#look at summary
summary(condense_both_removensint )
#test terms
Anova(condense_both_removensint , type=3)
#check overdispersion
overdisp.glmer(condense_both_removensint )
#look at residuals
plot(condense_both_removensint )
#test overall model
anova(condense_both_removensint , Null_model)

#remove 3 way interaction, condense cross, condense color patter
condense_both_remove3way <- glmer(propsolidblue~ swamp_other+Color_blue+water+swamp_other:Color_blue + swamp_other:water + Color_blue:water  + (1|sire/newdam), nAGQ = 0, family=binomial(link="logit"), weights = totalmaleswithcolor, data=census)
#look at summary
summary(condense_both_remove3way)
#test terms
Anova(condense_both_remove3way, type=3)
#check overdispersion
overdisp.glmer(condense_both_remove3way)
#look at residuals
plot(condense_both_remove3way)
#test overall model
anova(condense_both_remove3way, Null_model)

#only additive model
condense_additive_model <- glmer(propsolidblue~ swamp_other+Color_blue+water + (1|sire/newdam), nAGQ = 0, family=binomial(link="logit"), weights = totalmaleswithcolor, data=census)
#look at summary
summary(condense_additive_model)
#test terms
Anova(condense_additive_model, type=3)
#check overdispersion
overdisp.glmer(condense_additive_model)
#look at residuals
plot(condense_additive_model)
#test overall model
anova(condense_additive_model, Null_model)

#Determine which model is best
#compare models with AIC
library(AICcmodavg)

#set up the model list
cand.models <- list()

#assign the models
cand.models[[1]] <- orig_model
cand.models[[2]] <- condense_color_model 
cand.models[[3]] <- condense_both_full_model
cand.models[[4]]<-Null_model
cand.models[[5]]<- condense_cross_model
cand.models[[6]] <- condense_both_removensint 
cand.models[[7]] <- condense_both_remove3way
cand.models[[8]] <- condense_additive_model

#name the models
mod.names <-c("Original","Condense_blue","Condense blue and crosses","Null model","Condense crosses","Condense + remove ns int", "Condense + remove 3way", "Additive Only")
#examine the models
#The Table below is Supplemental Table 1a
#References to these results in third paragraph of results
aictab(cand.set = cand.models, modnames = mod.names)

#print the best one
#The table below is Table 2A in the main text
#References to these results in third paragraph of results
Anova(condense_both_removensint , type=3)
#check overdispersion
overdisp.glmer(condense_both_removensint )
#test overall model
anova(condense_both_removensint , Null_model)

#-----------------------------------------------------------------
#--------------We performed the same analysis on the proportion of males with 'any blue'--------------------------------
#--------------First, we look at the null model
Null_anyblue <- glmer(propanyblue~  (1|sire/newdam), nAGQ = 0, family=binomial(link="logit"), weights = totalmaleswithcolor, data=census)

#full original model
originalmodel_anyblue <- glmer(propanyblue~ cross*malecolorpattern*water + (1|sire/newdam), nAGQ = 0, family=binomial(link="logit"), weights = totalmaleswithcolor, data=census)
summary(originalmodel_anyblue) #look at the summary
Anova(originalmodel_anyblue, type=3) #look at the type 3 analysis
anova(originalmodel_anyblue,Null_anyblue) #look at the support for the full model
overdisp.glmer(originalmodel_anyblue) #check for overdispersion
plotresid(originalmodel_anyblue) #look at the residuals

#condense color pattern
#repeat the analysis but condense the color pattern to blue versus no-blue
cond_color_any <- glmer(propanyblue~ cross*Color_blue*water + (1|sire/dam), nAGQ = 0,family=binomial(link="logit"), weights = totalmaleswithcolor, data=census)
summary(cond_color_any) #look at the summary
Anova(cond_color_any, type=3) #look at the type 3 analysis
anova(cond_color_any,Null_anyblue) #look at the significance for the full model
overdisp.glmer(cond_color_any)  #check for overdispersion

#condense crosses to swamp versus no swamp
cond_cross_any <- glmer(propanyblue~ swamp_other*malecolorpattern*water + (1|sire/dam), nAGQ = 0,family=binomial(link="logit"), weights = totalmaleswithcolor, data=census)
summary(cond_cross_any) #look at the summary
Anova(cond_cross_any, type=3) #look at the type 3 analysis
anova(cond_cross_any,Null_anyblue) #look at the significance for the full model
overdisp.glmer(cond_cross_any) #check for overdispersion

#condense both color pattern and crosses
cond_cross_color_any <- glmer(propanyblue~ swamp_other*Color_blue*water + (1|sire/dam), nAGQ = 0,family=binomial(link="logit"), weights = totalmaleswithcolor, data=census)
summary(cond_cross_color_any) #check the summary
Anova(cond_cross_color_any, type=3) #check the type 3 anova
anova(cond_cross_color_any,Null_anyblue) #look at the significance for the full model
overdisp.glmer(cond_cross_color_any) #check for overdispersion

#Because these models converged with many significant interactions, we did not consider any further models.

#make the list of candidate models
cand.models_any <- list()
cand.models_any [[1]] <- Null_anyblue
cand.models_any [[2]] <- cond_cross_color_any 
cand.models_any [[3]]<- cond_cross_any 
cand.models_any[[4]]<- cond_color_any
cand.models_any[[5]] <- originalmodel_anyblue

#give the models names
mod.names_any <-c("Null model","condense cross and condense color", "condense cross","condense color","original")

#the data below can be found in supplemental table 1B
aictab(cand.set = cand.models_any , modnames = mod.names_any)

#print the best model 
#the table below is found in table 2B in the text
Anova(cond_cross_any, type=3) #look at the type 3 analysis
anova(cond_cross_any,Null_anyblue) #look at the significance for the full model
overdisp.glmer(cond_cross_any) #check for overdispersion

