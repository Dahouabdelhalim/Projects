#the analyses below address the effects of treatments and their 
#interactions on the proportion of males with solid red and solid yellow
#anal fins as well as males expressing any red or any yellow coloration
library(ggplot2) #graphics
library(dplyr) #summarizing data
library(lme4) #mixed models
library(car) #type 3 models, particularly glmers
library(lmerTest) #needed for some tests of lmers
options(contrasts = c("contr.sum","contr.poly"))  #set contrasts for type 3 analyses
library(RVAideMemoire) #needed for tests of overdispersion
library(emmeans) #needed for post-hoc tests
library(ggpubr) #needed to add correlations to the plot

#user sets the working directory
setwd()

#read in the data
#read in the data
census <- read.csv("Killifish_plasticity.csv")
census

#make the sires as factors (not numeric data)
census$sire <- as.factor(census$sire)

#make the dams unique
census$newdam <- paste(census$sire,census$dam, sep="")

#from the results "Across all families, there was a negative correlation 
#between the proportion of sons expressing any red versus any yellow 
#coloration (R = -0.95, P < 0.001, df = 122)
cor.test(census$propanyyellow, census$propanyred)
#the following graph and statistics are not in the manuscript
#these are the basis for the statement: 
#and the same pattern was upheld in all treatment combinations"

ggplot(census, aes(x=propanyred, y=propanyyellow)) + 
  geom_point() + 
  facet_grid(cross~water) + 
  xlab("proportion sons expressing any red") + 
  ylab("proportion sons expressing any yellow") + 
  geom_smooth(method="lm") + 
  stat_cor() +
  theme_bw()

#--------------Analysis of frequencies for the proportion of males with-------------- 
#---------------------------solid red anal fins------------------

#specify the null model
null_model_noeffect <- glmer(propsolidred ~ (1|sire/dam), family=binomial(link="logit"),  nAGQ=0, weights = totalmaleswithcolor, data=census)

#specifiy the full model
full_red <- glmer(propsolidred ~ water*cross*malecolorpattern + (1|sire/dam), family=binomial(link="logit"),  nAGQ=0, weights = totalmaleswithcolor, data=census)
summary(full_red) #summary
plot(full_red) #check the residuals

#line 1 of table 4 has the following statistics
Anova(full_red, type=3) #type 3 anova
anova(full_red, null_model_noeffect) #test of full model
overdisp.glmer(full_red) #check for over dispersion

#------Line 2 of table 4 - the analysis on Any red ---------
#specify the null model
null_model_noeffect <- glmer(propanyred ~ (1|sire/dam), family=binomial(link="logit"),  nAGQ=0, weights = totalmaleswithcolor, data=census)
#specify the full model
full_red_any <- glmer(propanyred ~ water*cross*malecolorpattern + (1|sire/dam), nAGQ=0,family=binomial(link="logit"), weights = totalmaleswithcolor, data=census) 
summary(full_red_any) #look at the summary
plot(full_red_any) #look at the residuals

#line 2 of table 4 has the following data
Anova(full_red_any, type=3) #test of individual effects
anova(full_red_any, null_model_noeffect) #test of the full model
overdisp.glmer(full_red_any) #test for overdispersion

#--------Line 3 of table 4 model for solid yellow ----------------------
#specify the null model
null_model_noeffect <- glmer(propsolidyellow ~ (1|sire/dam), family=binomial(link="logit"),  nAGQ=0, weights = totalmaleswithcolor, data=census)
#specify the full model
full_yel <- glmer(propsolidyellow ~ water*cross*malecolorpattern + (1|sire/dam), family=binomial(link="logit"),  nAGQ=0, weights = totalmaleswithcolor, data=census)
summary(full_yel) # look at the summary
#look at the residuals
plotresid(full_yel)

#line 3 of table 4 has the following statistics 
Anova(full_yel, type=3) #test of individual effects
anova(null_model_noeffect, full_yel ) #test of full model
overdisp.glmer(full_yel) #check for overdispersion

#--------Line 4 of table 4 --- analysis of males with any yellow --------------
#specify the null model
null_model_noeffect <- glmer(propanyyellow ~ (1|sire/dam), family=binomial(link="logit"),  nAGQ=0, weights = totalmaleswithcolor, data=census)
#specify the full model
any_yel <- glmer(propanyyellow ~ water*cross*malecolorpattern + (1|sire/dam), family=binomial(link="logit"),  nAGQ=0, weights = totalmaleswithcolor, data=census)
summary(any_yel) #look at the summary
plot(any_yel)#look at the residuals 

#line 4 of table 4 has the following statistics
Anova(any_yel, type=3) #individual tests of effects
anova(null_model_noeffect, any_yel) #test of the full model
overdisp.glmer(any_yel) #check for overdispersion

#-----------------------------------------------------------------------------------
#statistical support for the following statement: 
#The proportion of sons with red coloration was ~1.9X higher for red than for yellow sires.
any_means <- census %>%
  group_by(malecolorpattern) %>%
  summarise(mean_any_red = mean(propanyred),
            mean_solidred=mean(propsolidred),
            mean_any_yel=mean(propanyyellow),
            mean_solidyel = mean(propanyyellow))

#look at it
any_means

#comparison of effects
mean(0.571,0.618)/mean(0.293,0.319) #red statements red sons 1.9X higher for red sires
mean(0.613,0.665)/mean(0.388,0.355) #yellow statements yellow sons 1.6X higher for yellow sires

#support for the following statements: 
#In all cases, there were no significant differences between yellow sires 
#(e.g., Y-Y vs. Y-B) nor between red sires (e.g., R-R vs. R-B). In most cases, 
#there were significant differences between red and yellow sires (e.g., R-R vs. 
#Y-Y, R-R vs. Y-B, R-B vs. Y-Y, R-B vs. Y-B). 

emmeans(full_red, pairwise ~ 'malecolorpattern')
emmeans(full_red_any, pairwise ~ 'malecolorpattern')
emmeans(full_yel, pairwise ~ 'malecolorpattern')
emmeans(any_yel, pairwise ~ 'malecolorpattern')


#support for the following statements: 
#We also found that population of origin also had 
#substantial effects on the proportion of sons with solid 
#yellow, any yellow, and any red fins (Table 4, Fig. S3-6). 
#This finding was not predicted a priori but is consistent 
#with the hypothesis of a locus of large effect.  These 
#differences occurred because offspring from swamp crosses 
#were less likely to express solid yellow and any yellow 
#compared to offspring from spring crosses (P < 0.002 in all 
#post-hoc tests, Figs. S4,6). Likewise, offspring from spring 
#crosses were less likely to express any red compared to offspring 
#from swamp crosses (P < 0.015 post-hoc test, Fig. S4).  T
emmeans(full_red_any, pairwise ~ 'cross')
emmeans(full_yel, pairwise ~ 'cross')
emmeans(any_yel, pairwise ~ 'cross')


