library(car)
library(lsmeans)
library(MASS)
library(ggplot2)
library(lme4)
library(lmerTest)
library(scales)
library(MuMIn)
library(dplyr)

####MANOVA tests####
manova_flume_level<-manova(cbind(density, Richness, Biomass, OM, o2_L_H)~Flume_Treatment, data=Flume_Main_DB_23)
summary(manova_flume_level)
summary.aov(manova_flume_level)
manova_basketlevel<-manova(cbind(density, Richness, Biomass, OM, o2_L_H)~Basket_Treatment, data=Flume_Main_DB_23)
summary(manova_basketlevel)
summary.aov(manova_basketlevel)

###Flume Level Density comparing control to engineer generalized linear mixed model###
m1a <- glmer(density ~  Flume_Treatment + (1 | Flume), data = Flume_Main_DB_23, family = poisson)
##reduced
m1reduced <- glmer(density ~ 1 + (1 | Flume), data = Flume_Main_DB_23, family = poisson)
###Likleyhood ratio test comparing the full to reduced model for density at the flume level###
anova(m1a,m1reduced,test = 'LRT')

###Summarize model, but suboptimal###
summary(m1reduced)
summary(m1a)
###this give post-hoc comparison comparing the groups###
pairs(emmeans(m1a, ~ Flume_Treatment))


###LRT tests###
###anova between full and reduced model###
####correct order is important place full model first###
anova(m1a, m1reducedtest = 'LRT')
####
Anova(m1a,type = 3)
####
drop1(m1a,test="Chisq")

###Flume Level Richness comparing control to engineer generalized linear mixed model###
RM1 <- glmer(Richness ~  Flume_Treatment + (1 | Flume), data = Flume_Main_DB_23, family = poisson)
summary(RM1)
##reduced flume level richness model for LRT###
mrreduced <- glmer(Richness ~ 1 + (1 | Flume), data = Flume_Main_DB_23, family = poisson)
anova(RM1,mrreduced,test = 'LRT')
###
pairs(emmeans(RM1, ~ Flume_Treatment))

####flume level biomass comparing control to engineer linear mixed effect model###
bio = lmer(log(Biomass+0.01) ~ Flume_Treatment + (1|Flume) , data = Flume_Main_DB_23)
summary(bio)
###Flume level POM comparing control to engineer linear mixed effect model###
OM = lmer(log(OM) ~ Flume_Treatment + (1|Flume) , data = Flume_Main_DB_23)
summary(OM)
anova(OM, ddf = 'Kenward-Roger') 

###Flume level Respiration comparing control to engineer linear mixed effect model###
resp=lmer(o2_L_H ~ Flume_Treatment + (1|Flume) , data = Flume_Main_DB_23)
summary(resp)
anova(resp, ddf = 'Kenward-Roger') 
lsmeans(resp, pairwise ~ 'Flume_Treatment', adjust = 'none')
###density as a function of OM and treatment##
m2a <- glmer(density ~  OM+F_Treat + (1 | Flume), data = Flume_Main_DB_23, family = poisson)
###interaction###
m2a <- glmer(density ~  OM*F_Treat + (1 | Flume), data = Flume_Main_DB_23, family = poisson)
##get basket level comparisons for GLMM##
##density##
bden <- glmer(density ~  Basket_Treatment + (1 | Flume), data = Flume_Main_DB_23, family = poisson(link = "log"))
summary(bden)
##reduced basket level density model for LRT###
bdenreduced <- glmer(density ~ 1 + (1 | Flume), data = Flume_Main_DB_23, family = poisson)
anova(bden,bdenreduced,test = 'LRT')
pairs(emmeans(bden, ~ Basket_Treatment))
##richness##
brich <- glmer(Richness ~  Basket_Treatment + (1 | Flume), data = Flume_Main_DB_23, family = poisson)
summary(brich)
####Lrt for richness##
brichreduced <- glmer(Richness ~ 1 + (1 | Flume), data = Flume_Main_DB_23, family = poisson)
anova(brich,brichreduced,test = 'LRT')
pairs(emmeans(brich, ~ Basket_Treatment))
###basket level biomass, POM, ER, LMM###
bio = lmer(log(Biomass+0.01) ~ Flume_Treatment + (1|Flume) , data = Flume_Main_DB_23)
summary(bio)
anova(bio, ddf = 'Kenward-Roger') 
lsmeans(OM, pairwise ~ 'Flume_Treatment', adjust = 'none')
####basket level POM lmm###
organicmatter=lmer(log(OM) ~ Basket_Treatment + (1|Flume) , data = Flume_Main_DB_23)
summary(organicmatter)
##basket_resp##
resp=lmer(o2_L_H ~ Basket_Treatment + (1|Flume) , data = Flume_Main_DB_23)
summary(resp)
#### test of differences in standard deviation of density at flume level anova##

dendev <- lm(den_Out_dev~ Treatment, data =dev)
summary(dendev)
denaov<-aov(den_Out_dev~ Treatment, data =dev)
summary(denaov)
dendev
denaov
#### test of differences in standard deviation of biomass at flume level anova##
res.aov <- aov(bio_Out_dev~ Treatment, data =Flume_Standard_deviation_23)
summary(res.aov)
lm.bio <- lm(bio_Out_dev~ Treatment, data =Flume_Standard_deviation_23)
summary(lm.bio)
#### test of differences in standard deviation of richness at flume level anova##
richdev.aov <- aov(rich_out_Dev~ Treatment, data =Flume_Standard_deviation_23)
summary(richdev.aov)
lm.rich <- lm(rich_out_Dev~ Treatment, data =Flume_Standard_deviation_23)
summary(lm.rich)
###organic matter deviation###
om_aov<-aov(OM_out_Dev~Treatment, data=Flume_Standard_deviation_23)
summary(om_aov)
lm_ov<-lm(OM_out_Dev~Treatment, data=Flume_Standard_deviation_23)
summary(lm_ov)
resp_aov<-aov(Resp_out_Dev~Treatment, data=Flume_Standard_deviation_23)
summary(resp_aov)
lm_resp<-lm(Resp_out_Dev~Treatment, data=Flume_Standard_deviation_23)
summary(lm_resp)
##kenaward roger tests###
##flume level##
bio = lmer(log(Biomass+0.01) ~ Flume_Treatment + (1|Flume) , data = Flume_Main_DB_23)
summary(bio)
anova(bio, ddf = 'Kenward-Roger')
##patchlevel_biomass##
bio = lmer(log(Biomass+0.01) ~ Basket_Treatment + (1|Flume) , data = Flume_Main_DB_23)
summary(bio)
anova(bio, ddf = 'Kenward-Roger')
###Flume level POM kenward###
OM = lmer(log(OM) ~ Flume_Treatment + (1|Flume) , data = Flume_Main_DB_23)
summary(OM)
anova(OM, ddf = 'Kenward-Roger')
##bakset level om kenward##
OM = lmer(log(OM) ~ Basket_Treatment + (1|Flume) , data = Flume_Main_DB_23)
summary(OM)
anova(OM, ddf = 'Kenward-Roger')
###Flume level resp kenward roger###
resp=lmer(o2_L_H ~ Flume_Treatment + (1|Flume) , data = Flume_Main_DB_23)
summary(resp)
anova(resp, ddf = 'Kenward-Roger') 
##basket level resp kenward###
resp=lmer(o2_L_H ~ Basket_Treatment + (1|Flume) , data = Flume_Main_DB_23)
summary(resp)
anova(resp, ddf = 'Kenward-Roger') 
###Post-hoc kenward###
##biomass flume##
lsmeans(bio, pairwise ~ 'Flume_Treatment', adjust = 'none')
##POM flume###
lsmeans(resp, pairwise ~ 'Flume_Treatment', adjust = 'none')
###resp flume##
##post hoc basket kenward###
###biomass##
bio = lmer(log(Biomass+0.01) ~ Basket_Treatment + (1|Flume) , data = Flume_Main_DB_23)
summary(bio)
lsmeans(bio, pairwise ~ 'Basket_Treatment', adjust = 'none')
##POM basket posthoc##
OM = lmer(log(OM) ~ Basket_Treatment + (1|Flume) , data = Flume_Main_DB_23)
lsmeans(OM, pairwise ~ 'Basket_Treatment', adjust = 'none')
###Ecosystem resp posthoc basket##
resp=lmer(o2_L_H ~ Basket_Treatment + (1|Flume) , data = Flume_Main_DB_23)
lsmeans(resp, pairwise ~ 'Basket_Treatment', adjust = 'none')
##post-hoc GLM##
##for density##
pairs(emmeans(m1a, ~ Flume_Treatment))
##richness##
brich <- glmer(Richness ~  Flume_Treatment + (1 | Flume), data = Flume_Main_DB_23, family = poisson)
pairs(emmeans(brich, ~ Flume_Treatment))
##post-hoc for basket level GLM density###
pairs(emmeans(bden, ~ Basket_Treatment))
###post-hoc basket level glm richness##
pairs(emmeans(brich, ~ Basket_Treatment))

###GLM for reach level affect of treatment on POM interaction response of density####
m2a <- glmer(density ~  OM*Flume_Treatment + (1 | Flume), data = Flume_Main_DB_23, family = poisson)
summary(m2a)###interaction term is not significant p=0.149, therefore additive model is next appropriate step###

m2reduced<-glmer(density ~  OM+Flume_Treatment + (1 | Flume), data = Flume_Main_DB_23, family = poisson)
summary(m2reduced)
#estimating r2 for model##
##=0.5974269###
r.squaredGLMM(m2reduced)
###LRT to test if interaction model outperforms reduced additive model##
anova(m2a,m2reduced,test = 'LRT')###p-value=0.1508, therefore reduced model is supported here##
###reduce the model further remove treatment effect all together###
M2reduceX2<-glmer(density ~  OM + (1 | Flume), data = Flume_Main_DB_23, family = poisson)
summary(M2reduceX2)
anova(m2reduced,M2reduceX2,test = 'LRT')###signficant p-value confirming to have the additive effect of flume treatment in the model###
###GLM for reach level affect of treatment on POM interaction response of Richness####
mra <- glmer(Richness ~  OM*Flume_Treatment + (1 | Flume), data = Flume_Main_DB_23, family = poisson)
summary(mra)###no significant interaction effect p=0.914###
##therefore lets reduce it to an additive model###
mra_reduced<- glmer(Richness ~  OM+Flume_Treatment + (1 | Flume), data = Flume_Main_DB_23, family = poisson)
summary(mra_reduced)###signficant intercept effect of treatment P=0.00747###
###LRT test between interactive and additive model###
anova(mra,mra_reduced,test = 'LRT')###p=0.914 confirming reduced model is best choice###
##LRT between additive and model with no factor for treatment###
mra_reduced_2<- glmer(Richness ~  OM + (1 | Flume), data = Flume_Main_DB_23, family = poisson)
anova(mra_reduced,mra_reduced_2,test = 'LRT')##significant p, therefore keep the treatment effect in the model##
#estimating r2 for model##
##=0.3123734###
r.squaredGLMM(mra_reduced)
summary(mra_reduced)
#####GLM for patch level affect of basket treatment on POM interaction response of density####
##density##
mda <- glmer(density ~  OM*Basket_Treatment + (1 | Flume), data = Flume_Main_DB_23, family = poisson)
summary(mda)
###no interaction effect###
mda_reduced<-glmer(density ~  OM+Flume_Treatment + (1 | Flume), data = Flume_Main_DB_23, family = poisson)
summary(mda_reduced)
##significant intercept effect with caddisfly compared to control, P<0.0001##
##LRT between interaction and additive model##
anova(mda,mda_reduced, test= 'LRT')##P=0.04083 suggesting some support for the interactive model, but based on output from summary for glmer going with additive###
###LRT for reduced compared to no treamtent effect
mda_reduced_2<-glmer(density ~  OM+ (1 | Flume), data = Flume_Main_DB_23, family = poisson)
anova(mda_reduced,mda_reduced_2, test= 'LRT')###P=<0.0001, importance of treatment##
#####GLM for patch level affect of basket treatment on POM interaction response of richness####
mrab <- glmer(Richness ~  OM*Basket_Treatment + (1 | Flume), data = Flume_Main_DB_23, family = poisson)
summary(mrab)###no interaction effects###
mrab_reduced<- glmer(Richness ~  OM+Basket_Treatment + (1 | Flume), data = Flume_Main_DB_23, family = poisson)
summary(mrab_reduced)##signficant intercept effect of caddisfly treatment compared to the control##
##LRT###basket level richness###
anova(mrab,mrab_reduced, test= 'LRT')###p=0.2334 stick with the additive model###
##LRT to compare additive to model with no treatment factor###
mrab_red_2<-glmer(Richness ~  OM + (1 | Flume), data = Flume_Main_DB_23, family = poisson)
anova(mrab_reduced,mrab_red_2, test= 'LRT')
###importance of basket treatment##P=0.002817###

###LMER test of treatment OM effect###
##Flume level###
##biomass##
bio = lmer(log(Biomass+0.01) ~ OM*Flume_Treatment + (1|Flume) , data = Flume_Main_DB_23)
summary(bio)##no interaction effect###
bio_ad=lmer(log(Biomass+0.01) ~ OM+Flume_Treatment + (1|Flume) , data = Flume_Main_DB_23)
summary(bio_ad)
###respiration###
resp=lmer(o2_L_H ~ OM*Flume_Treatment + (1|Flume) , data =  Flume_Main_DB_23)
summary(resp)
###kenward and posthoc for biomass and resp##
anova(bio_ad, ddf = 'Kenward-Roger') 
anova(resp, ddf = 'Kenward-Roger')
lsmeans(resp, pairwise ~ 'Flume_Treatment', adjust = 'none')

lsmeans(bio_ad, pairwise ~ 'Flume_Treatment', adjust = 'none')
anova(bio_ad, ddf = 'Kenward-Roger') 
bio_ad=lmer(log(Biomass+0.01) ~ OM+Flume_Treatment + (1|Flume) , data = Flume_Main_DB_23)
summary(bio_ad)

bio_ad=lmer(log(Biomass+0.01) ~ OM+Basket_Treatment + (1|Flume) , data = Flume_Main_DB_23)
summary(bio_ad)
resp=lmer(o2_L_H ~ OM+Basket_Treatment * (1|Flume) , data =  Flume_Main_DB_23)
summary(resp)
resp=lmer(o2_L_H ~ OM*F_Treat + (1|Flume) , data = Flume_Main_DB_23)
summary(resp)

###nutrient model fit##
nut.fit = lm(Ammonium~Total_den*Treatment, data=Nutrient_POM_DB_23)
summary(nut.fit)
##

