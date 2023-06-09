#R script used in Bracic, Bohn et al (2021)
#author: Marko Bracic (contribution from Holger Schielzeth)

###################################################################################################
# Model for analyzing repeatability of CJB
###################################################################################################
# This code generates values from the statistical analysis of cognitive judgment bias test used in:
  # In-text values in the manuscript Results section
  # Figure 4: Repeatability of cognitive judgment bias
  # Supplementary Figure 3: Individual choice scores across repeated CJB test
  # Additional analysis 
    # model with all ambiguous cues
    # Repeatability of reference cues
    # random slope models
    # retest effect


### Load packages ###
library(lme4)#mixed models
library(lmerTest)#gives p-values from lme4 models
library(rptR)#extraction of variance components for estimating repeatabilities
library(ggplot2)#plotting
library(ggpubr)#plotting


### Loading the data ###

md<-read.table("CJB_repeatablity.txt", header = TRUE)


### Preparing the data ###

md$cue = factor(md$cue, c("P", "NP", "M", "NN", "N"))
  #setting up all cues as factors (Positive,Near Positive, Middle and Near Negative, and Negative)

md$cueC = relevel(md$cue, ref="M")
  # M cue set as reference cue (model values will not change)


### Models for repeatability calculation of the choice score ####

#Separate model for each ambiguous cue with random intercepts for individuals

#Near Positive cue
summary(lmer(score ~ test + (1|ind), data=subset(md, md$cue=="NP")))
  #checking the model

set.seed(5639);
R_NP = rptGaussian(score ~ test + (1|ind), grname="ind", 
                   data=subset(md, md$cue=="NP"), nboot=1000)
R_NP 
  #it provides repeatability score, p-values, and confidence intervals (if nboot is not = 0)
  #set.seed allows reproducible results: to get the same CI each time (if not used CI will slightly change for every run)


#Middle cue
summary(lmer(score ~ test + (1|ind), data=subset(md, md$cue=="M")))

set.seed(8546);
R_M = rptGaussian(score ~ test + (1|ind), grname="ind", 
                  data=subset(md, md$cue=="M"), nboot=1000)


#Near Neagative cue
summary(lmer(score ~ test + (1|ind), grname="ind", 
             data=subset(md, md$cue=="NN")))

set.seed(5938);
R_NN = rptGaussian(score ~ test + (1|ind), grname="ind", 
                   data=subset(md, md$cue=="NN"), nboot=1000)


#Table with R values, confidence intervals,and p-values (reported in manuscript in-text)

R_cjb = data.frame(factor(c("Near Positive", "Middle", "Near Negative"), 
                          levels=c("Near Positive", "Middle", "Near Negative")))
colnames(R_cjb)[1] = "Cue"
R_cjb$Repeatablility = unlist(lapply(c(R_NP$R$ind, R_M$R$ind, R_NN$R$ind), round,2))
R_cjb$CI_min = unlist(lapply(c(R_NP$CI_emp$`2.5%`,R_M$CI_emp$`2.5%`,R_NN$CI_emp$`2.5%`), round,2))
R_cjb$CI_max = unlist(lapply(c(R_NP$CI_emp$`97.5%`,R_M$CI_emp$`97.5%`,R_NN$CI_emp$`97.5%`), round,2))
R_cjb$p = unlist(lapply(c(R_NP$P$LRT_P,R_M$P$LRT_P,R_NN$P$LRT_P), format.pval,1, 0.001,0))


##########################################################################################
#                         Figure Repeatability of cognitive judgment bias             
##########################################################################################

#Preparing the plot from data obtained from above models

R_cjb = data.frame(factor(c("Near Positive", "Middle", "Near Negative"), 
                          levels=c("Near Positive", "Middle", "Near Negative")))
colnames(R_cjb)[1] = "Cue"
R_cjb$Repeatablility = unlist(c(R_NP$R$ind, R_M$R$ind, R_NN$R$ind))
R_cjb$CI_min = unlist(c(R_NP$CI_emp$`2.5%`,R_M$CI_emp$`2.5%`,R_NN$CI_emp$`2.5%`))
R_cjb$CI_max = unlist(c(R_NP$CI_emp$`97.5%`,R_M$CI_emp$`97.5%`,R_NN$CI_emp$`97.5%`))
  #same table as above but values not rounded

#plot of repeatability estimate and confidence intervals (for ambiguous cues only)
F4 = ggdotplot(R_cjb, x ="Cue", 
          y= "Repeatablility", 
          size =0, 
          add = "point", 
          add.params = list(size=2)) +
  geom_errorbar(aes(ymin=CI_min, ymax=CI_max, width = 0.1)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1)) + 
  annotate("text", x = 1, y = (R_NP$CI_emp$`97.5%`+ 0.031), 
           label = paste0("p = ",format.pval(R_NP$P$LRT_P, 1, 0.001,0)), size =3.5) +
  annotate("text", x = 2, y = (R_M$CI_emp$`97.5%`+ 0.031), 
           label = paste0("p = ",format.pval(R_M$P$LRT_P, 1, 0.001,0)), size =3.5) +
  annotate("text", x = 3, y = (R_NN$CI_emp$`97.5%`+ 0.031), 
           label = "p > 0.99", size =3.5) +
  scale_x_discrete(expand = expansion(mult = c(0, 0))) +
  xlab("\\n")+#empty space under x axis and removes axis name
  theme(axis.ticks.x=element_blank(), axis.line.x=element_blank(),
        axis.title=element_text(size=18), text=element_text(size=14)) +
  coord_cartesian(clip = 'off')#points are not hidden by axis

#tiff("Figure 4.tiff", res=350,width = 12.9, height = 12.8, units = "cm")
#F4
#dev.off()
#run the line above to save the figure


##########################################################################################
#             Additional analysis reported in Supplementary Data          
##########################################################################################

# Supp. Figure 3.Individual variation in choice score 

SF3= ggplot(md, aes(x=test, y=score)) +
  geom_line(aes(colour=ind)) +  
  scale_colour_grey(guide = "none") +
  ylab("Choice score") +
  xlab("CJB test order") +
  theme_classic() +
  theme(axis.title=element_text(size=20), axis.text  = element_text(size=14), 
        strip.text = element_text(size = 15), legend.position = "none") +
  ylim(-1,1) +
  facet_wrap(vars(cue), ncol=5, 
             labeller = labeller(cue = c("P" = "Positive", "NP" = "Near Positive",
                                         "M" = "Middle", "NN" = "Near Negative", "N" = "Negative")))
#ggsave("Supplementary Figure 3.tiff", width = 26, height = 15, units = "cm", dpi=350)



##########################################################################################
#             Additional analysis (not reported in the manuscript )             
##########################################################################################

### Preparing the data ###

md$cueNum = as.numeric(md$cue)
#cue as continues variable 

md$cueC = relevel(md$cue, ref="M")
# M cue set as reference cue (model values will not change)

md$cueselect = md$cue=="NP" | md$cue=="M" | md$cue=="NN"
#easy way to select which cues we want
#only ambiguous cues are used: Near Positive, Middle and Near Negative
#compares if value is one of the selected cues and gives back true or falls


### Setting up statistical model ###

#Centering variables -> increases model interpretation

md$testC = md$test - 2.5
#four repeated tests centered to mean
#mean = 2.5 = ((1+2+3+4)/4 tests)

md$cueNumC = md$cueNum - 3
# five cues centered to mean


### Additional Models for repeatability of the choice score during model selection ####

# Checking results for all ambiguous cues and reference cues

# For all ambiguous cues combined (cue*test interaction in the model)
R_AMB = rptGaussian(score ~ cue*testC + (1|ind), grname="ind", 
                    data=subset(md, md$cueselect), nboot=1000)
# also repeatable: lower then for NP and M but higher than for NN  --> support that CJB is repeatable

# Positive cue
R_P = rptGaussian(score ~ test + (1|ind), grname="ind", 
                  data=subset(md, md$cue=="P"), nboot=1000)
#There is an outlier in P cue (individual 23CLR) which inflates repeatability, remove?
R_P = rptGaussian(score ~ test + (1|ind), grname="ind", 
                  data=subset(md, md$cue=="P" & md$ind!="23CLR"), nboot=1000)
#23CLR removed -> much lower repeatability 

# Negative cue
R_N = rptGaussian(score ~ test + (1|ind), grname="ind", 
                  data=subset(md, md$cue=="N"), nboot=1000)


# Checking for individual specific response to cue and repeated tests ####

# Random slope model [individual-specific responses to cue] -> not needed
mod = lmer(score ~ cueNumC*test + (cueNumC|ind), 
           data=subset(md, md$cue=="NP" | md$cue=="M" | md$cue=="NN"))
modred = lmer(score ~ cueNumC*test + (1|ind), 
              data=subset(md, md$cue=="NP" | md$cue=="M" | md$cue=="NN"))
anova(mod, modred)
#comparing too models -> no sig. difference

# Random slope model with all cues [individual-specific time trends]-> not needed
mod = lmer(score ~ cueNumC*test + (test|ind), 
           data=subset(md, md$cue=="NP" | md$cue=="M" | md$cue=="NN"))
modred = lmer(score ~ cueNumC*test + (1|ind), 
              data=subset(md, md$cue=="NP" | md$cue=="M" | md$cue=="NN"))
anova(mod, modred)

# Random slope models[individual-specific time trends] -> not needed
mod = lmer(score ~ test + (test|ind), data=subset(md, md$cue=="NP"))
modred = lmer(score ~ test + (1|ind), data=subset(md, md$cue=="NP"))
anova(mod, modred)

mod = lmer(score ~ test + (test|ind), data=subset(md, md$cue=="M"))
modred = lmer(score ~ test + (1|ind), data=subset(md, md$cue=="M"))
anova(mod, modred)

mod = lmer(score ~ test + (test|ind), data=subset(md, md$cue=="NN"))
modred = lmer(score ~ test + (1|ind), data=subset(md, md$cue=="NN"))
anova(mod, modred)

#INDIVIDUALS DO NOT HAVE SIGNIFICANTLY DIFFERENT SLOPES IN EACH TEST OR CUE
# no need for random slope test|ind or cueNumC|ind in the model



### Centered linear mixed-effect model for retest effect on the choice score  ###

# Random intercept model with cue*test interaction
mod_RetestEffect = lmer(score ~ cueC*testC + (1|ind), 
                        data=subset(md, md$cue=="NP" | md$cue=="M" | md$cue=="NN"))

summary(mod_RetestEffect)  
# mice become more pessimistic over tests in middle cue, but more optimistic for the NP cue
  # significant but small effect