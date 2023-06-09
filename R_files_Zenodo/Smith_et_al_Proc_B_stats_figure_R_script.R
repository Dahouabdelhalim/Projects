R.Version()
require(dplyr)
require(tidyr)
require(lme4)
require(plyr)
require(readxl)
require(ggplot2)
require(multcomp)
require(lmerTest)
require(merTools)
require(boot)

#load the data and compress it to one row of responsiveness and learning per individual 
base <- "C:/Users/Andres/Box/Manuscripts/2019 Smith et al/drafts/ProcB/analysis/"

data_long <- read.csv(paste0(base, "01_data/02_for_submission/DS_dataset_long.csv"))
data_long$treatment<- factor(data_long$treatment, levels = c("control","pre-eclosion","post-eclosion","continual"))

#columns with .C are "centred" values. Not needed as we can do it in the analysis.
datasetIndividual <- unique(data_long[, !colnames(data_long) %in% 
                                        c("trial","PER_outcome","PER_response","cumsum"#,"temp_PER"
                                          #"Mcal.left.C","Lcal.left.C","Mcal.right.C","Lcal.right.C",
                                          #"Calyces.Total.C","lobes.left.C","lobes.right.C","MB.left.C", "MB.right.C"
                                          )])
#theme for plotting: 
plot.theme <- theme_classic(base_size = 16)+theme(strip.background = element_blank(),axis.title.x = element_blank())
ggplot(datasetIndividual, aes(x = treatment, y = ITD))+geom_violin()+geom_jitter(aes(color = colony))+plot.theme

#some summaries:
#number of bees
a = ddply(datasetIndividual, c("treatment", "age", "colony"), summarise,
      totalNo    = length(id)
)
#number of individuals that respond:
b = ddply(datasetIndividual, c("treatment", "age", "colony", "show_PER_in_response_test"), summarise,
      Responders    = length(id)
)

#number of individuals that learn:
c = ddply(datasetIndividual, c("treatment", "age", "colony", "learner"), summarise,
      Learners    = length(id)
)

#number of individuals that were scanned:
d = ddply(datasetIndividual[datasetIndividual$brain_vol_taken=="y",], c("treatment", "age", "colony", "learner"), summarise,
          Scanned    = length(id)
)

respond = merge(a,b, by = c("treatment", "age", "colony")); write.csv(respond, paste0(base,"/02_output/01_responsive.csv"), row.names = FALSE)
learn   = merge(a,c, by = c("treatment", "age", "colony")); write.csv(learn,   paste0(base,"/02_output/02_learner.csv"), row.names = FALSE)
scaned  = merge(a,d, by = c("treatment", "age", "colony")); write.csv(scaned,  paste0(base,"/02_output/03_scanned.csv"), row.names = FALSE)

#Responsiveness:
#same analysis in the master dataset
m2 <- glmer(show_PER_in_response_test ~ treatment * age + ITD  + (1|colony), data = datasetIndividual,   family = binomial)

#the mixed effects model is not required here as the fit is singular and the model 02_outputs are identical:
#my analysis
m3 <- glm  (show_PER_in_response_test ~ treatment * age + ITD              , data = datasetIndividual,   family = binomial)

anova(m2,m3)
summary(m3)

#irrespective of which test we use, if you have a look at the anova 02_output for this model you see that only treatment and 
#age have p<0.05 or  F values of >2 indicating a significant effect. 
#therefore remove the interaction from the model.
# anova(m2);
anova(m3, test = "Chisq")

#is the model structure we have specified the best fitting model? i.e are there interactive effects between age and treatment?

m4 <- glm  (show_PER_in_response_test ~ treatment + age + ITD              , data = datasetIndividual,   family = binomial)

#there does not appear to be an interaction between age and treatment
testInteraction <-
anova(m3,m4,test = "Chisq")
write.csv(as.data.frame(testInteraction), paste0(base,"/02_output/04_Responsiveness_ChiSqModelWithAndWithoutInteraction_GLMBinomial.csv") )

anova(m3,test = "Chisq")
anova(m4,test = "Chisq")

summary(m4); summary( glht(m4, mcp(treatment="Tukey")) )

write.csv(as.data.frame(summary(m4)$coef), paste0(base,"/02_output/04_Responsiveness_ModelSummary_GLMBinomial.csv") )

#plot the data:
newdata <- expand.grid(
  treatment = levels(datasetIndividual$treatment) ,
  age       = levels(datasetIndividual$age),
  ITD       = median(datasetIndividual$ITD)
)
newdata$fit    = predict(m4, newdata, se.fit = TRUE, interval = "predict", type = "response")$fit
newdata$se.fit = predict(m4, newdata, se.fit = TRUE, interval = "predict", type = "response")$se.fit

write.csv(newdata, paste0(base,"/02_output/05_AverageResponsiveness.csv") )

#get the raw proportion values from the data
a = datasetIndividual[!is.na(datasetIndividual$show_PER_in_response_test),]
b = ddply(a, c("treatment", "age","show_PER_in_response_test"), summarise, N    = length(id))
c = ddply(a, c("treatment", "age"), summarise, N    = length(id))
d = merge(b, c, by = c("treatment", "age"))
d$prop <- d$N.x/d$N.y
resp.e = d[d$show_PER_in_response_test =="Yes",]

p1 = ggplot(newdata, aes(x = age, y = fit ))+geom_point(size = 4)+
  geom_errorbar(data = newdata, aes(ymin=newdata$fit-(newdata$se.fit*1.96), ymax=newdata$fit+(newdata$se.fit*1.96)))+
  geom_point(data = resp.e, aes(x = age, y = prop),size = 4, color = "red", shape = 18)+
  facet_grid(~treatment)+plot.theme+
  plot.theme+ylab("Proportion responding");p1

#The summary ANOVA table tells us that at day three there is no difference in the proportion of individuals responding between control and 
#pre eclosion individuals, but there is a reduction in the proportion responding in the post-eclosion and continual treatments.

#there is an effect of age, older bees are more likely to respond.

#There is no interactive effect of treatment and age, indicating that the effect of treatment is the same in 3 and 12 day old individuals. 


#learners:

#the mixed effects model is not required here as the fit is singular and the model 02_outputs are identical:
m2 <- glmer(learner ~ treatment * age + ITD  + (1|colony), data = datasetIndividual, family = binomial(), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))
m3 <- glm(learner ~ treatment * age + ITD, data = datasetIndividual, family = binomial )
#summary(m2);
anova(m2,m3)
summary(m3)

#irrespective of which test we use, if you have a look at the anova 02_output for this model you see that only treatment and 
#ITD have p<0.05 or  F values of >2 indicating a significant effect. 
#therefore, interaction removed from the model.
anova(m3, test = "Chisq")

#is the model structure we've specified the best fitting model? i.e are there interactive effects between age and treatment?
m4 <- glm(learner ~ treatment + age + ITD, data = datasetIndividual, family = binomial ); summary(m4)

#there does not appear to be an interaction between age and treatment
anova(m3,m4,test = "Chisq")

testInteraction <-
anova(m3,m4,test = "Chisq")
write.csv(as.data.frame(testInteraction), paste0(base,"/02_output/07_Learning_ChiSqModelWithAndWithoutInteraction_GLMBinomial.csv") )

anova(m4,test = "Chisq")
summary(m4)

summary( glht(m4, mcp(treatment="Tukey")) )

write.csv(as.data.frame(summary(m4)$coef), paste0(base,"/02_output/06_Learning_ModelSummary_GLMBinomial.csv") )

#plot the data from the best model:

newdata3 <- expand.grid(
  treatment = levels(datasetIndividual$treatment) ,
  age       = levels(datasetIndividual$age),
  ITD       = median(datasetIndividual$ITD)
)
newdata3$fit    = predict(m4, newdata3, se.fit = TRUE, interval = "predict", type = "response")$fit
newdata3$se.fit = predict(m4, newdata3, se.fit = TRUE, interval = "predict", type = "response")$se.fit

write.csv(newdata, paste0(base,"/02_output/07_AverageLearning.csv") )

#get the raw proportion values from the data
a = datasetIndividual[!is.na(datasetIndividual$learner),]
b = ddply(a, c("treatment", "age","learner"), summarise, N    = length(id))
c = ddply(a, c("treatment", "age"), summarise, N    = length(id))
d = merge(b, c, by = c("treatment", "age"))
d$prop <- d$N.x/d$N.y
learn.e = d[d$learner =="Yes",]


p2 = 
  ggplot(newdata3, aes(x = age, y = fit))+geom_point(size = 4)+
  geom_errorbar(data = newdata3, aes(ymin=newdata3$fit-(newdata3$se.fit*1.96), ymax=newdata3$fit+(newdata3$se.fit*1.96)))+
  geom_point(data = learn.e, aes(x = age, y = prop), color = "red", size = 4, shape = 18)+
  facet_grid(~treatment)+plot.theme+ylab("Proportion learning");p2

require(cowplot)

plot_grid(p1,p2,labels = "AUTO", ncol = 1, align = 'v')+
ggsave( paste0(base,"/02_output/08_responsivenessAandlearningB.png"), height = 6, width = 9) 

#this tells us that of the individuals that are responsive, at day three each treatment group has a lower proportion of individuals
#that are able to learn relative to the control. There is no effect of age on this relationship

#additionally there is an effect of body size as ITD is positively associated with proportion that learnt.

#clean out the environment
rem <- as.character(c(ls()))
rem = rem[!rem %in% c("base", "datasetIndividual", "Volume_and_Learning", "plot.theme","data_long")]
rm(list =rem)

#the next analysis is the brain size analysis

brain_vol <- datasetIndividual[datasetIndividual$brain_vol_taken=="y",] #contains the untransformed brain volumes fot other structures


#different brain measurements:
# "Mcal.left"                                "Lcal.left"                                "Mcal.right"                              
# "Lcal.right"                               "Calyces.Total"                            "lobes.left"                              
# "lobes.right"                              "MB.left"                                  "MB.right"                                
# "Calyces.Total.C" should be just the two calyces and not the lobes ("stem") which is also part of the mushroom body.

# note that calyces was spelt wrong in the datasheet, and so cayleces was used in the model calls

#total mushroom body size (must have measurement for all three components on the left and the right)
totalMB <- brain_vol[  !is.na(brain_vol$Mcal.left) & !is.na(brain_vol$Lcal.left) & !is.na(brain_vol$lobes.left) &
            !is.na(brain_vol$Mcal.right) & !is.na(brain_vol$Lcal.right) & !is.na(brain_vol$lobes.right),]


#lobe size (mushroom body component - not thought to be part of learning)
totalLO <- totalMB[  !is.na(totalMB$lobes.left) &
                         !is.na(totalMB$lobes.right),]


#total Calyces size (mushroom body component - thought to be part of learning)
totalCayleces <- totalMB[  
  !is.na(totalMB$Mcal.left) & !is.na(totalMB$Lcal.left) & 
    !is.na(totalMB$Mcal.right) & !is.na(totalMB$Lcal.right),]


#add the volumes for the whole structures
totalMB$TotalMB = totalMB$Mcal.left + totalMB$Mcal.right+
                  totalMB$Lcal.left + totalMB$Lcal.right+
                  totalMB$lobes.left+ totalMB$lobes.right

totalLO$lobes   = totalLO$lobes.left + totalLO$lobes.right

totalCayleces$TotalCY = totalCayleces$Mcal.left + totalCayleces$Mcal.right+
                        totalCayleces$Lcal.left + totalCayleces$Lcal.right

#others (central body, antennal lobe, lobula & Medula ):
totalMB$antennal_lobe <- totalMB$left_antennal_lobe + totalMB$right_antennal_lobe
totalMB$lobula        <- totalMB$left_lobula+totalMB$right_lobula
totalMB$medulla       <- totalMB$left_medulla+totalMB$right_medulla

#scale by body size:
totalMB$scaledMB       <- totalMB$TotalMB/totalMB$ITD
totalMB$scaledAL       <- totalMB$antennal_lobe/totalMB$ITD
totalMB$scaledCB       <- totalMB$central_body/totalMB$ITD
totalMB$scaledLOBULA   <- totalMB$lobula/totalMB$ITD
totalMB$scaledMED      <- totalMB$medulla/totalMB$ITD
totalLO$scaledLO       <- totalLO$lobes  /totalLO$ITD
totalCayleces$scaledCY <- totalCayleces$TotalCY/totalCayleces$ITD

#is there any relationship between the treatment/age of an individual and the size of i) the full MB or ii) its component parts
#no interaction but significant treatment effect. Mushroom bodies are smaller in individuals exposed to pesticides as Adults.

#Total MB size.
MB1  <- lmer(scaledMB   ~ treatment * age + ITD + (1|colony), data = totalMB); anova(MB1);  summary(MB1)  #linear mixed-effects model
MB2  <- lmer(scaledMB   ~ treatment + age + ITD + (1|colony), data = totalMB); anova(MB2);  summary(MB2)  #linear mixed-effects model
anova(MB1,MB2,test = "Chisq")
summary(MB2)
summary( glht(MB2, mcp(treatment="Tukey")) )

#No treatment by age interaction. 
#Reduction in MBs of workers exposed as adults.
#significant reduction in relative MB size as ITD increases

write.csv(as.data.frame(summary(MB2)$coef), paste0(base,"/02_output/09_TotalMB_BrainSize_ModelSummary_LMER.csv") )

#just the lobes.
LO1  <- lmer( scaledLO ~ treatment * age + ITD + (1|colony), data = totalLO); anova(LO1);  summary(LO1)  #linear mixed-effects model
LO2  <- lmer( scaledLO ~ treatment + age + ITD + (1|colony), data = totalLO); anova(LO2);  summary(LO2)  #linear mixed-effects model
anova(LO1,LO2,test = "Chisq"); 
summary(LO2)
summary( glht(LO2, mcp(treatment="Tukey")) )

#No treatment or age effect.
#significant reduction in relative MB size as ITD increases

write.csv(as.data.frame(summary(LO2)$coef), paste0(base,"/02_output/10_Lobes_BrainSize_ModelSummary_LMER.csv") )

#just the Calyces:
MBC1  <- lmer( scaledCY   ~ treatment * age + ITD + (1|colony), data = totalCayleces); anova(MBC1);  summary(MBC1)  #linear mixed-effects model
MBC2  <- lmer( scaledCY   ~ treatment + age + ITD + (1|colony), data = totalCayleces); anova(MBC2);  summary(MBC2)  #linear mixed-effects model
anova(MB1,MB2,test = "Chisq"); 
summary(MBC2)
summary( glht(MBC2, mcp(treatment="Tukey")) )

#no treatment by age interaction. 
#Significant treatment effects - each pesticice treatment has smaller calyces than the control
#no effect of age
#significant reduction in relative MB size as ITD increases

write.csv(as.data.frame(summary(MBC2)$coef), paste0(base,"/02_output/11_Cayleces_BrainSize_ModelSummary_LMER.csv") )

#plot the data:
newdata <- data.frame(expand.grid(
  treatment = levels(totalMB$treatment) ,
  age       = levels(totalMB$age),
  ITD       = c(median(totalMB$ITD)) 
))

#Bootstrap the CIs fo the model predictions:
#get the bootstrapped CI for the model parameters
#example only, not required for CI around mean.

predFun       <- function(x) predict(x,newdata=newdata,re.form=NA)
newdata$MB    <- predFun(MB2)
bb            <- bootMer(MB2,FUN=predFun,nsim=200); bb_ci <- as.data.frame(t(apply(bb$t,2,quantile,c(0.025,0.975)))); names(bb_ci)  <- c("MB_lwr","MB_upr")
newdata       <- cbind(newdata,bb_ci)
newdata$LO    <- predFun(LO2)
bb            <- bootMer(LO2,FUN=predFun,nsim=200); bb_ci <- as.data.frame(t(apply(bb$t,2,quantile,c(0.025,0.975)))); names(bb_ci)  <- c("LO_lwr","LO_upr")
newdata       <- cbind(newdata,bb_ci)
newdata$CY    <- predFun(MBC2)
bb            <- bootMer(MBC2,FUN=predFun,nsim=200); bb_ci <- as.data.frame(t(apply(bb$t,2,quantile,c(0.025,0.975)))); names(bb_ci)  <- c("CY_lwr","CY_upr")
newdata       <- cbind(newdata,bb_ci)

write.csv(newdata, paste0(base,"/02_output/12_AverageBrainSize.csv") )


#calculate the raw means for mushroom bodies:
MBrawMean = ddply(totalMB, c("treatment", "age"), summarise, N = length(id), mean = mean(TotalMB/ITD),sd   = sd(TotalMB/ITD), ci   = sd / sqrt(N)*1.96)
#calculate the raw means for the lobes bodies:
LOrawMean = ddply(totalLO, c("treatment", "age"), summarise, N = length(id), mean = mean(lobes/ITD),sd   = sd(lobes/ITD), ci   = sd / sqrt(N)*1.96)
#calculate the raw means for Cayleces:
CYrawMean = ddply(totalCayleces, c("treatment", "age"), summarise, N = length(id), mean = mean(TotalCY/ITD),sd   = sd(TotalCY/ITD), ci   = sd / sqrt(N)*1.96)


MB = ggplot()+
  geom_point(   data = newdata, aes(x = age, y = MB), color = "black", size = 4)+
  geom_errorbar(data = newdata, aes(x = age, ymin=MB_lwr, ymax=MB_upr), color = "black")+  
  geom_point(   data = MBrawMean, aes(x = age, y=mean), color = "red", shape = 18, size = 4)+
  facet_grid(~treatment)+plot.theme+labs(x = "", y = expression(paste("Mushroom body  ", mm^3) ) ) 

LO = ggplot()+
  geom_point(   data = newdata, aes(x = age, y = LO), color = "black", size = 4)+
  geom_errorbar(data = newdata, aes(x = age, ymin=LO_lwr, ymax=LO_upr), color = "black")+  
  geom_point(   data = LOrawMean, aes(x = age, y=mean), color = "red", shape = 18, size = 4)+
  facet_grid(~treatment)+plot.theme+labs(x = "", y = expression(paste("Lobe  ", mm^3)) ) 

CY = ggplot()+
  geom_point(   data = newdata, aes(x = age, y = CY), color = "black", size = 4)+
  geom_errorbar(data = newdata, aes(x = age, ymin=CY_lwr, ymax=CY_upr), color = "black")+  
  geom_point(   data = CYrawMean, aes(x = age, y=mean), color = "red", shape = 18, size = 4)+
  facet_grid(~treatment)+plot.theme+labs(x = "", y = expression(paste("Calyces  ", mm^3)) ) 

plot_grid(MB,LO,CY, labels = "AUTO", ncol = 1, align = 'v')+ggsave(paste0(base,"/02_output/13_AverageMBSize.png"), height = 9, width = 9)


#test each brain part separately for effects of pesticide treatment (central body, antennal lobe, lobula & Medula ). 
#ignore the fact that there are slightly different individuals within each sample - is described in methods/text
#do an anova and extract model predictions

CB1     <- lmer((central_body/ITD)       ~ treatment + age + ITD + (1|colony), data = datasetIndividual); anova(CB1)     ; summary(CB1) #linear mixed-effects model
AL1     <- lmer((antLobe/ITD)       ~ treatment + age + ITD + (1|colony), data = datasetIndividual); anova(AL1)     ; summary(AL1) #linear mixed-effects model
LOBULA1 <- lmer((lobula/ITD)   ~ treatment + age + ITD + (1|colony), data = datasetIndividual); anova(LOBULA1) ; summary(LOBULA1) #linear mixed-effects model
MED1    <- lmer((medulla /ITD)     ~ treatment + age + ITD + (1|colony), data = datasetIndividual); anova(MED1)    ; summary(MED1) #linear mixed-effects model


plot(CB1)
plot(AL1)
plot(LOBULA1)
plot(MED1)

write.csv(as.data.frame(summary(CB1)$coef), paste0(base,"/02_output/14_CentralBodyBrainSize_ModelSummary_LMER.csv") )
write.csv(as.data.frame(summary(AL1)$coef), paste0(base,"/02_output/15_AntennalLobeBrainSize_ModelSummary_LMER.csv") )
write.csv(as.data.frame(summary(LOBULA1)$coef), paste0(base,"/02_output/16_LobulaBrainSize_ModelSummary_LMER.csv") )
write.csv(as.data.frame(summary(MED1)$coef), paste0(base,"/02_output/17_MedullaBodyBrainSize_ModelSummary_LMER.csv") )

structure_size <- ddply(datasetIndividual, c("treatment", "age"), summarise,
                        MB.N    = length(totalMB[!is.na(totalMB)]),
                        MB.mean = mean(totalMB, na.rm=TRUE),
                        MB.sd   = sd(totalMB, na.rm=TRUE),
                        MB.se   = MB.sd/sqrt(MB.N),
                        
                        CA.N    = length(calyces[!is.na(calyces)]),
                        CA.mean = mean(calyces, na.rm=TRUE),
                        CA.sd   = sd(calyces, na.rm=TRUE),
                        CA.se   = CA.sd/sqrt(CA.N),
                        
                        LO.N    = length(lobe[!is.na(lobe)]),
                        LO.mean = mean(lobe, na.rm=TRUE),
                        LO.sd   = sd(lobe, na.rm=TRUE),
                        LO.se   = LO.sd/sqrt(LO.N),
                        
                        CB.N    = length(central_body[!is.na(central_body)]),
                        CB.mean = mean(central_body, na.rm=TRUE),
                        CB.sd   = sd(central_body, na.rm=TRUE),
                        CB.se   = CB.sd/sqrt(CB.N),
                        #
                        AN.N    = length(antLobe[!is.na(antLobe)]),
                        AN.mean = mean(antLobe, na.rm=TRUE),
                        AN.sd   = sd(antLobe, na.rm=TRUE),
                        AN.se   = AN.sd/sqrt(AN.N),
                        #
                        LOBULA.N    = length(lobula[!is.na(lobula)]),
                        LOBULA.mean = mean(lobula, na.rm=TRUE),
                        LOBULA.sd   = sd(lobula, na.rm=TRUE),
                        LOBULA.se   = LOBULA.sd/sqrt(LOBULA.N),
                        #
                        MED.N    = length(medulla[!is.na(medulla)]),
                        MED.mean = mean(medulla, na.rm=TRUE),
                        MED.sd   = sd(medulla, na.rm=TRUE),
                        MED.se   = MED.sd/sqrt(MED.N)

)



write.csv(structure_size, paste0(base,"/02_output/18_RawBrainSizeSummary.csv") )


#clean out the environment
rem <- as.character(c(ls()))
rem = rem[!rem %in% c("base", "datasetIndividual", "Volume_and_Learning", "plot.theme","data_long","totalMB")]
rm(list =rem)


#Learning score by trial:

trial<- data_long[data_long$outcome =="Successfully completed PER test" & data_long$sum_score >0 ,]
summary(datasetIndividual[datasetIndividual$outcome =="Successfully completed PER test" & 
                            datasetIndividual$sum_score >0 ,]$treat)
trial$trial<- as.numeric(as.character(gsub("t","",trial$trial)))

trial$fail <- 10-trial$cumsum
trial$prop <- trial$cumsum/trial$fail

#the number of individuals in each treatment group is too small to do much with.
#aggregate the treatment and age groups to consist of exposed and control individials
#add a column to indicate control and exposed individuals (irrespective of treatment & age)
trial$agg_treat <-as.factor(ifelse(trial$treatment =="control", "control", "treatment"))

summary(trial$agg_treat)

#first, lets calculate the raw proportions to have an idea of what the data is doing

prop_respond <- ddply(trial, c("agg_treat","trial"), summarise,
                      N       = length(id),        #total tested
                      respond = sum(PER_response) #responded to PER
)

#adding trial and identifying the data as repeated measures causes the model to fail to converge
#adding colony as a nested random effect causes the model to fail to converge
#we can add individual as a random factor

trial$PER_response

# m0 = lmer(PER_response~agg_treat+trial+(1|id), data = trial); anova(m0) 
# m1 = lmer(PER_response~agg_treat+trial+ITD+(1|id), data = trial); anova(m1); anova(m0,m1)     #ITD doesnt improve the model fit
# m2 = lmer(PER_response~agg_treat*trial+(1|id), data = trial); anova(m2); anova(m1,m2)         #model better w/o the interaction
# m3 = lmer(PER_response~agg_treat+poly(trial,2)+(1|id), data = trial); anova(m3); anova(m1,m3) #model better with 2nd order poly
# m4 = lmer(PER_response~agg_treat+poly(trial,3)+(1|id), data = trial); anova(m4); anova(m3,m4) #model better w/o 3rd order poly
# m5 = lmer(PER_response~agg_treat*poly(trial,2)+(1|id), data = trial); anova(m5); anova(m4,m5) #model better w/o interaction between treat & 2rd order poly

#this needs to be a binomial response not gaussian:

m0 = glmer(PER_response~agg_treat+trial+(1|id),        family = "binomial", data = trial); anova(m0); summary(m0)
m0b = glm(PER_response~agg_treat+trial,        family = "binomial", data = trial); anova(m0,m0b, test = "Chisq") #random effects improves the fit
m1 = glmer(PER_response~agg_treat+trial+ITD+(1|id),    family = "binomial", data = trial); anova(m1); anova(m0,m1) #ITD doesnt improve the model fit
m2 = glmer(PER_response~agg_treat*trial+(1|id),        family = "binomial", data = trial); anova(m2); anova(m0,m2) #model better w/o the interaction
m3 = glmer(PER_response~agg_treat+poly(trial,2)+(1|id),family = "binomial", data = trial); anova(m3); anova(m0,m3) #model better with 2nd order poly
m4 = glmer(PER_response~agg_treat+poly(trial,3)+(1|id),family = "binomial", data = trial); anova(m4); anova(m3,m4) #model better w/o 3rd order poly
m5 = glmer(PER_response~agg_treat*poly(trial,2)+(1|id),family = "binomial", data = trial); anova(m5); anova(m3,m5) #model better w/o interaction between treat & 2rd order poly

summary(m3)

#plot the data:
newdata <- data.frame(expand.grid(
  agg_treat = levels(trial$agg_treat),
  trial     = unique(trial$trial)
))

predFun            <- function(x) predict(x,newdata=newdata,re.form=NA, type = "response")
newdata$respond    <- predFun(m3)
bb                 <- bootMer(m3,FUN=predFun,nsim=200); bb_ci <- as.data.frame(t(apply(bb$t,2,quantile,c(0.025,0.975)))); names(bb_ci)  <- c("resp_lwr","resp_upr")
newdata            <- cbind(newdata,bb_ci)

#
ggplot()+
  geom_point(data = prop_respond, aes(x = trial, y = respond/N, shape = agg_treat, color = agg_treat), size = 2)+ylim(0,1)+
  #geom_smooth(method = lm, formula = y ~ x + I(x^2),se = FALSE, aes(linetype=agg_treat))+
  geom_ribbon(data = newdata, aes(x = trial, ymin = resp_lwr, ymax = resp_upr, linetype=agg_treat, fill = agg_treat ),alpha = 0.05) + 
  geom_line(  data = newdata, aes(x = trial, y = respond, group = agg_treat, linetype=agg_treat, color = agg_treat) )+
  scale_x_continuous(breaks=seq(0,10,1))+theme_classic(base_size = 16)+theme(legend.title = element_blank())+
  labs(x = "Trial number", y = "Proportion exhibiting \\n learnt response" )+
  ggsave(paste0(base,"/02_output/19_responsiveness_by_treatment_or_control.png"), height = 6, width = 6)


write.csv(as.data.frame(summary(m3)$coef), paste0(base,"/02_output/20_responsiveness_by_treatment.csv") )

#the model is picking up the general increase in responsiveness observed as the number of trials increases
#the model is picking up that this increase is non linear
#the model shows pesticide treated individuals have a lower proportion of PER responses. 

#the correlation between between Calyces size & learning score:
PER<- datasetIndividual[datasetIndividual$outcome == "Successfully completed PER test" & 
                        datasetIndividual$brain_vol_taken =="y" & 
                       !is.na(datasetIndividual$calyces),] 

PER$agg_treat <-as.factor(ifelse(PER$treatment =="control", "control", "treatment"))
PER$scaledCY<- PER$calyces/PER$ITD

#however this analysis is not great as it is bounded by 0 and 10 so we should be using some kind of proportional model
#at the moment the model is predicting learning scores of >15 for the control treatment which is impossible.
#do a binomial with the cumulative score as successes and 10-cumscore

y = cbind(PER$sum_score, (10-PER$sum_score))

m0 = glm(y~scaledCY+agg_treat,          family = "binomial", data = PER)      ; anova(m0, test = "Chisq")   
m1 = glm(y~scaledCY*agg_treat,          family = "binomial", data = PER)      ; anova(m0,m1, test = "Chisq")  #interaction improves the model
m2 = glm(y~scaledCY*agg_treat+ITD,      family = "binomial", data = PER)      ; anova(m1,m2, test = "Chisq")  #ITD does not improve the model
m3 = glm(y~poly(scaledCY,2)*agg_treat,  family = "binomial", data = PER)      ; anova(m1,m3, test = "Chisq")  #a polynomial does not improve the model

summary(m1)

#plot the data:
newdata <- data.frame(expand.grid(
  agg_treat = levels(PER$agg_treat),
  scaledCY     = seq(0.01, 0.035, by=0.001)
))

newdata$fit    = predict(m1, newdata, se.fit = TRUE, interval = "predict", type = "response")$fit
newdata$se.fit = predict(m1, newdata, se.fit = TRUE, interval = "predict", type = "response")$se.fit

#
ggplot()+ 
  geom_point(data = PER, aes(x = scaledCY, y = sum_score/10, group = agg_treat, shape = agg_treat,color = agg_treat))+
  #ylim(0,1)+
  geom_ribbon(data = newdata, aes(x = scaledCY, ymin=newdata$fit-(newdata$se.fit*1.96), ymax=newdata$fit+(newdata$se.fit*1.96),
                                 group = agg_treat, fill = agg_treat),alpha = 0.05) +
  geom_line(  data = newdata, aes(x = scaledCY, y = fit, group = agg_treat, linetype=agg_treat, color = agg_treat) )+
  theme_classic(base_size = 16)+theme(legend.title = element_blank())+
  labs(x = "Scaled calyces size", y = "learning score" )+
  ggsave(paste0(base,"/02_output/22_Calyces_and_learning_score_.png"), height = 6, width = 8)


write.csv(as.data.frame(summary(m1)$coef), paste0(base,"/02_output/21_Calyces_and_learning_score_LM.csv") )