#Code to replicate analyses in Reichert et al. "No reproductive benefits of dear enemy recognition in a territorial songbird"
#Run in R version 4.0.2
####Load packages####

library(rptR) #Version 0.9.22
library(lme4) #Version 1.1-23
library(lmerTest) #Version 3.1-2
library(ordinal) #Version 2019.12-10
library(glmmTMB) #Version 1.0.2.1
library(censReg) #Version 0.5-30
library(tidyverse) #Version 1.3.0
library(magrittr) #Version 1.5
library(cowplot) #Version 1.0.0

####Load Data####

DryadData <- read_csv("FinalDryadData.csv") 
IndividualChickData <- read_csv("FinalIndividualChickData.csv") 

#Add on a convenience variable that can be used as proxy for individual ID (the nestbox + year of testing)
DryadData$BoxYear <- paste(DryadData$Box,DryadData$Year, sep = ".")

#Lazily adding on a variable for number of extrapair fledglings
FledgeData<-IndividualChickData %>%
  group_by(BoxYear) %>%
  summarize(NumberOwnFledged = sum(Fledged=="Y" & EPY == "No"), NumberEPY = sum(EPY=="Yes"), NumberNotEPY = sum(EPY=="No"), NumExtraPairFledglings = sum(EPY=="Yes" & Fledged =="Y"))

DryadData$NumExtraPairFledglings<-FledgeData$NumExtraPairFledglings[match(DryadData$BoxYear, FledgeData$BoxYear)]
DryadData$NumExtraPairFledglings[is.na(DryadData$NumExtraPairFledglings)]<-0

#Make the approach category variable into an ordered factor, so later models behave:
DryadData$ApproachCat<-as.factor(DryadData$ApproachCat)
DryadData$ApproachCat<-ordered(DryadData$ApproachCat,levels=c(5,10,15,20,25,30,"None"))
#Reordering the stimulus rate variable, so that the low rate will be the reference in models:
DryadData$StimRate <- as.factor(DryadData$StimRate)
DryadData$StimRate <-factor(DryadData$StimRate, levels(DryadData$StimRate)[c(2,1)])

####Table 1####
#For viewing any model results, use summary(ModelName)

ResponseDecrease<-glmer(ResponseBinary~TotalTrial+StimRate+(1|BoxYear),data=filter(DryadData,Phase=="H1"),family=binomial(link="logit")) 

ApproachDecrease<-clmm(ApproachCat~TotalTrial+StimRate+(1|BoxYear),data=filter(DryadData,Phase=="H1"))

SongDecrease<-glmmTMB(NumSongs~TotalTrial+StimRate+(1|BoxYear),data=filter(DryadData,Phase=="H1"), family="poisson", ziformula = ~1) 

####Figure 2####

Fig2A<-DryadData %>%
  filter(Phase=="H1" | Phase=="D1") %>% 
  group_by(TotalTrial) %>%
  summarize(PropResponding = mean(ResponseBinary), count = n()) %>%
  mutate(SE = sqrt(PropResponding*(1-PropResponding)/count), Day = c(rep(1:3,each=3),3)) %>%
  
  ggplot(aes(as.factor(TotalTrial),PropResponding,color=as.factor(Day))) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin=PropResponding-SE,ymax=PropResponding+SE),width=0.3, size=1) +
  scale_colour_manual(values=c("black","blue","red")) +
  labs(x = "Trial Number", y = "Proportion Responding") +
  scale_y_continuous(limits=c(0,1)) +
  scale_x_discrete(breaks = seq(1,10,by=1),labels=c("N1","N2","N3","N4","N5","N6","N7","N8","N9","S")) +
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(), panel.background=element_blank(),axis.text.x=element_text(colour='black',size=20), axis.text.y=element_text(colour='black',size=20), axis.title.x=element_text(colour='black',size=24,vjust=0.1,margin = margin(t=10,r=0,b=0,l=0)), axis.title.y=element_text(colour='black',size=24,vjust=0.3,margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.ticks=element_line(colour='black',size=1,linetype=1), axis.ticks.length=unit(0.3,"cm"),legend.position="none") +
  theme(axis.line.x=element_line(colour='black',size=0.5,linetype=1),axis.line.y=element_line(colour='black',size=0.5,linetype=1)) 

Fig2A

Fig2B<-DryadData %>%
  filter(Phase=="H1" | Phase=="D1") %>% 
  
  ggplot(aes(as.factor(TotalTrial),NumSongs,color=as.factor(Day))) + 
  geom_jitter(size=3,alpha=0.3,width=0.2,height=0) +
  stat_summary(fun = mean,fun.min=mean,fun.max=mean, geom ="crossbar",width = 0.5,color=c(rep("black",3),rep("blue",3),rep("red",4)))+
  scale_colour_manual(values=c("black","blue","red")) +
  labs(x = "Trial Number", y = "Number Songs") +
  scale_x_discrete(breaks = seq(1,10,by=1),labels=c("N1","N2","N3","N4","N5","N6","N7","N8","N9","S")) + 
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(), panel.background=element_blank(),axis.text.x=element_text(colour='black',size=20), axis.text.y=element_text(colour='black',size=20), axis.title.x=element_text(colour='black',size=24,vjust=0.1,margin = margin(t=10,r=0,b=0,l=0)), axis.title.y=element_text(colour='black',size=24,vjust=0.3,margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.ticks=element_line(colour='black',size=1,linetype=1), axis.ticks.length=unit(0.3,"cm"),legend.position="none") +
  theme(axis.line.x=element_line(colour='black',size=0.5,linetype=1),axis.line.y=element_line(colour='black',size=0.5,linetype=1)) 

Fig2B

Figure2<-cowplot::plot_grid(Fig2A, Fig2B, labels = c('A', 'B'), label_size = 24, nrow=1)

#ggsave(filename = "Figure2.pdf", plot = Figure2, device=cairo_pdf, width=20, height=8 )
Figure2

####Table 2####
ResponseLastHvsD<-glmer(ResponseBinary~as.factor(TotalTrial)+StimRate+(1|BoxYear),data=filter(DryadData, Phase=="H1"  & TotalTrial>=9| Phase == "D1" & TotalTrial>=9),family=binomial(link="logit"),glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))) 

ApproachLastHvsD<-clmm(ApproachCat~as.factor(TotalTrial)+StimRate+(1|BoxYear),data=filter(DryadData, Phase=="H1"  & TotalTrial>=9| Phase == "D1" & TotalTrial>=9)) 

NumSongsLastHvsD<-glmmTMB(NumSongs~as.factor(TotalTrial)+StimRate+(1|BoxYear),data=filter(DryadData, Phase=="H1"  & TotalTrial>=9| Phase == "D1" & TotalTrial>=9), family="poisson", ziformula = ~1) 

####Stats on who learned and who didn't####
DryadData %>%
  filter(Phase=="H1"  & TotalTrial>=9| Phase == "D1" & TotalTrial>=9) %>%
  group_by(BoxYear) %>%
  mutate(Learned = if_else(ResponseBinary[1]==0 & ResponseBinary[2]==1, "y", "n")) %>%
  ungroup() %>%
  filter(Phase=="H1") %>%
  summarize(Learners=sum(Learned=="y"), NotLearners = sum(Learned=="n")) 

####Testing other possible effects on learning####
#Make a few more variables to help out.
DryadData<-DryadData%>%
  group_by(BoxYear, Repetition) %>%
  mutate(RespondedNinth = ifelse(ResponseBinary[TotalTrial==9] == 1, "y", "n"), RespondedTenth = ifelse(ResponseBinary[TotalTrial==10] == 1, "y", "n")) 

DryadData<-DryadData%>%
  group_by(BoxYear, Repetition) %>%
  mutate(Nonresponders = ifelse(RespondedNinth=="n" & RespondedTenth=="n", "y", "n")) %>%
  ungroup()

DryadData$Learned<-as.factor(DryadData$Learned)

####Table S2####
#Standard criterion
LearningEffects<-glm(Learned ~ Year +  StartDateOfExpt + StimRate + ExptStartDateRelativeToLayDate, data = filter(DryadData,Phase=="H1", TotalTrial==9), family="binomial")
#Standard criterion without nonresponders
LearningEffectsNN<-glm(Learned ~ Year +  StartDateOfExpt + StimRate + ExptStartDateRelativeToLayDate, data = filter(DryadData,Phase=="H1", TotalTrial==9, Nonresponders=="n"), family="binomial")
#Response slope
LearningEffectsRS<-lm(ResponseSlope ~ Year +  StartDateOfExpt + StimRate + ExptStartDateRelativeToLayDate, data = filter(DryadData,Phase=="H1", TotalTrial==9 & Nonresponders=="n"))
#Song slope
LearningEffectsSS<-lm(SongSlope ~ Year +  StartDateOfExpt + StimRate + ExptStartDateRelativeToLayDate, data = filter(DryadData,Phase=="H1", TotalTrial==9 & Nonresponders=="n"))

#Effects of Age, reported in text
#Standard criterion
LearningPhysicalAgeOnly<-glm(Learned~age,data=filter(DryadData,Phase=="H1", TotalTrial==9),family=binomial(link="logit")) 
#Response slope
LearningPhysicalAgeOnlyRS<-lm(ResponseSlope~age,data=filter(DryadData,Phase=="H1", TotalTrial==9 & Nonresponders=="n"))
#Song slope
LearningPhysicalAgeOnlySS<-lm(SongSlope~age,data=filter(DryadData,Phase=="H1", TotalTrial==9 & Nonresponders=="n")) 


####Repeatability of learning####
#Standard criterion
LearningRepeatability <- rpt(as.factor(Learned) ~ Phase + StimRate + (1|BoxYear), grname = "BoxYear", data = filter(DryadData, BoxYear %in% filter(DryadData, Phase=="H2")$BoxYear & TotalTrial==9 & BoxYear!="CB1.2017"), datatype="Binary", nboot=1000,npermut=0) #As reported in the text, the individual that did not respond at all during the second repetition is removed (CB1.2017)

#Numbers that learned or didn't learn in both sets of playbacks.
table(H1Learned = DryadData$Learned[DryadData$Phase=="H1" & DryadData$TotalTrial==9 & DryadData$BoxYear %in% filter(DryadData, Phase=="H2")$BoxYear & DryadData$BoxYear!="CB1.2017"], H2Learned = DryadData$Learned[DryadData$Phase=="H2" & DryadData$TotalTrial==9 & DryadData$BoxYear %in% filter(DryadData, Phase=="H2")$BoxYear & DryadData$BoxYear!="CB1.2017"])

#Repeatabilities for the slope measures. As reported in the text, there aren't enough numbers of individuals here, and if you run this it's one long string of singularity errors. 
ResponseSlopeRepeatability <- rpt(ResponseSlope ~ Phase  + StimRate + (1|BoxYear), data=filter(DryadData, BoxYear %in% filter(DryadData, Phase=="H2")$BoxYear & TotalTrial==9 & BoxYear!="CB1.2017" & Nonresponders=="n"), grname = "BoxYear", datatype = "Gaussian", nboot = 1000, npermut = 0) 

SongSlopeRepeatability <- rpt(SongSlope~ Phase + StimRate + (1|BoxYear), data=filter(DryadData, BoxYear %in% filter(DryadData, Phase=="H2")$BoxYear & TotalTrial==9 & BoxYear!="CB1.2017" & Nonresponders=="n"), grname = "BoxYear", datatype = "Gaussian", nboot = 1000, npermut = 0) 
####Extrapair young####
LearnEPY<-glm(cbind(NumExtraPairFledglings,NumOwnFledged)~Learned,data=filter(DryadData,Phase=="H1", TotalTrial==9),family="binomial")

#Number with unknown paternity (see methods):
IndividualChickData %>%
  summarize(NumberFledged = sum(Fledged=="Y", na.rm=T), NumberFledgedUnknownPaternity = sum(Fledged=="Y" & KnownEPStatus=="N", na.rm=T), NumberWeighed = sum(!is.na(OffspringMass)), NumberWeighedUnknownPaternity = sum(!is.na(OffspringMass) & KnownEPStatus=="N"))

#Number that actually were extrapair
IndividualChickData$IncludeInFitnessAnalysis<-DryadData$IncludeInFitnessAnalysis[DryadData$Repetition==1][match(IndividualChickData$BoxYear, DryadData$BoxYear[DryadData$Repetition==1])]
IndividualChickData$Learned<-DryadData$Learned[DryadData$Repetition==1][match(IndividualChickData$BoxYear, DryadData$BoxYear[DryadData$Repetition==1])]
IndividualChickData$Nonresponders<-DryadData$Nonresponders[DryadData$Repetition==1][match(IndividualChickData$BoxYear, DryadData$BoxYear[DryadData$Repetition==1])]
IndividualChickData$SongSlope<-DryadData$SongSlope[DryadData$Repetition==1][match(IndividualChickData$BoxYear, DryadData$BoxYear[DryadData$Repetition==1])]
IndividualChickData$ResponseSlope<-DryadData$ResponseSlope[DryadData$Repetition==1][match(IndividualChickData$BoxYear, DryadData$BoxYear[DryadData$Repetition==1])]

IndividualChickData %>%
  summarize(NumberEPY = sum(EPY=="Yes" & IncludeInFitnessAnalysis=="y", na.rm=T))

####Fitness variables####
####Table S3####
#Standard criterion
ClutchLearned<-glm(ClutchSize~Learned,data=filter(DryadData,Phase=="H1", TotalTrial==9, IncludeInFitnessAnalysis == "y"),family="poisson") 

ZeroInflEPY<-glmmTMB(NumOwnFledged~Learned,data=filter(DryadData,Phase=="H1", TotalTrial==9, IncludeInFitnessAnalysis == "y"),family="poisson", ziformula = ~1) 

TotalFledgeMassTobitEPY<-censReg(TotalFledglingBiomassNoEPY~Learned,data=filter(DryadData,Phase=="H1", TotalTrial==9, IncludeInFitnessAnalysis == "y", BoxYear!="IN001.2018")) 

LearnedAvgMass<-glm(AvgFledglingWeightNoEPY~Learned+BroodSize, data=filter(DryadData,Phase=="H1", TotalTrial==9, IncludeInFitnessAnalysis == "y")) 

IndividualMass<-lmer(OffspringMass~Learned + BroodSize + (1|BoxYear),data=filter(IndividualChickData, IncludeInFitnessAnalysis == "y" & EPY == "No")) 

#Excluding those birds that did not respond to stranger. (N=29 remain, and N=19 didn't respond to stranger). 

ClutchLearnedNN<-glm(ClutchSize~Learned,data=filter(DryadData,Phase=="H1", TotalTrial==9, IncludeInFitnessAnalysis == "y", Nonresponders=="n"),family="poisson") 

ZeroInflEPYNN<-glmmTMB(NumOwnFledged~Learned,data=filter(DryadData,Phase=="H1", TotalTrial==9, IncludeInFitnessAnalysis == "y", Nonresponders=="n"),family="poisson", ziformula = ~1) 

TotalFledgeMassTobitEPYNN<-censReg(TotalFledglingBiomassNoEPY~Learned,data=filter(DryadData,Phase=="H1", TotalTrial==9, IncludeInFitnessAnalysis == "y", Nonresponders=="n", BoxYear!="IN001.2018")) 

LearnedAvgMassNN<-glm(AvgFledglingWeightNoEPY~Learned+BroodSize, data=filter(DryadData,Phase=="H1", TotalTrial==9, IncludeInFitnessAnalysis == "y",  Nonresponders=="n")) 

IndividualMassNN<-lmer(OffspringMass~Learned + BroodSize + (1|BoxYear),data=filter(IndividualChickData, IncludeInFitnessAnalysis == "y" & EPY == "No" & Nonresponders=="n")) 

#Response Slope

ClutchLearnedResponseSlope<-glm(ClutchSize~ResponseSlope,data=filter(DryadData,Phase=="H1", TotalTrial==9, IncludeInFitnessAnalysis == "y", Nonresponders=="n"),family="poisson") 

ZeroInflEPYResponseSlope<-glmmTMB(NumOwnFledged~ResponseSlope,data=filter(DryadData,Phase=="H1", TotalTrial==9, IncludeInFitnessAnalysis == "y", Nonresponders=="n"),family="poisson", ziformula = ~1) 

TotalFledgeMassTobitEPYResponseSlope<-censReg(TotalFledglingBiomassNoEPY~ResponseSlope,data=filter(DryadData,Phase=="H1", TotalTrial==9, IncludeInFitnessAnalysis == "y", Nonresponders=="n", BoxYear!="IN001.2018"))

LearnedAvgMassRS<-glm(AvgFledglingWeightNoEPY~ResponseSlope+BroodSize, data=filter(DryadData,Phase=="H1", TotalTrial==9, IncludeInFitnessAnalysis == "y", Nonresponders=="n")) 

IndividualMassRS<-lmer(OffspringMass~ResponseSlope + BroodSize + (1|BoxYear),data=filter(IndividualChickData, IncludeInFitnessAnalysis == "y" & EPY == "No" & Nonresponders=="n")) 

#Song Slope

ClutchLearnedSongSlope<-glm(ClutchSize~SongSlope,data=filter(DryadData,Phase=="H1", TotalTrial==9, IncludeInFitnessAnalysis == "y", Nonresponders=="n"),family="poisson") 

ZeroInflEPYSongSlope<-glmmTMB(NumOwnFledged~SongSlope,data=filter(DryadData,Phase=="H1", TotalTrial==9, IncludeInFitnessAnalysis == "y", Nonresponders=="n"),family="poisson", ziformula = ~1) 

TotalFledgeMassTobitEPYSongSlope<-censReg(TotalFledglingBiomassNoEPY~SongSlope,data=filter(DryadData,Phase=="H1", TotalTrial==9, IncludeInFitnessAnalysis == "y", Nonresponders=="n", BoxYear!="IN001.2018"))

LearnedAvgMassSS<-glm(AvgFledglingWeightNoEPY~SongSlope+BroodSize, data=filter(DryadData,Phase=="H1", TotalTrial==9, IncludeInFitnessAnalysis == "y", Nonresponders=="n")) 

IndividualMassSS<-lmer(OffspringMass~SongSlope + BroodSize + (1|BoxYear),data=filter(IndividualChickData, IncludeInFitnessAnalysis == "y" & EPY == "No" & Nonresponders=="n")) 

####Figure 3####
Fig3A <- ggplot(filter(DryadData,Phase=="H1", TotalTrial==9, IncludeInFitnessAnalysis == "y"),aes(x=as.character(Learned),y=ClutchSize)) +
  geom_jitter(width=0.1,height=0,size=3,alpha=0.3) +
  stat_summary(fun = median,fun.min=median, fun.max=median, geom ="crossbar",width = 0.3,color="red")+
  ylab("Clutch size") +
  xlab("")+
  scale_x_discrete(labels=c("Non-learners", "Learners"))+
  scale_y_continuous(limits=c(0,9), breaks=c(0,3,6,9)) +
  theme_classic()+
  theme(axis.text.x = element_text(size=24,color="black"),axis.text.y=element_text(size=24,color="black"),axis.title.y=element_text(size=28,margin = margin(t = 0, r = 20, b = 0, l = 0)))

Fig3B <- ggplot(filter(DryadData,Phase=="H1", TotalTrial==9, IncludeInFitnessAnalysis == "y"),aes(x=as.character(Learned),y=NumOwnFledged)) +
  geom_jitter(width=0.1,height=0,size=3,alpha=0.3) +
  stat_summary(fun = median,fun.min=median,fun.max=median, geom ="crossbar",width = 0.3,color="red")+
  ylab("Number fledged") +
  xlab("")+
  scale_x_discrete(labels=c("Non-learners", "Learners"))+
  theme_classic()+
  theme(axis.text.x = element_text(size=24,color="black"),axis.text.y=element_text(size=24,color="black"),axis.title.y=element_text(size=28,margin = margin(t = 0, r = 20, b = 0, l = 0)))

Fig3C<-ggplot(filter(DryadData,Phase=="H1", TotalTrial==9, IncludeInFitnessAnalysis == "y", BoxYear!="IN001.2018"),aes(x=as.character(Learned),y=TotalFledglingBiomassNoEPY)) +
  geom_jitter(width=0.1,height=0,size=3,alpha=0.3) +
  stat_summary(fun = median,fun.min=median,fun.max=median, geom ="crossbar",width = 0.3,color="red")+
  ylab("Total offspring mass (g)") +
  xlab("")+
  scale_x_discrete(labels=c("Non-learners", "Learners"))+
  theme_classic()+
  theme(axis.text.x = element_text(size=24,color="black"),axis.text.y=element_text(size=24,color="black"),axis.title.y=element_text(size=28,margin = margin(t = 0, r = 20, b = 0, l = 0)))

Fig3D<-DryadData %>%
  filter(IncludeInFitnessAnalysis == "y" & TotalTrial==9 & Repetition==1) %>%
  ggplot( aes(x=as.character(Learned), y=AvgFledglingWeightNoEPY)) +
  geom_jitter(width=0.1,height=0,size=3,alpha=0.3) +
  stat_summary(fun = median,fun.min=median, fun.max=median, geom ="crossbar",width = 0.3,color="red")+
  ylab("Average Mass (g)") +
  xlab("")+
  scale_x_discrete(labels=c("Non-learners", "Learners"))+
  scale_y_continuous(limits=c(0,20), breaks=c(0,5,10,15,20)) +
  theme_classic()+
  theme(axis.text.x = element_text(size=24,color="black"),axis.text.y=element_text(size=24,color="black"),axis.title.y=element_text(size=28,margin = margin(t = 0, r = 20, b = 0, l = 0)))

Figure3<-cowplot::plot_grid(Fig3A, Fig3B, Fig3C, Fig3D, labels = c('A', 'B', 'C', 'D'), label_size = 24, nrow=2)

#ggsave(filename = "Figure3.pdf", plot = Figure3, device=cairo_pdf, width=12, height=12 )
Figure3

####Table S1####
#Note that the stoat predated nests were just removed manually for the purposes of this table (the 3 DD boxes from 2018).
TableS1 <-DryadData %>%
  group_by(Site) %>%
  filter(Phase=="H1" & TotalTrial==9) %>%
  summarize(NumberLearned = sum(Learned=="y"), NumberTested = n(), MeanClutch = mean(ClutchSize), sdClutch = sd(ClutchSize), MeanNumOwnFledged = mean(NumOwnFledged),  sdNumOwnFledged = sd(NumOwnFledged), MeanTotalMass = mean(TotalFledglingBiomassNoEPY), sdTotalFledglingBiomassNoEPY = sd(TotalFledglingBiomassNoEPY), MeanAvgMass = mean(AvgFledglingWeightNoEPY, na.rm=TRUE), sdAvgMass = sd(AvgFledglingWeightNoEPY, na.rm=TRUE))

