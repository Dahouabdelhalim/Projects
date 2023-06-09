# Data analysis script for Nettle and Fraser 'Hunger affects cooperation in a Public Goods Game but not an Ultimatum Game'
# Daniel Nettle / Sam Fraser
# Version of April 2020

###Head#######
# Packages
library(dplyr)
library(ggplot2)
library(psych)
library(lmerTest)
library(cowplot)
library(tidyr)
library(reshape2)
library(metafor)
library(mosaic)
theme_set(theme_bw())

# Get experiment 1 data
d1<-read.csv("Experiment_1.csv")
d1$QuestionsCorrect=8*d1$QuestionsCorrect
d1$Condition.Name="Breakfast"
d1$Condition.Name[d1$Condition==0]="No Breakfast"
d1$Condition.Name=factor(d1$Condition.Name)

# Get experiment 2 data
d2=read.csv("Experiment_2.csv")
d2$UniqueID=factor(paste(d2$Experiment, d2$Participant.Number))
d2$Sample="Sample 1"
d2$Sample[d2$Experiment=="Replication"]="Sample 2"
d2$Previous.Rounds=d2$Round-1
ph=subset(d2, Punishment.Round=="Y")
ph$Sent.Punishment=(as.numeric(as.character(ph$Punishment.Sent))>0)
ph$MCO = (ph$Group.Contribution-ph$Contribution)/3
# t is the latency data for experiment 2
t=read.csv("Experiment_2_time.csv")
t$Condition.Name="Breakfast"
t$Condition.Name[t$Condition==0]="No breakfast"
t$Round[t$Trial.Number<11]=t$Trial.Number[t$Trial.Number<11]
t$Round[t$Trial.Number>10]=t$Trial.Number[t$Trial.Number>10]-10
t$Previous.Rounds=t$Round-1

#### Experiment 1 ####
# Descriptives
xtabs(~d1$Gender)
xtabs(~d1$Condition)

###Manipulation check#####
table(d1$BreakfastToday, d1$Condition.Name)
prop.table(table(d1$BreakfastToday, d1$Condition.Name), margin=2)
summary(lm(HowHungry~Condition.Name, data=d1))
describeBy(d1$HowHungry, d1$Condition.Name)
t.test(d1$HowHungry ~ d1$Condition.Name, var.equal=T)
describeBy(d1$ProposedAmount, d1$ProposerFirst)
describeBy(d1$LowestAcceptable, d1$ProposerFirst)
cor.test(d1$ProposedAmount, d1$LowestAcceptable)

### Experimental effects #####
# Proposer
summary(lm(ProposedAmount~Condition.Name, data=d1))
describeBy(d1$ProposedAmount, d1$Condition.Name)
# Dichotomous outcome
d1$Proposer.Cut=cut(d1$ProposedAmount, c(-0.5, 4.75, 10))
table(d1$ProposedAmount)/106
table(d1$Proposer.Cut)
table(d1$Proposer.Cut, d1$Condition.Name)
prop.table(table(d1$Proposer.Cut, d1$Condition.Name), margin=2)
chisq.test(d1$Proposer.Cut, d1$Condition.Name)

# Responder role
# Lowest acceptable
summary(lm(LowestAcceptable~Condition.Name, data=d1))
describeBy(d1$LowestAcceptable, d1$Condition.Name)
# Dichotomise into those who demanded at least 50% and those who did not
d1$Lowestcut=cut(d1$LowestAcceptable, c(-0.5, 4.75, 10))
prop.table(table(d1$Lowestcut, d1$Condition.Name), margin=2)
chisq.test(d1$Lowestcut, d1$Condition.Name)

####Figure 1#####
# Make summary by condition and amount
fig1.dat=summarise(group_by(d1, ProposedAmount, Condition.Name), n=n())
fig1.dat$key=paste(fig1.dat$Condition.Name, fig1.dat$ProposedAmount)
fig1.dat$key
# Make data frame of possible contributions
Amount=rep(seq(0, 10, by=0.5), times=2)
Condition=c(rep("Breakfast", times=21), rep("No Breakfast", times=21))
x=data.frame(Amount, Condition)
x$key=paste(x$Condition, x$Amount)
# Merge
fig1.dat=merge(x, fig1.dat, by="key", all.x=T)
fig1.dat$n[is.na(fig1.dat$n)]=0
fig1.dat$prop=fig1.dat$n/60
fig1.dat$prop[fig1.dat$Condition=="No Breakfast"]=fig1.dat$n[fig1.dat$Condition=="No Breakfast"]/46
# Now make figure
fig1=ggplot(fig1.dat, aes(x=Amount, y=prop, fill=Condition)) + 
  geom_bar(stat="identity", position=position_dodge()) + 
  scale_fill_manual(values=c("blue", "red")) + 
  guides(fill=FALSE) + 
  ylab("Proportion of proposers") + 
  xlab("Amount proposed")
fig1

#Lowest acceptable offer (fig1b)
# Make summary by condition and amount
fig1b.dat=summarise(group_by(d1, LowestAcceptable, Condition.Name), n=n())
fig1b.dat$key=paste(fig1b.dat$Condition.Name, fig1b.dat$LowestAcceptable)
# Make data frame of possible acceptables
Lowest=rep(seq(0, 10, by=0.5), times=2)
Condition=c(rep("Breakfast", times=21), rep("No Breakfast", times=21))
xb=data.frame(Lowest, Condition)
xb$key=paste(xb$Condition, xb$Lowest)
# Merge
fig1b.dat=merge(xb, fig1b.dat, by="key", all.x=T)
fig1b.dat$n[is.na(fig1b.dat$n)]=0
fig1b.dat$prop=fig1b.dat$n/60
fig1b.dat$prop[fig1b.dat$Condition=="No Breakfast"]=fig1b.dat$n[fig1b.dat$Condition=="No Breakfast"]/46
# Now make figure
fig1b=ggplot(fig1b.dat, aes(x=Lowest, y=prop, fill=Condition)) + 
  geom_bar(stat="identity", position=position_dodge()) + 
  scale_fill_manual(values=c("blue", "red")) + 
  ylab("Proportion of responders") + 
  xlab("Minumum acceptable") + 
  theme(legend.position = c(0.825, 0.5))
fig1b

png("fig1.png", res=300, , width=10, height=4, units="in")
plot_grid(fig1, fig1b, labels=c("A", "B"))
dev.off()


##### Decision latencies####
summary(lm(log(TimeProposer)~Condition.Name, data=d1))
describeBy(log(d1$TimeProposer), d1$Condition.Name)
summary(lm(log(TimeResponder)~Condition.Name, data=d1))
describeBy(log(d1$TimeResponder), d1$Condition.Name)
cor.test(log(d1$TimeProposer), d1$ProposedAmount)
cor.test(log(d1$TimeResponder), d1$LowestAcceptable)

###Cognitive reflection task#####
summary(lm(d1$QuestionsCorrect~d1$Condition.Name))
describeBy(d1$QuestionsCorrect, d1$Condition.Name)
summary(lm(log(d1$TotalQuestionTime)~d1$Condition.Name))
describeBy(log(d1$TotalQuestionTime), d1$Condition.Name)
cor.test(d1$QuestionsCorrect, log(d1$TotalQuestionTime))
cor.test(d1$QuestionsCorrect, d1$ProposedAmount)
cor.test(d1$QuestionsCorrect, d1$LowestAcceptable)

###Using hunger instead of treatment#####
# Proposer
summary(lm(ProposedAmount~HowHungry, data=d1))
summary(lm(LowestAcceptable~HowHungry, data=d1))

### Experiment 2#####
###Manipulation checks####
# Overall
mc.data=subset(d2, Round==1 & Punishment.Round=="N")
table(mc.data$BreakfastToday, mc.data$Condition.Name)
prop.table(table(mc.data$BreakfastToday, mc.data$Condition.Name), margin=2)
t.test(mc.data$HowHungry~mc.data$Condition.Name, var.equal=T)
describeBy(mc.data$HowHungry, mc.data$Condition.Name)

##Comprehension questions#####
comp.data=summarise(group_by(d2, UniqueID), Condition.Name=last(Condition.Name), Sample=last(Sample), comp1=6*mean(Non.Punishment.Questions.Correct), comp2 = 2*mean(Punishment.Questions.Correct), 
                    check1=sd(Non.Punishment.Questions.Correct), check2 = sd(Punishment.Questions.Correct),)
# Comparison by condition
psych::describeBy(comp.data$comp1, comp.data$Condition.Name)
t.test(comp.data$comp1~comp.data$Condition.Name, var.equal=T)
psych::describeBy(comp.data$comp2, comp.data$Condition.Name)
t.test(comp.data$comp2~comp.data$Condition.Name, var.equal=T)
# Effect of comphrension questions correct on contribution
mcomp=lmer(Contribution~Previous.Rounds*Punishment.Round + I(6*Non.Punishment.Questions.Correct) + (1|UniqueGroup/UniqueID), data=d2, REML=F)
summary(mcomp)
# Effect of comprehension of punishment on decision to punish
mpcomp=glmer(Sent.Punishment ~    
           scale(MCO) + scale(Punishment.Questions.Correct) +  
           (1|UniqueGroup/UniqueID), data=ph, family=binomial, na.action="na.fail")
summary(mpcomp)

###Overall analysis (model 1)#####
m1=lmer(Contribution~Previous.Rounds*Punishment.Round*Condition.Name + (1|UniqueGroup/UniqueID), data=d2, REML=F)
summary(m1)
describeBy(d2$Contribution[d2$Round==1], paste(d2$Condition.Name[d2$Round==1], d2$Punishment.Round[d2$Round==1]))
#  Figure 2a
round.summarya=summarise(group_by(d2, Punishment.Round, Round, Condition.Name, UniqueGroup), 
                         Contribution=mean(Contribution))
round.summary=summarise(group_by(round.summarya, Punishment.Round, Round, Condition.Name), 
                        Mean.Contribution=mean(Contribution), se=describe(Contribution)$se)
round.summary$Game="Punishment"
round.summary$Game[round.summary$Punishment.Round=="N"]="No punishment"
fig2a=ggplot(round.summary, aes(x=factor(Round), y=Mean.Contribution, colour=Condition.Name)) + 
  geom_point(aes(shape=Game)) + 
  geom_line(aes(x=Round, linetype=Game)) + 
  coord_cartesian(ylim=c(0, 20)) + 
  geom_errorbar(aes(ymax=Mean.Contribution+se, ymin=Mean.Contribution-se), width=0.1) + 
  xlab("Round") + 
  ylab("Mean Contribution") + 
  scale_colour_manual(name="Condition", values=c("blue", "red")) 
fig2a

### Model 2 (no punishment game, round 1)####
m2 = lm(Contribution~Condition.Name, data=subset(d2, Punishment.Round=="N" & Previous.Rounds==0))
summary(m2)
fig2bdata=summarise(group_by(subset(d2, Punishment.Round=="N" & Round==1), Condition.Name), mc=describe(Contribution)$mean, se=describe(Contribution)$se)
fig2b=ggplot(as.data.frame(fig2bdata), aes(x=Condition.Name, y=mc, fill=Condition.Name)) + 
  geom_bar(stat="identity") + 
  geom_errorbar(aes(ymax=mc+se, ymin=mc-se), width=0.25) + 
  ylab("Contribution") + 
  xlab("Condition") + 
  guides(fill=F) + 
  scale_fill_manual(values=c("blue", "red")) + 
  scale_x_discrete(labels=c("Brk.", "No Brk." ))
fig2b

#### Model 3 (no punishment game, after round 1)#####
m3=lmer(Contribution~Lagged.Contribution+Lagged.MCO*Condition.Name + (1|UniqueGroup/UniqueID), data=subset(d2, Punishment.Round=="N"), REML=F)
summary(m3)
# Figure 2c
fig2c=ggplot(data=subset(d2, Punishment.Round=="N"), aes(x=Lagged.MCO, y=Contribution, colour=Condition.Name)) + 
  geom_smooth(method="lm") + 
  xlab("Lagged contribution of others") + 
  scale_colour_manual(name="Condition", values=c("Blue", "Red")) + 
  geom_abline(intercept=0, slope=1, linetype="dotted")  +
  coord_cartesian(xlim=c(0, 20), ylim=c(0, 20)) + 
  guides(colour=FALSE)
fig2c

### Model 4 (punishment game, decision to punish)#### 
ph=subset(d2, Punishment.Round=="Y")
ph$Sent.Punishment=(as.numeric(as.character(ph$Punishment.Sent))>0)
table(ph$Sent.Punishment, ph$Condition.Name)
colSums(table(ph$Sent.Punishment, ph$Condition.Name))
prop.table(table(ph$Sent.Punishment, ph$Condition.Name), margin=2)
chisq.test(ph$Sent.Punishment, ph$Condition.Name)
ph$MCO = (ph$Group.Contribution-ph$Contribution)/3
m4=glmer(Sent.Punishment ~    
           scale(MCO)*Condition.Name +  
           (1|UniqueGroup/UniqueID), data=ph, family=binomial, na.action="na.fail")
summary(m4)
# Figure 2d
ph$MCO.bin=cut(ph$MCO, c(-1, 5, 10, 15, 20))
ph$Lagged.Contribution.Bin=cut(ph$Lagged.Contribution, c(-1, 5, 10, 15, 20))
fig2d.dat=summarise(group_by(ph, MCO.bin, Condition.Name), mean.cont=mean(Sent.Punishment))
fig2d.dat$midpoint=c(2.5, 2.5, 7.5, 7.5, 12.5, 12.5, 17.5, 17.5)
fig2d=ggplot(fig2d.dat, aes(x=midpoint, y=mean.cont, fill=Condition.Name)) + 
  geom_bar(stat="identity", position=position_dodge()) + 
  ylab("Probability of punishing") + 
  xlab("Mean contribution of others") + 
  scale_fill_manual(name="Condition", values=c("blue", "red"), labels=c("Breakfast", "No breakfast"))  +
  guides(fill=FALSE)
fig2d

### Model 5 (response to punishment)####
m5=lmer(Contribution~ Lagged.Contribution*Condition.Name*Lagged.Punished + (1|UniqueGroup/UniqueID), data=ph, REML=F)
summary(m5)


####Assemble figure 2 (first run previous sections)#####
bottom.part=plot_grid(fig2b, fig2c, fig2d, labels=c("B", "C", "D"), ncol=3, scale=0.9, rel_widths = c(0.2, 0.4, 0.4))
fig2=plot_grid(fig2a, bottom.part, ncol=1, rel_heights = c(0.6, 0.4), labels=c("A", ""))
fig2
png("figure2.png", res=300, width=8*300, height=8*300)
fig2
dev.off()
###Latency to make decisions (models 6 and 7)#####
# Latency to decide on a contribution
options(scipen=99)
m6 = lmer(log(TimeContribution)~Previous.Rounds*Punishment.Round*Condition.Name + (1|Group/Individual), data=t, REML=F)
summary(m6)
summarise(group_by(t, Round, Punishment.Round, Condition.Name), log.latency=mean(log(TimeContribution)), sd=describe(log(TimeContribution))$sd)
summarise(group_by(t, Punishment.Round, Condition.Name), log.latency=mean(log(TimeContribution)), sd=describe(log(TimeContribution))$sd)
summarise(group_by(t, Punishment.Round, Condition.Name), latency=mean(TimeContribution), sd=sd(TimeContribution))

# Does latency to decide predict how much a person contributes?
cor.test(log(t$TimeContribution[t$Punishment.Round=="N" & t$Round==1]), t$Contribution[t$Punishment.Round=="N" & t$Round==1])

# What about latency to decide whether to punish (model 7)?
m7 = lmer(log(TimePunishment)~Previous.Rounds*Condition.Name + (1|Group/Individual), data=subset(t, Punishment.Round=="Y"), REML=F)
summary(m7)

###Figure 3 (latency to decide) #####
round.summary.time.1=summarise(group_by(t, Punishment.Round, Round, Condition.Name, Group), 
                          Time=mean(log(TimeContribution)))
round.summary.time=summarise(group_by(round.summary.time.1, Punishment.Round, Round, Condition.Name), 
                          Mean.Time=mean(Time), se=describe(Time)$se)
round.summary.time$Game="Punishment"
round.summary.time$Game[round.summary.time$Punishment.Round=="N"]="No punishment"

fig3=ggplot(round.summary.time, aes(x=factor(Round), y=Mean.Time, colour=Condition.Name)) + 
  geom_point() + 
  geom_line(aes(x=Round)) + 
  geom_errorbar(aes(ymax=Mean.Time+se, ymin=Mean.Time-se), width=0.1) + 
  xlab("Round") + 
  ylab("Mean decision time (log s)") +
  scale_colour_manual(name="Condition", values=c("blue", "red")) + 
  facet_wrap( ~ Game)  
png("fig3.png",res=300, units="in", width=8, height=4)
fig3
dev.off()

###Comparison of two experiments#####
# First re-analyse the main results with contributions scaled
# 1. Proposer in exp 1
c1=lm(scale(ProposedAmount)~Condition.Name, data=d1)
summary(c1)
# 2. Responder in exp 1
c2=lm(scale(LowestAcceptable)~Condition.Name, data=d1)
# Now exp 2
# 3. Contribution in first no punishment game
c3=lm(scale(Contribution)~Condition.Name, data=subset(d2, Punishment.Round=="N" & Previous.Rounds==0))
# 4. Response to MCO, no punishment game
c4=lmer(scale(Contribution)~scale(Lagged.Contribution)+scale(Lagged.MCO)*Condition.Name + (1|UniqueGroup/UniqueID), data=subset(d2, Punishment.Round=="N"), REML=F)
summary(c4)
# 5. Willingness to punish
# Converted from a log odds ratio
table(ph$Condition.Name, ph$Sent.Punishment)
orrr(table(ph$Sent.Punishment, ph$Condition.Name), verbose=TRUE)
# Manually calculate odds ratio
oddspunbreakfast=717/643
oddspunnobreak = 575/705
or=oddspunnobreak/oddspunbreakfast
# So: the betas and ses are:
betas=c(summary(c1)$coefficients[2, 1],
summary(c3)$coefficients[2, 1],
summary(c4)$coefficients[5, 1],
summary(c2)$coefficients[2, 1], 
(log(or))/1.81)
# Ses
ses=c(summary(c1)$coefficients[2, 2],
        summary(c3)$coefficients[2, 2],
        summary(c4)$coefficients[5, 2],
        summary(c2)$coefficients[2, 2], 
      ((log(or))/1.81-(log(1/1.593))/1.81)/1.96)

re1=rma(yi=betas[1:2], sei=ses[1:2], method="REML")
re2=rma(yi=betas[4:5], sei=ses[4:5], method="REML")
summary(re2)
png("figure4.png", res=300, units="in", width=8, height=3)
par(mar=c(4,4,1,2))
forest.default(x=betas, sei=ses, psize=1, 
               slab=c("Exp. 1 proposer", 
                      "Exp. 2. First round, no pun.", 
                      "Exp. 2. Response to others, no pun.",
                      "Exp. 1 responder",
                      "Exp. 2 propensity to punish"), 
               ylim=c(-1, 12), 
               rows=c(9, 8, 5, 2, 1), 
               xlab="Effect size (95% confidence interval)", cex=0.9)
addpoly(re1, row=7, mlab=" ", cex=0.9)
addpoly(re2, row=0, mlab=" ", cex=0.9)
dev.off()

