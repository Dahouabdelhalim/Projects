# Analysis for Dreyer and Shingleton, 2019

library(readxl)
library(lme4)
library(emmeans)
library(ROCR)
library(ggplot2)
library(lmerTest)
library(reshape2)
library(plyr)
library(gridExtra)
library(ez)
library(car)
#test for overdispersion in GLMs with Poisson and binomial link functions.
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

# Set the working directory that contains all the data files.
setwd("~/XXXXXX")
# All male data
allmales <- read.csv("allmales.csv")
data<-allmales

# Figure 1A and Table S1
# compare trait size among genotype including GFP allele
traitsize<-lmer(genital~mIIS+mGFP+(1|mvial), data=data)
  plot(traitsize)
  dev.off()
  qqnorm(resid(traitsize))
  qqline(resid(traitsize))
  shapiro.test(resid(traitsize))
anova(traitsize)
emmeans(traitsize, list(pairwise ~ mIIS))
detach("package:lmerTest", unload=TRUE)
traitsize<-lmer(genital~mIIS+mGFP+(1|mvial), data=data)
means_preds = ezPredict(fit = traitsize, boot=TRUE, iterations=10000)
ezPlot2(preds=means_preds, x = mIIS, CI = (1-0.05/3), do_plot=FALSE)

library(lmerTest)
traitsize<-lmer(wing~mIIS+mGFP+(1|mvial), data=data)
  plot(traitsize)
  dev.off()
  qqnorm(resid(traitsize))
  qqline(resid(traitsize))
  shapiro.test(resid(traitsize))
anova(traitsize)
emmeans(traitsize, list(pairwise ~ mIIS), adjust = "tukey")
detach("package:lmerTest", unload=TRUE)
traitsize<-lmer(wing~mIIS+mGFP+(1|mvial), data=data)
means_preds = ezPredict(fit = traitsize, boot=TRUE, iterations=10000)
ezPlot2(preds=means_preds, x = mIIS, CI = (1-0.05/3), do_plot=FALSE)

library(lmerTest)
traitsize<-lmer(pupa~mIIS+mGFP+(1|mvial), data=data, REML=FALSE)
  plot(traitsize)
  dev.off()
  qqnorm(resid(traitsize))
  qqline(resid(traitsize))
  shapiro.test(resid(traitsize))
anova(traitsize)
meanstable<-(lsmeans(traitsize, list(pairwise ~ mIIS), adjust = "tukey"))[[1]]
meansonly<-(summary(meanstable))[1:2]
detach("package:lmerTest", unload=TRUE)
traitsize<-lmer(pupa~mIIS+mGFP+(1|mvial), data=data)
means_preds = ezPredict(fit = traitsize, boot=TRUE, iterations=10000)
ezPlot2(preds=means_preds, x = mIIS, CI = (1-0.05/3), do_plot=FALSE)


library(lmerTest)
genitalp <- ggplot(data, aes(x=mIIS, y=genital))+
  geom_boxplot(notch=TRUE) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.1)+
  scale_x_discrete(limits=c("dfoxo.3x","uas-gfp","inr.ca"))+
  ylim(7.5,8.5)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle=90, hjust=1))


wingp <- ggplot(data, aes(x=mIIS, y=wing))+
  geom_boxplot(notch=TRUE) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.1)+
  scale_x_discrete(limits=c("dfoxo.3x","uas-gfp","inr.ca"))+
 ylim(13.8,14.8)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle=90, hjust=1))

pupap <- ggplot(data, aes(x=mIIS, y=pupa))+
  geom_boxplot(notch=TRUE) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.1)+
  scale_x_discrete(limits=c("dfoxo.3x","uas-gfp","inr.ca"))+
  ylim(13.8,14.8)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle=90, hjust=1))

grid.arrange(genitalp, wingp, pupap, ncol=3)



#Table S2: Compare copulation success with trait size

No_Competition <- read.csv("No_Competition_Only.csv")
data<-No_Competition

matetest<-glm(mate_dig~genital+wing+pupa, family=binomial(link = "logit"), data=data)
overdisp_fun(matetest)
  summary(matetest)

#Figure S1

plot_range<-range(data$genital)

newdata <- with(data, data.frame(genital = seq(plot_range[[1]], plot_range[[2]], 0.01), wing=mean(wing, na.rm=TRUE), 
                                 pupa=mean(pupa, na.rm=TRUE)))
preds <- as.data.frame(predict(matetest, newdata, type="response", se.fit=TRUE))

newdata<-cbind(newdata, preds)

std <- qnorm(0.95 / 2 + 0.5)
newdata$ymin <- preds$fit - (qnorm(0.975)*preds$se.fit) 
newdata$ymax <- preds$fit + (qnorm(0.975)*preds$se.fit)

p <- ggplot(data, aes(x=genital, y=mate_dig)) 
p + geom_jitter(height=0.02, width=0, size=1, aes(color = factor(mIIS))) +
  geom_ribbon(data=newdata, aes(y=fit, ymin=ymin, ymax=ymax), alpha=0.5) + 
  geom_line(data=newdata, aes(y=fit)) + 
  labs(x="genital", y="P(mating)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


#Table S3: Compare offspring success with trait size
offspringtest<-glm(offspring_dig~genital+wing+pupa, family=binomial(link = "logit"), data=data)
overdisp_fun(offspringtest)
summary(offspringtest)


#Figure 1B

plot_range<-range(data$genital)

newdata <- with(data, data.frame(genital = seq(plot_range[[1]], plot_range[[2]], 0.01), wing=mean(wing, na.rm=TRUE), 
                                 pupa=mean(pupa, na.rm=TRUE)))
preds <- as.data.frame(predict(offspringtest, newdata, type="response", se.fit=TRUE))

newdata<-cbind(newdata, preds)

std <- qnorm(0.95 / 2 + 0.5)
newdata$ymin <- preds$fit - (qnorm(0.975)*preds$se.fit) 
newdata$ymax <- preds$fit + (qnorm(0.975)*preds$se.fit)

p <- ggplot(data, aes(x=genital, y=offspring_dig)) 
p + geom_jitter(height=0.02, width=0, size=1, aes(color = factor(mIIS))) +
  geom_ribbon(data=newdata, aes(y=fit, ymin=ymin, ymax=ymax), alpha=0.5) + 
  geom_line(data=newdata, aes(y=fit)) + 
  labs(x="genital", y="P(offspring)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))



#Figure S2: Compare probability of copulating/producing offspring with genotype when no competition

#convert logits to probs
logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}

matetest<-glmer(mate_dig~mIIS+(1|mvial), family=binomial(link = "logit"), data=data)
overdisp_fun(matetest)
anova(matetest)
emmeans(matetest, list(pairwise ~ mIIS), adjust = "tukey")
logits<-summary(emmeans(matetest, list(pairwise ~ mIIS), adjust = "tukey"))
mate_prob<-logit2prob(logits[[1]][,2])
mate_prob
mate_lower<-(logits[[1]][,5])
mate_lowerprob<-logit2prob(mate_lower)
mate_lowerCI<-mate_prob-mate_lowerprob
mate_lowerCI
mate_upper<-(logits[[1]][,6])
mate_upperprob<-logit2prob(mate_upper)
mate_upperCI<-mate_upperprob-mate_prob
mate_upperCI

offspringtest<-glmer(offspring_dig~mIIS+(1|mvial), family=binomial(link = "logit"), data=data)
overdisp_fun(offspringtest)
summary(offspringtest)
emmeans(offspringtest, list(pairwise ~ mIIS), adjust = "tukey")
logits<-summary(emmeans(offspringtest, list(pairwise ~ mIIS), adjust = "tukey"))
off_prob<-logit2prob(logits[[1]][,2])
off_prob
off_lower<-(logits[[1]][,5])
off_lowerprob<-logit2prob(off_lower)
off_lowerCI<-off_prob-off_lowerprob
off_lowerCI
off_upper<-(logits[[1]][,6])
off_upperprob<-logit2prob(off_upper)
off_upperCI<-off_upperprob-off_prob
off_upperCI

genotype<-c(rep("dFOXO",2),rep("control",2),rep("InR",2))
success<-rep(c("copulation", "offspring"),3)
value<-c(mate_prob[1], off_prob[1], mate_prob[3], off_prob[3],mate_prob[2], off_prob[2])
lowCI<-c(mate_lowerCI[1], off_lowerCI[1], mate_lowerCI[3], off_lowerCI[3],mate_lowerCI[2], off_lowerCI[2])
highCI<-c(mate_upperCI[1], off_upperCI[1], mate_upperCI[3], off_upperCI[3],mate_upperCI[2], off_upperCI[2])
chartdata<-data.frame(genotype, success, value, lowCI, highCI)

level_order<-factor(chartdata$genotype, level=c('dFOXO', 'control', "InR"))
ggplot(chartdata, aes(fill=success, y=value, x=level_order))+
  geom_bar(width=0.5, position = position_dodge(width=0.6), stat="identity") +
  scale_x_discrete(labels = c('>dFOXO.3X','>GFP','>InR.CA'))+
  geom_errorbar(aes(ymin=value-lowCI, ymax=value+highCI), width=.2, position=position_dodge(.6))+ 
  labs(x="genotype", y="P(copulation/offspring)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

#Figure S3: Compare number of offspring by genotype when males produces offspring
data<-subset(No_Competition, mate=="Yes")
tadultstest1<-glmer(tadults~mIIS+(1|mvial),family=poisson(link="log"), data=data)
overdisp_fun(tadultstest1)

# Significant overdisperssion. Use an observation-level random effects (OLRE) model.
data$OLRE <- seq_len(nrow(data))
tadultstest<-glmer(tadults~mIIS+(1|mvial)+(1|OLRE),family=poisson(link="log"), data=data)
overdisp_fun(tadultstest)
#Confirm OLRE model is better
AIC(tadultstest1,tadultstest)
summary(tadultstest)

emmeans(tadultstest, list(pairwise ~ mIIS))
chartdata<-summary(emmeans(tadultstest, list(pairwise ~ mIIS))[[1]])
level_order<-factor(chartdata$mIIS, level=c('dfoxo.3x', 'uas-gfp', "inr.ca"))
ggplot(chartdata, aes(y=emmean, x=level_order))+
  geom_bar(width=0.5, position = position_dodge(width=0.6), stat="identity")+
  scale_x_discrete(labels = c('>dFOXO.3X','>GFP','>InR.CA'))+
  geom_errorbar(aes(ymin=asymp.LCL, ymax=asymp.UCL), width=.2, position=position_dodge(.6))+ 
  labs(x="genotype", y="log(mean number of offspring") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_fill_grey(start = 0.5, end = .9)


#Figure 1C and Table S4: Compare the effect of relative trait size on probability of producing offspring between directly competing males
Direct_Competition_Only <-  read.csv("Direct_Competition_Only.csv")
data<-Direct_Competition_Only  


matetest<-glm(offspring_dig~genitalrel+wingrel+puparel+mGFP, family=binomial(link = "logit"), data=data)
overdisp_fun(matetest)
summary(matetest)
logits<-summary(emmeans(matetest, list(pairwise ~ genitalrel), adjust = "tukey"))
logits
# Calc prob
genprob<-logit2prob(logits[[1]][,2])
genprob
#Calc CI
genlower<-(logits[[1]][,5])
genlowerprob<-logit2prob(genlower)
genlowerCI<-genprob-genlowerprob
genlowerCI
genupper<-(logits[[1]][,6])
genupperprob<-logit2prob(genupper)
genupperCI<-genupperprob-genprob
genupperCI



logits<-summary(emmeans(matetest, list(pairwise ~ wingrel), adjust = "tukey"))
logits
# Calc prob
wingprob<-logit2prob(logits[[1]][,2])
wingprob
#Calc CI
winglower<-(logits[[1]][,5])
winglowerprob<-logit2prob(winglower)
winglowerCI<-wingprob-winglowerprob
winglowerCI
wingupper<-(logits[[1]][,6])
wingupperprob<-logit2prob(wingupper)
wingupperCI<-wingupperprob-wingprob
wingupperCI

trait<-c(rep("genital",2),rep("wing",2))
size<-rep(c("larger", "smaller"),2)
value<-c(genprob[1], genprob[2], wingprob[1], wingprob[2])
lowCI<-c(genlowerprob[1],genlowerprob[2],winglowerprob[1], winglowerprob[2])
highCI<-c(genupperprob[1], genupperprob[2], wingupperprob[1], wingupperprob[2])
chartdata<-data.frame(trait, size, value, lowCI, highCI)

level_order<-factor(chartdata$trait, level=c('genital', 'wing'))
ggplot(chartdata, aes(fill=size, y=value, x=level_order))+
  geom_bar(width=0.5, position = position_dodge(width=0.6), stat="identity") +
  scale_x_discrete(labels = c('genital','wing'))+
  geom_errorbar(aes(ymin=lowCI, ymax=highCI), width=.2, position=position_dodge(.6))+ 
  labs(x="trait", y="P(offspring)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_fill_grey(start = 0.5, end = .9)



#Table S5 and Figure S4
gsize<-lmer(genital~wing+pupa+mate+(1|fid)+mGFP, data=data)
  plot(gsize)
  dev.off()
  qqnorm(resid(gsize))
  qqline(resid(gsize))
  shapiro.test(resid(gsize))
summary(gsize)
library(lmerTest)
mean_genital<-lsmeansLT(gsize)[c(1:2),1]
lowCI_genital<-lsmeansLT(gsize)[c(1:2),5]
highCI_genital<-lsmeansLT(gsize)[c(1:2),6]

wsize<-lmer(wing~mate+genital+pupa+(1|fid)+mGFP, data=data)
  plot(wsize)
  dev.off()
  qqnorm(resid(wsize))
  qqline(resid(wsize))
  shapiro.test(resid(wsize))
summary(wsize)
mean_wing<-lsmeansLT(wsize)[c(1:2),1]
lowCI_wing<-lsmeansLT(wsize)[c(1:2),5]
highCI_wing<-lsmeansLT(wsize)[c(1:2),6]

psize<-lmer(pupa~mate+genital+wing+(1|fid)+mGFP, data=data)
  plot(psize)
  dev.off()
  qqnorm(resid(psize))
  qqline(resid(psize))
  shapiro.test(resid(psize))
summary(psize)
mean_pupa<-lsmeansLT(psize)[c(1:2),1]
lowCI_pupa<-lsmeansLT(psize)[c(1:2),5]
highCI_pupa<-lsmeansLT(psize)[c(1:2),6]

trait<-c(rep("genital",2),rep("wing",2),rep("pupa",2))
size<-rep(c("don't mate", "mate"),3)
value<-c(mean_genital[1], mean_genital[2], mean_wing[1], mean_wing[2], mean_pupa[1],mean_pupa[2])
lowCI<-c(lowCI_genital[1], lowCI_genital[2], lowCI_wing[1], lowCI_wing[2], lowCI_pupa[1],lowCI_pupa[2])
highCI<-c(highCI_genital[1], highCI_genital[2], highCI_wing[1], highCI_wing[2], highCI_pupa[1],highCI_pupa[2])
chartdata<-data.frame(trait, size, value, lowCI, highCI)

level_order<-factor(chartdata$trait, level=c('genital', 'wing','pupa'))
ggplot(chartdata, aes(fill=size, y=value, x=level_order))+
  geom_bar(width=0.5, position = position_dodge(width=0.6), stat="identity") +
  scale_x_discrete(labels = c('genital','wing','pupa'))+
  geom_errorbar(aes(ymin=lowCI, ymax=highCI), width=.2, position=position_dodge(.6))+ 
  labs(x="trait", y="trait size") +
  coord_cartesian(ylim=c(5,15))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_fill_grey(start = 0.5, end = .9)



#Table S6 Figure 1D
Indirect_Competition_Only <- read.csv("Indirect_Competition_Only.csv")

data<-subset(Indirect_Competition_Only, mate=="Yes")
tadultstest<-lmer(sqrt(tadults)~relgenital+relpupal+relwing+order+(1|fid)+mGFP, data=data)
  plot(tadultstest)
  dev.off()
  qqnorm(resid(tadultstest))
  qqline(resid(tadultstest))
  shapiro.test(resid(tadultstest))
summary(tadultstest)

emmeans(tadultstest, list(pairwise ~ relgenital), adjust = "tukey")
genmeanoff<-as.data.frame(summary(emmeans(tadultstest, list(pairwise ~ relgenital), adjust = "tukey"))[[1]])
genmeanoff
emmeans(tadultstest, list(pairwise ~ relwing), adjust = "tukey")
wingmeanoff<-as.data.frame(summary(emmeans(tadultstest, list(pairwise ~ relwing), adjust = "tukey"))[[1]])
wingmeanoff
emmeans(tadultstest, list(pairwise ~ relpupal), adjust = "tukey")
pupmeanoff<-as.data.frame(summary(emmeans(tadultstest, list(pairwise ~ relpupal), adjust = "tukey"))[[1]])
pupmeanoff

trait<-c(rep("genital",2),rep("wing",2),rep("pupal",2))
size<-rep(c("larger", "smaller"),3)
value<-c((genmeanoff[2])[1,1]^2, (genmeanoff[2])[2,1]^2,(wingmeanoff[2])[1,1]^2,(wingmeanoff[2])[2,1]^2, (pupmeanoff[2])[1,1]^2,(pupmeanoff[2])[2,1]^2)
lowCI<-c((genmeanoff[5])[1,1]^2, (genmeanoff[5])[2,1]^2,(wingmeanoff[5])[1,1]^2,(wingmeanoff[5])[2,1]^2, (pupmeanoff[5])[1,1]^2,exp(pupmeanoff[5])[2,1]^2)
highCI<-c((genmeanoff[6])[1,1]^2, (genmeanoff[6])[2,1]^2,(wingmeanoff[6])[1,1]^2,(wingmeanoff[6])[2,1]^2, (pupmeanoff[6])[1,1]^2,(pupmeanoff[6])[2,1]^2)
chartdata<-data.frame(trait, size, value, lowCI, highCI)

level_order<-factor(chartdata$trait, level=c('genital', 'wing', "pupal"))
ggplot(chartdata, aes(fill=size, y=value, x=level_order))+
  geom_bar(width=0.5, position = position_dodge(width=0.6), stat="identity")+
  scale_x_discrete(labels = c('genital','wing','pupal'))+
  geom_errorbar(aes(ymin=lowCI, ymax=highCI), width=.2, position=position_dodge(.6))+ 
  labs(x="trait", y="mean number of offspring") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_fill_grey(start = 0.5, end = .9)

#Table S7
data<-Indirect_Competition_Only
gsize<-lmer(genital~order+wing+pupa+offspring+(1|fid)+mGFP, data=data)
plot(gsize)
dev.off()
qqnorm(resid(gsize))
qqline(resid(gsize))
shapiro.test(resid(gsize))
summary(gsize)
library(lmerTest)
mean_genital<-lsmeansLT(gsize)[c(1:2),1]
lowCI_genital<-lsmeansLT(gsize)[c(1:2),5]
highCI_genital<-lsmeansLT(gsize)[c(1:2),6]

wsize<-lmer(wing~order+genital+pupa+offspring+(1|fid)+mGFP, data=data)
plot(wsize)
dev.off()
qqnorm(resid(wsize))
qqline(resid(wsize))
shapiro.test(resid(wsize))
summary(wsize)
mean_wing<-lsmeansLT(wsize)[c(1:2),1]
lowCI_wing<-lsmeansLT(wsize)[c(1:2),5]
highCI_wing<-lsmeansLT(wsize)[c(1:2),6]

psize<-lmer(pupa~order+genital+wing+offspring+(1|fid)+mGFP, data=data)
plot(psize)
dev.off()
qqnorm(resid(psize))
qqline(resid(psize))
shapiro.test(resid(psize))
summary(psize)
mean_pupa<-lsmeansLT(psize)[c(1:2),1]
lowCI_pupa<-lsmeansLT(psize)[c(1:2),5]
highCI_pupa<-lsmeansLT(psize)[c(1:2),6]


#Tables S8 and S9
data<-No_Competition

clatency<-lmer(log(clatency)~mIIS+(1|mvial),  data=data)
anova(clatency)
plot(clatency)
dev.off()
qqnorm(resid(clatency))
qqline(resid(clatency))
shapiro.test(resid(clatency))
emmeans(clatency, list(pairwise ~ mIIS))

cduration<-lmer(log(cduration)~mIIS+(1|mvial),  data=data)
anova(cduration)
plot(cduration)
dev.off()
qqnorm(resid(cduration))
qqline(resid(cduration))
shapiro.test(resid(cduration))
emmeans(cduration, list(pairwise ~ mIIS))

mlatency<-lmer(log(mlatency)~mIIS+(1|mvial),  data=data)
anova(mlatency)
plot(mlatency)
dev.off()
qqnorm(resid(mlatency))
qqline(resid(mlatency))
shapiro.test(resid(mlatency))
emmeans(mlatency, list(pairwise ~ mIIS))
#Because slight deviation from normality in residuals, repeat with bootstraping
detach("package:lmerTest", unload=TRUE)
mlatency<-lmer(wing~mIIS+(1|mvial), data=data)
means_preds = ezPredict(fit = mlatency, boot=TRUE, iterations=10000)
ezPlot2(preds=means_preds, x = mIIS, CI = (1-0.05/3), do_plot=FALSE)
attach("package:lmerTest")

mduration<-lmer(log(mduration)~mIIS+(1|mvial),  data=data)
anova(mduration)
plot(mduration)
dev.off()
qqnorm(resid(mduration))
qqline(resid(mduration))
shapiro.test(resid(mduration))
emmeans(mduration, list(pairwise ~ mIIS))

clatency<-lm(log(clatency)~genital+wing+pupa,  data=data)
Anova(clatency,type="III")
plot(clatency)
shapiro.test(resid(clatency))

cduration<-lm(log(cduration)~genital+wing+pupa,data=data)
Anova(cduration,type="III")
plot(cduration)
shapiro.test(resid(cduration))

mlatency<-lm(log(mlatency)~genital+wing+pupa,data=data)
Anova(mlatency,type="III")
plot(mlatency)
shapiro.test(resid(mlatency))

mduration<-lm(log(mduration)~genital+wing+pupa,data=data)
Anova(mduration,type="III")
plot(mlatency)
shapiro.test(resid(mduration))

#Figure S6
mlatency<-lmer(log(mlatency)~mIIS+(1|mvial),  data=data)
emmeans(mlatency, list(pairwise ~ mIIS))
chartdata<-summary(emmeans(mlatency, list(pairwise ~ mIIS))[[1]])
level_order<-factor(chartdata$mIIS, level=c('dfoxo.3x', 'uas-gfp', "inr.ca"))
ggplot(chartdata, aes(y=emmean, x=level_order))+
  geom_bar(width=0.5, position = position_dodge(width=0.6), stat="identity")+
  scale_x_discrete(labels = c('>dFOXO.3X','>GFP','>InR.CA'))+
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.2, position=position_dodge(.6))+ 
  labs(x="genotype", y="log(mating latency)(sec)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_fill_grey(start = 0.5, end = .9)
