###   SUPPLEMENTAL ANALYSES
###   Authors: M.C. Sabal, J.E. Merz, S.H. Alonzo, E.P. Palkovacs

##################################################################


# load required packages
library(ggplot2); library(coxme); library(lme4); library(gridExtra); library(car)


# load data for behavioral assay analysis - 63 salmon x 2 trials x 2 splits = 252
assay <- read.table("sabal_et_al_2020_behavioral_assay_data.txt", header = TRUE, sep = "")


#### APPENDIX S1 ####
## Additional analyses on physical traits of salmon groups (hatchery, wild-upstream, wild-downstream).

# Do salmon traits vary by salmon group?

# Get dataframe only of individual salmon and their physical traits
sam <- unique(assay[,c("date", "fish_no", "salmon_group", "FL", "WT", "Kbcond", "ATPresid")])

# Check if physical traits correlated. They are  not.
cordat<-sam[,c("FL", "Kbcond", "ATPresid")]
cordat<-na.omit(cordat)
cor(cordat)

# MANOVA - Do multiple physical traits vary by group? -Yes.
summary(manova(cbind(FL, ATPresid, Kbcond) ~ salmon_group, data=sam))

# Subsequent ANOVAs on each trait
summary(aov(Kbcond ~ salmon_group, data=sam)); TukeyHSD(aov(Kbcond ~ salmon_group, data=sam))
summary(aov(FL ~ salmon_group, data=sam)); TukeyHSD(aov(FL ~ salmon_group, data=sam))
summary(aov(ATPresid ~ salmon_group, data=sam)); TukeyHSD(aov(ATPresid ~ salmon_group, data=sam))


# Do salmon traits affect salmon's reaction to predation risk? - Table S1
# Mixed-effects cox model: look at only interactions of physical traits with predator treatment
coxme.traits<-coxme(Surv((assay$time.cox), assay$censor) ~ ATPresid*pred_treat + FL*pred_treat + Kbcond*pred_treat + (1|salmon_group/fish_no), data=assay)
Anova(coxme.traits, type="III") #Supplemental Table S1



#### APPENDIX S2 ####

## Calculate speed as body lengths per second (bl/s)
# split length is 0.9144 meters convereted to mm
assay$split.length.mm <- 0.9144 * 1000

# split length in body lengths for individual salmon
assay$split.length.bl <- assay$split.length.mm / assay$FL

# speed in bl/s
assay$speed.bls <- assay$split.length.bl / assay$time.secs


## Linear mixed-effects model on truncated dataset:
## does salmon speed vary by predator treatment, split, or salmon group?

### M/S ###
## Mixed-effects cox model (speed m/s)
lmer_mod_ms<-lmer(log(speed.ms) ~ salmon_group + pred_treat + salmon_group:pred_treat + split + salmon_group:split + split:pred_treat + (1|fish_no), data=assay)

## ANOVA
Anova(lmer_mod_ms, type="III") #Supplemental Table S2

## Linear contrasts - Table S3
# contrasts: for each salmon group, does predator treatment affect speed (m/s)?
lsmeans(lmer_mod_ms, pairwise ~ pred_treat | salmon_group)

### BL/S ###
## Mixed-effects cox model (speed m/s)
lmer_mod_bls<-lmer(log(speed.bls) ~ salmon_group + pred_treat + salmon_group:pred_treat + split + salmon_group:split + split:pred_treat + (1|fish_no), data=assay)

## ANOVA
Anova(lmer_mod_bls, type="III") #Supplemental Table S2

## Linear contrasts - Table S3
# contrasts: for each salmon group, does predator treatment affect speed (bl/s)?
lsmeans(lmer_mod_bls, pairwise ~ pred_treat | salmon_group)


#### APPENDIX S3 ####
# contrasts: for each salmon group, does speed (m/s) differ? - Table S4
lsmeans(lmer_mod_ms, pairwise ~ salmon_group | split)

# contrasts: for each salmon group, does speed (bl/s) differ? - Table S4
lsmeans(lmer_mod_bls, pairwise ~ salmon_group | split)


#### FIGURES ####
#### Supplemental Figure S2
options(scipen=999)

#change order of factor levels from low predator experience to high
assay$salmon_group<-factor(assay$salmon_group, levels=c("hatchery", "wild-upstream", "wild-downstream"))


#Body Condition plot
Kbcond.plot<-ggplot(data=assay, aes(x=salmon_group, y=Kbcond, fill=salmon_group)) + geom_boxplot(size=0.5, width=0.8, outlier.shape=NA) +
  scale_fill_manual(values=c("firebrick3", "steelblue1", "royalblue")) + theme_bw() +
  theme(axis.title.x = element_blank()) + theme(legend.position = "bottom") +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=12), 
        legend.text=element_text(size=12), strip.text = element_text(size=12)) +
  theme(legend.title=element_blank()) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position=c(.85, .85)) +  theme(legend.position="none") + ylab("Condition factor (K)") +
  annotate("text", x=c(1, 2, 3), y=c(0.00142, 0.00142, 0.00142), label = c("a", "b", "b"), size=4) +
  annotate("text", x=0.7, y=0.0014, label="A", fontface="bold", size=6) +
  theme(axis.line = element_line(colour = "black"),panel.border = element_blank(), panel.background = element_blank()) 
Kbcond.plot

#ATPase plot
atp.plot<-ggplot(data=assay, aes(x=salmon_group, y=ATPresid, fill=salmon_group)) + geom_boxplot(size=0.5, width=0.8, outlier.shape=NA) +
  scale_fill_manual(values=c("firebrick3", "steelblue1", "royalblue")) + theme_bw() +
  theme(axis.title.x = element_blank()) + theme(legend.position = "bottom") +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=12), 
        legend.text=element_text(size=12), strip.text = element_text(size=12)) +
  theme(legend.title=element_blank()) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position="none") + ylab("ATPase residual") + ylim(c(-2, 1.4)) +
  annotate("text", x=c(1, 2, 3), y=c(1.34, 1.34, 1.34), label = c("a", "a", "b"), size=4) +
  annotate("text", x=0.65, y=1.2, label="B", fontface="bold", size=6) +
  theme(axis.line = element_line(colour = "black"),panel.border = element_blank(), panel.background = element_blank()) 
atp.plot

wt.ln.plot<-ggplot(data=assay, aes(x=FL, y=WT)) + geom_point(pch=21, size=2, aes(fill=salmon_group)) + theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + scale_fill_manual(values=c("firebrick3", "steelblue1", "royalblue")) +
  theme(legend.title=element_blank()) + theme(legend.position=c(.7,.2)) +
  ylab("Weight (g)") + xlab("Fork length (mm)") +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=12), 
        legend.text=element_text(size=10)) +  theme(legend.position="none") +
  geom_smooth(method='loess', color="black", alpha=0.2, size=0.5) +
  annotate("text", x=60, y=15.5, label="C", fontface="bold", size=6) +
  theme(legend.position=c(.35, .8), legend.background = element_rect(colour = NA)) +
  theme(axis.line = element_line(colour = "black"),panel.border = element_blank(), panel.background = element_blank()) +
  annotate("text", x=95, y=3, label="all pair-wise  groups\\nsignificantly different", size=3.5)
wt.ln.plot

# Supplemental Figure S2
grid.arrange(Kbcond.plot, atp.plot, wt.ln.plot, ncol=2)


# Supplemental Figure S3
# make reaction norm
fmt_dcimals <- function(decimals=0){ function(x) as.character(round(x,decimals)) } # set to two decimal places

rxndat<-assay[,c("fish_no", "speed.bls","split", "salmon_group", "pred_treat")]
rxndat<-unique(rxndat)

### speed by salmon group + split + pred_treat
aggr.dat<-aggregate(speed.bls ~ salmon_group + pred_treat + split, data=rxndat, mean)
a<-aggregate(speed.bls ~ salmon_group + pred_treat + split, data=rxndat, sd); colnames(a)[4]<-"sd"
b<-aggregate(speed.bls ~ salmon_group + pred_treat + split, data=rxndat, length); colnames(b)[4]<-"n"
aggr.dat<-join(aggr.dat, a); aggr.dat<-join(aggr.dat, b)
aggr.dat$se<-aggr.dat$sd / sqrt(aggr.dat$n)

aggr.dat$pred_treat<-as.factor(aggr.dat$pred_treat)
aggr.dat$split<-as.factor(aggr.dat$split)

levels(aggr.dat$pred_treat)<-c("no predator", "predator")
levels(aggr.dat$split)<-c("A1 to A2", "A2 to A3")

figS3<-ggplot(data=subset(aggr.dat, split != "A1 to A3"), aes(x=pred_treat, y=speed.bls, fill=salmon_group)) +
  geom_rect(xmin = -Inf,xmax = 1.5, ymin = -Inf,ymax = Inf,alpha = 0.3, fill="gray92") +
  geom_line(aes(group=salmon_group, color=salmon_group), size=0.5, position=position_dodge(.3)) +
  geom_errorbar(aes(ymax=speed.bls + se, ymin=speed.bls - se, color=salmon_group), width=0, position=position_dodge(.3), size=0.5) +
  geom_point(pch=21, size=1.5, position=position_dodge(.3)) + facet_wrap(~split, ncol=2) +
  scale_fill_manual(values=c("firebrick3", "royalblue", "steelblue1")) + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x=element_blank()) +
  theme(strip.text.x = element_text(face="bold")) +
  theme(strip.background = element_rect(color="black", fill="white")) +
  scale_color_manual(values=c("firebrick3","royalblue", "steelblue1")) +
  ylab("Speed (bl/s)") + theme(legend.position = "bottom") +
  theme(legend.title = element_blank(), strip.background=element_blank()) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=12), 
        legend.text=element_text(size=10), strip.text = element_text(size=12)) +
  theme(legend.position=c(0.88,0.84), legend.background = element_rect(fill="transparent")) +
  ylim(c(0,1.3)) + theme(axis.line = element_line(colour = "black"), panel.border = element_blank(), panel.background = element_blank())

figS3
