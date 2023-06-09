###   BEHAVIORAL ASSAY ANALYSIS
###   HOW DOES PREDATION RISK AFFECT JUVENILE SALMON TRAVEL SPEED?
###   Authors: M.C. Sabal, J.E. Merz, S.H. Alonzo, E.P. Palkovacs

##################################################################


# load required packages
library(coxme); library(car); library(lsmeans); library(effsize); library(plyr); 
library(ggplot2); library(chron)


# load data for behavioral assay analysis - 63 salmon x 2 trials x 2 splits = 252
assay <- read.table("sabal_et_al_2020_behavioral_assay_data.txt", header = TRUE, sep = "")
assay$time_start <-times(assay$time_start) #make times into time format.

envi_dat<-read.table("sabal_et_al_2020_environmental_data.txt", header = TRUE, sep ="")


#### ANALYSIS ####
## mixed-effects cox model
coxme_mod<-coxme(Surv((assay$time.secs), assay$censor) ~ salmon_group + pred_treat + salmon_group:pred_treat + split + salmon_group:split + split:pred_treat + (1|fish_no), data=assay)

## ANOVA - Table #1
Anova(coxme_mod, type="III")

## Linear contrasts - Table #2
# contrasts: for each salmon group, does predator treatment affect speed?
lsmeans(coxme_mod, pairwise ~ pred_treat | salmon_group)

# contrasts: for each salmon group, does split affect speed?
#lsmeans(coxme_mod, pairwise ~ split | salmon_group)

# contrasts: for each salmon group, does speed differ?
#lsmeans(coxme_mod, pairwise ~ salmon_group)



## Effect sizes - Table #2
#salmon group and predator treatment
hatch_pred<-assay[assay$salmon_group == "hatchery" & assay$pred_treat == "pred", c("speed.ms")]
hatch_nopred<-assay[assay$salmon_group == "hatchery" & assay$pred_treat == "no_pred", c("speed.ms")]
cohen.d(hatch_pred, hatch_nopred, paired=FALSE, na.rm=TRUE, hedges.correction = TRUE)

wildup_pred<-assay[assay$salmon_group == "wild-upstream" & assay$pred_treat == "pred", c("speed.ms")]
wildup_nopred<-assay[assay$salmon_group == "wild-upstream" & assay$pred_treat == "no_pred", c("speed.ms")]
cohen.d(wildup_pred, wildup_nopred, paired=FALSE, na.rm=TRUE, hedges.correction = TRUE)

wilddown_pred<-assay[assay$salmon_group == "wild-downstream" & assay$pred_treat == "pred", c("speed.ms")]
wilddown_nopred<-assay[assay$salmon_group == "wild-downstream" & assay$pred_treat == "no_pred", c("speed.ms")]
cohen.d(wilddown_pred, wilddown_nopred, paired=FALSE, na.rm=TRUE, hedges.correction = TRUE)

#salmon group and split
#hatch12<-assay[assay$salmon_group == "hatchery" & assay$split == "A1_A2", c("speed.ms")]
#hatch23<-assay[assay$salmon_group == "hatchery" & assay$split == "A2_A3", c("speed.ms")]
#cohen.d(hatch12, hatch23, paired=FALSE, na.rm=TRUE, hedges.correction = TRUE)

#wildup12<-assay[assay$salmon_group == "wild-upstream" & assay$split == "A1_A2", c("speed.ms")]
#wildup23<-assay[assay$salmon_group == "wild-upstream" & assay$split == "A2_A3", c("speed.ms")]
#cohen.d(wildup12, wildup23, paired=FALSE, na.rm=TRUE, hedges.correction = TRUE)

#wilddown12<-assay[assay$salmon_group == "wild-downstream" & assay$split == "A1_A2", c("speed.ms")]
#wilddown23<-assay[assay$salmon_group == "wild-downstream" & assay$split == "A2_A3", c("speed.ms")]
#cohen.d(wilddown12, wilddown23, paired=FALSE, na.rm=TRUE, hedges.correction = TRUE)



## Non-target variables on reaction of salmon to predation risk.
# In the tests below, we are interested if there is a significant 
# interaction between the target variable and predator treatment on speed.
# In this case, we are not interested in main effects as these are 
# captured in the main analysis.

# time of day
timeofday_mod<-lm(log(speed.ms) ~ time_start*pred_treat, data=assay)
summary(timeofday_mod)

# order
order_mod<-lm(log(speed.ms) ~ order*pred_treat, data=assay)
summary(order_mod)

# time between trials
timebtw_fun<-function(x){
  out<-max(x$time_start) - min(x$time_start)
  out}  # custom function to calculate the time between predator and no predator trials for an individual salmon.

# apply function over dataframe to get times for each salmon.
out_timebtw<-ddply(assay, .(fish_no), timebtw_fun); colnames(out_timebtw)[2]<-"time_btw_trials"

# add back into main dataframe
assay<-join(out_timebtw, assay)

# analysis
timebtw_mod<-lm(log(speed.ms) ~ time_btw_trials*pred_treat, data=assay)
summary(timebtw_mod)

## Summarize environmental data
mean(envi_dat$airtempC)
sd(envi_dat$airtempC)

mean(envi_dat$watertempC)
sd(envi_dat$watertempC)

mean(envi_dat$turbidity, na.rm=T)
sd(envi_dat$turbidity, na.rm=T)

mean(envi_dat$watervelocity)
sd(envi_dat$watervelocity)

## Figure #3
# make reaction norm
fmt_dcimals <- function(decimals=0){ function(x) as.character(round(x,decimals)) } # set to two decimal places

rxndat<-assay[,c("fish_no", "speed.ms","split", "salmon_group", "pred_treat")]
rxndat<-unique(rxndat)

### speed by salmon group + split + pred_treat
aggr.dat<-aggregate(speed.ms ~ salmon_group + pred_treat + split, data=rxndat, mean)
a<-aggregate(speed.ms ~ salmon_group + pred_treat + split, data=rxndat, sd); colnames(a)[4]<-"sd"
b<-aggregate(speed.ms ~ salmon_group + pred_treat + split, data=rxndat, length); colnames(b)[4]<-"n"
aggr.dat<-join(aggr.dat, a); aggr.dat<-join(aggr.dat, b)
aggr.dat$se<-aggr.dat$sd / sqrt(aggr.dat$n)

aggr.dat$pred_treat<-as.factor(aggr.dat$pred_treat)
aggr.dat$split<-as.factor(aggr.dat$split)

levels(aggr.dat$pred_treat)<-c("no predator", "predator")
levels(aggr.dat$split)<-c("A1 to A2", "A2 to A3")

fig3.final<-ggplot(data=subset(aggr.dat, split != "A1 to A3"), aes(x=pred_treat, y=speed.ms, fill=salmon_group)) +
  geom_rect(xmin = -Inf,xmax = 1.5, ymin = -Inf,ymax = Inf,alpha = 0.3, fill="gray92") +
  geom_line(aes(group=salmon_group, color=salmon_group), size=0.5, position=position_dodge(.3)) +
  geom_errorbar(aes(ymax=speed.ms + se, ymin=speed.ms - se, color=salmon_group), width=0, position=position_dodge(.3), size=0.5) +
  geom_point(pch=21, size=1.5, position=position_dodge(.3)) + facet_wrap(~split, ncol=2) +
  scale_fill_manual(values=c("firebrick3", "royalblue", "steelblue1")) + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x=element_blank()) +
  theme(strip.text.x = element_text(face="bold")) +
  theme(strip.background = element_rect(color="black", fill="white")) +
  scale_color_manual(values=c("firebrick3","royalblue", "steelblue1")) +
  ylab("Speed (m/s)") + theme(legend.position = "bottom") +
  theme(legend.title = element_blank(), strip.background=element_blank()) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=12), 
        legend.text=element_text(size=10), strip.text = element_text(size=12)) +
  theme(legend.position=c(0.88,0.84), legend.background = element_rect(fill="transparent")) +
  ylim(c(0,0.115)) + theme(axis.line = element_line(colour = "black"), panel.border = element_blank(), panel.background = element_blank())

fig3.final
