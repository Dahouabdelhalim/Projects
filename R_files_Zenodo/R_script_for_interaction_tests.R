# R analysis for the voluntary Interaction tests for the contrast study data - run two separate analyses one that looks at
# the time points 1, 5 and 9 and one that looks at 19 and 27. 

# specify the data 
InteractionTests$Treatment <- as.factor(InteractionTests$Treatment)
InteractionTests$Day <-as.factor(InteractionTests$Day)
InteractionTests$PrePost <- as.factor(InteractionTests$PrePost)
InteractionTests$Cage <-as.factor(InteractionTests$Cage)


# next thing to do is to look at calculating the averages so then we can see what the data looks like.

library(plyr)

# by treatment 
int_treat <- ddply(InteractionTests, "Treatment", summarise, av=mean(VoluntaryInt), 
                   N = length(VoluntaryInt), SD = sd(VoluntaryInt), SE   = SD / sqrt(N))

# by day 
int_day <-ddply(InteractionTests, "Day", summarise, av=mean(VoluntaryInt), 
                N= length(VoluntaryInt), SD=sd(VoluntaryInt), SE = SD/sqrt(N))

# by pre post (time)
int_time <-ddply(InteractionTests, "PrePost", summarise, av=mean(VoluntaryInt), 
                N= length(VoluntaryInt), SD=sd(VoluntaryInt), SE = SD/sqrt(N))

# by time, treatment and prepost

int_treat_day_time <-ddply(InteractionTests, c("Treatment", "Day", "PrePost"), summarise, av=mean(VoluntaryInt),
                           N = length(VoluntaryInt), SD=sd(VoluntaryInt), SE=SD/sqrt(N))

# next thing to do is to plot these data so can see what they look like. 

library(ggplot2)

#plotting all factors - essentially this is what we want to look at. 
ggplot(InteractionTests, aes(x=Treatment, y=VoluntaryInt)) +geom_bar(stat = 'identity') +
  theme_bw() +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlab("Treatment") + ylab("Voluntary Interaction" ) + facet_grid(Day ~ PrePost)

# lets have a look at the raw data to see how it is distributed 

# need to split the data so we look at the days that were conducted in line with the handling phase (1,5 and 9) and then days that
# were conducted during the drinking experiments

HandlingPhaseIntdata <- subset(InteractionTests , Day%in%c(1,5,9))
DrinkingPhaseIntdata <- subset(InteractionTests, Day%in%c(19, 27))

HandlingPhaseIntdata <- droplevels(HandlingPhaseIntdata)
DrinkingPhaseIntdata <- droplevels(DrinkingPhaseIntdata)

# first will do the analysis for the interaction tests conducted in line with the handling phase 
library(ggpubr)

ggdensity(HandlingPhaseIntdata$VoluntaryInt, main = "density plot of Voluntary Interaction", xlab = "Voluntary Interaction")
ggqqplot(HandlingPhaseIntdata$VoluntaryInt)# definitely the data are non-normal when look at this. Relatively bimodal in
# terms of ditibution, looks like we have a peak at 0 and then one slightly higher up too. Therefore instead of trying to 
# transform this data I think it makes more sense to run a Generalised Linear Mixed Model (GLMM) where you have voluntary
# interaction as the dependent varaible, treatment day and time as  
# fixed factors and cage number as the random factor. I think fitting a gamma distribution might be best in this case. 

library(lme4)
library(optimx)

m1gam =glmer((VoluntaryInt) ~ Treatment + Day + PrePost + Treatment:Day + Treatment:PrePost + Day:PrePost + Treatment:Day: PrePost + (1|Cage), data=HandlingPhaseIntdata, 
             family= Gamma(link= "log"), REML= FALSE, control = glmerControl(optimizer = 'optimx', optCtrl = list(method='nlminb')))
# this uses a different optimiser to the bobyqa one as this still failed to converge this uses the optimx optimiser and the method nlminb 
m1gam
summary(m1gam)
# lets look at the fit of the model by looking at our residuals - seem like a pretty good fit. 

hist(resid(m1gam))
plot(m1gam)

# drop the three way interaction 

m1.1 =glmer((VoluntaryInt) ~ Treatment + Day + PrePost + Treatment:Day + Treatment:PrePost + Day:PrePost + (1|Cage), data=HandlingPhaseIntdata, 
             family= Gamma(link= "log"), REML= FALSE, control = glmerControl(optimizer = 'optimx', optCtrl = list(method='nlminb')))

anova(m1gam, m1.1) # looks like we have a significant three way interaction. X2=29.422, p=4.084e-7

# drop the treatment x day interaction 

m1.2 =glmer((VoluntaryInt) ~ Treatment + Day + PrePost + Treatment:PrePost + Day:PrePost + (1|Cage), data=HandlingPhaseIntdata, 
             family= Gamma(link= "log"), REML= FALSE, control = glmerControl(optimizer = 'optimx', optCtrl = list(method='nlminb')))

anova(m1gam, m1.2) # this shows a significant treatment x day interaction. X2=47.25, p=1.353e-9

# drop the treatment x prepost interaction 
m1.3 =glmer((VoluntaryInt) ~ Treatment + Day + PrePost + Treatment:Day + Day:PrePost + (1|Cage), data=HandlingPhaseIntdata, 
            family= Gamma(link= "log"), REML= FALSE, control = glmerControl(optimizer = 'optimx', optCtrl = list(method='nlminb')))

anova(m1gam, m1.3) # this shows a significant treatment x prepost interaction x2=29.664, p=1.624e-6

# drop the day x prepost interaction

m1.4 =glmer((VoluntaryInt) ~ Treatment + Day + PrePost + Treatment:Day + Treatment:PrePost + (1|Cage), data=HandlingPhaseIntdata, 
            family= Gamma(link= "log"), REML= FALSE, control = glmerControl(optimizer = 'optimx', optCtrl = list(method='nlminb')))

anova(m1gam, m1.4) # significant day x prepost interaction x2=49.886, p=3.813e-10

# next create a model with just the main effects in and then take each one out in turn to determine its effect. 

m1.5 =glmer((VoluntaryInt) ~ Treatment + Day + PrePost + (1|Cage), data=HandlingPhaseIntdata, 
            family= Gamma(link= "log"), REML= FALSE, control = glmerControl(optimizer = 'optimx', optCtrl = list(method='nlminb')))

# drop treatment 
m1.6 =glmer((VoluntaryInt) ~ Day + PrePost + (1|Cage), data=HandlingPhaseIntdata, 
            family= Gamma(link= "log"), REML= FALSE, control = glmerControl(optimizer = 'optimx', optCtrl = list(method='nlminb')))

anova(m1.5, m1.6) # main effect of treatment x2=62.132, p=3.21e-15

# drop day

m1.7 =glmer((VoluntaryInt) ~ Treatment + PrePost + (1|Cage), data=HandlingPhaseIntdata, 
            family= Gamma(link= "log"), REML= FALSE, control = glmerControl(optimizer = 'optimx', optCtrl = list(method='nlminb')))

anova(m1.5, m1.7) # main effect of day x2=25.306, p=3.198e-6

# drop pre post time 

m1.8 =glmer((VoluntaryInt) ~ Treatment + Day + (1|Cage), data=HandlingPhaseIntdata, 
            family= Gamma(link= "log"), REML= FALSE, control = glmerControl(optimizer = 'optimx', optCtrl = list(method='nlminb')))

anova(m1.5, m1.8) # non significant main effect of pre post time x2=1.386, p=0.2391

# so what I need to do next is look at the pairwise comparisons for the significant interactions... 
library(multcomp)

# lets look at the significant three way int (treatment, day, prepost)
HandlingPhaseIntdata$int1 <- with(HandlingPhaseIntdata, interaction(Treatment, Day, PrePost, sep = "x"))
m1int1 = glmer((VoluntaryInt) ~ int1 + 
                (1 |Cage), data=HandlingPhaseIntdata, family= Gamma(link= "log"), REML= FALSE, control = glmerControl(optimizer = 'optimx', optCtrl = list(method='nlminb')))

summary(glht(m1int1, mcp(int1="Tukey"))) # this shows that for the pre data for tail - decrease after day 1 post then changes
# across the post handling interaction across days. No differences in the tunnel for when comparing pre across days or post across days
# this seems to differ to what I had previously. 

# next lets look at the significant treatment x day interaction
HandlingPhaseIntdata$int2 <- with(HandlingPhaseIntdata, interaction(Treatment, Day, sep = "x"))
m1int2 = glmer((VoluntaryInt +1) ~ int2 + 
                 (1 |Cage), data=HandlingPhaseIntdata, family= Gamma(link= "log"), REML= FALSE, control = glmerControl(optimizer = 'optimx', optCtrl = list(method='nlminb')))

summary(glht(m1int2, mcp(int2="Tukey"))) # this shows that tunnel are always higher at each day than tail. No change across time with
# tunnel handled mice. Tail handled mice differ between all day combinations. 

# next lets look at the significant treatment x prepost interaction 

HandlingPhaseIntdata$int3 <- with(HandlingPhaseIntdata, interaction(Treatment, PrePost, sep = "x"))
m1int3 = glmer((VoluntaryInt +1) ~ int3 + 
                 (1 |Cage), data=HandlingPhaseIntdata, family= Gamma(link= "log"), REML= FALSE, control = glmerControl(optimizer = 'optimx', optCtrl = list(method='nlminb')))

summary(glht(m1int3, mcp(int3="Tukey"))) # this shows that for tail handled mice there is a differences between pre and post
# handling interaction tests whereas this is not the case for tunnel handled mice. 

# next lets look at the significant day x prepost interaction

HandlingPhaseIntdata$int4 <- with(HandlingPhaseIntdata, interaction(Day, PrePost, sep = "x"))
m1int4 = glmer((VoluntaryInt +1) ~ int4 + 
                 (1 |Cage), data=HandlingPhaseIntdata, family= Gamma(link= "log"), REML= FALSE, control = glmerControl(optimizer = 'optimx', optCtrl = list(method='nlminb')))

summary(glht(m1int4, mcp(int4="Tukey"))) # this shows that when comparing across post sessions (irrespective of treatment)
# that days 1 and 5 and 1 and 9 differ but the pre sessions dont differ. 

# we also have main effescts of treatment - tunnel handled higher overall but we also have a main effect of day so lets look at this

summary(glht(m1.5, mcp(Day="Tukey"))) # this shows that day 1 is signficantly different to both 5 and 9 but no difference
# between 5 and 9. 



# next lets do the same again looking at the interaction tests during the sucrose drinking experiments (days 19 and 27).

# lets look at the distrubution of the raw data

ggdensity(DrinkingPhaseIntdata$VoluntaryInt, main = "density plot of Voluntary Interaction", xlab = "Voluntary Interaction")
ggqqplot(DrinkingPhaseIntdata$VoluntaryInt) # very similar to the plot for the other days - data non-normal and looks like 
# we have a bimodal distiburion. Therefore rather than transformation will run a Generalised Linear mixed model as above with gamma distribution

m2gam =glmer((VoluntaryInt +1) ~ Treatment + Day + PrePost + Treatment:Day + Treatment:PrePost + Day:PrePost + Treatment:Day: PrePost + (1|Cage), data=DrinkingPhaseIntdata, 
             family= Gamma(link= "log"), REML= FALSE, control = glmerControl(optimizer = 'optimx', optCtrl = list(method='L-BFGS-B'))) # this uses a different
# optimiser need to have optimx loaded for this - the standard of using bobyqa didnt help the model converge whereas this
# optimiser seems to work fine and get no error messages. 

summary(m2gam)


# lets look at the fit of the model by looking at our residuals - doesnt seem like a bad fit. 

hist(resid(m2gam))
plot(m2gam)

# drop the three way interaction 

m2.1 =glmer((VoluntaryInt +1) ~ Treatment + Day + PrePost + Treatment:Day + Treatment:PrePost + Day:PrePost + (1|Cage), data=DrinkingPhaseIntdata, 
             family= Gamma(link= "log"), REML= FALSE, control = glmerControl(optimizer = 'optimx', optCtrl = list(method='L-BFGS-B')))
anova(m2gam, m2.1) # non-significant interaction remove this altogether from the model for below comparisons...x2=0.30, p=0.5835

# drop the treatment x day interaction 

m2.2 =glmer((VoluntaryInt +1) ~ Treatment + Day + PrePost + Treatment:PrePost + Day:PrePost + (1|Cage), data=DrinkingPhaseIntdata, 
            family= Gamma(link= "log"), REML= FALSE, control = glmerControl(optimizer = 'optimx', optCtrl = list(method='L-BFGS-B'))) # still fails to converge

anova(m2.1, m2.2) # this shows a non-significant treatment x day interaction. X2=1.06, p=0.3039

# drop the day x prepost interaction
m2.3 =glmer((VoluntaryInt +1) ~ Treatment + Day + PrePost + Treatment:Day + Treatment:PrePost + (1|Cage), data=DrinkingPhaseIntdata, 
            family= Gamma(link= "log"), REML= FALSE, control = glmerControl(optimizer = 'optimx', optCtrl = list(method='L-BFGS-B')))

anova(m2.1, m2.3) # non-significant day x prepost interaction x2=0.7247, p=0.3946

# drop the treatment x prepost interaction 
m2.4 =glmer((VoluntaryInt +1) ~ Treatment + Day + PrePost + Treatment:Day + Day:PrePost + (1|Cage), data=DrinkingPhaseIntdata, 
            family= Gamma(link= "log"), REML= FALSE, control = glmerControl(optimizer = 'optimx', optCtrl = list(method='L-BFGS-B')))

anova(m2.1, m2.4) # significant treatment x prepost effect x2=10.495, p=0.0012

# Just have model with main effects in 
m2.5 =glmer((VoluntaryInt +1) ~ Treatment + Day + PrePost + (1|Cage), data=DrinkingPhaseIntdata, 
            family= Gamma(link= "log"), REML= FALSE, control = glmerControl(optimizer = 'optimx', optCtrl = list(method='L-BFGS-B')))

# then remove treatment 
m2.6 =glmer((VoluntaryInt +1) ~ Day + PrePost + (1|Cage), data=DrinkingPhaseIntdata, 
            family= Gamma(link= "log"), REML= FALSE, control = glmerControl(optimizer = 'optimx', optCtrl = list(method='L-BFGS-B')))

anova(m2.5, m2.6) # this shows a main effect of treatment x2=52.199, p=5.014e-13

# drop day
m2.7 =glmer((VoluntaryInt +1) ~ Treatment + PrePost + (1|Cage), data=DrinkingPhaseIntdata, 
            family= Gamma(link= "log"), REML= FALSE, control = glmerControl(optimizer = 'optimx', optCtrl = list(method='L-BFGS-B')))

anova(m2.5, m2.7) # no main effect of day x2=1.3155, p=0.2514

# drop pre post time 
m2.8 =glmer((VoluntaryInt +1) ~ Treatment + Day + (1|Cage), data=DrinkingPhaseIntdata, 
            family= Gamma(link= "log"), REML= FALSE, control = glmerControl(optimizer = 'optimx', optCtrl = list(method='L-BFGS-B')))

anova(m2.5, m2.8) # main effect of prepost (just p=0.047) - x2=3.9416, p=0.047

# significant treatment x prepost interaction 
DrinkingPhaseIntdata$int1 <- with(DrinkingPhaseIntdata, interaction(Treatment, PrePost, sep = "x"))
m2int1 = glmer((VoluntaryInt +1) ~ int1 + 
                 (1 |Cage), data=DrinkingPhaseIntdata, family= Gamma(link= "log"), REML= FALSE, control = glmerControl(optimizer = 'optimx', optCtrl = list(method='L-BFGS-B')))
summary(glht(m2int1, mcp(int1="Tukey")))

# the next thing to do is look at where the signficant difference lie in terms of the significant main effects.
library(multcomp)

summary(glht(m2.5, mcp(Treatment="Tukey"))) # this basically just shows that overall tunnel handled mice have larger 
# lick cluster sizes. 

summary(glht(m2.5, mcp(PrePost="Tukey"))) # this just shows that overall pre and post are different. 
