# Read in data files.  (Files have same data in a long and a wide format).  
dat.long<-read.csv("long.csv")
dat.wide<-read.csv("wide.csv")

#get average time between first and second tests for each female
date.first<-as.Date(as.character(dat.wide$date_1), format = "%Y%m%d")
date.second<-as.Date(as.character(dat.wide$date_2), format = "%Y%m%d")
days.btw<-as.numeric(date.second - date.first)
dat.wide<-cbind.data.frame(dat.wide, days.btw)
min(dat.wide$days.btw) #min 13 days
mean(dat.wide$days.btw) #mean 27
sd(dat.wide$days.btw) #SD 18

#get condition of each female using Scaled Mass Index
library(smatr)
SMI = function(x, y, x.0 = mean(x, na.rm = T)) {
  B.sma <- coef(sma(log(y) ~ log(x), robust = F, na.action = na.omit))[2]
  result = y * (x.0 / x) ^ B.sma
  return(result)
}

dat.long$smi<-SMI(x = dat.long$svl, y = dat.long$mass)

#how many females were tested in both conditions?
test<-dat.wide[, c(3,9,13,21,24)]
test 
#2 females tested only in pure SM condition, 3 tested only in mixed chorus condition
#48 females tested in both conditions.



#ANALYZING TESTS IN WHICH FEMALES MADE CHOICES
#summarize choices and apply exact binomial tests:
#in the mixed chorus treatment
summary(dat.long[which(dat.long$chorus == "mixed"),13]) #15F, 29S
binom.test(x=15, n = (15+29)) #prob choosing fast = 0.34, sig diff from 50% at p = 0.049
#in the pureSM chorus treatment
summary(dat.long[which(dat.long$chorus == "pureSM"),13]) #22F, #18S
binom.test(x=22, n = (18+22)) #prob choosing fast = 0.55, not sig diff from 0 at p = 0.636

#Run glms to test whether chorus treatment predicts choices of F vs. S for females that made choices
#removing 'No Choice' trials:
dat.long.choices<-dat.long[dat.long$mate_choice!="NC",]
dat.long.choices<-droplevels(dat.long.choices)
dat.long.choices$mate_choice <- as.factor(dat.long.choices$mate_choice)

#get data to test for side effect (side of arena from which each stimulus was played was alternated
#between trials, but this will allow for an explicit test of side effect)
#how many trials was Fast=R vs. Fast = L?
table(dat.long.choices[,c(7,8)]) #fast on left = 21 + 22 = 43; fast on right = 20 + 21 = 41

#create a column for side of speaker playing fast call
fast.side <- rep(NA, times = nrow(dat.long.choices)) #blank column
fast.side <- factor(fast.side, levels = c("L", "R"))
dat.long.choices <- cbind.data.frame(dat.long.choices, fast.side)

for (i in 1:nrow(dat.long.choices)) {
  if(is.na(dat.long.choices[i,7])==T) {
    next
  }
  #where stimuli are FvS, leading speaker = fast speaker
  if(dat.long.choices[i,7]=="FvS") {
    dat.long.choices$fast.side[i] <-  dat.long.choices[i,8]
  }
  #where stimuli are SvF, leading speaker = slow speaker. 
    # where stimuli are SvF and leading speaker = L, fast speaker = R
  if(dat.long.choices[i, 7] == "SvF" & dat.long.choices[i,8] == "L") {
    dat.long.choices$fast.side[i] <- "R"
  }
  # where stimuli are SvF and leading speaker = R, fast speaker = L
  if(dat.long.choices[i, 7] == "SvF" & dat.long.choices[i,8] == "R") {
    dat.long.choices$fast.side[i] <- "L"
  }
}

library(lme4)
mod1<-glmer(mate_choice ~ chorus + (1|ID), data = dat.long.choices, family = "binomial") 
mod1a<-glmer(mate_choice ~ chorus + fast.side + (1|ID), data = dat.long.choices, family = "binomial") 
mod.null<-glmer(mate_choice ~ (1|ID), data = dat.long.choices, family = "binomial") 
anova(mod1, mod1a, mod.null) #p = 0.04037, chorus treatment sig. predicts choice
summary(mod1)
summary(mod1a)
#including speaker assignment (l vs r to fast call) does not improve model fit.

mod2<-glmer(mate_choice ~ svl + (1|ID), data = dat.long.choices, family = "binomial") 
mod2a<-glmer(mate_choice ~ svl + fast.side + (1|ID), data = dat.long.choices, family = "binomial") 
anova(mod2, mod2a, mod.null) # female size does not predict choice
#including speaker assignment does not sig. improve model fit
summary(mod2)
summary(mod2a)

mod3<-glmer(mate_choice ~ smi + (1|ID), data = dat.long.choices, family = "binomial") 
mod3a<-glmer(mate_choice ~ smi + fast.side + (1|ID), data = dat.long.choices, family = "binomial") 
anova(mod3, mod3a, mod.null) # female condition does not predict choice
anova(mod3a, mod.null)

#test for effect of treatment order on choice 
#(are females more likely to choose slow in pure-SM chorus if they've heard SB chorus previously)
#this is 'first exp' variable in dataframe.
mod.first.exp<-glmer(mate_choice ~ chorus + first_exp + (1|ID), data = dat.long.choices, family = "binomial")
mod.first.expa<-glmer(mate_choice ~ chorus + first_exp + fast.side + (1|ID), data = dat.long.choices, family = "binomial")
anova(mod1, mod.first.exp, mod.first.expa) # no sig. effect of treatment order or trt order + side

#reviewer question: do we have sufficient power to detect a treatment order effect
#conducting a power analysis for a glmer using simulations via simr
library(simr)
# NAs in the dataframe were causing errors with the simulation, re-fitting with na-omitted dataframe:
mod.first.exp2<-glmer(mate_choice ~ chorus + first_exp + (1|ID), data = na.omit(dat.long.choices), family = "binomial")
#the simulations below take a while to run
powerSim(mod.first.exp2, test = fixed("first_exp", method = "lr"), nsim = 1000, progress = F)
powerSim(mod.first.exp2, test = fixed("chorus", method = "lr"), nsim = 1000, progress = F)


#alternative approach for reviewers: backwards model selection
full.mod <- glmer(mate_choice ~ chorus + smi + svl + first_exp + fast.side + (1|ID), data = dat.long.choices, family = "binomial")
drop1(full.mod, test = "Chisq")
drop1(update(full.mod,  ~ . -svl), test = "Chisq")
drop1(update(full.mod,  ~ . -svl -smi), test = "Chisq")
drop1(update(full.mod,  ~ . -svl -smi -first_exp), test = "Chisq")
drop1(update(full.mod,  ~ . -svl -smi -first_exp -fast.side), test = "Chisq")
#best model includes only chorus background as significant predictor of choice.  


#post hoc tests using emmeans
library(emmeans)
emmeans(mod1, specs = "chorus", infer =c(T, T), type = "response")


#EXAMINING LATENCIES  

#identify distribution that best fits latency data
hist(dat.long.choices$latency) #looks like poisson, let's check fit quantitatively
library(fitdistrplus)
no.na.latency <- dat.long.choices$latency[is.na(dat.long.choices$latency)==F]
poissonfit<-fitdistrplus::fitdist(data = no.na.latency, distr = "pois")
nbfit<-fitdistrplus::fitdist(data = no.na.latency, distr = "nbinom")
qqcomp(list(poissonfit, nbfit), addlegend = T,
       xlegend = "topleft")
#negative binomial is best fit. 

#does latency differ by chorus treatment?  or female phenotype?
mod.latency<-glmer.nb(formula = latency ~ chorus + fast.side + (1|ID), data = dat.long.choices)
drop1(mod.latency, test = "Chisq") #latency does not differ between chorus treatments

mod2.latency<-glmer.nb(formula = latency ~ scale(svl) + fast.side + (1|ID), data = dat.long.choices,
                       control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=1e6))) #had some convergence issues with this one, so changed the optimizer settings
drop1(mod2.latency, test = "Chisq") #latency marginally (not sig.) higher for bigger females

mod3.latency<-glmer.nb(formula = latency ~ smi + fast.side + (1|ID), data = dat.long.choices)
drop1(mod3.latency, test = "Chisq")

#backwards model selection approach for reviewers:
mod.full.latency <- glmer.nb(formula = latency ~ chorus + scale(svl) + smi + fast.side + (1|ID), data = dat.long.choices)
drop1(mod.full.latency, test = "Chisq")
drop1(update(mod.full.latency,  ~ . -fast.side), test = "Chisq")
drop1(update(mod.full.latency,  ~ . -fast.side -chorus), test = "Chisq")
drop1(update(mod.full.latency,  ~ . -fast.side -chorus -smi), test = "Chisq")


#is latency associated with mate choice?
mod.choice.latency <- glmer(formula = mate_choice ~ latency + chorus + (1|ID), data = dat.long.choices, family = "binomial")
mod.null <- glmer(formula = mate_choice ~ (1|ID), data = dat.long.choices, family = "binomial")
anova(mod.choice.latency, mod.null, test = "Chisq")

plot(mate_choice ~ latency, data = dat.long.choices) 
#short and medium latencies are fairly evenly distributed among F vs S choices,
#but for latencies over 10-15 min, slow choices are more common.  

#make a good plot of how latencies differ by female choice (Fig. S5)
library(ggplot2)
plot.hist.latency <- ggplot(data = dat.long.choices, 
                         mapping = aes(x = latency, color = mate_choice,
                                       fill = mate_choice)) +
  geom_histogram(position = "dodge", bins = 10)+
  theme_classic()+
  theme(legend.position = "top")+
  scale_color_grey()+
  scale_fill_grey()+
  labs(fill = "Mate Choice")+
  labs(color = "Mate Choice")+
  xlab("Latency")+
  ylab("Count")
plot.hist.latency

#export the plot
svg(filename="supp.latency.hist.svg", 
    width=4, 
    height=3, 
    pointsize=12, 
    antialias = c("subpixel")
)
plot.hist.latency
dev.off()

# Additional QC details
#females were tested in 2 bouts due to an interruption in lab work
checkdates.short <- dat.wide[dat.wide$days.btw < 20,] #first bout
length(unique(checkdates.short$ID)) #34 females
unique(checkdates.short$days.btw) #13-14 days between testing

checkdates.long <- dat.wide[dat.wide$days.btw > 20,] #second bout
length(unique(checkdates.long$ID)) #19 females
unique(checkdates.long$days.btw) #49-51 days between testing

#test for time of testing effect
#first need to get day-of-testing on a numeric scale
dat.long.choices$date.formatted<-as.Date(as.character(dat.long.choices$date), format = "%Y%m%d")
unique(dat.long.choices$date.formatted) #7 unique testing dates
dat.long.choices$date.numeric <- dat.long.choices$date.formatted - as.Date(as.character(20200926), format = "%Y%m%d")
#now test for effect of date of testing on preferences:
mod.date <- glmer(mate_choice ~ date.numeric + (1|ID), data = dat.long.choices, family = "binomial")
mod.date.null <- glmer(mate_choice ~ (1|ID), data = dat.long.choices, family = "binomial")
anova(mod.date, mod.date.null)


#MAKING FIGURE 1
#Panel A
bar.counts <- c(22,18,15,29)
chorus <- c("Pure Species", "Pure Species", "Mixed Species", "Mixed Species")
chorus <- factor(chorus, levels = c("Pure Species", "Mixed Species"))
Stimulus <- c("Fast", "Slow","Fast", "Slow")
bar.df <- cbind.data.frame(bar.counts, chorus, Stimulus)

panelA <- ggplot(data = bar.df, aes(x=chorus, y=bar.counts, fill = Stimulus)) +
  geom_bar(stat="identity", color = "black", position = position_dodge())+
  theme_classic()+
  ylab("Number of females choosing stimulus")+
  ylim(0, 34)+
  scale_fill_manual(values = c("darkgrey", "white"))+
  xlab("Chorus background")+
  theme(legend.position="top")+
  theme(legend.title= element_blank())+
#annotate with stats
  geom_segment(x = 0.7, xend = 1.3, y = 25, yend = 25)+
  geom_segment(x = 1.7, xend = 2.3, y = 32, yend = 32)+
  annotate("text", label ="p = 0.64", x = 1, y = 27, size = 3.3, fontface = "italic")+
  annotate("text", label ="p = 0.049", x = 2, y = 34, size = 3.3, fontface = "italic")+
  theme(axis.title = element_text(size = 14))+
  theme(axis.text=element_text(size=12))



# panel B 
#remove NC data for plot (only want to plot choices)
dat.wide.plotting <- dat.wide
dat.wide.plotting[is.na(dat.wide.plotting$mate_choice_1)==F & dat.wide.plotting$mate_choice_1=="NC",13] <- NA
dat.wide.plotting[is.na(dat.wide.plotting$mate_choice_2)==F & dat.wide.plotting$mate_choice_2=="NC",24] <- NA
dat.wide.plotting <- droplevels(dat.wide.plotting)

#get mate choice in pure-SM and mixed chorus treatment in wide format data
mate.choice.pure <- rep(NA, times = nrow(dat.wide.plotting))
mate.choice.mixed <- rep(NA, times = nrow(dat.wide.plotting))
mate.choice.pure <- factor(mate.choice.pure, levels = c("S", "F"))
mate.choice.mixed <- factor(mate.choice.mixed, levels = c("S", "F"))

for(i in 1:nrow(dat.wide.plotting)) {
  if(dat.wide.plotting$chorus_1[i] == "pureSM") {
    mate.choice.pure[i] <- dat.wide.plotting$mate_choice_1[i]
    mate.choice.mixed[i] <- dat.wide.plotting$mate_choice_2[i]
  }
  if(dat.wide$chorus_1[i] == "mixed") {
    mate.choice.mixed[i] <- dat.wide.plotting$mate_choice_1[i]
    mate.choice.pure[i] <- dat.wide.plotting$mate_choice_2[i]
  }
}

dat.wide.plotting <- cbind.data.frame(dat.wide.plotting, mate.choice.pure, mate.choice.mixed)

#need to add some vertical jitter in the y-direction in order to see the different individual
#females' points/lines.  
#to make this visually clear, we are going to categorize the females by their choice patterns,
#then add the jitter sequentially alternating between different choice patterns.  
#start by getting females IDs grouped by their choice patterns:
#Fast in both:
f.both.ids <- droplevels(dat.wide.plotting[is.na(dat.wide.plotting$mate.choice.mixed)==F &
                                  dat.wide.plotting$mate.choice.mixed== "F" &
                                  is.na(dat.wide.plotting$mate.choice.pure)==F &
                                  dat.wide.plotting$mate.choice.pure=="F",3]) #8 females
#slow in both 
s.both.ids <- droplevels(dat.wide.plotting[is.na(dat.wide.plotting$mate.choice.mixed)==F &
                                             dat.wide.plotting$mate.choice.mixed== "S" &
                                             is.na(dat.wide.plotting$mate.choice.pure)==F &
                                             dat.wide.plotting$mate.choice.pure=="S",3]) #11 females
#'adaptive' choice pattern: F in pure, S in mixed:
adaptive.ids <- droplevels(dat.wide.plotting[is.na(dat.wide.plotting$mate.choice.mixed)==F &
                                             dat.wide.plotting$mate.choice.mixed== "S" &
                                             is.na(dat.wide.plotting$mate.choice.pure)==F &
                                             dat.wide.plotting$mate.choice.pure=="F",3]) #11 females

#'maladaptive' choice pattern: S in pure, F in mixed:
mal.ids <- droplevels(dat.wide.plotting[is.na(dat.wide.plotting$mate.choice.mixed)==F &
                                               dat.wide.plotting$mate.choice.mixed== "F" &
                                               is.na(dat.wide.plotting$mate.choice.pure)==F &
                                               dat.wide.plotting$mate.choice.pure=="S",3]) #3 females
#partial data: females that made a choice in only 1 chorus condition
par.ids <- droplevels(dat.wide.plotting[is.na(dat.wide.plotting$mate.choice.mixed)==T |
                                         is.na(dat.wide.plotting$mate.choice.pure)==T,3]) #20 females
#make a vector of female IDs that alternates between the categories of choice patterns. We will then 
#add a constant 0.01 unit jitter in the y direction in that sequence. this is purely for plotting aesthetics,
#so that the jitter maintains parallel lines between females (readability) and 
#maximizes space between points/lines showing females of the same choice category.  
#since there are different numbers of females in each choice category, we're going to (roughly) evenly 
#intersperse them.  
#I worked out an order (by hand on paper) for female ID (given in the vector below) that alternates 
#between  the choice categories. 
#again the jitter is arbitrary and this ordering is just to make things look tidy.
ids.order <- c("pfennfld10294", "7819_LF4", "pfennfld10111", "pfennfld11364", "pfennfld10289", 
               "pfennfld10293", "pfennfld10213", "pfennfld10116", "pfennfld10279", "pfennfld11363", 
               "11317 PFENNFLD", "1732_LF2", "10333 PFENNFLD", "pfennfld10291", "12162 PFENNFLD", 
               "11416 PFENNFLD", "1126_LF2", "pfennfld11422", "pfennfld10297", "1732_LF1", 
               "10316 PFENNFLD", "4116_LF2", "1124_LF2", "11316 PFENNFLD", "pfennfld10215", "6800_LF1",
               "12163 PFENNFLD", "pfennfld11401", "12021 PFENNFLD", "9801_LF4", "8901_LF2", "1505_LF3",
               "3500_RF3", "9907_LF3", "1832_LF4", "10290 PFENNFLD", "12165 PFENNFLD", "12164 PFENNFLD",
               "12102 PFENNFLD", "10207 PFENNFLD", "10319 PFENNFLD", "1504_LF4", "1504_LF1", 
               "10328 PFENNFLD", "12168 PFENNFLD", "1847_LF3", "10085 PFENNFLD", "10325 PFENNFLD", 
               "9901_LF3", "6800_LF2", "12020 PFENNFLD", "10205 PFENNFLD", "9912_LF3")
#empty vector for order # in the plotting dataframe:
dat.wide.plotting$jitter.order <- rep(NA, times = nrow(dat.wide.plotting))

#fill the jitter.order vector according to ids.order
for(i in 1:length(ids.order)){
  dat.wide.plotting[dat.wide.plotting$ID == ids.order[i], 32]<- i
}

#sort the dataframe by jitter.order
dat.wide.plotting <- dat.wide.plotting[order(dat.wide.plotting$jitter.order),]

#now to actually make and add the jitter:
nrow(dat.wide.plotting) #53
#make some offsets so that females are staggered in the y direction so each female is visible
yax.displace <- seq(from = -0.26, to = 0.26, by = 0.01) #53 displacements, one for each female

jittered.pure <- yax.displace + as.numeric(dat.wide.plotting$mate.choice.pure)
jittered.mixed <- yax.displace + as.numeric(dat.wide.plotting$mate.choice.mixed)
plot.data <- cbind.data.frame(dat.wide.plotting,jittered.pure, jittered.mixed)

panelB <- ggplot(data = plot.data, mapping = aes(x = 1, y = 1))+
  #add points showing females' choices
  geom_point(x = 0, y = jittered.pure, shape = 21, size = 3)+
  geom_point(x = 1, y = jittered.mixed, shape = 21, size = 3)+
  theme_classic()+
  scale_x_continuous(limits = c(-0.05, 1.2), breaks = c(0, 1), 
                     labels = c("Pure Species", "Mixed Species"), 
                     name = "Chorus Background")+
  scale_y_continuous(limits = c(0.7,2.3), breaks = c(1,2), 
                     labels = c("Slow", "Fast"), name = "Female choice")

#add lines connecting choices for each female's choices
panelB <- panelB + 
  geom_segment(data = plot.data, mapping = aes(x = 1, y = 1), 
               x = rep(0, times = 53), y = plot.data$jittered.pure, 
               xend = rep(1, times = 53), yend =plot.data$jittered.mixed, 
               size = 0.4)

panelB <- panelB + theme(axis.text=element_text(size=12)) +
  theme(axis.title = element_text(size = 14))

#make it a 2 panel plot
library(ggpubr)
fig1 <- ggarrange(panelA, 
                  ggplot() + theme_void(), #add some more space between plots A and B
                  panelB, nrow = 3, ncol = 1, heights = c(0.9, 0.1, 0.9),
                  labels = c("A", "", "B"))

#export the plot
svg(filename="fig1.svg", 
    width=4, 
    height=10, 
    pointsize=12, 
    antialias = c("subpixel")
)
fig1
dev.off()
  

