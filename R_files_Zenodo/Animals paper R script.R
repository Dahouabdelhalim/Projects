
# 2020 SFS experiment -  analysis of 12 wild-caught birds for Animals paper

# This script analyses the training data and the first 11 days of continuous foraging
# 
# ------------------------------------------------------------------------------------

# Melissa Bateson
# 15th March 2022


rm(list=ls())

# load necessary packages
library(tidyverse)
library(lme4)
library(lmerTest)
library(lubridate)
library(base)
library(irr)
library(gridExtra)




#### Autoshaping data ####


load("Training data.RData")


# Plot Figure 4
Figure4 = ggplot(data = t1, aes(x = session, y = good.trials, colour = bird)) +
  geom_point() +
  geom_line() +
  xlab("Session") +
  scale_x_continuous(breaks= c(1,2,3), labels = c(1,2,3)) +
  ylab("No. trials with a peck") +
  scale_colour_discrete(name = "Bird") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position=c(0.1,0.6)) +
  theme(legend.box.background = element_rect(color="black", size=0.5))
Figure4





#### Daily food consumption data #### 


load("Consumption data.RData")



# Basic descriptive stats
mean.eaten.room = f1 %>%
  group_by(room) %>%
  summarise(mean.eaten.room = mean(total.eaten, na.rm = T))

mean.eaten = mean.eaten.room %>%
  summarise(mean.eaten = mean(mean.eaten.room, na.rm = T),
            sd.eaten = sd(mean.eaten.room)) %>%
  mutate(se.eaten = sd.eaten/sqrt(2))

# Did the amount eaten per day change linearly over the 11 days?
m = lmer(total.eaten ~ day + (1|room), data = f1)
plot(m)
hist(resid(m))
summary(m)
anova(m)

# create a date.room variable for later data merging
f1$date.room = paste(f1$date, f1$room)
f1$date = NULL
f1$room = NULL





#### Key pecking data ####



load("Peck data.RData")




# Is the total amount eaten in a room on a day predicted by the number of reinforcements earned?

# Compute the number of pecks per aviary day
pecks.date = p %>%
  group_by(room, date.room) %>%
  summarise(total.pecks = n())

# Merge the consumption data with the pecking data
f1 = merge(f1, pecks.date, by = "date.room")


# Does the total pecks in a room per day predict the total eaten in that room that day?

m = lmer(total.eaten ~ total.pecks + (1|room), data = f1, na.action = na.exclude)
plot(m)
hist(resid(m))
summary(m) # the total amount eaten in a day is predicted by the pecks
anova(m)



pecks.hr = p %>%
  group_by(room,sex,bird,day,bin.1h) %>%
  summarise(total.pecks = n()) %>%
  mutate(bird.day.hour = paste(bird,day,bin.1h))

# note that the birds only had 15 mins to forage after 1700, therefore need to convert values to peck rates/hr
# values for bin.1hr = 17 need to be multiplied by 4
pecks.hr$peck.rate[pecks.hr$bin.1h<17] = pecks.hr$total.pecks[pecks.hr$bin.1h<17]
pecks.hr$peck.rate[pecks.hr$bin.1h==17] = pecks.hr$total.pecks[pecks.hr$bin.1h==17]*4


pecks.day = pecks.hr %>%
  group_by(room,sex,bird,day) %>%
  summarise(total.pecks.day = sum(total.pecks)) %>%
  ungroup() %>%
  mutate(bird.day = paste(bird,day))

pecks.bird = pecks.day %>%
  group_by(room,sex,bird) %>%
  summarise(total.pecks = sum(total.pecks.day),
            mean.pecks.day = mean(total.pecks.day, na.rm = TRUE),
            sd.pecks.day = sd(total.pecks.day, na.rm = TRUE)) %>%
  ungroup()



sum.pecks = sum(pecks.day$total.pecks.day) # 34,426
mean.pecks.per.day = mean(pecks.day$total.pecks.day)
sd.pecks.per.day = sd(pecks.day$total.pecks.day)
mean.pecks.per.hour = mean.pecks.per.day/9.25 # mean number of pecks/bird/hr = 28.19


# Did the rate of reinforcement change linearly over the 11 days of continuous foraging?

m = lmer(total.pecks.day ~ day + sex + (1|bird), data = pecks.day) # room explains none of variance so leave out
plot(m)
hist(resid(m))
summary(m)
anova(m)


# Did rate of pecking change as a function of time of day?

m = lmer((peck.rate) ~ bin.1h + sex + (1|bird), data = pecks.hr) # room explains none of variance so leave out
plot(m) # argument for transformation being necessary?
hist(resid(m))
summary(m) # peck rate declines as hour increases
anova(m)

# Plot Figure 5 - pecks per hour by bird
mean.pecks.hr = pecks.hr %>%
  group_by(room,sex,bird,bin.1h) %>%
  summarise(mean.peck.rate = mean(peck.rate, na.rm = TRUE),
            sd.peck.rate = sd(peck.rate, na.rm = TRUE),
            n = n()) %>%
  mutate(CI = 1.96*(sd.peck.rate/sqrt(n)),
         label = paste(room,bird))


Figure5 = ggplot(data = mean.pecks.hr, aes(x=bin.1h, y=mean.peck.rate, colour = sex)) +
  scale_colour_manual(values=c("magenta", "blue")) +
  scale_x_continuous(breaks = c(8,10,12,14,16)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin=(mean.peck.rate-CI), ymax=(mean.peck.rate+CI)), width=0.2) +
  facet_wrap(~label, ncol = 6) +
  ylab(expression(Reinforcements.hour^-1)) +
  xlab("Time") +
  guides(colour = "none") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
Figure5



# Is total pecks per day repeatable within bird?

pecks.per.day = select(pecks.day, "bird", "day", "total.pecks.day")
p2 = pivot_wider(data = pecks.per.day,names_from = "day",values_from = "total.pecks.day")
p2 = select(p2,"1","2","3","4","5","6","7","8","9","10","11")
icc(ratings = p2, model = "twoway", type = "agreement",r0 = 0, conf.level = 0.95)
# This gives an ICC of 0.85 with 95% CI of 0.725 - 0.944 (i.e. moderate-excellent reliability)







#### Mass data ####



load("Mass data.RData")


d=merge(d, f1, by="date.room", all.x=FALSE)


# Remove any masses that are biologically completely impossible and likely due to calibration problems
# Note that massess below 50 are not recorded, but very high massess occur occasionally
d = filter(d,mass.g<120)

# Total masses before filtering = 8857
# After filtering = 8841
# No > 120 g = 16
# = 0.18 % of masses excluded 


# sort the dataset by date and time
d=d[order(d$date.time), ]

# Number of masses
no.masses = d %>%
  group_by(room,sex,bird,day) %>%
  summarise(total.masses = n())

mass.stats = no.masses %>%
  group_by(room,sex,bird) %>%
  summarise(total = sum(total.masses, na.rm = TRUE), mean.per.day = mean(total.masses, na.rm = TRUE),
            sd = sd(total.masses, na.rm = TRUE))

mean.masses = sum(mass.stats$total)/132 # mean number of masses/bird/day = 66.98
mean.masses.per.hour = mean.masses/9.25 # mean number of masses/bird/hr = 7.24


d$row.number = seq(1:length(d$mass.g)) # Add a unique row number for later merging


#### Polynomial fitting for 12 birds in rooms 114 and 116 ####

# The aim of this part of the script is to fit the mass data obtained for each individual bird on each day of the experiment.
# The reason for this is to obtain a model of the mass change over the day that eliminates the noise in the measurements.
# We can then use the fitted values from the model in subsequent analyses (e.g. the maximum fitted mass in a day; the fitted mass at 1800 hrs)

# set the mass band either side of the fitted function for retaining masses
# I have finally settled on quite a high value that only excludes very extreme masses.
g.filter = 10

# make a new data frame to accumulate cleaned and fitted masses in
d.clean = data.frame(time.dec = numeric(0),time.s.clean = numeric(0), mass.g.clean= numeric(0), m.fit = numeric(0))
#make empty vectors for results
bird = NULL
day.treat = NULL
no.masses = NULL
no.masses.clean = NULL
AIC.linear = NULL
AIC.quadratic = NULL
AIC.cubic = NULL
linear.p = NULL
quadratic.p = NULL
cubic.p = NULL
r_sq.fit = NULL
AIC.fit = NULL
mass.earliest = NULL
mass.latest = NULL
mass.max = NULL
mass.min = NULL
mass.7 = NULL
mass.8 = NULL
mass.9 = NULL
mass.10 = NULL
mass.11 = NULL
mass.12 = NULL
mass.13 = NULL
mass.14 = NULL
mass.15 = NULL
mass.16 = NULL
mass.17 = NULL
time.max.mass = NULL
time.min.mass = NULL
time.earliest.mass = NULL
time.latest.mass = NULL

row = 1

for(i in c("P74",  "P79",  "P98",  "P75",  "P91",  "P99", "P77",  "P78",  "P80",  "P92",  "P93",  "P95")){ print(i)
    
    #if (i=="P74"|i=="P79"|i=="P98"|i=="P75"|i=="P91"|i=="P99")   {no.days= 11} 
    #if (i=="P77"|i=="P78"|i=="P80"|i=="P92"|i=="P93"|i=="P95")   {no.days= 12} 
  
  no.days = 11 # use 11 days since we have these for all birds
    
    for(k in 1:no.days){ print(k)
      
      bird[row] = i
      day.treat[row] = k
      
      s = subset(d, (bird == i & day == k))
      
      # stuff that goes inside the function
      
      no.masses[row] = length(s$mass.g) # save the number of masses that day
      
      # check to see that there is enough data to do the fitting
      # 'degree' must be less than number of unique points; need degree+1 points for fit that is not just the data
      
      if (length(s$mass.g)<10) {
        # Set fit parameters to NA
        AIC.linear[row] = NA
        AIC.quadratic[row] = NA
        AIC.cubic[row] = NA
        linear.p[row] = NA
        quadratic.p[row] = NA
        cubic.p[row] = NA
        r_sq.fit[row] = NA
        AIC.fit[row] = NA
        mass.earliest[row] = NA
        mass.latest[row] = NA
        mass.max[row] = NA
        mass.min[row] = NA
        mass.7[row] = NA
        mass.8[row] = NA
        mass.9[row] = NA
        mass.10[row] = NA
        mass.11[row] = NA
        mass.12[row] = NA
        mass.13[row] = NA
        mass.14[row] = NA
        mass.15[row] = NA
        mass.16[row] = NA
        mass.17[row] = NA
        time.max.mass[row] = NA
        time.min.mass[row] = NA
        time.earliest.mass[row] = NA
        time.latest.mass[row] = NA
        
        row = row + 1
      }
      
      if (length(s$mass.g) >=10) {
        
        # compare different degrees of polynomial
        ml = lm(mass.g ~ poly(time.dec,degree=1), data = s, na.action = na.exclude)
        mq = lm(mass.g ~ poly(time.dec,degree=2), data = s, na.action = na.exclude)
        mc = lm(mass.g ~ poly(time.dec,degree=3), data = s, na.action = na.exclude)
        
        AIC.linear[row] = AIC(ml)
        AIC.quadratic[row] = AIC(mq)
        AIC.cubic[row] = AIC(mc)
        
        # now fit the cubic
        
        m = lm(mass.g ~ poly(time.dec,degree=3), data = s) # NB this fits a cubic orthogonal polynomial
        
        s$m.fit.raw = predict(m)
        
        # calculate a fixed mass band for inclusion either side of the fitted slope
        s$m.lb=(s$m.fit.raw - g.filter)
        s$m.ub=(s$m.fit.raw + g.filter)
        # filter out any data outside these bounds
        s1 = subset(s, (mass.g<(m.fit.raw+g.filter) & mass.g>(m.fit.raw-g.filter)))
        s1$mass.g.clean = s1$mass.g
        s1$time.s.clean = s1$time.dec
        
        no.masses.clean[row] = length(s1$mass.g.clean) # save the number of masses retained in the final cleaned dataset
        
        # recheck to see that there is enough data to do the fitting after the cleaning
        
        if (length(s1$mass.g.clean)<10) {
          # Set fit parameters to NA
          AIC.linear[row] = NA
          AIC.quadratic[row] = NA
          AIC.cubic[row] = NA
          linear.p[row] = NA
          quadratic.p[row] = NA
          cubic.p[row] = NA
          r_sq.fit[row] = NA
          AIC.fit[row] = NA
          mass.earliest[row] = NA
          mass.latest[row] = NA
          mass.max[row] = NA
          mass.min[row] = NA
          mass.7[row] = NA
          mass.8[row] = NA
          mass.9[row] = NA
          mass.10[row] = NA
          mass.11[row] = NA
          mass.12[row] = NA
          mass.13[row] = NA
          mass.14[row] = NA
          mass.15[row] = NA
          mass.16[row] = NA
          mass.17[row] = NA
          time.max.mass[row] = NA
          time.min.mass[row] = NA
          time.earliest.mass[row] = NA
          time.latest.mass[row] = NA
          
          row = row + 1
        }
        
        if (length(s1$mass.g.clean) >=10) {
          
          # refit the cleaned dataset
          m1 = lm(mass.g.clean ~ poly(time.s.clean,degree=3), data = s1, na.action = na.exclude) 
          s1$m.fit = predict(m1) # save the fitted values from the curve fitting
          
          # extract coefficients and p-values from the fit
          linear.p[row] = summary(m1)$coefficients[2,4] # p-value for the linear term in the orthogonal fit
          quadratic.p[row] = summary(m1)$coefficients[3,4] # p-value for the quadratic term in the orthogonal fit
          cubic.p[row] = summary(m1)$coefficients[4,4] # p-value for the cubic term in the orthogonal fit
          
          # NB need to fit non-orthogonal polynomial to get the correct coefficients for the fitted line
          m2 = lm(mass.g.clean ~ time.s.clean + I(time.s.clean^2) + I(time.s.clean^3), data = s1, na.action = na.exclude)
          # eqn for cubic: y = a + bx + cx^2 + dx^3 where x is the minutes since midnight
          # extract parameter estimates b, c and d
          par.a = as.numeric(m2$coefficients[1])
          par.b = as.numeric(m2$coefficients[2])
          par.c = as.numeric(m2$coefficients[3])
          par.d = as.numeric(m2$coefficients[4])
          
          # calculate the fitted mass at each hour by plugging the relevant times into the fitted funtion
          mass.7[row] = par.a + (par.b*420) + (par.c*(420^2)) + (par.d*(420^3))
          mass.8[row] = par.a + (par.b*480) + (par.c*(480^2)) + (par.d*(480^3))
          mass.9[row] = par.a + (par.b*540) + (par.c*(540^2)) + (par.d*(540^3))
          mass.10[row] = par.a + (par.b*600) + (par.c*(600^2)) + (par.d*(600^3))
          mass.11[row] = par.a + (par.b*660) + (par.c*(660^2)) + (par.d*(660^3))
          mass.12[row] = par.a + (par.b*720) + (par.c*(720^2)) + (par.d*(720^3))
          mass.13[row] = par.a + (par.b*780) + (par.c*(780^2)) + (par.d*(780^3))
          mass.14[row] = par.a + (par.b*840) + (par.c*(840^2)) + (par.d*(840^3))
          mass.15[row] = par.a + (par.b*900) + (par.c*(900^2)) + (par.d*(900^3))
          mass.16[row] = par.a + (par.b*960) + (par.c*(960^2)) + (par.d*(960^3))
          mass.17[row] = par.a + (par.b*1020) + (par.c*(1020^2)) + (par.d*(1020^3))
          
          
          # Extract parameters from the fit
          r_sq.fit[row] = summary(m1)$r.squared # r-squared for the fit
          AIC.fit[row] = AIC(m1) # AIC for the final fit
          mass.earliest[row] = s1$m.fit[which(s1$time.s.clean==min(s1$time.s.clean))] # the earliest fitted mass
          mass.latest[row] = s1$m.fit[which(s1$time.s.clean==max(s1$time.s.clean))] # the latest fitted mass
          mass.max[row] = max(s1$m.fit) # the maximumn fitted mass
          mass.min[row] = min(s1$m.fit) # the minimum fitted mass
          time.max.mass[row] = median(s1$time.s.clean[which(s1$m.fit==max(s1$m.fit))]) # finds the time with the max value (median if more than one)
          time.min.mass[row] = median(s1$time.s.clean[which(s1$m.fit==min(s1$m.fit))]) # finds the time with the min value (median if more than one)
          time.earliest.mass[row] = min(s1$time.s.clean) # the time at which the first fitted mass occurred
          time.latest.mass[row] = max(s1$time.s.clean) # the time at which the last fitted mass occurred
          
          # Don't extrapolate beyond the data available each day 
          # exclude any fits for which the latest mass that day is more than one hour previously.
          if (time.latest.mass[row] < 720) {mass.13[row]=NA}
          if (time.latest.mass[row] < 780) {mass.14[row]=NA}
          if (time.latest.mass[row] < 840) {mass.15[row]=NA}
          if (time.latest.mass[row] < 900) {mass.16[row]=NA}
          if (time.latest.mass[row] < 960) {mass.17[row]=NA}
          
          # exclude any fits for which the earliest mass that day is more than one hour later.
          if (time.earliest.mass[row] > 540) {mass.8[row]=NA}
          if (time.earliest.mass[row] > 600) {mass.9[row]=NA}
          if (time.earliest.mass[row] > 660) {mass.10[row]=NA}
          if (time.earliest.mass[row] > 720) {mass.11[row]=NA}
        
          row = row + 1
          
          
          s1 = subset(s1, select = c("row.number","time.s.clean","mass.g.clean","m.fit"))
          d.clean = rbind.data.frame(d.clean, s1)
        }
        
      }
  }
}

d.new = merge(d,d.clean, by = "row.number", all.x = TRUE)

fit.results = data.frame(bird,day.treat, no.masses,no.masses.clean, AIC.linear, AIC.quadratic, AIC.cubic, linear.p, quadratic.p, cubic.p, r_sq.fit, AIC.fit, mass.earliest, 
                         mass.latest, mass.max, mass.min, mass.7, mass.8,
                         mass.9, mass.10, mass.11, mass.12, mass.13, mass.14, mass.15, mass.16, mass.17,
                         time.max.mass, time.min.mass, time.earliest.mass, time.latest.mass)

# create a room column

fit.results = fit.results %>% 
  mutate(sex = recode(bird,"P77"="Female", "P78"="Male", "P80"="Male", "P92"="Female", "P93"="Male", "P95"="Female",
                       "P74"="Male", "P75"="Female", "P79"="Male", "P91"="Female", "P98"="Male", "P99"="Female")) %>%
  mutate(room = recode(bird,"P77"=114, "P78"=114, "P80"=114, "P92"=114, "P93"=114, "P95"=114,
                     "P74"=116, "P75"=116, "P79"=116, "P91"=116, "P98"=116, "P99"=116))

fit.results$room = as.factor(fit.results$room)

fit.results$linear.p = round(fit.results$linear.p,4)
fit.results$quadratic.p = round(fit.results$quadratic.p,4)
fit.results$cubic.p = round(fit.results$cubic.p,4)






# Plot the raw mass data for each bird with the fits - Figure 3 and Figures S1-S11

FigureS1 = ggplot(data=filter(d.new,bird=="P74"), aes(x=time.only, y=mass.g)) +
  geom_point(colour = "red") +
  geom_point(aes(x=time.only, y=mass.g.clean), colour = "dark gray") +
  facet_wrap(~day) +
  geom_smooth(aes(x=time.only, y=m.fit), method = "loess", colour = "black", size = 0.5)  +
  ylab("Mass (g)") +
  xlab("Time") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
FigureS1

FigureS2 = ggplot(data=filter(d.new,bird=="P75"), aes(x=time.only, y=mass.g)) +
  geom_point(colour = "red") +
  geom_point(aes(x=time.only, y=mass.g.clean), colour = "dark gray") +
  facet_wrap(~day) +
  geom_smooth(aes(x=time.only, y=m.fit), method = "loess", colour = "black", size = 0.5)  +
  ylab("Mass (g)") +
  xlab("Time") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
FigureS2

Figure3 = ggplot(data=filter(d.new,bird=="P77"), aes(x=time.only, y=mass.g)) +
  geom_point(colour = "red", size = 1) +
  geom_point(aes(x=time.only, y=mass.g.clean), colour = "dark gray") +
  facet_wrap(~day) +
  geom_smooth(aes(x=time.only, y=m.fit), method = "loess", colour = "black", size = 0.5)  +
  ylab("Mass (g)") +
  xlab("Time") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
Figure3


FigureS3= ggplot(data=filter(d.new,bird=="P78"), aes(x=time.only, y=mass.g)) +
  geom_point(colour = "red") +
  geom_point(aes(x=time.only, y=mass.g.clean), colour = "dark gray") +
  facet_wrap(~day) +
  geom_smooth(aes(x=time.only, y=m.fit), method = "loess", colour = "black", size = 0.5)  +
  ylab("Mass (g)") +
  xlab("Time") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
FigureS3

FigureS4 = ggplot(data=filter(d.new,bird=="P79"), aes(x=time.only, y=mass.g)) +
  geom_point(colour = "red") +
  geom_point(aes(x=time.only, y=mass.g.clean), colour = "dark gray") +
  facet_wrap(~day) +
  geom_smooth(aes(x=time.only, y=m.fit), method = "loess", colour = "black", size = 0.5)  +
  ylab("Mass (g)") +
  xlab("Time") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
FigureS4

FigureS5 = ggplot(data=filter(d.new,bird=="P80"), aes(x=time.only, y=mass.g)) +
  geom_point(colour = "red") +
  geom_point(aes(x=time.only, y=mass.g.clean), colour = "dark gray") +
  facet_wrap(~day) +
  geom_smooth(aes(x=time.only, y=m.fit), method = "loess", colour = "black", size = 0.5)  +
  ylab("Mass (g)") +
  xlab("Time") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
FigureS5

FigureS6 = ggplot(data=filter(d.new,bird=="P91"), aes(x=time.only, y=mass.g)) +
  geom_point(colour = "red") +
  geom_point(aes(x=time.only, y=mass.g.clean), colour = "dark gray") +
  facet_wrap(~day) +
  geom_smooth(aes(x=time.only, y=m.fit), method = "loess", colour = "black", size = 0.5)  +
  ylab("Mass (g)") +
  xlab("Time") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
FigureS6

FigureS7 = ggplot(data=filter(d.new,bird=="P92"), aes(x=time.only, y=mass.g)) +
  geom_point(colour = "red") +
  geom_point(aes(x=time.only, y=mass.g.clean), colour = "dark gray") +
  facet_wrap(~day) +
  geom_smooth(aes(x=time.only, y=m.fit), method = "loess", colour = "black", size = 0.5)  +
  ylab("Mass (g)") +
  xlab("Time") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
FigureS7

FigureS8 = ggplot(data=filter(d.new,bird=="P93"), aes(x=time.only, y=mass.g)) +
  geom_point(colour = "red") +
  geom_point(aes(x=time.only, y=mass.g.clean), colour = "dark gray") +
  facet_wrap(~day) +
  geom_smooth(aes(x=time.only, y=m.fit), method = "loess", colour = "black", size = 0.5)  +
  ylab("Mass (g)") +
  xlab("Time") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
FigureS8

FigureS9 = ggplot(data=filter(d.new,bird=="P95"), aes(x=time.only, y=mass.g)) +
  geom_point(colour = "red") +
  geom_point(aes(x=time.only, y=mass.g.clean), colour = "dark gray") +
  facet_wrap(~day) +
  geom_smooth(aes(x=time.only, y=m.fit), method = "loess", colour = "black", size = 0.5)  +
  ylab("Mass (g)") +
  xlab("Time") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
FigureS9

FigureS10 = ggplot(data=filter(d.new,bird=="P98"), aes(x=time.only, y=mass.g)) +
  geom_point(colour = "red") +
  geom_point(aes(x=time.only, y=mass.g.clean), colour = "dark gray") +
  facet_wrap(~day) +
  geom_smooth(aes(x=time.only, y=m.fit), method = "loess", colour = "black", size = 0.5)  +
  ylab("Mass (g)") +
  xlab("Time") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
FigureS10

FigureS11 = ggplot(data=filter(d.new,bird=="P99"), aes(x=time.only, y=mass.g)) +
  geom_point(colour = "red") +
  geom_point(aes(x=time.only, y=mass.g.clean), colour = "dark gray") +
  facet_wrap(~day) +
  geom_smooth(aes(x=time.only, y=m.fit), method = "loess", colour = "black", size = 0.5)  +
  ylab("Mass (g)") +
  xlab("Time") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
FigureS11



# Plot the fit results for each bird across the hours of the day - Figure 7

fits.time = subset(fit.results, select = c("room","bird","sex","day.treat","mass.8","mass.9",
                                           "mass.10","mass.11","mass.12","mass.13","mass.14",
                                           "mass.15","mass.16"))
fits.time2 = pivot_longer(data = fits.time, cols = c("mass.8","mass.9",
                                                     "mass.10","mass.11","mass.12","mass.13","mass.14",
                                                     "mass.15","mass.16"),values_to = "mass", names_to = "fitted.time")

fits.time2$time[fits.time2$fitted.time=="mass.8"] = 8
fits.time2$time[fits.time2$fitted.time=="mass.9"] = 9
fits.time2$time[fits.time2$fitted.time=="mass.10"] = 10
fits.time2$time[fits.time2$fitted.time=="mass.11"] =11
fits.time2$time[fits.time2$fitted.time=="mass.12"] = 12
fits.time2$time[fits.time2$fitted.time=="mass.13"] = 13
fits.time2$time[fits.time2$fitted.time=="mass.14"] = 14
fits.time2$time[fits.time2$fitted.time=="mass.15"] = 15
fits.time2$time[fits.time2$fitted.time=="mass.16"] = 16


mean.mass.hr = fits.time2 %>%
  group_by(room,sex,bird,time) %>%
  summarise(mean.mass = mean(mass, na.rm = TRUE),
            sd.mass = sd(mass, na.rm = TRUE),
            n = n()) %>%
  mutate(CI = 1.96*(sd.mass/sqrt(n)),
         label = paste(room,bird))


Figure7 = ggplot(data = mean.mass.hr, aes(x=time, y=mean.mass, colour = sex)) +
  scale_colour_manual(values=c("blue", "magenta")) +
  scale_x_continuous(breaks = c(8,10,12,14,16)) +
  scale_y_continuous(breaks = c(70,75,80,85,90,95)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin=(mean.mass-CI), ymax=(mean.mass+CI)), width=0.2) +
  facet_wrap(~label, ncol = 6) +
  ylab("Estimated body mass (g)") +
  xlab("Time") +
  guides(colour = "none") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
Figure7




# How good are the polynomial fits?

# Which degree of polynomial has the lowest AIC?
fit.results$min.AIC = apply(fit.results[,5:7], 1, min)
fit.results$best.fit[fit.results$min.AIC==fit.results$AIC.linear] = "Linear"
fit.results$best.fit[fit.results$min.AIC==fit.results$AIC.quadratic] = "Quadratic"
fit.results$best.fit[fit.results$min.AIC==fit.results$AIC.cubic] = "Cubic"

xtabs(~fit.results$best.fit) # 35 linear, 39 quadratic and 57 cubic: corresponds to 26.7%, 29.8% and 43.4% respectively

fit.results$best.fit = factor(fit.results$best.fit, levels = c("Cubic","Quadratic","Linear")) # change order of levels

Figure6 = ggplot(fit.results, aes(x=bird, fill=best.fit)) +
  geom_bar(position="fill") +
  labs(fill="Best fitting model") +
  xlab("Bird ID") +
  ylab("Proportion of days") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
Figure6

# How good are the cubic fits?
mean(fit.results$r_sq.fit, na.rm = TRUE) # Mean = 0.54 for the fits
sd(fit.results$r_sq.fit, na.rm = TRUE) 




#### Descriptive stats for masses ####

# How many many mass measurements to we have and how many were excluded?

fit.totals = fit.results %>%
  group_by(bird) %>%
  summarise(total = sum(no.masses, na.rm = TRUE), 
            total.clean = sum(no.masses.clean, na.rm = TRUE),
            mean.clean = mean(no.masses.clean, na.rm = TRUE),
            sd.clean = sd(no.masses.clean, na.rm = TRUE),
            mean.r_squared = mean(r_sq.fit, na.rm = TRUE),
            sd.r_squared = sd(r_sq.fit, na.rm = TRUE)) %>%
  mutate(exclusions = total - total.clean)

sum(fit.totals$exclusions)/sum(fit.totals$total)*100  
# proportion excluded = 162/8841 = 1.83%

mean(fit.totals$mean.clean) # mean clean masses per day per bird
sd(fit.totals$mean.clean)

mean(fit.totals$total.clean) # total clean masses per bird
sd(fit.totals$total.clean)


# Daily mass gain and loss

# Compute daily gain
fit.results$daily.gain = fit.results$mass.16-fit.results$mass.8
options(digits = 3)

# Compute overnight mass loss (dusk.mass.yesterday - dawn.mass)
# loss should be related to EE (more specifically sleeping MR)
fit.results = fit.results%>%
  group_by(bird) %>%
  arrange(day.treat) %>%
  mutate(dusk.mass.yesterday = lag(mass.16)) %>%
  ungroup() %>%
  mutate(night.loss = dusk.mass.yesterday-mass.8) %>%
  mutate(bird.day = paste(bird,day.treat))

# Merge in the pecking data
pecks.day = select(pecks.day, c("bird.day","total.pecks.day"))
fit.results = merge(fit.results, pecks.day, by = "bird.day")
  
  


# Descriptive stats by bird
descriptives.bird = fit.results %>%
  group_by(room, sex, bird) %>%
  summarise(dawn.mass = mean(mass.8, na.rm = TRUE), sd.dawn = sd(mass.8, na.rm = TRUE),
            noon.mass = mean(mass.12, na.rm = TRUE), sd = sd(mass.12, na.rm = TRUE),
            dusk.mass = mean(mass.16, na.rm = TRUE), sd.dusk = sd(mass.16, na.rm = TRUE),
            gain = mean(daily.gain, na.rm = TRUE), sd.gain = sd(daily.gain, na.rm = TRUE),
            loss = mean(night.loss, na.rm = TRUE), sd.loss = sd(night.loss, na.rm = TRUE),
            feeds = mean(total.pecks.day, na.rm = TRUE), sd.feeds = sd(total.pecks.day, na.rm = TRUE)) %>%
  ungroup()

# Descriptive stats by sex
descriptives.sex = descriptives.bird %>%
  group_by(sex) %>%
  summarise(dawn = mean(dawn.mass, na.rm = TRUE), sd.dawn = sd(dawn.mass, na.rm = TRUE),
            noon = mean(noon.mass, na.rm = TRUE), sd = sd(noon.mass, na.rm = TRUE),
            dusk = mean(dusk.mass, na.rm = TRUE), sd.dusk = sd(dusk.mass, na.rm = TRUE),
            daily.gain = mean(gain, na.rm = TRUE), sd.gain = sd(gain, na.rm = TRUE),
            night.loss = mean(loss, na.rm = TRUE), sd.loss = sd(loss, na.rm = TRUE),
            pecks = mean(feeds, na.rm = TRUE), sd.pecks = sd(feeds, na.rm = TRUE)) %>%
  ungroup()


descriptives.overall = descriptives.bird %>%
  summarise(dawn = mean(dawn.mass, na.rm = TRUE), sd.dawn = sd(dawn.mass, na.rm = TRUE),
            noon = mean(noon.mass, na.rm = TRUE), sd = sd(noon.mass, na.rm = TRUE),
            dusk = mean(dusk.mass, na.rm = TRUE), sd.dusk = sd(dusk.mass, na.rm = TRUE),
            daily.gain = mean(gain, na.rm = TRUE), sd.gain = sd(gain, na.rm = TRUE),
            night.loss = mean(loss, na.rm = TRUE), sd.loss = sd(loss, na.rm = TRUE),
            pecks = mean(feeds, na.rm = TRUE), sd.pecks = sd(feeds, na.rm = TRUE))



fit.results = fit.results%>%
  rename(dawn.mass = mass.8) %>%
  rename(dusk.mass = mass.16) %>%
  rename(noon.mass = mass.12)



# read in the biometric data for each bird


load("Biometric data.RData")



fit.results = merge(fit.results, b1, by = "bird")

descriptives.bird = merge(descriptives.bird, b1, by = "bird")

descriptives.bird = descriptives.bird %>%
  arrange(room, desc(sex), bird)




# Do we see the expected sex difference in skeletal size (tarsus length)?
m = lm(tarsus.mean ~ sex, data = descriptives.bird)
summary(m) # yes - although not quite significant in this sample



# Reliability of mass measures 

# dawn mass
mass.dawn = select(fit.results, "bird", "day.treat", "dawn.mass")
mass.dawn.wide = pivot_wider(data = mass.dawn,names_from = "day.treat",values_from = "dawn.mass")
mass.dawn.wide = select(mass.dawn.wide,"1","2","3","4","5","6","7","8","9","10","11")
icc(ratings = mass.dawn.wide, model = "twoway", type = "agreement",r0 = 0, conf.level = 0.95)
# This gives an ICC of 0.90 with 95% CI of 0.79 - 0.97 (i.e. good-excellent reliability)


# noon mass
mass.noon = select(fit.results, "bird", "day.treat", "noon.mass")
mass.noon.wide = pivot_wider(data = mass.noon,names_from = "day.treat",values_from = "noon.mass")
mass.noon.wide = select(mass.noon.wide,"1","2","3","4","5","6","7","8","9","10","11")
icc(ratings = mass.noon.wide, model = "twoway", type = "agreement",r0 = 0, conf.level = 0.95)
# This gives an ICC of 0.91 with 95% CI of 0.82 - 0.97 (i.e. good-excellent reliability)

# dusk mass
mass.dusk = select(fit.results, "bird", "day.treat", "dusk.mass")
mass.dusk.wide = pivot_wider(data = mass.dusk,names_from = "day.treat",values_from = "dusk.mass")
mass.dusk.wide = select(mass.dusk.wide,"1","2","3","4","5","6","7","8","9","10","11")
icc(ratings = mass.dusk.wide, model = "twoway", type = "agreement",r0 = 0, conf.level = 0.95)
# This gives an ICC of 0.89 with 95% CI of 0.79 - 0.97 (i.e. good-excellent reliability)

# daily gain (a measure of daily energy expenditure, EE)
gain = select(fit.results, "bird", "day.treat", "daily.gain")
gain.wide = pivot_wider(data = gain,names_from = "day.treat",values_from = "daily.gain")
gain.wide = select(gain.wide,"1","2","3","4","5","6","7","8","9","10","11")
icc(ratings = gain.wide, model = "twoway", type = "agreement",r0 = 0, conf.level = 0.95)
# This gives an ICC of 0.16 with 95% CI of 0.01 - 0.44 (i.e. poor reliability)

# night loss (a measure of sleeping MR)
loss = fit.results %>%
  filter(day.treat>1) %>% # no loss for the first day
  select("bird", "day.treat", "night.loss")
loss.wide = pivot_wider(data = loss,names_from = "day.treat",values_from = "night.loss")
loss.wide = select(loss.wide,"2","3","4","5","6","7","8","9","10","11")
icc(ratings = loss.wide, model = "twoway", type = "agreement",r0 = 0, conf.level = 0.95)
# This gives an ICC of 0.16 with 95% CI of 0.04 - 0.45 (i.e. poor reliability)

# Note that noon mass is the most reliable measurement of body mass






#### Do birds that eat more during the day weigh more at dusk? ####

# Analysis by bird day
m = lmer(dusk.mass ~ total.pecks.day + sex + (1|bird), data = fit.results, na.action = na.exclude)
summary(m) # dusk mass is predicted by the number of reinforcements earned in a day
anova(m)
plot(m)
hist(resid(m))
fit.results$fitted = fitted(m)

fit.results$predicted.dusk.mass = napredict(omit = m,fitted(m)) # this code pads the predictions to account for those missing due to NAs

Figure8a = ggplot(data = fit.results, aes(x = total.pecks.day, y = dusk.mass, colour = bird)) +
  geom_point() +
  #stat_smooth(method = "lm", aes(colour = bird), se = FALSE) +
  geom_line(data = fit.results[!is.na(fit.results$dusk.mass),], aes(x = total.pecks.day, y=predicted.dusk.mass, colour = bird)) +
  guides(colour = "none") +
  ylim(70,102) +
  xlab(expression(Reinforcements.day^-1)) +
  ylab("Dusk mass (g)") +
  annotate("text", x = 100, y = 102, label = "(a)", size = 5) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
Figure8a

# Is this a within-subjects effect?
# Do within-subject centering (van der Pol and Wright method) to separate within-and between-subject effects

d.bird = fit.results %>%
  group_by(bird) %>%
  summarise(mean.pecks.bird = mean(total.pecks.day, na.rm = TRUE)) 

fit.results = merge(fit.results,d.bird, by = "bird")
fit.results$centered.pecks = fit.results$total.pecks.day - fit.results$mean.pecks.bird

# First fit eqn 2 from van de Pol and Wright to decompose within and between subject effects
m = lmer(dusk.mass ~ mean.pecks.bird + centered.pecks + sex + (1|bird), data = fit.results)
summary(m) # this models suggests that the effect is all within subjects
anova(m)

# now fit eqn 3 to test whether the difference between within and between parameter estimates is significant
m = lmer(dusk.mass ~ total.pecks.day + mean.pecks.bird + sex + (1|bird), data = fit.results)
summary(m) # this models suggests a marginally non-significant difference beween the parameter estimates for within and between subject effects
anova(m)

# plot the within-subjects effect
Figure8b = ggplot(data = fit.results, aes(x = centered.pecks, y = dusk.mass)) +
  geom_point(aes(colour = bird)) +
  stat_smooth(method = "lm", colour = 'black') +
  guides(colour = "none") +
  ylim(70,102) +
  xlab(expression(Centered~reinforcements.day^-1)) +
  ylab("Dusk mass (g)") +
  annotate("text", x = -250, y = 102, label = "(b)", size = 5) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
Figure8b

# plot the between-subjects effect
d.mass = fit.results %>%
  group_by(bird) %>%
  summarise(mean.mass.bird = mean(dusk.mass, na.rm = TRUE)) 

d.bird = merge(d.bird, d.mass, by = "bird")

Figure8c = ggplot(data = d.bird, aes(x = mean.pecks.bird, y = mean.mass.bird)) +
  geom_point(aes(colour = bird)) +
  stat_smooth(method = "lm", colour = 'black') +
  scale_colour_discrete(name  ="Bird ID") +
  ylim(70,102) +
  xlab(expression(Mean~reinforcements.day^-1)) +
  ylab("Mean dusk mass (g)") +
  annotate("text", x = 150, y = 102, label = "(c)", size = 5) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
Figure8c

grid.arrange(Figure8a, Figure8b, Figure8c, ncol = 3, nrow = 1, widths = c(1,1,1.3))



