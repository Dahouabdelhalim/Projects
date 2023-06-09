# R analysis of contrast study data - excluded 7 animals who fell off the maze before the end of the test.

# specify the data 

EPM$Treatment <- as.factor(EPM$Treatment)


# summarise the data - calculating the means for each factor so then can plot these and see what the data looks like and 
# what you might expect the analysis to look like.

library (plyr)

#by treatment for Open arm entries 
treat_entries <- ddply(EPM, "Treatment", summarise, av=mean(OpenArmEntries), 
                   N = length(OpenArmEntries), SD = sd(OpenArmEntries), SE   = SD / sqrt(N))


#by treatment for duration on open arms 
treat_duration <- ddply(EPM, "Treatment", summarise, av=mean(DurationOpenArms), 
                       N = length(DurationOpenArms), SD = sd(DurationOpenArms), SE   = SD / sqrt(N))

#by treatment for percentage on open arms 
treat_duration_percent <- ddply(EPM, "Treatment", summarise, av=mean(PercentOpenArms), 
                        N = length(PercentOpenArms), SD = sd(PercentOpenArms), SE   = SD / sqrt(N))

#by treatment for PSA
treat_PSA <- ddply(EPM, "Treatment", summarise, av=mean(PSA), 
                        N = length(PSA), SD = sd(PSA), SE   = SD / sqrt(N))

# plot the data to see what it looks like

library(ggplot2)

# open arm entries for each treatment 
ggplot(treat_entries, aes(x=Treatment, y=av)) +geom_bar(stat = 'identity') +
  theme_bw() +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlab("Treatment") + ylab("Open Arm Entries" ) + geom_errorbar(aes(ymin=av-SE, ymax=av+SE), width=.1)

# duration on open arms for each treatment 
ggplot(treat_duration, aes(x=Treatment, y=av)) +geom_bar(stat = 'identity') +
  theme_bw() +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlab("Treatment") + ylab("Duration on Open Arms" ) + geom_errorbar(aes(ymin=av-SE, ymax=av+SE), width=.1)

# percentage of time on the open arms for each treatment - looks identical to previous
ggplot(treat_duration_percent, aes(x=Treatment, y=av)) +geom_bar(stat = 'identity') +
  theme_bw() +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlab("Treatment") + ylab("Percentage of total time on open arms" ) + geom_errorbar(aes(ymin=av-SE, ymax=av+SE), width=.1)

# Number of protected stretch attend postures for both treatments 
ggplot(treat_PSA, aes(x=Treatment, y=av)) +geom_bar(stat = 'identity') +
  theme_bw() +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlab("Treatment") + ylab("Number of PSA" ) + geom_errorbar(aes(ymin=av-SE, ymax=av+SE), width=.1)

#lets look at the raw data to see how the distribution looks 

library(ggpubr)

# Open arm entries 
ggdensity(EPM$OpenArmEntries, main = "density plot of Open arm entries", xlab = "Open arm entries")
ggqqplot(EPM$OpenArmEntries) # distribution doesnt look too bad - QQ plot looks okay and the distrubution seems slightly 
# positively skewed. 

# Duration on open arms 
ggdensity(EPM$DurationOpenArms, main = "density plot of duration open arms", xlab = "Duration Open arms")
ggqqplot(EPM$DurationOpenArms) # again the distribution doesnt seem too bad or far from normality. slight deviation on qq plot
# at higher values

# percentage of time on open arms 
ggdensity(EPM$PercentOpenArms, main = "density plot of percentage open arms", xlab = "Percentage Open arms")
ggqqplot(EPM$PercentOpenArms) # again very similar to above - slight deviation at the higher values. 

# PSA 
ggdensity(EPM$PSA, main = "density plot of PSA", xlab = "PSA")
ggqqplot(EPM$PSA) # this data is higly positively skewed - will want to try and transform this to try and make it less skewed. 

# lets try a square root transformation

sqrt_PSA <- sqrt(EPM$PSA)
ggdensity(sqrt_PSA, main = "density plot of sqrt PSA", xlab = "Sqrt PSA") # looks better on the density plot 
ggqqplot(sqrt_PSA)

# The data seems relatively normally distributed across each of these DVs and therefore I think we should just do simple t tests
# make sure do the t test on the square root data for the PSA. 

# Independent Samples t test 

# Open Arm Entries 
t.test(OpenArmEntries~ Treatment, data = EPM) # significantly different 

#Welch Two Sample t-test
#
#data:  OpenArmEntries by Treatment
#t = -3.2104, df = 54.121, p-value = 0.002231
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -5.263510 -1.216785
#sample estimates:
#  mean in group Tail mean in group Tunnel 
#7.035714            10.275862 

# duration on open arms 
t.test(DurationOpenArms~Treatment, data =EPM) # significantly different

#Welch Two Sample t-test
#
#data:  DurationOpenArms by Treatment
#t = -3.6067, df = 53.6, p-value = 0.0006811
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -31.487617  -8.985659
#sample estimates:
#  mean in group Tail mean in group Tunnel 
#37.35750             57.59414 

# percentage of total time spent on open arms 
t.test(PercentOpenArms~Treatment, data = EPM) # significantly different 

#Welch Two Sample t-test
#
#data:  PercentOpenArms by Treatment
#t = -3.6068, df = 53.598, p-value = 0.000681
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -10.495257  -2.995112
#sample estimates:
#  mean in group Tail mean in group Tunnel 
#12.44964             19.19483 

# PSA - remember the stats are done on the squart root transformed data 
t.test(sqrt_PSA~Treatment, data=EPM) # significantly different 

#Welch Two Sample t-test

#data:  sqrt_PSA by Treatment
#t = 2.0574, df = 44.412, p-value = 0.04555
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  0.01032385 0.98773297
#sample estimates:
#  mean in group Tail mean in group Tunnel 
#2.214479             1.715451 
