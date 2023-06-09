# R analysis for the open field data collected for the contrast study. 

# specify the data 

OF_data <- as.factor(OF_data$Treatment)

# next want to calculate the means in order to see what the data looks likw ehen plotted. 

library (plyr)

#by treatment for duration spent in the centre 
treat_Dur <- ddply(OF_data, "Treatment", summarise, av=mean(DurationCentre), 
                       N = length(DurationCentre), SD = sd(DurationCentre), SE   = SD / sqrt(N))

# by treatment for frequency in the centre 
treat_freq <- ddply(OF_data, "Treatment", summarise, av=mean(FreqCentre), 
                   N = length(FreqCentre), SD = sd(FreqCentre), SE   = SD / sqrt(N))

# percentage of time spent in the centre 
treat_perc <- ddply(OF_data, "Treatment", summarise, av=mean(PercentageCentre), 
                    N = length(PercentageCentre), SD = sd(PercentageCentre), SE   = SD / sqrt(N))

# movement 
treat_move <- ddply(OF_data, "Treatment", summarise, av=mean(Movement), 
                    N = length(Movement), SD = sd(Movement), SE   = SD / sqrt(N))

# velocity
treat_velocity <- ddply(OF_data, "Treatment", summarise, av=mean(Velocity), 
                                   N = length(Velocity), SD = sd(Velocity), SE   = SD / sqrt(N))

# distance
treat_distance <- ddply(OF_data, "Treatment", summarise, av=mean(Distance), 
                                   N = length(Distance), SD = sd(Distance), SE   = SD / sqrt(N))

# next lets plot the data to see what we expect from our analysis 

library(ggplot2)

# duration in centre 
ggplot(treat_Dur, aes(x=Treatment, y=av)) +geom_bar(stat = 'identity') +
  theme_bw() +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlab("Treatment") + ylab("Duration in centre" ) + geom_errorbar(aes(ymin=av-SE, ymax=av+SE), width=.1)

# frequency in centre 
ggplot(treat_freq, aes(x=Treatment, y=av)) +geom_bar(stat = 'identity') +
  theme_bw() +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlab("Treatment") + ylab("Frequency in centre" ) + geom_errorbar(aes(ymin=av-SE, ymax=av+SE), width=.1)

# percentage of time spent in the centre 
ggplot(treat_perc, aes(x=Treatment, y=av)) +geom_bar(stat = 'identity') +
  theme_bw() +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlab("Treatment") + ylab("Percentage of time in centre" ) + geom_errorbar(aes(ymin=av-SE, ymax=av+SE), width=.1)

# movement 
ggplot(treat_move, aes(x=Treatment, y=av)) +geom_bar(stat = 'identity') +
  theme_bw() +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlab("Treatment") + ylab("Total time spent moving" ) + geom_errorbar(aes(ymin=av-SE, ymax=av+SE), width=.1)

# velocity
ggplot(treat_velocity, aes(x=Treatment, y=av)) +geom_bar(stat = 'identity') +
  theme_bw() +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlab("Treatment") + ylab("Mean velocity (cm/s)" ) + geom_errorbar(aes(ymin=av-SE, ymax=av+SE), width=.1)

# distance travelled
ggplot(treat_distance, aes(x=Treatment, y=av)) +geom_bar(stat = 'identity') +
  theme_bw() +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlab("Treatment") + ylab("Distance Travelled" ) + geom_errorbar(aes(ymin=av-SE, ymax=av+SE), width=.1)

#lets look at the raw data to see how the distribution looks 

library(ggpubr)

# Duration in the centre
ggdensity(OF_data$DurationCentre, main = "density plot of Duration in centre", xlab = "Duration in centre")
ggqqplot(OF_data$DurationCentre)  # doesnt look too bad. 

# freq in centre 
ggdensity(OF_data$FreqCentre, main = "density plot of Freq in centre", xlab = "Freq in centre")
ggqqplot(OF_data$FreqCentre) # again doesnt look too bad 

# percentage of time in centre 
ggdensity(OF_data$PercentageCentre, main = "density plot of percentage of time in centre", xlab = "Percentage in centre")
ggqqplot(OF_data$PercentageCentre) # not too bad 

# movement 
ggdensity(OF_data$Movement, main = "density plot of time spent moving", xlab = "Movement")
ggqqplot(OF_data$Movement)

# velocity 
ggdensity(OF_data$Velocity, main = "density plot of velocity", xlab = "Velocity")
ggqqplot(OF_data$Velocity)

#distance 
ggdensity(OF_data$Distance, main = "density plot of Distance travelled", xlab = "Distance")
ggqqplot(OF_data$Distance)

# I think all this data looks relatively normally distributed therefore I think it makes sense to run independent samples
# t tests. 

# Duration in centre 
t.test(DurationCentre~ Treatment, data = OF_data) # significantly different 

#Welch Two Sample t-test

#data:  DurationCentre by Treatment
#t = -2.9442, df = 58.671, p-value = 0.004635
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -38.642155  -7.367845
#sample estimates:
#  mean in group Tail mean in group Tunnel 
#99.23094            122.23594

# frequency in the centre 
t.test(FreqCentre~ Treatment, data = OF_data) # sig different

#Welch Two Sample t-test

#data:  FreqCentre by Treatment
#t = -2.3628, df = 61.554, p-value = 0.02131
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -20.942095  -1.745405
#sample estimates:
#  mean in group Tail mean in group Tunnel 
#76.75000             88.09375 

# percentage of time in the centre 
t.test(PercentageCentre~Treatment, data=OF_data)# significantly different 

#Welch Two Sample t-test

#data:  PercentageCentre by Treatment
#t = -2.9162, df = 58.66, p-value = 0.005013
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -6.390918 -1.189082
#sample estimates:
#  mean in group Tail mean in group Tunnel 
#16.67312             20.46313

# movement 
t.test(Movement~Treatment, data=OF_data) # not different 

#Welch Two Sample t-test

#data:  Movement by Treatment
#t = 0.11338, df = 61.969, p-value = 0.9101
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -13.57034  15.20221
#sample estimates:
#  mean in group Tail mean in group Tunnel 
#260.7044             259.8884 

# velocity 
t.test(Velocity~Treatment, data=OF_data) # not significant (but relatively marginal)

#Welch Two Sample t-test

#data:  Velocity by Treatment
#t = -1.7576, df = 60.528, p-value = 0.08387
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -1.716985  0.110735
#sample estimates:
#  mean in group Tail mean in group Tunnel 
#10.86063             11.66375 

# distance 
t.test(Distance~Treatment, data=OF_data) # non significant but marginal

#Welch Two Sample t-test

#data:  Distance by Treatment
#t = -1.8565, df = 60.824, p-value = 0.06822
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -1056.21719    39.22156
#sample estimates:
#  mean in group Tail mean in group Tunnel 
#6457.169             6965.667 



