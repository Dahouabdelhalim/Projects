# R script for physiological data for the contrast study - have bodyweight, adrenal and thymus weight and then calculated
# relative to their bodyweight (as a percentage)

# specify the data 

Physiology$Treatment <- as.factor(Physiology$Treatment)

# next thing to do is to work out the means so then can plot these to visualise the data.

library(plyr)

#for the bodyweight 
treat_bodyweight <- ddply(Physiology, "Treatment", summarise, av=mean(Bodyweight), 
                   N = length(Bodyweight), SD = sd(Bodyweight), SE   = SD / sqrt(N))

# for adrenal weight 
treat_adrenal <- ddply(Physiology, "Treatment", summarise, av=mean(AdrenalWeight), 
                          N = length(AdrenalWeight), SD = sd(AdrenalWeight), SE   = SD / sqrt(N))

# for thymus weight 
treat_thymus <- ddply(Physiology, "Treatment", summarise, av=mean(ThymusWeight), 
                          N = length(ThymusWeight), SD = sd(ThymusWeight), SE   = SD / sqrt(N))

# adrenal relative to bodyweight
treat_Adr_bodyweight <- ddply(Physiology, "Treatment", summarise, av=mean(AdrenalBodyweight), 
                          N = length(AdrenalBodyweight), SD = sd(AdrenalBodyweight), SE   = SD / sqrt(N))

# Thymus relative to bodyweight
treat_thy_bodyweight <- ddply(Physiology, "Treatment", summarise, av=mean(ThyBodyweight), 
                              N = length(ThyBodyweight), SD = sd(ThyBodyweight), SE   = SD / sqrt(N))

# next lets plot the data to see what we expect from our analysis 

library(ggplot2)

# bodyweight
ggplot(treat_bodyweight, aes(x=Treatment, y=av)) +geom_bar(stat = 'identity') +
  theme_bw() +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlab("Treatment") + ylab("Bodyweight" ) + geom_errorbar(aes(ymin=av-SE, ymax=av+SE), width=.1)

#Adrenal weight 
ggplot(treat_adrenal, aes(x=Treatment, y=av)) +geom_bar(stat = 'identity') +
  theme_bw() +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlab("Treatment") + ylab("Adrenal Weight" ) + geom_errorbar(aes(ymin=av-SE, ymax=av+SE), width=.1)

# thymus weight 
ggplot(treat_thymus, aes(x=Treatment, y=av)) +geom_bar(stat = 'identity') +
  theme_bw() +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlab("Treatment") + ylab("Thymus Weight" ) + geom_errorbar(aes(ymin=av-SE, ymax=av+SE), width=.1)

# adrenal relative to bodyweight 
ggplot(treat_Adr_bodyweight, aes(x=Treatment, y=av)) +geom_bar(stat = 'identity') +
  theme_bw() +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlab("Treatment") + ylab("Adrenal Weight/Bodyweight" ) + geom_errorbar(aes(ymin=av-SE, ymax=av+SE), width=.1)

# thymus relative to bodyweight 
ggplot(treat_thy_bodyweight, aes(x=Treatment, y=av)) +geom_bar(stat = 'identity') +
  theme_bw() +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlab("Treatment") + ylab("Thymus Weight/Bodyweight" ) + geom_errorbar(aes(ymin=av-SE, ymax=av+SE), width=.1)

#lets look at the raw data to see how the distribution looks 

library(ggpubr)

# Bodyweight
ggdensity(Physiology$Bodyweight, main = "density plot of Bodyweight", xlab = "Bodyweight")
ggqqplot(Physiology$Bodyweight) 

# adrenal weight 
ggdensity(Physiology$AdrenalWeight, main = "density plot of Adrenal weight", xlab = "Adrenal Weight")
ggqqplot(Physiology$AdrenalWeight) # looks good. - therefore for publication I think we should present data for this. 

# thymus weight 
ggdensity(Physiology$ThymusWeight, main = "density plot of Thymus weight", xlab = "Thymus Weight")
ggqqplot(Physiology$ThymusWeight) # looks okay

# adrenal relative to bodyweight 
ggdensity(Physiology$AdrenalBodyweight, main = "density plot of Adrenal weight/Bodyweight", xlab = "Adrenal Weight/Bodyweight")
ggqqplot(Physiology$AdrenalBodyweight) # has clear points - biphasic peaks 

# thymus relative to bodyweight 
ggdensity(Physiology$ThyBodyweight, main = "density plot of Thymus weight/Bodyweight", xlab = "Thymus Weight/Bodyweight")
ggqqplot(Physiology$ThyBodyweight) # looks ok 

# independent samples t tests 

#Bodyweight
t.test(Bodyweight~Treatment, data=Physiology) # no difference

#Welch Two Sample t-test

#data:  Bodyweight by Treatment
#t = -1.1909, df = 29.758, p-value = 0.2431
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -1.8159647  0.4784647
#sample estimates:
#  mean in group Tail mean in group Tunnel 
#25.80625             26.47500 

# adrenal weight 
t.test(AdrenalWeight~Treatment, data=Physiology) # this is significantly different 

#Welch Two Sample t-test

#data:  AdrenalWeight by Treatment
#t = 2.3653, df = 28.002, p-value = 0.02517
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
 # 0.181732 2.530768
#sample estimates:
#  mean in group Tail mean in group Tunnel 
#7.32500              5.96875 

# thymus weight 
t.test(ThymusWeight~Treatment, data=Physiology) # non significant 

#Welch Two Sample t-test

#data:  ThymusWeight by Treatment
#t = 1.3432, df = 29.991, p-value = 0.1893
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -3.3443 16.1943
#sample estimates:
#  mean in group Tail mean in group Tunnel 
#53.10625             46.68125 

# adrenal weight realtive to bodyweight 
t.test(AdrenalBodyweight~Treatment, data=Physiology) # significantly different 

#Welch Two Sample t-test

#data:  AdrenalBodyweight by Treatment
#t = 2.4974, df = 27.642, p-value = 0.01875
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  0.001232768 0.012517232
#sample estimates:
#  mean in group Tail mean in group Tunnel 
#0.028750             0.021875

# thymus relative to bodyweight 
t.test(ThyBodyweight~Treatment, data=Physiology) # non significant 

#Welch Two Sample t-test

#data:  ThyBodyweight by Treatment
#t = 1.4644, df = 29.744, p-value = 0.1536
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.01086597  0.06586597
#sample estimates:
#  mean in group Tail mean in group Tunnel 
#0.205625             0.178125 



