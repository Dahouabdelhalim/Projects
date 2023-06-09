# R analysis for the preshift data for both consumption and lick cluster size.... for each concentration of sucrose. 

# want to do a 2x2 anova looking at the difference between sucrose concentration and handling method. 

# specify the data
PreshiftCons$Treatment <- as.factor(PreshiftCons$Treatment)
PreshiftCons$Conc <- as.factor(PreshiftCons$Conc)
PreshiftCons$Consumption <- as.numeric(PreshiftCons$Consumption)

# consumption data first

library(lme4)

m0 <- lm((Consumption) ~ Treatment + Conc + Treatment: Conc, data = PreshiftCons)
m0
summary(m0)

hist(resid(m0)) # looking at the fit of the residuals this seems like a nice fit.
plot(m0)

# remove the interaction 

m1 <- lm((Consumption) ~ Treatment + Conc, data = PreshiftCons)
anova(m0, m1) #F=0.144, p=0.7061

# remove treatment

m2 <- lm((Consumption) ~ Conc, data = PreshiftCons)
anova(m1, m2) #F=1.6263, p=0.207

#remove concentration 
m3 <- lm((Consumption) ~ Treatment, data = PreshiftCons)
anova(m1, m3) # F=4.3867, p=0.04038

# lets plot the data .... 
library(ggplot2)

#by handling x conc
mean_Cons_Hand_Conc <- ddply(PreshiftCons, c("Treatment", "Conc"), summarise, av=mean(Consumption, na.rm = TRUE),
                             N = length(Consumption-1), SD=sd(Consumption, na.rm = TRUE), SE=SD/sqrt(N))

ggplot(data=mean_Cons_Hand_Conc, aes(x=Conc, y=av, group=Treatment)) +
  xlab("Conc") + ylab("Consumption" ) +
  geom_errorbar(aes(ymin=av-SE, ymax=av+SE), width=.1) +
  geom_line(aes(color=Treatment))+
  geom_point(aes(color=Treatment))

# same for lick cluster size 250

m4 <- lm((LCS250) ~ Treatment + Conc + Treatment: Conc, data = PreshiftLCS)
m4
summary(m4)

hist(resid(m4)) # looking at the fit of the residuals this seems like a nice fit.
plot(m4)

# remove the interaction 

m5 <- lm((LCS250) ~ Treatment + Conc, data = PreshiftLCS)
anova(m4, m5) #F=0.8291, p=0.3662

# remove treatment 

m6 <- lm((LCS250) ~ Conc, data = PreshiftLCS)
anova(m5, m6) # marginal F=3.1361, p=0.08157

# remove conc 

m7 <- lm((LCS250) ~ Treatment, data = PreshiftLCS)
anova(m5, m7) # significant F=4.2818, p=0.04277

# same for lick cluster size 500

m8 <- lm((LCS500) ~ Treatment + Conc + Treatment: Conc, data = PreshiftLCS)
m8
summary(m8)

hist(resid(m8)) # looking at the fit of the residuals this seems like a nice fit.
plot(m8)

# remove the interaction 

m9 <- lm((LCS500) ~ Treatment + Conc, data = PreshiftLCS)
anova(m8, m9) #F=1.3891, p=0.2432

# remove treatment 

m10 <- lm((LCS500) ~ Conc, data = PreshiftLCS)
anova(m9, m10) # marginal F=6.7744, p=0.01159

# remove conc 

m11 <- lm((LCS500) ~ Treatment, data = PreshiftLCS)
anova(m9, m11) # significant F=4.0706, p=0.04804

# same for lick cluster size 1000

m12 <- lm((LCS1000) ~ Treatment + Conc + Treatment: Conc, data = PreshiftLCS)
m12
summary(m12)

hist(resid(m12)) # looking at the fit of the residuals this seems like a nice fit.
plot(m12)

# remove the interaction 

m13 <- lm((LCS1000) ~ Treatment + Conc, data = PreshiftLCS)
anova(m12, m13) #F=2.1889, p=0.1442

# remove treatment 

m14 <- lm((LCS1000) ~ Conc, data = PreshiftLCS)
anova(m13, m14) # marginal F=9.6494, p=0.002874

# remove conc 

m15 <- lm((LCS1000) ~ Treatment, data = PreshiftLCS)
anova(m13, m15) # significant F=2.0205 p=0.1603
