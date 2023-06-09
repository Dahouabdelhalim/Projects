# analysis looking at the control mice consumption data for the postshift phase across the two 
# postshiftphases (in line with reviewer comments)
# want to do a generalised linear mixed model to determined whether we have a handling effect and whether the tail handled mice 
# are anhedonic across both postshift phases.

# need to remmeber this data is just for the control mice not the mice undergoing contrast effects 
# so we have a n=8 for each group 

# datafile = ConsumptionControlsPostshift.csv

# specify the data
ConsumptionControlsPostshift$Handling.Method <- as.factor(ConsumptionControlsPostshift$Handling.Method)
ConsumptionControlsPostshift$Mouse.ID <- as.factor(ConsumptionControlsPostshift$Mouse.ID)
ConsumptionControlsPostshift$SucroseConc <- as.factor(ConsumptionControlsPostshift$SucroseConc)
ConsumptionControlsPostshift$PostshiftPhase <- as.factor(ConsumptionControlsPostshift$PostshiftPhase)

# lets run a genearlised linear mixed model 

library(lme4)

m1.0 <- lmer((Consumption) ~ Handling.Method + SucroseConc + PostshiftPhase + Handling.Method: SucroseConc + 
               + Handling.Method:PostshiftPhase + PostshiftPhase:SucroseConc + Handling.Method: SucroseConc: PostshiftPhase
             + (1|Mouse.ID), data=ConsumptionControlsPostshift)

summary(m1.0)

hist(resid(m1.0)) # looking at the fit of the residuals this seems like a good fit.
plot(m1.0) 

# lets remove each variable to determine their effect.

#lets remove the three way interaction

m2.0 <- lmer((Consumption) ~ Handling.Method + SucroseConc + PostshiftPhase + Handling.Method: SucroseConc + 
               + Handling.Method:PostshiftPhase + PostshiftPhase:SucroseConc 
             + (1|Mouse.ID), data=ConsumptionControlsPostshift)
anova(m1.0, m2.0)

# Df     AIC     BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# m2.0  9 -223.37 -191.50 120.69  -241.37                         
# m1.0 10 -221.51 -186.09 120.75  -241.51 0.1333      1      0.715

# lets remove the Handling method x conc interaction
m3.0 <- lmer((Consumption) ~ Handling.Method + SucroseConc + PostshiftPhase + 
               + Handling.Method:PostshiftPhase + PostshiftPhase:SucroseConc 
             + (1|Mouse.ID), data=ConsumptionControlsPostshift)
anova(m1.0, m3.0)

# Df     AIC     BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# m3.0  8 -224.13 -195.80 120.06  -240.13                         
# m1.0 10 -221.51 -186.09 120.75  -241.51 1.3771      2     0.5023

# lets remove the handling method x postshift phase interaction

m4.0 <- lmer((Consumption) ~ Handling.Method + SucroseConc + PostshiftPhase + 
               + PostshiftPhase:SucroseConc + Handling.Method:SucroseConc 
             + (1|Mouse.ID), data=ConsumptionControlsPostshift)
anova(m1.0, m4.0)

# Df     AIC     BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# m4.0  8 -225.16 -196.84 120.58  -241.16                         
# m1.0 10 -221.51 -186.09 120.75  -241.51 0.3409      2     0.8433

# lets remove the postshiftphase x sucrose conc interaction 
m5.0 <- lmer((Consumption) ~ Handling.Method + SucroseConc + PostshiftPhase + 
               + Handling.Method:PostshiftPhase + Handling.Method:SucroseConc 
             + (1|Mouse.ID), data=ConsumptionControlsPostshift)
anova(m1.0, m5.0)

# Df     AIC     BIC logLik deviance  Chisq Chi Df Pr(>Chisq)
# m5.0  8 -225.07 -196.74 120.54  -241.07                         
# m1.0 10 -221.51 -186.09 120.75  -241.51 0.4334      2     0.8052

# lets remove main effects - handling method - lets create a model with just the main effects included
m6.0 <- lmer((Consumption) ~ Handling.Method + SucroseConc + PostshiftPhase + (1|Mouse.ID), data=ConsumptionControlsPostshift)

# lets remove handling method 
m6.1 <- lmer((Consumption) ~ SucroseConc + PostshiftPhase + (1|Mouse.ID), data=ConsumptionControlsPostshift)
anova(m6.0, m6.1)

# Df     AIC     BIC logLik deviance  Chisq Chi Df Pr(>Chisq)  
# m6.1  5 -224.18 -206.48 117.09  -234.18                           
# m6.0  6 -227.63 -206.38 119.81  -239.63 5.4443      1    0.01963 *

# lets remove sucrose concentration 
m6.2 <- lmer((Consumption) ~ Handling.Method + PostshiftPhase + (1|Mouse.ID), data=ConsumptionControlsPostshift)
anova(m6.0, m6.2)

# Df     AIC     BIC logLik deviance  Chisq Chi Df Pr(>Chisq)   
# m6.2  5 -221.24 -203.54 115.62  -231.24                            
# m6.0  6 -227.63 -206.38 119.81  -239.63 8.3853      1   0.003783 **

# lets remove postshift phase
m6.3 <- lmer((Consumption) ~ Handling.Method + SucroseConc + (1|Mouse.ID), data=ConsumptionControlsPostshift)
anova(m6.0, m6.3)

# Df     AIC     BIC logLik deviance  Chisq Chi Df Pr(>Chisq)  
# m6.3  5 -226.54 -208.84 118.27  -236.54                           
# m6.0  6 -227.63 -206.38 119.81  -239.63 3.0833      1     0.0791 .

# next lets try and visualise the data to see what it looks like ... 

library(ggplot2)

#by time and treatment 
mean_Cons_Hand_Conc_PostshiftPhase <- ddply(ConsumptionControlsPostshift, c("Handling.Method", "SucroseConc","PostshiftPhase"), summarise, av=mean(Consumption, na.rm = TRUE),
                             N = length(Consumption-1), SD=sd(Consumption, na.rm = TRUE), SE=SD/sqrt(N))


ggplot(data=mean_Cons_Hand_Conc_PostshiftPhase, aes(x=PostshiftPhase, y=av, group=SucroseConc)) +
  xlab("Postshift Phase") + ylab("Consumption" ) +
  geom_errorbar(aes(ymin=av-SE, ymax=av+SE), width=.1) +
  geom_line(aes(color=SucroseConc))+
  geom_point(aes(color=SucroseConc)) +
  facet_wrap(~Handling.Method)
  