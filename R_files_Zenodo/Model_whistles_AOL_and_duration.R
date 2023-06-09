## Manuscript: "Anthropogenic noise impairs cooperation in bottlenose dolphins"
## by SÃ¸rensen, Haddock, Guarino, Jaakkola, McMullen, Jensen, Tyack and King
##
##
##
## Citations:
## Florian Hartig (2022). DHARMa: Residual Diagnostics for Hierarchical
## (Multi-Level / Mixed) Regression Models. R package version 0.4.5.
## http://florianhartig.github.io/DHARMa/
##   
#### Last updated: 14-11-2022
################################################################################
# Clear environment
rm(list = ls())

# Load packages
library(glmmTMB) 
library(ggdist)
library(DHARMa) 
library(lmerTest) 
library(plyr)
library(effects)
library(ggplot2)

#### Load data ####
my_data <- read.csv("Data_Whistle_AOL_and_duration.csv", header=T) # Load data
attach(my_data)
str(my_data)
my_data$treatment <- as.factor(my_data$treatment)                                                     

## Subset data for Delta and Reese:
Reese = my_data[animal == 1, ]                     # Reese: 150 whistles
Delta = my_data[animal == 2, ]                     # Delta: 101 whistles

########################################################################################################
#### Delta whistle AOL ~ NL (as measured just before whistle emission) ####
########################################################################################################
# Run LMM with session number (1-20) and trial (accounting for whistles being repeated within a trial) as
# random effects:
Delta_AOL_full <- lmer(AOL_across_harmonics ~ NL + (1|session) + (1|trial), data = Delta) 
Delta_AOL_null <- lmer(AOL_across_harmonics ~  1 + (1|session) + (1|trial), data = Delta)
anova(Delta_AOL_full, Delta_AOL_null) 
summary(Delta_AOL_full) 
confint(Delta_AOL_full)

# Check model fit:
simulationOutput<-simulateResiduals(fittedModel = Delta_AOL_full, plot=T)

# Check diagnostic plots: 
par(mfrow=c(1,3))
hist(residuals(Delta_AOL_full))
qqnorm(residuals(Delta_AOL_full)) 
qqline(residuals(Delta_AOL_full))
plot(residuals(Delta_AOL_full) ~ fitted(Delta_AOL_full))

########################################################################################################
#### Reese whistle AOL ~ NL (as measured just before whistle emission) ####
########################################################################################################
# Run LMM with session number (1-20) and trial (accounting for the emission of more whistles within a trial) as
# random effects:
Reese_AOL_full <- lmer(AOL_across_harmonics ~ NL + (1|session) + (1|trial), data = Reese)
Reese_AOL_null <- lmer(AOL_across_harmonics ~  1 + (1|session) + (1|trial), data = Reese)
anova(Reese_AOL_full, Reese_AOL_null)
summary(Reese_AOL_full) 
confint(Reese_AOL_full)

# Checking modelfit: 
simulationOutput<-simulateResiduals(fittedModel = Reese_AOL_full, plot=T) 

# Diagnostic plots: 
par(mfrow=c(1,3))
hist(residuals(Reese_AOL_full))
qqnorm(residuals(Reese_AOL_full)) 
qqline(residuals(Reese_AOL_full))
plot(residuals(Reese_AOL_full) ~ fitted(Reese_AOL_null))

# The qq diagnostic plot can be improved by removing outliers identified by Cook's distance, which 
# results in an even more significant relationship between NL and whistle AOL.
# However, given that the model fit using DHARMa simulateResiduals looks good, we have presented the
# results for Reese_AOL_full in the manuscript (i.e. without removing outliers).

################################################################################
#### Reese whistle duration ~ NL (as measured just before whistle emission) ####
################################################################################
# Run LMM with session number (1-20) and trial (accounting for whistles being repeated within a trial) as
# random effects:

# For visualising in plot:
Reese_dur_full_not_log <- lmer(whistle_dur ~ NL + (1|session) + (1|trial), data = Reese) 

# log transformed whistle duration as normality assumption was violated
Reese2 <- Reese 
Reese2$whistle_dur = log(Reese$whistle_dur)
Reese_dur_full <- lmer(whistle_dur ~ NL + (1|session) + (1|trial), data = Reese2)
Reese_dur_null <- lmer(whistle_dur ~ 1  + (1|session) + (1|trial), data = Reese2) 
anova(Reese_dur_null, Reese_dur_full) ; 
summary(Reese_dur_full) ; 
confint(Reese_dur_full) ;

# Check modelfit
simulationOutput<-simulateResiduals(fittedModel = Reese_dur_full, plot=T) 

# Diagnositic plots
par(mfrow=c(1,3))
hist(residuals(Reese_dur_full))
qqnorm(residuals(Reese_dur_full)) 
qqline(residuals(Reese_dur_full))
plot(residuals(Reese_dur_full) ~ fitted(Reese_dur_full))

################################################################################
#### Delta whistle duration ~ NL (as measured just before whistle emission) ####
################################################################################
# Run LMM with session number (1-20) and trial (accounting for whistles being repeated within a trial) as
# random effects:

# For visualising in plot:
Delta_dur_full_not_log <- lmer(whistle_dur ~ NL + (1|session) + (1|trial), data = Delta)

# log transformed whistle duration as normality assumption was violated
Delta2 <- Delta 
Delta2$whistle_dur <- log(Delta$whistle_dur)
Delta_dur_full     <- lmer(whistle_dur ~ NL + (1|session) + (1|trial), data = Delta2)
Delta_dur_null     <- lmer(whistle_dur ~  1 + (1|session) + (1|trial), data = Delta2) 
anova(Delta_dur_null, Delta_dur_full) ; 
summary(Delta_dur_full) ; 
confint(Delta_dur_full) ;

# Check modelfit
simulationOutput<-simulateResiduals(fittedModel = Delta_dur_full, plot=T) 

# Diagnostic plots
par(mfrow=c(1,3))
hist(residuals(Delta_dur_full))
qqnorm(residuals(Delta_dur_full)) 
qqline(residuals(Delta_dur_full))
plot(residuals(Delta_dur_full) ~ fitted(Delta_dur_full))

################################################################################
#### Settings for plotting AOL figures ####
################################################################################
## Reese
R_effects_NL_pred   <- effects::predictorEffect(predictor = "NL", mod = Reese_AOL_full)
summary(R_effects_NL_pred) ; 
Reese_NL_pred    <- as.data.frame(R_effects_NL_pred)

## Delta
D_effects_NL_pred   <- effects::predictorEffect(predictor = "NL", mod = Delta_AOL_full)
summary(D_effects_NL_pred) ;
Delta_NL_pred    <- as.data.frame(D_effects_NL_pred)

# Set colours
my_colors <- RColorBrewer::brewer.pal(6, "Blues")[2:7]

################################################################################
#### Scatter-plot: Reese AOL vs NL ####
################################################################################
pp1 <- ggplot() + geom_point(data = Reese, aes(NL, AOL_across_harmonics, group = as.factor(treatment), color=as.factor(treatment)), size=5, alpha=1) +
  
  scale_colour_manual(values=my_colors) +
  
  geom_point() + xlim(95, 160) + ylim(120,172) + 
  
  geom_line(data=Reese_NL_pred, aes(x=NL, y=fit), color="#5B8FA8FF",lwd=2) +
  
  geom_ribbon(data= Reese_NL_pred, aes(x=NL, ymin=lower, ymax=upper), alpha= 0.3, fill="#5B8FA8FF") +
  
  labs(x="", y="") + 
  
  theme(legend.position = "none", axis.text = element_text(size = 25), panel.grid = element_blank(), panel.background = element_blank()
        ,axis.line = element_line(colour = "black"), axis.title.x = element_blank(), text = element_text(size=28), panel.border = element_rect(fill=NA)) +
  
#  scale_color_brewer(palette = "Blues", direction = 1) +

  scale_y_continuous(expand = c(0,0), limits = c(120, 172)) 
pp1
################################################################################
#### Scatter-plot: Delta AOL vs NL ####
################################################################################
pp2 <- ggplot() + geom_point(data = Delta, aes(NL, AOL_across_harmonics, group=as.factor(treatment), color=as.factor(treatment)), size=5, alpha=1) +
  
  scale_colour_manual(values=my_colors) +
  
  geom_point() + xlim(95, 160) + geom_point() + ylim(120,172) +
  
  geom_line(data = Delta_NL_pred, aes(x=NL, y=fit), color="#5B8FA8FF",lwd=2) +
  
  geom_ribbon(data = Delta_NL_pred , aes(x=NL, ymin=lower, ymax=upper), alpha= 0.3, fill="#5B8FA8FF") +
  
    labs(x="ambient NL (dB re 1 ?Pa, rms)", y="") + 
  
  theme(legend.position = "none", axis.text = element_text(size = 25), panel.grid = element_blank(), panel.background = element_blank()
        ,axis.line = element_line(colour = "black"), axis.title.x = element_blank(), text = element_text(size=28), panel.border = element_rect(fill=NA)) +
  
  # scale_color_brewer(palette = "Blues", direction = 1) +
  
  scale_y_continuous(expand = c(0,0), limits = c(120, 172)) 
pp2
################################################################################
#### Settings for plotting duration figures ####
################################################################################
# Reese
R_effects_dur_NL_pred <- effects::predictorEffect(predictor = "NL", mod = Reese_dur_full_not_log) ;
summary(R_effects_dur_NL_pred) ;
Reese_dur_NL_pred     <- as.data.frame(R_effects_dur_NL_pred) ;

# Delta
D_effects_dur_NL_pred <- effects::predictorEffect(predictor = "NL", mod = Delta_dur_full_not_log) ;
summary(D_effects_dur_NL_pred) ;
Delta_dur_NL_pred     <- as.data.frame(D_effects_dur_NL_pred) ;

################################################################################
## SCATTER PLOT:  Reese whistle duration vs NL 
################################################################################
dd1 <- ggplot() + geom_point(data = Reese, aes(NL, whistle_dur, group=as.factor(treatment), color = as.factor(treatment)), size = 5, alpha = 1) +
  
  geom_point() + xlim(95, 160) + geom_point() + ylim(0, 3.5) +
  
  geom_line(data = Reese_dur_NL_pred, aes(x=NL, y=fit), color="#5B8FA8FF",lwd=2) +
  
  geom_ribbon(data = Reese_dur_NL_pred, aes(x=NL, ymin=lower, ymax=upper), alpha= 0.3, fill="#5B8FA8FF") +
  
  labs(x="", y="") +
  
  theme(legend.position = "none", axis.text = element_text(size = 25), panel.grid = element_blank(), panel.background = element_blank()
        , axis.line = element_line(colour = "black"), text = element_text(size=28), panel.border = element_rect(fill=NA)) +
  
  scale_color_brewer(palette = "Blues", direction = 1) +
  
  scale_y_continuous(expand = c(0,0), limits = c(0, 3.5)) 

dd1 

################################################################################
## SCATTER PLOT:  Delta whistle duration vs NL
################################################################################
dd2 <- ggplot() + geom_point(data = Delta, aes(NL, whistle_dur, group = as.factor(treatment), color = as.factor(treatment)), size = 5, alpha = 1) +
  
  geom_point() + xlim(95, 160) + geom_point() + ylim(0, 3.5) +
  
  geom_line(data = Delta_dur_NL_pred, aes(x=NL, y=fit), color="#5B8FA8FF",lwd=2) + 
  
  geom_ribbon(data = Delta_dur_NL_pred, aes(x=NL, ymin=lower, ymax=upper), alpha = 0.3, fill="#5B8FA8FF") +
  
  labs(x="", y="") +
  
  theme(legend.position = "none", axis.text = element_text(size = 25), panel.grid = element_blank(), panel.background = element_blank()
        , axis.line = element_line(colour = "black"), text = element_text(size=28), panel.border = element_rect(fill=NA)) +
  
  scale_color_brewer(palette = "Blues", direction = 1) +
  
  scale_y_continuous(expand = c(0,0), limits = c(0, 3.5))

dd2


################################ END OF SCRIPT #################################

