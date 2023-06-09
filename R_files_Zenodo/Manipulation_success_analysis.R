####----------------------------------------------------------------------------------------------
# Density manipulation 
####----------------------------------------------------------------------------------------------



# Packages 
####-------------------------


library(data.table)
library(ggplot2)
library(patchwork)
library(dplyr)
library(effects)
library(glmmTMB)
library(emmeans)

# Load data
load("data_density.RData")
# period: pre=prior to density manipulation, during=during density manipulation, exp_day: experimental day, 
# N_rec: Number of recordings, N_ids: Number of individuals visited, total_rec: total number of recordings to both loggers, prop: proportion of total recordings


#   1)  Number of recordings 
####-------------------------

m1 <- glmmTMB(cbind(N_rec, total_rec-N_rec) ~ period*treatment + exp_day + (1|Site) + ar1(factor(exp_day)+0|logger), 
              data=data_density, family = 'binomial') 

# Post-hoc comparisons between combinations of period and manipulation 
emmeans(m1, list(pairwise ~ period*treatment), adjust = "tukey", type = "response")    

# Extract effects
e <- allEffects(m1)$`period:treatment` %>% data.frame %>% data.table

# Plot
col.treatment <- c("high" = "orange", "low" = "darkslategrey")
(plot1 <- ggplot() +
    geom_jitter(data=data_density, aes(y = prop, x = period, colour=treatment, group=treatment), alpha=.15,width = 0.1) +
    geom_point(data = e, aes(y = fit, x = period, colour=treatment, group=treatment), size=2,position=position_dodge(width = 0.2)) + 
    geom_line(data = e, aes(y = fit, x = period, colour=treatment, group=treatment), position=position_dodge(width = 0.2)) +
    geom_errorbar(data = e, aes(x = period , ymin=(lower), ymax=(upper), width=.1, colour=treatment), position=position_dodge(width = 0.2)) + 
    theme_classic() + ylab("Proportion of recordings") + scale_color_manual(values=col.treatment) + 
    theme(axis.text.x = element_text( color="black", size=15), 
          axis.text.y = element_text( color="black", size=15), 
          axis.title.x = element_blank(), axis.title.y = element_text(size=16), legend.position = "none"))



#   2)  Number of individuals 
####-------------------------

m2 <- glmmTMB(N_ids ~ period*treatment + exp_day + (1|Site) + ar1(factor(exp_day)+0|logger), 
              data=data_density, family = 'poisson') 

# Post-hoc comparisons between combinations of period and manipulation 
emmeans(m2, list(pairwise ~ period*treatment), adjust = "tukey", type = "response")  

# Extract effects
e <- allEffects(m2)$`period:treatment` %>% data.frame %>% data.table

# Plot
col.treatment <- c("high" = "orange", "low" = "darkslategrey")
(plot2 <- ggplot() +
    geom_jitter(data=data_density, aes(y = N_ids, x = period, colour=treatment, group=treatment), alpha=.15,width = 0.1) +
    geom_point(data = e, aes(y = fit, x = period, colour=treatment, group=treatment), size=2,position=position_dodge(width = 0.2)) + 
    geom_line(data = e, aes(y = fit, x = period, colour=treatment, group=treatment), position=position_dodge(width = 0.2)) +
    geom_errorbar(data = e, aes(x = period , ymin=(lower), ymax=(upper), width=.1, colour=treatment), position=position_dodge(width = 0.2)) + 
    theme_classic() + ylab("Number of individuals") + scale_color_manual(values=col.treatment) + 
    theme(axis.text.x = element_text( color="black", size=15), 
          axis.text.y = element_text( color="black", size=15), 
          axis.title.x = element_blank(), axis.title.y = element_text(size=16),
          strip.text = element_blank(), legend.position = "none"))

plot1+plot2






