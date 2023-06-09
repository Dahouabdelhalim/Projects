source("EvolExp_results_dataclean.R")

library(Rmisc)

#### FECUNDITY (f) ####
# Subsetting and summarizing data for analysis 
# Ignoring "Control" replicates for all following analyses, denoted by "C" in Treatment names
f_nocontrol <- subset(f.reps, Treatment %in% c("H", "L"))
# "H" (meant for 'high frequency') used during data collection in raw file denotes "FAST";
# "L" (meant for 'low frequency') used during data collection in raw file denotes "SLOW"

# Summary of data without "Control"'s
# NB: "Source_pop" denotes sourced stock pop container from pre-exp rearing, not used in analysis (just used for extra measure of replicability)
f_nocontrol.summary <- summarySE(f_nocontrol, measurevar = "f", groupvars = c("Source_pop", "Treatment", "rep"), na.rm = T)
f_nocontrol.summary$var <- (f_nocontrol.summary$sd)^2
f_nocontrol_popvars <- f_nocontrol.summary[, c("Source_pop", "Treatment", "rep", "var")]
f_nocontrol_popvars_summ <- summarySE(f_nocontrol_popvars, measurevar = "var", groupvars = c("Treatment"))

#### MATURATION (mu) ####
# Subsetting and summarizing data for analysis 
# Ignoring "Control" replicates for all following analyses, denoted by "C" in Treatment names
mu_nocontrol <- subset(mu.reps, Treatment %in% c("H", "L"))
# "H" (meant for 'high frequency') used during data collection in raw file denotes "FAST";
# "L" (meant for 'low frequency') used during data collection in raw file denotes "SLOW"

# Summary of data without "Control"'s
# NB: "Source_pop" denotes sourced stock pop container from pre-exp rearing, not used in analysis (just used for extra measure of replicability)
mu_nocontrol.summary <- summarySE(mu_nocontrol, measurevar = "mu", groupvars = c("Source_pop", "Treatment", "rep"), na.rm = T)
mu_nocontrol.summary$var <- (mu_nocontrol.summary$sd)^2
mu_nocontrol_popvars <- mu_nocontrol.summary[, c("Source_pop", "Treatment", "rep", "var")]
mu_nocontrol_popvars_summ <- summarySE(mu_nocontrol_popvars, measurevar = "var", groupvars = c("Treatment"))
  