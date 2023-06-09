source("EvolExp_results_dataclean.R")

library(ggplot2)

#### FECUNDITY (f) ####

# Cleaning f intra-population variance by treatment and replicate pop
# Ignoring "Control" replicates for all following analyses, denoted by "C" in Treatment names
f_nocontrol <- subset(f.reps, Treatment %in% c("H", "L", "S"))
# "H" (meant for 'high frequency') used during data collection in raw file denotes "FAST";
# "L" (meant for 'low frequency') used during data collection in raw file denotes "SLOW"
# "S" denotes "STOCHASTIC"
# NB: "Source_pop" denotes sourced stock pop container from pre-exp rearing, not used in analysis (just used for extra measure of replicability)

f_nocontrol$Treatment <- gsub("H", "Fast", f_nocontrol$Treatment)
f_nocontrol$Treatment <- gsub("L", "Slow", f_nocontrol$Treatment)
f_nocontrol$rep_label <- paste(substr(f_nocontrol$Source_pop,1,1), f_nocontrol$Treatment, f_nocontrol$rep, sep = "_")
f_nocontrol <- f_nocontrol[,c("Treatment", "rep_label", "f")]

# Plot individual f & population-level means for the three treatments
f.mean_expand.p <- ggplot(transform(f_nocontrol, Treatment=factor(Treatment, levels=c("Fast", "Slow", "S"))), 
                           aes(rep_label, f)) + 
  geom_jitter(aes(col = Treatment), 
              position = position_jitterdodge(jitter.width = 1, dodge.width = 0.5),
              size = 1.2, alpha = 0.4) +
  stat_summary(fun.data="mean_sdl",  fun.args = list(mult=1), 
               geom = "pointrange",  size = 0.4,
               position = position_dodge(0), na.rm = T) +
  scale_color_manual(values = c("#EB95FF", "#5AFFC3", "#36F2F9")) +
  facet_wrap(~Treatment, scales = "free_x",
             labeller = labeller(Treatment = c("Fast" = "Fast", "Slow" = "Slow", "S" = "Slow Stochastic")),
             strip.position = "bottom") +
  xlab("Treatment x replicate populations") +
  ylab(expression(paste(italic("f"), " of individuals (Mean ± SD of population)"))) +
  ggtitle(expression(paste("Mean (", italic("f"), ") (fecundity) in experimental populations"))) +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, family = "Avenir", size = 13),
        panel.grid = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x =  element_blank(),
        strip.text = element_text(family = "Avenir", size = 11)) 

# Summarizing f intra-population variance by treatment and replicate pop
f_nocontrol.summary <- summarySE(f_nocontrol, measurevar = "f", groupvars = c("Treatment", "rep_label"), na.rm = T)
f_popmeans_detslow.vec <- f_nocontrol.summary$f[f_nocontrol.summary$Treatment == "Slow"] # Deterministic Slow
f_popmeans_detfast.vec <- f_nocontrol.summary$f[f_nocontrol.summary$Treatment == "Fast"] # Deterministic Fast
f_popmeans_detstoch.vec <- f_nocontrol.summary$f[f_nocontrol.summary$Treatment == "S"] # Stochastic Slow

# Mann-Whitney tests
wilcox.test(f_popmeans_detslow.vec, f_popmeans_detfast.vec, alternative = "two.sided")
wilcox.test(f_popmeans_detslow.vec, f_popmeans_detstoch.vec, alternative = "two.sided")


#### MATURATION (mu) ####

# Summarize mu intra-population variance by treatment and replicate pop
# Ignoring "Control" replicates for all following analyses, denoted by "C" in Treatment names
mu_nocontrol <- subset(mu.reps, Treatment %in% c("H", "L", "S"))
# "H" (meant for 'high frequency') used during data collection in raw file denotes "FAST";
# "L" (meant for 'low frequency') used during data collection in raw file denotes "SLOW"
# "S" denotes "STOCHASTIC"
# NB: "Source_pop" denotes sourced stock pop container from pre-exp rearing, not used in analysis (just used for extra measure of replicability)

mu_nocontrol$Treatment <- gsub("H", "Fast", mu_nocontrol$Treatment)
mu_nocontrol$Treatment <- gsub("L", "Slow", mu_nocontrol$Treatment)
mu_nocontrol$rep_label <- paste(substr(mu_nocontrol$Source_pop,1,1), mu_nocontrol$Treatment, mu_nocontrol$rep, sep = "_")
mu_nocontrol <- mu_nocontrol[,c("Treatment", "rep_label", "mu")]

# Plot individual mu & population-level means for the three treatments
mu.mean_expand.p <- ggplot(transform(mu_nocontrol, Treatment=factor(Treatment, levels=c("Fast", "Slow", "S"))), 
                                     aes(rep_label, mu)) + 
  geom_jitter(aes(col = Treatment), 
              position = position_jitterdodge(jitter.width = 1, dodge.width = 0.5),
              size = 1.2) +
  stat_summary(fun.data="mean_sdl",  fun.args = list(mult=1), 
               geom = "pointrange",  size = 0.4,
               position = position_dodge(0), na.rm = T) +
  scale_color_manual(values = c("#EB95FF", "#5AFFC3", "#36F2F9")) +
  facet_wrap(~Treatment, scales = "free_x",
             labeller = labeller(Treatment = c("Fast" = "Fast", "Slow" = "Slow", "S" = "Slow Stochastic")),
             strip.position = "bottom") +
  xlab("Treatment x replicate populations") +
  ylab(expression(paste(italic("\\u03BC"), " of individuals (Mean ± SD of population)"))) +
  ggtitle(expression(paste("Mean (", italic("\\u03BC"), ") (maturation rate) in experimental populations"))) +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, family = "Avenir", size = 13),
        panel.grid = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x =  element_blank(),
        strip.text = element_text(family = "Avenir", size = 11)) 

# Summarizing mu intra-population variance by treatment and replicate pop
mu_nocontrol.summary <- summarySE(mu_nocontrol, measurevar = "mu", groupvars = c("Treatment", "rep_label"), na.rm = T)
mu_popmeans_detslow.vec <- mu_nocontrol.summary$mu[mu_nocontrol.summary$Treatment == "Slow"]
mu_popmeans_detfast.vec <- mu_nocontrol.summary$mu[mu_nocontrol.summary$Treatment == "Fast"]
mu_popmeans_detstoch.vec <- mu_nocontrol.summary$mu[mu_nocontrol.summary$Treatment == "S"]

# Mann-Whitney tests
wilcox.test(mu_popmeans_detslow.vec, mu_popmeans_detfast.vec, alternative = "two.sided")
wilcox.test(mu_popmeans_detslow.vec, mu_popmeans_detstoch.vec, alternative = "two.sided")


