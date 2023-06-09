library(plyr)
library(reshape2)
library(Rmisc)
library(ggplot2)

source("EvolExp_results_dataclean.R")

#### Initial data subsetting and summarizing prior to bootstrap analysis ####

# Summarize f intra-population variance by treatment
f_nocontrol <- subset(f.reps, Treatment %in% c("H", "L", "S"))
f_nocontrol$Treatment <- gsub("H", "Fast", f_nocontrol$Treatment)
f_nocontrol$Treatment <- gsub("L", "Slow", f_nocontrol$Treatment)
# "H" (meant for 'high frequency') used during data collection in raw file denotes "FAST";
# "L" (meant for 'low frequency') used during data collection in raw file denotes "SLOW"
# "S" denotes "STOCHASTIC"
f_nocontrol$rep_label <- paste(substr(f_nocontrol$Source_pop,1,1), f_nocontrol$rep, sep = "_")
f_nocontrol <- f_nocontrol[,c("Treatment", "rep_label", "f")]

f_nocontrol.summary <- summarySE(f_nocontrol, measurevar = "f", groupvars = c("Treatment", "rep_label"), na.rm = T)
f_nocontrol.summary$var <- (f_nocontrol.summary$sd)^2
f_nocontrol_popvars <- f_nocontrol.summary[, c("Treatment", "var")]
f_nocontrol_popvars$Treatment <- factor(f_nocontrol_popvars$Treatment)


# Summarize mu intra-population variance by source_pop and treatment
mu_nocontrol <- subset(mu.reps, Treatment %in% c("H", "L", "S"))
mu_nocontrol$Treatment <- gsub("H", "Fast", mu_nocontrol$Treatment)
mu_nocontrol$Treatment <- gsub("L", "Slow", mu_nocontrol$Treatment)
# "H" (meant for 'high frequency') used during data collection in raw file denotes "FAST";
# "L" (meant for 'low frequency') used during data collection in raw file denotes "SLOW"
# "S" denotes "STOCHASTIC"
mu_nocontrol$rep_label <- paste(substr(mu_nocontrol$Source_pop,1,1), mu_nocontrol$rep, sep = "_")
mu_nocontrol <- mu_nocontrol[,c("Treatment", "rep_label", "mu")]

mu_nocontrol.summary <- summarySE(mu_nocontrol, measurevar = "mu", groupvars = c("Treatment", "rep_label"), na.rm = T)
mu_nocontrol.summary$var <- (mu_nocontrol.summary$sd)^2
mu_nocontrol_popvars <- mu_nocontrol.summary[, c("Treatment", "var")]
mu_nocontrol_popvars$Treatment <- factor(mu_nocontrol_popvars$Treatment)


#### Building Bootstrap function ####

# Number of bootstrap resampling
bootstrap_num = 50000

# Bootstrap p-val function using all raw trait data (before calculating variance per replicate pop)
# This function is for comparing Slow and Fast Treatments
bootstrap_slowfast_pval <- function(raw.df, var_summ){
  # vector of all raw measurements across source pops and deterministic treatments
  raw_slowfast <- subset(raw.df[,4], raw.df$Treatment == "H" | raw.df$Treatment == "L")
  # total number of slow and fast treatment replicate populations
  # NB: Again, note differences between treatment code used during data collection and analysis (H = Fast; L = Slow)
  slow_resample_num <- sum(var_summ$Treatment == "Slow")
  fast_resample_num <- sum(var_summ$Treatment == "Fast")
  total_resample_num <- slow_resample_num + fast_resample_num
  
  # random sample of 30 from all raw *individual-level* measurements
  boot_mean_diff <-c()
  univ_random_sample <- replicate(total_resample_num, sample(raw_slowfast, 30, replace = T))
  univ_variance_pmf <- apply(univ_random_sample, 2, sd)^2
  
  for(i in 1:bootstrap_num){
    bootvars_slow <- mean(sample(univ_variance_pmf, slow_resample_num, replace = T))
    bootvars_fast <- mean(sample(univ_variance_pmf, fast_resample_num, replace = T))
    boot_mean_diff[i] <- bootvars_slow - bootvars_fast
  }
  
  fast_observed_var <- var_summ[which(var_summ$Treatment == "Fast"), "var"]
  slow_observed_var <- var_summ[which(var_summ$Treatment == "Slow"), "var"]
  observed_diff <- mean(slow_observed_var) - mean(fast_observed_var)
  
  pval = sum(boot_mean_diff >= observed_diff) / length(boot_mean_diff)
  
  output <- list(pval, observed_diff, boot_mean_diff)
  return(output)
  
}

# seed
addTaskCallback(function(...) {set.seed(321);TRUE})

# Bootstrap (using above function) of difference in MU between Deterministic Slow vs. Deterministic Fast treatments
mu_detslowfast_boot_output <- bootstrap_slowfast_pval(mu.reps, mu_nocontrol_popvars)
# p-value of that treatment-level difference
mu_detslowfast_boot_output[[1]]

# MU Permutation distribution plot; compare with observed difference
# Deterministic Slow vs. Deterministic Fast
ggplot() +
  geom_histogram(aes(mu_detslowfast_boot_output[[3]]), bins = 1000, col = "dodgerblue4") +
  scale_y_continuous(expand = c(0, 0)) +
  geom_vline(xintercept = mu_detslowfast_boot_output[[2]], linetype = "dashed", size = 1, col = "red") +
  geom_rect(aes(xmin = -Inf,
                xmax = mu_detslowfast_boot_output[[2]],
                ymin = -Inf,
                ymax = Inf), fill = "gray20", alpha = 0.3) +
  annotate("text", x = mu_detslowfast_boot_output[[2]]+0.0000025, y = 150, 
           label = paste("Observed difference \\n =", 
                         round(mu_detslowfast_boot_output[[2]], 8)),
           family = "Avenir", col = "red") +
  annotate("text", x = mu_detslowfast_boot_output[[2]]+0.0000025, y = 75,
           label = paste("p =", round(mu_detslowfast_boot_output[[1]], 3)),
            col = "red", fontface = "italic") +
  xlab("Permuted calculations of differences in intrapopulation variance") +
  ylab("Frequency") +
  ggtitle("Permutation test of difference in Var(μ) (maturation rate) between Fast and Slow treatments")+
  theme(plot.title = element_text(hjust = 0.5, family = "Avenir"),
        axis.text = element_text(family = "Avenir"),
        axis.title = element_text(family = "Avenir"))


# Bootstrap (using above function) of difference in F between Deterministic Slow vs. Deterministic Fast treatments
f_detslowfast_boot_output <- bootstrap_slowfast_pval(f.reps, f_nocontrol_popvars)
# p-value of that treatment-level difference
f_detslowfast_boot_output[[1]]

# F Permutation distribution plot; compare with observed difference
# Deterministic Slow vs. Deterministic Fast
ggplot() +
  geom_histogram(aes(f_detslowfast_boot_output[[3]]), bins = 1000, col = "dodgerblue4") +
  scale_y_continuous(expand = c(0, 0)) +
  geom_vline(xintercept = f_detslowfast_boot_output[[2]], linetype = "dashed", size = 1, col = "red") +
  geom_rect(aes(xmin = -Inf,
                xmax = f_detslowfast_boot_output[[2]],
                ymin = -Inf,
                ymax = Inf), fill = "gray20", alpha = 0.3) +
  annotate("text", x = f_detslowfast_boot_output[[2]]+0.375, y = 150, 
           label = paste("Observed difference \\n =", round(f_detslowfast_boot_output[[2]], 2)), 
           family = "Avenir", col = "red") +
  annotate("text", x = f_detslowfast_boot_output[[2]]+0.375, y = 75,
           label = paste("p =", round(f_detslowfast_boot_output[[1]], 3)),
           col = "red", fontface = "italic") +
  xlab("Permuted calculations of differences in intrapopulation variance") +
  ylab("Frequency") +
  ggtitle("Permutation test of difference in Var(f) (fecundity) between Fast and Slow treatments")+
  theme(plot.title = element_text(hjust = 0.5, family = "Avenir"),
        axis.text = element_text(family = "Avenir"),
        axis.title = element_text(family = "Avenir"))





# This function is for comparing Deterministic Slow and STOCHASTIC SLOW Treatments
bootstrap_detstoch_pval <- function(raw.df, var_summ){
  # vector of all raw mu measurements across source pops and deterministic treatments
  raw_slowstoch <- subset(raw.df[,4], raw.df$Treatment == "L" | raw.df$Treatment == "S") #L for low (i.e. slow), and S for stochastic
  # total number of slow and fast treatment replicate populations
  slow_resample_num <- sum(var_summ$Treatment == "Slow")
  stoch_resample_num <- sum(var_summ$Treatment == "S") #stochastic
  total_resample_num <- slow_resample_num + stoch_resample_num
  
  # random sample of size 10 populations from all raw mu measurements
  
  boot_mean_diff <-c()
  univ_random_sample <- replicate(total_resample_num, sample(raw_slowstoch, 30, replace = T))
  univ_variance_pmf <- apply(univ_random_sample, 2, sd)^2
  
  for(i in 1:bootstrap_num){
    bootvars_slow <- mean(sample(univ_variance_pmf, slow_resample_num, replace = T))
    bootvars_stoch <- mean(sample(univ_variance_pmf, stoch_resample_num, replace = T))
    boot_mean_diff[i] <- bootvars_stoch - bootvars_slow
  }
  
  slow_observed_var <- var_summ[which(var_summ$Treatment == "Slow"), "var"]
  stoch_observed_var <- var_summ[which(var_summ$Treatment == "S"), "var"]
  observed_diff <- mean(stoch_observed_var, na.rm = T) - mean(slow_observed_var, na.rm = T)
  
  pval = sum(boot_mean_diff >= observed_diff) / length(boot_mean_diff)
  
  output <- list(pval, observed_diff, boot_mean_diff)
  return(output)
  
  
}

addTaskCallback(function(...) {set.seed(321);TRUE})

# Bootstrap (using above function) of difference in MU between Deterministic Slow vs. Stochastic Slow treatments
mu_detstoch_boot_output <- bootstrap_detstoch_pval(mu.reps, mu_nocontrol_popvars)
# p-value of that treatment-level difference
mu_detstoch_boot_output[[1]]

# MU Permutation distribution plot; compare with observed difference
# Deterministic Slow vs. Stochastic Slow
ggplot() +
  geom_histogram(aes(mu_detstoch_boot_output[[3]]), bins = 1000, col = "dodgerblue4") +
  scale_y_continuous(expand = c(0, 0)) +
  geom_vline(xintercept = mu_detstoch_boot_output[[2]], linetype = "dashed", size = 1, col = "red") +
  geom_rect(aes(xmin = -Inf,
                xmax = mu_detstoch_boot_output[[2]],
                ymin = -Inf,
                ymax = Inf), fill = "gray20", alpha = 0.3) +
  annotate("text", x = mu_detstoch_boot_output[[2]]+0.00001, y = 150, 
           label = paste("Observed difference \\n =", 
                         round(mu_detstoch_boot_output[[2]], 8)),
           family = "Avenir", col = "red") +
  annotate("text", x = mu_detstoch_boot_output[[2]]+0.00001, y = 75,
           label = paste("p =", round(mu_detstoch_boot_output[[1]], 3)),
           col = "red", fontface = "italic") +
  xlab("Permuted calculations of differences in intrapopulation variance") +
  ylab("Frequency") +
  ggtitle("Permutation test of difference in Var(μ) (maturation rate) between \\nDeterministic and Stochastic Slow treatments")+
  theme(plot.title = element_text(hjust = 0.5, family = "Avenir"),
        axis.text = element_text(family = "Avenir"),
        axis.title = element_text(family = "Avenir"))

# Bootstrap (using above function) of difference in F between Deterministic Slow vs. Stochastic Slow treatments
f_detstoch_boot_output <- bootstrap_detstoch_pval(f.reps, f_nocontrol_popvars)
# p-value of that treatment-level difference
f_detstoch_boot_output[[1]]

# F Permutation distribution plot; compare with observed difference
# Deterministic Slow vs. Stochastic Slow
ggplot() +
  geom_histogram(aes(f_detstoch_boot_output[[3]]), bins = 1000, col = "dodgerblue4") +
  scale_y_continuous(expand = c(0, 0)) +
  geom_vline(xintercept = f_detstoch_boot_output[[2]], linetype = "dashed", size = 1, col = "red") +
  geom_rect(aes(xmin = -Inf,
                xmax = f_detstoch_boot_output[[2]],
                ymin = -Inf,
                ymax = Inf), fill = "gray20", alpha = 0.3) +
  annotate("text", x = f_detstoch_boot_output[[2]]+2, y = 150, 
           label = paste("Observed difference \\n =", 
                         round(f_detstoch_boot_output[[2]], 8)),
           family = "Avenir", col = "red") +
  annotate("text", x = f_detstoch_boot_output[[2]]+2, y = 75,
           label = paste("p =", round(f_detstoch_boot_output[[1]], 3)),
           col = "red", fontface = "italic") +
  xlab("Permuted calculations of differences in intrapopulation variance") +
  ylab("Frequency") +
  ggtitle("Permutation test of difference in Var(f) (fecundity) between \\nDeterministic and Stochastic Slow treatments")+
  theme(plot.title = element_text(hjust = 0.5, family = "Avenir"),
        axis.text = element_text(family = "Avenir"),
        axis.title = element_text(family = "Avenir"))
