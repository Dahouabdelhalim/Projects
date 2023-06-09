####
### addessi_etal_rcode.R
### Created by Jeffrey R. Stevens on 28 Sept 2013 (jeffrey.r.stevens@gmail.com)
### Finalized on: 2020-07-16
### Summary: This script calculates descriptive and inferential statistics, 
###     and generates figures for the analysis of capuchin intertemporal choice data.
### Instructions: Place this file and the data files (addessi_etal_data.csv)
### 	in the same directory.  Create a folder called "figures". Set the R
### 	working directory to this directory.  At the R command prompt, type 
### 	> source("addessi_etal_rcode.R")
### 	This will run the script, adding all of the calculated variables to the
### 	workspace and saving PNG versions of the figures in the figures directory.
### Uses: This script can be reproduced and modified for personal and scientific use.
### Data files: Description of the data columns:
###  addessi_etal_data.csv--intertemporal choice data
###   subject - subject name
###   session - session number
###   date - date
###   condition - experimental condition (High cost, Low cost same, Low cost different--see manuscript for explanation)
###   trial_num - trial number
###   ll_side - side on which the larger, later option was placed
###   choice - choice between larger, later (1) and smaller, sooner (0) option
###   last_choice - choice in previous trial (0 = smaller, sooner; 1 = larger, later)
####

####
# Load libraries and define functions -------------------------------------
####

library(BayesFactor)  # needed to calculate Bayes factors
library(bayestestR)  # needed for estimating Bayes factors for models using BICs
library(brms)  # needed to fit Bayesian models
library(emmeans)  # needed to calculate estimated marginal means
library(foreach)  # needed for iterations
library(here)		# needed for here
library(lme4)	# needed for GLMMs
library(papaja)    # needed for within-subjects confidence intervals
library(tidyverse)	# needed for tidyverse

  
###
## Create themes for plots
###		
theme_plots <- function () { 
  theme_bw(base_size=20) %+replace% 
    theme(
      panel.grid = element_blank(),
      legend.position = "none",
      
    )
}

theme_legend <- function () { 
  theme_bw(base_size=30) %+replace% 
    theme(
      panel.grid = element_blank(),
      legend.title = element_blank()  # remove legend title
    )
}


####
# Input and prepare data -------------------------------------
####

## Input data
data <- read_csv(here("data/addessi_etal_data.csv")) %>% 	# input data file
  mutate(condition = factor(condition, levels = c("High cost", "Low cost different", "Low cost same"))) %>% 	# assign levels of condition
  arrange(subject, session)  # sort rows

## Summarize data over sessions
condition_subject <- data %>% 
  group_by(subject, condition) %>%  # for each subject and session
  summarize(num_sessions = max(session), # calculate the number of sessions
            min_session = num_sessions - 5, # find the minmum number of sessions for stable data
            mean_choice = mean(choice, na.rm = TRUE))  # calculate the mean choice proportions

## Extract only the stable data
# Initiate data frame
stable_data <- data[1, ]  # get first row of data
stable_data$end_session <- NA  # add end_session column
stable_data <- stable_data[-1, ]  # remove data

# Extract stable data for each subject and condition
foreach(current_subject = unique(data$subject)) %do% {  # for each subject
  foreach(current_condition = unique(data$condition)) %do% {  # for each condition
    current_data <- filter(data, subject == current_subject & condition == current_condition & session > condition_subject$min_session[which(condition_subject$subject == current_subject & condition_subject$condition == current_condition)])  # filter the last five sessions for current subject and condition
    current_data$end_session <- current_data$session - max(current_data$session) + 5  # renumber sessions to be 1-5
    stable_data <- bind_rows(stable_data, current_data)  # append data
  }
}

## Summarize data over sessions for stable data
condition_subject_stable <- stable_data %>% 
  group_by(subject, condition) %>%  # for each subject and condition
  summarize(num_sessions = max(session), # calculate the total number of sessions
            mean_choice = mean(choice, na.rm = TRUE))  # calculate the mean choice proportions

subject_means <- condition_subject_stable %>% 
  group_by(subject) %>%  # for each subject
  summarize(mean_choice = mean(mean_choice))  # calculate the mean choice proportions
  
# Calculate number of standard deviations from the overall mean for Robot
overall_mean <- mean(subject_means$mean_choice)  # calculate overall mean
overall_sd <- sd(subject_means$mean_choice)  # calculate overall SD
robot_sds <- (filter(condition_subject_stable, subject == "Robot" & condition == "Low cost same")$mean_choice - overall_mean) / overall_sd  # calculate departure from overall mean for Robot's Low cost same response


####
# Descriptive statistics --------------------------------------------------
####

## Calculate means and within-subjects 95% CIs for all data
condition_means2 <- condition_subject_stable %>% 
  group_by(condition) %>%  # for each condition
  summarize(choice_mean = mean(mean_choice)  # calculate mean choice proportions
  )
condition_means2$choice_ci <- wsci(condition_subject_stable, id = "subject", factors = "condition", dv = "mean_choice")$mean_choice  # calculate within-subjects 95% confidence intervals for choice proportions
condition_means2$choice_lowerci <- condition_means2$choice_mean - condition_means2$choice_ci  # calculate lower CI
condition_means2$choice_upperci <- condition_means2$choice_mean + condition_means2$choice_ci  # calculate upper CI

## Calculate means and within-subjects 95% CIs for data without Robot
condition_subject_stable_trimmed <- condition_subject_stable %>% 
  filter(subject != "Robot")  # remove Robot
condition_means2_trimmed <- condition_subject_stable_trimmed %>% 
  group_by(condition) %>%   # for each condition
  summarize(choice_mean = mean(mean_choice)  # calculate mean choice proportions
  )
condition_means2_trimmed$choice_ci <- wsci(condition_subject_stable_trimmed, id = "subject", factors = "condition", dv = "mean_choice")$mean_choice   # calculate within-subjects 95% confidence intervals for choice proportions
condition_means2_trimmed$choice_lowerci <- condition_means2_trimmed$choice_mean - condition_means2_trimmed$choice_ci  # calculate lower CI
condition_means2_trimmed$choice_upperci <- condition_means2_trimmed$choice_mean + condition_means2_trimmed$choice_ci  # calculate upper CI


####
# Inferential statistics --------------------------------------------------
####

# _ Choice preferences differ from chance --------------------------------------------------------

## Create data frame for each condition
high_cost <- filter(condition_subject_stable, condition == "High cost")
low_cost_different <- filter(condition_subject_stable, condition == "Low cost different")
low_cost_same <- filter(condition_subject_stable, condition == "Low cost same")

## Check normality
shapiro.test(high_cost$mean_choice)
shapiro.test(low_cost_different$mean_choice)
shapiro.test(low_cost_same$mean_choice)

## Conduct one-sample Wilcoxon tests
wilcox.test(high_cost$mean_choice, mu = 0.5)
wilcox.test(low_cost_different$mean_choice, mu = 0.5)
wilcox.test(low_cost_same$mean_choice, mu = 0.5)

## Conduct Bayes factor t-tests
ttestBF(high_cost$mean_choice, mu = 0.5)  
ttestBF(low_cost_different$mean_choice, mu = 0.5)
ttestBF(low_cost_same$mean_choice, mu = 0.5)

# _ Models all data --------------------------------------------------------

# __ Frequentist analysis --------------------------------------------------------
## Random effects 
# Models
glmer_intercept <- glm(choice ~ 1, data = stable_data, family = binomial(link="logit"))  # model with intercept only (no random effects)
glmer_rand <- glmer(choice ~ (1 | subject), data = stable_data, family = binomial(link="logit"))  # model with random effect of subject
anova(glmer_rand, glmer_intercept)  # likelihood ratio test of random against intercept only
  # Results: glmer_rand has lower AIC than glmer_intercept

# Estimated Bayes factor
bf_models(glmer_rand, denominator = glmer_intercept)  # BF for random against intercept only
  # Results: BF = 1.080e+11

## Fixed effects 
# Models
glmer_condition <- glmer(choice ~ condition +  (1 | subject), data = stable_data, family = binomial(link="logit"))  # model with condition
glmer_trial_num <- glmer(choice ~ trial_num + (1 | subject), data = stable_data, family = binomial(link="logit"))  # model with trial number
glmer_full <- glmer(choice ~ condition + trial_num + (1 | subject), data = stable_data, family = binomial(link="logit"))  # model with condition and trial number

# Backwards selection
summary(glmer_full)  # test full model
  # Results: effect of condition
summary(glmer_condition)  # drop trial_num
  # Results: effect of condition

# Forward selection
anova(glmer_rand, glmer_condition)  # likelihood ratio test condition against random
  # Results: glmer_condition has lower AIC than glmer_rand
anova(glmer_rand, glmer_trial_num)  # likelihood ratio test trial_num against random
  # Results: glmer_trial_num does not have lower AIC than glmer_rand
anova(glmer_condition, glmer_full)  # likelihood ratio test condition and trial_num against condition
  # Results: glmer_full does not have lower AIC than glmer_condition

# Estimated Bayes factors
bf_models(glmer_condition, denominator = glmer_full)  # BF for condition against condition and trial_num
  # Results: BF = 28.90
bf_models(glmer_condition, denominator = glmer_rand)  # BF for condition against random
  # Results: BF = 1.348e+10

## Condition contrasts
summary(glmer_condition)  # summarize best model
condition_emmeans <- emmeans(glmer_condition, ~ condition)  # compute estimated marginal means
pairs(condition_emmeans)  # compute pairwise contrasts
# Results:
# contrast                           estimate    SE  df z.ratio p.value
# High cost - Low cost different       -2.214 0.334 Inf -6.628  <.0001 
# High cost - Low cost same            -1.606 0.342 Inf -4.693  <.0001 
# Low cost different - Low cost same    0.608 0.222 Inf  2.741  0.0169 

## Exploratory test of effect of previous choices
# Models
glmer_last_choice_intercept <- glmer(choice ~ 1 + (1 | subject), data = filter(stable_data, !is.na(last_choice)), family = binomial(link="logit"))  # model with intercept only
glmer_last_choice <- glmer(choice ~ last_choice + (1 | subject), data = filter(stable_data, !is.na(last_choice)), family = binomial(link="logit"))  # model with last_choice
glmer_last_choice_condition <- glmer(choice ~ condition + (1 | subject), data = filter(stable_data, !is.na(last_choice)), family = binomial(link="logit"))  # model with condition
glmer_last_choice_full <- glmer(choice ~ condition + last_choice + (1 | subject), data = filter(stable_data, !is.na(last_choice)), family = binomial(link="logit"))  # model with condition and last_choice

# Backwards selection
summary(glmer_last_choice_full)
  # Results: Effects of both condition and last_choice

# Forwards selection
anova(glmer_last_choice_intercept, glmer_last_choice_condition, glmer_last_choice)  # likelihood ratio test condition and last_choice against intercept
  # Results: glmer_last_choice_condition has lowest AIC
anova(glmer_last_choice_condition, glmer_last_choice_full)  # likelihood ratio test condition and last_choice against condition
  # Results: glmer_last_choice_full has lowest AIC

# Estimated Bayes factor
bf_models(glmer_last_choice_full, denominator = glmer_last_choice_condition)  # BF for condition and last_choice against condition
  # Results: BF = 7.134

# __ Bayesian analysis --------------------------------------------------------
## Prepare for analysis
## Note that the Bayesian analyses require a lot of computation and can take a lot of time to compute
ncores = parallel::detectCores()  # find numbers of cores to parallelize analysis
options(contrasts = c('contr.bayes', 'contr.poly'), mc.cores = parallel::detectCores())  # set options
fixed_priors <- set_prior("normal(0, 1)", class = "b")  # set fixed priors
intercept_priors <- set_prior("student_t(3, 0, 2.5)", class = "Intercept")  # set intercept priors

## Random effects
# Models
bglmer_intercept <- brm(choice ~ 1, data = stable_data, family = bernoulli(link = "logit"), prior = intercept_priors, iter = 20000, cores = ncores, save_all_pars = TRUE)  # model with intercept only
bglmer_rand <- brm(choice ~ 1 + (1|subject), data = stable_data, family = bernoulli(link = "logit"), prior = intercept_priors, iter = 20000, cores = ncores, save_all_pars = TRUE)  # model with random effect of subject

# Bayes factor
bglmer_intercept_rand_bf <- bf_models(bglmer_rand, denominator = bglmer_intercept)  # BF for random against intercept
  # BF = 1.95e+12

## Fixed effects
# Models
bglmer_condition <- brm(choice ~ condition +  (1|subject), data = stable_data, family = bernoulli(link = "logit"), prior = fixed_priors, iter = 20000, cores = ncores, save_all_pars = TRUE)  # model with condition
bglmer_trial_num <- brm(choice ~ trial_num + (1|subject), data = stable_data, family = bernoulli(link = "logit"), prior = fixed_priors, iter = 20000, cores = ncores, save_all_pars = TRUE)  # model with condition and trial number
bglmer_full <- brm(choice ~ condition + trial_num + (1|subject), data = stable_data, family = bernoulli(link = "logit"), prior = fixed_priors, iter = 20000, cores = ncores, save_all_pars = TRUE)  # model with condition and trial number

# Bayes factors
# bglmer_rand_condition_bf <- bf_models(bglmer_condition, denominator = bglmer_rand)  # BF for condition against random
  # Fails due to error (Failed to initialize module pointer: no such symbol _rcpp_module_boot_stan_fit4model431b2c110d60_dd4649e402fd4afe76a096455fc34056_mod in package /tmp/RtmpUg1j1b/file431b29cf7668.so)
bglmer_full_condition_bf <- bf_models(bglmer_condition, denominator = bglmer_full)  # BF for condition against condition and trial_num
  # Results: BF = 15.926

## Condition contrasts
condition_emmeans_b <- emmeans(bglmer_condition, ~ condition)  # compute estimated marginal means
condition_pairs_b <- pairs(condition_emmeans_b)  # compute pairwise contrasts
contrasts_bf <- bayesfactor_parameters(condition_pairs_b, prior = bglmer_condition)  # compute pairwise contrasts
# Results:
# High cost - Low cost different     | 4.901e+06
# High cost - Low cost same          |  3613.912
# Low cost different - Low cost same |     6.757

## Exploratory  test of effect of previous choices
# Models
bglmer_last_choice_condition <- brm(choice ~ condition + (1 | subject), data = filter(stable_data, !is.na(last_choice)), family = bernoulli(link = "logit"), prior = fixed_priors, iter = 20000, cores = ncores, save_all_pars = TRUE)  # model with condition
bglmer_last_choice_full <- brm(choice ~ condition + last_choice + (1 | subject), data = filter(stable_data, !is.na(last_choice)), family = bernoulli(link = "logit"), prior = fixed_priors, iter = 20000, cores = ncores, save_all_pars = TRUE)  # model with condition and last_choice

# Bayes factor
bglmer_last_choice_condition_bf <- bf_models(bglmer_last_choice_full, denominator = bglmer_last_choice_condition)  # BF for condition and last_choice against condition
  # Results: BF = 44.154


# _ Models without Robot's data --------------------------------------------------------
## Remove Robot's data
stable_data_trimmed <- filter(stable_data, subject != "Robot")

# __ Frequentist analysis --------------------------------------------------------
## Random effects 
# Models
glmer_intercept_trimmed <- glm(choice ~ 1, data = stable_data_trimmed, family = binomial(link="logit"))  # model with intercept only (no random effects)
glmer_rand_trimmed <- glmer(choice ~ (1 | subject), data = stable_data_trimmed, family = binomial(link="logit"))  # model with random effect of subject
anova(glmer_rand_trimmed, glmer_intercept_trimmed)  # likelihood ratio test of random against intercept only
  # Results: glmer_rand_trimmed has lower AIC

# Estimated Bayes factor
bf_models(glmer_rand_trimmed, denominator = glmer_intercept_trimmed)  # BF for random against intercept only
# Results: BF = 26136.862

## Fixed effects 
# Models
glmer_condition_trimmed <- glmer(choice ~ condition +  (1 | subject), data = stable_data_trimmed, family = binomial(link="logit"))  # model with condition
glmer_trial_num_trimmed <- glmer(choice ~ trial_num +  (1 | subject), data = stable_data_trimmed, family = binomial(link="logit"))  # model with trial_num
glmer_full_trimmed <- glmer(choice ~ condition + trial_num + (1 | subject), data = stable_data_trimmed, family = binomial(link="logit"))  # model with condition and trial number

# Backwards selection
summary(glmer_full_trimmed)  # test condition and trial_num model
  # Results: effect of condition
summary(glmer_condition_trimmed)  # drop trial_num
  # Results: effect of condition

# Forward selection
anova(glmer_rand_trimmed, glmer_condition_trimmed)  # likelihood ratio test condition against random
  # Results: glmer_condition_trimmed has lower AIC than glmer_rand_trimmed
anova(glmer_rand_trimmed, glmer_trial_num_trimmed)  # likelihood ratio test trial_num against random
  # Results: glmer_trial_num_trimmed does not have lower AIC than glmer_rand_trimmed
anova(glmer_condition_trimmed, glmer_full_trimmed)  # likelihood ratio test condition and trial_num against condition
  # Results: glmer_full_trimmed does not have lower AIC than glmer_condition_trimmed

# Estimated Bayes factors
bf_models(glmer_condition_trimmed, denominator = glmer_full_trimmed)  # BF for condition against condition and trial_num
  # Results: BF = 23.20
bf_models(glmer_condition_trimmed, denominator = glmer_rand_trimmed)  # BF for condition against random
  # Results: BF = 1.401e+17

## Condition contrasts
summary(glmer_condition_trimmed)  # summarize best model
condition_emmeans2_trimmed <- emmeans(glmer_condition_trimmed, ~ condition)  # compute estimated marginal means
pairs(condition_emmeans2_trimmed)  # compute pairwise contrasts
# Results:
# contrast                           estimate   SE  df z.ratio p.value
# High cost - Low cost different        -2.95 0.44 Inf -6.697  <.0001 
# High cost - Low cost same             -1.18 0.48 Inf -2.450  0.0380 
# Low cost different - Low cost same     1.77 0.29 Inf  6.106  <.0001 

## Exploratory test of effect of previous choices
# Models
glmer_last_choice_intercept_trimmed <- glmer(choice ~ 1 + (1 | subject), data = filter(stable_data_trimmed, !is.na(last_choice)), family = binomial(link="logit"))  # model with intercept only
glmer_last_choice_trimmed <- glmer(choice ~ last_choice + (1 | subject), data = filter(stable_data_trimmed, !is.na(last_choice)), family = binomial(link="logit"))  # model with last_choice
glmer_last_choice_condition_trimmed <- glmer(choice ~ condition + (1 | subject), data = filter(stable_data_trimmed, !is.na(last_choice)), family = binomial(link="logit"))  # model with condition
glmer_last_choice_full_trimmed <- glmer(choice ~ condition + last_choice + (1 | subject), data = filter(stable_data_trimmed, !is.na(last_choice)), family = binomial(link="logit"))  # model with condition and last_choice

# Backwards selection
summary(glmer_last_choice_full_trimmed)
# Results: Effects of condition but not last_choice

# Forwards selection
anova(glmer_last_choice_intercept_trimmed, glmer_last_choice_condition_trimmed, glmer_last_choice_trimmed)  # likelihood ratio test condition and last_choice against intercept
# Results: glmer_last_choice_condition_trimmed has lower AIC than glmer_last_choice_intercept_trimmed
anova(glmer_last_choice_condition_trimmed, glmer_last_choice_full_trimmed)  # likelihood ratio test condition and last_choice against condition
# Results: glmer_last_choice_full does not have lower AIC than glmer_last_choice_condition_trimmed

# Estimated Bayes factor
bf_models(glmer_last_choice_full_trimmed, denominator = glmer_last_choice_condition_trimmed)  # BF for condition and trial_num against condition
# Results: BF = 0.143

# __ Bayesian analysis --------------------------------------------------------
## Random effects
# Models
bglmer_intercept_trimmed <- brm(choice ~ 1, data = stable_data_trimmed, family = bernoulli(link = "logit"), prior = intercept_priors, iter = 20000, cores = ncores, save_all_pars = TRUE)  # model with random effect of subject
bglmer_rand_trimmed <- brm(choice ~ 1 + (1|subject), data = stable_data_trimmed, family = bernoulli(link = "logit"), prior = intercept_priors, iter = 20000, cores = ncores, save_all_pars = TRUE)  # model with random effect of subject

# Bayes factor
bglmer_intercept_rand_trimmed_bf <- bf_models(bglmer_rand_trimmed, denominator = bglmer_intercept_trimmed)  # BF for random against intercept
  # Results: BF = 4.200e+05

## Fixed effects
# Models
bglmer_condition_trimmed <- brm(choice ~ condition +  (1|subject), data = stable_data_trimmed, family = bernoulli(link = "logit"), prior = fixed_priors, iter = 20000, cores = ncores, save_all_pars = TRUE)  # model with condition
bglmer_trial_num_trimmed <- brm(choice ~ trial_num + (1|subject), data = stable_data_trimmed, family = bernoulli(link = "logit"), prior = fixed_priors, iter = 20000, cores = ncores, save_all_pars = TRUE)  # model with condition and trial number
bglmer_full_trimmed <- brm(choice ~ condition + trial_num + (1|subject), data = stable_data_trimmed, family = bernoulli(link = "logit"), prior = fixed_priors, iter = 20000, cores = ncores, save_all_pars = TRUE)  # model with condition and trial number

# Bayes factors
# bglmer_rand_condition_trial_num_trimmed_bf <- bf_models(bglmer_condition_trimmed, bglmer_trial_num_trimmed, denominator = bglmer_rand_trimmed)  # BF for condition and trial_num against random
  # Results: BF(condition) = 8.379e+17, BF(trial_num) = 0.078 (sometimes fails to initialize module pointer: no such symbol _rcpp_module_boot_stan_fit4model431b2c110d60_dd4649e402fd4afe76a096455fc34056_mod in package /tmp/RtmpUg1j1b/file431b29cf7668.so)
bglmer_full_condition_trimmed_bf <- bf_models(bglmer_condition_trimmed, denominator = bglmer_full_trimmed)  # BF for condition against condition and trial_num
  # Results: BF = 11.92

## Condition contrasts
condition_emmeans_trimmed <- emmeans(bglmer_condition_trimmed, ~ condition)  # compute estimated marginal means
condition_pairs_trimmed <- pairs(condition_emmeans_trimmed)  # compute pairwise contrasts
contrasts_trimmed_bf <- bayesfactor_parameters(condition_pairs_trimmed, prior = bglmer_condition_trimmed)  # compute pairwise contrasts
# High cost - Low cost different     | 2.559e+07
# High cost - Low cost same          |     6.301
# Low cost different - Low cost same | 2.136e+06

## Exploratory test of effect of previous choices
# Models
bglmer_last_choice_condition_trimmed <- brm(choice ~ condition + (1 | subject), data = filter(stable_data_trimmed, !is.na(last_choice)), family = bernoulli(link = "logit"), prior = fixed_priors, iter = 20000, cores = ncores, save_all_pars = TRUE)  # model with condition
bglmer_last_choice_full_trimmed <- brm(choice ~ condition + last_choice + (1 | subject), data = filter(stable_data_trimmed, !is.na(last_choice)), family = bernoulli(link = "logit"), prior = fixed_priors, iter = 20000, cores = ncores, save_all_pars = TRUE)  # model with condition and last choice

# Bayes factor
bglmer_last_choice_condition_trimmed_bf <- bf_models(bglmer_last_choice_full_trimmed, denominator = bglmer_last_choice_condition_trimmed)  # BF for condition and last_choice against condition
  # Results: BF = 0.844


####
# Plot figures ------------------------------------------------------------
####

# _ Session means per subject ----------------

# Calculate mean choices per subject, session, and condition
session_means2 <- data %>% 
  group_by(subject, session, condition) %>%  # for each subject, session, and condition
  summarize(choice_means = mean(choice, na.rm = TRUE))   # calculate mean choice proportion

# Plot mean choices per subject, session, and condition
ggplot(session_means2, aes(x = session, y = choice_means, group = condition, col = condition)) +
  geom_line() +  # plot lines
  facet_wrap(~ subject) +  # facet by subject
  labs(x = "Sessions", y = "Proportion choosing larger, later") +  # add axis labels
  theme_legend() +  # use custom theme
  theme(legend.position = c(0.75, 0.1))  # position legend
ggsave(here("figures/sessions_per_subject.png"), width = 12, height = 12)  # save figure

# _ Condition effects on choice -------------------------------------------

# Plot choice proportions per condition and subject
ggplot(condition_subject_stable, aes(x = condition, y = mean_choice, group = subject)) +
  geom_line(color = "grey70") +  # plot lines
  geom_point(data = condition_means2, aes(x = condition, y = choice_mean), size = 4, inherit.aes = FALSE) +  # plot means as points
  geom_linerange(data = condition_means2, aes(x = condition, ymin = choice_lowerci, ymax = choice_upperci), size = 1, inherit.aes = FALSE) +  # plot within-subjects 95% CIs
  labs(x = "Condition", y = "Proportion choosing larger, later") +  # add axis labels
  theme_plots()  # use custom theme
ggsave(here("figures/choice_by_condition_subject.png"), width = 6.75, height = 6)  # save figure

# Plot choice proportions per condition and subject without Robot
condition_subject_stable %>% 
  filter(subject != "Robot") %>%  # remove subject Robot
  ggplot(aes(x = condition, y = mean_choice, group = subject)) +
  geom_line(color = "grey70") +  # plot lines
  geom_point(data = condition_means2_trimmed, aes(x = condition, y = choice_mean), size = 4, inherit.aes = FALSE) +  # plot means as points
  geom_linerange(data = condition_means2_trimmed, aes(x = condition, ymin = choice_lowerci, ymax = choice_upperci), size = 1, inherit.aes = FALSE) +  # plot within-subjects 95% CIs
  labs(x = "Condition", y = "Proportion choosing larger, later") +  # add axis labels
  theme_plots()  # use custom theme
ggsave(here("figures/choice_by_condition_subject_norobot.png"), width = 6.75, height = 6)  # save figure

