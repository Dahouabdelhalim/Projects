###################################################
### Juan F. Duque (jfduque89@gmail.com) & Jeffrey R. Stevens
###	finalized on 2018-03-29
### Summary: This script takes our data and extracts the relevant information 
###   needed to run analyses script, "Duque_et_al_2018_rcode.R"
### Instructions: First, download this script, "Duque_et_al_2018_rcode.R".  
###   Second, make two folders, labeled "data" and "figures", 
###   both lowercase.  Download the "all_data.csv" data file into 
###   the data directory and set the working directory to where the scripts are. 
###   At the R command prompt, type source('Duque_et_al_2018_rcode.R').  This will run the 
###   script and output analyses and figures relevant for Duque et al "Mesotocin influences pinyon jay prosociality". 
###  **********
### Description of the relevant columns for "all_data.csv":
###
###  dataset - either 1 (for Experiment 1, no hormones) or 2 (for Experiment 2, with hormonal manipulation) 
###  donor - a unique numeric identifier for each subject in the study (donor=subject)
###  date - date that session in corresponding row was conducted
###  session - sequential numbering for each subject's session
###  trial - specific trial number within a session
###  partner - location of partner (either N- no partner, L-left sidecage, or R-right sidecage)
###  recipient - if present, a unique numeric identifier for each recipient in the study 
###  payoff_condition - trial type; either Pre-Test (attention trials), Side bias (bias trials), or Prosocial/Altruistic
###  partner_present - either 'Partner absent' or 'Partner present'
###  hormone - hormone condition: NA for Experiment 1; A, B, C, or SAL for Experiment 2 (experimenters were blind to condition)
###  right.left - subject's choice on that trial, either right or left wire
###  response - FOR PRE-TEST (i.e. attention trials): 1 if subject correctly chose side with food, 0 if incorrect
###  response - FOR SIDE BIAS (i.e. bias trials): 1 if subject chose LEFT, 0 if RIGHT
###  response - FOR PROSOCIAL/ALTRUISTIC: IF partner absent, 1 if subject chose LEFT, 0 if RIGHT
###  response - FOR PROSOCIAL/ALTRUISTIC: IF partner present, 1 if subject chose same side as partner, 0 if not.

###################################################
################
# Clear variables and load libraries
################
rm(list=ls()) 
library(BayesFactor)
library(car)
library(papaja)
library(tidyverse)

################
# Functions
################

## Calculates absolute prosocial and altruistic tendency per session
absoluteTendency <- function(df) {
  # Separate side bias and prosocial/altruism trials
  df_side <- subset(df, payoff_condition == "Side bias" & partner_present == "Partner present")   # subset side bias trials
  df_alt_pro <- subset(df, partner_present == "Partner present" & (payoff_condition == "Prosocial" | payoff_condition == "Altruism"))  # subset prosocial and altruism trials
  
  # Calculate matching for side bias trials
  df_side$p_match <- ifelse((df_side$partner == "R" & df_side$response == 1) | (df_side$partner == "L" & df_side$response == 0), 1, ifelse(!is.na(df_side$partner) & !is.na(df_side$response), 0, NA)) # create column in which choice matches partner side
  side_match <- df_side %>% group_by(donor, session, recipient) %>% summarize(bias = mean(p_match, na.rm = TRUE)) # calculate proportion choice of side matching partner for side bias trials
  
  # Calculate matching for prosocial/altruism trials
  alt_pro_match <- df_alt_pro %>% group_by(donor, session, recipient, payoff_condition) %>% summarize(response = mean(response, na.rm = TRUE))  # calculate proportion choice of side matching partner for prosocial/altruism trials
  
  # Calculate difference between prosocial/altruism matching and side bias matching (absolute tendency)
  alt_pro_match_wide <- spread(alt_pro_match, payoff_condition, response) # spread data to wide format
  alt_pro_side_match <- inner_join(alt_pro_match_wide, side_match, by = c("donor", "session", "recipient")) # join prosocial, altruism, and side bias data
  alt_pro_side_match <- alt_pro_side_match %>% 
    mutate(bias_altruism = ifelse(is.na(Altruism), NA, bias),  # create bias column with NAs
      prosocial_absolute = Prosocial - bias,  # calculate difference between prosocial and side bias matching
      altruism_absolute = Altruism - bias  # calculate difference between altruism and side bias matching
    )
  alt_pro_side_match  # return data frame
}

## Convert model BIC values to Bayes factors
bic_bf <- function(null_model, alternative_model) {
  new_bf <- exp((BIC(null_model) - BIC(alternative_model)) / 2) # convert BICs to Bayes factor
  names(new_bf) <- NULL   # remove BIC label
  return(new_bf)  # return Bayes factor
}

################
# Prepare data
################

all_data <- read_csv("data/all_data.csv") %>%   # read in data
  mutate(payoff_condition = Recode(payoff_condition, "'Altruistic'='Altruism'"),  # recode altruistic to altruism
    payoff_condition = factor(payoff_condition, levels = c("Pre-test", "Side bias", "Prosocial", "Altruism")),  # reorder payoff conditions
    hormone = Recode(hormone, "'A'='High MT';'B'='Low MT';'C'='Saline'") # recode hormone conditions
  )

experiment1data <- subset(all_data, dataset == 1)  # extract experiment 1 data
experiment2data <- subset(all_data, dataset == 2)  # extract experiment 2 data
experiment2data_high <- subset(experiment2data, hormone=="High MT")  # extract high MT data
experiment2data_low <- subset(experiment2data, hormone=="Low MT")  # extract low MT data
experiment2data_saline <- subset(experiment2data, hormone=="Saline")  # extract saline data

col.blind.colors <- c("#D55E00", "#0072B2", "#009E73", "#E69F00", "#F0E442", "#CC79A7", "#56B4E9", "black", "grey50")  # assign colorblind safe colors

################
# Experiment 1: Behavioral data
################

## Prepare data
tendency_wide_all1 <- absoluteTendency(experiment1data) # calculate absolute tendency
tendency_wide1 <- tendency_wide_all1 %>% 
  group_by(donor) %>%  # for each donor
  summarize(bias_prosocial = mean(bias, na.rm = TRUE),  # calculate mean bias matching
    bias_altruism = mean(bias_altruism, na.rm = TRUE), 
    Prosocial_Absolute = mean(prosocial_absolute, na.rm = TRUE),  # calculate mean prosocial absolute tendency
    Altruism_Absolute = mean(altruism_absolute, na.rm = TRUE),  # calculate mean altruism absolute tendency
    Prosocial_Weighted = ifelse(Prosocial_Absolute > 0, Prosocial_Absolute / (1 - bias_prosocial), ifelse(Prosocial_Absolute < 0, Prosocial_Absolute / bias_prosocial, 0)),   # calculate prosocial weighted tendency
    Altruism_Weighted = ifelse(Altruism_Absolute > 0, Altruism_Absolute / (1 - bias_altruism), ifelse(Altruism_Absolute < 0, Altruism_Absolute / bias_altruism, 0)))  # calculate altruism weighted tendency
tendency_long1 <- tendency_wide1 %>% 
  gather(key = condition, value = tendency, Prosocial_Absolute, Altruism_Absolute, Prosocial_Weighted, Altruism_Weighted) %>%  # convert wide data frame to long
  separate(condition, c("payoff_condition", "measure")) %>% # separate payoff condition and measure
  mutate(payoff_condition = factor(payoff_condition, levels = c("Prosocial", "Altruism"))) # reorder payoff condition levels

# Absolute tendency
absolute_tendency1 <- filter(tendency_long1, measure == "Absolute")  # extract absolute tendency data
absolute_tendency_summary1 <- absolute_tendency1 %>% 
  group_by(payoff_condition) %>% # for each payoff condition
  summarize(tend = mean(tendency)) # calculate mean tendency
absolute_tendency_summary1$ci <- wsci(data = filter(tendency_long1, measure == "Absolute"), id = "donor", factors = c("payoff_condition"), dv = "tendency")$tendency # calculate within-subjects 95% CI per payoff condition
absolute_tendency_summary1$lower <- absolute_tendency_summary1$tend - absolute_tendency_summary1$ci  # create column of upper CIs
absolute_tendency_summary1$upper <- absolute_tendency_summary1$tend + absolute_tendency_summary1$ci  # create column of lower CIs

# Weighted tendency
weighted_tendency1 <- filter(tendency_long1, measure == "Weighted")  # extract weighted tendency data
weighted_tendency_summary1 <- weighted_tendency1 %>% 
  group_by(payoff_condition) %>%  # for each payoff condition
  summarize(tend = mean(tendency)) # calculate mean tendency
weighted_tendency_summary1$ci <- wsci(data = filter(tendency_long1, measure == "Weighted"), id = "donor", factors = c("payoff_condition"), dv = "tendency")$tendency # calculate within-subjects 95% CI per payoff condition
weighted_tendency_summary1$lower <- weighted_tendency_summary1$tend - weighted_tendency_summary1$ci  # create column of upper CIs
weighted_tendency_summary1$upper <- weighted_tendency_summary1$tend + weighted_tendency_summary1$ci  # create column of lower CIs

## Compare means to zero
prosocial_absolute_t1 <- t.test(tendency_wide1$Prosocial_Absolute, mu = 0)   # frequentist t-test greater than 0
prosocial_absolute_bft1 <- ttestBF(tendency_wide1$Prosocial_Absolute, mu = 0)   # Bayesian t-test greater than 0
prosocial_absolute_bf1 <- round(extractBF(prosocial_absolute_bft1)$bf, 1) # extract Bayes factor

altruism_absolute_t1 <- t.test(tendency_wide1$Altruism_Absolute, mu = 0)   # frequentist t-test greater than 0
altruism_absolute_bft1 <- ttestBF(tendency_wide1$Altruism_Absolute, mu = 0)   # Bayesian t-test greater than 0
altruism_absolute_bf1 <- round(extractBF(altruism_absolute_bft1)$bf, 1) # extract Bayes factor

# Compare means to zero
prosocial_weighted_t1 <- t.test(tendency_wide1$Prosocial_Weighted, mu = 0)   # frequentist t-test greater than 0
prosocial_weighted_bft1 <- ttestBF(tendency_wide1$Prosocial_Weighted, mu = 0)   # Bayesian t-test greater than 0
prosocial_weighted_bf1 <- round(extractBF(prosocial_weighted_bft1)$bf, 1) # extract Bayes factor

altruism_weighted_t1 <- t.test(tendency_wide1$Altruism_Weighted, mu = 0)   # frequentist t-test greater than 0
altruism_weighted_bft1 <- ttestBF(tendency_wide1$Altruism_Weighted, mu = 0)   # Bayesian t-test greater than 0
altruism_weighted_bf1 <- round(extractBF(altruism_weighted_bft1)$bf, 1) # extract Bayes factor

## Plot absolute tendency
ggplot(absolute_tendency1, aes(x = payoff_condition, y = tendency)) +
  geom_jitter(shape = 1, size = 8, width = 0.015, height = 0.02) +   # plot subject data points
  geom_point(data = absolute_tendency_summary1, aes(x = payoff_condition, y = tend), size = 12, shape = 18) + # plot mean per condition
  geom_linerange(data = absolute_tendency_summary1, aes(x = payoff_condition, ymin = lower, ymax = upper), inherit.aes = FALSE, size = 2) +  # plot CIs per condition
  geom_hline(yintercept = 0) + # plot chance line
  geom_text(aes(x = 1, y = -0.1, label = paste("BF =", prosocial_absolute_bf1)), size = 15) + # print Bayes factors
  geom_text(aes(x = 2, y = -0.1, label = paste("BF =", altruism_absolute_bf1)), size = 15) + # print Bayes factors
  labs(x = "Payoff condition", y = "Absolute tendency") + # label axes
  theme_classic() + # use classic theme
  theme(axis.title = element_text(size = 55), axis.text = element_text(size = 40), plot.margin = unit(c(5, 1, 3, 1), "mm")) # resize fonts
ggsave("figures/absolute_tendency1.pdf", width = 12, height = 12)  # create PDF file

## Plot weighted tendency
ggplot(weighted_tendency1, aes(x = payoff_condition, y = tendency)) +
  geom_jitter(shape = 1, size = 8, width = 0.01, height = 0.025) +   # plot subject data points
  geom_point(data = weighted_tendency_summary1, aes(x = payoff_condition, y = tend), size = 12, shape = 18) + # plot mean per condition
  geom_linerange(data = weighted_tendency_summary1, aes(x = payoff_condition, ymin = lower, ymax = upper), inherit.aes = FALSE, size = 2) +  # plot CIs per condition
  geom_hline(yintercept = 0) + # plot chance line
  geom_text(aes(x = 1, y = -0.2, label = paste("BF =", prosocial_weighted_bf1)), size = 15) + # print Bayes factors
  geom_text(aes(x = 2, y = -0.2, label = paste("BF =", altruism_weighted_bf1)), size = 15) + # print Bayes factors
  labs(x = "Payoff condition", y = "Weighted tendency") + # label axes
  theme_classic() + # use classic theme
  theme(axis.title = element_text(size = 55), axis.text = element_text(size = 40), plot.margin = unit(c(5, 1, 3, 1), "mm")) # resize fonts
ggsave("figures/weighted_tendency1.pdf", width = 12, height = 12)  # create PDF file

################
# Experiment 2: Hormone manipulation
# Same Trial Types (bias, altruism, and prosocial)
# But each subject received one of three hormone administrations.
# Three hormone conditions: High MT (high dose of mesotocin given to subject), Low MT, or Saline
################

## Prepare data
# High MT data
tendency_high_wide_all2 <- absoluteTendency(experiment2data_high) # calculate prosocial tendency correction
tendency_high_wide2 <- tendency_high_wide_all2 %>% 
  group_by(donor) %>%  # for each donor
  summarize(bias_prosocial = mean(bias, na.rm = TRUE),  # calculate mean bias matching
    bias_altruism = mean(bias_altruism, na.rm = TRUE),   # calculate mean altruism bias matching
    Prosocial_Absolute = mean(prosocial_absolute, na.rm = TRUE),  # calculate mean prosocial absolute tendency
    Altruism_Absolute = mean(altruism_absolute, na.rm = TRUE),  # calculate mean altruism absolute tendency
    Prosocial_Weighted = ifelse(Prosocial_Absolute > 0, Prosocial_Absolute / (1 - bias_prosocial), ifelse(Prosocial_Absolute < 0, Prosocial_Absolute / bias_prosocial, 0)),   # calculate prosocial weighted tendency
    Altruism_Weighted = ifelse(Altruism_Absolute > 0, Altruism_Absolute / (1 - bias_altruism), ifelse(Altruism_Absolute < 0, Altruism_Absolute / bias_altruism, 0)))  # calculate altruism weighted tendency
tendency_high_long2 <- tendency_high_wide2 %>% 
  gather(key = condition, value = tendency, Prosocial_Absolute, Altruism_Absolute, Prosocial_Weighted, Altruism_Weighted) %>%  # convert wide data frame to long
  separate(condition, c("payoff_condition", "measure")) %>% # separate payoff condition and measure
  mutate(payoff_condition = factor(payoff_condition, levels = c("Prosocial", "Altruism"))) # reorder payoff condition levels

# Low MT data
tendency_low_wide_all2 <- absoluteTendency(experiment2data_low) # calculate prosocial tendency correction
tendency_low_wide2 <- tendency_low_wide_all2 %>% 
  group_by(donor) %>%  # for each donor
  summarize(bias_prosocial = mean(bias, na.rm = TRUE),  # calculate mean bias matching
    bias_altruism = mean(bias_altruism, na.rm = TRUE),   # calculate mean altruism bias matching
    Prosocial_Absolute = mean(prosocial_absolute, na.rm = TRUE),  # calculate mean prosocial absolute tendency
    Altruism_Absolute = mean(altruism_absolute, na.rm = TRUE),  # calculate mean altruism absolute tendency
    Prosocial_Weighted = ifelse(Prosocial_Absolute > 0, Prosocial_Absolute / (1 - bias_prosocial), ifelse(Prosocial_Absolute < 0, Prosocial_Absolute / bias_prosocial, 0)),   # calculate prosocial weighted tendency
    Altruism_Weighted = ifelse(Altruism_Absolute > 0, Altruism_Absolute / (1 - bias_altruism), ifelse(Altruism_Absolute < 0, Altruism_Absolute / bias_altruism, 0)))  # calculate altruism weighted tendency
tendency_low_long2 <- tendency_low_wide2 %>% 
  gather(key = condition, value = tendency, Prosocial_Absolute, Altruism_Absolute, Prosocial_Weighted, Altruism_Weighted) %>%  # convert wide data frame to long
  separate(condition, c("payoff_condition", "measure")) %>% # separate payoff condition and measure
  mutate(payoff_condition = factor(payoff_condition, levels = c("Prosocial", "Altruism"))) # reorder payoff condition levels

# Saline data
tendency_saline_wide_all2 <- absoluteTendency(experiment2data_saline) # calculate prosocial tendency correction
tendency_saline_wide2 <- tendency_saline_wide_all2 %>% 
  group_by(donor) %>%  # for each donor
  summarize(bias_prosocial = mean(bias, na.rm = TRUE),  # calculate mean bias matching
    bias_altruism = mean(bias_altruism, na.rm = TRUE),   # calculate mean altruism bias matching
    Prosocial_Absolute = mean(prosocial_absolute, na.rm = TRUE),  # calculate mean prosocial absolute tendency
    Altruism_Absolute = mean(altruism_absolute, na.rm = TRUE),  # calculate mean altruism absolute tendency
    Prosocial_Weighted = ifelse(Prosocial_Absolute > 0, Prosocial_Absolute / (1 - bias_prosocial), ifelse(Prosocial_Absolute < 0, Prosocial_Absolute / bias_prosocial, 0)),   # calculate prosocial weighted tendency
    Altruism_Weighted = ifelse(Altruism_Absolute > 0, Altruism_Absolute / (1 - bias_altruism), ifelse(Altruism_Absolute < 0, Altruism_Absolute / bias_altruism, 0)))  # calculate altruism weighted tendency
tendency_saline_long2 <- tendency_saline_wide2 %>% 
  gather(key = condition, value = tendency, Prosocial_Absolute, Altruism_Absolute, Prosocial_Weighted, Altruism_Weighted) %>%  # convert wide data frame to long
  separate(condition, c("payoff_condition", "measure")) %>% # separate payoff condition and measure
  mutate(payoff_condition = factor(payoff_condition, levels = c("Prosocial", "Altruism"))) # reorder payoff condition levels

# Combine data
tendency_long2 <- bind_rows(tendency_high_long2, tendency_low_long2, tendency_saline_long2) %>%   # combine hormone data
mutate(payoff_condition = factor(payoff_condition, levels = c("Prosocial", "Altruism")), # reorder payoff condition levels
    hormone = rep(c("High MT", "Low MT", "Saline"), each = dim(tendency_high_long2)[1])
  )
       
# Absolute tendency
absolute_tendency2 <- filter(tendency_long2, measure == "Absolute")  # extract absolute tendency dat
absolute_tendency_summary2 <- absolute_tendency2 %>% 
  group_by(payoff_condition, hormone) %>% # for each payoff condition
  summarize(tend = mean(tendency)) # calculate mean tendency
absolute_tendency_summary2$ci <- wsci(data = filter(tendency_long2, measure == "Absolute"), id = "donor", factors = c("payoff_condition", "hormone"), dv = "tendency")$tendency # calculate within-subjects 95% CI per payoff condition
absolute_tendency_summary2$lower <- absolute_tendency_summary2$tend - absolute_tendency_summary2$ci  # create column of upper CIs
absolute_tendency_summary2$upper <- absolute_tendency_summary2$tend + absolute_tendency_summary2$ci  # create column of lower CIs

# Weighted tendency
weighted_tendency2 <- filter(tendency_long2, measure == "Weighted")  # extract weighted tendency data
weighted_tendency_summary2 <- weighted_tendency2 %>% 
  group_by(payoff_condition, hormone) %>%  # for each payoff condition
  summarize(tend = mean(tendency)) # calculate mean tendency
weighted_tendency_summary2$ci <- wsci(data = filter(tendency_long2, measure == "Weighted"), id = "donor", factors = c("payoff_condition", "hormone"), dv = "tendency")$tendency # calculate within-subjects 95% CI per payoff condition
weighted_tendency_summary2$lower <- weighted_tendency_summary2$tend - weighted_tendency_summary2$ci  # create column of upper CIs
weighted_tendency_summary2$upper <- weighted_tendency_summary2$tend + weighted_tendency_summary2$ci  # create column of lower CIs

## Compare means to zero
# High MT
prosocial_absolute_high_t2 <- t.test(tendency_high_wide2$Prosocial_Absolute, mu = 0)   # frequentist t-test greater than 0
prosocial_absolute_high_bft2 <- ttestBF(tendency_high_wide2$Prosocial_Absolute, mu = 0)   # Bayesian t-test greater than 0
prosocial_absolute_high_bf2 <- round(extractBF(prosocial_absolute_high_bft2)$bf, 1) # extract Bayes factor

altruism_absolute_high_t2 <- t.test(tendency_high_wide2$Altruism_Absolute, mu = 0)   # frequentist t-test greater than 0
altruism_absolute_high_bft2 <- ttestBF(tendency_high_wide2$Altruism_Absolute, mu = 0)   # Bayesian t-test greater than 0
altruism_absolute_high_bf2 <- round(extractBF(altruism_absolute_high_bft2)$bf, 1) # extract Bayes factor

prosocial_weighted_high_t2 <- t.test(tendency_high_wide2$Prosocial_Weighted, mu = 0)   # frequentist t-test greater than 0
prosocial_weighted_high_bft2 <- ttestBF(tendency_high_wide2$Prosocial_Weighted, mu = 0)   # Bayesian t-test greater than 0
prosocial_weighted_high_bf2 <- round(extractBF(prosocial_weighted_high_bft2)$bf, 1) # extract Bayes factor

altruism_weighted_high_t2 <- t.test(tendency_high_wide2$Altruism_Weighted, mu = 0)   # frequentist t-test greater than 0
altruism_weighted_high_bft2 <- ttestBF(tendency_high_wide2$Altruism_Weighted, mu = 0)   # Bayesian t-test greater than 0
altruism_weighted_high_bf2 <- round(extractBF(altruism_weighted_high_bft2)$bf, 1) # extract Bayes factor

# Low MT
prosocial_absolute_low_t2 <- t.test(tendency_low_wide2$Prosocial_Absolute, mu = 0)   # frequentist t-test greater than 0
prosocial_absolute_low_bft2 <- ttestBF(tendency_low_wide2$Prosocial_Absolute, mu = 0)   # Bayesian t-test greater than 0
prosocial_absolute_low_bf2 <- round(extractBF(prosocial_absolute_low_bft2)$bf, 1) # extract Bayes factor

altruism_absolute_low_t2 <- t.test(tendency_low_wide2$Altruism_Absolute, mu = 0)   # frequentist t-test greater than 0
altruism_absolute_low_bft2 <- ttestBF(tendency_low_wide2$Altruism_Absolute, mu = 0)   # Bayesian t-test greater than 0
altruism_absolute_low_bf2 <- round(extractBF(altruism_absolute_low_bft2)$bf, 1) # extract Bayes factor

prosocial_weighted_low_t2 <- t.test(tendency_low_wide2$Prosocial_Weighted, mu = 0)   # frequentist t-test greater than 0
prosocial_weighted_low_bft2 <- ttestBF(tendency_low_wide2$Prosocial_Weighted, mu = 0)   # Bayesian t-test greater than 0
prosocial_weighted_low_bf2 <- round(extractBF(prosocial_weighted_low_bft2)$bf, 1) # extract Bayes factor

altruism_weighted_low_t2 <- t.test(tendency_low_wide2$Altruism_Weighted, mu = 0)   # frequentist t-test greater than 0
altruism_weighted_low_bft2 <- ttestBF(tendency_low_wide2$Altruism_Weighted, mu = 0)   # Bayesian t-test greater than 0
altruism_weighted_low_bf2 <- round(extractBF(altruism_weighted_low_bft2)$bf, 1) # extract Bayes factor

# Saline  
prosocial_absolute_saline_t2 <- t.test(tendency_saline_wide2$Prosocial_Absolute, mu = 0)   # frequentist t-test greater than 0
prosocial_absolute_saline_bft2 <- ttestBF(tendency_saline_wide2$Prosocial_Absolute, mu = 0)   # Bayesian t-test greater than 0
prosocial_absolute_saline_bf2 <- round(extractBF(prosocial_absolute_saline_bft2)$bf, 1) # extract Bayes factor

altruism_absolute_saline_t2 <- t.test(tendency_saline_wide2$Altruism_Absolute, mu = 0)   # frequentist t-test greater than 0
altruism_absolute_saline_bft2 <- ttestBF(tendency_saline_wide2$Altruism_Absolute, mu = 0)   # Bayesian t-test greater than 0
altruism_absolute_saline_bf2 <- round(extractBF(altruism_absolute_saline_bft2)$bf, 1) # extract Bayes factor

prosocial_weighted_saline_t2 <- t.test(tendency_saline_wide2$Prosocial_Weighted, mu = 0)   # frequentist t-test greater than 0
prosocial_weighted_saline_bft2 <- ttestBF(tendency_saline_wide2$Prosocial_Weighted, mu = 0)   # Bayesian t-test greater than 0
prosocial_weighted_saline_bf2 <- round(extractBF(prosocial_weighted_saline_bft2)$bf, 1) # extract Bayes factor

altruism_weighted_saline_t2 <- t.test(tendency_saline_wide2$Altruism_Weighted, mu = 0)   # frequentist t-test greater than 0
altruism_weighted_saline_bft2 <- ttestBF(tendency_saline_wide2$Altruism_Weighted, mu = 0)   # Bayesian t-test greater than 0
altruism_weighted_saline_bf2 <- round(extractBF(altruism_weighted_saline_bft2)$bf, 1) # extract Bayes factor

# Add Bayes factors to summary tables
absolute_tendency_summary2$bf <- c(prosocial_absolute_high_bf2, prosocial_absolute_low_bf2, prosocial_absolute_saline_bf2, altruism_absolute_high_bf2, altruism_absolute_low_bf2, altruism_absolute_saline_bf2)
absolute_tendency_summary2$bayes <- paste("BF =", absolute_tendency_summary2$bf)
weighted_tendency_summary2$bf <- c(prosocial_weighted_high_bf2, prosocial_weighted_low_bf2, prosocial_weighted_saline_bf2, altruism_weighted_high_bf2, altruism_weighted_low_bf2, altruism_weighted_saline_bf2)
weighted_tendency_summary2$bayes <- paste("BF =", weighted_tendency_summary2$bf)

## Plot absolute tendency
ggplot(absolute_tendency2, aes(x = payoff_condition, y = tendency)) +
  facet_wrap(~hormone) +  # separate panels for each payoff condition
  geom_jitter(shape = 1, size = 8, width = 0.015, height = 0.025) +   # plot subject data points
  geom_point(data = absolute_tendency_summary2, aes(x = payoff_condition, y = tend), size = 12, shape = 18) + # plot mean per condition
  geom_linerange(data = absolute_tendency_summary2, aes(x = payoff_condition, ymin = lower, ymax = upper), inherit.aes = FALSE, size = 2) +  # plot CIs per condition
  geom_hline(yintercept = 0) + # plot chance line
  geom_text(data = absolute_tendency_summary2, aes(x = payoff_condition, y = -0.4, label = bayes), size = 12, color = "black") +
  labs(x = "Payoff condition", y = "Absolute tendency") + # label axes
  scale_color_manual(name="", values = col.blind.colors) +  # define line colors 
  theme_classic() + # use classic theme
  theme(axis.title = element_text(size = 55), axis.text = element_text(size = 40), strip.text.x = element_text(size = 40, margin = margin(3,0,4,0, "mm")), plot.margin = unit(c(5, 1, 3, 1), "mm"), legend.position = "null") # resize fonts
ggsave("figures/absolute_tendency2.pdf", width = 20, height = 12)  # create PDF file

## Plot weighted tendency
ggplot(weighted_tendency2, aes(x = payoff_condition, y = tendency)) +
  facet_wrap(~hormone) +  # separate panels for each payoff condition
  geom_jitter(shape = 1, size = 8, width = 0.02, height = 0.035) +   # plot subject data points
  geom_point(data = weighted_tendency_summary2, aes(x = payoff_condition, y = tend), size = 12, shape = 18) + # plot mean per condition
  geom_linerange(data = weighted_tendency_summary2, aes(x = payoff_condition, ymin = lower, ymax = upper), inherit.aes = FALSE, size = 2) +  # plot CIs per condition
  geom_hline(yintercept = 0) + # plot chance line
  geom_text(data = weighted_tendency_summary2, aes(x = payoff_condition, y = -0.8, label = bayes), size = 12, color = "black") +
  labs(x = "Payoff condition", y = "Weighted tendency") + # label axes
  scale_color_manual(name="", values = col.blind.colors) +  # define line colors 
  theme_classic() + # use classic theme
  theme(axis.title = element_text(size = 55), axis.text = element_text(size = 40), strip.text.x = element_text(size = 40, margin = margin(3,0,4,0, "mm")), plot.margin = unit(c(5, 1, 3, 1), "mm"), legend.position = "null") # resize fonts
ggsave("figures/weighted_tendency2.pdf", width = 20, height = 12)  # create PDF file

##################################
# Supplementary
##################################

## Experiment 1 
# Matching values for each subject (donor)
matching_subjects1 <- tendency_wide_all1 %>% 
  group_by(donor) %>%  # for each donor
  summarize(Bias = round(mean(bias), 2), Prosocial = round(mean(Prosocial, na.rm = TRUE), 2), Altruism = round(mean(Altruism, na.rm = TRUE), 2)) %>%   # calculate mean matching for bias, prosocial, and altruism trials
  rename(Donor = donor)  # rename donor column

write.table(matching_subjects1, file = "figures/Subject_RawMatching1.txt", sep = ",", quote = FALSE, row.names = FALSE)  # export table

matching_subjects_long1 <- matching_subjects1 %>% 
  gather(key = payoff_condition, value = tendency, Prosocial, Altruism, Bias) %>%  # convert from wide to long
  mutate(Donor = factor(Donor),   # convert to factor
    payoff_condition = factor(payoff_condition, levels = c("Prosocial", "Bias", "Altruism"))  # reorder factor levels
  )

ggplot(matching_subjects_long1, aes(payoff_condition, tendency, group = Donor)) +  # plot matching values per condition
  geom_line(aes(col = Donor), size = 1) + # plot separate lines for donors
  stat_summary(aes(payoff_condition, tendency), fun.data = mean_cl_boot, inherit.aes = FALSE, size = 2.5, shape = 18)+
  labs(x = "Payoff condition", y = "Matching") + # label axes
  theme_classic() + # use classic theme
  theme(axis.title = element_text(size = 55), axis.text = element_text(size = 40), strip.text.x = element_text(size = 40, margin = margin(3, 0, 4, 0, "mm")), plot.margin = unit(c(5, 1, 3, 1), "mm"), legend.position = "null") # resize fonts
ggsave("figures/matching1.pdf", width = 12, height = 12) #Create PDF file

# Tendency values for each subject (donor)
tendency_subjects1 <- tendency_wide1 %>% 
  select(-bias_prosocial, -bias_altruism) %>% # remove bias column
  rename(Donor = donor) %>%  # rename donor column
  mutate(Prosocial_Absolute = round(Prosocial_Absolute, 2), Altruism_Absolute = round(Altruism_Absolute, 2), Prosocial_Weighted = round(Prosocial_Weighted, 2), Altruism_Weighted = round(Altruism_Weighted, 2))  # round values

write.table(tendency_subjects1, file = "figures/Subject_Tendency1.txt", sep = ",", quote = FALSE, row.names = FALSE) # export table

# Combine matching and tendency
matching_tendency1 <- matching_subjects1 %>% 
  inner_join(tendency_subjects1, by = "Donor")  # join matching and tendency data

write.table(matching_tendency1, file = "figures/Subject_RawMatching_Tendency1.txt", sep = ",", quote = FALSE, row.names = FALSE)  # export table

## Experiment 2 
# Matching values for each subject (donor)
matching_subjects2_high <- tendency_high_wide_all2 %>% 
  group_by(donor) %>%   # for each donor
  summarize(Bias = round(mean(bias), 2), Prosocial = round(mean(Prosocial, na.rm = TRUE), 2), Altruism = round(mean(Altruism, na.rm = TRUE), 2)) %>%   # calculate mean matching for bias, prosocial, and altruism trials
  mutate(Hormone = "High MT")  # add hormone column

matching_subjects2_low <- tendency_low_wide_all2 %>% 
  group_by(donor) %>%   # for each donor
  summarize(Bias = round(mean(bias), 2), Prosocial = round(mean(Prosocial, na.rm = TRUE), 2), Altruism = round(mean(Altruism, na.rm = TRUE), 2)) %>%   # calculate mean matching for bias, prosocial, and altruism trials
  mutate(Hormone = "Low MT")  # add hormone column

matching_subjects2_sal <- tendency_saline_wide_all2 %>% 
  group_by(donor) %>%   # for each donor
  summarize(Bias = round(mean(bias), 2), Prosocial = round(mean(Prosocial, na.rm = TRUE), 2), Altruism = round(mean(Altruism, na.rm = TRUE), 2)) %>%   # calculate mean matching for bias, prosocial, and altruism trials
  mutate(Hormone = "Saline")  # add hormone column

matching_subjects2 <- bind_rows(matching_subjects2_high, matching_subjects2_low, matching_subjects2_sal) %>%  # combine data
  arrange(donor) %>% # arrange by donor
  rename(Donor = donor)  # rename donor column

write.table(matching_subjects2, file = "figures/Subject_RawMatching2.txt", sep = ",", quote = FALSE, row.names = FALSE)  # export data

matching_subjects_long2 <- matching_subjects2 %>% 
  gather(key = payoff_condition, value = tendency, Prosocial, Altruism, Bias) %>%  # convert from wide to long
  mutate(Donor = factor(Donor),   # convert to factor
    payoff_condition = factor(payoff_condition, levels = c("Prosocial", "Bias", "Altruism"))  # reorder factor levels
  )

ggplot(matching_subjects_long2, aes(payoff_condition, tendency, group = Donor)) +  # plot matching values per condition
  geom_line(aes(col = Donor), size = 1) + # plot separate lines for donors
  stat_summary(aes(payoff_condition, tendency), fun.data = mean_cl_boot, inherit.aes = FALSE, size = 2.5, shape = 18)+
  facet_wrap(~Hormone) +  # facet by hormone
  labs(x = "Payoff condition", y = "Matching") + # label axes
  theme_classic() + # use classic theme
  theme(axis.title = element_text(size = 55), axis.text = element_text(size = 40), strip.text.x = element_text(size = 40, margin = margin(3, 0, 4, 0, "mm")), plot.margin = unit(c(5, 1, 3, 1), "mm"), legend.position = "null") # resize fonts
ggsave("figures/matching2.pdf", width = 20, height = 12) # create PDF file

# Tendency values for each subject (donor)
tendency_subjects2 <- bind_rows(tendency_high_wide2, tendency_low_wide2, tendency_saline_wide2) %>% 
  select(-bias_prosocial, -bias_altruism) %>%  # remove bias column
  mutate(Hormone = rep(c("High MT", "Low MT", "Saline"), each = dim(tendency_high_wide2)[1])) %>% 
  rename(Donor = donor) %>%  # rename donor column
  mutate(Prosocial_Absolute = round(Prosocial_Absolute, 2), Altruism_Absolute = round(Altruism_Absolute, 2), Prosocial_Weighted = round(Prosocial_Weighted, 2), Altruism_Weighted = round(Altruism_Weighted, 2))  # round values

write.table(tendency_subjects2, file = "figures/Subject_Tendency2.txt", sep = ",", quote = FALSE, row.names = FALSE)  # export data

# Combine matching and tendency
matching_tendency2 <- matching_subjects2 %>% 
  inner_join(tendency_subjects2, by = c("Donor", "Hormone")) %>%  # join matching and tendency data
  select(Donor, Hormone, everything())  # reorder columns

write.table(matching_tendency2, file = "figures/Subject_RawMatching_Tendency2.txt", sep = ",", quote = FALSE, row.names = FALSE)  # export data


## Matching values received by partners in Experiment 1
matching_recipients1 <- tendency_wide_all1 %>% 
  group_by(recipient) %>%  # for each recipient
  summarize(Prosocial_Absolute = round(mean(prosocial_absolute), 2), Altruistic_Absolute = round(mean(altruism_absolute, na.rm = TRUE), 2))  # calculate mean absolute tendency

write.table(matching_recipients1, file = "figures/Recipient_Matching.txt", sep = ",", quote = FALSE, row.names = FALSE)
