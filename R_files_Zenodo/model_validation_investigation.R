
library(tidyverse)
library(mgcv)
library(rlang)
library(furrr)

# Setup parallel with furrr
options(mc.cores = 16)
future::plan(multisession)

load(file = "Data/Diuron_photo_model_and_data2.RDA") # Diuron photosynthesis models and data
load(file = "Data/Diuron_growth_model_and_data2.RDA") # Diuron growth models and data

source("Scripts/model_validation_functions.R")

# Remove the missing block so we have the same number of reps as the real data
Diuron_growth_data <- Diuron_growth_data %>% filter(block != "1")


model_validation_sim <- function(model, error_str, interaction_str, return_data = FALSE){
  
  # ==========================================
  # This function generates a single simulation for a particular endpoint (model) with a 
  # specified amount of error (error_str) with a certain interaction strength (interaction_str).
  # Within this simulation a fake dataset with a random component is generated and
  # we see how well the modelling approach and the raw statistic perform at capturing
  # the known interaction.
  # 
  # return_data = TRUE: returns the a list with the dataset plus the computed statistics (used for visual investigation)
  #
  # ==========================================
  
  # Expression to generate the simulation for the photosynthesis endpoint
  gen_photo_dat_expr <- expr({

    fake_block <- rnorm(n = 4, sd = error_str, mean = 0) #0.05, 0.07, 0.03
    sample_error <- rnorm(300, mean = 0, sd = error_str) #0.05, 0.07, 0.03
    mean_t0 <-  mean(Diuron_photo_data$t0)
    
    # Generate the fake data
    dat_fake <- within(Diuron_photo_data, {
      
      # With random and fixed effects
      Yield <- t0 + hours/720*2.5 -(80-Light_num)/(800)*2.5 - Diuron_num/30*2.5 + 
        ifelse((80 - Light_num)*Diuron_num == 0, 0, 1 )*interaction_str + ifelse(sign(interaction_str) == -1, abs(interaction_str), 0)  +
        fake_block[as.numeric(block)]  + sample_error
        
      # With fixed effects only
      True_val <- mean_t0 + hours/720*2.5 -(80-Light_num)/(800)*2.5 - Diuron_num/30*2.5 + 
        ifelse((80 - Light_num)*Diuron_num == 0, 0, 1 )*interaction_str + ifelse(sign(interaction_str) == -1, abs(interaction_str), 0)  
    }) %>% ungroup()
    
    # Run the gam model
    m1_diuron_yield <- gam(Yield ~ hours_fact +
                             offset(t0)+
                             s(Diuron_num, Light_num, k = 15, by = hours_fact) +
                             #s(sample_id, bs = "re")+
                             #s(Light_num, k = 3, by = hours_fact) +
                             #s(Diuron_num, k = 5, by =hours_fact) +
                             s(block, bs = "re"),
                           data = dat_fake)
    
    # Compute the CIs using empirical bayes
    Diuron_CItrend3 <- compute_IR_and_CIs(dat = filter(dat_fake, block == 2) %>% mutate(t0 = mean_t0), 
                                          model = m1_diuron_yield, response = "photosynthesis", 
                                          variable = "Diuron") %>% 
      select(X2.5., X50., X97.5., hours, Light_num, Diuron_num)
    
  })
  
  
  # Expression to generate a dataset for the growth endpoint
  gen_growth_dat_expr <- expr({
    
    fake_block <- rnorm(n = 4, sd = error_str, mean = 0)
    sample_error <- rnorm(225, sd = error_str)
    mean_t0 <-  mean(Diuron_growth_data$t0)
    
    # Generate the fake data
    dat_fake <- within(Diuron_growth_data, {
      
      # With random and fixed effects
      celld <- t0 * exp(hours/72/2 - (80 -Light_num)/(80*2) - Diuron_num/(3*2) + 
                          ifelse((80 - Light_num)*Diuron_num == 0, 0, 1 )*interaction_str + 
                          fake_block[as.numeric(block)] + sample_error)
      # With fixed effects only
      True_val <- mean_t0 * exp(hours/72/2 - (80 -Light_num)/(80*2) - Diuron_num/(3*2) + 
                                  ifelse((80 - Light_num)*Diuron_num == 0, 0, 1 )*interaction_str)
    })
    
    # Run the gam model
    m1_diuron_light <- gam(log(celld) ~ 0 + hours_fact + offset(log(t0)) +
                             s(Diuron_num, Light_num, k = 15, by = hours_fact) +
                             #s(sample_id, bs = "re")+
                             #s(Light_num, k = 3, by = hours_fact) +
                             #s(Diuron_num, k = 5, by =hours_fact) +
                             s(block, bs = "re"),
                           data = dat_fake)
    
    # Compute the CIs using empirical bayes
    Diuron_CItrend3 <- compute_IR_and_CIs(dat = filter(dat_fake, block == 2) %>% mutate(t0 = mean_t0), 
                                          model = m1_diuron_light, response = "growth", 
                                          variable = "Diuron") %>% 
      select(X2.5., X50., X97.5., hours, Light_num, Diuron_num)
  })
  
  # Select the expression for the desired endpoint
  data_eval <- switch(model, "photosynthesis" = gen_photo_dat_expr, "growth" = gen_growth_dat_expr)
  
  # Evaluate the endpoint expression
  eval(data_eval)


  # Compute the real statistic based on the response with no random terms
  true_statistic <- dat_fake %>% filter(block == 2) %>% 
    mutate(Treatment = case_when(Light_num == 80 & Diuron_num == 0 ~ "Control",
                                 Light_num == 80  ~ "Diuron_only",
                                 Diuron_num == 0  ~ "Light_only",
                                 TRUE ~ "Light_diuron")) %>% 
    group_by(Light_num, hours) %>% mutate(Light_only = True_val[Diuron_num == 0]) %>% # Handy trick similar to pivot wider
    group_by(Diuron_num, hours) %>% mutate(Diuron_only = True_val[Light_num == 80]) %>%
    group_by(hours) %>% 
    mutate(Control = True_val[Light_num == 80 & Diuron_num == 0]) %>%
    ungroup() %>% filter(Treatment == "Light_diuron") %>% within({
      
      metric <-  log((True_val/Control)/(Diuron_only*Light_only/Control^2)) # Compute the statistic
      
    })
  
  # Get the response based on the endpoint
  resp_var <- switch(model, "growth" = expr(celld), "photosynthesis" = expr(Yield))
  
  # Compute the raw statistic based on the response with the random and fixed effects
  eval(expr({naive_statistic <- dat_fake %>% 
    mutate(Treatment = case_when(Light_num == 80 & Diuron_num == 0 ~ "Control",
                                 Light_num == 80  ~ "Diuron_only",
                                 Diuron_num == 0  ~ "Light_only",
                                 TRUE ~ "Light_diuron")) %>% 
    group_by(Light_num, hours, block) %>% mutate(Light_only = `[`(!!resp_var, Diuron_num == 0)) %>% 
    group_by(Diuron_num, hours, block) %>% mutate(Diuron_only = `[`(!!resp_var, Light_num == 80)) %>% 
    group_by(hours, block) %>% 
    mutate(Control = `[`(!!resp_var, Light_num == 80 & Diuron_num == 0)) %>% 
    ungroup() %>% filter(Treatment == "Light_diuron") %>% within({
      
      metric <-  log((!!resp_var/Control)/(Diuron_only*Light_only/Control^2))
      
    }) %>% group_by(Light_num, hours, Diuron_num) %>% 
    # Compute the approx 95% CIs
    summarise(naive_metric = mean(metric), SE = sd(metric)/sqrt(length(metric)), 
              naive_lwr = naive_metric - 1.96*SE, naive_upr = naive_metric + 1.96*SE) %>% 
    select(-SE)
  }))
  
  # Compute classification metrics for the modelling approach and the raw statistic 
  test_dat <- left_join(true_statistic, Diuron_CItrend3, by = c("hours", "Light_num", "Diuron_num")) %>% 
    left_join( naive_statistic, by = c("hours", "Light_num", "Diuron_num")) %>% 
    within( { 
      is_captured_model <- metric >= X2.5. & metric <= X97.5.
      correctly_model <- sign(X50.) == sign(metric) & sign(X2.5.) == sign(X97.5.)
      incorrectly_model <- sign(X50.) != sign(metric) & sign(X2.5.) == sign(X97.5.)
      is_captured_naive <- metric >= naive_lwr & metric <= naive_upr
      correctly_naive <- sign(naive_metric) == sign(metric) & sign(naive_lwr) == sign(naive_upr)
      incorrectly_naive <- sign(naive_metric) != sign(metric) & sign(naive_lwr) == sign(naive_upr)
    
    })
  
  # Format results
  summary_data <- test_dat %>% select(is_captured_model, correctly_model, incorrectly_model, is_captured_naive, correctly_naive, incorrectly_naive) %>% 
    colSums(na.rm = TRUE) %>% t() %>% data.frame(., error_str, interaction_str, model)
  
  # Return what the user wants                                               
  if(return_data){
    
    return(list(summary_data, test_dat))
    
  } else {
    
  return(summary_data) 
    
}
  
}


# Photosynthesis model validation -----------------------------------------

# Specify the levels of natural variation, interaction strengths, and number of reps for the photosynthesis endpoint
n_reps <- 200
photo_sims_params <- expand_grid(error_str = c(0.01, 0.025, 0.05, 0.1), 
                                 interaction_str = c(-0.5, -0.25, -0.1, 0.1, 0.25, 0.5), 
                                 i = seq_len(n_reps), model = "photosynthesis")
photo_sims_params$i <- NULL  

# set the seed and run the simulations
set.seed(101)
photo_sims <- future_pmap(photo_sims_params, model_validation_sim) %>% do.call("rbind", .)

# Process the simulations so we have proportions
photo_sims_summary <- photo_sims %>% group_by(error_str, interaction_str, model) %>% 
  summarise(p_captured_model = sum(is_captured_model)/(n_reps*40), p_correctly_model = sum(correctly_model)/(n_reps*40),  
            p_incorrectly_model = sum(incorrectly_model)/(n_reps*40), p_captured_naive = sum(is_captured_naive)/(n_reps*40),  
            p_correctly_naive = sum(correctly_naive)/(n_reps*40),  p_incorrectly_naive = sum(incorrectly_naive)/(n_reps*40)) %>% 
  pivot_longer(cols = p_captured_model:p_incorrectly_naive, names_to = "eval_metric", values_to = "value") %>% 
  mutate(`Interaction Type` = factor(case_when(interaction_str == -0.5 ~ "Strong synergism",
                                               interaction_str == -0.25 ~ "Moderate synergism",
                                               interaction_str == -0.1 ~ "Weak synergism",
                                               interaction_str == 0.1 ~ "Weak antagonism",
                                               interaction_str == 0.25 ~ "Moderate antagonism",
                                               interaction_str == 0.5 ~ "Strong antagonism"), 
                                     levels = c("Weak synergism", "Weak antagonism",
                                                "Moderate synergism", "Moderate antagonism",
                                                "Strong synergism", "Strong antagonism")),
         Method = ifelse(grepl("_model", eval_metric), "Model", "Raw statistic"),
         
  )

## Plot to access ability to correctly classify the interaction
ggplot(filter(photo_sims_summary, eval_metric %in% c("p_correctly_model", "p_correctly_naive")),
       aes( x = factor(round(error_str, 2)), y = value, color = Method)) + ylab("Proportion of interactions correctly classified") + 
  xlab("Additive Error") +
  geom_point(position = position_dodge(0.5)) + facet_wrap(. ~ `Interaction Type`) + theme_bw()

ggsave("Outputs/Model validation/proportion_correctly_classified_photo.png", width = 25, height = 15, dpi = 300, units = "cm")
##

## Plot to access the rate of misclassifications
Error_rate_summary <- photo_sims_summary %>% group_by(error_str, `Interaction Type`, Method) %>% 
  summarise(Error_rate = value[grepl("_incorrectly", eval_metric)]/(value[grepl("_correctly", eval_metric)] + value[grepl("_incorrectly", eval_metric)])) %>% 
  ungroup()

ggplot(Error_rate_summary,
       aes( x = factor(error_str), y = Error_rate, color = Method)) +  ylab("Proportion of classifications incorrect") + 
  xlab("Additive Error") +
  geom_point(position = position_dodge(0.5)) + facet_wrap(. ~ `Interaction Type`) + theme_bw()

ggsave("Outputs/Model validation/error_rate_photo.png", width = 25, height = 15, dpi = 300, units = "cm")
##

## Plot to access the effectiveness of the credible intervals at capturing the true statistic
ggplot(filter(photo_sims_summary, (eval_metric %in% c("p_captured_model"))),
       aes( x = factor(round(error_str, 2)), y = value, color = `Interaction Type`)) + 
  geom_abline(slope = 0, intercept = 0.95, linetype = "dashed") +
  geom_point(position = position_dodge(0.35))  + theme_bw() + ylab("Proportion of True Values Captured by the 95CI") + 
  xlab("Additive Error") + ylim(c(0, 1))

ggsave("Outputs/Model validation/photo_95CI_capture.png", width = 17, height = 12, dpi = 300, units = "cm")
##

# Growth model validation -----------------------------------------

# Specify the levels of natural variation, interaction strengths, and number of reps for the growth endpoint
n_reps <- 200
growth_sims_params <- expand_grid(error_str = c(0.01, 0.025, 0.05, 0.1),
                                  interaction_str = c(-0.8,-0.5, -0.1, 0.1, 0.5, 0.8), 
                                  i = seq_len(n_reps), model = "growth")
growth_sims_params$i <- NULL  


# set the seed and run the simulations
set.seed(102)
growth_sims <- future_pmap(growth_sims_params, model_validation_sim) %>% do.call("rbind", .)

# Process the simulations so we have proportions
growth_sims_summary <- growth_sims %>% group_by(error_str, interaction_str, model) %>% 
  summarise(p_captured_model = sum(is_captured_model)/(n_reps*40), p_correctly_model = sum(correctly_model)/(n_reps*40),  
            p_incorrectly_model = sum(incorrectly_model)/(n_reps*40), p_captured_naive = sum(is_captured_naive)/(n_reps*40),  
            p_correctly_naive = sum(correctly_naive)/(n_reps*40),  p_incorrectly_naive = sum(incorrectly_naive)/(n_reps*40)) %>% 
  pivot_longer(cols = p_captured_model:p_incorrectly_naive, names_to = "eval_metric", values_to = "value") %>% 
  mutate(`Interaction Type` = factor(case_when(interaction_str == -0.8 ~ "Strong synergism",
                                        interaction_str == -0.5 ~ "Moderate synergism",
                                        interaction_str == -0.1 ~ "Weak synergism",
                                        interaction_str == 0.1 ~ "Weak antagonism",
                                        interaction_str == 0.5 ~ "Moderate antagonism",
                                        interaction_str == 0.8 ~ "Strong antagonism"), 
                                     levels = c("Weak synergism", "Weak antagonism",
                                                "Moderate synergism", "Moderate antagonism",
                                                "Strong synergism", "Strong antagonism")),
         Method = ifelse(grepl("_model", eval_metric), "Model", "Raw statistic"),
         
         )

## Plot to access ability to correctly classify the interaction
ggplot(filter(growth_sims_summary, eval_metric %in% c("p_correctly_model", "p_correctly_naive")),
       aes( x = factor(error_str), y = value, color = Method)) + ylab("Proportion of interactions correctly classified") + 
  xlab("Multiplicative Error") +
  geom_point(position = position_dodge(0.5)) + facet_wrap(. ~ `Interaction Type`) + theme_bw()

ggsave("Outputs/Model validation/proportion_correctly_classified_growth.png", width = 25, height = 15, dpi = 300, units = "cm")
##

## Plot to access the rate of misclassifications
Error_rate_summary <- growth_sims_summary %>% group_by(error_str, `Interaction Type`, Method) %>% 
  summarise(Error_rate = value[grepl("_incorrectly", eval_metric)]/(value[grepl("_correctly", eval_metric)] + value[grepl("_incorrectly", eval_metric)])) %>% 
  ungroup()
  
ggplot(Error_rate_summary,
       aes( x = factor(error_str), y = Error_rate, color = Method)) +  ylab("Proportion of classifications incorrect") + 
  xlab("Multiplicative Error") +
  geom_point(position = position_dodge(0.5)) + facet_wrap(. ~ `Interaction Type`) + theme_bw()

ggsave("Outputs/Model validation/error_rate_growth.png", width = 25, height = 15, dpi = 300, units = "cm")
##

## Plot to access the effectiveness of the credible intervals at capturing the true statistic
ggplot(filter(growth_sims_summary, (eval_metric %in% c("p_captured_model"))),
       aes( x = factor(error_str), y = value, color = `Interaction Type`)) + 
  geom_abline(slope = 0, intercept = 0.95, linetype = "dashed") +
  geom_point(position = position_dodge(0.35))  + theme_bw() + ylab("Proportion of True Values Captured by the 95CI") + 
  xlab("Multiplicative Error") + ylim(c(0, 1))

ggsave("Outputs/Model validation/growth_95CI_capture.png", width = 17, height = 12, dpi = 300, units = "cm")
##


# Test a single model simulation ------------------------------------------

# COMMENTED: this code allows us to visually investigate what is going on in a single model simulation

# check <- model_validation_sim(model = "growth", error_str = 0.1, interaction_str = 0.5, return_data = TRUE))
# 
# check[[1]]
# 
# ggplot(check[[2]]) + aes(x = factor(hours), y = X50., colour = paste(Light_num, Diuron_num), fill = paste(Light_num, Diuron_num)) +
#   #geom_line(aes(x = as.integer(hours_fact))) +
#   geom_hline(yintercept = 0, lty= "dashed") + theme_classic() + xlab("Hours") +
#   geom_point(position = position_dodge2(0.5)) +
#   # geom_point(aes(x = hours_fact, y = X97.5., colour = Treatment, alpha = (X97.5. == 5.5)),
#   #            position = position_dodge2(0.5), shape = 2) +
#   # geom_point(aes(x = hours_fact, y = X2.5., colour = Treatment, alpha = (X2.5. == -5.5)),
#   #            position = position_dodge2(0.5), shape = 6) +
#   geom_linerange(mapping = aes(ymin = X2.5., ymax = X97.5.), position = position_dodge2(0.5)) +
#   geom_point(aes(y = naive_metric), position = position_dodge2(0.57), colour = "darkgrey") +
#   geom_linerange(mapping = aes(ymin = naive_lwr, ymax = naive_upr), colour = "darkgrey", position = position_dodge2(0.57)) +
#   geom_point(aes(y = metric), position = position_dodge2(0.5), colour = "black") +
#   ylab(expression("I"[R])) + guides(alpha = "none")
