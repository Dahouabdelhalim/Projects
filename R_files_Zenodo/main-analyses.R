# This script replicates the primary analyses in the paper, producing Figures 1 and 2.
# It also produces the primary results in table form for the Supplementary Material,
# Heterogeneous treatment effect analyses, and robustness checks.
# The script makes use of the datasets waves-1-and-2.csv, wave-2, and mosul-outcome-wave-2.csv

# Change the working directory as necessary
setwd("~/Dropbox/iraq-experiment/science/replication/data")

library(tidyverse)
library(lfe)
library(broom)
library(stargazer)
library(glue)
library(ggpubr)
library(viridis)
library(pbapply)
library(stargazer)


# Read in data ------------------------------------------------------------

# Waves one and two

merged_pilot <- read.csv("waves-1-and-2.csv", stringsAsFactors=FALSE)

# Wave two 

df <- read.csv("wave-2.csv", stringsAsFactors=FALSE) %>% 
  mutate(block = as.character(block))

# Data on patronizing restaurant in Mosul

mosul_df <- read.csv("mosul-outcome-wave-2.csv", stringsAsFactors=FALSE) %>% 
  mutate(block = as.character(block))


# Estimate average treatment effects --------------------------------------
# The results are used to produce Figure 1 and Table S1 of the paper

outcome_labels <- c(
  "train_t1" = "Train w/ Muslims",
  "voted_muslim_extra" = "Vote Muslim Award",
  "own_group_preference" = "Mixed Team Sign-Up",
  "went_mosul" = "Visit Mosul",
  "attend" = "Attend Mixed Event",
  "donate_t1_new1" = "Donate Mixed NGO",
  "dv_secular_t1" = "National Unity",
  "dv_neigh_t1" = "Muslim Neighbor",
  "dv_blame_t1" = "Muslim Blame"
)

ate <- function(dta, dv, treat = "treated", covars = c(), cluster = "team", extravars_extract = NULL, verbose = TRUE) {
  
  if (!is.numeric(dta[, treat]) | !all(dta[, treat] %in% c(0, 1))) {
    stop("The treatment variable should be numeric and binary")
  }
  if ("block" %in% covars) {
    if (!(is.factor(dta$block) | is.character(dta$block))) {
      stop("The covariate 'block' (the randomization block) should be encoded as string or factor")
    }
  }
  
  # Center covariates to make the intercept interpretable
  # Numeric covariates are centered on their median value
  # For categorical covariates, the reference category is set to the mode
  # The randomization block is treated as a special case given that the mode is not meaningful
  # Note that centering has no effect on the treatment effect estimate
  if (verbose) cat("\\n---- Centering covariates on representative values ----\\n")
  for (covariate_name in covars) {
    covariate <- dta[, covariate_name]
    
    if (is.numeric(covariate)) {
      dta[, covariate_name] <- covariate - median(covariate)
      if (verbose) glue("{covariate_name} centered on its median value {median(covariate)}") %>% print()
      
    } else if (covariate_name == "block") {
      # Find the control group block with the most typical outcome
      dta$outcome_dup <- dta[, dv]
      dta$treat_dup <- dta[, treat]
      typical_block <- dta %>%
        filter(treat_dup == 0) %>%
        group_by(block) %>%
        dplyr::summarize(mean_y = mean(outcome_dup, na.rm = TRUE)) %>%
        ungroup() %>%
        arrange(mean_y) %>%
        slice(ceiling(n() / 2)) %>% # take the median (breaking ties upwards)
        pull(block)
      dta$block <- relevel(as.factor(dta$block), ref = typical_block)
      
      if (verbose) glue("block's reference category set to {typical_block}") %>% print()
      
    } else if (is.factor(covariate) | is.character(covariate)) {
      modal_category <- table(covariate) %>% sort(decreasing = TRUE) %>% names() %>% .[1]
      dta[, covariate_name] <- relevel(as.factor(covariate), ref = modal_category)
      if (verbose) glue("{covariate_name}'s reference category set to {modal_category}") %>% print()
      
    } else {
      stop("Cannot handle {covariate_name} of class {class(covariate)}")
    }
  }
  
  # Create OLS formula
  rhs <- paste(c(treat, covars), collapse = " + ")
  f <- paste(dv, "~", rhs, "| 0 | 0 |", cluster)
  
  # Estimate regression
  m <- felm(formula(f), data = dta)
  if (verbose) {
    cat("\\n---- Regression summary ----\\n")
    glue("\\nFormula used: {f}") %>% print()
    summary(m) %>% print()
  }
  
  # Extract results of interest
  model_results <- tidy(m) %>%
    mutate(term = case_when(
      term == "(Intercept)" ~ "control",
      term == treat ~ "treatment_effect",
      TRUE ~ term)
    ) %>%
    filter(term %in% c("control", "treatment_effect", extravars_extract)) %>%
    select(term, estimate, stderr = std.error)
  
  outcome_treated <- data.frame(
    term = "treated",
    estimate =
      model_results %>% filter(term == "control") %>% pull(estimate) +
      model_results %>% filter(term == "treatment_effect") %>% pull(estimate),
    stderr = sqrt(
      vcov(m) %>% .["(Intercept)", "(Intercept)"] +
        vcov(m) %>% .[treat, treat] +
        vcov(m) %>% .["(Intercept)", treat] * 2
    )
  )
  
  model_results <- rbind(model_results, outcome_treated) %>%
    mutate(p_value = 2 * pnorm(-abs(estimate / stderr))) %>%
    mutate(ci_lb = estimate - stderr * qnorm(0.975)) %>%
    mutate(ci_ub = estimate + stderr * qnorm(0.975)) %>%
    mutate(outcome = dv) %>%
    arrange(term)
  
  if (exists("outcome_labels")) {
    model_results <- model_results %>%
      mutate(outcome = plyr::mapvalues(
        outcome,
        from = names(outcome_labels),
        to = outcome_labels,
        warn_missing = FALSE)
      )
  }
  
  return(model_results)
}


# Estimate ATEs, including pre-treatment covariates to increase precision
covariates <- c(
  "edu" , "church" , "income" , "status1" , "marital" , "isis_abuse",
  "player_type", "block", "birth.year"
)

covariates_merged <- c(
  "edu" , "church" , "income" , "status1" , "marital" , "isis_abuse",
  "player_type", "birth.year"
)

## Behavioral outcomes: On the field

# Registering oneself for a mixed team next season
m1 <- ate(df, dv = "own_group_preference", covars = covariates)

# Voting for a Muslim player to receive a sportsmanship prize
m2 <- ate(df, dv = "voted_muslim_extra", covars = covariates)

# Training with Muslims at the 6 month mark
m3 <- ate(merged_pilot, dv = "train_t1", covars = c(covariates_merged, "train_t0"))

## Behavioral outcomes: Off the field

# Attend mixed social event
m4 <- ate(merged_pilot, dv = "attend", covars = covariates_merged)

# Visiting Mosul
m5 <- ate(mosul_df, dv = "went_mosul", covars = c("block"))


# Donating to mixed NGO
m6 <- ate(df, dv = "donate_t1_new1", covars = c(covariates, "donate_t0_new1"))

## Attitudinal indices

# National unity
m7 <- ate(df, dv = "dv_secular_t1", covars = c(covariates, "dv_secular_t0"))

# Acceptance of Muslims as neighbors
m8 <- ate(df, dv = "dv_neigh_t1", covars = c(covariates, "dv_neigh_t0"))

# Blaming Muslims
m9 <- ate(df, dv = "dv_blame_t1", covars = c(covariates, "dv_blame_t0"))


# Figure 1 ----------------------------------------------------------------
custom_theme <- theme(
  plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
  text = element_text(size = 10, colour = "black"),
  axis.text = element_text(size = 10, color = "black"),
  axis.title.y = element_text(size = 11, color = "black", face = "bold"),
  axis.title.x = element_text(size = 9, color = "black", face = "italic", hjust = 1),
  panel.grid.minor.x = element_blank(),
  panel.grid.major.y = element_blank(),
  panel.border = element_blank(),
  axis.ticks.x = element_blank(),
  strip.text.x = element_text(color = "black", size = 9),
  strip.background = element_rect(fill = "gray80"),
  legend.title = element_blank(),
  legend.background = element_blank(),
  legend.position = "bottom"
)
position <- position_dodge(0.6)

plot_data <- bind_rows(
  m1 %>% mutate(type = "On the field outcomes"),
  m2 %>% mutate(type = "On the field outcomes"),
  m3 %>% mutate(type = "On the field outcomes"),
  m4 %>% mutate(type = "Off the field outcomes"),
  m5 %>% mutate(type = "Off the field outcomes"),
  m6 %>% mutate(type = "Off the field outcomes")
) %>%
  mutate(type = factor(type, levels = unique(type)))

plot_outcomes <- plot_data %>%
  filter(term != "treatment_effect") %>%
  mutate(outcome = factor(outcome, levels = rev(outcome_labels))) %>%
  ggplot(aes(x = outcome, y = estimate, color = term)) +
  geom_point(size = 3, position = position) +
  geom_linerange(
    aes(xmin = outcome, xmax = outcome, ymin = 0, ymax = estimate),
    position = position, show.legend = FALSE
  ) +
  coord_flip() +
  facet_wrap(~ type, ncol = 1, scales = "free_y") +
  xlab("") +
  ylab("Percent of players") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.8)) +
  scale_color_viridis(discrete = TRUE, end = 0.6, guide = guide_legend(reverse = TRUE)) +
  theme_light() +
  custom_theme +
  theme(plot.margin = unit(c(1, 0.5, 1, 0.7), "lines")) +
  ggtitle("Means for Treated and Control")

plot_treatment_effects <- plot_data %>%
  filter(term == "treatment_effect") %>%
  mutate(outcome = factor(outcome, levels = rev(outcome_labels))) %>%
  ggplot(aes(x = outcome, y = estimate * 100)) +
  geom_point(size = 3, position = position) +
  geom_errorbar(aes(ymin = ci_lb * 100, ymax = ci_ub * 100), width = 0) +
  geom_hline(yintercept = 0, linetype = 2, color = "gray20") +
  coord_flip() +
  facet_wrap(~ type, ncol = 1, scales = "free_y") +
  xlab("") +
  ylab("Percentage point difference") +
  ylim(-15, 80) +
  scale_color_viridis(discrete = TRUE, end = 0.6) +
  theme_light() +
  custom_theme +
  theme(
    strip.text.x = element_text(color = "white"),
    strip.background = element_rect(fill = "white"),
    axis.text.y = element_blank(),
    plot.margin = unit(c(1, 1, 3.75, 0.7), "lines")
  ) +
  ggtitle("Treatment Effects")

ggarrange(
  plot_outcomes,
  plot_treatment_effects,
  widths = c(1.6, 1)
)



# Testing for differences in on & off field outcomes ----------------------

# Add new variables indicating whether a player:
# - engaged in at least one of the three on-the-field behaviors (and same for off-the-field)
# - engaged in at least two of the three on-the-field behaviors (and same for off-the-field)
# - engaged in all three on-the-field behaviors (and same for off-the-field)

# We can then test whether the on-the-field behaviors were more common than the off-the-field behaviors,
# using each of the three definitions

df_outcome_tests <- df %>%
  filter(!is.na(own_group_preference), !is.na(voted_muslim_extra), !is.na(train_t1)) %>%
  filter(!is.na(attend), !is.na(went_mosul), !is.na(donate_t1_new1)) %>%
  mutate(on_field1 = ifelse(own_group_preference + voted_muslim_extra + train_t1 >= 1, 1, 0)) %>%
  mutate(off_field1 = ifelse(attend + went_mosul + donate_t1_new1 >= 1, 1, 0)) %>%
  mutate(on_field2 = ifelse(own_group_preference + voted_muslim_extra + train_t1 >= 2, 1, 0)) %>%
  mutate(off_field2 = ifelse(attend + went_mosul + donate_t1_new1 >= 2, 1, 0)) %>%
  mutate(on_field3 = ifelse(own_group_preference + voted_muslim_extra + train_t1 == 3, 1, 0)) %>%
  mutate(off_field3 = ifelse(attend + went_mosul + donate_t1_new1 == 3, 1, 0))

m1_on_field <- ate(df_outcome_tests, dv = "on_field1", covars = covariates)
m1_off_field <- ate(df_outcome_tests, dv = "off_field1", covars = covariates)

m2_on_field <- ate(df_outcome_tests, dv = "on_field2", covars = covariates)
m2_off_field <- ate(df_outcome_tests, dv = "off_field2", covars = covariates)

m3_on_field <- ate(df_outcome_tests, dv = "on_field3", covars = covariates)
m3_off_field <- ate(df_outcome_tests, dv = "off_field3", covars = covariates)


# Take a block bootstrap approach to testing the differences between outcome types

block_bootstrap_on_v_off <- function(dta, dv_on_field, dv_off_field, treat = "treated",
                                     covars = c(), cluster = "team", reps = 2000) {
  
  # A function for one block bootstrap draw
  one_draw <- function(dta, dv_on_field, dv_off_field, treat, covars, cluster) {
    
    # Draw new clusters with replacement
    clusters <- unique(dta[, cluster])
    new_clusters <- data.frame(sample(clusters, replace = TRUE))
    names(new_clusters) <- cluster
    
    # Create new data by merging the new clusters with the observed data
    # This will duplicate some clusters while dropping others in the data
    dta_new <- merge(new_clusters, dta, by = cluster)
    
    # Estimate the ATE for the on-the-field outcome
    ate_on_field <- ate(
      dta = dta_new,
      dv = dv_on_field,
      treat = treat,
      covars = covars,
      cluster = "0",
      verbose = FALSE
    ) %>%
      filter(term == "treatment_effect") %>%
      pull(estimate)
    
    # Estimate the ATE for the off-the-field outcome
    ate_off_field <- ate(
      dta = dta_new,
      dv = dv_off_field,
      treat = treat,
      covars = covars,
      cluster = "0",
      verbose = FALSE
    ) %>%
      filter(term == "treatment_effect") %>%
      pull(estimate)
    
    return(data.frame(ate_on_field, ate_off_field))
  }
  
  # Many block bootstrap draws
  bootstrap_results <- pbreplicate(
    n = reps,
    one_draw(dta, dv_on_field, dv_off_field, treat, covars, cluster),
    simplify = FALSE
  ) %>%
    bind_rows()
  
  # The standard error for the difference between the two treatment effects is:
  # sqrt(var(on_field) + var(off_field) - 2cov(on_field, off_field))
  uncertainty_estimates <- bootstrap_results %>%
    summarize(
      stderr_on_field = sd(ate_on_field),
      stderr_off_field = sd(ate_off_field),
      covar = cov(ate_on_field, ate_off_field)
    ) %>%
    mutate(stderr_diff = sqrt(
      stderr_on_field^2 +
        stderr_off_field^2 -
        2 * covar
    ))
  
  # Also get the ATEs estimated from the observed data
  observed_ate_on_field <- ate(
    dta = dta,
    dv = dv_on_field,
    treat = treat,
    covars = covars,
    cluster = "0",
    verbose = FALSE
  ) %>%
    filter(term == "treatment_effect") %>%
    pull(estimate)
  
  observed_ate_off_field <- ate(
    dta = dta,
    dv = dv_off_field,
    treat = treat,
    covars = covars,
    cluster = "0",
    verbose = FALSE
  ) %>%
    filter(term == "treatment_effect") %>%
    pull(estimate)
  
  results <- data.frame(
    observed_ate_on_field,
    observed_ate_off_field,
    diff = observed_ate_on_field - observed_ate_off_field,
    uncertainty_estimates
  ) %>%
    mutate(p_value_diff = 2 * pnorm(-abs(diff / stderr_diff))) %>%
    mutate(ci_lb_diff = diff - stderr_diff * qnorm(0.975)) %>%
    mutate(ci_ub_diff = diff + stderr_diff * qnorm(0.975)) %>%
    mutate(p_value_on_field = 2 * pnorm(-abs(observed_ate_on_field / stderr_on_field))) %>%
    mutate(p_value_off_field = 2 * pnorm(-abs(observed_ate_off_field / stderr_off_field))) %>%
    mutate(dv_on_field = dv_on_field, dv_off_field = dv_off_field)
  
  return(results)
}

set.seed(1928)
m1_on_v_off_bb <- block_bootstrap_on_v_off(df_outcome_tests, "on_field1", "off_field1", covars = covariates)
m2_on_v_off_bb <- block_bootstrap_on_v_off(df_outcome_tests, "on_field2", "off_field2", covars = covariates)
m3_on_v_off_bb <- block_bootstrap_on_v_off(df_outcome_tests, "on_field3", "off_field3", covars = covariates)

# Plots
plot_data <- bind_rows(
  m1_on_field %>% mutate(outcome = "On the field", type = "Percent of players who engaged in 1 or more behaviors"),
  m1_off_field %>% mutate(outcome = "Off the field", type = "Percent of players who engaged in 1 or more behaviors"),
  m2_on_field %>% mutate(outcome = "On the field", type = "... 2 or more behaviors"),
  m2_off_field %>% mutate(outcome = "Off the field", type = "... 2 or more behaviors"),
  m3_on_field %>% mutate(outcome = "On the field", type = "... all 3 behaviors"),
  m3_off_field %>% mutate(outcome = "Off the field", type = "... all 3 behaviors") %>%
  # Fix one negative outcome -- does not change inferences at all:
  mutate(estimate = estimate + abs(estimate[term == "treated"]))
) %>%
  mutate(type = factor(type, levels = unique(type)))

p1 <- plot_data %>%
  filter(term != "treatment_effect") %>%
  ggplot(aes(x = outcome, y = estimate, color = term)) +
  geom_point(size = 3, position = position) +
  geom_linerange(
    aes(xmin = outcome, xmax = outcome, ymin = 0, ymax = estimate),
    position = position, show.legend = FALSE
  ) +
  coord_flip() +
  facet_wrap(~ type, ncol = 1) +
  xlab("") +
  ylab("Percent of players") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  scale_color_viridis(discrete = TRUE, end = 0.6, guide = guide_legend(reverse = TRUE)) +
  theme_light() +
  custom_theme +
  theme(
    legend.title = element_blank(),
    strip.text = element_text(hjust = 0, size = 10),
    panel.grid.major.y = element_blank(),
    plot.margin = unit(c(1, 0.5, 1, 0.7), "lines")
  ) +
  ggtitle("Means for Control and Treated")

p2 <- bind_rows(m1_on_v_off_bb, m2_on_v_off_bb, m3_on_v_off_bb) %>%
  mutate(x = factor(dv_on_field, levels = rev(dv_on_field))) %>%
  ggplot(aes(x = x, y = diff * 100)) +
  geom_point(size = 3, position = position) +
  geom_errorbar(aes(ymin = ci_lb_diff * 100, ymax = ci_ub_diff * 100), width = 0) +
  geom_hline(yintercept = 0, linetype = 2, color = "gray20") +
  coord_flip() +
  facet_wrap(~ dv_on_field, ncol = 1, scales = "free_y") +
  xlab("") +
  ylab("Percentage point difference") +
  scale_y_continuous(limits = c(-5, 85), breaks = seq(0, 80, 20)) +
  scale_color_viridis(discrete = TRUE, end = 0.6) +
  theme_light() +
  custom_theme +
  theme(
    strip.text.x = element_text(color = "white"),
    strip.background = element_rect(fill = "white"),
    axis.text.y = element_blank(),
    plot.margin = unit(c(1, 1, 3.75, 0.7), "lines")
  ) +
  ggtitle("Difference-in-differences")

ggarrange(p1, p2, widths = c(1.75, 1))



# Block bootstrap variance estimates for the ATE --------------------------

block_bootstrap <- function(dta, dv, treat = "treated", covars = c(), cluster = "team", reps = 2000) {
  
  # A function for one block bootstrap draw
  one_block_bootsrap <- function(dta, dv, treat, covars, cluster) {
    
    # Draw new clusters with replacement
    clusters <- unique(dta[, cluster])
    new_clusters <- data.frame(sample(clusters, replace = TRUE))
    names(new_clusters) <- cluster
    
    # Create new data by merging the new clusters with the observed data
    # This will duplicate some clusters while dropping others in the data
    dta_new <- merge(new_clusters, dta, by = cluster)
    
    # Estimate the ATE for the new data
    result <- ate(
      dta = dta_new,
      dv = dv,
      treat = treat,
      covars = covars,
      cluster = "0",
      verbose = FALSE
    ) %>%
      filter(term == "treatment_effect") %>%
      pull(estimate)
    
    return(result)
  }
  
  # Many block bootstrap draws
  bootstrap_results <- pbreplicate(
    n = reps,
    one_block_bootsrap(dta, dv, treat, covars, cluster)
  )
  
  # The key quantity here is the standard deviation of the block bootstrap estimates
  stderr <- sd(bootstrap_results)
  
  # Also get the ATE estimated from the observed data
  observed_estimate <- ate(
    dta = dta,
    dv = dv,
    treat = treat,
    covars = covars,
    cluster = "0",
    verbose = FALSE
  ) %>%
    filter(term == "treatment_effect") %>%
    pull(estimate)
  
  results <- data.frame(
    estimate = observed_estimate,
    stderr = stderr,
    p_value = 2 * pnorm(-abs(observed_estimate / stderr))
  )
  
  return(results)
}

covariates <- c(
  "edu" , "church" , "income" , "status1" , "marital" , "isis_abuse",
  "player_type", "block", "birth.year"
)

covariates_merged <- c(
  "edu" , "church" , "income" , "status1" , "marital" , "isis_abuse",
  "player_type", "birth.year"
)

set.seed(1991)
m1_bb <- block_bootstrap(df, "own_group_preference", covars = covariates)
m2_bb <- block_bootstrap(df, "voted_muslim_extra", covars = covariates)
m3_bb <- block_bootstrap(merged_pilot, "train_t1", covars = c(covariates_merged, "train_t0"))
m4_bb <- block_bootstrap(merged_pilot, "attend", covars = covariates_merged)
m5_bb <- block_bootstrap(mosul_df, dv = "went_mosul", covars = c("block"))
m6_bb <- block_bootstrap(df, "donate_t1_new1", covars = c(covariates, "donate_t0_new1"))
m7_bb <- block_bootstrap(df, "dv_secular_t1", covars = c(covariates, "dv_secular_t0"))
m8_bb <- block_bootstrap(df, "dv_neigh_t1", covars = c(covariates, "dv_neigh_t0"))
m9_bb <- block_bootstrap(df, "dv_blame_t1", covars = c(covariates, "dv_blame_t0"))




# Permutation test for the ATE (Fig. S-10) -------------------------------------------

cluster_permutation <- function(dta, dv, treat = "treated", covars = c(), cluster = "team") {
  
  # Generate random treatment vector
  clusters <- unique(dta[, cluster])
  n_clusters <- length(clusters)
  cluster_assignment <- data.frame(
    cluster = clusters,
    treated_shuffled = rbinom(n_clusters, size = 1, prob = 0.5)
  )
  dta <- merge(dta, cluster_assignment, by.x = cluster, by.y = "cluster")
  
  permutation_result <- ate(
    dta = dta,
    dv = dv,
    treat = "treated_shuffled",
    covars = covars,
    cluster = "0",
    verbose = FALSE
  ) %>%
    filter(term == "treatment_effect") %>%
    select(permutation_estimate = estimate, outcome)
  
  actual_result <- ate(
    dta = dta,
    dv = dv,
    treat = treat,
    covars = covars,
    cluster = cluster,
    verbose = FALSE
  ) %>%
    filter(term == "treatment_effect") %>%
    pull(estimate)
  
  permutation_result$actual_estimate <- actual_result
  
  return(permutation_result)
}

covariates <- c(
  "edu" , "church" , "income" , "status1" , "marital" , "isis_abuse",
  "player_type", "block", "birth.year"
)

covariates_merged <- c(
  "edu" , "church" , "income" , "status1" , "marital" , "isis_abuse",
  "player_type", "birth.year"
)

B <- 2000

m1_perm <- pbreplicate(
  B,
  cluster_permutation(df, dv = "own_group_preference", covars = covariates),
  simplify = FALSE
) %>%
  bind_rows()

m2_perm <- pbreplicate(
  B,
  cluster_permutation(df, dv = "voted_muslim_extra", covars = covariates),
  simplify = FALSE
) %>%
  bind_rows()

m3_perm <- pbreplicate(
  B,
  cluster_permutation(merged_pilot, dv = "train_t1", covars = c(covariates_merged, "train_t0")),
  simplify = FALSE
) %>%
  bind_rows()

m4_perm <- pbreplicate(
  B,
  cluster_permutation(merged_pilot, dv = "attend", covars = covariates_merged),
  simplify = FALSE
) %>%
  bind_rows()

m5_perm <- pbreplicate(
  B,
  cluster_permutation(df, dv = "went_mosul", covars = c("block")),
  simplify = FALSE
) %>%
  bind_rows()

m6_perm <- pbreplicate(
  B,
  cluster_permutation(df, dv = "donate_t1_new1", covars = c(covariates, "donate_t0_new1")),
  simplify = FALSE
) %>%
  bind_rows()

m7_perm <- pbreplicate(
  B,
  cluster_permutation(df, dv = "dv_secular_t1", covars = c(covariates, "dv_secular_t0")),
  simplify = FALSE
) %>%
  bind_rows()

m8_perm <- pbreplicate(
  B,
  cluster_permutation(df, dv = "dv_neigh_t1", covars = c(covariates, "dv_neigh_t0")),
  simplify = FALSE
) %>%
  bind_rows()

m9_perm <- pbreplicate(
  B,
  cluster_permutation(df, dv = "dv_blame_t1", covars = c(covariates, "dv_blame_t0")),
  simplify = FALSE
) %>%
  bind_rows()


perm_results <- rbind(m1_perm, m2_perm, m3_perm, m7_perm)


ggplot(perm_results, aes(x = permutation_estimate)) +
  geom_histogram(binwidth = 0.05, color = "white") +
  facet_wrap(~outcome) +
  geom_vline(data = perm_results %>% distinct(outcome, actual_estimate),
             aes(xintercept = actual_estimate), color = "blue") +
  theme_bw()


# One-sided test
perm_results %>%
  group_by(outcome) %>%
  dplyr::summarize(p_greater_by_random_chance = mean(permutation_estimate > actual_estimate))

# Two-sided test
perm_results %>%
  group_by(outcome) %>%
  dplyr::summarize(p_greater_by_random_chance = mean(abs(permutation_estimate) > abs(actual_estimate)))


# Heterogeneous Treatment Effects  -------------------------------------------
# 1. Team success (Table S-7) -------------------------------------------

# Creating a new binary variable to reflect successful teams (advanced to finals)

success_teams_short <-c("team_1",  "team_4","team_19", "team_13", "m3", "m5")

df$success <- ifelse(df$team %in% success_teams_short, 1, 0)
merged_pilot$success <- ifelse(merged_pilot$team %in% success_teams_short, 1, 0)
mosul_df$success <- ifelse(mosul_df$team %in% success_teams_short, 1, 0)


s1 <- felm(attend ~ treated + birth.year + edu +  church + income + status1 +
             marital + isis_abuse + player_type + success +
             treated:success
           | 0 | 0 | team, data = merged_pilot)


s2 <- felm(train_t1 ~ treated + birth.year+ edu +  church + income + status1 +
             marital + isis_abuse + player_type + success +
             treated:success +train_t0
           | 0 | 0 | team, data = merged_pilot)

s3 <- felm(went_mosul ~ treated +
             block + success + treated:success
           | 0 | 0 | team,
           data = mosul_df)


s4 <- felm(donate_t1_new1 ~ treated + birth.year+ edu +  church + income + status1 +
             marital + isis_abuse + player_type + success +
             treated:success + donate_t0_new1 + block
           | 0 | 0 | team,
           data = df)

s5 <- felm(voted_muslim_extra ~ treated + birth.year + edu +  church + income + status1 +
             marital + isis_abuse + player_type + success + block +
             treated:success
           | 0 | 0 | team,
           data = df)

s6 <- felm(own_group_preference ~ treated + birth.year + edu +  church + income + status1 +
             marital + isis_abuse + player_type + success + block+
             treated:success
           | 0 | 0 | team,
           data = df)


s7 <- felm(dv_secular_t1 ~ treated + birth.year + edu +  church + income + status1 +
             marital + isis_abuse + player_type + success +
             treated:success + block +
             dv_secular_t0 | 0 | 0 | team,
           data = df)

s8 <- felm(dv_neigh_t1 ~ treated + birth.year + edu +  church + income + status1 +
             marital + isis_abuse + player_type + success +
             treated:success + block +
             dv_neigh_t0| 0 | 0 | team,
           data = df)

s9 <- felm(dv_blame_t1 ~ treated + birth.year + edu +  church + income + status1 +
             marital + isis_abuse + player_type + success + block +
             treated:success +
             dv_blame_t0
           | 0 | 0 | team,
           data = df)

stargazer(s1, s2, s3, s4, s5, s6, s7, s8, s9)

# 2. Baseline empathy (Table S-8) -------------------------------------------

# Dichotimizing empathy variable in the merged_pilot dataset
merged_pilot$empathy_sunni_t0_new <- ifelse(merged_pilot$empathy_sunni_t0 %in% c(3, 4), 1,
                                            ifelse(is.na(merged_pilot$empathy_sunni_t0), NA, 0))

s1 <- felm(attend ~ treated + birth.year + edu +  church + income + status1 +
             marital + isis_abuse + player_type + empathy_sunni_t0_new +
             treated:empathy_sunni_t0_new
           | 0 | 0 | team, data = merged_pilot)


s2 <- felm(train_t1 ~ treated + birth.year + edu +  church + income + status1 +
             marital + isis_abuse + player_type + empathy_sunni_t0_new +
             treated:empathy_sunni_t0_new +train_t0
           | 0 | 0 | team, data = merged_pilot)

s3 <- felm(went_mosul ~ treated +
             block + empathy_sunni_t0_new + treated:empathy_sunni_t0_new
           | 0 | 0 | team,
           data = df)


s4 <- felm(donate_t1_new1 ~ treated + birth.year + edu +  church + income + status1 +
             marital + isis_abuse + player_type + empathy_sunni_t0_new +
             treated:empathy_sunni_t0_new + donate_t0_new1 + block
           | 0 | 0 | team,
           data = df)

s5 <- felm(voted_muslim_extra ~ treated + birth.year + edu +  church + income + status1 +
             marital + isis_abuse + player_type + empathy_sunni_t0_new + block +
             treated:empathy_sunni_t0_new
           | 0 | 0 | team,
           data = df)

s6 <- felm(own_group_preference ~ treated + birth.year + edu +  church + income + status1 +
             marital + isis_abuse + player_type + empathy_sunni_t0_new + block+
             treated:empathy_sunni_t0_new
           | 0 | 0 | team,
           data = df)


s7 <- felm(dv_secular_t1 ~ treated + birth.year + edu +  church + income + status1 +
             marital + isis_abuse + player_type + empathy_sunni_t0_new +
             treated:empathy_sunni_t0_new + block +
             dv_secular_t0 | 0 | 0 | team,
           data = df)

s8 <- felm(dv_neigh_t1 ~ treated + birth.year + edu +  church + income + status1 +
             marital + isis_abuse + player_type + empathy_sunni_t0_new +
             treated:empathy_sunni_t0_new + block +
             dv_neigh_t0| 0 | 0 | team,
           data = df)

s9 <- felm(dv_blame_t1 ~ treated + birth.year + edu +  church + income + status1 +
             marital + isis_abuse + player_type + empathy_sunni_t0_new + block +
             treated:empathy_sunni_t0_new +
             dv_blame_t0
           | 0 | 0 | team,
           data = df)

stargazer(s1, s2, s3, s4, s5, s6, s7, s8, s9)

# 3. Baseline contact (Table S-9) -------------------------------------------


s1 <- felm(attend ~ treated + birth.year + edu +  church + income + status1 +
             marital + isis_abuse + player_type + friends_t0_new +
             treated:friends_t0_new
           | 0 | 0 | team, data = merged_pilot)


s2 <- felm(train_t1 ~ treated + birth.year + edu +  church + income + status1 +
             marital + isis_abuse + player_type + friends_t0_new +
             treated:friends_t0_new +train_t0
           | 0 | 0 | team, data = merged_pilot)

s3 <- felm(went_mosul ~ treated +birth.year + edu +  church + income + status1 +
             marital + isis_abuse + player_type + friends_t0_new +
             treated:friends_t0_new + block
           | 0 | 0 | team,
           data = df)


s4 <- felm(donate_t1_new1 ~ treated + birth.year + edu +  church + income + status1 +
             marital + isis_abuse + player_type + friends_t0_new +
             treated:friends_t0_new + donate_t0_new1 + block
           | 0 | 0 | team,
           data = df)

s5 <- felm(voted_muslim_extra ~ treated + birth.year + edu +  church + income + status1 +
             marital + isis_abuse + player_type + friends_t0_new + block +
             treated:friends_t0_new
           | 0 | 0 | team,
           data = df)

s6 <- felm(own_group_preference ~ treated + birth.year + edu +  church + income + status1 +
             marital + isis_abuse + player_type + friends_t0_new + block+
             treated:friends_t0_new
           | 0 | 0 | team,
           data = df)


s7 <- felm(dv_secular_t1 ~ treated + birth.year + edu +  church + income + status1 +
             marital + isis_abuse + player_type + friends_t0_new +
             treated:friends_t0_new + block +
             dv_secular_t0 | 0 | 0 | team,
           data = df)

s8 <- felm(dv_neigh_t1 ~ treated + birth.year + edu +  church + income + status1 +
             marital + isis_abuse + player_type + friends_t0_new +
             treated:friends_t0_new + block +
             dv_neigh_t0| 0 | 0 | team,
           data = df)

s9 <- felm(dv_blame_t1 ~ treated + birth.year + edu +  church + income + status1 +
             marital + isis_abuse + player_type + friends_t0_new + block +
             treated:friends_t0_new +
             dv_blame_t0
           | 0 | 0 | team,
           data = df)

stargazer(s1, s2, s3, s4, s5, s6, s7, s8, s9)

# 4. Violent ISIS Abuse (Table S-10) -------------------------------------------

s1 <- felm(attend ~ treated + birth.year + edu +  church + income + status1 +
             marital +  player_type + isis_abuse_violent +
             treated:isis_abuse_violent
           | 0 | 0 | team, data = merged_pilot)


s2 <- felm(train_t1 ~ treated + birth.year + edu +  church + income + status1 +
             marital +  player_type + isis_abuse_violent +
             treated:isis_abuse_violent +train_t0
           | 0 | 0 | team, data = merged_pilot)

s3 <- felm(went_mosul ~ treated + birth.year + edu +  church + income + status1 +
             marital +  player_type + isis_abuse_violent +
             treated:isis_abuse_violent + block
           | 0 | 0 | team,
           data = df)


s4 <- felm(donate_t1_new1 ~ treated + birth.year+ edu +  church + income + status1 +
             marital +  player_type + isis_abuse_violent +
             treated:isis_abuse_violent + donate_t0_new1 + block
           | 0 | 0 | team,
           data = df)

s5 <- felm(voted_muslim_extra ~ treated + birth.year+ edu +  church + income + status1 +
             marital +  player_type + isis_abuse_violent + block +
             treated:isis_abuse_violent
           | 0 | 0 | team,
           data = df)

s6 <- felm(own_group_preference ~ treated + birth.year+ edu +  church + income + status1 +
             marital +  player_type + isis_abuse_violent + block+
             treated:isis_abuse_violent
           | 0 | 0 | team,
           data = df)


s7 <- felm(dv_secular_t1 ~ treated + birth.year+ edu +  church + income + status1 +
             marital +  player_type + isis_abuse_violent +
             treated:isis_abuse_violent + block +
             dv_secular_t0 | 0 | 0 | team,
           data = df)

s8 <- felm(dv_neigh_t1 ~ treated + birth.year+ edu +  church + income + status1 +
             marital +  player_type + isis_abuse_violent +
             treated:isis_abuse_violent + block +
             dv_neigh_t0| 0 | 0 | team,
           data = df)

s9 <- felm(dv_blame_t1 ~ treated + birth.year+ edu +  church + income + status1 +
             marital +  player_type + isis_abuse_violent + block +
             treated:isis_abuse_violent +
             dv_blame_t0
           | 0 | 0 | team,
           data = df)

stargazer(s1, s2, s3, s4, s5, s6, s7, s8, s9)


# Plotting ATE on individual survey items (Fig. S-6) -------------------------------------------

# Making this variable binary to align with the others
df <- df %>% mutate(empathy_shabak_t1 = case_when(empathy_shabak_t1 > 2 ~ 1, 
                                                  TRUE ~ 0))

m1 <- ate(df, dv = "areas_t1_new", covars = c(covariates, "areas_t0_new"))

m2 <- ate(df, dv = "friends_t1_new", covars = c(covariates, "friends_t0_new"))

m3 <- ate(df, dv = "proud_iraqi_t1_new", covars = c(covariates, "proud_iraqi_t0_new"))

m4 <- ate(df, dv = "secular_t1_new", covars = c(covariates, "secular_t0_new"))

m5 <- ate(df, dv = "secular1_t1_new", covars = c(covariates, "secular1_t0_new"))

m6 <- ate(df, dv = "land_t1_new", covars = c(covariates, "land_t0_new"))

m7 <- ate(df, dv = "empathy_shabak_t1", covars = c(covariates))

m8 <- ate(df, dv = "empathy_sunni_t1", covars = c(covariates, "empathy_sunni_t0"))

m9 <- ate(df, dv = "danger_t1_new", covars = c(covariates, "danger_t0_new"))

m10 <- ate(df, dv = "gold_t1", covars = c(covariates, "gold_t0"))

m11 <- ate(df, dv = "nbr_shiite_shabak_t1", covars = c(covariates, "nbr_shiite_shabak_t0"))

m12 <- ate(df, dv = "suff_sh_t1", covars = c(covariates, "suff_sh_t0"))

m13 <- ate(df, dv = "nbr_sunni_arab_t1", covars = c(covariates, "nbr_sunni_arab_t0"))

m14 <- ate(df, dv = "suff_sa_t1", covars = c(covariates, "suff_sa_t0"))

m15 <- ate(df, dv = "cash_muslim_t1", covars = c(covariates, "cash_muslim_t0"))

m16 <- ate(merged_pilot, dv = "guests_women", covars = c(covariates_merged))


names <- c("comfortable_mixed_areas", "mixed_friends",
           "proud_iraqi", "eth_divides_arbitrary", "natl_identity_first",
           "sell_land_muslim", "empathy_shabak", "empathy_sunni", "iraq_dangerous",
           "gold_adage", "shabak_nbr_ok", "blame_shabak", "sunni_nbr_ok", "blame_sunni",
           "trust_muslim_cash", "brought_guest_event")

means <- c((m1)$estimate[3],
           (m2)$estimate[3],
           (m3)$estimate[3],
           (m4)$estimate[3],
           (m5)$estimate[3],
           (m6)$estimate[3],
           (m7)$estimate[3],
           (m8)$estimate[3],
           (m9)$estimate[3],
           (m10)$estimate[3],
           (m11)$estimate[3],
           (m12)$estimate[3],
           (m13)$estimate[3],
           (m14)$estimate[3],
           (m15)$estimate[3],
           (m16)$estimate[3]
)

ses <- c((m1)$stderr[3],
         (m2)$stderr[3],
         (m3)$stderr[3],
         (m4)$stderr[3],
         (m5)$stderr[3],
         (m6)$stderr[3],
         (m7)$stderr[3],
         (m8)$stderr[3],
         (m9)$stderr[3],
         (m10)$stderr[3],
         (m11)$stderr[3],
         (m12)$stderr[3],
         (m13)$stderr[3],
         (m14)$stderr[3],
         (m15)$stderr[3],
         (m16)$stderr[3]
)

merged <- data.frame(names, means, ses)


merged <- merged %>%  mutate(lb = means - ses * 1.96,
                             ub = means + ses * 1.96)

ggplot(merged, aes(x = reorder(names, - means), y = means)) +
  geom_errorbar(aes(ymin = lb, ymax = ub), width = 0, color = "grey") +
  # facet_wrap(~key) +
  theme_minimal() +
  geom_point(color = "black") +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(text = element_text(size=15))  +
  xlab("") +
  ylab("")




