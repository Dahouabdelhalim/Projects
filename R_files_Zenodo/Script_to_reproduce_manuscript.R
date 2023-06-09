# Hi! This script will allow you to replicate the stats in Levy et al 2000,
  # 'Higher dominance rank is associated with lower glucocorticoids in wild female baboons: A rank metric comparison' in Hormones & Behavior
# If you have any questions, feel free to contact me at emily.j.levy@duke.edu or ejlevy91@gmail.com
# I ran two informal re-analyses of other papers - Beehner et al 2006 and Gesquiere et al 2011 - that I mention in the discussion. 
  # Because these were on previously published datasets, I did not include the code or datasets to re-run those stats.
# This script requires two datasets:
  # 1. Dryad_dataset, which will run all the analyses except for Figure 2
  # 2. Rank coefficients for fig 2, which will allow you to recreate Figure 2
# Enjoy!



# Load packages -----------------------------------------------------------

library(readxl) # to load dataset
library(faraway) # for VIF score
library(glmmTMB) # for models
library(performance) # for R2 calculations
library(ggplot2) # plots
library(forcats) # For Fig 2 (reorders variables)
library(dplyr) # for data cleaning
library(broom) # data cleaning in power analysis
library(writexl) # save power analysis datasets



# Load and prep dataset ---------------------------------------------------

df <- read_xlsx('Dryad_dataset.xlsx', sheet = 1) # Might be re-named, insert name here!
str(df)

# Convert variables to correct class.
  # id, social_group, and hydroyear should be factors (as should all variables that are obviously characters, such as season)
  # All others should be numeric
df$id <- as.factor(df$id)
df$social_group <- as.factor(df$social_group)
df$hydroyear <- as.factor(df$hydroyear)
df$season <- as.factor(df$season)

# Re-level season so wet is base
df$season <- relevel(df$season, ref = "Wet")

# Ready to go!



# VIF scores --------------------------------------------------------------

# Ordinal rank and group size are necessarily correlated. We z-scored group size to avoid variance inflation.

# Not transformed
vif(lm(fgc_concentration_log10 ~ adult_female_groupsize + I(adult_female_groupsize^2) + age_years + reproductive_status + season + years_collected_to_meth + years_meth_to_assay, data = df)) 
vif(lm(fgc_concentration_log10 ~ ordinal_rank + adult_female_groupsize + I(adult_female_groupsize^2) + age_years + reproductive_status + season + years_collected_to_meth + years_meth_to_assay, data = df)) 
vif(lm(fgc_concentration_log10 ~ proportional_rank + adult_female_groupsize + I(adult_female_groupsize^2) + age_years + reproductive_status + season + years_collected_to_meth + years_meth_to_assay, data = df)) 
vif(lm(fgc_concentration_log10 ~ hi_mid_low_rank + adult_female_groupsize + I(adult_female_groupsize^2) + age_years + reproductive_status + season + years_collected_to_meth + years_meth_to_assay, data = df)) 
vif(lm(fgc_concentration_log10 ~ alpha_or_not + adult_female_groupsize + I(adult_female_groupsize^2) + age_years + reproductive_status + season + years_collected_to_meth + years_meth_to_assay, data = df)) 
vif(lm(fgc_concentration_log10 ~ elo_standardized + adult_female_groupsize + I(adult_female_groupsize^2) + age_years + reproductive_status + season + years_collected_to_meth + years_meth_to_assay, data = df)) 
vif(lm(fgc_concentration_log10 ~ ordinal_rank+alpha_or_not + adult_female_groupsize + I(adult_female_groupsize^2) + age_years + reproductive_status + season + years_collected_to_meth + years_meth_to_assay, data = df)) 
vif(lm(fgc_concentration_log10 ~ proportional_rank+alpha_or_not + adult_female_groupsize + I(adult_female_groupsize^2) + age_years + reproductive_status + season + years_collected_to_meth + years_meth_to_assay, data = df)) 

# Transformed
vif(lm(fgc_concentration_log10 ~ adult_female_groupsize_standardized + I(adult_female_groupsize_standardized^2) + age_years + reproductive_status + season + years_collected_to_meth + years_meth_to_assay, data = df)) 
vif(lm(fgc_concentration_log10 ~ ordinal_rank + adult_female_groupsize_standardized + I(adult_female_groupsize_standardized^2) + age_years + reproductive_status + season + years_collected_to_meth + years_meth_to_assay, data = df)) 
vif(lm(fgc_concentration_log10 ~ proportional_rank + adult_female_groupsize_standardized + I(adult_female_groupsize_standardized^2) + age_years + reproductive_status + season + years_collected_to_meth + years_meth_to_assay, data = df)) 
vif(lm(fgc_concentration_log10 ~ hi_mid_low_rank + adult_female_groupsize_standardized + I(adult_female_groupsize_standardized^2) + age_years + reproductive_status + season + years_collected_to_meth + years_meth_to_assay, data = df)) 
vif(lm(fgc_concentration_log10 ~ alpha_or_not + adult_female_groupsize_standardized + I(adult_female_groupsize_standardized^2) + age_years + reproductive_status + season + years_collected_to_meth + years_meth_to_assay, data = df)) 
vif(lm(fgc_concentration_log10 ~ elo_standardized + adult_female_groupsize_standardized + I(adult_female_groupsize_standardized^2) + age_years + reproductive_status + season + years_collected_to_meth + years_meth_to_assay, data = df)) 
vif(lm(fgc_concentration_log10 ~ ordinal_rank+alpha_or_not + adult_female_groupsize_standardized + I(adult_female_groupsize_standardized^2) + age_years + reproductive_status + season + years_collected_to_meth + years_meth_to_assay, data = df)) 
vif(lm(fgc_concentration_log10 ~ proportional_rank+alpha_or_not + adult_female_groupsize_standardized + I(adult_female_groupsize_standardized^2) + age_years + reproductive_status + season + years_collected_to_meth + years_meth_to_assay, data = df)) 



# Models ------------------------------------------------------------------

# These are the 8 main models

null <- glmmTMB(fgc_concentration_log10 ~ age_years + reproductive_status + season + 
                  adult_female_groupsize_standardized + I(adult_female_groupsize_standardized^2) + 
                  years_collected_to_meth + years_meth_to_assay +
                  (age_years|id) + (1|social_group) + (1|hydroyear), data = df, family = "gaussian")
proprank <- glmmTMB(fgc_concentration_log10 ~ proportional_rank + age_years + reproductive_status + season + 
                      adult_female_groupsize_standardized + I(adult_female_groupsize_standardized^2) + 
                      years_collected_to_meth + years_meth_to_assay +
                      (age_years|id) + (1|social_group) + (1|hydroyear), data = df, family = "gaussian")
alpha.proprank <- glmmTMB(fgc_concentration_log10 ~ alpha_or_not + proportional_rank + age_years + reproductive_status + season + 
                            adult_female_groupsize_standardized + I(adult_female_groupsize_standardized^2) + 
                            years_collected_to_meth + years_meth_to_assay + 
                            (age_years|id) + (1|social_group) + (1|hydroyear), data = df, family = "gaussian")
ordrank <- glmmTMB(fgc_concentration_log10 ~ ordinal_rank + age_years + reproductive_status + season + 
                     adult_female_groupsize_standardized + I(adult_female_groupsize_standardized^2) + 
                     years_collected_to_meth + years_meth_to_assay +
                     (age_years|id) + (1|social_group) + (1|hydroyear), data = df, family = "gaussian")
alpha.ordrank <- glmmTMB(fgc_concentration_log10 ~ alpha_or_not + ordinal_rank + age_years + reproductive_status + season + 
                           adult_female_groupsize_standardized + I(adult_female_groupsize_standardized^2) + 
                           years_collected_to_meth + years_meth_to_assay +
                           (age_years|id) + (1|social_group) + (1|hydroyear), data = df, family = "gaussian")
elo.score <- glmmTMB(fgc_concentration_log10 ~ elo_standardized + age_years + reproductive_status + season + 
                       adult_female_groupsize_standardized + I(adult_female_groupsize_standardized^2) + 
                       years_collected_to_meth + years_meth_to_assay +
                       (age_years|id) + (1|social_group) + (1|hydroyear), data = df, family = "gaussian")
alpha <- glmmTMB(fgc_concentration_log10 ~ alpha_or_not + age_years + reproductive_status + season + 
                   adult_female_groupsize_standardized + I(adult_female_groupsize_standardized^2) + 
                   years_collected_to_meth + years_meth_to_assay +
                   (age_years|id) + (1|social_group) + (1|hydroyear), data = df, family = "gaussian")
himidlow <- glmmTMB(fgc_concentration_log10 ~ hi_mid_low_rank + age_years + reproductive_status + season + 
                      adult_female_groupsize_standardized + I(adult_female_groupsize_standardized^2) + 
                      years_collected_to_meth + years_meth_to_assay +
                      (age_years|id) + (1|social_group) + (1|hydroyear), data = df, family = "gaussian")

# Model coefficients and standard errors
  # Coefficients were converted to % change in fGC for every 1 unit change in the predictor by antilogging, subtracting by 1, and then multiplying by 100
  # Standard errors were used to calculate 95% CI by calculating [coefficient +/- 1.96*SE]
summary(alpha)
summary(alpha.proprank)
summary(alpha.ordrank)
summary(proprank)
summary(null)
summary(elo.score)
summary(ordrank)
summary(himidlow)

# Model AIC scores
AIC(alpha)
AIC(alpha.proprank)
AIC(alpha.ordrank)
AIC(proprank)
AIC(null)
AIC(elo.score)
AIC(ordrank)
AIC(himidlow)



# Normality and homoscedasticity ------------------------------------------

# Function
assumptions <- function(variable){
  p1 <- smoothScatter(resid(variable) ~ fitted(variable), xlab = "Fitted Values", ylab = "Residuals")
  abline(h=0, col = "blue")
  p2 <- qqnorm(resid(variable))
  qqline(resid(variable))
}

# Insert model names to run function for each model
assumptions(alpha)
assumptions(alpha.proprank)
assumptions(alpha.ordrank)
assumptions(proprank)
assumptions(null)
assumptions(elo.score)
assumptions(ordrank)
assumptions(himidlow)



# Models with hybrid score ------------------------------------------------

# This is in the supplement
hybrid <- subset(df, !is.na(hybridscore))

length(unique(hybrid$id)) # 192 females
length(unique(hybrid$social_group)) # 15 groups

# Models
h.null <- glmmTMB(fgc_concentration_log10 ~ hybridscore + age_years + reproductive_status + season + adult_female_groupsize_standardized + I(adult_female_groupsize_standardized^2) + years_collected_to_meth + years_meth_to_assay +
                    (age_years|id) + (1|social_group) + (1|hydroyear), data = hybrid, family = "gaussian")
h.proprank <- glmmTMB(fgc_concentration_log10 ~ hybridscore + proportional_rank + age_years + reproductive_status + season + adult_female_groupsize_standardized + I(adult_female_groupsize_standardized^2) + years_collected_to_meth + years_meth_to_assay +
                        (age_years|id) + (1|social_group) + (1|hydroyear), data = hybrid, family = "gaussian")
h.alpha.proprank <- glmmTMB(fgc_concentration_log10 ~ hybridscore + alpha_or_not + proportional_rank + age_years + reproductive_status + season + adult_female_groupsize_standardized + I(adult_female_groupsize_standardized^2) + years_collected_to_meth + years_meth_to_assay +
                              (age_years|id) + (1|social_group) + (1|hydroyear), data = hybrid, family = "gaussian")
h.ordrank <- glmmTMB(fgc_concentration_log10 ~ hybridscore + ordinal_rank + age_years + reproductive_status + season + adult_female_groupsize_standardized + I(adult_female_groupsize_standardized^2) + years_collected_to_meth + years_meth_to_assay +
                       (age_years|id) + (1|social_group) + (1|hydroyear), data = hybrid, family = "gaussian")
h.alpha.ordrank <- glmmTMB(fgc_concentration_log10 ~ hybridscore + alpha_or_not + ordinal_rank + age_years + reproductive_status + season + adult_female_groupsize_standardized + I(adult_female_groupsize_standardized^2) + years_collected_to_meth + years_meth_to_assay +
                             (age_years|id) + (1|social_group) + (1|hydroyear), data = hybrid, family = "gaussian")
h.elo.score <- glmmTMB(fgc_concentration_log10 ~ hybridscore + elo_standardized + age_years + reproductive_status + season + adult_female_groupsize_standardized + I(adult_female_groupsize_standardized^2) + years_collected_to_meth + years_meth_to_assay +
                         (age_years|id) + (1|social_group) + (1|hydroyear), data = hybrid, family = "gaussian")
h.alpha <- glmmTMB(fgc_concentration_log10 ~ hybridscore + alpha_or_not + age_years + reproductive_status + season + adult_female_groupsize_standardized + I(adult_female_groupsize_standardized^2) + years_collected_to_meth + years_meth_to_assay +
                     (age_years|id) + (1|social_group) + (1|hydroyear), data = hybrid, family = "gaussian")
h.null.alpha <- glmmTMB(fgc_concentration_log10 ~ alpha_or_not + age_years + reproductive_status + season + adult_female_groupsize_standardized + I(adult_female_groupsize_standardized^2) + years_collected_to_meth + years_meth_to_assay +
                          (age_years|id) + (1|social_group) + (1|hydroyear), data = hybrid, family = "gaussian")
h.himidlow <- glmmTMB(fgc_concentration_log10 ~ hybridscore + hi_mid_low_rank + age_years + reproductive_status + season + adult_female_groupsize_standardized + I(adult_female_groupsize_standardized^2) + years_collected_to_meth + years_meth_to_assay +
                        (age_years|id) + (1|social_group) + (1|hydroyear), data = hybrid, family = "gaussian")

# Model coefficients
summary(h.alpha)
summary(h.alpha.proprank)
summary(h.alpha.ordrank)
summary(h.proprank)
summary(h.elo.score)
summary(h.null)
summary(h.ordrank)
summary(h.himidlow)

# deltaAIC of the alpha model with hybrid score vs the alpha model without hybrid score 
  # ie, does hybrid score improve model fit? Nope.
AIC(h.null.alpha)-AIC(h.alpha)



# R^2 values --------------------------------------------------------------

# We performed two different R^2 calculations: 
  # 1. For each main model
  # 2. In the alpha model, but removing one fixed effect at a time

# 1.
  # We already ran the models above, just need to calculate R^2
r2(alpha)
r2(alpha.proprank)
r2(alpha.ordrank)
r2(proprank)
r2(elo.score)
r2(null)
r2(ordrank)
r2(himidlow)

# 2.
  # Need to run the models first
no.rank <- glmmTMB(fgc_concentration_log10 ~ age_years + reproductive_status + season + adult_female_groupsize_standardized + I(adult_female_groupsize_standardized^2) + years_collected_to_meth + years_meth_to_assay + 
                     (age_years|id) + (1|social_group) + (1|hydroyear), data = df, family = "gaussian")
no.age_years <- glmmTMB(fgc_concentration_log10 ~ alpha_or_not + reproductive_status + season + adult_female_groupsize_standardized + I(adult_female_groupsize_standardized^2) + years_collected_to_meth + years_meth_to_assay + 
                    (1|id) + (1|social_group) + (1|hydroyear), data = df, family = "gaussian")
no.repro <- glmmTMB(fgc_concentration_log10 ~ alpha_or_not + age_years + season + adult_female_groupsize_standardized + I(adult_female_groupsize_standardized^2) + years_collected_to_meth + years_meth_to_assay + 
                      (age_years|id) + (1|social_group) + (1|hydroyear), data = df, family = "gaussian")
no.season <- glmmTMB(fgc_concentration_log10 ~ alpha_or_not + age_years + reproductive_status + adult_female_groupsize_standardized + I(adult_female_groupsize_standardized^2) + years_collected_to_meth + years_meth_to_assay + 
                       (age_years|id) + (1|social_group) + (1|hydroyear), data = df, family = "gaussian")
no.grpsize <- glmmTMB(fgc_concentration_log10 ~ alpha_or_not + age_years + reproductive_status + season + years_collected_to_meth + years_meth_to_assay + 
                        (age_years|id) + (1|social_group) + (1|hydroyear), data = df, family = "gaussian")
no.years_collected_to_meth <- glmmTMB(fgc_concentration_log10 ~ alpha_or_not + age_years + reproductive_status + season + adult_female_groupsize_standardized + I(adult_female_groupsize_standardized^2) + years_meth_to_assay + 
                                        (age_years|id) + (1|social_group) + (1|hydroyear), data = df, family = "gaussian")
no.years_meth_to_assay <- glmmTMB(fgc_concentration_log10 ~ alpha_or_not + age_years + reproductive_status + season + adult_female_groupsize_standardized + I(adult_female_groupsize_standardized^2) + years_collected_to_meth +
                                    (age_years|id) + (1|social_group) + (1|hydroyear), data = df, family = "gaussian")

# And now calculate R^2 values
r2(alpha)
r2(no.rank)
r2(no.age_years)
r2(no.repro)
r2(no.season)
r2(no.grpsize)
r2(no.years_collected_to_meth)
r2(no.years_meth_to_assay)



# Fig 2: Forest Plot Comparing Rank Coefficients --------------------------------------

# I took the coefficients and SEs for the rank variables in each main model (see 'Models' above). 
  # I put them in a new Excel spreadsheet to make life a little easier
ranks <- read_xlsx('Rank coefficients for fig 2.xlsx', sheet=1)

ggplot(ranks, aes(x=fct_reorder(Model_Variable, -AIC), y=Converted_Estimate_absolute, color=Model)) + 
  theme(text=element_text(size=12), axis.text=element_text(color='black'), 
        panel.background = element_rect(fill = NA),
        axis.line=element_line()) +
  geom_rect(xmin=0, xmax=4.5, 
            ymin=-4, ymax=14, size=0, fill='grey80', color='grey80') +
  geom_rect(xmin=4.5, xmax=5.5, 
            ymin=-4, ymax=14, size=0, fill='grey90', color='grey90') +
  geom_hline(yintercept=0, linetype=2) +   
  geom_point(size=3) +
  geom_segment(aes(yend=High_95CI_absolute, y=Low_95CI_absolute, xend=Model_Variable), size=1.35) +
  coord_flip() +
  labs(x='Model and Variable', y='Rank Estimate Converted to % Change in fGC') +
  scale_color_brewer(breaks=levels(with(ranks, reorder(Model, AIC))),palette='Dark2') +
  theme(text=element_text(size=12), axis.text=element_text(color='black'), 
        panel.background = element_rect(fill = NA)) 
ggsave('Rank effect sizes.pdf', width = 10, height = 3)
# I then added labels in Illustrator



# Supplemental Figure: Group size histogram ----------------------------------------------------

ggplot(df, aes(x=adult_female_groupsize)) + geom_histogram(bins=32, stat='count') + 
  theme_minimal() + xlim(0,32) + labs(x='Adult Female Group Size', y='Count')
ggsave('Group size histogram.pdf', width = 4, height = 3)



# Removing random effect of ID ----------------------

# This pops up quickly in the discussion when discussing reasons why we might see inconsistent results across studies of female baboons

ordrank.no.id <- glmmTMB(fgc_concentration_log10 ~ ordinal_rank + age_years + reproductive_status + season + adult_female_groupsize_standardized + I(adult_female_groupsize_standardized^2) + 
                           years_collected_to_meth + years_meth_to_assay +
                           (1|social_group) + (1|hydroyear), data = df, family = "gaussian")

summary(ordrank.no.id)

# Coefficient & SE for ordinal rank
or.no.id <- 0.0008495
or.no.id.se <- 0.0002692

# Anti-logging coefficient to get percent change in fGC for every 1-unit change in ordinal rank
10^or.no.id
# Anti-logging and calculate 95% confidence interval
10^(or.no.id - 1.96*or.no.id.se)
10^(or.no.id + 1.96*or.no.id.se) 



# Power Simulation/Calculation --------------------------------------------------------------

# This is mentioned in the discussion, and the figures and methods are described in the supplement
# Ali Galezo was super helpful in writing this!
# Keep in mind that every time you do this, the results will be slightly different

# This function specifices nfemales and nsamples - this is what we're using to simulate the Weingrill and Wittig studies.
simulation <- function(nfemales, nsamples){
  females_sub <- df %>% group_by(id) %>% dplyr::summarise(n=n()) %>% filter(n >= nsamples/nfemales) # Removes errors where it selects a female with a couple of samples and then can't randomly pick enough rows of data.
  females <- sample(unique(females_sub$id), nfemales, replace=FALSE) # Returns nfemales different ids.
  females_data <- df %>% filter(id %in% females) # Dataset with all the rows of data from each of those n females.
  split <- split(females_data, females_data$id) # Turns females_data into a list.
  one_row_list <- lapply(split, function(x){
    x[sample(nrow(x), 1), ]
  })                                #Pick one random row from each of these females. 
  one_row_df <- do.call("rbind", one_row_list) # Binds the list together (should be nfemales rows)
  rows_left <- anti_join(females_data, one_row_df, 
                         by=c("id", "years_collected_to_meth", "years_meth_to_assay", "age_years", "fgc_concentration", "reproductive_status", 
                              "social_group", "ordinal_rank", "proportional_rank", "adult_female_groupsize", 
                              "hydroyear", "season", "fgc_concentration_log10", "adult_female_groupsize_standardized", "hi_mid_low_rank",
                              "alpha_or_not", "elo_standardized", "hybridscore")) # Removes samples that were used to populate before. Warning message is fine.
  remainder_rows <- rows_left[sample(nrow(rows_left), (nsamples-nfemales)), ] # Should be df with rows=nsamples-nfemales
  final <- rbind(remainder_rows, one_row_df)
  model.proportional_rank <- glmmTMB(fgc_concentration_log10 ~ proportional_rank + age_years + reproductive_status + season + adult_female_groupsize_standardized + I(adult_female_groupsize_standardized^2) + years_collected_to_meth + years_meth_to_assay +
                           (1|id) + (1|social_group) + (1|hydroyear), data = final, family='gaussian') # Note that I removed random slope of age due to convergence problems
  tibble.proportional_rank <- broom.mixed::tidy(model.proportional_rank, effects='fixed', data=final)
  pval.proportional_rank <- as.numeric(tibble.proportional_rank[2,7])
  return(pval.proportional_rank)
}

## Sample size from Weingrill et al 2004

# Run the function 10000 times
rep1 <- replicate(10000, simulation(nfemales=10, nsamples=260))
# Transform to data frame
rep1_df <- as.data.frame(rep1)
# Remove NAs due to model convergence problems
rep1_df_no_na <- subset(rep1_df, !is.na(rep1)) 
# Calculate the percent of 'significant' p-values.
p1 <- nrow(subset(rep1_df_no_na, rep1 <= 0.05))/nrow(rep1_df_no_na)*100 
# Plot
ggplot(rep1_df, aes(x=rep1)) + 
  geom_histogram(binwidth=0.01) + 
  geom_vline(aes(xintercept=0.05), color='red') + 
  labs(x='p-value', y='Count', title='Simulating Weingrill et al 2004 Sample Size', subtitle='260 fecal samples from 10 females') +
  theme_minimal()
# Save plot
ggsave('Weingrill Simulation.pdf', width = 5, height = 3)
# Save dataset
write_xlsx(rep1_df, 'WeingrillSimulation.xlsx')

## Sample size from Wittig et al 2008

# Run the function 10000 times
rep2 <- replicate(10000, simulation(nfemales=22, nsamples=532))
# Transform to data frame
rep2_df <- as.data.frame(rep2)
# Remove NAs due to model convergence problems
rep2_df_no_na <- subset(rep2_df, !is.na(rep2)) 
# Calculate the percent of 'significant' p-values.
p2 <- nrow(subset(rep2_df_no_na, rep2 <= 0.05))/nrow(rep2_df_no_na)*100 
# Plot
ggplot(rep2_df, aes(x=rep2)) + 
  geom_histogram(binwidth=0.01) + 
  geom_vline(aes(xintercept=0.05), color='red') + 
  labs(x='p-value', y='Count', title='Simulating Wittig et al 2008 Sample Size', subtitle='532 fecal samples from 22 females') +
  theme_minimal()
# Save plot
ggsave('Wittig Simulation.pdf', width = 5, height = 3)
# Save dataset
write_xlsx(rep2_df, 'WittigSimulation.xlsx')

