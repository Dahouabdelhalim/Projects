######### Code for richnes X fire severity X productivity model for Ecology Report Brodie et al., 2021
######### Brodie, E. G., J. E. D. Miller, and H. D. Safford. 2021. Productivity modifies the effects of fire severity on understory diversity. Ecology.





# Load data and libraries
library(tidyverse)
library(brms)
library(tidybayes)
d <- read_csv("manuscript_data_code/sev_div_prod_data.csv") # will need to modify to direct to file on your computer

# scale and center continuous predictors for analysis
d$hli_st <- scale(d$hli)
d$tsf_st <- scale(d$tsf)
d$p_st <- scale(d$average_ndvi)
d$sev_int <- as.integer(d$fire_sev) # ordered categorical variables should be fed in as integers or ordered factors


##### RUN MODEL ---------------------------------------------------
# intercept-only model first
fit0 <- brm(formula = sp_rich ~ (1|fire_id), 
            data = d, family = negbinomial())

# full model
# model takes about 3 minutes to run so saving it as a file on your computer with the 'file' argument is recommended
fit1 <- brm(formula = sp_rich ~ mo(sev_int)*p_st + tsf_st + hli_st + (1|fire_id), 
                     data = d, family = negbinomial(), 
                     cores = 4, chains = 4, file = "manuscript_data_code/fit_sevdivprod.rds") # you will need to change the file argument here to save model to your computer


##### CHECK MODEL ---------------------------------------------------
plot(fit1) # check for mixing of chains
# visual posterior predictive check for model fit
pp_check(fit1) # all data together
pp_check(fit1,type = "intervals_grouped", group = "sev_int")#grouped by severity

# model summary and marginal effects
summary(fit1, prob = .95)
conditional_effects(fit1, method = "fitted")
conditional_effects(fit1, effects = "sev_fac:p_st", method = "fitted")


# compare fit of intercept-only and full model
loo(fit0, fit1)
# Model comparisons:
#  elpd_diff se_diff
# fit1   0.0       0.0  
# fit0 -55.7       8.7 
# full model posterior has a significantly higher estimated pointwise predictive density than intercept-only


##### VISUALIZE MODEL ---------------------------------------------------
# get values for productivity quantiles to produce marginal effects for
quantile(d$p_st, c(.15,.5,.85))
#  15%         50%         85% 
# -0.98571038 -0.05613479  1.09723592 

# draw and summarize posterior samples for all combinations of fire severity and our three prodcutivity quantiles with other predictors at their median values (i.e. interaction effect of fire severity and productivity conditional on all other predictors at median, "conditional effects")
post_pred1 <- d %>%
  modelr::data_grid(sev_int = unique(sev_int), tsf_st = median(d$tsf_st), hli_st = median(d$hli_st), p_st = c(-0.98571038,-0.05613479, 1.09723592)) %>% 
  add_fitted_draws(fit1,re_formula = ~0, scale = 'response', n = 1000) %>% 
  mean_hdci(.value)

# plot marginal effects of interaction between fire severity and productivity
post_pred1$p_st_fac <- as.factor(post_pred1$p_st)
ggplot( data = post_pred1, aes(x = sev_int, y = .value, group = p_st_fac, color = p_st_fac)) +
  geom_pointinterval(aes(ymin = .lower, ymax = .upper), point_size = .6, size = .4, position = position_dodge(width = .42)) +
  xlab("Fire severity class") +
  ylab("Species richness") +
  ylim(3,28) +
  scale_x_discrete(labels = c("Unburned", "Low", "Low-moderate", "Moderate", "High-moderate", "High")) +
  scale_color_manual(values = c("#CBB9BD","#957D37", "#046211"),labels = c( "Low",  "Medium", "High"), name = "Productivity") +
  theme_bw() +
  theme(text = element_text(size = 7),
        axis.text.x = element_text(angle = 30, vjust = .8, size = 6),
        axis.text.y = element_text(size = 6),
        legend.text = element_text(size = 6))
  



##### MODEL CONTRASTS ---------------------------------------------------
# draw and summarize posterior samples for all fire severity classes
post_pred2 <- d %>%
  modelr::data_grid(sev_int, tsf_st = median(d$tsf_st), p_st = median(d$p_st), hli_st = median(d$hli_st)) %>% 
  add_fitted_draws(fit1,re_formula = ~0, scale = 'response') 

# PAIRWISE CONTRASTS FOR FIRE SEVERITY CLASSES
# create pairs
for_pairs <- c("0", "1", "2", "3", "4", "5")
pairs <- combn(for_pairs,2)
print(pairs) # pairwise combinations of classes


bonf_val <- (1-.05/15) #  correction for 15 tests
contrasts <-  list()

for(i in 1:15) {
  spdiff <-  (post_pred2[post_pred2$sev_int == pairs[1,i],] - post_pred2[post_pred2$sev_int == pairs[2,i],])
  contrasts[[i]] <- hdci(spdiff$.value, bonf_val)
}

#put contrasts into a readable dataframe
contrasts_df <- data.frame(matrix(unlist(contrasts), ncol = max(lengths(contrasts)), byrow = TRUE)) %>% 
  cbind(t(pairs))
contrasts_df <- contrasts_df %>%
  mutate(sig = sign(contrasts_df$X1) == sign(contrasts_df$X2)) %>% 
  rename("sev1" = "1", "sev2" = "2", "lower_ci" = X1, "upper_ci" = X2) %>% 
  dplyr::select("sev1", "sev2", "lower_ci", "upper_ci", "sig")

head(contrasts_df)







