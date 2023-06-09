######### Code for species trait X fire severity X productivity models for Ecology Report Brodie et al., 2021
######### Brodie, E. G., J. E. D. Miller, and H. D. Safford. 2021. Productivity modifies the effects of fire severity on understory diversity. Ecology.


# Load data and libraries
library(tidyverse)
library(brms)
library(tidybayes)
d <- read_csv("manuscript_data_code/subalpine_traits_data.csv") # will need to modify to direct to file on your computer

# scale and center continuous predictors for analysis
d$tsf_st <- scale(d$tsf)
d$p_st <- scale(d$average_ndvi)
d$sev_int <- as.integer(d$fire_sev)



#################################### **LIFESPAN MODEL** ##########################################
##### RUN MODEL ---------------------------------------------------
# takes about 6 hours to run so would recommend saving the output as shown with the 'file' argument. You will need to change the file in order to save it to your computer
fit_lifespan <- brm(value ~ mo(sev_int)*lifespan*p_st + (lifespan | plot) + (mo(sev_int) | uu), 
             data = d, family = 'bernoulli', iter = 4000, chains = 4, cores = 4, file = "manuscript_data_code/fit_lifespan.rds")


##### CHECK MODEL ---------------------------------------------------
plot(fit_lifespan) # check for mixing of chains
# visual posterior predictive check for model fit
pp_check(fit_lifespan, nsamples = 100) # all data together

# model summary and marginal effects
summary(fit_lifespan)


##### MODEL VISUALIZATION ---------------------------------------------------
# get values for productivity quantiles to produce marginal effects for
quantile(d$p_st, c(.15,.5,.85))
# 15%         50%         85% 
# -1.01221138  0.01229895  1.13537307 

# draw and summarize posterior samples for all combinations of fire severity, lifespan, and our three prodcutivity quantiles
post_pred_lifespan <- d %>%
  modelr::data_grid(sev_int = unique(sev_int),  lifespan = unique(lifespan), p_st = c(-1.01221138, 0.01229895, 1.13537307)) %>% 
  add_fitted_draws(fit_lifespan, re_formula = ~0, scale = 'response') 
unique(post_pred_lifespan$p_st)

#plot predictions in a tryptic by creating three separate plots for low, medium, and high productivity
low <- post_pred_lifespan %>% 
  filter(p_st == -1.01221138) %>% 
  mutate(lifespan = factor(lifespan, levels = c("an_bien", "short_peren", "long_peren"))) %>% 
  ggplot(., aes(y = .value, x = as.factor(sev_int), group = lifespan, color = lifespan)) +
  stat_pointinterval(.width = .95, point_size = 1.7, size = 1.2, position = position_dodge(w=0.35)) +
  labs(x = "", y = "Probability of occurrence", title = "Low productivity", tag = "a") +
  theme_bw() +
  scale_color_manual(name = "Life history", values = c("#2a9d8f","#f4a261", "#e76f51"), guide = guide_legend(override.aes = list(color = "white"))) +
  theme(
    legend.text = element_text(color = "white"),
    legend.title = element_text(color = "white"),
    axis.title.y = element_text(colour = "white")
  ) +
  coord_cartesian(ylim = c(0,0.35)) +
  scale_x_discrete(labels = c("", "", "", "", "", "")) 



med <- post_pred_lifespan %>% 
  filter(p_st == 0.01229895) %>%
  mutate(lifespan = factor(lifespan, levels = c("an_bien", "short_peren", "long_peren"))) %>% 
  ggplot(., aes(y = .value, x = as.factor(sev_int), group = lifespan, color = lifespan)) +
  stat_pointinterval(.width = .95, point_size = 1.7, size = 1.2, position = position_dodge(w=0.35)) +
  labs( title = "Medium productivity",
        x = "",
        y="Probability of occurrence",
        tag = "b") +
  theme_bw() +
  coord_cartesian(ylim = c(0,0.35)) +
  scale_color_manual(name = "Life history", values = c("#2a9d8f","#f4a261", "#e76f51"), labels = c("Annual/\\nbiennial", "Short-lived \\nperennial", "Long-lived \\nperennial")) + # labels = c(labels)
  scale_x_discrete(labels = c("", "", "", "", "", ""))


high <- post_pred_lifespan %>% 
  filter(p_st == 1.13537307) %>%
  mutate(lifespan = factor(lifespan, levels = c("an_bien", "short_peren", "long_peren"))) %>% 
  ggplot(., aes(y = .value, x =as.factor(sev_int), group = lifespan, color = lifespan)) +
  stat_pointinterval(.width = .95, point_size = 1.7, size = 1.2, position = position_dodge(w=0.35)) +
  labs( title = "High productivity",
        x="Fire Severity", 
        y="Probability of occurrence",
        tag = "c") +
  theme_bw() +
  scale_color_manual(name = "Life history", values = c("#2a9d8f","#f4a261", "#e76f51"), guide = guide_legend(override.aes = list(color = "white"))) +
  theme(
    legend.text = element_text(color = "white"),
    legend.title = element_text(color = "white"),
    axis.title.y = element_text(colour = "white")
  ) +
  coord_cartesian(ylim = c(0,0.35)) +
  scale_x_discrete(labels = c("Unburned", "Low", "Low-moderate", "Moderate", "High-moderate", "High"))


tryptic_lifespan <- gridExtra::grid.arrange( low, med, high, nrow = 3 )




##### MODEL CONTRASTS ---------------------------------------------------

# Difference between long-lived and annual/biennial species at low, medium, and high productivity
# at low productivity
lprod_anbien <- post_pred_lifespan %>% 
  filter(p_st == -1.01221138, lifespan == "an_bien") 
lprod_long <- post_pred_lifespan %>% 
  filter(p_st == -1.01221138, lifespan == "long_peren") 
lprod_diff <- (lprod_anbien[,".value"] - lprod_long[,".value"]) %>% 
  cbind(lprod_anbien[,"sev_int"]) %>% 
  group_by(sev_int) %>% 
  mean_hdci(.value, .width = c(0.95, 0.90))
lprod_diff %>% 
  ggplot(aes(x = sev_int, y = .value, ymin = .lower, ymax = .upper)) +
  geom_pointinterval() +
  geom_hline(yintercept = 0)

# at medium productivity
mprod_anbien <- post_pred_lifespan %>% 
  filter(p_st == 0.01229895, lifespan == "an_bien")
mprod_long <- post_pred_lifespan %>% 
  filter(p_st == 0.01229895, lifespan == "long_peren")
mprod_diff <- (mprod_anbien[,".value"] - mprod_long[,".value"]) %>% 
  cbind(mprod_anbien[,"sev_int"]) %>% 
  group_by(sev_int) %>% 
  mean_hdci(.value, .width = c(0.95, 0.90))
mprod_diff %>% 
  ggplot(aes(x = sev_int, y = .value, ymin = .lower, ymax = .upper)) +
  geom_pointinterval() +
  geom_hline(yintercept = 0)


# at high productivity
hprod_anbien <- post_pred_lifespan %>% 
  filter(p_st == 1.13537307, lifespan == "an_bien")
hprod_long <- post_pred_lifespan %>% 
  filter(p_st == 1.13537307, lifespan == "long_peren")
hprod_diff <- (hprod_anbien[,".value"] - hprod_long[,".value"]) %>% 
  cbind(hprod_anbien[,"sev_int"]) %>% 
  group_by(sev_int) %>% 
  mean_hdci(.value, .width = c(0.95, 0.90))
hprod_diff %>% 
  ggplot(aes(x = sev_int, y = .value, ymin = .lower, ymax = .upper)) +
  geom_pointinterval() +
  geom_hline(yintercept = 0) 







#################################### **DISPERSAL MODEL** ##########################################
##### RUN MODEL ---------------------------------------------------
# model takes a few hours to run so saving as .rds with 'file' argument is recommended
fit_dispersal <- brm(value ~ mo(sev_int)*far*p_st + (far | plot) + (mo(sev_int) | uu), 
            data = d, family = 'bernoulli', iter = 3000, chains = 4, cores = 4, 
            file = "manuscript_data_code/fit_dispersal.rds") # change file argument to save to your computer

##### CHECK MODEL ---------------------------------------------------
plot( fit_dispersal) # check for mixing of chains
# visual posterior predictive check for model fit
pp_check(fit_dispersal, nsamples = 50) # all data together

# model summary and marginal effects
summary(fit_dispersal)


##### MODEL VISUALIZATION ---------------------------------------------------
# get values for productivity quantiles to produce marginal effects for
quantile(d$p_st, c(.15,.5,.85))
# 15%         50%         85% 
# -1.01221138  0.01229895  1.13537307 

# draw and summarize posterior samples for all combinations of fire severity, lifespan, and our three prodcutivity quantiles
post_pred_dispersal <- d %>%
  modelr::data_grid(sev_int = unique(sev_int),  far = unique(far), p_st = c(-1.01221138, 0.01229895, 1.13537307)) %>% 
  add_fitted_draws(fit_dispersal, re_formula = ~0, scale = 'response') 

#plot predictions in a tryptic by creating three separate plots for low, medium, and high productivity
low <- post_pred_dispersal %>% 
  filter(p_st == -1.01221138) %>% 
  ggplot(., aes(y = .value, x = as.factor(sev_int), group = as.factor(far), color = as.factor(far))) +
  stat_pointinterval(.width = .95, point_size = 1.7, size = 1.2, position = position_dodge(w=0.35)) +
  scale_x_discrete(labels = c("Unburned", "Low", "Low-moderate", "Moderate", "High-moderate", "High")) +
  labs(x = NULL, y = "Probability of occurrence", title = "Low productivity", tag = "a") +
  ylim(0,0.35) +
  theme_bw() +
  scale_color_viridis_d(name = "Dispersal \\n strategy", begin = .1, end = .7, guide = guide_legend(override.aes = list(color = "white"))) +
  theme(
    legend.text = element_text(color = "white"),
    legend.title = element_text(color = "white"),
    axis.title.y = element_text(colour = "white")
  ) +
  scale_x_discrete(labels = c("", "", "", "", "", "")) 


med <- post_pred_dispersal %>% 
  filter(p_st == 0.01229895) %>% 
  ggplot(., aes(y = .value, x = as.factor(sev_int), group = as.factor(far), color = as.factor(far))) +
  stat_pointinterval(.width = .95, point_size = 1.7, size = 1.2, position = position_dodge(w=0.35)) +
  labs( title = "Medium productivity",
        x=NULL, 
        y="Probability of occurrence",
        tag = "b") +
  ylim(0,0.35) +
  theme_bw() +
  scale_color_viridis_d(name = "Dispersal \\n strategy", begin = .1, end = .7, labels = c("Near", "Far")) +
  scale_x_discrete(labels = c("", "", "", "", "", "")) 


high <- post_pred_dispersal %>% 
  filter(p_st == 1.13537307) %>%  
  ggplot(., aes(y = .value, x = as.factor(sev_int), group = as.factor(far), color = as.factor(far))) +
  stat_pointinterval(.width = .95, point_size = 1.7, size = 1.2, position = position_dodge(w=0.35)) +
  labs( title = "High productivity",
        x="Fire Severity", 
        y="Probability of occurrence",
        tag = "c") +
  scale_color_viridis_d(name = "Dispersal \\n strategy", begin = .1, end = .7, guide = guide_legend(override.aes = list(color = "white"))) +
  scale_x_discrete(labels = c("Unburned", "Low", "Low-moderate", "Moderate", "High-moderate", "High")) +
  theme_bw() +
  theme(
    legend.text = element_text(color = "white"),
    legend.title = element_text(color = "white"),
    axis.title.y = element_text(colour = "white")
  ) 


tryptic <- grid.arrange( low, med, high, nrow = 3, ncol = 1 )



##### MODEL CONTRASTS ---------------------------------------------------

# Difference between far-dispersed and near-dispersed species at low, medium, and high productivity
# difference for low productivity
lprod_far <- post_pred_dispersal %>% 
  filter(far == 1, p_st == -1.01221138) 
lprod_near <- post_pred_dispersal %>% 
  filter(far == 0, p_st == -1.01221138) 
lprod_diff <- (lprod_far[,".value"] - lprod_near[,".value"]) %>% 
  cbind(lprod_far[,"sev_int"]) %>% 
  group_by(sev_int) %>% 
  mean_hdci(.value, .width = c(0.95, 0.90))
lprod_diff %>% 
  ggplot(aes(x = sev_int, y = .value, ymin = .lower, ymax = .upper)) +
  geom_pointinterval() +
  geom_hline(yintercept = 0)

# difference for medium productivity
mprod_far <- post_pred_dispersal %>% 
  filter(far == 1, p_st == 0.01229895) 
mprod_near <- post_pred_dispersal %>% 
  filter(far == 0, p_st == 0.01229895) 
mprod_diff <- (mprod_far[,".value"] - mprod_near[,".value"]) %>% 
  cbind(mprod_far[,"sev_int"]) %>% 
  group_by(sev_int) %>% 
  mean_hdci(.value, .width = c(0.95, 0.90))
mprod_diff %>% 
  ggplot(aes(x = sev_int, y = .value, ymin = .lower, ymax = .upper)) +
  geom_pointinterval() +
  geom_hline(yintercept = 0)


# difference for high productivity
hprod_far <- post_pred_dispersal %>% 
  filter(far == 1, p_st == 1.13537307) 
hprod_near <- post_pred_dispersal %>% 
  filter(far == 0, p_st == 1.13537307) 
hprod_diff <- (hprod_far[,".value"] - hprod_near[,".value"]) %>% 
  cbind(hprod_far[,"sev_int"]) %>% 
  group_by(sev_int) %>% 
  mean_hdci(.value, .width = c(0.95, 0.90))
hprod_diff %>% 
  ggplot(aes(x = sev_int, y = .value, ymin = .lower, ymax = .upper)) +
  geom_pointinterval() +
  geom_hline(yintercept = 0) 





#################################### **GEOPHYTE MODEL** ##########################################
##### RUN MODEL ---------------------------------------------------
# model takes a few hours to run so saving as .rds with 'file' argument is recommended
fit_geo <- brm(value ~ mo(sev_int)*geo*p_st + (geo | plot) + (mo(sev_int) | uu), 
    data = d, family = 'bernoulli', iter = 3000, chains = 4, cores = 4, file = "manuscript_data_code/fit_geo.rds") # change file argument to your computer

##### CHECK MODEL ---------------------------------------------------
plot( fit_geo) # check for mixing of chains
# visual posterior predictive check for model fit
pp_check(fit_geo, nsamples = 50) # all data together

# model summary and marginal effects
summary(fit_geo)







