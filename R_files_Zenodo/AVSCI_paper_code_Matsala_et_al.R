library(tidyverse)
library(lme4)
library(DHARMa)
library(effects)
library(vegan)
library(gridExtra)
library(grid)
library(ggfortify)
library(lmtest)
library(betareg)
library(jtools)

# read data from GitHUB repo (or use downloaded .csv data): 
df.plots <- read_csv('https://raw.githubusercontent.com/Janzeero-PhD/Chernobyl_young_forests/master/df_plots.csv')
df.all.added <- read_csv('https://raw.githubusercontent.com/Janzeero-PhD/Chernobyl_young_forests/master/df_added.csv')
df.plots$cont_code <- as.factor(df.plots$cont_code)

# update adjusted distances and shares of conifers:
df.plots.upd <- read_csv('https://raw.githubusercontent.com/Janzeero-PhD/Chernobyl_young_forests/master/plots_upd.csv')
df.plots$dist <- df.plots.upd$dist_update[match(df.plots$id,
                                                df.plots.upd$id
                                                )]
df.all.added$dist_stand <- df.plots.upd$dist_update[match(df.all.added$IDPlots,
                                                df.plots.upd$id
)]
df.plots$prop <- df.plots.upd$prop_upd[match(df.plots$id,
                                                df.plots.upd$id
)]
df.all.added$prop_stand <- df.plots.upd$prop_upd[match(df.all.added$IDPlots,
                                                          df.plots.upd$id
)]
# dist - distance to the closest stand where trees are in age to be able to spread seeds
# prop - proportion of conifer tree species in such stand
# elevation - elevation from SRTM DEM (m)
# dist_water - distance to the closest river or channel (vector layers from OSM Rivers)
# cont_code (or Cont) - relative level of soil radioactive contamination (low, medium, and high)

#### Composition ####

### Scots pine basal area proportion

df.plots <- df.plots %>%
  mutate(response = (prop_BA_pine)) %>%
  mutate(dist_stand_scaled = scale(dist),
         prop_stand_scaled = scale(prop),
         elevation_scaled = scale(elevation),
         dist_water_scaled = scale(dist_water))

df.plots %>%
  group_by(cont_code) %>%
  summarise(mean = mean(response),
            se = sd(response) / sqrt(n()))

df.plots <- df.plots %>%
  mutate(.,  response = case_when(
    response == 0 ~ 0.001,
    response == 1 ~ 0.999,
    TRUE ~ response
  ))

fit_null <- betareg(response ~  dist_stand_scaled + prop_stand_scaled + elevation_scaled + dist_water_scaled, data = df.plots)

fit_BA_pine <- betareg(response ~ cont_code + dist_stand_scaled + prop_stand_scaled + elevation_scaled + dist_water_scaled, data = df.plots)


hist(residuals(fit_BA_pine))

summary(fit_BA_pine)

lrtest(fit_null, fit_BA_pine)

effects <- effect("cont_code", fit_BA_pine) %>%
  as.data.frame(.)

### Silver birch basal area proportion

df.plots <- df.plots %>%
  mutate(response = (prop_BA_birch)) %>%
  mutate(dist_stand_scaled = scale(dist),
         prop_stand_scaled = scale(prop),
         elevation_scaled = scale(elevation),
         dist_water_scaled = scale(dist_water))

ggplot(df.plots) +
  geom_histogram(aes(x = response))

df.plots <- df.plots %>%
  mutate(.,  response = case_when(
    response == 0 ~ 0.001,
    response == 1 ~ 0.999,
    TRUE ~ response
  ))

fit_null <- betareg(response ~  dist_stand_scaled + prop_stand_scaled + elevation_scaled + dist_water_scaled, data = df.plots)

fit_BA_birch <- betareg(response ~ cont_code + dist_stand_scaled + prop_stand_scaled + elevation_scaled + dist_water_scaled, data = df.plots)

hist(residuals(fit_BA_birch))

summary(fit_BA_birch)

lrtest(fit_null, fit_BA_birch)
effects <- effect("cont_code", fit_BA_birch) %>%
  as.data.frame(.)

### Shannon diversity

df.plots <- df.plots %>%
  mutate(response = (shannon_exp)) %>%
  mutate(dist_stand_scaled = scale(dist),
         prop_stand_scaled = scale(prop),
         elevation_scaled = scale(elevation),
         dist_water_scaled = scale(dist_water))

ggplot(df.plots) +
  geom_histogram(aes(x = response))

fit_null <- MASS::glm.nb(response ~  dist_stand_scaled + prop_stand_scaled + elevation_scaled + dist_water_scaled, data = df.plots)

fit_shannon_div <- MASS::glm.nb(response ~ cont_code + dist_stand_scaled + prop_stand_scaled + elevation_scaled + dist_water_scaled, data = df.plots)

hist(residuals(fit_shannon_div))

summary(fit_shannon_div)
summ(fit_shannon_div)
lrtest(fit_null, fit_shannon_div)

effects <- effect("cont_code", fit_shannon_div) %>%
  as.data.frame(.)
#### Establishment #####

### Stocking density

df.plots <- df.plots %>%
  mutate(response = (n_large)) %>%
  mutate(dist_stand_scaled = scale(dist),
         prop_stand_scaled = scale(prop),
         elevation_scaled = scale(elevation),
         dist_water_scaled = scale(dist_water))

hist(df.plots$response)

fit_null <- glmer(response ~  dist_stand_scaled + prop_stand_scaled + elevation_scaled + dist_water_scaled + (1 | id), 
                data = df.plots, 
                family = poisson(link = "log"))

fit_deform<- glmer(response ~  cont_code + dist_stand_scaled + prop_stand_scaled + elevation_scaled + dist_water_scaled + (1 | id), 
                 data = df.plots, 
                 family = poisson(link = "log"))

simulationOutput <- simulateResiduals(fittedModel = fit_deform, n = 250)
summary(fit_deform)
summ(fit_deform)
plot(simulationOutput)

lrtest(fit_null, fit_deform)
effects <- effect("cont_code", fit_deform) %>%
  as.data.frame(.)
### Total deformation

df.all.added.deformation <- df.all.added %>%
  filter(alive == 'yes') %>%
  group_by(IDPlots, Cont) %>%
  summarize(deformated = sum(deformated == "yes"),
            total = n(),
            dist_stand = mean(dist_stand, na.rm = TRUE),
            prop_stand = mean(prop_stand, na.rm = TRUE),
            elevation = mean(elevation, na.rm = TRUE),
            dist_water = mean(dist_water, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(dist_stand_scaled = scale(dist_stand),
         prop_stand_scaled = scale(prop_stand),
         elevation_scaled = scale(elevation),
         dist_water_scaled = scale(dist_water))

fit_null <- glmer(cbind(deformated, total - deformated) ~  dist_stand_scaled + prop_stand_scaled + elevation_scaled + dist_water_scaled + (1 | IDPlots), 
                  data = df.all.added.deformation, 
                  family = binomial(link = "logit"))

fit_deform<- glmer(cbind(deformated, total - deformated) ~ Cont + dist_stand_scaled + prop_stand_scaled + elevation_scaled + dist_water_scaled + (1 | IDPlots), 
                   data = df.all.added.deformation, 
                   family = binomial(link = "logit"))

simulationOutput <- simulateResiduals(fittedModel = fit_deform, n = 250)
summary(fit_deform)
summ(fit_deform)
plot(simulationOutput)

lrtest(fit_null, fit_deform)
effects <- effect("Cont", fit_deform) %>%
  as.data.frame(.)
### Scots pine deformation

df.all.added.PISY <- df.all.added %>%
  filter(Species_n == 'PISY') %>%
  filter(alive == 'yes') %>%
  group_by(IDPlots, Cont) %>%
  summarize(deformated = sum(deformated == "yes"),
            total = n(),
            dist_stand = mean(dist_stand, na.rm = TRUE),
            prop_stand = mean(prop_stand, na.rm = TRUE),
            elevation = mean(elevation, na.rm = TRUE),
            dist_water = mean(dist_water, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(dist_stand_scaled = scale(dist_stand),
         prop_stand_scaled = scale(prop_stand),
         elevation_scaled = scale(elevation),
         dist_water_scaled = scale(dist_water))

fit_null <- glmer(cbind(deformated, total - deformated) ~  dist_stand_scaled + prop_stand_scaled + elevation_scaled + dist_water_scaled + (1 | IDPlots), 
                data = df.all.added.PISY, 
                family = binomial(link = "logit"))

fit_deform_pine <- glmer(cbind(deformated, total - deformated) ~ Cont + dist_stand_scaled + prop_stand_scaled + elevation_scaled + dist_water_scaled + (1 | IDPlots), 
                       data = df.all.added.PISY, 
                       family = binomial(link = "logit"))

simulationOutput <- simulateResiduals(fittedModel = fit_deform_pine, n = 250)

plot(simulationOutput)

lrtest(fit_null, fit_deform_pine)
anova(fit_null, fit_deform_pine)
summary(fit_deform_pine)
summ(fit_deform_pine)
effects <- effect("Cont", fit_deform_pine) %>%
  as.data.frame(.)

p <- ggplot(effects, aes(x = factor(Cont, levels = c("Low", "Medium", "High")),
                    y = fit * 100)) +
  geom_point(data = df.all.added.PISY, 
             aes(x = factor(Cont, levels = c("Low", "Medium", "High")), 
                 y = deformated / total * 100),
             position = position_jitter(width = 0.05), col = '#bdbdbd') +
  geom_point() +
  geom_errorbar(aes(ymin = (fit - se) * 100, ymax = (fit + se) * 100), width = 0.1) +
  scale_x_discrete(breaks = c('Low', 'Medium', 'High'),
                   labels = c('Low', 'Medium', 'High')) +
  theme(legend.position = 'none') +
  labs(y = "Proportion of Scots pines deformed (%)",
       x = "Soil contamination level")

ggsave("effect_scottspine_deformation.jpg", dpi = 600, p, width = 3.5, height = 3.5)

### Silver birch deformation

df.all.added.BEPE <- df.all.added %>%
  filter(Species_n == 'BEPE') %>%
  filter(alive == 'yes') %>%
  group_by(IDPlots, Cont) %>%
  summarize(deformated = sum(deformated == "yes"),
            total = n(),
            dist_stand = mean(dist_stand, na.rm = TRUE),
            prop_stand = mean(prop_stand, na.rm = TRUE),
            elevation = mean(elevation, na.rm = TRUE),
            dist_water = mean(dist_water, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(dist_stand_scaled = scale(dist_stand),
         prop_stand_scaled = scale(prop_stand),
         elevation_scaled = scale(elevation),
         dist_water_scaled = scale(dist_water))

fit_null <- glmer(cbind(deformated, total - deformated) ~  dist_stand_scaled + prop_stand_scaled + elevation_scaled + dist_water_scaled + (1 | IDPlots), 
                  data = df.all.added.BEPE, 
                  family = binomial(link = "logit"))

fit_deform_birch <- glmer(cbind(deformated, total - deformated) ~ Cont + dist_stand_scaled + prop_stand_scaled + elevation_scaled + dist_water_scaled + (1 | IDPlots), 
                         data = df.all.added.BEPE, 
                         family = binomial(link = "logit"))

simulationOutput <- simulateResiduals(fittedModel = fit_deform_birch, n = 250)

plot(simulationOutput)
summary(fit_deform_birch)
summ(fit_deform_birch)
lrtest(fit_null, fit_deform_birch)

effects <- effect("Cont", fit_deform_birch) %>%
  as.data.frame(.)
#### Structure ####

### Presence of understory

df.plots <- df.plots %>%
  mutate(response = n_micro > 0) %>%
  mutate(dist_stand_scaled = scale(dist),
         prop_stand_scaled = scale(prop),
         elevation_scaled = scale(elevation),
         dist_water_scaled = scale(dist_water))

fit_null <- glm(response ~  dist_stand_scaled + prop_stand_scaled + elevation_scaled + dist_water_scaled , 
                data = df.plots, family = 'binomial')

fit_under_stocked <- glm(response ~ cont_code + dist_stand_scaled + prop_stand_scaled + elevation_scaled + dist_water_scaled, 
                         data = df.plots, family = 'binomial')
fit_under_stocked_water <- glm(response ~ cont_code + dist_stand_scaled + prop_stand_scaled + elevation_scaled, 
                         data = df.plots, family = 'binomial')
fit_under_stocked_elev <- glm(response ~ cont_code + dist_stand_scaled + prop_stand_scaled + dist_water_scaled, 
                         data = df.plots, family = 'binomial')
summary(fit_under_stocked)
summ(fit_under_stocked)
simulationOutput <- simulateResiduals(fittedModel = fit_under_stocked, n = 250)
summ(glm(response ~ cont_code, data = df.plots, family = 'binomial'))
plot(simulationOutput)

lrtest(fit_null, fit_under_stocked)
lrtest(fit_under_stocked, fit_under_stocked_water)
lrtest(fit_under_stocked, fit_under_stocked_elev)

effects <- effect("cont_code", fit_under_stocked) %>%
  as.data.frame(.)
effects_water <- effect("dist_water_scaled", fit_under_stocked) %>%
  as.data.frame(.)

p <- ggplot(effects, aes(x = factor(cont_code),
                    y = fit)) +
  geom_point() +
  geom_errorbar(aes(ymin = fit - se, ymax = fit + se), width = 0.1) +
  scale_x_discrete(breaks = c('0', '1', '2'),
                   labels = c('Low', 'Medium', 'High')) +
  labs(y = "Probability of tree regeneration",
       x = "Soil contamination level")

ggsave("effect_understory_presence.jpg", dpi = 600, p_ua, width = 4.5, height = 4.5)

### DBH

df.all.added <- df.all.added %>%
  mutate(response = log(DBH)) %>%
  mutate(dist_stand_scaled = scale(dist_stand),
         prop_stand_scaled = scale(prop_stand),
         elevation_scaled = scale(elevation),
         dist_water_scaled = scale(dist_water))

fit_null <- lmer(response ~  dist_stand_scaled + prop_stand_scaled + elevation_scaled + dist_water_scaled + (1 | IDPlots), 
                 data = df.all.added %>% filter(., DBH > 0), REML = FALSE)

fit_dbh <- lmer(response ~ Cont + dist_stand_scaled + prop_stand_scaled + elevation_scaled + dist_water_scaled + (1 | IDPlots), 
                data = df.all.added %>% filter(., DBH > 0), REML = FALSE)

simulationOutput <- simulateResiduals(fittedModel = fit_dbh, n = 250)
summary(fit_dbh)
summ(fit_dbh)
plot(simulationOutput)

hist(residuals(fit_dbh))

lrtest(fit_null, fit_dbh)

effects <- effect("Cont", fit_dbh) %>%
  as.data.frame(.)
### Height

df.all.added.top90height <- df.all.added %>%
  filter(Height_m > quantile(df.all.added$Height_m, 0.9, na.rm = TRUE)) %>%
  mutate(response = log(Height_m)) %>%
  mutate(dist_stand_scaled = scale(dist_stand),
         prop_stand_scaled = scale(prop_stand),
         elevation_scaled = scale(elevation),
         dist_water_scaled = scale(dist_water))

fit_null <- lmer(response ~  dist_stand_scaled + prop_stand_scaled + elevation_scaled + dist_water_scaled + (1 | IDPlots), 
                 data = df.all.added.top90height %>% filter(., DBH > 0), REML = FALSE)

fit_height <- lmer(response ~ Cont + dist_stand_scaled + prop_stand_scaled + elevation_scaled + dist_water_scaled + (1 | IDPlots), 
                data = df.all.added.top90height %>% filter(., DBH > 0), REML = FALSE)

simulationOutput <- simulateResiduals(fittedModel = fit_height, n = 250)
summary(fit_height)
summ(fit_height)
plot(simulationOutput)

lrtest(fit_null, fit_height)

effects <- effect("Cont", fit_height) %>%
  as.data.frame(.)
### Mortality

df.all.added.mortality <- df.all.added %>%
  group_by(IDPlots, Cont) %>%
  summarize(dead = sum(alive == "no"),
            total = n(),
            dist_stand = mean(dist_stand, na.rm = TRUE),
            prop_stand = mean(prop_stand, na.rm = TRUE),
            elevation = mean(elevation, na.rm = TRUE),
            dist_water = mean(dist_water, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(dist_stand_scaled = scale(dist_stand),
         prop_stand_scaled = scale(prop_stand),
         elevation_scaled = scale(elevation),
         dist_water_scaled = scale(dist_water))

fit_null <- glmer(cbind(dead, total - dead) ~  dist_stand_scaled + prop_stand_scaled + elevation_scaled + dist_water_scaled + (1 | IDPlots), 
                  data = df.all.added.mortality, 
                  family = binomial(link = "logit"))

fit_dead <- glmer(cbind(dead, total - dead) ~ Cont + dist_stand_scaled + prop_stand_scaled + elevation_scaled + dist_water_scaled + (1 | IDPlots), 
                  data = df.all.added.mortality, 
                  family = binomial(link = "logit"))

simulationOutput <- simulateResiduals(fittedModel = fit_dead, n = 250)
summary(fit_dead)
summ(fit_dead)
plot(simulationOutput)

lrtest(fit_null, fit_dead)

effects <- effect("Cont", fit_dead) %>%
  as.data.frame(.)
### Structural diversity

df.plots <- df.plots %>%
  mutate(response = (shannon_diam)) %>%
  mutate(dist_stand_scaled = scale(dist),
         prop_stand_scaled = scale(prop),
         elevation_scaled = scale(elevation),
         dist_water_scaled = scale(dist_water))

ggplot(df.plots) +
  geom_histogram(aes(x = response))

fit_null <- lm(response ~  dist_stand_scaled + prop_stand_scaled + elevation_scaled + dist_water_scaled, data = df.plots)

fit_shannon_div <- lm(response ~ cont_code + dist_stand_scaled + prop_stand_scaled + elevation_scaled + dist_water_scaled, data = df.plots)

hist(residuals(fit_shannon_div))

summary(fit_shannon_div)
summ(fit_shannon_div)
lrtest(fit_null, fit_shannon_div)

effects <- effect("cont_code", fit_shannon_div) %>%
  as.data.frame(.)

# for paper cover image:
theme_set(theme_classic())
data.paper <- data.frame(Type = c('def','under','def','under','def','under'),
                         Level = c('1','1','2','2','3','3'),
                         Value = c(5,30.5,24.6,54.4,20.4,18.2),
                         se = c(2.1,8.7,7.4,10.9,6.7,7.2))
colors = c('#440154','#20908d','#fde725')

ggplot(data.paper, aes(Type,Value))+
  geom_col(aes(fill=Level), position='dodge')+
  scale_fill_manual(values = colors)+
  theme(legend.position='none',
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank())+
  labs(y = 'Predicted mean value, %')+
  #scale_x_discrete(breaks = c('def', 'under'),
  #                 labels = c('Scots pine stem deformation', 'Presense of understory in a plot'))
  scale_y_continuous(breaks = c(0,10,20,30,40,50,60))+
  geom_errorbar(aes(ymin=Value-se, ymax=Value+se, group=Level),position=position_dodge(.9),width=0.2)+
  annotate("text", label = c(expression(paste("Mean values of model"))), 
           x = 1, y = 60, size = 4)+
  annotate("text", label = c(expression(paste("contamination effects:"))), 
           x = 1, y = 55, size = 4)+
  annotate("text", label = c(expression(paste("a - Scots pine stem deformation (%)"))), 
           x = 1.15, y = 50, size = 4)+
  annotate("text", label = c(expression(paste("b - Presence of understory"))), 
           x = 0.98, y = 45, size = 4)+
  annotate("text", label = c(expression(paste("trees on the plot (%)"))), 
           x = 1, y = 40, size = 4)+
  annotate("text", label = c(expression(paste("a"))), 
           x = 1, y = -2, size = 4)+
  annotate("text", label = c(expression(paste("b"))), 
           x = 2, y = -2, size = 4)
ggsave("paper.jpg", dpi = 600, width = 4.5, height = 3.5)
