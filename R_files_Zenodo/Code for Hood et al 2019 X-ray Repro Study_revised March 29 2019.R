#### Hood et al 2019 " Prior reproduction alters how mitochondria respond to an oxidative event" ####
#### Correspondence for data analysis and code : Ryan J. Weaver (RyanJWeaver1@gmail.com) ####

library(cowplot)
library(tidyverse)
library(stringr)
library(broom)
library(psych)
library(forcats)
library(nlme)

dat = read_csv("start.dat 3-29-19.csv")
# get data out of column names and into cells
cleandat <- dat %>% 
  # first get all columns in the form organ.complex.metric
  #### gather is a function that takes the names from columns, here 5-50
  ## and puts them in 1 column, with the first argument, here y.group, 
  ## then takes all the entries below each original column, and 
  ## puts it into another column, with the second argument, here value
  gather(y.group, value, 5:50) 

cleandat   

cleandat %>% 
  filter(!str_detect(y.group, "po"),
         !str_detect(y.group, "c2."),
         !str_detect(y.group, "s2"),
         y.group != "brain.cs",
                y.group != "brain.hne",
                y.group != "brain.pc",
                y.group != "hrt.cs",
                y.group != "hrt.hne",
                y.group != "hrt.pc") -> cleandat

write.csv(cleandat, file = "raw data March 29.csv", row.names = F)

#### getting data together####

cleandat %>% 
  group_by(y.group, x.ray, repro) %>% 
  summarise(average = mean(value, na.rm = T),
            sd = sd(value, na.rm = T),
            n = length(value)) -> summary.dat

write.csv(summary.dat, file = "summary repro dat March 29.csv", row.names = F)
#In excel category, measure, part etc were added manually####


summary.dat <- read.csv("summary repro dat March 29.csv")


#### Standard linear model ####

####  Writing the function to run 1 model for all values ####

fit_mod = function(x){
  lm(value ~ repro*x.ray, data = x)
}

#### Mapping the function across all samples and Putting model outputs in a dataframe

modest<- cleandat %>% 
  group_by(y.group) %>% 
  nest(repro, x.ray, value) %>% 
  mutate(mod = map(data, fit_mod),
         tmod = map(mod, tidy)) %>% 
  select(y.group, tmod) %>% 
  unnest() %>%
  mutate(x.ray = ifelse(term == "(Intercept)" |
                          term == "repror", "con", "x"),
         repro = ifelse(str_detect(term, "repro"), "r", "nr"),
         p.value = round(p.value, 5))
modest

modest%>% 
  mutate(termnew = case_when(
    term == "(Intercept)" ~ "Virgin - Control",
    term == "repror" ~ "Reproduced - Control",
    term == "x.rayx" ~ "Virgin - X-ray",
    term == "repror:x.rayx" ~ "Reproduced - X-ray"),
    termnew = fct_relevel(termnew, "Virgin - Control", 
                          "Reproduced - Control", "Virgin - X-ray",
                          "Reproduced - X-ray"))-> estplot

write.csv(modest, file = "Linear Model Estimates March 29.csv", row.names = F)

mod1 = lm(value~group, data = cleandat)

rstandard(mod1, infl = lm.influence(mod1, do.coef = FALSE),
          sd = sqrt(deviance(mod1)/df.residual(mod1)),
          type = c("sd.1", "predictive"))

#### extract standarized residuals function ####

lm_resid = function (x) {
  rstandard(x, infl = lm.influence(x, do.coef = FALSE),
          sd = sqrt(deviance(x)/df.residual(x)),
          type = c("sd.1", "predictive"))
}


#### Extracting model fit values

mod_resids <- cleandat %>% 
  group_by(y.group) %>% 
  nest(repro, x.ray, value) %>% 
  mutate(mod = map(data, fit_mod),
         resids = map(mod, lm_resid),
         fitted = map(mod, "fitted.values")) %>% 
  select(y.group, resids, fitted) %>% 
  unnest() %>% 
  mutate(st.resids = (resids-mean(resids))/sd(resids)) 
mod_resids

jpeg(file = "March 29 Linear Model resid check.jpg", units = "in", width = 14, height = 8.5, res = 500)

mod_resids %>% 
  filter(y.group != "brain.cs",
         y.group != "brain.hne",
         y.group != "brain.pc",
         y.group != "hrt.cs",
         y.group != "hrt.hne",
         y.group != "hrt.pc") %>% 
ggplot( aes(x = fitted, y = resids)) +
  geom_point(size = 2, alpha = 0.9) +
  geom_hline(yintercept= 0) +
  facet_wrap(~y.group, scales = "free")
dev.off()




#### Trying a different model to better fit the data ####
##### Testing for homogeneity of variance among groups


cleandat$group <-as.factor(cleandat$group)

fligner.test(value ~ group, data = cleandat)

# Fligner-Killeen test of homogeneity of variances
# 
# data:  value by group
# Fligner-Killeen:med chi-squared = 16.962, df = 3, p-value = 0.0007197


#### The results above suggest there is unequal variances among groups. ####

######################################
#### GLS Model with unequal variances ####
######################################


####  Writing the function to run 1 model for all values ####

fit_mod_lme = function(x){
  summary(gls(value ~ repro*x.ray, weights = varIdent(form=~ 1|group), data = x, na.action = na.omit))
}

#### Function to extract the results fromt each gls model ####

tidy_gls = function(x){
  tibble(term = names(x$coefficients),
         estimate = x$coefficients,
         std.error = x$tTable[,2],
         statistic = x$tTable[,3],
         p.value = x$tTable[,4])
}

#### Mapping the function across all samples and Putting model outputs in a dataframe

modest_lme<- cleandat %>% 
  group_by(y.group) %>% 
  nest(repro, x.ray, value, group) %>% 
  mutate(mod = map(data, fit_mod_lme),
         tmod = map(mod, tidy_gls)) %>% 
  select(y.group, tmod) %>% 
  unnest() %>% 
  mutate(x.rayx = ifelse(term == "(Intercept)" |
                           term == "repror", "con", "x"),
         repror = ifelse(str_detect(term, "repro"), "r", "nr"),
         p.value = round(p.value, 5))
modest_lme



modest_lme%>% 
  mutate(termnew = case_when(
    term == "(Intercept)" ~ "Virgin - Control",
    term == "repror" ~ "Reproduced - Control",
    term == "x.rayx" ~ "Virgin - X-ray",
    term == "repror:x.rayx" ~ "Reproduced - X-ray"),
    termnew = fct_relevel(termnew, "Virgin - Control", 
                          "Reproduced - Control", "Virgin - X-ray",
                          "Reproduced - X-ray"))-> estplot_lme

write.csv(modest_lme, file = "GLS Model Estimates March 29.csv", row.names = F)


#### Extracting pearson standardized residuals from the gls model ####
# This requires changing the gls function to not report the summary automatically as was done above ###

resid_fit_mod_lme = function(x){
  gls(value ~ repro*x.ray, weights = varIdent(form=~ 1|group), data = x, na.action = na.omit)
}

#### This function extracts the residuals from the gls object in the pearson standarized form 
#### the standardized residuals (raw residuals divided by the corresponding standard errors)

lme_resid = function(x){
  residuals(x, type = "p")
}


#### Extracting model fit values and standardized residuals ####

lme_mod_st_resids <- cleandat %>% 
  group_by(y.group) %>% 
  nest(repro, x.ray,group,value) %>% 
  mutate(mod = map(data, resid_fit_mod_lme),
         resids = map(mod, lme_resid),
         fitted = map(mod, "fitted")) %>% 
  select(y.group, resids, fitted) %>% 
  unnest() 
  
lme_mod_st_resids

jpeg(file = "March 29 GLS Model resid check.jpg", units = "in", width = 14, height = 8.5, res = 500)

lme_mod_st_resids %>% 
  filter(y.group != "brain.cs",
         y.group != "brain.hne",
         y.group != "brain.pc",
         y.group != "hrt.cs",
         y.group != "hrt.hne",
         y.group != "hrt.pc") %>% 
ggplot( aes(x = fitted, y = resids)) +
  geom_point(size = 2, alpha = 0.9) +
  geom_hline(yintercept= 0) +
  facet_wrap(~y.group, scales = "free")
dev.off()


#### Analysis of litter size at weaning and x-ray effect on mito measures ####

pupdat <- cleandat %>% 
  filter(repro == "r")

write.csv(pupdat, file = "Repro only dat March 29.csv", row.names = F)
#### Manually added pup data for each individual

pupdat <- read.csv("Repro only dat March 29.csv")

#### Discretizing reproductive output to normal and high ####
pupdat %>%
  mutate(rep.out = ifelse(pup >10, "high", "normal")) -> pupdat

pupdat %>% 
  group_by(id, repro, x.ray, group, pup, rep.out) %>% 
  summarize(pup.wean = mean(pup))->pupdat2

jpeg(file = "Pup histogram by xray March 29.jpg", units = "in", width = 8, height = 5.5, res = 500)

ggplot(pupdat2, aes(x = pup, fill = x.ray))+
  geom_histogram(binwidth = 1, position = "dodge")+
  scale_fill_manual(values = c("gray30","red"),
                     labels = c("Control", "X-ray"),
                     name = "") 
dev.off()

#### Is there a difference in pup number between xray and control groups?####

pupmod <-glm(pup ~ x.ray, family = "poisson", data = pupdat2)

summary(pupmod)

# Call:
#   glm(formula = pup ~ x.ray, family = "poisson", data = pupdat2)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -2.2391  -1.0447   0.2009   0.8936   1.6305  
# 
# Coefficients:
#               Estimate Std. Error z value Pr(>|z|)    
#   (Intercept)  2.54553    0.09901  25.709   <2e-16 ***
#   x.rayx      -0.07743    0.13520  -0.573    0.567    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for poisson family taken to be 1)
# 
# Null deviance: 24.471  on 17  degrees of freedom
# Residual deviance: 24.144  on 16  degrees of freedom
# AIC: 105.48
# 
# Number of Fisher Scoring iterations: 4

### No evidence of a difference, makes sense because they were irradiated after they gave birth ####
pupdat %>% 
  filter(y.group == "sm.c1.rcr") %>% 
  mutate(scalepup = (pup-mean(pup))/sd(pup) )->testpup
testpup

library(lme4)

lmemod <- gls(value ~  pup*x.ray, data = testpup, na.action = na.omit,weights = varIdent (form = ~1|x.ray))
summary(lmemod, verbose = TRUE)

# Generalized least squares fit by REML
# Model: value ~ pup * x.ray 
# Data: testpup 
# AIC      BIC    logLik
# 247.3115 251.1458 -117.6557
# 
# Variance function:
#   Structure: Different standard deviations per stratum
# Formula: ~1 | x.ray 
# Parameter estimates:
#   con       x 
# 1.00000 2.40374 
# Coefficients:
#                 Value Std.Error   t-value p-value  Interpretation
# (Intercept)  3177.393  440.6841  7.210137  0.0000   Mean CS value when of zero pups in the control group
# pup           -83.982   32.8011 -2.560350  0.0227   For every 1 pup, CS decreased by 83 points in the control goup
# x.rayx      -2812.527 1083.1396 -2.596643  0.0211   X-ray decreases mean CS activity by 2812 points on averageif they had zero pups
# pup:x.rayx    115.375   86.4071  1.335245  0.2031   The number of pups doesn't really affect CS activity in the X-ray group
                                                      #115.37 - 83.98 = 32 so for every 1 pup in the xray group, CS increased by 32



#### Model function 
fit_mod_pup = function(x){
  summary(gls(value~pup*x.ray,weights = varIdent (form = ~1|x.ray),
              na.action = na.omit ,data = x))
}

tidy_gls = function(x){
  tibble(term = names(x$coefficients),
         estimate = x$coefficients,
         std.error = x$tTable[,2],
         statistic = x$tTable[,3],
         p.value = x$tTable[,4])
}

#### Mapping the function across all samples and Putting model outputs in a dataframe

modest.pup<- pupdat %>%
  filter(y.group != "brain.cs",
         y.group != "brain.hne",
         y.group != "brain.pc",
         y.group != "hrt.cs",
         y.group != "hrt.hne",
         y.group != "hrt.pc") %>%
  group_by(y.group) %>% 
  nest(value, pup, x.ray) %>% 
  mutate(mod = map(data, fit_mod_pup),
         tmod = map(mod, tidy_gls)) %>% 
  select(y.group, tmod) %>% 
  unnest() %>% 
  mutate(p.value = round(p.value, 5))

modest.pup

write.csv(modest.pup, file = "GLS Pup random Model X on Repro Results March 29.csv", row.names = F)

#### Extracting pearson standardized residuals from the gls model ####
# This requires changing the gls function to not report the summary automatically as was done above ###

resid_fit_mod_lme = function(x){
  gls(value~pup*x.ray,weights = varIdent (form = ~1|x.ray),
      na.action = na.omit ,data = x)
}

#### This function extracts the residuals from the gls object in the pearson standarized form 
#### the standardized residuals (raw residuals divided by the corresponding standard errors)

lme_resid = function(x){
  residuals(x, type = "p")
}


#### Extracting model fit values and standardized residuals ####

lme_mod_st_resids <- pupdat %>% 
  filter(y.group != "brain.cs",
         y.group != "brain.hne",
         y.group != "brain.pc",
         y.group != "hrt.cs",
         y.group != "hrt.hne",
         y.group != "hrt.pc") %>% 
  group_by(y.group) %>% 
  nest(value, pup, x.ray) %>% 
  mutate(mod = map(data, resid_fit_mod_lme),
         resids = map(mod, lme_resid),
         fitted = map(mod, "fitted")) %>% 
  select(y.group, resids, fitted) %>% 
  unnest() 

lme_mod_st_resids

jpeg(file = "March 29 GLS Pup Model resid check.jpg", units = "in", width = 14, height = 8.5, res = 500)

ggplot(lme_mod_st_resids, aes(x = fitted, y = resids)) +
  geom_point(size = 2, alpha = 0.9) +
  geom_hline(yintercept= 0) +
  facet_wrap(~y.group, scales = "free")
dev.off()




#### Plotting Effect of pups on mito measures



fmt_dcimals <- function(decimals=0){
  function(x) format(x,nsmall = decimals,scientific = FALSE)
}

#### The effect of pups on Muscle RCR Figure

jpeg(file = "April 3 Pup Muscle RCR effect.jpg", units = "in", width = 8.5, height = 5.5, res = 500)

pupdat %>% 
  filter(y.group == "sm.c1.rcr") %>%
  ggplot(aes(x = pup,y = value, col = x.ray)) +
  geom_smooth(aes(col = x.ray, lty = x.ray), method="lm", alpha = 0.2, se = FALSE) +
  geom_point(aes(col =x.ray), size = 6, alpha =0.7) +
    scale_color_manual(values = c("gray30","red"),
                     labels = c("Control", "X-ray"),
                     name = "") +
  scale_linetype_manual(values = c(1, 3),
                        labels = c("Control", "X-ray"),
                        name = "") +
  xlim(5,20)+
    xlab("Pups Weaned") +
  ylab("Respiratory Control Ratio")+
  scale_y_continuous(limits = c(2,10),
    labels = fmt_dcimals(1))+
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 22))+
  theme(legend.position = "")
dev.off()

#### The effect of pups on Muscle RCR Figure



jpeg(file = "April 3 Pup Liver RCR effect.jpg", units = "in", width = 8.5, height = 5.5, res = 500)

pupdat %>% 
  filter(y.group == "liv.c1.rcr") %>%
  ggplot(aes(x = pup,y = value, col = x.ray)) +
  geom_smooth(aes(col = x.ray, lty = x.ray), method="lm", alpha = 0.2, se = FALSE) +
    geom_point(aes(col =x.ray), size = 6, alpha =0.7) +
  scale_color_manual(values = c("gray30","red"),
                     labels = c("Control", "X-ray"),
                     name = "") +
    scale_linetype_manual(values = c(1, 3),
                           labels = c("Control", "X-ray"),
                           name = "") +
  xlim(5,20)+
  xlab("Pups Weaned") +
  ylab("Respiratory Control Ratio")+
  scale_y_continuous(limits = c(3,10),
    labels = fmt_dcimals(1))+
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 22))+
  theme(legend.position = "top")
dev.off()





#### THe effect of pups and x-ray on Muscle CS Figure ####

jpeg(file = "April 3 Pup Muscle CS effect.jpg", units = "in", width = 8.5, height = 5.5, res = 500)

pupdat %>% 
  filter(y.group == "sm.cs") %>%
  ggplot(aes(x = pup,y = value, col = x.ray)) +
  geom_smooth(aes(col = x.ray, lty = x.ray), method="lm", alpha = 0.2, se = FALSE) +
   geom_point(aes(col =x.ray), size = 6, alpha =0.7) +
  scale_color_manual(values = c("gray30","red"),
                     labels = c("Control", "X-ray"),
                     name = "") +
  scale_linetype_manual(values = c(1, 3),
                        labels = c("Control", "X-ray"),
                        name = "") +
  xlim(5,20)+
  xlab("Pups Weaned") +
  ylab(expression(paste("CS activity" ~~(nmol~min^-1~mg~protein^-1))))+
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 22))+
  theme(legend.position = "top")
dev.off()

#### The effect of pups on Liver H2 production Figure ####

jpeg(file = "April 3 Pup Liver H2 effect.jpg", units = "in", width = 8.5, height = 5.5, res = 500)

 pupdat %>% 
  filter(y.group == "liv.h2") %>%
  ggplot(aes(x = pup,y = value, col = x.ray)) +
  geom_smooth(aes(col = x.ray, lty = x.ray), method="lm", alpha = 0.2, se = FALSE) +
  geom_point(aes(col =x.ray), size = 6, alpha =0.7) +
  scale_color_manual(values = c("gray30","red"),
                     labels = c("Control", "X-ray"),
                     name = "") +
  scale_linetype_manual(values = c(1, 3),
                        labels = c("Control", "X-ray"),
                        name = "") +
  xlim(5,20)+
  xlab("Pups Weaned") +
  ylab(expression(paste("H2O2 production" ~(pmol~min^-1~mg~protein^-1))))+
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 22))+
  theme(legend.position = "")
dev.off()


#### Putting all 4 plots together ####

comboplot<- plot_grid(pupfig1,pupfig2,pupfig3,pupfig4, align = "hv", 
                      labels = c("A", "B","C","D"), scale = 0.9)
save_plot(comboplot, file = "Pup Effect Plot combo.jpg", base_width = 14, base_height = 8.5)

jpeg(file = "March 29 Pup Xray Interaction Figures.jpg", units = "in", width = 14, height = 8.5, res = 500)

pupdat %>% 
  filter(y.group != "brain.cs",
         y.group != "brain.hne",
         y.group != "brain.pc",
         y.group != "hrt.cs",
         y.group != "hrt.hne",
         y.group != "hrt.pc") %>% 
  ggplot(aes(x = pup,y = value , col = x.ray)) +
  #geom_smooth(method="gam", aes(fill = x.ray), alpha = 0.2)+
  geom_point( size = 2) +
  scale_color_manual(values = c("gray30","red"),
                     labels = c("Control", "X-ray"),
                     name = "") +
   # scale_fill_manual(values = c("gray30","red"),
   #                   labels = c("Control", "X-ray"),
   #                   name = "") +
  xlim(5,20)+
  xlab("Pups Weaned") +
  ylab("Response")+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16))+
  theme(legend.position = "top")+
  facet_wrap(~y.group, scales = "free")
dev.off()


#### Calculate Hedges' g and 95% CI ####


summary.dat %>% 
  mutate(repro = ifelse(grepl("nr", repro), "Virgin", "Reproduced"),
         x.ray= ifelse(grepl("con", x.ray), "Control", "X-ray")) -> summary.dat

#### Function for Calculating Hedges's g and the associated variables #####


cohensd = function(x){
  ctrl = filter(x, level == 0)

  x %>% 
    mutate(var.pooled = ((n - 1)*sd^2 + (ctrl$n - 1)*ctrl$sd^2)/(n + ctrl$n - 2),
           s.pooled = sqrt(var.pooled),
           n.pooled= n + ctrl$n,
           j = 1-(3/((4*n.pooled)-9)),
           g = ((average - ctrl$average)/s.pooled)*j) 
  
}

extract.lci = function(g, n.pooled){
  cohen.d.ci(g, n.pooled)[1]
}

extract.uci = function(g, n.pooled){
  cohen.d.ci(d = g, n = n.pooled)[3]
}

summary.dat %>% 
  as.tibble() %>% 
  group_by(y.group) %>% 
  nest(average, sd, n, level, x.ray, repro, measure, part, category) %>%
  mutate(g = map(data, cohensd)) %>%
  unnest(g) %>% 
  mutate(lci = map2_dbl(g, n.pooled, extract.lci), 
         uci = map2_dbl(g, n.pooled, extract.uci)) -> newplotdatg



newplotdatg %>% 
  mutate(group = paste(repro,x.ray, sep = "-")) -> newplotdatg

newplotdatg$measure <- factor(newplotdatg$measure, levels = c("CS","RCR",
                          "S3","S4", "C1" ,"C2","C3","C4","H2","HNE","PC"))

write.csv(newplotdatg, file = "Final Data Summary_g March 29.csv", row.names = F)


#### Liver Summary Figures ####


newplotdatg$group <-factor(newplotdatg$group, levels = c("Virgin-Control", "Reproduced-Control", "Virgin-X-ray", "Reproduced-X-ray"))



# Complex Activity 

jpeg(file = "March 29 Complex Activity Liver Summary.jpg", units = "in", width = 8.5, height = 5.5, res = 500)

newplotdatg %>% 
  filter(category == "Activity") %>% 
  filter(measure != "CS",
         level != 0,
         part == "liv") %>% 
  #group != "Reproduced-X-ray") %>%
  ggplot(aes(x = group, col = measure)) +
  geom_point(aes(y = g), size = 8, alpha = 0.7,
             position = position_dodge(width = 0.7)) +
  geom_hline(yintercept = 0, lty =2)+
  geom_errorbar( aes(ymin = lci, ymax = uci),
                 position = position_dodge(width = 0.7), width = 0.2) +
  scale_color_manual(values = c( "seagreen3", "darkgoldenrod2", "sienna4","gray20"),
                     labels = c("Complex I", "Complex II", "Complex III", "Complex IV"),
                     name = "") +
  scale_x_discrete(labels = c('Reproduction', 'X-ray', "X-ray + Reproduction"))+
  xlab("") +
  ylab("Hedges' g") +
  ylim(-1.5,2.5)+
  #ggtitle("Liver Complex Activity") + theme(plot.title = element_text(hjust=1))+
  theme(legend.position = "top")+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))
dev.off()


# Damage Markers

jpeg(file = "March 29 Damage Marker Liver Summary .jpg", units = "in", width = 8.5, height = 5.5, res = 500)

newplotdatg %>% 
  filter(category == "Damage") %>% 
  filter(level != 0,
         part == "liv" )%>% 
  #group != "Reproduced-X-ray") %>% 
  ggplot(aes(x = group, col = measure)) +
  geom_point(aes(y = g, shape = measure, color = measure), size = 8, alpha = 0.7,
             position = position_dodge(width = 0.7)) +
  scale_shape_manual(values=c(15, 17, 17),
                     labels = c("H2O2", "4-HNE", "Protein Carbonyl"),
                     name = "") +
  geom_hline(yintercept = 0, lty =2)+
  geom_errorbar( aes(ymin = lci, ymax = uci),
                 position = position_dodge(width = 0.7), width = 0.2) +
  scale_color_manual(values = c("darkgoldenrod2", "sienna4","gray20"),
                     labels = c("H2O2", "4-HNE", "Protein Carbonyl"),
                     name = "") +
  scale_x_discrete(labels = c('Reproduction', 'X-ray', "X-ray + Reproduction"))+
  xlab("") +
  ylab("Hedges' g") +
  ylim(-2.5,2)+
  theme(legend.position = "top")+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))
dev.off()

# Respiration Measures 

jpeg(file = "March 29 Respiration Measures Liver Summary.jpg", units = "in", width = 8.5, height = 5.5, res = 500)

newplotdatg %>% 
  filter(category == "Respiration" | measure == "CS" ) %>% 
  filter(level != 0,
         part == "liv",
         !grepl("c2", y.group),
         #group != "Reproduced-X-ray",
         measure != "PO") %>% 
  ggplot(aes(x = group, col = measure)) +
  geom_point(aes(y = g, shape = measure, col = measure), size = 8, alpha = 0.7,
             position = position_dodge(width = 0.7)) +
  scale_shape_manual(values=c(15,16, 16, 16,17),
                     labels = c( "CS","RCR", "state 3", "state 4"),
                     name = "")+
  geom_hline(yintercept = 0, lty =2)+
  geom_errorbar( aes(ymin = lci, ymax = uci),
                 position = position_dodge(width = 0.7), width = 0.2) +
  scale_color_manual(values = c("slate blue","chartreuse4", "palegreen3", "paleturquoise4"),
                     labels = c( "CS","RCR", "state 3", "state 4"),
                     name = "") +
  scale_x_discrete(labels = c('Reproduction', 'X-ray', "X-ray + Reproduction"))+
  xlab("") +
  ylab("Hedges' g") +
  theme(legend.position = "top")+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))
dev.off()

#### Skeletal Muscle Summary Figures ####

# Complex Activity #

jpeg(file = "March 29 Complex Activity Muscle Summary.jpg", units = "in", width = 8.5, height = 5.5, res = 500)

newplotdatg %>% 
  filter(category == "Activity") %>% 
  filter(measure != "CS",
         level != 0,
         part == "sm" )%>% 
  #group != "Reproduced-X-ray") %>% 
  ggplot(aes(x = group, col = measure)) +
  geom_point(aes(y = g), size = 8, alpha = 0.7,
             position = position_dodge(width = 0.7)) +
  geom_hline(yintercept = 0, lty =2)+
  geom_errorbar( aes(ymin = lci, ymax = uci),
                 position = position_dodge(width = 0.7), width = 0.2) +
  scale_color_manual(values = c( "seagreen3", "darkgoldenrod2", "sienna4","gray20"),
                     labels = c("Complex I", "Complex II", "Complex III", "Complex IV"),
                     name = "") +
  scale_x_discrete(labels = c('Reproduction', 'X-ray', "X-ray + Reproduction"))+
  xlab("") +
  ylim(-4,2)+
  ylab("Hedges' g") +
  theme(legend.position = "top")+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))

dev.off()

# Damage Markers 

jpeg(file = "March 29 Damage Marker Muscle Summary.jpg", units = "in", width = 8.5, height = 5.5, res = 500)

newplotdatg %>% 
  filter(category == "Damage") %>% 
  filter(level != 0,
         part == "sm" )%>% 
  #group != "Reproduced-X-ray") %>% 
  ggplot(aes(x = group, col = measure)) +
  geom_point(aes(y = g, shape = measure, color = measure), size = 8, alpha = 0.7,
             position = position_dodge(width = 0.7)) +
  scale_shape_manual(values=c(15, 17, 17),
                     labels = c("H2O2", "4-HNE", "Protein Carbonyl"),
                     name = "") +
  geom_hline(yintercept = 0, lty =2)+
  geom_errorbar( aes(ymin = lci, ymax = uci),
                 position = position_dodge(width = 0.7), width = 0.2) +
  scale_color_manual(values = c("darkgoldenrod2", "sienna4","gray20"),
                     labels = c("H2O2", "4-HNE", "Protein Carbonyl"),
                     name = "") +
  scale_x_discrete(labels = c('Reproduction', 'X-ray', "X-ray + Reproduction"))+
  xlab("") +
  ylab("Hedges' g") +
  ylim(-1.7,1.5)+
  theme(legend.position = "top")+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))
dev.off()

# Respiration Measures (no complex 2)

jpeg(file = "March 29 Respiration Measures Muscle Summary.jpg", units = "in", width = 8.5, height = 5.5, res = 500)

newplotdatg %>% 
  filter(category == "Respiration" | measure == "CS" ) %>% 
  filter(level != 0,
         part == "sm",
         !grepl("c2", y.group),
         #group != "Reproduced-X-ray",
         measure != "PO") %>% 
  ggplot(aes(x = group, col = measure)) +
  geom_point(aes(y = g, shape = measure, col = measure), size = 8, alpha = 0.7,
             position = position_dodge(width = 0.7)) +
  scale_shape_manual(values=c(15,16, 16, 16,17),
                     labels = c( "CS","RCR", "state 3", "state 4"),
                     name = "")+
  geom_hline(yintercept = 0, lty =2)+
  geom_errorbar( aes(ymin = lci, ymax = uci),
                 position = position_dodge(width = 0.7), width = 0.2) +
  scale_color_manual(values = c("slate blue","chartreuse4", "palegreen3", "paleturquoise4"),
                     labels = c( "CS","RCR", "state 3", "state 4"),
                     name = "") +
  scale_x_discrete(labels = c('Reproduction', 'X-ray', "X-ray + Reproduction"))+
  xlab("") +
  ylab("Hedges' g") +
  ylim(-2.5,2)+
  theme(legend.position = "top")+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))
dev.off()



#################################
################################

#### Interaction plots ####




dat %>% mutate(repro = ifelse(grepl("nr", repro), "Virgin", "Reproduced"),
               x.ray= ifelse(grepl("con", x.ray), "Control", "X-ray")) -> dat
dat$repro <- as.factor(dat$repro)
dat$repro <- relevel(dat$repro, ref ="Virgin")
dat$x.ray <- as.factor(dat$x.ray)
dat$x.ray <- relevel(dat$x.ray, ref ="Control")

newplotdatg$repro <- as.factor(newplotdatg$repro)
newplotdatg$repro <- relevel(newplotdatg$repro, ref ="Virgin")
newplotdatg$x.ray <- as.factor(newplotdatg$x.ray)
newplotdatg$x.ray <- relevel(newplotdatg$x.ray, ref ="Control")


#### Extracting model estimated standard error ####
modest_lme %>%
  filter(y.group =="sm.c1.s3")->mod.lme.se

newplotdatg %>% 
  filter(y.group =="sm.c1.s3") %>% 
  mutate(se = mod.lme.se$std.error)%>% 
  mutate(average = average*1000,
         se = se*1000,
         sd = sd*1000) -> filtdat_lme

filtdat_lme$repro <- as.factor(filtdat_lme$repro)
filtdat_lme$repro <- relevel(filtdat_lme$repro, ref ="Virgin")
filtdat_lme$x.ray <- as.factor(filtdat_lme$x.ray)
filtdat_lme$x.ray <- relevel(filtdat_lme$x.ray, ref ="Control")



jpeg(file = "March 29 LME Muscle S3 Respiration Interaction.jpg", units = "in", width = 5.5, height = 5.5, res = 500)


filtdat_lme %>% 
ggplot(aes(x = x.ray, col = repro)) +
  geom_point(data = dat, aes(y = sm.c1.s3*1000, col = repro),size = 4, 
             alpha =0.7,position = position_dodge(width = 0.4))+
  geom_point( aes(y = average, fill =x.ray), size = 8, shape =15, alpha = 0.7,
              position = position_dodge(width = 0.4))+
  stat_summary(fun.y = mean, geom="line", aes(x = filtdat_lme$x.ray, y = filtdat_lme$average,
                                              group = factor(repro)), 
               position = position_dodge(width = 0.4)) +
  geom_errorbar( aes(ymin = average - se, ymax = average + se),
                 position = position_dodge(width = 0.4), width = 0.2) +
  scale_color_manual(values = c("gray30", "red"),
                     labels = c("Virgin", "Reproduced"),
                     name = "") +
  scale_fill_manual(values = c("gray30", "red"),
                    labels = c("Virgin", "Reproduced"),
                    name = "") +
  xlab("Exposure") +
  ylab(expression(paste("O"[2]*~"consumption" ~~(pmol~min^-1~CS~activity^-1))))+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) +
  theme(legend.position = "top")

dev.off()

#### Muscle Complex I activity Figure ####

modest_lme %>%
  filter(y.group =="sm.c1")->mod.lme.se

newplotdatg %>% 
  filter(y.group =="sm.c1") %>% 
  mutate(se = mod.lme.se$std.error) -> filtdat_lme

filtdat_lme$repro <- as.factor(filtdat_lme$repro)
filtdat_lme$repro <- relevel(filtdat_lme$repro, ref ="Virgin")
filtdat_lme$x.ray <- as.factor(filtdat_lme$x.ray)
filtdat_lme$x.ray <- relevel(filtdat_lme$x.ray, ref ="Control")

jpeg(file = "March 29 LME Muscle Complex I Activity Interaction.jpg", units = "in", width = 5.5, height = 5.5, res = 500)


filtdat_lme %>% 
  ggplot(aes(x = x.ray, col = repro)) +
  geom_point(data = dat, aes(y = sm.c1, col = repro),size = 4, 
             alpha =0.7,position = position_dodge(width = 0.4))+
  geom_point( aes(y = average, fill =x.ray), size = 8, shape =15, alpha = 0.7,
              position = position_dodge(width = 0.4))+
  stat_summary(fun.y = mean, geom="line", aes(x = filtdat_lme$x.ray, y = filtdat_lme$average,
                                              group = factor(repro)), 
               position = position_dodge(width = 0.4)) +
  geom_errorbar( aes(ymin = average - se, ymax = average + se),
                 position = position_dodge(width = 0.4), width = 0.2) +
  scale_color_manual(values = c("gray30", "red"),
                     labels = c("Virgin", "Reproduced"),
                     name = "") +
  scale_fill_manual(values = c("gray30", "red"),
                    labels = c("Virgin", "Reproduced"),
                    name = "") +
  xlab("Exposure") +
  ylab(expression(paste("Complex I Activity" ~~(nmol~min^-1~mg~protein^-1))))+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) +
  theme(legend.position = "top")

dev.off()

#### Muscle CS activity Figure ####

modest_lme %>%
  filter(y.group =="sm.cs")->mod.lme.se

newplotdatg %>% 
  filter(y.group =="sm.cs") %>% 
  mutate(se = mod.lme.se$std.error) -> filtdat_lme

filtdat_lme$repro <- as.factor(filtdat_lme$repro)
filtdat_lme$repro <- relevel(filtdat_lme$repro, ref ="Virgin")
filtdat_lme$x.ray <- as.factor(filtdat_lme$x.ray)
filtdat_lme$x.ray <- relevel(filtdat_lme$x.ray, ref ="Control")

jpeg(file = "March 29 LME Muscle CS Activity Interaction.jpg", units = "in", width = 5.5, height = 5.5, res = 500)


filtdat_lme %>% 
  ggplot(aes(x = x.ray, col = repro)) +
  geom_point(data = dat, aes(y = sm.cs, col = repro),size = 4, 
             alpha =0.7,position = position_dodge(width = 0.4))+
  geom_point( aes(y = average, fill =x.ray), size = 8, shape =15, alpha = 0.7,
              position = position_dodge(width = 0.4))+
  stat_summary(fun.y = mean, geom="line", aes(x = filtdat_lme$x.ray, y = filtdat_lme$average,
                                              group = factor(repro)), 
               position = position_dodge(width = 0.4)) +
  geom_errorbar( aes(ymin = average - se, ymax = average + se),
                 position = position_dodge(width = 0.4), width = 0.2) +
  scale_color_manual(values = c("gray30", "red"),
                     labels = c("Virgin", "Reproduced"),
                     name = "") +
  scale_fill_manual(values = c("gray30", "red"),
                    labels = c("Virgin", "Reproduced"),
                    name = "") +
  xlab("Exposure") +
  ylab(expression(paste("Citrate Synthase Activity" ~~(nmol~min^-1~mg~protein^-1))))+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) +
  theme(legend.position = "top")

dev.off()

#### Liver CS activity Figure ####

modest_lme %>%
  filter(y.group =="liv.cs")->mod.lme.se

newplotdatg %>% 
  filter(y.group =="liv.cs") %>% 
  mutate(se = mod.lme.se$std.error) -> filtdat_lme

filtdat_lme$repro <- as.factor(filtdat_lme$repro)
filtdat_lme$repro <- relevel(filtdat_lme$repro, ref ="Virgin")
filtdat_lme$x.ray <- as.factor(filtdat_lme$x.ray)
filtdat_lme$x.ray <- relevel(filtdat_lme$x.ray, ref ="Control")

jpeg(file = "March 29 LME Liver CS Activity Interaction.jpg", units = "in", width = 5.5, height = 5.5, res = 500)


filtdat_lme %>% 
  ggplot(aes(x = x.ray, col = repro)) +
  geom_point(data = dat, aes(y = liv.cs, col = repro),size = 4, 
             alpha =0.7,position = position_dodge(width = 0.4))+
  geom_point( aes(y = average, fill =x.ray), size = 8, shape =15, alpha = 0.7,
              position = position_dodge(width = 0.4))+
  stat_summary(fun.y = mean, geom="line", aes(x = filtdat_lme$x.ray, y = filtdat_lme$average,
                                              group = factor(repro)), 
               position = position_dodge(width = 0.4)) +
  geom_errorbar( aes(ymin = average - se, ymax = average + se),
                 position = position_dodge(width = 0.4), width = 0.2) +
  scale_color_manual(values = c("gray30", "red"),
                     labels = c("Virgin", "Reproduced"),
                     name = "") +
  scale_fill_manual(values = c("gray30", "red"),
                    labels = c("Virgin", "Reproduced"),
                    name = "") +
  xlab("Exposure") +
  ylab(expression(paste("Citrate Synthase Activity" ~~(nmol~min^-1~mg~protein^-1))))+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) +
  theme(legend.position = "top")

dev.off()


#### Treating each group as a level ####

######################################
#### GLS Model with unequal variances ####
######################################

cleandat %>% 
  mutate(treatgroup = case_when(repro == "nr" & x.ray == "con" ~ "nr.con", 
                                repro == "nr" & x.ray == "x" ~ "nr.xray",
                                repro == "r" & x.ray == "con" ~ "repro.con",
                                repro == "r" & x.ray == "x" ~ "repro.xray")) ->cleandat2

####  Writing the function to run 1 model for all values ####

fit_mod_lme2 = function(x){
  summary(gls(value ~ treatgroup, weights = varIdent(form=~ 1|group), data = x, na.action = na.omit))
}

#### Function to extract the results fromt each gls model ####

tidy_gls = function(x){
  tibble(term = names(x$coefficients),
         estimate = x$coefficients,
         std.error = x$tTable[,2],
         statistic = x$tTable[,3],
         p.value = x$tTable[,4])
}

#### Mapping the function across all samples and Putting model outputs in a dataframe

modest_lme2<- cleandat2 %>% 
  group_by(y.group) %>% 
  nest(treatgroup,value, group) %>% 
  mutate(mod = map(data, fit_mod_lme2),
         tmod = map(mod, tidy_gls)) %>% 
  select(y.group, tmod) %>% 
  unnest() %>% 
  mutate(p.value = round(p.value, 7))

modest_lme2



modest_lme2%>% 
  mutate(termnew = case_when(
    term == "(Intercept)" ~ "Virgin - Control",
    term == "treatgroupnr.xray" ~ "Virgin - X-ray",
    term == "treatgrouprepro.con" ~ "Reproduced - Control",
    term == "treatgrouprepro.xray" ~ "Reproduced - X-ray"),
    termnew = fct_relevel(termnew, "Virgin - Control", 
                          "Reproduced - Control", "Virgin - X-ray",
                          "Reproduced - X-ray"))-> estplot_lme2

write.csv(estplot_lme2, file = "4 group GLS Model Estimates April 1.csv", row.names = F)


