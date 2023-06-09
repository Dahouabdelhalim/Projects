## regression models explaining sign of total detritus TCAP within sites

library(tidyverse)
library(vegan)
library(car)
'%nin%' <- function(x,y)!('%in%'(x,y))

#input files here

syncsa_summary_list <-readRDS("syncsa_summary_sept2021.rds") #from 06_script
site_covariates<-read_csv("site_covariates.csv")
allaxes1234 <- read_csv("allaxes1234.csv")

###############################
#                             #
#     total detritus          #
#                             #
###############################

det_covariates <- site_covariates %>% 
  filter(include.detritus==1)

#combine all explanatory variables

detritus_reg <- syncsa_summary_list$single_total_det.biomass.log.axis %>% 
  rownames_to_column(var = "sites") %>% 
  left_join(det_covariates)

det_covariates <- detritus_reg %>% 
  select(sites, n:mean_logdetritus) %>% 
  mutate(include.detritus = 1)


m1 <- glm(log(TCAP.p) ~ n + mean_logdetritus+ BC2 + BC4 + BC15 + BC17 + chao + andes, family = gaussian, data = detritus_reg)

Anova(m1, type = 2) #negative effect of n and andes north on log p (i.e. n and andes north associted with significance), 
#positive effect of cv env on log p (i.e. cv detritus associated with non-sig)

m2 <- glm(TCAP.rho~ n +mean_logdetritus + BC2 + BC4 + BC15 + BC17 + chao + andes, family = gaussian, data = detritus_reg)

Anova(m2, type = 2) #positive effect of andes north on rho


#######################
#                     #
# water bromeliads    #
#                     #
#######################

wat_covariates <- site_covariates %>% 
  filter(include.water==1)

#combine all explanatory variables

water_reg <- syncsa_summary_list$single_actualw.biomass.log.axis %>% 
  rownames_to_column(var = "sites") %>% 
  left_join(wat_covariates)

m1w <- glm(log(TCAP.p) ~ n +mean_logwater+ +BC2 + BC4 + BC15 + BC17 + chao + andes, family = gaussian, data = water_reg)

Anova(m1w, type = 2) #negative effect of n and andes north on log p (i.e. n and andes north associted with significance), 
#positive effect of cv env on log p (i.e. cv detritus associated with non-sig)

m2w <- glm(TCAP.rho ~ n + mean_logwater+BC2 + BC4 + BC15 + BC17 + chao + andes, family = gaussian, data = water_reg)

Anova(m2w, type = 2) #NEGATIVE effect of andes north on rho
summary(m2w)

# -----------------------------------------------

#figure of TCAP.rho (water) ~ BC15

# general theme for ggplot2 

general_theme <- theme(panel.border = element_blank(),
                       panel.background = element_blank(),
                       plot.title=element_text(size=12),
                       axis.title = element_text(face = "bold", size = 14, colour = "black"),
                       axis.text = element_text(size = 10, colour = "black"), 
                       axis.line.x = element_line(colour = "black", size = 0.25),
                       axis.line.y = element_line(colour = "black", size = 0.25),
                       panel.grid.major = element_blank(),
                       #panel.grid.major = element_line(colour = "grey70", size = 0.2),
                       legend.background = element_blank(),
                       legend.key= element_blank(),
                       strip.background = element_rect(colour="white", fill="white"),
                       strip.text.x = element_text(size = 10),
                       axis.ticks = element_line(colour = "grey70", size = 0.2),
                       axis.title.x = element_text(margin = margin(10, 0, 0, 0)),
                       axis.title.y = element_text(margin = margin(0, 10, 0, 0)))

plot1<-ggplot(water_reg, aes(x=BC15, y=TCAP.rho)) + 
  geom_point(size=3, alpha = 0.75, shape = 19) +
  geom_smooth(method = "lm", se = TRUE, size = 1)+
  labs(y = "TCAP Procrustes correlation",
       x = "Precipitation annual seasonality") +
  general_theme

plot1

ggsave(filename = "figures/water_reg.pdf",
       plot = plot1,
       width = 183, height = 120, units = "mm", dpi = 300) #font problem

ggsave(filename = "figures/water_reg.jpeg",
       plot = plot1,
       width = 183, height = 120, units = "mm", dpi = 600)

