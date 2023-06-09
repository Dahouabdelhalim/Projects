### RRN and RRT over time - 50d incubation

# This script works with unordered data (incub_resp_unordered).
# It creates Figure S4.

library(lme4)
library(lmerTest)
library(multcomp)
library(multcompView)
library(emmeans)

library(tidyverse)
library(RColorBrewer)
library(ggpmisc)
library(ggpubr)
#--------------------------

# load script, which loads needed data
source("load_data_script.R")

#### trends per treatment over incubation - incub_resp ####

model.resp.rate.1 <- lmer(respiration_ngC_g_h_n_at_20.C ~ treatment * d_since_start_num + (1|sample_ID), incub_resp_unordered, REML = FALSE)
model.resp.rate.2 <- lmer(respiration_ngC_g_h_n_at_20.C ~ treatment + d_since_start_num + (1|sample_ID), incub_resp_unordered, REML = FALSE)
anova(model.resp.rate.1, model.resp.rate.2)
# only sign. difference if log-transformed, and than only slightly lower AIC in model.resp.rate.1

model.resp.rate <- lmer(log(cumul_resp_total_time) ~ landuse + treatment + d_since_start_num + (1|site) + (1|sample_ID), incub_resp_unordered)
summary(model.resp.rate)
anova(model.resp.rate)
plot(model.resp.rate)
plot(predict(model.resp.rate))

lmer.emm <- emmeans(model.resp.rate,~ treatment|landuse)
plot(lmer.emm, comparisons = TRUE)
multcomp::cld(lmer.emm, alpha = 0.05, Letters = letters, reversed = TRUE)


#### calculate RRT and RRN over time ####

T10_incub_resp_unordered <- incub_resp_unordered[incub_resp_unordered$treatment=="T10",c(2:4,21,38,40)]
colnames(T10_incub_resp_unordered) <- c("sample_ID", "site", "landuse", "d_of_incub", "T10_resp_rate", "T10_cumul_resp")
T20_incub_resp_unordered <- incub_resp_unordered[incub_resp_unordered$treatment=="T20",c(2:4,21,38,40)]
colnames(T20_incub_resp_unordered) <- c("sample_ID", "site", "landuse", "d_of_incub", "T20_resp_rate", "T20_cumul_resp")
T20N_incub_resp_unordered <- incub_resp_unordered[incub_resp_unordered$treatment=="T20+N",c(2:4,21,38,40)]
colnames(T20N_incub_resp_unordered) <- c("sample_ID", "site", "landuse", "d_of_incub", "T20+N_resp_rate", "T20+N_cumul_resp")

RRT_incub_resp_unordered <- merge(T10_incub_resp_unordered, T20_incub_resp_unordered, by=1:4)
RRTN_incub_resp_unordered <- merge(T20N_incub_resp_unordered, RRT_incub_resp_unordered, by=1:4)


RRTN_incub_resp_unordered$RRT_cumul_resp <- RRTN_incub_resp_unordered$T20_cumul_resp/RRTN_incub_resp_unordered$T10_cumul_resp
RRTN_incub_resp_unordered$RRN_cumul_resp <- RRTN_incub_resp_unordered$`T20+N_cumul_resp`/RRTN_incub_resp_unordered$T20_cumul_resp
RRTN_incub_resp_unordered$RRT_resp_rate <- RRTN_incub_resp_unordered$T20_resp_rate/RRTN_incub_resp_unordered$T10_resp_rate
RRTN_incub_resp_unordered$RRN_resp_rate <- RRTN_incub_resp_unordered$`T20+N_resp_rate`/RRTN_incub_resp_unordered$T20_resp_rate


#### create mean RRT and RRN of 3 replicates ####

RRNmeanresp <- RRTN_incub_resp_unordered %>%
  group_by(site, landuse, d_of_incub) %>%
  summarise(mean_RR_cumul_resp = mean(RRN_cumul_resp),
            sd_RR_cumul_resp = sd(RRN_cumul_resp),
            mean_RR_resp_rate = mean(RRN_resp_rate),
            sd_RR_resp_rate = sd(RRN_resp_rate)) %>%
  ungroup()
RRNmeanresp$RR <- c("N")

RRTmeanresp <- RRTN_incub_resp_unordered %>%
  group_by(site, landuse, d_of_incub) %>%
  summarise(mean_RR_cumul_resp = mean(RRT_cumul_resp),
            sd_RR_cumul_resp = sd(RRT_cumul_resp),
            mean_RR_resp_rate = mean(RRT_resp_rate),
            sd_RR_resp_rate = sd(RRT_resp_rate)) %>%
  ungroup()
RRTmeanresp$RR <- c("T")

RRTNmeanresp <- rbind(RRTmeanresp, RRNmeanresp)
RRTNmeanresp$landuse <- ordered(RRTNmeanresp$landuse, levels = c("Forest", "Grassland", "Cropland"))


#### Figure S4 - RRT of respiration rate over time ####

ggplot(data=RRTmeanresp, aes(x=as.numeric(d_of_incub), y=mean_RR_resp_rate, colour = landuse)) +
  geom_point(aes(shape = landuse), size = 2) +
  theme_bw(base_size = 13) + 
  theme(strip.text.y = element_text(angle = 45)) +
  theme(strip.background = element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "bottom") +
  xlab(expression(days~of~incubation)) +
  ylab(expression(RR[T]~respiration~rate)) +
  scale_color_manual(values = c("#d3b023", "#8F7700FF","#EFC000FF")) + # landuse yellow
  scale_shape_manual(values = c(15,17,16)) +
  geom_hline(yintercept = 1, colour = "#999999", linetype = "dashed") +
  facet_grid(site~., scales = "free")



#### plot RRN of respiration rate over time ####

ggplot(data=RRNmeanresp, aes(x=as.numeric(d_of_incub), y=mean_RR_resp_rate, colour = landuse)) +
  geom_point(aes(shape = landuse), size = 2) +
  theme_bw(base_size = 13) + 
  theme(strip.text.y = element_text(angle = 45)) +
  theme(strip.background = element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "bottom") +
  xlab(expression(days~of~incubation)) +
  ylab(expression(RR[N]~respiration~rate)) +
  scale_color_manual(values = c("#0073C2FF", "#7AA6DCFF","#003C67FF")) + # landuse blue
  scale_shape_manual(values = c(15,17,16)) +
  geom_hline(yintercept = 1, colour = "#999999", linetype = "dashed") +
  facet_grid(site~., scales = "free")
