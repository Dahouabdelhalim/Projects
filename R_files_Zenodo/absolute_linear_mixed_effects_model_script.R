### linear mixed-effect models - effect of land use on absolute values

# Linear mixed-effects modelling works with unordered data (sample_data_unordered).

library(lme4)
library(lmerTest)
library(multcomp)
library(multcompView)
library(emmeans)
#--------------------------


#### CUE ####

# Test if influence of landuse on CUE differs depending on site (treatment is repetition)
mixed.lmer.1 <- lmer(CUE ~ landuse + site + (1|treatment), data = sample_data_unordered, REML = F)
mixed.lmer.2 <- lmer(CUE ~ landuse * site + (1|treatment), data = sample_data_unordered, REML = F)
anova(mixed.lmer.1,mixed.lmer.2)
# significant interactive effect ***

mixed.lmer <- lmer(CUE ~ site * landuse + (1|treatment), data = sample_data_unordered)
aov_CUE_interactive <- anova(mixed.lmer)
aov_CUE_interactive
summary(mixed.lmer)
plot(mixed.lmer)
lmer.emm <- emmeans(mixed.lmer, pairwise ~ landuse|site, type = "response")
plot(lmer.emm, comparisons = TRUE)
letters.CUE.interactive <- multcomp::cld(lmer.emm, alpha = 0.05, Letters = letters, reversed = TRUE)
letters.CUE.interactive$group <- gsub(" ", "", letters.CUE.interactive$.group)
letters.CUE.interactive
# --> landuse effect varies between sites, so there is no landuse specific effect per se on CUE

# Is there a significant difference in CUE between land uses over all three sites?
# Test if influence of landuse on CUE (site, treatment is repetition)
mixed.lmer <- lmer(CUE ~ landuse + (1|site) + (1|treatment), data = sample_data_unordered)
aov_CUE_lu <- anova(mixed.lmer)
aov_CUE_lu
summary(mixed.lmer)
plot(mixed.lmer)
lmer.emm <- emmeans(mixed.lmer, ~ landuse)
plot(lmer.emm, comparisons = TRUE)
letters.CUE.lu <- multcomp::cld(lmer.emm, alpha = 0.05, Letters = letters, reversed = TRUE)
letters.CUE.lu$group <- gsub(" ", "", letters.CUE.lu$.group)
letters.CUE.lu
# --> no lu effect over all three sites 



#### cumulative respiration ####

# Test if influence of landuse on cumulresp differs depending on site (treatment is repetition)
mixed.lmer.1 <- lmer(log(cumul_resp_incub_ngC_g_total_time) ~ site + landuse + (1|treatment), data = sample_data_unordered, REML = F)
mixed.lmer.2 <- lmer(log(cumul_resp_incub_ngC_g_total_time) ~ site * landuse + (1|treatment), data = sample_data_unordered, REML = F)
anova(mixed.lmer.1,mixed.lmer.2)
# significant interactive effect *** (#choose the model with lower AIC)

mixed.lmer <- lmer(log(cumul_resp_incub_ngC_g_total_time) ~ site * landuse + (1|treatment), data = sample_data_unordered)
aov_cumulresp_interactive <- anova(mixed.lmer)
aov_cumulresp_interactive
summary(mixed.lmer)
plot(mixed.lmer)
lmer.emm <- emmeans(mixed.lmer, ~ landuse|site)
plot(lmer.emm, comparisons = TRUE)
letters.cumulresp.interactive <- multcomp::cld(lmer.emm, alpha = 0.05, Letters = letters, reversed = TRUE)
letters.cumulresp.interactive$group <- gsub(" ", "", letters.cumulresp.interactive$.group)
letters.cumulresp.interactive
# --> as seen before there is a significant interaction between landuse and site


# Test if influence of landuse on cumulresp Corg differs depending on site (treatment is repetition)
mixed.lmer <- lmer(log(cumul_resp_Corg) ~ site * landuse + (1|treatment), data = sample_data_unordered)
aov_cumulrespCorg_interactive <- anova(mixed.lmer)
aov_cumulrespCorg_interactive
summary(mixed.lmer)
plot(mixed.lmer)
lmer.emm <- emmeans(mixed.lmer, ~ landuse|site)
plot(lmer.emm, comparisons = TRUE)
letters.cumulrespCorg.interactive <- multcomp::cld(lmer.emm, alpha = 0.05, Letters = letters, reversed = TRUE)
letters.cumulrespCorg.interactive$group <- gsub(" ", "", letters.cumulrespCorg.interactive$.group)
letters.cumulrespCorg.interactive
# --> as seen before there is a significant interaction between landuse and site


# Is there a significant difference in cumulresp between land uses over all three sites?
# test if influence of landuse on cumulresp (site, treatment is repetition)
mixed.lmer <- lmer(log(cumul_resp_incub_ngC_g_total_time) ~ landuse + (1|site) + (1|treatment), data = sample_data_unordered)
aov_cumulresp_lu <- anova(mixed.lmer)
aov_cumulresp_lu
summary(mixed.lmer)
plot(mixed.lmer)
lmer.emm <- emmeans(mixed.lmer, ~ landuse)
plot(lmer.emm, comparisons = TRUE)
letters.cumulresp.lu <- multcomp::cld(lmer.emm, alpha = 0.05, Letters = letters, reversed = TRUE)
letters.cumulresp.lu$group <- gsub(" ", "", letters.cumulresp.lu$.group)
#letters.cumulresp.lu
# --> general lu effect over all three sites (forest highest cumulresp)



#### Cmic/Corg ####

# Test if Cmic over Corg differs depending on site and landuse
aov_CmicCorg_interactive <- aov(Cmic_Corg_perc ~ site * landuse, data = gensoil)
summary(aov_CmicCorg_interactive)
lmer.emm <- emmeans(aov_CmicCorg_interactive, ~ landuse|site)
plot(lmer.emm, comparisons = TRUE)
letters.CmicCorg.interactive <- multcomp::cld(lmer.emm, alpha = 0.05, Letters = letters, reversed = TRUE)
letters.CmicCorg.interactive$group <- gsub(" ", "", letters.CmicCorg.interactive$.group)
letters.CmicCorg.interactive
# --> only dependent on site, no interaction


#### Cgrowth ####

# Is there a significant difference in Cgrowth between land uses over all three sites?
mixed.lmer <- lmer(log(Cgrowth_ngC_gDW_h) ~ site * landuse + (1|treatment), data = sample_data_unordered)
aov_Cgrowth_interactive <- anova(mixed.lmer)
aov_Cgrowth_interactive
summary(mixed.lmer)
plot(mixed.lmer)
lmer.emm <- emmeans(mixed.lmer, ~ landuse|site)
plot(lmer.emm, comparisons = TRUE)
letters.Cgrowth.interactive <- multcomp::cld(lmer.emm, alpha = 0.05, Letters = letters, reversed = TRUE)
letters.Cgrowth.interactive$group <- gsub(" ", "", letters.Cgrowth.interactive$.group)
letters.Cgrowth.interactive
# --> there is a significant interaction between landuse and site