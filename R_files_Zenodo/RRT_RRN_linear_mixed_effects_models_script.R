### linear mixed-effect models - effect of land use on response ratios

# Linear mixed-effects modelling works with unordered data.

library(lme4)
library(lmerTest)
library(multcomp)
library(multcompView)
library(emmeans)
#--------------------------

# load script, which loads needed data
source("load_data_script.R")

# load script, which loads needed data
source("RRT_RRN_calculation_script.R")



#### cumul resp ####

# linear mixed-effect model RRT - cumulresp
mixed.lmer <- lmer(cumul_resp_incub_ngC_g_total_time ~ landuse + (1|site), data = multipleRRT)
aov.RRT_cumulresp <- anova(mixed.lmer)
aov.RRT_cumulresp
summary(mixed.lmer)
plot(mixed.lmer)
lmer.emm <- emmeans(mixed.lmer, ~ landuse)
plot(lmer.emm, comparisons = TRUE)
letters.cld.cumulresp.RRT <- multcomp::cld(lmer.emm, alpha = 0.05, Letters = letters, reversed = TRUE)
letters.cld.cumulresp.RRT$group <- gsub(" ", "", letters.cld.cumulresp.RRT$.group)
letters.cld.cumulresp.RRT$RR <- c("T")

# linear mixed-effect model RRN - cumulresp
mixed.lmer <- lmer(cumul_resp_incub_ngC_g_total_time ~ landuse + (1|site), data = multipleRRN)
aov.RRN_cumulresp <- anova(mixed.lmer)
aov.RRN_cumulresp
summary(mixed.lmer)
plot(mixed.lmer)
lmer.emm <- emmeans(mixed.lmer, ~ landuse)
plot(lmer.emm, comparisons = TRUE)
letters.cld.cumulresp.RRN <- multcomp::cld(lmer.emm, alpha = 0.05, Letters = letters, reversed = TRUE)
letters.cld.cumulresp.RRN$group <- gsub(" ", "", letters.cld.cumulresp.RRN$.group)
letters.cld.cumulresp.RRN$RR <- c("N")

letters.cld.cumulresp <- rbind(letters.cld.cumulresp.RRT, letters.cld.cumulresp.RRN)
letters.cld.cumulresp



#### cumul resp Corg ####

# linear mixed-effect model RRT - cumulresp
mixed.lmer <- lmer(cumul_resp_Corg ~ landuse + (1|site), data = multipleRRT)
aov.RRT_cumulrespCorg <- anova(mixed.lmer)
aov.RRT_cumulrespCorg
summary(mixed.lmer)
plot(mixed.lmer)
lmer.emm <- emmeans(mixed.lmer, ~ landuse)
plot(lmer.emm, comparisons = TRUE)
letters.cld.cumulrespCorg.RRT <- multcomp::cld(lmer.emm, alpha = 0.05, Letters = letters, reversed = TRUE)
letters.cld.cumulrespCorg.RRT$group <- gsub(" ", "", letters.cld.cumulrespCorg.RRT$.group)
letters.cld.cumulrespCorg.RRT$RR <- c("T")

# linear mixed-effect model RRN - cumulresp
mixed.lmer <- lmer(cumul_resp_Corg ~ landuse + (1|site), data = multipleRRN)
aov.RRN_cumulrespCorg <- anova(mixed.lmer)
aov.RRN_cumulrespCorg
summary(mixed.lmer)
plot(mixed.lmer)
lmer.emm <- emmeans(mixed.lmer, ~ landuse)
plot(lmer.emm, comparisons = TRUE)
letters.cld.cumulrespCorg.RRN <- multcomp::cld(lmer.emm, alpha = 0.05, Letters = letters, reversed = TRUE)
letters.cld.cumulrespCorg.RRN$group <- gsub(" ", "", letters.cld.cumulrespCorg.RRN$.group)
letters.cld.cumulrespCorg.RRN$RR <- c("N")

letters.cld.cumulrespCorg <- rbind(letters.cld.cumulrespCorg.RRT, letters.cld.cumulrespCorg.RRN)
letters.cld.cumulrespCorg


#### Cresp ####

# linear mixed-effect model RRT - Cresp
mixed.lmer <- lmer(Cresp_CUE_ngC_gDW_h ~ landuse + (1|site), data = multipleRRT)
aov.RRT_Cresp <- anova(mixed.lmer)
aov.RRT_Cresp
summary(mixed.lmer)
plot(mixed.lmer)
lmer.emm <- emmeans(mixed.lmer, ~ landuse)
plot(lmer.emm, comparisons = TRUE)
letters.cld.Cresp.RRT <- multcomp::cld(lmer.emm, alpha = 0.05, Letters = letters, reversed = TRUE)
letters.cld.Cresp.RRT$group <- gsub(" ", "", letters.cld.Cresp.RRT$.group)
letters.cld.Cresp.RRT$RR <- c("T")

# linear mixed-effect model RRN - Cresp
mixed.lmer <- lmer(Cresp_CUE_ngC_gDW_h ~ landuse + (1|site), data = multipleRRN)
aov.RRN_Cresp <- anova(mixed.lmer)
aov.RRN_Cresp
summary(mixed.lmer)
plot(mixed.lmer)
lmer.emm <- emmeans(mixed.lmer, ~ landuse)
plot(lmer.emm, comparisons = TRUE)
letters.cld.Cresp.RRN <- multcomp::cld(lmer.emm, alpha = 0.05, Letters = letters, reversed = TRUE)
letters.cld.Cresp.RRN$group <- gsub(" ", "", letters.cld.Cresp.RRN$.group)
letters.cld.Cresp.RRN$RR <- c("N")

letters.cld.Cresp <- rbind(letters.cld.Cresp.RRT, letters.cld.Cresp.RRN)
letters.cld.Cresp



#### CUE ####

# linear mixed-effect model RRT - CUE
mixed.lmer <- lmer(CUE ~ landuse + (1|site), data = multipleRRT)
aov.RRT_CUE <- anova(mixed.lmer)
aov.RRT_CUE
summary(mixed.lmer)
plot(mixed.lmer)
lmer.emm <- emmeans(mixed.lmer, ~ landuse)
plot(lmer.emm, comparisons = TRUE)
letters.cld.CUE.RRT <- multcomp::cld(lmer.emm, alpha = 0.05, Letters = letters, reversed = TRUE)
letters.cld.CUE.RRT$group <- gsub(" ", "", letters.cld.CUE.RRT$.group)
letters.cld.CUE.RRT$RR <- c("T")

# linear mixed-effect model RRN - CUE
mixed.lmer <- lmer(CUE ~ landuse + (1|site), data = multipleRRN)
aov.RRN_CUE <- anova(mixed.lmer)
aov.RRN_CUE
summary(mixed.lmer)
plot(mixed.lmer)
lmer.emm <- emmeans(mixed.lmer, ~ landuse)
plot(lmer.emm, comparisons = TRUE)
letters.cld.CUE.RRN <- multcomp::cld(lmer.emm, alpha = 0.05, Letters = letters, reversed = TRUE)
letters.cld.CUE.RRN$group <- gsub(" ", "", letters.cld.CUE.RRN$.group)
letters.cld.CUE.RRN$RR <- c("N")

letters.cld.CUE <- rbind(letters.cld.CUE.RRT, letters.cld.CUE.RRN)
letters.cld.CUE


#### Cgrowth ####

# linear mixed-effect model RRT - Cgrowth
mixed.lmer <- lmer(Cgrowth_ngC_gDW_h ~ landuse + (1|site), data = multipleRRT)
aov.RRT_Cgrowth <- anova(mixed.lmer)
aov.RRT_Cgrowth
summary(mixed.lmer)
plot(mixed.lmer)
lmer.emm <- emmeans(mixed.lmer, ~ landuse)
plot(lmer.emm, comparisons = TRUE)
letters.cld.Cgrowth.RRT <- multcomp::cld(lmer.emm, alpha = 0.05, Letters = letters, reversed = TRUE)
letters.cld.Cgrowth.RRT$group <- gsub(" ", "", letters.cld.Cgrowth.RRT$.group)
letters.cld.Cgrowth.RRT$RR <- c("T")

# linear mixed-effect model RRN - Cgrowth
#mu <- multipleRRN %>% filter(site!="NO")

mixed.lmer <- lmer(Cgrowth_ngC_gDW_h ~ landuse + (1|site), data = multipleRRN)
aov.RRN_Cgrowth <- anova(mixed.lmer)
aov.RRN_Cgrowth
summary(mixed.lmer)
plot(mixed.lmer)
lmer.emm <- emmeans(mixed.lmer, ~ landuse)
plot(lmer.emm, comparisons = TRUE)
letters.cld.Cgrowth.RRN <- multcomp::cld(lmer.emm, alpha = 0.05, Letters = letters, reversed = TRUE)
letters.cld.Cgrowth.RRN$group <- gsub(" ", "", letters.cld.Cgrowth.RRN$.group)
letters.cld.Cgrowth.RRN$RR <- c("N")

letters.cld.Cgrowth <- rbind(letters.cld.Cgrowth.RRT, letters.cld.Cgrowth.RRN)
letters.cld.Cgrowth



#### Cmic at the end of incubation ####

# linear mixed-effect model RRT - Cmicend
mixed.lmer <- lmer(Cmic_ugC_gDW ~ landuse + (1|site), data = multipleRRT)
aov.RRT_Cmicend <- anova(mixed.lmer)
aov.RRT_Cmicend
summary(mixed.lmer)
plot(mixed.lmer)
lmer.emm <- emmeans(mixed.lmer, ~ landuse)
plot(lmer.emm, comparisons = TRUE)
letters.cld.Cmicend.RRT <- multcomp::cld(lmer.emm, alpha = 0.05, Letters = letters, reversed = TRUE)
letters.cld.Cmicend.RRT$group <- gsub(" ", "", letters.cld.Cmicend.RRT$.group)
letters.cld.Cmicend.RRT$RR <- c("T")

# linear mixed-effect model RRN - Cmicend
mixed.lmer <- lmer(Cmic_ugC_gDW ~ landuse + (1|site), data = multipleRRN)
aov.RRN_Cmicend <- anova(mixed.lmer)
aov.RRN_Cmicend
summary(mixed.lmer)
plot(mixed.lmer)
lmer.emm <- emmeans(mixed.lmer, ~ landuse)
plot(lmer.emm, comparisons = TRUE)
letters.cld.Cmicend.RRN <- multcomp::cld(lmer.emm, alpha = 0.05, Letters = letters, reversed = TRUE)
letters.cld.Cmicend.RRN$group <- gsub(" ", "", letters.cld.Cmicend.RRN$.group)
letters.cld.Cmicend.RRN$RR <- c("N")

letters.cld.Cmicend <- rbind(letters.cld.Cmicend.RRT, letters.cld.Cmicend.RRN)
letters.cld.Cmicend



#### Nmic at the end of incubation ####   

# linear mixed-effect model RRT - Nmicend
mixed.lmer <- lmer(Nmic_ugN_gDW ~ landuse + (1|site), data = multipleRRT)
aov.RRT_Nmicend <- anova(mixed.lmer)
aov.RRT_Nmicend
summary(mixed.lmer)
plot(mixed.lmer)
lmer.emm <- emmeans(mixed.lmer, ~ landuse)
plot(lmer.emm, comparisons = TRUE)
letters.cld.Nmicend.RRT <- multcomp::cld(lmer.emm, alpha = 0.05, Letters = letters, reversed = TRUE)
letters.cld.Nmicend.RRT$group <- gsub(" ", "", letters.cld.Nmicend.RRT$.group)
letters.cld.Nmicend.RRT$RR <- c("T")

# linear mixed-effect model RRN - Nmicend
mixed.lmer <- lmer(Nmic_ugN_gDW ~ landuse + (1|site), data = multipleRRN)
aov.RRN_Nmicend <- anova(mixed.lmer)
aov.RRN_Nmicend
summary(mixed.lmer)
plot(mixed.lmer)
lmer.emm <- emmeans(mixed.lmer, ~ landuse)
plot(lmer.emm, comparisons = TRUE)
letters.cld.Nmicend.RRN <- multcomp::cld(lmer.emm, alpha = 0.05, Letters = letters, reversed = TRUE)
letters.cld.Nmicend.RRN$group <- gsub(" ", "", letters.cld.Nmicend.RRN$.group)
letters.cld.Nmicend.RRN$RR <- c("N")

letters.cld.Nmicend <- rbind(letters.cld.Nmicend.RRT, letters.cld.Nmicend.RRN)
letters.cld.Nmicend



# linear mixed-effect model RRN - Nmicend - test if response driven by CD
multipleRRN.wo.CD <- multipleRRN[multipleRRN$site!="CD",]

mixed.lmer <- lmer(Nmic_ugN_gDW ~ landuse + (1|site), data = multipleRRN.wo.CD)
aov.RRN_Nmicend.wo.CD <- anova(mixed.lmer)
aov.RRN_Nmicend.wo.CD
summary(mixed.lmer)
plot(mixed.lmer)
lmer.emm <- emmeans(mixed.lmer, ~ landuse)
plot(lmer.emm, comparisons = TRUE)
letters.cld.Nmicend.RRN.wo.CD <- multcomp::cld(lmer.emm, alpha = 0.05, Letters = letters, reversed = TRUE)
letters.cld.Nmicend.RRN.wo.CD$group <- gsub(" ", "", letters.cld.Nmicend.RRN.wo.CD$.group)
letters.cld.Nmicend.RRN.wo.CD$RR <- c("N")
letters.cld.Nmicend.RRN.wo.CD
# no significant difference between land use types without CD site --> can't include CD site, 
# as N data in CFE samples from +N treatment are bad

# linear mixed-effect model RRN - Cgrowth - test if response driven by CD
multipleRRN.wo.CD <- multipleRRN[multipleRRN$site!="CD",]

mixed.lmer <- lmer(Cgrowth_ngC_gDW_h ~ landuse + (1|site), data = multipleRRN.wo.CD)
aov.RRN_Cgrowth.wo.CD <- anova(mixed.lmer)
aov.RRN_Cgrowth.wo.CD
summary(mixed.lmer)
plot(mixed.lmer)
lmer.emm <- emmeans(mixed.lmer, ~ landuse)
plot(lmer.emm, comparisons = TRUE)
letters.cld.Cgrowth.RRN.wo.CD <- multcomp::cld(lmer.emm, alpha = 0.05, Letters = letters, reversed = TRUE)
letters.cld.Cgrowth.RRN.wo.CD$group <- gsub(" ", "", letters.cld.Cgrowth.RRN.wo.CD$.group)
letters.cld.Cgrowth.RRN.wo.CD$RR <- c("N")
letters.cld.Cgrowth.RRN.wo.CD
# no significant difference between land use types without CD site --> can't include CD site, 
# as N data in CFE samples from +N treatment are bad


#### mass specific growth rate ####  

# linear mixed-effect model RRT - massspec
mixed.lmer <- lmer(mass_specific_growth_rate_1perd ~ landuse + (1|site), data = multipleRRT)
aov.RRT_massspec <- anova(mixed.lmer)
aov.RRT_massspec
summary(mixed.lmer)
plot(mixed.lmer)
lmer.emm <- emmeans(mixed.lmer, ~ landuse)
plot(lmer.emm, comparisons = TRUE)
letters.cld.massspec.RRT <- multcomp::cld(lmer.emm, alpha = 0.05, Letters = letters, reversed = TRUE)
letters.cld.massspec.RRT$group <- gsub(" ", "", letters.cld.massspec.RRT$.group)
letters.cld.massspec.RRT$RR <- c("T")

# linear mixed-effect model RRN - massspec
mixed.lmer <- lmer(mass_specific_growth_rate_1perd ~ landuse + (1|site), data = multipleRRN)
aov.RRN_massspec <- anova(mixed.lmer)
aov.RRN_massspec
summary(mixed.lmer)
plot(mixed.lmer)
lmer.emm <- emmeans(mixed.lmer, ~ landuse)
plot(lmer.emm, comparisons = TRUE)
letters.cld.massspec.RRN <- multcomp::cld(lmer.emm, alpha = 0.05, Letters = letters, reversed = TRUE)
letters.cld.massspec.RRN$group <- gsub(" ", "", letters.cld.massspec.RRN$.group)
letters.cld.massspec.RRN$RR <- c("N")

letters.cld.massspec <- rbind(letters.cld.massspec.RRT, letters.cld.massspec.RRN)
letters.cld.massspec



#### turnover ####

# linear mixed-effect model RRT - turnover
mixed.lmer <- lmer(turnover_d ~ landuse + (1|site), data = multipleRRT)
aov.RRT_turnover <- anova(mixed.lmer)
aov.RRT_turnover
summary(mixed.lmer)
plot(mixed.lmer)
lmer.emm <- emmeans(mixed.lmer, ~ landuse)
plot(lmer.emm, comparisons = TRUE)
letters.cld.turnover.RRT <- multcomp::cld(lmer.emm, alpha = 0.05, Letters = letters, reversed = TRUE)
letters.cld.turnover.RRT$group <- gsub(" ", "", letters.cld.turnover.RRT$.group)
letters.cld.turnover.RRT$RR <- c("T")

# linear mixed-effect model RRN - turnover
mixed.lmer <- lmer(turnover_d ~ landuse + (1|site), data = multipleRRN)
aov.RRN_turnover <- anova(mixed.lmer)
aov.RRN_turnover
summary(mixed.lmer)
plot(mixed.lmer)
lmer.emm <- emmeans(mixed.lmer, ~ landuse)
plot(lmer.emm, comparisons = TRUE)
letters.cld.turnover.RRN <- multcomp::cld(lmer.emm, alpha = 0.05, Letters = letters, reversed = TRUE)
letters.cld.turnover.RRN$group <- gsub(" ", "", letters.cld.turnover.RRN$.group)
letters.cld.turnover.RRN$RR <- c("N")

letters.cld.turnover <- rbind(letters.cld.turnover.RRT, letters.cld.turnover.RRN)
letters.cld.turnover



#### extractableCnF ####

# linear mixed-effect model RRT - extractableCnF
mixed.lmer <- lmer(nF_ugC_gDW ~ landuse + (1|site), data = multipleRRT)
aov.RRT_extractableCnF <- anova(mixed.lmer)
aov.RRT_extractableCnF
summary(mixed.lmer)
plot(mixed.lmer)
lmer.emm <- emmeans(mixed.lmer, ~ landuse)
plot(lmer.emm, comparisons = TRUE)
letters.cld.extractableCnF.RRT <- multcomp::cld(lmer.emm, alpha = 0.05, Letters = letters, reversed = TRUE)
letters.cld.extractableCnF.RRT$group <- gsub(" ", "", letters.cld.extractableCnF.RRT$.group)
letters.cld.extractableCnF.RRT$RR <- c("T")

# linear mixed-effect model RRN - extractableCnF
mixed.lmer <- lmer(nF_ugC_gDW ~ landuse + (1|site), data = multipleRRN)
aov.RRN_extractableCnF <- anova(mixed.lmer)
aov.RRN_extractableCnF
summary(mixed.lmer)
plot(mixed.lmer)
lmer.emm <- emmeans(mixed.lmer, ~ landuse)
plot(lmer.emm, comparisons = TRUE)
letters.cld.extractableCnF.RRN <- multcomp::cld(lmer.emm, alpha = 0.05, Letters = letters, reversed = TRUE)
letters.cld.extractableCnF.RRN$group <- gsub(" ", "", letters.cld.extractableCnF.RRN$.group)
letters.cld.extractableCnF.RRN$RR <- c("N")

letters.cld.extractableCnF <- rbind(letters.cld.extractableCnF.RRT, letters.cld.extractableCnF.RRN)
letters.cld.extractableCnF



#### extractableNnF ####

# linear mixed-effect model RRT - extractableNnF
mixed.lmer <- lmer(nF_ugN_gDW ~ landuse + (1|site), data = multipleRRT)
aov.RRT_extractableNnF <- anova(mixed.lmer)
aov.RRT_extractableNnF
summary(mixed.lmer)
plot(mixed.lmer)
lmer.emm <- emmeans(mixed.lmer, ~ landuse)
plot(lmer.emm, comparisons = TRUE)
letters.cld.extractableNnF.RRT <- multcomp::cld(lmer.emm, alpha = 0.05, Letters = letters, reversed = TRUE)
letters.cld.extractableNnF.RRT$group <- gsub(" ", "", letters.cld.extractableNnF.RRT$.group)
letters.cld.extractableNnF.RRT$RR <- c("T")

# linear mixed-effect model RRN - extractableNnF
mixed.lmer <- lmer(nF_ugN_gDW ~ landuse + (1|site), data = multipleRRN)
aov.RRN_extractableNnF <- anova(mixed.lmer)
aov.RRN_extractableNnF
summary(mixed.lmer)
plot(mixed.lmer)
lmer.emm <- emmeans(mixed.lmer, ~ landuse)
plot(lmer.emm, comparisons = TRUE)
letters.cld.extractableNnF.RRN <- multcomp::cld(lmer.emm, alpha = 0.05, Letters = letters, reversed = TRUE)
letters.cld.extractableNnF.RRN$group <- gsub(" ", "", letters.cld.extractableNnF.RRN$.group)
letters.cld.extractableNnF.RRN$RR <- c("N")

letters.cld.extractableNnF <- rbind(letters.cld.extractableNnF.RRT, letters.cld.extractableNnF.RRN)
letters.cld.extractableNnF



#### Cmic:Nmic end ####

# linear mixed-effect model RRT - Cmic:Nmic end
mixed.lmer <- lmer(CN_mic ~ landuse + (1|site), data = multipleRRT)
aov.RRT_CNmic <- anova(mixed.lmer)
aov.RRT_CNmic
summary(mixed.lmer)
plot(mixed.lmer)
lmer.emm <- emmeans(mixed.lmer, ~ landuse)
plot(lmer.emm, comparisons = TRUE)
letters.cld.CNmic.RRT <- multcomp::cld(lmer.emm, alpha = 0.05, Letters = letters, reversed = TRUE)
letters.cld.CNmic.RRT$group <- gsub(" ", "", letters.cld.CNmic.RRT$.group)
letters.cld.CNmic.RRT$RR <- c("T")

# linear mixed-effect model RRN - Cmic:Nmic end
mixed.lmer <- lmer(CN_mic ~ landuse + (1|site), data = multipleRRN)
aov.RRN_CNmic <- anova(mixed.lmer)
aov.RRN_CNmic
summary(mixed.lmer)
plot(mixed.lmer)
lmer.emm <- emmeans(mixed.lmer, ~ landuse)
plot(lmer.emm, comparisons = TRUE)
letters.cld.CNmic.RRN <- multcomp::cld(lmer.emm, alpha = 0.05, Letters = letters, reversed = TRUE)
letters.cld.CNmic.RRN$group <- gsub(" ", "", letters.cld.CNmic.RRN$.group)
letters.cld.CNmic.RRN$RR <- c("N")

letters.cld.CNmic <- rbind(letters.cld.CNmic.RRT, letters.cld.CNmic.RRN)
letters.cld.CNmic