library(lme4)
library(lmerTest)
library(emmeans)
library(DHARMa)


### 01. import plant resistance arsenal data -----------------------------------
setwd("~/Desktop/DRYAD/")
dat <- read.csv(file = "02_Rotter_Mimgut_chemistry_data.csv", header = T, stringsAsFactors = F )

# add unique identifier for maternal line
dat$Family <- paste(dat$Population, dat$Family, sep = "-")
dat$Family <- as.factor(dat$Family)

# remove rows (plants) that were not sampled for PPGs
dat <- dat[!is.na(dat$UNK_10), ]
nrow(dat) # 483 plants sampled

# summarize number of samples
table(dat$Population)
table(dat$Family)

# analyze chemistry data
      # PPGs - unk10, calcA, conan, verb, calcB, mimulo, unk16
      # models should include fixed effects for: Subregion and Population
      # and random effects for: Family (maternal line)

names(dat)[7:13] # PPGS of interest

# explore overall variation of PPGs
max(dat$UNK_10) / min(dat$UNK_10[dat$UNK_10 > 0]) # 3339.636
max(dat$VERB) / min(dat$VERB[dat$VERB > 0]) # 3606.153
max(dat$CALC_B) / min(dat$CALC_B[dat$CALC_B > 0]) # 144428.1
max(dat$MIMULO) / min(dat$MIMULO[dat$MIMULO > 0]) # 1676.599
max(dat$UNK_16) / min(dat$UNK_16[dat$UNK_16 > 0]) # 8286.347
max(dat$CALC_A) / min(dat$CALC_A[dat$CALC_A > 0]) # 2211.691
max(dat$CONAN) / min(dat$CONAN[dat$CONAN > 0]) # ~2995.171


# to facilitate log-transformation, add constant to all PPG values -------------
  # smallest observed concentration value for that particular PPG

dat$UNK_10 <- dat$UNK_10 + min(dat$UNK_10[dat$UNK_10 > 0])
dat$VERB <- dat$VERB + min(dat$VERB[dat$VERB > 0])
dat$CALC_B <- dat$CALC_B + min(dat$CALC_B[dat$CALC_B > 0])
dat$MIMULO <- dat$MIMULO + min(dat$MIMULO[dat$MIMULO > 0])
dat$UNK_16 <- dat$UNK_16 + min(dat$UNK_16[dat$UNK_16 > 0])

min(dat$UNK_10[dat$UNK_10 > 0]) # 0.00723909
  # min(dat$CALC_A[dat$CALC_A > 0]) # 0.04301171 (no need, use negative binomial)
  # min(dat$CONAN[dat$CONAN > 0]) # 0.07283934 (no need, use negative binomial)
min(dat$VERB[dat$VERB > 0]) # 0.01494847
min(dat$CALC_B[dat$CALC_B > 0]) # 0.000281202
min(dat$MIMULO[dat$MIMULO > 0]) # 0.004471893
min(dat$UNK_16[dat$UNK_16 > 0]) # 0.004901809


# add column with log-transformed PPG data -------------------------------------
dat$UNK_10_log <- log(dat$UNK_10)
dat$CALC_A_log <- log(dat$CALC_A)
dat$CONAN_log <- log(dat$CONAN)
dat$VERB_log <- log(dat$VERB)
dat$CALC_B_log <- log(dat$CALC_B)
dat$MIMULO_log <- log(dat$MIMULO)
dat$UNK_16_log <- log(dat$UNK_16)


### 03. model PPGs -------------------------------------------------------------
names(dat)[7:13] # PPGS of interest


# UNK_10 -----------------------------------------------------------------------

# view data distribution
hist(dat$UNK_10_log, breaks = 50)

mod_UNK10  <-  lmer(UNK_10_log ~ Subregion + Population + (1|Family), 
      data = dat)

summary(mod_UNK10)
simulateResiduals(fittedModel = mod_UNK10, plot = T)


# CALC_A -----------------------------------------------------------------------

# view data distribution
hist(dat$CALC_A, breaks = 50)

m_nb_CALCA <- glmer.nb(CALC_A ~ Subregion + Population + (1|Family), data = dat)

summary(m_nb_CALCA)
simulateResiduals(fittedModel = m_nb_CALCA, plot = T)


# CONAN ------------------------------------------------------------------------

# view data distribution
hist(dat$CONAN, breaks = 50)

m_nb_CONAN <- glmer.nb(CONAN ~ Subregion + Population + (1|Family), data = dat)

summary(m_nb_CONAN)
simulateResiduals(fittedModel = m_nb_CONAN, plot = T)


# VERB ------------------------------------------------------------------------

# view data distribution
hist(dat$VERB_log, breaks = 50)

mod_VERB  <-  lmer(VERB_log ~ Subregion + Population + (1|Family), 
                    data = dat)

summary(mod_VERB)
simulateResiduals(fittedModel = mod_VERB, plot = T)


# CALC_B -----------------------------------------------------------------------

# view data distribution
hist(dat$CALC_B_log, breaks = 50)

mod_CALCB  <-  lmer(CALC_B_log ~ Subregion + Population + (1|Family), 
                   data = dat)

summary(mod_CALCB)
simulateResiduals(fittedModel = mod_CALCB, plot = T)


# MIMULO -----------------------------------------------------------------------

# view data distribution
hist(dat$MIMULO_log, breaks = 50)

mod_MIMULO  <-  lmer(MIMULO_log ~ Subregion + Population + (1|Family), 
                    data = dat)

summary(mod_MIMULO)
simulateResiduals(fittedModel = mod_MIMULO, plot = T)


# UNK_16 -----------------------------------------------------------------------

# view data distribution
hist(dat$UNK_16_log, breaks = 50)

mod_UNK16  <-  lmer(UNK_16_log ~ Subregion + Population + (1|Family), 
                    data = dat)

summary(mod_UNK16)
simulateResiduals(fittedModel = mod_UNK16, plot = T)



### 04. generate model predictions ---------------------------------------------

UNK10_preds <-  as.data.frame(emmeans(mod_UNK10, spec = "Population", type = "response"))
UNK10_preds$predicted <- exp(UNK10_preds$emmean) - min(dat$UNK_10[dat$UNK_10 > 0])

CALCA_preds <-  as.data.frame(emmeans(m_nb_CALCA, spec = "Population", type = "response"))
CONAN_preds <-  as.data.frame(emmeans(m_nb_CONAN, spec = "Population", type = "response"))

VERB_preds <-  as.data.frame(emmeans(mod_VERB, spec = "Population", type = "response"))
VERB_preds$predicted <- exp(VERB_preds$emmean) - min(dat$VERB[dat$VERB > 0])

CALCB_preds <-  as.data.frame(emmeans(mod_CALCB, spec = "Population", type = "response"))
CALCB_preds$predicted <- exp(CALCB_preds$emmean) - min(dat$CALC_B[dat$CALC_B > 0])

MIMULO_preds <-  as.data.frame(emmeans(mod_MIMULO, spec = "Population", type = "response"))
MIMULO_preds$predicted <- exp(MIMULO_preds$emmean) - min(dat$MIMULO[dat$MIMULO > 0])

UNK16_preds <-  as.data.frame(emmeans(mod_UNK16, spec = "Population", type = "response"))
UNK16_preds$predicted <- exp(UNK16_preds$emmean) - min(dat$UNK_16[dat$UNK_16 > 0])


### 05. save trait data --------------------------------------------------------
trait_dat <- data.frame(
  Population = UNK10_preds$Population,
  UNK_10 = UNK10_preds$predicted,
  CALC_A = CALCA_preds$response,
  CONAN = CONAN_preds$response,
  VERB = VERB_preds$predicted,
  CALC_B = CALCB_preds$predicted,
  MIMULO = MIMULO_preds$predicted,
  UNK_16 = UNK16_preds$predicted)

# write.csv(trait_dat, file = "predicted_trait_values_01_07_2022.csv", row.names = F)

# explore population level variation
max(trait_dat$UNK_10) / min(trait_dat$UNK_10[trait_dat$UNK_10 > 0]) # 21x
max(trait_dat$VERB) / min(trait_dat$VERB[trait_dat$VERB > 0]) # 7x
max(trait_dat$CALC_B) / min(trait_dat$CALC_B[trait_dat$CALC_B > 0]) # 48x
max(trait_dat$MIMULO) / min(trait_dat$MIMULO[trait_dat$MIMULO > 0]) # 10x
max(trait_dat$UNK_16) / min(trait_dat$UNK_16[trait_dat$UNK_16 > 0]) # 20x
max(trait_dat$CALC_A) / min(trait_dat$CALC_A[trait_dat$CALC_A > 0]) # 5 x
max(trait_dat$CONAN) / min(trait_dat$CONAN[trait_dat$CONAN > 0]) # 7x


### 05b. build figure showing average concentration of compounds by region
      # use model predictions --------------------------------------------------

location_dat <- read.csv(file = "01_Rotter_Mimgut_locations.csv"   , header = T, stringsAsFactors = F)
c_dat <- merge(trait_dat, location_dat, by.x = "Population", by.y = "Population")

# format for ggplot
p1 = c_dat[, c("Subregion", "UNK_10")]
p2 = c_dat[, c("Subregion", "CALC_A")]
p3 = c_dat[, c("Subregion", "CONAN")]
p4 = c_dat[, c("Subregion", "VERB")]
p5 = c_dat[, c("Subregion", "CALC_B")]
p6 = c_dat[, c("Subregion", "MIMULO")]
p7 = c_dat[, c("Subregion", "UNK_16")]

names(p1) <- c("Subregion", "Concentration")
names(p2) <- c("Subregion", "Concentration")
names(p3) <- c("Subregion", "Concentration")
names(p4) <- c("Subregion", "Concentration")
names(p5) <- c("Subregion", "Concentration")
names(p6) <- c("Subregion", "Concentration")
names(p7) <- c("Subregion", "Concentration")

all_chem_dat <- rbind(p1,p2,p3,p4,p5,p6,p7)

all_chem_dat$PPG <- rep(c("unknown PPG 10", "calceolarioside A", "conandroside", "verbascoside", "calceolarioside B", "mimuloside", "unknown PPG 16"), each = 41)

my_colors <- c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3")
      # alpha, coastal, cordilleran, ena, northern, southern, uk

library(ggplot2)
dev.new()

ggplot(data = all_chem_dat, aes(x = Subregion, y = Concentration)) +
      geom_boxplot(fill = rep(my_colors, times = 7)) +
      facet_wrap(~PPG, scales = "free") +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
      theme(axis.title.x = element_blank()) +
      labs(y = "PPG concentration") +
      theme(axis.title.y = element_text(face = "bold", size = 14)) 
# quartz.save(file = "Figure_S1_modeled_PPG_conc_across_regions.jpg", type = "jpg", dpi = 300)


### 06. variance partitioning --------------------------------------------------
library(partR2)
library(rptR)

R2_UNK10_mar <- partR2(mod_UNK10, partvars = c("Subregion", "Population"), 
              R2_type = "marginal", nboot = 100)
summary(R2_UNK10_mar)


R2_VERB_mar <- partR2(mod_VERB, partvars = c("Subregion", "Population"), 
                       R2_type = "marginal", nboot = 1000)
summary(R2_VERB_mar)


R2_CALCB_mar <- partR2(mod_CALCB, partvars = c("Subregion", "Population"), 
                      R2_type = "marginal", nboot = 1000)
summary(R2_CALCB_mar)


R2_MIMULO_mar <- partR2(mod_MIMULO, partvars = c("Subregion", "Population"), 
                       R2_type = "marginal", nboot = 1000)
summary(R2_MIMULO_mar)


R2_UNK16_mar <- partR2(mod_UNK16, partvars = c("Subregion", "Population"), 
                        R2_type = "marginal", nboot = 1000)
summary(R2_UNK16_mar)


### random effects variance partitioning in rptR ###

# "In combination with ratio = TRUE, grname = "Fixed" will return 
# the fraction of the total variance that is explained by variance 
# in the linear predictor, i.e. by fixed effects
# Incidently, the combination grname = "Fixed", ratio = TRUE 
# estimates a form of the coefficient of determination R2 
# that we have called the marginal R2 of mixed effects models 
# (Nakagawa & Schielzeth 2013)"

names(dat)[7:13] # PPGS of interest 
  # CALC_A and CONAN used negative binomial models, which don't work with rptR

rpt_UNK10 <- rpt(UNK_10_log ~ Subregion + Population + (1|Family), 
                                      grname = c("Family", "Fixed"), 
                                      data = dat, 
                                      datatype = "Gaussian", 
                                      nboot = 1000, 
                                      npermut = 1000,
                                      adjusted = FALSE,
                                      ratio = TRUE)

rpt_VERB <- rpt(VERB_log ~ Subregion + Population + (1|Family), 
                                      grname = c("Family", "Fixed"), 
                                      data = dat, 
                                      datatype = "Gaussian", 
                                      nboot = 1000, 
                                      npermut = 1000,
                                      adjusted = FALSE,
                                      ratio = TRUE)

rpt_CALCB <- rpt(CALC_B_log ~ Subregion + Population + (1|Family), 
                grname = c("Family", "Fixed"), 
                data = dat, 
                datatype = "Gaussian", 
                nboot = 1000, 
                npermut = 1000,
                adjusted = FALSE,
                ratio = TRUE)

rpt_MIMULO <- rpt(MIMULO_log ~ Subregion + Population + (1|Family), 
                 grname = c("Family", "Fixed"), 
                 data = dat, 
                 datatype = "Gaussian", 
                 nboot = 1000, 
                 npermut = 1000,
                 adjusted = FALSE,
                 ratio = TRUE)

rpt_UNK16 <- rpt(UNK_16_log ~ Subregion + Population + (1|Family), 
                  grname = c("Family", "Fixed"), 
                  data = dat, 
                  datatype = "Gaussian", 
                  nboot = 1000, 
                  npermut = 1000,
                  adjusted = FALSE,
                  ratio = TRUE)



### build forest plot ----------------------------------------------------------
library(ggplot2)

# import .csv containing variance estimates from above 
dat <- read.csv(file = "02_PPG_variance_effects.csv", header = T, stringsAsFactors = F)

# add complete PPG names
dat$PPG2 <- rep(
      c("unknown PPG 10",
        "verbascoside",
        "calceolarioside B",
        "mimuloside",
        "unknown PPG 16"), each = 4)

dev.new(width = 8.5, height = 5, noRStudioGD = TRUE)

ggplot(data = dat, aes(y = term, x = R2, xmin = lower_CI, xmax = upper_CI)) +
      geom_errorbarh(height = 0.25, size = 1.25, col = "darkgray") +
      geom_point(size = 2.5, col = dat$plot_color) + 
      facet_wrap(~ PPG2) +
      labs(y = "") +
      labs(x = "Variance explained") +
      theme(strip.text = element_text(size = 12, face = "bold")) +
      theme(axis.title.x = element_text(size = 12, face = "bold")) +
      scale_y_discrete(labels = c("Maternal Family", "Region + Population", "Population", "Biogeographic Region")) +
      theme_bw() 

# quartz.save(file = "Figure2_variance_partitioning_09_06_2022.jpg", type = "jpg", dpi = 300)


