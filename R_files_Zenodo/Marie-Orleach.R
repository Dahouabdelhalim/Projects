## Main R script associated to the manuscript entitled 
## "The repeatable opportunity for selection differs between pre- and postcopulatory fitness components" 
## by L. Marie-Orleach, N. Vellnow, and L. Schärer

library(rptR) #for rptR::rpt
library(hablar) #for hablar::s
library(lmerTest) #for lmerTest::lmer

source(".../Marie-Orleach.et.al_R.SCRIPT_FUNCTIONS.R")
iteration <- 10000

#### 1. variance and covariance in mRS and all four fitness components ####
dataset_focal <- read.delim(".../Marie-Orleach.et.al_DATASET.FOCAL.txt")
dataset_focal <- FITCOMP (dataset_focal)

#### 1.1 standard variance in mRS and all four fitness components ####
stdvar_mRS <- STDVAR ("mRS", dataset_focal)
stdvar_F   <- STDVAR ("F"  , dataset_focal)
stdvar_MS  <- STDVAR ("MS" , dataset_focal)
stdvar_STE <- STDVAR ("STE", dataset_focal)
stdvar_SFE <- STDVAR ("SFE", dataset_focal)
 
#### 1.2 binomial sampling error in STE and SFE ####
bse_STE <- BSE ("STE", dataset_focal)
bse_SFE <- BSE ("SFE", dataset_focal)

#### 1.3 covariances ####
stdcov_F.MS    <- STDCOV ("F",   "MS",  dataset_focal)
stdcov_F.STE   <- STDCOV ("F",   "STE", dataset_focal)
stdcov_F.SFE   <- STDCOV ("F",   "SFE", dataset_focal)
stdcov_MS.STE  <- STDCOV ("MS",  "STE", dataset_focal)
stdcov_MS.SFE  <- STDCOV ("MS",  "SFE", dataset_focal)
stdcov_STE.SFE <- STDCOV ("STE", "SFE", dataset_focal)

#### 1.4 total variance ####
totvar <- TOTVAR (dataset_focal)

#### 1.5 variance and covariance bootstrap (95% CIs) ####
boot_varcovar_outcome <- BOOT_VARCOVAR (dataset_focal, iteration)

#### 1.6 variance pairwise comparisons (P values) ####
pvar_F.MS    <- PWCOMP_VAR ("F",   "MS",  dataset_focal, iteration)
pvar_F.STE   <- PWCOMP_VAR ("F",   "STE", dataset_focal, iteration)
pvar_F.SFE   <- PWCOMP_VAR ("F",   "SFE", dataset_focal, iteration)
pvar_MS.STE  <- PWCOMP_VAR ("MS",  "STE", dataset_focal, iteration)
pvar_MS.SFE  <- PWCOMP_VAR ("MS",  "SFE", dataset_focal, iteration)
pvar_STE.SFE <- PWCOMP_VAR ("STE", "SFE", dataset_focal, iteration)



#### 2. repeatability in mRS and all four fitness components ####
dataset_group <- read.delim(".../Marie-Orleach.et.al_DATASET.MATING.GROUP.txt")
dataset_group <- FITCOMP (dataset_group) 

## /!\\ warnings 'boundary (singular) fit: see ?isSingular' arise when the random effect (1|focal) explains 0 variance.
## /!\\ These warnings arise for F, and for the permutation tests of all fitness components (which is expected).
## /!\\ These warnings appear in all steps in which repeatability is assessed (3.2 & 3.3)
rpt_mRS <- RPT ("mRS", TRANS_mRS, dataset_group, iteration) 
rpt_F   <- RPT ("F",   TRANS_F,   dataset_group, iteration)
rpt_MS  <- RPT ("MS",  TRANS_MS,  dataset_group, iteration)
rpt_STE <- RPT ("STE", TRANS_STE, dataset_group, iteration)
rpt_SFE <- RPT ("SFE", TRANS_SFE, dataset_group, iteration)



#### 3. repeatable variance in mRS and all four fitness components ####
dataset_focal <- read.delim(".../Marie-Orleach.et.al_DATASET.FOCAL.txt")
dataset_focal <- FITCOMP (dataset_focal)
dataset_group <- read.delim(".../Marie-Orleach.et.al_DATASET.MATING.GROUP.txt")
dataset_group <- FITCOMP (dataset_group) 

#### 3.1 repeatable variance ####
## /!\\ missing data are explained in the Methods of the article
rptvar_mRS <- RPTVAR ("mRS", TRANS_mRS, dataset_focal, dataset_group)
rptvar_F   <- RPTVAR ("F",   TRANS_F,   dataset_focal, dataset_group)
rptvar_MS  <- RPTVAR ("MS",  TRANS_MS,  dataset_focal, dataset_group)
rptvar_STE <- RPTVAR ("STE", TRANS_STE, dataset_focal, dataset_group)
rptvar_SFE <- RPTVAR ("SFE", TRANS_SFE, dataset_focal, dataset_group)

#### 3.2 repeatable variance bootstrap (95% CIs) ####
boot_rptvar_outcome <- BOOT_RPTVAR (dataset_focal, dataset_group, iteration)

#### 3.3 pairwise comparisons in repeatable variance ####
prptvar_F.MS    <- PWCOMP_RPTVAR ("F",   "MS",  TRANS_F,   TRANS_MS,  dataset_focal, dataset_group, iteration)
prptvar_F.STE   <- PWCOMP_RPTVAR ("F",   "STE", TRANS_F,   TRANS_STE, dataset_focal, dataset_group, iteration)
prptvar_F.SFE   <- PWCOMP_RPTVAR ("F",   "SFE", TRANS_F,   TRANS_SFE, dataset_focal, dataset_group, iteration)
prptvar_MS.STE  <- PWCOMP_RPTVAR ("MS",  "STE", TRANS_MS,  TRANS_STE, dataset_focal, dataset_group, iteration)
prptvar_MS.SFE  <- PWCOMP_RPTVAR ("MS",  "SFE", TRANS_MS,  TRANS_SFE, dataset_focal, dataset_group, iteration)
prptvar_STE.SFE <- PWCOMP_RPTVAR ("STE", "SFE", TRANS_STE, TRANS_SFE, dataset_focal, dataset_group, iteration)



#### 4. group and batch effects ####
dataset_group <- read.delim(".../Marie-Orleach.et.al_DATASET.MATING.GROUP.txt")
dataset_group <- FITCOMP (dataset_group) 

m_mRS <- lmer(REL(TRANS_mRS(mRS)) ~  (1|focal) + as.factor(mating_group) + as.factor(batch) + as.factor(mating_group)*as.factor(batch), data=dataset_group)
m_F   <- lmer(REL(TRANS_F(F))     ~  (1|focal) + as.factor(mating_group) + as.factor(batch) + as.factor(mating_group)*as.factor(batch), data=dataset_group)
m_MS  <- lmer(REL(TRANS_MS(MS))   ~  (1|focal) + as.factor(mating_group) + as.factor(batch) + as.factor(mating_group)*as.factor(batch), data=dataset_group)
m_STE <- lmer(REL(TRANS_STE(STE)) ~  (1|focal) + as.factor(mating_group) + as.factor(batch) + as.factor(mating_group)*as.factor(batch), data=dataset_group)
m_SFE <- lmer(REL(TRANS_SFE(SFE)) ~  (1|focal) + as.factor(mating_group) + as.factor(batch) + as.factor(mating_group)*as.factor(batch), na.action=na.omit, data=dataset_group)

anova(m_mRS)
anova(m_F)
anova(m_MS)
anova(m_STE)
anova(m_SFE)