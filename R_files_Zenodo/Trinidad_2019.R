
#---------------------------------------------------------------------#
# R script associated with the paper "Effect of heterospecific and    #
# conspecific competition on individual differences                   #
# in tadpole behavior"                                                #
# Author: Cammy Beyts                                                 #
# ORCID ID: 0000-0002-4729-2982                                       #
# Email: cammy.beyts@ed.ac.uk                                         #
# Date: First created April 2021.  Last modified 21 Oct 2022          #
#---------------------------------------------------------------------#

#script naming convention for assays and treatment names:
##ACT = activity assay
##EXP = exploration assay
##PRED = predation risk taking assay
##nocomp = no competition treatment
##consp = conspecific treatment
##hetero = heterospecific treatment

####load packages####
library(ggplot2)
library(dplyr) 
library(tidyverse)
library(brms)
library(rstanarm)
library(shinystan)
library(bayesplot)
library(kableExtra)
library(gridExtra)
library(emmeans)
library(tidybayes)
library(coda)

####ggplot themes #####

#ggplot theme for plotting figures
theme_catplots <- function() {
  theme_classic() + # use a white background
    theme(
      axis.title.x = element_text(face = "bold", size =14, colour = "black"),
      axis.title.y = element_text(face = "bold", size = 20, colour = "black"),
      axis.text.y = element_text(face = "bold", size = 14),
      axis.text.x = element_text(face = "bold", size = 14),
      strip.text.x = element_text(face = "bold", size=12),
      strip.background = element_rect(colour=FALSE, fill=FALSE),
      legend.position="none"
    )
}


####custom function to calculate treatment differences ####

#returns mean and 95% CI of treatment difference
treat_compare <- function(treat1, treat2){
  contr <- abs(treat1 - treat2) 
  contr <- round(data.frame(mean(contr), HPDinterval(as.mcmc(contr))), digits = 3) 
  colnames(contr) <- c("Mean", "2.5", "97.5")
  return(contr)
}

####Get data into R####
read.data <- read.delim("Trinidad_2019_ACT_EXP_PRED.txt", header = TRUE)

#copy treatment column and rename
read.data <- read.data %>%
  mutate(NameTreatment=Treatment)

#convert order values 
read.data <- read.data %>%
  mutate(Treatment=recode(Treatment, NoCompetition = 0, Conspecific = 1, Heterospecific= 2)) %>% #treatments converted into 1, 2 and 3
  mutate(Rep=recode(Rep, 0, 1, 2, 3, 4, 5)) %>%  #reps 123456 recoded to 012345 %>%
  mutate(Scale_SVL = scale(SVL))

#set up data types
full.data1 <- read.data %>%
  mutate(Treatment=as.factor(Treatment)) %>%
  mutate(Set=as.factor(Set)) %>% 
  mutate(RepCon=as.numeric(Rep))

#for activity data: change values of 0 into the lowest value recorded for lognormal distribution model 
full.data1$ACT.Dist.Pixels[full.data1$ACT.Dist.Pixels %in% 0] <- 1.41

####brms chains####
n_thin <- 10
burnin <- 1000
pd_size <- 3000
n_chains <- 4
n_iter <- round(pd_size / n_chains) * n_thin + burnin


####brms - body size random slopes model ####
SVL_brms_2022c <- brm(
  SVL ~ Treatment + RepCon + (1|FocalNest) + (1 + RepCon|a|gr(TadpoleID, by = Treatment)),
  data = full.data1,
  family = gaussian(),
  chains = n_chains, iter = n_iter, warmup = burnin, thin = n_thin
)
#save(SVL_brms_2022c, file = "SVL_brms_2022c.rda")
load(file = "SVL_brms_2022c.rda")
summary(SVL_brms_2022c)

#plot of SVL raw values#
svl_raw_plot <- ggplot() +
  geom_boxplot(data = full.data1, aes(x = Treatment, y = SVL)) +
  labs(
    x = "\\n Treatment",
    y = "SVL (mm) \\n"
  ) +
  scale_x_discrete(labels = c("No Competition", "Conspecific", "Heterospecific")) +
  theme_catplots()
svl_raw_plot
ggsave(filename = "brms_tables_2022/Raw_SVL_plot.png", svl_raw_plot, height = 8, width = 8)
  
  
  

#### body size random slopes model - posterior model estimates #####

#fixed effects
fixef_mod_svl <- fixef(SVL_brms_2022c, summary = FALSE)
#variance for tadpole ID
v_id_svl <- (VarCorr(SVL_brms_2022c, summary = FALSE)$TadpoleID$sd)^2
#variances for eggmass ID
v_eggmass_svl <- (VarCorr(SVL_brms_2022c, summary = FALSE)$FocalNest$sd)^2

column_names_svl <- c("Mean", 
                      "2.5", 
                      "97.5"
)

row_names_svl <- c("No competition",
                   "Conspecific",
                   "Heterospecific",
                   "Trial",
                   "Egg Mass ID",
                   "No competition ",
                   "Conspecific ",
                   "Heterospecific ",
                   "No competition  ",
                   "Conspecific  ",
                   "Heterospecific  "
)

#get summary of fixed effects
#population level treatment effects - remove intercept (no competition) from other treatment groups
fixef_mod_SVl_nocomp <- fixef_mod_svl[,1] 
fixef_mod_SVl_consp <- fixef_mod_svl[,1] + fixef_mod_svl[,2] #remove intercept
fixef_mod_SVl_hetero <- fixef_mod_svl[,1] + fixef_mod_svl[,3] #remove intercept
#trial number
fixef_mod_SVL_rep <- fixef_mod_svl[,4]
#combine fixed effects into one df
fixef_mod_SVL <- cbind(fixef_mod_SVl_nocomp, fixef_mod_SVl_consp, fixef_mod_SVl_hetero, fixef_mod_SVL_rep)
fixef_mod_SVL <- round(data.frame(colMeans(fixef_mod_SVL), HPDinterval(as.mcmc(fixef_mod_SVL))), digits = 3) 
colnames(fixef_mod_SVL) <- column_names_svl

#get summary of eggmass
v_eggmass_SVL <- v_eggmass_svl[,1]
v_eggmass_SVL <- round(data.frame(mean(v_eggmass_SVL), HPDinterval(as.mcmc(v_eggmass_SVL))), digits = 3)
colnames(v_eggmass_SVL) <- column_names_svl

#get summary of variances in intial body size (variance in body size in first trial)
v_id_SVL_int <- v_id_svl[,c(1,3,5)]
v_id_SVL_int <- round(data.frame(colMeans(v_id_SVL_int), HPDinterval(as.mcmc(v_id_SVL_int))), digits = 3) 
colnames(v_id_SVL_int) <- column_names_svl

#get summary of variances in change in body size (variance in change in body size from first to last trial)
v_id_SVL_plast <- v_id_svl[,c(2,4,6)]
v_id_SVL_plast <- round(data.frame(colMeans(v_id_SVL_plast), HPDinterval(as.mcmc(v_id_SVL_plast))), digits = 3) 
colnames(v_id_SVL_plast) <- column_names_svl

#combine fixed effects, EggMassID and ID variance into one data frame
SVL_df <- rbind(fixef_mod_SVL, v_eggmass_SVL, v_id_SVL_int, v_id_SVL_plast)
rownames(SVL_df) <- row_names_svl

#kable table of body size results
SVL_kable <- SVL_df  %>%
  knitr::kable(
    format = "html",
    escape = TRUE
  ) %>%
  kable_styling(bootstrap_options = c("hover", "condensed"), full_width = F) %>% 
  add_header_above(c(" " = 2, "95% CI" = 2)) %>%
  add_header_above(c(" " = 1, "Body Size" = 3)) %>%
  pack_rows("Population means", 1, 4) %>%
  pack_rows("Variance among egg masses", 5, 5) %>%
  pack_rows("Variance among individuals in initial body size", 6, 8) %>%
  pack_rows("Variance among individuals in change in body size", 9, 11)
SVL_kable

#### body size random slopes model - treatment comparisons using posterior estimates #####

#get population treatment comparisons
fixef_mod_SVl_nocomp_vs_consp <- treat_compare(fixef_mod_SVl_nocomp, fixef_mod_SVl_consp)
fixef_mod_SVl_nocomp_vs_hetero <- treat_compare(fixef_mod_SVl_nocomp, fixef_mod_SVl_hetero)
fixef_mod_SVl_consp_vs_hetero <- treat_compare(fixef_mod_SVl_consp, fixef_mod_SVl_hetero)
SVL_fe_diff <- rbind(fixef_mod_SVl_nocomp_vs_consp, fixef_mod_SVl_nocomp_vs_hetero, fixef_mod_SVl_consp_vs_hetero)

#get variance among individuals in initial body size - treatment comparisons
v_id_SVL_int_nocomp <- v_id_svl[,c(1)]
v_id_SVL_int_consp <- v_id_svl[,c(3)]
v_id_SVL_int_hetero <- v_id_svl[,c(5)]
v_id_SVL_int_nocomp_vs_consp <- treat_compare(v_id_SVL_int_nocomp, v_id_SVL_int_consp)
v_id_SVL_int_nocomp_vs_hetero <- treat_compare(v_id_SVL_int_nocomp, v_id_SVL_int_hetero)
v_id_SVL_int_consp_vs_hetero <- treat_compare(v_id_SVL_int_consp, v_id_SVL_int_hetero)
v_id_SVL_int_diff <- rbind(v_id_SVL_int_nocomp_vs_consp, v_id_SVL_int_nocomp_vs_hetero, v_id_SVL_int_consp_vs_hetero)

#get variance among individuals in change in body size - treatment comparisons
v_id_SVL_plast_nocomp <- v_id_svl[,c(2)]
v_id_SVL_plast_consp <- v_id_svl[,c(4)]
v_id_SVL_plast_hetero <- v_id_svl[,c(6)]
v_id_SVL_plast_nocomp_vs_consp <- treat_compare(v_id_SVL_plast_nocomp, v_id_SVL_plast_consp)
v_id_SVL_plast_nocomp_vs_hetero <- treat_compare(v_id_SVL_plast_nocomp, v_id_SVL_plast_hetero)
v_id_SVL_plast_consp_vs_hetero <- treat_compare(v_id_SVL_plast_consp, v_id_SVL_plast_hetero)
v_id_SVL_plast_diff <- rbind(v_id_SVL_plast_nocomp_vs_consp, v_id_SVL_plast_nocomp_vs_hetero, v_id_SVL_plast_consp_vs_hetero)

SVL_diff <- rbind(SVL_fe_diff, v_id_SVL_int_diff, v_id_SVL_plast_diff)

rownames(SVL_diff) <- c("No Competition - Conspecific",
                        "No Competition - Heterospecific",
                        "Conspecific - Heterospecific",
                        "No Competition - Conspecific ",
                        "No Competition - Heterospecific ",
                        "Conspecific - Heterospecific ",
                        "No Competition - Conspecific  ",
                        "No Competition - Heterospecific  ",
                        "Conspecific - Heterospecific  "
                        )

#kable table of treatment differences in body size results
SVL_kable_diff <- SVL_diff  %>%
  knitr::kable(
    format = "html",
    escape = TRUE
  ) %>%
  kable_styling(bootstrap_options = c("hover", "condensed"), full_width = F) %>% 
  add_header_above(c(" " = 2, "95% CI" = 2)) %>%
  add_header_above(c(" " = 1, "Body Size" = 3)) %>%
  pack_rows("Population means", 1, 3) %>%
  pack_rows("Variance among individuals in initial body size", 4, 6) %>%
  pack_rows("Variance among individuals in change in body size", 7, 9)
SVL_kable_diff

####Plot SVL results####

#dataframe of treatment effects on tadpole SVL
fixef_mod_SVL_plot <- cbind(fixef_mod_SVl_nocomp, fixef_mod_SVl_consp, fixef_mod_SVl_hetero)
fixef_mod_SVL_plot <- round(data.frame(colMeans(fixef_mod_SVL_plot), HPDinterval(as.mcmc(fixef_mod_SVL_plot))), digits = 3) 
colnames(fixef_mod_SVL_plot) <- column_names_svl

#Effect of treatment on population estimates of tadpole body size
SVL_fe_mean_plot <- 
  ggplot(
    (rbind(fixef_mod_SVL_plot, SVL_fe_diff) %>%
       mutate(
         trait = rep(c("SVL"),each = 6),
         treat = rep(c("6.No Competition", "5.Conspecific", "4.Heterospecific", "3.NC-Cdiff", "2.NC-Hdiff", "1.C-Hdiff"),1)
       ) %>%
       rename(lower = "2.5",
              upper = "97.5")
    ),
    aes(y = treat, color = treat)
  )+
  geom_point(aes(x = Mean), size = 5, alpha = 0.7) +
  geom_linerange(aes(xmin = lower, xmax = upper)) +
  labs(
    title = " ",
    x = "Posterior Estimate with 95% Credible Intervals",
    y = "Treatment" 
  ) +
  scale_y_discrete(labels = c(bquote(Delta ~ "Conspecific-\\n Heterospecific"), bquote(Delta ~ "No Competition-\\n Heterospecific"), bquote(Delta ~ "No Competition-\\n Conspecific"), "Heterospecific", "Conspecific", "No Competition")) +
  scale_colour_manual(values = c("#56B4E9", "#009E73", "#CC79A7", "#D55E00", "#E69F00", "#000000")) +
  geom_vline(xintercept = c(0, 0), linetype = "dotted") +
  theme_catplots() +
  theme(axis.text.x = element_text(face="bold", size=8))
SVL_fe_mean_plot


#Effect of treatment on variance among individuals in initial tadpole body size
SVL_v_id_plot <- 
  ggplot(
    (rbind(v_id_SVL_int, v_id_SVL_int_diff) %>%
       mutate(
         trait = rep(c("SVL"),each = 6),
         treat = rep(c("6.No Competition", "5.Conspecific", "4.Heterospecific", "3.NC-Cdiff", "2.NC-Hdiff", "1.C-Hdiff"),1)
       ) %>%
       rename(lower = "2.5",
              upper = "97.5")
    ),
    aes(y = treat, color = treat)
  )+
  geom_point(aes(x = Mean), size = 5, alpha = 0.7) +
  geom_linerange(aes(xmin = lower, xmax = upper)) +
  labs(
    title = " ",
    x = "Posterior Estimate with 95% Credible Intervals",
    y = "Treatment" 
  ) +
  #facet_grid(. ~ trait, scales = "free") +
  scale_y_discrete(labels = c(bquote(Delta ~ "Conspecific-\\n Heterospecific"), bquote(Delta ~ "No Competition-\\n Heterospecific"), bquote(Delta ~ "No Competition-\\n Conspecific"), "Heterospecific", "Conspecific", "No Competition")) +
  scale_colour_manual(values = c("#56B4E9", "#009E73", "#CC79A7", "#D55E00", "#E69F00", "#000000")) +
  geom_vline(xintercept = c(0, 0), linetype = "dotted") +
  theme_catplots() +
  theme(axis.text.x = element_text(face="bold", size=10))
SVL_v_id_plot

#Effect of treatment on variance among individuals in tadpole growth rates
SVL_v_id_plast_plot <- 
  ggplot(
    (rbind(v_id_SVL_plast, v_id_SVL_plast_diff) %>%
       mutate(
         trait = rep(c("SVL"),each = 6),
         treat = rep(c("6.No Competition", "5.Conspecific", "4.Heterospecific", "3.NC-Cdiff", "2.NC-Hdiff", "1.C-Hdiff"),1)
       ) %>%
       rename(lower = "2.5",
              upper = "97.5")
    ),
    aes(y = treat, color = treat)
  )+
  geom_point(aes(x = Mean), size = 5, alpha = 0.7) +
  geom_linerange(aes(xmin = lower, xmax = upper)) +
  labs(
    title = " ",
    x = "Posterior Estimate with 95% Credible Intervals",
    y = "Treatment" 
  ) +
  #facet_grid(. ~ trait, scales = "free") +
  scale_y_discrete(labels = c(bquote(Delta ~ "Conspecific-\\n Heterospecific"), bquote(Delta ~ "No Competition-\\n Heterospecific"), bquote(Delta ~ "No Competition-\\n Conspecific"), "Heterospecific", "Conspecific", "No Competition")) +
  scale_colour_manual(values = c("#56B4E9", "#009E73", "#CC79A7", "#D55E00", "#E69F00", "#000000")) +
  geom_vline(xintercept = c(0, 0), linetype = "dotted") +
  theme_catplots() +
  theme(axis.text.x = element_text(face="bold", size=10))
SVL_v_id_plast_plot


##### mulitvariate BRMS model ####

#activity assay model
bf_ACT <- bf(
  ACT.Dist.Pixels ~ Treatment + Scale_SVL + RepCon + (1|FocalNest) + (1 |a| gr(TadpoleID, by = Treatment))) +
  lf(sigma ~ Treatment) + lognormal()

#exploration assay model 
bf_EXP <- bf(
  EXP.Dist.Pixels ~ Treatment + Scale_SVL + RepCon + (1|FocalNest) + (1 |a| gr(TadpoleID, by = Treatment))) +
  lf(hu ~ Treatment, sigma ~ Treatment) + hurdle_lognormal()

#predation-risk assay model
bf_PRED <- bf(
  PRED.Dist.Pixels ~ Treatment + Scale_SVL + RepCon + (1|FocalNest) + (1 |a| gr(TadpoleID, by = Treatment))) +
  lf(hu ~ Treatment, sigma ~ Treatment) + hurdle_lognormal()

#combine activity, exploration and predation-risk assay models into one multivariate model
ACT_EXP_PRED_brms_2022a <- brm(
  bf_ACT + bf_EXP + bf_PRED + set_rescor(FALSE), 
  chains = n_chains, iter = n_iter, warmup = burnin, thin = n_thin,
  data = full.data1
)
summary(ACT_EXP_PRED_brms_2022a)
#save(ACT_EXP_PRED_brms_2022a, file = "ACT_EXP_PRED_brms_2022a.rda")
load(file = "ACT_EXP_PRED_brms_2022a.rda")

#plot model outputs
ACT_EXP_PRED_brms_2022a_areas <- mcmc_areas(ACT_EXP_PRED_brms_2022a, 
                                         regex_pars = c("^b_"), prob = 0.80, #use "^sd_" to look at mcmc_areas of variances among individuals, "^b_" for fixed effects "^cor_" for correlations
                                         prob_outer = 0.95, # 95%
                                         point_est = "mean")
ACT_EXP_PRED_brms_2022a_areas


#####get multivariate model results - posterior model estimates####

##select (non treatment) fixed effects
fixef_mod <- fixef(ACT_EXP_PRED_brms_2022a, pars = c("ACTDistPixels_Scale_SVL",
                                                     "ACTDistPixels_RepCon",
                                                     "EXPDistPixels_Scale_SVL",
                                                     "EXPDistPixels_RepCon",
                                                     "PREDDistPixels_Scale_SVL",
                                                     "PREDDistPixels_RepCon"
                                                     ),
                   summary = F)

#select population treatment fixed effects
fixef_mod_treat <- fixef(ACT_EXP_PRED_brms_2022a, pars = c(
                                                     "ACTDistPixels_Intercept",
                                                     "ACTDistPixels_Treatment1",
                                                     "ACTDistPixels_Treatment2",
                                                     "EXPDistPixels_Intercept",
                                                     "EXPDistPixels_Treatment1",
                                                     "EXPDistPixels_Treatment2",
                                                     "PREDDistPixels_Intercept",
                                                     "PREDDistPixels_Treatment1",
                                                     "PREDDistPixels_Treatment2"
                                                     ),
                         summary = F)


##select variance in Egg Mass ID
#posterier estimates are converted from standard deviation to variance estimates
v_eggmass <- (VarCorr(ACT_EXP_PRED_brms_2022a,
                  summary = F)$FocalNest$sd)^2

#select sd among individuals
v_id <- (VarCorr(ACT_EXP_PRED_brms_2022a, summary = F)$TadpoleID$sd)

#select sd within individuals
v_sigma <- fixef(ACT_EXP_PRED_brms_2022a, pars = c("sigma_ACTDistPixels_Intercept", 
                                                           "sigma_EXPDistPixels_Intercept",
                                                           "sigma_PREDDistPixels_Intercept",
                                                           "sigma_ACTDistPixels_Treatment1",
                                                           "sigma_ACTDistPixels_Treatment2",
                                                           "sigma_EXPDistPixels_Treatment1",
                                                           "sigma_EXPDistPixels_Treatment2",
                                                           "sigma_PREDDistPixels_Treatment1",
                                                           "sigma_PREDDistPixels_Treatment2"
                                                           ), summary = F)

#select fixed effects from hurdle model 
fixed_mod_hurdle <- fixef(ACT_EXP_PRED_brms_2022a, pars = c("hu_EXPDistPixels_Intercept",
                                                             "hu_PREDDistPixels_Intercept",
                                                             "hu_EXPDistPixels_Treatment1",
                                                             "hu_EXPDistPixels_Treatment2",
                                                             "hu_PREDDistPixels_Treatment1",
                                                             "hu_PREDDistPixels_Treatment2"
                                                            ), summary = F)

#select individual level correlations among individuals
cor_id_list <- VarCorr(ACT_EXP_PRED_brms_2022a, summary = F)$TadpoleID$cor


column_names <- c("Mean", 
                  "2.5", 
                  "97.5"
)


row_names_ACTdummy <- c("Hurdle_Intercept",
                        "Hurdle_Conspecific",
                        "Hurdle_Heterospecific"
)

row_names_rep <- c("No Competition", "Conspecific", "Heterospecific")


#get summary of (non treatment) fixed effects
fixef_mod_ACT <- fixef_mod[,c(1:2)]
fixef_mod_ACT <- round(data.frame(colMeans(fixef_mod_ACT), HPDinterval(as.mcmc(fixef_mod_ACT))), digits = 3) 
fixef_mod_EXP <- fixef_mod[,c(3:4)]
fixef_mod_EXP <- round(data.frame(colMeans(fixef_mod_EXP), HPDinterval(as.mcmc(fixef_mod_EXP))), digits = 3) 
fixef_mod_PRED <- fixef_mod[,c(5:6)]
fixef_mod_PRED <- round(data.frame(colMeans(fixef_mod_PRED), HPDinterval(as.mcmc(fixef_mod_PRED))), digits = 3) 
colnames(fixef_mod_ACT) <- column_names
colnames(fixef_mod_EXP) <- column_names
colnames(fixef_mod_PRED) <- column_names

#get summary of treatment fixed effects
#activity assay
fixef_mod_treat_ACT_all <- fixef_mod_treat[,c(1, 4:5)]
fixef_mod_treat_ACT_nocomp <- fixef_mod_treat_ACT_all[,1]
fixef_mod_treat_ACT_consp <- fixef_mod_treat_ACT_all[,1] + fixef_mod_treat_ACT_all[,2]
fixef_mod_treat_ACT_hetero <- fixef_mod_treat_ACT_all[,1] + fixef_mod_treat_ACT_all[,3]
#combine treatment fixed effects into one df
fixef_mod_treat_ACT <- cbind(fixef_mod_treat_ACT_nocomp, fixef_mod_treat_ACT_consp, fixef_mod_treat_ACT_hetero)
fixef_mod_treat_ACT <- round(data.frame(colMeans(fixef_mod_treat_ACT), HPDinterval(as.mcmc(fixef_mod_treat_ACT))), digits = 3) 
colnames(fixef_mod_treat_ACT) <- column_names

#exploration assay
fixef_mod_treat_EXP_all <- fixef_mod_treat[,c(2, 6:7)]
fixef_mod_treat_EXP_nocomp <- fixef_mod_treat_EXP_all[,1]
fixef_mod_treat_EXP_consp <- fixef_mod_treat_EXP_all[,1] + fixef_mod_treat_EXP_all[,2]
fixef_mod_treat_EXP_hetero <- fixef_mod_treat_EXP_all[,1] + fixef_mod_treat_EXP_all[,3]
#combine treatment fixed effects into one df
fixef_mod_treat_EXP <- cbind(fixef_mod_treat_EXP_nocomp, fixef_mod_treat_EXP_consp, fixef_mod_treat_EXP_hetero)
fixef_mod_treat_EXP <- round(data.frame(colMeans(fixef_mod_treat_EXP), HPDinterval(as.mcmc(fixef_mod_treat_EXP))), digits = 3) 
colnames(fixef_mod_treat_EXP) <- column_names

#predation assay
fixef_mod_treat_PRED_all <- fixef_mod_treat[,c(3, 8:9)]
fixef_mod_treat_PRED_nocomp <- fixef_mod_treat_PRED_all[,1]
fixef_mod_treat_PRED_consp <- fixef_mod_treat_PRED_all[,1] + fixef_mod_treat_PRED_all[,2]
fixef_mod_treat_PRED_hetero <- fixef_mod_treat_PRED_all[,1] + fixef_mod_treat_PRED_all[,3]
#combine treatment fixed effects into one df
fixef_mod_treat_PRED <- cbind(fixef_mod_treat_PRED_nocomp, fixef_mod_treat_PRED_consp, fixef_mod_treat_PRED_hetero)
fixef_mod_treat_PRED <- round(data.frame(colMeans(fixef_mod_treat_PRED), HPDinterval(as.mcmc(fixef_mod_treat_PRED))), digits = 3) 
colnames(fixef_mod_treat_PRED) <- column_names

#get summary of eggmass
#Activity assay
v_eggmass_ACT <- v_eggmass[,1]
v_eggmass_ACT <- round(data.frame(mean(v_eggmass_ACT), HPDinterval(as.mcmc(v_eggmass_ACT))), digits = 3)
colnames(v_eggmass_ACT) <- column_names
#Exploration assay
v_eggmass_EXP <- v_eggmass[,2]
v_eggmass_EXP <- round(data.frame(mean(v_eggmass_EXP), HPDinterval(as.mcmc(v_eggmass_EXP))), digits = 3)
colnames(v_eggmass_EXP) <- column_names
#Predation assay
v_eggmass_PRED <- v_eggmass[,3]
v_eggmass_PRED <- round(data.frame(mean(v_eggmass_PRED), HPDinterval(as.mcmc(v_eggmass_PRED))), digits = 3)
colnames(v_eggmass_PRED) <- column_names

#get summary of SD among individuals and convert to variance
#Activity assay
v_id_ACT_all <- v_id[,c(1, 4, 7)]
v_id_ACT_nocomp <- v_id_ACT_all[,1]
v_id_ACT_nocomp <- (v_id_ACT_nocomp)^2
v_id_ACT_consp <- v_id_ACT_all[,2]^2
v_id_ACT_hetero <- v_id_ACT_all[,3]^2
v_id_ACT <- cbind(v_id_ACT_nocomp, v_id_ACT_consp, v_id_ACT_hetero)
v_id_ACT <- round(data.frame(colMeans(v_id_ACT), HPDinterval(as.mcmc(v_id_ACT, prob=0.95))), digits = 3) 
colnames(v_id_ACT) <- column_names

#Exploration assay
v_id_EXP_all <- v_id[,c(2, 5, 8)]
v_id_EXP_nocomp <- v_id_EXP_all[,1]
v_id_EXP_nocomp <- (v_id_EXP_nocomp)^2
v_id_EXP_consp <- v_id_EXP_all[,2]^2
v_id_EXP_hetero <- v_id_EXP_all[,3]^2
v_id_EXP <- cbind(v_id_EXP_nocomp, v_id_EXP_consp, v_id_EXP_hetero)
v_id_EXP <- round(data.frame(colMeans(v_id_EXP), HPDinterval(as.mcmc(v_id_EXP))), digits = 3) 
colnames(v_id_EXP) <- column_names

#Predation assay
v_id_PRED_all <- v_id[,c(3, 6, 9)]
v_id_PRED_nocomp <- v_id_PRED_all[,1]
v_id_PRED_nocomp <- (v_id_PRED_nocomp)^2
v_id_PRED_consp <- v_id_PRED_all[,2]^2
v_id_PRED_hetero <- v_id_PRED_all[,3]^2
v_id_PRED <- cbind(v_id_PRED_nocomp, v_id_PRED_consp, v_id_PRED_hetero)
v_id_PRED <- round(data.frame(colMeans(v_id_PRED), HPDinterval(as.mcmc(v_id_PRED))), digits = 3) 
colnames(v_id_PRED) <- column_names

#get summary of SD within individuals and convert to variance
v_sigma_ACT_all <- v_sigma[,c(1, 4, 5)]
#Activity assay
v_sigma_ACT_nocomp <- v_sigma_ACT_all[,1]
v_sigma_ACT_nocomp <- (exp(v_sigma_ACT_nocomp))^2
v_sigma_ACT_consp <- v_sigma_ACT_all[,1] + v_sigma_ACT_all[,2]
v_sigma_ACT_consp <- (exp(v_sigma_ACT_consp))^2
v_sigma_ACT_hetero <- v_sigma_ACT_all[,1] + v_sigma_ACT_all[,3]
v_sigma_ACT_hetero <- (exp(v_sigma_ACT_hetero))^2
v_sigma_ACT <- cbind(v_sigma_ACT_nocomp, v_sigma_ACT_consp, v_sigma_ACT_hetero)
v_sigma_ACT <- round(data.frame(colMeans(v_sigma_ACT), HPDinterval(as.mcmc(v_sigma_ACT, prob=0.95))), digits = 3) 
colnames(v_sigma_ACT) <- column_names

#Exploration assay
v_sigma_EXP_all <- v_sigma[,c(2, 6, 7)]
v_sigma_EXP_nocomp <- v_sigma_EXP_all[,1]
v_sigma_EXP_nocomp <- (exp(v_sigma_EXP_nocomp))^2
v_sigma_EXP_consp <- v_sigma_EXP_all[,1] + v_sigma_EXP_all[,2]
v_sigma_EXP_consp <- (exp(v_sigma_EXP_consp))^2
v_sigma_EXP_hetero <- v_sigma_EXP_all[,1] + v_sigma_EXP_all[,3]
v_sigma_EXP_hetero <- (exp(v_sigma_EXP_hetero))^2
v_sigma_EXP <- cbind(v_sigma_EXP_nocomp, v_sigma_EXP_consp, v_sigma_EXP_hetero)
v_sigma_EXP <- round(data.frame(colMeans(v_sigma_EXP), HPDinterval(as.mcmc(v_sigma_EXP, prob=0.95))), digits = 3) 
colnames(v_sigma_EXP) <- column_names

#Predation assay
v_sigma_PRED_all <- v_sigma[,c(3, 8, 9)]
v_sigma_PRED_nocomp <- v_sigma_PRED_all[,1]
v_sigma_PRED_nocomp <- (exp(v_sigma_PRED_nocomp))^2
v_sigma_PRED_consp <- v_sigma_PRED_all[,1] + v_sigma_PRED_all[,2]
v_sigma_PRED_consp <- (exp(v_sigma_PRED_consp))^2
v_sigma_PRED_hetero <- v_sigma_PRED_all[,1] + v_sigma_PRED_all[,3]
v_sigma_PRED_hetero <- (exp(v_sigma_PRED_hetero))^2
v_sigma_PRED <- cbind(v_sigma_PRED_nocomp, v_sigma_PRED_consp, v_sigma_PRED_hetero)
v_sigma_PRED <- round(data.frame(colMeans(v_sigma_PRED), HPDinterval(as.mcmc(v_sigma_PRED, prob=0.95))), digits = 3) 
colnames(v_sigma_PRED) <- column_names

#get summary of hurdle model
#Exploration assay
fixed_mod_hurdle_EXP_all <- fixed_mod_hurdle[,c(1, 3, 4)]
fixed_mod_hurdle_EXP_nocomp <- fixed_mod_hurdle_EXP_all[,1]
fixed_mod_hurdle_EXP_nocomp <- exp(fixed_mod_hurdle_EXP_nocomp)
#convert estimates from logit scale to probability
fixed_mod_hurdle_EXP_nocomp<- fixed_mod_hurdle_EXP_nocomp/(1+fixed_mod_hurdle_EXP_nocomp)
fixed_mod_hurdle_EXP_consp <- fixed_mod_hurdle_EXP_all[,1] +  fixed_mod_hurdle_EXP_all[,2]
fixed_mod_hurdle_EXP_consp <- exp(fixed_mod_hurdle_EXP_consp)
#convert estimates from logit scale to probability
fixed_mod_hurdle_EXP_consp <- fixed_mod_hurdle_EXP_consp/(1+fixed_mod_hurdle_EXP_consp)
fixed_mod_hurdle_EXP_hetero <- fixed_mod_hurdle_EXP_all[,1] +  fixed_mod_hurdle_EXP_all[,3]
fixed_mod_hurdle_EXP_hetero <- exp(fixed_mod_hurdle_EXP_hetero)
#convert estimates from logit scale to probability
fixed_mod_hurdle_EXP_hetero <- fixed_mod_hurdle_EXP_hetero/(1+fixed_mod_hurdle_EXP_hetero)
fixed_mod_hurdle_EXP <- cbind(fixed_mod_hurdle_EXP_nocomp, fixed_mod_hurdle_EXP_consp, fixed_mod_hurdle_EXP_hetero)
fixed_mod_hurdle_EXP <- round(data.frame(colMeans(fixed_mod_hurdle_EXP), HPDinterval(as.mcmc(fixed_mod_hurdle_EXP, prob=0.95))), digits = 3) 
colnames(fixed_mod_hurdle_EXP) <- column_names

#Predation assay
fixed_mod_hurdle_PRED_all <- fixed_mod_hurdle[,c(2, 5, 6)]
fixed_mod_hurdle_PRED_nocomp <- fixed_mod_hurdle_PRED_all[,1]
fixed_mod_hurdle_PRED_nocomp <- exp(fixed_mod_hurdle_PRED_nocomp)
#convert estimates from logit scale to probability
fixed_mod_hurdle_PRED_nocomp <- fixed_mod_hurdle_PRED_nocomp/(1+fixed_mod_hurdle_PRED_nocomp)
fixed_mod_hurdle_PRED_consp <- fixed_mod_hurdle_PRED_all[,1] +  fixed_mod_hurdle_PRED_all[,2]
fixed_mod_hurdle_PRED_consp <- exp(fixed_mod_hurdle_PRED_consp)
#convert estimates from logit scale to probability
fixed_mod_hurdle_PRED_consp <- fixed_mod_hurdle_PRED_consp/(1+fixed_mod_hurdle_PRED_consp)
fixed_mod_hurdle_PRED_hetero <- fixed_mod_hurdle_PRED_all[,1] +  fixed_mod_hurdle_PRED_all[,3]
fixed_mod_hurdle_PRED_hetero  <- exp(fixed_mod_hurdle_PRED_hetero)
#convert estimates from logit scale to probability
fixed_mod_hurdle_PRED_hetero <- fixed_mod_hurdle_PRED_hetero/(1+fixed_mod_hurdle_PRED_hetero)
fixed_mod_hurdle_PRED <- cbind(fixed_mod_hurdle_PRED_nocomp, fixed_mod_hurdle_PRED_consp, fixed_mod_hurdle_PRED_hetero)
fixed_mod_hurdle_PRED <- round(data.frame(colMeans(fixed_mod_hurdle_PRED), HPDinterval(as.mcmc(fixed_mod_hurdle_PRED, prob=0.95))), digits = 3) 
colnames(fixed_mod_hurdle_PRED) <- column_names


#combine results of fixed effects and variances into one dataframe for each assay
EXP_df <- rbind(fixef_mod_treat_EXP, fixef_mod_EXP, v_eggmass_EXP, v_id_EXP, v_sigma_EXP, fixed_mod_hurdle_EXP)
PRED_df <- rbind(fixef_mod_treat_PRED, fixef_mod_PRED, v_eggmass_PRED, v_id_PRED, v_sigma_PRED, fixed_mod_hurdle_PRED)
ACT_df <- rbind(fixef_mod_treat_ACT, fixef_mod_ACT, v_eggmass_ACT, v_id_ACT, v_sigma_ACT)


#activity assay does not have hurdle part of model, NA need to be added into the dataframe for these
ACT_dum <- data.frame(matrix(NA, nrow = 3, ncol = 3,
                             dimnames = list(NULL, paste0("ColumnName_", 1:3))) )
colnames(ACT_dum) <- column_names
rownames(ACT_dum) <- row_names_ACTdummy
ACT_df <- rbind(ACT_df,ACT_dum)

#combine activity, exploration and predation dataframes together and label the row names. 
ACT_EXP_PRED_df <- cbind(ACT_df, EXP_df, PRED_df)
rownames(ACT_EXP_PRED_df) <- c("No Competition",
                               "Conspecific",
                               "Heterospecific",
                               "SVL",
                               "Trial",
                               "Egg Mass ID",
                               "No Competition ",
                               "Conspecific ",
                               "Heterospecific ",
                               "No Competition  ",
                               "Conspecific  ",
                               "Heterospecific  ",
                               "No Competition   ",
                               "Conspecific   ",
                               "Heterospecific   "
                               )
  
  
#create kable table of the final dataframe for the multivariate model
ACT_EXP_PRED_kable <- ACT_EXP_PRED_df  %>%
  knitr::kable(
    format = "html",
    escape = TRUE
  ) %>%
  kable_styling(bootstrap_options = c("hover", "condensed"), full_width = F) %>% 
  add_header_above(c(" " = 2, "95% CI" = 2, " " = 1, "95% CI" = 2, " " = 1, "95% CI" = 2)) %>%
  add_header_above(c(" " = 1, "Activity" = 3, "Exploration" = 3, "Predation" = 3)) %>%
  pack_rows("Population means", 1, 5) %>%
  pack_rows("Variance among egg masses", 6, 6) %>%
  pack_rows("Variance among individuals", 7, 9) %>%
  pack_rows("Variance within individuals", 10, 12) %>%
  pack_rows("Probability remained in acclimation zone", 13, 15)
ACT_EXP_PRED_kable


#####get treatment difference in multivariate model results####

#Population mean treatment differences
#Activity assay
fixef_mod_treat_ACT_nocomp_vs_consp <- treat_compare(fixef_mod_treat_ACT_nocomp, fixef_mod_treat_ACT_consp)
fixef_mod_treat_ACT_nocomp_vs_hetero <- treat_compare(fixef_mod_treat_ACT_nocomp, fixef_mod_treat_ACT_hetero)
fixef_mod_treat_ACT_consp_vs_hetero <- treat_compare(fixef_mod_treat_ACT_consp, fixef_mod_treat_ACT_hetero)
#combine treatment effects into one df
fixef_mod_treat_ACT_diff <- rbind(fixef_mod_treat_ACT_nocomp_vs_consp, 
                                   fixef_mod_treat_ACT_nocomp_vs_hetero, 
                                   fixef_mod_treat_ACT_consp_vs_hetero)

#Exploration assay
fixef_mod_treat_EXP_nocomp_vs_consp <- treat_compare(fixef_mod_treat_EXP_nocomp, fixef_mod_treat_EXP_consp)
fixef_mod_treat_EXP_nocomp_vs_hetero <- treat_compare(fixef_mod_treat_EXP_nocomp, fixef_mod_treat_EXP_hetero)
fixef_mod_treat_EXP_consp_vs_hetero <- treat_compare(fixef_mod_treat_EXP_consp, fixef_mod_treat_EXP_hetero)
#combine treatment effects into one df
fixef_mod_treat_EXP_diff <- rbind(fixef_mod_treat_EXP_nocomp_vs_consp, 
                                   fixef_mod_treat_EXP_nocomp_vs_hetero, 
                                   fixef_mod_treat_EXP_consp_vs_hetero)

#Predation assay
fixef_mod_treat_PRED_nocomp_vs_consp <- treat_compare(fixef_mod_treat_PRED_nocomp, fixef_mod_treat_PRED_consp)
fixef_mod_treat_PRED_nocomp_vs_hetero <- treat_compare(fixef_mod_treat_PRED_nocomp, fixef_mod_treat_PRED_hetero)
fixef_mod_treat_PRED_consp_vs_hetero <- treat_compare(fixef_mod_treat_PRED_consp, fixef_mod_treat_PRED_hetero)
#combine treatment effects into one df
fixef_mod_treat_PRED_diff <- rbind(fixef_mod_treat_PRED_nocomp_vs_consp, 
                              fixef_mod_treat_PRED_nocomp_vs_hetero, 
                              fixef_mod_treat_PRED_consp_vs_hetero)


#Treatment comparisons for variance among individuals
#Activity assay
v_id_ACT_nocomp_vs_consp <- treat_compare(v_id_ACT_nocomp, v_id_ACT_consp)
v_id_ACT_nocomp_vs_hetero <- treat_compare(v_id_ACT_nocomp, v_id_ACT_hetero)
v_id_ACT_consp_vs_hetero <- treat_compare(v_id_ACT_consp, v_id_ACT_hetero)
#combine treatment effects into one df
v_id_ACT_diff <- rbind(v_id_ACT_nocomp_vs_consp, 
                        v_id_ACT_nocomp_vs_hetero, 
                        v_id_ACT_consp_vs_hetero)

#Exploration assay
v_id_EXP_nocomp_vs_consp <- treat_compare(v_id_EXP_nocomp, v_id_EXP_consp)
v_id_EXP_nocomp_vs_hetero <- treat_compare(v_id_EXP_nocomp, v_id_EXP_hetero)
v_id_EXP_consp_vs_hetero <- treat_compare(v_id_EXP_consp, v_id_EXP_hetero)
#combine treatment effects into one df
v_id_EXP_diff <- rbind(v_id_EXP_nocomp_vs_consp, 
                        v_id_EXP_nocomp_vs_hetero, 
                        v_id_EXP_consp_vs_hetero)

#Predation assay
v_id_PRED_nocomp_vs_consp <- treat_compare(v_id_PRED_nocomp, v_id_PRED_consp)
v_id_PRED_nocomp_vs_hetero <- treat_compare(v_id_PRED_nocomp, v_id_PRED_hetero)
v_id_PRED_consp_vs_hetero <- treat_compare(v_id_PRED_consp, v_id_PRED_hetero)
#combine treatment effects into one df
v_id_PRED_diff <- rbind(v_id_PRED_nocomp_vs_consp, 
                                   v_id_PRED_nocomp_vs_hetero, 
                                   v_id_PRED_consp_vs_hetero)

#variance within individuals treatment differences
#Activity assay
v_sigma_ACT_nocomp_vs_consp <- treat_compare(v_sigma_ACT_nocomp, v_sigma_ACT_consp)
v_sigma_ACT_nocomp_vs_hetero <- treat_compare(v_sigma_ACT_nocomp, v_sigma_ACT_hetero)
v_sigma_ACT_consp_vs_hetero <- treat_compare(v_sigma_ACT_consp, v_sigma_ACT_hetero)
#combine treatment effects into one df
v_sigma_ACT_diff <- rbind(v_sigma_ACT_nocomp_vs_consp, 
                           v_sigma_ACT_nocomp_vs_hetero, 
                           v_sigma_ACT_consp_vs_hetero)

#Exploration assay
v_sigma_EXP_nocomp_vs_consp <- treat_compare(v_sigma_EXP_nocomp, v_sigma_EXP_consp)
v_sigma_EXP_nocomp_vs_hetero <- treat_compare(v_sigma_EXP_nocomp, v_sigma_EXP_hetero)
v_sigma_EXP_consp_vs_hetero <- treat_compare(v_sigma_EXP_consp, v_sigma_EXP_hetero)
#combine treatment effects into one df
v_sigma_EXP_diff <- rbind(v_sigma_EXP_nocomp_vs_consp, 
                           v_sigma_EXP_nocomp_vs_hetero, 
                           v_sigma_EXP_consp_vs_hetero)

#Predation assay
v_sigma_PRED_nocomp_vs_consp <- treat_compare(v_sigma_PRED_nocomp, v_sigma_PRED_consp)
v_sigma_PRED_nocomp_vs_hetero <- treat_compare(v_sigma_PRED_nocomp, v_sigma_PRED_hetero)
v_sigma_PRED_consp_vs_hetero <- treat_compare(v_sigma_PRED_consp, v_sigma_PRED_hetero)
#combine treatment effects into one df
v_sigma_PRED_diff <- rbind(v_sigma_PRED_nocomp_vs_consp, 
                        v_sigma_PRED_nocomp_vs_hetero, 
                        v_sigma_PRED_consp_vs_hetero)


#hurdle model treatment differences
#Exploration assay
fixed_mod_hurdle_EXP_nocomp_vs_consp <- treat_compare(fixed_mod_hurdle_EXP_nocomp, fixed_mod_hurdle_EXP_consp)
fixed_mod_hurdle_EXP_nocomp_vs_hetero <- treat_compare(fixed_mod_hurdle_EXP_nocomp, fixed_mod_hurdle_EXP_hetero)
fixed_mod_hurdle_EXP_consp_vs_hetero <- treat_compare(fixed_mod_hurdle_EXP_consp, fixed_mod_hurdle_EXP_hetero)
#combine treatment effects into one df
fixed_mod_hurdle_EXP_diff <- rbind(fixed_mod_hurdle_EXP_nocomp_vs_consp, 
                                    fixed_mod_hurdle_EXP_nocomp_vs_hetero, 
                                    fixed_mod_hurdle_EXP_consp_vs_hetero)

#predation assay
fixed_mod_hurdle_PRED_nocomp_vs_consp <- treat_compare(fixed_mod_hurdle_PRED_nocomp, fixed_mod_hurdle_PRED_consp)
fixed_mod_hurdle_PRED_nocomp_vs_hetero <- treat_compare(fixed_mod_hurdle_PRED_nocomp, fixed_mod_hurdle_PRED_hetero)
fixed_mod_hurdle_PRED_consp_vs_hetero <- treat_compare(fixed_mod_hurdle_PRED_consp, fixed_mod_hurdle_PRED_hetero)
#combine treatment effects into one df
fixed_mod_hurdle_PRED_diff <- rbind(fixed_mod_hurdle_PRED_nocomp_vs_consp, 
                           fixed_mod_hurdle_PRED_nocomp_vs_hetero, 
                           fixed_mod_hurdle_PRED_consp_vs_hetero)

#activity assay does not have hurdle part of model, NA need to be added into the dataframe for these
ACT_diff_dum <- data.frame(matrix(NA, nrow = 3, ncol = 3,
                             dimnames = list(NULL, paste0("ColumnName_", 1:3))) )
colnames(ACT_diff_dum) <- column_names

#combine results of treatment differences in the mean population, variance among individuals, variance within individuals and hurdle model into one dataframe for each assay
EXP_diff_df <- rbind(fixef_mod_treat_EXP_diff, v_id_EXP_diff, v_sigma_EXP_diff, fixed_mod_hurdle_EXP_diff)
PRED_diff_df <- rbind(fixef_mod_treat_PRED_diff, v_id_PRED_diff, v_sigma_PRED_diff, fixed_mod_hurdle_PRED_diff)
ACT_diff_df <- rbind(fixef_mod_treat_ACT_diff, v_id_ACT_diff, v_sigma_ACT_diff, ACT_diff_dum)

#combine activity, exploration and predation dataframes together and label the row names. 
ACT_EXP_PRED_diff_df <- cbind(ACT_diff_df, EXP_diff_df, PRED_diff_df)
rownames(ACT_EXP_PRED_diff_df) <- c("No Competition - Conspecific",
                               "No Competition - Heterospecific",
                               "Conspecific - Heterospecific",
                               "No Competition - Conspecific ",
                               "No Competition - Heterospecific ",
                               "Conspecific - Heterospecific ",
                               "No Competition - Conspecific  ",
                               "No Competition - Heterospecific  ",
                               "Conspecific - Heterospecific  ",
                               "No Competition - Conspecific    ",
                               "No Competition - Heterospecific    ",
                               "Conspecific - Heterospecific    "
)

#create kable table of the final dataframe
ACT_EXP_PRED_diff_kable <- ACT_EXP_PRED_diff_df  %>%
  knitr::kable(
    format = "html",
    escape = TRUE
  ) %>%
  kable_styling(bootstrap_options = c("hover", "condensed"), full_width = F) %>% 
  add_header_above(c(" " = 2, "95% CI" = 2, " " = 1, "95% CI" = 2, " " = 1, "95% CI" = 2)) %>%
  add_header_above(c(" " = 1, "Activity" = 3, "Exploration" = 3, "Predation" = 3)) %>%
  pack_rows("Population means", 1, 3) %>%
  pack_rows("Variance among individuals", 4, 6) %>%
  pack_rows("Variance within individuals", 7, 9) %>%
  pack_rows("Probability remained in acclimation zone", 10, 12)
ACT_EXP_PRED_diff_kable


#####get multivariate model results - posterior repeatability estimates ####

#Calculate adjusted estimates of repeatability from multivariate model posterior estimates
#Activity assay
rep_ACT_nocomp <- v_id_ACT_nocomp/(v_id_ACT_nocomp+v_sigma_ACT_nocomp)
rep_ACT_consp <- v_id_ACT_consp/(v_id_ACT_consp+v_sigma_ACT_consp)
rep_ACT_hetero <- v_id_ACT_hetero/(v_id_ACT_hetero+v_sigma_ACT_hetero)
#combine treatment effects into single df
rep_ACT <- cbind(rep_ACT_nocomp, rep_ACT_consp, rep_ACT_hetero)
rep_ACT <- round(data.frame(colMeans(rep_ACT), HPDinterval(as.mcmc(rep_ACT))), digits = 3) 
colnames(rep_ACT) <- column_names

#Exploration assay
rep_EXP_nocomp <- v_id_EXP_nocomp/(v_id_EXP_nocomp+v_sigma_EXP_nocomp)
rep_EXP_consp <- v_id_EXP_consp/(v_id_EXP_consp+v_sigma_EXP_consp)
rep_EXP_hetero <- v_id_EXP_hetero/(v_id_EXP_hetero+v_sigma_EXP_hetero)
#combine treatment effects into single df
rep_EXP <- cbind(rep_EXP_nocomp, rep_EXP_consp, rep_EXP_hetero)
rep_EXP <- round(data.frame(colMeans(rep_EXP), HPDinterval(as.mcmc(rep_EXP))), digits = 3) 
colnames(rep_EXP) <- column_names

#Predation assay
rep_PRED_nocomp <- v_id_PRED_nocomp/(v_id_PRED_nocomp+v_sigma_PRED_nocomp)
rep_PRED_consp <- v_id_PRED_consp/(v_id_PRED_consp+v_sigma_PRED_consp)
rep_PRED_hetero <- v_id_PRED_hetero/(v_id_PRED_hetero+v_sigma_PRED_hetero)
#combine treatment effects into single df
rep_PRED <- cbind(rep_PRED_nocomp, rep_PRED_consp, rep_PRED_hetero)
rep_PRED <- round(data.frame(colMeans(rep_PRED), HPDinterval(as.mcmc(rep_PRED))), digits = 3) 
colnames(rep_PRED) <- column_names

#combine repeatability estimates from each assay into a single data frame
ACT_EXP_PRED_rep_df <- cbind(rep_ACT, rep_EXP, rep_PRED)
rownames(ACT_EXP_PRED_rep_df) <- c("No Competition",
                                   "Conspecific",
                                   "Heterospecific"
)

#create kable table of the repeatability results
ACT_EXP_PRED_rep_df_kable <- ACT_EXP_PRED_rep_df  %>%
  knitr::kable(
    format = "html",
    escape = TRUE
  ) %>%
  kable_styling(bootstrap_options = c("hover", "condensed"), full_width = F) %>% 
  add_header_above(c(" " = 2, "95% CI" = 2, " " = 1, "95% CI" = 2, " " = 1, "95% CI" = 2)) %>%
  add_header_above(c(" " = 1, "Activity" = 3, "Exploration" = 3, "Predation" = 3))
ACT_EXP_PRED_rep_df_kable

#####get treatment differences between repeatability estimates####
#Activity assay
rep_ACT_nocomp_vs_consp <- treat_compare(rep_ACT_nocomp, rep_ACT_consp)
rep_ACT_nocomp_vs_hetero <- treat_compare(rep_ACT_nocomp, rep_ACT_hetero)
rep_ACT_consp_vs_hetero <- treat_compare(rep_ACT_consp, rep_ACT_hetero)
#combine treatment differences into single df
rep_ACT_diff <- rbind(rep_ACT_nocomp_vs_consp, rep_ACT_nocomp_vs_hetero, rep_ACT_consp_vs_hetero)

#Exploration assay
rep_EXP_nocomp_vs_consp <- treat_compare(rep_EXP_nocomp, rep_EXP_consp)
rep_EXP_nocomp_vs_hetero <- treat_compare(rep_EXP_nocomp, rep_EXP_hetero)
rep_EXP_consp_vs_hetero <- treat_compare(rep_EXP_consp, rep_EXP_hetero)
#combine treatment differences into single df
rep_EXP_diff <- rbind(rep_EXP_nocomp_vs_consp, rep_EXP_nocomp_vs_hetero, rep_EXP_consp_vs_hetero)

#Predation assay
rep_PRED_nocomp_vs_consp <- treat_compare(rep_PRED_nocomp, rep_PRED_consp)
rep_PRED_nocomp_vs_hetero <- treat_compare(rep_PRED_nocomp, rep_PRED_hetero)
rep_PRED_consp_vs_hetero <- treat_compare(rep_PRED_consp, rep_PRED_hetero)
#combine treatment differences into single df
rep_PRED_diff <- rbind(rep_PRED_nocomp_vs_consp, rep_PRED_nocomp_vs_hetero, rep_PRED_consp_vs_hetero)

#combine treatment differences of repeatability estimates from each assay into a single data frame
ACT_EXP_PRED_rep_diff_df <- cbind(rep_ACT_diff, rep_EXP_diff, rep_PRED_diff)
rownames(ACT_EXP_PRED_rep_diff_df) <- c("No Competition - Conspecific",
                                        "No Competition - Heterospecific",
                                        "Conspecific - Heterospecific"
)

#create kable table of the final dataframe
ACT_EXP_PRED_rep_diff_df_kable <- ACT_EXP_PRED_rep_diff_df  %>%
  knitr::kable(
    format = "html",
    escape = TRUE
  ) %>%
  kable_styling(bootstrap_options = c("hover", "condensed"), full_width = F) %>% 
  add_header_above(c(" " = 2, "95% CI" = 2, " " = 1, "95% CI" = 2, " " = 1, "95% CI" = 2)) %>%
  add_header_above(c(" " = 1, "Activity" = 3, "Exploration" = 3, "Predation" = 3))
ACT_EXP_PRED_rep_diff_df_kable


##### get multivariate model results - posterior correlation estimates #####
row_names_corr <-  c(
  paste(c("ActVI", "ExpVI", "PredVI"), rep("ActVI", 3), sep = "_"),
  paste(c("ActVI", "ExpVI", "PredVI"), rep("ExpVI", 3), sep = "_"),
  paste(c("ActVI", "ExpVI", "PredVI"), rep("PredVI", 3), sep = "_")
)


#obtain pairwise among individual correlations between Activity, Exploration and Predation assays
cor_id_list <- VarCorr(ACT_EXP_PRED_brms_2022a, summary = FALSE)$TadpoleID$cor
cor_id_list_nocomp <- cor_id_list[, 1:3, 1:3] %>%
  as.data.frame.array()
cor_id_list_consp <- cor_id_list[, 4:6, 4:6] %>%
  as.data.frame.array()
cor_id_list_hetero <- cor_id_list[, 7:9, 7:9] %>%
  as.data.frame.array()

#obtain mean estimate and 95% CI of pairwise correalations between assays for each treatment
#no competition treatment
cor_id_list_nocomp <-  round(data.frame(colMeans(cor_id_list_nocomp), HPDinterval(as.mcmc(cor_id_list_nocomp))), digits = 3) %>%
  `rownames<-`(row_names_corr) %>%
  `colnames<-`(column_names) %>%
  slice(-c(1,4,5,7,8,9))
#conspeific treatment
cor_id_list_consp  <-  round(data.frame(colMeans(cor_id_list_consp), HPDinterval(as.mcmc(cor_id_list_consp))), digits = 3)  %>%
  `rownames<-`(row_names_corr) %>%
  `colnames<-`(column_names) %>%
  slice(-c(1,4,5,7,8,9))
#heterospeific treatment
cor_id_list_hetero  <-  round(data.frame(colMeans(cor_id_list_hetero), HPDinterval(as.mcmc(cor_id_list_hetero))), digits = 3)  %>%
  `rownames<-`(row_names_corr) %>%
  `colnames<-`(column_names) %>%
  slice(-c(1,4,5,7,8,9))


#create single dataframe of pairwise correalations for all treatments
cor_nocomp_consp_hetero <- cbind(cor_id_list_nocomp, cor_id_list_consp, cor_id_list_hetero) 

#create table of pairwise correalations for all treatments
cor_nocomp_consp_hetero_kable <- cor_nocomp_consp_hetero  %>%
  knitr::kable(
    format = "html",
    escape = TRUE
  ) %>%
  kable_styling(bootstrap_options = c("hover", "condensed"), full_width = F) %>% 
  add_header_above(c(" " = 2, "95% CI" = 2, " " = 1, "95% CI" = 2, " " = 1, "95% CI" = 2)) %>%
  add_header_above(c(" " = 1, "No Competition" = 3, "Conspecific" = 3, "Heterospecific" = 3))
cor_nocomp_consp_hetero_kable


#####get treatment differences between among individual correlation estimates####
#pairewise correlations between assays, treatment differences
cor_id_list_nocomp_consp <- cor_id_list_nocomp - cor_id_list_consp
cor_id_list_nocomp_hetero <- cor_id_list_nocomp - cor_id_list_hetero
cor_id_list_consp_hetero <- cor_id_list_consp - cor_id_list_hetero

#obtain mean estimate and 95% CI of pairwise correalations between assays - treatment comparisons
#no competition - conspecific treatment
cor_id_list_nocomp_consp <-  round(data.frame(colMeans(cor_id_list_nocomp_consp), HPDinterval(as.mcmc(cor_id_list_nocomp_consp))), digits = 3) %>%
  `rownames<-`(row_names_corr) %>%
  `colnames<-`(column_names) %>%
  slice(-c(1,4,5,7,8,9))
#no competition - heterospeific treatment
cor_id_list_nocomp_hetero  <-  round(data.frame(colMeans(cor_id_list_nocomp_hetero), HPDinterval(as.mcmc(cor_id_list_nocomp_hetero))), digits = 3)  %>%
  `rownames<-`(row_names_corr) %>%
  `colnames<-`(column_names) %>%
  slice(-c(1,4,5,7,8,9))
#conspecific - heterospeific treatment
cor_id_list_consp_hetero  <-  round(data.frame(colMeans(cor_id_list_consp_hetero), HPDinterval(as.mcmc(cor_id_list_consp_hetero))), digits = 3)  %>%
  `rownames<-`(row_names_corr) %>%
  `colnames<-`(column_names) %>%
  slice(-c(1,4,5,7,8,9))


#create single dataframe of pairwise correalations for all treatment comparisons
cor_nocomp_consp_nocomp_hetero_consp_hetero <- cbind(cor_id_list_nocomp_consp, cor_id_list_nocomp_hetero, cor_id_list_consp_hetero)


#create table of pairwise correlations for all treatment comparisons
cor_nocomp_consp_nocomp_hetero_consp_hetero_kable <- cor_nocomp_consp_nocomp_hetero_consp_hetero  %>%
  knitr::kable(
    format = "html",
    escape = TRUE
  ) %>%
  kable_styling(bootstrap_options = c("hover", "condensed"), full_width = T) %>% 
  add_header_above(c(" " = 2, "95% CI" = 2, " " = 1, "95% CI" = 2, " " = 1, "95% CI" = 2)) %>%
  add_header_above(c(" " = 1, "No Competition-Conspecific" = 3, "No Competition-Heterospecific" = 3, "Conspecific-Heterospecific" = 3))
cor_nocomp_consp_nocomp_hetero_consp_hetero_kable


#### population, variance among, variance within and repeatability plots ######

#Effect of treatment on population estimates of activity, exploration and predatory risk taking assays

ACT_EXP_PRED_mean_plot <- 
  ggplot(
    (rbind(fixef_mod_treat_ACT, fixef_mod_treat_ACT_diff, fixef_mod_treat_EXP, fixef_mod_treat_EXP_diff, fixef_mod_treat_PRED, fixef_mod_treat_PRED_diff) %>%
       mutate(
         trait = rep(c("Activity", "Exploration", "Predation"),each = 6),
         treat = rep(c("6.No Competition", "5.Conspecific", "4.Heterospecific", "3.NC-Cdiff", "2.NC-Hdiff", "1.C-Hdiff"),3)
       ) %>%
       rename(lower = "2.5",
              upper = "97.5")
    ),
    aes(y = treat, color = treat)
  )+
  geom_point(aes(x = Mean), size = 5, alpha = 0.7) +
  geom_linerange(aes(xmin = lower, xmax = upper)) +
  labs(
    x = "\\n Posterior Estimate with 95% Credible Intervals",
    y = "Treatment"
  ) +
  facet_grid(. ~ trait, scales = "free") +
  scale_y_discrete(labels = c(bquote(Delta ~ "Conspecific-\\n Heterospecific"), bquote(Delta ~ "No Competition-\\n Heterospecific"), bquote(Delta ~ "No Competition-\\n Conspecific"), "Heterospecific", "Conspecific", "No Competition")) +
  scale_colour_manual(values = c("#56B4E9", "#009E73", "#CC79A7", "#D55E00", "#E69F00", "#000000")) +
  geom_vline(xintercept = c(0, 0), linetype = "dotted") +
  theme_catplots() +
  theme(axis.text.x = element_text(face="bold", size=10))
ACT_EXP_PRED_mean_plot


#Effect of treatment on among individual variance in activity, exploration and predatory risk taking assays
ACT_EXP_PRED_var_id_plot <- 
  ggplot(
    (rbind(v_id_ACT, v_id_ACT_diff, v_id_EXP, v_id_EXP_diff, v_id_PRED, v_id_PRED_diff) %>%
  mutate(
    trait = rep(c("Activity", "Exploration", "Predation"),each = 6),
    treat = rep(c("6.No Competition", "5.Conspecific", "4.Heterospecific", "3.NC-Cdiff", "2.NC-Hdiff", "1.C-Hdiff"),3)
  ) %>%
    rename(lower = "2.5",
           upper = "97.5")
  ),
  aes(y = treat, color = treat)
  )+
  geom_point(aes(x = Mean), size = 5, alpha = 0.7) +
  geom_linerange(aes(xmin = lower, xmax = upper)) +
  labs(
    x = "\\n Posterior Estimate with 95% Credible Intervals",
    y = "Treatment"
  ) +
  facet_grid(. ~ trait, scales = "free") +
  scale_y_discrete(labels = c(bquote(Delta ~ "Conspecific-\\n Heterospecific"), bquote(Delta ~ "No Competition-\\n Heterospecific"), bquote(Delta ~ "No Competition-\\n Conspecific"), "Heterospecific", "Conspecific", "No Competition")) +
  scale_colour_manual(values = c("#56B4E9", "#009E73", "#CC79A7", "#D55E00", "#E69F00", "#000000")) +
  geom_vline(xintercept = c(0, 0), linetype = "dotted") +
  theme_catplots() +
  theme(axis.text.x = element_text(face="bold", size=10))
ACT_EXP_PRED_var_id_plot


#Effect of treatment on within individual variance in activity, exploration and predatory risk taking assays
ACT_EXP_PRED_var_sigma_plot <- 
  ggplot(
    (rbind(v_sigma_ACT, v_sigma_ACT_diff, v_sigma_EXP, v_sigma_EXP_diff, v_sigma_PRED, v_sigma_PRED_diff) %>%
       mutate(
         trait = rep(c("Activity", "Exploration", "Predation"),each = 6),
         treat = rep(c("6.No Competition", "5.Conspecific", "4.Heterospecific", "3.NC-Cdiff", "2.NC-Hdiff", "1.C-Hdiff"),3)
       ) %>%
       rename(lower = "2.5",
              upper = "97.5")
    ),
    aes(y = treat, color = treat)
  ) +
  geom_point(aes(x = Mean), size = 5, alpha = 0.7) +
  geom_linerange(aes(xmin = lower, xmax = upper)) +
  labs(
    x = " \\n Posterior Estimate with 95% Credible Intervals", 
    y = "Treatment" 
  ) +
  facet_grid(. ~ trait, scales = "free") +
  scale_y_discrete(labels = c(bquote(Delta ~ "Conspecific-\\n Heterospecific"), bquote(Delta ~ "No Competition-\\n Heterospecific"), bquote(Delta ~ "No Competition-\\n Conspecific"), "Heterospecific", "Conspecific", "No Competition")) +
  scale_colour_manual(values = c("#56B4E9", "#009E73", "#CC79A7", "#D55E00", "#E69F00", "#000000")) +
  geom_vline(xintercept = c(0, 0), linetype = "dotted") +
  theme_catplots() +
  theme(axis.text.x = element_text(face="bold", size=10))
ACT_EXP_PRED_var_sigma_plot

#Effect of treatment on behavioral repeatability in activity, exploration and predatory risk taking assays
ACT_EXP_PRED_repeat_plot <- 
  ggplot(
    (rbind(rep_ACT, rep_ACT_diff, rep_EXP, rep_EXP_diff, rep_PRED, rep_PRED_diff) %>%
       mutate(
         trait = rep(c("Activity", "Exploration", "Predation"),each = 6),
         treat = rep(c("6.No Competition", "5.Conspecific", "4.Heterospecific", "3.NC-Cdiff", "2.NC-Hdiff", "1.C-Hdiff"),3)
       ) %>%
       rename(lower = "2.5",
              upper = "97.5")
    ),
    aes(y = treat, color = treat)
  ) +
  geom_point(aes(x = Mean), size = 5, alpha = 0.7) +
  geom_linerange(aes(xmin = lower, xmax = upper)) +
  labs(
    x = " \\n Posterior Estimate with 95% Credible Intervals", 
    y = "Treatment" 
  ) +
  facet_grid(. ~ trait, scales = "free") +
  scale_y_discrete(labels = c(bquote(Delta ~ "Conspecific-\\n Heterospecific"), bquote(Delta ~ "No Competition-\\n Heterospecific"), bquote(Delta ~ "No Competition-\\n Conspecific"), "Heterospecific", "Conspecific", "No Competition")) +
  scale_colour_manual(values = c("#56B4E9", "#009E73", "#CC79A7", "#D55E00", "#E69F00", "#000000")) +
  geom_vline(xintercept = c(0, 0), linetype = "dotted") +
  theme_catplots() +
  theme(axis.text.x = element_text(face="bold", size=12))
ACT_EXP_PRED_repeat_plot

#Effect of treatment on the proportion of tadpoles remaining in the acclimation zone in activity, exploration and predatory risk taking assays
ACT_EXP_PRED_hurdle_plot <- 
  ggplot(
    (rbind(fixed_mod_hurdle_EXP, fixed_mod_hurdle_EXP_diff, fixed_mod_hurdle_PRED, fixed_mod_hurdle_PRED_diff) %>%
       mutate(
         trait = rep(c("Exploration", "Predation"),each = 6),
         treat = rep(c("6.No Competition", "5.Conspecific", "4.Heterospecific", "3.NC-Cdiff", "2.NC-Hdiff", "1.C-Hdiff"),2)
       ) %>%
       rename(lower = "2.5",
              upper = "97.5")
    ),
    aes(y = treat, color = treat)
  ) +
  geom_point(aes(x = Mean), size = 5, alpha = 0.7) +
  geom_linerange(aes(xmin = lower, xmax = upper)) +
  labs(
    x = " \\n Posterior Estimate with 95% Credible Intervals", 
    y = "Treatment" 
  ) +
  facet_grid(. ~ trait, scales = "free") +
  scale_y_discrete(labels = c(bquote(Delta ~ "Conspecific-\\n Heterospecific"), bquote(Delta ~ "No Competition-\\n Heterospecific"), bquote(Delta ~ "No Competition-\\n Conspecific"), "Heterospecific", "Conspecific", "No Competition")) +
  scale_colour_manual(values = c("#56B4E9", "#009E73", "#CC79A7", "#D55E00", "#E69F00", "#000000")) +
  geom_vline(xintercept = c(0, 0), linetype = "dotted") +
  theme_catplots() +
  theme(axis.text.x = element_text(face="bold", size=12))
ACT_EXP_PRED_hurdle_plot
