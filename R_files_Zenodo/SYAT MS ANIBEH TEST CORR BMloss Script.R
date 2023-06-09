#................................................................

## Analyses for manuscript: 
## Haemosporidian infections influence risk-taking behaviours in young male blackcaps Sylvia atricapilla

#................................................................

#Influence of behaviours associated with infection on host energy reserves ####

## SUMMARY OF THE SCRIPT:----

# **Cautious note**: YOU MUST FIRST RUN DATA SCRIPT
# TO GENERATE ADEQUATED TRANSFORMED DATA FOR ANALYSES

##................................................................

library(dplyr)
library(tidyverse)
library(purrr)
library(vegan)
library(ggplot2)
library(reshape2)
library(psych)
library(nFactors)
library(lme4)
library(lmerTest)
library(rstatix)
library(car)
library(PerformanceAnalytics )
library(ggpubr)

library(lsmeans)
library(afex)

library(Hmisc)

    set.seed (123)
#................................................................

    # (1) MERGE DATABASE FROM EACH BEHAVIOUR ####
    #.......First, include only variable of interest and rename ====
    
# We included only behavioural variables significantly related 
# with haemosporidian infection and percentage of weight loss 
# during the experiment 

    ##..............EXPLORATION ====

    EXPL_BMloss_red <- MS_SYAT_BH_EXPL[,c("ID",
                                      "T_INITIAL_log",
                                      "T_CROSS_log")] %>%
    dplyr::rename(T_INITIAL_log_EXPL = T_INITIAL_log,
                  T_CROSS_log_EXPL = T_CROSS_log)

    ##..............FORAGING ==== 
    # To control the effect on body mass loss of birds that ate or not during 
    # experiment
    FORAG_BIOL_DB_red <- FORAG_BIOL_DB[, c("ID","EAT_BIN_f")]

    ## ..............PREDATOR PRESENCE IMMEDIATE BEHAVIOUR (Pi)====

    names(MS_SYAT_BH_Pi)

    PRED_i_vars <- c("ID","STIMULUS",
                 "KEY_BEHAV", "RBT")

    PRED_i_BMloss_red <- MS_SYAT_BH_Pi[,PRED_i_vars]

    PRED_i_BMloss_red$KEY_BEHAV <- as.numeric(paste(PRED_i_BMloss_red$KEY_BEHAV))

    ### From long to wide format
    names(PRED_i_BMloss_red)

    PRED_i_BMloss_red_w  <- reshape(PRED_i_BMloss_red,
                                idvar = "ID",
                                timevar = "STIMULUS",
                                direction = "wide"
    )

    ### Differences between type of stimulus

    PRED_i_BMloss_red_w <- mutate(PRED_i_BMloss_red_w,
                              Crown_dif= (KEY_BEHAV.PRED-
                                              KEY_BEHAV.CTRL),
                              Crown_match = ifelse(Crown_dif==1,1,0),
                              RBT_PRED_dif=(RBT.PRED-
                                                RBT.CTRL))

    ## ..............ALARM CALL ====        

    Acall_BMloss_red <- MS_SYAT_BH_Acall[,c("ID","Acall_BIN")]

    ## ..............BODY MASS LOSS ====
    names(MS_SYAT_BH_BIOL)

    BMloss_red <- MS_SYAT_BH_BIOL[,c("ID","BM_Loss")]

#.......Second, MERGE DATABASES ====

PARAS_BMloss_DB_ALL <- Reduce(function(x,y) merge(x,y,
                                                  by="ID",
                                                  all=TRUE),
                              list(EXPL_BMloss_red,
                                   PRED_i_BMloss_red_w,
                                   Acall_BMloss_red,
                                   BMloss_red,
                                   FORAG_BIOL_DB_red))
    
PARAS_BMloss_DB_ALL_f <- merge (PARAS_BMloss_DB_ALL, 
                                    MS_SYAT_BH_BIOL,
                                    by.x ="ID", 
                                    by.y="ID",
                                    all.x = TRUE,
                                    no.dups=FALSE)
    
PARAS_BMloss_DB_ALL_f <- mutate(PARAS_BMloss_DB_ALL_f,
                                    INF_BIN_f = as.factor(INF_BIN),
                                    BM_Loss = BM_Loss.x)
#................................................................

# (2) CORRELATIONS BETWEEN VARIABLES OF INTEREST AND BODY MASS LOSS ####

# .......DIFFERENCES IN BODY MASS WITH FORAGING ####

BMloss_EAT.lm <- lm (BM_Loss ~EAT_BIN_f,
                   data = PARAS_BMloss_DB_ALL_f)
Anova(BMloss_EAT.lm
  ,type = 3)

BMloss_EAT <- residuals(BMloss_EAT.lm)

PARAS_BMloss_DB_ALL_f$BMloss_EAT <- BMloss_EAT

BMloss_EAT_meansd <- PARAS_BMloss_DB_ALL_f %>% 
  group_by(EAT_BIN_f) %>%
  summarise(BMloss_mean = mean(BM_Loss, na.rm = TRUE),
            BMloss_sd = sd (BM_Loss, na.rm = TRUE))

BMloss_EAT_meansd

#.....................................................................
    # .....Differences of continuous variables ====  
    #.......BM_Loss and T_INITIAL_log_EXPL ====

    #........ STATUS OF INFECTION ===

    Bmass_Texpl_glm.0 <- lm(BMloss_EAT ~ T_INITIAL_log_EXPL*INF_BIN_f ,
                      data=PARAS_BMloss_DB_ALL_f )

    par(mfrow = c(2, 3))
    plot( Bmass_Texpl_glm.0,
      ask = F, which = 1:6) 

    Anova(Bmass_Texpl_glm.0, type =3)

    # Not significant interaction

    Bmass_Texpl_glm.1 <- lm(BMloss_EAT ~ T_INITIAL_log_EXPL+INF_BIN_f,
                        data=PARAS_BMloss_DB_ALL_f )

    Anova(Bmass_Texpl_glm.1, type =3)
    summary(Bmass_Texpl_glm.1)
    
    # We removed status of infection
    
    Bmass_Texpl_glm <- lm(BMloss_EAT ~ T_INITIAL_log_EXPL,
                            data=PARAS_BMloss_DB_ALL_f )
    
    Anova(Bmass_Texpl_glm, type =3)
    summary(Bmass_Texpl_glm)
    Bmass_Texpl_glm.anova <- Anova(Bmass_Texpl_glm, type =3)

#.....................................................................

    #.......BM_Loss and RBT_PRED_dif====
    
    #........ STATUS OF INFECTION ===

    Bmass_RBT_glm.0 <- lm(BMloss_EAT ~ RBT_PRED_dif*INF_BIN_f,
                      data=PARAS_BMloss_DB_ALL_f )

    Anova(Bmass_RBT_glm.0, type =3)

    # Not significant interaction

    Bmass_RBT_glm.1 <- lm(BMloss_EAT ~ RBT_PRED_dif+INF_BIN_f,
                    data=PARAS_BMloss_DB_ALL_f )

    Anova(Bmass_RBT_glm.1, type =3)

    # We removed status of infection
    
    Bmass_RBT_glm <- lm(BMloss_EAT ~ RBT_PRED_dif,
                          data=PARAS_BMloss_DB_ALL_f )
    
    Anova(Bmass_RBT_glm, type =3)
    summary(Bmass_RBT_glm)
    
    Bmass_RBT_glm.anova <- Anova(Bmass_RBT_glm, type =3)
#.................................................................
    # .....Differences of dichotomous variables ====          
    # Convert to factor

    PARAS_BMloss_DB_ALL_f <- mutate(PARAS_BMloss_DB_ALL_f,
                              Crown_match_f = factor(Crown_match),
                              Acall_BIN_f = factor(Acall_BIN)
    )

#.................................................................
    #.......BM_Loss and Crown performance according to risk ====

    #........ STATUS OF INFECTION ===
    
    Bmass_Crown_glm.0 <- lm(BMloss_EAT ~ Crown_match_f*INF_BIN_f,
                      data=PARAS_BMloss_DB_ALL_f)

    Anova(Bmass_Crown_glm.0,
      type=3)
    
    # interaction not significant
    
    Bmass_Crown_glm.1 <- lm(BMloss_EAT ~ Crown_match_f+INF_BIN_f,
                            data=PARAS_BMloss_DB_ALL_f)
    
    Anova(Bmass_Crown_glm.1,
          type=3)
    
    # We removed status of infection
    
    Bmass_Crown_glm <- lm(BMloss_EAT ~ Crown_match_f,
                          data=PARAS_BMloss_DB_ALL_f)
    par(mfrow = c(2, 3))
    plot(Bmass_Crown_glm,
         ask = F, which = 1:6)
    
    Anova(Bmass_Crown_glm,
          type=3)
    summary(Bmass_Crown_glm)
    
    #.................................................................    
    #.......BM_Loss and Alarm call====
    
    #........ STATUS OF INFECTION ===

    Bmass_Acall_G_glm.0 <- lm(BMloss_EAT ~ Acall_BIN_f*INF_BIN_f,
                      data=PARAS_BMloss_DB_ALL_f)

    Anova(Bmass_Acall_G_glm.0,
      type=3) 
    
    # interaction not significant

    Bmass_Acall_G_glm.1 <- lm(BMloss_EAT ~ Acall_BIN_f+INF_BIN_f,
                          data=PARAS_BMloss_DB_ALL_f)

    Anova(Bmass_Acall_G_glm.1,
      type=3) 
    
    # We removed status of infection

    Bmass_Acall_G_glm.1 <- lm(BMloss_EAT ~ Acall_BIN_f,
                              data=PARAS_BMloss_DB_ALL_f)
    par(mfrow = c(2, 3))
    plot(Bmass_Acall_G_glm,
         ask = F, which = 1:6)
    
    Anova(Bmass_Acall_G_glm.1,
          type=3)
    summary(Bmass_Acall_G_glm.1)

  #.........................................................
  # CORRELATIONS: r parameter and associated P values  ####
    
    library(correlation)
    
    PARAS_BMloss_corrDB <- PARAS_BMloss_DB_ALL_f[,
                                                 c("T_INITIAL_log_EXPL",
                                                   "RBT_PRED_dif",
                                                   "Crown_match_f",
                                                   "Acall_BIN_f",
                                                   "BMloss_EAT")]
    # EXPLORATION ASSESSMENT TEST
    
    BM_Loss_EXPL_cor <- correlation(PARAS_BMloss_DB_ALL_f[,
                                      c("T_INITIAL_log_EXPL",
                                      "BM_Loss")],
                             include_factors = TRUE,
                             method = "pearson",
                             p_adjust = "none"
    )
    
    BM_Loss_EXPL_cor$p
 
    # RISK-TAKING ASSESSMENT TESTS
    
    #....RBT_PRED_dif
    
    PARAS_BMloss_RISKT_DB <- na.omit(PARAS_BMloss_DB_ALL_f[,
                                                           c("RBT_PRED_dif",
                                                             "BMloss_EAT")])
    
    BM_Loss_RBT_PRED_cor <- correlation(PARAS_BMloss_RISKT_DB,
                                        include_factors = TRUE,
                                        method = "pearson",
                                        p_adjust = "none"
    )
   
    #....Crown_match_f
    
    BMloss_Crown_match_f.lm <- lm(BMloss_EAT ~ Crown_match_f, 
                               data=PARAS_BMloss_DB_ALL_f)
    
    BM_Loss_Crown_matchdif <- Anova(BMloss_Crown_match_f.lm, type=3)
    
    #....Acall_BIN_f
    
    BMloss_Acall_BIN_f.lm <- lm(BMloss_EAT ~ Acall_BIN_f, 
                                  data=PARAS_BMloss_DB_ALL_f)
    
    BM_Loss_Acalldif <- Anova(BMloss_Acall_BIN_f.lm, type=3)
    
    #(3) FDR CORRECTION BY MULTIPLE TESTING ####
    
    FDR_BM_Loss_BEHAV_effect <- data.frame (BEHAV_Parameter=c("T_INITIAL_log_EXPL",
                                                        "RBT_PRED_dif",
                                                        "Crown_match_f",
                                                        "Acall_BIN_f"))
    FDR_BM_Loss_BEHAV_effect$p<- NA
    
    # We add values
    
    FDR_BM_Loss_BEHAV_effect[FDR_BM_Loss_BEHAV_effect$BEHAV_Parameter=="T_INITIAL_log_EXPL",
                       "p"] <- Bmass_Texpl_glm.anova[row.names(Bmass_Texpl_glm.anova)=="T_INITIAL_log_EXPL",
                                                     "Pr(>F)"]
    FDR_BM_Loss_BEHAV_effect[FDR_BM_Loss_BEHAV_effect$BEHAV_Parameter=="RBT_PRED_dif",
                       "p"] <- Bmass_RBT_glm.anova[row.names(Bmass_RBT_glm.anova)=="RBT_PRED_dif",
                                                     "Pr(>F)"]
    
    FDR_BM_Loss_BEHAV_effect[FDR_BM_Loss_BEHAV_effect$BEHAV_Parameter=="Crown_match_f",
                       "p"] <- BM_Loss_Crown_matchdif["Crown_match_f","Pr(>F)"]
    
    FDR_BM_Loss_BEHAV_effect[FDR_BM_Loss_BEHAV_effect$BEHAV_Parameter=="Acall_BIN_f",
                       "p"] <- BM_Loss_Acalldif["Acall_BIN_f","Pr(>F)"]
    
    FDR_BM_Loss_BEHAV_effect <- arrange(FDR_BM_Loss_BEHAV_effect, p)

    FDR_BM_Loss_BEHAV_effect$p_BDRadjust <-p.adjust(FDR_BM_Loss_BEHAV_effect$p,
                                              method="fdr")    
    
    library(kableExtra)
    FDR_BM_Loss_BEHAV_effect %>%
      kbl(digits = 4) %>%
      kable_styling()
    
    
    