#................................................................

## Analyses for manuscript: 
## Haemosporidian infections influence risk-taking behaviours in young male blackcaps Sylvia atricapilla

#................................................................

#................................................................
        # SUMMARY OF THE SCRIPT:====

        # **Cautious note**: YOU MUST FIRST RUN DATA SCRIPT
        # TO GENERATE ADEQUATED TRANSFORMED DATA FOR ANALYSES

# ...For FORAGING analyses use = MS_SYAT_BH_FORAG database ====
#................................................................

library(dplyr)
library(vegan)
library(ggplot2)
library(reshape2)
library(psych)
library(nFactors)
library(lme4)
library(lmerTest)
library(rstatix)
library(elrm)
library(PerformanceAnalytics )
library(ggpubr)
library(car)

library(lsmeans)
library(afex)
library(jtools)
library(multcomp)

set.seed(123)    
#................................................................

    ##...Foraging behaviour ANALYSES ----
##....................................................................
  
### Delete individuals without foraging behaviour

    FORAG.DATA <- !is.na(MS_SYAT_BH_FORAG$EAT_BIN)

    MS_SYAT_BH_FORAG_f <- MS_SYAT_BH_FORAG[FORAG.DATA,]

    ### Merge Behaviour and Biological database

    FORAG_BIOL_DB <- merge (MS_SYAT_BH_FORAG_f, 
                       MS_SYAT_BH_BIOL,
                       by.x ="ID", 
                       by.y="ID",
                       all.x = TRUE,
                       no.dups=FALSE)
    str(FORAG_BIOL_DB)

    ### Variable structure for analyses

    FORAG_BIOL_DB <- mutate(FORAG_BIOL_DB,
                       INF_BIN_f = as.factor(INF_BIN),
                       MULT_INF_f = as.factor(MULT_INF),
                       EAT_BIN_f = as.factor(EAT_BIN))

    ### Subset of dataframe for individuals eating and infected
    
    FORAG_BIOL_DB_EAT <- filter (FORAG_BIOL_DB,
                                   EAT_BIN ==1)%>%
        mutate(T_1_EAT_stdz = scale(T_1_EAT,
                               center=TRUE, scale =TRUE),
               N_PECKS_log_stdz = scale(N_PECKS_log,
                                  center=TRUE, scale =TRUE))

##....................................................................    
    ###...(1) Eating (YES/NO)..."EAT_BIN_f" ####

    # How many birds ate? 
    
    sum(FORAG_BIOL_DB$EAT_BIN==1) #25 birds

    # How long does it take start eating?
    
    options(digits=10)
    mean(FORAG_BIOL_DB_EAT$T_1_EAT)
    round(sd(FORAG_BIOL_DB_EAT$T_1_EAT),digits=2)
    
    # Differences between eating and not eating birds with BC
    
    BC_EAT.glm = glm(EAT_BIN_f  ~ BC, 
                     family=binomial,
                     data=FORAG_BIOL_DB)
    
    # overdispersion test
    phi <- sum((residuals(BC_EAT.glm, type="pearson"))^2)/BC_EAT.glm$df.residual
    print(c("Pearson overdispersion =", round(phi, 3)),quote=FALSE)
    
    drop1(BC_EAT.glm,~.,
          test="Chisq") # no differences
    
    #........STATUS OF INFECTION  ====
    
    E1 <- xtabs(~EAT_BIN_f+interaction(INF_BIN_f,TREAT)
                ,data=FORAG_BIOL_DB)
    
    Eat_G_glm.0 = glm(EAT_BIN_f  ~ INF_BIN_f +
                        TREAT +
                        INF_BIN_f:TREAT, 
                      family=binomial,
                      data=FORAG_BIOL_DB)
    
    # overdispersion test
    phi <- sum((residuals(Eat_G_glm.0, type="pearson"))^2)/Eat_G_glm.0$df.residual
    print(c("Pearson overdispersion =", round(phi, 3)),quote=FALSE)
    
    drop1(Eat_G_glm.0,~.,
          test="Chisq") # no significant interaction
    
    # Model without interaction
    
    Eat_G_glm = glm(EAT_BIN_f  ~ INF_BIN_f +
                      TREAT , 
                    family=binomial,
                    data=FORAG_BIOL_DB)
    
    Eat_G_LRT_db <- drop1(Eat_G_glm,
                        test="Chisq")
    Eat_G_LRT_db
    
    summary(Eat_G_glm)
    
    # overdispersion test
    phi <- sum((residuals(Eat_G_glm, type="pearson"))^2)/Eat_G_glm$df.residual
    print(c("Pearson overdispersion =", round(phi, 3)), quote=FALSE)
    
    #........MULTIPLE STATUS OF INFECTION ====
    
    ME1 <- xtabs(~EAT_BIN_f+interaction(MULT_INF_f,TREAT)
                ,data=FORAG_BIOL_DB)
    
    
    Eat_MLG_glm.0 = glm(EAT_BIN_f ~ MULT_INF_f +
                          TREAT +
                          MULT_INF_f:TREAT , 
                        family=binomial,
                        data=FORAG_BIOL_DB)
    
    # overdispersion test
    phi <- sum((residuals(Eat_MLG_glm.0, type="pearson"))^2)/Eat_MLG_glm.0$df.residual
    print(c("Pearson overdispersion =", round(phi, 3)),quote=FALSE)
    
    drop1(Eat_MLG_glm.0,~.,
          test="Chisq") 
    
    # Interaction not significant 
    
    Eat_MLG_glm = glm(EAT_BIN_f~ MULT_INF_f +
                        TREAT , 
                      family=binomial,
                      data=FORAG_BIOL_DB)
    
    Eat_MLG_LRT_db <- drop1(Eat_MLG_glm,
                           test="Chisq")
    Eat_MLG_LRT_db
    
    summary (Eat_MLG_glm)
    
    phi <- sum((residuals(Eat_MLG_glm, type="pearson"))^2)/Eat_MLG_glm$df.residual
    print(c("Pearson overdispersion =", round(phi, 3)), quote=FALSE)

## **SUBSET OF EATING BIRDS** ####

    # **Caution note**: n=25 Interaction between infection and treatment not tested  
    # because of low sample size
    
##....................................................................
    ###...(2) Latency to eat "T_1_EAT" ####
    # Sample size
    
    nrow(FORAG_BIOL_DB_EAT[!is.na(FORAG_BIOL_DB_EAT$T_1_EAT),]) #25
    
    # Correlation between BC and latency to eat
    
    T_1_EAT_BC_lm = lm(T_1_EAT_stdz  ~ BC,
                      data=FORAG_BIOL_DB_EAT)
    
    Anova(T_1_EAT_BC_lm, type=3)
    summary(T_1_EAT_BC_lm)# not significant
    
    #........STATUS OF INFECTION  ====
    
    T_1_EAT_G_lm = lm(T_1_EAT_stdz  ~ INF_BIN_f +
                            TREAT ,
                        data=FORAG_BIOL_DB_EAT)
    
    par(mfrow=c(2,3))
    plot(T_1_EAT_G_lm,
         ask = F, which = 1:6)
    dfbetasPlots(T_1_EAT_G_lm)
    
    Latency_to_eat_G_F_db <- Anova(T_1_EAT_G_lm,
                                 type=3)
    Latency_to_eat_G_F_db
    
    summary(T_1_EAT_G_lm)
    confint(T_1_EAT_G_lm)
   
    #........MULTIPLE STATUS OF INFECTION ====
    
    T_1_EAT_MLG_lm = lm(T_1_EAT_stdz ~ MULT_INF_f +
                                   TREAT , 
                               data=FORAG_BIOL_DB_EAT)
    
    Latency_to_eat_MLG_F_db <- Anova(T_1_EAT_MLG_lm,
                                 type=3)
    
    Latency_to_eat_MLG_F_db
    
    plot(T_1_EAT_MLG_lm,
         ask = F, which = 1:6)
    dfbetasPlots(T_1_EAT_MLG_lm)

    # Estimates for the plot
    # POST-HOC for testing differences between groups
    # All pairwise comparison between levels of infection
    
    confint(T_1_EAT_MLG_lm) # with no corrections
    
    T_1_EAT_MLG.Tuk_stdz <-lsmeans(T_1_EAT_MLG_lm, 
                              pairwise~MULT_INF_f,
                              level = .95, 
                              adjust = "tukey",
                              infer = c(TRUE,TRUE) )
    
    T_1_EAT_MLG.Tuk_stdz
    
    ###...(3) Number of pecks "N_PECKS_log" ####
    # Sample size
    
    nrow(FORAG_BIOL_DB_EAT[!is.na(FORAG_BIOL_DB_EAT$N_PECKS_log),])# 24
    
    xtabs(~TREAT+MULT_INF_f,
          FORAG_BIOL_DB_EAT)
    
    # Correlation between BC and number of pecks
    
    N_PECKS_EAT_BC_lm = lm(N_PECKS_log_stdz  ~ BC,
                       data=FORAG_BIOL_DB_EAT)
    
    Anova(N_PECKS_EAT_BC_lm, type=3)
    summary(N_PECKS_EAT_BC_lm)# almost not significant (P = 0.06)
    
    #........STATUS OF INFECTION  ====

    N_PECKS_log_G_lm = lm(N_PECKS_log_stdz ~ INF_BIN_f +
                            TREAT ,
                        data=FORAG_BIOL_DB_EAT)
    
    par(mfrow = c(2, 3))
    plot(N_PECKS_log_G_lm,
         ask = F, which = 1:6)
    dfbetasPlots(N_PECKS_log_G_lm)
    
    N_pecks_log_G_F_db <- Anova(N_PECKS_log_G_lm,
                                 type=3)
    N_pecks_log_G_F_db
    
    summary(N_PECKS_log_G_lm)
    confint(N_PECKS_log_G_lm)
    
    #........MULTIPLE STATUS OF INFECTION ====

    N_PECKS_log_MLG_lm = lm(N_PECKS_log_stdz~ MULT_INF_f +
                                   TREAT ,
                               data=FORAG_BIOL_DB_EAT)
    par(mfrow = c(2, 3))
    plot(N_PECKS_log_MLG_lm,
         ask = F, which = 1:6)
    dfbetasPlots(N_PECKS_log_MLG_lm)
    
    N_pecks_log_MLG_F_db <- Anova(N_PECKS_log_MLG_lm,
                              type=3) 
   
    N_pecks_log_MLG_F_db
    
    # Estimates for the plot
    # POST-HOC for testing differences between groups
    # All pairwise comparison between levels of infection
    
    confint(N_PECKS_log_MLG_lm) # with no corrections
    
    N_PECKS_log_MLG.Tuk_stdz <-lsmeans(N_PECKS_log_MLG_lm, 
                                pairwise~MULT_INF_f,
                                level = .95, 
                                adjust = "tukey",
                                infer = c(TRUE,TRUE) )
    
    N_PECKS_log_MLG.Tuk_stdz
    
    ###...FDR CORRECTION BY MULTIPLE TESTING OF FORAGING VARIABLES  ----
    
    Eat_G_LRT_db$FORAG_var <- "EAT_BIN"
    Eat_G_LRT_db$Predictor_var <- row.names(Eat_G_LRT_db)
    Eat_MLG_LRT_db$FORAG_var <- "EAT_BIN"
    Eat_MLG_LRT_db$Predictor_var <- row.names(Eat_MLG_LRT_db)
    
    Latency_to_eat_G_F_db$FORAG_var <- "Latency_to_eat"
    Latency_to_eat_G_F_db$Predictor_var <- row.names(Latency_to_eat_G_F_db)
    Latency_to_eat_MLG_F_db$FORAG_var <- "Latency_to_eat"
    Latency_to_eat_MLG_F_db$Predictor_var <- row.names(Latency_to_eat_MLG_F_db)
    
    N_pecks_log_G_F_db$FORAG_var <- "N_pecks_log"
    N_pecks_log_G_F_db$Predictor_var <- row.names(N_pecks_log_G_F_db)
    N_pecks_log_MLG_F_db$FORAG_var <- "N_pecks_log"
    N_pecks_log_MLG_F_db$Predictor_var <- row.names(N_pecks_log_MLG_F_db)
    
    
    FDR_FORAG_tests_p <- data.frame( FORAG_var = 
                                     c(rep("EAT_BIN",4),
                                       rep("Latency_to_eat",4),
                                       rep("N_pecks_log",4)),
                                     Predictor_var =
                                     c(rep(c("INF_BIN_f", "TREAT.1",
                                       "MULT_INF_f","TREAT.2"),3))
    )
    
    FDR_FORAG_tests_p$Predictor_var_FORAG_var <- paste(FDR_FORAG_tests_p$FORAG_var,
                                                       FDR_FORAG_tests_p$Predictor_var,
                                                       sep="_")
    unique(FDR_FORAG_tests_p$Predictor_var_FORAG_var)
    FDR_FORAG_tests_p$p <- NA
    
    # for binomial eat
    FDR_FORAG_tests_p[FDR_FORAG_tests_p$Predictor_var_FORAG_var== "EAT_BIN_INF_BIN_f",
                      "p"]<- Eat_G_LRT_db[Eat_G_LRT_db$Predictor_var=="INF_BIN_f",
                                        "Pr(>Chi)"]
    FDR_FORAG_tests_p[FDR_FORAG_tests_p$Predictor_var_FORAG_var== "EAT_BIN_TREAT.1",
                      "p"]<- Eat_G_LRT_db[Eat_G_LRT_db$Predictor_var=="TREAT",
                                        "Pr(>Chi)"]
    FDR_FORAG_tests_p[FDR_FORAG_tests_p$Predictor_var_FORAG_var== "EAT_BIN_MULT_INF_f",
                      "p"]<- Eat_MLG_LRT_db[Eat_MLG_LRT_db$Predictor_var=="MULT_INF_f",
                                        "Pr(>Chi)"]
    FDR_FORAG_tests_p[FDR_FORAG_tests_p$Predictor_var_FORAG_var== "EAT_BIN_TREAT.2",
                      "p"]<- Eat_MLG_LRT_db[Eat_MLG_LRT_db$Predictor_var=="TREAT",
                                        "Pr(>Chi)"]
    # for latency to eat
    FDR_FORAG_tests_p[FDR_FORAG_tests_p$Predictor_var_FORAG_var== "Latency_to_eat_INF_BIN_f",
                      "p"]<- Latency_to_eat_G_F_db[Latency_to_eat_G_F_db$Predictor_var=="INF_BIN_f",
                                        "Pr(>F)"]
    FDR_FORAG_tests_p[FDR_FORAG_tests_p$Predictor_var_FORAG_var== "Latency_to_eat_TREAT.1",
                      "p"]<- Latency_to_eat_G_F_db[Latency_to_eat_G_F_db$Predictor_var=="TREAT",
                                        "Pr(>F)"]
    FDR_FORAG_tests_p[FDR_FORAG_tests_p$Predictor_var_FORAG_var== "Latency_to_eat_MULT_INF_f",
                      "p"]<- Latency_to_eat_MLG_F_db[Latency_to_eat_MLG_F_db$Predictor_var=="MULT_INF_f",
                                           "Pr(>F)"]
    FDR_FORAG_tests_p[FDR_FORAG_tests_p$Predictor_var_FORAG_var== "Latency_to_eat_TREAT.2",
                      "p"]<- Latency_to_eat_MLG_F_db[Latency_to_eat_MLG_F_db$Predictor_var=="TREAT",
                                           "Pr(>F)"]
    
    # for number of pecks
    
    FDR_FORAG_tests_p[FDR_FORAG_tests_p$Predictor_var_FORAG_var== "N_pecks_log_INF_BIN_f",
                      "p"]<- N_pecks_log_G_F_db[N_pecks_log_G_F_db$Predictor_var=="INF_BIN_f",
                                                 "Pr(>F)"]
    FDR_FORAG_tests_p[FDR_FORAG_tests_p$Predictor_var_FORAG_var== "N_pecks_log_TREAT.1",
                      "p"]<- N_pecks_log_G_F_db[N_pecks_log_G_F_db$Predictor_var=="TREAT",
                                                 "Pr(>F)"]
    FDR_FORAG_tests_p[FDR_FORAG_tests_p$Predictor_var_FORAG_var== "N_pecks_log_MULT_INF_f",
                      "p"]<- N_pecks_log_MLG_F_db[N_pecks_log_MLG_F_db$Predictor_var=="MULT_INF_f",
                                                    "Pr(>F)"]
    FDR_FORAG_tests_p[FDR_FORAG_tests_p$Predictor_var_FORAG_var== "N_pecks_log_TREAT.2",
                      "p"]<- N_pecks_log_MLG_F_db[N_pecks_log_MLG_F_db$Predictor_var=="TREAT",
                                                    "Pr(>F)"]
    
    
    # Arrange database by p values
    
    FDR_FORAG_tests_p <- arrange(FDR_FORAG_tests_p, 
                                 p)
    
    # FDR adjustment of p values
    
    FDR_FORAG_tests_p$FDR_p_adjust <- p.adjust(FDR_FORAG_tests_p$p,
                                             method="fdr")
    
    FDR_FORAG_tests_p$bonferroni_p_adjust <- p.adjust(FDR_FORAG_tests_p$p,
                                               method="bonferroni")
  