#................................................................

## Analyses for manuscript: 
## Haemosporidian infections influence risk-taking behaviours in young male blackcaps Sylvia atricapilla

#................................................................
        # SUMMARY OF THE SCRIPT: ====

    # **Cautious note**: YOU MUST FIRST RUN DATA SCRIPT
# TO GENERATE ADEQUATED TRANSFORMED DATA FOR ANALYSES

# ...For EXPLORATION analyses use = MS_SYAT_BH_EXPL database ====
#................................................................

library(dplyr)
library(vegan)
library(ggplot2)
library(officer)
library(Rmisc)
library(rvg)
library(tidyverse)
library(here)
library(reshape2)
library(psych)
library(nFactors)
library(lme4)
library(lmerTest)
library(rstatix)
library(PerformanceAnalytics )
library(ggpubr)
library(emmeans)
library(afex)

# robust estimates
        library(car)
        library(MASS)
        library(robustbase)

set.seed(123)     
#................................................................

##...Exploratory behaviour ANALYSES ----
#.......................................

    ### Delete individuals without exploratory behaviour measured

    EXPL.DATA <- !is.na(MS_SYAT_BH_EXPL$T_INITIAL)
    
    MS_SYAT_BH_EXPL_f <- MS_SYAT_BH_EXPL[EXPL.DATA,]

    ### Merge Behaviour and Biological database

    EXPL_BIOL_DB <- merge (MS_SYAT_BH_EXPL_f, 
                           MS_SYAT_BH_BIOL,
                           by.x ="ID", 
                           by.y="ID",
                           all.x = TRUE,
                           no.dups=FALSE)
    str(EXPL_BIOL_DB)

    ### Variable structure for analyses
    
    EXPL_BIOL_DB <- mutate(EXPL_BIOL_DB,
                           INF_BIN_f = as.factor(INF_BIN),
                           MULT_INF = as.factor(MULT_INF))
    
    ### Exclusion of one individual that did not performed any movement
    ### during exploration.
    ### Results did not change qualitatively but model residuals clearly improved
    
    Bird.not.move <-EXPL_BIOL_DB[EXPL_BIOL_DB$N_MOV==0, "ID"]# Bird31
    
    EXPL_BIOL_DB_out <- filter(EXPL_BIOL_DB,
                               ID != Bird.not.move)
    
    ### Standardized variables for standardized estimates

    EXPL_BIOL_DB_out <- mutate(EXPL_BIOL_DB_out,
                               T_INITIAL_log_stdz = as.numeric(scale(T_INITIAL_log,
                                                          center=TRUE,
                                                          scale= TRUE)),
                               T_CROSS_log_stdz = as.numeric(scale(T_CROSS_log,
                                                        center=TRUE,
                                                        scale= TRUE)),
                               Activity_stdz = as.numeric(scale(Activity,
                                                     center=TRUE,
                                                     scale = TRUE)),
                               RBT_stdz = as.numeric(scale(RBT,
                                                center = TRUE,
                                                scale = TRUE)))

    # Correlation between BC and exploratory behavioural variables
    
    names(EXPL_BIOL_DB_out)
    
    chart.Correlation(EXPL_BIOL_DB_out[, c("BC",
                                           "T_INITIAL_log_stdz",
                                           "T_CROSS_log_stdz",  
                                           "Activity_stdz",
                                           "RBT_stdz")])
    corr.test(EXPL_BIOL_DB_out[, c("BC",
                                   "T_INITIAL_log_stdz",
                                   "T_CROSS_log_stdz",  
                                   "Activity_stdz",
                                   "RBT_stdz")])
    
    ###...(1) Activity level..."Activity_stdz" ####
    
    # Previous normality test
    EXPL_BIOL_DB_out %>%
      group_by(INF_BIN_f) %>%
      shapiro_test(Activity_stdz) 
    
    ggqqplot(EXPL_BIOL_DB_out, "Activity_stdz")
    # Levene's test
    
    EXPL_BIOL_DB_out %>%
      levene_test(Activity_stdz ~ INF_BIN_f,
                  center=mean) 
    
    EXPL_BIOL_DB_out %>%
      levene_test(Activity_stdz ~ TREAT,
                  center=mean) 
    
    EXPL_BIOL_DB_out %>%
      levene_test(Activity_stdz ~ TREAT*INF_BIN_f,
                  center=mean) 
    
    EXPL_BIOL_DB_out %>%
      levene_test(Activity_stdz ~ MULT_INF,
                  center=mean) 
    
    #........ STATUS OF INFECTION  ====
    
    Activ_G_EXPL_lm.0 <- lm(Activity_stdz ~ INF_BIN_f*TREAT,
                          data=EXPL_BIOL_DB_out)
    
    par(mfrow = c(2, 3))
    plot(Activ_G_EXPL_lm.0,
         ask = F, which = 1:6)
    par(mfrow = c(1, 1))
    dfbetasPlots(Activ_G_EXPL_lm.0)
    
    Anova(Activ_G_EXPL_lm.0,
          type=3)
    
    # Interaction not significant
    
    Activ_G_EXPL_lm <- lm(Activity_stdz ~ INF_BIN_f+TREAT,
                        data=EXPL_BIOL_DB_out)
    par(mfrow = c(2, 3))
    plot(Activ_G_EXPL_lm,
         ask = F, which = 1:6)
    dfbetasPlots(Activ_G_EXPL_lm)
    
    Activ_G_EXPL_F <- Anova(Activ_G_EXPL_lm,
                          type=3)
    
    Activ_G_EXPL_F_db <- data.frame(Activ_G_EXPL_F)  # database with F and P for FDR correction
    
    summary(Activ_G_EXPL_lm,
            type=3)
    confint(Activ_G_EXPL_lm)
 
    #........ MULTIPLE STATUS OF INFECTION  ====
    
    Activ_MLG_EXPL_lm_0 <- lm(Activity_stdz ~ MULT_INF*TREAT,
                             data=EXPL_BIOL_DB_out)
    par(mfrow = c(2, 3))
    plot(Activ_MLG_EXPL_lm_0,
         ask = F, which = 1:6) 
    dfbetasPlots(Activ_MLG_EXPL_lm_0)
    
    Anova(Activ_MLG_EXPL_lm_0,
          type=3)
    
    # Interaction not significant
    
    Activ_MLG_EXPL_lm <- lm(Activity_stdz ~ MULT_INF+TREAT,
                           data=EXPL_BIOL_DB_out)
    par(mfrow = c(2, 3))
    plot(Activ_MLG_EXPL_lm,
         ask = F, which = 1:6)
    dfbetasPlots(Activ_MLG_EXPL_lm)
    
    Activ_MLG_EXPL_F<- Anova(Activ_MLG_EXPL_lm,
                            type=3)
    
    Activ_MLG_EXPL_F_db <- data.frame(Activ_MLG_EXPL_F)  # database with F and p for FDR correction
    
    summary (Activ_MLG_EXPL_lm)
    
    # Estimates of the plot in standard units
    
    confint(Activ_MLG_EXPL_lm) # with no corrections
    
    Activ_MLG_EXPL.Tuk_stdz <-lsmeans(Activ_MLG_EXPL_lm, 
                                        pairwise~MULT_INF,
                                        level = .95, 
                                        adjust = "tukey",
                                        infer = c(TRUE,TRUE) )
    
    Activ_MLG_EXPL.Tuk_stdz 
    
    ###...(2) Relative body turns..."RBT_stdz"####
    
    # Previous normality test
    
    EXPL_BIOL_DB_out %>%
      group_by(INF_BIN_f) %>%
      shapiro_test(RBT_stdz) 
    
    ggqqplot(EXPL_BIOL_DB_out, "RBT_stdz")
    
    # Levene's test
    
    EXPL_BIOL_DB_out %>%
      levene_test(RBT_stdz ~ INF_BIN_f,
                  center=mean) 
    
    EXPL_BIOL_DB_out %>%
      levene_test(RBT_stdz ~ TREAT,
                  center=mean) 
    
    EXPL_BIOL_DB_out %>%
      levene_test(RBT_stdz~ TREAT*INF_BIN_f,
                  center=mean) 
    
    EXPL_BIOL_DB_out %>%
      levene_test(RBT_stdz~ MULT_INF,
                  center=mean) 
    
    #........ STATUS OF INFECTION  ====
    
    RBT_G_EXPL_lm.0 <- lm(RBT_stdz~ INF_BIN_f*TREAT,
                          data=EXPL_BIOL_DB_out)
    
    par(mfrow = c(2, 3))
    plot( RBT_G_EXPL_lm.0,
          ask = F, which = 1:6)
    dfbetasPlots(RBT_G_EXPL_lm.0)
    
    Anova(RBT_G_EXPL_lm.0,
          type=3)
    
    # Interaction not significant
    
    RBT_G_EXPL_lm <- lm(RBT_stdz~ INF_BIN_f+TREAT,
                      data=EXPL_BIOL_DB_out)
    par(mfrow = c(2, 3))
    plot(RBT_G_EXPL_lm,
         ask = F, which = 1:6)
    
    RBT_G_EXPL_F <- Anova(RBT_G_EXPL_lm,
                        type=3)
    
    RBT_G_EXPL_F_db <- data.frame(RBT_G_EXPL_F)  # database with F and P for FDR correction
    
    summary(RBT_G_EXPL_lm)
    confint(RBT_G_EXPL_lm)
    
    #........ MULTIPLE STATUS OF INFECTION ====
    
    RBT_MLG_EXPL_lm.0 <- lm(RBT_stdz~ MULT_INF*TREAT,
                           data= EXPL_BIOL_DB_out)
    
    par(mfrow = c(2, 3))
    plot(RBT_MLG_EXPL_lm.0,
         ask = F, which = 1:6)
    dfbetasPlots(RBT_MLG_EXPL_lm.0)
    
    Anova(RBT_MLG_EXPL_lm.0,
          type=3)
    
    # Interaction not significant
    
    RBT_MLG_EXPL_lm <- lm(RBT_stdz~ MULT_INF+TREAT,
                         data=EXPL_BIOL_DB_out)
    
    par(mfrow = c(2, 3))
    plot(RBT_MLG_EXPL_lm,
         ask = F, which = 1:6) 
    dfbetasPlots(RBT_MLG_EXPL_lm)
    
    RBT_MLG_EXPL_F <- Anova(RBT_MLG_EXPL_lm,
                           type=3)
    
    RBT_MLG_EXPL_F_db <- data.frame(RBT_MLG_EXPL_F)  # database with F and p for FDR correction
    
    ## Estimates of the plot in standard units
    
    summary(RBT_MLG_EXPL_lm,
            type=3)
    
    confint(RBT_MLG_EXPL_lm) # with no corrections
    
    RBT_MLG_EXPL.Tuk_stdz <-lsmeans(RBT_MLG_EXPL_lm, 
                                     pairwise~MULT_INF,
                                     level = .95, 
                                     adjust = "tukey",
                                     infer = c(TRUE,TRUE) )
    
    RBT_MLG_EXPL.Tuk_stdz
    
    ###...(3) Initial latency..."T_INITIAL_log_stdz" ####
    
    # Previous normality test
    
    EXPL_BIOL_DB_out %>%
        group_by(INF_BIN_f) %>%
        shapiro_test(T_INITIAL_log_stdz) 
    
    ggqqplot(EXPL_BIOL_DB_out, "T_INITIAL_log_stdz")

    # Levene's test
    
    EXPL_BIOL_DB_out %>%
        levene_test(T_INITIAL_log_stdz ~ INF_BIN_f) 
    
    EXPL_BIOL_DB_out %>%
        levene_test(T_INITIAL_log_stdz ~ TREAT,
                    center=mean) 
    
    EXPL_BIOL_DB_out %>%
        levene_test(T_INITIAL_log_stdz ~ TREAT*INF_BIN_f,
                    center=mean) 
    
    EXPL_BIOL_DB_out %>%
        levene_test(T_INITIAL_log_stdz ~ MULT_INF,
                    center=mean) 
    
    #........ STATUS OF INFECTION ====
    
    T1_G_EXPL_lm.0 <- lm(T_INITIAL_log_stdz ~ INF_BIN_f*TREAT,
                             data=EXPL_BIOL_DB_out)
    Anova(T1_G_EXPL_lm.0,
          type=3)
    
    par(mfrow = c(2, 3))
    plot(T1_G_EXPL_lm.0,
         ask = F, which = 1:6)
    dfbetasPlots(T1_G_EXPL_lm.0)
    
    # Not significant interaction
    
    T1_G_EXPL_lm <- lm(T_INITIAL_log_stdz ~ INF_BIN_f+TREAT,
                           data=EXPL_BIOL_DB_out)
    par(mfrow = c(2, 3))
    plot(T1_G_EXPL_lm,
         ask = F, which = 1:6)
    dfbetasPlots(T1_G_EXPL_lm)
    
    T1_G_EXPL_F <- Anova(T1_G_EXPL_lm,
                       type=3)
    
    T1_G_EXPL_F_db <- data.frame(T1_G_EXPL_F)  # database with F and P for FDR correction
    
    summary(T1_G_EXPL_lm)
    confint(T1_G_EXPL_lm)
    
    # Plot of the effect
    
    ggplot(data= EXPL_BIOL_DB_out,
           aes(y=T_INITIAL_log_stdz, 
               x=INF_BIN_f ))  +
        geom_boxplot(aes(fill=INF_BIN_f,
                         alpha=0.5)) + 
        geom_point(aes(alpha=0.5)) +
        theme_classic()+
        theme ( legend.position = "none")

    #........ MULTIPLE STATUS OF INFECTION  ====

    T1_EXPL_MLG_lm.0 <- lm(T_INITIAL_log_stdz ~ MULT_INF*TREAT,
                                data=EXPL_BIOL_DB_out)
    
    Anova(T1_EXPL_MLG_lm.0,
          type=3)
    
    par(mfrow = c(2, 3))
    plot(T1_EXPL_MLG_lm.0,
         ask = F, which = 1:6)
    dfbetasPlots(T1_EXPL_MLG_lm.0)
    
    # Not significant interaction
  
    T1_EXPL_MLG_lm <- lm(T_INITIAL_log_stdz ~ MULT_INF+TREAT,
                              data=EXPL_BIOL_DB_out)
   par(mfrow = c(2, 3))
    plot(T1_EXPL_MLG_lm,
         ask = F, which = 1:6)
    dfbetasPlots(T1_EXPL_MLG_lm)
    
    T1_EXPL_MLG_F <- Anova(T1_EXPL_MLG_lm,
                           type=3)
    
    T1_EXPL_MLG_F_db <- data.frame(T1_EXPL_MLG_F)  # database with F and p for FDR correction
   
    # POST-HOC for testing differences between groups
    # Pairwise comparison within levels of infection
    
    confint(T1_EXPL_MLG_lm) # with no corrections
   
    T1_EXPL_MLG.Tuk_stdz <-lsmeans(T1_EXPL_MLG_lm, 
                                 pairwise~MULT_INF,
                                 level = .95, 
                                 adjust = "tukey",
                                 infer = c(TRUE,TRUE) )
    
    T1_EXPL_MLG.Tuk_stdz

    # Plot with the effects
    
    ggplot(data= filter(EXPL_BIOL_DB_out,
                        !is.na(MULT_INF)),
           aes(y=T_INITIAL_log_stdz, 
               x=MULT_INF))  +
        geom_boxplot(aes(fill=MULT_INF,
                         alpha=0.5)) + 
        geom_point(aes(alpha=0.5)) +
        theme_classic()+
        theme ( legend.position = "none")

    ###...(4) Latency to cross..."T_CROSS_log_stdz"####

    #     Previous normality test
    
    EXPL_BIOL_DB_out %>%
        group_by(INF_BIN_f) %>%
        shapiro_test(T_CROSS_log_stdz) # not normal for infected individuals. 
    
    ggqqplot(EXPL_BIOL_DB_out, "T_CROSS_log_stdz")# some individuals never crossed
    
    # frequency of individuals that did not crossed
    
    length(EXPL_BIOL_DB_out[EXPL_BIOL_DB_out$T_CROSS ==300,
                            "ID"]) # 9 out of 40
    # **Caution note**: we will check carefully normality of residuals in
    # the models
    
    # Levene's test
    
    EXPL_BIOL_DB_out %>%
        levene_test(T_CROSS_log_stdz ~ INF_BIN_f,
                    center=mean) 
    
    EXPL_BIOL_DB_out %>%
        levene_test(T_CROSS_log_stdz ~ TREAT,
                    center=mean) 
    
    EXPL_BIOL_DB_out %>%
        levene_test(T_CROSS_log_stdz ~ TREAT*INF_BIN_f,
                    center=mean) 

    EXPL_BIOL_DB_out %>%
        levene_test(T_CROSS_log_stdz ~ MULT_INF,
                    center=mean) 
    
    #........ STATUS OF INFECTION  ====
    
    Lat_X_G_EXPL_lm.0 <- lm(T_CROSS_log_stdz ~ INF_BIN_f*TREAT,
                            data=EXPL_BIOL_DB_out)
    
    par(mfrow = c(2, 3))
    plot( Lat_X_G_EXPL_lm.0,
          ask = F, which = 1:6) 
    dfbetasPlots(Lat_X_G_EXPL_lm.0)
    
    Anova(Lat_X_G_EXPL_lm.0,
          type=3)
    
    # Interaction not significant
    
    Lat_X_G_EXPL_lm <- lm(T_CROSS_log_stdz ~ INF_BIN_f+TREAT,
                          data=EXPL_BIOL_DB_out)
    par(mfrow = c(2, 3))
    plot(Lat_X_G_EXPL_lm,
         ask = F, which = 1:6)
    dfbetasPlots(Lat_X_G_EXPL_lm)
    
    Lat_X_G_EXPL_F <- Anova(Lat_X_G_EXPL_lm,
                          type=3)
   
    Lat_X_G_EXPL_F_db <- data.frame(Lat_X_G_EXPL_F)  # database with F and P for FDR correction
    
    summary(Lat_X_G_EXPL_lm)
    confint(Lat_X_G_EXPL_lm)
    
    #........ MULTIPLE STATUS OF INFECTION  ====
    
    Lat_X_MLG_EXPL_lm.0 = lm(T_CROSS_log_stdz  ~ MULT_INF*TREAT,
                             data=EXPL_BIOL_DB_out)
    par(mfrow = c(2, 3))
    plot(Lat_X_MLG_EXPL_lm.0,
         ask = F, which = 1:6)
    dfbetasPlots(Lat_X_MLG_EXPL_lm.0)
    
    Anova(Lat_X_MLG_EXPL_lm.0,
          type=3)
    
    # Interaction not significant
    
    Lat_X_MLG_EXPL_lm = lm(T_CROSS_log_stdz ~ MULT_INF +
                               TREAT ,
                           data=EXPL_BIOL_DB_out)
    par(mfrow = c(2, 3))
    plot(Lat_X_MLG_EXPL_lm,
         ask = F, which = 1:6)
    dfbetasPlots(Lat_X_MLG_EXPL_lm)
    
    Lat_X_MLG_EXPL_F <- Anova(Lat_X_MLG_EXPL_lm,
                             type=3)
    
    Lat_X_MLG_EXPL_F_db <- data.frame(Lat_X_MLG_EXPL_F)  # database with F and P for FDR correction
    
    ## Estimates of the plot in standard units
    
    summary(Lat_X_MLG_EXPL_lm,
            type=3)
    
    confint(Lat_X_MLG_EXPL_lm) # with no corrections
    
    Lat_X_MLG_EXPL.Tuk_stdz <-lsmeans(Lat_X_MLG_EXPL_lm, 
                                   pairwise~MULT_INF,
                                   level = .95, 
                                   adjust = "tukey",
                                   infer = c(TRUE,TRUE) )
    
    Lat_X_MLG_EXPL.Tuk_stdz

    ###...(5) SPEED OF EXPLORATION: latency to displace to each ####
    # new cage position      *####

    Long_T6_EXPL_var <- list(c("T_INITIAL_log",
                             "T_SECOND_log",
                             "T_THIRD_log",
                             "T_FOURTH_log",
                             "T_FIFTH_log",
                             "T_SIXTH_log"))
    
    # Transformation from wide to long format:

    EXPL_BIOL_DB_out_Long <-reshape(EXPL_BIOL_DB_out, 
                                 direction='long', 
                                 varying=Long_T6_EXPL_var, 
                                 timevar="Position",
                                 times=c("T_1", "T_2", 
                                         "T_3", "T_4", 
                                         "T_5", "T_6"),
                                 v.names="T_EXPL_log",
                                 idvar="ID") 
    
    #     Previous normality test
    EXPL_BIOL_DB_out_Long <- mutate(EXPL_BIOL_DB_out_Long,
                                    Position_f = as.factor(Position)
    )
    
    EXPL_BIOL_DB_out_Long %>%
        group_by(Position_f) %>%
        shapiro_test(T_EXPL_log)# bad for T5 and T6,many individuals
                # did not reach 6 positions
    
    xtabs(~T_EXPL_log + Position_f, 
          EXPL_BIOL_DB_out_Long  )
    
    ggqqplot(EXPL_BIOL_DB_out_Long, 
             "T_EXPL_log",
             facet.by = "Position")
    
    # Repeated measures analyses
    
    EXPL_BIOL_DB_out_Long <- mutate (EXPL_BIOL_DB_out_Long,
                                     T_EXPL_log_stdz = scale(T_EXPL_log,
                                                             center=TRUE, scale=TRUE))
    
    #........ STATUS OF INFECTION  ====
    
    Lat_pos_G_EXPL_0 <-  aov_car(T_EXPL_log_stdz~ Position_f*TREAT*INF_BIN_f+
                                        Error(ID/Position_f),
                                    data=EXPL_BIOL_DB_out_Long,
                                    type=3)
    Lat_pos_G_EXPL_0
   
    # Interaction not significant
    
    Lat_pos_G_EXPL <- aov_car(T_EXPL_log_stdz~ Position_f +
                                        TREAT+
                                        INF_BIN_f+
                                        Error(ID/Position_f),
                                    data=EXPL_BIOL_DB_out_Long,
                                    type=3)
    Lat_pos_G_EXPL
    
    Lat_pos_G_EXPL_F <- Lat_pos_G_EXPL$anova_table
    
    Lat_pos_G_EXPL_F_db <- data.frame(Lat_pos_G_EXPL_F)  # database with F and p for FDR correction
    
    summary(Lat_pos_G_EXPL)
    
    # Residuals plot
    plots <- list() 
    
    for (i in c(1:6)){
        p1 = ggqqplot(as.numeric(residuals(Lat_pos_G_EXPL$lm)[,i]), 
                       title=paste(colnames(residuals(Lat_pos_G_EXPL$lm))[i]))
        plots[[i]] <- p1 
    }
    
    multiplot(plotlist = plots, 
              cols = 3)
    
    # POST-HOC 
    
    # Within position between level of infection
    
    Lat_pos_G_EXPL.Bonf <-lsmeans(Lat_pos_G_EXPL, 
                                 pairwise~INF_BIN_f|Position_f ,
                                 level = .95, 
                                 adjust = "bonferroni",
                                 infer = c(TRUE,TRUE) )
    
    #adjusting for all tests together
    summary(Lat_pos_G_EXPL.Bonf,
            by=NULL,adjust="Bonferroni")
    
     # Plot of the results   
    Lat_pos_G_plot.1 <- afex_plot(Lat_pos_G_EXPL, x="Position_f",
                                    trace="INF_BIN_f",
                                    data_arg= list(cex=2.5, pch=21),
                                    mapping = c("shape", "fill"),
                                within_vars="ID")+
        theme_classic() + 
        theme(legend.position="bottom")
    Lat_pos_G_plot.1
    
    EXPL_T1.T6.PLOT <- ggplot(data=EXPL_BIOL_DB_out_Long, 
                                aes(y=T_EXPL_log,
                                    x=Position_f,
                                    group=ID))+ 
        theme_classic()  +
        geom_point(shape=16, size=2,
                   alpha=0.8, aes(colour=factor(INF_BIN_f), 
                                  fill = factor(INF_BIN_f)))+
      scale_fill_manual(values=c("#332288","#DDCC77"))+
      scale_colour_manual(values=c("#332288","#DDCC77"))+
        geom_line(size=1,
                  aes(col=INF_BIN_f,
                      alpha=0.7))+
      annotate("text", x = 1, y = 2.7, 
               label = "**", size=7)+
      annotate("text", x = 2, y = 2.7, 
               label = "**", size=7)+
      xlab ("Positions")+
      ylab ("Latency time (log seconds)")+
      theme(axis.text.x=element_text(size=14,
                                     vjust = -0.5,
                                     colour="black"),
            axis.text.y=element_text(size=14,
                                     hjust = -0.5,
                                     colour="black"),
            axis.title.x=element_text(size=14,
                                      vjust=-3,
                                      colour="black"),
            axis.title.y = element_text(size=14,
                                        vjust=3,
                                        colour="black"),
            axis.ticks = element_line(colour="black"))+
      scale_y_continuous(breaks=c(seq(0,2.5, by=0.5)))+
      scale_x_discrete(labels=c("T_1" = "T1", 
                                "T_2" = "T2",
                                "T_3" = "T3",
                                "T_4" = "T4",
                                "T_5" = "T5",
                                "T_6" = "T6"))
    
    Fig3 <- EXPL_T1.T6.PLOT
    Fig3
    
    Fig3.nolegend <- Fig3+
      theme(legend.position = "none",
            strip.background = element_blank(),
            strip.text.x = element_blank())
    Fig3.nolegend
    
    # Code to export figure 3 in pptx format:
    
    Fig3_dml <- rvg::dml(ggobj = Fig3.nolegend)
  
      # initialize PowerPoint slide
    officer::read_pptx() %>%
      # add slide
    officer::add_slide() %>%
      # specify object and location of object
    officer::ph_with(Fig3_dml,
                     ph_location(height = 4.2, width = 5)) %>%
      # export slide
    base::print(
      target = here::here(
        "Fig3.pptx"
      )
    )
    
    #........ MULTIPLE STATUS OF INFECTION ====
    
    Lat_pos_MLG_EXPL_0 <-  aov_car(T_EXPL_log~ Position_f*TREAT*MULT_INF+
                                     Error(ID/Position_f),
                                 data=EXPL_BIOL_DB_out_Long,
                                 type=3)
    Lat_pos_MLG_EXPL_0
    
    # Interaction not significant
    
    Lat_pos_MLG_EXPL <- aov_car(T_EXPL_log~ Position_f +
                                TREAT+
                                MULT_INF+
                                Error(ID/Position_f),
                              data=EXPL_BIOL_DB_out_Long,
                              type=3)
    Lat_pos_MLG_EXPL
    
    Lat_pos_MLG_EXPL_F <- Lat_pos_MLG_EXPL$anova_table
    
    Lat_pos_MLG_EXPL_F_db <- data.frame(Lat_pos_MLG_EXPL_F) # database with F and P for FDR correction
    
    summary(Lat_pos_MLG_EXPL)
    
    # Residuals plot
    
    plots <- list()  
    for (i in c(1:6)){
        p1 = ggqqplot(as.numeric(residuals(Lat_pos_MLG_EXPL$lm)[,i]), 
                       title=paste(colnames(residuals(Lat_pos_MLG_EXPL$lm))[i]))
        plots[[i]] <- p1 
    }
    
    multiplot(plotlist = plots, 
              cols = 3)
    
    # POST-HOC 
    
    # Within position between level of infection
    
    Lat_pos_MLG_EXPL.Bonf <-lsmeans(Lat_pos_MLG_EXPL, 
                                  pairwise~MULT_INF|Position_f ,
                                  level = .95, 
                                  adjust = "Bonferroni",
                                  infer = c(TRUE,TRUE) )
    Lat_pos_MLG_EXPL.Bonf
    
    #adjusting for all tests together
    summary(Lat_pos_MLG_EXPL.Bonf,
            by=NULL,adjust="bonferroni") # correction by 18 tests

    # Plot of the results   
    Lat_pos_MLG_plot.1 <- afex_plot(Lat_pos_MLG_EXPL, x="Position_f",
                                trace="MULT_INF",
                                data_arg= list(cex=2.5, pch=21),
                                mapping = c("shape", "fill"))+
        theme_classic() + 
        theme(legend.position="bottom")
    Lat_pos_MLG_plot.1

    
    ###...FDR CORRECTION BY MULTIPLE TESTING DURING THE EXPLORATORY TEST  ----
    
    # (1) lm results
    T1_G_EXPL_F_db$EXPL_var <- "Initial_Latency_log"
    T1_G_EXPL_F_db$Predictor_var <- row.names(T1_G_EXPL_F_db)
    T1_EXPL_MLG_F_db$EXPL_var <- "Initial_Latency_log"
    T1_EXPL_MLG_F_db$Predictor_var <- row.names(T1_EXPL_MLG_F_db)
    
    Lat_X_G_EXPL_F_db$EXPL_var <- "Latency_to_cross_log"
    Lat_X_G_EXPL_F_db$Predictor_var <- row.names(Lat_X_G_EXPL_F_db)
    Lat_X_MLG_EXPL_F_db$EXPL_var <- "Latency_to_cross_log"
    Lat_X_MLG_EXPL_F_db$Predictor_var <- row.names(Lat_X_MLG_EXPL_F_db)
    
    Activ_G_EXPL_F_db$EXPL_var <- "Activity"
    Activ_G_EXPL_F_db$Predictor_var <- row.names(Activ_G_EXPL_F_db)
    Activ_MLG_EXPL_F_db$EXPL_var <- "Activity"
    Activ_MLG_EXPL_F_db$Predictor_var <- row.names(Activ_MLG_EXPL_F_db)
    
    RBT_G_EXPL_F_db$EXPL_var <- "Relative_body_turns"
    RBT_G_EXPL_F_db$Predictor_var <- row.names(RBT_G_EXPL_F_db)
    RBT_MLG_EXPL_F_db$EXPL_var <- "Relative_body_turns"
    RBT_MLG_EXPL_F_db$Predictor_var <- row.names(RBT_MLG_EXPL_F_db)
    
    # Merge all results from linear model testing
    
    EXPL_lm_allvars.0 <- rbind(T1_G_EXPL_F_db,
                             T1_EXPL_MLG_F_db,
                             Lat_X_G_EXPL_F_db,
                             Lat_X_MLG_EXPL_F_db,
                             Activ_G_EXPL_F_db,
                             Activ_MLG_EXPL_F_db,
                             RBT_G_EXPL_F_db,
                             RBT_MLG_EXPL_F_db)
    
    unique(EXPL_lm_allvars.0$Predictor_var)
    
    EXPL_lm_allvars <- filter(EXPL_lm_allvars.0,
                               Predictor_var %in% c("INF_BIN_f",
                                                    "TREAT",
                                                    "MULT_INF"))
    EXPL_lm_allvars <- EXPL_lm_allvars[,c("Pr..F.",
                         "EXPL_var", "Predictor_var")]
    
    # (2) repeated measurement analyses
    
    # Now, we prepare data from repeated measurement analyses of "latency to reach
    # new positions" to join with the above database
    
    Lat_pos_G_EXPL_F_db$EXPL_var <- "Latency_to_reach_new_positions_log"
    Lat_pos_G_EXPL_F_db$Predictor_var <- row.names(Lat_pos_G_EXPL_F_db)  
    
    Lat_pos_MLG_EXPL_F_db$EXPL_var <- "Latency_to_reach_new_positions_log"
    Lat_pos_MLG_EXPL_F_db$Predictor_var <- row.names(Lat_pos_MLG_EXPL_F_db) 
    
    # Merge both results from "latency to reach new positions"
    
    Lat_pos_EXPL_F_db_all.0 <- rbind(Lat_pos_G_EXPL_F_db,
                                   Lat_pos_MLG_EXPL_F_db)%>%
      filter(Predictor_var %in% c("TREAT:Position_f",
                                  "INF_BIN_f:Position_f",
                                  "MULT_INF:Position_f"))
    
    # We are only interested in the repeated measurement test
    
    Lat_pos_EXPL_F_db_all <- Lat_pos_EXPL_F_db_all.0 [,c("Pr..F.",
                                                      "EXPL_var", 
                                                      "Predictor_var")]
    
    # (3) Merge and FDR correction
    
    FDR_EXPL_tests_p <- rbind(Lat_pos_EXPL_F_db_all,
                          EXPL_lm_allvars)
    
    FDR_EXPL_tests_p <- arrange(FDR_EXPL_tests_p, Pr..F.)
    
    FDR_EXPL_tests_p$FDR_p_adjust <-p.adjust(FDR_EXPL_tests_p$Pr..F.,
                                                  method="fdr")

    