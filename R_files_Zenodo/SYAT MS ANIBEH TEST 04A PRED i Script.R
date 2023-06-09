#.........................................................
## Analyses for manuscript: 
## Haemosporidian infections influence risk-taking behaviours in young male blackcaps Sylvia atricapilla
#.........................................................

    ## SUMMARY OF THE SCRIPT:====

# **Cautious note**: YOU MUST FIRST RUN DATA SCRIPT
# TO GENERATE ADEQUATED TRANSFORMED DATA FOR ANALYSES

# ...For PREDATOR EXPOSURE analyses use: 

# (A)MS_SYAT_BH_Pi = For immediate behaviour to predator presence

# In each database and for each variable two type of test:
    # ... EFFICACY OF THE TEST (APPENDIX A3):
    #for "order of the STIMULUS_f playback" and "type of STIMULUS_f"
    # ... PARASITISM EFFECT when the behavioural variable passed ok
    # the "efficacy of the test".

#...................................................................

library(dplyr)
library(vegan)
library(ggplot2)
library(reshape2)
library(psych)
library(nFactors)
library(lme4)
library(lmerTest)
library(rstatix)
library(lsmeans)
library(afex)
library(blmeco)
library(LMERConvenienceFunctions)
library(elrm)

library(PerformanceAnalytics )
library(ggpubr)
library(effects)
library(gridExtra)
library(Rmisc)
library(officer)
library(rvg)
library(tidyverse)

library(broom)
library(jtools)
    
    set.seed(123)
#...................................................................
    ## PRELIMINARY PREPARATION OF DATABASES FOR ANALYSES ----
#...................................................................

    #### PREDATOR EXPOSURE IMMEDIATE BEHAVIOUR (Pi)
    
    ### Merge Behaviour (IMMEDIATE) and Biological database

    PRED_i_BIOL_DB <- merge (MS_SYAT_BH_Pi, 
                                MS_SYAT_BH_BIOL,
                                by.x ="ID", 
                                by.y="ID",
                                all.x = TRUE,
                                no.dups=FALSE)
    str(PRED_i_BIOL_DB)

    ### Variable structure for analyses

    PRED_i_BIOL_DB <- mutate(PRED_i_BIOL_DB,
                               INF_BIN_f = factor(INF_BIN),
                               MULT_INF_f = factor(MULT_INF),
                               STIM_ORDER_f = factor(STIM_ORDER),
                               STIMULUS_f =factor(STIMULUS),
                               KEY_BEHAV = factor(KEY_BEHAV),
                             
                             #standardized variables for the plots
                             
                             Activity_stdz = as.numeric(paste(scale(Activity,
                                                   center=TRUE, scale=TRUE))),
                             RBT_stdz = as.numeric(paste(scale (RBT,
                                               center=TRUE, scale=TRUE))),
                             T_INITIAL_log_stdz = as.numeric(paste(scale (T_INITIAL_log,
                                              center=TRUE, scale=TRUE)))
                             )

#...................................................................
    ##...PREDATOR EXPOSURE analyses ----
#...................................................................
    # **Caution note**: When we had warnings in glmer models because
    # model failed to converged we include the following method (bobyqa)
    # to improve convergence controlling errors
    
    mi.control <- glmerControl(optimizer="bobyqa", 
                               optCtrl=list(maxfun=100000))
#.................................................................
    #### IMMEDIATE BEHAVIOUR TO PREDATOR EXPOSURE ####
#.................................................................
    
    #............................................................... 
    ###.... (1) Coherent crown raising response ####
    
    ###........ EFFICACY OF THE TEST ####
    
    xtabs(~KEY_BEHAV+STIMULUS_f
        ,data=PRED_i_BIOL_DB)

    ### We classified individuals in two categories:
    
        ###   (1) when individuals raised head crown only when the
        ### predator was present 
        ###   (0) Other alternatives

        # Transform from long to wide format:
      
    PRED_i_BIOL_DB_w <-reshape(PRED_i_BIOL_DB, 
                                   direction='wide',
                                   timevar="STIMULUS",
                                   v.names= "KEY_BEHAV",
                                   idvar="ID") 

    # We created new dichotomous variable:
    
    PRED_i_BIOL_DB_w <- mutate(PRED_i_BIOL_DB_w,
                                    CROWN_ID = ifelse(KEY_BEHAV.PRED==1 &
                                                          KEY_BEHAV.CTRL==0,
                                                      1,0))
        
    # Test only with the order of the stimulus
    
    xtabs(~CROWN_ID+STIM_ID_ORDER_f
          ,data=PRED_i_BIOL_DB_w)
    
    CR_ORS_PRED_glm = glm(CROWN_ID ~ STIM_ID_ORDER_f, 
                                 family=binomial,
                                 data=PRED_i_BIOL_DB_w)
    
    drop1(CR_ORS_PRED_glm,
          test="Chisq")
    
    # The order of the STIMULUS_f did not influence likely of head crown rising
    # according to risk
    
    xtabs(~KEY_BEHAV + STIMULUS_f,
          data=PRED_i_BIOL_DB)
    
    # ... Effect of the type of stimulus alone
    
    CR_STIM_PRED_glmer = glmer(KEY_BEHAV ~ STIMULUS_f +
                                   (1 | ID), 
                               family=binomial,
                               data=PRED_i_BIOL_DB)
    
    drop1(CR_STIM_PRED_glmer,
          test="Chisq")
    
     plot(allEffects (CR_STIM_PRED_glmer))
     
    # Birds erected head crown more likely with predator than with control 
    # (LRT=79.23, P < 0.001, 30 out of 42 individuals erected head crown only with predator)
     
     # Correlation with body condition (BC)
     
     CR_BC_PRED_glm = glm(CROWN_ID ~ BC, 
                          family=binomial,
                          data=PRED_i_BIOL_DB_w)
     
     drop1(CR_BC_PRED_glm,
           test="Chisq") # not significant
     # overdispersion
     phi <- sum((residuals(CR_BC_PRED_glm, 
                           type="pearson"))^2)/CR_BC_PRED_glm$df.residual
     print(c("Pearson overdispersion =", 
             round(phi, 3)), quote=FALSE)
    
    ###........ PARASITISM EFFECT: ####
    
    #........ STATUS OF INFECTION ====
     
     CRX1 <- xtabs(~CROWN_ID+interaction(INF_BIN_f,TREAT)
           ,data=PRED_i_BIOL_DB_w) 

     CRX1dat <- data.frame(INF_BIN_f = rep(0:1, 2), TREAT = rep(c("PQ","Water"),
                                                                     each=2), 
                                CROWN_ID = c(CRX1[2, ]),
                                n = colSums(CRX1))
     CRX1dat  # view collapsed data set
     
     # Exact logistic regression because some cells were empty
     set.seed(123)
     CRX1exact = elrm(CROWN_ID/n ~ INF_BIN_f +
                      TREAT +
                      INF_BIN_f:TREAT, 
                    interest = ~INF_BIN_f +
                      TREAT +
                      INF_BIN_f:TREAT,
                    iter= 1e+05, 
                     burnIn=2000, data=CRX1dat,
                    r=2)
     
     summary(CRX1exact) # interaction not significant
     # plot(CRX1exact)

     # Testing simple effects without interaction
     set.seed(123)
     CRX1exact2 = elrm(CROWN_ID/n ~ INF_BIN_f +
                      TREAT, 
                    interest = ~INF_BIN_f +
                      TREAT,
                    iter= 1e+05, 
                    burnIn=2000, data=CRX1dat,
                    r=2)
     
     summary(CRX1exact2) 
     
     # Extract P values for FDR correction
     
     CROWN_INF_BIN_p <- data.frame(p =CRX1exact2$p.values)
     CROWN_INF_BIN_p$Predictor_var <- rownames(CROWN_INF_BIN_p)
     CROWN_INF_BIN_p$PRED_var <- "CROWN_rise"
     
     # plot(CRX1exact2)
    
    #........ MULTIPLE STATUS OF INFECTION ====

     CRMLGX1 <- xtabs(~CROWN_ID+interaction(MULT_INF_f,TREAT)
                 ,data=PRED_i_BIOL_DB_w)
     
     CRMLGX1dat <- data.frame(MULT_INF_f = rep(0:2, 2), TREAT = rep(c("PQ","Water"),
                                                              each = 3), 
                         CROWN_ID = CRMLGX1[2, ], n = colSums(CRMLGX1))
     CRMLGX1dat  # view collapsed data set
     
     # Exact logistic regression because some cells were empty
     set.seed(123)
     CRMLGX1exact1 = elrm(CROWN_ID/n ~ MULT_INF_f +
                      TREAT +
                       MULT_INF_f:TREAT, 
                    interest = ~MULT_INF_f +
                      TREAT +
                      MULT_INF_f:TREAT,
                    iter= 1e+05, 
                    burnIn=2000, data=CRMLGX1dat,
                    r=2)
     CRMLGX1exact2 <- update(CRMLGX1exact1, iter = 1e+05, 
                        burnIn = 5000)
     
     summary(CRMLGX1exact2) # interaction significant
     # plot(CRMLGX1exact2)
     
     # Extract P values for FDR correction
     
     CROWN_MLGINF_p <- data.frame(p =CRMLGX1exact2$p.values)
     CROWN_MLGINF_p$Predictor_var <- rownames(CROWN_MLGINF_p)
     CROWN_MLGINF_p$PRED_var <- "CROWN_rise"
     # only interaction for the FDR correction
     CROWN_MLGINF_int_p <- filter(CROWN_MLGINF_p,
                                 Predictor_var=="MULT_INF_f:TREATWater") 
    # POST-HOC
    
    # For SINGLE INFECTED 
     
     MLGX1_sINF <- xtabs(~CROWN_ID + TREAT
                  ,data=filter(PRED_i_BIOL_DB_w,
                               MULT_INF_f==1))
     
     MLGX1_sINFglm = glm(CROWN_ID ~ TREAT, 
                           family=binomial,
                       data=filter(PRED_i_BIOL_DB_w,
                                   MULT_INF_f==1))
     
     drop1(MLGX1_sINFglm,
           test="Chisq") # not significant
     
     # overdispersion test
     phi <- sum((residuals(MLGX1_sINFglm, type="pearson"))^2)/MLGX1_sINFglm$df.residual
     print(c("Pearson overdispersion =", round(phi, 3)), quote=FALSE)
    
    # For MULTIPLE INFECTED 
    
     MLGX1_mINF <- xtabs(~CROWN_ID + TREAT
                       ,data=filter(PRED_i_BIOL_DB_w,
                                    MULT_INF_f==2))
     
     MLGX1_mINFglm = glm(CROWN_ID ~ TREAT, 
                       family=binomial,
                       data=filter(PRED_i_BIOL_DB_w,
                                   MULT_INF_f==2))
     
     drop1(MLGX1_mINFglm,
           test="Chisq") # not significant
     summary(MLGX1_mINFglm)
     
     # overdispersion test
     phi <- sum((residuals(MLGX1_mINFglm, type="pearson"))^2)/MLGX1_mINFglm$df.residual
     print(c("Pearson overdispersion =", round(phi, 3)), quote=FALSE)

     # For NOT INFECTED
    
     MLGX1_noINF <- xtabs(~CROWN_ID + TREAT
                       ,data=filter(PRED_i_BIOL_DB_w,
                                    MULT_INF_f==0))
     
     MLGX1_noINFdat <- data.frame(TREAT = rep(c("PQ","Water")), 
                               CROWN_ID = MLGX1_noINF[2, ], 
                               n = colSums(MLGX1_noINF))
     MLGX1_noINFdat  # view collapsed data set
     
     # Exact logistic regression because some cells were empty
     set.seed(123)
     MLGX1_noINFexact = elrm(CROWN_ID/n ~ TREAT, 
                          interest = ~TREAT,
                          iter= 1e+05, 
                          burnIn=2000, data=MLGX1_noINFdat,
                          r=2)
     summary(MLGX1_noINFexact) # not significant
     # plot(MLGX1_noINFexact)
     
     # Multiple-test correction
     
     p.adjust( c(0.06, 0.19,0.55), method = "bonferroni")

     #............................................................... 
     ###.... (2) Activity level ####
     
     # ...Previous tests
     
     # Identify outliers
     
     PRED_i_BIOL_DB %>%
       group_by(STIMULUS_f) %>%
       identify_outliers(Activity_stdz)
     
     # Normality tests
     
     PRED_i_BIOL_DB %>%
       group_by(INF_BIN_f, TREAT, 
                STIMULUS_f) %>%
       shapiro_test(Activity_stdz)
     
     ggqqplot(PRED_i_BIOL_DB, "Activity_stdz",
              facet.by = c("STIMULUS_f"))
     
     ggplot (PRED_i_BIOL_DB,
             aes(x=Activity_stdz))+
       geom_histogram()+
       facet_wrap(vars(STIMULUS_f))
     
     # Check normality also in residuals of the models
     
     # Levene's tests
     
     PRED_i_BIOL_DB %>%
       group_by(STIMULUS_f) %>%
       levene_test(Activity_stdz ~ STIM_ID_ORDER_f,
                   center=mean)
     
     PRED_i_BIOL_DB %>%
       group_by(STIMULUS_f) %>%
       levene_test(Activity_stdz ~ TREAT*INF_BIN_f,
                   center=mean) 
     
     PRED_i_BIOL_DB %>%
       group_by(STIMULUS_f) %>%
       levene_test(Activity_stdz ~ TREAT,
                   center=mean) 
     
     PRED_i_BIOL_DB %>%
       group_by(STIMULUS_f) %>%
       levene_test(Activity_stdz ~ INF_BIN_f,
                   center=mean) 
     
     PRED_i_BIOL_DB %>%
       group_by(STIMULUS_f) %>%
       levene_test(Activity_stdz ~ MULT_INF_f,
                   center=mean) 
     
     # Homogeneity of covariances assumption
     # the Box's M is highly sensitive, so unless p < 0.001 and your sample sizes 
     # are unequal, we ignored it.
     
     box_m(PRED_i_BIOL_DB[,
                          "Activity_stdz",
                          drop = FALSE],PRED_i_BIOL_DB$TREAT) 
     
     PRED_i_BIOL_DB_MULT_G <- PRED_i_BIOL_DB %>%
       filter(!is.na(MULT_INF_f))# select not NA values
     
     box_m(PRED_i_BIOL_DB_MULT_G[,
                                 "Activity_stdz",
                                 drop = FALSE], 
           PRED_i_BIOL_DB_MULT_G$MULT_INF_f) 
     
     box_m(PRED_i_BIOL_DB[,
                          "Activity_stdz",
                          drop = FALSE], 
           PRED_i_BIOL_DB$STIM_ID_ORDER_f)
     
     ###........ EFFICACY OF THE TEST #### 
     #........ Previous visualization====
     
     #........Differences between groups
     
     ggplot(data= filter(PRED_i_BIOL_DB),
            aes(y=Activity_stdz, 
                x= STIMULUS_f)) +
       geom_boxplot(aes(fill=as.factor(STIM_ID_ORDER_f))) + 
       geom_point() +
       theme_classic()+
       facet_wrap(vars(STIM_ID_ORDER_f))+
       ggtitle("PREDATOR ORDER") +
       theme(axis.text.x = element_text(angle = 90)) 
     
     #........Individual response
     
     ggplot(data=PRED_i_BIOL_DB, 
            aes(x=ID,
                y=Activity_stdz,
                group=ID)) +
       theme_classic()+
       geom_point(aes(color=STIMULUS_f, size=1.5,alpha=0.5))+
       geom_line()+
       theme(legend.position = "none",
             axis.text.x = element_text(angle = 90))
     
     #........ Repeated measures test ====
     
     Act_OR_PRED_0 =aov_car(Activity_stdz ~ STIMULUS_f*STIM_ID_ORDER_f +
                              Error(ID/STIMULUS_f),
                            data=PRED_i_BIOL_DB,
                            type=3)
     
     summary(Act_OR_PRED_0) # Not significant
     
     # Test effect
     
     PRED_i_BIOL_DB %>%
       pairwise_t_test(Activity_stdz ~ STIMULUS_f,
                       paired = TRUE)
     
     # Bird activity was higher with the sparrowhawk than with the bottle
     
     # Correlation with body condition (BC)
     PRED_i_BIOL_DB <-mutate(PRED_i_BIOL_DB,
                             BC_stdz = scale(BC))
     
     Act_BC_PRED =aov_car(Activity_stdz ~ STIMULUS_f*BC_stdz +
                            Error(ID/STIMULUS_f),
                          data=PRED_i_BIOL_DB,
                          factorize = FALSE,
                          type=3)
     
     summary(Act_BC_PRED)
     
     ###........ PARASITISM EFFECT: ####
     #........ STATUS OF INFECTION  ====
     
     Activ_G_PRED.0 = aov_car(Activity_stdz~ INF_BIN_f*TREAT*
                                STIMULUS_f +  
                                Error(ID/STIMULUS_f), 
                              data=PRED_i_BIOL_DB,
                              type=3)
     summary(Activ_G_PRED.0)
     Activ_G_PRED.0 
     
     # Interaction not significant
     
     Activ_G_PRED = aov_car(Activity_stdz~ INF_BIN_f +
                              TREAT+
                              STIMULUS_f +
                              INF_BIN_f:STIMULUS_f +
                              TREAT:STIMULUS_f +  
                              Error(ID/STIMULUS_f), 
                            data=PRED_i_BIOL_DB,
                            type=3)
     summary(Activ_G_PRED)
     Activ_G_PRED # Not significant
     
     Activ_G_PRED_F <- Activ_G_PRED$anova_table
     
     Activ_G_PRED_F_db <- data.frame(Activ_G_PRED_F)  # database with F and P for FDR correction
     
     # Model residuals
     
     plots <- list()  
     for (i in c(1:2)){
       
       p1 = ggqqplot(as.numeric(residuals(Activ_G_PRED$lm)[,i]), 
                     title=paste(colnames(residuals(Activ_G_PRED$lm))[i]))
       plots[[i]] <- p1 
     }
     
     multiplot(plotlist = plots, 
               cols = 2)
     
     # POST-HOC-Estimates (for plotting)
     
     (ls.Activ.G.1 <- lsmeans(Activ_G_PRED,
                            c("STIMULUS_f"),
                            by="INF_BIN_f"))
     
     # between levels of infection
     Estim_PostHoc_Activ_G <-pairs(update(pairs(ls.Activ.G.1), by=NULL))
     
     confint(Estim_PostHoc_Activ_G)
     
     #........ MULTIPLE STATUS OF INFECTION ====
     
     Activ_MLG_PRED.0 = aov_car(Activity_stdz~ MULT_INF_f*TREAT*
                                  STIMULUS_f +  
                                  Error(ID/STIMULUS_f), 
                                data=PRED_i_BIOL_DB,
                                type=3)
     summary(Activ_MLG_PRED.0)
     Activ_MLG_PRED.0 
     
     # Interaction not significant
     
     Activ_MLG_PRED = aov_car(Activity_stdz~ MULT_INF_f +
                                TREAT+
                                STIMULUS_f +
                                MULT_INF_f:STIMULUS_f +
                                TREAT:STIMULUS_f +  
                                Error(ID/STIMULUS_f), 
                              data=PRED_i_BIOL_DB,
                              type=3)
     summary(Activ_MLG_PRED)
     Activ_MLG_PRED 
     
     Activ_MLG_PRED_F <- Activ_MLG_PRED$anova_table
     
     Activ_MLG_PRED_F_db <- data.frame(Activ_MLG_PRED_F)  # database with F and P for FDR correction
     
     # Model residuals
     
     plots <- list()  
     for (i in c(1:2)){
       
       p1 = ggqqplot(as.numeric(residuals(Activ_MLG_PRED$lm)[,i]), 
                     title=paste(colnames(residuals(Activ_MLG_PRED$lm))[i]))
       plots[[i]] <- p1 
     }
     
     multiplot(plotlist = plots, 
               cols = 2)
     
     # POST-HOC TREAT
     
     (ls.MLG.TREAT <- lsmeans(Activ_MLG_PRED,
                              c("STIMULUS_f"),
                              by="TREAT"))
     
     update(pairs(ls.MLG.TREAT), by=NULL, adjust = "Bonferroni")#
     
     # POST-HOC-Estimates (for plotting) WITH INFECTION
     
     (ls.Activ.MLG.1 <- lsmeans(Activ_MLG_PRED,
                              c("STIMULUS_f"),
                              by="MULT_INF_f"))
     
     # between levels of infection
     Estim_PostHoc_Activ_MLG <-pairs(update(pairs(ls.Activ.MLG.1), by=NULL))
     
     confint(Estim_PostHoc_Activ_MLG, adjust="none")
     confint(Estim_PostHoc_Activ_MLG, adjust="tukey")
     
     #............................................................... 
     ###.... (3) Relative body turns (RBT) ####
     
     # ...Previous tests
     
     # Identify outliers
     
     PRED_i_BIOL_DB %>%
       group_by(STIMULUS_f) %>%
       identify_outliers(RBT_stdz)
     
     # Normality tests
     
     PRED_i_BIOL_DB %>%
       group_by(INF_BIN_f, TREAT, 
                STIMULUS_f) %>%
       shapiro_test(RBT_stdz)
     
     ggqqplot(PRED_i_BIOL_DB, "RBT_stdz",
              facet.by = c("STIMULUS_f"))
     
     ggplot (PRED_i_BIOL_DB,
             aes(x=RBT_stdz))+
       geom_histogram()+
       facet_wrap(vars(STIMULUS_f))
     
     # Check normality also in residuals of the models
     
     # Levene's tests
     
     PRED_i_BIOL_DB %>%
       group_by(STIMULUS_f) %>%
       levene_test(RBT_stdz ~ STIM_ID_ORDER_f,
                   center=mean) 
     
     PRED_i_BIOL_DB %>%
       group_by(STIMULUS_f) %>%
       levene_test(RBT_stdz ~ TREAT*INF_BIN_f,
                   center=mean) 
  
     PRED_i_BIOL_DB %>%
       group_by(STIMULUS_f) %>%
       levene_test(RBT_stdz ~ TREAT,
                   center=mean) 
     
     PRED_i_BIOL_DB %>%
       group_by(STIMULUS_f) %>%
       levene_test(RBT_stdz ~ INF_BIN_f,
                   center=mean) 
     
     PRED_i_BIOL_DB %>%
       group_by(STIMULUS_f) %>%
       levene_test(RBT_stdz ~ MULT_INF_f,
                   center=mean) 
     
     # Homogeneity of covariances assumption
     # the Box's M is highly sensitive, so unless p < 0.001 and your sample sizes 
     # are unequal, we ignored it.
     
     box_m(PRED_i_BIOL_DB[,
                          "RBT_stdz",
                          drop = FALSE],PRED_i_BIOL_DB$TREAT) 
     
     PRED_i_BIOL_DB_MULT_G <- PRED_i_BIOL_DB %>%
       filter(!is.na(MULT_INF_f))# select not NA values
     
     box_m(PRED_i_BIOL_DB_MULT_G[,
                                 "RBT_stdz",
                                 drop = FALSE], 
           PRED_i_BIOL_DB_MULT_G$MULT_INF_f) 
     
     box_m(PRED_i_BIOL_DB[,
                          "RBT_stdz",
                          drop = FALSE], PRED_i_BIOL_DB$STIM_ID_ORDER_f)#ok
     ###........ EFFICACY OF THE TEST #### 
     #........ Previous visualization====
     
     #........Differences between groups
     
     ggplot(data= filter(PRED_i_BIOL_DB),
            aes(y=RBT_stdz, 
                x= STIMULUS_f)) +
       geom_boxplot(aes(fill=as.factor(STIM_ID_ORDER_f))) + 
       geom_point() +
       theme_classic()+
       facet_wrap(vars(STIM_ID_ORDER_f))+
       ggtitle("PREDATOR ORDER") +
       theme(axis.text.x = element_text(angle = 90)) 
     
     #........Individual response
     
     ggplot(data=PRED_i_BIOL_DB, 
            aes(x=ID,
                y=RBT_stdz,
                group=ID)) +
       theme_classic()+
       geom_point(aes(color=STIMULUS_f, size=1.5,alpha=0.5))+
       geom_line()+
       theme(legend.position = "none",
             axis.text.x = element_text(angle = 90))
     
     #........ Repeated measures analyses====
     
     RBT_OR_PRED_0 =aov_car(RBT_stdz ~ STIMULUS_f*STIM_ID_ORDER_f +
                              Error(ID/STIMULUS_f),
                            data=PRED_i_BIOL_DB,
                            type=3)
     
     summary(RBT_OR_PRED_0) 
     
     # Interaction not significant. Order of the stimulus did not affect
     # RBT_stdz response
     
     RBT_OR_PRED =aov_car(RBT_stdz ~ STIMULUS_f +
                            Error(ID/STIMULUS_f),
                          data=PRED_i_BIOL_DB,
                          type=3)
     
     summary(RBT_OR_PRED) 
     
     PRED_i_BIOL_DB %>%
       pairwise_t_test(RBT_stdz ~ STIMULUS_f,
                       paired = TRUE)
     
     # RBT was higher with the sparrowhawk than with the bottle
     
     # Correlation with body condition (BC)
     
     RBT_BC_PRED =aov_car(RBT_stdz ~ STIMULUS_f*BC_stdz +
                            Error(ID/STIMULUS_f),
                          data=PRED_i_BIOL_DB,
                          factorize = FALSE,
                          type=3)
     
     summary(RBT_BC_PRED)
     
     ###........ PARASITISM EFFECT: ####
     #........ STATUS OF INFECTION  ====
     
     RBT_G_PRED.0 = aov_car(RBT_stdz~ INF_BIN_f*TREAT*
                              STIMULUS_f +  
                              Error(ID/STIMULUS_f), 
                            data=PRED_i_BIOL_DB,
                            type=3)
     summary(RBT_G_PRED.0)
     RBT_G_PRED.0 
     
     # Interaction not significant
     
     RBT_G_PRED = aov_car(RBT_stdz~ INF_BIN_f +
                            TREAT+
                            STIMULUS_f +
                            INF_BIN_f:STIMULUS_f +
                            TREAT:STIMULUS_f +  
                            Error(ID/STIMULUS_f), 
                          data=PRED_i_BIOL_DB,
                          type=3)
     summary(RBT_G_PRED)
     RBT_G_PRED 
     
     RBT_G_PRED_F <- RBT_G_PRED$anova_table
     
     RBT_G_PRED_F_db <- data.frame(RBT_G_PRED_F)  # database with F and P for FDR correction
     
     # Model residuals
     
     plots <- list()  
     for (i in c(1:2)){
       
       p1 = ggqqplot(as.numeric(residuals(RBT_G_PRED$lm)[,i]), 
                     title=paste(colnames(residuals(RBT_G_PRED$lm))[i]))
       plots[[i]] <- p1 
     }
     
     multiplot(plotlist = plots, 
               cols = 2)
     
     # POST-HOC-Estimates
     
     (ls.RBT.G.1 <- lsmeans(RBT_G_PRED,
                            c("STIMULUS_f"),
                            by="INF_BIN_f"))
     
     # within individuals
     update(pairs(ls.RBT.G.1), by=NULL, adjust = "none")#  
     update(pairs(ls.RBT.G.1), by=NULL, adjust = "Bonferroni")#  
     
     # between levels of infection
     Estim_PostHoc_RBT_G <-pairs(update(pairs(ls.RBT.G.1), by=NULL))
    
     confint(Estim_PostHoc_RBT_G)

     # Plot with the effect of STATUS OF INFECTION 
     
     ggplot(data= PRED_i_BIOL_DB,
            aes(y=RBT_stdz, 
                x= STIMULUS_f,
                fill=INF_BIN_f))  +
       geom_boxplot() + 
       geom_point() +
       theme_classic()+
       facet_wrap(vars(INF_BIN_f))
     #........ MULTIPLE STATUS OF INFECTION ====
     
     RBT_MLG_PRED.0 = aov_car(RBT_stdz~ MULT_INF_f*TREAT*
                                STIMULUS_f +  
                                Error(ID/STIMULUS_f), 
                              data=PRED_i_BIOL_DB,
                              type=3)
     summary(RBT_MLG_PRED.0)
     RBT_MLG_PRED.0 
     
     # Interaction not significant
     
     RBT_MLG_PRED = aov_car(RBT_stdz~ MULT_INF_f +
                              TREAT+
                              STIMULUS_f +
                              MULT_INF_f:STIMULUS_f +
                              TREAT:STIMULUS_f +  
                              Error(ID/STIMULUS_f), 
                            data=PRED_i_BIOL_DB,
                            type=3)
     summary(RBT_MLG_PRED)
     RBT_MLG_PRED 
     
     RBT_MLG_PRED_F <- RBT_MLG_PRED$anova_table
     
     RBT_MLG_PRED_db <- data.frame(RBT_MLG_PRED_F)  # database with F and P for FDR correction
     
     # Model residuals
     
     plots <- list()  
     for (i in c(1:2)){
       
       p1 = ggqqplot(as.numeric(residuals(RBT_MLG_PRED$lm)[,i]), 
                     title=paste(colnames(residuals(RBT_MLG_PRED$lm))[i]))
       plots[[i]] <- p1 
     }
     
     multiplot(plotlist = plots, 
               cols = 2)
     
     # POST-HOC-Estimates
     
     (ls.RBT.MLG.1 <- lsmeans(RBT_MLG_PRED,
                              c("STIMULUS_f"),
                              by="MULT_INF_f"))
     
     # within individuals
     update(pairs(ls.RBT.MLG.1), by=NULL, adjust = "none")
     update(pairs(ls.RBT.MLG.1), by=NULL, adjust = "Bonferroni")
     
     # between levels of infection
     Estim_PostHoc_RBT_ML <-pairs(update(pairs(ls.RBT.MLG.1), by=NULL))
     
     confint(Estim_PostHoc_RBT_ML, adjust="none")
     confint(Estim_PostHoc_RBT_ML, adjust="tukey")

     # Plot with the effect of STATUS OF INFECTION 
     
     # function for computing mean, 95%CI, max and min values
     
     min.mean.sd.max <- function(x) {
       r <- c(mean(x,na.rm = T) - qt(.975, 
                                     length(!is.na(x) - 1))*(sd(x,na.rm = T)/sqrt(length(!is.na(x)))), 
              mean(x,na.rm = T) - sd(x,na.rm = T)/sqrt(length(!is.na(x))),
              mean(x,na.rm = T), 
              mean(x,na.rm = T) + sd(x,na.rm = T)/sqrt(length(!is.na(x))),
              mean(x,na.rm = T) + qt(.975, 
                                     length(!is.na(x) - 1))*(sd(x,na.rm = T)/sqrt(length(!is.na(x))))
       )
       names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
       r
     }
     
     PRED_i_BIOL_MULT_DB <- mutate(PRED_i_BIOL_DB,
                                   STIMULUS_f_long = ifelse(STIMULUS_f ==
                                                              "CTRL","Control",
                                                            ifelse(STIMULUS_f=="PRED",
                                                                   "Predator",NA))) 
     
     Fig4 <- ggplot(data= filter(PRED_i_BIOL_MULT_DB,
                                 !is.na(MULT_INF_f)),
                    aes(y=RBT, 
                        x= STIMULUS_f_long,
                        fill=MULT_INF_f,
                        alpha=0.5))
     set.seed(123)
     
     Fig4 <- Fig4 + stat_summary(fun.data = min.mean.sd.max, 
                                 geom = "boxplot")+
       scale_fill_manual(values=c("#332288","#DDCC77","#661100"))+
       geom_jitter(position=position_jitter(width=.2), size=2,
                   shape=16,
                   aes(alpha=0.5)) +
       theme_classic()+
       facet_wrap(vars(MULT_INF_f))+
       ylab("Relative body turns")+ 
       xlab("") +
       theme(axis.text.x = element_text(angle = 45,
                                        vjust = 0.5,
                                        hjust=0.5, 
                                        size=14,
                                        colour="black"),
             axis.text.y = element_text (size = 14,
                                         hjust = 0.5,
                                         colour="black"),
             axis.title.y = element_text(size = 14,
                                         vjust=3,
                                         colour="black"),
             axis.ticks = element_line(colour="black"))+
       scale_y_continuous(breaks=c(-1.5,-1,-0.5,0,0.5,1,1.5))
     
     Fig4.nolegend <- Fig4 +
       theme(legend.position = "none",
             strip.background = element_blank(),
             strip.text.x = element_blank())
     
     Fig4.nolegend
     
     # Code to export figure 4 in pptx format:
     
     Fig4_dml <- rvg::dml(ggobj = Fig4.nolegend)
     
     # initialize PowerPoint slide 
     officer::read_pptx() %>%
       # add slide 
       officer::add_slide() %>%
       # specify object and location of object 
       officer::ph_with(Fig4_dml,
                        ph_location(height = 5, width = 5.5)) %>%
       # export slide 
       base::print(
         target = here::here(
           "Fig4.pptx"
         )
       )
     
     #............................................................... 
    ###.... (4) Initial latency ####

    # ...Previous tests
    
    # Identify outliers
    
    PRED_i_BIOL_DB %>%
        group_by(STIMULUS_f) %>%
        identify_outliers(T_INITIAL_log_stdz)
    
    # Normality tests
    
    PRED_i_BIOL_DB %>%
        group_by(INF_BIN_f, TREAT, 
                 STIMULUS_f) %>%
        shapiro_test(T_INITIAL_log_stdz)
    
    ggqqplot(PRED_i_BIOL_DB, "T_INITIAL_log",
             facet.by = c("STIMULUS_f"))
    
    ggplot (PRED_i_BIOL_DB,
            aes(x=T_INITIAL_log_stdz))+
        geom_histogram()+
        facet_wrap(vars(STIMULUS_f))
    
    # Check normality also in residuals of the models
    
    # Levene's tests
 
    PRED_i_BIOL_DB %>%
        group_by(STIMULUS_f) %>%
        levene_test(T_INITIAL_log_stdz ~ STIM_ID_ORDER_f,
                    center=mean) 
    
    PRED_i_BIOL_DB %>%
        group_by(STIMULUS_f) %>%
        levene_test(T_INITIAL_log_stdz ~ TREAT*INF_BIN_f,
                    center=mean) 
    
    PRED_i_BIOL_DB %>%
        group_by(STIMULUS_f) %>%
        levene_test(T_INITIAL_log_stdz ~ TREAT,
                    center=mean) 
    
    PRED_i_BIOL_DB %>%
        group_by(STIMULUS_f) %>%
        levene_test(T_INITIAL_log_stdz ~ INF_BIN_f,
                    center=mean) 
    
    PRED_i_BIOL_DB %>%
        group_by(STIMULUS_f) %>%
        levene_test(T_INITIAL_log_stdz ~ MULT_INF_f,
                    center=mean)
    
    # Homogeneity of covariances assumption
    # the Box's M is highly sensitive, so unless p < 0.001 and your sample sizes 
    # are unequal, we ignored it.
    
    box_m(PRED_i_BIOL_DB[,
                            "T_INITIAL_log_stdz",
                            drop = FALSE],PRED_i_BIOL_DB$TREAT) 
    
    PRED_i_BIOL_DB_MULT_G <- PRED_i_BIOL_DB %>%
        filter(!is.na(MULT_INF_f))# select not NA values
    
    box_m(PRED_i_BIOL_DB_MULT_G[,
                                   "T_INITIAL_log_stdz",
                                   drop = FALSE], 
          PRED_i_BIOL_DB_MULT_G$MULT_INF_f) 
    
    box_m(PRED_i_BIOL_DB[,
                            "T_INITIAL_log_stdz",
                            drop = FALSE], 
          PRED_i_BIOL_DB$STIM_ID_ORDER_f)
    
    ###........ EFFICACY OF THE TEST #### 
    
    #........ Previous visualization====
    
    #........Differences between groups
    
    ggplot(data= filter(PRED_i_BIOL_DB),
           aes(y=T_INITIAL_log_stdz, 
               x= STIMULUS_f))  +
        geom_boxplot(aes(fill=as.factor(STIM_ID_ORDER_f))) + 
        geom_point() +
        theme_classic()+
        facet_wrap(vars(STIM_ID_ORDER_f))+
        ggtitle("PREDATOR ORDER") +
        theme(axis.text.x = element_text(angle = 90)) 
    
    #........Individual response
    
    ggplot(data=PRED_i_BIOL_DB, 
           aes(x=ID,
               y=T_INITIAL_log_stdz,
               group=ID)) +
        theme_classic()+
        geom_point(aes(color=STIMULUS_f, size=1.5,alpha=0.5))+
        geom_line()+
        theme(legend.position = "none",
              axis.text.x = element_text(angle = 90))
    
    #........ Repeated measures test ====
    
    Lat_OR_PRED_0 =aov_car(T_INITIAL_log_stdz ~ STIMULUS_f*STIM_ID_ORDER_f +
                               Error(ID/STIMULUS_f),
                               data=PRED_i_BIOL_DB,
                               type=3)
    
    summary(Lat_OR_PRED_0)
    
    # Post-hoc analyses to explore further effect of the order of the conspecific
    # stimulus
   
    Lat.two.way <- PRED_i_BIOL_DB %>%
        dplyr::select(STIM_ID_ORDER_f,T_INITIAL_log_stdz,
                      STIMULUS_f, ID) %>%
        group_by(STIM_ID_ORDER_f) %>%
        anova_test(dv = T_INITIAL_log_stdz , 
                   wid = ID, 
                   within = c(STIMULUS_f),
                   type=3)
    Lat.two.way
    
    p.adjust(c(0.002, 0.541), method="bonferroni")
    
    # ORDER 1 = PRED first # not significant
    # ORDER 2 = PRED second # Significant
    
# **When the sparrow hack stimulus was presented in the second place **
# **latency differs between type of stimulus**.
    
    # Sampling subset that includes the individuals in which sparrowhawk
    # was shown in second place.
    
    PRED_Lat_sel_02=filter(PRED_i_BIOL_DB,
                           STIM_ID_ORDER_f==2) %>%
        mutate(ID=factor(ID))
    
    length(unique(PRED_Lat_sel_02$ID)) # 21 individuals
    
    PRED_Lat_sel_02 %>%
        pairwise_t_test(T_INITIAL_log_stdz ~ STIMULUS_f,
                        paired = TRUE)
    
    # Correlation with body condition (BC)
    PRED_Lat_sel_02 <-mutate(PRED_Lat_sel_02,
                             BC_stdz = scale(BC))
    
    Lat_BC_PRED =aov_car(T_INITIAL_log_stdz ~ STIMULUS_f*BC_stdz +
                           Error(ID/STIMULUS_f),
                         data=PRED_Lat_sel_02,
                         factorize = FALSE,
                         type=3)
    
    summary(Lat_BC_PRED)# not significant
    
    ###........ PARASITISM EFFECT ####
    
    # Not tested interactions due to low sample size.
    
    #........ STATUS OF INFECTION ====
    
    Lat_G_PRED = aov_car(T_INITIAL_log_stdz~ INF_BIN_f +
                              TREAT+
                              STIMULUS_f +
                              INF_BIN_f:STIMULUS_f +
                              TREAT:STIMULUS_f +
                              Error(ID/STIMULUS_f), 
                            data=PRED_Lat_sel_02,
                            type=3)
    summary(Lat_G_PRED)
    Lat_G_PRED # no effect
    
    Lat_G_PRED_F <- Lat_G_PRED$anova_table
    
    Lat_G_PRED_F_db <- data.frame(Lat_G_PRED_F)  # database with F and P for FDR correction
    
    
    # Model residuals

    plots <- list()  
    for (i in c(1:2)){
        
        p1 = ggqqplot(as.numeric(residuals(Lat_G_PRED$lm)[,i]), 
                      title=paste(colnames(residuals(Lat_G_PRED$lm))[i]))
        plots[[i]] <- p1 
    }
    
    multiplot(plotlist = plots, 
              cols = 2)
    
    # POST-HOC-Estimates (for plotting)
    
    (ls.Lat_G.1 <- lsmeans(Lat_G_PRED,
                           c("STIMULUS_f"),
                           by="INF_BIN_f"))
    
    # between levels of infection
    Estim_PostHoc_Lat_G <-pairs(update(pairs(ls.Lat_G.1), by=NULL))
    
    confint(Estim_PostHoc_Lat_G)

    #........ MULTIPLE STATUS OF INFECTION ====
    
    PRED_Lat_sel_02_MLG <- filter(PRED_Lat_sel_02,
                                  !is.na(MULT_INF_f))
    
    Lat_MLG_PRED = aov_car(T_INITIAL_log_stdz~ MULT_INF_f +
                             TREAT+
                             STIMULUS_f +
                             MULT_INF_f:STIMULUS_f +
                             TREAT:STIMULUS_f +
                             Error(ID/STIMULUS_f), 
                         data=PRED_Lat_sel_02_MLG,
                         type=3)
    summary(Lat_MLG_PRED)
    Lat_MLG_PRED # no effect
    
    Lat_MLG_PRED_F <- Lat_MLG_PRED$anova_table
    
    Lat_MLG_PRED_F_db <- data.frame(Lat_MLG_PRED_F)  # database with F and P for FDR correction
    
    # Model residuals

    plots <- list()  
    for (i in c(1:2)){
        
        p1 = ggqqplot(as.numeric(residuals(Lat_MLG_PRED$lm)[,i]), 
                      title=paste(colnames(residuals(Lat_MLG_PRED$lm))[i]))
        plots[[i]] <- p1 
    }
    
    multiplot(plotlist = plots, 
              cols = 2)
    
    # POST-HOC-Estimates (for plotting)
    
    (ls.Lat.MLG.1 <- lsmeans(Lat_MLG_PRED,
                             c("STIMULUS_f"),
                             by="MULT_INF_f"))
    
    # between levels of infection
    Estim_PostHoc_Lat_ML <-pairs(update(pairs(ls.Lat.MLG.1), by=NULL))
    
    confint(Estim_PostHoc_Lat_ML, adjust="none")
    confint(Estim_PostHoc_Lat_ML, adjust="tukey")
    
    #............................................................... 
    ###...(5) FDR CORRECTION BY MULTIPLE TESTING DURING THE EXPLORATORY TEST####

    Lat_G_PRED_F_db$PRED_var <- "Latency_log"
    Lat_G_PRED_F_db$Predictor_var <- row.names(Lat_G_PRED_F_db)
    Lat_MLG_PRED_F_db$PRED_var <- "Latency_log"
    Lat_MLG_PRED_F_db$Predictor_var <- row.names(Lat_MLG_PRED_F_db)
    
    Activ_G_PRED_F_db$PRED_var <- "Activity"
    Activ_G_PRED_F_db$Predictor_var <- row.names(Activ_G_PRED_F_db)
    Activ_MLG_PRED_F_db$PRED_var <- "Activity"
    Activ_MLG_PRED_F_db$Predictor_var <- row.names(Activ_MLG_PRED_F_db)
    
    RBT_G_PRED_F_db$PRED_var <- "Relative_body_turns"
    RBT_G_PRED_F_db$Predictor_var <- row.names(RBT_G_PRED_F_db)
    RBT_MLG_PRED_db$PRED_var <- "Relative_body_turns"
    RBT_MLG_PRED_db$Predictor_var <- row.names(RBT_MLG_PRED_db)
    
    # Merge all results from anova repeated measurement analyses
    
    PRED_i_p.0 <- rbind(Lat_G_PRED_F_db,
                      Lat_MLG_PRED_F_db,
                      Activ_G_PRED_F_db,
                      Activ_MLG_PRED_F_db,
                      RBT_G_PRED_F_db,
                      RBT_MLG_PRED_db)
    names(PRED_i_p.0)
    PRED_i_p <- filter(PRED_i_p.0,
                       Predictor_var %in% c("INF_BIN_f:STIMULUS_f",
                                            "TREAT:STIMULUS_f",
                                            "MULT_INF_f:STIMULUS_f"))
    PRED_i_p <- PRED_i_p[,c("Pr..F.",
                            "PRED_var", "Predictor_var")]
    
    # Merge all results from glm
    
    CROWN_all_p <- rbind(CROWN_INF_BIN_p,
                         CROWN_MLGINF_int_p)
    
    
    CROWN_all_p <- CROWN_all_p %>% 
      dplyr::select(p, PRED_var, Predictor_var)%>%
      dplyr::rename(Pr..F.= p) %>%
      dplyr::filter(!Predictor_var=="joint")
    
    #  # Merge both database together  and FDR correction
    
    FDR_PREDi_allp <- rbind(CROWN_all_p,
                        PRED_i_p)
    
    FDR_PREDi_allp <- arrange(FDR_PREDi_allp, Pr..F.)
    
    FDR_PREDi_allp$FDR_p_adjust <-p.adjust(FDR_PREDi_allp$Pr..F.,
                                             method="fdr")
 